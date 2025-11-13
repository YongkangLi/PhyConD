package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.model.MetropolisHastingsParticle;
import com.yongkangl.phylogenetics.util.Mutation;
import com.yongkangl.phylogenetics.util.Taxon;
import com.yongkangl.phylogenetics.util.Utils;
import jeigen.DenseMatrix;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;

public class MetropolisHastingsPath extends MetropolisHastingsParticle {
    private final Taxon initial;
    private final Taxon terminal;
    private final double T;
    private SiteModel proposal;
    private SiteModel target;
    private final List<Mutation> mutations;

    // These fields are used to store intermediate results
    // to compute the Metropolis-Hastings ratio
    private int targetSite;
    private double jointProbability;
    private List<Mutation> neighbors;
    private List<Mutation> neighborsOfNeighbors;
    private List<Mutation> mutationsToEvaluate;
    private List<Mutation> mutationsProposed;
    private DenseMatrix[] forwardMatrices;
    private DenseMatrix[] backwardMatrices;
    private boolean dirty = true;
    private double logPiOverQ;


    // These fields are used to diagnose the convergence
    private List<Double> firstJumpTimes;
    private List<Integer> numberOfMutations;
    private List<Double> logPis;
    private List<Integer> orders;

    public MetropolisHastingsPath(Taxon initial, Taxon terminal, double T, IndependentSiteModel proposal, SiteModel target) {
        this.initial = new Taxon(initial);
        this.terminal = new Taxon(terminal);
        this.T = T;
        this.proposal = proposal;
        this.target = target;

        int contextLength = target.getContextLength();
        if (proposal.getContextLength() < contextLength) {
            mutations = initialize(initial, terminal, T, false, contextLength);
        } else {
            mutations = initialize(initial, terminal, T, false);
        }
        this.targetSite = target.getContextLength() - 1;

        this.firstJumpTimes = new ArrayList<>();
        this.numberOfMutations = new ArrayList<>();
        this.logPis = new ArrayList<>();
        this.orders = new ArrayList<>();
    }

    public MetropolisHastingsPath(Taxon initial, Taxon terminal, double T, SiteModel proposal, SiteModel target, List<Mutation> mutations, int targetSite, List<Mutation> mutationsToEvaluate, List<Mutation> neighbors, List<Mutation> neighborsOfNeighbors) {
        this.initial = new Taxon(initial);
        this.terminal = new Taxon(terminal);
        this.T = T;
        this.proposal = proposal;
        this.target = target;
        this.mutations = mutations;
        this.targetSite = targetSite;
        this.mutationsToEvaluate = mutationsToEvaluate;
        this.neighbors = neighbors;
        this.neighborsOfNeighbors = neighborsOfNeighbors;
        this.dirty = true;
    }

    public void setTarget(SiteModel target) {
        this.target = target;
    }

    public void step(int site, int steps) {
        for (int step = 0; step < steps; step++) {
            MetropolisHastingsParticle proposal = hobolthPropose(site);
            double numerator = proposal.logPiOverQ();
            double denominator = logPiOverQ();
            double metropolisHastingsRatio = Math.exp(numerator - denominator);
            double u = new RandomDataGenerator().nextUniform(0.0, 1.0);
            if (u < metropolisHastingsRatio) {
                copyFrom(proposal);
            }
            log();
        }
    }

    public void copyFrom(MetropolisHastingsParticle from) {
        if (from instanceof MetropolisHastingsPath) {
            MetropolisHastingsPath path = (MetropolisHastingsPath) from;
            initial.copySequence(path.initial);
            terminal.copySequence(path.terminal);
            mutations.clear();
            mutations.addAll(path.mutations);
            logPiOverQ = path.logPiOverQ;
            dirty = false;
            neighbors = null;
            neighborsOfNeighbors = null;
            mutationsToEvaluate = null;
            mutationsProposed = null;
            forwardMatrices = null;
            backwardMatrices = null;
        }
    }

    public MetropolisHastingsParticle propose() {
        return hobolthPropose();
    }

    public MetropolisHastingsParticle independentPropose() {
        int contextLength = target.getContextLength();
        targetSite = (targetSite - contextLength + 1) % (initial.len() - 2 * target.getContextLength()) + contextLength;

        mutationsToEvaluate = new ArrayList<>();
        neighbors = new ArrayList<>();
        neighborsOfNeighbors = new ArrayList<>();
        for (int i = 0; i < mutations.size(); i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            if (currentSite == targetSite) {
                mutationsToEvaluate.add(mutation);
                mutations.remove(i);
                i--;
            } else if (targetSite - contextLength <= currentSite && currentSite <= targetSite + contextLength) {
                neighbors.add(mutation);
                neighborsOfNeighbors.add(mutation);
            } else if (targetSite - 2 * contextLength <= currentSite && currentSite <= targetSite + 2 * contextLength) {
                neighborsOfNeighbors.add(mutation);
            }
        }

        dirty = true;

        SubstitutionModel proposalModel = new SubstitutionModel(proposal, initial);
        List<Mutation> newMutations = proposalModel.endpointedSimulation(targetSite, terminal.get(targetSite), 0, T);

        MetropolisHastingsParticle proposed = new MetropolisHastingsPath(initial, terminal, T, proposal, target, Mutation.mergeMutations(mutations, newMutations), targetSite, newMutations, neighbors, neighborsOfNeighbors);
        Mutation.addAllMutations(mutations, mutationsToEvaluate);
        return proposed;
    }

    public DenseMatrix getLastForwardMatrice() {
        return forwardMatrices[forwardMatrices.length - 1];
    }

    public DenseMatrix getFirstBackwardMatrice() {
        return backwardMatrices[0];
    }

    public DenseMatrix[] getBackwardMatrices() {
        return backwardMatrices;
    }

    public static int hobolthEndpointedOnChildren(List<DenseMatrix> childrenBackwards, List<Integer> childrenStates) {
        List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
        double[] probabilities = new double[4];
        for (int state = 0; state < 4; state++) {
            probabilities[state] = 1.0;
            for (int child = 0; child < childrenBackwards.size(); child++) {
                probabilities[state] *= childrenBackwards.get(child).get(state, childrenStates.get(child));
            }
            possibleStates.add(new Pair<>(state, probabilities[state]));
        }
        EnumeratedDistribution<Integer> stateDistribution = new EnumeratedDistribution<>(possibleStates);
        int endState = stateDistribution.sample();
        return endState;
    }

    public double hobolthEndpointedOnChildren(int site, DenseMatrix forward, List<DenseMatrix> childrenBackwards, List<Integer> childrenStates) {
        List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
        double[] probabilities = new double[4];
        for (int state = 0; state < 4; state++) {
            probabilities[state] = 1.0;
            for (int child = 0; child < childrenBackwards.size(); child++) {
                probabilities[state] *= childrenBackwards.get(child).get(state, childrenStates.get(child));
            }
            probabilities[state] *= forward.get(initial.get(site), state);
            possibleStates.add(new Pair<>(state, probabilities[state]));
        }
        EnumeratedDistribution<Integer> stateDistribution = new EnumeratedDistribution<>(possibleStates);
        int endState = stateDistribution.sample();
        setTerminal(site, endState);
        return hobolthEndpointed(site);
    }

    public double hobolthForward(int site) {
        int contextLength = target.getContextLength();
        targetSite = site;
        // update the target site to form a deterministic scan
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        mutationsToEvaluate = new ArrayList<>();
        neighbors = new ArrayList<>();
        neighborsOfNeighbors = new ArrayList<>();
        for (int i = 0; i < mutations.size(); i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            if (currentSite == site) {
                mutationsToEvaluate.add(mutation);
                mutations.remove(i);
                i--;
            } else if (site - contextLength <= currentSite && currentSite <= site + contextLength) {
                neighbors.add(mutation);
                neighborsOfNeighbors.add(mutation);
            } else if (site - 2 * contextLength <= currentSite && currentSite <= site + 2 * contextLength) {
                neighborsOfNeighbors.add(mutation);
            }
        }

        // To signal that the logPiOverQ needs to be re-computed!!!
        dirty = true;
        substitutionModel.reset();
        forwardMatrices = substitutionModel.forward(site, neighbors, T);
        return logHobolthPiOverQ(mutationsToEvaluate);
    }

    public double hobolthBackward(int site) {
        int contextLength = target.getContextLength();
        targetSite = site;
        // update the target site to form a deterministic scan
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        mutationsToEvaluate = new ArrayList<>();
        neighbors = new ArrayList<>();
        neighborsOfNeighbors = new ArrayList<>();
        for (int i = 0; i < mutations.size(); i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            if (currentSite == site) {
                mutationsToEvaluate.add(mutation);
                mutations.remove(i);
                i--;
            } else if (site - contextLength <= currentSite && currentSite <= site + contextLength) {
                neighbors.add(mutation);
                neighborsOfNeighbors.add(mutation);
            } else if (site - 2 * contextLength <= currentSite && currentSite <= site + 2 * contextLength) {
                neighborsOfNeighbors.add(mutation);
            }
        }

        // To signal that the logPiOverQ needs to be re-computed!!!
        dirty = true;
        substitutionModel.reset();
        backwardMatrices = substitutionModel.backward(site, neighbors, T);
        return logHobolthPiOverQ(mutationsToEvaluate);
    }

    public double hobolthEndpointed(int site) {
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        return hobolthEndpointed(site, substitutionModel.backward(site, neighbors, T));
    }

    public double hobolthEndpointed(int site, DenseMatrix[] matrices) {
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        int endState = terminal.get(site);
        jointProbability = 1.0;
        double t = 0.0;
        mutationsProposed = new ArrayList<>();
        for (int i = 0; i <= neighbors.size(); i++) {
            double endTime;
            int nextEndState;
            if (i < neighbors.size()) {
                endTime = neighbors.get(i).getTime();
                double[] distribution = substitutionModel.intermediateStateDistribution(site, endTime - t);
                List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
                double normalizingConstant = 0.0;
                double[] unnormalizedProbabilities = new double[4];
                for (int state = 0; state < 4; state++) {
                    unnormalizedProbabilities[state] = distribution[state] * matrices[i+1].get(state, endState);
                    possibleStates.add(new Pair<>(state, unnormalizedProbabilities[state]));
                    normalizingConstant += unnormalizedProbabilities[state];
                }
                EnumeratedDistribution<Integer> stateDistribution = new EnumeratedDistribution<>(possibleStates);
                nextEndState = stateDistribution.sample();
                jointProbability *= unnormalizedProbabilities[nextEndState] / normalizingConstant;
            } else {
                endTime = T;
                nextEndState = endState;
            }
            mutationsProposed.addAll(substitutionModel.endpointedSimulation(site, nextEndState, t, endTime));
            substitutionModel.mutate(site, nextEndState);
            if (i != neighbors.size()) {
                substitutionModel.mutate(neighbors.get(i));
            }
            t = endTime;
        }
        return logHobolthPiOverQ(mutationsProposed);
    }

    public void accept(boolean accept) {
        if (accept) {
            Mutation.addAllMutations(mutations, mutationsProposed);
        } else {
            Mutation.addAllMutations(mutations, mutationsToEvaluate);
        }
    }

    public MetropolisHastingsParticle hobolthPropose(int site) {
        int contextLength = target.getContextLength();
        targetSite = site;
        // update the target site to form a deterministic scan
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        mutationsToEvaluate = new ArrayList<>();
        neighbors = new ArrayList<>();
        neighborsOfNeighbors = new ArrayList<>();
        int endState = terminal.get(site);
        for (int i = 0; i < mutations.size(); i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            if (currentSite == site) {
                mutationsToEvaluate.add(mutation);
                mutations.remove(i);
                i--;
            } else if (site - contextLength <= currentSite && currentSite <= site + contextLength) {
                neighbors.add(mutation);
                neighborsOfNeighbors.add(mutation);
            } else if (site - 2 * contextLength <= currentSite && currentSite <= site + 2 * contextLength) {
                neighborsOfNeighbors.add(mutation);
            }
        }

        // To signal that the logPiOverQ needs to be re-computed!!!
        dirty = true;

        double t = 0.0;
        substitutionModel.reset();

        List<Mutation> newMutations = new ArrayList<>();
        DenseMatrix[] matrices = substitutionModel.backward(site, neighbors, T);
        jointProbability = 1.0;
        for (int i = 0; i <= neighbors.size(); i++) {
            double endTime;
            int nextEndState;
            if (i < neighbors.size()) {
                endTime = neighbors.get(i).getTime();
                double[] distribution = substitutionModel.intermediateStateDistribution(site, endTime - t);
                List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
                double normalizingConstant = 0.0;
                double unnormalizedProbabilities[] = new double[4];
                for (int state = 0; state < 4; state++) {
                    unnormalizedProbabilities[state] = distribution[state] * matrices[i+1].get(state, endState);
                    possibleStates.add(new Pair<>(state, unnormalizedProbabilities[state]));
                    normalizingConstant += unnormalizedProbabilities[state];
                }
                EnumeratedDistribution<Integer> stateDistribution = new EnumeratedDistribution<>(possibleStates);
                nextEndState = stateDistribution.sample();
                jointProbability *= unnormalizedProbabilities[nextEndState] / normalizingConstant;
            } else {
                endTime = T;
                nextEndState = endState;
            }
            newMutations.addAll(substitutionModel.endpointedSimulation(site, nextEndState, t, endTime));
            substitutionModel.mutate(site, nextEndState);
            if (i != neighbors.size()) {
                substitutionModel.mutate(neighbors.get(i));
            }
            t = endTime;
        }

        MetropolisHastingsParticle proposed = new MetropolisHastingsPath(initial, terminal, T, proposal, target, Mutation.mergeMutations(mutations, newMutations), targetSite, newMutations, neighbors, neighborsOfNeighbors);
        Mutation.addAllMutations(mutations, mutationsToEvaluate);
        return proposed;
    }

    public MetropolisHastingsParticle hobolthPropose() {
        int contextLength = target.getContextLength();
        // update the target site to form a deterministic scan
        targetSite = (targetSite - contextLength + 1) % (initial.len() - 2 * target.getContextLength()) + contextLength;
        return hobolthPropose(targetSite);
    }

    public double logPiOverQ() {
        return logHobolthPiOverQ();
    }

    public double logIndependentPiOverQ() {
        return logPi() - logIndependentQ();
    }

    public double logHobolthPiOverQ(List<Mutation> toEvaluate) {
        double ll = 0.0;
        List<Mutation> keyMutations = Mutation.mergeMutations(neighborsOfNeighbors, toEvaluate);
        int size = keyMutations.size();
        int len = initial.len();
        int contextLength = target.getContextLength();
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        double t = 0.0;
        for (int i = 0; i < size; i++) {
            Mutation mutation = keyMutations.get(i);
            int currentSite = mutation.getSite();
            double rateSum = 0.0;
            for (int site = targetSite - contextLength; site <= targetSite + contextLength; site++) {
                if (site != targetSite && contextLength <= site && site < len - contextLength) {
                    rateSum += substitutionModel.getMutability(site);
                }
            }
            ll -= (mutation.getTime() - t) * rateSum;
            double mutationRate = substitutionModel.mutate(currentSite, mutation.getTo());
            if (mutationRate <= 0) {
                System.out.println("Mutation rate is negative: " + mutationRate);
            }
            if (targetSite - contextLength <= currentSite && currentSite <= targetSite + contextLength && currentSite != targetSite) {
                ll += Math.log(mutationRate);
            }
            t = mutation.getTime();
        }
        double rateSum = 0.0;
        for (int site = targetSite - contextLength; site <= targetSite + contextLength; site++) {
            if (site != targetSite && contextLength <= site && site < len - contextLength) {
                rateSum += substitutionModel.getMutability(site);
            }
        }
        ll -= (T - t) * rateSum;

        return ll;
    }

    public double logHobolthPiOverQ() {
        if (dirty) {
            logPiOverQ = logHobolthPiOverQ(mutationsToEvaluate);
            dirty = false;
            neighbors = null;
            neighborsOfNeighbors = null;
            mutationsToEvaluate = null;
        }
        return logPiOverQ;
    }

    public double logPi() {
        double ll = 0.0;
        List<Mutation> keyMutations = Mutation.mergeMutations(neighborsOfNeighbors, mutationsToEvaluate);
        int size = keyMutations.size();
        int len = initial.len();
        int contextLength = target.getContextLength();
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        double t = 0.0;
        for (int i = 0; i < size; i++) {
            Mutation mutation = keyMutations.get(i);
            int currentSite = mutation.getSite();
            double rateSum = 0.0;
            for (int site = targetSite - contextLength; site <= targetSite + contextLength; site++) {
                if (contextLength <= site && site < len - contextLength) {
                    rateSum += substitutionModel.getMutability(site);
                }
            }
            ll -= (mutation.getTime() - t) * rateSum;
            double mutationRate = substitutionModel.mutate(currentSite, mutation.getTo());
            if (targetSite - contextLength <= currentSite && currentSite <= targetSite + contextLength) {
                ll += Math.log(mutationRate);
            }
            t = mutation.getTime();
        }
        double rateSum = 0.0;
        for (int site = targetSite - contextLength; site <= targetSite + contextLength; site++) {
            if (contextLength <= site && site < len - contextLength) {
                rateSum += substitutionModel.getMutability(site);
            }
        }
        ll -= (T - t) * rateSum;
        return ll;
    }

    public double logFullPi() {
        double ll = 0.0;
        int size = mutations.size();
        int len = initial.len();
        int contextLength = target.getContextLength();
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        double t = 0.0;
        for (int i = 0; i < size; i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            double rateSum = 0.0;
            for (int site = contextLength; site < len - contextLength; site++) {
                rateSum += substitutionModel.getMutability(site);
            }
            ll -= (mutation.getTime() - t) * rateSum;
            double mutationRate = substitutionModel.mutate(currentSite, mutation.getTo());
            ll += Math.log(mutationRate);
            t = mutation.getTime();
        }
        double rateSum = 0.0;
        for (int site = contextLength; site < len - contextLength; site++) {
            rateSum += substitutionModel.getMutability(site);
        }
        ll -= (T - t) * rateSum;
        return ll;
    }

    public double logIndependentQ() {
        LogLikelihoodDelegateImpl logLikelihoodDelegate = new LogLikelihoodDelegateImpl(initial, proposal);
        return logLikelihoodDelegate.logLikelihoodIndependentSingleSite(targetSite, mutationsToEvaluate, T);
    }

    public double logHolbothQ() {
        double ll = 0.0;
        List<Mutation> mutations = Mutation.mergeMutations(neighbors, mutationsToEvaluate);
        int size = mutations.size();
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        double t = 0.0;
        for (int i = 0; i < size; i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            ll -= (mutation.getTime() - t) * substitutionModel.getMutability(targetSite);
            double mutationRate = substitutionModel.mutate(currentSite, mutation.getTo());
            if (currentSite == targetSite) {
                ll += Math.log(mutationRate);
            }
            t = mutation.getTime();
        }
        ll -= (T - t) * substitutionModel.getMutability(targetSite);
        return ll;
    }

    public double logImportanceWeight (SiteModel denominator, SiteModel numerator) {
        int len = initial.len();
        int contextLength = numerator.getContextLength();
        double ll = 0.0;
        int size = mutations.size();
        SubstitutionModel denominatorModel = new SubstitutionModel(denominator, initial);
        SubstitutionModel numeratorModel = new SubstitutionModel(numerator, initial);
        double t = 0.0;
        for (int i = 0; i < size; i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            double rateSum = 0.0;
            for (int site = contextLength; site < len - contextLength; site++) {
                    rateSum += numeratorModel.getMutability(site) - denominatorModel.getMutability(site);
            }
            ll -= (mutation.getTime() - t) * rateSum;
            ll += Math.log(numeratorModel.mutate(currentSite, mutation.getTo())) - Math.log(denominatorModel.mutate(currentSite, mutation.getTo()));
            t = mutation.getTime();
        }
        double rateSum = 0.0;
        for (int site = contextLength; site < initial.len() - contextLength; site++) {
            rateSum += numeratorModel.getMutability(site) - denominatorModel.getMutability(site);
        }
        ll -= (T - t) * rateSum;
        return ll;
    }

    private List<Mutation> initialize(Taxon initial, Taxon terminal, double T, boolean hobolth) {
        return initialize(initial, terminal, T, hobolth, proposal.getContextLength());
    }

    private List<Mutation> initialize(Taxon initial, Taxon terminal, double T, boolean hobolth, int contextLength) {
        int len = initial.len();
        List<Mutation>[] paths = new List[len - 2 * contextLength];

        for (int i = contextLength; i < initial.len() - contextLength; i++) {
            double[][] kernel = BlockwiseSiteModel.getKernel(1, new int[0], new int[0], proposal);
            BlockPathSampler path = new BlockPathSampler(i, initial.get(i), terminal.get(i), T, new BlockwiseSiteModel(kernel), hobolth);
            paths[i - contextLength] = path.getMutations();
        }
        return Mutation.mergePaths(paths);
    }

    public String toString() {
        StringBuilder time = new StringBuilder();
        StringBuilder site = new StringBuilder();
        StringBuilder to = new StringBuilder();
        time.append("[");
        site.append("[");
        to.append("[");
        if (!mutations.isEmpty()) {
            Mutation mutation = mutations.get(0);
            time.append(mutation.getTime());
            site.append(mutation.getSite());
            to.append(mutation.getTo());
        }
        for (int i = 1; i < mutations.size(); i++){
            Mutation mutation = mutations.get(i);
            time.append(", ").append(mutation.getTime());
            site.append(", ").append(mutation.getSite());
            to.append(", ").append(mutation.getTo());
        }
        time.append("]");
        site.append("]");
        to.append("]");
        return "(" + time + ", " + site + ", " + to + ")";
    }

    public Mutation getMutation(int i) {
        return mutations.get(i);
    }

    public int getNumberOfMutations() {
        return mutations.size();
    }

    public void setInitial(int site, int state) {
        initial.set(site, state);
    }

    public void setTerminal(int site, int state) {
        terminal.set(site, state);
    }

    public int getInitialState(int site) {
        return initial.get(site);
    }

    public int getTerminalState(int site) {
        return terminal.get(site);
    }

    public double getFirstJumpTime(int step) {
        return firstJumpTimes.get(step);
    }

    public int getNumberOfMutations(int step) {
        return numberOfMutations.get(step);
    }

    public double getLogPi(int step) {
        return logPis.get(step);
    }

    public double getPi(int step) {
        return Math.exp(logPis.get(step));
    }

    public int getOrder(int step) {
        return orders.get(step);
    }

    public double[] calculateAutoCovariance(int nSteps, int endLag, int burnIn) {
        double[] autoCovariances = new double[endLag];
        double mean = 0.0;
        for (int step = burnIn; step < nSteps; step++) {
            mean += getFirstJumpTime(step);
        }
        mean /= (nSteps - burnIn);
        for (int lag = 0; lag < endLag; lag++) {
            autoCovariances[lag] = 0.0;
            for (int step = burnIn; step < nSteps - lag; step++) {
                autoCovariances[lag] = (getFirstJumpTime(step) - mean) * (getFirstJumpTime(step + lag) - mean);
            }
            autoCovariances[lag] /= (nSteps - lag - burnIn);
        }
        return autoCovariances;
    }

    public double calculateAutoCorrelation(int nSteps, int lag, int burnIn) {
        int n = nSteps - lag - burnIn;
        double[] array1 = new double[n];
        double[] array2 = new double[n];
        for (int i = 0; i < n; i++) {
            array1[i] = getOrder(burnIn + i);
            array2[i] = getOrder(burnIn + i + lag);
        }
        double mean1 = Utils.sum(array1) / n;
        double mean2 = Utils.sum(array2) / n;
        double autoCovariance = 0.0;
        for (int i = 0; i < n; i++) {
            double dev1 = array1[i] - mean1;
            double dev2 = array2[i] - mean2;
            autoCovariance += dev1 * dev2;
            array1[i] = dev1 * dev1;
            array2[i] = dev2 * dev2;
        }
        autoCovariance /= n;
        double std1 = Math.sqrt(Utils.sum(array1) / n);
        double std2 = Math.sqrt(Utils.sum(array2) / n);
        return autoCovariance / (std1 * std2);
    }

    public void log() {
    }

    public static class K80CpGSufficientStatistics {
        private double n_ts;
        private double n_tv;
        private double n_CpG;
        private double count;

        public K80CpGSufficientStatistics(double n_ts, double n_tv, double n_CpG, double count) {
            this.n_ts = n_ts;
            this.n_tv = n_tv;
            this.n_CpG = n_CpG;
            this.count = count;
        }

        public double getNts() {
            return n_ts;
        }

        public double getNtv() {
            return n_tv;
        }

        public double getNCpG() {
            return n_CpG;
        }

        public double getCount() {
            return count;
        }
    }

    public K80CpGSufficientStatistics getK80CpGSufficnetStatistics() {
        int len = target.getContextLength();

        int n_ts = 0;
        int n_tv = 0;
        int n_CpG = 0;
        double weighted_n_CpG = 0.0;
        int count = 0;

        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);

        for (int site = len; site < initial.getSequence().length - 1; site++) {
            if (initial.get(site -1) == Nucleotides.C_STATE && initial.get(site) == Nucleotides.G_STATE) {
                n_CpG++;
            }
        }

        double t = 0.0;

        for (Mutation mutation : mutations) {
            int site = mutation.getSite();
            int fromState = substitutionModel.getCurrent().get(site);
            int toState = mutation.getTo();
            double time = mutation.getTime();

            weighted_n_CpG += n_CpG * (time - t);

            // Check if it is a transition or transversion
            if (isTransition(fromState, toState)) {
                n_ts++;
            } else if (isTransversion(fromState, toState)) {
                n_tv++;
            }

            substitutionModel.mutate(site, toState);

            int leftNeighbor = substitutionModel.getCurrent().get(site - 1);
            int rightNeighbor = substitutionModel.getCurrent().get(site + 1);
            int oldLocalNcpg = ((leftNeighbor == Nucleotides.C_STATE && fromState == Nucleotides.G_STATE) || (fromState == Nucleotides.C_STATE && rightNeighbor == Nucleotides.G_STATE)) ? 1 : 0;
            int newLocalNcpg = ((leftNeighbor == Nucleotides.C_STATE && toState == Nucleotides.G_STATE) || (toState == Nucleotides.C_STATE && rightNeighbor == Nucleotides.G_STATE)) ? 1 : 0;
            n_CpG = n_CpG - oldLocalNcpg + newLocalNcpg;

            count += oldLocalNcpg;

            t = time;
        }
        weighted_n_CpG += n_CpG * (T - t);

        return new K80CpGSufficientStatistics(n_ts, n_tv, weighted_n_CpG, count);
    }

    private static boolean isTransition(int fromState, int toState) {
        return (fromState == Nucleotides.A_STATE && toState == Nucleotides.G_STATE) ||
                (fromState == Nucleotides.G_STATE && toState == Nucleotides.A_STATE) ||
                (fromState == Nucleotides.C_STATE && toState == Nucleotides.UT_STATE) ||
                (fromState == Nucleotides.UT_STATE && toState == Nucleotides.C_STATE);
    }

    private static boolean isTransversion(int fromState, int toState) {
        return (fromState == Nucleotides.A_STATE && (toState == Nucleotides.C_STATE || toState == Nucleotides.UT_STATE)) ||
                (fromState == Nucleotides.G_STATE && (toState == Nucleotides.C_STATE || toState == Nucleotides.UT_STATE)) ||
                (fromState == Nucleotides.C_STATE && (toState == Nucleotides.A_STATE || toState == Nucleotides.G_STATE)) ||
                (fromState == Nucleotides.UT_STATE && (toState == Nucleotides.A_STATE || toState == Nucleotides.G_STATE));
    }
}
