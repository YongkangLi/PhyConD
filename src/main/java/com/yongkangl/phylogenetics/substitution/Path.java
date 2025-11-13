package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.model.SMCParticle;
import com.yongkangl.phylogenetics.util.Mutation;
import com.yongkangl.phylogenetics.util.Taxon;
import jeigen.DenseMatrix;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;

public class Path extends SMCParticle {
    private final Taxon initial;
    private final Taxon terminal;
    private final double T;
    private final int blockSize;
    private final SiteModel target;
//    private final ComplexBlockwiseSubstitutionModel complexBlockwiseSubstitutionModel;
    private BlockwiseSubstitutionModel blockwiseSubstitutionModel;
    private final List<Mutation> mutations;
    private final int mutationSteps;

    public Path(int nSteps, Taxon initial, Taxon terminal, double T, int blockSize, BlockwiseSubstitutionModel blockwiseSubstitutionModel,  SiteModel target, int mutationSteps) {
        super(nSteps);
        this.initial = initial;
        this.terminal = terminal;
        this.T = T;
        this.blockSize = blockSize;

        this.blockwiseSubstitutionModel = new BlockwiseSubstitutionModel(blockwiseSubstitutionModel);

        this.target = target;
        mutations = sample(false);
        this.mutationSteps = mutationSteps;
    }

    public void copyFrom(Path from) {
        super.copyFrom(from);
        mutations.clear();
        mutations.addAll(from.mutations);
    }

    public double logLikelihood() {
        LogLikelihoodDelegate delegate = new LogLikelihoodDelegateImpl(initial, target);
        return delegate.logLikelihoodFullPath(mutations, T);
    }

    public double logMetropolisHastingsRatio(List<Mutation> neighborsOfNeighbors, List<Mutation> neighbors, List<Mutation> newMutations, List<Mutation> originalMutations, int targetSite) {
        LogLikelihoodDelegateImpl targetImpl = new LogLikelihoodDelegateImpl(initial, target);
//        LogLikelihoodDelegateImpl proposalImpl = new LogLikelihoodDelegateImpl(initial, proposal);
        double temperature = getCurrentTemperature();
        return targetImpl.logGeometricallyTemperedLikelihood(blockwiseSubstitutionModel, Mutation.mergeMutations(neighborsOfNeighbors, newMutations), T, targetSite, temperature)
//             - proposalImpl.logLikelihoodIndependentSingleSite(targetSite, newMutations, T)
             - targetImpl.logLikelihoodDependentSingleSite(targetSite, Mutation.mergeMutations(neighbors, newMutations), T)
             - targetImpl.logGeometricallyTemperedLikelihood(blockwiseSubstitutionModel, Mutation.mergeMutations(neighborsOfNeighbors, originalMutations), T, targetSite, temperature)
//             + proposalImpl.logLikelihoodIndependentSingleSite(targetSite, originalMutations, T);
             + targetImpl.logLikelihoodDependentSingleSite(targetSite, Mutation.mergeMutations(neighbors, originalMutations), T);
//        return 0.0;
    }

    private List<Mutation> sample(boolean hobolth) {
        return sample(hobolth, target.getContextLength());
    }

    private List<Mutation> sample(boolean hobolth, int contextLength) {
        int length = initial.len();
        List<Mutation>[] paths = new List[length - 2 * contextLength];

        int[] boundaries = blockwiseSubstitutionModel.getBoundaries();
        for (int i = 0; i < blockwiseSubstitutionModel.getEffectiveLength(); i++) {
            int start = boundaries[i];
            int end = boundaries[i + 1];
            BlockPathSampler blockPath = new BlockPathSampler(start, initial.getEncodedBlock(start, end), terminal.getEncodedBlock(start, end), T, blockwiseSubstitutionModel.getBlockwiseSiteModel(i), hobolth);
            List<Mutation>[] sitePaths = blockPath.getSitePaths();
            System.arraycopy(sitePaths, 0, paths, start - contextLength, end - start);
        }
        return Mutation.mergePaths(paths);
    }

    private double logLikelihood(double temperature) {
//        LogLikelihoodDelegate delegate = new LogLikelihoodDelegateImpl(initial, proposal);
//        double ll = delegate.logLikelihoodFullPath(mutations, T, target.getContextLength());
        return temperature * logLikelihood() + (1 - temperature) * blockwiseSubstitutionModel.logLikelihoodFullPath(mutations, T);
    }

    public double calculateStepLogWeight(int step) {
        return 0.0;
    }

    public void mutate() {
        determinsticScan();
    }

    public void determinsticScan() {
        int contextLength = target.getContextLength();
        if (contextLength > initial.len() - contextLength - 1) {
            throw new IllegalArgumentException("Sequence length is too short");
        }
        for (int i = 0; i < mutationSteps; i++) {
            for (int site = contextLength; site < initial.len() - contextLength; site++) {
                mutate(site);
            }
        }
    }

    public void mutate(int site) {
        SubstitutionModel substitutionModel = new SubstitutionModel(target, initial);
        int contextLength = target.getContextLength();
        ArrayList<Mutation> originalMutations = new ArrayList<>();
        ArrayList<Mutation> neighbors = new ArrayList<>();
        ArrayList<Mutation> neighborsOfNeighbors = new ArrayList<>();
        int endState = initial.get(site);
        for (int i = 0; i < mutations.size(); i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            if (currentSite == site) {
                endState = mutation.getTo();
                originalMutations.add(mutation);
                mutations.remove(i);
                i--;
            } else if (site - contextLength <= currentSite && currentSite <= site + contextLength) {
                neighbors.add(mutation);
                neighborsOfNeighbors.add(mutation);
            } else if (site - 2 * contextLength <= currentSite && currentSite <= site + 2 * contextLength) {
                neighborsOfNeighbors.add(mutation);
            }
        }
        double t = 0.0;
        substitutionModel.reset();

//        SubstitutionModel proposalModel = new SubstitutionModel(proposal, initial);
//        List<Mutation> newMutations = proposalModel.samplePath(site, endState, 0, T);

        List<Mutation> newMutations = new ArrayList<>();
        DenseMatrix[] matrices = substitutionModel.backward(site, neighbors, T);
        for (int i = 0; i <= neighbors.size(); i++) {
            double endTime;
            int nextEndState;
            if (i < neighbors.size()) {
                endTime = neighbors.get(i).getTime();
                double[] distribution = substitutionModel.intermediateStateDistribution(site, endTime - t);
                List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
                for (int state = 0; state < 4; state++) {
                    possibleStates.add(new Pair<>(state, distribution[state] * matrices[i + 1].get(endState, state)));
                }
                EnumeratedDistribution<Integer> stateDistribution = new EnumeratedDistribution<>(possibleStates);
                nextEndState = stateDistribution.sample();
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

        double u = new RandomDataGenerator().nextUniform(0, 1);
        double metropolisHastingsRatio = Math.exp(logMetropolisHastingsRatio(neighborsOfNeighbors, neighbors, newMutations, originalMutations, site));
        if (u < metropolisHastingsRatio) {
            Mutation.addAllMutations(mutations, newMutations);
        } else {
            Mutation.addAllMutations(mutations, originalMutations);
        }
    }

    public String toString() {
        return "(" + Mutation.toString(mutations) + ", " + getStepLogWeight(false) + ")";
    }

    public int numberOfMutations() {
        return mutations.size();
    }
}
