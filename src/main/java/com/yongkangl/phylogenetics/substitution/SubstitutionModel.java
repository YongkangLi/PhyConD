package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Mutation;
import com.yongkangl.phylogenetics.util.Taxon;
import jeigen.DenseMatrix;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SubstitutionModel {
    private final SiteModel siteModel;
    private final Taxon initial;
    private Taxon current;
    private final double[] mutabilities;
    private final List<Integer> masks;

    public SubstitutionModel(SubstitutionModel copyFrom) {
        initial = new Taxon(copyFrom.getCurrent());
        siteModel = copyFrom.getSiteModel();
        current = new Taxon(initial.len());
        mutabilities = new double[initial.len()];
        masks = new ArrayList<>();
        reset();
    }

    public SubstitutionModel(SiteModel siteModel, Taxon initial) {
        this.siteModel = siteModel;
        this.initial = new Taxon(initial.len());
        this.initial.copySequence(initial);

        current = new Taxon(initial.len());
        mutabilities = new double[initial.len()];

        masks = new ArrayList<>();

        reset();
    }

    public void addMask(int mask) {
        masks.add(mask);
        mutabilities[mask] = 0.0;
    }

    public Taxon getInitialSequence() {
        return initial;
    }

    public SiteModel getSiteModel() {
        return siteModel;
    }

    public void reset() {
        current.copySequence(initial);

        int contextLength = siteModel.getContextLength();

        for (int i = 0; i < contextLength; i++) {
            mutabilities[i] = 0.0;
        }

        for (int i = contextLength; i < initial.len() - contextLength; i++) {
            mutabilities[i] = getMutability(i);
        }

        for (int i = initial.len() - contextLength; i < initial.len(); i++) {
            mutabilities[i] = 0.0;
        }

        for (int mask: masks) {
            mutabilities[mask] = 0.0;
        }
    }

    private int getContext(int site) {
        int contextLength = siteModel.getContextLength();
        assert site >= contextLength && site < current.len() - contextLength;
        int encoded =  0;
        for (int i = 0; i < 2 * contextLength + 1; i++) {
            encoded <<= 2;
            encoded |= current.get(site - contextLength + i);
        }
        return encoded;
    }

    public double getMutability(int site) {
        int encoded = getContext(site);
        return siteModel.getMutability(encoded);
    }

//    public double geometricPathMutability(BlockwiseSubstitutionModel blockwiseSubstitutionModel, int site, double temperature) {
//        double mutability = 0.0;
//        int current = getCurrent().get(site);
//        for (int to = 0; to < 4; to++) {
//            if (to != current) {
//                mutability += Math.pow(blockwiseSubstitutionModel.compute(current, to), 1 - temperature) * Math.pow(getRate(site, to), temperature);
//            }
//        }
//        return mutability;
//    }

    public double getRate(int site, int to) {
        int encoded = getContext(site);
        encoded <<= 2;
        encoded |= to;
        return siteModel.getRate(encoded);
    }

    public double getRate(int site, int to, double temperature) {
        return Math.pow(getRate(site, to), temperature);
    }

    public void set(int site, int to) {
        current.set(site, to);
        int len = siteModel.getContextLength();
        for (int i = site < len ? len : site - len; i <= site + len && i < current.len() - len; i++) {
            updateMutability(i);
        }
    }

    public double mutate(int site, int to) {
        double ret = getRate(site, to);
        set(site, to);
        return ret;
    }

    public double mutate(int site, int to, double temperature) {
        return temperature * Math.log(mutate(site, to));
    }

    public double mutate(Mutation mutation) {
        int site = mutation.getSite();
        int to = mutation.getTo();
        return mutate(site, to);
    }

    public void updateMutability(int site) {
        int len = siteModel.getContextLength();
        if (site >= len && site < current.len() - len) {
            mutabilities[site] = getMutability(site);
        }
        for (int mask: masks) {
            mutabilities[mask] = 0.0;
        }
    }

    public double getMutabilitiesSum(int exclusion) {
        double sum = 0.0;
        for (int i = exclusion; i < current.len() - exclusion; i++) {
            sum += mutabilities[i];
        }
        return sum;
    }
    public double getMutabilitiesSum() {
        return Arrays.stream(mutabilities).sum();
    }

    public int getContextLength() {
        return siteModel.getContextLength();
    }

    private DenseMatrix getTransitionMatrix(int context, double T) {
        context = context & (~ (0x3 << 2 * siteModel.getContextLength()));
        double[][] Q = new double[4][4];
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                int encoded = context | (i << (2 * siteModel.getContextLength()));
                Q[i][j] = siteModel.getRate((encoded << 2) | j);
            }
        }
        IndependentLogLikelihood logLikelihood = new IndependentLogLikelihood(Q);
        return logLikelihood.setTime(T);
    }

    public DenseMatrix[] forward(int site, List<Mutation> mutations, double T) {
        int[] contexts = new int[mutations.size() + 1];
        contexts[0] = getContext(site) & (~(0x3 << 2 * siteModel.getContextLength()));
        for (int i = 1; i < mutations.size() + 1; i++) {
            Mutation mutation = mutations.get(i - 1);
            int s = site + siteModel.getContextLength() - mutation.getSite();
            contexts[i] = (contexts[i - 1] & (~(0x3 << 2 * s))) | (mutation.getTo() << (2 * s));
        }
        DenseMatrix[] matrices = new DenseMatrix[mutations.size() + 1];
        if (mutations.size() > 0) {
            matrices[0] = getTransitionMatrix(contexts[0], mutations.get(0).getTime());
            for (int i = 1; i < mutations.size(); i++) {
                matrices[i] = matrices[i - 1].mmul(getTransitionMatrix(contexts[i], mutations.get(i).getTime() - mutations.get(i - 1).getTime()));
            }
            matrices[mutations.size()] = matrices[mutations.size() - 1].mmul(getTransitionMatrix(contexts[mutations.size()], T - mutations.get(mutations.size() - 1).getTime()));
        } else {
            matrices[0] = getTransitionMatrix(contexts[0], T);
        }
        return matrices;
    }

    // TODO: Examine whether the non-reversibility of the model is a problem
    public DenseMatrix[] backward(int site, List<Mutation> mutations, double T) {
        int[] contexts = new int[mutations.size() + 1];
        contexts[0] = getContext(site) & (~ (0x3 << 2 * siteModel.getContextLength()));
        for (int i = 1; i < mutations.size() + 1; i++) {
            Mutation mutation = mutations.get(i - 1);
            int s = site + siteModel.getContextLength() - mutation.getSite();
            contexts[i] = (contexts[i - 1] & (~ (0x3 << 2 * s))) | (mutation.getTo() << (2 * s));
        }
        DenseMatrix[] matrices = new DenseMatrix[mutations.size() + 1];
        if (!mutations.isEmpty()) {
            matrices[mutations.size()] = getTransitionMatrix(contexts[mutations.size()], T - mutations.get(mutations.size() - 1).getTime());
            for (int i = mutations.size() - 1; i > 0; i--) {
//            matrices[i] = matrices[i + 1].mmul(getTransitionMatrix(contexts[i + 1], mutations.get(i + 1).getTime() - mutations.get(i).getTime()));
                matrices[i] = getTransitionMatrix(contexts[i], mutations.get(i).getTime() - mutations.get(i - 1).getTime()).mmul(matrices[i + 1]);
            }
            matrices[0] = getTransitionMatrix(contexts[0], mutations.get(0).getTime()).mmul(matrices[1]);
        } else {
            matrices[0] = getTransitionMatrix(contexts[0], T);
        }
        return matrices;
    }

    public double[] intermediateStateDistribution(int site, double T) {
        int encoded = getContext(site) & (~ (0x3 << 2 * siteModel.getContextLength()));
        DenseMatrix matrix = getTransitionMatrix(encoded, T);
        double[] distribution = new double[4];
        for (int i = 0; i < 4; i++) {
            distribution[i] = matrix.get(current.get(site), i);
        }
        return distribution;
    }

    public int sampleNextState(int site) {
        int state = current.get(site);
        List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
        for (int i = 0; i < 4; i++) {
            if (i != state) {
                possibleStates.add(new Pair<>(i, getRate(site, i)));
            }
        }
        return new EnumeratedDistribution<>(possibleStates).sample();
    }

    public Mutation sampleNextMutation(int site, double t) {
        double dt = new ExponentialDistribution(1.0 / getMutability(site)).sample();
        int to = sampleNextState(site);
        return new Mutation(t + dt, site, to);
    }

    public Mutation sampleNextMutationJointly(double t) {
        int len = siteModel.getContextLength();
        List<Pair<Integer, Double>> possibleSites = new ArrayList<>();
        double mutabilitySum = 0;
        for (int i = len; i <current.len() - len; i++) {
            double mutability = mutabilities[i];
            possibleSites.add(new Pair<>(i,mutability));
            mutabilitySum += mutability;
        }
        int site = new EnumeratedDistribution<>(possibleSites).sample();
        double dt = new ExponentialDistribution(1.0 / mutabilitySum).sample();
        int to = sampleNextState(site);
        return new Mutation(t + dt, site, to);
    }

    public List<Mutation> forwardSimulation(int site, double startTime, double endTime) {
        double t = startTime;
        ArrayList<Mutation> mutations = new ArrayList<>();
        while (t < endTime) {
            Mutation mutation = sampleNextMutation(site, t);
            if (mutation.getTime() < endTime) {
                mutate(mutation.getSite(), mutation.getTo());
                t = mutation.getTime();
                mutations.add(mutation);
            } else {
                break;
            }
        }
        return mutations;
    }

    public List<Mutation> endpointedSimulation(int site, int nextState, double startTime, double endTime) {
        int initialState = current.get(site);
        ArrayList<Mutation> mutations;
        do {
            double t = startTime;
            mutate(site, initialState);
            mutations = new ArrayList<>();
            while (t < endTime) {
                Mutation mutation = sampleNextMutation(site, t);
                if (mutation.getTime() < endTime) {
                    mutate(mutation.getSite(), mutation.getTo());
                    t = mutation.getTime();
                    mutations.add(mutation);
                } else {
                    break;
                }
            }
        } while (current.get(site) != nextState);
        return mutations;
    }

    private void simulateIndependently(int site, double T) {
        double t = 0.0;
        while (t < T) {
            Mutation mutation = sampleNextMutation(site, t);
            if (mutation.getTime() < T) {
                mutate(mutation.getSite(), mutation.getTo());
                t = mutation.getTime();
            } else {
                break;
            }
        }
    }

    public Taxon simulateIndependently(double T) {
        int l = getContextLength();
        for (int site = l; site < initial.len() - l; site++) {
            simulateIndependently(site, T);
        }
        Taxon terminal = current;
        current = new Taxon(initial.len());
        current.copySequence(initial);
        return terminal;
    }

    public Taxon simulateDependently(double T) {
        double t = 0.0;
        while (true) {
            Mutation mutation = sampleNextMutationJointly(t);
            t = mutation.getTime();
            if (t < T) {
                mutate(mutation);
            }
            else {
                break;
            }
        }
        Taxon terminal = current;
        current = new Taxon(initial.len());
        current.copySequence(initial);
        return terminal;
    }

    public Taxon getCurrent() {
        return current;
    }
}
