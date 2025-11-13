package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Mutation;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Yongkang Li
 * @version $Id$
 */
public class BlockPathSampler {
    /**
     * Single Block Path Sampler
     * Sample a series of Mutations that occur within a single block of a given Sequence.
     * This sampler is based on a block-wise site model.
     */

    private final int startSite;
    private final int startState;
    private final int endState;
    private final double T;
    private final int blockSize;
    private final BlockwiseSiteModel blockwiseSiteModel;
    private final List<Mutation> mutations;
    private final List<Mutation>[] sitePaths;

    public BlockPathSampler(int site, int start, int end, double T, BlockwiseSiteModel blockwiseSiteModel, boolean hobolth) {
        this.startSite = site;
        startState = start;
        endState = end;
        this.T = T;
        this.blockwiseSiteModel = blockwiseSiteModel;
        this.blockSize = blockwiseSiteModel.getBlockSize();

        if (startState != endState) {
            // If start != end, draw the first mutation conditional on this information
            // in order to speed up the accept/reject algorithm
            mutations = hobolth ? hobolthRejectionSampling() : rejectionSampling();
        } else {
            mutations = rejectionSampling();
        }
        sitePaths = splitSitePaths();
    }

    private ArrayList<Mutation> rejectionSampling() {
        ArrayList<Mutation> potential;
        do {
            potential = forwardSampling(startState, 0, T);
        } while (potential.isEmpty() ? (startState != endState) : (potential.get(potential.size() - 1).getTo() != endState));
        return potential;
    }

    private ArrayList<Mutation> hobolthRejectionSampling() {
        // Draw initial time using inverse CDF method (see eqn (2.1) of Hobolth paper)
        boolean success = false;
        ArrayList<Mutation> potential = new ArrayList<>();
        do {
            double rate = blockwiseSiteModel.getRate(startSite, startState);
            double t = -Math.log(1 - (new UniformRealDistribution(0.0, 1.0).sample() * (1 - Math.exp(T * rate)))) / (-rate);
            int sampledStart = sampleNextState(startState);
            Mutation sampledMutation = new Mutation(t, startSite, sampledStart);

            int maxTrials = 100000;
            for (int trials = 0; trials < maxTrials; trials++) {
                potential = new ArrayList<>();
                potential.add(sampledMutation);
                potential.addAll(forwardSampling(sampledStart, t, T));
                if (potential.get(potential.size() - 1).getTo() == endState) {
                    success = true;
                    break;
                }
            }
        } while (!success);
        return potential;
    }

    private ArrayList<Mutation> forwardSampling(int start, double t, double T) {
        int state = start;
        double time = t;
        ArrayList<Mutation> potential = new ArrayList<>();
        while (time < T) {
            double dt = new ExponentialDistribution(- 1 / blockwiseSiteModel.getRate(state, state)).sample();
            time += dt;
            if (time < T) {
                state = sampleNextState(state);
                potential.add(new Mutation(time, startSite, state));
            }
        }
        return potential;
    }

    private int sampleNextState(int state) {
        List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
        for (int i = 0; i < (1 << (2 * blockSize)); i++) {
            if (i != state) {
                possibleStates.add(new Pair<>(i, - blockwiseSiteModel.getRate(state, i) / blockwiseSiteModel.getRate(state, state)));
            }
        }
        return new EnumeratedDistribution<>(possibleStates).sample();
    }

    public Mutation get(int i) {
        return mutations.get(i);
    }

    public List<Mutation> getMutations() {
        return mutations;
    }

    public List<Mutation>[] getSitePaths() {
        return sitePaths;
    }

    private List<Mutation>[] splitSitePaths() {
        List<Mutation>[] paths = new ArrayList[blockSize];
        for (int i = 0; i < blockSize; i++) {
            paths[i] = new ArrayList<>();
        }

        int[] currentState = BlockwiseSiteModel.decode(blockSize, startState);
        for (Mutation mutation : mutations) {
            int[] nextState = BlockwiseSiteModel.decode(blockSize, mutation.getTo());
            for (int offset = 0; offset < blockSize; offset++) {
                if (currentState[offset] != nextState[offset]) {
                    paths[offset].add(new Mutation(mutation.getTime(), startSite + offset, nextState[offset]));
                }
            }
            currentState = nextState;
        }
        return paths;
    }

    public static boolean isEmpty(ArrayList<Mutation>[] mutations) {
        for (ArrayList<Mutation> m : mutations) {
            if (!m.isEmpty()) {
                return false;
            }
        }
        return true;
    }
}
