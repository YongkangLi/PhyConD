package com.yongkangl.phylogenetics.model;

import com.yongkangl.phylogenetics.substitution.*;
import com.yongkangl.phylogenetics.util.Taxon;
import com.yongkangl.phylogenetics.util.Utils;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class PairwiseSequenceSMC implements Runnable {
    private ExecutorService executorService;
    private final Taxon initial;
    private final Taxon terminal;
    private final double T;
    private final int blockSize;
    private final BlockwiseSubstitutionModel blockwiseSubstitutionModel;
    private final SiteModel target;
    private final int nSteps;
    private final int nParticles;
    private final int mutationSteps;
    private final Path[] paths;
    private final double[] stepLogWeights;
    private double ess;
    private final double threshold;

    public PairwiseSequenceSMC(Taxon initial, Taxon terminal, double T, int blockSize, double mutabilityThreshold, SiteModel proposal, SiteModel target, int nSteps, int nParticles, int mutationSteps, double threshold) {
        this.initial = initial;
        this.terminal = terminal;
        this.T = T;
        this.blockSize = blockSize;
        this.blockwiseSubstitutionModel = new BlockwiseSubstitutionModel(proposal, target.getContextLength(), blockSize, mutabilityThreshold, initial);
        this.target = target;
        this.nSteps = nSteps;
        this.nParticles = nParticles;
        this.mutationSteps = mutationSteps;
        this.threshold = threshold;

        paths = new Path[nParticles];
        stepLogWeights = new double[nSteps];
    }

    public PairwiseSequenceSMC(Taxon initial, Taxon terminal, double T, int blockSize, double mutabilityThreshold, SiteModel target, int nSteps, int nParticles, int mutationSteps, double threshold) {
        this(initial, terminal, T, blockSize, mutabilityThreshold, target, target, nSteps, nParticles, mutationSteps, threshold);
    }

    public double getLogWeight() {
        return Utils.sum(stepLogWeights);
    }

    public double getESS() {
        return ess;
    }

    public Pair<Double, Double> importanceSampling() {
        boolean toShutdown = false;
        if (executorService == null) {
            toShutdown = true;
            executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        }
        for (int i = 0; i < nParticles; i++) {
            paths[i] = new Path(nSteps, initial, terminal, T, blockSize, blockwiseSubstitutionModel, target, mutationSteps);
        }
        double[] logWeights = new double[nParticles];
        calculateStepLogWeights();
        if (toShutdown) {
            executorService.shutdown();
        }
        for (int j = 0; j < nParticles; j++) {
            logWeights[j] = paths[j].getStepLogWeight(true);
        }
        double logWeight = Utils.logSum(logWeights) - Math.log(nParticles);
        double ESS = Utils.ESS(logWeights);
        return new Pair<>(logWeight, ESS);
    }

    public double getProposalLogLikelihood() {
        return blockwiseSubstitutionModel.logLikelihood(initial, terminal, T);
    }

    public void run() {
        boolean toShutdown = false;
        if (executorService == null) {
            toShutdown = true;
            executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        }

        for (int i = 0; i < nParticles; i++) {
            paths[i] = new Path(nSteps, initial, terminal, T, blockSize, blockwiseSubstitutionModel, target, mutationSteps);
        }
        double[] logWeights = new double[nParticles];

        for (int i = 0; i < nSteps; i++) {
            calculateStepLogWeights();
            for (int j = 0; j < nParticles; j++) {
                logWeights[j] = paths[j].getStepLogWeight(false);
            }
            stepLogWeights[i] = Utils.logSum(logWeights) - Math.log(nParticles);
            if (i != nSteps - 1) {
                resample();
                mutate();
            }
        }

        ess = Utils.ESS(logWeights);

        if (toShutdown) {
            executorService.shutdown();
        }
    }

    private void mutate() {
        Future<?>[] futures = new Future[nParticles];
        for (int i = 0; i < nParticles; i++) {
            Path path = paths[i];
            futures[i] = executorService.submit(path::determinsticScan);
        }
        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private void calculateStepLogWeights() {
        Future<?>[] futures = new Future[nParticles];
        for (int i = 0; i < nParticles; i++) {
            Path path = paths[i];
            futures[i] = executorService.submit(() -> path.calculateStepLogWeight());
        }
        for (int i = 0; i < nParticles; i++) {
            try {
                futures[i].get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private void resample() {
        List<Pair<Integer, Double>> possibleParticles = new ArrayList<>();
        for (int i = 0; i < nParticles; i++) {
            possibleParticles.add(new Pair<>(i, Math.exp(paths[i].getStepLogWeight(true))));
        }
        EnumeratedDistribution<Integer> dist = new EnumeratedDistribution<>(possibleParticles);
        int[] sampledCount = new int[nParticles];
        for (int i = 0; i < nParticles; i++) {
            sampledCount[i] = 0;
        }
        for (int i = 0; i < nParticles; i++) {
            sampledCount[dist.sample()]++;
        }
        List<Integer> toCopy = new ArrayList<>();
        List<Integer> toPaste = new ArrayList<>();
        for (int i = 0; i < nParticles; i++) {
            if (sampledCount[i] == 0) {
                toPaste.add(i);
            } else if (sampledCount[i] > 1) {
                for (int j = 0; j < sampledCount[i] - 1; j++) {
                    toCopy.add(i);
                }
            }
        }
        for (int i = 0; i < toCopy.size(); i++) {
            paths[toPaste.get(i)].copyFrom(paths[toCopy.get(i)]);
        }
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < nParticles; i++) {
            s.append(paths[i].toString()).append("\n");
        }
        return s.toString();
    }

    public double summaryStatistics() {
        int count = 0;
        for (int i = 0; i < nParticles; i++) {
            count += paths[i].numberOfMutations();
        }
        return ((double) count) / nParticles;
    }
}
