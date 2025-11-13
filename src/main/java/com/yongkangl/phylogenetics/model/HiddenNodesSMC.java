package com.yongkangl.phylogenetics.model;

import com.yongkangl.phylogenetics.io.TreeNode;
import com.yongkangl.phylogenetics.sampler.HiddenNodes;
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

public class HiddenNodesSMC implements Runnable {
    private ExecutorService executorService;
    private final TreeNode root;
    private final SiteModel[] modelSequence;
    private final int nSteps;
    private final int nParticles;
    private final int mutationSteps;
    private final HiddenNodes[] hiddenNodes;
    private final double[] stepLogWeights;
    private double ess;

    public HiddenNodesSMC(TreeNode root, SiteModel proposal, SiteModel target, int nSteps, int nParticles, int mutationSteps) {
        this.root = root;
        this.modelSequence = SiteModel.temper(proposal, target, nSteps);
        this.nSteps = nSteps;
        this.nParticles = nParticles;
        this.mutationSteps = mutationSteps;

        hiddenNodes = new HiddenNodes[nParticles];
        stepLogWeights = new double[nSteps];
    }

    public double getLogWeight() {
        return Utils.sum(stepLogWeights);
    }

    public double getESS() {
        return ess;
    }

    public void run() {
        boolean toShutdown = false;
        if (executorService == null) {
            toShutdown = true;
            executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        }

        Future<HiddenNodes>[] futures = new Future[nParticles];

        for (int i = 0; i < nParticles; i++) {
            futures[i] = executorService.submit(() -> new HiddenNodes(nSteps, root, mutationSteps, modelSequence));
        }

        for (int i = 0; i < nParticles; i++) {
            try {
                hiddenNodes[i] = futures[i].get();
                System.out.println(hiddenNodes[i].mutations());
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        double[] logWeights = new double[nParticles];

        for (int i = 0; i < nSteps; i++) {
            calculateStepLogWeights();
            for (int j = 0; j < nParticles; j++) {
                logWeights[j] = hiddenNodes[j].getStepLogWeight(false);
            }
            stepLogWeights[i] = Utils.logSum(logWeights) - Math.log(nParticles);
            System.out.println(stepLogWeights[i] + " (ESS: " + Utils.ESS(logWeights) + ")");
            if (i != nSteps - 1) {
                // Advanced in ``resample''
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
            HiddenNodes hiddenNodeSet = hiddenNodes[i];
            futures[i] = executorService.submit(hiddenNodeSet::mutate);
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
            HiddenNodes hiddenNodeSet = hiddenNodes[i];
            futures[i] = executorService.submit(() -> hiddenNodeSet.calculateStepLogWeight());
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
            possibleParticles.add(new Pair<>(i, Math.exp(hiddenNodes[i].getStepLogWeight(true))));
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
            hiddenNodes[toPaste.get(i)].copyFrom(hiddenNodes[toCopy.get(i)]);
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (HiddenNodes hiddenNode : hiddenNodes) {
            sb.append(hiddenNode.toString()).append("\n");
        }
        return sb.toString();
    }
}
