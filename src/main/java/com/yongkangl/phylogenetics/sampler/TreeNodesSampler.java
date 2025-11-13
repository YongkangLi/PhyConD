package com.yongkangl.phylogenetics.sampler;

import com.yongkangl.phylogenetics.io.TreeNode;
import com.yongkangl.phylogenetics.substitution.*;
import com.yongkangl.phylogenetics.util.Taxon;
import com.yongkangl.phylogenetics.util.Utils;

import java.util.ArrayList;
import java.util.List;

public class TreeNodesSampler implements Runnable {
    private final SiteModel proposal;
    private final SiteModel target;
    private final int nSteps;
    private final int nParticles;
    private final double mutabilityThreshold;
    private final double threshold;

    private double logTreeWeight;
    private double[] ESSs;

    private Taxon sampledUCA;

    private class Lineage {
        private final Taxon initial;
        private final Taxon terminal;
        private final double T;
        private double logWeight;

        public Lineage(Taxon initial, Taxon terminal, double T) {
            this.initial = initial;
            this.terminal = terminal;
            this.T = T;
        }

        public void setLogWeight(double logWeight) {
            this.logWeight = logWeight;
        }
    }

    private void sampleLineages(TreeNode root) {
        root.sampleClone(proposal);
    }

    private void sampleLineages(TreeNode root, boolean fixed, boolean mixUCA) {
//        sampledUCA = root.getProposalTaxonSampler().sample();
        if (!fixed) {
//            if (mixUCA) {
//                Taxon CDR3Sample = root.sample(190, 235);
//                sampledUCA.set(190, CDR3Sample.subSequence(190, 235));
//                logSampleWeight += root.evaluateTargetLogLikelihood(sampledUCA, 190, 235) - root.evaluateProposalLogLikelihood(sampledUCA, 190, 235);
//            }
            root.sample();
            sampledUCA = root.getFixedTaxon();
            sampleLineages(root);
        }
    }

    public TreeNodesSampler(TreeNode root, SiteModel proposal, SiteModel target, int nSteps, int nParticles, double mutabilityThreshold, double threshold, boolean fixed, boolean mixUCA) {
        this.proposal = proposal;
        this.target = target;
        this.nSteps = nSteps;
        this.nParticles = nParticles;
        this.mutabilityThreshold = mutabilityThreshold;
        this.threshold = threshold;

        sampleLineages(root, fixed, mixUCA);
    }

    @Override
    public void run() {
//        int size = lineages.size();
//        SequenceMonteCarlo[] sequenceMonteCarlos = new SequenceMonteCarlo[size];
//        double[] logWeights = new double[size];
//        ESSs = new double[size];
//
//        for (int i = 0; i < size; i++) {
//            Lineage lineage = lineages.get(i);
//            sequenceMonteCarlos[i] = new SequenceMonteCarlo(lineage.initial, lineage.terminal, lineage.T, 1, mutabilityThreshold, proposal, target, 1, nParticles, 0, threshold);
//            Pair<Double, Double> results = sequenceMonteCarlos[i].importanceSampling();
//            BlockwiseSubstitutionModel ism = new BlockwiseSubstitutionModel(proposal, target.getContextLength(), 1, Double.POSITIVE_INFINITY, lineage.initial);
//            logWeights[i] = results.getFirst() + sequenceMonteCarlos[i].getProposalLogLikelihood() - ism.logLikelihood(lineage.initial, lineage.terminal, lineage.T);
//            ESSs[i] = results.getSecond();
//            lineage.setLogWeight(logWeights[i]);
//            System.out.println(lineage.initial + " -> " + lineage.terminal + ": " + logWeights[i] + " (ESS: " + ESSs[i] + ")");
//        }
//        logTreeWeight = Utils.sum(logWeights);
    }

//    @Override
//    public void run() {
//        int size = lineages.size();
//        double[] logWeights = new double[size];
//        ESSs = new double[size];
//        IndependentLogLikelihood ill = new IndependentLogLikelihood(BlockwiseSiteModel.getKernel(1, new int[]{0}, new int[]{0}, pathProposal));
//
//        for (int i = 0; i < size; i++) {
//            Lineage lineage = lineages.get(i);
//            double logDSM = DependentLogLikelihood.logLikelihood(target, lineage.initial, lineage.terminal, lineage.T);
//            double logISM = ill.logLikelihood(lineage.initial, lineage.terminal, lineage.T, target.getContextLength());
//            logWeights[i] = logDSM - logISM;
//            lineage.setLogWeight(logWeights[i]);
//            System.out.println(lineage.initial + " -> " + lineage.terminal + ": " + logWeights[i]);
//            ESSs[i] = 1024.0;
//        }
//        logTreeWeight = Utils.sum(logWeights);
//    }

    public double getLogWeight() {
        return logTreeWeight;
    }

    public double getMinESS() {
        return Utils.min(ESSs);
    }

    public double[] getESSs() {
        return ESSs;
    }

    @Override
    public String toString() {
//        StringBuilder sb = new StringBuilder();
//        for (Lineage lineage : lineages) {
//            sb.append(lineage.initial).append(" -> ").append(lineage.terminal).append(": ").append(lineage.logWeight).append("\n");
//        }
//        return sb.toString();
        return sampledUCA.toString();
    }
}
