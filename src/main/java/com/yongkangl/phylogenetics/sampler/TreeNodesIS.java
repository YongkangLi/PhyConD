package com.yongkangl.phylogenetics.sampler;


import com.yongkangl.phylogenetics.io.TreeNode;
import com.yongkangl.phylogenetics.substitution.SiteModel;
import com.yongkangl.phylogenetics.util.Utils;
import org.apache.commons.math3.util.Pair;

import java.util.Arrays;

public class TreeNodesIS {
    private final TreeNode tree;
    private final SiteModel proposal;
    private final SiteModel target;
    private final int nTrees;
    private final int nSteps;
    private final int nParticles;
    private final double mutabilityThreshold;
    private final double threshold;
    private final boolean fixInternalNodes;
    private final boolean mixUCA;
    private final TreeNodesSampler[] samplers;
    private double ess;
    private double[] subESSs;
    private double[][] subESSss;

    public TreeNodesIS(TreeNode root, SiteModel proposal, SiteModel target, int nTrees, int nSteps, int nParticles, double mutabilityThreshold, double threshold, boolean fixInternalNodes, boolean mixUCA) {
        tree = root;
        this.proposal = proposal;
        this.target = target;


        this.nTrees = nTrees;
        this.nSteps = nSteps;
        this.mutabilityThreshold = mutabilityThreshold;
        this.threshold = threshold;
        this.nParticles = nParticles;

        this.fixInternalNodes = fixInternalNodes;
        this.mixUCA = mixUCA;

        samplers = new TreeNodesSampler[nTrees];

        subESSs = new double[nTrees];
        subESSss = new double[nTrees][];
    }

    public double run() {
//        System.out.println("Log Tree Weight, Log Sample Weight");
        double[] logWeights = new double[nTrees];
        for (int i = 0; i < nTrees; i++) {
            samplers[i] = new TreeNodesSampler(tree, proposal, target, nSteps, nParticles, mutabilityThreshold, threshold, fixInternalNodes, mixUCA);
//            samplers[i].run();
//            logWeights[i] = samplers[i].getLogWeight();
//            System.out.println(logWeights[i]);
//            subESSs[i] = samplers[i].getMinESS();
//            subESSss[i] = samplers[i].getESSs();
        }
//        double logSumWeights = Utils.logSum(logWeights);
//        ess =  Math.exp(2 * logSumWeights - Utils.logSquareSum(logWeights));
//        return logSumWeights - Math.log(nTrees);
        return 0.0;
    }

    public double getESS() {
        return ess;
    }

    public double[] getSubESS() {
        return subESSs;
    }

    public String getSubESSs() {
        String s = "{\n";
        for (int i = 0; i < subESSss.length; i++) {
            s += (Arrays.toString(subESSss[i]) + "\n");
        }
        s += "}";
        return s;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < nTrees; i++) {
            sb.append(samplers[i].toString()).append("\n");
        }
        return sb.toString();
    }
}
