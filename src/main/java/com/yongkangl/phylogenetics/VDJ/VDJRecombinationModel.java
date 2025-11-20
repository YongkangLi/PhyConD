package com.yongkangl.phylogenetics.VDJ;

import com.yongkangl.phylogenetics.substitution.Nucleotides;
import com.yongkangl.phylogenetics.util.Taxon;
import com.yongkangl.phylogenetics.util.Utils;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class VDJRecombinationModel {
    private boolean initialized = false;
    private VGeneSegmentModel vGeneSegmentModel;
    private DGeneSegmentModel dGeneSegmentModel;
    private JGeneSegmentModel jGeneSegmentModel;

    private int vIdx;
    private String vGene;
    private int[] l1s;
    private double[] logPvs;

    private int jIdx;
    private String jGene;
    private int[] l4s;
    private double[] logPjs;

    private Map<Integer, Double> logP2s;
    private Map<Integer, Double> logP3s;

    private List<Triple<Integer, Integer, Double>> logPs;

    private double logP;

    public VDJRecombinationModel(String folder) {
        File test = new File(folder);
        if (test.exists() && test.isDirectory()) {
            try {
                this.vGeneSegmentModel = new VGeneSegmentModel(folder + "/vgene");
                this.dGeneSegmentModel = new DGeneSegmentModel(folder + "/dgene");
                this.jGeneSegmentModel = new JGeneSegmentModel(folder + "/jgene");
                initialized = true;
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }

    public double computeVDJProb(String sequence) {
        if (initialized) {
            List<Pair<Integer, Pair<int[], double[]>>> vs = vGeneSegmentModel.filterV(sequence);
            if (vs.isEmpty()) {
                return Double.NEGATIVE_INFINITY;
            }
            Pair<Integer, Pair<int[], double[]>> v = vs.get(0);
            vIdx = v.getKey();
            vGene = vGeneSegmentModel.translateIdToName(vIdx);
            l1s = v.getValue().getFirst();
            logPvs = v.getValue().getSecond();

            List<Pair<Integer, Pair<int[], double[]>>> js = jGeneSegmentModel.filterJ(vIdx, sequence);
            if (js.isEmpty()) {
                return Double.NEGATIVE_INFINITY;
            }
            Pair<Integer, Pair<int[], double[]>> j = js.get(0);
            jIdx = j.getKey();
            jGene = jGeneSegmentModel.translateIdToName(jIdx);
            l4s = j.getValue().getFirst();
            logPjs = j.getValue().getSecond();

            logP2s = new ConcurrentHashMap<>();
            for (int l2 = l1s[0] + 1; l2 < l4s[l4s.length - 1] - 1; l2++) {
                double logP2 = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < l1s.length; i++) {
                    int l1 = l1s[i];
                    if (l1 < l2) {
                        logP2 = Utils.logAdd(logP2, VDInsertionSegmentModel.getInsertionSequenceProb(sequence.substring(l1, l2), sequence.charAt(l1 - 1)) + logPvs[i]);
                    }
                }
                logP2s.put(l2, logP2);
            }

            logP3s = new ConcurrentHashMap<>();
            for (int l3 = l1s[0] + 1; l3 < l4s[l4s.length - 1]; l3++) {
                double logP3 = Double.NEGATIVE_INFINITY;
                for (int i = 0; i < l4s.length; i++) {
                    int l4 = l4s[i];
                    if (l3 < l4) {
                        logP3 = Utils.logAdd(logP3, DJInsertionSegmentModel.getInsertionSequenceProb(sequence.substring(l3, l4), sequence.charAt(l4)) + logPjs[i]);
                    }
                }
                logP3s.put(l3, logP3);
            }

            logPs = new ArrayList<>();
            logP = Double.NEGATIVE_INFINITY;
            List<Triple<Integer, Integer, Double>> ds = dGeneSegmentModel.analyzeD(vIdx, jIdx, sequence, l1s[0] + 1, l4s[l4s.length - 1] - 1);
            for (Triple<Integer, Integer, Double> d : ds) {
                int l2 = d.getLeft();
                int l3 = d.getMiddle();
                double logProb = d.getRight() + logP2s.get(l2) + logP3s.get(l3);
                logPs.add(Triple.of(l2, l3, logProb));
                logP = Utils.logAdd(logP, logProb);
            }
            logPs.sort(Comparator.comparing(Triple::getRight));

            return logP;
        } else {
            double[] stationary = new double[]{0.19610813492990845, 0.32047847148641945, 0.2074584913818569, 0.27595490220181523};
            logP = 0.0;
            Taxon taxon = new Taxon(sequence);
            for (int i = 0; i < sequence.length(); i++) {
                logP += Math.log(stationary[taxon.get(i)]);
            }
            return logP;
        }
    }
}