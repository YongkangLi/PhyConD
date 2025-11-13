package com.yongkangl.phylogenetics.VDJ;

import java.util.*;
import java.io.*;
import java.util.concurrent.ConcurrentHashMap;

import com.fasterxml.jackson.databind.*;
import com.fasterxml.jackson.core.type.TypeReference;
import org.apache.commons.math3.util.Pair;
import org.apache.commons.lang3.tuple.Triple;

public class DGeneSegmentModel {

    private double[][][] dChoiceProbs;
    private double[][] d5DelProbs;
    private double[][][] d3DelProbs;
    private Map<String, List<Triple<Integer, Integer, Integer>>> affixMap;
    private Map<Integer, String> idToName;

    public DGeneSegmentModel(String jsonDirPath) throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        dChoiceProbs = mapper.readValue(new File(jsonDirPath + "/d_choice_probs.json"), double[][][].class);
        d5DelProbs = mapper.readValue(new File(jsonDirPath + "/d_5_del_probs.json"), double[][].class);
        d3DelProbs = mapper.readValue(new File(jsonDirPath + "/d_3_del_probs.json"), double[][][].class);

        Map<String, List<List<Integer>>> rawAffixMap =
            mapper.readValue(new File(jsonDirPath + "/d_affix_source_map.json"),
                             new TypeReference<Map<String, List<List<Integer>>>>() {});

        affixMap = new HashMap<>();
        for (Map.Entry<String, List<List<Integer>>> e : rawAffixMap.entrySet()) {
            List<Triple<Integer, Integer, Integer>> triples = new ArrayList<>();
            for (List<Integer> triplet : e.getValue()) {
                triples.add(Triple.of(triplet.get(0), triplet.get(1), triplet.get(2)));
            }
            affixMap.put(e.getKey(), triples);
        }

        Map<String, String> nameMap =
            mapper.readValue(new File(jsonDirPath + "/d_gene_map.json"),
                             new TypeReference<Map<String, String>>() {});

        idToName = new HashMap<>();
        for (Map.Entry<String, String> e : nameMap.entrySet()) {
            idToName.put(Integer.parseInt(e.getKey()), e.getValue());
        }
    }

    public List<Pair<Integer, Triple<int[], int[], double[]>>> filterD (int vIdx, int jIdx, String sequence) {
        int L = sequence.length();
        Map<Integer, List<Triple<Integer, Integer, Double>>> map = new ConcurrentHashMap<>();
        for (int l2 = 290-27; l2 < L - 21; l2++) {
            for (int l3 = l2 + 1; l3 <= L - 21; l3++) {
                Pair<Double, List<Map.Entry<Integer, Double>>> pair = getProbabilityAndPosterior(vIdx, jIdx, sequence.substring(l2, l3));
                double prob = pair.getFirst();
                if (prob > 0) {
                    for (Map.Entry<Integer, Double> entry : pair.getSecond()) {
                        double val = prob * entry.getValue();
                        if (map.containsKey(entry.getKey())) {
                            map.get(entry.getKey()).add(Triple.of(l2, l3, val));
                        } else {
                            map.put(entry.getKey(), new ArrayList<>(Collections.singletonList(Triple.of(l2, l3, val))));
                        }
                    }
                }
            }
        }
        List<Pair<Integer, Triple<int[], int[], double[]>>> pairs = new ArrayList<>();
        for (Map.Entry<Integer, List<Triple<Integer, Integer, Double>>> entry : map.entrySet()) {
            List<Triple<Integer, Integer, Double>> list = entry.getValue();
            int size = list.size();
            int[] s = new int[size];
            int[] e = new int[size];
            double[] p = new double[size];
            list.sort(Comparator.comparing(Triple::getLeft));
            for (int i = 0; i < size; i++) {
                s[i] = list.get(i).getLeft();
                e[i] = list.get(i).getMiddle();
                p[i] = list.get(i).getRight();
            }
            pairs.add(new Pair<>(entry.getKey(), Triple.of(s, e, p)));
        }
        return pairs;
    }

    public List<Triple<Integer, Integer, Double>> analyzeD(int vIdx, int jIdx, String sequence, int minl2, int maxl3) {
        int L = sequence.length();
        List<Triple<Integer, Integer, Double>> list = new ArrayList<>();
        for (int l2 = minl2; l2 < maxl3; l2++) {
            for (int l3 = l2 + 1; l3 <= maxl3; l3++) {
                Pair<Double, List<Map.Entry<Integer, Double>>> pair = getProbabilityAndPosterior(vIdx, jIdx, sequence.substring(l2, l3));
                double prob = pair.getFirst();
                if (prob > 0) {
                    list.add(Triple.of(l2, l3, Math.log(prob)));
                }
            }
        }
        list.sort(Comparator.comparing(Triple::getRight));
        return list;
    }

    public Pair<Double, List<Map.Entry<Integer, Double>>> getProbabilityAndPosterior(int vIdx, int jIdx, String subsequence) {
        List<Triple<Integer, Integer, Integer>> sources = affixMap.getOrDefault(subsequence, Collections.emptyList());
        double total = 0.0;
        Map<Integer, Double> weightMap = new HashMap<>();

        for (Triple<Integer, Integer, Integer> t : sources) {
            int dIdx = t.getLeft();
            int d5 = t.getMiddle();
            int d3 = t.getRight();

            if (vIdx >= dChoiceProbs.length || jIdx >= dChoiceProbs[0].length || dIdx >= dChoiceProbs[0][0].length)
                continue;
            if (d5 >= d5DelProbs[0].length || d3 >= d3DelProbs[0][0].length) continue;

            double p = dChoiceProbs[vIdx][jIdx][dIdx] * d5DelProbs[dIdx][d5] * d3DelProbs[dIdx][d5][d3];
            total += p;
            weightMap.put(dIdx, weightMap.getOrDefault(dIdx, 0.0) + p);
        }

        List<Map.Entry<Integer, Double>> posterior = new ArrayList<>();
        for (Map.Entry<Integer, Double> e : weightMap.entrySet()) {
            posterior.add(new AbstractMap.SimpleEntry<>(e.getKey(), e.getValue() / total));
        }

        return new Pair<>(total, posterior);
    }
}
