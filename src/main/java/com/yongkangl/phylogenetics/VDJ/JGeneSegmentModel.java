package com.yongkangl.phylogenetics.VDJ;

import java.util.*;
import java.io.*;
import java.util.concurrent.ConcurrentHashMap;

import com.fasterxml.jackson.databind.*;
import com.fasterxml.jackson.core.type.TypeReference;
import org.apache.commons.math3.util.Pair;

public class JGeneSegmentModel {

    private double[][] jChoiceProbs;
    private double[][] jDelProbs;
    private Map<String, List<Pair<Integer, Integer>>> affixMap;
    private Map<Integer, String> idToName;

    public JGeneSegmentModel(String jsonDirPath) throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        jChoiceProbs = mapper.readValue(new File(jsonDirPath + "/j_choice_probs.json"), double[][].class);
        jDelProbs = mapper.readValue(new File(jsonDirPath + "/j_5_del_probs.json"), double[][].class);

        Map<String, List<List<Integer>>> rawAffixMap =
            mapper.readValue(new File(jsonDirPath + "/j_affix_source_map.json"),
                             new TypeReference<Map<String, List<List<Integer>>>>() {});

        affixMap = new HashMap<>();
        for (Map.Entry<String, List<List<Integer>>> e : rawAffixMap.entrySet()) {
            List<Pair<Integer, Integer>> entries = new ArrayList<>();
            for (List<Integer> pair : e.getValue()) {
                entries.add(new Pair<>(pair.get(0), pair.get(1)));
            }
            affixMap.put(e.getKey(), entries);
        }

        Map<String, String> nameMap =
            mapper.readValue(new File(jsonDirPath + "/j_gene_map.json"),
                             new TypeReference<Map<String, String>>() {});

        idToName = new HashMap<>();
        for (Map.Entry<String, String> e : nameMap.entrySet()) {
            idToName.put(Integer.parseInt(e.getKey()), e.getValue());
        }
    }

    public List<Pair<Integer, Pair<int[], double[]>>> filterJ (int vIdx, String sequence) {
        int L = sequence.length();
        Map<Integer, List<Pair<Integer, Double>>> map = new ConcurrentHashMap<>();
        for (int l4 = L - 63; l4 <= L - 21; l4++) {
            Pair<Double, List<Map.Entry<Integer, Double>>> pair = getProbabilityAndPosterior(vIdx, sequence.substring(l4, L));
            double prob = pair.getFirst();
            if (prob > 0) {
                double logProb = Math.log(prob);
                for (Map.Entry<Integer, Double> entry : pair.getSecond()) {
                    double val = logProb + Math.log(entry.getValue());
                    if (map.containsKey(entry.getKey())) {
                        map.get(entry.getKey()).add(new Pair<>(l4, val));
                    } else {
                        map.put(entry.getKey(), new ArrayList<>(Collections.singletonList(new Pair<>(l4, val))));
                    }
                }
            }
        }
        List<Pair<Integer, Pair<int[], double[]>>> pairs = new ArrayList<>();
        for (Map.Entry<Integer, List<Pair<Integer, Double>>> entry : map.entrySet()) {
            List<Pair<Integer, Double>> list = entry.getValue();
            int size = list.size();
            int[] l = new int[size];
            double[] logP = new double[size];
            list.sort(Comparator.comparing(Pair::getFirst));
            for (int i = 0; i < size; i++) {
                l[i] = list.get(i).getFirst();
                logP[i] = list.get(i).getSecond();
            }
            pairs.add(new Pair<>(entry.getKey(), new Pair<>(l, logP)));
        }
        return pairs;
    }

    public Pair<Double, List<Map.Entry<Integer, Double>>> getProbabilityAndPosterior(int vIdx, String affix) {
        List<Pair<Integer, Integer>> sources = affixMap.getOrDefault(affix, Collections.emptyList());
        double total = 0.0;
        Map<Integer, Double> weighted = new HashMap<>();

        for (Pair<Integer, Integer> pair : sources) {
            int jIdx = pair.getFirst();
            int d5 = pair.getSecond();
            if (jIdx >= jChoiceProbs[0].length || d5 >= jDelProbs[0].length) continue;

            double p = jChoiceProbs[vIdx][jIdx] * jDelProbs[jIdx][d5];
            total += p;
            weighted.put(jIdx, weighted.getOrDefault(jIdx, 0.0) + p);
        }

        List<Map.Entry<Integer, Double>> posterior = new ArrayList<>();
        for (Map.Entry<Integer, Double> e : weighted.entrySet()) {
            posterior.add(new AbstractMap.SimpleEntry<>(e.getKey(), e.getValue() / total));
        }

        return new Pair<>(total, posterior);
    }

    public String translateIdToName(int id) {
        return idToName.getOrDefault(id, "Unknown");
    }
}
