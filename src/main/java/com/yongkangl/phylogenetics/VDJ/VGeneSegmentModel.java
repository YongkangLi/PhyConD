package com.yongkangl.phylogenetics.VDJ;

import java.util.*;
import java.io.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

import com.fasterxml.jackson.databind.*;
import com.fasterxml.jackson.core.type.TypeReference;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.util.Pair;

public class VGeneSegmentModel {

    private double[] vChoiceProbs;
    private double[][] vDelProbs;
    private Map<String, List<Pair<Integer, Integer>>> prefixMap;
    private Map<Integer, String> idToName;

    public VGeneSegmentModel(String jsonDirPath) throws IOException {
        ObjectMapper mapper = new ObjectMapper();
        vChoiceProbs = mapper.readValue(new File(jsonDirPath + "/v_choice_probs.json"), double[].class);
        vDelProbs = mapper.readValue(new File(jsonDirPath + "/v_del_probs.json"), double[][].class);

        Map<String, List<List<Integer>>> rawPrefixMap =
            mapper.readValue(new File(jsonDirPath + "/v_prefix_source_map.json"),
                             new TypeReference<Map<String, List<List<Integer>>>>() {});

        prefixMap = new HashMap<>();
        for (Map.Entry<String, List<List<Integer>>> e : rawPrefixMap.entrySet()) {
            List<Pair<Integer, Integer>> pairs = new ArrayList<>();
            for (List<Integer> entry : e.getValue()) {
                pairs.add(new Pair<>(entry.get(0), entry.get(1)));
            }
            prefixMap.put(e.getKey(), pairs);
        }

        Map<String, String> nameMap =
            mapper.readValue(new File(jsonDirPath + "/v_gene_map.json"),
                             new TypeReference<Map<String, String>>() {});

        idToName = new HashMap<>();
        for (Map.Entry<String, String> e : nameMap.entrySet()) {
            idToName.put(Integer.parseInt(e.getKey()), e.getValue());
        }
    }

    public List<Pair<Integer, Pair<int[], double[]>>> filterV (String sequence) {
        Map<Integer, List<Pair<Integer, Double>>> map = new ConcurrentHashMap<>();
        for (int l1 = 290 - 27; l1 <= 305; l1++) {
            Pair<Double, List<Map.Entry<Integer, Double>>> pair = getProbabilityAndPosterior(sequence.substring(0, l1));
            double prob = pair.getFirst();
            if (prob > 0) {
                double logProb = Math.log(prob);
                for (Map.Entry<Integer, Double> entry : pair.getSecond()) {
                    double val = logProb + Math.log(entry.getValue());
                    if (map.containsKey(entry.getKey())) {
                        map.get(entry.getKey()).add(new Pair<>(l1, val));
                    } else {
                        map.put(entry.getKey(), new ArrayList<>(Collections.singletonList(new Pair<>(l1, val))));
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

    public Pair<Double, List<Map.Entry<Integer, Double>>> getProbabilityAndPosterior(String sub) {
        List<Pair<Integer, Integer>> sources = prefixMap.getOrDefault(sub, Collections.emptyList());
        double total = 0.0;
        Map<Integer, Double> weighted = new HashMap<>();

        for (Pair<Integer, Integer> pair : sources) {
            int vIdx = pair.getFirst();
            int d3 = pair.getSecond();
            if (vIdx >= vChoiceProbs.length || d3 >= vDelProbs[0].length) continue;

            double p = vChoiceProbs[vIdx] * vDelProbs[vIdx][d3];
            total += p;
            weighted.put(vIdx, weighted.getOrDefault(vIdx, 0.0) + p);
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

    public void sample(Random random) {
        EnumeratedIntegerDistribution choiceDistribution = new EnumeratedIntegerDistribution(
            IntStream.range(0, vChoiceProbs.length).toArray(),
            vChoiceProbs
        );
        int vIdx = choiceDistribution.sample();
        EnumeratedIntegerDistribution delDistribution = new EnumeratedIntegerDistribution(
            IntStream.range(0, vDelProbs[vIdx].length).toArray(),
            vDelProbs[vIdx]
        );

        // To return the subsequence
    }
}
