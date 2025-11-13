package com.yongkangl.phylogenetics.io;

import com.yongkangl.phylogenetics.util.Taxon;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class FastaParser {
    private final Map<String, Taxon> taxons;
    private int alignmentLen;
    public FastaParser (String filePath) {
        taxons = new ConcurrentHashMap<>();
        alignmentLen = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line;
            String header = null;
            StringBuilder sequence = new StringBuilder();
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    // Save previous sequence
                    if (header != null) {
                        Taxon taxon = new Taxon(sequence.toString());
                        taxons.put(header, taxon);
                        if (alignmentLen == 0){
                            alignmentLen = taxon.len();
                        } else if (alignmentLen != taxon.len()) {
                            System.err.println("Variable Sequence Lengths!");
                        }
                    }
                    // Start new sequence
                    header = line.substring(1).trim();
                    sequence = new StringBuilder();
                } else {
                    sequence.append(line.trim());
                }
            }
            if (header != null) {
                Taxon taxon = new Taxon(sequence.toString());
                taxons.put(header, taxon);
                if (alignmentLen == 0){
                    alignmentLen = taxon.len();
                } else if (alignmentLen != taxon.len()) {
                    System.err.println("Variable Sequence Lengths!");
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public Taxon queryTaxon(String name) {
        return taxons.get(name);
    }

    public double[] queryPartial(String name) {
        Taxon taxon = queryTaxon(name);
        double[] partial = new double[alignmentLen << 2];
        for (int i = 0; i < alignmentLen; i++) {
            for (int base = 0; base < 4; base++) {
                partial[(i << 2) + base] = 0.0;
            }
            partial[(i << 2) + taxon.get(i)] = 1.0;
        }
        return partial;
    }

    public int getAlignmentLength() {
        return alignmentLen;
    }
}