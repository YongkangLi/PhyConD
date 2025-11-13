package com.yongkangl.phylogenetics.VDJ;

import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;

import java.util.Arrays;
import java.util.stream.IntStream;

public class VDInsertionSegmentModel {

    private static final double[] INSERTION_PROB = new double[] {
        0.0156841, 0.0133492, 0.0318305, 0.0536327, 0.0717323, 0.0693499, 0.0708693, 0.069049, 0.0662165, 0.066167, 0.0614392, 0.0567546, 0.048248, 0.0434769, 0.0372504, 0.0319681, 0.0268793, 0.0238947, 0.0174149, 0.0164607, 0.013425, 0.0118015, 0.0103139, 0.00842441, 0.00742351, 0.00613087, 0.00593826, 0.00488709, 0.00476727, 0.00383626, 0.00369003, 0.0032656, 0.00290255, 0.00231115, 0.00221705, 0.00199928, 0.00161274, 0.00166865, 0.00121428, 0.00123824, 0.00112886, 0.00116173, 0.000625371, 0.000682141, 0.00100915, 0.000449861, 0.000618188, 0.000637948, 0.000438992, 0.000137678, 0.000358883, 0.000192473, 0.000448392, 9.58974e-05, 0.000217938, 0.000195527, 0.000250542, 4.3086e-05, 0.00013789, 0.00010244, 0.000332023
    };

    private static final double[][] DINUCLEOTIDE_PROB = new double[][] {
        {0.253311, 0.240039, 0.302671, 0.20398},
        {0.155895, 0.405713, 0.237885, 0.200507},
        {0.221125, 0.169374, 0.424277, 0.185225},
        {0.237648, 0.263564, 0.236513, 0.262274}
    };

    private static double getInsertionLengthProb(int length) {
        if (length < 0 || length >= INSERTION_PROB.length) return 0.0;
        return INSERTION_PROB[length];
    }

    private static double getDinucleotideProb(int fromIndex, int toIndex) {
        if (fromIndex < 0 || fromIndex >= 4 || toIndex < 0 || toIndex >= 4) return 0.0;
        return DINUCLEOTIDE_PROB[fromIndex][toIndex];
    }

    public static double getInsertionSequenceProb(String seq, char nucleotide) {
        int length = seq.length();
        if (length == 0 || length >= INSERTION_PROB.length) return Double.NEGATIVE_INFINITY;

        double logProb = Math.log(INSERTION_PROB[length]);
        int first = ntIndex(seq.charAt(0));
        logProb += Math.log(DINUCLEOTIDE_PROB[ntIndex(nucleotide)][first]);
        for (int i = 0; i < length - 1; i++) {
            int from = ntIndex(seq.charAt(i));
            int to = ntIndex(seq.charAt(i + 1));
            logProb += Math.log(DINUCLEOTIDE_PROB[to][from]);
        }
        return logProb;
    }

    private static int ntIndex(char c) {
        switch (Character.toUpperCase(c)) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return -1;
        }
    }

    private String sample(char nucleotide) {
        EnumeratedIntegerDistribution lengthDist = new EnumeratedIntegerDistribution(
            IntStream.range(0, INSERTION_PROB.length).toArray(),
            INSERTION_PROB
        );
        int length = lengthDist.sample();
        StringBuilder sb = new StringBuilder(length);
        if (length > 0) {
            EnumeratedIntegerDistribution dinucleotideDist = new EnumeratedIntegerDistribution(
                IntStream.range(0, 4).toArray(),
                DINUCLEOTIDE_PROB[ntIndex(nucleotide)]
            );
            sb.append(dinucleotideDist.sample());
        }
        for (int i = 1; i < length; i++) {
            EnumeratedIntegerDistribution dinucleotideDist = new EnumeratedIntegerDistribution(
                    IntStream.range(0, 4).toArray(),
                    DINUCLEOTIDE_PROB[ntIndex(sb.charAt(i - 1))]
            );
            sb.append(dinucleotideDist.sample());
        }
        return sb.toString();
    }
}
