package com.yongkangl.phylogenetics.VDJ;

public class DJInsertionSegmentModel {

    private static final double[] INSERTION_PROB = new double[] {
        0.0172366, 0.0200147, 0.0386513, 0.065385, 0.0735531, 0.0714316, 0.0604028, 0.0661654, 0.0623203, 0.0546426, 0.0535173, 0.0459073, 0.0430245, 0.0392942, 0.0327065, 0.0294925, 0.024265, 0.023849, 0.0191015, 0.0171603, 0.0139567, 0.0119792, 0.0114401, 0.0115652, 0.00812883, 0.00909221, 0.00674935, 0.00693261, 0.00682681, 0.0051028, 0.00533472, 0.00511431, 0.00438476, 0.00320017, 0.00356963, 0.00353649, 0.00222131, 0.00285969, 0.002356, 0.00179071, 0.00198818, 0.00119387, 0.00111915, 0.00190812, 0.00112761, 0.00103636, 0.00105627, 0.000817446, 0.000566862, 0.000859937, 0.000797913, 0.000560903, 0.00018577, 0.000383531, 0.000290172, 0.00025551, 0.000234736, 0.000260275, 0.000293701, 4.18954e-05, 0.000758746
    };

    private static final double[][] DINUCLEOTIDE_PROB = new double[][] {
        {0.274652, 0.217093, 0.290202, 0.218054},
        {0.206493, 0.387621, 0.196134, 0.209753},
        {0.18501, 0.215686, 0.41834, 0.180964},
        {0.196224, 0.260825, 0.258822, 0.284129}
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
        int first = ntIndex(seq.charAt(length - 1));
        logProb += Math.log(DINUCLEOTIDE_PROB[ntIndex(nucleotide)][first]);
        for (int i = length - 1; i > 0; i--) {
            int from = ntIndex(seq.charAt(i));
            int to = ntIndex(seq.charAt(i - 1));
            logProb += Math.log(DINUCLEOTIDE_PROB[from][to]);
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
}
