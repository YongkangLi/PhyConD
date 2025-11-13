package com.yongkangl.phylogenetics.util;

import com.yongkangl.phylogenetics.substitution.BlockwiseSiteModel;
import com.yongkangl.phylogenetics.substitution.Nucleotides;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;

public class Taxon {
    private final int length;
    private final int[] sequence;
    private String name;

    public Taxon(String string) {
        length = string.length();
        sequence = new int[length];
        for (int i = 0; i < length; i++) {
            switch (string.charAt(i)) {
                case 'A':
                    sequence[i] = Nucleotides.A_STATE;
                    break;
                case 'C':
                    sequence[i] = Nucleotides.C_STATE;
                    break;
                case 'G':
                    sequence[i] = Nucleotides.G_STATE;
                    break;
                case 'U':
                case 'T':
                    sequence[i] = Nucleotides.UT_STATE;
                    break;
                case '-':
                    sequence[i] = Nucleotides.GAP_STATE;
                    break;
                default:
                    sequence[i] = Nucleotides.UNKNOWN_STATE;
            }
        }
    }

    public Taxon(Taxon copyFrom) {
        this.name = copyFrom.name;
        this.length = copyFrom.len();
        this.sequence = new int[length];
        copySequence(copyFrom);
    }

    public Taxon(int length) {
        this.length = length;
        sequence = new int[length];
    }

    public Taxon(int[] sequence) {
        this.length = sequence.length;
        this.sequence = sequence;
    }

    public Taxon(int[] sequence, int length) {
        this.length = length;
        this.sequence = sequence;
    }

    public Taxon(double[] partials) {
        this.length = partials.length >> 2;
        this.sequence = new int[this.length];
        for (int i = 0; i < this.length; i++) {
            List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
            for (int j = 0; j < 4; j++) {
                possibleStates.add(new Pair<>(j, partials[i * 4 + j]));
            }
            EnumeratedDistribution<Integer> distribution = new EnumeratedDistribution<>(possibleStates);
            this.sequence[i] = distribution.sample();
        }
    }

    public void copySequence(Taxon copyFrom) {
        if (length >= 0) {
            System.arraycopy(copyFrom.getSequence(), 0, this.sequence, 0, length);
        }
    }

    public int len() {
        return length;
    }

    public int get(int i) {
        if (i < 0 || i >= length) {
            System.out.println("Invalid index: " + i + " for " + length);
        }
        return sequence[i];
    }

    public int getEncodedBlock(int site, int end) {
        return BlockwiseSiteModel.encode(subSequence(site, end));
    }

    public int[] getSequence() {
        return sequence;
    }

    public int[] subSequence(int start, int end) {
        int[] subSequence = new int[end - start];
        System.arraycopy(sequence, start, subSequence, 0, end - start);
        return subSequence;
    }

    public void set(int i, int state) {
        sequence[i] = state;
    }

    public void set(int offset, int[] states) {
        System.arraycopy(states, 0, sequence, offset, states.length);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < length; i++) {
            switch (sequence[i]) {
                case Nucleotides.A_STATE:
                    sb.append('A');
                    break;
                case Nucleotides.C_STATE:
                    sb.append('C');
                    break;
                case Nucleotides.G_STATE:
                    sb.append('G');
                    break;
                case Nucleotides.UT_STATE:
                    sb.append('T');
                    break;
                case Nucleotides.GAP_STATE:
                    sb.append('-');
                    break;
                default:
                    sb.append('?');
            }
        }
        return sb.toString();
    }

    public static int HammingDistance(Taxon a, Taxon b) {
        int distance = 0;
        for (int i = 0; i < a.len(); i++) {
            if (a.get(i) != b.get(i)) {
                distance++;
            }
        }
        return distance;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
}
