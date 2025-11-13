package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Taxon;
import jeigen.DenseMatrix;

import java.util.Arrays;

public class BlockwiseSiteModel {
    private final int blockSize;
    private final double[][] Q;

    public BlockwiseSiteModel(double[][] Q) {
//        setContextLength(0);
        this.blockSize = Integer.numberOfTrailingZeros(Q.length) / 2;
        this.Q = Q;
    }

    public int getBlockSize() {
        return blockSize;
    }

    public double getRate(int from, int to) {
        return Q[from][to];
    }

    public double getMutability(int from) {
        return -Q[from][from];
    }

    public static int encode (int[] states) {
        int encoded = 0;
        for (int state : states) {
            encoded = encoded << 2;
            encoded += state;
        }
        return encoded;
    }

    public static int[] decode (int blockSize, int encoded) {
        int[] states = new int[blockSize];
        for (int i = blockSize - 1; i >= 0; i--) {
            states[i] = encoded & 3;
            encoded = encoded >> 2;
        }
        return states;
    }

    public static double[][] getKernel(int subLength, int[] leftNeighbors, int[] rightNeighbors, SiteModel siteModel) {
        int contextLength = leftNeighbors.length;
        int n = 1 << (2 * subLength);
        double[][] infinitesimal = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                infinitesimal[i][j] = 0.0;
            }
        }
        for (int i = 0; i < n; i++) {
            int[] states = decode(subLength, i);
            int[] statesWithNeighbors = new int[subLength + 2 * contextLength];
            System.arraycopy(leftNeighbors, 0, statesWithNeighbors, 0, contextLength);
            System.arraycopy(states, 0, statesWithNeighbors, contextLength, subLength);
            System.arraycopy(rightNeighbors, 0, statesWithNeighbors, contextLength + subLength, contextLength);
            Taxon taxon = new Taxon(statesWithNeighbors);
            SubstitutionModel substitutionModel = new SubstitutionModel(siteModel, taxon);
            for (int site = 0; site < subLength; site++) {
                int originalState = states[site];
                for (int state = 0; state < 4; state++) {
                    if (state != originalState) {
                        double rate = substitutionModel.getRate(contextLength + site, state);
                        states[site] = state;
                        int j = encode(states);
                        infinitesimal[i][j] = rate;
                    }
                }
                states[site] = originalState;
            }
            infinitesimal[i][i] = -Arrays.stream(infinitesimal[i]).sum();
        }
        return infinitesimal;
    }

    public double logLikelihood(int from, int to, double T) {
        DenseMatrix kernel = new DenseMatrix(this.Q).mul(T);
        DenseMatrix transitionMatrix = kernel.mexp();
        return Math.log(transitionMatrix.get(from, to));
    }
}
