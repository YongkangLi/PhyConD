package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Taxon;
import jeigen.DenseMatrix;

public class IndependentLogLikelihood {
    private final DenseMatrix kernel;

    public IndependentLogLikelihood(double[][] Q) {
        kernel = new DenseMatrix(Q);
    }

    public DenseMatrix setTime(double T) {
        return kernel.mul(T).mexp();
    }

    public double logLikelihood(Taxon initial, Taxon terminal, double T, int contextLength) {
        DenseMatrix P = setTime(T);
        double ll = 0;
        for (int i = contextLength; i < initial.len() - contextLength; i++) {
            ll += Math.log(P.get(initial.get(i), terminal.get(i)));
        }
        return ll;
    }
}
