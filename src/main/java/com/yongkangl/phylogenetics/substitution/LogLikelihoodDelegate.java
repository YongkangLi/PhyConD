package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Mutation;

import java.util.List;

public interface LogLikelihoodDelegate {
    public double logLikelihoodFullPath(List<Mutation> mutations, double T);
    public double logLikelihoodFullPath(List<Mutation> mutations, double T, int excluded);
}
