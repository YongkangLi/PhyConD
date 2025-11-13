package com.yongkangl.phylogenetics.substitution;


public class IndependentSiteModel extends SiteModel {
    public IndependentSiteModel(double[][] Q) {
        setContextLength(0);
        setContextIndependentRates(Q);
        setMutabilities();
    }

    public int getContextLength() {
        return 0;
    }

    @Override
    public double getContextDependentRate(int encoded) {
        return 1.0;
    }
}
