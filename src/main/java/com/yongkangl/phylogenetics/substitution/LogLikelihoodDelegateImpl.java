package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Mutation;
import com.yongkangl.phylogenetics.util.Taxon;

import java.util.List;

public class LogLikelihoodDelegateImpl implements LogLikelihoodDelegate {
    private final SubstitutionModel substitutionModel;

    public LogLikelihoodDelegateImpl(Taxon initial, SiteModel siteModel) {
        this.substitutionModel = new SubstitutionModel(siteModel, initial);
    }

    public double logLikelihoodIndependentSingleSite(int site, List<Mutation> mutations, double T) {
        substitutionModel.reset();
        double ll = 0;
        double time = 0;
        for (Mutation mutation : mutations) {
            ll -= (mutation.getTime() - time) * substitutionModel.getMutability(site);
            ll += Math.log(substitutionModel.mutate(site, mutation.getTo()));
            time = mutation.getTime();
        }
        ll -= (T - time) * substitutionModel.getMutability(site);
        return ll;
    }

    public double logLikelihoodDependentSingleSite(int targetSite, List<Mutation> mutations, double T) {
        double ll = 0.0;
        int size = mutations.size();
        substitutionModel.reset();
        double t = 0.0;
        for (int i = 0; i < size; i++) {
            Mutation mutation = mutations.get(i);
            int currentSite = mutation.getSite();
            ll -= (mutation.getTime() - t) * substitutionModel.getMutability(targetSite);
            double mutationRate = substitutionModel.mutate(currentSite, mutation.getTo());
            if (currentSite == targetSite) {
                ll += Math.log(mutationRate);
            }
            t = mutation.getTime();
        }
        ll -= (T - t) * substitutionModel.getMutability(targetSite);
        return ll;
    }

    public double logLikelihoodFullPath(List<Mutation> mutations, double T, int excluded) {
        // THIS IS VERY IMPORTANT!!!
        substitutionModel.reset();
        double ll = 0;
        double time = 0;
        for (Mutation mutation : mutations) {
            int site = mutation.getSite();
            ll -= (mutation.getTime() - time) * substitutionModel.getMutabilitiesSum(excluded);
            ll += Math.log(substitutionModel.mutate(site, mutation.getTo()));
            time = mutation.getTime();
        }
        ll -= (T - time) * substitutionModel.getMutabilitiesSum(excluded);
        return ll;
    }

    public double logLikelihoodFullPath(List<Mutation> mutations, double T) {
        // THIS IS VERY IMPORTANT!!!
        substitutionModel.reset();
        double ll = 0;
        double time = 0;
        for (Mutation mutation : mutations) {
            int site = mutation.getSite();
            ll -= (mutation.getTime() - time) * substitutionModel.getMutabilitiesSum();
            ll += Math.log(substitutionModel.mutate(site, mutation.getTo()));
            time = mutation.getTime();
        }
        ll -= (T - time) * substitutionModel.getMutabilitiesSum();
        return ll;
    }

    public double logGeometricallyTemperedLikelihood(BlockwiseSubstitutionModel blockwiseSubstitutionModel, List<Mutation> mutations, double T, int targetSite, double temperature) {
        int contextLength = substitutionModel.getContextLength();
        substitutionModel.reset();
        double lr = 0;
        double time = 0;
        for (Mutation mutation : mutations) {
            int currentSite = mutation.getSite();
            double dt = mutation.getTime() - time;
            for (int site = targetSite - contextLength; site <= targetSite + contextLength; site++) {
                if (contextLength <= site && site < substitutionModel.getCurrent().len() - contextLength) {
                    lr -= dt * blockwiseSubstitutionModel.geometricPathMutability(substitutionModel, site, temperature);
                }
            }

            int from = substitutionModel.getCurrent().get(currentSite);
            int to = mutation.getTo();
            double mutationRate = substitutionModel.mutate(currentSite, mutation.getTo());
            if (targetSite - contextLength <= currentSite && currentSite <= targetSite + contextLength) {
                lr += (1 - temperature) * Math.log(blockwiseSubstitutionModel.getRate(currentSite, to)) + temperature * Math.log(mutationRate);
            }
            time = mutation.getTime();
        }
        double dt = T - time;
        for (int site = targetSite - contextLength; site <= targetSite + contextLength; site++) {
            if (contextLength <= site && site < substitutionModel.getCurrent().len() - contextLength) {
                lr -= dt * blockwiseSubstitutionModel.geometricPathMutability(substitutionModel, site, temperature);
            }
        }
        return lr;
    }
}
