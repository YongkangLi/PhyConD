package com.yongkangl.phylogenetics.substitution;

import java.util.Arrays;

/**
 * @author Yongkang Li
 * @version $Id$
 */
public class SiteModel {
    private int contextLength;
    private double[][] contextIndependentRates;
    private double[] contextDependentRates;
    private double[] mutabilities;

    public static SiteModel temper (SiteModel proposal, SiteModel target, double temperature) {
        int contextLength = target.getContextLength();
        double[][] independentRates = new double[4][4];
        for (int from = 0; from < 4; from++) {
            double sum = 0.0;
            for (int to = 0; to < 4; to++) {
                independentRates[from][to] = Math.exp((1 - temperature) * Math.log(proposal.getContextIndependentRate(from, to)) + temperature * Math.log(target.getContextIndependentRate(from, to)));
                if (from != to) {
                    sum += independentRates[from][to];
                }
            }
            independentRates[from][from] = -sum;
        }
        int neighborPossibilities = 1 << (2 * contextLength);
        int totalPossibilities = 1 << (4 * contextLength + 4);
        double[] dependentRates = new double[totalPossibilities];
        for (int leftNeighbor = 0; leftNeighbor < neighborPossibilities; leftNeighbor++) {
            for (int rightNeighbor = 0; rightNeighbor < neighborPossibilities; rightNeighbor++) {
                for (int from = 0; from < 4; from++) {
                    int encoded = (leftNeighbor << (2 * contextLength + 2)) | (from << (2 * contextLength)) | (rightNeighbor);
                    double sum = 0.0;
                    for (int to = 0; to < 4; to++) {
                        int full = (encoded << 2) | to;
                        int adjusted = (from << 2) | to;
                        dependentRates[full] = Math.exp((1 - temperature) * Math.log(proposal.getContextDependentRate(adjusted)) + temperature * Math.log(target.getContextDependentRate(full)));
                        if (from != to) {
                            sum += dependentRates[full];
                        }
                    }
                    dependentRates[(encoded << 2) | from] = -sum;
                }
            }
        }
        SiteModel model = new SiteModel();
        model.setContextLength(contextLength);
        model.setContextIndependentRates(independentRates);
        model.setContextDependentRates(dependentRates);
        model.setMutabilities();
        return model;
    }

    public static SiteModel[] temper(SiteModel proposal, SiteModel target, double[] temperatures) {
        int n = temperatures.length;
        SiteModel[] models = new SiteModel[n+2];
        models[0] = proposal;
        for (int i = 1; i <= n; i++) {
            models[i] = temper(proposal, target, temperatures[i-1]);
        }
        models[n+1] = target;
        return models;
    }

    public static SiteModel[] temper(SiteModel proposal, SiteModel target, int n) {
        double[] temperatures = new double[n-1];
        for (int i = 0; i < n-1; i++) {
//            temperatures[i] = Math.pow(((double) (i+1) / ((double) n)), 1.25);
//            temperatures[i] = 1.0 / (1.0 + Math.exp(- 5.0 * (((double) (i+1) / ((double) n)) - 0.5)));
//            temperatures[i] = 0.5 * (1.0 + Math.sin(Math.PI * (((double) (i+1) / ((double) n)) - 0.5)));
            double k = 1.0;
            double x = (((double) (i+1)) / ((double) n));
            temperatures[i] = 0.5 + 0.5 * Math.atan(k * (x - 0.5)) / Math.atan(k / 2.0);
        }
        return temper(proposal, target, temperatures);
    }

    public void setContextLength(int contextLength) {
        this.contextLength = contextLength;
    }

    public void setContextIndependentRates(double[][] contextIndependentRates) {
        this.contextIndependentRates = contextIndependentRates;
    }

    public void setContextDependentRates(double[] contextDependentRates) {
        this.contextDependentRates = contextDependentRates;
    }

    public void setMutabilities() {
        int len = (int) (Math.pow(4, 2 * contextLength + 1));
        mutabilities = new double[len];
        for (int encoded = 0; encoded < len; encoded++) {
            mutabilities[encoded] = 0.0;
            int from = (encoded >> (2 * contextLength)) & 0x3;
            for (int state = 0; state < 4; state++) {
                if (state != from) {
                    mutabilities[encoded] += getContextDependentRate((encoded << 2) | state) * getContextIndependentRate(from, state);
                }
            }
        }
    }

    public int getContextLength() {
        return contextLength;
    }

    public double getContextIndependentRate(int from, int to) {
        return contextIndependentRates[from][to];
    }

    // the encode is left-context | site | right-context | to
    public double getContextDependentRate(int encoded) {
        return contextDependentRates[encoded];
    }

    public double getRate(int encoded) {
        int from = (encoded >> (2 * getContextLength() + 2)) & 0x3;
        int to = encoded & 0x3;
        return getContextDependentRate(encoded) * getContextIndependentRate(from, to);
    }

    // the encode is left-context | site | right-context
    public double getMutability(int encoded) {
        return mutabilities[encoded];
    }

    public void scaleContextDependentRates(double scale) {
        for (int i = 0; i < contextDependentRates.length; i++) {
            contextDependentRates[i] *= scale;
        }
        for (int i = 0; i < mutabilities.length; i++) {
            mutabilities[i] *= scale;
        }
    }
}
