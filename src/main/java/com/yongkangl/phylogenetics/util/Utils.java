package com.yongkangl.phylogenetics.util;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.ToDoubleFunction;

public class Utils {
    public static double logAdd(double logA, double logB) {
        if (logA == Double.NEGATIVE_INFINITY && logB == Double.NEGATIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        } else if (logA < logB) {
                return logB + Math.log(1 + Math.exp(logA - logB));
        } else {
            return logA + Math.log(1 + Math.exp(logB - logA));
        }
    }

    private static double logSum(double[] logValues, int start, int end) {
        double ret;
        if (end - start == 1) {
            ret = logValues[start];
        } else if (end - start == 2) {
            ret = logAdd(logValues[start], logValues[start + 1]);
        } else {
            int mid = (start + end) / 2;
            ret = logAdd(logSum(logValues, start, mid), logSum(logValues, mid, end));
        }
        return ret;
    }

    public static double logSum(double[] logValues) {
        double ret;
        if (logValues.length == 0) {
            ret = Double.NEGATIVE_INFINITY;
        } else if (logValues.length == 1) {
            ret = logValues[0];
        } else {
            ret = logSum(logValues, 0, logValues.length);
        }
        return ret;
    }

    public static double logSquareSum(double[] logValues) {
        double newLogValues[] = new double[logValues.length];
        for (int i = 0; i < logValues.length; i++) {
            newLogValues[i] = 2 * logValues[i];
        }
        return logSum(newLogValues);
    }

    public static double ESS(double[] logWeights) {
        return Math.exp(2 * Utils.logSum(logWeights) - Utils.logSquareSum(logWeights));
    }

    private static double min(double[] values, int start, int end) {
        if (end - start == 1) {
            return values[start];
        } else if (end - start == 2) {
            return Math.min(values[start], values[start + 1]);
        } else {
            int mid = (start + end) / 2;
            return Math.min(min(values, start, mid), min(values, mid, end));
        }
    }
    public static double min(double[] values) {
        return min(values, 0, values.length);
    }

    private static double sum(double[] values, int start, int end) {
        if (end - start == 1) {
            return values[start];
        } else if (end - start == 2) {
            return values[start] + values[start + 1];
        } else {
            int mid = (start + end) / 2;
            return sum(values, start, mid) + sum(values, mid, end);
        }
    }
    public static double sum(double[] values) {
        return sum(values, 0, values.length);
    }
}
