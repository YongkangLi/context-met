package com.yongkangl.phylogenetics.substitution;

import jeigen.DenseMatrix;


public class IndependentSiteModel extends SiteModel {
    private final IndependentLogLikelihood logLikelihood;

    public IndependentSiteModel(double[][] Q) {
        setContextLength(0);
        setContextIndependentRates(Q);
        logLikelihood = new IndependentLogLikelihood(Q);
        setMutabilities();
    }

    public int getContextLength() {
        return 0;
    }

    @Override
    public double getContextDependentRate(int encoded) {
        return 1.0;
    }

    public double compute(int from, int to) {
        return getRate(((from << 2) | to));
    }

    public DenseMatrix getTransitionProbability (double T) {
        return logLikelihood.setTime(T);
    }
}
