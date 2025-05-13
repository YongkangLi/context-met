package com.yongkangl.phylogenetics.substitution;

public class K80CpGSiteModel extends CpGSiteModel {
    public K80CpGSiteModel(double lambda, double alpha, double beta) {
        super(new double[][]{
                {-0.25 * (alpha + 2 * beta), 0.25 * beta, 0.25 * alpha, 0.25 * beta},
                {0.25 * beta, -0.25 * (alpha + 2 * beta), 0.25 * beta, 0.25 * alpha},
                {0.25 * alpha, 0.25 * beta, -0.25 * (alpha + 2 * beta), 0.25 * beta},
                {0.25 * beta, 0.25 * alpha, 0.25 * beta, -0.25 * (alpha + 2 * beta)}
        }, 1.0 / lambda);
    }
}