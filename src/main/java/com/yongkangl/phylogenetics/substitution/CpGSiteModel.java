package com.yongkangl.phylogenetics.substitution;

import static org.apache.commons.math3.util.ArithmeticUtils.pow;

public class CpGSiteModel extends SiteModel {
    public int getContextLength() {
        return 1;
    }

    public CpGSiteModel(double[][] Q, double phi) {
        setContextLength(1);
        double[][] contextIndependentRates = new double[4][4];
        for (int i = 0; i < 4; i++) {
            System.arraycopy(Q[i], 0, contextIndependentRates[i], 0, 4);
        }
        int len = pow(4, 4);
        double[] contextDependentRates = new double[len];
        for (int encoded = 0; encoded < len; encoded++) {
            int left = encoded >> 6;
            int from = (encoded >> 4) & 0x3;
            int right = (encoded >> 2) & 0x3;
            if (left == Nucleotides.C_STATE && from == Nucleotides.G_STATE) {
                contextDependentRates[encoded] = phi;
            } else if (from == Nucleotides.C_STATE && right == Nucleotides.G_STATE) {
                contextDependentRates[encoded] = phi;
            } else {
                contextDependentRates[encoded] = 1.0;
            }
        }
        setContextDependentRates(contextDependentRates);
        setContextIndependentRates(contextIndependentRates);
        setMutabilities();
    }
}
