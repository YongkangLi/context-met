package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Taxon;
import jeigen.DenseMatrix;

public class DependentLogLikelihood {
    public static double logLikelihood(SiteModel siteModel, Taxon initial, Taxon terminal, double T) {
        int len = initial.len();
        int contextLength = siteModel.getContextLength();
        int[] leftNeighbors = initial.subSequence(0, contextLength);
        int[] rightNeighbors = initial.subSequence(len - contextLength, len);
        double[][] kernel = BlockwiseSiteModel.getKernel(len - 2 * contextLength, leftNeighbors, rightNeighbors, siteModel);
        DenseMatrix P = new DenseMatrix(kernel).mul(T).mexp();
        return Math.log(P.get(BlockwiseSiteModel.encode(initial.subSequence(contextLength, len - contextLength)), BlockwiseSiteModel.encode(terminal.subSequence(contextLength, len - contextLength))));
    }
}
