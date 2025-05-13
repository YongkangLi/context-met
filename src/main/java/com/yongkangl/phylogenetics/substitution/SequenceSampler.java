package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Taxon;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.List;

public class SequenceSampler implements Sampler {
    private final int length;
    private final EnumeratedDistribution<Integer>[] enumeratedSiteDistributions;

    public SequenceSampler(double[] partials, int sequenceLength) {
        this.length = sequenceLength;
        enumeratedSiteDistributions = new EnumeratedDistribution[length];
        int index = 0;
        for (int site = 0; site < length; site++) {
            List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
            for (int state = 0; state < 4; state++) {
                possibleStates.add(new Pair<>(state, partials[index++]));
            }
            enumeratedSiteDistributions[site] = new EnumeratedDistribution<>(possibleStates);
        }
    }

    public Taxon sample() {
        int[] states = new int[length];
        for (int site = 0; site < length; site++) {
            states[site] = enumeratedSiteDistributions[site].sample();
        }
        return new Taxon(states);
    }
}
