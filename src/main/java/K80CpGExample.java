import com.yongkangl.phylogenetics.model.MarkovChain;
import com.yongkangl.phylogenetics.substitution.*;
import com.yongkangl.phylogenetics.util.Taxon;
import com.yongkangl.phylogenetics.util.Utils;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class K80CpGExample {
    private static final char[] nucleotides = {'A', 'C', 'G', 'T'};
    private static final Map<Character, Integer> nucleotideDict = new HashMap<>();

    static {
        for (int i = 0; i < nucleotides.length; i++) {
            nucleotideDict.put(nucleotides[i], i);
        }
    }

    public static Taxon sampleFromStationary(double lambda_CG, double[] pi, double[] proposal, int sequenceLength, int numSamples) {
        String[] samples = generateSequences(numSamples, sequenceLength, proposal);
        double[] logWeights = new double[numSamples];

        for (int i = 0; i < numSamples; i++) {
            logWeights[i] = logP(samples[i], lambda_CG, pi) - logQ(samples[i], proposal);
        }

        double logWeightsSum = Utils.logSum(logWeights);
        double[] normalizedWeights = new double[numSamples];
        for (int i = 0; i < numSamples; i++) {
            normalizedWeights[i] = Math.exp(logWeights[i] - logWeightsSum);
        }

        int weightedSampleIndex = sampleIndex(normalizedWeights);
        String sample = samples[weightedSampleIndex];

        int count = 0;
        for (int j = 0; j < sample.length() - 1; j++) {
            if (sample.charAt(j) == 'C' && sample.charAt(j + 1) == 'G') {
                count++;
            }
        }
        return new Taxon(sample);
    }

    private static double logP(String x, double lambda_CG, double[] pi) {
        double productTerm = 0.0;
        int n = x.length();
        for (int j = 0; j < n - 2; j++) {
            productTerm += (x.charAt(j) == 'C' && x.charAt(j + 1) == 'G' ? Math.log(lambda_CG) : 0) + Math.log(pi[nucleotideDict.get(x.charAt(j + 1))]);
        }
        productTerm += (x.charAt(n - 2) == 'C' && x.charAt(n - 1) == 'G' ? Math.log(lambda_CG) : 0);
        return productTerm;
    }

    private static double logQ(String x, double[] proposal) {
        double productTerm = 0.0;
        int n = x.length();
        for (int j = 1; j < n; j++) {
            productTerm += Math.log(proposal[nucleotideDict.get(x.charAt(j))]);
        }
        return productTerm;
    }

    private static String[] generateSequences(int numSamples, int length, double[] probabilities) {
        Random random = new Random();
        String[] sequences = new String[numSamples];
        for (int i = 0; i < numSamples; i++) {
            StringBuilder sb = new StringBuilder(length);
            for (int j = 0; j < length; j++) {
                double p = random.nextDouble();
                double cumulativeProbability = 0.0;
                for (int k = 0; k < probabilities.length; k++) {
                    cumulativeProbability += probabilities[k];
                    if (p <= cumulativeProbability) {
                        sb.append(nucleotides[k]);
                        break;
                    }
                }
            }
            sequences[i] = sb.toString();
        }
        return sequences;
    }

    private static int sampleIndex(double[] probabilities) {
        Random random = new Random();
        double p = random.nextDouble();
        double cumulativeProbability = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            cumulativeProbability += probabilities[i];
            if (p <= cumulativeProbability) {
                return i;
            }
        }
        return probabilities.length - 1; // Fallback
    }

    public static void main(String[] args) {
        boolean toCorrect = true;
        int proposalChoice = 0;
        double[][] Q = {{-0.2, 0.05, 0.1, 0.05}, {0.05, -0.2, 0.05, 0.1}, {0.1, 0.05, -0.2, 0.05}, {0.05, 0.1, 0.05, -0.2}};
        double trueLambda = 0.15;
        double trueAlpha = 0.4;
        double trueBeta = 0.2;
        double[] pi = {0.25, 0.25, 0.25, 0.25};
        double[] pr = {0.4, 0.1, 0.1, 0.4};
        double T = 1.0;

        IndependentSiteModel proposal = new IndependentSiteModel(Q);
        K80CpGSiteModel trueTarget = new K80CpGSiteModel(trueLambda, trueAlpha, trueBeta);

        int sequenceLength = 750;
        int numPairs = 3;
        Taxon[] initials = new Taxon[numPairs];
        Taxon[] terminals = new Taxon[numPairs];
        int[] xCpGs = new int[numPairs];
        for (int i = 0; i < numPairs; i++) {
            initials[i] = sampleFromStationary(trueLambda, pi, pr, sequenceLength + 2 * trueTarget.getContextLength(), 10000);
            SubstitutionModel trueSubstitutionModel = new SubstitutionModel(trueTarget, initials[i]);
            terminals[i] = trueSubstitutionModel.simulateDependently(T);
            for (int j = 0; j < initials[i].len() - 1; j++) {
                if (initials[i].get(j) == Nucleotides.C_STATE && initials[i].get(j+1) == Nucleotides.G_STATE) {
                    xCpGs[i]++;
                }
            }
            System.out.println(xCpGs[i]);
        }
        int xCpG_sum = Arrays.stream(xCpGs).sum();

        int nIterations = 16;
        int[] sampleSizes = {10, 10, 10, 10, 50, 50, 50, 50, 200, 200, 200, 200, 500, 500, 500, 500};

        double lambda = 0.08;
        double alpha = 0.32;
        double beta = 0.16;
        for (int iteration = 0; iteration < nIterations; iteration++) {
            K80CpGSiteModel target = new K80CpGSiteModel(lambda, alpha, beta);
            int nParticles = sampleSizes[iteration];
            int nSteps = sequenceLength * 256;
            MetropolisHastingsPath[] paths = new MetropolisHastingsPath[nParticles];
            double[] n_ts = new double[numPairs];
            double[] n_tv = new double[numPairs];
            double[] n_CpG = new double[numPairs];
            double[] count = new double[numPairs];
            for (int i = 0; i < numPairs; i++) {
                for (int j = 0; j < nParticles; j++) {
                    paths[j] = new MetropolisHastingsPath(toCorrect, proposalChoice, initials[i], terminals[i], T, proposal, target);
                }
                MarkovChain<MetropolisHastingsPath> metropolisHastings = new MarkovChain<>(paths);
                metropolisHastings.run(nSteps);

                // E-step
                MetropolisHastingsPath.K80CpGSufficientStatistics[] statistics = new MetropolisHastingsPath.K80CpGSufficientStatistics[nParticles];
                for (int j = 0; j < nParticles; j++) {
                    statistics[j] = paths[j].getK80CpGSufficnetStatistics();
                }
                MetropolisHastingsPath.K80CpGSufficientStatistics sampleMean = calculateSampleMeans(statistics);
                n_ts[i] = sampleMean.getNts();
                n_tv[i] = sampleMean.getNtv();
                n_CpG[i] = sampleMean.getNCpG();
                count[i] = sampleMean.getCount();
            }

            double n_CpG_sum = Utils.sum(n_CpG);
            double n_ts_sum = Utils.sum(n_ts);
            double n_tv_sum = Utils.sum(n_tv);
            double count_sum = Utils.sum(count);

            // M-step
            // numerically update lambda
            lambda = estimateLambda(numPairs, lambda, alpha, beta, sequenceLength, xCpG_sum, count_sum, n_CpG_sum, 0.0001);

            // analytically update alpha and beta
            double denominator = 2 * n_CpG_sum / lambda + numPairs * sequenceLength * T - 2 * n_CpG_sum;
            alpha = 4 * n_ts_sum / denominator;
            beta = 4 * n_tv_sum / (2 * denominator);

            System.out.println("Iteration " + iteration + ": lambda=" + lambda + ", alpha=" + alpha + ", beta=" + beta);
        }
    }

    public static MetropolisHastingsPath.K80CpGSufficientStatistics calculateSampleMeans(MetropolisHastingsPath.K80CpGSufficientStatistics[] statistics) {
        int size = statistics.length;
        double meanNts = 0.0;
        double meanNtv = 0.0;
        double meanNCpG = 0.0;
        double meanCount = 0.0;
        for (MetropolisHastingsPath.K80CpGSufficientStatistics statistic : statistics) {
            meanNts += ((double) statistic.getNts());
            meanNtv += ((double) statistic.getNtv());
            meanNCpG += statistic.getNCpG();
            meanCount += ((double) statistic.getCount());
        }
        meanNts /= size;
        meanNtv /= size;
        meanNCpG /= size;
        meanCount /= size;
        return new MetropolisHastingsPath.K80CpGSufficientStatistics(meanNts, meanNtv, meanNCpG, meanCount);
    }

    public static double estimateLambda(int numPairs, double lambda, double alpha, double beta, int sequenceLen, int xCpG, double count, double meanNCpG, double step) {
        double h = 0.000000001;
        double derivative = (logL(numPairs, lambda + h, alpha, beta, sequenceLen, xCpG, count, meanNCpG) - logL(numPairs, lambda, alpha, beta, sequenceLen, xCpG, count, meanNCpG)) / h;
        return lambda + step * derivative;
    }

    public static double logL(int numPairs, double lambda, double alpha, double beta, int sequenceLen, int xCpG, double count, double meanNCpG) {
        return - numPairs * logZ(lambda, sequenceLen) + (xCpG - count) * Math.log(lambda) - (0.5 * alpha + beta) * meanNCpG / lambda;
    }

    public static double logZ(double lambda, int sequenceLen) {
        double[][] A = new double[4][4];

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                A[i][j] = 0.25;
            }
        }
        A[Nucleotides.C_STATE][Nucleotides.G_STATE] *= lambda;

        RealMatrix matrix = MatrixUtils.createRealMatrix(A);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);
        double[] realEigenvalues = eigenDecomposition.getRealEigenvalues();
        int maxID = 0;
        double maxEigenVal = realEigenvalues[0];
        for (int i = 1; i < 4; i++) {
            if (realEigenvalues[i] > maxEigenVal) {
                maxID = i;
                maxEigenVal = realEigenvalues[i];
            }
        }

        RealMatrix leftEigenvectorsMatrix = MatrixUtils.inverse(eigenDecomposition.getVT());
        RealVector vector = leftEigenvectorsMatrix.getColumnVector(maxID);
        double norm = vector.getNorm();
        double l1Sum = 0.0;
        for (int i = 0; i < vector.getDimension(); i++) {
            l1Sum += vector.getEntry(i) / norm;
        }
        l1Sum = Math.abs(l1Sum);

//        double t1 = Math.sqrt(0.75 + 0.25 * lambda);
//        double maxEigenVal = 0.5 + 0.5 * t1;
//        double denominator = 0.875 + t1 + 0.125 * lambda;
//        double l1Sum = 1.0 + 0.25 * (2.5 + 3.0 * t1 + 1.5 * lambda + t1 * lambda) / denominator + 2 * (-0.875 - t1 + 0.75 * lambda + t1 * lambda + 0.125 * lambda * lambda) / (denominator * (-1.0 + lambda));

        return sequenceLen * Math.log(maxEigenVal) + Math.log(l1Sum);
    }
}
