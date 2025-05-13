import com.yongkangl.phylogenetics.model.PairwiseSequenceSMC;
import com.yongkangl.phylogenetics.substitution.*;
import com.yongkangl.phylogenetics.util.Taxon;

import java.util.Arrays;

public class CpGCounterExample {
    public static void main(String[] args) {
        double lambda = 0.1;
        double[][] Q = {{-1.0, 1.0 / 3, 1.0 / 3, 1.0 / 3}, {1.0 / 3, -1.0, 1.0 / 3, 1.0 / 3}, {1.0 / 3, 1.0 / 3, -1.0, 1.0 / 3}, {1.0 / 3, 1.0 / 3, 1.0 / 3, -1.0}};
        double[] Ts = {0.3, 1.0, 2.0, 5.0};

        Taxon initial = new Taxon("ACGTA");
        Taxon terminal = new Taxon("ACTTA");

        int nParticles = 2048;
        int nSteps = 16;
        int mutationStep = 16;
        int N = 1;
        for(double T: Ts) {
            double[] transitionProbabilities = new double[N];
            for(int i = 0; i < N; i++) {
                transitionProbabilities[i] = SMC(initial, terminal, T, Q, lambda, nParticles, nSteps, mutationStep);
            }
            String arrayString = Arrays.toString(transitionProbabilities);
            System.out.println(T + ": " + arrayString.substring(1, arrayString.length() - 1));
        }

//        double T = 0.3;
//        double lambda = 0.05;
//        SMC(initial, terminal, T, Q, lambda, 2048, 16, 16);

//        double[] lambda_list = {0.02, 0.05, 0.1, 0.2, 0.5, 0.9, 0.95, 1.05, 1.1, 1.2, 2.0, 5.0};
//        double[] T_list = {0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 2.0, 5.0};
//        int[] m_list = {1, 2, 4, 6, 8, 10, 12, 16, 20, 24, 32, 40, 48, 64, 80, 96, 128, 160, 192, 256};

//        int n = 16;
//        double[] estimations = new double[n];
//        for (int m : m_list) {
//            double mean = 0.0;
//            for (int i = 0; i < n; i++) {
//                estimations[i] = SMC(initial, terminal, 0.3, Q, 0.1, 2048, 16, m);
//                mean += estimations[i];
//            }
//            mean /= n;
//            double std = 0.0;
//            for (int i = 0; i < n; i++) {
//                double error = estimations[i] - mean;
//                std += error * error;
//            }
//            std = Math.sqrt(std / n);
//            System.out.println("m=" + m + ": mean=" + mean + ", std=" + std);
//        }
    }

    public static double SMC(Taxon initial, Taxon terminal, double T, double[][] Q, double lambda, int nParticles, int nSteps, int mutationSteps) {
        IndependentSiteModel proposal = new IndependentSiteModel(Q);
        SiteModel target = new CpGSiteModel(Q, 1.0 / lambda);

        PairwiseSequenceSMC sequenceMonteCarlo = new PairwiseSequenceSMC(initial, terminal, T, 1, Double.POSITIVE_INFINITY, proposal, target, nSteps, nParticles, mutationSteps, 1.0);
        sequenceMonteCarlo.run();

        IndependentLogLikelihood denominator = new IndependentLogLikelihood(Q);

        return Math.exp(sequenceMonteCarlo.getLogWeight() + denominator.logLikelihood(initial, terminal, T, target.getContextLength()));
    }

    public static void IS(Taxon initial, Taxon terminal, double T, double[][] Q, double lambda, int nParticles) {
        IndependentSiteModel proposal = new IndependentSiteModel(Q);
        SiteModel target = new CpGSiteModel(Q, 1.0 / lambda);

        PairwiseSequenceSMC sequenceMonteCarlo = new PairwiseSequenceSMC(initial, terminal, T, 1, Double.POSITIVE_INFINITY, target, 1, nParticles, 0, 1.0);
        sequenceMonteCarlo.run();

        IndependentLogLikelihood denominator = new IndependentLogLikelihood(Q);

        System.out.println("T=" + T + ", lambda=" + lambda + ": " + Math.exp(sequenceMonteCarlo.getLogWeight() + denominator.logLikelihood(initial, terminal, T, target.getContextLength())));
    }
}