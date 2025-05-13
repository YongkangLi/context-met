import com.yongkangl.phylogenetics.model.MarkovChain;
import com.yongkangl.phylogenetics.model.Particle;
import com.yongkangl.phylogenetics.substitution.*;
import com.yongkangl.phylogenetics.util.Taxon;
import org.apache.commons.cli.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.stream.IntStream;

public class MetropolisHastingsExample {
    public static void main(String[] args) throws IOException {
        Options options = new Options();
        options.addOption("f", "phi", true, "interaction parameter");
        options.addOption("s", "steps", true, "number of steps");
        options.addOption("p", "particles", true, "number of particles");
        options.addOption("c", "correct", true, "to perform metropolis correction or not");
        options.addOption("m", "metropolis", true, "proposal distribution for the metropolis chain");
        options.addOption("x", "initial", true, "initial sequence");
        options.addOption("y", "terminal", true, "terminal sequence");
        options.addOption("t", "T", true, "time interval");
        options.addOption("o", "output", true, "output path");

        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("Error parsing command line: " + e.getMessage());
            System.exit(1);
        }

        double phi = 10.0;
        int nSteps = 16384;
        int nParticles = 8192;

        boolean toCorrect = true;
        int proposalChoice = 1;
        Taxon initial = new Taxon("ACGAA");
        Taxon terminal = new Taxon("AAAAA");
        double T = 0.02;

        String outputPath = "example_output/";

        if (cmd.hasOption("phi")) {
            try {
                phi = Double.parseDouble(cmd.getOptionValue("phi"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid interaction parameter phi");
                System.exit(1);
            }
        }
        if (cmd.hasOption("steps")) {
            try {
                nSteps = Integer.parseInt(cmd.getOptionValue("steps"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number for steps");
                System.exit(1);
            }
        }
        if (cmd.hasOption("particles")) {
            try {
                nParticles = Integer.parseInt(cmd.getOptionValue("particles"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number for particles");
                System.exit(1);
            }
        }
        if (cmd.hasOption("correct")) {
            try {
                toCorrect = Boolean.parseBoolean(cmd.getOptionValue("correct"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid boolean for correction");
                System.exit(1);
            }
        }
        if (cmd.hasOption("m")) {
            try {
                proposalChoice = Integer.parseInt(cmd.getOptionValue("m"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid number for proposal");
                System.exit(1);
            }
        }
        if (cmd.hasOption("x")) {
            try {
                initial = new Taxon(cmd.getOptionValue("x"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid initial sequence");
                System.exit(1);
            }
        }
        if (cmd.hasOption("y")) {
            try {
                terminal = new Taxon(cmd.getOptionValue("y"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid terminal sequence");
                System.exit(1);
            }
        }
        if (cmd.hasOption("t")) {
            try {
                T = Double.parseDouble(cmd.getOptionValue("t"));
            } catch (NumberFormatException e) {
                System.err.println("Invalid time interval");
                System.exit(1);
            }
        }
        if (cmd.hasOption("o")) {
            outputPath = cmd.getOptionValue("o");
        }

        if (!outputPath.endsWith(File.separator)) {
            outputPath += File.separator;
        }

        double[][] Q = {{-1.0, 1.0 / 3, 1.0 / 3, 1.0 / 3}, {1.0 / 3, -1.0, 1.0 / 3, 1.0 / 3}, {1.0 / 3, 1.0 / 3, -1.0, 1.0 / 3}, {1.0 / 3, 1.0 / 3, 1.0 / 3, -1.0}};

        IndependentSiteModel proposal = new IndependentSiteModel(Q);
        SiteModel target = new CpGSiteModel(Q, phi);

        MetropolisHastingsPath[] paths = new MetropolisHastingsPath[nParticles];
        for (int i = 0; i < nParticles; i++) {
            paths[i] = new MetropolisHastingsPath(toCorrect, proposalChoice, initial, terminal, T, proposal, target);
        }
        MarkovChain<MetropolisHastingsPath> metropolisHastings = new MarkovChain<>(paths);
        metropolisHastings.run(nSteps);

        String prefix = toCorrect ? (proposalChoice == 1 ? "mh_hobolth_" : "mh_indpt_") : "uncorrected_";
        FileWriter numberOfMutationsFile = null;
        FileWriter firstJumpTimeFile = null;
        FileWriter gelmanRubinFile = null;
        FileWriter autoCorrelationFile = null;
        numberOfMutationsFile = new FileWriter(outputPath + prefix +"number_of_mutations.txt");
        firstJumpTimeFile = new FileWriter(outputPath + prefix + "first_jump_time.txt");
        gelmanRubinFile = new FileWriter(outputPath + prefix + "gelman_rubin.txt");
        autoCorrelationFile = new FileWriter(outputPath + prefix + "auto_correlation.txt");

        for (int i = 0; i < nParticles; i++) {
            numberOfMutationsFile.write(String.valueOf(paths[i].getNumberOfMutations(nSteps - 1)));
            numberOfMutationsFile.write("\n");
            firstJumpTimeFile.write(String.valueOf(paths[i].getFirstJumpTime(nSteps - 1)));
            firstJumpTimeFile.write("\n");
        }
        numberOfMutationsFile.close();
        firstJumpTimeFile.close();

        double acceptanceRatio = IntStream.range(0, nParticles).mapToDouble(i -> paths[i].getAccepted()).average().getAsDouble();
        acceptanceRatio /= nSteps;
        System.out.println("Acceptance ratio: " + acceptanceRatio);

        if (Taxon.HammingDistance(initial, terminal) == 2) {
            int reverses = 0;
            for (int i = 0; i < nParticles; i++) {
                for (int step = 1; step < nSteps; step++) {
                    if (paths[i].getOrder(step - 1) != paths[i].getOrder(step)) {
                        reverses++;
                    }
                }
            }
            System.out.println("Reverses: " + ((double) reverses) / ((((double) nParticles)) * ((double) nSteps)));
        }

        calculateGelmanRubin(gelmanRubinFile, paths, nSteps, nParticles, 0);
        gelmanRubinFile.close();
        calculateAutoCorrelations(autoCorrelationFile, paths, nSteps, nParticles, 3072, 8192);
        autoCorrelationFile.close();
    }

    public static void calculateGelmanRubin(FileWriter fileWriter, MetropolisHastingsPath[] paths, int nSteps, int nParticles, int burnIn) throws IOException {
        double[] partialSums = IntStream.range(0, nParticles).mapToDouble(chain -> 0).toArray();
        double[] chainMeans = new double[nParticles];

        for (int step = burnIn + 1; step < nSteps; step++) {
            double L = (double) (step - burnIn);
            double chainMeanSum = 0.0;
            for (int chain = 0; chain < nParticles; chain++) {
                partialSums[chain] += paths[chain].getFirstJumpTime(step);
                chainMeans[chain] = partialSums[chain] / L;
                chainMeanSum += chainMeans[chain];
            }
            double grandMean = chainMeanSum / nParticles;
            double squareError = 0.0;
            for (int chain = 0; chain < nParticles; chain++) {
                double error = chainMeans[chain] - grandMean;
                squareError += error * error;
            }
            double B = (L / (nParticles - 1)) * squareError;
            double W = 0.0;
            double[] chainVariances = IntStream.range(0, nParticles).mapToDouble(chain -> 0).toArray();
            for (int chain = 0; chain < nParticles; chain++) {
                for (int i = burnIn + 1; i <= step; i++) {
                    double error = paths[chain].getFirstJumpTime(i) - chainMeans[chain];
                    chainVariances[chain] += error * error;
                }
                chainVariances[chain] /= (L - 1);
                W += chainVariances[chain];
            }
            W /= nParticles;
            double ratio = ((L - 1) / L) + (B / (L * W));
            if (step - burnIn >= 3) {
                fileWriter.write(step + ", " + ratio);
                fileWriter.write("\n");
            }
        }
    }

    public static void calculateAutoCorrelations(FileWriter fileWriter, MetropolisHastingsPath[] paths, int nSteps, int nParticles, int endLag, int burnIn) throws IOException {
        double[] averageAutoCorrelations = new double[endLag / 3];
//        for (int lag = 0; lag < endLag; lag++) {
//            averageAutoCorrelations[lag] = 0.0;
//        }
//        for (MetropolisHastingsPath path: paths) {
//            double[] autoCovariances = path.calculateAutoCovariance(nSteps, endLag, burnIn);
//            for (int lag = 1; lag < endLag; lag++) {
//                averageAutoCorrelations[lag] += autoCovariances[lag] / autoCovariances[0];
//            }
//        }
//        for (int lag = 1; lag < endLag; lag++) {
//            averageAutoCorrelations[lag] /= nParticles;
//        }
//
//        for (int lag = 3; lag < endLag; lag += 3) {
//            System.out.println(lag / 3 + ", " + averageAutoCorrelations[lag]);
//        }

        for (int j = 0; j < endLag / 3; j++) {
            averageAutoCorrelations[j] = 0.0;
        }

        for (int i = 0; i < nParticles; i++) {
            for (int j = 0; j < endLag / 3; j++) {
                averageAutoCorrelations[j] = paths[i].calculateAutoCorrelation(nSteps, (j + 1) * 3, burnIn);
            }
        }

        for (int j = 1; j < endLag/ 3; j++) {
            fileWriter.write((j + 1) + ", " + Math.abs(averageAutoCorrelations[j] / nParticles));
            fileWriter.write("\n");
        }
    }
}
