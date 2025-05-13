# On Gibbs Sampling for Endpoint-Conditioned Neighbor-Dependent Sequence Evolution Models

## Description

This repository contains the source code accompanying our paper *On Gibbs Sampling for Endpoint-Conditioned Neighbor-Dependent Sequence Evolution Models* [[1]](#1).
The code implements corrected Markov chain Monte Carlo (MCMC) methods for sampling paths in models of DNA sequence evolution with neighbor-dependent substitution rates.
In particular, it includes a Metropolis-Hastings correction to the previously proposed algorithm by Hobolth (2008) [[2]](#2), which we show does not target the correct full conditional distribution.

In our paper [[1]](#1), we analyze the target distribution of the Markov chain in [[2]](#2) theoretically, demonstrate its bias empirically, and introduce two alternative Metropolis-Hastings samplers—one using the Hobolth algorithm as a proposal and another based on an independent-site model. 
We provide empirical comparisons of the two methods in terms of convergence and efficiency. Our results show that the independent-site proposal leads to better mixing.
This repository enables full reproduction of the main figures (particularly Figure 2 and Figure 3), and is released under the Apache-2.0 license for non-commercial use.

## Instructions

### Dependencies
The code is implemented in Java and plotting scripts are implemented in Python. Therefore, you need to have the following dependencies installed:
- Java 1.8
- Python 3.9 (including `numpy`, `matplotlib`, `seaborn`, `pandas`, `scipy`, and `itertools`)
 
### Repository Structure
    
    context_met/
    ├── src/main/java
    │   ├── MetropolisHastingsExample.java
    │   ├── K80CpGExample.java
    │   └──CpGCounterExample.java
    ├── scripts
    │   ├── generate_true_samples.py
    │   ├── plot_densities.py
    │   ├── plot_histograms.py
    │   ├── plot_auto_correlation.py
    │   └── plot_gelman_rubin.py
    └── executables/
        └── metropolis-hastings.jar

### Arguments

| shorthand | full form    | comments              | example/default value |
|-----------|--------------|-----------------------|-----------------------|
| -f        | --phi        | interaction parameter | 10                    |
| -s        | --steps      | number of steps       | 16384                 |
| -p        | --particles  | number of particles   | 8192                  |
| -c        | --correction | whether to correction | true                  |
| -m        | --metropolis | metropolis proposal   | 1                     |
| -x        | --initial    | initial sequence      | ACGAA                 |
| -y        | --terminal   | terminal sequence     | AAAAA                 |
| -t        | --T          | time interval         | 0.02                  |
| -o        | --output     | output path           | example_output        |

By running the executable file, the acceptance ratio and frequency of reverses will appear in the command line output while all the other outputs will be saved in the output path specified in the command line arguments (`example_output` as the default value). The output files will contain the following information:
- `*_first_jump_time.txt`: The time of the first jump in each sample.
- `*_number_of_mutations.txt`: The number of mutations in each sample.
- `*_auto_correlation.txt`: The auto-correlation of the first jump site.
- `*_gelman_rubin.txt`: Gelman-Rubin diagnostic statistic of the first jump time.

### Reproduction of Figures

1. Produce results from the Hobolth algorithm (The name of the corresponding output files starts with `uncorrected_`)
   
   `java -jar executables/metropolis-hastings.jar -f 10 -s 16384 -p 8192 -c false -m 1 -x ACGAA -y AAAAA -t 0.02 -o example_output`

2. Produce results from metropolis algorithm with the Hobolth algorithm as a proposal (The name of the corresponding output files starts with `mh_hobolth_`)

   `java -jar executables/metropolis-hastings.jar -f 10 -s 16384 -p 8192 -c true -m 1 -x ACGAA -y AAAAA -t 0.02 -o example_output`

3. Produce results from metropolis algorithm with the independent model as a proposal (The name of the corresponding output files starts with `mh_indpt_`)

   `java -jar executables/metropolis-hastings.jar -f 10 -s 16384 -p 8192 -c true -m 0 -x ACGAA -y AAAAA -t 0.02 -o example_output`

4. Generate samples from the exact distribution

   `python scripts/generate_true_samples.py`

5. Plot the densities of the first jump time (`example_output/first_jump_time.pdf`)

   `python scripts/plot_densities.py`
 
6. Plot the histograms of the number of mutations (`example_output/number_of_mutations.pdf`)

   `python scripts/plot_histograms.py`

7. Plot the auto-correlation of the first jump site (`example_output/auto_correlation.pdf`)

    `python scripts/plot_auto_correlation.py`
 
8. Plot the Gelman-Rubin diagnostic statistic of the first jump time (`example_output/gelman_rubin.pdf`)

   `python scripts/plot_gelman_rubin.py`

## References
<a id="1">[1]</a> *Li, Y., Mathews, J., and Schmidler, S. C.* (2025). **On Gibbs Sampling for Endpoint-Conditioned Neighbor-Dependent Sequence Evolution Models**. Journal of Computational and Graphical Statistics.

<a id="2">[2]</a> *Hobolth, A.* (2008). **A Markov chain Monte Carlo Expectation Maximization algorithm for statistical analysis of DNA sequence evolution with neighbor-dependent substitution rates**. Journal of Computational and Graphical Statistics 17, 138–162. [![DOI:](https://zenodo.org/badge/DOI/10.1198/106186008X289010.svg)](https://doi.org/10.1198/106186008X289010)

> ⚠️ **Warning: Research Prototype**  
> This repository contains research code and is provided for academic purposes only.  
> It is not intended or guaranteed to be production-ready. Use at your own risk.
