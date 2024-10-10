# Variational Inference in Location-Scale Families

This directory contains supplemental material for the paper "Variational Inference in Location-Scale Families: Exact Recovery of the Mean and Correlation Matrix". The code can be used to reproduce the experimental results and figures of the paper.

The VI and MCMC methods are implemented in Stan. All the relevant stan files are in the `model` folder and all the data file are in the `data` files.

There are two R scripts to run the experiments:

*  `examples_illustration.R` generates Figures 1, 2a, 2b, 3, 4, 5, 6 and 9, for Sections 2.3, 3.2 and 4.2, and for Appendix B.2.
* `experiments.R` generates Figure 7 and the results for Table 1 for Section 5.
