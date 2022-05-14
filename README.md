# Generalized Root Causal Inference

This repository contains R code for the Generalized Root Causal Inference (GRCI) algorithm for identifying patient/sample-specific root causes of a binary target using the heteroscedastic noise model (HNM). HNM models both the conditional expectation m(X) and the conditional mean absolute deviation s(X) (similar to the conditional standard deviation) using non-linear functions: Y = m(X) + s(X)E.

The Experiments folder contains code needed to replicate the experimental results.

# Installation

Click the green 'Code' button up top and download the .zip file. Then:

> library(devtools)

Extract 'GRCI-master.zip' and then extract 'RANN-master.zip.' Install RANN:

> install("Directory_to.../RANN-master")

Then install GRCI:

> install("Directory_to.../GRCI-master")

> library(GRCI)

# Run GRCI

> G = generate_DAG_HNM(p=10,en=2)

> X = sample_DAG_HNM(nsamps = 1000,G)

> out = GRCI(X$data[,-G$Y],X$data[,G$Y])

> print(out$order); print(out$scores[1:5,]) # print reverse partial order and the corresponding error terms


