# Generalized Root Causal Inference

This repository contains code for the Generalized Root Causal Inference (GRCI) algorithm for identify patient/sample-specific root causes of a binary variable using the heteroscedastic noise model (HNM). HNM models both the conditional expectation and conditional mean absolute deviation (similar to the conditional variance) using non-linear functions.

The Experiments folder contains code needed to replicate the synthetic data results.

# Installation

> library(devtools)

> install_github("ericstrobl/GRCI")

> library(GRCI)

# Run GRCI on the Data

> G = generate_DAG_HNM(p=10,en=2)

> X = sample_DAG_HNM(nsamps = 1000,G)

> out = GRCI(X$data[,-G$Y],X$data[,G$Y])

> print(out$order); print(out$scores[1:5,]) # print reverse partial order and their corresponding error terms


