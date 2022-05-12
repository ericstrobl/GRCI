# Generalized Root Causal Inference

This repository contains code for the Generalized Root Causal Inference (GRCI) algorithm for identify patient/sample-specific root causes of a binary variable using the heteroscedastic noise model (HNM). HNM models both the conditional expectation and conditional mean absolute deviation (similar to conditional variance) using non-linear functions.

The Experiments folder contains code needed to replicate the synthetic data results.

# Installation

> library(devtools)

> install_github("ericstrobl/GRCI")

> library(GRCI)

# Run GRCI on the Data

> suffStat = list(); suffStat$data = synth_data[,resort_p];

> out = CIM(suffStat, GCM_wrap, alpha=0.01, p=ncol(suffStat$data), waves=waves) # run CIM

> print(out$f_star) # print F*


