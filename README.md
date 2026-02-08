# Bayesian Canonical Correlation Analysis for Multi-Omics Data

## Overview
This repository provides an implementation of a Bayesian Sparse Canonical Correlation Analysis (CCA) framework for integrative analysis of multi-omics data. The method is designed to identify shared latent structure across multiple high-dimensional biological data modalities.

The repository includes core analysis functions, a reproducible breast cancer multi-omics demonstration, and a rendered analysis report.

---

## Repository Contents

| File | Description |
|------|-------------|
| `Bayesian_Summary_Data_Analysis.R` | Core functions for Bayesian CCA using summary data |
| `Breast_Cancer_Analysis_Demonstration.Rmd` | Reproducible analysis workflow |
| `Breast_Cancer_Analysis_Demonstration.pdf` | Rendered demonstration report |
| `Breast Cancer Data Demonstration.Rdata` | Example dataset used in the analysis |

Several other files contains competing method functions
---


## Method Summary

Canonical Correlation Analysis (CCA) is a classical statistical technique used to identify relationships between two sets of variables. However, traditional CCA is often unstable in:

- High-dimensional settings
- Noisy biological data
- Multi-omics integration problems

The Bayesian CCA framework implemented here:

- Introduces prior-based shrinkage
- Stabilizes estimation in high-dimensional regimes
- Enables inference using summary-level data
- Provides posterior uncertainty quantification

This approach is particularly suited for integrative biomedical and multi-omics studies.

---

## Requirements

### Software
- R (version ≥ 4.0 recommended)

### Required R Packages

Install the required packages before running the code:

```r
install.packages(c(
  "MASS",
  "mvtnorm",
  "ggplot2",
  "dplyr",
  "rmarkdown"
))

