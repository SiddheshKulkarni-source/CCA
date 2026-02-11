---

# Bayesian Sparse Canonical Correlation Analysis

## Multiplicative Half-Cauchy Factor Model (MHCFM)

---

## Overview

This repository provides an implementation of the **Multiplicative Half-Cauchy Factor Model (MHCFM)** for sparse Canonical Correlation Analysis (CCA).

The method is developed in:

> Kulkarni, S., Pal, S., & Gaskins, J. (2026).
> *A Bayesian Methodology for Estimation of Sparse Canonical Correlation.*

The MHCFM framework is designed for:

* High-dimensional multi-omics integration
* Sparse cross-view signal detection
* Posterior uncertainty quantification
* Robust inference under p > n regimes

The repository includes:

* Core MHCFM implementation
* Competing Bayesian methods used in simulation studies
* Reproducible breast cancer multi-omics analysis
* Posterior summary and inference tools

---

## Repository Structure

| File                                       | Description                                                        |
| ------------------------------------------ | ------------------------------------------------------------------ |
| `MHCFM.R`                                  | Core implementation of the Multiplicative Half-Cauchy Factor Model |
| `Bayesian_Summary_Data_Analysis.R`         | Posterior summary, credible intervals, shrinkage diagnostics       |
| `Breast_Cancer_Analysis_Demonstration.Rmd` | Reproducible real-data analysis                                    |
| `Breast_Cancer_Analysis_Demonstration.pdf` | Rendered demonstration report                                      |
| `Breast Cancer Data Demonstration.Rdata`   | Processed DNA/RNA dataset used in the manuscript                   |
| `Simulation_*.R`                           | Scripts for generating simulation datasets                         |
| `Competitor_*.R`                           | Code for competing CCA methods used in benchmarking                |

**Note:** Earlier development versions included alternative shrinkage models.
The repository now contains only the finalized MHCFM implementation to ensure consistency with the manuscript.

---

## Method Summary

Classical CCA becomes unstable in:

* High-dimensional settings
* p >> n regimes
* Noisy biological data
* Multi-omics integration problems

The MHCFM addresses these issues by:

* Modeling shared latent structure via Bayesian factor models
* Applying multiplicative Half-Cauchy shrinkage priors
* Encouraging structured sparsity in loading matrices
* Estimating canonical correlations through posterior sampling
* Providing credible intervals for canonical correlations and direction vectors

Unlike frequentist sparse CCA methods, MHCFM:

* Provides full posterior uncertainty
* Incorporates hierarchical shrinkage
* Remains stable in high-dimensional regimes

---

## Computational Considerations

MHCFM uses MCMC sampling.

Typical settings used in the manuscript:

burn_iter <- 25000
mcmc_iter <- 75000
thin <- 15
d <- 20

### Approximate Runtime

Runtime depends on:

* Sample size (n)
* Dimensionality (p₁, p₂)
* Number of factors (d)
* MCMC iterations

Representative runtime ranges:

* Moderate dimensional (n > p): several minutes
* High dimensional (p > n): tens of minutes
* Frequentist competitors: typically seconds

Bayesian methods are computationally intensive due to:

* Latent variable sampling
* Hierarchical shrinkage updates
* Graphical horseshoe updates

Users are encouraged to reduce `mcmc_iter` for testing purposes.

---

## Installation Requirements

### R Version

R ≥ 4.0 recommended

### Required Packages

Install required packages before running:

install.packages(c(
"MASS",
"mvtnorm",
"psych",
"ggplot2",
"dplyr",
"rmarkdown",
"coda"
))

---

## Reproducing the Breast Cancer Analysis

1. Download all files into a working directory.
2. Open `Breast_Cancer_Analysis_Demonstration.Rmd`.
3. Ensure the following files are present:

   * `MHCFM.R`
   * `Bayesian_Summary_Data_Analysis.R`
   * `Breast Cancer Data Demonstration.Rdata`
4. Knit the Rmd file.

If precomputed output is unavailable, run:

MHCFM_Output <- MHCFM(
d = 20,
burn_iter = 25000,
mcmc_iter = 75000,
thin = 15,
X_1 = X_1,
X_2 = X_2,
CCA_select = 10
)

---

## Output Objects

The `MHCFM()` function returns:

| Object                         | Description                                     |
| ------------------------------ | ----------------------------------------------- |
| `A_1_MCMC`                     | Loading matrices for View 1                     |
| `A_2_MCMC`                     | Loading matrices for View 2                     |
| `Mu_1_MCMC`, `Mu_2_MCMC`       | Mean vectors                                    |
| `Omega_1_MCMC`, `Omega_2_MCMC` | Precision matrices                              |
| `CCA_MCMC`                     | Posterior samples of canonical correlations     |
| `Direction_CCA_Vec1_MCMC`      | Posterior samples of direction vectors (View 1) |
| `Direction_CCA_Vec2_MCMC`      | Posterior samples of direction vectors (View 2) |
| `log_det_MCMC`                 | Log-determinant of covariance                   |
| `Log_likelihood_MCMC_Grand`    | Posterior log-likelihood                        |

Posterior summaries are generated using:

Analysis_Summary()

This function provides:

* Posterior means
* Credible intervals
* Effective sample sizes
* Shrinkage diagnostics
* Convergence traceplots

---

## Overshrinkage Diagnostic

The summary function computes:

Shrinkage_calculator

This represents the posterior probability that the first canonical correlation is below a user-specified threshold.

High values may indicate excessive shrinkage.

---

## Simulation Reproducibility

Simulation code used in the manuscript is included in this repository.

To reproduce simulation results:

1. Generate datasets using provided simulation scripts.
2. Fit MHCFM and competitor methods.
3. Summarize posterior results using summary functions.

Note: Full simulation replication requires substantial computation time.

---

## Notes for Reviewers and Users

* All manuscript results can be reproduced using the provided code and data.
* Sensitivity analyses with respect to factor dimension (d) and shrinkage parameter (ζ) are described in the manuscript and Supplementary Materials.
* The repository has been aligned with the final manuscript version.

---

## Contact

Siddhesh Kulkarni
Bayesian Machine Learning Scientist
Translational Informatics & Predictive Sciences

For questions regarding implementation, please open an issue in this repository.

---


