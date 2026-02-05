# Bayesian Canonical Correlation Analysis using DFSM

## Overview
This repository contains R implementations of a Bayesian Canonical Correlation Analysis framework based on the Dynamic Factor Structure Model (DFSM). The project demonstrates applications to breast cancer multi-omics data.

## Contents

- DFSM.R  
  Implementation of Bayesian Dynamic Factor Structure Model

- Bayesian_Summary_Data_Analysis.R  
  Functions for summary data based inference

- Breast_Cancer_Analysis_CCA.Rmd  
  Full analysis workflow

- Breast_Cancer_Analysis_Demonstration.Rmd  
  Reproducible demonstration

## Requirements

- R (>= 4.0)
- Packages:
  - MASS
  - mvtnorm
  - ggplot2
  - dplyr
  - rmarkdown

## Running the Analysis

```r
source("DFSM.R")
