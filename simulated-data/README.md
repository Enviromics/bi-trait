````markdown
# df-reml

Derivative-free REML implementation for fitting bi-trait genomic--enviromic mixed models.

This folder contains the R scripts used to fit the bi-trait mixed model contrasting **EXP** and **CST**. In the manuscript, these two contexts are referred to as **BTY** and **CTY**, respectively: **EXP corresponds to Breeding Trial Yield (BTY)**, and **CST corresponds to Cultivar Target Yield (CTY)**.

## Overview

The scripts in this folder implement a derivative-free restricted maximum likelihood (DF-REML) approach for estimating variance--covariance components and predicting genotype performance across two correlated contexts.

The model was developed for methodological testing with simulated EXP/CST data and later adapted conceptually to the empirical BTY/CTY structure used in the manuscript.

The main goal is to estimate how performance in one context can inform prediction in the other, while accounting for:

- genomic relationships among genotypes;
- genetic covariance between EXP and CST;
- context-specific residual variances;
- genotype-specific responses to environmental covariates.

## Model context

In the original implementation:

- `EXP` represents experimental trial observations;
- `CST` represents commercial stand observations.

In the manuscript terminology:

- `EXP` is interpreted as `BTY`, Breeding Trial Yield;
- `CST` is interpreted as `CTY`, Cultivar Target Yield.

Thus, the EXP/CST implementation provides the computational basis for the BTY-to-CTY modeling framework described in the article.

## Main script

The main script is:

```text
DF-REML_bi-trait.R
````

This script fits a bi-trait genomic--enviromic mixed model using a derivative-free optimizer.

## Input files

The script expects two main input files:

```text
phenos_and_ECs.txt
M.txt
```

where:

* `phenos_and_ECs.txt` contains phenotypic records, genotype identifiers, context labels (`EXP` or `CST`), locations, coordinates, response values, and environmental covariates;
* `M.txt` contains the SNP marker matrix used to compute the genomic relationship matrix.

The expected phenotype file includes columns such as:

```text
Lon, Lat, LOC, TYP, GCD, Y, EC1, EC2, ...
```

where:

* `LOC` is the location identifier;
* `TYP` is the context indicator (`EXP` or `CST`);
* `GCD` is the genotype or cultivar identifier;
* `Y` is the phenotypic response;
* `EC1`, `EC2`, etc. are environmental covariates.

## Model structure

The fitted model includes:

* fixed intercepts for EXP and CST;
* genotype-specific random intercepts;
* genotype-specific random slopes for environmental covariates;
* a 2 × 2 genetic variance--covariance matrix linking EXP and CST;
* context-specific residual variances;
* a genomic relationship matrix derived from SNP markers.

The genomic relationship matrix is used to connect related genotypes and to support prediction across partially observed contexts.

## Parameter estimation

Variance components are estimated by minimizing the negative REML log-likelihood using a derivative-free optimizer.

The implementation uses a bounded parameterization for:

* genetic variances for EXP and CST;
* the genetic correlation between EXP and CST;
* residual variances for EXP and CST.

The genetic covariance is reconstructed from the estimated genetic variances and genetic correlation.

## Outputs

The script returns or reports:

* estimated genetic variance--covariance matrix;
* estimated residual variances;
* genetic correlation between EXP and CST;
* BLUPs for genotype intercepts and environmental slopes;
* predicted values;
* fitted vs. observed summaries;
* computational diagnostics, including log-likelihood, number of iterations, and elapsed time.

## Notes

This folder is primarily focused on the DF-REML implementation. The empirical CassavaBase case study has its own folder and uses the article terminology BTY/CTY directly.

The EXP/CST notation is kept here to preserve consistency with the original simulation and implementation scripts.

```
```
