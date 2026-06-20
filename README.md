````markdown
# bi-trait

R code and processed input files supporting the manuscript:

**How to enviromically predict breeding genotypes as if they were commercial cultivars?**

This repository supports a bi-trait genomic--enviromic framework for predicting breeding-stage genotypes on a cultivar-target scale. The central idea is to treat Breeding Trial Yield (BTY) and Cultivar Target Yield (CTY) as distinct but related trait-stage expressions of yield, and to model them jointly using genomic relationships, BTY--CTY covariance, and environmental information.

## Overview

The repository contains two complementary components:

1. **Illustrative simulation**  
   A controlled simulated example used to demonstrate the bi-trait genomic--enviromic DF-REML implementation, inspect variance--covariance estimates, compare identity and SNP-based relationship structures, and generate spatial recommendation maps.

2. **Empirical CassavaBase case study**  
   A real-data case study using processed CassavaBase records from the IITA and NaCRRI breeding programs. This analysis evaluates BTY-to-CTY prediction under an operational, unbalanced breeding-data setting.

The empirical analysis compares three candidate models:

- `H0`: genomic BTY-only baseline;
- `H1`: bi-trait genomic model using BTY and CTY jointly;
- `H2`: bi-trait genomic--enviromic model using BTY, CTY, and environmental principal components.

Prediction performance is evaluated using bridge-genotype validation, location-centered evaluation, Leave-One-CTY-Location-Out (LOLO) validation, and LOLO-strict validation.

## Repository structure

```text
bi-trait/
│
├── CassavaBase/
│   ├── processed empirical input files
│   ├── phenotypic responses
│   ├── SNP marker matrix
│   ├── environmental covariates
│   └── supporting documentation for the empirical case study
│
├── df-reml/
│   └── derivative-free REML implementation scripts
│
├── simulation scripts and example inputs
│
└── README.md
````

## CassavaBase folder

The `CassavaBase` folder contains the processed input files used in the empirical case study. These files correspond only to the cleaned subset used in the manuscript, not to a redistribution of the full CassavaBase database.

The folder includes analysis-ready files for:

* processed phenotypic responses;
* genotype and location identifiers;
* SNP marker data;
* environmental covariates;
* basic file descriptions and consistency checks.

The empirical cassava data originate from [CassavaBase](https://cassavabase.org/). Users of these processed files should cite CassavaBase as the original data source.

## Simulation component

The simulation component provides a controlled example for testing the bi-trait genomic--enviromic model. It generates paired BTY and CTY data with shared genetic material, SNP markers, environmental covariates, and spatially structured target environments.

The simulation is not intended to replace the empirical validation. In the manuscript, it is used as an implementation example, while the CassavaBase case study provides the main empirical evaluation.

## Model implementation

The model is fitted using a derivative-free restricted maximum likelihood (DF-REML) approach. The implementation estimates genetic and residual variance--covariance components between BTY and CTY and allows genotype-specific responses to environmental covariates or environmental principal components.

The genomic relationship matrix is computed from SNP markers and used to connect related genotypes across trait stages.

## Data availability

The simulated data and the processed input files used in the CassavaBase case study are provided in this repository.

The empirical cassava breeding data originate from CassavaBase. The files made available here correspond only to the processed subset used in the manuscript, including phenotypic summaries, genotype and location identifiers, SNP markers, and extracted environmental covariates.

## Citation

If you use this repository, please cite the manuscript associated with this code and data:

Associated preprint: Resende, R. T. (2025). *How to enviromically predict breeding genotypes as if they were commercial cultivars?* bioRxiv. https://doi.org/10.1101/2025.05.28.656616

Please also cite CassavaBase when using the empirical cassava data.
```
```
