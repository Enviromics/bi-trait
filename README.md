# bi-trait
R code for bi-trait enviromics data simulation and mixed-model inference, contrasting experimental trials (EXP) and commercial stands (CST) in reproducible methodological studies.

## 1. Overview

This repository provides R code supporting a methodological framework for bi-trait enviromics analyses, where the two traits represent distinct data-generation contexts: experimental trials (EXP) and commercial stands (CST). The objective is to study how information generated under controlled experimental conditions can be statistically integrated with large-scale, heterogeneous commercial data using mixed models with environmental covariates.

The repository includes (i) a spatially explicit simulator to generate paired EXP and CST datasets under controlled imbalance, shared genetic material, and structured enviromic gradients; and (ii) a derivative-free REML (DF-REML) implementation for fitting bi-context mixed models, with emphasis on estimating genetic and residual variance–covariance components between EXP and CST.

The code is intended for reproducible methodological investigations in plant breeding and quantitative genetics.

## 2. Simulation of EXP/CST enviromic data

Simulation is implemented in the function `simulate_exp_cst()` and targets paired datasets for EXP (experimental trials) and CST (commercial stands). Design choices are intentionally strict: a raster-like target population of environments (TPE), a fixed pixel grid, a one-to-one LOC–pixel mapping, and environmental covariates (ECs) that are constant within each location.

### Core design principles

- **Spatial structure**  
  The TPE is defined as a square grid (`n_grid × n_grid`) on continuous coordinates within `[1,10] × [1,10]`. No jittering or free coordinates are allowed; each location is assigned to exactly one raster pixel.

- **Enviromics**  
  Latent spatial factors are constructed from separable one-dimensional random-walk surfaces along longitude and latitude. Environmental covariates (EC1…ECp) are derived from these latent fields through simple non-linear transformations and scaling, remaining constant within locations.

- **Genetic material**  
  The simulated population includes breeding genotypes (GEN) and cultivars (CUL, used as checks). SNP markers are generated using AlphaSimR, yielding a marker matrix (M) and a genomic relationship matrix (G) computed via the VanRaden method.

- **Phenotype generation**  
  Phenotypes combine context-specific means (EXP, CST), location effects, enviromic contributions, genetic effects with cross-context covariance, and residual noise.

- **Sampling and imbalance**  
  EXP data follow a BLUE-like structure with at most one observation per genotype–location combination. Genotype presence per location is probabilistic, subject to connectivity constraints. A subset of cultivars acts as core checks present in all EXP locations, with additional cultivars entering variably. CST data represent commercial stands, with one cultivar per location and enforced coverage across cultivars.

### Key simulation controls

Frequently modified parameters include:
- Spatial resolution and environments: `n_grid`, `n_loc_exp`, `n_loc_cst`
- Enviromic structure: `q_latent`, `p_ec`
- Genetic design: `n_gen`, `n_cul`, `n_snp`
- EXP imbalance and checks: `pi_gen_exp`, `n_cul_core`, `pi_cul_var_exp`
- Signal levels: `mu_exp`, `mu_cst`, `vg_exp`, `vg_cst`, `cov_g`, `ve_exp`, `ve_cst`

### Outputs

The function returns `list(pheno, M, G, tpe)`, where:

**pheno** is the observation-level table used for downstream inference, with one row per evaluated material × environment. Columns are:
- `Lon`, `Lat`: continuous spatial coordinates of the location, inherited from the raster pixel.
- `LOC`: environment (location) identifier.
- `TYP`: data-generation context (`EXP` = experimental trial; `CST` = commercial stand).
- `CHK`: check indicator (`YES`/`NO` for EXP; `NA` for CST).
- `GCD`: genetic material identifier (breeding genotype `GEN*` or cultivar `CUL*`).
- `Y`: simulated response variable.
- `EC1` … `ECp`: environmental covariates, constant within each `LOC`.

**M** is the marker matrix (one row per genetic material, one column per SNP), coded as 0/1/2 allele counts and used to represent the genomic information underlying genetic effects.

**G** is the genomic relationship matrix computed from `M` (VanRaden formulation), capturing realized genetic relatedness among all GEN and CUL entries and used directly in the bi-trait mixed model.

**tpe** is an optional square raster-level table describing the target population of environments, containing pixel coordinates and the latent spatial factors and environmental covariates used to generate the location-level enviromic structure.

## 3. Bi-trait enviromics model (EXP vs CST)

This repository implements a bi-trait enviromics mixed model contrasting EXP and CST within a unified REML framework. The two contexts are treated as correlated observation systems of the same biological entities, enabling joint inference while preserving their distinct data-generation processes.

Model fitting is performed via a **Derivative-Free REML (DF-REML)** algorithm and follows a random regression formulation on environmental covariates, allowing genotype-specific environmental responses to differ between EXP and CST.

### Model structure

The model includes:
- Context-specific fixed intercepts (EXP, CST).
- Genotype-level random regression coefficients on ECs, nested within context.
- A 2×2 genetic variance–covariance matrix (G₀) linking EXP and CST.
- Context-specific residual variances (R₀).

Genetic relationships are modeled using a genomic relationship matrix derived from SNP markers, ensuring consistency between simulation and inference.

### Input data files

The DF-REML analysis operates on:
- `phenos_and_ECs.txt`: phenotypic data with genotype identifiers, context (EXP/CST), locations, coordinates, responses, and environmental covariates.
- `M.txt`: SNP marker matrix used to construct the genomic relationship matrix.

These files are directly compatible with the simulation output and allow full reproduction of the analyses.

### Outputs

The DF-REML script returns variance–covariance estimates (G₀, R₀), genetic correlations between EXP and CST, BLUPs for genotype intercepts and EC slopes, predicted values, and diagnostic summaries. The framework supports direct evaluation of how genotypes assessed in EXP are expected to perform when deployed as CST across heterogeneous environments.
