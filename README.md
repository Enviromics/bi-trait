# bi-trait
R code for bi-trait enviromics data simulation and mixed-model inference, contrasting experimental trials (EXP) and commercial stands (CST) in reproducible methodological studies.

## 1. Overview

This repository provides R code supporting a methodological framework for bi-trait enviromics analyses, where the two traits represent distinct data-generation contexts: experimental trials (EXP) and commercial stands (CST). The objective is to study how information generated under controlled experimental conditions can be statistically integrated with large-scale, heterogeneous commercial data using mixed models enriched with environmental covariates.

The repository includes (i) a spatially explicit simulator designed to generate EXP and CST datasets under controlled imbalance, shared genetic material, and structured enviromic gradients; and (ii) a derivative-free REML (DF-REML) implementation for fitting bi-context mixed models, focusing on the estimation of genetic and residual variance–covariance components between EXP and CST.

The code is intended for reproducible methodological investigations in plant breeding and quantitative genetics, particularly for assessing the feasibility and limits of linking experimental and commercial performance through enviromics. The repository accompanies a manuscript under review and is not designed as a production-ready analysis pipeline.

## 2. Simulation of EXP/CST enviromic data

Simulation is implemented in `simulate_exp_cst()` and targets a paired dataset with two contexts: EXP (experimental trials) and CST (commercial stands). Design choices are strict on purpose: a raster-like TPE with a fixed pixel grid, a strict LOC↔pixel mapping, and ECs copied from pixel to LOC (constant within LOC). :contentReference[oaicite:0]{index=0}

### Core design (what is enforced)
2.1. **TPE and coordinates**
   - Square grid (`n_grid × n_grid`) on continuous coordinates in `[1,10] × [1,10]`.
   - No jitter and no “free” coordinates; each LOC is assigned to exactly one pixel.

2.2. **Enviromics (ECs)**
   - Latent spatial factors (`q_latent`: F1..Fq) built from separable 1D multiplicative random-walk surfaces along Lon/Lat, then smoothed and scaled.
   - ECs (`p_ec`: EC1..ECp) derived from subsets of latent factors using mild non-linear transforms and small noise; EC values are constant within LOC.

2.3. **Genetics**
   - Materials = GEN (breeding genotypes) + CUL (cultivars/checks).
   - SNPs simulated with AlphaSimR (`n_snp`), returning `M` (0/1/2) and `G` (VanRaden).

2.4. **Phenotypes**
   - Context-specific means (`mu_exp`, `mu_cst`), LOC random effect (`sd_loc`), EC contribution (context-specific regression; `sd_beta`), genetic effects with cross-context covariance (`vg_exp`, `vg_cst`, `cov_g`), and residuals (`ve_exp`, `ve_cst`).

2.5. **Sampling rules and imbalance**
   - **EXP:** BLUE-like (≤ 1 obs per (GCD, LOC)).
     - GEN presence per LOC_EXP is Bernoulli (`pi_gen_exp`), with constraints: each GEN appears ≥ 1 time and each LOC_EXP contains ≥ 1 GEN.
     - CUL checks: core checks (`n_cul_core`) appear in all EXP LOCs; remaining cultivars enter by LOC with `pi_cul_var_exp`.
   - **CST:** one cultivar per LOC; enforced coverage so every CUL appears at least once overall (`n_loc_cst ≥ n_cul`).

### Minimal parameter map (what you actually change most often)
- **Reproducibility:** `seed`
- **Spatial resolution / number of environments:** `n_grid`, `n_loc_exp`, `n_loc_cst`, `overlap_locs`
- **Enviromics:** `q_latent`, `p_ec`
- **Genetics:** `n_gen`, `n_cul`, `n_snp`
- **EXP imbalance / checks:** `pi_gen_exp`, `n_cul_core`, `pi_cul_var_exp`
- **Signal levels:** `mu_exp`, `mu_cst`, `vg_exp`, `vg_cst`, `cov_g`, `ve_exp`, `ve_cst`, `sd_loc`, `sd_beta`

### Outputs
`simulate_exp_cst()` returns `list(pheno, M, G, tpe)` where `tpe` is optional (pixel-level table with latent factors and ECs).

## 3. Bi-trait enviromics model (EXP vs CST)

This repository implements a bi-trait enviromics mixed model contrasting **experimental trials (EXP)** and **commercial stands (CST)** within a unified REML framework. EXP and CST are treated as two correlated observation contexts of the same biological system, enabling joint inference while preserving their distinct data-generating mechanisms.

Model fitting is performed via a **Derivative-Free REML (DF-REML)** algorithm. The formulation follows a random regression structure on environmental covariates (ECs), allowing genotype-specific environmental responses to be estimated separately for EXP and CST.

### Model specification

For each observation, the response is modeled using:
- Context-specific fixed intercepts (EXP, CST).
- Genotype-level random regression coefficients on ECs, nested within context.
- A 2×2 genetic variance–covariance matrix (G₀) linking EXP and CST.
- Context-specific residual variances (R₀).

Genetic relationships among genotypes are modeled using a genomic relationship matrix derived from marker data (VanRaden formulation), ensuring consistency between simulation and inference.

### Input data files

The DF-REML analysis operates on two input files:

- `phenos_and_ECs.txt`  
  Phenotypic dataset containing genotype identifiers, context labels (EXP/CST), locations, spatial coordinates, response values, and associated environmental covariates.

- `M.txt`  
  SNP marker matrix used to construct the genomic relationship matrix.

These files are either produced by the simulation module or provided directly for reproducing the analyses reported in the associated manuscript.

### Outputs

The DF-REML script returns:
- Estimates of genetic and residual variance–covariance components (G₀ and R₀).
- Genetic correlation between EXP and CST.
- BLUPs for genotype intercepts and EC slopes by context.
- Predicted values for observed data and for genotypes projected as CST.
- Diagnostic summaries and visualization-ready objects.

This framework enables direct evaluation of how genotypes evaluated under EXP conditions are expected to perform when deployed as CST across heterogeneous environments.

