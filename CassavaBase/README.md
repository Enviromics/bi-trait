# Cassava BTY--CTY Data

This folder contains the cleaned input data used in the CassavaBase BTY--CTY case study.

## Files

| File       | Description                                                           |
| ---------- | --------------------------------------------------------------------- |
| `phe.csv`  | Phenotypic pseudo-responses at the `Loc × TS × Gen` level.            |
| `ec.csv`   | Environmental covariates for the 26 locations used in the case study. |
| `M10K.csv` | SNP marker matrix with 10,000 markers for the 378 genotypes.          |

## Phenotypic data

`phe.csv` contains the response variable used in the case study.

Main columns:

* `Program`: breeding program source;
* `Loc`: location code;
* `TS`: trait-stage code;
* `Gen`: genotype identifier;
* `n_obs`: number of original records used;
* `y_mean`: mean of original observations;
* `y_asterisk`: final pseudo-response used in the analyses.

Trait-stage codes:

| Code  | Definition                                                                         |
| ----- | ---------------------------------------------------------------------------------- |
| `BTY` | Breeding Trial Yield, derived from `FYLD` in Preliminary Yield Trials.             |
| `CTY` | Cultivar Target Yield, derived from `DYLD` in Variety Release and Regional Trials. |

## Environmental covariates

`ec.csv` contains location identifiers, coordinates, and environmental covariates.

Main columns:

* `Loc`
* `Loc_name`
* `Lon`
* `Lat`
* environmental covariates used as ECs.

## SNP markers

`M10K.csv` contains one row per genotype. The first column is `Gen`; the remaining columns are SNP markers coded numerically.

The 10,000 SNPs were selected without phenotype information, after basic marker filtering and even spacing among eligible markers.

## Data dimensions

| File       | Rows | Columns |
| ---------- | ---: | ------: |
| `phe.csv`  |  867 |       7 |
| `ec.csv`   |   26 |     387 |
| `M10K.csv` |  378 |   10001 |

## Basic checks

```r
phe <- data.table::fread("phe.csv")
ec  <- data.table::fread("ec.csv")
M   <- data.table::fread("M10K.csv")

setequal(unique(phe$Gen), M$Gen)
setequal(unique(phe$Loc), ec$Loc)
```

Expected counts:

```r
data.table::uniqueN(phe$Gen) # 378
data.table::uniqueN(phe$Loc) # 26
```
