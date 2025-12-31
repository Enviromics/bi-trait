# ===============================
# 1) Simulation function
# ===============================
# This function simulates an EXP (experimental trials) + CST (commercial stands) dataset with:
# (i) A square raster-like TPE (n_grid x n_grid) on continuous coordinates in [1,10] x [1,10] with no jitter.
# (ii) A strict LOC ↔ pixel mapping: each LOC is assigned to exactly one raster pixel; no free coordinates.
# (iii) Environmental covariates (EC1..ECp) generated at pixel level and then copied to LOC level (EC constant within LOC).
# (iv) A latent-factor field (F1..Fq) built from separable 1D multiplicative random-walk surfaces along Lon/Lat,
#      then smoothed (loess), rescaled, optionally warped, and perturbed by Gaussian bumps; each latent field is scaled.
# (v) ECs derived from random subsets of latent fields with mild non-linear transforms and small additive noise, then scaled.
# (vi) Materials = GEN (breeding genotypes) + CUL (cultivars). SNPs are generated via AlphaSimR; M (0/1/2) and G (VanRaden) are returned.
# (vii) Bi-trait-like genetic structure for "EXP vs CST" is imposed by sampling correlated genetic effects using G (among materials) and G0 (2x2).
# (viii) Phenotypes are generated as: Y = mu_TYP + LOC_RE + EC_effect_TYP + genetic_effect_TYP + residual_TYP.
# (ix) EXP imbalance is "BLUE-like": for TYP == "EXP", there is at most one observation per (GCD, LOC).
#      - GEN presence is Bernoulli per (GEN, LOC_EXP) with pi_gen_exp, but forced constraints ensure:
#        each GEN appears at least once in EXP, and each LOC_EXP contains at least one GEN.
#      - CUL checks in EXP: core checks are deterministic (CULA, CULB, ...) up to n_cul_core and appear in all LOC_EXP;
#        remaining cultivars enter each LOC_EXP with pi_cul_var_exp.
# (x) CST is one cultivar per LOC (with replacement), but with a forced coverage constraint so that all CUL appear at least once overall.
# Notes on robustness:
# - Multiple numerical safeguards are included (symmetrization, diagonal jitter, safe chol) to avoid failures in G/G0 sampling.
# - Parameter validation enforces probability bounds and core-check bounds, and enforces n_loc_cst >= n_cul for CST full coverage.

rm(list = ls())
suppressPackageStartupMessages({
  library(AlphaSimR)
  library(AGHmatrix)
})

simulate_exp_cst <- function(par = list()) {

  .stopif <- function(cond, msg) if (isTRUE(cond)) stop(msg, call. = FALSE)
  .sym <- function(A) (A + t(A)) / 2

  .safe_chol <- function(A, eps0 = 1e-8, max_tries = 12) {
    A <- .sym(A)
    eps <- eps0 * max(1, mean(diag(A)))
    for (i in seq_len(max_tries)) {
      out <- try(chol(A + diag(eps, nrow(A))), silent = TRUE)
      if (!inherits(out, "try-error")) return(out)
      eps <- eps * 10
    }
    stop("Cholesky failed.", call. = FALSE)
  }

  .loess_smooth <- function(x, span = 0.60) {
    n <- length(x)
    if (n < 4) return(as.numeric(x))
    idx <- seq_len(n)
    fit <- try(stats::loess(x ~ idx, span = span), silent = TRUE)
    if (inherits(fit, "try-error")) return(as.numeric(x))
    as.numeric(stats::predict(fit, idx))
  }

  .rescale_range <- function(x, to = c(0.85, 1.15)) {
    rx <- range(x, finite = TRUE)
    if (!all(is.finite(rx)) || diff(rx) < 1e-12) return(rep(mean(to), length(x)))
    (x - rx[1]) / diff(rx) * (to[2] - to[1]) + to[1]
  }

  .rw_mult <- function(n, sd = 0.05) {
    x <- numeric(n); x[1] <- 1
    if (n > 1) for (i in 2:n) x[i] <- x[i - 1] * rnorm(1, mean = 1, sd = sd)
    x
  }

  .loc_ids <- function(n) sprintf("LOC%02d", seq_len(n))

  .gen_ids <- function(n, width = 0L) {
    n <- as.integer(n); width <- as.integer(width)
    if (width > 0L) sprintf(paste0("GEN%0", width, "d"), seq_len(n)) else paste0("GEN", seq_len(n))
  }

  .cul_ids <- function(n) {
    L <- toupper(letters)
    tag <- function(i) {
      s <- ""
      while (i > 0) {
        r <- (i - 1) %% 26
        s <- paste0(L[r + 1], s)
        i <- (i - 1) %/% 26
      }
      s
    }
    paste0("CUL", vapply(seq_len(as.integer(n)), tag, character(1)))
  }

  .scale_vec <- function(x) as.numeric(scale(x))

  .vanraden_manual <- function(M) {
    p <- colMeans(M) / 2
    W <- sweep(M, 2, 2 * p, "-")
    denom <- 2 * sum(p * (1 - p))
    .stopif(!is.finite(denom) || denom <= 0, "VanRaden denominator <= 0.")
    G <- (W %*% t(W)) / denom
    G <- .sym(G); diag(G) <- diag(G) + 1e-8
    G
  }

  .center01 <- function(v) (v - mean(v)) / max(1e-12, stats::sd(v))

  .gauss_bump <- function(n, cx, cy, sx, sy) {
    xs <- seq_len(n); ys <- seq_len(n)
    X <- outer(xs, rep(1, n)); Y <- outer(rep(1, n), ys)
    exp(-0.5 * ((X - cx)^2 / sx^2 + (Y - cy)^2 / sy^2))
  }

  par0 <- list(
    seed = 11L,
    n_grid = 20L,
    q_latent = 4L, p_ec = 5L, return_tpe = TRUE,
    n_loc_exp = 5L, n_loc_cst = 19L, overlap_locs = FALSE,
    n_gen = 5L, n_cul = 7L, gen_id_width = 0L,
    n_gen_exp_range = c(5L, 5L), n_chk_exp = 7L,
    n_cul_core = 3L, pi_gen_exp = 0.50, pi_cul_var_exp = 0.60,
    n_snp = 50L, n_chr = NULL, n_founders = 120L, seg_sites = 1400L,
    min_maf_try = c(0.05, 0.02, 0.01), n_progeny = 3L, n_cross_min = 60L,
    mu_exp = 250, mu_cst = 235, vg_exp = 80, vg_cst = 35, cov_g = 25, ve_exp = 80, ve_cst = 35,
    sd_loc = 6.0, sd_beta = 2.2,
    .loess_span = 0.95, .axis_range = c(0.85, 1.15), .sd_rw_mult = 0.05,
    .sd_latent_pix_mult = 0.005, .sd_latent_pix_add = 0.00, .sd_ec_add = 0.05,
    .warp_sd = 0.12, .n_bumps = 2L, .bump_amp = 0.25, .bump_scale = c(6, 18),
    y_digits = 1L
  )

  par <- utils::modifyList(par0, par)
  set.seed(as.integer(par$seed))

  # --- Parameter checks for EXP imbalance and CST coverage ---
  .stopif(length(par$n_cul_core) != 1L || !is.finite(par$n_cul_core), "n_cul_core must be finite (length 1).")
  .stopif(par$n_cul_core < 0 || par$n_cul_core > par$n_cul, "n_cul_core must be in [0, n_cul].")
  .stopif(length(par$pi_gen_exp) != 1L || !is.finite(par$pi_gen_exp) || par$pi_gen_exp < 0 || par$pi_gen_exp > 1,
          "pi_gen_exp must be in [0,1].")
  .stopif(length(par$pi_cul_var_exp) != 1L || !is.finite(par$pi_cul_var_exp) || par$pi_cul_var_exp < 0 || par$pi_cul_var_exp > 1,
          "pi_cul_var_exp must be in [0,1].")

  ng_in <- par$n_grid
  .stopif(length(ng_in) != 1 || !is.finite(ng_in) || ng_in < 2 || abs(ng_in - round(ng_in)) > 1e-12,
          "n_grid must be an integer >= 2.")
  ng <- as.integer(ng_in)

  lon_seq <- seq(1, 10, length.out = ng)
  lat_seq <- seq(1, 10, length.out = ng)
  tpe0 <- expand.grid(Lon = lon_seq, Lat = lat_seq)
  n_pix <- nrow(tpe0)

  q <- as.integer(par$q_latent)
  .stopif(q < 1, "q_latent must be >= 1.")
  F <- matrix(0, n_pix, q)

  axis_c <- .center01(seq_len(ng))

  for (j in seq_len(q)) {
    x_lon <- .rw_mult(ng, sd = par$.sd_rw_mult)
    x_lat <- .rw_mult(ng, sd = par$.sd_rw_mult)
    x_lon <- .rescale_range(.loess_smooth(x_lon, span = par$.loess_span), to = par$.axis_range)
    x_lat <- .rescale_range(.loess_smooth(x_lat, span = par$.loess_span), to = par$.axis_range)

    Z <- outer(x_lon, x_lat, FUN = "*")

    if (is.finite(par$.warp_sd) && par$.warp_sd > 0) {
      a <- rnorm(1, 0, par$.warp_sd); b <- rnorm(1, 0, par$.warp_sd)
      Z <- Z * exp(outer(axis_c, rep(1, ng)) * a + outer(rep(1, ng), axis_c) * b)
    }

    nb <- as.integer(par$.n_bumps)
    if (is.finite(nb) && nb > 0) {
      for (h in seq_len(nb)) {
        cx <- runif(1, 1, ng); cy <- runif(1, 1, ng)
        sx <- runif(1, par$.bump_scale[1], par$.bump_scale[2])
        sy <- runif(1, par$.bump_scale[1], par$.bump_scale[2])
        Z <- Z + rnorm(1, 0, par$.bump_amp) * .gauss_bump(ng, cx, cy, sx, sy)
      }
    }

    if (is.finite(par$.sd_latent_pix_mult) && par$.sd_latent_pix_mult > 0) {
      Z <- Z * matrix(rnorm(ng * ng, mean = 1, sd = par$.sd_latent_pix_mult), ng, ng)
    }
    if (is.finite(par$.sd_latent_pix_add) && par$.sd_latent_pix_add > 0) {
      Z <- Z + matrix(rnorm(ng * ng, mean = 0, sd = par$.sd_latent_pix_add), ng, ng)
    }

    F[, j] <- .scale_vec(as.vector(Z))
  }
  colnames(F) <- paste0("F", seq_len(q))

  p <- as.integer(par$p_ec)
  .stopif(p < 1, "p_ec must be >= 1.")
  EC <- matrix(0, n_pix, p)

  for (k in seq_len(p)) {
    m <- sample(seq_len(min(3L, q)), 1)
    idx <- sample(seq_len(q), m, replace = FALSE)
    base <- .scale_vec(F[, idx, drop = FALSE] %*% rnorm(m))

    tr <- sample(1:4, 1)
    z <- base
    if (tr == 2) z <- z + 0.20 * z^2
    if (tr == 3) z <- tanh(1.0 * z)
    if (tr == 4 && m >= 2) z <- z + 0.25 * (z * F[, idx[2]])

    if (is.finite(par$.sd_ec_add) && par$.sd_ec_add > 0) z <- z + rnorm(n_pix, 0, par$.sd_ec_add)
    EC[, k] <- .scale_vec(z)
  }
  colnames(EC) <- paste0("EC", seq_len(p))
  tpe <- if (isTRUE(par$return_tpe)) cbind(tpe0, F, EC) else NULL

  n_exp <- as.integer(par$n_loc_exp); n_cst <- as.integer(par$n_loc_cst)
  .stopif(n_exp < 1 || n_cst < 1, "n_loc_exp and n_loc_cst must be >= 1.")
  .stopif(!isTRUE(par$overlap_locs) && (n_exp + n_cst > n_pix), "Not enough pixels for non-overlapping LOCs.")

  pix_exp <- sample(seq_len(n_pix), n_exp, replace = FALSE)
  pix_cst <- if (isTRUE(par$overlap_locs)) sample(seq_len(n_pix), n_cst, replace = TRUE) else sample(setdiff(seq_len(n_pix), pix_exp), n_cst, replace = FALSE)

  loc_tbl <- data.frame(
    LOC = .loc_ids(n_exp + n_cst),
    TYP = c(rep("EXP", n_exp), rep("CST", n_cst)),
    pix = c(pix_exp, pix_cst),
    stringsAsFactors = FALSE
  )
  loc_tbl$Lon <- tpe0$Lon[loc_tbl$pix]
  loc_tbl$Lat <- tpe0$Lat[loc_tbl$pix]
  loc_tbl <- cbind(loc_tbl, EC[loc_tbl$pix, , drop = FALSE])

  gen_id <- .gen_ids(par$n_gen, width = par$gen_id_width)
  .stopif(length(gen_id) < 1L, "n_gen must be >= 1.")
  cul_id <- .cul_ids(par$n_cul)
  gcd_id <- c(gen_id, cul_id)

  .stopif(n_cst < length(cul_id), "n_loc_cst must be >= n_cul to force all CUL to appear in CST (1 CUL per LOC).")

  n_snp <- as.integer(par$n_snp)
  .stopif(n_snp < 5, "n_snp must be >= 5.")

  n_chr <- par$n_chr
  if (is.null(n_chr)) n_chr <- min(10L, n_snp)
  n_chr <- as.integer(n_chr)
  .stopif(n_chr < 1, "n_chr must be >= 1.")

  snp_per_chr <- rep(floor(n_snp / n_chr), n_chr)
  rem <- n_snp - sum(snp_per_chr)
  if (rem > 0) snp_per_chr[seq_len(rem)] <- snp_per_chr[seq_len(rem)] + 1L

  founders <- try(runMacs(
    nInd = as.integer(par$n_founders),
    nChr = n_chr,
    segSites = as.integer(par$seg_sites),
    inbred = FALSE,
    species = "GENERIC"
  ), silent = TRUE)
  if (inherits(founders, "try-error")) {
    founders <- quickHaplo(
      nInd = as.integer(par$n_founders),
      nChr = n_chr,
      segSites = as.integer(par$seg_sites),
      inbred = FALSE
    )
  }

  SP <- SimParam$new(founders)
  ok <- FALSE
  for (maf in par$min_maf_try) {
    tmp <- try(SP$addSnpChip(nSnpPerChr = snp_per_chr, minSnpFreq = maf), silent = TRUE)
    if (!inherits(tmp, "try-error")) { ok <- TRUE; break }
  }
  .stopif(!ok, "Failed to set SNP chip. Try increasing seg_sites or relaxing min_maf_try.")

  base_pop <- newPop(founders, simParam = SP)
  need1 <- length(cul_id); need2 <- length(gen_id)
  n_progeny <- as.integer(par$n_progeny)
  n_cross1 <- max(as.integer(par$n_cross_min), ceiling(need1 / n_progeny))
  n_cross2 <- max(as.integer(par$n_cross_min), ceiling(need2 / n_progeny))

  pop1 <- randCross(base_pop, nCrosses = n_cross1, nProgeny = n_progeny, simParam = SP)
  pop2 <- randCross(pop1,     nCrosses = n_cross2, nProgeny = n_progeny, simParam = SP)

  M1 <- pullSnpGeno(pop1, simParam = SP); storage.mode(M1) <- "integer"
  M2 <- pullSnpGeno(pop2, simParam = SP); storage.mode(M2) <- "integer"

  .stopif(nrow(M1) < need1, "Not enough phase-1 individuals for cultivars. Increase n_cross_min or n_progeny.")
  .stopif(nrow(M2) < need2, "Not enough phase-2 individuals for breeding genotypes. Increase n_cross_min or n_progeny.")

  idx_cul <- sample(seq_len(nrow(M1)), need1, replace = FALSE)
  idx_gen <- sample(seq_len(nrow(M2)), need2, replace = FALSE)

  M <- rbind(M2[idx_gen, , drop = FALSE], M1[idx_cul, , drop = FALSE])
  rownames(M) <- gcd_id
  colnames(M) <- sprintf("snp%02d", seq_len(ncol(M)))

  G <- try(AGHmatrix::Gmatrix(M, method = "VanRaden"), silent = TRUE)
  if (inherits(G, "try-error")) G <- .vanraden_manual(M)
  G <- .sym(G); diag(G) <- diag(G) + 1e-8
  rownames(G) <- colnames(G) <- gcd_id

  M_out <- data.frame(GCD = rownames(M), M, row.names = NULL, check.names = FALSE)

  # --- EXP imbalance (BLUE-like): <= 1 obs per (GCD, LOC) ---
  # Checks: core checks are deterministic (CULA, CULB, ...), appear in all EXP LOCs
  # Other cultivars enter irregularly by LOC with pi_cul_var_exp
  # Genotypes: presence/absence by (GEN, LOC_EXP); each GEN appears at least once; each LOC_EXP has >= 1 GEN

  n_cul_core <- as.integer(par$n_cul_core)
  pi_gen_exp <- as.numeric(par$pi_gen_exp)
  pi_cul_var_exp <- as.numeric(par$pi_cul_var_exp)

  cul_core <- if (n_cul_core > 0L) cul_id[seq_len(n_cul_core)] else character(0)
  cul_var <- setdiff(cul_id, cul_core)

  pres_gen <- matrix(stats::rbinom(n_exp * length(gen_id), 1L, pi_gen_exp) == 1L, nrow = n_exp)
  cg <- colSums(pres_gen)
  if (any(cg == 0L)) {
    z <- which(cg == 0L)
    for (k in z) pres_gen[sample.int(n_exp, 1), k] <- TRUE
  }
  rg <- rowSums(pres_gen)
  if (any(rg == 0L)) {
    z <- which(rg == 0L)
    for (j in z) pres_gen[j, sample.int(length(gen_id), 1)] <- TRUE
  }

  pres_cul <- NULL
  if (length(cul_var) > 0L && pi_cul_var_exp > 0) {
    pres_cul <- matrix(stats::rbinom(n_exp * length(cul_var), 1L, pi_cul_var_exp) == 1L, nrow = n_exp)
  }

  exp_rows <- vector("list", n_exp)
  for (j in seq_len(n_exp)) {
    g_sel <- gen_id[pres_gen[j, ]]
    c_sel <- cul_core
    if (!is.null(pres_cul)) c_sel <- c(c_sel, cul_var[pres_cul[j, ]])
    c_sel <- unique(c_sel)

    exp_rows[[j]] <- data.frame(
      Lon = loc_tbl$Lon[j], Lat = loc_tbl$Lat[j], LOC = loc_tbl$LOC[j],
      TYP = "EXP",
      CHK = c(rep("NO", length(g_sel)), rep("YES", length(c_sel))),
      GCD = c(g_sel, c_sel),
      stringsAsFactors = FALSE
    )
  }
  exp_df <- do.call(rbind, exp_rows)
  .stopif(any(duplicated(paste(exp_df$GCD, exp_df$LOC))), "EXP duplicated GCD x LOC.")

  # --- CST: 1 cultivar per LOC (with replacement), forcing all CUL to appear at least once overall ---
  cst_idx <- which(loc_tbl$TYP == "CST")
  n_cst_loc <- length(cst_idx)
  cst_draw <- sample(cul_id, n_cst_loc, replace = TRUE)
  missing_cul <- setdiff(cul_id, cst_draw)
  if (length(missing_cul) > 0L) {
    idx_replace <- sample.int(n_cst_loc, length(missing_cul), replace = FALSE)
    cst_draw[idx_replace] <- missing_cul
  }

  cst_df <- data.frame(
    Lon = loc_tbl$Lon[cst_idx], Lat = loc_tbl$Lat[cst_idx], LOC = loc_tbl$LOC[cst_idx],
    TYP = "CST", CHK = NA_character_,
    GCD = cst_draw,
    stringsAsFactors = FALSE
  )

  df <- rbind(exp_df, cst_df)
  loc_map <- match(df$LOC, loc_tbl$LOC)
  df <- cbind(df, loc_tbl[loc_map, colnames(EC), drop = FALSE])

  G0 <- matrix(c(par$vg_exp, par$cov_g, par$cov_g, par$vg_cst), 2, 2)
  U <- (t(.safe_chol(G)) %*% matrix(rnorm(length(gcd_id) * 2), nrow = length(gcd_id), ncol = 2) %*% .safe_chol(G0))
  rownames(U) <- gcd_id; colnames(U) <- c("EXP", "CST")

  loc_re <- rnorm(nrow(loc_tbl), 0, par$sd_loc)
  beta_exp <- rnorm(p, 0, par$sd_beta)
  beta_cst <- rnorm(p, 0, par$sd_beta)

  idx_exp_obs <- df$TYP == "EXP"
  EC_obs <- as.matrix(df[, colnames(EC), drop = FALSE])

  env <- numeric(nrow(df))
  env[idx_exp_obs]  <- as.vector(EC_obs[idx_exp_obs,  , drop = FALSE] %*% beta_exp)
  env[!idx_exp_obs] <- as.vector(EC_obs[!idx_exp_obs, , drop = FALSE] %*% beta_cst)

  typ_idx <- ifelse(idx_exp_obs, 1L, 2L)
  u_i <- U[df$GCD, , drop = FALSE][cbind(seq_len(nrow(df)), typ_idx)]

  mu <- ifelse(idx_exp_obs, par$mu_exp, par$mu_cst)
  e_sd <- ifelse(idx_exp_obs, sqrt(par$ve_exp), sqrt(par$ve_cst))

  df$Y <- round(mu + loc_re[loc_map] + env + u_i + rnorm(nrow(df), 0, e_sd), as.integer(par$y_digits))

  pheno <- df[, c("Lon", "Lat", "LOC", "TYP", "CHK", "GCD", "Y", colnames(EC)), drop = FALSE]
  out <- list(pheno = pheno, M = M_out, G = G)
  if (isTRUE(par$return_tpe)) out$tpe <- tpe
  out
}

# ===============================
# 2) Apply simulation (numeric)
# ===============================
# This block defines a 'par' list and runs the simulation, then prints minimal diagnostics to verify core constraints.
# Key parameters:
# seed: RNG seed for full reproducibility of raster, SNPs, and phenotypes.
# n_grid: raster resolution (n_grid x n_grid) used to define the TPE and the pool of available pixels for LOC sampling.
# q_latent: number of latent spatial factors (F1..Fq) used to construct ECs.
# p_ec: number of environmental covariates EC1..ECp returned in pheno and (optionally) in tpe.
# return_tpe: if TRUE, returns the full raster table (Lon, Lat, F's, EC's) as sim$tpe.
# n_loc_exp: number of EXP environments (LOC01..LOCn), sampled from distinct pixels when overlap_locs = FALSE.
# n_loc_cst: number of CST environments; CST pixels are sampled from remaining pixels (no overlap) or with replacement if overlap_locs = TRUE.
# overlap_locs: if TRUE, EXP and CST can share pixels; if FALSE, EXP and CST pixels are disjoint.
# n_gen: number of breeding genotypes (GEN1..GENn).
# n_cul: number of cultivars (CULA..); these are checks in EXP and cultivars in CST.
# n_cul_core: number of core checks in EXP; deterministic selection as the first cultivars (CULA..).
# pi_gen_exp: Bernoulli probability that a given GEN is observed in a given LOC_EXP (presence/absence imbalance driver).
# pi_cul_var_exp: Bernoulli probability that a non-core cultivar is included as a check in a given LOC_EXP.
# n_snp: number of SNP markers; controls M dimension and the resolution of G.
# Post-run checks printed:
# - Basic counts and inclusion check: all observed GCD are present in M.
# - TYP x CHK contingency to confirm encoding (EXP has YES/NO; CST has NA).
# - No duplicated (GCD, LOC) pairs in EXP.
# - Core checks present in all EXP LOCs.
# - All cultivars appear at least once in CST, and CST contains only CUL.
# - ECs constant within LOC (by checking each EC has one unique value per LOC).

par <- list(
  seed = 2026,
  n_grid = 20,
  q_latent = 4,
  p_ec = 5,
  return_tpe = TRUE,
  n_loc_exp = 4,
  n_loc_cst = 20,
  overlap_locs = FALSE,
  n_gen = 5,
  n_cul = 7,
  n_cul_core = 2,
  pi_gen_exp = 0.30,
  pi_cul_var_exp = 0.30,
  n_snp = 50
)

sim <- simulate_exp_cst(par)
pheno <- sim$pheno

cat("\n--- pheno (head) ---\n")
print(head(pheno, 10))

cat("\n--- Counts ---\n")
cat("Rows:", nrow(pheno), "\n")
cat("LOCs:", length(unique(pheno$LOC)), "\n")
cat("GCDs observed:", length(unique(pheno$GCD)), "\n")
cat("All GCD in M:", all(unique(pheno$GCD) %in% sim$M$GCD), "\n")
cat("M rows (GEN+CUL):", nrow(sim$M), "\n")
cat("G dim:", paste(dim(sim$G), collapse = " x "), "\n")

cat("\n--- TYP x CHK ---\n")
print(with(pheno, table(TYP, CHK, useNA = "ifany")))

cat("\n--- EXP duplicated GCD x LOC ---\n")
exp_df <- pheno[pheno$TYP == "EXP", , drop = FALSE]
dup_exp <- any(duplicated(paste(exp_df$GCD, exp_df$LOC)))
print(!dup_exp)

cul_all <- sim$M$GCD[grepl("^CUL", sim$M$GCD)]
core_exp <- if (par$n_cul_core > 0L) cul_all[seq_len(par$n_cul_core)] else character(0)

cat("\n--- Core checks present in all EXP LOCs ---\n")
ok_core <- TRUE
if (length(core_exp) > 0L) {
  locs_exp <- unique(exp_df$LOC)
  ok_core <- all(vapply(locs_exp, function(L) all(core_exp %in% exp_df$GCD[exp_df$LOC == L]), logical(1)))
}
print(ok_core)

cat("\n--- All cultivars appear at least once in CST ---\n")
cst_df <- pheno[pheno$TYP == "CST", , drop = FALSE]
ok_cst_cov <- all(cul_all %in% cst_df$GCD)
print(ok_cst_cov)

cat("\n--- CST contains only CUL ---\n")
print(all(grepl("^CUL", cst_df$GCD)))

cat("\n--- Y summary ---\n")
print(summary(pheno$Y))

ec_cols <- grep("^EC\\d+$", names(pheno), value = TRUE)
ok_ec <- TRUE
for (k in ec_cols) {
  v <- tapply(pheno[[k]], pheno$LOC, function(x) length(unique(round(x, 10))) == 1L)
  ok_ec <- ok_ec && all(unlist(v))
}
cat("\n--- ECs constant within LOC ---\n")
print(ok_ec)

# ===============================
# 3) Visualization
# ===============================
# This block provides visual inspection of: (i) the MET-like response table, (ii) spatial layout, (iii) raster EC fields,
# (iv) an illustrative association between Y and EC1, and (v) a publication-quality genomic relationship heatmap.
# Plots produced:
# (A) Yield (Y) across environments (LOC) for each material:
#     - Lines are grouped by material x type (interaction(GEN, TYP)) and linetype encodes TYP (EXP solid, CST dashed).
#     - Shapes encode material type globally: GEN uses circles (16), CUL uses triangles (17), for all LOC (EXP and CST).
#     - Text labels are placed once per (material, TYP) to keep the plot readable while preserving identity.
# (B) LOC map (Lon vs Lat):
#     - Shows where each LOC is located on the [1,10]x[1,10] coordinate system to confirm the raster sampling logic.
# (C) Raster EC maps (facet over EC1..ECp):
#     - Visualizes the generated EC surfaces over the full TPE grid; useful to verify spatial patterns and scaling.
# (D) Raster EC1 + sampled LOC points:
#     - Overlays EXP/CST LOC sampling locations on top of the EC1 raster to confirm strict pixel mapping.
# (E) Quick EC1 association diagnostic:
#     - Plots a TYP-adjusted Y vs EC1 scatter and overlays a simple fitted line from Y ~ EC1 + TYP (visual sanity check).
# (F) G matrix heatmap with dendrogram (publication-oriented):
#     - Computes a clustering order from a distance derived from G, then draws:
#       * a left dendrogram (ggdendro) aligned to the heatmap rows,
#       * a heatmap of G in the clustered order,
#       * diagonal rectangles highlighting k clusters (k_clusters),
#       * an optional vertical dashed cut-height indicator in the dendrogram panel.
#     - Intended to inspect relatedness structure among GEN and CUL and provide a figure-ready summary of G.


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  library(ggdendro)
  library(patchwork)
})

dados <- pheno |>
  transmute(
    TYP = as.character(TYP),
    LOC = as.character(LOC),
    CHK = as.character(CHK),
    CHK_plot = ifelse(is.na(CHK), "CST", CHK),
    GEN = as.character(GCD),
    Y = as.numeric(Y),
    Lon = as.numeric(Lon),
    Lat = as.numeric(Lat)
  )

coords <- dados |>
  distinct(LOC, Lon, Lat) |>
  arrange(LOC)

dados$LOC <- factor(dados$LOC, levels = sort(unique(dados$LOC)))
labels_duplos <- dados |>
  arrange(GEN, TYP, LOC) |>
  group_by(GEN, TYP) |>
  slice(1) |>
  ungroup()

ggplot(dados, aes(LOC, Y, group = interaction(GEN, TYP), color = GEN)) +
  geom_line(aes(linetype = TYP), linewidth = 0.9) +
  geom_point(aes(shape = ifelse(grepl("^GEN", GEN), "GEN", "CUL")), size = 2.0) +
  geom_text(
    data = labels_duplos, aes(label = GEN),
    hjust = 1.15, vjust = 0.5, size = 2.8
  ) +
  scale_linetype_manual(values = c(EXP = "solid", CST = "dashed")) +
  scale_shape_manual(values = c(GEN = 16, CUL = 17), breaks = NULL) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.01))) +
  labs(x = "Environments (LOC)", y = "Yield (Y)") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggplot(coords, aes(Lon, Lat, label = LOC)) +
  geom_point(size = 2.6) +
  geom_text(vjust = -0.6, size = 3) +
  coord_fixed() +
  scale_x_continuous(limits = c(1, 10), breaks = 1:10, expand = c(0, 0)) +
  scale_y_reverse(limits = c(10, 1), breaks = 1:10, expand = c(0, 0)) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude")

if (!is.null(sim$tpe)) {

  EC <- sim$tpe |>
    select(Lon, Lat, starts_with("EC"))

  # --- Pixel size in plotting coordinates (works for any n_grid, avoids needing n_grid) ---
  ux <- sort(unique(EC$Lon))
  uy <- sort(unique(EC$Lat))
  dx <- if (length(ux) > 1) min(diff(ux)) else 0
  dy <- if (length(uy) > 1) min(diff(uy)) else 0

  # Expand limits by half a pixel so tiles at 1 and 10 are not clipped
  lims_x <- c(1 - dx/2, 10 + dx/2)
  lims_y <- c(10 + dy/2, 1 - dy/2)

  EC_long <- EC |>
    pivot_longer(starts_with("EC"), names_to = "EC", values_to = "Value")

  print(
    ggplot(EC_long, aes(Lon, Lat, fill = Value)) +
      geom_tile(width = dx, height = dy) +
      facet_wrap(~ EC, ncol = min(5, length(unique(EC_long$EC)))) +
      scale_fill_gradientn(colors = viridis(100)) +
      coord_fixed() +
      scale_x_continuous(limits = lims_x, breaks = 1:10, expand = c(0, 0)) +
      scale_y_reverse(limits = lims_y, breaks = 1:10, expand = c(0, 0)) +
      labs(x = "Longitude", y = "Latitude", fill = "Value") +
      theme_bw() +
      theme(
        strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 8)
      )
  )

  locs <- dados |>
    distinct(Lon, Lat, LOC, TYP) |>
    arrange(LOC)

  print(
    ggplot() +
      geom_tile(data = EC, aes(x = Lon, y = Lat, fill = EC1), width = dx, height = dy) +
      geom_point(data = locs, aes(x = Lon, y = Lat, shape = TYP),
                 size = 2.2, color = "black") +
      coord_fixed() +
      scale_fill_gradientn(colors = viridis(100)) +
      scale_x_continuous(limits = lims_x, breaks = 1:10, expand = c(0, 0)) +
      scale_y_reverse(limits = lims_y, breaks = 1:10, expand = c(0, 0)) +
      theme_bw() +
      labs(x = "Longitude", y = "Latitude", fill = "EC1")
  )

  dat <- left_join(dados, EC, by = c("Lon", "Lat"))
  fit0 <- lm(Y ~ TYP, data = dat)
  y_adj <- residuals(fit0) + coef(fit0)[1]
  fit1 <- lm(Y ~ EC1 + TYP, data = dat)

  plot(dat$EC1, y_adj, xlab = "EC1", ylab = "Y (TYP-adjusted + intercept)")
  abline(coef(fit1)[c("(Intercept)", "EC1")], lty = 2, col = "blue")
}


# --- G heatmap + left dendrogram + cluster boxes (ggplot2) ---
# Assumes you already have: G <- sim$G  (numeric matrix with rownames/colnames = GCD)

G <- sim$G
stopifnot(is.matrix(G), nrow(G) == ncol(G), !is.null(rownames(G)))

# --- Hierarchical clustering order (distance from G; always non-negative) ---
D  <- as.dist(max(G, na.rm = TRUE) - G)
hc <- hclust(D, method = "average")
dd <- ggdendro::dendro_data(hc, type = "rectangle")

ord   <- dd$labels$label[order(dd$labels$x)]
G_ord <- G[ord, ord]
n     <- length(ord)

# --- Heatmap data (numeric axes, y reversed so ord[1] appears at the top) ---
df <- expand.grid(Row = seq_len(n), Col = seq_len(n))
df$Value <- G_ord[cbind(df$Row, df$Col)]
df$x <- df$Col
df$y <- n - df$Row + 1

# --- Cluster boxes along the diagonal (choose k as you like) ---
k_clusters <- 4L
cl <- cutree(hc, k = k_clusters)[ord]
r  <- rle(cl)
ends   <- cumsum(r$lengths)
starts <- c(1L, head(ends, -1L) + 1L)

rects <- data.frame(
  xmin = starts - 0.5,
  xmax = ends + 0.5,
  ymin = (n - ends + 1) - 0.5,
  ymax = (n - starts + 1) + 0.5
)

# --- Dendrogram segments mapped to the same y scale as the heatmap ---
seg <- dd$segments
seg$pos    <- n - seg$x + 1
seg$posend <- n - seg$xend + 1

# --- Optional cut height indicator (vertical dashed red line) ---
h_cut <- 0.60 * max(hc$height)

p_dend <- ggplot(seg) +
  geom_segment(aes(x = y, y = pos, xend = yend, yend = posend),
               linewidth = 0.45, lineend = "round") +
  geom_vline(xintercept = h_cut, linetype = "dashed", linewidth = 0.7, color = "red") +
  scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
  scale_x_reverse(expand = expansion(mult = c(0.02, 0.02))) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 5))

lim <- range(df$Value, finite = TRUE)

p_heat <- ggplot(df, aes(x = x, y = y, fill = Value)) +
  geom_tile() +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA, color = "black", linewidth = 0.5
  ) +
  coord_fixed(expand = FALSE) +
  scale_x_continuous(
    breaks = seq_len(n),
    labels = ord,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = seq_len(n),
    labels = rev(ord),
    position = "right",
    expand = c(0, 0)
  ) +
  scale_fill_gradient2(
    low = "#b3a2ff", mid = "white", high = "#e41a1c",
    midpoint = 0, limits = lim,
    name = "Genomic\nRelationship"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 5, 5, 0),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

p <- p_dend + p_heat + plot_layout(widths = c(1.2, 4.6), guides = "collect") &
  theme(legend.position = "right")

print(p)

# Optional save (high-res)
# ggsave("G_heatmap_dendrogram.png", p, width = 12, height = 8, dpi = 300)
