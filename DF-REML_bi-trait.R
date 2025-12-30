# ============================================================
# 1) DF-REML bi-trait (EXP vs CST) with random regression on ECs
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(nloptr)
  library(AGHmatrix)
})

.sym <- function(A) (A + t(A)) / 2

.read_M <- function(file) {
  m1 <- try(read.table(file, header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE), silent = TRUE)
  if (!inherits(m1, "try-error")) {
    M <- as.matrix(m1)
    storage.mode(M) <- "integer"
    return(M)
  }
  m2 <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
  M <- as.matrix(m2)
  storage.mode(M) <- "integer"
  M
}

.vanraden_manual <- function(M) {
  p <- colMeans(M) / 2
  W <- sweep(M, 2, 2 * p, "-")
  denom <- 2 * sum(p * (1 - p))
  if (!is.finite(denom) || denom <= 0) stop("VanRaden denominator <= 0.", call. = FALSE)
  G <- (W %*% t(W)) / denom
  G <- .sym(G); diag(G) <- diag(G) + 1e-8
  G
}

.build_Z <- function(dat, ec_cols) {
  N <- nrow(dat)
  G <- nlevels(dat$GEN)
  Tn <- nlevels(dat$TYP)
  if (Tn != 2L) stop("TYP must have exactly 2 levels.", call. = FALSE)

  W <- cbind(1, as.matrix(dat[, ec_cols, drop = FALSE]))
  n_re <- ncol(W)

  gen_i <- as.integer(dat$GEN)
  typ_i <- as.integer(dat$TYP)

  col0 <- (gen_i - 1L) * (Tn * n_re) + (typ_i - 1L) * n_re
  Z <- matrix(0, nrow = N, ncol = Tn * n_re * G)

  rr <- seq_len(N)
  for (k in seq_len(n_re)) Z[cbind(rr, col0 + k)] <- W[, k]

  list(Z = Z, typ_i = typ_i, n_re = n_re)
}

.fit_dfreml_bitra <- function(dat, M_mat) {

  need <- c("Y", "TYP", "LOC", "Lon", "Lat")
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop(paste("Missing columns:", paste(miss, collapse = ", ")), call. = FALSE)

  if (!("GEN" %in% names(dat))) {
    if ("GCD" %in% names(dat)) dat$GEN <- dat$GCD else stop("Missing GEN (or GCD).", call. = FALSE)
  }

  ec_cols <- grep("^EC[0-9]+$", names(dat), value = TRUE)
  if (!length(ec_cols)) stop("No EC columns found (EC1..ECp).", call. = FALSE)

  dat$Y <- suppressWarnings(as.numeric(dat$Y))
  if (!any(is.finite(dat$Y))) stop("Y has no finite values.", call. = FALSE)

  for (nm in c("Lon", "Lat", ec_cols)) dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))

  # Enforce stable TYP semantics and ordering (EXP first, then CST)
  ty <- unique(as.character(dat$TYP))
  if (!("EXP" %in% ty)) stop("TYP must contain 'EXP'.", call. = FALSE)
  ty_other <- setdiff(ty, "EXP")
  if (length(ty_other) != 1L) stop("TYP must have exactly 2 levels: EXP and CST.", call. = FALSE)
  type_exp <- "EXP"
  type_cst <- ty_other[1]

  dat$GEN <- factor(dat$GEN)
  dat$LOC <- factor(dat$LOC, levels = sort(unique(as.character(dat$LOC))))
  dat$TYP <- factor(as.character(dat$TYP), levels = c(type_exp, type_cst))

  dat <- dat[order(dat$GEN, dat$LOC), ]
  dat <- droplevels(dat)

  ids <- levels(dat$GEN)

  if (is.null(rownames(M_mat))) {
    if (nrow(M_mat) == length(ids)) rownames(M_mat) <- ids
  }
  if (is.null(rownames(M_mat))) stop("M.txt must provide rownames (GCD) or have nrow(M) == nlevels(GEN).", call. = FALSE)
  if (!all(ids %in% rownames(M_mat))) stop("Some GEN levels are missing in M.txt rownames.", call. = FALSE)

  M_mat <- M_mat[ids, , drop = FALSE]
  storage.mode(M_mat) <- "integer"

  A <- try(AGHmatrix::Gmatrix(SNPmatrix = M_mat, method = "VanRaden", ploidy = 2), silent = TRUE)
  if (inherits(A, "try-error")) A <- .vanraden_manual(M_mat)
  A <- .sym(A); diag(A) <- diag(A) + 1e-8
  rownames(A) <- colnames(A) <- ids

  X <- model.matrix(~ 0 + TYP, data = dat)
  colnames(X) <- levels(dat$TYP)

  Zobj <- .build_Z(dat, ec_cols)
  Z <- Zobj$Z
  typ_i <- Zobj$typ_i
  n_re <- Zobj$n_re

  Y <- as.numeric(dat$Y)
  N <- length(Y)
  G <- nlevels(dat$GEN)

  # Initial values from GEN means by TYP (robust to imbalance)
  tab <- aggregate(Y ~ GEN + TYP, data = dat, FUN = mean)
  wide <- reshape(tab, timevar = "TYP", idvar = "GEN", direction = "wide")
  colnames(wide) <- sub("^Y\\.", "", colnames(wide))
  rn <- as.character(wide$GEN); wide$GEN <- NULL; rownames(wide) <- rn

  tlev <- levels(dat$TYP)
  y1 <- wide[, tlev[1], drop = TRUE]
  y2 <- wide[, tlev[2], drop = TRUE]

  s2_1 <- var(y1, use = "complete.obs"); if (!is.finite(s2_1) || s2_1 <= 0) s2_1 <- var(tab$Y)
  s2_2 <- var(y2, use = "complete.obs"); if (!is.finite(s2_2) || s2_2 <= 0) s2_2 <- var(tab$Y)
  s12  <- cov(y1, y2, use = "pairwise.complete.obs"); if (!is.finite(s12)) s12 <- 0
  r0   <- s12 / sqrt(max(1e-12, s2_1 * s2_2)); r0 <- max(min(r0, 0.99), -0.99)

  init_theta <- c(log(s2_1 / 2), log(s2_2 / 2), atanh(r0), log(s2_1 / 2), log(s2_2 / 2))

  reml_negloglik <- function(theta) {
    g11 <- exp(theta[1]); g22 <- exp(theta[2])
    r_gen <- tanh(theta[3])
    g12 <- r_gen * sqrt(g11 * g22)
    r11 <- exp(theta[4]); r22 <- exp(theta[5])

    Ire <- diag(n_re)
    Sigma <- rbind(cbind(g11 * Ire, g12 * Ire),
                   cbind(g12 * Ire, g22 * Ire))
    K <- kronecker(A, Sigma)

    V <- Z %*% K %*% t(Z)
    V[cbind(seq_len(N), seq_len(N))] <- V[cbind(seq_len(N), seq_len(N))] + ifelse(typ_i == 1L, r11, r22) + 1e-10

    cholV <- try(chol(V), silent = TRUE)
    if (inherits(cholV, "try-error")) return(1e10)

    logdetV <- 2 * sum(log(diag(cholV)))

    ViX <- backsolve(cholV, forwardsolve(t(cholV), X))
    ViY <- backsolve(cholV, forwardsolve(t(cholV), Y))

    C <- crossprod(X, ViX)
    cholC <- try(chol(C), silent = TRUE)
    if (inherits(cholC, "try-error")) return(1e10)

    logdetC <- 2 * sum(log(diag(cholC)))
    beta <- try(solve(C, crossprod(X, ViY)), silent = TRUE)
    if (inherits(beta, "try-error")) return(1e10)

    yVy <- sum(Y * ViY)
    bCb <- sum(beta * crossprod(X, ViY))

    ll <- 0.5 * (-logdetV - logdetC - yVy + bCb)
    -ll
  }

  t0 <- Sys.time()
  opt <- nloptr(
    x0 = init_theta,
    eval_f = reml_negloglik,
    lb = c(-20, -20, -10, -20, -20),
    ub = c( 20,  20,  10,  20,  20),
    opts = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = 1e-8, maxeval = 5000)
  )
  t1 <- Sys.time()

  th <- opt$solution
  g11 <- exp(th[1]); g22 <- exp(th[2])
  r_gen <- tanh(th[3]); g12 <- r_gen * sqrt(g11 * g22)
  r11 <- exp(th[4]); r22 <- exp(th[5])

  G0 <- matrix(c(g11, g12, g12, g22), 2, 2, byrow = TRUE)
  R0 <- diag(c(r11, r22)); rownames(R0) <- colnames(R0) <- levels(dat$TYP)

  Ire <- diag(n_re)
  Sigma <- rbind(cbind(g11 * Ire, g12 * Ire),
                 cbind(g12 * Ire, g22 * Ire))
  K <- kronecker(A, Sigma)

  V <- Z %*% K %*% t(Z)
  V[cbind(seq_len(N), seq_len(N))] <- V[cbind(seq_len(N), seq_len(N))] + ifelse(typ_i == 1L, r11, r22) + 1e-10
  cholV <- chol(V)

  ViY <- backsolve(cholV, forwardsolve(t(cholV), Y))
  ViX <- backsolve(cholV, forwardsolve(t(cholV), X))

  XtViX <- crossprod(X, ViX)
  beta <- as.numeric(solve(XtViX, crossprod(X, ViY)))
  names(beta) <- colnames(X)

  resid <- Y - as.vector(X %*% beta)
  Vinv_resid <- backsolve(cholV, forwardsolve(t(cholV), resid))
  u_hat <- as.vector(K %*% crossprod(Z, Vinv_resid))
  y_hat <- as.vector(X %*% beta + Z %*% u_hat)

  u_arr <- array(u_hat, dim = c(n_re, 2L, G))
  blup_df <- data.frame(Gen = rep(levels(dat$GEN), each = 2),
                        Type = rep(levels(dat$TYP), times = G),
                        int = NA_real_, stringsAsFactors = FALSE)
  for (nm in ec_cols) blup_df[[nm]] <- NA_real_
  for (g in seq_len(G)) {
    for (tt in 1:2) {
      rr <- (g - 1L) * 2L + tt
      blup_df$int[rr] <- u_arr[1, tt, g]
      if (length(ec_cols)) blup_df[rr, ec_cols] <- u_arr[2:n_re, tt, g]
    }
  }

  coords <- unique(dat[, c("LOC", "Lon", "Lat")])
  coords <- coords[order(coords$LOC), ]
  rownames(coords) <- NULL

  list(
    dat = dat, ec_cols = ec_cols,
    type_exp = type_exp, type_cst = type_cst,
    M = M_mat, A = A, X = X, Z = Z,
    theta = th, G0 = G0, R0 = R0, r_gen = r_gen,
    beta = beta, u_hat = u_hat, blup_df = blup_df, y_hat = y_hat,
    logLik_REML = -opt$objective,
    elapsed_sec = as.numeric(difftime(t1, t0, units = "secs")),
    iterations = if (!is.null(opt$iterations)) opt$iterations else NA,
    coords = coords
  )
}

# ============================================================
# 2) Application (inputs: phenos_and_ECs.txt, M.txt)
# ============================================================

dat <- read.table("phenos_and_ECs.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
M_mat <- .read_M("M.txt")

fit <- .fit_dfreml_bitra(dat, M_mat)

cat("Elapsed time (sec):", round(fit$elapsed_sec, 3), "\n")
cat("REML log-likelihood:", fit$logLik_REML, "\n")
cat("Iterations:", fit$iterations, "\n\n")

cat("Estimated G0:\n"); print(fit$G0)
cat("\nEstimated R0:\n"); print(fit$R0)
cat("\nGenetic correlation:", round(fit$r_gen, 6), "\n\n")

yhat <- as.numeric(fit$y_hat)
Yobs <- as.numeric(fit$dat$Y)

cor_obs_pred <- stats::cor(Yobs, yhat, use = "complete.obs")
rmse_obs_pred <- sqrt(mean((Yobs - yhat)^2, na.rm = TRUE))
cat("Obs vs Pred (all data): r =", round(cor_obs_pred, 4), "| RMSE =", round(rmse_obs_pred, 4), "\n\n")

tmp <- fit$dat
tmp$y_hat <- yhat
pred_tab <- aggregate(y_hat ~ GEN + TYP, data = tmp, FUN = mean)
pred_wide <- reshape(pred_tab, timevar = "TYP", idvar = "GEN", direction = "wide")
colnames(pred_wide) <- sub("^y_hat\\.", "", colnames(pred_wide))
rownames(pred_wide) <- as.character(pred_wide$GEN)
pred_wide$GEN <- NULL

cat("Predicted means by GEN and TYP:\n"); print(pred_wide)
cat("\nCovariance of predicted means:\n"); print(cov(pred_wide, use = "complete.obs"))
cat("\nCorrelation of predicted means:\n"); print(cor(pred_wide, use = "complete.obs"))

# ============================================================
# 3) Visualization (EXP first, then CST; GEN predicted as CST)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(viridis)
  library(ggdendro)
  library(patchwork)
})

A <- fit$A
dat2 <- fit$dat
coords <- fit$coords
ec_cols <- fit$ec_cols
tlev <- levels(dat2$TYP)

# GEN color palette (match OLD.R figure A/B)
gen_colors <- c(
  "GEN1" = "#00CFFF",
  "GEN2" = "#4B83F5",
  "GEN3" = "#C67BF7",
  "GEN4" = "#FF79C6",
  "GEN5" = "#FF5E7C"
)

# Keep EXP first, then CST/CSL in all faceting/plots
if ("EXP" %in% tlev) {
  other <- setdiff(tlev, "EXP")
  if (length(other) == 1L) dat2$TYP <- factor(as.character(dat2$TYP), levels = c("EXP", other))
}
tlev <- levels(dat2$TYP)

# --- A heatmap + dendrogram ---
d <- as.dist(1 - A)
hc <- hclust(d, method = "ward.D2")
ord <- hc$order
row_order <- rownames(A)[ord]
col_order <- colnames(A)[ord]

dend <- ggdendro::dendro_data(hc)
p_dendro <- ggplot(ggdendro::segment(dend)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  theme_void()

A_melt <- as.data.frame(as.table(A))
names(A_melt) <- c("Genotype1", "Genotype2", "Value")

p_heatmap <- ggplot(
  A_melt,
  aes(x = factor(Genotype1, levels = col_order),
      y = factor(Genotype2, levels = row_order),
      fill = Value)
) +
  geom_tile(color = "grey90") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "Genomic\nRelationship") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.title = element_blank())

print(p_dendro / p_heatmap + plot_layout(heights = c(1, 4)))

# --- Observed vs predicted (facet order follows TYP levels: EXP first) ---
obs_pred <- data.frame(
  GEN = dat2$GEN, TYP = dat2$TYP,
  Y = as.numeric(dat2$Y), y_hat = yhat,
  stringsAsFactors = FALSE
)
obs_pred$GROUP <- ifelse(grepl("^GEN", as.character(obs_pred$GEN)), "Breeding Genotype", "Commercial Cultivar")

cors <- tapply(seq_len(nrow(obs_pred)), obs_pred$TYP,
               function(ii) cor(obs_pred$Y[ii], obs_pred$y_hat[ii], use = "complete.obs"))
cors_df <- data.frame(TYP = names(cors), r = as.numeric(cors), stringsAsFactors = FALSE)
cors_df$TYP <- factor(cors_df$TYP, levels = tlev)
cors_df$r_label <- paste0("r = ", round(cors_df$r, 3))

print(
  ggplot(obs_pred, aes(x = y_hat, y = Y)) +
    geom_smooth(method = "lm", color = "black", se = FALSE, linetype = 2) +
    geom_point(aes(color = GEN, shape = GROUP), size = 2) +
    scale_color_manual(values = gen_colors, na.value = "grey60") +
    scale_shape_manual(values = c("Breeding Genotype" = 16, "Commercial Cultivar" = 17)) +
    facet_wrap(~ TYP, scales = "free_x") +
    geom_text(data = cors_df, aes(x = -Inf, y = Inf, label = r_label),
              hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 4) +
    labs(x = "Predicted", y = "Observed", shape = NULL) +
    theme_bw()
)

# --- GEN predicted as CST (never EXP) across all LOC ---
env_cov <- unique(dat2[, c("LOC", ec_cols), drop = FALSE])
env_cov <- env_cov[order(env_cov$LOC), ]

type_stand <- fit$type_cst
beta_stand <- as.numeric(fit$beta[type_stand])

gen_focus <- grep("^GEN", levels(dat2$GEN), value = TRUE)
if (!length(gen_focus)) gen_focus <- levels(dat2$GEN)
if (length(gen_focus) > 5L) gen_focus <- gen_focus[1:5]

# Restrict palette to the GENs actually present (keeps mapping stable)
gen_colors_use <- gen_colors[intersect(names(gen_colors), gen_focus)]
if (!length(gen_colors_use)) gen_colors_use <- gen_colors

key_all <- paste(fit$blup_df$Gen, fit$blup_df$Type)
idx_bt <- match(paste(gen_focus, type_stand), key_all)
blup_stand <- fit$blup_df[idx_bt, , drop = FALSE]
rownames(blup_stand) <- blup_stand$Gen

loc_levels <- levels(dat2$LOC)
pred_grid <- expand.grid(GEN = gen_focus, LOC = loc_levels, stringsAsFactors = FALSE)
pred_grid <- merge(pred_grid, env_cov, by = "LOC", all.x = TRUE, sort = FALSE)
pred_grid$LOC <- factor(pred_grid$LOC, levels = loc_levels)
pred_grid <- pred_grid[order(pred_grid$GEN, pred_grid$LOC), ]

Bint <- blup_stand[as.character(pred_grid$GEN), "int", drop = TRUE]
Bslp <- as.matrix(blup_stand[as.character(pred_grid$GEN), ec_cols, drop = FALSE])
Ecv  <- as.matrix(pred_grid[, ec_cols, drop = FALSE])
pred_grid$y_hat <- beta_stand + Bint + rowSums(Bslp * Ecv)

label_start <- pred_grid[pred_grid$LOC == loc_levels[1], ]
label_end   <- pred_grid[pred_grid$LOC == loc_levels[length(loc_levels)], ]

print(
  ggplot(pred_grid, aes(x = LOC, y = y_hat, group = GEN, color = GEN)) +
    geom_point(size = 2) +
    geom_line() +
    geom_text(data = label_start, aes(label = GEN), hjust = 1.2, vjust = 0.5, size = 3.5, show.legend = FALSE) +
    geom_text(data = label_end,   aes(label = GEN), hjust = -0.2, vjust = 0.5, size = 3.5, show.legend = FALSE) +
    scale_color_manual(values = gen_colors_use, na.value = "grey60") +
    scale_x_discrete(expand = expansion(mult = c(0.125, 0.125))) +
    labs(x = "Location (LOC)", y = paste0("Predicted as ", type_stand)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.minor = element_blank(),
          legend.position = "none")
)

# --- LOC coordinates ---
print(
  ggplot(coords, aes(x = Lon, y = Lat, label = LOC)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, size = 3) +
    scale_x_continuous(breaks = 1:10, limits = c(0.5, 10.5)) +
    scale_y_continuous(breaks = 1:10, limits = c(0.5, 10.5)) +
    coord_fixed() +
    theme_bw()
)

# --- PC1 map of EC space ---
coords_env <- merge(coords, env_cov, by = "LOC", all.x = TRUE, sort = FALSE)
coords_env <- coords_env[order(coords_env$LOC), ]

pca <- prcomp(coords_env[, ec_cols, drop = FALSE], scale. = TRUE)
coords_env$PC1 <- pca$x[, 1]

lo_pc1 <- loess(PC1 ~ Lon * Lat, data = coords_env, span = 0.75)
grid <- expand.grid(
  Lon = seq(min(coords_env$Lon), max(coords_env$Lon), by = 0.25),
  Lat = seq(min(coords_env$Lat), max(coords_env$Lat), by = 0.25)
)
grid$PC1_pred <- as.numeric(predict(lo_pc1, newdata = grid))
grid <- grid[is.finite(grid$PC1_pred), , drop = FALSE]

print(
  ggplot(grid, aes(x = Lon, y = Lat, fill = PC1_pred)) +
    geom_raster() +
    geom_contour(aes(z = PC1_pred), bins = 10, color = "white") +
    geom_point(data = coords_env, aes(x = Lon, y = Lat), color = "black", size = 2, inherit.aes = FALSE) +
    geom_text(data = coords_env, aes(x = Lon, y = Lat, label = LOC), vjust = -0.5, size = 3, inherit.aes = FALSE) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = 1:10) +
    scale_fill_viridis_c(name = "PC1") +
    coord_fixed() +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw()
)

# --- EC maps (loess on coords) ---
gridEC <- grid[, c("Lon", "Lat"), drop = FALSE]
for (ec in ec_cols) {
  lo <- loess(as.formula(paste(ec, "~ Lon * Lat")), data = coords_env, span = 0.75)
  gridEC[[paste0(ec, "_pred")]] <- as.numeric(predict(lo, newdata = gridEC))
}

grid_long <- do.call(rbind, lapply(ec_cols, function(ec) {
  data.frame(Lon = gridEC$Lon, Lat = gridEC$Lat,
             Panel = paste0(ec, "_pred"),
             pred = gridEC[[paste0(ec, "_pred")]],
             stringsAsFactors = FALSE)
}))
empty_panel <- data.frame(Lon = coords_env$Lon, Lat = coords_env$Lat,
                          Panel = "LOC", pred = NA_real_, stringsAsFactors = FALSE)

plot_data <- rbind(empty_panel, grid_long)
plot_data$Panel <- factor(plot_data$Panel, levels = c("LOC", paste0(ec_cols, "_pred")))

print(
  ggplot(plot_data, aes(x = Lon, y = Lat)) +
    facet_wrap(~ Panel, ncol = 3) +
    geom_raster(data = subset(plot_data, Panel != "LOC"), aes(fill = pred)) +
    geom_point(data = coords_env, aes(x = Lon, y = Lat), color = "black", size = 1, inherit.aes = FALSE) +
    geom_text(data = coords_env, aes(x = Lon, y = Lat, label = LOC), vjust = 0, size = 2, inherit.aes = FALSE) +
    scale_fill_viridis_c(na.value = "transparent", name = "EC value") +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = 1:10) +
    coord_fixed() +
    labs(x = "Longitude", y = "Latitude") +
    theme_bw()
)

# --- Recommendation map (winner GEN) for each TYP; facet order EXP then CST ---
ec_pred_cols <- paste0(ec_cols, "_pred")

pred_idx <- expand.grid(
  Gen = gen_focus,
  Type = levels(dat2$TYP),
  idx = seq_len(nrow(gridEC)),
  stringsAsFactors = FALSE
)
pred_data <- cbind(pred_idx, gridEC[pred_idx$idx, c("Lon", "Lat", ec_pred_cols), drop = FALSE])

key_sub <- paste(fit$blup_df$Gen, fit$blup_df$Type)
mpos <- match(paste(pred_data$Gen, pred_data$Type), key_sub)

coef_int <- fit$blup_df$int[mpos]
coef_slp <- as.matrix(fit$blup_df[mpos, ec_cols, drop = FALSE])

Epred <- as.matrix(pred_data[, ec_pred_cols, drop = FALSE])
colnames(Epred) <- sub("_pred$", "", colnames(Epred))
Epred <- Epred[, ec_cols, drop = FALSE]

pred_data$beta <- as.numeric(fit$beta[pred_data$Type])
pred_data$y_hat <- pred_data$beta + coef_int + rowSums(coef_slp * Epred)

grp <- interaction(pred_data$Type, pred_data$idx, drop = TRUE)
imax <- tapply(seq_len(nrow(pred_data)), grp, function(ii) ii[which.max(pred_data$y_hat[ii])])
winner <- pred_data[as.integer(imax), c("Type", "Lon", "Lat", "Gen", "y_hat")]
winner$Type <- factor(winner$Type, levels = levels(dat2$TYP))

print(
  ggplot(winner, aes(x = Lon, y = Lat, fill = Gen)) +
    geom_tile() +
    facet_wrap(~ Type) +
    scale_fill_manual(values = gen_colors_use, na.value = "grey80") +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = 1:10) +
    coord_fixed() +
    labs(x = "Longitude", y = "Latitude", fill = "Winning GEN") +
    theme_bw()
)
