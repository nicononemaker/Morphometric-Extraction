# =============================================================================
# PERMANOVA  —  PART 1  (single file: analysis, all sections)
#
# Self-contained. Config + Sections 1-7 in order, single-threaded.
# Trait set uses stipe_to_length_ratio (logged); CO2 predictor is spco2_clim.
# Runs region PERMANOVA/PERMDISP/pairwise (1), environmental PERMANOVA (2),
# region-vs-env varpart (3), nested ordination (4), residual ordination (5),
# R2-by-scale trend (6), and distance-decay/Mantel (7). Writes all *.rds caches
# that permanova_2_figures.R reads.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(vegan) })

# ---------------------------------------------------------------- config ------
# --- locate and source config.R (works with the .Rproj open, via source(), or
# from the console) --------------------------------------------------------
.find_config <- function() {
  for (p in c("config.R", "../config.R")) if (file.exists(p)) return(p)
  # last resort: derive from this script's own path (RStudio or source())
  sp <- if (requireNamespace("rstudioapi", quietly = TRUE))
    tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "") else ""
  if (!nzchar(sp)) sp <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) "")
  if (nzchar(sp)) {
    cand <- file.path(dirname(dirname(sp)), "config.R")   # repo root from analysis/
    if (file.exists(cand)) return(cand)
  }
  stop("config.R not found. Open Morphometric-Extraction.Rproj, or set your ",
       "working directory to the repo root, then re-run.")
}
source(.find_config())
IN    <- AGG_CSV
OUT   <- OUT_DIR
CACHE <- file.path(OUT, "cache");        dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

TRAITS <- c("length_frond_mm","max_width_mm","stipe_to_length_ratio",
            "solidity","peak_position_frac","branched_fraction","circularity")
LOG    <- c("length_frond_mm","max_width_mm","stipe_to_length_ratio")

PREDICTORS <- c("SST_clim","nitrate_clim","phosphate_clim","salinity_clim",
                "wave.height","kd490_clim","PAR_mean","iron_clim","spco2_clim")

REGIONS <- c("REALM","PROVINC","ECOREGION")
MIN_N   <- 10
NPERM   <- 9999

# Parallel permutations for varpart adonis2 fits. Leaves 2 logical cores free
# for the OS. Set to 1 to force single-threaded (exact reproducibility of the
# non-parallel seed stream). On Windows this uses a socket cluster: real but
# sub-linear speedup (per-fit overhead + data copy to workers).
NCORES  <- max(1L, parallel::detectCores() - 2L)

LAT_CANDS <- c("lat","latitude","decimalLatitude","Latitude","LAT","y")
LON_CANDS <- c("lon","long","longitude","decimalLongitude","Longitude","LON","x")

say <- function(...) cat(..., "\n")

# ---------------------------------------------------------------- load --------
rec <- read.csv(IN, stringsAsFactors = FALSE, check.names = FALSE)
TRAITS     <- intersect(TRAITS, names(rec))
PREDICTORS <- intersect(PREDICTORS, names(rec))
REGIONS    <- intersect(REGIONS, names(rec))
LAT <- LAT_CANDS[LAT_CANDS %in% names(rec)][1]
LON <- LON_CANDS[LON_CANDS %in% names(rec)][1]
say("Traits    :", paste(TRAITS, collapse = ", "))
say("Predictors:", paste(PREDICTORS, collapse = ", "))
say("Coords    :", if (!is.na(LAT) && !is.na(LON)) paste(LON, "/", LAT) else "none found")

# shared collinearity pruning (|r|>0.9, priority order)
prune_predictors <- function(dat, preds) {
  cm <- cor(dat[, preds, drop = FALSE], use = "pairwise.complete.obs")
  kept <- character(0)
  for (p in preds) {
    if (!length(kept)) { kept <- p; next }
    if (max(abs(cm[p, kept])) > 0.9) {
      hit <- kept[which.max(abs(cm[p, kept]))]
      say(sprintf("  drop %-16s (|r|=%.2f with %s)", p, max(abs(cm[p, kept])), hit))
    } else kept <- c(kept, p)
  }
  kept
}

# centroid distance matrix D[C] in the space of the dissimilarity measure
# (Anderson et al. §5.1 back-transform via the averaged Gower matrix).
# Works for any dissimilarity; for Euclidean it reduces to distances among means.
centroid_dist <- function(D, g) {
  g   <- factor(g)
  lev <- levels(g)
  # Gower-centered inner-product matrix G from D
  A  <- -0.5 * as.matrix(D)^2
  n  <- nrow(A)
  J  <- diag(n) - matrix(1/n, n, n)
  G  <- J %*% A %*% J
  # averaged (g x g) Gower matrix
  Gc <- matrix(NA_real_, length(lev), length(lev), dimnames = list(lev, lev))
  for (a in seq_along(lev)) for (b in seq_along(lev)) {
    ia <- which(g == lev[a]); ib <- which(g == lev[b])
    Gc[a, b] <- mean(G[ia, ib])
  }
  # back-transform to distances among centroids
  DC <- matrix(0, length(lev), length(lev), dimnames = list(lev, lev))
  for (a in seq_along(lev)) for (b in seq_along(lev)) {
    val <- Gc[a, a] - 2 * Gc[a, b] + Gc[b, b]
    DC[a, b] <- sqrt(max(val, 0))
  }
  as.dist(DC)
}

# =============================================================================
# SECTION 1  region PERMANOVA / PERMDISP / pairwise  (cached per scale)
# =============================================================================
analyze_region <- function(GROUP) {
  say(sprintf("\n########## ANALYZE REGION  %s ##########", GROUP))
  d <- rec
  for (t in LOG) d[[t]] <- log(d[[t]])
  d[[GROUP]][trimws(tolower(d[[GROUP]])) %in% c("other","","na","unknown")] <- NA
  cc <- complete.cases(d[, TRAITS], d[[GROUP]]); d <- d[cc, ]
  keep <- names(which(table(d[[GROUP]]) >= MIN_N)); d <- d[d[[GROUP]] %in% keep, ]
  if (length(keep) < 3 || nrow(d) < 30) { say("  too few; skipped"); return(invisible()) }

  g <- factor(d[[GROUP]])
  D <- dist(scale(d[, TRAITS]), method = "euclidean")
  say(sprintf("  %d records, %d groups", nrow(d), nlevels(g)))

  set.seed(1); ad <- adonis2(D ~ g, permutations = NPERM)
  bd_test <- betadisper(D, g, bias.adjust = TRUE)
  set.seed(1); pd <- permutest(bd_test, permutations = NPERM, pairwise = TRUE)

  # pairwise PERMANOVA (Holm)
  M <- as.matrix(D); combs <- combn(levels(g), 2, simplify = FALSE)
  pw <- do.call(rbind, lapply(combs, function(p) {
    sel <- g %in% p; set.seed(1)
    a <- adonis2(as.dist(M[sel, sel]) ~ g[sel], permutations = NPERM)
    data.frame(r1 = p[1], r2 = p[2], R2 = round(a$R2[1], 4), P = a$`Pr(>F)`[1])
  }))
  pw$P_adj <- p.adjust(pw$P, "holm"); pw$sig <- pw$P_adj < 0.05

  # centroid coords (trait-space PCoA) via plain betadisper
  bd  <- betadisper(D, g)
  cen <- as.data.frame(bd$centroids[, 1:2]); names(cen) <- c("Ax1","Ax2")
  cen$grp <- rownames(cen); cen$n <- as.integer(table(g)[cen$grp])
  eig <- bd$eig; ve <- 100 * eig[1:2] / sum(eig[eig > 0])
  pts <- as.data.frame(bd$vectors[, 1:2]); names(pts) <- c("Ax1","Ax2"); pts$grp <- g

  # centroid distance matrix D[C] (§5.1) + its mMDS coords for the inset
  DC <- centroid_dist(D, g)
  mds <- tryCatch(vegan::metaMDS(DC, k = 2, trace = 0, autotransform = FALSE),
                  error = function(e) NULL)
  mds_pts <- if (!is.null(mds)) {
    x <- as.data.frame(vegan::scores(mds, display = "sites"))
    names(x) <- c("MDS1","MDS2"); x$grp <- rownames(x); x
  } else NULL
  mds_stress <- if (!is.null(mds)) mds$stress else NA_real_

  # geographic centroids for the map (if coords exist)
  geo <- NULL; specimens_xy <- NULL
  if (!is.na(LAT) && !is.na(LON)) {
    dd <- d %>% mutate(.grp = g,
                       .lat = suppressWarnings(as.numeric(.data[[LAT]])),
                       .lon = suppressWarnings(as.numeric(.data[[LON]]))) %>%
      filter(is.finite(.lat), is.finite(.lon))
    geo <- dd %>% group_by(.grp) %>%
      summarise(lon = mean(.lon), lat = mean(.lat), n = dplyr::n(), .groups = "drop")
    specimens_xy <- dd %>% dplyr::select(.grp, .lon, .lat)
  }

  res <- list(
    GROUP = GROUP, n = nrow(d), levels = levels(g),
    adonis = ad, permdisp = pd, pairwise = pw,
    centroids_trait = cen, specimens_trait = pts, var_expl = ve,
    DC = DC, mds_pts = mds_pts, mds_stress = mds_stress,
    geo = geo, specimens_xy = specimens_xy,
    summary = data.frame(grouping = GROUP, n = nrow(d), groups = nlevels(g),
                         R2 = round(ad$R2[1], 4), F = round(ad$F[1], 2),
                         adonis_p = ad$`Pr(>F)`[1],
                         permdisp_p = round(pd$tab$`Pr(>F)`[1], 4),
                         pairs_sig = sum(pw$sig), pairs_total = nrow(pw)))
  saveRDS(res, file.path(CACHE, sprintf("region_%s.rds", tolower(GROUP))))
  say(sprintf("  cached region_%s.rds  (pairwise %d/%d sig; mMDS stress %.3f)",
              tolower(GROUP), sum(pw$sig), nrow(pw), mds_stress))
  res$summary
}

region_summ <- bind_rows(lapply(REGIONS, analyze_region))
saveRDS(region_summ, file.path(CACHE, "region_summary.rds"))


# =============================================================================
# SECTION 2  environmental PERMANOVA (marginal)  — cached
# =============================================================================
say("\n########## ANALYZE ENVIRONMENT ##########")
dE <- rec
for (t in LOG) dE[[t]] <- log(dE[[t]])
ccE <- complete.cases(dE[, TRAITS], dE[, PREDICTORS]); dE <- dE[ccE, ]
say(sprintf("  complete-case records: %d", nrow(dE)))
keepE <- prune_predictors(dE, PREDICTORS)
say("  surviving:", paste(keepE, collapse = ", "))
DE <- dist(scale(dE[, TRAITS]), method = "euclidean")
fE <- as.formula(paste("DE ~", paste(keepE, collapse = " + ")))

# marginal (unique) contribution of each predictor
set.seed(1); adE <- adonis2(fE, data = dE, by = "margin", permutations = NPERM)

# #3  TOTAL model R2 (whole formula) so we can show the shared gap.
#     total R2 - sum(marginal R2) = variance shared among collinear predictors,
#     which the marginal bars cannot attribute to any single predictor.
set.seed(1); adE_all <- adonis2(fE, data = dE, permutations = NPERM)   # by="terms" default
rn_all <- rownames(adE_all)
res_row <- grep("^Residual$", rn_all, ignore.case = TRUE)
R2_total <- if (length(res_row) == 1 && "R2" %in% names(adE_all)) {
  1 - adE_all$R2[res_row]                      # total explained = 1 - residual R2
} else NA_real_
R2_marg_sum <- sum(as.data.frame(adE)$R2[rownames(adE) %in% keepE])
R2_shared  <- if (is.finite(R2_total)) max(R2_total - R2_marg_sum, 0) else NA_real_
say(sprintf("  total model R2 = %s | sum of marginal = %.3f | shared = %s",
            ifelse(is.finite(R2_total), sprintf("%.3f", R2_total), "NA"),
            R2_marg_sum,
            ifelse(is.finite(R2_shared), sprintf("%.3f", R2_shared), "NA")))

# #5  bootstrap CIs on each predictor's MARGINAL R2, via DIRECT matrix algebra.
#     (adonis2(by="margin", permutations=0) per rep is pathologically slow; this
#     computes marginal R2 straight from the Gower-centered matrix instead.)
#       G = -1/2 J D^2 J ; SST = tr(G) ; H = X(X'X)^-1 X' ; SSR = tr(HGH)
#       marginal R2_j = (SSR_full - SSR_without_j) / SST
marg_R2 <- function(Dmat, Xfull) {
  n  <- attr(Dmat, "Size"); Dm <- as.matrix(Dmat)
  J  <- diag(n) - 1/n
  G  <- -0.5 * J %*% (Dm^2) %*% J
  SST <- sum(diag(G))
  # tr(HGH) = tr(GH) = tr( (X'G X)(X'X)^-1 )  -- only k x k matrices, no n x n hat
  trGH <- function(X) {
    X <- cbind(1, X)
    sum(diag(crossprod(X, G %*% X) %*% solve(crossprod(X))))
  }
  ssr_full <- trGH(Xfull)
  vapply(seq_len(ncol(Xfull)), function(j)
    (ssr_full - trGH(Xfull[, -j, drop = FALSE])) / SST, numeric(1))
}
NBOOT <- 999
say(sprintf("  bootstrapping marginal-R2 CIs (%d reps, matrix method)...", NBOOT))
boot_mat <- matrix(NA_real_, nrow = NBOOT, ncol = length(keepE),
                   dimnames = list(NULL, keepE))
set.seed(1); t0 <- Sys.time()
for (b in seq_len(NBOOT)) {
  idx <- sample.int(nrow(dE), replace = TRUE)
  db  <- dE[idx, , drop = FALSE]
  Db  <- dist(scale(db[, TRAITS]), method = "euclidean")
  Xb  <- as.matrix(db[, keepE, drop = FALSE])
  boot_mat[b, ] <- tryCatch(marg_R2(Db, Xb), error = function(e) rep(NA_real_, length(keepE)))
}
say(sprintf("  bootstrap CIs done in %.1f s",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
ci <- t(apply(boot_mat, 2, function(x)
  quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)))
colnames(ci) <- c("lo", "med", "hi")
ci_df <- data.frame(term = keepE, ci, row.names = NULL)

env_res <- list(kept = keepE, adonis = adE,
                tab = { t <- as.data.frame(adE); t$term <- rownames(t); t },
                R2_total = R2_total, R2_marg_sum = R2_marg_sum, R2_shared = R2_shared,
                ci = ci_df, nboot = NBOOT)
saveRDS(env_res, file.path(CACHE, "env.rds"))
say("  cached env.rds (with total R2, shared gap, bootstrap CIs)")

# sequential (incremental) R2 decomposition for Option 3 (deterministic)
{
  fitR2 <- function(Dmat, X) {
    n <- attr(Dmat, "Size"); Dm <- as.matrix(Dmat)
    J <- diag(n) - 1/n; G <- -0.5 * J %*% (Dm^2) %*% J
    SST <- sum(diag(G)); Xi <- cbind(1, X)
    sum(diag(crossprod(Xi, G %*% Xi) %*% solve(crossprod(Xi)))) / SST
  }
  marg <- as.data.frame(adE); marg$term <- rownames(marg)
  ord  <- marg$term[order(-marg$R2)]; ord <- ord[ord %in% keepE]
  cumv <- numeric(length(ord))
  for (i in seq_along(ord)) cumv[i] <- fitR2(DE, as.matrix(dE[, ord[1:i], drop = FALSE]))
  sq <- data.frame(term = ord, seq_R2 = round(c(cumv[1], diff(cumv)), 4),
                   cum_R2 = round(cumv, 4))
  saveRDS(sq, file.path(CACHE, "env_sequential.rds"))
  say("  cached env_sequential.rds"); print(sq, row.names = FALSE)
}


# =============================================================================
# SECTION 3  variance partitioning region vs environment  — cached per scale
# =============================================================================
analyze_varpart <- function(GROUP) {
  say(sprintf("\n########## ANALYZE VARPART  %s ##########", GROUP))
  d <- rec
  for (t in LOG) d[[t]] <- log(d[[t]])
  d[[GROUP]][trimws(tolower(d[[GROUP]])) %in% c("other","","na","unknown")] <- NA
  cc <- complete.cases(d[, TRAITS], d[, PREDICTORS]) & !is.na(d[[GROUP]]); d <- d[cc, ]
  keep <- names(which(table(d[[GROUP]]) >= MIN_N)); d <- d[d[[GROUP]] %in% keep, ]
  d[[GROUP]] <- factor(d[[GROUP]])
  if (nrow(d) < 30 || nlevels(d[[GROUP]]) < 3) { say("  too few; skipped"); return(invisible()) }
  kept <- prune_predictors(d, PREDICTORS)

  D <- dist(scale(d[, TRAITS]), method = "euclidean")
  env_rhs <- paste(kept, collapse = " + ")
  f_E  <- as.formula(paste("D ~", env_rhs))
  f_R  <- as.formula(paste("D ~", GROUP))
  f_ER <- as.formula(paste("D ~", env_rhs, "+", GROUP))
  f_RE <- as.formula(paste("D ~", GROUP, "+", env_rhs))

  set.seed(1); aE  <- adonis2(f_E,  data = d, by = "terms",  permutations = NPERM)
  set.seed(1); aR  <- adonis2(f_R,  data = d, by = "terms",  permutations = NPERM)
  set.seed(1); aRE <- adonis2(f_RE, data = d, by = "terms",  permutations = NPERM)
  set.seed(1); aER <- adonis2(f_ER, data = d, by = "terms",  permutations = NPERM)
  set.seed(1); aRE_margin <- adonis2(f_RE, data = d, by = "margin", permutations = NPERM)

  R2_E        <- sum(aE$R2[rownames(aE) %in% kept])
  R2_R        <- aR$R2[rownames(aR) == GROUP]
  R2_R_give_E <- aER$R2[rownames(aER) == GROUP]
  R2_E_give_R <- sum(aRE$R2[rownames(aRE) %in% kept])
  shared      <- R2_R - R2_R_give_E
  resid_R     <- R2_R_give_E
  pure_E      <- R2_E_give_R
  total_expl  <- pure_E + shared + resid_R
  p_R_give_E  <- aRE_margin$`Pr(>F)`[rownames(aRE_margin) == GROUP]
  frac_abs    <- if (R2_R > 0) shared / R2_R else NA

  part <- data.frame(
    scale = GROUP, region_alone = round(R2_R,4), env_alone = round(R2_E,4),
    pure_env = round(pure_E,4), shared = round(shared,4),
    pure_region = round(resid_R,4), region_given_env_p = p_R_give_E,
    total_explained = round(total_expl,4), frac_region_absorbed = round(frac_abs,3))
  res <- list(GROUP = GROUP, kept = kept,
              R2_E = R2_E, R2_R = R2_R, pure_E = pure_E, shared = shared,
              resid_R = resid_R, total_expl = total_expl,
              p_R_give_E = p_R_give_E, frac_absorbed = frac_abs, part = part)
  saveRDS(res, file.path(CACHE, sprintf("varpart_%s.rds", tolower(GROUP))))
  say(sprintf("  cached varpart_%s.rds  (%.0f%% of region absorbed by env)",
              tolower(GROUP), 100 * frac_abs))
  part
}

varpart_summ <- bind_rows(lapply(REGIONS, analyze_varpart))
saveRDS(varpart_summ, file.path(CACHE, "varpart_summary.rds"))


# =============================================================================
# ||  EXPLORATORY — NOT USED IN THE MANUSCRIPT  ||
# ||  Sections 4-7 below (nested ordination, residual ordination, R2-by-scale
# ||  trend, distance-decay/Mantel) were exploratory. They are NOT the source
# ||  of any manuscript figure or table and each writes its own separate .rds
# ||  cache. To trim the repo to paper-only outputs, delete Sections 4-7 here
# ||  AND the matching figure blocks in 04_permanova_figures.R that read
# ||  nested_ordination.rds / residual_ordination.rds / r2_trend.rds /
# ||  distdecay_*.rds. Left runnable so nothing breaks until you remove both.
# =============================================================================

# =============================================================================
# SECTION 4  NESTED main-effects ordination (Anderson Image-2 style)
#   All region scales' centroids embedded in ONE shared trait-space mMDS, with
#   parent<-child nesting links (ecoregion -> province -> realm). Uses the SAME
#   complete-case frame across scales so the embedding is common.
# =============================================================================
say("\n########## NESTED MAIN-EFFECTS ORDINATION ##########")
build_nested <- function() {
  # need all three scales present and a nesting relationship in the data
  scales <- intersect(c("REALM","PROVINC","ECOREGION"), REGIONS)
  if (length(scales) < 2) { say("  <2 scales available; skipped"); return(NULL) }

  d <- rec
  for (t in LOG) d[[t]] <- log(d[[t]])
  # one shared complete-case frame: complete traits AND all scale labels present
  ok <- complete.cases(d[, TRAITS])
  for (s in scales) {
    d[[s]][trimws(tolower(d[[s]])) %in% c("other","","na","unknown")] <- NA
    ok <- ok & !is.na(d[[s]])
  }
  d <- d[ok, ]
  if (nrow(d) < 30) { say("  too few shared complete cases; skipped"); return(NULL) }
  say(sprintf("  shared frame: %d records across %d scales", nrow(d), length(scales)))

  D <- dist(scale(d[, TRAITS]), method = "euclidean")

  # centroid coords per scale, each in its OWN D[C] space, then embed ALL
  # centroids jointly: build one combined centroid-distance matrix across every
  # level of every scale so they live in a common ordination.
  # Approach: stack a grouping that is unique per (scale,level), get D[C] once.
  combo <- factor(do.call(paste, c(lapply(scales, function(s) paste0(substr(s,1,3), ":", d[[s]])),
                                   sep = "___")))
  # but we want per-level centroids, not per-combination; instead compute each
  # scale's centroids in the shared G and assemble.
  A  <- -0.5 * as.matrix(D)^2; n <- nrow(A)
  J  <- diag(n) - matrix(1/n, n, n); G <- J %*% A %*% J

  cen_all <- list(); map <- list()
  for (s in scales) {
    gs <- factor(d[[s]]); lev <- levels(gs)
    # centroid inner products within shared G
    Gc <- matrix(NA_real_, length(lev), length(lev), dimnames = list(lev, lev))
    for (a in seq_along(lev)) for (b in seq_along(lev)) {
      ia <- which(gs == lev[a]); ib <- which(gs == lev[b])
      Gc[a, b] <- mean(G[ia, ib])
    }
    cen_all[[s]] <- Gc
  }

  # Assemble a global among-centroid distance matrix over ALL levels of ALL
  # scales using shared-G inner products between every centroid pair.
  all_levels <- do.call(rbind, lapply(scales, function(s) {
    data.frame(scale = s, level = levels(factor(d[[s]])), stringsAsFactors = FALSE)
  }))
  K <- nrow(all_levels)
  # centroid inner product between (scale a, level i) and (scale b, level j)
  cip <- function(sa, la, sb, lb) {
    ia <- which(d[[sa]] == la); ib <- which(d[[sb]] == lb)
    mean(G[ia, ib])
  }
  Gg <- matrix(NA_real_, K, K)
  for (i in seq_len(K)) for (j in seq_len(K))
    Gg[i, j] <- cip(all_levels$scale[i], all_levels$level[i],
                    all_levels$scale[j], all_levels$level[j])
  DCg <- matrix(0, K, K)
  for (i in seq_len(K)) for (j in seq_len(K))
    DCg[i, j] <- sqrt(max(Gg[i,i] - 2*Gg[i,j] + Gg[j,j], 0))
  mds <- tryCatch(vegan::metaMDS(as.dist(DCg), k = 2, trace = 0, autotransform = FALSE),
                  error = function(e) NULL)
  if (is.null(mds)) { say("  mMDS failed; skipped"); return(NULL) }
  xy <- as.data.frame(vegan::scores(mds, display = "sites"))
  names(xy) <- c("MDS1","MDS2")
  coords <- cbind(all_levels, xy)
  coords$n <- vapply(seq_len(K), function(i)
    sum(d[[all_levels$scale[i]]] == all_levels$level[i]), integer(1))

  # nesting links: for each finer-scale level, find its (modal) parent level
  # at the next-coarser scale, so we can draw child->parent segments.
  order_coarse_to_fine <- intersect(c("REALM","PROVINC","ECOREGION"), scales)
  links <- list()
  if (length(order_coarse_to_fine) >= 2) {
    for (k in 2:length(order_coarse_to_fine)) {
      child <- order_coarse_to_fine[k]; parent <- order_coarse_to_fine[k-1]
      for (lv in unique(d[[child]])) {
        rows <- d[[parent]][d[[child]] == lv]
        par_lv <- names(sort(table(rows), decreasing = TRUE))[1]
        cc <- coords[coords$scale == child  & coords$level == lv, ]
        pc <- coords[coords$scale == parent & coords$level == par_lv, ]
        if (nrow(cc) && nrow(pc))
          links[[length(links)+1]] <- data.frame(
            x = cc$MDS1, y = cc$MDS2, xend = pc$MDS1, yend = pc$MDS2,
            child_scale = child)
      }
    }
  }
  links_df <- if (length(links)) do.call(rbind, links) else NULL

  list(coords = coords, links = links_df, stress = mds$stress,
       scales = order_coarse_to_fine, n = nrow(d))
}
nested <- build_nested()
if (!is.null(nested)) {
  saveRDS(nested, file.path(CACHE, "nested_ordination.rds"))
  say(sprintf("  cached nested_ordination.rds (stress %.3f, %d centroids)",
              nested$stress, nrow(nested$coords)))
}


# =============================================================================
# SECTION 5  RESIDUAL ordination (Anderson §5.2)
#   Remove the environment effect from G via the residual-maker (I-H)G(I-H),
#   then ordinate region centroids in the residualized space. Shows whether
#   pure-region structure remains once environment is partialled out.
# =============================================================================
say("\n########## RESIDUAL ORDINATION (region | environment) ##########")
build_residual <- function(GROUP) {
  d <- rec
  for (t in LOG) d[[t]] <- log(d[[t]])
  d[[GROUP]][trimws(tolower(d[[GROUP]])) %in% c("other","","na","unknown")] <- NA
  cc <- complete.cases(d[, TRAITS], d[, PREDICTORS]) & !is.na(d[[GROUP]]); d <- d[cc, ]
  keep <- names(which(table(d[[GROUP]]) >= MIN_N)); d <- d[d[[GROUP]] %in% keep, ]
  d[[GROUP]] <- factor(d[[GROUP]])
  if (nrow(d) < 30 || nlevels(d[[GROUP]]) < 3) { say(sprintf("  %s too few; skipped", GROUP)); return(NULL) }
  kept <- prune_predictors(d, PREDICTORS)

  D <- dist(scale(d[, TRAITS]), method = "euclidean")
  A <- -0.5 * as.matrix(D)^2; n <- nrow(A)
  J <- diag(n) - matrix(1/n, n, n); G <- J %*% A %*% J

  # environment design + hat matrix H; residual-maker (I-H)
  X <- model.matrix(as.formula(paste("~", paste(kept, collapse = " + "))), data = d)
  H <- X %*% solve(crossprod(X), t(X))
  I <- diag(n)
  GR <- (I - H) %*% G %*% (I - H)         # environment-residualized Gower

  # centroids of GROUP in residualized space
  g <- d[[GROUP]]; lev <- levels(g)
  Gc <- matrix(NA_real_, length(lev), length(lev), dimnames = list(lev, lev))
  for (a in seq_along(lev)) for (b in seq_along(lev)) {
    ia <- which(g == lev[a]); ib <- which(g == lev[b]); Gc[a,b] <- mean(GR[ia, ib])
  }
  DCr <- matrix(0, length(lev), length(lev), dimnames = list(lev, lev))
  for (a in seq_along(lev)) for (b in seq_along(lev))
    DCr[a,b] <- sqrt(max(Gc[a,a] - 2*Gc[a,b] + Gc[b,b], 0))
  mds <- tryCatch(vegan::metaMDS(as.dist(DCr), k = 2, trace = 0, autotransform = FALSE),
                  error = function(e) NULL)
  if (is.null(mds)) { say(sprintf("  %s mMDS failed; skipped", GROUP)); return(NULL) }
  xy <- as.data.frame(vegan::scores(mds, display = "sites")); names(xy) <- c("MDS1","MDS2")
  xy$grp <- rownames(xy); xy$n <- as.integer(table(g)[xy$grp])
  saveRDS(list(GROUP = GROUP, coords = xy, stress = mds$stress, kept = kept),
          file.path(CACHE, sprintf("residual_%s.rds", tolower(GROUP))))
  say(sprintf("  cached residual_%s.rds (stress %.3f)", tolower(GROUP), mds$stress))
}
invisible(lapply(REGIONS, build_residual))


# =============================================================================
# SECTION 6  R2-by-scale trend (assembled from region + varpart summaries)
# =============================================================================
r2_trend <- region_summ %>%
  transmute(scale = grouping, groups, permanova_R2 = R2, adonis_p) %>%
  left_join(varpart_summ %>% transmute(scale, pure_region, frac_region_absorbed),
            by = "scale")
ord_lv <- intersect(c("REALM","PROVINC","ECOREGION"), r2_trend$scale)
r2_trend$scale <- factor(r2_trend$scale, levels = ord_lv)
r2_trend <- r2_trend[order(r2_trend$scale), ]
saveRDS(r2_trend, file.path(CACHE, "r2_by_scale.rds"))
say("\nR2-by-scale trend cached:"); print(r2_trend, row.names = FALSE)


# =============================================================================
# SECTION 7  DISTANCE-DECAY / MANTEL
#   Tests "nearer regions -> more likely NON-significantly different". Per scale:
#   pairwise great-circle geo distance vs trait D[C]; simple + partial Mantel
#   (control env), and logistic P(significant) ~ geo distance (perm p).
# =============================================================================
say("\n########## DISTANCE-DECAY / MANTEL ##########")
have_geosphere <- requireNamespace("geosphere", quietly = TRUE)
gc_km <- function(lon1, lat1, lon2, lat2) {
  if (have_geosphere) geosphere::distHaversine(c(lon1,lat1), c(lon2,lat2))/1000
  else { R <- 6371; toR <- pi/180; dlat <- (lat2-lat1)*toR; dlon <- (lon2-lon1)*toR
    a <- sin(dlat/2)^2 + cos(lat1*toR)*cos(lat2*toR)*sin(dlon/2)^2
    R*2*atan2(sqrt(a), sqrt(1-a)) }
}
env_centroid_dist <- function(GROUP) {
  d <- rec; for (t in LOG) d[[t]] <- log(d[[t]])
  d[[GROUP]][trimws(tolower(d[[GROUP]])) %in% c("other","","na","unknown")] <- NA
  cc <- complete.cases(d[, PREDICTORS]) & !is.na(d[[GROUP]]); d <- d[cc, ]
  keep <- names(which(table(d[[GROUP]]) >= MIN_N)); d <- d[d[[GROUP]] %in% keep, ]
  if (!nrow(d)) return(NULL)
  ag <- aggregate(scale(d[, PREDICTORS]), list(grp = d[[GROUP]]), mean)
  rownames(ag) <- ag$grp; ag$grp <- NULL
  as.matrix(dist(ag, method = "euclidean"))
}
build_distdecay <- function(GROUP) {
  tag <- tolower(GROUP)
  rf <- file.path(CACHE, sprintf("region_%s.rds", tag))
  if (!file.exists(rf)) return(invisible())
  R <- readRDS(rf); if (is.null(R$geo)) return(invisible())
  geo <- R$geo; geo$.grp <- as.character(geo$.grp)
  DCt <- as.matrix(R$DC); pw <- R$pairwise
  # FULL set (trait-complete, matches heatmap) for geo/trait/sig/simple-Mantel/logistic
  lev <- intersect(rownames(DCt), geo$.grp); if (length(lev) < 4) return(invisible())
  DCt <- DCt[lev, lev]
  gm <- matrix(0, length(lev), length(lev), dimnames = list(lev, lev))
  for (i in seq_along(lev)) for (j in seq_along(lev)) if (i < j) {
    gi <- geo[geo$.grp == lev[i], ]; gj <- geo[geo$.grp == lev[j], ]
    gm[i,j] <- gm[j,i] <- gc_km(gi$lon, gi$lat, gj$lon, gj$lat)
  }
  ut <- which(upper.tri(gm), arr.ind = TRUE)
  pairs <- data.frame(r1 = lev[ut[,1]], r2 = lev[ut[,2]],
                      geo_km = gm[ut], traitDC = DCt[ut])
  k <- function(a,b) paste(pmin(a,b), pmax(a,b))
  pw$key <- k(pw$r1, pw$r2); pairs$key <- k(pairs$r1, pairs$r2)
  pairs$sig <- pw$sig[match(pairs$key, pw$key)]
  # ENV subset (env-complete) for partial Mantel only
  envM_full <- env_centroid_dist(GROUP)
  env_lev <- if (is.null(envM_full)) character(0) else intersect(lev, rownames(envM_full))
  dropped <- setdiff(lev, env_lev)
  if (length(dropped))
    say(sprintf("  [%s] partial-Mantel drops: %s", tag, paste(dropped, collapse = ", ")))
  em <- if (is.null(envM_full)) NULL else as.matrix(envM_full)
  pairs$envDC <- if (is.null(em)) NA_real_ else
    mapply(function(a,b) if (a %in% env_lev && b %in% env_lev) em[a,b] else NA_real_,
           pairs$r1, pairs$r2)
  set.seed(1); mant  <- vegan::mantel(as.dist(DCt), as.dist(gm), permutations = 9999)
  pmant <- NULL
  if (length(env_lev) >= 4) {
    set.seed(1); pmant <- tryCatch(
      vegan::mantel.partial(as.dist(DCt[env_lev, env_lev]), as.dist(gm[env_lev, env_lev]),
                            as.dist(em[env_lev, env_lev]), permutations = 9999),
      error = function(e) NULL)
  }
  glm_fit <- tryCatch(glm(as.integer(sig %in% TRUE) ~ geo_km, data = pairs, family = binomial()),
                      error = function(e) NULL)
  perm_p <- NA_real_
  if (!is.null(glm_fit)) {
    obs <- coef(glm_fit)[["geo_km"]]; set.seed(1)
    nul <- replicate(4999, { pp <- pairs; pp$sig <- sample(pp$sig)
      f <- tryCatch(glm(as.integer(sig %in% TRUE) ~ geo_km, data = pp, family = binomial()),
                    error = function(e) NULL)
      if (is.null(f)) NA_real_ else coef(f)[["geo_km"]] })
    perm_p <- mean(abs(nul) >= abs(obs), na.rm = TRUE)
  }
  saveRDS(list(GROUP = GROUP, pairs = pairs, mantel = mant, pmantel = pmant,
               glm = glm_fit, glm_perm_p = perm_p,
               n_regions_full = length(lev), n_regions_env = length(env_lev),
               dropped_env = dropped),
          file.path(CACHE, sprintf("distdecay_%s.rds", tag)))
  say(sprintf("  cached distdecay_%s.rds (full %d realms/%d pairs; env %d; Mantel r=%.2f, partial r=%s)",
              tag, length(lev), nrow(pairs), length(env_lev), mant$statistic,
              if (!is.null(pmant)) sprintf("%.2f", pmant$statistic) else "NA"))
}
invisible(lapply(REGIONS, build_distdecay))

# also stash run-level metadata the figure script needs
saveRDS(list(REGIONS = REGIONS, TRAITS = TRAITS, PREDICTORS = PREDICTORS,
             LAT = LAT, LON = LON, NPERM = NPERM, MIN_N = MIN_N),
        file.path(CACHE, "meta.rds"))

say("\nANALYSIS DONE. Cache written to:", CACHE)
say("Next: run permanova_2_figures.R (fast; no permutations).")

