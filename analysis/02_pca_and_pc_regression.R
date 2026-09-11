# =============================================================================
# Regional and environmental analysis
#
# This script:
#   - identifies and corrects predictor and trait skew (log transforms)
#   - builds the PCA on the trait subset and interprets the axes
#   - tests for latitudinal gradients in the first two PCs
#       (OLS slope + R2, Spearman rho; abs and signed latitude)
#   - tests for differences in the distribution of the first two PCs across
#       MEOW regions (assumption-gated: Welch's ANOVA + Games-Howell, or
#       Kruskal-Wallis + Dunn with BH adjustment; compact letter display)
#   - tests for differences in the distribution of the first two PCs across
#       climate zones (tropical vs temperate by region majority; assumption-
#       gated Welch's t-test or Mann-Whitney U, with effect size)
#   - tests for predictive power of environmental variables on the PCs
#       (single-variable regression, then OLS multiple regression, upgraded
#       to GLS with an exponential spatial correlation structure where Moran's
#       I flags residual spatial autocorrelation)
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(car); library(ape); library(nlme)
})

# =============================================================================
# SETUP: config, paths, columns, palettes
# =============================================================================

## --- locate and source config.R (works with the .Rproj open, via source()) ----
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
IN  <- AGG_CSV
OUT <- OUT_DIR

## --- trait and predictor definitions ----
TRAITS <- c("length_frond_mm","max_width_mm","stipe_to_length_ratio",
            "solidity","peak_position_frac","branched_fraction","circularity")
LOG    <- c("length_frond_mm","max_width_mm", "stipe_to_length_ratio")

# env predictors, priority order (pruning honors this)
# NOTE: abs_latitude REMOVED from the multiple-regression predictor set --
# it is a coordinate proxy that competed with SST for temperature variance
# (dropped SST's VIF from ~10 to ~5). Latitude is handled separately in the
# dedicated latitude-gradient section as a descriptive framework.
PREDICTORS <- c("SST_clim","nitrate_clim","phosphate_clim",
                "salinity_clim","wave.height", "kd490_clim", "PAR_mean", "iron_clim", "spco2_clim")

## --- per-predictor LOG flags (applied upstream; regression + figures share them) ----
LOG_PREDICTORS <- c("nitrate_clim","phosphate_clim","iron_clim",
                    "kd490_clim","wave.height")
LOG_OFFSET_FRAC <- 0.01   # for vars with zeros: add 1% of the positive min before log

## --- coordinate + realm columns ----
LATCOL <- "Latitude"; LONCOL <- "Longitude"
REALMCOL <- "realm"   # realm column for point coloring (auto-detected below)

## --- load data + resolve which columns are actually present ----
rec <- read.csv(IN, stringsAsFactors = FALSE, check.names = FALSE)
TRAITS     <- intersect(TRAITS, names(rec))
PREDICTORS <- intersect(PREDICTORS, names(rec))
stopifnot(all(c(LATCOL, LONCOL) %in% names(rec)))

realm_candidates <- names(rec)[grepl("realm", names(rec), ignore.case = TRUE)]
if (length(realm_candidates)) REALMCOL <- realm_candidates[1] else REALMCOL <- NA

## --- complete cases on TRAITS ONLY: this is the full PCA population ----
#  The PCA is computed on every trait-complete specimen. Specimens missing
#  environmental data are NOT dropped here -- they are filtered out later as a
#  SUBSAMPLE, only where env analysis is explicitly needed. PC scores are
#  computed once on the full population and never recomputed.
d <- rec
for (t in LOG) d[[t]] <- log(d[[t]])
cc <- complete.cases(d[, TRAITS])
d  <- d[cc, ]
cat(sprintf("Trait-complete records (full PCA population): %d\n", nrow(d)))

## --- realm color palette (tropical -> reds, temperate -> blues) ----
realm_cols <- c(
  # Tropical realms -> shades of red
  "Atlantic Warm Water"          = "#67000d",
  "Indo-Pacific Warm Water"      = "#a50f15",
  "Tropical Atlantic"            = "#cb181d",
  "Tropical Eastern Pacific"     = "#ef3b2c",
  "Western Indo-Pacific"         = "#fb6a4a",
  "Central Indo-Pacific"         = "#fc9272",
  "Eastern Indo-Pacific"         = "#fcbba1",
  # Temperate realms -> shades of blue
  "Temperate Northern Atlantic"  = "#08306b",
  "Temperate Northern Pacific"   = "#08519c",
  "Temperate South America"      = "#2171b5",
  "Temperate Australasia"        = "#4292c6"
)

## --- clean display names for env predictors (used by coef + bubble plots) ----
# keyed by BOTH raw and log_ names so lookup works before/after log transform.
DISPLAY_NAME <- c(
  "SST_clim"             = "Sea surface temperature",
  "nitrate_clim"         = "Nitrate",
  "log_nitrate_clim"     = "Nitrate (log)",
  "phosphate_clim"       = "Phosphate",
  "log_phosphate_clim"   = "Phosphate (log)",
  "salinity_clim"        = "Salinity",
  "wave.height"          = "Wave height",
  "log_wave.height"      = "Wave height (log)",
  "kd490_clim"           = "Turbidity (Kd490)",
  "log_kd490_clim"       = "Turbidity (Kd490, log)",
  "PAR_mean"             = "Photosynthetically available radiation",
  "iron_clim"            = "Iron",
  "log_iron_clim"        = "Iron (log)",
  "abs_latitude"         = "Absolute latitude",
  "spco2_clim"           = "Surface ocean pCO\u2082",
  "log_spco2_clim"       = "Surface ocean pCO\u2082 (log)"
)

## --- clean display names for traits (used in biplot) ----
TRAIT_DISPLAY <- c(
  "length_frond_mm"        = "Frond length",
  "max_width_mm"           = "Max width",
  "stipe_to_length_ratio"  = "Stipe width:Length",
  "solidity"               = "Solidity",
  "peak_position_frac"     = "Peak position",
  "branched_fraction"      = "Branched fraction",
  "circularity"            = "Circularity"
)

## --- apply predictor log transforms UPSTREAM (shared by regression + figures) ----
LOG_PREDICTORS <- intersect(LOG_PREDICTORS, PREDICTORS)
pretty_name <- setNames(PREDICTORS, PREDICTORS)   # raw -> label, updated below

if (length(LOG_PREDICTORS)) {
  cat("\nLog-transforming predictors:", paste(LOG_PREDICTORS, collapse=", "), "\n")
  for (p in LOG_PREDICTORS) {
    x <- d[[p]]
    pos <- x[is.finite(x) & x > 0]
    off <- if (any(x <= 0, na.rm = TRUE) && length(pos))
      LOG_OFFSET_FRAC * min(pos) else 0
    if (off > 0) cat(sprintf("  %s has <=0 values; log(x + %.3g)\n", p, off))
    newname <- paste0("log_", p)
    d[[newname]] <- log(x + off)
    PREDICTORS[PREDICTORS == p] <- newname
    pretty_name[[newname]] <- paste0("log ", p)
  }
}

## --- helper: raw-or-log term -> clean label (falls back to pretty_name, then raw) ----
clean_lab <- function(x) {
  out <- ifelse(x %in% names(DISPLAY_NAME), DISPLAY_NAME[x],
                ifelse(x %in% names(pretty_name), pretty_name[x], x))
  unname(out)
}

## --- significance-star helper (shared across all tests below) ----
star <- function(p) ifelse(is.na(p), "",
                           ifelse(p < .001, "***",
                                  ifelse(p < .01,  "**",
                                         ifelse(p < .05,  "*",
                                                ifelse(p < .1,   ".", "")))))

# =============================================================================
# PCA: build morphospace on the trait subset
# =============================================================================

## --- fit the PCA (single fit; PC scores reused everywhere) ----
Xtr <- scale(d[, TRAITS])
pca <- prcomp(Xtr, center = FALSE, scale. = FALSE)
ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
d$PC1 <- pca$x[,1]; d$PC2 <- pca$x[,2]

cat(sprintf("\nPC1 %.1f%%  PC2 %.1f%%  (cumulative %.1f%%)\n", ve[1], ve[2], ve[1]+ve[2]))
cat("\nLoadings PC1/PC2 (interpret what each axis means):\n")
load <- round(pca$rotation[, 1:2], 3)
print(load)
write.csv(load, file.path(OUT, "pc_loadings.csv"))

## --- full loadings for EVERY PC axis ----
load_all <- round(pca$rotation, 3)
colnames(load_all) <- paste0("PC", seq_len(ncol(load_all)))
cat("\nFull loadings (all PC axes):\n")
print(load_all)
write.csv(load_all, file.path(OUT, "pc_loadings_all.csv"))

## --- scree / cumulative-variance table (all PCs) ----
# ve is % variance per PC (defined above). Save the full table now; the
# figure is drawn later (save_scree) once ggplot is loaded.
scree_df <- data.frame(
  PC = factor(paste0("PC", seq_along(ve)), levels = paste0("PC", seq_along(ve))),
  variance_pct = round(ve, 3),
  cumulative_pct = round(cumsum(ve), 3)
)
write.csv(scree_df, file.path(OUT, "pc_scree_variance.csv"), row.names = FALSE)
cat("\nVariance explained per PC (with cumulative):\n")
print(scree_df, row.names = FALSE)

# =============================================================================
# SUBSAMPLING: env-complete subset + plotting pools
# =============================================================================

## --- env subsample: specimens with complete predictors + coords ----
#   The PCA above was fit on the FULL trait-complete population. This is a
#   SUBSET of that population, used ONLY by the environment-dependent analyses
#   (collinearity pruning, OLS/GLS regression, scatter grids, coefficient plot,
#   bubble plot). PC1/PC2 are carried over UNCHANGED -- the PCA is never
#   recomputed, so PC scores mean the same thing in every figure. Trait-only
#   figures (biplot, region boxplots) keep the full population `d`.
env_cc <- complete.cases(d[, PREDICTORS], d[, c(LATCOL, LONCOL)])
d_env  <- d[env_cc, ]
cat(sprintf("Env-complete subsample: %d of %d trait-complete specimens\n",
            nrow(d_env), nrow(d)))

## --- plotting pools (drop NA-realm rows for FIGURES ONLY) ----
#   dp      = ENV pool (d_env) -> env figures: scatter grids, coef, bubble
#   dp_full = FULL trait pool (d) -> trait figures: biplot, region boxplots
#   Both share the same NA-realm drop; PC scores are identical (one PCA).
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

drop_na_realm <- function(x, tag) {
  if (!is.na(REALMCOL) && REALMCOL %in% names(x)) {
    realm_raw <- x[[REALMCOL]]
    keep_plot <- !is.na(realm_raw) & !(trimws(as.character(realm_raw)) %in% c("NA",""))
    n_drop <- sum(!keep_plot)
    if (n_drop > 0)
      cat(sprintf("[plots:%s] dropping %d NA-realm specimen(s) from figures (kept in all stats)\n",
                  tag, n_drop))
    out <- x[keep_plot, ]
    out$.realm <- droplevels(as.factor(out[[REALMCOL]]))
  } else {
    cat(sprintf("[plots:%s] no realm column found; plotting all points in one color.\n", tag))
    out <- x
    out$.realm <- factor("all")
  }
  out
}

dp      <- drop_na_realm(d_env, "env")    # env figures
dp_full <- drop_na_realm(d,     "full")   # trait figures (biplot, region boxplots)

realm_levels <- union(levels(dp$.realm), levels(dp_full$.realm))
missing_lv <- setdiff(realm_levels, names(realm_cols))
if (length(missing_lv)) {
  cat("[plots] realms not in reference palette (grey):", paste(missing_lv, collapse=", "), "\n")
  realm_cols[missing_lv] <- "#BBBBBB"
}

# =============================================================================
# PCA FIGURES: biplot + scree
# =============================================================================

## --- PCA biplot (pure morphospace -> full trait-complete population) ----
save_biplot <- function() {
  scores <- data.frame(PC1 = dp_full$PC1, PC2 = dp_full$PC2, realm = dp_full$.realm)
  ld <- as.data.frame(pca$rotation[, 1:2]); names(ld) <- c("PC1", "PC2")
  ld$trait <- ifelse(rownames(ld) %in% names(TRAIT_DISPLAY),
                     TRAIT_DISPLAY[rownames(ld)], rownames(ld))
  score_rad <- min(max(abs(scores$PC1)), max(abs(scores$PC2)))
  load_rad  <- max(sqrt(ld$PC1^2 + ld$PC2^2))
  arrow_scale <- 0.8 * score_rad / load_rad
  ld$xend <- ld$PC1 * arrow_scale; ld$yend <- ld$PC2 * arrow_scale
  ld$lx <- ld$xend * 1.08; ld$ly <- ld$yend * 1.08
  xlab <- sprintf("PC1 (%.1f%%)", ve[1]); ylab <- sprintf("PC2 (%.1f%%)", ve[2])
  bp <- ggplot(scores, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_point(aes(color = realm), size = 1.6, alpha = 0.8) +
    scale_color_manual(values = realm_cols, drop = FALSE, name = "REALM") +
    geom_segment(data = ld, inherit.aes = FALSE, aes(x = 0, y = 0, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.18, "cm")), color = "black", linewidth = 0.5) +
    geom_text(data = ld, inherit.aes = FALSE, aes(x = lx, y = ly, label = trait),
              size = 3, fontface = "bold", color = "black") +
    labs(x = xlab, y = ylab, title = "PCA biplot",
         subtitle = "Shades of blue = temperate realm; Shades of red = tropical realm") +
    coord_equal() + theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 11), legend.key.size = unit(0.4, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  ggsave(file.path(OUT, "pca_biplot_linear.png"), bp, width = 8, height = 6, dpi = 150)
  cat("[biplot] wrote pca_biplot_linear.png\n")
  lg <- ggplot(dp_full, aes(x = PC1, y = PC2, color = .realm)) + geom_point() +
    scale_color_manual(values = realm_cols, drop = FALSE, name = "REALM") +
    theme_bw(base_size = 10) + guides(color = guide_legend(override.aes = list(size = 3)))
  ggsave(file.path(OUT, "scatter_realm_legend.png"), lg, width = 5, height = 5, dpi = 150)
}

## --- scree plot + cumulative variance ----
#   Bars = variance explained per PC; line+points = running cumulative %.
#   Reference lines at 70% and 90% mark common retention thresholds.
save_scree <- function() {
  sd <- scree_df
  sd$cum <- sd$cumulative_pct
  ymax <- max(100, max(sd$cum))
  g <- ggplot(sd, aes(x = PC)) +
    geom_col(aes(y = variance_pct), fill = "#4292c6", color = "grey30", width = 0.7) +
    geom_line(aes(y = cum, group = 1), color = "#c0392b", linewidth = 0.8) +
    geom_point(aes(y = cum), color = "#c0392b", size = 2) +
    geom_text(aes(y = variance_pct, label = sprintf("%.1f%%", variance_pct)),
              vjust = -0.4, size = 2.8, color = "grey20") +
    geom_text(aes(y = cum, label = sprintf("%.0f%%", cum)),
              vjust = -0.8, size = 2.8, color = "#c0392b") +
    geom_hline(yintercept = c(70, 90), linetype = "dashed", color = "grey60", linewidth = 0.3) +
    annotate("text", x = length(ve), y = 71, label = "70%", hjust = 1, size = 2.6, color = "grey50") +
    annotate("text", x = length(ve), y = 91, label = "90%", hjust = 1, size = 2.6, color = "grey50") +
    scale_y_continuous(limits = c(0, ymax + 5), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Variance explained (%)",
         title = "PCA scree + cumulative variance",
         subtitle = "Bars = per-PC variance; red line = cumulative") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
  ggsave(file.path(OUT, "pca_scree.png"), g, width = 7, height = 5, dpi = 150)
  cat("[scree] wrote pca_scree.png + pc_scree_variance.csv\n")
}

save_biplot()
save_scree()
cat("[biplot] biplot + realm legend written.\n")

# =============================================================================
# LATITUDINAL GRADIENT: PC1/PC2 vs latitude
# =============================================================================
#   Latitude is a coordinate proxy, not a mechanism, so it is analyzed on its
#   own here rather than competing with SST inside the multiple model. For BOTH
#   PC1 and PC2 we fit against:
#     - abs_latitude  (distance from equator; the classic diversity/size axis)
#     - Latitude      (signed; exposes hemisphere asymmetry if present)
#   Reports OLS slope, R2, and Spearman rho; draws a 4-panel figure (2 PCs x 2
#   latitude flavors) with realm-colored points and a fit line drawn only where
#   the slope is significant. Runs on the env subsample `dp` (coordinates
#   required); PC scores are the same single-PCA scores.

save_latitude_gradient <- function() {
  ## --- build abs_latitude on the plotting frame if not already present ----
  if (!"abs_latitude" %in% names(dp) && LATCOL %in% names(dp))
    dp$abs_latitude <- abs(dp[[LATCOL]])
  
  LAT_VARS <- intersect(c("abs_latitude", LATCOL), names(dp))
  if (!length(LAT_VARS)) { cat("[lat] no latitude column available; skipping.\n"); return(invisible()) }
  
  LAT_LABEL <- c(abs_latitude = "Absolute latitude (\u00b0)",
                 Latitude      = "Latitude (\u00b0, signed)")
  
  ## --- stats table across PC x latitude-flavor ----
  rows <- list(); panels <- list()
  for (PC in c("PC1","PC2")) {
    for (LV in LAT_VARS) {
      dd <- data.frame(y = dp[[PC]], x = dp[[LV]], realm = dp$.realm)
      dd <- dd[is.finite(dd$x) & is.finite(dd$y), ]
      if (nrow(dd) < 3) next
      m  <- lm(y ~ x, data = dd); sm <- summary(m)
      slope   <- coef(m)[2]
      slope_p <- if (nrow(sm$coefficients) >= 2) sm$coefficients[2,4] else NA
      r2      <- sm$r.squared
      sp      <- suppressWarnings(cor.test(dd$x, dd$y, method = "spearman"))
      gate    <- is.finite(slope_p) && slope_p < .05
      
      rows[[length(rows)+1]] <- data.frame(
        PC = PC, latitude = LV,
        slope = round(slope,4), slope_p = signif(slope_p,3),
        R2 = round(r2,4),
        spearman_rho = round(unname(sp$estimate),3),
        spearman_p = signif(sp$p.value,3),
        n = nrow(dd), sig = star(slope_p),
        row.names = NULL, stringsAsFactors = FALSE)
      
      lab <- sprintf("R2 = %.3f\nrho = %.3f\nslope p = %s\nn = %d",
                     r2, unname(sp$estimate),
                     if (is.na(slope_p)) "NA" else format.pval(slope_p, digits = 2, eps = 1e-3),
                     nrow(dd))
      xlab <- if (LV %in% names(LAT_LABEL)) LAT_LABEL[[LV]] else LV
      g <- ggplot(dd, aes(x = x, y = y)) +
        geom_point(aes(color = realm), size = 1.5, alpha = 0.8, show.legend = FALSE)
      if (gate)
        g <- g + geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                             color = "black", linewidth = 0.7)
      g <- g +
        scale_color_manual(values = realm_cols, drop = FALSE) +
        annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05,
                 size = 2.7, label = lab) +
        labs(x = xlab, y = PC, title = sprintf("%s vs %s", PC, xlab)) +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(size = 9, hjust = 0.5))
      panels[[length(panels)+1]] <- g
    }
  }
  
  lat_tab <- do.call(rbind, rows)
  write.csv(lat_tab, file.path(OUT, "latitude_gradient_stats.csv"), row.names = FALSE)
  cat("\n=== Latitudinal gradient (standalone) ===\n")
  print(lat_tab[, c("PC","latitude","slope","slope_p","R2","spearman_rho","sig")],
        row.names = FALSE)
  
  ## --- assemble panel figure ----
  ncol <- length(LAT_VARS)   # abs_lat + signed side by side; PCs stack in rows
  ttl  <- "Latitudinal gradients in morphological axes (standalone; not in mult. regression)"
  built <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    built <- Reduce(`+`, panels) + plot_layout(ncol = ncol) +
      plot_annotation(title = ttl) & theme(plot.title = element_text(size = 12))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    built <- gridExtra::arrangeGrob(grobs = panels, ncol = ncol, top = ttl)
  } else {
    for (i in seq_along(panels))
      ggsave(file.path(OUT, sprintf("latitude_gradient_panel_%d.png", i)),
             panels[[i]], width = 5, height = 4, dpi = 150)
    cat("[lat] install patchwork/gridExtra for a combined panel.\n"); return(invisible())
  }
  w <- 4.4*ncol; h <- 3.8*2
  ggsave(file.path(OUT, "latitude_gradient.png"), built, width = w, height = h, dpi = 150, limitsize = FALSE)
  cat("[lat] wrote latitude_gradient.png + latitude_gradient_stats.csv\n")
}
save_latitude_gradient()

# =============================================================================
# SPATIAL ANALYSIS - MEOW REGIONS
# =============================================================================
#   Tests for differences in the distribution of PC1 and PC2 across MEOW
#   region scales (realm / province / ecoregion). Each test is assumption-
#   gated: Welch's ANOVA + Games-Howell (parametric) or Kruskal-Wallis + Dunn
#   with BH adjustment (non-parametric), with a compact letter display (CLD)
#   via multcompView. Runs on the FULL trait-complete population (dp_full),
#   NOT the env subsample -- PC scores come from the single PCA (never
#   recomputed), so no second PCA and no sign-alignment are needed. This
#   recovers ecoregions that were only dropped from env figures for missing
#   environmental data, matching the PERMANOVA pool.

## --- packages for post-hoc + CLD ----
suppressPackageStartupMessages({
  library(ggplot2)
  have_mc  <- requireNamespace("multcompView", quietly = TRUE)
  have_rst <- requireNamespace("rstatix",      quietly = TRUE)
})

## --- switch region pool to the full trait-complete population ----
dp <- dp_full
cat(sprintf("region morphology pool (full trait-complete): %d specimens\n",
            nrow(dp)))

## --- detect MEOW scale columns present in the data (realm/province/ecoregion) ----
# NOTE: the data contains TWO province-like columns -- the correct MEOW
# "PROVINCE" (from the ecoregion polygon join) and a stale free-text
# "State.Province" (Darwin Core locality: states/islands/parishes, ~200
# blanks). We must pick the MEOW one. Strategy: for province, prefer an EXACT
# "province"/"provinc" match and explicitly exclude "state.province".
MEOW_PATTERNS <- c(realm = "realm", province = "provinc", ecoregion = "ecoregion")
meow_cols <- sapply(names(MEOW_PATTERNS), function(nm) {
  pat <- MEOW_PATTERNS[[nm]]
  hit <- names(dp)[grepl(pat, names(dp), ignore.case = TRUE)]
  if (nm == "province") {
    hit <- hit[!grepl("state.?provinc", hit, ignore.case = TRUE)]     # drop free-text State.Province
    exact <- hit[grepl("^provinc(e)?$", hit, ignore.case = TRUE)]     # prefer exact MEOW PROVINC/PROVINCE
    if (length(exact)) return(exact[1])
  }
  if (length(hit)) hit[1] else NA_character_
})
names(meow_cols) <- names(MEOW_PATTERNS)
meow_cols <- meow_cols[!is.na(meow_cols)]
cat("[MEOW] scales found:",
    paste(sprintf("%s=%s", names(meow_cols), meow_cols), collapse = ", "), "\n")

## --- clean display names for each MEOW scale (used in titles) ----
SCALE_DISPLAY <- c(realm = "MEOW Realm", province = "MEOW Province", ecoregion = "MEOW Ecoregion")
scale_disp <- function(s) if (s %in% names(SCALE_DISPLAY)) unname(SCALE_DISPLAY[s]) else
  paste0(toupper(substring(s,1,1)), substring(s,2))

## --- accumulator for the publication-stats block (filled by region_test_plot) ----
REGION_STATS <- list()

## --- realm-inherited shading helper (provinces/ecoregions inherit parent hue) ----
#  Each province/ecoregion inherits its parent realm's thermal class from the
#  data, then gets a graded shade within that class (ordered by the region's
#  median latitude so the ramp reads geographically).

# tropical vs temperate realm sets (matches the realm_cols split upstream)
TROPICAL_REALMS <- c("Atlantic Warm Water","Indo-Pacific Warm Water",
                     "Tropical Atlantic","Tropical Eastern Pacific",
                     "Western Indo-Pacific","Central Indo-Pacific",
                     "Eastern Indo-Pacific")
TEMPERATE_REALMS <- c("Temperate Northern Atlantic","Temperate Northern Pacific",
                      "Temperate South America","Temperate Australasia",
                      "Southern Ocean")

# sequential ramps (dark -> light), same hues as realm_cols
RED_RAMP  <- c("#67000d","#a50f15","#cb181d","#ef3b2c","#fb6a4a","#fc9272","#fcbba1")
BLUE_RAMP <- c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1")

REALMCOL_MEOW <- if ("REALM" %in% names(dp)) "REALM" else
{ hit <- names(dp)[grepl("realm", names(dp), ignore.case = TRUE)]
if (length(hit)) hit[1] else NA_character_ }

# build a named color vector for the regions present in `dat`
build_scale_cols <- function(dat, scale_col) {
  regs <- levels(dat$region)
  if (is.na(REALMCOL_MEOW)) {  # no parent realm available -> grey
    return(setNames(rep("#BBBBBB", length(regs)), regs))
  }
  # parent realm = the modal realm among specimens of each region
  parent <- sapply(regs, function(r) {
    rl <- dp[[REALMCOL_MEOW]][ as.character(dp[[scale_col]]) == r ]
    rl <- rl[!is.na(rl) & !(trimws(as.character(rl)) %in% c("NA",""))]
    if (!length(rl)) return(NA_character_)
    names(sort(table(rl), decreasing = TRUE))[1]
  })
  # region median latitude, for ordering shades within a thermal class
  latc <- if ("Latitude" %in% names(dp)) "Latitude" else
  { h <- names(dp)[grepl("lat", names(dp), ignore.case = TRUE)]; if (length(h)) h[1] else NA }
  reg_lat <- sapply(regs, function(r) {
    if (is.na(latc)) return(NA_real_)
    median(abs(dp[[latc]][ as.character(dp[[scale_col]]) == r ]), na.rm = TRUE)
  })
  
  cls <- ifelse(parent %in% TROPICAL_REALMS, "trop",
                ifelse(parent %in% TEMPERATE_REALMS, "temp", "other"))
  
  out <- setNames(rep("#BBBBBB", length(regs)), regs)
  assign_ramp <- function(members, ramp, ascending = TRUE) {
    if (!length(members)) return(invisible())
    o <- order(reg_lat[members], decreasing = !ascending, na.last = TRUE)
    members <- members[o]
    idx <- if (length(members) == 1) 1 else
      round(seq(1, length(ramp), length.out = length(members)))
    out[members] <<- ramp[idx]
  }
  # tropical: low latitude -> darkest red; temperate: high latitude -> darkest blue
  assign_ramp(regs[cls == "trop"], RED_RAMP,  ascending = TRUE)
  assign_ramp(regs[cls == "temp"], BLUE_RAMP, ascending = FALSE)
  out
}

## --- per-scale, per-PC region test + CLD boxplot ----
region_test_plot <- function(PC, scale_name, scale_col) {
  dat <- data.frame(y = dp[[PC]], region = as.character(dp[[scale_col]]),
                    stringsAsFactors = FALSE)
  dat <- dat[is.finite(dat$y) & !is.na(dat$region) &
               !(trimws(dat$region) %in% c("NA","")), ]
  dat$region <- droplevels(as.factor(dat$region))
  
  tab <- table(dat$region)
  keep_reg <- names(tab)[tab >= 10]
  dat <- dat[dat$region %in% keep_reg, ]
  dat$region <- droplevels(dat$region)
  ng <- nlevels(dat$region)
  if (ng < 2) {
    cat(sprintf("[%s:%s] <2 groups with n>=2; skipping.\n", scale_name, PC))
    return(invisible())
  }
  
  pmat_df <- NULL   # guard: always defined before the CLD block
  tag <- sprintf("%s | %s (%d groups)", PC, scale_name, ng)
  
  # ---- assumption checks ----
  aov_fit <- aov(y ~ region, data = dat)
  res <- residuals(aov_fit)
  sh_x <- if (length(res) > 5000) sample(res, 5000) else res
  shap_p <- tryCatch(shapiro.test(sh_x)$p.value, error = function(e) NA)
  lev_p  <- tryCatch(car::leveneTest(y ~ region, data = dat)[1, "Pr(>F)"],
                     error = function(e) NA)
  parametric <- is.finite(shap_p) && shap_p >= 0.05 &&
    is.finite(lev_p)  && lev_p  >= 0.05
  
  cat(sprintf("\n===== %s =====\n", tag))
  cat(sprintf("  Shapiro p = %s | Levene p = %s\n",
              format.pval(shap_p, digits = 3, eps = 1e-4),
              format.pval(lev_p,  digits = 3, eps = 1e-4)))
  
  # ---- omnibus + post-hoc (parametric vs non-parametric) ----
  if (parametric) {
    method_lab <- "Welch's ANOVA + Games-Howell"
    om <- oneway.test(y ~ region, data = dat, var.equal = FALSE)
    cat(sprintf("  -> parametric. Welch F = %.3f, p = %s\n", om$statistic,
                format.pval(om$p.value, digits = 3, eps = 1e-4)))
    REGION_STATS[[length(REGION_STATS)+1]] <<- data.frame(
      PC = PC, scale = scale_name, groups = ng, n = nrow(dat),
      test = "Welch ANOVA", statistic = round(unname(om$statistic),3),
      df1 = round(unname(om$parameter[1]),1), df2 = round(unname(om$parameter[2]),1),
      p = signif(om$p.value,3),
      shapiro_p = signif(shap_p,3), levene_p = signif(lev_p,3),
      row.names = NULL, stringsAsFactors = FALSE)
    if (have_rst) {
      ph <- tryCatch(rstatix::games_howell_test(dat, y ~ region),
                     error = function(e) { cat("  [posthoc error]", conditionMessage(e), "\n"); NULL })
      if (!is.null(ph)) pmat_df <- as.data.frame(ph[, c("group1","group2","p.adj")])
    } else {
      pt <- pairwise.t.test(dat$y, dat$region, p.adjust.method = "holm", pool.sd = FALSE)
      pm <- pt$p.value
      pmat_df <- do.call(rbind, lapply(rownames(pm), function(r)
        do.call(rbind, lapply(colnames(pm), function(c)
          if (!is.na(pm[r,c])) data.frame(group1=c, group2=r, p.adj=pm[r,c]) else NULL))))
    }
  } else {
    method_lab <- "Kruskal-Wallis + Dunn (BH)"
    kw <- kruskal.test(y ~ region, data = dat)
    cat(sprintf("  -> non-parametric. KW chi2 = %.3f, df = %d, p = %s\n",
                kw$statistic, kw$parameter,
                format.pval(kw$p.value, digits = 3, eps = 1e-4)))
    REGION_STATS[[length(REGION_STATS)+1]] <<- data.frame(
      PC = PC, scale = scale_name, groups = ng, n = nrow(dat),
      test = "Kruskal-Wallis", statistic = round(unname(kw$statistic),3),
      df1 = unname(kw$parameter), df2 = NA,
      p = signif(kw$p.value,3),
      shapiro_p = signif(shap_p,3), levene_p = signif(lev_p,3),
      row.names = NULL, stringsAsFactors = FALSE)
    if (have_rst) {
      ph <- tryCatch(rstatix::dunn_test(dat, y ~ region, p.adjust.method = "BH"),
                     error = function(e) { cat("  [posthoc error]", conditionMessage(e), "\n"); NULL })
      if (!is.null(ph)) pmat_df <- as.data.frame(ph[, c("group1","group2","p.adj")])
    } else {
      pw <- pairwise.wilcox.test(dat$y, dat$region, p.adjust.method = "BH")
      pm <- pw$p.value
      pmat_df <- do.call(rbind, lapply(rownames(pm), function(r)
        do.call(rbind, lapply(colnames(pm), function(c)
          if (!is.na(pm[r,c])) data.frame(group1=c, group2=r, p.adj=pm[r,c]) else NULL))))
    }
  }
  
  # ---- count significant pairwise comparisons ----
  if (!is.null(pmat_df) && nrow(pmat_df)) {
    n_sig   <- sum(pmat_df$p.adj < 0.05, na.rm = TRUE)
    n_total <- sum(!is.na(pmat_df$p.adj))
    cat(sprintf("[%s %s] significant pairs: %d of %d (%.0f%%)\n",
                scale_name, PC, n_sig, n_total,
                100 * n_sig / max(n_total, 1)))
    write.csv(
      data.frame(scale = scale_name, PC = PC,
                 n_sig = n_sig, n_total = n_total,
                 prop = round(n_sig / max(n_total, 1), 3),
                 test = method_lab),
      file.path(OUT, sprintf("region_%s_%s_pairsig.csv", scale_name, PC)),
      row.names = FALSE)
  } else {
    cat(sprintf("[%s %s] no pairwise matrix available; count skipped\n",
                scale_name, PC))
  }
  
  # ---- compact letter display (matrix input -> hyphen-safe) ----
  cld_map <- NULL
  if (have_mc && !is.null(pmat_df) && nrow(pmat_df)) {
    grps <- levels(dat$region)
    pmat <- matrix(NA_real_, length(grps), length(grps),
                   dimnames = list(grps, grps))
    for (i in seq_len(nrow(pmat_df))) {
      a <- as.character(pmat_df$group1[i]); b <- as.character(pmat_df$group2[i])
      if (a %in% grps && b %in% grps) { pmat[a,b] <- pmat_df$p.adj[i]; pmat[b,a] <- pmat_df$p.adj[i] }
    }
    diag(pmat) <- 1
    cld <- tryCatch(multcompView::multcompLetters(pmat, threshold = 0.05,
                                                  compare = "<")$Letters,
                    error = function(e) { cat("  [CLD error]", conditionMessage(e), "\n"); NULL })
    if (!is.null(cld))
      cld_map <- data.frame(region = names(cld), cld = unname(cld),
                            stringsAsFactors = FALSE)
  }
  if (is.null(cld_map)) {
    cat("  [note] CLD unavailable; no letters.\n")
    cld_map <- data.frame(region = levels(dat$region), cld = "",
                          stringsAsFactors = FALSE)
  }
  
  ytop <- tapply(dat$y, dat$region, max, na.rm = TRUE)
  cld_map$y <- ytop[cld_map$region] + 0.05 * diff(range(dat$y, na.rm = TRUE))
  cld_map$n <- as.integer(table(dat$region)[cld_map$region])
  
  med <- tapply(dat$y, dat$region, median, na.rm = TRUE)   # PC median (pairwise + reporting)
  
  # ---- x-axis order = each region's median |latitude| (tropical -> temperate) ----
  # Ordering is by geography (median abs latitude), NOT by PC value, so every
  # boxplot reads equator-to-pole left-to-right and the ramp shading lines up.
  # abslat is looked up on dp (full frame, has LATCOL) for the same regions.
  if (LATCOL %in% names(dp)) {
    reg_abslat <- tapply(abs(dp[[LATCOL]]),
                         as.character(dp[[scale_col]]), median, na.rm = TRUE)
    lev_ord <- levels(dat$region)
    lev_ord <- lev_ord[order(reg_abslat[lev_ord], na.last = TRUE)]
  } else {
    lev_ord <- names(sort(med))   # fallback: PC median if no latitude column
  }
  dat$region     <- factor(dat$region,     levels = lev_ord)
  cld_map$region <- factor(cld_map$region, levels = lev_ord)
  
  write.csv(cld_map[, c("region","cld","n")],
            file.path(OUT, sprintf("region_%s_%s_cld.csv", scale_name, PC)),
            row.names = FALSE)
  
  # ---- pairwise p-values (drive the CLD) + per-pair effect size ----
  # Rank-biserial r per pair from the group values (direction: group1 - group2).
  # Exported so specific comparisons are citable, e.g. "highest vs next (p=..)".
  rank_biserial <- function(a, b) {
    xa <- dat$y[dat$region == a]; xb <- dat$y[dat$region == b]
    if (length(xa) < 1 || length(xb) < 1) return(NA_real_)
    w <- suppressWarnings(wilcox.test(xa, xb))
    1 - (2 * unname(w$statistic)) / (length(xa) * length(xb))   # >0: a tends higher
  }
  if (!is.null(pmat_df) && nrow(pmat_df)) {
    pp <- pmat_df
    pp$median1 <- round(med[as.character(pp$group1)], 3)
    pp$median2 <- round(med[as.character(pp$group2)], 3)
    pp$rank_biserial_r <- round(mapply(rank_biserial,
                                       as.character(pp$group1), as.character(pp$group2)), 3)
    pp$p.adj <- signif(pp$p.adj, 3)
    pp$sig   <- star(pp$p.adj)
    pp <- pp[order(pp$p.adj), ]
    write.csv(pp, file.path(OUT, sprintf("region_%s_%s_pairwise.csv", scale_name, PC)),
              row.names = FALSE)
    
    # ---- focused: most-temperate group (highest |lat|) vs the next ----
    # lev_ord is now ascending in median |latitude|, so the last entry is the
    # highest-latitude (most temperate) region.
    hi   <- lev_ord[length(lev_ord)]          # highest median |latitude|
    nxt  <- lev_ord[length(lev_ord) - 1]
    hit  <- pp[(pp$group1 == hi  & pp$group2 == nxt) |
                 (pp$group1 == nxt & pp$group2 == hi), ]
    cat(sprintf("  [most temperate vs next] highest-|lat| = %s (%s median %.3f); next = %s (%s median %.3f)\n",
                hi, PC, med[hi], nxt, PC, med[nxt]))
    if (nrow(hit))
      cat(sprintf("     %s vs %s: adjusted p = %s, rank-biserial r = %.3f%s\n",
                  hi, nxt, format.pval(hit$p.adj[1], digits = 3, eps = 1e-4),
                  hit$rank_biserial_r[1],
                  if (hit$p.adj[1] < 0.05) " (significant)" else " (n.s.)"))
    else
      cat("     [top vs next comparison not in pairwise table]\n")
  }
  
  # ---- boxplot with CLD ----
  # ecoregion / province have many groups -> widen and shrink text
  wide <- ng > 15
  # realm scale: exact-name palette; finer scales: realm-inherited shades
  fill_cols <- if (scale_name == "realm") realm_cols else build_scale_cols(dat, scale_col)
  
  g <- ggplot(dat, aes(x = region, y = y)) +
    geom_boxplot(aes(fill = region), outlier.size = 0.5, width = 0.65,
                 color = "grey25", linewidth = 0.35, show.legend = FALSE) +
    geom_text(data = cld_map, aes(x = region, y = y, label = cld),
              inherit.aes = FALSE, size = if (wide) 2.6 else 3.5, fontface = "bold") +
    geom_text(data = cld_map, aes(x = region, y = -Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -0.6, size = if (wide) 1.9 else 2.5,
              color = "grey40") +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    labs(x = NULL, y = PC,
         title = sprintf("%s by %s", PC, scale_disp(scale_name)),
         subtitle = sprintf("%s; shared letter = not different (p<.05)", method_lab)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1,
                                     size = if (wide) 6 else 9),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
  
  w <- if (wide) max(10, ng * 0.42) else 8
  ggsave(file.path(OUT, sprintf("region_boxplot_%s_%s_cld.png", scale_name, PC)),
         g, width = w, height = 6, dpi = 150, limitsize = FALSE)
  cat(sprintf("[%s] wrote region_boxplot_%s_%s_cld.png + csv\n",
              scale_name, scale_name, PC))
}

for (sc in names(meow_cols))
  for (PC in c("PC1","PC2"))
    region_test_plot(PC, sc, meow_cols[[sc]])

## --- paired PC1/PC2 region figure (faceted, stacked) ----
#   Both axes shown per region: PC1 panel on top, PC2 below, sharing ONE
#   x-order set by region median |latitude|. Each PC keeps its own y-scale
#   (orthogonal axes, not comparable magnitudes) via facet free-y. CLD letters
#   are computed independently per PC and per panel:
#     PC1 -> CAPITAL letters ; PC2 -> lowercase letters
#   so a reader never cross-reads a PC1 group against a PC2 group.

# core test for ONE PC on a given scale; returns dat + cld (no plotting).
region_cld_only <- function(PC, scale_col, upper = TRUE) {
  dat <- data.frame(y = dp[[PC]], region = as.character(dp[[scale_col]]),
                    stringsAsFactors = FALSE)
  dat <- dat[is.finite(dat$y) & !is.na(dat$region) &
               !(trimws(dat$region) %in% c("NA","")), ]
  dat$region <- droplevels(as.factor(dat$region))
  keep_reg <- names(table(dat$region))[table(dat$region) >= 10]
  dat <- dat[dat$region %in% keep_reg, ]; dat$region <- droplevels(dat$region)
  if (nlevels(dat$region) < 2) return(NULL)
  
  aov_fit <- aov(y ~ region, data = dat); res <- residuals(aov_fit)
  sh_x <- if (length(res) > 5000) sample(res, 5000) else res
  shap_p <- tryCatch(shapiro.test(sh_x)$p.value, error = function(e) NA)
  lev_p  <- tryCatch(car::leveneTest(y ~ region, data = dat)[1,"Pr(>F)"], error = function(e) NA)
  parametric <- is.finite(shap_p) && shap_p >= 0.05 && is.finite(lev_p) && lev_p >= 0.05
  
  pmat_df <- NULL
  if (parametric) {
    method_lab <- "Welch's ANOVA + Games-Howell"
    if (have_rst) {
      ph <- tryCatch(rstatix::games_howell_test(dat, y ~ region), error = function(e) NULL)
      if (!is.null(ph)) pmat_df <- as.data.frame(ph[, c("group1","group2","p.adj")])
    }
  } else {
    method_lab <- "Kruskal-Wallis + Dunn (BH)"
    if (have_rst) {
      ph <- tryCatch(rstatix::dunn_test(dat, y ~ region, p.adjust.method = "BH"), error = function(e) NULL)
      if (!is.null(ph)) pmat_df <- as.data.frame(ph[, c("group1","group2","p.adj")])
    }
  }
  
  cld_map <- NULL
  if (have_mc && !is.null(pmat_df) && nrow(pmat_df)) {
    grps <- levels(dat$region)
    pmat <- matrix(NA_real_, length(grps), length(grps), dimnames = list(grps, grps))
    for (i in seq_len(nrow(pmat_df))) {
      a <- as.character(pmat_df$group1[i]); b <- as.character(pmat_df$group2[i])
      if (a %in% grps && b %in% grps) { pmat[a,b] <- pmat_df$p.adj[i]; pmat[b,a] <- pmat_df$p.adj[i] }
    }
    diag(pmat) <- 1
    cld <- tryCatch(multcompView::multcompLetters(pmat, threshold = 0.05, compare = "<")$Letters,
                    error = function(e) NULL)
    if (!is.null(cld)) cld_map <- data.frame(region = names(cld), cld = unname(cld), stringsAsFactors = FALSE)
  }
  if (is.null(cld_map)) cld_map <- data.frame(region = levels(dat$region), cld = "", stringsAsFactors = FALSE)
  
  # PC1 -> UPPER, PC2 -> lower
  cld_map$cld <- if (upper) toupper(cld_map$cld) else tolower(cld_map$cld)
  cld_map$n <- as.integer(table(dat$region)[cld_map$region])
  list(dat = dat, cld = cld_map, method = method_lab)
}

paired_region_plot <- function(scale_name, scale_col) {
  r1 <- region_cld_only("PC1", scale_col, upper = TRUE)
  r2 <- region_cld_only("PC2", scale_col, upper = FALSE)
  if (is.null(r1) || is.null(r2)) {
    cat(sprintf("[paired:%s] insufficient groups; skipping.\n", scale_name)); return(invisible())
  }
  
  # shared x-order = region median |latitude| (tropical -> temperate), so the
  # paired panels read equator-to-pole and match the single boxplots.
  if (LATCOL %in% names(dp)) {
    reg_abslat <- tapply(abs(dp[[LATCOL]]),
                         as.character(dp[[scale_col]]), median, na.rm = TRUE)
    ord <- levels(r1$dat$region)
    ord <- ord[order(reg_abslat[ord], na.last = TRUE)]
  } else {
    med1 <- tapply(r1$dat$y, r1$dat$region, median, na.rm = TRUE)   # fallback: PC1 median
    ord  <- names(sort(med1))
  }
  
  # keep only regions present in BOTH PCs (should be identical, but guard)
  ord <- ord[ord %in% levels(r2$dat$region)]
  
  d1 <- r1$dat; d1$PC <- "PC1"
  d2 <- r2$dat; d2$PC <- "PC2"
  both <- rbind(d1, d2)
  both <- both[both$region %in% ord, ]
  both$region <- factor(both$region, levels = ord)
  both$PC <- factor(both$PC, levels = c("PC1","PC2"))
  
  c1 <- r1$cld; c1$PC <- "PC1"
  c2 <- r2$cld; c2$PC <- "PC2"
  cld <- rbind(c1, c2)
  cld <- cld[cld$region %in% ord, ]
  cld$region <- factor(cld$region, levels = ord)
  cld$PC <- factor(cld$PC, levels = c("PC1","PC2"))
  # letter y-position = per-PC panel max + headroom
  ytop <- tapply(both$y, list(both$region, both$PC), max, na.rm = TRUE)
  rng  <- tapply(both$y, both$PC, function(v) diff(range(v, na.rm = TRUE)))
  cld$y <- mapply(function(rg, pc) ytop[rg, pc] + 0.05 * rng[pc],
                  as.character(cld$region), as.character(cld$PC))
  
  wide <- length(ord) > 15
  fill_cols <- if (scale_name == "realm") realm_cols else build_scale_cols(both, scale_col)
  
  g <- ggplot(both, aes(x = region, y = y)) +
    geom_boxplot(aes(fill = region), outlier.size = 0.5, width = 0.65,
                 color = "grey25", linewidth = 0.35, show.legend = FALSE) +
    geom_text(data = cld, aes(x = region, y = y, label = cld),
              inherit.aes = FALSE, size = if (wide) 2.6 else 3.3, fontface = "bold") +
    geom_text(data = cld, aes(x = region, y = -Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -0.6, size = if (wide) 1.9 else 2.4, color = "grey40") +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    facet_grid(PC ~ ., scales = "free_y", switch = "y") +
    labs(x = NULL, y = NULL,
         title = sprintf("PC1 & PC2 by %s", scale_disp(scale_name)),
         subtitle = sprintf("x-order by region median |lat|; PC1 = CAPITAL letters, PC2 = lowercase (letters compared WITHIN each panel only)\nPC1: %s | PC2: %s",
                            r1$method, r2$method)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = if (wide) 6 else 9),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8, color = "grey30"),
          strip.placement = "outside", strip.text.y = element_text(face = "bold"))
  
  w <- if (wide) max(10, length(ord) * 0.42) else 8
  ggsave(file.path(OUT, sprintf("region_boxplot_%s_PC1_PC2_paired.png", scale_name)),
         g, width = w, height = 9, dpi = 150, limitsize = FALSE)
  cat(sprintf("[paired] wrote region_boxplot_%s_PC1_PC2_paired.png\n", scale_name))
}

for (sc in names(meow_cols))
  paired_region_plot(sc, meow_cols[[sc]])

# =============================================================================
# SPATIAL ANALYSIS - CLIMATE ZONES
# =============================================================================
#   Tests for differences in the distribution of PC1 and PC2 across climate
#   zones (tropical vs temperate). A region is labeled "Tropical" if >50% of
#   its specimens fall below TROPIC_BOUND (23.5 deg) |latitude|, else
#   "Temperate"; specimens are then POOLED by their region's label and the two
#   zones compared. Assumption-gated: Welch's t-test (parametric) or
#   Mann-Whitney U (non-parametric), with an effect size (Cohen's d or
#   rank-biserial r). First as stats-only across all MEOW scales, then as a
#   plotted CLD boxplot at the ecoregion scale.

TROPIC_BOUND <- 23.5

## --- stats-only: tropical vs temperate by region-majority, per MEOW scale ----
thermal_region_test <- function(PC, scale_name, scale_col) {
  if (!(LATCOL %in% names(dp))) { cat("[thermal] no latitude column; skipping.\n"); return(invisible()) }
  base <- data.frame(y = dp[[PC]], region = as.character(dp[[scale_col]]),
                     abslat = abs(dp[[LATCOL]]), stringsAsFactors = FALSE)
  base <- base[is.finite(base$y) & is.finite(base$abslat) &
                 !is.na(base$region) & !(trimws(base$region) %in% c("NA","")), ]
  # restrict to the same n>=10 regions the boxplots use
  keep <- names(table(base$region))[table(base$region) >= 10]
  base <- base[base$region %in% keep, ]
  if (!nrow(base)) { cat(sprintf("[thermal:%s:%s] no regions >=10; skipping.\n", scale_name, PC)); return(invisible()) }
  
  # per-region majority label: >50% of specimens below the tropic line
  frac_trop <- tapply(base$abslat < TROPIC_BOUND, base$region, mean)
  reg_label <- ifelse(frac_trop > 0.5, "Tropical", "Temperate")
  base$class <- factor(reg_label[base$region], levels = c("Tropical","Temperate"))
  
  n_reg_trop <- sum(reg_label == "Tropical"); n_reg_temp <- sum(reg_label == "Temperate")
  if (nlevels(droplevels(base$class)) < 2) {
    cat(sprintf("[thermal:%s:%s] all regions one class (%d trop / %d temp); skipping.\n",
                scale_name, PC, n_reg_trop, n_reg_temp)); return(invisible())
  }
  
  yt <- base$y[base$class == "Tropical"]; ye <- base$y[base$class == "Temperate"]
  nt <- length(yt); ne <- length(ye)
  
  # ---- assumption checks ----
  aov_fit <- aov(y ~ class, data = base); res <- residuals(aov_fit)
  sh_x <- if (length(res) > 5000) sample(res, 5000) else res
  shap_p <- tryCatch(shapiro.test(sh_x)$p.value, error = function(e) NA)
  lev_p  <- tryCatch(car::leveneTest(y ~ class, data = base)[1,"Pr(>F)"], error = function(e) NA)
  parametric <- is.finite(shap_p) && shap_p >= 0.05 && is.finite(lev_p) && lev_p >= 0.05
  
  cat(sprintf("\n===== %s | %s: tropical vs temperate (by region majority) =====\n", PC, scale_name))
  cat(sprintf("  regions: %d tropical, %d temperate | specimens: %d tropical, %d temperate\n",
              n_reg_trop, n_reg_temp, nt, ne))
  cat(sprintf("  Shapiro p = %s | Levene p = %s\n",
              format.pval(shap_p, digits = 3, eps = 1e-4),
              format.pval(lev_p,  digits = 3, eps = 1e-4)))
  
  # ---- two-group test: Welch t (parametric) / Mann-Whitney U (non-para) ----
  if (parametric) {
    method_lab <- "Welch's t-test"
    tt <- t.test(y ~ class, data = base, var.equal = FALSE)
    stat <- unname(tt$statistic); dfr <- unname(tt$parameter); pval <- tt$p.value
    sp <- sqrt(((nt-1)*var(yt) + (ne-1)*var(ye)) / (nt + ne - 2))
    eff <- (mean(yt) - mean(ye)) / sp; eff_lab <- "Cohen's d"
    cat(sprintf("  -> parametric. Welch t = %.3f, df = %.1f, p = %s | Cohen's d = %.3f\n",
                stat, dfr, format.pval(pval, digits = 3, eps = 1e-4), eff))
  } else {
    method_lab <- "Mann-Whitney U"
    wt <- suppressWarnings(wilcox.test(y ~ class, data = base))
    stat <- unname(wt$statistic); dfr <- NA; pval <- wt$p.value
    eff <- 1 - (2 * stat) / (nt * ne); eff_lab <- "rank-biserial r"
    cat(sprintf("  -> non-parametric. Mann-Whitney U = %.0f, p = %s | rank-biserial r = %.3f\n",
                stat, format.pval(pval, digits = 3, eps = 1e-4), eff))
  }
  
  REGION_STATS[[length(REGION_STATS)+1]] <<- data.frame(
    PC = PC, scale = paste0("thermal_", scale_name), groups = 2, n = nrow(base),
    test = method_lab, statistic = round(stat,3),
    df1 = if (is.na(dfr)) NA else round(dfr,1), df2 = NA,
    p = signif(pval,3), shapiro_p = signif(shap_p,3), levene_p = signif(lev_p,3),
    row.names = NULL, stringsAsFactors = FALSE)
  
  THERMAL_EFFECTS[[length(THERMAL_EFFECTS)+1]] <<- data.frame(
    PC = PC, scale = scale_name,
    regions_tropical = n_reg_trop, regions_temperate = n_reg_temp,
    n_tropical = nt, n_temperate = ne,
    median_tropical = round(median(yt),3), median_temperate = round(median(ye),3),
    test = method_lab, effect_metric = eff_lab, effect_size = round(eff,3),
    p = signif(pval,3), row.names = NULL, stringsAsFactors = FALSE)
}

THERMAL_EFFECTS <- list()
for (sc in names(meow_cols))
  for (PC in c("PC1","PC2"))
    thermal_region_test(PC, sc, meow_cols[[sc]])

## --- region-majority thermal classification table (which region -> which label) ----
thermal_labels <- do.call(rbind, lapply(names(meow_cols), function(sc) {
  scale_col <- meow_cols[[sc]]
  b <- data.frame(region = as.character(dp[[scale_col]]), abslat = abs(dp[[LATCOL]]),
                  stringsAsFactors = FALSE)
  b <- b[is.finite(b$abslat) & !is.na(b$region) & !(trimws(b$region) %in% c("NA","")), ]
  keep <- names(table(b$region))[table(b$region) >= 10]
  b <- b[b$region %in% keep, ]
  if (!nrow(b)) return(NULL)
  ft <- tapply(b$abslat < TROPIC_BOUND, b$region, mean)
  nn <- table(b$region)
  data.frame(scale = sc, region = names(ft),
             n = as.integer(nn[names(ft)]),
             pct_below_tropic = round(100*ft,1),
             label = ifelse(ft > 0.5, "Tropical", "Temperate"),
             row.names = NULL, stringsAsFactors = FALSE)
}))
if (!is.null(thermal_labels)) {
  write.csv(thermal_labels, file.path(OUT, "thermal_region_labels.csv"), row.names = FALSE)
  cat("\n=== Region thermal labels (majority <23.5 deg) ===\n")
  print(thermal_labels, row.names = FALSE)
}

thermal_eff_tab <- do.call(rbind, THERMAL_EFFECTS)
if (!is.null(thermal_eff_tab) && nrow(thermal_eff_tab)) {
  write.csv(thermal_eff_tab, file.path(OUT, "thermal_class_effects.csv"), row.names = FALSE)
  cat("\n=== Tropical vs temperate tests (per scale) ===\n")
  print(thermal_eff_tab, row.names = FALSE)
}

## --- plotted companion: climate-zone CLD boxplot (ecoregion majority) ----
#   The plotted companion to the ecoregion thermal test. Climate zone is a
#   PER-ECOREGION label (>50% of specimens below TROPIC_BOUND); specimens are
#   pooled by their ecoregion's zone and PC1/PC2 compared across the two zones
#   in the same CLD-boxplot style as the region section (assumption-gated
#   Welch t / Mann-Whitney; shared-letter CLD). Boxes ordered tropical (low
#   |lat|) left, temperate (high |lat|) right.
climate_zone_boxplot <- function(PC, scale_col = meow_cols[["ecoregion"]]) {
  if (is.null(scale_col) || is.na(scale_col)) {
    cat("[climzone] no ecoregion column; skipping.\n"); return(invisible())
  }
  if (!(LATCOL %in% names(dp))) { cat("[climzone] no latitude column; skipping.\n"); return(invisible()) }
  
  base <- data.frame(y = dp[[PC]], region = as.character(dp[[scale_col]]),
                     abslat = abs(dp[[LATCOL]]), stringsAsFactors = FALSE)
  base <- base[is.finite(base$y) & is.finite(base$abslat) &
                 !is.na(base$region) & !(trimws(base$region) %in% c("NA","")), ]
  # same n>=10 ecoregion set the boxplots + thermal tests use
  keep <- names(table(base$region))[table(base$region) >= 10]
  base <- base[base$region %in% keep, ]
  if (!nrow(base)) { cat(sprintf("[climzone:%s] no ecoregions >=10; skipping.\n", PC)); return(invisible()) }
  
  # per-ecoregion majority label, then pool specimens by their region's zone
  frac_trop <- tapply(base$abslat < TROPIC_BOUND, base$region, mean)
  reg_label <- ifelse(frac_trop > 0.5, "Tropical", "Temperate")
  base$zone <- factor(reg_label[base$region], levels = c("Tropical","Temperate"))
  base$zone <- droplevels(base$zone)
  if (nlevels(base$zone) < 2) {
    cat(sprintf("[climzone:%s] only one zone present; skipping.\n", PC)); return(invisible())
  }
  
  # ---- assumption checks (mirror the region section) ----
  aov_fit <- aov(y ~ zone, data = base); res <- residuals(aov_fit)
  sh_x <- if (length(res) > 5000) sample(res, 5000) else res
  shap_p <- tryCatch(shapiro.test(sh_x)$p.value, error = function(e) NA)
  lev_p  <- tryCatch(car::leveneTest(y ~ zone, data = base)[1,"Pr(>F)"], error = function(e) NA)
  parametric <- is.finite(shap_p) && shap_p >= 0.05 && is.finite(lev_p) && lev_p >= 0.05
  
  yt <- base$y[base$zone == "Tropical"]; ye <- base$y[base$zone == "Temperate"]
  nt <- length(yt); ne <- length(ye)
  
  cat(sprintf("\n===== %s | climate zone (ecoregion majority) =====\n", PC))
  cat(sprintf("  specimens: %d tropical, %d temperate | Shapiro p = %s | Levene p = %s\n",
              nt, ne, format.pval(shap_p, digits = 3, eps = 1e-4),
              format.pval(lev_p, digits = 3, eps = 1e-4)))
  
  # ---- two-group test: Welch t (parametric) / Mann-Whitney U (non-para) ----
  if (parametric) {
    method_lab <- "Welch's t-test"
    tt <- t.test(y ~ zone, data = base, var.equal = FALSE)
    pval <- tt$p.value
    cat(sprintf("  -> parametric. Welch t = %.3f, df = %.1f, p = %s\n",
                unname(tt$statistic), unname(tt$parameter),
                format.pval(pval, digits = 3, eps = 1e-4)))
  } else {
    method_lab <- "Mann-Whitney U"
    wt <- suppressWarnings(wilcox.test(y ~ zone, data = base))
    pval <- wt$p.value
    cat(sprintf("  -> non-parametric. Mann-Whitney U = %.0f, p = %s\n",
                unname(wt$statistic), format.pval(pval, digits = 3, eps = 1e-4)))
  }
  
  # ---- CLD: two groups -> shared letter if n.s., distinct if p<.05 ----
  if (is.finite(pval) && pval < 0.05) {
    cld_map <- data.frame(zone = c("Tropical","Temperate"), cld = c("a","b"),
                          stringsAsFactors = FALSE)
  } else {
    cld_map <- data.frame(zone = c("Tropical","Temperate"), cld = c("a","a"),
                          stringsAsFactors = FALSE)
  }
  
  # ---- order the two boxes by median |latitude| of their pooled specimens ----
  zone_medlat <- tapply(base$abslat, base$zone, median, na.rm = TRUE)
  zord <- names(sort(zone_medlat))           # low |lat| (tropical) -> high (temperate)
  base$zone     <- factor(base$zone,     levels = zord)
  cld_map <- cld_map[match(zord, cld_map$zone), ]
  cld_map$zone <- factor(cld_map$zone, levels = zord)
  
  cld_map$y <- tapply(base$y, base$zone, max, na.rm = TRUE)[cld_map$zone] +
    0.05 * diff(range(base$y, na.rm = TRUE))
  cld_map$n <- as.integer(table(base$zone)[cld_map$zone])
  
  write.csv(
    data.frame(zone = as.character(cld_map$zone), cld = cld_map$cld, n = cld_map$n,
               median_abslat = round(unname(zone_medlat[as.character(cld_map$zone)]), 2),
               PC = PC, test = method_lab, p = signif(pval, 3)),
    file.path(OUT, sprintf("climate_zone_ecoregion_%s_cld.csv", PC)), row.names = FALSE)
  
  # ---- boxplot (darkest tropical-red / temperate-blue from the region ramps) ----
  zone_cols <- c(Tropical = RED_RAMP[1], Temperate = BLUE_RAMP[1])
  
  g <- ggplot(base, aes(x = zone, y = y)) +
    geom_boxplot(aes(fill = zone), outlier.size = 0.5, width = 0.55,
                 color = "grey25", linewidth = 0.35, show.legend = FALSE) +
    geom_text(data = cld_map, aes(x = zone, y = y, label = cld),
              inherit.aes = FALSE, size = 4.5, fontface = "bold") +
    geom_text(data = cld_map, aes(x = zone, y = -Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -0.6, size = 3, color = "grey40") +
    scale_fill_manual(values = zone_cols, drop = FALSE) +
    labs(x = NULL, y = PC,
         title = sprintf("%s by climate zone (MEOW ecoregion majority)", PC),
         subtitle = sprintf("%s; zone = >50%% of ecoregion's specimens below %.1f\u00b0 |lat|; shared letter = not different (p<.05)",
                            method_lab, TROPIC_BOUND)) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(size = 11, face = "bold"),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 8.5, color = "grey30"))
  
  ggsave(file.path(OUT, sprintf("climate_zone_ecoregion_boxplot_%s.png", PC)),
         g, width = 5.5, height = 6, dpi = 150)
  cat(sprintf("[climzone] wrote climate_zone_ecoregion_boxplot_%s.png + cld.csv\n", PC))
}

for (PC in c("PC1","PC2"))
  climate_zone_boxplot(PC)

# =============================================================================
# ENVIRONMENTAL ANALYSIS - predictor collinearity setup
# =============================================================================
#   Diagnoses and prunes the collinearity structure among env predictors
#   before regression. A correlation heatmap documents which predictors are
#   redundant; priority-ordered pruning drops any predictor with |r| > 0.9
#   against an already-kept predictor. Runs on the env subsample (predictors
#   are only defined there), matching the regression population exactly.

## --- environmental predictor correlation heatmap ----
#   Pairwise correlation matrix ACROSS the env predictors themselves -- the
#   collinearity structure behind the pruning step. Upper triangle blanked;
#   fill = Pearson r on a fixed [-1,1] diverging scale; |r|>0.9 pairs flagged
#   (these get pruned). Predictors on the LOG-ADJUSTED scale fed to the model,
#   with clean display names.
save_env_corr_heatmap <- function() {
  # Use full PREDICTORS (pre-pruning) so the reader sees why each was pruned.
  P  <- PREDICTORS
  M  <- as.matrix(d_env[, P])
  cmx <- cor(M, use = "pairwise.complete.obs")           # Pearson r
  labs_clean <- clean_lab(P)
  dimnames(cmx) <- list(labs_clean, labs_clean)
  
  # order variables by hierarchical clustering on (1 - |r|) so collinear
  # blocks sit together -- makes the pruned clusters visually obvious.
  dord <- tryCatch(hclust(as.dist(1 - abs(cmx)))$order, error = function(e) seq_along(P))
  ord_labs <- labs_clean[dord]
  
  # long form, keep lower triangle + diagonal only (blank the mirror image)
  long <- expand.grid(row = ord_labs, col = ord_labs,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  long$r <- mapply(function(a, b) cmx[a, b], long$row, long$col)
  ri <- match(long$row, ord_labs); ci <- match(long$col, ord_labs)
  long$r[ci > ri] <- NA                                   # blank upper triangle
  long$row <- factor(long$row, levels = ord_labs)
  long$col <- factor(long$col, levels = ord_labs)
  
  # label every visible tile; flag the |r|>0.9 pairs that pruning targets
  long$tile_lab <- ifelse(is.na(long$r), "", sprintf("%.2f", long$r))
  long$flag <- !is.na(long$r) & long$row != long$col & abs(long$r) > 0.9
  long$txt_col <- ifelse(!is.na(long$r) & abs(long$r) > 0.6, "white", "grey10")
  
  g <- ggplot(long, aes(x = col, y = row, fill = r)) +
    geom_tile(color = "grey30", linewidth = 0.4) +
    geom_tile(data = subset(long, flag), color = "black", linewidth = 1.1, fill = NA) +
    geom_text(aes(label = tile_lab, color = txt_col), size = 3, fontface = "bold") +
    scale_fill_gradient2(low = "powderblue", mid = "white", high = "firebrick1",
                         midpoint = 0, limits = c(-1, 1), na.value = "white",
                         name = "Pearson r") +
    scale_color_identity() +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(ord_labs)) +            # diagonal top-left -> bottom-right
    coord_equal() +
    labs(x = NULL, y = NULL,
         title = "Environmental predictor correlation structure",
         subtitle = "Log-adjusted predictors; clustered by |r|") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          panel.grid = element_blank(),
          legend.key.height = unit(1.1, "cm"))
  
  n <- length(P); sz <- max(6, 0.7 * n + 3)
  ggsave(file.path(OUT, "env_predictor_corr_heatmap.png"), g, width = sz, height = sz, dpi = 150)
  
  # full (untriangulated) matrix as CSV, in the clustered order, for the record
  write.csv(round(cmx[dord, dord], 4),
            file.path(OUT, "env_predictor_corr_matrix.csv"))
  cat("[corr] wrote env_predictor_corr_heatmap.png + env_predictor_corr_matrix.csv\n")
}
save_env_corr_heatmap()

## --- collinearity pruning of predictors (|r| > 0.9, priority order) ----
cm <- cor(d_env[, PREDICTORS], use = "pairwise.complete.obs")
kept <- character(0)
for (p in PREDICTORS) {
  if (!length(kept)) { kept <- p; next }
  if (max(abs(cm[p, kept])) > 0.9) cat(sprintf("  prune %s (collinear)\n", p))
  else kept <- c(kept, p)
}
cat("Env predictors:", paste(kept, collapse=", "), "\n")

## --- z-score the kept predictors + build the regression RHS ----
for (p in kept) d_env[[paste0(p,".z")]] <- as.numeric(scale(d_env[[p]]))
zp <- paste0(kept, ".z")
rhs <- paste(zp, collapse = " + ")

# =============================================================================
# ENVIRONMENTAL ANALYSIS - single-variable regression
# =============================================================================
#   Each predictor regressed against each PC ALONE (OLS slope, R2, Spearman
#   rho), then a scatter grid of PC vs each predictor with a fit line drawn
#   only where the slope is significant. Shows each predictor's raw marginal
#   signal before the multiple model partials it against correlated predictors.

## --- univariate regression table ----
cat("\nUnivariate regression run.\n")
univariate_tab <- function(PC) {
  rows <- lapply(kept, function(p) {
    z  <- paste0(p, ".z")
    m  <- lm(as.formula(paste(PC, "~", z)), data = d_env)
    s  <- summary(m)
    co <- s$coefficients
    sp <- suppressWarnings(cor.test(d_env[[p]], d_env[[PC]], method = "spearman"))
    data.frame(
      PC = PC,
      predictor = if (p %in% names(pretty_name)) pretty_name[p] else p,
      slope = round(co[2,1],4), se = round(co[2,2],4),
      t = round(co[2,3],2), p = signif(co[2,4],3), sig = star(co[2,4]),
      R2 = round(s$r.squared,4),
      spearman_rho = round(unname(sp$estimate),3),
      spearman_p = signif(sp$p.value,3),
      n = length(m$residuals),
      row.names = NULL, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
uni <- rbind(univariate_tab("PC1"), univariate_tab("PC2"))
uni <- do.call(rbind, lapply(split(uni, uni$PC), function(x) x[order(-x$R2), ]))
write.csv(uni, file.path(OUT, "pc_env_regression_univariate.csv"), row.names = FALSE)

cat("\n=== Univariate regression (each predictor regressed alone) ===\n")
print(uni[, c("PC","predictor","slope","se","t","p","sig","R2","spearman_rho")], row.names = FALSE)

## --- univariate scatter grids: PC1 & PC2 vs each env predictor (linear fit) ----
panel_fit <- function(PC, pred) {
  dd <- data.frame(y = dp[[PC]], x = dp[[pred]], realm = dp$.realm)
  dd <- dd[is.finite(dd$x) & is.finite(dd$y), ]
  n  <- nrow(dd)
  sp  <- suppressWarnings(cor.test(dd$x, dd$y, method = "spearman"))
  lin <- lm(y ~ x, data = dd); sm <- summary(lin)
  r2  <- sm$r.squared
  slope_p <- if (nrow(sm$coefficients) >= 2) sm$coefficients[2, 4] else NA
  gate <- is.finite(slope_p) && slope_p < .05
  list(dd = dd, n = n, r2 = r2, rho = unname(sp$estimate), sp_p = sp$p.value,
       slope_p = slope_p, gate = gate)
}

make_panel <- function(PC, pred, ranked = NULL) {
  ft <- if (is.null(ranked)) panel_fit(PC, pred) else ranked
  dd <- ft$dd
  line_tag <- if (ft$gate) "fit" else "no fit (ns)"
  lab <- sprintf("R2 = %.3f\nr = %.3f\np = %s\nslope p = %s\nn = %d | %s",
                 ft$r2, ft$rho, format.pval(ft$sp_p, digits = 2, eps = 1e-3),
                 if (is.na(ft$slope_p)) "NA" else format.pval(ft$slope_p, digits = 2, eps = 1e-3),
                 ft$n, line_tag)
  g <- ggplot(dd, aes(x = x, y = y)) +
    geom_point(aes(color = realm), size = 1.5, alpha = 0.8, show.legend = FALSE)
  if (ft$gate)
    g <- g + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "black", linewidth = 0.7)
  g +
    scale_color_manual(values = realm_cols, drop = FALSE) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size = 2.7, label = lab) +
    labs(x = clean_lab(pred), title = clean_lab(pred)) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 9, hjust = 0.5), axis.title.x = element_blank())
}

save_grid <- function(PC) {
  fits <- lapply(kept, function(p) panel_fit(PC, p))
  ord  <- order(vapply(fits, function(f) f$r2, numeric(1)), decreasing = TRUE)
  kk   <- kept[ord]; fits <- fits[ord]
  diag_tab <- data.frame(
    PC = PC, predictor = vapply(kk, function(p) clean_lab(p), character(1)),
    n   = vapply(fits, `[[`, numeric(1), "n"),
    R2  = round(vapply(fits, `[[`, numeric(1), "r2"), 4),
    spearman_r = round(vapply(fits, `[[`, numeric(1), "rho"), 3),
    spearman_p = signif(vapply(fits, `[[`, numeric(1), "sp_p"), 3),
    slope_p    = signif(vapply(fits, `[[`, numeric(1), "slope_p"), 3),
    line_drawn = vapply(fits, `[[`, logical(1), "gate"),
    row.names = NULL)
  write.csv(diag_tab, file.path(OUT, sprintf("scatter_%s_linear_diagnostics.csv", PC)), row.names = FALSE)
  panels <- Map(function(p, f) make_panel(PC, p, ranked = f), kk, fits)
  ncol   <- min(3, length(panels)); nrow <- ceiling(length(panels)/ncol)
  ttl    <- sprintf("%s vs environment (sorted by R2; linear fit, log-adj predictors)", PC)
  built <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    built <- Reduce(`+`, panels) + plot_layout(ncol = ncol) +
      plot_annotation(title = ttl) & theme(plot.title = element_text(size = 12))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    built <- gridExtra::arrangeGrob(grobs = panels, ncol = ncol, top = ttl)
  } else {
    for (i in seq_along(panels))
      ggsave(file.path(OUT, sprintf("scatter_%s_linear_%s.png", PC, kk[i])), panels[[i]], width = 5, height = 4, dpi = 150)
    cat("[scatter] install patchwork/gridExtra for combined panels.\n"); return(invisible())
  }
  w <- 4.2*ncol; h <- 3.6*nrow
  ggsave(file.path(OUT, sprintf("scatter_%s_vs_env_linear.png", PC)), built, width = w, height = h, dpi = 150, limitsize = FALSE)
  cat(sprintf("[scatter] wrote scatter_%s_vs_env_linear.png + diagnostics.csv\n", PC))
}
save_grid("PC1")
save_grid("PC2")

# =============================================================================
# ENVIRONMENTAL ANALYSIS - OLS multiple regression
# =============================================================================
#   Full multiple regression of each PC on all kept (z-scored) predictors,
#   with VIF collinearity diagnostics and Moran's I on the residuals to test
#   for residual spatial autocorrelation. A significant Moran's I flags that
#   OLS standard errors are anticonservative and triggers the GLS upgrade.

## --- fit + diagnostics per PC (VIF + Moran's I) ----
fit_report <- function(PC) {
  cat(sprintf("\n===================  %s ~ environment (OLS)  ===================\n", PC))
  f <- as.formula(paste(PC, "~", rhs))
  m <- lm(f, data = d_env)
  sm <- summary(m)
  # compact fit line (full coefficient table prints in the univariate section)
  cat(sprintf("  R2 = %.3f | adj R2 = %.3f | F(%d,%d) = %.1f | p = %s | n = %d\n",
              sm$r.squared, sm$adj.r.squared, sm$fstatistic[2], sm$fstatistic[3],
              sm$fstatistic[1],
              format.pval(pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3],
                             lower.tail = FALSE)), nrow(d_env)))
  
  cat("\nVIF (collinearity; >5 shaky, >10 bad):\n")
  v <- tryCatch(vif(m), error = function(e) setNames(rep(NA, length(zp)), zp))
  print(round(v, 2))
  
  coords <- as.matrix(d_env[, c(LONCOL, LATCOL)])
  dm <- as.matrix(dist(coords))
  diag(dm) <- 0
  w <- 1 / dm; diag(w) <- 0
  w[!is.finite(w)] <- 0
  mi <- tryCatch(Moran.I(residuals(m), w),
                 error = function(e) list(observed=NA, expected=NA, p.value=NA))
  cat(sprintf("\nMoran's I on residuals: I = %.4f (expected %.4f), p = %s\n",
              mi$observed, mi$expected, format.pval(mi$p.value)))
  if (!is.na(mi$p.value) && mi$p.value < 0.05)
    cat("  -> residuals spatially autocorrelated; OLS SEs are anticonservative.\n",
        "     GLS spatial upgrade fit below.\n")
  else
    cat("  -> no significant residual autocorrelation; OLS inference stands.\n")
  
  list(model = m, moran = mi, formula = f)
}

res1 <- fit_report("PC1")
res2 <- fit_report("PC2")

## --- tidy OLS multiple-regression table ----
tidy_ols <- function(res, PC) {
  s   <- summary(res$model)$coefficients
  sm  <- summary(res$model)
  terms <- rownames(s)
  pred  <- sub("\\.z$", "", terms)
  pretty <- ifelse(pred %in% names(pretty_name), pretty_name[pred], pred)
  data.frame(
    PC = PC, term = terms, predictor = pretty,
    estimate = round(s[,1],4), se = round(s[,2],4),
    t = round(s[,3],2), p = signif(s[,4],3), sig = star(s[,4]),
    model_R2 = round(sm$r.squared,4),
    model_adjR2 = round(sm$adj.r.squared,4),
    model_F = round(sm$fstatistic[1],2),
    model_df = paste(sm$fstatistic[2], sm$fstatistic[3], sep="/"),
    n = length(res$model$residuals),
    moran_p = signif(res$moran$p.value,3),
    row.names = NULL, stringsAsFactors = FALSE)
}

mult <- rbind(tidy_ols(res1,"PC1"), tidy_ols(res2,"PC2"))
cat("\n=== OLS multiple regression (no spatial autocorrelation correction) ===\n")
print(mult[, c("PC","predictor","estimate","se","t","p","sig","model_adjR2")], row.names = FALSE)

# =============================================================================
# ENVIRONMENTAL ANALYSIS - GLS multiple regression
# =============================================================================
#   Spatially-corrected refit, run ONLY where Moran's I flagged residual
#   autocorrelation. Fits GLS with an exponential correlation structure
#   (corExp, nugget = TRUE) on jittered coordinates, from several starting
#   ranges, keeping the highest-logLik fit whose range did not collapse to
#   zero (which would silently reduce GLS to OLS). Reports the fitted range /
#   nugget and warns if the GLS SEs are indistinguishable from OLS.

gls_upgrade <- function(res, PC) {
  mi <- res$moran
  if (is.na(mi$p.value) || mi$p.value >= 0.05) return(NULL)
  cat(sprintf("\n[GLS spatial fit] %s: ", PC))
  set.seed(42)                                    # reproducible jitter
  d_env$.jit_lon <- d_env[[LONCOL]] + rnorm(nrow(d_env), 0, 1e-4)
  d_env$.jit_lat <- d_env[[LATCOL]] + rnorm(nrow(d_env), 0, 1e-4)
  
  # Multi-start: the corExp range can collapse to ~0 (a bad local optimum) from
  # the default starting value, silently reducing GLS to OLS. Fit from several
  # starting ranges (degrees) and keep the fit with the best (highest) logLik,
  # among those whose range didn't collapse.
  start_ranges <- c(1, 5, 10, 25, 50)
  fit_one <- function(r0) tryCatch(
    gls(res$formula, data = d_env,
        correlation = corExp(value = c(r0, 0.1), form = ~ .jit_lon + .jit_lat,
                             nugget = TRUE),
        method = "REML"),
    error = function(e) NULL)
  fits <- Filter(Negate(is.null), lapply(start_ranges, fit_one))
  # keep fits whose fitted range is not essentially zero (> 1e-3 deg)
  good <- Filter(function(g) {
    cs <- tryCatch(coef(g$modelStruct$corStruct, unconstrained = FALSE),
                   error = function(e) NULL)
    !is.null(cs) && "range" %in% names(cs) && cs[["range"]] > 1e-3
  }, fits)
  pick <- if (length(good)) good else fits          # fall back to any fit
  if (!length(pick)) { cat("  GLS failed: no starting value converged\n"); return(NULL) }
  g <- pick[[ which.max(vapply(pick, logLik, numeric(1))) ]]   # best likelihood
  
  # ---- spatial-structure diagnostic ----
  # If corExp fails to fit (range -> 0), GLS silently collapses to OLS. Extract
  # the fitted range + nugget and flag when the SEs are unchanged from OLS.
  cs <- tryCatch(coef(g$modelStruct$corStruct, unconstrained = FALSE),
                 error = function(e) NULL)
  rng <- if (!is.null(cs) && "range"  %in% names(cs)) unname(cs["range"])  else NA
  nug <- if (!is.null(cs) && "nugget" %in% names(cs)) unname(cs["nugget"]) else NA
  ols_se  <- summary(res$model)$coefficients[, "Std. Error"]
  gls_se  <- summary(g)$tTable[, "Std.Error"]
  se_delta <- max(abs(gls_se - ols_se[names(gls_se)]), na.rm = TRUE)
  cat(sprintf("  [corExp] fitted range = %.4g deg, nugget = %.4g | max |GLS-OLS SE| = %.3g\n",
              rng, nug, se_delta))
  if (is.finite(se_delta) && se_delta < 1e-6)
    cat("  [WARN] GLS SEs are identical to OLS -> spatial structure did NOT fit",
        "(range collapsed). These p-values are effectively OLS, not spatially corrected.\n")
  
  tt <- summary(g)$tTable   # kept for the returned data.frame (printed later in the GLS table)
  data.frame(PC = PC, term = rownames(tt),
             estimate = tt[,"Value"], se = tt[,"Std.Error"],
             t = tt[,"t-value"], p = tt[,"p-value"],
             row.names = NULL, stringsAsFactors = FALSE)
}

gls1 <- gls_upgrade(res1, "PC1")
gls2 <- gls_upgrade(res2, "PC2")

## --- tidy GLS (spatially-corrected) table ----
gls_tidy <- function(gls_df, PC) {
  if (is.null(gls_df)) return(NULL)
  pred   <- sub("\\.z$", "", gls_df$term)
  pretty <- ifelse(pred %in% names(pretty_name), pretty_name[pred], pred)
  data.frame(PC = PC, term = gls_df$term, predictor = pretty,
             estimate = round(gls_df$estimate,4), se = round(gls_df$se,4),
             t = round(gls_df$t,2), p = signif(gls_df$p,3), sig = star(gls_df$p),
             row.names = NULL, stringsAsFactors = FALSE)
}

gls_out <- do.call(rbind, list(gls_tidy(gls1,"PC1"), gls_tidy(gls2,"PC2")))
if (!is.null(gls_out)) {
  cat("\n=== GLS multiple regression (spatial autocorrelation correction) ===\n")
  print(gls_out[, c("PC","predictor","estimate","se","t","p","sig")], row.names = FALSE)
} else {
  cat("No GLS table: residuals were spatially clean, OLS inference stands.\n")
}

# =============================================================================
# ENVIRONMENTAL ANALYSIS - table exports and extra figures
# =============================================================================
#   Writes the regression tables to CSV, then draws the two summary figures:
#   a coefficient plot (standardized betas + 95% CI, GLS-or-OLS source of
#   record) and a Spearman bubble plot (raw bivariate rho), both sharing one
#   fixed predictor y-order so a reader can track each variable straight across
#   from raw correlation to partialled coefficient.

## --- regression table exports ----
write.csv(mult, file.path(OUT, "pc_env_regression_multiple.csv"), row.names = FALSE)
write.csv(mult[, c("PC","term","estimate","se","t","p","moran_p")],
          file.path(OUT, "pc_env_regression.csv"), row.names = FALSE)
cat("\nCompare univariate vs OLS multiple: a predictor strong ALONE but weak in the",
    "\nFULL model is sharing its signal with a correlated predictor.\n")

if (!is.null(gls_out)) {
  write.csv(gls_out, file.path(OUT, "pc_env_regression_gls.csv"), row.names = FALSE)
}
cat("\nWrote: pc_env_regression_univariate.csv, pc_env_regression_multiple.csv,",
    "pc_env_regression.csv, pc_env_regression_gls.csv, pc_loadings.csv ->", OUT, "\n")

## --- coefficient source for plotting: GLS where available, else OLS ----
coef_source <- function(PC) {
  gls_df <- if (PC == "PC1") gls1 else gls2
  if (!is.null(gls_df)) {
    src <- "GLS (spatial)"
    df  <- data.frame(term = gls_df$term, estimate = gls_df$estimate,
                      se = gls_df$se, p = gls_df$p)
  } else {
    res <- if (PC == "PC1") res1 else res2
    s   <- summary(res$model)$coefficients
    src <- "OLS"
    df  <- data.frame(term = rownames(s), estimate = s[,1], se = s[,2], p = s[,4])
  }
  df <- df[df$term != "(Intercept)", , drop = FALSE]
  df$PC <- PC; df$source <- src
  df$lo <- df$estimate - 1.96 * df$se
  df$hi <- df$estimate + 1.96 * df$se
  df$sig <- df$p < 0.05
  df
}

## --- shared predictor y-order (built from PC1 |beta|; used by both figures) ----
cf_all <- rbind(coef_source("PC1"), coef_source("PC2"))
cf_all$label <- clean_lab(sub("\\.z$", "", cf_all$term))
ord_src <- cf_all[cf_all$PC == "PC1", ]
SHARED_Y_LEVELS <- ord_src$label[order(abs(ord_src$estimate))]   # ascending -> largest at top

## --- coefficient plot (standardized beta +/- 95% CI, GLS-or-OLS) ----
#    - point = standardized beta, whiskers = 95% CI (est +/- 1.96 SE)
#    - significant (p<.05) in color (blue neg / red pos), non-significant grey
#    - SHARED fixed y-order + clean names, matching the bubble plot below
save_coef_plot <- function() {
  cf <- cf_all
  cf$PC <- factor(cf$PC, levels = c("PC1","PC2"))
  cf$label <- factor(cf$label, levels = SHARED_Y_LEVELS)
  cf$dir <- ifelse(!cf$sig, "n.s.",
                   ifelse(cf$estimate < 0, "neg", "pos"))
  cf$dir <- factor(cf$dir, levels = c("neg","n.s.","pos"))
  src_lab <- paste(unique(cf$source), collapse = " / ")
  g <- ggplot(cf, aes(x = estimate, y = label, color = dir)) +
    geom_vline(xintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.6) +
    geom_point(aes(size = sig)) +
    scale_color_manual(
      values = c(neg = "#3182bd", `n.s.` = "grey70", pos = "#c0392b"),
      labels = c(neg = "negative (p<.05)", `n.s.` = "n.s.", pos = "positive (p<.05)"),
      name = NULL, drop = FALSE) +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 2), guide = "none") +
    facet_wrap(~ PC) +
    labs(x = "Standardized coefficient (95% CI)", y = NULL,
         title = "Environmental drivers of morphological axes",
         subtitle = sprintf("%s multiple regression; predictors z-scored", src_lab)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT, "coef_plot_PC1_PC2.png"), g, width = 9, height = 8, dpi = 150)
  cat("[coef] wrote coef_plot_PC1_PC2.png\n")
}
save_coef_plot()

## --- Spearman bubble plot (raw bivariate rho; same y-order as coef plot) ----
#    - raw bivariate rho per predictor; bubble AREA = |rho|; stars = sig
#    - SAME shared y-order and clean names as the coef plot, so a reader can
#      track each variable straight across from raw -> controlled.
save_bubble_plot <- function() {
  make_rho <- function(PC) {
    do.call(rbind, lapply(kept, function(p) {
      sp <- suppressWarnings(cor.test(d_env[[p]], d_env[[PC]], method = "spearman"))
      data.frame(PC = PC, label = clean_lab(p),
                 rho = unname(sp$estimate), p = sp$p.value,
                 row.names = NULL, stringsAsFactors = FALSE)
    }))
  }
  rr <- rbind(make_rho("PC1"), make_rho("PC2"))
  rr$stars <- ifelse(rr$p < .001, "***", ifelse(rr$p < .01, "**",
                                                ifelse(rr$p < .05, "*", ifelse(rr$p < .1, ".", ""))))
  rr$PC <- factor(rr$PC, levels = c("PC1","PC2"))
  # SHARED y-order (same as coef plot). Any predictor not in the shared list
  # (shouldn't happen) falls to the end.
  lev <- SHARED_Y_LEVELS[SHARED_Y_LEVELS %in% rr$label]
  lev <- c(lev, setdiff(unique(rr$label), lev))
  rr$label <- factor(rr$label, levels = lev)
  rr$star_x <- rr$rho + ifelse(rr$rho >= 0, 0.04, -0.04)
  rr$star_h <- ifelse(rr$rho >= 0, 0, 1)
  g <- ggplot(rr, aes(x = rho, y = label)) +
    geom_vline(xintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(size = abs(rho), fill = rho), shape = 21, color = "grey20", stroke = 0.4) +
    geom_text(aes(x = star_x, label = stars, hjust = star_h), size = 4, vjust = 0.8) +
    scale_size_area(max_size = 9, guide = "none") +
    scale_fill_gradient2(low = "#3182bd", mid = "grey90", high = "#c0392b",
                         midpoint = 0, name = expression(Spearman~rho)) +
    facet_wrap(~ PC) +
    labs(x = expression(Spearman~rho), y = NULL,
         title = "Raw environmental correlations with morphological axes",
         subtitle = "Bubble area = |rho|; *** p<.001  ** p<.01  * p<.05  . p<.1") +
    xlim(min(rr$rho) - 0.12, max(rr$rho) + 0.12) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          strip.text = element_text(face = "bold"), legend.position = "right")
  ggsave(file.path(OUT, "spearman_bubble_PC1_PC2.png"), g, width = 9, height = 8, dpi = 150)
  cat("[bubble] wrote spearman_bubble_PC1_PC2.png\n")
}
save_bubble_plot()

cat("\n[done] coefficient plot + Spearman bubble plot written to", OUT, "\n")

# =============================================================================
# SUMMARY TEXT EXPORT
# =============================================================================
#   Consolidated, citable numbers for a methods/results section: sample sizes,
#   PCA variance + loadings, regional omnibus tests, tropical/temperate
#   contrasts, multiple-regression fit + spatial diagnostics, per-predictor
#   coefficients (GLS-or-OLS source of record), and reproducibility notes.
#   Written to publication_stats_summary.txt + companion CSVs.

sink_file <- file.path(OUT, "publication_stats_summary.txt")
sink(sink_file, split = FALSE)  # write to file only (split=TRUE also floods console)

cat("=====================================================================\n")
cat(" PUBLICATION STATISTICS SUMMARY\n")
cat(" Global Asparagopsis taxiformis morphometrics -> environment\n")
cat("=====================================================================\n")

cat("\n--- 1. SAMPLE SIZES ---\n")
cat(sprintf("  Traits analyzed (%d): %s\n", length(TRAITS), paste(TRAITS, collapse = ", ")))
cat(sprintf("  Full trait-complete population (PCA): n = %d\n", nrow(d)))
cat(sprintf("  Env-complete subsample (regression): n = %d (%.1f%% of trait pool)\n",
            nrow(d_env), 100*nrow(d_env)/nrow(d)))
cat(sprintf("  Dropped for missing environmental data: n = %d\n", nrow(d) - nrow(d_env)))

cat("\n--- 2. PCA VARIANCE STRUCTURE ---\n")
print(scree_df, row.names = FALSE)
cat(sprintf("  PC1 + PC2 retained: %.1f%% cumulative variance\n", ve[1] + ve[2]))
cat("\n  PC1/PC2 loadings:\n")
print(load)

cat("\n--- 3. REGIONAL DIFFERENTIATION (omnibus tests per MEOW scale) ---\n")
region_tab <- do.call(rbind, REGION_STATS)
if (!is.null(region_tab) && nrow(region_tab)) {
  print(region_tab, row.names = FALSE)
  write.csv(region_tab, file.path(OUT, "publication_regional_tests.csv"), row.names = FALSE)
} else cat("  (no regional tests recorded)\n")

cat("\n--- 3b. TROPICAL vs TEMPERATE (region-majority, per scale, with effect size) ---\n")
if (exists("thermal_eff_tab") && !is.null(thermal_eff_tab) && nrow(thermal_eff_tab)) {
  print(thermal_eff_tab, row.names = FALSE)
  cat("  (regions labeled tropical if >50% of specimens below 23.5 deg; specimens pooled by label)\n")
  cat("  (effect: |d| ~ 0.2 small, 0.5 medium, 0.8 large; rank-biserial |r| similar scale)\n")
} else cat("  (no thermal contrast recorded)\n")

cat("\n--- 4. MULTIPLE REGRESSION MODEL FIT (trait morphology ~ environment) ---\n")
for (r in list(res1, res2)) {
  sm <- summary(r$model); PC <- all.vars(r$formula)[1]
  cat(sprintf("  %s: R2 = %.3f, adj R2 = %.3f, F(%d,%d) = %.2f, p %s | n = %d\n",
              PC, sm$r.squared, sm$adj.r.squared,
              sm$fstatistic[2], sm$fstatistic[3], sm$fstatistic[1],
              format.pval(pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3],
                             lower.tail = FALSE), digits = 3, eps = 2.2e-16),
              length(r$model$residuals)))
  cat(sprintf("       Moran's I on residuals = %.4f, p = %s -> %s\n",
              r$moran$observed, format.pval(r$moran$p.value, digits = 3, eps = 1e-16),
              if (!is.na(r$moran$p.value) && r$moran$p.value < 0.05)
                "spatially autocorrelated; GLS reported" else "spatially clean; OLS reported"))
}

cat("\n--- 5. PREDICTOR COEFFICIENTS (source-of-record) ---\n")
cat("  (GLS spatial estimates where residuals autocorrelated, else OLS; predictors z-scored)\n")
pub_coef <- do.call(rbind, lapply(c("PC1","PC2"), function(PC) {
  cs <- coef_source(PC)             # already GLS-or-OLS, intercept dropped
  cs$predictor <- clean_lab(sub("\\.z$", "", cs$term))
  data.frame(PC = PC, predictor = cs$predictor,
             beta = round(cs$estimate,4), se = round(cs$se,4),
             p = signif(cs$p,3), sig = star(cs$p), source = cs$source,
             row.names = NULL, stringsAsFactors = FALSE)
}))
print(pub_coef, row.names = FALSE)
write.csv(pub_coef, file.path(OUT, "publication_coefficients.csv"), row.names = FALSE)

cat("\n--- 6. REPRODUCIBILITY NOTES ---\n")
cat(sprintf("  R %s; predictor log-transforms: %s\n",
            getRversion(), paste(LOG_PREDICTORS, collapse = ", ")))
cat(sprintf("  Collinearity pruning |r|>0.9; predictors retained: %s\n", paste(kept, collapse = ", ")))
cat(sprintf("  GLS spatial structure: exponential (corExp), nugget = TRUE, seed = 42\n"))
cat(sprintf("  Post-hoc adjustment: BH (non-parametric) / Games-Howell (parametric)\n"))

cat("\n=====================================================================\n")
sink()