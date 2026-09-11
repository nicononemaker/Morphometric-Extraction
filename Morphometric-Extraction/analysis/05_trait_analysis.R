# =============================================================================
# TRAIT-FACING ANALYSIS TEMPLATE
#
# For each trait you list in TRAITS, this runs the same pipeline the PC script
# runs on PC1/PC2, but with the RAW trait as the response variable:
#
#   (A) Raw bivariate correlations  : Spearman rho + linear R2 vs each predictor
#   (B) Multiple regression         : OLS (+VIF) with Moran's I spatial check,
#                                     GLS spatial upgrade if residuals autocorr.
#   (C) Regional comparison         : per MEOW scale (realm/province/ecoregion),
#                                     assumption-checked Welch-ANOVA+Games-Howell
#                                     (parametric) OR Kruskal-Wallis+Dunn/BH
#                                     (non-parametric), boxplot with CLD letters.
#
# USAGE: set TRAITS <- c("trait_a","trait_b",...). Every trait gets the full
# pipeline. Predictor set, log flags, and coloring mirror the PC script so
# outputs are directly comparable.
#
# OUTPUTS: everything writes to OUT/Trait_Analysis/<trait>/ so each trait's
# figures and tables are self-contained.
#
# DESIGN NOTES (mirrors PC_1_2_full_analysis.R):
#   - Predictor log transforms applied ONCE upstream; regression + figures share.
#   - Traits with heavy positive skew can be log-transformed via TRAIT_LOG.
#   - Multiple regression / correlations use the ENV subsample (predictors +
#     coords complete). Regional comparison uses the FULL trait-complete pool,
#     so it recovers regions dropped only for missing environmental data.
#   - abs_latitude is NOT a mult-regression predictor (coordinate proxy); a
#     standalone latitude-gradient block is provided per trait.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(car); library(ape); library(nlme)
  library(ggplot2); library(scales)
})

# =====================================================================
# USER INPUTS -- edit these
# =====================================================================
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
IN  <- AGG_CSV
OUT <- file.path(OUT_DIR, "Trait_Analysis")   # <-- subfolder
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# >>> LIST THE TRAITS TO ANALYZE HERE <<<
TRAITS <- c("length_frond_mm","max_width_mm","stipe_width_mm",
            "solidity","peak_position_frac","branched_fraction","circularity", "aspect_ratio", "area_mm2","stipe_to_length_ratio")

# traits to log-transform before analysis (heavy positive skew). Set c() for none.
TRAIT_LOG <- c("length_frond_mm","max_width_mm","stipe_width_mm", "aspect_ratio", "area_mm2")

# env predictors for the multiple regression (abs_latitude intentionally excluded)
PREDICTORS <- c("SST_clim","nitrate_clim","phosphate_clim",
                "salinity_clim","wave.height","kd490_clim",
                "PAR_mean","iron_clim","spco2_clim")

# predictors log-transformed upstream (skew-positive fields)
LOG_PREDICTORS  <- c("nitrate_clim","phosphate_clim","iron_clim",
                     "kd490_clim","wave.height", "spco2_clim")
LOG_OFFSET_FRAC <- 0.01

LATCOL <- "Latitude"; LONCOL <- "Longitude"
REALMCOL <- "realm"
MIN_N  <- 10          # min specimens/region for the regional comparison

# =====================================================================
# SETUP (mirrors the PC script)
# =====================================================================
rec <- read.csv(IN, stringsAsFactors = FALSE, check.names = FALSE)
TRAITS     <- intersect(TRAITS, names(rec))
PREDICTORS <- intersect(PREDICTORS, names(rec))
stopifnot(length(TRAITS) >= 1)
stopifnot(all(c(LATCOL, LONCOL) %in% names(rec)))

realm_candidates <- names(rec)[grepl("realm", names(rec), ignore.case = TRUE)]
REALMCOL <- if (length(realm_candidates)) realm_candidates[1] else NA

# ---- realm palette (same as PC script) ----
realm_cols <- c(
  "Atlantic Warm Water"          = "#67000d",
  "Indo-Pacific Warm Water"      = "#a50f15",
  "Tropical Atlantic"            = "#cb181d",
  "Tropical Eastern Pacific"     = "#ef3b2c",
  "Western Indo-Pacific"         = "#fb6a4a",
  "Central Indo-Pacific"         = "#fc9272",
  "Eastern Indo-Pacific"         = "#fcbba1",
  "Temperate Northern Atlantic"  = "#08306b",
  "Temperate Northern Pacific"   = "#08519c",
  "Temperate South America"      = "#2171b5",
  "Temperate Australasia"        = "#4292c6"
)

DISPLAY_NAME <- c(
  "SST_clim"            = "Sea surface temperature",
  "nitrate_clim"        = "Nitrate", "log_nitrate_clim"   = "Nitrate (log)",
  "phosphate_clim"      = "Phosphate","log_phosphate_clim"= "Phosphate (log)",
  "salinity_clim"       = "Salinity",
  "wave.height"         = "Wave height","log_wave.height" = "Wave height (log)",
  "kd490_clim"     = "Turbidity (Kd490)","log_kd490_clim" = "Turbidity (Kd490, log)",
  "PAR_mean"            = "Photosynthetically available radiation",
  "iron_clim"           = "Iron", "log_iron_clim"        = "Iron (log)",
  "spco2_clim"      = "Surface ocean pCO\u2082",
  "log_spco2_clim"  = "Surface ocean pCO\u2082 (log)",
  "abs_latitude"        = "Absolute latitude"
)

TRAIT_DISPLAY <- c(
  "length_frond_mm"    = "Frond length",  "max_width_mm"      = "Max width",
  "stipe_width_mm"     = "Stipe width",   "solidity"          = "Solidity",
  "peak_position_frac" = "Peak position", "branched_fraction" = "Branched fraction",
  "circularity"        = "Circularity"
)
trait_lab <- function(t) if (t %in% names(TRAIT_DISPLAY)) unname(TRAIT_DISPLAY[t]) else t

star <- function(p) ifelse(is.na(p), "",
                    ifelse(p < .001, "***", ifelse(p < .01, "**",
                    ifelse(p < .05, "*", ifelse(p < .1, ".", "")))))

# ---- trait-complete pool = FULL population (regional comparison uses this) ---
d <- rec
for (t in intersect(TRAIT_LOG, TRAITS)) d[[t]] <- log(d[[t]])
cc <- complete.cases(d[, TRAITS])
d  <- d[cc, ]
cat(sprintf("Trait-complete records (full pool): %d\n", nrow(d)))

# ---- predictor log transforms upstream ----
LOG_PREDICTORS <- intersect(LOG_PREDICTORS, PREDICTORS)
pretty_name <- setNames(PREDICTORS, PREDICTORS)
if (length(LOG_PREDICTORS)) {
  for (p in LOG_PREDICTORS) {
    x <- d[[p]]; pos <- x[is.finite(x) & x > 0]
    off <- if (any(x <= 0, na.rm = TRUE) && length(pos)) LOG_OFFSET_FRAC * min(pos) else 0
    newname <- paste0("log_", p)
    d[[newname]] <- log(x + off)
    PREDICTORS[PREDICTORS == p] <- newname
    pretty_name[[newname]] <- paste0("log ", p)
  }
}
clean_lab <- function(x) {
  out <- ifelse(x %in% names(DISPLAY_NAME), DISPLAY_NAME[x],
         ifelse(x %in% names(pretty_name), pretty_name[x], x))
  unname(out)
}

# ---- ENV subsample (predictors + coords complete): used by (A) and (B) ------
env_cc <- complete.cases(d[, PREDICTORS], d[, c(LATCOL, LONCOL)])
d_env  <- d[env_cc, ]
cat(sprintf("Env-complete subsample: %d of %d\n", nrow(d_env), nrow(d)))

# ---- collinearity pruning (|r|>0.9, priority order) -- shared by all traits --
cm <- cor(d_env[, PREDICTORS], use = "pairwise.complete.obs")
kept <- character(0)
for (p in PREDICTORS) {
  if (!length(kept)) { kept <- p; next }
  if (max(abs(cm[p, kept])) > 0.9) cat(sprintf("  prune %s (collinear)\n", p)) else kept <- c(kept, p)
}
cat("Env predictors kept:", paste(kept, collapse = ", "), "\n")
for (p in kept) d_env[[paste0(p,".z")]] <- as.numeric(scale(d_env[[p]]))
zp  <- paste0(kept, ".z")
rhs <- paste(zp, collapse = " + ")

# ---- realm coloring frames ----
add_realm <- function(x) {
  if (!is.na(REALMCOL) && REALMCOL %in% names(x)) {
    keep <- !is.na(x[[REALMCOL]]) & !(trimws(as.character(x[[REALMCOL]])) %in% c("NA",""))
    x <- x[keep, ]; x$.realm <- droplevels(as.factor(x[[REALMCOL]]))
  } else x$.realm <- factor("all")
  x
}
dp_env  <- add_realm(d_env)   # env figures (correlations)
dp_full <- add_realm(d)       # regional comparison
missing_lv <- setdiff(union(levels(dp_env$.realm), levels(dp_full$.realm)), names(realm_cols))
if (length(missing_lv)) realm_cols[missing_lv] <- "#BBBBBB"

have_mc  <- requireNamespace("multcompView", quietly = TRUE)
have_rst <- requireNamespace("rstatix",      quietly = TRUE)
have_pw  <- requireNamespace("patchwork",    quietly = TRUE)

# =====================================================================
# (A) RAW bivariate correlations: trait vs each predictor
# =====================================================================
trait_correlations <- function(TRAIT, odir) {
  panel_fit <- function(pred) {
    dd <- data.frame(y = dp_env[[TRAIT]], x = dp_env[[pred]], realm = dp_env$.realm)
    dd <- dd[is.finite(dd$x) & is.finite(dd$y), ]
    sp  <- suppressWarnings(cor.test(dd$x, dd$y, method = "spearman"))
    lin <- lm(y ~ x, data = dd); sm <- summary(lin)
    slope_p <- if (nrow(sm$coefficients) >= 2) sm$coefficients[2,4] else NA
    list(dd = dd, n = nrow(dd), r2 = sm$r.squared, rho = unname(sp$estimate),
         sp_p = sp$p.value, slope_p = slope_p,
         gate = is.finite(slope_p) && slope_p < .05)
  }
  fits <- lapply(kept, panel_fit)
  ord  <- order(vapply(fits, function(f) f$r2, numeric(1)), decreasing = TRUE)
  kk   <- kept[ord]; fits <- fits[ord]

  diag_tab <- data.frame(
    trait = TRAIT, predictor = vapply(kk, clean_lab, character(1)),
    n = vapply(fits, `[[`, numeric(1), "n"),
    R2 = round(vapply(fits, `[[`, numeric(1), "r2"), 4),
    spearman_rho = round(vapply(fits, `[[`, numeric(1), "rho"), 3),
    spearman_p   = signif(vapply(fits, `[[`, numeric(1), "sp_p"), 3),
    slope_p      = signif(vapply(fits, `[[`, numeric(1), "slope_p"), 3),
    line_drawn   = vapply(fits, `[[`, logical(1), "gate"),
    row.names = NULL)
  write.csv(diag_tab, file.path(odir, sprintf("%s_correlations.csv", TRAIT)), row.names = FALSE)

  panels <- Map(function(p, ft) {
    dd <- ft$dd
    lab <- sprintf("R2 = %.3f\nrho = %.3f\nslope p = %s\nn = %d",
                   ft$r2, ft$rho,
                   if (is.na(ft$slope_p)) "NA" else format.pval(ft$slope_p, digits = 2, eps = 1e-3),
                   ft$n)
    g <- ggplot(dd, aes(x = x, y = y)) +
      geom_point(aes(color = realm), size = 1.4, alpha = 0.8, show.legend = FALSE)
    if (ft$gate) g <- g + geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                                      color = "black", linewidth = 0.7)
    g + scale_color_manual(values = realm_cols, drop = FALSE) +
      annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size = 2.6, label = lab) +
      labs(x = clean_lab(p), y = trait_lab(TRAIT), title = clean_lab(p)) +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(size = 9, hjust = 0.5), axis.title.x = element_blank())
  }, kk, fits)

  ncol <- min(3, length(panels)); nrow <- ceiling(length(panels)/ncol)
  ttl  <- sprintf("%s vs environment (sorted by R2)", trait_lab(TRAIT))
  if (have_pw) {
    library(patchwork)
    built <- Reduce(`+`, panels) + plot_layout(ncol = ncol) +
      plot_annotation(title = ttl) & theme(plot.title = element_text(size = 12))
    ggsave(file.path(odir, sprintf("%s_scatter_vs_env.png", TRAIT)), built,
           width = 4.2*ncol, height = 3.6*nrow, dpi = 150, limitsize = FALSE)
  } else {
    for (i in seq_along(panels))
      ggsave(file.path(odir, sprintf("%s_scatter_%s.png", TRAIT, kk[i])),
             panels[[i]], width = 5, height = 4, dpi = 150)
  }
  cat(sprintf("  [A] correlations written for %s\n", TRAIT))
}

# =====================================================================
# (B) MULTIPLE REGRESSION: trait ~ environment  (OLS + VIF + Moran + GLS)
# =====================================================================
trait_regression <- function(TRAIT, odir) {
  f <- as.formula(paste(TRAIT, "~", rhs))
  m <- lm(f, data = d_env); sm <- summary(m)

  v <- tryCatch(vif(m), error = function(e) setNames(rep(NA, length(zp)), zp))

  coords <- as.matrix(d_env[, c(LONCOL, LATCOL)])
  dm <- as.matrix(dist(coords)); diag(dm) <- 0
  w <- 1/dm; diag(w) <- 0; w[!is.finite(w)] <- 0
  mi <- tryCatch(Moran.I(residuals(m), w),
                 error = function(e) list(observed=NA, expected=NA, p.value=NA))

  # OLS coefficient table
  s <- sm$coefficients
  ols <- data.frame(
    trait = TRAIT, term = rownames(s),
    predictor = ifelse(sub("\\.z$","",rownames(s)) %in% names(pretty_name),
                       pretty_name[sub("\\.z$","",rownames(s))], sub("\\.z$","",rownames(s))),
    estimate = round(s[,1],4), se = round(s[,2],4), t = round(s[,3],2),
    p = signif(s[,4],3), sig = star(s[,4]),
    vif = round(c(NA, v[match(rownames(s)[-1], names(v))]), 2),
    model_R2 = round(sm$r.squared,4), model_adjR2 = round(sm$adj.r.squared,4),
    moran_p = signif(mi$p.value,3), n = length(m$residuals),
    row.names = NULL, stringsAsFactors = FALSE)
  write.csv(ols, file.path(odir, sprintf("%s_regression_OLS.csv", TRAIT)), row.names = FALSE)

  # GLS spatial upgrade if residuals autocorrelated
  coef_df <- data.frame(term = rownames(s)[-1], estimate = s[-1,1], se = s[-1,2],
                        p = s[-1,4], source = "OLS")
  if (!is.na(mi$p.value) && mi$p.value < 0.05) {
    set.seed(42)
    dg <- d_env
    dg$.jit_lon <- dg[[LONCOL]] + rnorm(nrow(dg), 0, 1e-4)
    dg$.jit_lat <- dg[[LATCOL]] + rnorm(nrow(dg), 0, 1e-4)
    g <- tryCatch(gls(f, data = dg,
                      correlation = corExp(form = ~ .jit_lon + .jit_lat, nugget = TRUE),
                      method = "REML"),
                  error = function(e) { cat("   GLS failed:", conditionMessage(e), "\n"); NULL })
    if (!is.null(g)) {
      tt <- summary(g)$tTable
      gls_tab <- data.frame(
        trait = TRAIT, term = rownames(tt),
        predictor = ifelse(sub("\\.z$","",rownames(tt)) %in% names(pretty_name),
                           pretty_name[sub("\\.z$","",rownames(tt))], sub("\\.z$","",rownames(tt))),
        estimate = round(tt[,"Value"],4), se = round(tt[,"Std.Error"],4),
        t = round(tt[,"t-value"],2), p = signif(tt[,"p-value"],3), sig = star(tt[,"p-value"]),
        row.names = NULL, stringsAsFactors = FALSE)
      write.csv(gls_tab, file.path(odir, sprintf("%s_regression_GLS.csv", TRAIT)), row.names = FALSE)
      coef_df <- data.frame(term = rownames(tt)[-1], estimate = tt[-1,"Value"],
                            se = tt[-1,"Std.Error"], p = tt[-1,"p-value"], source = "GLS (spatial)")
    }
  }

  # coefficient plot (single trait)
  coef_df$label <- clean_lab(sub("\\.z$","", coef_df$term))
  coef_df$label <- factor(coef_df$label, levels = coef_df$label[order(abs(coef_df$estimate))])
  coef_df$lo <- coef_df$estimate - 1.96*coef_df$se
  coef_df$hi <- coef_df$estimate + 1.96*coef_df$se
  coef_df$dir <- ifelse(coef_df$p >= .05, "n.s.", ifelse(coef_df$estimate < 0, "neg","pos"))
  coef_df$dir <- factor(coef_df$dir, levels = c("neg","n.s.","pos"))
  gcf <- ggplot(coef_df, aes(x = estimate, y = label, color = dir)) +
    geom_vline(xintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = c(neg="#3182bd", `n.s.`="grey70", pos="#c0392b"),
                       labels = c(neg="negative (p<.05)", `n.s.`="n.s.", pos="positive (p<.05)"),
                       name = NULL, drop = FALSE) +
    labs(x = "Standardized coefficient (95% CI)", y = NULL,
         title = sprintf("Environmental drivers of %s", trait_lab(TRAIT)),
         subtitle = sprintf("%s; predictors z-scored", coef_df$source[1])) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
  ggsave(file.path(odir, sprintf("%s_coef_plot.png", TRAIT)), gcf, width = 7, height = 6, dpi = 150)
  cat(sprintf("  [B] regression (%s) written for %s\n", coef_df$source[1], TRAIT))
}

# =====================================================================
# (C) REGIONAL COMPARISON across MEOW scales (full trait pool)
# =====================================================================
MEOW_PATTERNS <- c(realm = "realm", province = "province", ecoregion = "ecoregion")
meow_cols_map <- sapply(MEOW_PATTERNS, function(pat) {
  hit <- names(dp_full)[grepl(pat, names(dp_full), ignore.case = TRUE)]
  if (length(hit)) hit[1] else NA_character_
})
meow_cols_map <- meow_cols_map[!is.na(meow_cols_map)]

TROPICAL_REALMS <- c("Atlantic Warm Water","Indo-Pacific Warm Water","Tropical Atlantic",
                     "Tropical Eastern Pacific","Western Indo-Pacific","Central Indo-Pacific",
                     "Eastern Indo-Pacific")
TEMPERATE_REALMS <- c("Temperate Northern Atlantic","Temperate Northern Pacific",
                      "Temperate South America","Temperate Australasia","Southern Ocean")
RED_RAMP  <- c("#67000d","#a50f15","#cb181d","#ef3b2c","#fb6a4a","#fc9272","#fcbba1")
BLUE_RAMP <- c("#08306b","#08519c","#2171b5","#4292c6","#6baed6","#9ecae1")
REALMCOL_MEOW <- { h <- names(dp_full)[grepl("realm", names(dp_full), ignore.case = TRUE)]
                   if (length(h)) h[1] else NA_character_ }

build_scale_cols <- function(dat, scale_col) {
  regs <- levels(dat$region)
  if (is.na(REALMCOL_MEOW)) return(setNames(rep("#BBBBBB", length(regs)), regs))
  parent <- sapply(regs, function(r) {
    rl <- dp_full[[REALMCOL_MEOW]][ as.character(dp_full[[scale_col]]) == r ]
    rl <- rl[!is.na(rl) & !(trimws(as.character(rl)) %in% c("NA",""))]
    if (!length(rl)) return(NA_character_); names(sort(table(rl), decreasing = TRUE))[1]
  })
  latc <- if ("Latitude" %in% names(dp_full)) "Latitude" else NA
  reg_lat <- sapply(regs, function(r) if (is.na(latc)) NA_real_ else
    median(abs(dp_full[[latc]][ as.character(dp_full[[scale_col]]) == r ]), na.rm = TRUE))
  cls <- ifelse(parent %in% TROPICAL_REALMS, "trop", ifelse(parent %in% TEMPERATE_REALMS, "temp","other"))
  out <- setNames(rep("#BBBBBB", length(regs)), regs)
  assign_ramp <- function(members, ramp, ascending = TRUE) {
    if (!length(members)) return(invisible())
    members <- members[order(reg_lat[members], decreasing = !ascending, na.last = TRUE)]
    idx <- if (length(members) == 1) 1 else round(seq(1, length(ramp), length.out = length(members)))
    out[members] <<- ramp[idx]
  }
  assign_ramp(regs[cls=="trop"], RED_RAMP, TRUE)
  assign_ramp(regs[cls=="temp"], BLUE_RAMP, FALSE)
  out
}

region_test_plot <- function(TRAIT, scale_name, scale_col, odir) {
  dat <- data.frame(y = dp_full[[TRAIT]], region = as.character(dp_full[[scale_col]]),
                    stringsAsFactors = FALSE)
  dat <- dat[is.finite(dat$y) & !is.na(dat$region) & !(trimws(dat$region) %in% c("NA","")), ]
  dat$region <- droplevels(as.factor(dat$region))
  keep_reg <- names(table(dat$region))[table(dat$region) >= MIN_N]
  dat <- dat[dat$region %in% keep_reg, ]; dat$region <- droplevels(dat$region)
  ng <- nlevels(dat$region)
  if (ng < 2) { cat(sprintf("   [%s:%s] <2 groups; skip.\n", scale_name, TRAIT)); return(invisible()) }

  aov_fit <- aov(y ~ region, data = dat); res <- residuals(aov_fit)
  sh_x <- if (length(res) > 5000) sample(res, 5000) else res
  shap_p <- tryCatch(shapiro.test(sh_x)$p.value, error = function(e) NA)
  lev_p  <- tryCatch(car::leveneTest(y ~ region, data = dat)[1,"Pr(>F)"], error = function(e) NA)
  parametric <- is.finite(shap_p) && shap_p >= 0.05 && is.finite(lev_p) && lev_p >= 0.05

  pmat_df <- NULL
  if (parametric) {
    method_lab <- "Welch's ANOVA + Games-Howell"
    om <- oneway.test(y ~ region, data = dat, var.equal = FALSE)
    omni_p <- om$p.value
    if (have_rst) {
      ph <- tryCatch(rstatix::games_howell_test(dat, y ~ region), error = function(e) NULL)
      if (!is.null(ph)) pmat_df <- as.data.frame(ph[, c("group1","group2","p.adj")])
    } else {
      pt <- pairwise.t.test(dat$y, dat$region, p.adjust.method = "holm", pool.sd = FALSE)
      pm <- pt$p.value
      pmat_df <- do.call(rbind, lapply(rownames(pm), function(r) do.call(rbind, lapply(colnames(pm),
        function(c) if (!is.na(pm[r,c])) data.frame(group1=c, group2=r, p.adj=pm[r,c]) else NULL))))
    }
  } else {
    method_lab <- "Kruskal-Wallis + Dunn (BH)"
    kw <- kruskal.test(y ~ region, data = dat); omni_p <- kw$p.value
    if (have_rst) {
      ph <- tryCatch(rstatix::dunn_test(dat, y ~ region, p.adjust.method = "BH"), error = function(e) NULL)
      if (!is.null(ph)) pmat_df <- as.data.frame(ph[, c("group1","group2","p.adj")])
    } else {
      pw <- pairwise.wilcox.test(dat$y, dat$region, p.adjust.method = "BH")
      pm <- pw$p.value
      pmat_df <- do.call(rbind, lapply(rownames(pm), function(r) do.call(rbind, lapply(colnames(pm),
        function(c) if (!is.na(pm[r,c])) data.frame(group1=c, group2=r, p.adj=pm[r,c]) else NULL))))
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

  ytop <- tapply(dat$y, dat$region, max, na.rm = TRUE)
  cld_map$y <- ytop[cld_map$region] + 0.05 * diff(range(dat$y, na.rm = TRUE))
  cld_map$n <- as.integer(table(dat$region)[cld_map$region])
  lev_ord <- names(sort(tapply(dat$y, dat$region, median, na.rm = TRUE)))
  dat$region     <- factor(dat$region,     levels = lev_ord)
  cld_map$region <- factor(cld_map$region, levels = lev_ord)

  write.csv(cld_map[, c("region","cld","n")],
            file.path(odir, sprintf("%s_region_%s_cld.csv", TRAIT, scale_name)), row.names = FALSE)

  wide <- ng > 15
  fill_cols <- if (scale_name == "realm") realm_cols else build_scale_cols(dat, scale_col)
  g <- ggplot(dat, aes(x = region, y = y)) +
    geom_boxplot(aes(fill = region), outlier.size = 0.5, width = 0.65,
                 color = "grey25", linewidth = 0.35, show.legend = FALSE) +
    geom_text(data = cld_map, aes(x = region, y = y, label = cld), inherit.aes = FALSE,
              size = if (wide) 2.6 else 3.5, fontface = "bold") +
    geom_text(data = cld_map, aes(x = region, y = -Inf, label = paste0("n=", n)),
              inherit.aes = FALSE, vjust = -0.6, size = if (wide) 1.9 else 2.5, color = "grey40") +
    scale_fill_manual(values = fill_cols, drop = FALSE) +
    labs(x = NULL, y = trait_lab(TRAIT),
         title = sprintf("%s by %s", trait_lab(TRAIT), scale_name),
         subtitle = sprintf("%s (omnibus p = %s); shared letter = not different (p<.05)",
                            method_lab, format.pval(omni_p, digits = 2, eps = 1e-4))) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = if (wide) 6 else 9),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
  w <- if (wide) max(10, ng * 0.42) else 8
  ggsave(file.path(odir, sprintf("%s_region_boxplot_%s.png", TRAIT, scale_name)),
         g, width = w, height = 6, dpi = 150, limitsize = FALSE)
  cat(sprintf("  [C] %s | %s (%d groups, %s)\n", TRAIT, scale_name, ng, method_lab))
}

# =====================================================================
# Trait-level spatial autocorrelation check (Moran's I on OLS residuals)
#   Mirrors the PC1/PC2 framework: same inverse-distance weights, same
#   Moran.I call, same p<0.05 GLS-upgrade rule. One row per trait so you
#   can report which trait models needed spatial correction.
#
#   Assumes (as in the PC script):
#     - d_env : env-complete subsample with coords + predictors
#     - kept  : retained predictor names ; zp <- paste0(kept,".z") already built
#     - rhs   <- paste(zp, collapse = " + ")
#     - TRAITS, LOG, LATCOL, LONCOL defined
#   Requires: library(ape) for Moran.I ; library(nlme) for the GLS upgrade.
# =====================================================================

# ---- shared inverse-distance weight matrix (identical to sec 3) ------------
coords <- as.matrix(d_env[, c(LONCOL, LATCOL)])
dm <- as.matrix(dist(coords)); diag(dm) <- 0
W  <- 1 / dm; diag(W) <- 0; W[!is.finite(W)] <- 0

# ---- response transform: log the same traits the PCA logged --------------
trait_response <- function(tr) {
  x <- d_env[[tr]]
  if (tr %in% LOG) log(x) else x
}

# ---- per-trait OLS + Moran's I on residuals ------------------------------
trait_moran <- function(tr) {
  d_env$.y <- trait_response(tr)
  m  <- lm(as.formula(paste(".y ~", rhs)), data = d_env)
  mi <- tryCatch(ape::Moran.I(residuals(m), W),
                 error = function(e) list(observed=NA, expected=NA,
                                          sd=NA, p.value=NA))
  data.frame(
    trait     = tr,
    logged    = tr %in% LOG,
    moran_I   = round(mi$observed, 4),
    moran_exp = round(mi$expected, 4),
    moran_p   = signif(mi$p.value, 3),
    autocor   = !is.na(mi$p.value) && mi$p.value < 0.05,
    model     = ifelse(!is.na(mi$p.value) && mi$p.value < 0.05,
                       "GLS (spatially corrected)", "OLS"),
    n         = length(residuals(m)),
    row.names = NULL, stringsAsFactors = FALSE)
}

trait_moran_tab <- do.call(rbind, lapply(TRAITS, trait_moran))

cat("\n=== Trait-level spatial autocorrelation (Moran's I on OLS residuals) ===\n")
print(trait_moran_tab, row.names = FALSE)
cat(sprintf("\n%d of %d trait models show residual autocorrelation (p<0.05) -> GLS reported.\n",
            sum(trait_moran_tab$autocor), nrow(trait_moran_tab)))
write.csv(trait_moran_tab,
          file.path(OUT, "trait_moran_autocor.csv"), row.names = FALSE)

# ---- OPTIONAL: GLS upgrade for the flagged traits ------------------------
# Same corExp(nugget=TRUE) structure as the PC models. Uncomment to fit and
# export spatially-corrected trait coefficients for the autocorrelated traits.
#
# set.seed(42)
# d_env$.jit_lon <- d_env[[LONCOL]] + rnorm(nrow(d_env), 0, 1e-4)
# d_env$.jit_lat <- d_env[[LATCOL]] + rnorm(nrow(d_env), 0, 1e-4)
# gls_trait <- function(tr) {
#   d_env$.y <- trait_response(tr)
#   g <- tryCatch(
#     nlme::gls(as.formula(paste(".y ~", rhs)), data = d_env,
#               correlation = nlme::corExp(form = ~ .jit_lon + .jit_lat, nugget = TRUE),
#               method = "REML"),
#     error = function(e) { cat("  GLS failed for", tr, ":", conditionMessage(e), "\n"); NULL })
#   if (is.null(g)) return(NULL)
#   tt <- summary(g)$tTable
#   data.frame(trait = tr, term = rownames(tt),
#              estimate = tt[,"Value"], se = tt[,"Std.Error"],
#              t = tt[,"t-value"], p = tt[,"p-value"],
#              row.names = NULL, stringsAsFactors = FALSE)
# }
# flagged <- trait_moran_tab$trait[trait_moran_tab$autocor]
# gls_trait_tab <- do.call(rbind, lapply(flagged, gls_trait))
# write.csv(gls_trait_tab, file.path(OUT, "trait_gls_coefficients.csv"), row.names = FALSE)

# =====================================================================
# DRIVER LOOP: run the full pipeline for every trait
# =====================================================================
for (TRAIT in TRAITS) {
  cat(sprintf("\n########## TRAIT: %s (%s) ##########\n", TRAIT, trait_lab(TRAIT)))
  odir <- file.path(OUT, TRAIT)
  dir.create(odir, showWarnings = FALSE, recursive = TRUE)

  trait_correlations(TRAIT, odir)                 # (A) raw correlations
  trait_regression(TRAIT, odir)                   # (B) multiple regression
  for (sc in names(meow_cols_map))                # (C) regional comparison
    region_test_plot(TRAIT, sc, meow_cols_map[[sc]], odir)
}

cat("\n[done] Trait_Analysis pipeline complete. Outputs in", OUT, "\n")
cat("       one subfolder per trait, each with correlations + regression + regional figures.\n")
