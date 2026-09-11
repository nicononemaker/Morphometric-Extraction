# =============================================================================
# 07_sort_outputs.R  —  organize analysis_out/ into manuscript vs exploratory,
# with themed subfolders. Run AFTER 02-06.
#
#   analysis_out/
#   ├── manuscript_results/
#   │   ├── manuscript_figures/       (paper Figs 1-5, S1-S7)
#   │   └── manuscript_tables/        (paper Table 1, S2-S9 + summary)
#   ├── exploratory_figures/
#   │   ├── regional/                 (region boxplots, paired, climate-zone)
#   │   ├── latitude/
#   │   ├── permanova/                (permanova/varpart/ordination figures)
#   │   ├── distance_decay/
#   │   └── diagnostics/              (scatter grids, coef, bubble, corr, scree)
#   └── exploratory_tables/
#       ├── regional/                 (region cld/pairwise, thermal)
#       ├── latitude/
#       ├── permanova/                (permanova/varpart tables)
#       └── diagnostics/              (scatter diagnostics, moran, correlations)
#
# Only files written by the CURRENT 02-06 pipeline are matched and moved. Stale
# CSVs from abandoned analyses (01_*, morph_*, omega_*, wr_*, stipe_*, etc.) are
# NOT matched and are LEFT IN PLACE, listed at the end. Safe to re-run.
# =============================================================================
if (file.exists("config.R")) source("config.R") else
if (file.exists("../config.R")) source("../config.R") else
  stop("config.R not found — set your working directory to the repo root.")

ROOTS <- c(OUT_DIR, file.path(OUT_DIR, "Trait_Analysis"))

# ---- destination folders ----------------------------------------------------
MR   <- file.path(OUT_DIR, "manuscript_results")
FIGm <- file.path(MR, "manuscript_figures")
TABm <- file.path(MR, "manuscript_tables")
EF   <- file.path(OUT_DIR, "exploratory_figures")
ET   <- file.path(OUT_DIR, "exploratory_tables")
subs <- c("regional","latitude","permanova","distance_decay","diagnostics")
dirs <- c(FIGm, TABm,
          file.path(EF, subs),
          file.path(ET, c("regional","latitude","permanova","diagnostics")))
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- routing rules: (regex, destination). First match wins. ----------------
# order matters: manuscript patterns first, then themed exploratory.
RULES <- list(
  # ---- MANUSCRIPT FIGURES ----
  list("^pca_biplot_linear\\.png$",                         FIGm),
  list("^pca_scree\\.png$",                                 FIGm),
  list("^region_boxplot_ecoregion_PC1_PC2_paired\\.png$",   FIGm),
  list("^region_boxplot_(realm|province)_PC1_PC2_paired\\.png$", FIGm),
  list("^climate_zone_ecoregion_boxplot_PC[12]\\.png$",     FIGm),
  list("^scatter_PC[12]_vs_env_linear\\.png$",              FIGm),
  list("^env_predictor_corr_heatmap\\.png$",                FIGm),
  list("^permanova_pairwise_(realm|province|ecoregion)\\.png$", FIGm),
  list("^varpart_venn_(realm|province|ecoregion)\\.png$",   FIGm),
  # ---- MANUSCRIPT TABLES ----
  list("^publication_.*\\.(csv|txt)$",                      TABm),
  list("^pc_loadings(_all)?\\.csv$",                        TABm),
  list("^pc_scree_variance\\.csv$",                         TABm),
  list("^env_predictor_corr_matrix\\.csv$",                 TABm),
  list("^pc_env_regression.*\\.csv$",                       TABm),
  list("^TABLE_(env_permanova|varpart_stats|env_R2_by_frame)\\.csv$", TABm),
  list("^Trait_Supp_(OLS|GLS)\\.csv$",                      TABm),
  # ---- EXPLORATORY FIGURES ----
  list("^region_boxplot_.*_cld\\.png$",             file.path(EF,"regional")),
  list("^climate_zone_.*\\.png$",                   file.path(EF,"regional")),
  list("^latitude_gradient.*\\.png$",               file.path(EF,"latitude")),
  list("^distdecay_.*\\.png$",                      file.path(EF,"distance_decay")),
  list("^(permanova_|varpart_|nested_ordination|residual_ordination|r2_by_scale|env_permanova|env_sequential|env_stack).*\\.png$",
                                                    file.path(EF,"permanova")),
  list("^TABLE_env_permanova\\.png$",               file.path(EF,"permanova")),
  list("^(scatter_.*_linear.*|coef_plot|spearman_bubble|scatter_realm_legend|.*_coef_plot|.*_scatter.*)\\.png$",
                                                    file.path(EF,"diagnostics")),
  # ---- EXPLORATORY TABLES ----
  list("^region_.*_(cld|pairwise|pairsig)\\.csv$",  file.path(ET,"regional")),
  list("^thermal_(class_effects|region_labels)\\.csv$", file.path(ET,"regional")),
  list("^climate_zone_.*_cld\\.csv$",               file.path(ET,"regional")),
  list("^pairwise_zone_summary_.*\\.csv$",          file.path(ET,"regional")),
  list("^latitude_gradient_stats\\.csv$",           file.path(ET,"latitude")),
  list("^(env_permanova_marginal|varpart_.*|permanova_.*)\\.csv$", file.path(ET,"permanova")),
  list("^(scatter_.*_diagnostics|trait_moran_autocor|trait_gls_coefficients|.*_correlations|.*_regression_(OLS|GLS)|.*_regression_.*)\\.csv$",
                                                    file.path(ET,"diagnostics"))
)

route <- function(fname) {
  for (r in RULES) if (grepl(r[[1]], fname)) return(r[[2]])
  NA_character_   # unmatched -> leave in place (stale)
}

moved <- 0; left <- character(0)
for (root in ROOTS) {
  if (!dir.exists(root)) next
  files <- list.files(root, pattern = "\\.(png|csv|txt)$", full.names = FALSE)
  for (f in files) {
    dest <- route(f)
    if (is.na(dest)) { left <- c(left, f); next }
    if (file.rename(file.path(root, f), file.path(dest, f))) moved <- moved + 1
    else left <- c(left, f)
  }
}

cat(sprintf("Sorted %d files.\n  manuscript: %s\n  exploratory figs: %s\n  exploratory tables: %s\n",
            moved, MR, EF, ET))
if (length(left)) {
  cat(sprintf("\nLeft in place (%d — stale/unmatched, not from current pipeline):\n", length(left)))
  cat(paste0("  ", left, collapse = "\n"), "\n")
}
