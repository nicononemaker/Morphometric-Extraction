library(dplyr); library(readr)

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
TA_DIR <- file.path(OUT_DIR, "Trait_Analysis")   # read per-trait tables from here
SUPP_OUT <- TA_DIR                                # and write the supp tables here

PCA_TRAITS <- c("length_frond_mm","max_width_mm","stipe_to_length_ratio",
                "solidity","peak_position_frac","branched_fraction","circularity")

TRAIT_LABEL <- c(length_frond_mm="Frond length", max_width_mm="Max width",
                 stipe_to_length_ratio="Stipe width:leng", solidity="Solidity",
                 peak_position_frac="Peak position", branched_fraction="Branched fraction",
                 circularity="Circularity")
PRED_LABEL <- c(SST_clim="SST", log_nitrate_clim="Nitrate",
                log_phosphate_clim="Phosphate", salinity_clim="Salinity",
                log_wave.height="Wave height", log_kd490_clim="Turbidity",
                PAR_mean="PAR", log_iron_clim="Iron", spco2_clim="pCO\u2082")

collect <- function(kind) {
  rows <- lapply(PCA_TRAITS, function(tr) {
    f <- file.path(TA_DIR, tr, sprintf("%s_regression_%s.csv", tr, kind))
    if (!file.exists(f)) { message("missing: ", f); return(NULL) }
    read_csv(f, show_col_types = FALSE)
  })
  out <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
  out <- out[out$predictor != "(Intercept)", ]
  out$Trait     <- TRAIT_LABEL[out$trait]
  out$Predictor <- ifelse(out$predictor %in% names(PRED_LABEL),
                          PRED_LABEL[out$predictor], out$predictor)
  keep <- intersect(c("Trait","Predictor","estimate","se","t","p","sig","vif",
                      "model_R2","model_adjR2","moran_p","n"), names(out))
  out <- out[, keep]
  nn <- intersect(c("estimate","se","t","vif","model_R2","model_adjR2"), names(out))
  out[nn] <- lapply(out[nn], round, 3)
  out$p <- signif(out$p, 3)
  out
}

ols_tab <- collect("OLS")
gls_tab <- collect("GLS")

write.csv(ols_tab, file.path(SUPP_OUT, "Trait_Supp_OLS.csv"), row.names = FALSE)
write.csv(gls_tab, file.path(SUPP_OUT, "Trait_Supp_GLS.csv"), row.names = FALSE)
cat("wrote Trait_Supp_OLS.csv and Trait_Supp_GLS.csv to", SUPP_OUT, "\n")
cat("OLS rows:", nrow(ols_tab), " GLS rows:", nrow(gls_tab), "\n")