# =============================================================================
# PERMANOVA MASTER  —  PART 2 of 2 : FIGURES (fast; no permutations)
#
# Reads the cache written by permanova_1_analysis.R and draws every figure.
# Re-run this as often as you like to iterate on aesthetics — it never permutes.
#
# Figures per region scale:
#   A  pairwise trait-divergence heatmap
#   B  trait-space PCoA centroids            (realm palette)
#   C  geographic centroid map               (realm palette)
#   C2 geographic EDGE map: region-to-region trait difference
#        edge WIDTH = centroid distance D[C] (§5.1)
#        edge COLOR = pairwise significance (Holm)
#        edges filtered to geographically proximate pairs (Delaunay) by default
#   E  trait-space mMDS "main-effects" inset (§5.2)  — separate figure
# Plus:
#   env marginal-R2 bar
#   varpart component bar + absorption bar   (per scale)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
})
have_maps      <- requireNamespace("maps", quietly = TRUE)
have_deldir    <- requireNamespace("deldir", quietly = TRUE)   # for Delaunay edge filter
have_newscale  <- requireNamespace("ggnewscale", quietly = TRUE) # dual color legends on edge map

# ------------------------------------------------------------ config / cache --
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
OUT   <- OUT_DIR
CACHE <- file.path(OUT, "cache")
stopifnot(dir.exists(CACHE))
say <- function(...) cat(..., "\n")

meta <- readRDS(file.path(CACHE, "meta.rds"))
REGIONS <- meta$REGIONS

# --- edge map design (Figure C2) ---------------------------------------------
# The edge map now always draws: (1) a BACKBONE of geographically adjacent realm
# pairs (Delaunay), straight lines, width = D[C], colored red (sig) / grey (n.s.);
# plus (2) curved blue arcs for ALL realm pairs that are NOT significantly
# different — the "surprising" convergences. EDGE_MODE below is retained only for
# reference and no longer gates the C2 layers.
EDGE_MODE <- "delaunay_sig"

realm_cols <- c(
  "Atlantic Warm Water"          = "red2",
  "Indo-Pacific Warm Water"      = "darkorange",
  "Tropical Atlantic"            = "gold",
  "Tropical Eastern Pacific"     = "pink",
  "Western Indo-Pacific"         = "deeppink",
  "Central Indo-Pacific"         = "sienna4",
  "Eastern Indo-Pacific"         = "lightcoral",
  "Temperate Northern Atlantic"  = "powderblue",
  "Temperate Northern Pacific"   = "midnightblue",
  "Temperate South America"      = "slateblue",
  "Temperate Australasia"        = "royalblue",
  "Southern Ocean"               = "darkcyan"
)
scale_for <- function(GROUP, present_levels, aes = "colour") {
  if (toupper(GROUP) == "REALM" && all(present_levels %in% names(realm_cols))) {
    vals <- realm_cols[present_levels]
    if (aes == "fill") scale_fill_manual(values = vals, name = NULL)
    else               scale_colour_manual(values = vals, name = NULL)
  } else {
    if (aes == "fill") scale_fill_discrete(name = NULL)
    else               scale_colour_discrete(name = NULL)
  }
}

# clean display name for a scale label (for titles/subtitles/legends)
nice_scale <- function(GROUP) {
  m <- c(REALM = "Realm", PROVINC = "Province", PROVINCE = "Province",
         ECOREGION = "Ecoregion")
  out <- m[toupper(GROUP)]
  if (is.na(out)) return(tools::toTitleCase(tolower(GROUP)))
  unname(out)
}

# Delaunay neighbor pairs among centroid coords (returns char matrix of grp pairs)
delaunay_pairs <- function(coords) {  # coords: data.frame(grp, x, y)
  if (!have_deldir || nrow(coords) < 3) return(NULL)
  dl <- deldir::deldir(coords$x, coords$y, suppressMsge = TRUE)
  seg <- dl$delsgs  # ind1, ind2 index into coords
  data.frame(r1 = coords$grp[seg$ind1], r2 = coords$grp[seg$ind2],
             stringsAsFactors = FALSE)
}

# ---- two-circle vegan-style Venn (region vs environment) --------------------
# left-only = pure region [R|E]; overlap = shared; right-only = pure env [E|R];
# a labeled "Residuals" note sits below. Circles are equal-size overlapping
# disks (classic vegan varpart look), NOT area-proportional.
venn_region_env <- function(pure_region, shared, pure_env, residual,
                            title, subtitle,
                            left_fill = "#B2182B", mid_fill = "#4393C3",
                            right_fill = "#92C5DE") {
  th <- seq(0, 2*pi, length.out = 200)
  r  <- 1.15; dx <- 0.72                     # radius and center offset
  circ <- function(cx) data.frame(x = cx + r*cos(th), y = r*sin(th))
  L <- circ(-dx); R <- circ(dx)
  fmt <- function(v) if (is.finite(v) && v > 0.0005) sprintf("%.3f", v) else "~0"

  ggplot() +
    geom_polygon(data = L, aes(x, y), fill = left_fill,  alpha = .55, colour = "grey30") +
    geom_polygon(data = R, aes(x, y), fill = right_fill, alpha = .55, colour = "grey30") +
    # value labels: pure values inside each circle's non-overlapping crescent,
    # shared value in the central overlap
    annotate("text", x = -dx-0.42, y = 0.12, label = "Pure\nregion", fontface = "bold", size = 3.1) +
    annotate("text", x = -dx-0.42, y = -0.28, label = fmt(pure_region), size = 3.4) +
    annotate("text", x = 0,        y = 0.12, label = "Shared", fontface = "bold", size = 3.1) +
    annotate("text", x = 0,        y = -0.28, label = fmt(shared), size = 3.4) +
    annotate("text", x =  dx+0.42, y = 0.12, label = "Pure\nenv.", fontface = "bold", size = 3.1) +
    annotate("text", x =  dx+0.42, y = -0.28, label = fmt(pure_env), size = 3.4) +
    # circle identity labels above
    annotate("text", x = -dx, y = r+0.22, label = "Region", fontface = "bold", size = 3.6, colour = left_fill) +
    annotate("text", x =  dx, y = r+0.22, label = "Environment", fontface = "bold", size = 3.6, colour = "#2166AC") +
    annotate("text", x = 0, y = -r-0.35,
             label = sprintf("Residuals = %s", fmt(residual)), size = 3, colour = "grey40") +
    coord_equal(clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"),
          plot.margin = margin(10, 30, 10, 30))
}

# =============================================================================
# per-scale region figures
# =============================================================================
draw_region <- function(GROUP) {
  tag <- tolower(GROUP)
  fp  <- file.path(CACHE, sprintf("region_%s.rds", tag))
  if (!file.exists(fp)) { say("  no cache for", GROUP); return(invisible()) }
  R <- readRDS(fp)
  lv <- R$levels
  say(sprintf("\n## FIGURES  %s ##", GROUP))

  # -------- FIGURE A  pairwise R2 heatmap --------
  pw <- R$pairwise

  # ---- classify + order realms by mean absolute latitude ----
  # tropical if mean |lat| < 23.5 (Tropic lines); order within group by |lat|.
  # mean |lat| computed from cached per-specimen coords when available.
  abslat <- setNames(rep(NA_real_, length(lv)), lv)
  if (!is.null(R$specimens_xy)) {
    sx <- R$specimens_xy
    ml <- tapply(abs(sx$.lat), sx$.grp, mean, na.rm = TRUE)
    abslat[names(ml)[names(ml) %in% lv]] <- ml[names(ml) %in% lv]
  }
  # fallback: mean of signed geo lat if specimen coords missing
  if (all(is.na(abslat)) && !is.null(R$geo)) {
    g2 <- R$geo; abslat[as.character(g2$.grp)] <- abs(g2$lat)
  }
  zone <- ifelse(is.na(abslat), "Unclassified",
                 ifelse(abslat < 23.5, "Tropical", "Temperate"))
  # order by ascending |latitude| overall: tropical realms (low |lat|) come first,
  # temperate (high |lat|) second. With a non-reversed y-axis this puts
  # tropical x tropical in the bottom-left (smallest |lat|) and temperate x
  # temperate top-right (largest |lat|), so |lat| increases monotonically along
  # both axes and the zone grouping agrees with the ordering.
  zlev <- c("Tropical", "Temperate", "Unclassified")
  ord_df <- data.frame(realm = lv, abslat = abslat, zone = factor(zone, levels = zlev),
                       stringsAsFactors = FALSE)
  ord_df <- ord_df[order(ord_df$zone, ord_df$abslat), ]
  ord <- ord_df$realm                        # new axis order

  # -------- within/between-zone pairwise divergence summary --------
  # Uses the SAME `pw` table and `zone` classification as the heatmap below,
  # so the reported numbers match figure A exactly. Buckets each region pair as
  # within-tropical, within-temperate, or between-zone, and reports how many
  # were significant (Holm P_adj < 0.05) and the mean/median R2 per bucket.
  {
    zvec <- setNames(zone, lv)            # region -> zone, same rule as heatmap
    pw_z <- pw
    pw_z$z1 <- zvec[as.character(pw_z$r1)]
    pw_z$z2 <- zvec[as.character(pw_z$r2)]

    # drop pairs touching an Unclassified region (they aren't in either zone)
    keep <- pw_z$z1 %in% c("Tropical","Temperate") &
            pw_z$z2 %in% c("Tropical","Temperate")
    pw_z <- pw_z[keep, ]

    if (nrow(pw_z)) {
      pw_z$bucket <- ifelse(pw_z$z1 != pw_z$z2, "Between-zone",
                     ifelse(pw_z$z1 == "Tropical", "Within-tropical",
                                                   "Within-temperate"))

      # significance flag: prefer Holm P_adj (what the figure uses); fall back to `sig`
      sig_flag <- if ("P_adj" %in% names(pw_z)) pw_z$P_adj < 0.05 else pw_z$sig %in% TRUE

      zsum <- do.call(rbind, lapply(
        split(seq_len(nrow(pw_z)), pw_z$bucket), function(ix) {
          data.frame(
            scale     = nice_scale(GROUP),
            bucket    = pw_z$bucket[ix[1]],
            n_pairs   = length(ix),
            n_sig     = sum(sig_flag[ix], na.rm = TRUE),
            pct_sig   = round(100 * sum(sig_flag[ix], na.rm = TRUE) / length(ix), 1),
            mean_R2   = round(mean(pw_z$R2[ix], na.rm = TRUE), 3),
            median_R2 = round(median(pw_z$R2[ix], na.rm = TRUE), 3),
            stringsAsFactors = FALSE)
        }))

      # overall within (trop+temp combined) vs between
      within_ix  <- which(pw_z$bucket %in% c("Within-tropical","Within-temperate"))
      between_ix <- which(pw_z$bucket == "Between-zone")
      zsum2 <- data.frame(
        scale   = nice_scale(GROUP),
        bucket  = c("Within (all)", "Between-zone (all)"),
        n_pairs = c(length(within_ix), length(between_ix)),
        n_sig   = c(sum(sig_flag[within_ix], na.rm=TRUE), sum(sig_flag[between_ix], na.rm=TRUE)),
        pct_sig = c(round(100*sum(sig_flag[within_ix],na.rm=TRUE)/max(length(within_ix),1),1),
                    round(100*sum(sig_flag[between_ix],na.rm=TRUE)/max(length(between_ix),1),1)),
        mean_R2 = c(round(mean(pw_z$R2[within_ix], na.rm=TRUE),3),
                    round(mean(pw_z$R2[between_ix], na.rm=TRUE),3)),
        median_R2 = c(round(median(pw_z$R2[within_ix], na.rm=TRUE),3),
                      round(median(pw_z$R2[between_ix], na.rm=TRUE),3)),
        stringsAsFactors = FALSE)

      zsum_all <- rbind(zsum, zsum2)
      cat(sprintf("\n--- Pairwise divergence by zone bucket (%s) ---\n", nice_scale(GROUP)))
      print(zsum_all, row.names = FALSE)
      write.csv(zsum_all,
                file.path(OUT, sprintf("pairwise_zone_summary_%s.csv", tag)),
                row.names = FALSE)
    } else {
      cat(sprintf("\n[zone summary] no classifiable pairs at %s scale; skipped\n",
                  nice_scale(GROUP)))
    }
  }

  # ---- long symmetric table with the new ordering ----
  sym <- bind_rows(pw[, c("r1","r2","R2","sig")],
                   setNames(pw[, c("r2","r1","R2","sig")], c("r1","r2","R2","sig")))
  # x = ord (temperate block first => leftmost). y = ord too (NOT reversed), so the
  # bottom-left corner is (temperate, temperate) = within-temperate comparisons.
  sym$r1 <- factor(sym$r1, levels = ord); sym$r2 <- factor(sym$r2, levels = ord)
  # fill value only for significant cells; n.s. -> NA (drawn white)
  sym$R2_sig <- ifelse(sym$sig %in% TRUE, sym$R2, NA_real_)
  sym$lab    <- ifelse(sym$sig %in% TRUE, sprintf("%.2f", sym$R2), "n.s.")

  # ---- zone side-strips: one tile-row/col strip flagging tropical/temperate ----
  zone_cols <- c(Tropical = "firebrick2", Temperate = "cornflowerblue", Unclassified = "grey70")
  strip_x <- data.frame(r = factor(ord, levels = ord),
                        zone = factor(ord_df$zone, levels = zlev))
  strip_y <- data.frame(r = factor(ord, levels = ord),
                        zone = factor(ord_df$zone[match(ord, ord_df$realm)], levels = zlev))
  # group label positions (midpoint of each zone block along the axis)
  zpos <- do.call(rbind, lapply(split(seq_along(ord), ord_df$zone[match(ord, ord_df$realm)]),
                   function(idx) if (length(idx)) data.frame(mid = mean(idx),
                                    lo = min(idx), hi = max(idx)) else NULL))
  zpos$zone <- rownames(zpos)

  sym$txt_col <- ifelse(sym$sig %in% TRUE, "grey10", "grey60")
  fA <- ggplot(sym, aes(r1, r2)) +
    geom_tile(aes(fill = R2_sig), colour = "grey55") +
    geom_text(aes(label = lab, fontface = ifelse(sig %in% TRUE, "bold", "plain"),
                  colour = txt_col), size = 2.7, show.legend = FALSE) +
    scale_colour_identity() +
    scale_fill_gradient(low = "#E5F5E0", high = "#00441B", na.value = "white",
                        name = expression(R^2), trans = "sqrt")

  # side-strips: a thin colored band beneath the x-axis and left of the y-axis.
  # placed just outside the tile grid using continuous coordinates on discrete axes.
  nL <- length(ord)
  if (have_newscale) {
    fA <- fA + ggnewscale::new_scale_fill() +
      # bottom strip (x zones), placed just below row 1
      geom_tile(data = strip_x, aes(x = r, y = 0.3, fill = zone),
                height = 0.4, inherit.aes = FALSE) +
      # left strip (y zones), placed just left of column 1
      geom_tile(data = strip_y, aes(x = 0.3, y = r, fill = zone),
                width = 0.4, inherit.aes = FALSE) +
      scale_fill_manual(values = zone_cols, name = "Zone", drop = TRUE)
  }
  # divider lines splitting the matrix into trop×trop / temp×temp / cross blocks.
  # tropical-first ordering => boundary sits after the tropical block on both axes.
  n_trop <- sum(ord_df$zone[match(ord, ord_df$realm)] == "Tropical")
  if (n_trop > 0 && n_trop < nL) {
    xb <- n_trop + 0.5                 # vertical divider (tropical | temperate on x)
    yb <- n_trop + 0.5                 # horizontal divider (tropical | temperate on y)
    fA <- fA +
      geom_vline(xintercept = xb, colour = "grey20", linewidth = 0.8) +
      geom_hline(yintercept = yb, colour = "grey20", linewidth = 0.8)
  }
  fA <- fA +
    coord_cartesian(clip = "off", xlim = c(0.1, nL + 0.5),
                    ylim = c(0.1, nL + 0.5)) +
    labs(title = sprintf("Pairwise trait divergence between %ss", tolower(nice_scale(GROUP))),
         subtitle = sprintf("PERMANOVA R\u00b2; %d of %d pairs significant (Holm p<0.05); white = n.s.; grouped tropical/temperate, ordered by |latitude|",
                            sum(pw$sig), nrow(pw)), x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
          plot.margin = margin(20, 12, 8, 12))
  sz <- max(7.5, 0.6 * length(lv))
  ggsave(file.path(OUT, sprintf("permanova_pairwise_%s.png", tag)), fA,
         width = sz + 1, height = sz, dpi = 140, bg = "white")

  # ============================================================================
  # FIGURE A2 / A3  p-value heatmaps (reuse ordering, strips, dividers from A)
  # ============================================================================
  # need the adjusted p in a symmetric table; P_adj is what `sig` was based on.
  if ("P_adj" %in% names(pw)) {
    symp <- bind_rows(pw[, c("r1","r2","P_adj","sig")],
                      setNames(pw[, c("r2","r1","P_adj","sig")], c("r1","r2","P_adj","sig")))
    symp$r1 <- factor(symp$r1, levels = ord); symp$r2 <- factor(symp$r2, levels = ord)

    # significance TIER (binned), only for significant cells; n.s. -> NA (white)
    tier <- cut(symp$P_adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                labels = c("p<0.001", "p<0.01", "p<0.05", "n.s."), right = TRUE)
    symp$tier <- ifelse(symp$sig %in% TRUE, as.character(tier), NA_character_)
    symp$tier <- factor(symp$tier, levels = c("p<0.05", "p<0.01", "p<0.001"))
    # cell label: the p-value if significant, else "n.s."
    symp$plab <- ifelse(symp$sig %in% TRUE,
                        ifelse(symp$P_adj < 0.001, "<0.001", sprintf("%.3f", symp$P_adj)),
                        "n.s.")
    symp$txt_col <- ifelse(symp$sig %in% TRUE, "grey10", "grey60")
    # binary green label
    symp$binlab <- ifelse(symp$sig %in% TRUE, "sig", "n.s.")
    symp$binfill <- ifelse(symp$sig %in% TRUE, "sig", NA_character_)

    # helper to append the shared strips + dividers + theme to a base heatmap
    add_frame <- function(g, ttl, sub) {
      if (have_newscale) {
        g <- g + ggnewscale::new_scale_fill() +
          geom_tile(data = strip_x, aes(x = r, y = 0.3, fill = zone),
                    height = 0.4, inherit.aes = FALSE) +
          geom_tile(data = strip_y, aes(x = 0.3, y = r, fill = zone),
                    width = 0.4, inherit.aes = FALSE) +
          scale_fill_manual(values = zone_cols, name = "Zone", drop = TRUE)
      }
      if (n_trop > 0 && n_trop < nL)
        g <- g + geom_vline(xintercept = n_trop + 0.5, colour = "grey20", linewidth = 0.8) +
                 geom_hline(yintercept = n_trop + 0.5, colour = "grey20", linewidth = 0.8)
      g + coord_cartesian(clip = "off", xlim = c(0.1, nL + 0.5), ylim = c(0.1, nL + 0.5)) +
        labs(title = ttl, subtitle = sub, x = NULL, y = NULL) +
        theme_bw(base_size = 9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              panel.grid = element_blank(), plot.margin = margin(20, 12, 8, 12))
    }

    # ---- FIGURE A2: tiered-green by significance level (3 discrete shades) ----
    tier_cols <- c("p<0.05" = "#C7E9C0", "p<0.01" = "#74C476", "p<0.001" = "#238B45")
    fA2 <- ggplot(symp, aes(r1, r2)) +
      geom_tile(aes(fill = tier), colour = "grey55") +
      geom_text(aes(label = plab, colour = txt_col,
                    fontface = ifelse(sig %in% TRUE, "bold", "plain")),
                size = 2.5, show.legend = FALSE) +
      scale_colour_identity() +
      scale_fill_manual(values = tier_cols, na.value = "white",
                        name = "Significance", drop = FALSE)
    fA2 <- add_frame(fA2,
      sprintf("Pairwise significance between %ss", tolower(nice_scale(GROUP))),
      sprintf("Holm-adjusted p; darker green = smaller p (more significant); white = n.s. (%d of %d sig)",
              sum(pw$sig), nrow(pw)))
    ggsave(file.path(OUT, sprintf("permanova_pairwise_pval_%s.png", tag)), fA2,
           width = sz + 1, height = sz, dpi = 140, bg = "white")

    # ---- FIGURE A3: binary green (significant vs n.s.), single shade ----
    # single green for all significant cells, but print the (adjusted) p-value.
    fA3 <- ggplot(symp, aes(r1, r2)) +
      geom_tile(aes(fill = binfill), colour = "grey55") +
      geom_text(aes(label = plab, colour = txt_col,
                    fontface = ifelse(sig %in% TRUE, "bold", "plain")),
                size = 2.5, show.legend = FALSE) +
      scale_colour_identity() +
      scale_fill_manual(values = c("sig" = "#41AB5D"), na.value = "white",
                        name = NULL, labels = "significant", drop = FALSE)
    fA3 <- add_frame(fA3,
      sprintf("Pairwise significant differences between %ss", tolower(nice_scale(GROUP))),
      sprintf("green = significant (Holm p<0.05, value shown); white = n.s. (%d of %d sig)",
              sum(pw$sig), nrow(pw)))
    ggsave(file.path(OUT, sprintf("permanova_pairwise_binary_%s.png", tag)), fA3,
           width = sz + 1, height = sz, dpi = 140, bg = "white")
    say(sprintf("  wrote p-value + binary heatmaps for %s", tag))
  }

  # -------- FIGURE B  trait-space PCoA centroids --------
  cen <- R$centroids_trait; pts <- R$specimens_trait; ve <- R$var_expl
  fB <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_point(data = pts, aes(Ax1, Ax2, colour = grp), size = .7, alpha = .2) +
    geom_point(data = cen, aes(Ax1, Ax2, colour = grp, size = n)) +
    geom_text_repel(data = cen, aes(Ax1, Ax2, label = grp, colour = grp),
                    size = 3, fontface = "bold", max.overlaps = 30, show.legend = FALSE) +
    scale_size(range = c(3, 9), name = "n") +
    scale_for(GROUP, lv, "colour") +
    labs(title = sprintf("%s centroids in trait space (PCoA)", nice_scale(GROUP)),
         subtitle = "large points = centroids; faded points = specimens",
         x = sprintf("PCoA 1 (%.1f%%)", ve[1]), y = sprintf("PCoA 2 (%.1f%%)", ve[2])) +
    guides(colour = "none") + theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("permanova_centroids_%s.png", tag)), fB,
         width = 9, height = 7, dpi = 140, bg = "white")

  # -------- FIGURE C  geographic centroid map --------
  if (!is.null(R$geo) && have_maps) {
    geo <- R$geo; sx <- R$specimens_xy
    world <- ggplot2::map_data("world")
    fC <- ggplot() +
      geom_polygon(data = world, aes(long, lat, group = group),
                   fill = "grey92", colour = "grey80", linewidth = 0.2) +
      geom_point(data = sx, aes(.lon, .lat, colour = .grp), size = .8, alpha = .25) +
      geom_point(data = geo, aes(lon, lat, colour = .grp, size = n)) +
      geom_text_repel(data = geo, aes(lon, lat, label = .grp, colour = .grp),
                      size = 3, fontface = "bold", max.overlaps = 30, show.legend = FALSE) +
      scale_size(range = c(3, 10), name = "n") +
      scale_for(GROUP, lv, "colour") +
      coord_quickmap(xlim = range(geo$lon) + c(-15, 15),
                     ylim = range(geo$lat) + c(-15, 15)) +
      labs(title = sprintf("%s centroids in geographic space", nice_scale(GROUP)),
           subtitle = "large points = region centroids (sized by n); faded points = specimens",
           x = NULL, y = NULL) +
      guides(colour = "none") + theme_minimal(base_size = 11) +
      theme(panel.grid = element_line(colour = "grey95"))
    ggsave(file.path(OUT, sprintf("permanova_geomap_%s.png", tag)), fC,
           width = 11, height = 6.5, dpi = 140, bg = "white")

    # -------- FIGURE C2  geographic EDGE map (trait difference) --------
    DCm <- as.matrix(R$DC)
    # full pairwise edge table with centroid distance D[C] and coords
    coords <- geo %>% transmute(grp = as.character(.grp), x = lon, y = lat)
    all_edges <- pw %>% transmute(r1, r2, sig,
                                  DC = mapply(function(a, b) DCm[a, b], r1, r2)) %>%
      left_join(coords, by = c("r1" = "grp")) %>% rename(x1 = x, y1 = y) %>%
      left_join(coords, by = c("r2" = "grp")) %>% rename(x2 = x, y2 = y) %>%
      filter(is.finite(x1), is.finite(x2))

    key <- function(a, b) paste(pmin(a, b), pmax(a, b))

    # which pairs are geographically adjacent (Delaunay)?
    dp <- delaunay_pairs(coords)
    near <- if (!is.null(dp)) key(dp$r1, dp$r2) else character(0)
    if (is.null(dp)) say("  (deldir not installed; backbone falls back to all pairs)")
    all_edges$adjacent <- key(all_edges$r1, all_edges$r2) %in% near

    # ---- layer 1: BACKBONE = adjacent pairs, colored by significance ----
    backbone <- if (length(near)) all_edges %>% filter(adjacent) else all_edges
    backbone$sig_lab <- factor(ifelse(backbone$sig %in% TRUE, "sig (Holm p<0.05)", "n.s."),
                               levels = c("sig (Holm p<0.05)", "n.s."))

    # ---- layer 2: SURPRISING = ALL non-significant pairs, as curved arcs ----
    #   (n.s. = realms NOT significantly different; the exception worth showing)
    surprising <- all_edges %>% filter(!(sig %in% TRUE))
    n_sig <- sum(all_edges$sig %in% TRUE); n_ns <- nrow(surprising)

    fC2 <- ggplot() +
      geom_polygon(data = world, aes(long, lat, group = group),
                   fill = "grey93", colour = "grey82", linewidth = 0.2) +
      # backbone straight edges
      geom_segment(data = backbone,
                   aes(x = x1, y = y1, xend = x2, yend = y2,
                       linewidth = DC, colour = sig_lab)) +
      scale_linewidth(range = c(0.3, 3.2), name = "Centroid\ndistance D[C]") +
      scale_colour_manual(values = c("sig (Holm p<0.05)" = "#B2182B", "n.s." = "grey65"),
                          drop = FALSE, name = NULL)
    # curved arcs for non-significant pairs (only if any exist)
    if (nrow(surprising)) {
      if (have_newscale) {
        fC2 <- fC2 +
          ggnewscale::new_scale_colour() +
          geom_curve(data = surprising,
                     aes(x = x1, y = y1, xend = x2, yend = y2),
                     curvature = 0.25, linewidth = 0.9, colour = "#2166AC",
                     alpha = 0.9, lineend = "round") +
          geom_point(data = data.frame(x = NA, y = NA), aes(x, y, colour = "not sig. different")) +
          scale_colour_manual(values = c("not sig. different" = "#2166AC"), name = NULL,
                              guide = guide_legend(override.aes = list(shape = NA, linetype = 1)))
      } else {
        # no ggnewscale: draw dashed blue arcs, describe them in the subtitle instead
        say("  (ggnewscale not installed; n.s. arcs drawn dashed, no separate legend key)")
        fC2 <- fC2 +
          geom_curve(data = surprising,
                     aes(x = x1, y = y1, xend = x2, yend = y2),
                     curvature = 0.25, linewidth = 0.9, colour = "#2166AC",
                     linetype = "dashed", alpha = 0.9, lineend = "round")
      }
    }
    fC2 <- fC2 +
      geom_point(data = geo, aes(lon, lat), size = 2.4, colour = "grey20") +
      geom_text_repel(data = geo, aes(lon, lat, label = .grp),
                      size = 2.8, fontface = "bold", max.overlaps = 40, colour = "grey15") +
      coord_quickmap(xlim = range(geo$lon) + c(-15, 15),
                     ylim = range(geo$lat) + c(-15, 15)) +
      labs(title = sprintf("Region-to-region trait divergence in geographic space (%s)", nice_scale(GROUP)),
           subtitle = sprintf("straight edges = neighboring %ss (width = D[C]; red = sig, grey = n.s.);  blue arcs = ALL %s pairs NOT significantly different  (%d sig, %d n.s.)",
                              tolower(nice_scale(GROUP)), tolower(nice_scale(GROUP)), n_sig, n_ns),
           x = NULL, y = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_line(colour = "grey96"))
    ggsave(file.path(OUT, sprintf("permanova_edgemap_%s.png", tag)), fC2,
           width = 11, height = 6.5, dpi = 140, bg = "white")
    say(sprintf("  wrote geomap + edgemap for %s (%d backbone, %d n.s. arcs)",
                tag, nrow(backbone), nrow(surprising)))
  } else if (!is.null(R$geo) && !have_maps) {
    say("  install.packages('maps') to enable geographic maps.")
  }

  # -------- FIGURE E  trait-space mMDS main-effects inset (§5.2), separate --------
  if (!is.null(R$mds_pts)) {
    m <- R$mds_pts
    fE <- ggplot(m, aes(MDS1, MDS2, colour = grp)) +
      geom_hline(yintercept = 0, colour = "grey90") +
      geom_vline(xintercept = 0, colour = "grey90") +
      geom_point(size = 5) +
      geom_text_repel(aes(label = grp), size = 3.2, fontface = "bold",
                      max.overlaps = 40, show.legend = FALSE) +
      scale_for(GROUP, lv, "colour") +
      labs(title = sprintf("%s main-effects ordination (mMDS on D[C])", nice_scale(GROUP)),
           subtitle = sprintf("distances among region centroids in trait space; stress = %.3f",
                              R$mds_stress),
           x = "mMDS 1", y = "mMDS 2") +
      guides(colour = "none") + theme_bw(base_size = 11)
    ggsave(file.path(OUT, sprintf("permanova_mmds_maineffects_%s.png", tag)), fE,
           width = 7.5, height = 6.5, dpi = 140, bg = "white")
    say(sprintf("  wrote mMDS main-effects inset for %s", tag))
  }
  invisible(TRUE)
}
invisible(lapply(REGIONS, draw_region))

# =============================================================================
# environment marginal-R2 bar
# =============================================================================
env <- readRDS(file.path(CACHE, "env.rds"))
envbar <- env$tab %>% filter(term %in% env$kept) %>%
  mutate(sig = .data[["Pr(>F)"]] < 0.05)
# join bootstrap CIs if present (older caches won't have them)
if (!is.null(env$ci)) envbar <- envbar %>% left_join(env$ci, by = "term")
envbar <- envbar %>% mutate(term = reorder(term, R2))

fEnv <- ggplot(envbar, aes(term, R2, fill = sig)) +
  geom_col(width = .7)
# #5  bootstrap 95% CI whiskers on each marginal R2
if (!is.null(env$ci)) {
  fEnv <- fEnv +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .18, colour = "grey25", linewidth = .4)
}
fEnv <- fEnv +
  geom_text(aes(label = sprintf("%.3f%s", R2, ifelse(sig, "*", "")),
                y = if (!is.null(env$ci)) hi else R2), hjust = -0.2, size = 3.2) +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "grey70"),
                    labels = c("n.s.", "p<0.05"), name = NULL) +
  coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, .22))) +
  labs(title = "Environmental drivers of trait variation",
       # #3  report total model R2 and the shared (unattributable) fraction
       subtitle = if (!is.null(env$R2_total) && is.finite(env$R2_total)) {
         sprintf("Marginal PERMANOVA R\u00b2 (unique contribution); whiskers = 95%% bootstrap CI (%d reps); * = p<0.05\nTotal model R\u00b2 = %.3f  |  sum of marginal = %.3f  |  shared among collinear predictors = %.3f",
                 env$nboot, env$R2_total, env$R2_marg_sum, env$R2_shared)
       } else "Marginal PERMANOVA R\u00b2 (by = \"margin\"); * = p<0.05",
       x = NULL, y = expression("Marginal " * R^2)) +
  theme_bw(base_size = 11) + theme(panel.grid.major.y = element_blank())
ggsave(file.path(OUT, "env_permanova_marginal.png"), fEnv, width = 8.5, height = 4.8, dpi = 140, bg = "white")
say("\nwrote env_permanova_marginal.png")

# #3 (companion) small stacked bar: total explained = attributed + shared
if (!is.null(env$R2_total) && is.finite(env$R2_total)) {
  shar <- data.frame(
    part = factor(c("Uniquely attributed (sum of marginals)","Shared among collinear predictors"),
                  levels = c("Shared among collinear predictors","Uniquely attributed (sum of marginals)")),
    R2 = c(env$R2_marg_sum, env$R2_shared))
  fShare <- ggplot(shar, aes(x = "", y = R2, fill = part)) +
    geom_col(width = .5, colour = "white") +
    geom_text(aes(label = sprintf("%.3f", R2)),
              position = position_stack(vjust = .5), size = 4, colour = "white", fontface = "bold") +
    scale_fill_manual(values = c("Uniquely attributed (sum of marginals)" = "#2166AC",
                                 "Shared among collinear predictors"      = "#A6CEE3"), name = NULL) +
    coord_flip() +
    labs(title = "Total environmental variance: attributed vs. shared",
         subtitle = sprintf("Total model R\u00b2 = %.3f; the shared slice cannot be assigned to any single predictor",
                            env$R2_total),
         x = NULL, y = expression(R^2)) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
  ggsave(file.path(OUT, "env_permanova_shared.png"), fShare, width = 8.5, height = 2.6, dpi = 140, bg = "white")
  say("wrote env_permanova_shared.png")
}

# =============================================================================
# OPTION 1  single stacked "total variance" bar:
#   nine unique marginal segments + one aggregate SHARED block + residual.
#   Two versions: full (residual shown) and zoomed (explained-only).
# =============================================================================
if (!is.null(env$R2_total) && is.finite(env$R2_total)) {
  # pull p-values from the marginal table (column is "Pr(>F)")
  pcol <- if ("Pr(>F)" %in% names(env$tab)) "Pr(>F)" else
          grep("^Pr", names(env$tab), value = TRUE)[1]
  uni <- env$tab %>% filter(term %in% env$kept) %>%
    transmute(seg = term, R2, p = .data[[pcol]], kind = "unique") %>% arrange(desc(R2))
  shared_row <- data.frame(seg = "Shared (collinear)", R2 = env$R2_shared, p = NA, kind = "shared")
  resid_row  <- data.frame(seg = "Unexplained residual",
                           R2 = max(1 - env$R2_total, 0), p = NA, kind = "residual")

  # ordered factor: unique (largest->smallest), then shared, then residual
  stack_full <- rbind(uni, shared_row, resid_row)
  stack_full$seg <- factor(stack_full$seg, levels = rev(stack_full$seg))

  # DISTINCT categorical color per predictor (not a single-hue ramp)
  n_uni <- nrow(uni)
  cat_pal <- if (n_uni <= 9) {
    c("#4E79A7","#F28E2B","#59A14F","#E15759","#B07AA1",
      "#76B7B2","#EDC948","#FF9DA7","#9C755F")[seq_len(n_uni)]
  } else grDevices::hcl.colors(n_uni, "Dark 3")
  uni_cols <- setNames(cat_pal, uni$seg)
  fill_map <- c(uni_cols,
                "Shared (collinear)"   = "#BDBDBD",
                "Unexplained residual" = "#F0F0F0")

  # annotated legend labels: "predictor  (R2=0.016, p=0.016*)"
  star <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else
                      if (p < 0.01) "**" else if (p < 0.05) "*" else ""
  lab_for <- function(seg, R2, p) {
    if (seg == "Shared (collinear)")   return(sprintf("Shared among predictors  (R\u00b2=%.3f)", R2))
    if (seg == "Unexplained residual") return(sprintf("Unexplained residual  (R\u00b2=%.3f)", R2))
    sprintf("%s  (R\u00b2=%.3f, p=%.3f%s)", seg, R2, p, star(p))
  }
  lab_map <- setNames(
    mapply(lab_for, as.character(stack_full$seg), stack_full$R2, stack_full$p),
    as.character(stack_full$seg))

  base_stack <- function(dat, ttl, sub) {
    ggplot(dat, aes(x = "", y = R2, fill = seg)) +
      geom_col(width = .55, colour = "white", linewidth = .3) +
      geom_text(aes(label = ifelse(R2 > 0.006, sprintf("%.3f", R2), "")),
                position = position_stack(vjust = .5), size = 3, colour = "grey15") +
      scale_fill_manual(values = fill_map, labels = lab_map, name = NULL,
                        guide = guide_legend(reverse = TRUE, ncol = 1)) +
      coord_flip() +
      labs(title = ttl, subtitle = sub, x = NULL, y = expression(R^2)) +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            legend.text = element_text(size = 8.5),
            legend.key.size = unit(0.9, "lines"))
  }

  # full: includes residual (scaled to total trait variance = 1)
  fFull <- base_stack(stack_full,
    "Total trait variance: unique + shared + unexplained",
    sprintf("environment explains R\u00b2 = %.3f total; unique parts sum to %.3f, shared = %.3f, residual = %.3f\nlegend shows each predictor's unique marginal R\u00b2 and p-value (* p<0.05, ** p<0.01, *** p<0.001)",
            env$R2_total, env$R2_marg_sum, env$R2_shared, max(1 - env$R2_total, 0)))
  ggsave(file.path(OUT, "env_stack_full.png"), fFull, width = 10, height = 3.6, dpi = 140, bg = "white")

  # zoomed: explained variance only (drop residual so unique-vs-shared is legible)
  stack_expl <- rbind(uni, shared_row)
  stack_expl$seg <- factor(stack_expl$seg, levels = rev(stack_expl$seg))
  fZoom <- base_stack(stack_expl,
    "Explained trait variance only: unique + shared",
    sprintf("of the %.3f environment explains, %.0f%% is uniquely attributable and %.0f%% is shared among collinear predictors\nlegend shows each predictor's unique marginal R\u00b2 and p-value",
            env$R2_total, 100*env$R2_marg_sum/env$R2_total, 100*env$R2_shared/env$R2_total))
  ggsave(file.path(OUT, "env_stack_explained.png"), fZoom, width = 10, height = 3.6, dpi = 140, bg = "white")
  say("wrote env_stack_full.png and env_stack_explained.png")
}

# =============================================================================
# OPTION 3  incremental / sequential bar:
#   cumulative R2 as predictors enter one at a time (by="terms" order). Each
#   segment = variance ADDED given the ones before it; overlap shows up as each
#   predictor's sequential piece being smaller than its standalone marginal.
#   Requires the sequential decomposition, which we recompute deterministically
#   here (no permutations) from the cached kept-set on the real data.
# =============================================================================
seq_fp <- file.path(CACHE, "env_sequential.rds")
if (file.exists(seq_fp)) {
  sq <- readRDS(seq_fp)
  sq$term <- factor(sq$term, levels = sq$term)
  fSeq <- ggplot(sq, aes(x = reorder(term, cum_R2), y = seq_R2)) +
    geom_col(aes(fill = seq_R2), width = .7) +
    geom_text(aes(label = sprintf("+%.3f", seq_R2)), hjust = -0.15, size = 3) +
    scale_fill_gradient(low = "#C6DBEF", high = "#08519C", guide = "none") +
    coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, .2))) +
    labs(title = "Sequential (incremental) environmental contribution",
         subtitle = sprintf("variance ADDED by each predictor given those already entered; total = %.3f\n(compare each to its standalone marginal in the main bar \u2014 the gap is shared variance)",
                            sum(sq$seq_R2, na.rm = TRUE)),
         x = NULL, y = expression("Sequential " * R^2)) +
    theme_bw(base_size = 11) + theme(panel.grid.major.y = element_blank())
  ggsave(file.path(OUT, "env_sequential.png"), fSeq, width = 8.5, height = 4.5, dpi = 140, bg = "white")
  say("wrote env_sequential.png")
} else {
  say("(env_sequential.rds not cached; skipping Option 3 \u2014 see analysis/one-off to generate it)")
}

# =============================================================================
# varpart component + absorption bars (per scale)
# =============================================================================
draw_varpart <- function(GROUP) {
  tag <- tolower(GROUP)
  fp  <- file.path(CACHE, sprintf("varpart_%s.rds", tag))
  if (!file.exists(fp)) return(invisible())
  V <- readRDS(fp)

  comp <- data.frame(
    component = factor(c("Pure environment","Shared (env-mediated region)","Pure region (identity)"),
                       levels = c("Pure region (identity)","Shared (env-mediated region)","Pure environment")),
    R2 = c(V$pure_E, V$shared, V$resid_R))
  fD1 <- ggplot(comp, aes(x = "", y = R2, fill = component)) +
    geom_col(width = .6, colour = "white") +
    geom_text(aes(label = ifelse(R2 > 0.004, sprintf("%.3f", R2), "")),
              position = position_stack(vjust = .5), size = 3.4, colour = "white", fontface = "bold") +
    scale_fill_manual(values = c("Pure region (identity)"       = "#B2182B",
                                 "Shared (env-mediated region)" = "#4393C3",
                                 "Pure environment"             = "#92C5DE"), name = NULL) +
    coord_flip() +
    labs(title = sprintf("Variance partition of trait space (%s)", nice_scale(GROUP)),
         subtitle = sprintf("Total explained R\u00b2 = %.3f", V$total_expl),
         x = NULL, y = expression(R^2)) +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
  ggsave(file.path(OUT, sprintf("varpart_components_%s.png", tag)), fD1,
         width = 9, height = 3.4, dpi = 140, bg = "white")

  absorb <- data.frame(
    part = factor(c("Absorbed by environment","Remaining region explained variance"),
                  levels = c("Remaining region explained variance","Absorbed by environment")),
    frac = c(V$frac_absorbed, 1 - V$frac_absorbed))
  fD2 <- ggplot(absorb, aes(x = "", y = frac, fill = part)) +
    geom_col(width = .5, colour = "white") +
    geom_text(aes(label = sprintf("%.0f%%", 100 * frac)),
              position = position_stack(vjust = .5), size = 5, colour = "white", fontface = "bold") +
    scale_fill_manual(values = c("Absorbed by environment"  = "#4393C3",
                                 "Remaining region explained variance" = "#B2182B"), name = NULL) +
    coord_flip() +
    labs(title = sprintf("Where does the regional signal go? (%s)", nice_scale(GROUP)),
         subtitle = sprintf("Region attributed R\u00b2 = %.3f  \u2192  %.0f%% cross explained by environment;  regional difference given env. p = %s",
                            V$R2_R, 100 * V$frac_absorbed, format.pval(V$p_R_give_E)),
         x = NULL, y = "Fraction of region-alone R\u00b2") +
    theme_bw(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")
  ggsave(file.path(OUT, sprintf("varpart_absorption_%s.png", tag)), fD2,
         width = 9, height = 2.8, dpi = 140, bg = "white")

  # ---- VENN 1: total variance, region vs environment (per scale) ----
  residual <- max(1 - V$total_expl, 0)
  vennA <- venn_region_env(
    pure_region = V$resid_R, shared = V$shared, pure_env = V$pure_E, residual = residual,
    title = sprintf("Region \u2229 Environment partition of trait variance: %s", nice_scale(GROUP)),
    subtitle = sprintf("total explained R\u00b2 = %.3f; values are fractions of total trait variance", V$total_expl))
  ggsave(file.path(OUT, sprintf("varpart_venn_%s.png", tag)), vennA,
         width = 7, height = 5, dpi = 140, bg = "white")

  # ---- VENN 2: where the region-explained variance goes (per scale) ----
  # same circles, but framed on the region signal: of Region-alone R2, the
  # overlap is env-absorbed and the left-only is region-specific identity.
  vennB <- venn_region_env(
    pure_region = V$resid_R, shared = V$shared, pure_env = V$pure_E, residual = residual,
    title = sprintf("Where the region-explained variance goes (%s)", nice_scale(GROUP)),
    subtitle = sprintf("Region based R\u00b2 = %.3f = pure region (%.3f) + shared with env (%.3f);  %.0f%% is env-absorbed",
                       V$R2_R, V$resid_R, V$shared, 100 * V$frac_absorbed))
  ggsave(file.path(OUT, sprintf("varpart_venn_regiongoes_%s.png", tag)), vennB,
         width = 7, height = 5, dpi = 140, bg = "white")

  say(sprintf("wrote varpart figures + 2 Venns for %s", tag))
}
invisible(lapply(REGIONS, draw_varpart))

# =============================================================================
# NESTED main-effects ordination (Image-2 style): all scales in one trait space
# =============================================================================
nfp <- file.path(CACHE, "nested_ordination.rds")
if (file.exists(nfp)) {
  N <- readRDS(nfp)
  co <- N$coords
  co$scale <- factor(co$scale, levels = intersect(c("REALM","PROVINC","ECOREGION"), unique(co$scale)))
  shp <- c(REALM = 17, PROVINC = 16, ECOREGION = 15)   # triangle / circle / square
  sz  <- c(REALM = 6,  PROVINC = 4,  ECOREGION = 2.6)
  realm_pts <- subset(co, scale == "REALM")
  finer_pts <- subset(co, scale != "REALM")

  fN <- ggplot()
  if (!is.null(N$links))
    fN <- fN + geom_segment(data = N$links,
                            aes(x = x, y = y, xend = xend, yend = yend),
                            colour = "grey78", linewidth = 0.3)
  # finer-scale centroids in grey, shaped by scale
  if (nrow(finer_pts))
    fN <- fN + geom_point(data = finer_pts,
                          aes(MDS1, MDS2, shape = scale, size = scale),
                          colour = "grey45")
  # realm centroids colored by palette
  fN <- fN +
    geom_point(data = realm_pts, aes(MDS1, MDS2, colour = level, size = scale), shape = 17) +
    ggrepel::geom_text_repel(data = realm_pts, aes(MDS1, MDS2, label = level),
                             size = 2.8, fontface = "bold", max.overlaps = 40, colour = "grey15") +
    scale_shape_manual(values = shp, name = "Scale",
                       breaks = intersect(c("PROVINC","ECOREGION"), levels(co$scale))) +
    scale_size_manual(values = sz, guide = "none") +
    scale_colour_manual(values = realm_cols, na.value = "grey45", guide = "none") +
    labs(title = "Nested biogeographic structure in trait space",
         subtitle = sprintf("mMDS of among-centroid distances D[C]; realms colored/labeled (triangles), provinces (circles) & ecoregions (squares) in grey; lines link each unit to its parent; stress = %.3f",
                            N$stress),
         x = "mMDS 1", y = "mMDS 2") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT, "nested_ordination.png"), fN, width = 9.5, height = 7.5, dpi = 140, bg = "white")
  say("wrote nested_ordination.png")
}

# =============================================================================
# VARPART — statistical summary table
#   For each scale: how region's explained variance decomposes into the part
#   ENVIRONMENT also accounts for (shared / co-explained) vs. the part unique to
#   region (pure region, i.e. region | environment), plus the marginal test of
#   region after environment is partialled out.
# =============================================================================
meta <- readRDS(file.path(CACHE, "meta.rds"))
fmtp <- function(p) ifelse(is.na(p), "—",
                           ifelse(p < .001, "<0.001", sprintf("%.3f", p)))

vp_rows <- list()
for (g in meta$REGIONS) {
  fp <- file.path(CACHE, sprintf("varpart_%s.rds", tolower(g)))
  if (!file.exists(fp)) next
  v <- readRDS(fp)
  vp_rows[[length(vp_rows)+1]] <- data.frame(
    Scale                       = g,
    `Region alone (R2_R)`       = round(v$R2_R, 4),
    `Env alone (R2_E)`          = round(v$R2_E, 4),
    `Shared / co-explained`     = round(v$shared, 4),
    `Pure region (R|E)`         = round(v$resid_R, 4),
    `Pure env (E|R)`            = round(v$pure_E, 4),
    `% of region absorbed by env` = sprintf("%.1f%%", 100 * v$frac_absorbed),
    `Total explained`           = round(v$total_expl, 4),
    `p(region | env)`           = fmtp(v$p_R_give_E),
    check.names = FALSE, stringsAsFactors = FALSE)
}
vp_table <- do.call(rbind, vp_rows)

# order coarse -> fine
ord_lv <- intersect(c("REALM","PROVINC","ECOREGION"), vp_table$Scale)
vp_table <- vp_table[match(ord_lv, vp_table$Scale), ]

write.csv(vp_table, file.path(OUT, "TABLE_varpart_stats.csv"), row.names = FALSE)
say("\nwrote TABLE_varpart_stats.csv"); print(vp_table, row.names = FALSE)

if (requireNamespace("flextable", quietly = TRUE)) {
  library(flextable)
  ft <- flextable(vp_table) %>%
    add_header_lines(
      "Variance partitioning of trait dissimilarity: region vs. environment (PERMANOVA, by-terms/by-margin)") %>%
    bold(part = "header") %>% autofit()
  save_as_docx(ft, path = file.path(OUT, "TABLE_varpart_stats.docx"))
  say("wrote TABLE_varpart_stats.docx")
}

# =============================================================================
# RESIDUAL ordination (region | environment) per scale
# =============================================================================
draw_residual <- function(GROUP) {
  tag <- tolower(GROUP)
  fp  <- file.path(CACHE, sprintf("residual_%s.rds", tag))
  if (!file.exists(fp)) return(invisible())
  Rr <- readRDS(fp); co <- Rr$coords
  fR <- ggplot(co, aes(MDS1, MDS2, colour = grp)) +
    geom_hline(yintercept = 0, colour = "grey90") +
    geom_vline(xintercept = 0, colour = "grey90") +
    geom_point(aes(size = n)) +
    ggrepel::geom_text_repel(aes(label = grp), size = 3, fontface = "bold",
                             max.overlaps = 40, show.legend = FALSE) +
    scale_size(range = c(3, 9), name = "n") +
    scale_for(GROUP, co$grp, "colour") +
    labs(title = sprintf("Residual %s structure after removing environment", tolower(nice_scale(GROUP))),
         subtitle = sprintf("mMDS of region centroids in the environment-residualized space (I\u2212H)G(I\u2212H); stress = %.3f\nspread here = trait differences among regions NOT explained by environment",
                            Rr$stress),
         x = "mMDS 1", y = "mMDS 2") +
    guides(colour = "none") + theme_bw(base_size = 11)
  ggsave(file.path(OUT, sprintf("residual_ordination_%s.png", tag)), fR,
         width = 9, height = 7, dpi = 140, bg = "white")
  say(sprintf("wrote residual_ordination_%s.png", tag))
}
invisible(lapply(REGIONS, draw_residual))

# =============================================================================
# R2-by-scale trend: how explanatory power changes as grouping gets finer
# =============================================================================
tfp <- file.path(CACHE, "r2_by_scale.rds")
if (file.exists(tfp)) {
  tr <- readRDS(tfp)
  tr$scale <- factor(tr$scale, levels = intersect(c("REALM","PROVINC","ECOREGION"), tr$scale))
  lab <- sprintf("%s\n(%d groups)", vapply(as.character(tr$scale), nice_scale, character(1)), tr$groups)
  tr$xlab <- factor(lab, levels = lab[order(tr$scale)])
  tr$sig <- ifelse(!is.na(tr$adonis_p) & tr$adonis_p < 0.05, "p<0.05", "n.s.")

  long <- rbind(
    data.frame(xlab = tr$xlab, metric = "Total region R\u00b2", value = tr$permanova_R2),
    data.frame(xlab = tr$xlab, metric = "Pure region (| env)", value = tr$pure_region))
  fT <- ggplot(long, aes(xlab, value, colour = metric, group = metric)) +
    geom_line(linewidth = 1) + geom_point(size = 3.4) +
    geom_text(data = tr, aes(xlab, permanova_R2,
              label = sprintf("R\u00b2=%.3f\n%s", permanova_R2, sig)),
              inherit.aes = FALSE, vjust = -0.6, size = 3) +
    scale_colour_manual(values = c("Total region R\u00b2" = "#B2182B",
                                   "Pure region (| env)" = "#2166AC"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(title = "How regional structure changes across nesting scales",
         subtitle = "coarse \u2192 fine; total region R\u00b2 vs. the pure-region fraction after removing environment",
         x = NULL, y = expression(R^2)) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(OUT, "r2_by_scale.png"), fT, width = 8.5, height = 5.5, dpi = 140, bg = "white")
  say("wrote r2_by_scale.png")
}

# =============================================================================
# DISTANCE-DECAY figures per scale:
#   (1) trait divergence D[C] vs geographic distance, colored sig/n.s., with
#       LOESS/linear fit and simple + partial Mantel r annotations
#   (2) logistic P(significant) vs geographic distance
# =============================================================================
draw_distdecay <- function(GROUP) {
  tag <- tolower(GROUP)
  fp  <- file.path(CACHE, sprintf("distdecay_%s.rds", tag))
  if (!file.exists(fp)) return(invisible())
  DD <- readRDS(fp); p <- DD$pairs
  p$sig_lab <- factor(ifelse(p$sig %in% TRUE, "sig (Holm p<0.05)", "n.s."),
                      levels = c("sig (Holm p<0.05)", "n.s."))

  nfull <- if (!is.null(DD$n_regions_full)) DD$n_regions_full else length(unique(c(p$r1, p$r2)))
  nenv  <- if (!is.null(DD$n_regions_env)) DD$n_regions_env else NA
  n_sig <- sum(p$sig %in% TRUE); n_ns <- sum(!(p$sig %in% TRUE))

  mr  <- if (!is.null(DD$mantel))  sprintf("simple Mantel r = %.2f, p = %s (%d %ss, all pairs)",
             DD$mantel$statistic, format.pval(DD$mantel$signif), nfull, tolower(nice_scale(GROUP))) else "simple Mantel: NA"
  pmr <- if (!is.null(DD$pmantel)) sprintf("partial Mantel (| env) r = %.2f, p = %s (%d env-complete %ss)",
             DD$pmantel$statistic, format.pval(DD$pmantel$signif), nenv, tolower(nice_scale(GROUP))) else
             sprintf("partial Mantel: not computed (only %s env-complete %ss)", nenv, tolower(nice_scale(GROUP)))

  # (1) distance-decay scatter  (FULL trait-complete set)
  f1 <- ggplot(p, aes(geo_km, traitDC)) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey30", fill = "grey85", linewidth = .6) +
    geom_point(aes(colour = sig_lab), size = 2.6, alpha = .85) +
    scale_colour_manual(values = c("sig (Holm p<0.05)" = "#B2182B", "n.s." = "#4393C3"),
                        drop = FALSE, name = NULL) +
    labs(title = sprintf("Morphological divergence vs. geographic distance (%s)", nice_scale(GROUP)),
         subtitle = sprintf("%d pairs (%d sig, %d n.s.);  %s\n%s",
                            nrow(p), n_sig, n_ns, mr, pmr),
         x = "Great-circle distance between region centroids (km)",
         y = "Trait-space centroid distance D[C]") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(OUT, sprintf("distdecay_scatter_%s.png", tag)), f1,
         width = 8.5, height = 6, dpi = 140, bg = "white")

  say(sprintf("wrote distance-decay figures for %s", tag))
}
invisible(lapply(REGIONS, draw_distdecay))



# =============================================================================
# ENVIRONMENTAL MODEL — publication tables
#   (A) Marginal PERMANOVA of the standalone environmental model, one row per
#       retained predictor + a Shared/Total/Residual summary block.
#   (B) Reconciliation table: environmental R2 across frames (standalone vs.
#       each scale's varpart), with N and retained predictor set, since these
#       differ by complete-case filter + MIN_N + per-subset pruning.
#   Writes CSVs + a formatted flextable (docx/png).
# =============================================================================
env  <- readRDS(file.path(CACHE, "env.rds"))
meta <- readRDS(file.path(CACHE, "meta.rds"))

pcol <- if ("Pr(>F)" %in% names(env$tab)) "Pr(>F)" else grep("^Pr", names(env$tab), value = TRUE)[1]

stars <- function(p) ifelse(is.na(p), "",
                            ifelse(p < .001, "***", ifelse(p < .01, "**",
                                                           ifelse(p < .05, "*", "ns"))))
fmtp  <- function(p) ifelse(is.na(p), "—",
                            ifelse(p < .001, "<0.001", sprintf("%.3f", p)))

pretty <- c(SST_clim = "Sea surface temp.", nitrate_clim = "Nitrate",
            phosphate_clim = "Phosphate", salinity_clim = "Salinity",
            wave.height = "Wave height", kd490_clim = "Turbidity (Kd490)",
            PAR_mean = "PAR", iron_clim = "Iron", spco2_clim = "Surface ocean pCO\u2082")
lab <- function(x) ifelse(x %in% names(pretty), pretty[x], x)

# --- (A) per-predictor marginal rows -----------------------------------------
marg <- env$tab %>%
  filter(term %in% env$kept) %>%
  transmute(
    Predictor  = lab(term),
    Df, SumOfSqs = round(SumOfSqs, 3),
    `Pseudo-F` = round(F, 3),
    R2         = round(R2, 4),
    P = .data[[pcol]], sig = stars(.data[[pcol]])) %>%
  arrange(desc(R2))
if (!is.null(env$ci)) {
  marg <- marg %>% left_join(
    env$ci %>% transmute(Predictor = lab(term),
                         `95% CI` = sprintf("[%.3f, %.3f]", lo, hi)),
    by = "Predictor")
}
marg$P <- fmtp(marg$P)

# --- (A) summary block --------------------------------------------------------
resid_R2 <- max(1 - env$R2_total, 0)
summ <- data.frame(
  Predictor = c("Sum of marginal (unique)", "Shared (collinear)",
                "Total model", "Residual"),
  Df = NA, SumOfSqs = NA, `Pseudo-F` = NA,
  R2 = round(c(env$R2_marg_sum, env$R2_shared, env$R2_total, resid_R2), 4),
  P = "—", sig = "", `95% CI` = "", check.names = FALSE)
if (is.null(env$ci)) { marg$`95% CI` <- NULL; summ$`95% CI` <- NULL }

env_table <- bind_rows(marg, summ)
write.csv(env_table, file.path(OUT, "TABLE_env_permanova.csv"), row.names = FALSE)
say("\nwrote TABLE_env_permanova.csv"); print(env_table, row.names = FALSE)

# --- (A) formatted flextable (docx + png) ------------------------------------
env_n <- sum(env$adonis$Df) + 1L
if (requireNamespace("flextable", quietly = TRUE)) {
  library(flextable)
  n_marg <- nrow(marg)
  ft <- flextable(env_table) %>%
    add_header_lines(sprintf(
      "Marginal PERMANOVA, standalone environmental model (n = %d; %d permutations; Euclidean on 7 scaled traits)",
      env_n, meta$NPERM)) %>%
    bold(part = "header") %>%
    hline(i = n_marg, border = officer::fp_border(width = 1.2)) %>%
    align(j = c("Df","SumOfSqs","Pseudo-F","R2","P"), align = "right", part = "all") %>%
    colformat_double(j = "R2", digits = 4) %>%
    autofit()
  save_as_docx(ft, path = file.path(OUT, "TABLE_env_permanova.docx"))
  tryCatch(save_as_image(ft, path = file.path(OUT, "TABLE_env_permanova.png")),
           error = function(e) say("  (png export needs webshot2/magick; docx written)"))
  say("wrote TABLE_env_permanova.docx")
} else {
  say("  (install 'flextable' for formatted docx/png output)")
}

# --- (B) environmental R2 across frames ---------------------------------------
rows <- list(
  data.frame(
    Frame        = "Standalone env. model",
    Scale        = "—  (no region filter)",
    N            = env_n,
    k_predictors = length(env$kept),
    env_R2       = round(env$R2_total, 4),
    Retained     = paste(env$kept, collapse = ", "),
    stringsAsFactors = FALSE)
)
for (g in meta$REGIONS) {
  fp <- file.path(CACHE, sprintf("varpart_%s.rds", tolower(g)))
  if (!file.exists(fp)) next
  vp <- readRDS(fp)
  d <- rec
  for (t in LOG) if (t %in% names(d)) d[[t]] <- log(d[[t]])
  d[[g]][trimws(tolower(d[[g]])) %in% c("other","","na","unknown")] <- NA
  cc <- complete.cases(d[, meta$TRAITS], d[, meta$PREDICTORS]) & !is.na(d[[g]])
  d <- d[cc, ]
  keep_grp <- names(which(table(d[[g]]) >= meta$MIN_N))
  n_g <- sum(d[[g]] %in% keep_grp)
  rows[[length(rows)+1]] <- data.frame(
    Frame        = "Varpart (env alone)",
    Scale        = g,
    N            = n_g,
    k_predictors = length(vp$kept),
    env_R2       = round(vp$R2_E, 4),
    Retained     = paste(vp$kept, collapse = ", "),
    stringsAsFactors = FALSE)
}
env_frames <- do.call(rbind, rows)
write.csv(env_frames, file.path(OUT, "TABLE_env_R2_by_frame.csv"), row.names = FALSE)
say("\nwrote TABLE_env_R2_by_frame.csv"); print(env_frames, row.names = FALSE)

say("\nALL FIGURES DONE. PNGs in:", OUT)
