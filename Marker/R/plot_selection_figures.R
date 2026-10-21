# =============================================================================
# plot_selection_figures.R
# Survivor-selection falsification: main + supplementary figures
# Style follows plot_density_main_figure.R (PI-approved): master palette,
# donor-level violin + jitter, clean stars, legends inside panels.
# =============================================================================

library(tidyverse)
library(patchwork)

source("colors_master.R")

theme_pub <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    plot.tag = element_text(size = 20, face = "bold"),
    plot.margin = margin(8, 12, 8, 8),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank()
  )

p_to_stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                                    p < 0.05 ~ "*", TRUE ~ "")

N_BINS <- 10

# =============================================================================
# DATA PREP — recomputes everything from Data/, so this file is standalone
# =============================================================================

prepare_selection_data <- function(n_bins = N_BINS) {
  read_any <- function(stem) {
    for (p in sprintf("Data/%s.%s", stem, c("rds", "csv"))) {
      if (file.exists(p)) return(if (grepl("rds$", p)) readRDS(p)
                                 else readr::read_csv(p, show_col_types = FALSE))
    }
    stop("cannot find Data/", stem)
  }
  norm_img <- function(x) as.character(suppressWarnings(as.numeric(as.character(x))))

  comp <- read_any("ratio_data") %>%
    mutate(Diagnosis  = factor(Diagnosis, levels = c("ND", "T1D")),
           Donor      = as.character(Donor),
           ImageID    = norm_img(ImageID),
           alpha_frac = Count.glu / (Count.ins + Count.glu)) %>%
    filter(is.finite(alpha_frac))

  dens <- read_any("density_data") %>%
    mutate(Donor = as.character(Donor), ImageID = norm_img(ImageID)) %>%
    dplyr::select(Donor, ImageID, full_tissue_area) %>%
    distinct()

  donor_area <- comp %>%
    distinct(Donor, Diagnosis, ImageID) %>%
    inner_join(dens, by = c("Donor", "ImageID")) %>%
    group_by(Donor, Diagnosis) %>%
    summarise(area_mm2 = sum(full_tissue_area), .groups = "drop")

  edges <- seq(0, 1, length.out = n_bins + 1)

  counts_by_bin <- comp %>%
    inner_join(dens, by = c("Donor", "ImageID")) %>%
    mutate(bin = cut(alpha_frac, edges, include.lowest = TRUE)) %>%
    count(Donor, Diagnosis, bin, name = "n_islets") %>%
    complete(nesting(Donor, Diagnosis), bin, fill = list(n_islets = 0)) %>%
    inner_join(donor_area, by = c("Donor", "Diagnosis")) %>%
    mutate(islets_per_mm2 = n_islets / area_mm2)

  # Bin midpoints give a continuous x-axis, which reads better than factor
  # labels when the panels are stacked.
  mids <- (edges[-length(edges)] + edges[-1]) / 2
  bin_lookup <- tibble(bin = levels(counts_by_bin$bin), mid = mids)

  ineq <- counts_by_bin %>%
    group_by(Diagnosis, bin) %>%
    summarise(mean_density = mean(islets_per_mm2), .groups = "drop") %>%
    mutate(bin = as.character(bin)) %>%
    left_join(bin_lookup, by = "bin")

  ineq_wide <- ineq %>%
    dplyr::select(bin, mid, Diagnosis, mean_density) %>%
    pivot_wider(names_from = Diagnosis, values_from = mean_density) %>%
    mutate(ratio = T1D / ND)

  shape <- comp %>%
    mutate(bin = as.character(cut(alpha_frac, edges, include.lowest = TRUE))) %>%
    count(Diagnosis, bin, name = "n") %>%
    group_by(Diagnosis) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    left_join(bin_lookup, by = "bin")

  burden <- comp %>%
    inner_join(dens, by = c("Donor", "ImageID")) %>%
    group_by(Donor, Diagnosis) %>%
    summarise(islets = n(), alpha = sum(Count.glu), beta = sum(Count.ins),
              endocrine = sum(Count.ins + Count.glu), .groups = "drop") %>%
    inner_join(donor_area, by = c("Donor", "Diagnosis")) %>%
    mutate(islets_per_mm2 = islets / area_mm2,
           alpha_per_mm2  = alpha / area_mm2,
           beta_per_mm2   = beta / area_mm2,
           endo_per_mm2   = endocrine / area_mm2)

  # Section area differs by diagnosis (exocrine atrophy), which inflates every
  # T1D density. Carried through so panel D can state the correction.
  ta <- comp %>%
    distinct(Donor, Diagnosis, ImageID) %>%
    inner_join(dens, by = c("Donor", "ImageID")) %>%
    group_by(Donor, Diagnosis) %>%
    summarise(mean_area = mean(full_tissue_area), .groups = "drop")
  ar <- ta %>% group_by(Diagnosis) %>%
    summarise(m = median(mean_area), .groups = "drop")
  area_fold <- ar$m[ar$Diagnosis == "T1D"] / ar$m[ar$Diagnosis == "ND"]

  # Every number that appears in a figure subtitle is computed here, so a
  # subtitle cannot drift out of step with the data behind it.
  nd_af <- comp$alpha_frac[comp$Diagnosis == "ND"]
  t1_af <- comp$alpha_frac[comp$Diagnosis == "T1D"]
  med <- function(v, g) median(v[burden$Diagnosis == g])
  stats <- list(
    pct_nd_above_t1d_median = 100 * mean(nd_af >= median(t1_af)),
    alpha_fold  = med(burden$alpha_per_mm2, "T1D") / med(burden$alpha_per_mm2, "ND"),
    beta_fold   = med(burden$beta_per_mm2,  "T1D") / med(burden$beta_per_mm2,  "ND"),
    islet_fold  = med(burden$islets_per_mm2, "T1D") / med(burden$islets_per_mm2, "ND"),
    endo_fold   = med(burden$endo_per_mm2,  "T1D") / med(burden$endo_per_mm2,  "ND"),
    max_ratio   = max(ineq_wide$ratio[is.finite(ineq_wide$ratio)]))

  list(comp = comp, ineq = ineq, ineq_wide = ineq_wide, shape = shape,
       burden = burden, area_fold = area_fold, edges = edges, stats = stats)
}

# =============================================================================
# MAIN FIGURE
# =============================================================================

# ── Panel A: composition shape ───────────────────────────────────────────────

plot_sel_a <- function(d) {
  ggplot(d$shape, aes(x = mid, y = 100 * prop, fill = Diagnosis)) +
    geom_col(position = position_dodge(0.08), width = 0.075, alpha = 0.85) +
    scale_fill_manual(values = colors_diag) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(x = "Alpha-cell fraction of islet", y = "Islets (% within group)",
         tag = "A", title = "Composition shifts almost completely",
         subtitle = sprintf(
           "Only %.2f%% of ND islets are as alpha-rich as the median T1D islet",
           d$stats$pct_nd_above_t1d_median),
         fill = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.25, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# ── Panel B: the density inequality ──────────────────────────────────────────

plot_sel_b <- function(d) {
  # Log y is essential: 0.015 vs 1.251 islets/mm2 in the same bin is invisible
  # on a linear axis.
  ggplot(d$ineq, aes(x = mid, y = mean_density, color = Diagnosis)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    scale_color_manual(values = colors_diag) +
    scale_y_log10(labels = scales::label_number(accuracy = 0.001)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(x = "Alpha-cell fraction of islet", y = "Islets per mm\u00B2 (log scale)",
         tag = "B", title = "Absolute density of alpha-rich islets rises",
         subtitle = "Selection can only remove islets, so T1D cannot exceed ND",
         color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.25, 0.20),
          legend.background = element_rect(fill = "white", color = NA))
}

# ── Panel C: the impossibility ───────────────────────────────────────────────

plot_sel_c <- function(d) {
  dd <- d$ineq_wide %>% filter(is.finite(ratio), ratio > 0)
  ggplot(dd, aes(x = mid, y = ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5,
               color = "gray30") +
    geom_line(linewidth = 0.9, color = "black") +
    geom_point(aes(color = ratio > 1), size = 2.6) +
    scale_color_manual(values = c("FALSE" = "gray55", "TRUE" = col_t1d),
                       guide = "none") +
    scale_y_log10(labels = scales::label_number(accuracy = 0.01)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    annotate("text", x = 0.95, y = max(dd$ratio) * 0.55,
             label = sprintf("%.0f\u00D7", max(dd$ratio)),
             size = 4, fontface = "bold", color = col_t1d) +
    annotate("text", x = 0.30, y = 1.45, label = "survival probability = 1",
             size = 3, color = "gray30") +
    labs(x = "Alpha-cell fraction of islet", y = "T1D / ND density ratio",
         tag = "C", title = "No survival probability can produce this",
         subtitle = "Points above the dashed line require probabilities > 1") +
    theme_pub
}

# ── Panel D: absolute cell burden ────────────────────────────────────────────

plot_sel_d <- function(d) {
  bl <- d$burden %>%
    dplyr::select(Donor, Diagnosis, beta_per_mm2, alpha_per_mm2, endo_per_mm2) %>%
    pivot_longer(-c(Donor, Diagnosis), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric,
      levels = c("beta_per_mm2", "alpha_per_mm2", "endo_per_mm2"),
      labels = c("Beta cells", "Alpha cells", "Total endocrine")))

  stars <- bl %>%
    group_by(metric) %>%
    summarise(
      p = suppressWarnings(wilcox.test(value[Diagnosis == "T1D"],
                                       value[Diagnosis == "ND"])$p.value),
      y = max(value, na.rm = TRUE) * 0.95, .groups = "drop") %>%
    mutate(lab = p_to_stars(p))

  ggplot(bl, aes(x = Diagnosis, y = value)) +
    geom_violin(aes(fill = Diagnosis), alpha = 0.25, color = NA, width = 0.7) +
    geom_jitter(aes(color = Diagnosis), width = 0.15, size = 1.5, alpha = 0.6) +
    geom_text(data = stars, aes(x = 1.5, y = y, label = lab), size = 8,
              inherit.aes = FALSE) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = NULL, y = "Cells per mm\u00B2 of tissue", tag = "D",
         title = "Alpha cells were created, not revealed",
         subtitle = sprintf(
           "Alpha %.2f\u00D7 higher in T1D (%.2f\u00D7 after correcting for %.0f%% smaller T1D sections)",
           d$stats$alpha_fold, d$stats$alpha_fold / d$area_fold,
           100 * (1 - d$area_fold))) +
    theme_pub
}

# =============================================================================
# SUPPLEMENTARY FIGURE
# =============================================================================

# ── S1: matched alpha-rich islets ────────────────────────────────────────────

plot_sel_s1 <- function(d) {
  t1 <- d$comp$alpha_frac[d$comp$Diagnosis == "T1D"]
  lo <- quantile(t1, .25); hi <- quantile(t1, .75)
  mt <- d$comp %>%
    filter(alpha_frac >= lo, alpha_frac <= hi) %>%
    mutate(total = Count.ins + Count.glu)

  p <- suppressWarnings(wilcox.test(mt$total[mt$Diagnosis == "T1D"],
                                    mt$total[mt$Diagnosis == "ND"])$p.value)

  ggplot(mt, aes(x = Diagnosis, y = total)) +
    geom_violin(aes(fill = Diagnosis), alpha = 0.25, color = NA, width = 0.7) +
    stat_summary(fun = median, geom = "crossbar", width = 0.35,
                 linewidth = 0.4, color = "black") +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_y_log10() +
    annotate("text", x = 1.5, y = max(mt$total) * 0.7, label = p_to_stars(p),
             size = 8) +
    labs(x = NULL, y = "Endocrine cells per islet (log scale)", tag = "A",
         title = "Alpha-rich ND islets vs T1D islets",
         subtitle = sprintf("Matched on alpha fraction %.2f\u2013%.2f; ND n = %d, T1D n = %d",
                            lo, hi, sum(mt$Diagnosis == "ND"),
                            sum(mt$Diagnosis == "T1D"))) +
    theme_pub
}

# ── S2: composition vs disease duration ──────────────────────────────────────

plot_sel_s2 <- function(d) {
  dd <- d$comp %>%
    filter(Diagnosis == "T1D") %>%
    group_by(Donor, Disease.Duration) %>%
    summarise(alpha_frac = mean(alpha_frac), islets = n(), .groups = "drop") %>%
    mutate(duration_bin = cut(Disease.Duration, breaks = c(0, 5, 15, Inf),
                              labels = c("0\u20135 yr", "5\u201315 yr", ">15 yr")))
  ct <- suppressWarnings(cor.test(dd$Disease.Duration, dd$alpha_frac,
                                  method = "spearman"))

  ggplot(dd, aes(x = Disease.Duration, y = alpha_frac)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7,
                alpha = 0.15) +
    geom_point(aes(color = duration_bin, size = islets), alpha = 0.8) +
    scale_color_manual(values = colors_duration, name = "Disease\nduration") +
    scale_size_continuous(range = c(1.5, 4.5), guide = "none") +
    labs(x = "Disease duration (years)", y = "Mean alpha-cell fraction",
         tag = "B", title = "Composition tracks duration",
         subtitle = sprintf(
           "Spearman rho = %.2f, p = %.1e \u2014 but a ratio rises if beta falls, so this does not isolate conversion",
           ct$estimate, ct$p.value)) +
    theme_pub +
    theme(legend.position = c(0.85, 0.35),
          legend.background = element_rect(fill = "white", color = NA))
}

# ── S3: islet density, the other half of "preserved cellularity" ────────────

plot_sel_s3 <- function(d) {
  p <- suppressWarnings(wilcox.test(
    d$burden$islets_per_mm2[d$burden$Diagnosis == "T1D"],
    d$burden$islets_per_mm2[d$burden$Diagnosis == "ND"])$p.value)

  ggplot(d$burden, aes(x = Diagnosis, y = islets_per_mm2)) +
    geom_violin(aes(fill = Diagnosis), alpha = 0.25, color = NA, width = 0.7) +
    geom_jitter(aes(color = Diagnosis), width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    annotate("text", x = 1.5, y = max(d$burden$islets_per_mm2) * 0.95,
             label = p_to_stars(p), size = 8) +
    labs(x = NULL, y = "Islets per mm\u00B2", tag = "C",
         title = "Cellularity is preserved per islet, not per organ",
         subtitle = sprintf(
           "Islet density falls %.0f%%; total endocrine cells per mm\u00B2 fall %.0f%%",
           100 * (1 - d$stats$islet_fold), 100 * (1 - d$stats$endo_fold))) +
    theme_pub
}

# =============================================================================
# GENERATE
# =============================================================================

generate_selection_main <- function(d = NULL, save = TRUE,
                                    output_dir = "Selection/Graphs") {
  if (is.null(d)) { cat("Preparing data...\n"); d <- prepare_selection_data() }
  cat("Generating selection main figure...\n")
  p_a <- plot_sel_a(d); cat("  A\n")
  p_b <- plot_sel_b(d); cat("  B\n")
  p_c <- plot_sel_c(d); cat("  C\n")
  p_d <- plot_sel_d(d); cat("  D\n")

  fig <- (p_a | p_b) / (p_c | p_d) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

  if (save) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    for (ext in c("pdf", "png", "tiff")) {
      ggsave(file.path(output_dir, paste0("selection_main_figure.", ext)), fig,
             width = 13, height = 9, dpi = 300,
             bg = if (ext == "pdf") NA else "white")
    }
    cat("Saved to", output_dir, "\n")
  }
  list(fig = fig, panels = list(a = p_a, b = p_b, c = p_c, d = p_d), data = d)
}

generate_selection_supp <- function(d = NULL, save = TRUE,
                                    output_dir = "Selection/Graphs") {
  if (is.null(d)) { cat("Preparing data...\n"); d <- prepare_selection_data() }
  cat("Generating selection supplementary figure...\n")
  s1 <- plot_sel_s1(d); cat("  A\n")
  s2 <- plot_sel_s2(d); cat("  B\n")
  s3 <- plot_sel_s3(d); cat("  C\n")

  fig <- (s1 | s2 | s3) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

  if (save) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    for (ext in c("pdf", "png", "tiff")) {
      ggsave(file.path(output_dir, paste0("selection_supp_figure.", ext)), fig,
             width = 15, height = 5, dpi = 300,
             bg = if (ext == "pdf") NA else "white")
    }
    cat("Saved to", output_dir, "\n")
  }
  list(fig = fig, panels = list(a = s1, b = s2, c = s3), data = d)
}

# Usage:
#   source("plot_selection_figures.R")
#   d    <- prepare_selection_data()
#   main <- generate_selection_main(d)
#   supp <- generate_selection_supp(d)
