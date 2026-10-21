# =============================================================================
# plot_marker_main_figure.R
# Main figure for the six per-cell rate models.
# Style follows plot_density_main_figure.R (PI-approved).
#
# A  Forest plot of all six rate ratios (the headline)
# B  Ki-67 rates by diagnosis, faceted by cell type
# C  CD3 rates by diagnosis, faceted by cell type
# D  The substitution arithmetic that motivates the whole analysis
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(patchwork)

source("colors_master.R")
source("Marker/R/marker_common.R")   # registry, data loading, model resolution
source("Marker/R/marker_effects.R")  # verified marginal effects

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

# =============================================================================
# DATA + MODEL COLLECTION
# =============================================================================

collect_marker_fits <- function(model_dir = MARKER_MODEL_DIR,
                                stain_handling = "double_only") {
  raw  <- load_marker_data(stain_handling)
  fits <- list()
  for (def in marker_defs) {
    full <- raw[[sub("_.*", "", def$key)]]
    r    <- load_rate_model(def, model_dir)
    if (is.null(r)) {
      cat("  no cached model for", def$label, "\n")
    } else {
      cat(sprintf("  %-12s %-6s (%s)\n", def$label, r$family, r$how))
    }
    fits[[def$label]] <- list(def = def, data = marker_subset(full, def),
                              full = full,
                              model = if (is.null(r)) NULL else r$model,
                              family = if (is.null(r)) NA_character_ else r$family)
  }
  list(fits = fits, raw = raw)
}

MARKER_FIGURE_VERSION <- "marginal_diagnosis (verified)"

# Delegates to marker_common.R so this figure and the all-comparisons plots
# cannot report different stars for the same comparison. They did: panel B
# showed ** for Ki-67 alpha from the coefficient table (p = 0.0063, the effect
# at Head/Female only) while 01_main_diagnosis showed *** from the emmeans
# contrast averaged over Region and Sex (p < 0.001).
rate_ratio_ci <- function(model, data = NULL) {
  marginal_diagnosis(model)
}

marker_rates_table <- function(coll) {
  bind_rows(lapply(coll$fits, function(f) {
    rr <- rate_ratio_ci(f$model, f$data)
    if (is.null(rr)) return(NULL)
    obs <- f$data %>%
      group_by(Diagnosis) %>%
      summarise(rate = PER * sum(.data[[f$def$num]]) / sum(.data[[f$def$den]]),
                .groups = "drop")
    rr %>% mutate(marker = f$def$marker, cell = f$def$cell,
                  label = f$def$label,
                  rate_ND  = obs$rate[obs$Diagnosis == "ND"],
                  rate_T1D = obs$rate[obs$Diagnosis == "T1D"],
                  donors = dplyr::n_distinct(f$data$Donor), .before = 1)
  }))
}

# =============================================================================
# PANELS
# =============================================================================

plot_marker_a <- function(tab) {
  d <- tab %>%
    mutate(cell = factor(cell, levels = c("acinar", "alpha", "beta")),
           marker = factor(marker, levels = c("Ki-67", "CD3")),
           sig = p_to_stars(p))

  need <- c("ratio", "lwr", "upr", "p")
  miss <- setdiff(need, names(tab))
  if (length(miss) > 0) {
    stop("marker_rates_table() produced no ", paste(miss, collapse = "/"),
         " column.\n  If you are seeing this straight after an update, an older ",
         "copy of\n  plot_marker_main_figure.R is still being sourced. Check ",
         "getwd() and that\n  Marker/R/ holds the current file (it should define ",
         "MARKER_FIGURE_VERSION).")
  }

  if (any(is.na(d$lwr))) {
    warning("no confidence interval for: ",
            paste(d$label[is.na(d$lwr)], collapse = ", "),
            " — plotted as points without ranges")
  }

  ggplot(d, aes(x = ratio, y = cell, color = marker)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5,
               color = "gray30") +
    # Point for every model; range only where an interval was obtained, so a
    # model with no CI still appears rather than being dropped.
    geom_point(size = 2.2, position = position_dodge(0.5)) +
    geom_pointrange(aes(xmin = lwr, xmax = upr), linewidth = 0.7, size = 0.5,
                    position = position_dodge(0.5), na.rm = TRUE) +
    geom_text(aes(label = sig), position = position_dodge(0.5),
              vjust = -0.9, size = 5, show.legend = FALSE) +
    scale_x_log10() +
    scale_color_manual(values = c("Ki-67" = col_head, "CD3" = col_tail),
                       name = NULL) +
    labs(x = "Rate ratio, T1D vs ND (log scale)", y = NULL, tag = "A",
         title = "Per-cell rates: what changes and what does not",
         subtitle = "Each point is a rate per cell of that type; offset removes the denominator") +
    theme_pub +
    theme(legend.position = c(0.88, 0.25),
          legend.background = element_rect(fill = "white", color = NA))
}

# cells: which bins to show. The acinar bin is excluded from the main figure by
#   default — for Ki-67 it was uninformative (p = 0.12, and the crude, donor-
#   level and model estimates pointed in three different directions), and for
#   both markers Count.Acinar is the weakest of the three denominators: 35% of
#   islets have none at all, and numerator-exceeds-denominator runs 1.8-2.2%
#   against 0.02% for insulin. It is still fitted, still in the Excel workbooks
#   and still in All_Comparisons/.
# scales: "free_y" lets each cell type use its own range, which matters because
#   the rates differ severalfold between bins. Use "fixed" to compare bins
#   directly at the cost of squashing the smaller one.
plot_marker_bc <- function(coll, which_marker, tag, title, subtitle = NULL,
                           cells = c("beta", "alpha"), scales = "free_y") {
  keep <- Filter(function(f) f$def$marker == which_marker &&
                   f$def$cell %in% cells && !is.null(f$model),
                 coll$fits)
  if (length(keep) == 0) return(NULL)

  dr <- bind_rows(lapply(keep, function(f) {
    f$data %>%
      group_by(Donor, Diagnosis) %>%
      summarise(num = sum(.data[[f$def$num]]), den = sum(.data[[f$def$den]]),
                .groups = "drop") %>%
      filter(den > 0) %>%
      mutate(rate = PER * num / den, cell = f$def$cell)
  })) %>%
    mutate(cell = factor(cell, levels = cells))

  stars <- bind_rows(lapply(keep, function(f) {
    rr <- rate_ratio_ci(f$model, f$data)
    if (is.null(rr)) return(NULL)
    tibble(cell = f$def$cell, lab = p_to_stars(rr$p))
  })) %>% mutate(cell = factor(cell, levels = levels(dr$cell)))

  stars <- stars %>%
    left_join(dr %>% group_by(cell) %>%
                summarise(y = max(rate) * 0.95, .groups = "drop"), by = "cell")

  ggplot(dr, aes(x = Diagnosis, y = rate)) +
    geom_violin(aes(fill = Diagnosis), alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(aes(color = Diagnosis), width = 0.15, size = 1.4, alpha = 0.6) +
    geom_text(data = stars, aes(x = 1.5, y = y, label = lab), size = 7,
              inherit.aes = FALSE) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    facet_wrap(~ cell, scales = scales,
               labeller = as_labeller(c(beta = "\u03b2-cells", alpha = "\u03b1-cells",
                                        acinar = "Acinar cells"))) +
    labs(x = NULL, y = sprintf("%s+ rate (%% of cells)", which_marker),
         tag = tag, title = title,
         subtitle = subtitle) +   # what the points and stars show belongs in the legend
    theme_pub
}

`%||%` <- function(a, b) if (is.null(a)) b else a

plot_marker_d <- function(coll) {
  comp <- bind_rows(lapply(names(coll$raw), function(tag) {
    coll$raw[[tag]] %>%
      group_by(Diagnosis) %>%
      summarise(Beta = mean(Count.ins), Alpha = mean(Count.glu),
                .groups = "drop") %>%
      mutate(marker = tag)
  })) %>%
    pivot_longer(c(Beta, Alpha), names_to = "cell", values_to = "cells") %>%
    mutate(cell = factor(cell, levels = c("Alpha", "Beta")),
           marker = recode(marker, ki67 = "Ki-67 set", cd3 = "CD3 set"))

  totals <- comp %>% group_by(marker, Diagnosis) %>%
    summarise(total = sum(cells), .groups = "drop")

  ggplot(comp, aes(x = Diagnosis, y = cells, fill = cell)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_text(data = totals, aes(x = Diagnosis, y = total, label = round(total, 1)),
              vjust = -0.5, size = 3.2, inherit.aes = FALSE) +
    scale_fill_manual(values = c("Beta" = col_female, "Alpha" = col_male),
                      name = NULL) +
    facet_wrap(~ marker) +
    expand_limits(y = max(totals$total) * 1.12) +
    labs(x = NULL, y = "Endocrine cells per islet", tag = "D",
         title = "Each lost beta cell is matched by about one gained alpha cell",
         subtitle = "Total height is near-identical; composition inverts") +
    theme_pub +
    theme(legend.position = c(0.5, 0.88),
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# GENERATE
# =============================================================================

# Two panels: Ki-67 and CD3, each faceted by cell type.
#
# The forest plot and the beta/alpha substitution panel are no longer composed
# into this figure — the forest plot restated what these two panels already
# show, and the substitution arithmetic belongs with the total cell-count
# figure, where it is not artificially split by marker. plot_marker_a() and
# plot_marker_d() are still defined above if either is wanted for a
# supplementary panel.
generate_marker_main_figure <- function(coll = NULL, save = TRUE,
                                        output_dir = "Marker/Graphs",
                                        cells = c("beta", "alpha"),
                                        scales = "free_y") {
  if (is.null(coll)) { cat("Loading models and data...\n"); coll <- collect_marker_fits() }
  tab <- marker_rates_table(coll)

  cat("Generating marker main figure...\n")
  p_a <- plot_marker_bc(coll, "Ki-67", "A", "Ki-67 proliferation rate",
                        cells = cells, scales = scales); cat("  A (Ki-67)\n")
  p_b <- plot_marker_bc(coll, "CD3",   "B", "CD3 association rate",
                        cells = cells, scales = scales); cat("  B (CD3)\n")
  if (is.null(p_a) || is.null(p_b)) {
    stop("no models available for one of the markers — run marker_preflight()")
  }

  fig <- (p_a | p_b) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))

  if (save) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    for (ext in c("pdf", "png", "tiff")) {
      # cairo_pdf draws the Greek facet labels; the default pdf() cannot.
      ggsave(file.path(output_dir, paste0("marker_main_figure.", ext)), fig,
             width = 11, height = 5, dpi = 300,
             device = if (ext == "pdf") cairo_pdf else NULL,
             bg = if (ext != "pdf") "white" else NULL)
    }
    cat("Saved to", output_dir, "\n")
    write_csv(tab, file.path(output_dir, "marker_rate_ratios.csv"))
  }

  cat("\n=== rate ratios (plot_marker_main_figure.R:", MARKER_FIGURE_VERSION,
      ") ===\n", sep = "")
  cat("Computed by marginal_diagnosis() — exact linear contrast, no emmeans.\n")
  cat("Run validate_marker_rates() before quoting any of these.\n\n")
  print(as.data.frame(tab %>% dplyr::select(marker, cell, rate_ND, rate_T1D,
                                            ratio, lwr, upr, p,
                                            marginal_ratio, donors)),
        digits = 4)

  cat("\nFigure shows:", paste(cells, collapse = " + "),
      "— acinar is fitted and reported but not plotted here.\n")

  list(fig = fig, table = tab, panels = list(a = p_a, b = p_b), coll = coll)
}
