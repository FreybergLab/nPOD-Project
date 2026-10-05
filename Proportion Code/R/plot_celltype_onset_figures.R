# =============================================================================
# plot_celltype_onset_figures.R
# T1D-only onset/duration analysis — cell-type composition
# Style: matches density main figure panels B & C
# =============================================================================
#
# Two panels per cell type:
#   Panel A: % Cell vs Age at Onset (marginalized over Region, Sex, Islet Size)
#            Points colored by disease duration bin
#   Panel B: % Cell vs Disease Duration (marginalized over Region, Sex, Islet Size)
#            Points colored by onset age bin
#
# Combined figure: 4 rows × 2 columns (one row per cell type)
#
# =============================================================================

library(tidyverse)
library(brms)
library(coda)
library(patchwork)
library(cowplot)

source("colors_master.R")

# =============================================================================
# THEME (matching density figure style)
# =============================================================================

theme_pub <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.tag = element_text(size = 20, face = "bold"),
    plot.margin = margin(8, 12, 8, 8),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank()
  )

# --- Onset age bins (for coloring duration panel) ---
colors_onset_bin <- c("0\u20135 yr" = "#FDAE61",
                      "5\u201315 yr" = "#F46D43",
                      ">15 yr"      = "#A50026")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

get_posterior_fitted <- function(model, newdata, ndraws = NULL) {
  fitted(model, newdata = newdata, re_formula = NA,
         summary = FALSE, ndraws = ndraws)
}

summarize_posterior <- function(x, prob = 0.95) {
  hpd <- coda::HPDinterval(coda::as.mcmc(x), prob = prob)
  tibble(
    estimate = median(x),
    lower.HPD = hpd[1],
    upper.HPD = hpd[2]
  )
}

prob_to_stars <- function(prob) {
  case_when(
    prob > 0.9995 ~ "***",
    prob > 0.995  ~ "**",
    prob > 0.975  ~ "*",
    TRUE ~ ""
  )
}

# =============================================================================
# MARGINALIZED PREDICTIONS
# =============================================================================
#
# For a focal continuous variable (age_at_onset_c or Disease.Duration_c),
# compute predictions marginalized over Region, Sex, and the other continuous
# covariates (set to 0 = their mean since centered).
#

get_onset_predictions <- function(model, data, focal_var, focal_range, n_points = 50) {
  
  focal_seq <- seq(focal_range[1], focal_range[2], length.out = n_points)
  
  # Expand over Region and Sex for marginalization
  grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    log_Islet.Cells_c = 0
  )
  
  all_preds <- list()
  
  for (i in seq_len(nrow(grid))) {
    nd <- data.frame(
      age_at_onset_c = if (focal_var == "age_at_onset_c") focal_seq else 0,
      Disease.Duration_c = if (focal_var == "Disease.Duration_c") focal_seq else 0,
      Region = factor(grid$Region[i], levels = levels(data$Region)),
      Sex = factor(grid$Sex[i], levels = levels(data$Sex)),
      log_Islet.Cells_c = 0,
      Donor = NA,
      ImageID = NA
    )
    
    post <- get_posterior_fitted(model, nd)
    all_preds[[i]] <- post
  }
  
  # Average across Region × Sex combinations (marginalize)
  n_draws <- nrow(all_preds[[1]])
  n_pts <- ncol(all_preds[[1]])
  avg_post <- matrix(0, nrow = n_draws, ncol = n_pts)
  for (i in seq_along(all_preds)) {
    avg_post <- avg_post + all_preds[[i]]
  }
  avg_post <- avg_post / length(all_preds)
  
  # Summarize to median + HPD
  results <- tibble(focal_val = focal_seq) %>%
    mutate(
      post_summary = map(seq_len(n_pts), function(j) {
        draws <- avg_post[, j]
        hpd <- coda::HPDinterval(coda::as.mcmc(draws), prob = 0.95)
        tibble(
          estimate = median(draws) * 100,
          lower = hpd[1] * 100,
          upper = hpd[2] * 100
        )
      })
    ) %>%
    unnest(post_summary)
  
  return(results)
}

# =============================================================================
# SLOPE SIGNIFICANCE (posterior probability of direction)
# =============================================================================

get_slope_pd <- function(model, focal_var, data) {
  
  # Two-point finite difference marginalized over Region × Sex
  grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    log_Islet.Cells_c = 0
  )
  
  all_slopes <- list()
  
  for (i in seq_len(nrow(grid))) {
    nd_low <- data.frame(
      age_at_onset_c = if (focal_var == "age_at_onset_c") -1 else 0,
      Disease.Duration_c = if (focal_var == "Disease.Duration_c") -1 else 0,
      Region = factor(grid$Region[i], levels = levels(data$Region)),
      Sex = factor(grid$Sex[i], levels = levels(data$Sex)),
      log_Islet.Cells_c = 0,
      Donor = NA, ImageID = NA
    )
    nd_high <- nd_low
    if (focal_var == "age_at_onset_c") {
      nd_high$age_at_onset_c <- 1
    } else {
      nd_high$Disease.Duration_c <- 1
    }
    
    post_low <- get_posterior_fitted(model, nd_low)
    post_high <- get_posterior_fitted(model, nd_high)
    all_slopes[[i]] <- (post_high[, 1] - post_low[, 1]) / 2
  }
  
  # Average over Region × Sex
  avg_slope <- Reduce("+", all_slopes) / length(all_slopes)
  pd <- pmax(mean(avg_slope > 0), mean(avg_slope < 0))
  
  return(list(pd = pd, stars = prob_to_stars(pd), draws = avg_slope))
}

# =============================================================================
# PANEL A: % Cell vs Age at Onset
# =============================================================================

plot_onset_panel <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  onset_mean <- mean(data$age_at_onset, na.rm = TRUE)
  onset_range <- range(data$age_at_onset_c, na.rm = TRUE)
  
  # Model predictions
  preds <- get_onset_predictions(model, data, "age_at_onset_c", onset_range)
  preds <- preds %>% mutate(age_at_onset = focal_val + onset_mean)
  
  # Donor-level means
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  donor_means <- data %>%
    group_by(Donor, age_at_onset, Disease.Duration) %>%
    summarise(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(duration_bin = cut(Disease.Duration, breaks = c(0, 5, 15, Inf),
                              labels = c("0\u20135 yr", "5\u201315 yr", ">15 yr")))
  
  if (is.null(y_max)) y_max <- max(c(donor_means$outcome, preds$upper), na.rm = TRUE) * 1.1
  
  # Significance
  sig <- get_slope_pd(model, "age_at_onset_c", data)
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = age_at_onset, ymin = lower, ymax = upper),
                fill = "gray50", alpha = 0.2) +
    geom_line(data = preds, aes(x = age_at_onset, y = estimate),
              color = "black", linewidth = 0.8) +
    geom_point(data = donor_means, aes(x = age_at_onset, y = outcome),
               size = 2, alpha = 0.7) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Age at onset (years)",
         y = paste0("% ", cell_label, " Cells\n(donor mean)"),
         title = paste0(cell_label, " Cell: onset age")) +
    theme_pub
  
  if (nchar(sig$stars) > 0) {
    plt <- plt + annotate("text",
      x = mean(range(donor_means$age_at_onset, na.rm = TRUE)),
      y = 95,
      label = sig$stars, size = 8)
  }
  
  return(plt)
}

# =============================================================================
# PANEL B: % Cell vs Disease Duration
# =============================================================================

plot_duration_panel <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  dur_mean <- mean(data$Disease.Duration, na.rm = TRUE)
  dur_range <- range(data$Disease.Duration_c, na.rm = TRUE)
  
  # Model predictions
  preds <- get_onset_predictions(model, data, "Disease.Duration_c", dur_range)
  preds <- preds %>% mutate(Disease.Duration = focal_val + dur_mean)
  
  # Donor-level means
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  donor_means <- data %>%
    group_by(Donor, Disease.Duration, age_at_onset) %>%
    summarise(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
    mutate(onset_bin = cut(age_at_onset, breaks = c(-Inf, 5, 15, Inf),
                           labels = c("0\u20135 yr", "5\u201315 yr", ">15 yr")))
  
  if (is.null(y_max)) y_max <- max(c(donor_means$outcome, preds$upper), na.rm = TRUE) * 1.1
  
  # Significance
  sig <- get_slope_pd(model, "Disease.Duration_c", data)
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Disease.Duration, ymin = lower, ymax = upper),
                fill = "gray50", alpha = 0.2) +
    geom_line(data = preds, aes(x = Disease.Duration, y = estimate),
              color = "black", linewidth = 0.8) +
    geom_point(data = donor_means, aes(x = Disease.Duration, y = outcome),
               size = 2, alpha = 0.7) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = "Disease duration (years)",
         y = paste0("% ", cell_label, " Cells\n(donor mean)"),
         title = paste0(cell_label, " Cell: disease duration")) +
    theme_pub
  
  if (nchar(sig$stars) > 0) {
    plt <- plt + annotate("text",
      x = mean(range(donor_means$Disease.Duration, na.rm = TRUE)),
      y = 95,
      label = sig$stars, size = 8)
  }
  
  return(plt)
}

# =============================================================================
# MASTER FUNCTION
# =============================================================================

generate_onset_figures <- function(t1d_data, models, output_dir = "Proportion/Graphs/Onset",
                                    save_plots = TRUE) {
  
  cat("\n=============================================================================\n")
  cat("GENERATING T1D ONSET / DURATION FIGURES\n")
  cat("=============================================================================\n")
  
  # models should be a named list: list(beta = ins_onset, alpha = glu_onset, ...)
  cell_specs <- list(
    beta  = list(var = "Percent.ins",  label = "\u03b2"),
    alpha = list(var = "Percent.glu",  label = "\u03b1"),
    delta = list(var = "Percent.soma", label = "\u03b4"),
    pp    = list(var = "Percent.PP",   label = "PP")
  )
  
  if (save_plots) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  results <- list()
  onset_panels <- list()
  duration_panels <- list()
  
  for (cell_type in names(cell_specs)) {
    if (!cell_type %in% names(models)) {
      cat("  Skipping", cell_type, "— no model provided\n")
      next
    }
    
    spec <- cell_specs[[cell_type]]
    model <- models[[cell_type]]
    
    cat("\n--- Processing", spec$label, "---\n")
    
    cat("  Creating onset panel...\n")
    p_onset <- plot_onset_panel(t1d_data, model, spec$var, spec$label)
    
    cat("  Creating duration panel...\n")
    p_duration <- plot_duration_panel(t1d_data, model, spec$var, spec$label)
    
    onset_panels[[cell_type]] <- p_onset
    duration_panels[[cell_type]] <- p_duration
    
    results[[cell_type]] <- list(onset = p_onset, duration = p_duration)
    
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(cell_type, "_onset.png")),
             p_onset, width = 6, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_duration.png")),
             p_duration, width = 6, height = 5, dpi = 300, bg = "white")
      
      # Paired onset + duration side by side
      p_paired <- p_onset + p_duration +
        plot_layout(ncol = 2) +
        plot_annotation(
          title = paste0("T1D ", spec$label, " Cell: Onset Age & Disease Duration"),
          theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
        )
      ggsave(file.path(output_dir, paste0(cell_type, "_onset_duration_paired.png")),
             p_paired, width = 12, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_onset_duration_paired.pdf")),
             p_paired, width = 12, height = 5, bg = "white")
      
      results[[cell_type]]$paired <- p_paired
      cat("  Saved 4 figures\n")
    }
  }
  
  # Combined figure: 4 rows × 2 columns
  cat("\nAssembling combined figure...\n")
  
  # Remove legends from individual panels for combined view, add shared legend
  combined <- wrap_plots(
    onset_panels$beta + theme(legend.position = "none") +
      labs(tag = "A", title = "\u03b2-Cell: onset age"),
    duration_panels$beta + theme(legend.position = "none") +
      labs(title = "\u03b2-Cell: disease duration"),
    onset_panels$alpha + theme(legend.position = "none") +
      labs(tag = "B", title = "\u03b1-Cell: onset age"),
    duration_panels$alpha + theme(legend.position = "none") +
      labs(title = "\u03b1-Cell: disease duration"),
    onset_panels$delta + theme(legend.position = "none") +
      labs(tag = "C", title = "\u03b4-Cell: onset age"),
    duration_panels$delta + theme(legend.position = "none") +
      labs(title = "\u03b4-Cell: disease duration"),
    onset_panels$pp + theme(legend.position = "none") +
      labs(tag = "D", title = "PP-Cell: onset age"),
    duration_panels$pp + theme(legend.position = "none") +
      labs(title = "PP-Cell: disease duration"),
    ncol = 2
  )
  
  results$combined <- combined
  
  if (save_plots) {
    ggsave(file.path(output_dir, "all_onset_duration_combined.png"),
           combined, width = 12, height = 16, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "all_onset_duration_combined.pdf"),
           combined, width = 12, height = 16, bg = "white")
    cat("Saved combined figure\n")
  }
  
  cat("\n=============================================================================\n")
  cat("COMPLETE!\n")
  cat("=============================================================================\n")
  
  return(results)
}

# =============================================================================
# USAGE
# =============================================================================

if (interactive()) {
  cat("\nTo generate onset/duration figures:\n")
  cat("  # Load T1D-only proportion models individually\n")
  cat("  onset_models <- list(\n")
  cat("    beta  = readRDS('Proportion/Models/ins_ao_model.rds'),\n")
  cat("    alpha = readRDS('Proportion/Models/glu_ao_model.rds'),\n")
  cat("    delta = readRDS('Proportion/Models/soma_ao_model.rds'),\n")
  cat("    pp    = readRDS('Proportion/Models/pp_ao_model.rds')\n")
  cat("  )\n")
  cat("  \n")
  cat("  # Prepare T1D data\n")
  cat("  t1d_data <- readRDS('Data/quad_data_t1d.rds') %>%\n")
  cat("    mutate(Region = factor(Region, levels = c('Head','Body','Tail')),\n")
  cat("           Sex = factor(Sex, levels = c('Female','Male')),\n")
  cat("           Donor = factor(Donor), ImageID = factor(ImageID)) %>%\n")
  cat("    droplevels()\n")
  cat("  \n")
  cat("  figures <- generate_onset_figures(t1d_data, onset_models)\n")
}
