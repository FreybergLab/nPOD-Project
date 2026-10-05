# =============================================================================
# plot_main_figures.R
# Main Paper Figures from *_powered Models
# Updated: density-matching style — stars only, no annotation boxes
# =============================================================================
#
# For each 3-way interaction, creates TWO complementary views:
#
# DIAGNOSIS × REGION × ISLET SIZE:
#   View A: Facet by Region (3 panels), Diagnosis overlaid
#           → Shows ND vs T1D contrast within each Region
#   View B: Facet by Diagnosis (2 panels), Region overlaid
#           → Shows Region contrasts within each Diagnosis
#
# DIAGNOSIS × REGION × AGE:
#   View A: Facet by Region (3 panels), Diagnosis overlaid
#           → Shows ND vs T1D contrast within each Region
#   View B: Facet by Diagnosis (2 panels), Region overlaid
#           → Shows Region contrasts within each Diagnosis
#
# Annotations: significance stars only (pd > 0.975), no boxes or slope text
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
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank()
  )

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
    prob > 0.995 ~ "**",
    prob > 0.975 ~ "*",
    TRUE ~ ""
  )
}

# =============================================================================
# SLOPE COMPUTATION: DIAGNOSIS × REGION
# =============================================================================

compute_slopes_by_diag_region <- function(model, continuous_var) {
  
  if (continuous_var == "Age_c") {
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"),
      Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"),
      Age_c = c(-1, 1),
      log_Islet.Cells_c = 0
    )
  } else {
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"),
      Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"),
      Age_c = 0,
      log_Islet.Cells_c = c(-1, 1)
    )
  }
  
  post <- get_posterior_fitted(model, grid)
  
  slopes <- list()
  slope_draws <- list()
  
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      cell_slopes <- list()
      
      for (s in c("Female", "Male")) {
        if (continuous_var == "Age_c") {
          low_idx <- which(grid$Diagnosis == diag & grid$Region == reg & 
                             grid$Sex == s & grid$Age_c == -1)
          high_idx <- which(grid$Diagnosis == diag & grid$Region == reg & 
                              grid$Sex == s & grid$Age_c == 1)
        } else {
          low_idx <- which(grid$Diagnosis == diag & grid$Region == reg & 
                             grid$Sex == s & grid$log_Islet.Cells_c == -1)
          high_idx <- which(grid$Diagnosis == diag & grid$Region == reg & 
                              grid$Sex == s & grid$log_Islet.Cells_c == 1)
        }
        cell_slopes[[s]] <- (post[, high_idx] - post[, low_idx]) / 2
      }
      
      # Average over Sex
      draws <- (cell_slopes[["Female"]] + cell_slopes[["Male"]]) / 2
      
      summ <- summarize_posterior(draws)
      summ$Diagnosis <- diag
      summ$Region <- reg
      summ$pd <- pmax(mean(draws > 0), mean(draws < 0))
      summ$stars <- prob_to_stars(summ$pd)
      
      key <- paste(diag, reg, sep = "_")
      slopes[[key]] <- summ
      slope_draws[[key]] <- draws
    }
  }
  
  slopes_df <- bind_rows(slopes)
  
  return(list(slopes = slopes_df, draws = slope_draws))
}

# =============================================================================
# MODEL PREDICTIONS
# =============================================================================

get_model_predictions <- function(model, continuous_var, cont_range, n_points = 50, 
                                   by_vars = c("Diagnosis", "Region")) {
  
  cont_seq <- seq(cont_range[1], cont_range[2], length.out = n_points)
  
  if (continuous_var == "Age_c") {
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"),
      Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"),
      Age_c = cont_seq,
      log_Islet.Cells_c = 0
    )
  } else {
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"),
      Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"),
      Age_c = 0,
      log_Islet.Cells_c = cont_seq
    )
  }
  
  post <- get_posterior_fitted(model, grid)
  
  results <- grid %>%
    mutate(cell_id = row_number()) %>%
    group_by(across(all_of(c(by_vars, continuous_var)))) %>%
    summarize(cell_ids = list(cell_id), .groups = "drop") %>%
    mutate(
      post_summary = map(cell_ids, function(ids) {
        draws <- rowMeans(post[, ids, drop = FALSE])
        hpd <- coda::HPDinterval(coda::as.mcmc(draws), prob = 0.95)
        tibble(
          estimate = median(draws) * 100,
          lower = hpd[1] * 100,
          upper = hpd[2] * 100
        )
      })
    ) %>%
    unnest(post_summary) %>%
    dplyr::select(-cell_ids)
  
  return(results)
}

# =============================================================================
# VIEW A: FACET BY REGION, DIAGNOSIS OVERLAID (ISLET SIZE)
# Shows: T1D vs ND contrast within each Region
# =============================================================================

plot_islet_by_region <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  slopes_info <- compute_slopes_by_diag_region(model, "log_Islet.Cells_c")
  
  islet_mean <- mean(data$log_Islet.Cells)
  islet_range <- range(data$log_Islet.Cells_c)
  
  model_preds <- get_model_predictions(model, "log_Islet.Cells_c", islet_range, 
                                        n_points = 50, by_vars = c("Diagnosis", "Region"))
  model_preds <- model_preds %>%
    mutate(log_Islet.Cells = log_Islet.Cells_c + islet_mean)
  
  data <- data %>%
    mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region) %>%
    summarize(
      outcome = mean(outcome_pct, na.rm = TRUE),
      log_islet_mean = mean(log_Islet.Cells, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (is.null(y_max)) {
    y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.2
  }
  
  # Compute contrast stars (T1D - ND) within each Region for significance
  star_df <- tibble(Region = c("Head", "Body", "Tail")) %>%
    mutate(
      star_info = map(Region, function(r) {
        t1d_draws <- slopes_info$draws[[paste("T1D", r, sep = "_")]]
        nd_draws <- slopes_info$draws[[paste("ND", r, sep = "_")]]
        diff <- t1d_draws - nd_draws
        pd <- pmax(mean(diff > 0), mean(diff < 0))
        tibble(stars = prob_to_stars(pd))
      })
    ) %>%
    unnest(star_info) %>%
    filter(nchar(stars) > 0)
  
  x_breaks <- seq(2.5, 5.5, by = 0.5)
  x_labels <- round(exp(x_breaks))
  
  p <- ggplot() +
    geom_ribbon(
      data = model_preds,
      aes(x = log_Islet.Cells, ymin = lower, ymax = upper, fill = Diagnosis),
      alpha = 0.25
    ) +
    geom_line(
      data = model_preds,
      aes(x = log_Islet.Cells, y = estimate, color = Diagnosis),
      linewidth = 1.3
    ) +
    geom_point(
      data = donor_means,
      aes(x = log_islet_mean, y = outcome, color = Diagnosis),
      alpha = 0.5, size = 2
    ) +
    facet_wrap(~ Region, ncol = 3) +
    scale_color_manual(values = colors_diag) +
    scale_fill_manual(values = colors_diag) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(
      title = paste0(cell_label, ": Diagnosis \u00D7 Region \u00D7 Islet Size"),
      x = "Islet Size (cells)",
      y = paste0("% ", cell_label, " Cells"),
      color = "Diagnosis", fill = "Diagnosis"
    ) +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = mean(x_breaks), y = y_max * 0.95, label = stars),
      size = 8, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# VIEW B: FACET BY DIAGNOSIS, REGION OVERLAID (ISLET SIZE)
# Shows: Region contrasts within each Diagnosis
# =============================================================================

plot_islet_by_diagnosis <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  slopes_info <- compute_slopes_by_diag_region(model, "log_Islet.Cells_c")
  
  islet_mean <- mean(data$log_Islet.Cells)
  islet_range <- range(data$log_Islet.Cells_c)
  
  model_preds <- get_model_predictions(model, "log_Islet.Cells_c", islet_range, 
                                        n_points = 50, by_vars = c("Diagnosis", "Region"))
  model_preds <- model_preds %>%
    mutate(log_Islet.Cells = log_Islet.Cells_c + islet_mean)
  
  data <- data %>%
    mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region) %>%
    summarize(
      outcome = mean(outcome_pct, na.rm = TRUE),
      log_islet_mean = mean(log_Islet.Cells, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (is.null(y_max)) {
    y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.25
  }
  
  # Compute region contrast stars within each Diagnosis
  star_df <- tibble(Diagnosis = c("ND", "T1D")) %>%
    mutate(
      star_info = map(Diagnosis, function(d) {
        h_draws <- slopes_info$draws[[paste(d, "Head", sep = "_")]]
        b_draws <- slopes_info$draws[[paste(d, "Body", sep = "_")]]
        t_draws <- slopes_info$draws[[paste(d, "Tail", sep = "_")]]
        
        hb <- h_draws - b_draws
        ht <- h_draws - t_draws
        bt <- b_draws - t_draws
        
        tibble(
          comparison = c("Head - Body", "Head - Tail", "Body - Tail"),
          pd = c(pmax(mean(hb > 0), mean(hb < 0)),
                 pmax(mean(ht > 0), mean(ht < 0)),
                 pmax(mean(bt > 0), mean(bt < 0)))
        ) %>%
          mutate(stars = prob_to_stars(pd)) %>%
          filter(nchar(stars) > 0) %>%
          mutate(label = paste(comparison, stars))
      })
    ) %>%
    unnest(star_info) %>%
    group_by(Diagnosis) %>%
    summarize(combined_stars = paste(label, collapse = "\n"), .groups = "drop") %>%
    filter(nchar(combined_stars) > 0)
  
  x_breaks <- seq(2.5, 5.5, by = 0.5)
  x_labels <- round(exp(x_breaks))
  
  p <- ggplot() +
    geom_ribbon(
      data = model_preds,
      aes(x = log_Islet.Cells, ymin = lower, ymax = upper, fill = Region),
      alpha = 0.2
    ) +
    geom_line(
      data = model_preds,
      aes(x = log_Islet.Cells, y = estimate, color = Region),
      linewidth = 1.3
    ) +
    geom_point(
      data = donor_means,
      aes(x = log_islet_mean, y = outcome, color = Region),
      alpha = 0.5, size = 2
    ) +
    facet_wrap(~ Diagnosis, ncol = 2) +
    scale_color_manual(values = colors_region) +
    scale_fill_manual(values = colors_region) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(
      title = paste0(cell_label, ": Diagnosis \u00D7 Region \u00D7 Islet Size"),
      x = "Islet Size (cells)",
      y = paste0("% ", cell_label, " Cells"),
      color = "Region", fill = "Region"
    ) +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = mean(x_breaks), y = y_max * 0.95, label = combined_stars),
      size = 3.5, lineheight = 0.85, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# VIEW A: FACET BY REGION, DIAGNOSIS OVERLAID (AGE)
# =============================================================================

plot_age_by_region <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  slopes_info <- compute_slopes_by_diag_region(model, "Age_c")
  
  age_mean <- mean(data$Age)
  age_range <- range(data$Age_c)
  
  model_preds <- get_model_predictions(model, "Age_c", age_range, 
                                        n_points = 50, by_vars = c("Diagnosis", "Region"))
  model_preds <- model_preds %>%
    mutate(Age = Age_c + age_mean)
  
  data <- data %>%
    mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region, Age) %>%
    summarize(
      outcome = mean(outcome_pct, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (is.null(y_max)) {
    y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.2
  }
  
  # Contrast stars (T1D - ND) within each Region
  star_df <- tibble(Region = c("Head", "Body", "Tail")) %>%
    mutate(
      star_info = map(Region, function(r) {
        t1d_draws <- slopes_info$draws[[paste("T1D", r, sep = "_")]]
        nd_draws <- slopes_info$draws[[paste("ND", r, sep = "_")]]
        diff <- t1d_draws - nd_draws
        pd <- pmax(mean(diff > 0), mean(diff < 0))
        tibble(stars = prob_to_stars(pd))
      })
    ) %>%
    unnest(star_info) %>%
    filter(nchar(stars) > 0)
  
  p <- ggplot() +
    geom_ribbon(
      data = model_preds,
      aes(x = Age, ymin = lower, ymax = upper, fill = Diagnosis),
      alpha = 0.25
    ) +
    geom_line(
      data = model_preds,
      aes(x = Age, y = estimate, color = Diagnosis),
      linewidth = 1.3
    ) +
    geom_point(
      data = donor_means,
      aes(x = Age, y = outcome, color = Diagnosis),
      alpha = 0.5, size = 2
    ) +
    facet_wrap(~ Region, ncol = 3) +
    scale_color_manual(values = colors_diag) +
    scale_fill_manual(values = colors_diag) +
    scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(
      title = paste0(cell_label, ": Diagnosis \u00D7 Region \u00D7 Age"),
      x = "Age (years)",
      y = paste0("% ", cell_label, " Cells"),
      color = "Diagnosis", fill = "Diagnosis"
    ) +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = 21, y = y_max * 0.95, label = stars),
      size = 8, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# VIEW B: FACET BY DIAGNOSIS, REGION OVERLAID (AGE)
# =============================================================================

plot_age_by_diagnosis <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  slopes_info <- compute_slopes_by_diag_region(model, "Age_c")
  
  age_mean <- mean(data$Age)
  age_range <- range(data$Age_c)
  
  model_preds <- get_model_predictions(model, "Age_c", age_range, 
                                        n_points = 50, by_vars = c("Diagnosis", "Region"))
  model_preds <- model_preds %>%
    mutate(Age = Age_c + age_mean)
  
  data <- data %>%
    mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region, Age) %>%
    summarize(
      outcome = mean(outcome_pct, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (is.null(y_max)) {
    y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.25
  }
  
  # Region contrast stars within each Diagnosis
  star_df <- tibble(Diagnosis = c("ND", "T1D")) %>%
    mutate(
      star_info = map(Diagnosis, function(d) {
        h_draws <- slopes_info$draws[[paste(d, "Head", sep = "_")]]
        b_draws <- slopes_info$draws[[paste(d, "Body", sep = "_")]]
        t_draws <- slopes_info$draws[[paste(d, "Tail", sep = "_")]]
        
        hb <- h_draws - b_draws
        ht <- h_draws - t_draws
        bt <- b_draws - t_draws
        
        tibble(
          comparison = c("Head - Body", "Head - Tail", "Body - Tail"),
          pd = c(pmax(mean(hb > 0), mean(hb < 0)),
                 pmax(mean(ht > 0), mean(ht < 0)),
                 pmax(mean(bt > 0), mean(bt < 0)))
        ) %>%
          mutate(stars = prob_to_stars(pd)) %>%
          filter(nchar(stars) > 0) %>%
          mutate(label = paste(comparison, stars))
      })
    ) %>%
    unnest(star_info) %>%
    group_by(Diagnosis) %>%
    summarize(combined_stars = paste(label, collapse = "\n"), .groups = "drop") %>%
    filter(nchar(combined_stars) > 0)
  
  p <- ggplot() +
    geom_ribbon(
      data = model_preds,
      aes(x = Age, ymin = lower, ymax = upper, fill = Region),
      alpha = 0.2
    ) +
    geom_line(
      data = model_preds,
      aes(x = Age, y = estimate, color = Region),
      linewidth = 1.3
    ) +
    geom_point(
      data = donor_means,
      aes(x = Age, y = outcome, color = Region),
      alpha = 0.5, size = 2
    ) +
    facet_wrap(~ Diagnosis, ncol = 2) +
    scale_color_manual(values = colors_region) +
    scale_fill_manual(values = colors_region) +
    scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(
      title = paste0(cell_label, ": Diagnosis \u00D7 Region \u00D7 Age"),
      x = "Age (years)",
      y = paste0("% ", cell_label, " Cells"),
      color = "Region", fill = "Region"
    ) +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = 21, y = y_max * 0.95, label = combined_stars),
      size = 3.5, lineheight = 0.85, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# MASTER FUNCTION
# =============================================================================

generate_main_3way_figures <- function(data, model_dir = "Proportion/Models",
                                       output_dir = "Proportion/Graphs/Main", save_plots = TRUE) {
  
  cat("\n=============================================================================\n")
  cat("GENERATING MAIN PAPER 3-WAY INTERACTION FIGURES\n")
  cat("Two views per interaction: contrast within facet variable\n")
  cat("=============================================================================\n")
  
  cell_specs <- list(
    beta  = list(file = "ins_powered_model.rds",  var = "Percent.ins",  label = "Beta Cell",  y_max = 100),
    alpha = list(file = "glu_powered_model.rds",  var = "Percent.glu",  label = "Alpha Cell", y_max = 100),
    delta = list(file = "soma_powered_model.rds", var = "Percent.soma", label = "Delta Cell", y_max = 60),
    pp    = list(file = "pp_powered_model.rds",   var = "Percent.PP",   label = "PP Cell",    y_max = 90)
  )
  
  if (save_plots) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  results <- list()
  
  for (cell_type in names(cell_specs)) {
    spec <- cell_specs[[cell_type]]
    
    cat("\n--- Processing", spec$label, "---\n")
    model_path <- file.path(model_dir, spec$file)
    cat("  Loading", model_path, "...\n")
    model <- readRDS(model_path)
    
    # Islet Size: Two views
    cat("  Creating Islet Size figures (2 views)...\n")
    p_islet_A <- plot_islet_by_region(data, model, spec$var, spec$label, spec$y_max)
    p_islet_B <- plot_islet_by_diagnosis(data, model, spec$var, spec$label, spec$y_max)
    
    # Age: Two views
    cat("  Creating Age figures (2 views)...\n")
    p_age_A <- plot_age_by_region(data, model, spec$var, spec$label, spec$y_max)
    p_age_B <- plot_age_by_diagnosis(data, model, spec$var, spec$label, spec$y_max)
    
    rm(model); gc()
    
    results[[cell_type]] <- list(
      islet_by_region = p_islet_A,
      islet_by_diagnosis = p_islet_B,
      age_by_region = p_age_A,
      age_by_diagnosis = p_age_B
    )
    
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(cell_type, "_islet_by_region.png")),
             p_islet_A, width = 14, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_islet_by_diagnosis.png")),
             p_islet_B, width = 12, height = 6, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_age_by_region.png")),
             p_age_A, width = 14, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_age_by_diagnosis.png")),
             p_age_B, width = 12, height = 6, dpi = 300, bg = "white")
      cat("  Saved 4 figures to", output_dir, "\n")
    }
  }
  
  cat("\n=============================================================================\n")
  cat("COMPLETE! Generated", length(cell_specs) * 4, "figures.\n")
  cat("=============================================================================\n")
  
  return(results)
}

if (interactive()) {
  cat("\nTo generate all main paper 3-way figures:\n")
  cat("  figures <- generate_main_3way_figures(quad_data)\n")
  cat("  # Models loaded individually from Proportion/Models/\n")
}
