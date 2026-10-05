# =============================================================================
# plot_supplementary_figures.R
# Supplementary Figures from *_mid Models
# Updated: density-matching style — stars only, no annotation boxes
# =============================================================================
#
# For each interaction, creates complementary views showing different contrasts:
#
# DIAGNOSIS × REGION × SEX (categorical):
#   View A: Facet by Region → Shows Diagnosis contrasts
#   View B: Facet by Diagnosis → Shows Region contrasts
#
# DIAGNOSIS × SEX × AGE/ISLET:
#   View A: Facet by Sex, Diagnosis overlaid → Diag contrasts within Sex
#   View B: Facet by Diagnosis, Sex overlaid → Sex contrasts within Diag
#
# REGION × SEX × AGE/ISLET:
#   View A: Facet by Sex, Region overlaid → Region contrasts within Sex
#   View B: Facet by Region, Sex overlaid → Sex contrasts within Region
#
# 4-WAY: 3×2 grid (Region × Sex), Diagnosis overlaid
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
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.margin = margin(8, 12, 8, 8),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )

# =============================================================================
# FACTOR LEVELS
# =============================================================================

region_levels <- c("Head", "Body", "Tail")

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
# MODEL PREDICTIONS (GENERAL)
# =============================================================================

get_model_predictions_general <- function(model, continuous_var, cont_range, n_points = 50,
                                           by_vars = c("Diagnosis", "Sex")) {
  
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
  
  if ("Region" %in% names(results)) {
    results$Region <- factor(results$Region, levels = region_levels)
  }
  
  return(results)
}

# =============================================================================
# SLOPE COMPUTATION (GENERAL)
# =============================================================================

compute_slopes_general <- function(model, continuous_var, by_vars) {
  
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
  
  combos <- unique(grid[, by_vars, drop = FALSE])
  
  slopes <- list()
  slope_draws <- list()
  
  for (i in 1:nrow(combos)) {
    combo <- combos[i, ]
    
    if (continuous_var == "Age_c") {
      match_low <- grid$Age_c == -1
      match_high <- grid$Age_c == 1
    } else {
      match_low <- grid$log_Islet.Cells_c == -1
      match_high <- grid$log_Islet.Cells_c == 1
    }
    
    for (v in by_vars) {
      match_low <- match_low & grid[[v]] == combo[[v]]
      match_high <- match_high & grid[[v]] == combo[[v]]
    }
    
    low_idx <- which(match_low)
    high_idx <- which(match_high)
    
    draws <- (rowMeans(post[, high_idx, drop = FALSE]) - 
                rowMeans(post[, low_idx, drop = FALSE])) / 2
    
    summ <- summarize_posterior(draws)
    for (v in by_vars) summ[[v]] <- combo[[v]]
    summ$pd <- pmax(mean(draws > 0), mean(draws < 0))
    summ$stars <- prob_to_stars(summ$pd)
    
    # Create key in consistent order based on by_vars
    key_parts <- sapply(by_vars, function(v) as.character(combo[[v]]))
    key <- paste(key_parts, collapse = "_")
    slopes[[key]] <- summ
    slope_draws[[key]] <- draws
  }
  
  slopes_df <- bind_rows(slopes)
  
  return(list(slopes = slopes_df, draws = slope_draws))
}

# =============================================================================
# DIAGNOSIS × REGION × SEX (CATEGORICAL)
# View A: Facet by Region (shows Diag contrasts within Region × Sex)
# =============================================================================

plot_diag_region_sex_A <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region, Sex) %>%
    summarize(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  
  if (is.null(y_max)) y_max <- max(donor_means$outcome, na.rm = TRUE) * 1.15
  
  cat_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), Age_c = 0, log_Islet.Cells_c = 0
  )
  post <- get_posterior_fitted(model, cat_grid)
  
  # Compute stars for T1D-ND contrast within each Region × Sex
  star_df <- expand.grid(Region = c("Head", "Body", "Tail"), 
                          Sex = c("Female", "Male"), stringsAsFactors = FALSE) %>%
    mutate(Region = factor(Region, levels = region_levels)) %>%
    mutate(
      star_info = map2(Region, Sex, function(r, s) {
        nd_idx <- which(cat_grid$Diagnosis == "ND" & cat_grid$Region == r & cat_grid$Sex == s)
        t1d_idx <- which(cat_grid$Diagnosis == "T1D" & cat_grid$Region == r & cat_grid$Sex == s)
        diff <- post[, t1d_idx] - post[, nd_idx]
        pd <- pmax(mean(diff > 0), mean(diff < 0))
        tibble(stars = prob_to_stars(pd))
      })
    ) %>%
    unnest(star_info) %>%
    filter(nchar(stars) > 0)
  
  p <- ggplot(donor_means, aes(x = Sex, y = outcome, fill = Diagnosis)) +
    geom_violin(position = position_dodge(0.9), alpha = 0.3, color = NA, scale = "width") +
    geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.9),
               size = 1.2, alpha = 0.5, aes(color = Diagnosis)) +
    facet_wrap(~ Region, ncol = 3) +
    scale_fill_manual(values = colors_diag) +
    scale_color_manual(values = colors_diag, guide = "none") +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = paste0(cell_label, ": Diagnosis \u00D7 Region \u00D7 Sex"),
         x = "", y = paste0("% ", cell_label, " Cells"), fill = "Diagnosis") +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = Sex, y = y_max * 0.95, label = stars),
      inherit.aes = FALSE, size = 6
    )
  }
  
  return(p)
}

# =============================================================================
# DIAGNOSIS × REGION × SEX (CATEGORICAL)
# View B: Facet by Diagnosis (shows Region contrasts within Diag × Sex)
# =============================================================================

plot_diag_region_sex_B <- function(data, model, outcome_var, cell_label, y_max = NULL) {
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region, Sex) %>%
    summarize(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  
  if (is.null(y_max)) y_max <- max(donor_means$outcome, na.rm = TRUE) * 1.2
  
  cat_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), Age_c = 0, log_Islet.Cells_c = 0
  )
  post <- get_posterior_fitted(model, cat_grid)
  
  # Compute region contrast stars within each Diagnosis × Sex
  star_df <- expand.grid(Diagnosis = c("ND", "T1D"), 
                          Sex = c("Female", "Male"), stringsAsFactors = FALSE) %>%
    mutate(
      star_info = map2(Diagnosis, Sex, function(d, s) {
        h_idx <- which(cat_grid$Diagnosis == d & cat_grid$Region == "Head" & cat_grid$Sex == s)
        b_idx <- which(cat_grid$Diagnosis == d & cat_grid$Region == "Body" & cat_grid$Sex == s)
        t_idx <- which(cat_grid$Diagnosis == d & cat_grid$Region == "Tail" & cat_grid$Sex == s)
        
        hb <- (post[, h_idx] - post[, b_idx])
        ht <- (post[, h_idx] - post[, t_idx])
        
        tibble(
          comparison = c("H-B", "H-T"),
          pd = c(pmax(mean(hb > 0), mean(hb < 0)),
                 pmax(mean(ht > 0), mean(ht < 0)))
        ) %>%
          mutate(stars = prob_to_stars(pd)) %>%
          filter(nchar(stars) > 0) %>%
          mutate(label = paste(comparison, stars))
      })
    ) %>%
    unnest(star_info) %>%
    group_by(Diagnosis, Sex) %>%
    summarize(combined_stars = paste(label, collapse = "\n"), .groups = "drop") %>%
    filter(nchar(combined_stars) > 0)
  
  p <- ggplot(donor_means, aes(x = Sex, y = outcome, fill = Region)) +
    geom_violin(position = position_dodge(0.9), alpha = 0.3, color = NA, scale = "width") +
    geom_point(position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.9),
               size = 1.2, alpha = 0.5, aes(color = Region)) +
    facet_wrap(~ Diagnosis, ncol = 2) +
    scale_fill_manual(values = colors_region) +
    scale_color_manual(values = colors_region, guide = "none") +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = paste0(cell_label, ": Diagnosis \u00D7 Region \u00D7 Sex"),
         x = "", y = paste0("% ", cell_label, " Cells"), fill = "Region") +
    theme_pub +
    theme(legend.position = "bottom")
  
  return(p)
}

# =============================================================================
# DIAGNOSIS × SEX × CONTINUOUS
# View A: Facet by Sex, Diagnosis overlaid (shows Diag contrast)
# =============================================================================

plot_diag_sex_cont_A <- function(data, model, outcome_var, cell_label, continuous_var, y_max = NULL) {
  
  is_age <- continuous_var == "Age_c"
  
  if (is_age) {
    cont_mean <- mean(data$Age)
    cont_range <- range(data$Age_c)
  } else {
    cont_mean <- mean(data$log_Islet.Cells)
    cont_range <- range(data$log_Islet.Cells_c)
  }
  
  model_preds <- get_model_predictions_general(model, continuous_var, cont_range,
                                                by_vars = c("Diagnosis", "Sex"))
  if (is_age) {
    model_preds <- model_preds %>% mutate(cont_val = Age_c + cont_mean)
  } else {
    model_preds <- model_preds %>% mutate(cont_val = log_Islet.Cells_c + cont_mean)
  }
  
  slopes_info <- compute_slopes_general(model, continuous_var, c("Diagnosis", "Sex"))
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  if (is_age) {
    donor_means <- data %>%
      group_by(Donor, Diagnosis, Sex, Age) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
      rename(cont_val = Age)
  } else {
    donor_means <- data %>%
      group_by(Donor, Diagnosis, Sex) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE),
                cont_val = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  
  if (is.null(y_max)) y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.2
  
  # Contrast stars (T1D - ND) within each Sex
  star_df <- tibble(Sex = c("Female", "Male")) %>%
    mutate(
      star_info = map(Sex, function(s) {
        t1d_draws <- slopes_info$draws[[paste("T1D", s, sep = "_")]]
        nd_draws <- slopes_info$draws[[paste("ND", s, sep = "_")]]
        diff <- t1d_draws - nd_draws
        pd <- pmax(mean(diff > 0), mean(diff < 0))
        tibble(stars = prob_to_stars(pd))
      })
    ) %>%
    unnest(star_info) %>%
    filter(nchar(stars) > 0)
  
  p <- ggplot() +
    geom_ribbon(data = model_preds, aes(x = cont_val, ymin = lower, ymax = upper, fill = Diagnosis), alpha = 0.25) +
    geom_line(data = model_preds, aes(x = cont_val, y = estimate, color = Diagnosis), linewidth = 1.3) +
    facet_wrap(~ Sex, ncol = 2) +
    scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
    scale_y_continuous(limits = c(0, y_max))
  
  # Only show donor mean points for age plots (islet size means cluster at same x)
  if (is_age) {
    p <- p + geom_point(data = donor_means, aes(x = cont_val, y = outcome, color = Diagnosis), alpha = 0.5, size = 2)
  }

  if (is_age) {
    p <- p + scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) + labs(x = "Age (years)")
    ann_x <- 21
  } else {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks))) + labs(x = "Islet Size (cells)")
    ann_x <- 4.0
  }
  
  p <- p +
    labs(title = paste0(cell_label, ": Diagnosis \u00D7 Sex \u00D7 ", ifelse(is_age, "Age", "Islet Size")),
         y = paste0("% ", cell_label, " Cells"), color = "Diagnosis", fill = "Diagnosis") +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = ann_x, y = y_max * 0.95, label = stars),
      size = 8, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# DIAGNOSIS × SEX × CONTINUOUS
# View B: Facet by Diagnosis, Sex overlaid (shows Sex contrast)
# =============================================================================

plot_diag_sex_cont_B <- function(data, model, outcome_var, cell_label, continuous_var, y_max = NULL) {
  
  is_age <- continuous_var == "Age_c"
  
  if (is_age) {
    cont_mean <- mean(data$Age)
    cont_range <- range(data$Age_c)
  } else {
    cont_mean <- mean(data$log_Islet.Cells)
    cont_range <- range(data$log_Islet.Cells_c)
  }
  
  model_preds <- get_model_predictions_general(model, continuous_var, cont_range,
                                                by_vars = c("Diagnosis", "Sex"))
  if (is_age) {
    model_preds <- model_preds %>% mutate(cont_val = Age_c + cont_mean)
  } else {
    model_preds <- model_preds %>% mutate(cont_val = log_Islet.Cells_c + cont_mean)
  }
  
  slopes_info <- compute_slopes_general(model, continuous_var, c("Diagnosis", "Sex"))
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  if (is_age) {
    donor_means <- data %>%
      group_by(Donor, Diagnosis, Sex, Age) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
      rename(cont_val = Age)
  } else {
    donor_means <- data %>%
      group_by(Donor, Diagnosis, Sex) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE),
                cont_val = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  
  if (is.null(y_max)) y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.2
  
  # Contrast stars (F - M) within each Diagnosis
  star_df <- tibble(Diagnosis = c("ND", "T1D")) %>%
    mutate(
      star_info = map(Diagnosis, function(d) {
        f_draws <- slopes_info$draws[[paste(d, "Female", sep = "_")]]
        m_draws <- slopes_info$draws[[paste(d, "Male", sep = "_")]]
        diff <- f_draws - m_draws
        pd <- pmax(mean(diff > 0), mean(diff < 0))
        tibble(stars = prob_to_stars(pd))
      })
    ) %>%
    unnest(star_info) %>%
    filter(nchar(stars) > 0)
  
  p <- ggplot() +
    geom_ribbon(data = model_preds, aes(x = cont_val, ymin = lower, ymax = upper, fill = Sex), alpha = 0.25) +
    geom_line(data = model_preds, aes(x = cont_val, y = estimate, color = Sex), linewidth = 1.3) +
    facet_wrap(~ Diagnosis, ncol = 2) +
    scale_color_manual(values = colors_sex) + scale_fill_manual(values = colors_sex) +
    scale_y_continuous(limits = c(0, y_max))
  
  # Only show donor mean points for age plots (islet size means cluster at same x)
  if (is_age) {
    p <- p + geom_point(data = donor_means, aes(x = cont_val, y = outcome, color = Sex), alpha = 0.5, size = 2)
  }

  if (is_age) {
    p <- p + scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) + labs(x = "Age (years)")
    ann_x <- 21
  } else {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks))) + labs(x = "Islet Size (cells)")
    ann_x <- 4.0
  }
  
  p <- p +
    labs(title = paste0(cell_label, ": Diagnosis \u00D7 Sex \u00D7 ", ifelse(is_age, "Age", "Islet Size")),
         y = paste0("% ", cell_label, " Cells"), color = "Sex", fill = "Sex") +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = ann_x, y = y_max * 0.95, label = stars),
      size = 8, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# REGION × SEX × CONTINUOUS
# View A: Facet by Sex, Region overlaid (shows Region contrast)
# =============================================================================

plot_region_sex_cont_A <- function(data, model, outcome_var, cell_label, continuous_var, y_max = NULL) {
  
  is_age <- continuous_var == "Age_c"
  
  if (is_age) {
    cont_mean <- mean(data$Age)
    cont_range <- range(data$Age_c)
  } else {
    cont_mean <- mean(data$log_Islet.Cells)
    cont_range <- range(data$log_Islet.Cells_c)
  }
  
  model_preds <- get_model_predictions_general(model, continuous_var, cont_range,
                                                by_vars = c("Region", "Sex"))
  if (is_age) {
    model_preds <- model_preds %>% mutate(cont_val = Age_c + cont_mean)
  } else {
    model_preds <- model_preds %>% mutate(cont_val = log_Islet.Cells_c + cont_mean)
  }
  
  slopes_info <- compute_slopes_general(model, continuous_var, c("Region", "Sex"))
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  if (is_age) {
    donor_means <- data %>%
      group_by(Donor, Region, Sex, Age) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
      rename(cont_val = Age)
  } else {
    donor_means <- data %>%
      group_by(Donor, Region, Sex) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE),
                cont_val = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  
  if (is.null(y_max)) y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.25
  
  # Contrast stars: Head-Tail within each Sex
  star_df <- tibble(Sex = c("Female", "Male")) %>%
    mutate(
      star_info = map(Sex, function(s) {
        h_draws <- slopes_info$draws[[paste("Head", s, sep = "_")]]
        b_draws <- slopes_info$draws[[paste("Body", s, sep = "_")]]
        t_draws <- slopes_info$draws[[paste("Tail", s, sep = "_")]]
        
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
    group_by(Sex) %>%
    summarize(combined_stars = paste(label, collapse = "\n"), .groups = "drop") %>%
    filter(nchar(combined_stars) > 0)
  
  if (is_age) {
    ann_x <- 21
  } else {
    ann_x <- 4.0
  }
  
  p <- ggplot() +
    geom_ribbon(data = model_preds, aes(x = cont_val, ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
    geom_line(data = model_preds, aes(x = cont_val, y = estimate, color = Region), linewidth = 1.3) +
    facet_wrap(~ Sex, ncol = 2) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    scale_y_continuous(limits = c(0, y_max))
  
  # Only show donor mean points for age plots (islet size means cluster at same x)
  if (is_age) {
    p <- p + geom_point(data = donor_means, aes(x = cont_val, y = outcome, color = Region), alpha = 0.5, size = 2)
  }

  if (is_age) {
    p <- p + scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) + labs(x = "Age (years)")
  } else {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks))) + labs(x = "Islet Size (cells)")
  }
  
  p <- p +
    labs(title = paste0(cell_label, ": Region \u00D7 Sex \u00D7 ", ifelse(is_age, "Age", "Islet Size")),
         y = paste0("% ", cell_label, " Cells"), color = "Region", fill = "Region") +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = ann_x, y = y_max * 0.95, label = combined_stars),
      size = 3.5, lineheight = 0.85, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# REGION × SEX × CONTINUOUS
# View B: Facet by Region, Sex overlaid (shows Sex contrast)
# =============================================================================

plot_region_sex_cont_B <- function(data, model, outcome_var, cell_label, continuous_var, y_max = NULL) {
  
  is_age <- continuous_var == "Age_c"
  
  if (is_age) {
    cont_mean <- mean(data$Age)
    cont_range <- range(data$Age_c)
  } else {
    cont_mean <- mean(data$log_Islet.Cells)
    cont_range <- range(data$log_Islet.Cells_c)
  }
  
  model_preds <- get_model_predictions_general(model, continuous_var, cont_range,
                                                by_vars = c("Region", "Sex"))
  if (is_age) {
    model_preds <- model_preds %>% mutate(cont_val = Age_c + cont_mean)
  } else {
    model_preds <- model_preds %>% mutate(cont_val = log_Islet.Cells_c + cont_mean)
  }
  
  slopes_info <- compute_slopes_general(model, continuous_var, c("Region", "Sex"))
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  if (is_age) {
    donor_means <- data %>%
      group_by(Donor, Region, Sex, Age) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
      rename(cont_val = Age)
  } else {
    donor_means <- data %>%
      group_by(Donor, Region, Sex) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE),
                cont_val = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  
  if (is.null(y_max)) y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.2
  
  # Contrast stars: F-M within each Region
  star_df <- tibble(Region = factor(c("Head", "Body", "Tail"), levels = region_levels)) %>%
    mutate(
      star_info = map(Region, function(r) {
        f_draws <- slopes_info$draws[[paste(r, "Female", sep = "_")]]
        m_draws <- slopes_info$draws[[paste(r, "Male", sep = "_")]]
        diff <- f_draws - m_draws
        pd <- pmax(mean(diff > 0), mean(diff < 0))
        tibble(stars = prob_to_stars(pd))
      })
    ) %>%
    unnest(star_info) %>%
    filter(nchar(stars) > 0)
  
  if (is_age) {
    ann_x <- 21
  } else {
    ann_x <- 4.0
  }
  
  p <- ggplot() +
    geom_ribbon(data = model_preds, aes(x = cont_val, ymin = lower, ymax = upper, fill = Sex), alpha = 0.25) +
    geom_line(data = model_preds, aes(x = cont_val, y = estimate, color = Sex), linewidth = 1.3) +
    facet_wrap(~ Region, ncol = 3) +
    scale_color_manual(values = colors_sex) + scale_fill_manual(values = colors_sex) +
    scale_y_continuous(limits = c(0, y_max))
  
  # Only show donor mean points for age plots (islet size means cluster at same x)
  if (is_age) {
    p <- p + geom_point(data = donor_means, aes(x = cont_val, y = outcome, color = Sex), alpha = 0.5, size = 2)
  }

  if (is_age) {
    p <- p + scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) + labs(x = "Age (years)")
  } else {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks))) + labs(x = "Islet Size (cells)")
  }
  
  p <- p +
    labs(title = paste0(cell_label, ": Region \u00D7 Sex \u00D7 ", ifelse(is_age, "Age", "Islet Size")),
         y = paste0("% ", cell_label, " Cells"), color = "Sex", fill = "Sex") +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = ann_x, y = y_max * 0.95, label = stars),
      size = 8, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# 4-WAY: DIAGNOSIS × REGION × SEX × CONTINUOUS
# =============================================================================

plot_4way <- function(data, model, outcome_var, cell_label, continuous_var, y_max = NULL) {
  
  is_age <- continuous_var == "Age_c"
  
  if (is_age) {
    cont_mean <- mean(data$Age)
    cont_range <- range(data$Age_c)
  } else {
    cont_mean <- mean(data$log_Islet.Cells)
    cont_range <- range(data$log_Islet.Cells_c)
  }
  
  model_preds <- get_model_predictions_general(model, continuous_var, cont_range,
                                                by_vars = c("Diagnosis", "Region", "Sex"))
  if (is_age) {
    model_preds <- model_preds %>% mutate(cont_val = Age_c + cont_mean)
  } else {
    model_preds <- model_preds %>% mutate(cont_val = log_Islet.Cells_c + cont_mean)
  }
  
  slopes_info <- compute_slopes_general(model, continuous_var, c("Diagnosis", "Region", "Sex"))
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  
  if (is_age) {
    donor_means <- data %>%
      group_by(Donor, Diagnosis, Region, Sex, Age) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
      rename(cont_val = Age)
  } else {
    donor_means <- data %>%
      group_by(Donor, Diagnosis, Region, Sex) %>%
      summarize(outcome = mean(outcome_pct, na.rm = TRUE),
                cont_val = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  
  if (is.null(y_max)) y_max <- max(c(donor_means$outcome, model_preds$upper), na.rm = TRUE) * 1.2
  
  # Contrast stars: T1D-ND within each Region × Sex
  star_df <- expand.grid(Region = c("Head", "Body", "Tail"), 
                          Sex = c("Female", "Male"), stringsAsFactors = FALSE) %>%
    mutate(Region = factor(Region, levels = region_levels)) %>%
    mutate(
      star_info = map2(Region, Sex, function(r, s) {
        nd_key <- paste("ND", r, s, sep = "_")
        t1d_key <- paste("T1D", r, s, sep = "_")
        nd_draws <- slopes_info$draws[[nd_key]]
        t1d_draws <- slopes_info$draws[[t1d_key]]
        diff <- t1d_draws - nd_draws
        pd <- pmax(mean(diff > 0), mean(diff < 0))
        tibble(stars = prob_to_stars(pd))
      })
    ) %>%
    unnest(star_info) %>%
    filter(nchar(stars) > 0)
  
  p <- ggplot() +
    geom_ribbon(data = model_preds, aes(x = cont_val, ymin = lower, ymax = upper, fill = Diagnosis), alpha = 0.25) +
    geom_line(data = model_preds, aes(x = cont_val, y = estimate, color = Diagnosis), linewidth = 1.2) +
    facet_grid(Region ~ Sex) +
    scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
    scale_y_continuous(limits = c(0, y_max))
  
  # Only show donor mean points for age plots (islet size means cluster at same x)
  if (is_age) {
    p <- p + geom_point(data = donor_means, aes(x = cont_val, y = outcome, color = Diagnosis), alpha = 0.4, size = 1.5)
  }

  if (is_age) {
    p <- p + scale_x_continuous(limits = c(0, 42), breaks = seq(0, 40, 10)) + labs(x = "Age (years)")
    ann_x <- 35
  } else {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks))) + labs(x = "Islet Size (cells)")
    ann_x <- 5.0
  }
  
  p <- p +
    labs(title = paste0(cell_label, ": Diagnosis \u00D7 Region \u00D7 Sex \u00D7 ", ifelse(is_age, "Age", "Islet Size"), " (4-way)"),
         y = paste0("% ", cell_label, " Cells"), color = "Diagnosis", fill = "Diagnosis") +
    theme_pub +
    theme(legend.position = "bottom")
  
  if (nrow(star_df) > 0) {
    p <- p + geom_text(
      data = star_df,
      aes(x = ann_x, y = y_max * 0.90, label = stars),
      size = 6, inherit.aes = FALSE
    )
  }
  
  return(p)
}

# =============================================================================
# MASTER FUNCTION
# =============================================================================

generate_supplementary_figures <- function(data, model_dir = "Proportion/Models",
                                            output_dir = "Proportion/Graphs/Supplementary", save_plots = TRUE) {
  
  cat("\n=============================================================================\n")
  cat("GENERATING SUPPLEMENTARY FIGURES FROM *_mid MODELS\n")
  cat("Two views per 3-way interaction + 4-way figures\n")
  cat("=============================================================================\n")
  
  # Ensure consistent Region ordering across all plots
  data$Region <- factor(data$Region, levels = region_levels)
  
  cell_specs <- list(
    beta  = list(file = "ins_mid_model.rds",  var = "Percent.ins",  label = "Beta Cell",  y_max = 100),
    alpha = list(file = "glu_mid_model.rds",  var = "Percent.glu",  label = "Alpha Cell", y_max = 100),
    delta = list(file = "soma_mid_model.rds", var = "Percent.soma", label = "Delta Cell", y_max = 60),
    pp    = list(file = "pp_mid_model.rds",   var = "Percent.PP",   label = "PP Cell",    y_max = 90)
  )
  
  if (save_plots) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  results <- list()
  
  for (cell_type in names(cell_specs)) {
    spec <- cell_specs[[cell_type]]
    cat("\n--- Processing", spec$label, "---\n")
    model_path <- file.path(model_dir, spec$file)
    cat("  Loading", model_path, "...\n")
    model <- readRDS(model_path)
    
    cat("  Creating Diag \u00D7 Region \u00D7 Sex (2 views)...\n")
    p1a <- plot_diag_region_sex_A(data, model, spec$var, spec$label, spec$y_max)
    p1b <- plot_diag_region_sex_B(data, model, spec$var, spec$label, spec$y_max)
    
    cat("  Creating Diag \u00D7 Sex \u00D7 Age (2 views)...\n")
    p2a <- plot_diag_sex_cont_A(data, model, spec$var, spec$label, "Age_c", spec$y_max)
    p2b <- plot_diag_sex_cont_B(data, model, spec$var, spec$label, "Age_c", spec$y_max)
    
    cat("  Creating Diag \u00D7 Sex \u00D7 Islet (2 views)...\n")
    p3a <- plot_diag_sex_cont_A(data, model, spec$var, spec$label, "log_Islet.Cells_c", spec$y_max)
    p3b <- plot_diag_sex_cont_B(data, model, spec$var, spec$label, "log_Islet.Cells_c", spec$y_max)
    
    cat("  Creating Region \u00D7 Sex \u00D7 Age (2 views)...\n")
    p4a <- plot_region_sex_cont_A(data, model, spec$var, spec$label, "Age_c", spec$y_max)
    p4b <- plot_region_sex_cont_B(data, model, spec$var, spec$label, "Age_c", spec$y_max)
    
    cat("  Creating Region \u00D7 Sex \u00D7 Islet (2 views)...\n")
    p5a <- plot_region_sex_cont_A(data, model, spec$var, spec$label, "log_Islet.Cells_c", spec$y_max)
    p5b <- plot_region_sex_cont_B(data, model, spec$var, spec$label, "log_Islet.Cells_c", spec$y_max)
    
    cat("  Creating 4-way Age...\n")
    p6 <- plot_4way(data, model, spec$var, spec$label, "Age_c", spec$y_max)
    
    cat("  Creating 4-way Islet...\n")
    p7 <- plot_4way(data, model, spec$var, spec$label, "log_Islet.Cells_c", spec$y_max)
    
    rm(model); gc()
    
    results[[cell_type]] <- list(
      diag_region_sex_A = p1a, diag_region_sex_B = p1b,
      diag_sex_age_A = p2a, diag_sex_age_B = p2b,
      diag_sex_islet_A = p3a, diag_sex_islet_B = p3b,
      region_sex_age_A = p4a, region_sex_age_B = p4b,
      region_sex_islet_A = p5a, region_sex_islet_B = p5b,
      fourway_age = p6, fourway_islet = p7
    )
    
    if (save_plots) {
      ggsave(file.path(output_dir, paste0(cell_type, "_diag_region_sex_byRegion.png")), p1a, width = 14, height = 6, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_diag_region_sex_byDiag.png")), p1b, width = 12, height = 6, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_diag_sex_age_bySex.png")), p2a, width = 12, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_diag_sex_age_byDiag.png")), p2b, width = 12, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_diag_sex_islet_bySex.png")), p3a, width = 12, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_diag_sex_islet_byDiag.png")), p3b, width = 12, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_region_sex_age_bySex.png")), p4a, width = 12, height = 6, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_region_sex_age_byRegion.png")), p4b, width = 14, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_region_sex_islet_bySex.png")), p5a, width = 12, height = 6, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_region_sex_islet_byRegion.png")), p5b, width = 14, height = 5, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_4way_age.png")), p6, width = 10, height = 10, dpi = 300, bg = "white")
      ggsave(file.path(output_dir, paste0(cell_type, "_4way_islet.png")), p7, width = 10, height = 10, dpi = 300, bg = "white")
      cat("  Saved 12 figures to", output_dir, "\n")
    }
  }
  
  cat("\n=============================================================================\n")
  cat("COMPLETE! Generated", length(cell_specs) * 12, "supplementary figures.\n")
  cat("=============================================================================\n")
  
  return(results)
}

if (interactive()) {
  cat("\nTo generate all supplementary figures:\n")
  cat("  figures <- generate_supplementary_figures(quad_data)\n")
  cat("  # Models loaded individually from Proportion/Models/\n")
}
