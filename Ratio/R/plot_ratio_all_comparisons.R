# =============================================================================
# plot_ratio_all_comparisons.R
# All possible post-hoc comparison plots for Ratio ordbetareg models
# =============================================================================
#
# Outcome: beta_alpha_ratio = beta / (alpha + beta)
# Data: ratio_data (cutoff_data, Stain != Single, Donor != 6473, ins+glu > 0)
#
# DIAGNOSIS MODELS (powered/mid/full):
#   Same formula structure as proportions — reuses proportion plot functions.
#   Powered: (Diag + Region + Sex + Age_c + log_IC_c)^2 + Diag:Region:{Age,Islet}
#   Mid:     + DRS, DS{Age,Islet}, RS{Age,Islet}, 4-way
#   Full:    (5 vars)^5
#   → 35 plots
#
# T1D MODEL (AO parameterization only):
#   DD + Age and DD + AO are algebraically equivalent (Age = DD + AO).
#   Only the AO model is fitted; DD slopes extracted from it directly.
#   (DD_c + Region + Sex + AO_c + log_IC_c)^2 +
#   DD:Region:log_IC_c + DD:Region:AO_c
#   → 19 plots (DD, AO, and Islet slopes all from one model)
#
# Total: ~54 plots
#
# =============================================================================

library(tidyverse)
library(brms)
library(coda)
library(patchwork)
library(ggtext)

source("colors_master.R")

# Source proportion helpers (posterior functions, diagnosis-tier plot functions)
source("Proportion/R/plot_proportion_all_comparisons.R")

# Also source supplementary functions for mid model
source("Proportion/R/plot_celltype_supplementary_figures_no_stars.R")

# =============================================================================
# SHARED SETTINGS
# =============================================================================

ratio_label <- "\u03B2/(\u03B1+\u03B2) Ratio"
ratio_var   <- "beta_alpha_ratio"
ratio_y_max <- 100

# Bold ratio label for titles (rendered via ggtext::element_markdown)
ratio_label_bold <- "**\u03B2/(\u03B1+\u03B2) Ratio**"

# Theme addition to render markdown in plot titles
theme_bold_title <- theme(plot.title = element_markdown())


# #############################################################################
#
#  T1D MODEL HELPERS (DD and AO parameterizations)
#
# #############################################################################

# Generic posterior means for T1D categorical grid (no Diagnosis factor)
get_t1d_categorical_means <- function(model, data, by_vars, cont_vars) {
  # cont_vars: named list of continuous variable names and their zero values
  grid <- expand.grid(
    Region = factor(c("Head", "Body", "Tail"), levels = region_levels),
    Sex = c("Female", "Male"),
    stringsAsFactors = FALSE
  )
  for (cv in names(cont_vars)) grid[[cv]] <- cont_vars[[cv]]
  
  post <- get_posterior_fitted(model, grid)
  combos <- distinct(grid[, by_vars, drop = FALSE])
  results <- list()
  for (i in seq_len(nrow(combos))) {
    mask <- rep(TRUE, nrow(grid))
    for (v in by_vars) mask <- mask & grid[[v]] == combos[[v]][i]
    draws <- rowMeans(post[, which(mask), drop = FALSE]) * 100
    summ <- summarize_posterior(draws)
    for (v in by_vars) summ[[v]] <- combos[[v]][i]
    results[[i]] <- summ
  }
  bind_rows(results) %>%
    mutate(estimate = estimate, lower = lower.HPD, upper = upper.HPD) %>%
    select(all_of(by_vars), estimate, lower, upper)
}

# Posterior pairwise contrast pd for T1D binary factor
get_t1d_binary_pd <- function(model, data, focal_var, cond_var, cont_vars) {
  grid <- expand.grid(
    Region = factor(c("Head", "Body", "Tail"), levels = region_levels),
    Sex = c("Female", "Male"), stringsAsFactors = FALSE
  )
  for (cv in names(cont_vars)) grid[[cv]] <- cont_vars[[cv]]
  post <- get_posterior_fitted(model, grid)
  focal_levels <- unique(grid[[focal_var]])
  
  if (is.null(cond_var)) {
    idx_a <- which(grid[[focal_var]] == focal_levels[1])
    idx_b <- which(grid[[focal_var]] == focal_levels[2])
    diff <- rowMeans(post[, idx_a, drop = FALSE]) - rowMeans(post[, idx_b, drop = FALSE])
    pd <- pmax(mean(diff > 0), mean(diff < 0))
    return(tibble(contrast = paste(focal_levels[1], "-", focal_levels[2]),
                  pd = pd, stars = prob_to_stars(pd)))
  }
  
  cond_levels <- unique(grid[[cond_var]])
  results <- list()
  for (cl in cond_levels) {
    mask <- grid[[cond_var]] == cl
    idx_a <- which(mask & grid[[focal_var]] == focal_levels[1])
    idx_b <- which(mask & grid[[focal_var]] == focal_levels[2])
    diff <- rowMeans(post[, idx_a, drop = FALSE]) - rowMeans(post[, idx_b, drop = FALSE])
    pd <- pmax(mean(diff > 0), mean(diff < 0))
    results[[length(results) + 1]] <- tibble(
      !!cond_var := cl, contrast = paste(focal_levels[1], "-", focal_levels[2]),
      pd = pd, stars = prob_to_stars(pd))
  }
  bind_rows(results)
}

# Region pairwise pd for T1D
get_t1d_region_pd <- function(model, data, cond_var, cont_vars) {
  grid <- expand.grid(
    Region = factor(c("Head", "Body", "Tail"), levels = region_levels),
    Sex = c("Female", "Male"), stringsAsFactors = FALSE
  )
  for (cv in names(cont_vars)) grid[[cv]] <- cont_vars[[cv]]
  post <- get_posterior_fitted(model, grid)
  pairs <- list(c("Head", "Body"), c("Head", "Tail"), c("Body", "Tail"))
  
  if (is.null(cond_var)) {
    results <- list()
    for (pr in pairs) {
      diff <- rowMeans(post[, which(grid$Region == pr[1]), drop = FALSE]) -
        rowMeans(post[, which(grid$Region == pr[2]), drop = FALSE])
      pd <- pmax(mean(diff > 0), mean(diff < 0))
      results[[length(results) + 1]] <- tibble(
        contrast = paste(pr[1], "-", pr[2]), pd = pd, stars = prob_to_stars(pd))
    }
    return(bind_rows(results))
  }
  
  cond_levels <- unique(grid[[cond_var]])
  results <- list()
  for (cl in cond_levels) {
    mask <- grid[[cond_var]] == cl
    for (pr in pairs) {
      diff <- rowMeans(post[, which(mask & grid$Region == pr[1]), drop = FALSE]) -
        rowMeans(post[, which(mask & grid$Region == pr[2]), drop = FALSE])
      pd <- pmax(mean(diff > 0), mean(diff < 0))
      results[[length(results) + 1]] <- tibble(
        !!cond_var := cl, contrast = paste(pr[1], "-", pr[2]),
        pd = pd, stars = prob_to_stars(pd))
    }
  }
  bind_rows(results)
}

# T1D continuous predictions — marginalized over non-focal categoricals
get_t1d_cont_preds <- function(model, data, focal_var, by_vars,
                               all_cont_vars, n_points = 50) {
  focal_range <- range(data[[focal_var]], na.rm = TRUE)
  focal_seq <- seq(focal_range[1], focal_range[2], length.out = n_points)
  
  grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), stringsAsFactors = FALSE
  )
  for (cv in all_cont_vars) grid[[cv]] <- 0
  # Overwrite focal with sequence
  grid <- grid[rep(seq_len(nrow(grid)), each = n_points), ]
  grid[[focal_var]] <- rep(focal_seq, times = nrow(grid) / n_points)
  
  post <- get_posterior_fitted(model, grid)
  
  grid %>%
    mutate(cell_id = row_number()) %>%
    group_by(across(all_of(c(by_vars, focal_var)))) %>%
    summarize(cell_ids = list(cell_id), .groups = "drop") %>%
    mutate(post_summary = map(cell_ids, function(ids) {
      draws <- rowMeans(post[, ids, drop = FALSE])
      hpd <- coda::HPDinterval(coda::as.mcmc(draws), prob = 0.95)
      tibble(estimate = median(draws) * 100, lower = hpd[1] * 100, upper = hpd[2] * 100)
    })) %>% unnest(post_summary) %>% select(-cell_ids)
}

# T1D slope computation via finite differences
compute_t1d_slopes <- function(model, data, focal_var, by_vars, all_cont_vars) {
  grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), stringsAsFactors = FALSE
  )
  for (cv in all_cont_vars) grid[[cv]] <- 0
  grid_lo <- grid; grid_hi <- grid
  grid_lo[[focal_var]] <- -1; grid_hi[[focal_var]] <- 1
  grid_both <- bind_rows(grid_lo %>% mutate(.level = "lo"),
                         grid_hi %>% mutate(.level = "hi"))
  post <- get_posterior_fitted(model, grid_both %>% select(-.level))
  
  n_half <- nrow(grid_lo)
  combos <- if (length(by_vars) > 0) distinct(grid[, by_vars, drop = FALSE]) else tibble(.dummy = 1)
  
  slopes <- list(); slope_draws <- list()
  for (i in seq_len(nrow(combos))) {
    mask <- rep(TRUE, n_half)
    if (length(by_vars) > 0)
      for (v in by_vars) mask <- mask & grid[[v]] == combos[[v]][i]
    lo_idx <- which(mask)
    hi_idx <- lo_idx + n_half
    draws <- (rowMeans(post[, hi_idx, drop = FALSE]) -
                rowMeans(post[, lo_idx, drop = FALSE])) / 2
    summ <- summarize_posterior(draws)
    if (length(by_vars) > 0) for (v in by_vars) summ[[v]] <- combos[[v]][i]
    summ$pd <- pmax(mean(draws > 0), mean(draws < 0))
    summ$stars <- prob_to_stars(summ$pd)
    key <- if (length(by_vars) > 0) paste(sapply(by_vars, function(v) combos[[v]][i]), collapse = "_") else "overall"
    slopes[[key]] <- summ; slope_draws[[key]] <- draws
  }
  list(slopes = bind_rows(slopes), draws = slope_draws)
}


# #############################################################################
#
#  T1D PLOT FUNCTIONS
#
# #############################################################################

# --- Categorical: Region ---
plot_t1d_region <- function(data, model, cont_vars, y_lab) {
  dm <- data %>% group_by(Donor, Region) %>%
    summarise(md = mean(beta_alpha_ratio * 100, na.rm = TRUE), .groups = "drop")
  emm <- get_t1d_categorical_means(model, data, "Region", cont_vars)
  ctr <- get_t1d_region_pd(model, data, NULL, cont_vars)
  sig <- ctr %>% filter(nchar(stars) > 0) %>% mutate(label = paste(contrast, stars))
  sig_text <- paste(sig$label, collapse = "\n")
  
  ggplot() +
    geom_violin(data = dm, aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm, aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    { if (nchar(sig_text) > 0) annotate("text", x = 2, y = max(dm$md) * 0.95,
                                        label = sig_text, size = 4, lineheight = 0.9) } +
    labs(x = "Pancreatic region", y = y_lab, title = "Region (T1D)") +
    theme_pub
}

# --- Categorical: Sex ---
plot_t1d_sex <- function(data, model, cont_vars, y_lab) {
  dm <- data %>% group_by(Donor, Sex) %>%
    summarise(md = mean(beta_alpha_ratio * 100, na.rm = TRUE), .groups = "drop")
  emm <- get_t1d_categorical_means(model, data, "Sex", cont_vars)
  ctr <- get_t1d_binary_pd(model, data, "Sex", NULL, cont_vars)
  
  ggplot() +
    geom_violin(data = dm, aes(x = Sex, y = md, fill = Sex),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm, aes(x = Sex, y = md, color = Sex),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_sex, guide = "none") +
    scale_color_manual(values = colors_sex, guide = "none") +
    { if (nchar(ctr$stars) > 0) annotate("text", x = 1.5, y = max(dm$md) * 0.95,
                                         label = ctr$stars, size = 8) } +
    labs(x = NULL, y = y_lab, title = "Sex (T1D)") +
    theme_pub
}

# --- Region | Sex and Sex | Region ---
plot_t1d_region_by_sex <- function(data, model, cont_vars, y_lab) {
  dm <- data %>% group_by(Donor, Region, Sex) %>%
    summarise(md = mean(beta_alpha_ratio * 100, na.rm = TRUE), .groups = "drop")
  emm <- get_t1d_categorical_means(model, data, c("Region", "Sex"), cont_vars)
  ctr <- get_t1d_region_pd(model, data, "Sex", cont_vars)
  sig_labels <- ctr %>% filter(nchar(stars) > 0) %>%
    mutate(label = paste(contrast, stars)) %>%
    group_by(Sex) %>% summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_sex <- tibble(Sex = c("Female", "Male"))
  sig_labels <- left_join(all_sex, sig_labels, by = "Sex") %>% mutate(label = replace_na(label, ""))
  y_max <- max(dm$md, na.rm = TRUE) * 1.05
  
  ggplot() +
    geom_violin(data = dm, aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm, aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ Sex) +
    geom_text(data = sig_labels, aes(x = 2, y = y_max * 0.95, label = label),
              size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = y_lab,
         title = "Region \u00D7 Sex (T1D)", subtitle = "Region contrasts within each Sex") +
    theme_pub
}

plot_t1d_sex_by_region <- function(data, model, cont_vars, y_lab) {
  dm <- data %>% group_by(Donor, Region, Sex) %>%
    summarise(md = mean(beta_alpha_ratio * 100, na.rm = TRUE), .groups = "drop")
  emm <- get_t1d_categorical_means(model, data, c("Region", "Sex"), cont_vars)
  ctr <- get_t1d_binary_pd(model, data, "Sex", "Region", cont_vars)
  y_max <- max(dm$md, na.rm = TRUE) * 1.05
  
  ggplot() +
    geom_violin(data = dm, aes(x = Region, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm, aes(x = Region, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) + scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:3, y = y_max * c(0.98, 0.93, 0.88), label = ctr$stars, size = 8) +
    labs(x = "Pancreatic region", y = y_lab,
         title = "Sex \u00D7 Region (T1D)", subtitle = "Sex contrast within each Region",
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# --- Generic T1D continuous slope plot ---
plot_t1d_slope <- function(data, model, focal_var, raw_var, by_var, colors,
                           all_cont_vars, y_lab, x_lab, title_text) {
  preds <- get_t1d_cont_preds(model, data, focal_var,
                              if (is.null(by_var)) character(0) else by_var,
                              all_cont_vars)
  mean_offset <- mean(data[[raw_var]], na.rm = TRUE)
  preds$x <- preds[[focal_var]] + mean_offset
  
  dm <- data %>% group_by(Donor, !!!syms(if (is.null(by_var)) character(0) else by_var)) %>%
    summarise(md = mean(beta_alpha_ratio * 100, na.rm = TRUE),
              x = mean(.data[[raw_var]], na.rm = TRUE), .groups = "drop")
  
  if (is.null(by_var)) {
    slope_info <- compute_t1d_slopes(model, data, focal_var, character(0), all_cont_vars)
    stars <- slope_info$slopes$stars[1]
    
    plt <- ggplot() +
      geom_ribbon(data = preds, aes(x = x, ymin = lower, ymax = upper),
                  fill = "grey60", alpha = 0.3) +
      geom_line(data = preds, aes(x = x, y = estimate), linewidth = 0.9) +
      geom_point(data = dm, aes(x = x, y = md), color = col_t1d, alpha = 0.6, size = 2) +
      labs(x = x_lab, y = y_lab, title = title_text) + theme_pub
    if (nchar(stars) > 0)
      plt <- plt + annotate("text", x = min(preds$x) + 1, y = max(dm$md) * 0.95,
                            label = stars, size = 8, hjust = 0)
    return(plt)
  }
  
  slope_info <- compute_t1d_slopes(model, data, focal_var, by_var, all_cont_vars)
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = x, ymin = lower, ymax = upper,
                                  fill = .data[[by_var]]), alpha = 0.2) +
    geom_line(data = preds, aes(x = x, y = estimate, color = .data[[by_var]]),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = x, y = md, color = .data[[by_var]]),
               alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
    labs(x = x_lab, y = y_lab, title = title_text) +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  
  # Annotate significant pairwise contrasts
  if (by_var == "Region") {
    pairs_df <- slope_info$slopes
    # Done via annotation of overall text
    sig_pairs <- list()
    for (pr in list(c("Head","Body"), c("Head","Tail"), c("Body","Tail"))) {
      d1 <- slope_info$draws[[pr[1]]]; d2 <- slope_info$draws[[pr[2]]]
      if (!is.null(d1) && !is.null(d2)) {
        diff <- d1 - d2; pd <- pmax(mean(diff > 0), mean(diff < 0))
        s <- prob_to_stars(pd)
        if (nchar(s) > 0) sig_pairs[[length(sig_pairs)+1]] <- paste(pr[1], "-", pr[2], s)
      }
    }
    if (length(sig_pairs) > 0)
      plt <- plt + annotate("text", x = min(preds$x) + 1, y = max(dm$md) * 0.95,
                            label = paste(sig_pairs, collapse = "\n"),
                            size = 4, hjust = 0, lineheight = 0.9)
  } else {
    s <- slope_info$slopes
    if (nrow(s) == 2) {
      d1 <- slope_info$draws[[s[[by_var]][1]]]; d2 <- slope_info$draws[[s[[by_var]][2]]]
      if (!is.null(d1) && !is.null(d2)) {
        diff <- d1 - d2; pd <- pmax(mean(diff > 0), mean(diff < 0))
        star <- prob_to_stars(pd)
        if (nchar(star) > 0)
          plt <- plt + annotate("text", x = min(preds$x) + 1, y = max(dm$md) * 0.95,
                                label = star, size = 8, hjust = 0)
      }
    }
  }
  plt
}

# --- T1D 3-way: continuous × Region, faceted at ±1 SD of another continuous ---
plot_t1d_3way <- function(data, model, focal_var, raw_focal, cond_var,
                          all_cont_vars, y_lab, x_lab, title_text) {
  focal_range <- range(data[[focal_var]], na.rm = TRUE)
  focal_seq <- seq(focal_range[1], focal_range[2], length.out = 50)
  mean_offset <- mean(data[[raw_focal]], na.rm = TRUE)
  cond_sd <- sd(data[[cond_var]], na.rm = TRUE)
  
  preds_list <- list()
  for (lev in c(-1, 1)) {
    grid <- expand.grid(
      Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"), stringsAsFactors = FALSE
    )
    for (cv in all_cont_vars) grid[[cv]] <- 0
    grid <- grid[rep(seq_len(nrow(grid)), each = length(focal_seq)), ]
    grid[[focal_var]] <- rep(focal_seq, times = nrow(grid) / length(focal_seq))
    grid[[cond_var]] <- lev * cond_sd
    
    post <- get_posterior_fitted(model, grid)
    result <- grid %>% mutate(cell_id = row_number()) %>%
      group_by(Region, !!sym(focal_var)) %>%
      summarize(cell_ids = list(cell_id), .groups = "drop") %>%
      mutate(post_summary = map(cell_ids, function(ids) {
        draws <- rowMeans(post[, ids, drop = FALSE])
        hpd <- coda::HPDinterval(coda::as.mcmc(draws), prob = 0.95)
        tibble(estimate = median(draws) * 100, lower = hpd[1] * 100, upper = hpd[2] * 100)
      })) %>% unnest(post_summary) %>% select(-cell_ids)
    
    cond_label_short <- sub("_c$", "", cond_var)
    result$cond_level <- factor(
      ifelse(lev == -1, paste0(cond_label_short, " = mean \u2212 1 SD"),
             paste0(cond_label_short, " = mean + 1 SD")))
    preds_list[[length(preds_list) + 1]] <- result
  }
  preds <- bind_rows(preds_list)
  preds$x <- preds[[focal_var]] + mean_offset
  preds$Region <- factor(preds$Region, levels = region_levels)
  y_max <- max(preds$upper, na.rm = TRUE) * 1.1
  
  ggplot(preds) +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
    geom_line(aes(x = x, y = estimate, color = Region), linewidth = 0.9) +
    facet_wrap(~ cond_level) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(x = x_lab, y = y_lab, title = title_text) +
    theme_pub + theme(legend.position = "bottom")
}


# #############################################################################
#
#  MASTER EXECUTION
#
# #############################################################################

generate_all_ratio_plots <- function(model_dir = "Ratio/Models",
                                     output_dir = "Ratio/Graphs/All_Comparisons") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING ALL RATIO POST-HOC COMPARISON PLOTS\n")
  cat("=============================================================================\n\n")
  
  # --- Load data ---
  cat("Loading data...\n")
  ratio_data <- readRDS("Data/ratio_data.rds") %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
           Region = factor(Region, levels = region_levels),
           Sex = factor(Sex, levels = c("Female", "Male")),
           Donor = factor(Donor), ImageID = factor(ImageID))
  
  ratio_data_t1d <- readRDS("Data/ratio_data_t1d.rds") %>%
    mutate(Region = factor(Region, levels = region_levels),
           Sex = factor(Sex, levels = c("Female", "Male")),
           Donor = factor(Donor), ImageID = factor(ImageID)) %>%
    droplevels()
  
  cat("  ratio_data:", nrow(ratio_data), "rows,", n_distinct(ratio_data$Donor), "donors\n")
  cat("  ratio_t1d: ", nrow(ratio_data_t1d), "rows,", n_distinct(ratio_data_t1d$Donor), "donors\n\n")
  
  y_lab <- paste0(ratio_label, " (%)")
  
  # =========================================================================
  # DIAGNOSIS TIER (powered / mid / full) — reuse proportion functions
  # =========================================================================
  
  dir_pow  <- file.path(output_dir, "Diagnosis", "Powered")
  dir_mid  <- file.path(output_dir, "Diagnosis", "Mid")
  dir_full <- file.path(output_dir, "Diagnosis", "Full")
  
  # --- Powered ---
  cat("--- Diagnosis Powered (21 plots) ---\n")
  m <- readRDS(file.path(model_dir, "ratio_powered_model.rds"))
  
  save_plot(plot_main_diagnosis(ratio_data, m, ratio_var, ratio_label), "01_main_diagnosis", dir_pow)
  save_plot(plot_main_region(ratio_data, m, ratio_var, ratio_label), "02_main_region", dir_pow)
  save_plot(plot_main_sex(ratio_data, m, ratio_var, ratio_label), "03_main_sex", dir_pow)
  save_plot(plot_age_overall(ratio_data, m, ratio_var, ratio_label), "04_age_overall", dir_pow)
  save_plot(plot_islet_overall(ratio_data, m, ratio_var, ratio_label), "05_islet_overall", dir_pow)
  save_plot(plot_diagnosis_by_region(ratio_data, m, ratio_var, ratio_label), "06_diagnosis_by_region", dir_pow, 8)
  save_plot(plot_region_by_diagnosis(ratio_data, m, ratio_var, ratio_label), "07_region_by_diagnosis", dir_pow, 9, 5)
  save_plot(plot_diagnosis_by_sex(ratio_data, m, ratio_var, ratio_label), "08_diagnosis_by_sex", dir_pow)
  save_plot(plot_sex_by_diagnosis(ratio_data, m, ratio_var, ratio_label), "09_sex_by_diagnosis", dir_pow)
  save_plot(plot_region_by_sex(ratio_data, m, ratio_var, ratio_label), "10_region_by_sex", dir_pow, 9, 5)
  save_plot(plot_sex_by_region(ratio_data, m, ratio_var, ratio_label), "11_sex_by_region", dir_pow, 8)
  save_plot(plot_cont_by_cat(ratio_data, m, ratio_var, ratio_label, "Age_c", "Diagnosis", colors_diag, "Diagnosis"),
            "12_age_by_diagnosis", dir_pow)
  save_plot(plot_cont_by_cat(ratio_data, m, ratio_var, ratio_label, "Age_c", "Region", colors_region, "Region"),
            "13_age_by_region", dir_pow)
  save_plot(plot_cont_by_cat(ratio_data, m, ratio_var, ratio_label, "Age_c", "Sex", colors_sex, "Sex"),
            "14_age_by_sex", dir_pow)
  save_plot(plot_cont_by_cat(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", "Diagnosis", colors_diag, "Diagnosis"),
            "15_islet_by_diagnosis", dir_pow)
  save_plot(plot_cont_by_cat(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", "Region", colors_region, "Region"),
            "16_islet_by_region", dir_pow)
  save_plot(plot_cont_by_cat(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", "Sex", colors_sex, "Sex"),
            "17_islet_by_sex", dir_pow)
  save_plot(plot_3way_facet_region(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", ratio_y_max),
            "18_3way_islet_facet_region", dir_pow, 14, 5)
  save_plot(plot_3way_facet_diagnosis(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", ratio_y_max),
            "19_3way_islet_facet_diagnosis", dir_pow, 12, 6)
  save_plot(plot_3way_facet_region(ratio_data, m, ratio_var, ratio_label, "Age_c", ratio_y_max),
            "20_3way_age_facet_region", dir_pow, 14, 5)
  save_plot(plot_3way_facet_diagnosis(ratio_data, m, ratio_var, ratio_label, "Age_c", ratio_y_max),
            "21_3way_age_facet_diagnosis", dir_pow, 12, 6)
  rm(m); gc()
  
  # --- Mid ---
  cat("\n--- Diagnosis Mid (12 supplementary plots) ---\n")
  m <- readRDS(file.path(model_dir, "ratio_mid_model.rds"))
  save_plot(plot_diag_region_sex_A(ratio_data, m, ratio_var, ratio_label, ratio_y_max),
            "22_diag_region_sex_byRegion", dir_mid, 14, 6)
  save_plot(plot_diag_region_sex_B(ratio_data, m, ratio_var, ratio_label, ratio_y_max),
            "23_diag_region_sex_byDiag", dir_mid, 12, 6)
  save_plot(plot_diag_sex_cont_A(ratio_data, m, ratio_var, ratio_label, "Age_c", ratio_y_max),
            "24_diag_sex_age_bySex", dir_mid, 12, 5)
  save_plot(plot_diag_sex_cont_B(ratio_data, m, ratio_var, ratio_label, "Age_c", ratio_y_max),
            "25_diag_sex_age_byDiag", dir_mid, 12, 5)
  save_plot(plot_diag_sex_cont_A(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", ratio_y_max),
            "26_diag_sex_islet_bySex", dir_mid, 12, 5)
  save_plot(plot_diag_sex_cont_B(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", ratio_y_max),
            "27_diag_sex_islet_byDiag", dir_mid, 12, 5)
  save_plot(plot_region_sex_cont_A(ratio_data, m, ratio_var, ratio_label, "Age_c", ratio_y_max),
            "28_region_sex_age_bySex", dir_mid, 12, 6)
  save_plot(plot_region_sex_cont_B(ratio_data, m, ratio_var, ratio_label, "Age_c", ratio_y_max),
            "29_region_sex_age_byRegion", dir_mid, 14, 5)
  save_plot(plot_region_sex_cont_A(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", ratio_y_max),
            "30_region_sex_islet_bySex", dir_mid, 12, 6)
  save_plot(plot_region_sex_cont_B(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", ratio_y_max),
            "31_region_sex_islet_byRegion", dir_mid, 14, 5)
  save_plot(plot_4way(ratio_data, m, ratio_var, ratio_label, "Age_c", ratio_y_max),
            "32_4way_age", dir_mid, 10, 10)
  save_plot(plot_4way(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", ratio_y_max),
            "33_4way_islet", dir_mid, 10, 10)
  rm(m); gc()
  
  # --- Full ---
  cat("\n--- Diagnosis Full (2 proof-of-concept plots) ---\n")
  m <- readRDS(file.path(model_dir, "ratio_full_model.rds"))
  save_plot(plot_5way(ratio_data, m, ratio_var, ratio_label, "Age_c", "log_Islet.Cells_c", ratio_y_max),
            "34_5way_age_x_islet_cond", dir_full, 16, 10)
  save_plot(plot_5way(ratio_data, m, ratio_var, ratio_label, "log_Islet.Cells_c", "Age_c", ratio_y_max),
            "35_5way_islet_x_age_cond", dir_full, 16, 10)
  rm(m); gc()
  
  
  # =========================================================================
  # T1D MODEL (AO parameterization — DD slopes extracted from same model)
  # Since Age = DD + AO, the DD and AO models are algebraically equivalent.
  # Only the AO model is fitted; DD slopes come out of it directly.
  # =========================================================================
  
  dir_t1d <- file.path(output_dir, "T1D")
  cat("\n--- T1D Model (21 plots from AO parameterization) ---\n")
  m <- readRDS(file.path(model_dir, "ratio_ao_model.rds"))
  t1d_cont <- c("Disease.Duration_c" = 0, "age_at_onset_c" = 0, "log_Islet.Cells_c" = 0)
  t1d_all <- c("Disease.Duration_c", "age_at_onset_c", "log_Islet.Cells_c")
  
  # Categorical
  save_plot(plot_t1d_region(ratio_data_t1d, m, t1d_cont, y_lab), "01_region", dir_t1d)
  save_plot(plot_t1d_sex(ratio_data_t1d, m, t1d_cont, y_lab), "02_sex", dir_t1d)
  save_plot(plot_t1d_region_by_sex(ratio_data_t1d, m, t1d_cont, y_lab), "03_region_by_sex", dir_t1d, 9, 5)
  save_plot(plot_t1d_sex_by_region(ratio_data_t1d, m, t1d_cont, y_lab), "04_sex_by_region", dir_t1d, 8)
  
  # DD slopes (extracted from AO model)
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "Disease.Duration_c", "Disease.Duration", NULL, NULL,
                           t1d_all, y_lab, "Disease duration (years)", "DD overall (T1D)"),
            "05_dd_overall", dir_t1d)
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "Disease.Duration_c", "Disease.Duration", "Region", colors_region,
                           t1d_all, y_lab, "Disease duration (years)", "DD \u00D7 Region (T1D)"),
            "06_dd_by_region", dir_t1d)
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "Disease.Duration_c", "Disease.Duration", "Sex", colors_sex,
                           t1d_all, y_lab, "Disease duration (years)", "DD \u00D7 Sex (T1D)"),
            "07_dd_by_sex", dir_t1d)
  
  # AO slopes
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "age_at_onset_c", "age_at_onset", NULL, NULL,
                           t1d_all, y_lab, "Age at onset (years)", "AO overall (T1D)"),
            "08_ao_overall", dir_t1d)
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "age_at_onset_c", "age_at_onset", "Region", colors_region,
                           t1d_all, y_lab, "Age at onset (years)", "AO \u00D7 Region (T1D)"),
            "09_ao_by_region", dir_t1d)
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "age_at_onset_c", "age_at_onset", "Sex", colors_sex,
                           t1d_all, y_lab, "Age at onset (years)", "AO \u00D7 Sex (T1D)"),
            "10_ao_by_sex", dir_t1d)
  
  # Islet slopes
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "log_Islet.Cells_c", "log_Islet.Cells", NULL, NULL,
                           t1d_all, y_lab, "Islet size (cells)", "Islet overall (T1D)"),
            "11_islet_overall", dir_t1d)
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "log_Islet.Cells_c", "log_Islet.Cells", "Region", colors_region,
                           t1d_all, y_lab, "Islet size (cells)", "Islet \u00D7 Region (T1D)"),
            "12_islet_by_region", dir_t1d)
  save_plot(plot_t1d_slope(ratio_data_t1d, m, "log_Islet.Cells_c", "log_Islet.Cells", "Sex", colors_sex,
                           t1d_all, y_lab, "Islet size (cells)", "Islet \u00D7 Sex (T1D)"),
            "13_islet_by_sex", dir_t1d)
  
  # 3-way: DD × Region × AO
  save_plot(plot_t1d_3way(ratio_data_t1d, m, "Disease.Duration_c", "Disease.Duration", "age_at_onset_c",
                          t1d_all, y_lab, "Disease duration (years)",
                          "DD \u00D7 Region \u00D7 AO (T1D)"),
            "14_3way_dd_region_ao", dir_t1d, 10, 5)
  save_plot(plot_t1d_3way(ratio_data_t1d, m, "age_at_onset_c", "age_at_onset", "Disease.Duration_c",
                          t1d_all, y_lab, "Age at onset (years)",
                          "AO \u00D7 Region \u00D7 DD (T1D)"),
            "15_3way_ao_region_dd", dir_t1d, 10, 5)
  
  # 3-way: DD × Region × Islet
  save_plot(plot_t1d_3way(ratio_data_t1d, m, "Disease.Duration_c", "Disease.Duration", "log_Islet.Cells_c",
                          t1d_all, y_lab, "Disease duration (years)",
                          "DD \u00D7 Region \u00D7 Islet (T1D)"),
            "16_3way_dd_region_islet", dir_t1d, 10, 5)
  save_plot(plot_t1d_3way(ratio_data_t1d, m, "log_Islet.Cells_c", "log_Islet.Cells", "Disease.Duration_c",
                          t1d_all, y_lab, "Islet size (cells)",
                          "Islet \u00D7 Region \u00D7 DD (T1D)"),
            "17_3way_islet_region_dd", dir_t1d, 10, 5)
  
  # 3-way: AO × Region × Islet
  save_plot(plot_t1d_3way(ratio_data_t1d, m, "age_at_onset_c", "age_at_onset", "log_Islet.Cells_c",
                          t1d_all, y_lab, "Age at onset (years)",
                          "AO \u00D7 Region \u00D7 Islet (T1D)"),
            "18_3way_ao_region_islet", dir_t1d, 10, 5)
  save_plot(plot_t1d_3way(ratio_data_t1d, m, "log_Islet.Cells_c", "log_Islet.Cells", "age_at_onset_c",
                          t1d_all, y_lab, "Islet size (cells)",
                          "Islet \u00D7 Region \u00D7 AO (T1D)"),
            "19_3way_islet_region_ao", dir_t1d, 10, 5)
  rm(m); gc()
  
  
  cat("\n=============================================================================\n")
  cat("COMPLETE! ~54 plots saved to", output_dir, "\n")
  cat("  Diagnosis/Powered (21), Mid (12), Full (2)\n")
  cat("  T1D (19 — DD & AO slopes from single AO model)\n")
  cat("=============================================================================\n")
}


# #############################################################################
#
#  SUPPLEMENTARY FIGURE S10 — β/(α+β) Ratio overview
#
#  Changes from earlier version:
#   - Titles all use Greek letters and are bolded (via ggtext::element_markdown)
#   - Subtitles removed from panels A and B
#   - Legend repositioned on panels C and D
#   - Points removed from panel C (Diagnosis × Islet Size)
#
# #############################################################################

generate_supplementary_figure_S10 <- function(
    model_dir  = "Ratio/Models",
    output_file = "Ratio/Graphs/Supplementary_Figure_S10_-_Beta-Alpha-cell_ratios.pdf") {
  
  cat("\n--- Generating Supplementary Figure S10 ---\n")
  
  # --- Load data ---
  ratio_data <- readRDS("Data/ratio_data.rds") %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
           Region = factor(Region, levels = region_levels),
           Sex = factor(Sex, levels = c("Female", "Male")),
           Donor = factor(Donor), ImageID = factor(ImageID))
  
  ratio_data_t1d <- readRDS("Data/ratio_data_t1d.rds") %>%
    mutate(Region = factor(Region, levels = region_levels),
           Sex = factor(Sex, levels = c("Female", "Male")),
           Donor = factor(Donor), ImageID = factor(ImageID)) %>%
    droplevels()
  
  y_lab <- paste0(ratio_label, " (%)")
  
  # Helper: bold ratio title with suffix
  bold_title <- function(suffix) paste0(ratio_label_bold, " \u2014 ", suffix)
  
  # --- Load models ---
  m_pow <- readRDS(file.path(model_dir, "ratio_powered_model.rds"))
  m_t1d <- readRDS(file.path(model_dir, "ratio_ao_model.rds"))
  t1d_all <- c("Disease.Duration_c", "age_at_onset_c", "log_Islet.Cells_c")
  
  # ==========================================================================
  # Panel A: Diagnosis × Region — remove subtitle
  # ==========================================================================
  pA <- plot_diagnosis_by_region(ratio_data, m_pow, ratio_var, ratio_label) +
    labs(title = bold_title("Diagnosis \u00D7 Region"), subtitle = NULL) +
    theme_bold_title
  
  # ==========================================================================
  # Panel B: Diagnosis × Sex — remove subtitle
  # ==========================================================================
  pB <- plot_diagnosis_by_sex(ratio_data, m_pow, ratio_var, ratio_label) +
    labs(title = bold_title("Diagnosis \u00D7 Sex"), subtitle = NULL) +
    theme_bold_title
  
  # ==========================================================================
  # Panel C: Diagnosis × Islet Size — move legend, remove points
  # ==========================================================================
  pC <- plot_cont_by_cat(ratio_data, m_pow, ratio_var, ratio_label,
                         "log_Islet.Cells_c", "Diagnosis", colors_diag, "Diagnosis")
  # Remove geom_point layers
  pC$layers <- pC$layers[!sapply(pC$layers, function(l) inherits(l$geom, "GeomPoint"))]
  pC <- pC +
    labs(title = bold_title("Diagnosis \u00D7 Islet Size")) +
    theme_bold_title +
    theme(legend.position = c(0.85, 0.50),
          legend.background = element_rect(fill = "white", color = NA))
  
  # ==========================================================================
  # Panel D: Diagnosis × Age — move legend
  # ==========================================================================
  pD <- plot_cont_by_cat(ratio_data, m_pow, ratio_var, ratio_label,
                         "Age_c", "Diagnosis", colors_diag, "Diagnosis") +
    labs(title = bold_title("Diagnosis \u00D7 Age")) +
    theme_bold_title +
    theme(legend.position = c(0.85, 0.50),
          legend.background = element_rect(fill = "white", color = NA))
  
  # ==========================================================================
  # Panel E: Disease Duration (T1D only)
  # ==========================================================================
  pE <- plot_t1d_slope(ratio_data_t1d, m_t1d, "Disease.Duration_c", "Disease.Duration",
                       NULL, NULL, t1d_all, y_lab, "Disease duration (years)",
                       bold_title("Disease Duration")) +
    theme_bold_title
  
  # ==========================================================================
  # Panel F: Age of Disease Onset (T1D only)
  # ==========================================================================
  pF <- plot_t1d_slope(ratio_data_t1d, m_t1d, "age_at_onset_c", "age_at_onset",
                       NULL, NULL, t1d_all, y_lab, "Age at onset (years)",
                       bold_title("Age of Disease Onset")) +
    theme_bold_title
  
  rm(m_pow, m_t1d); gc()
  
  # ==========================================================================
  # Assemble with patchwork
  # ==========================================================================
  fig <- (pA | pB) / (pC | pD) / (pE | pF) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(output_file, fig, width = 16, height = 18, dpi = 300)
  cat("Saved Supplementary Figure S10 to:", output_file, "\n")
}


# --- Run ---
if (interactive()) {
  cat("\nTo generate all ratio comparison plots, run:\n")
  cat("  generate_all_ratio_plots()\n")
  cat("\nTo generate Supplementary Figure S10, run:\n")
  cat("  generate_supplementary_figure_S10()\n")
}