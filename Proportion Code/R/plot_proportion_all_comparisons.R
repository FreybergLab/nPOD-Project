# =============================================================================
# plot_proportion_all_comparisons.R
# All possible post-hoc comparison plots for Proportion ordbetareg models
# =============================================================================
#
# Three model tiers per cell type (ins, glu, soma, pp):
#
#   POWERED: (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
#            Diagnosis:Region:log_Islet.Cells_c + Diagnosis:Region:Age_c
#            → Main body of paper
#
#   MID:     powered + Diagnosis:Region:Sex + Diagnosis:Sex:{Age,Islet} +
#            Region:Sex:{Age,Islet} + 4-way {Age,Islet}
#            → Supplementary (only interactions NOT in powered)
#
#   FULL:    (5 vars)^5 → Proof of concept (5-way grid)
#
# PLOTS FROM POWERED MODEL (per cell type):
#   Main effects:
#     1.  Diagnosis (violin)
#     2.  Region (violin)
#     3.  Sex (violin)
#     4.  Age (overall slope)
#     5.  Islet size (overall slope)
#   Two-way categorical (6 views):
#     6.  Diagnosis | Region
#     7.  Region | Diagnosis
#     8.  Diagnosis | Sex
#     9.  Sex | Diagnosis
#    10.  Region | Sex
#    11.  Sex | Region
#   Two-way continuous × categorical (6):
#    12.  Age × Diagnosis
#    13.  Age × Region
#    14.  Age × Sex
#    15.  Islet × Diagnosis
#    16.  Islet × Region
#    17.  Islet × Sex
#   Three-way from powered (4 = 2 interactions × 2 views):
#    18.  Diag × Region × Islet — facet Region, color Diag
#    19.  Diag × Region × Islet — facet Diag, color Region
#    20.  Diag × Region × Age — facet Region, color Diag
#    21.  Diag × Region × Age — facet Diag, color Region
#
# PLOTS FROM MID MODEL (per cell type, only NEW interactions):
#    22-23. Diag × Region × Sex (view A + B)
#    24-25. Diag × Sex × Age (view A + B)
#    26-27. Diag × Sex × Islet (view A + B)
#    28-29. Region × Sex × Age (view A + B)
#    30-31. Region × Sex × Islet (view A + B)
#    32.    4-way Diag × Region × Sex × Age
#    33.    4-way Diag × Region × Sex × Islet
#
# PLOTS FROM FULL MODEL (per cell type):
#    34.    5-way: faceted Region × Sex grid, Diag colored, Age on x, at ±1 SD Islet
#    35.    5-way: faceted Region × Sex grid, Diag colored, Islet on x, at ±1 SD Age
#
# Total: 35 plots × 4 cell types = 140 plots
#
# =============================================================================

library(tidyverse)
library(brms)
library(coda)
library(patchwork)
library(cowplot)

source("colors_master.R")

# =============================================================================
# THEME
# =============================================================================

theme_pub <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title       = element_text(size = 11, face = "bold"),
    axis.text        = element_text(size = 9, color = "black"),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 10, color = "gray40"),
    plot.tag         = element_text(size = 18, face = "bold"),
    plot.margin      = margin(8, 12, 8, 8),
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 10),
    strip.text       = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )

region_levels <- c("Head", "Body", "Tail")


# =============================================================================
# CORE POSTERIOR HELPERS
# =============================================================================

get_posterior_fitted <- function(model, newdata, ndraws = NULL) {
  fitted(model, newdata = newdata, re_formula = NA,
         summary = FALSE, ndraws = ndraws)
}

summarize_posterior <- function(x, prob = 0.95) {
  hpd <- coda::HPDinterval(coda::as.mcmc(x), prob = prob)
  tibble(estimate = median(x), lower.HPD = hpd[1], upper.HPD = hpd[2])
}

prob_to_stars <- function(prob) {
  case_when(prob > 0.9995 ~ "***", prob > 0.995 ~ "**", prob > 0.975 ~ "*", TRUE ~ "")
}

# Posterior means for a categorical grid — averaged over non-focal factors
get_categorical_means <- function(model, data, by_vars) {
  grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = factor(c("Head", "Body", "Tail"), levels = region_levels),
    Sex = c("Female", "Male"),
    Age_c = 0,
    log_Islet.Cells_c = 0,
    stringsAsFactors = FALSE
  )
  post <- get_posterior_fitted(model, grid)
  
  # For each unique combo of by_vars, average posterior draws over non-focal
  combos <- distinct(grid[, by_vars, drop = FALSE])
  results <- list()
  for (i in seq_len(nrow(combos))) {
    mask <- rep(TRUE, nrow(grid))
    for (v in by_vars) mask <- mask & grid[[v]] == combos[[v]][i]
    idx <- which(mask)
    draws <- rowMeans(post[, idx, drop = FALSE]) * 100
    summ <- summarize_posterior(draws)
    for (v in by_vars) summ[[v]] <- combos[[v]][i]
    results[[i]] <- summ
  }
  bind_rows(results) %>%
    mutate(estimate = estimate, lower = lower.HPD, upper = upper.HPD) %>%
    select(all_of(by_vars), estimate, lower, upper)
}

# Posterior pairwise contrasts for a 2-level factor | conditioning factor
get_binary_contrast_pd <- function(model, data, focal_var, cond_var = NULL) {
  grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = factor(c("Head", "Body", "Tail"), levels = region_levels),
    Sex = c("Female", "Male"),
    Age_c = 0, log_Islet.Cells_c = 0, stringsAsFactors = FALSE
  )
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
      !!cond_var := cl,
      contrast = paste(focal_levels[1], "-", focal_levels[2]),
      pd = pd, stars = prob_to_stars(pd)
    )
  }
  bind_rows(results)
}

# Region pairwise contrasts | conditioning factor
get_region_contrast_pd <- function(model, data, cond_var = NULL) {
  grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = factor(c("Head", "Body", "Tail"), levels = region_levels),
    Sex = c("Female", "Male"),
    Age_c = 0, log_Islet.Cells_c = 0, stringsAsFactors = FALSE
  )
  post <- get_posterior_fitted(model, grid)
  pairs <- list(c("Head", "Body"), c("Head", "Tail"), c("Body", "Tail"))
  
  if (is.null(cond_var)) {
    results <- list()
    for (pr in pairs) {
      idx_a <- which(grid$Region == pr[1])
      idx_b <- which(grid$Region == pr[2])
      diff <- rowMeans(post[, idx_a, drop = FALSE]) - rowMeans(post[, idx_b, drop = FALSE])
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
      idx_a <- which(mask & grid$Region == pr[1])
      idx_b <- which(mask & grid$Region == pr[2])
      diff <- rowMeans(post[, idx_a, drop = FALSE]) - rowMeans(post[, idx_b, drop = FALSE])
      pd <- pmax(mean(diff > 0), mean(diff < 0))
      results[[length(results) + 1]] <- tibble(
        !!cond_var := cl,
        contrast = paste(pr[1], "-", pr[2]), pd = pd, stars = prob_to_stars(pd))
    }
  }
  bind_rows(results)
}

# Marginalized model predictions for continuous variable
get_continuous_preds <- function(model, data, continuous_var, by_vars, n_points = 50) {
  if (continuous_var == "Age_c") {
    cont_range <- range(data$Age_c)
    cont_seq <- seq(cont_range[1], cont_range[2], length.out = n_points)
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"), Age_c = cont_seq, log_Islet.Cells_c = 0,
      stringsAsFactors = FALSE)
  } else {
    cont_range <- range(data$log_Islet.Cells_c)
    cont_seq <- seq(cont_range[1], cont_range[2], length.out = n_points)
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"), Age_c = 0, log_Islet.Cells_c = cont_seq,
      stringsAsFactors = FALSE)
  }
  post <- get_posterior_fitted(model, grid)
  
  grid %>%
    mutate(cell_id = row_number()) %>%
    group_by(across(all_of(c(by_vars, continuous_var)))) %>%
    summarize(cell_ids = list(cell_id), .groups = "drop") %>%
    mutate(
      post_summary = map(cell_ids, function(ids) {
        draws <- rowMeans(post[, ids, drop = FALSE])
        hpd <- coda::HPDinterval(coda::as.mcmc(draws), prob = 0.95)
        tibble(estimate = median(draws) * 100, lower = hpd[1] * 100, upper = hpd[2] * 100)
      })
    ) %>% unnest(post_summary) %>% select(-cell_ids)
}

# Slope computation: finite differences over ±1 unit, marginalized
compute_slopes <- function(model, continuous_var, by_vars) {
  if (continuous_var == "Age_c") {
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"), Age_c = c(-1, 1), log_Islet.Cells_c = 0,
      stringsAsFactors = FALSE)
  } else {
    grid <- expand.grid(
      Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
      Sex = c("Female", "Male"), Age_c = 0, log_Islet.Cells_c = c(-1, 1),
      stringsAsFactors = FALSE)
  }
  post <- get_posterior_fitted(model, grid)
  combos <- distinct(grid[, by_vars, drop = FALSE])
  
  slopes <- list(); slope_draws <- list()
  for (i in seq_len(nrow(combos))) {
    if (continuous_var == "Age_c") {
      match_lo <- grid$Age_c == -1; match_hi <- grid$Age_c == 1
    } else {
      match_lo <- grid$log_Islet.Cells_c == -1; match_hi <- grid$log_Islet.Cells_c == 1
    }
    for (v in by_vars) { match_lo <- match_lo & grid[[v]] == combos[[v]][i]
    match_hi <- match_hi & grid[[v]] == combos[[v]][i] }
    draws <- (rowMeans(post[, which(match_hi), drop = FALSE]) -
                rowMeans(post[, which(match_lo), drop = FALSE])) / 2
    summ <- summarize_posterior(draws)
    for (v in by_vars) summ[[v]] <- combos[[v]][i]
    summ$pd <- pmax(mean(draws > 0), mean(draws < 0))
    summ$stars <- prob_to_stars(summ$pd)
    key <- paste(sapply(by_vars, function(v) combos[[v]][i]), collapse = "_")
    slopes[[key]] <- summ; slope_draws[[key]] <- draws
  }
  list(slopes = bind_rows(slopes), draws = slope_draws)
}


# =============================================================================
# SAVE HELPER
# =============================================================================

save_plot <- function(plt, filename, output_dir, width = 7, height = 5) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("pdf", "png", "tiff")) {
    ggsave(file.path(output_dir, paste0(filename, ".", ext)), plt,
           width = width, height = height, dpi = 300,
           bg = if (ext != "pdf") "white" else NULL)
  }
  cat("    Saved:", filename, "\n")
}


# #############################################################################
#
#  POWERED MODEL PLOTS  (21 per cell type)
#
# #############################################################################


# --- 1. Main effect: Diagnosis ---
plot_main_diagnosis <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Diagnosis) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, "Diagnosis")
  ctr <- get_binary_contrast_pd(model, data, "Diagnosis")
  
  ggplot() +
    geom_violin(data = dm, aes(x = Diagnosis, y = md, fill = Diagnosis),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm, aes(x = Diagnosis, y = md, color = Diagnosis),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    { if (nchar(ctr$stars) > 0) annotate("text", x = 1.5, y = max(dm$md) * 0.95,
                                         label = ctr$stars, size = 8) } +
    labs(x = NULL, y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Diagnosis")) +
    theme_pub
}

# --- 2. Main effect: Region ---
plot_main_region <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Region) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, "Region")
  ctr <- get_region_contrast_pd(model, data)
  sig <- ctr %>% filter(nchar(stars) > 0) %>%
    mutate(label = paste(contrast, stars))
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
    labs(x = "Pancreatic region", y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Region")) +
    theme_pub
}

# --- 3. Main effect: Sex ---
plot_main_sex <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Sex) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, "Sex")
  ctr <- get_binary_contrast_pd(model, data, "Sex")
  
  ggplot() +
    geom_violin(data = dm, aes(x = Sex, y = md, fill = Sex),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm, aes(x = Sex, y = md, color = Sex),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_sex, guide = "none") +
    scale_color_manual(values = colors_sex, guide = "none") +
    { if (nchar(ctr$stars) > 0) annotate("text", x = 1.5, y = max(dm$md) * 0.95,
                                         label = ctr$stars, size = 8) } +
    labs(x = NULL, y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Sex")) +
    theme_pub
}

# --- 4. Age overall slope ---
plot_age_overall <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  age_mean <- mean(data$Age)
  preds <- get_continuous_preds(model, data, "Age_c", character(0)) %>%
    mutate(Age = Age_c + age_mean)
  dm <- data %>% group_by(Donor, Age) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  
  slope_info <- compute_slopes(model, "Age_c", character(0))
  stars <- slope_info$slopes$stars[1]
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper),
                fill = "grey60", alpha = 0.3) +
    geom_line(data = preds, aes(x = Age, y = estimate), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md), alpha = 0.5, size = 1.8) +
    labs(x = "Age (years)", y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Age (overall)")) +
    theme_pub
  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = min(data$Age) + 2, y = max(dm$md) * 0.95,
                          label = stars, size = 8, hjust = 0)
  plt
}

# --- 5. Islet size overall slope ---
plot_islet_overall <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  islet_mean <- mean(data$log_Islet.Cells)
  preds <- get_continuous_preds(model, data, "log_Islet.Cells_c", character(0)) %>%
    mutate(log_Islet.Cells = log_Islet.Cells_c + islet_mean)
  dm <- data %>% group_by(Donor) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE),
              lic = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  
  slope_info <- compute_slopes(model, "log_Islet.Cells_c", character(0))
  stars <- slope_info$slopes$stars[1]
  
  x_breaks <- seq(2.5, 5.5, by = 0.5)
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = log_Islet.Cells, ymin = lower, ymax = upper),
                fill = "grey60", alpha = 0.3) +
    geom_line(data = preds, aes(x = log_Islet.Cells, y = estimate), linewidth = 0.9) +
    geom_point(data = dm, aes(x = lic, y = md), alpha = 0.5, size = 1.8) +
    scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks))) +
    labs(x = "Islet size (cells)", y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Islet size (overall)")) +
    theme_pub
  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = x_breaks[2], y = max(dm$md) * 0.95,
                          label = stars, size = 8, hjust = 0)
  plt
}


# --- 6. Diagnosis | Region ---
plot_diagnosis_by_region <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Diagnosis, Region) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, c("Diagnosis", "Region"))
  ctr <- get_binary_contrast_pd(model, data, "Diagnosis", "Region")
  y_max <- max(dm$md, na.rm = TRUE) * 1.05
  
  ggplot() +
    geom_violin(data = dm, aes(x = Region, y = md, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm, aes(x = Region, y = md, color = Diagnosis, group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) + scale_color_manual(values = colors_diag) +
    { if (any(nchar(ctr$stars) > 0))
      annotate("text", x = 1:3, y = y_max * c(0.98, 0.93, 0.88),
               label = ctr$stars, size = 8) } +
    labs(x = "Pancreatic region", y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Diagnosis \u00D7 Region"),
         subtitle = "Diagnosis contrast within each Region",
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# --- 7. Region | Diagnosis ---
plot_region_by_diagnosis <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Diagnosis, Region) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, c("Diagnosis", "Region"))
  ctr <- get_region_contrast_pd(model, data, "Diagnosis")
  sig_labels <- ctr %>%
    filter(nchar(stars) > 0) %>%
    mutate(label = paste(contrast, stars)) %>%
    group_by(Diagnosis) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_diag <- tibble(Diagnosis = c("ND", "T1D"))
  sig_labels <- left_join(all_diag, sig_labels, by = "Diagnosis") %>%
    mutate(label = replace_na(label, ""))
  y_max <- max(dm$md, na.rm = TRUE) * 1.05
  
  ggplot() +
    geom_violin(data = dm, aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm, aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ Diagnosis) +
    geom_text(data = sig_labels, aes(x = 2, y = y_max * 0.95, label = label),
              size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Region \u00D7 Diagnosis"),
         subtitle = "Region contrasts within each Diagnosis") +
    theme_pub
}

# --- 8. Diagnosis | Sex ---
plot_diagnosis_by_sex <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Diagnosis, Sex) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, c("Diagnosis", "Sex"))
  ctr <- get_binary_contrast_pd(model, data, "Diagnosis", "Sex")
  y_max <- max(dm$md, na.rm = TRUE) * 1.05
  
  ggplot() +
    geom_violin(data = dm, aes(x = Sex, y = md, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm, aes(x = Sex, y = md, color = Diagnosis, group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) + scale_color_manual(values = colors_diag) +
    annotate("text", x = 1:2, y = y_max * c(0.98, 0.93), label = ctr$stars, size = 8) +
    labs(x = NULL, y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Diagnosis \u00D7 Sex"),
         subtitle = "Diagnosis contrast within each Sex",
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# --- 9. Sex | Diagnosis ---
plot_sex_by_diagnosis <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Diagnosis, Sex) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, c("Diagnosis", "Sex"))
  ctr <- get_binary_contrast_pd(model, data, "Sex", "Diagnosis")
  y_max <- max(dm$md, na.rm = TRUE) * 1.05
  
  ggplot() +
    geom_violin(data = dm, aes(x = Diagnosis, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm, aes(x = Diagnosis, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) + scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:2, y = y_max * c(0.98, 0.90), label = ctr$stars, size = 8) +
    labs(x = NULL, y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Sex \u00D7 Diagnosis"),
         subtitle = "Sex contrast within each Diagnosis",
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# --- 10. Region | Sex ---
plot_region_by_sex <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Region, Sex) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, c("Region", "Sex"))
  ctr <- get_region_contrast_pd(model, data, "Sex")
  sig_labels <- ctr %>%
    filter(nchar(stars) > 0) %>%
    mutate(label = paste(contrast, stars)) %>%
    group_by(Sex) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_sex <- tibble(Sex = c("Female", "Male"))
  sig_labels <- left_join(all_sex, sig_labels, by = "Sex") %>%
    mutate(label = replace_na(label, ""))
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
    labs(x = "Pancreatic region", y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Region \u00D7 Sex"),
         subtitle = "Region contrasts within each Sex") +
    theme_pub
}

# --- 11. Sex | Region ---
plot_sex_by_region <- function(data, model, outcome_var, cell_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  dm <- data %>% group_by(Donor, Region, Sex) %>%
    summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop")
  emm <- get_categorical_means(model, data, c("Region", "Sex"))
  ctr <- get_binary_contrast_pd(model, data, "Sex", "Region")
  y_max <- max(dm$md, na.rm = TRUE) * 1.05
  
  ggplot() +
    geom_violin(data = dm, aes(x = Region, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm, aes(x = Region, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) + scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:3, y = y_max * c(0.98, 0.93, 0.88), label = ctr$stars, size = 8) +
    labs(x = "Pancreatic region", y = paste0("% ", cell_label, " cells"),
         title = paste(cell_label, "— Sex \u00D7 Region"),
         subtitle = "Sex contrast within each Region",
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}


# --- 12-14. Continuous × single categorical (Age) ---
plot_cont_by_cat <- function(data, model, outcome_var, cell_label,
                             cont_var, cat_var, colors, cat_label) {
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  is_age <- cont_var == "Age_c"
  cont_mean <- if (is_age) mean(data$Age) else mean(data$log_Islet.Cells)
  
  preds <- get_continuous_preds(model, data, cont_var, cat_var)
  if (is_age) preds <- preds %>% mutate(x = Age_c + cont_mean)
  else preds <- preds %>% mutate(x = log_Islet.Cells_c + cont_mean)
  
  if (is_age) {
    dm <- data %>% group_by(Donor, !!sym(cat_var), Age) %>%
      summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
      rename(x = Age)
  } else {
    dm <- data %>% group_by(Donor, !!sym(cat_var)) %>%
      summarise(md = mean(outcome_pct, na.rm = TRUE),
                x = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  
  slope_info <- compute_slopes(model, cont_var, cat_var)
  
  x_lab <- if (is_age) "Age (years)" else "Islet size (cells)"
  cont_label <- if (is_age) "Age" else "Islet Size"
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = x, ymin = lower, ymax = upper,
                                  fill = .data[[cat_var]]), alpha = 0.2) +
    geom_line(data = preds, aes(x = x, y = estimate, color = .data[[cat_var]]),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = x, y = md, color = .data[[cat_var]]),
               alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors) + scale_fill_manual(values = colors) +
    labs(x = x_lab, y = paste0("% ", cell_label, " cells"),
         title = paste0(cell_label, " — ", cat_label, " \u00D7 ", cont_label)) +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  
  if (!is_age) {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    plt <- plt + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks)))
  }
  plt
}


# --- 18-21. Three-way: Diag × Region × Continuous (existing pattern) ---
# These reuse the existing plot_celltype_main_figures.R functions.
# We include them here for completeness with the same interface.

plot_3way_facet_region <- function(data, model, outcome_var, cell_label, cont_var, y_max = NULL) {
  is_age <- cont_var == "Age_c"
  cont_mean <- if (is_age) mean(data$Age) else mean(data$log_Islet.Cells)
  cont_range <- if (is_age) range(data$Age_c) else range(data$log_Islet.Cells_c)
  preds <- get_continuous_preds(model, data, cont_var, c("Diagnosis", "Region"))
  if (is_age) preds <- preds %>% mutate(x = Age_c + cont_mean)
  else preds <- preds %>% mutate(x = log_Islet.Cells_c + cont_mean)
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  if (is_age) {
    dm <- data %>% group_by(Donor, Diagnosis, Region, Age) %>%
      summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>%
      rename(x = Age)
  } else {
    dm <- data %>% group_by(Donor, Diagnosis, Region) %>%
      summarise(md = mean(outcome_pct, na.rm = TRUE),
                x = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  if (is.null(y_max)) y_max <- max(c(dm$md, preds$upper), na.rm = TRUE) * 1.2
  
  slope_info <- compute_slopes(model, cont_var, c("Diagnosis", "Region"))
  star_df <- tibble(Region = factor(c("Head", "Body", "Tail"), levels = region_levels)) %>%
    mutate(star_info = map(Region, function(r) {
      t1d_d <- slope_info$draws[[paste("T1D", r, sep = "_")]]
      nd_d <- slope_info$draws[[paste("ND", r, sep = "_")]]
      diff <- t1d_d - nd_d
      pd <- pmax(mean(diff > 0), mean(diff < 0))
      tibble(stars = prob_to_stars(pd))
    })) %>% unnest(star_info) %>% filter(nchar(stars) > 0)
  
  cont_label <- if (is_age) "Age" else "Islet Size"
  x_lab <- if (is_age) "Age (years)" else "Islet size (cells)"
  
  p <- ggplot() +
    geom_ribbon(data = preds, aes(x = x, ymin = lower, ymax = upper, fill = Diagnosis), alpha = 0.25) +
    geom_line(data = preds, aes(x = x, y = estimate, color = Diagnosis), linewidth = 1.3) +
    geom_point(data = dm, aes(x = x, y = md, color = Diagnosis), alpha = 0.5, size = 2) +
    facet_wrap(~ Region, ncol = 3) +
    scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = paste0(cell_label, ": Diag \u00D7 Region \u00D7 ", cont_label),
         subtitle = "Facet by Region, Diagnosis overlaid",
         x = x_lab, y = paste0("% ", cell_label, " cells")) +
    theme_pub + theme(legend.position = "bottom")
  
  if (!is_age) {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks)))
  }
  if (nrow(star_df) > 0)
    p <- p + geom_text(data = star_df, aes(x = if (is_age) 21 else mean(c(2.5,5.5)),
                                           y = y_max * 0.95, label = stars),
                       size = 8, inherit.aes = FALSE)
  p
}

plot_3way_facet_diagnosis <- function(data, model, outcome_var, cell_label, cont_var, y_max = NULL) {
  is_age <- cont_var == "Age_c"
  cont_mean <- if (is_age) mean(data$Age) else mean(data$log_Islet.Cells)
  preds <- get_continuous_preds(model, data, cont_var, c("Diagnosis", "Region"))
  if (is_age) preds <- preds %>% mutate(x = Age_c + cont_mean)
  else preds <- preds %>% mutate(x = log_Islet.Cells_c + cont_mean)
  
  data <- data %>% mutate(outcome_pct = .data[[outcome_var]] * 100)
  if (is_age) {
    dm <- data %>% group_by(Donor, Diagnosis, Region, Age) %>%
      summarise(md = mean(outcome_pct, na.rm = TRUE), .groups = "drop") %>% rename(x = Age)
  } else {
    dm <- data %>% group_by(Donor, Diagnosis, Region) %>%
      summarise(md = mean(outcome_pct, na.rm = TRUE),
                x = mean(log_Islet.Cells, na.rm = TRUE), .groups = "drop")
  }
  if (is.null(y_max)) y_max <- max(c(dm$md, preds$upper), na.rm = TRUE) * 1.25
  
  slope_info <- compute_slopes(model, cont_var, c("Diagnosis", "Region"))
  star_df <- tibble(Diagnosis = c("ND", "T1D")) %>%
    mutate(star_info = map(Diagnosis, function(d) {
      h <- slope_info$draws[[paste(d, "Head", sep = "_")]]
      b <- slope_info$draws[[paste(d, "Body", sep = "_")]]
      t <- slope_info$draws[[paste(d, "Tail", sep = "_")]]
      tibble(comparison = c("H-B", "H-T", "B-T"),
             pd = c(pmax(mean((h-b)>0), mean((h-b)<0)),
                    pmax(mean((h-t)>0), mean((h-t)<0)),
                    pmax(mean((b-t)>0), mean((b-t)<0)))) %>%
        mutate(stars = prob_to_stars(pd)) %>%
        filter(nchar(stars) > 0) %>%
        mutate(label = paste(comparison, stars))
    })) %>% unnest(star_info) %>%
    group_by(Diagnosis) %>%
    summarize(combined = paste(label, collapse = "\n"), .groups = "drop") %>%
    filter(nchar(combined) > 0)
  
  cont_label <- if (is_age) "Age" else "Islet Size"
  x_lab <- if (is_age) "Age (years)" else "Islet size (cells)"
  
  p <- ggplot() +
    geom_ribbon(data = preds, aes(x = x, ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
    geom_line(data = preds, aes(x = x, y = estimate, color = Region), linewidth = 1.3) +
    geom_point(data = dm, aes(x = x, y = md, color = Region), alpha = 0.5, size = 2) +
    facet_wrap(~ Diagnosis, ncol = 2) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = paste0(cell_label, ": Diag \u00D7 Region \u00D7 ", cont_label),
         subtitle = "Facet by Diagnosis, Region overlaid",
         x = x_lab, y = paste0("% ", cell_label, " cells")) +
    theme_pub + theme(legend.position = "bottom")
  
  if (!is_age) {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks)))
  }
  if (nrow(star_df) > 0)
    p <- p + geom_text(data = star_df,
                       aes(x = if (is_age) 21 else mean(c(2.5,5.5)),
                           y = y_max * 0.95, label = combined),
                       size = 3.5, lineheight = 0.85, inherit.aes = FALSE)
  p
}


# #############################################################################
#
#  MID MODEL SUPPLEMENTARY PLOTS  (12 per cell type)
#  These are the EXTRA interactions beyond the powered model.
#  Uses the existing supplementary figure patterns.
#
# #############################################################################

# These are imported directly from plot_celltype_supplementary_figures_no_stars.R
# Functions: plot_diag_region_sex_A/B, plot_diag_sex_cont_A/B,
#            plot_region_sex_cont_A/B, plot_4way
# Already defined in that script — source it or inline below.
# For self-contained use, we source the supplementary script:

# source("plot_celltype_supplementary_figures_no_stars.R")
# The user already has that script; the master function below references it.


# #############################################################################
#
#  FULL MODEL PLOTS  (2 per cell type: 5-way at ±1 SD)
#
# #############################################################################

plot_5way <- function(data, model, outcome_var, cell_label, x_cont_var,
                      cond_cont_var, y_max = NULL) {
  is_age_x <- x_cont_var == "Age_c"
  is_age_cond <- cond_cont_var == "Age_c"
  cont_mean_x <- if (is_age_x) mean(data$Age) else mean(data$log_Islet.Cells)
  cond_sd <- sd(data[[cond_cont_var]], na.rm = TRUE)
  
  preds_list <- list()
  for (lev in c(-1, 1)) {
    x_range <- if (is_age_x) range(data$Age_c) else range(data$log_Islet.Cells_c)
    x_seq <- seq(x_range[1], x_range[2], length.out = 40)
    if (is_age_x) {
      grid <- expand.grid(
        Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
        Sex = c("Female", "Male"), Age_c = x_seq, stringsAsFactors = FALSE)
      grid$log_Islet.Cells_c <- lev * cond_sd
    } else {
      grid <- expand.grid(
        Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
        Sex = c("Female", "Male"), log_Islet.Cells_c = x_seq, stringsAsFactors = FALSE)
      grid$Age_c <- lev * cond_sd
    }
    post <- get_posterior_fitted(model, grid)
    
    # Summarize per Diagnosis × Region × Sex × x
    result <- grid %>%
      mutate(cell_id = row_number()) %>%
      group_by(Diagnosis, Region, Sex, !!sym(x_cont_var)) %>%
      summarize(cell_ids = list(cell_id), .groups = "drop") %>%
      mutate(post_summary = map(cell_ids, function(ids) {
        draws <- if (length(ids) > 1) rowMeans(post[, ids, drop = FALSE]) else post[, ids]
        hpd <- coda::HPDinterval(coda::as.mcmc(draws), prob = 0.95)
        tibble(estimate = median(draws) * 100, lower = hpd[1] * 100, upper = hpd[2] * 100)
      })) %>% unnest(post_summary) %>% select(-cell_ids)
    
    cond_label <- if (is_age_cond) "Age" else "Islet"
    result$cond_level <- factor(
      ifelse(lev == -1,
             paste0(cond_label, " = mean \u2212 1 SD"),
             paste0(cond_label, " = mean + 1 SD")),
      levels = c(paste0(cond_label, " = mean \u2212 1 SD"),
                 paste0(cond_label, " = mean + 1 SD")))
    preds_list[[length(preds_list) + 1]] <- result
  }
  preds <- bind_rows(preds_list)
  if (is_age_x) preds$x <- preds$Age_c + cont_mean_x
  else preds$x <- preds$log_Islet.Cells_c + cont_mean_x
  
  if (is.null(y_max)) y_max <- max(preds$upper, na.rm = TRUE) * 1.1
  x_lab <- if (is_age_x) "Age (years)" else "Islet size (cells)"
  
  p <- ggplot(preds) +
    geom_ribbon(aes(x = x, ymin = lower, ymax = upper, fill = Diagnosis), alpha = 0.2) +
    geom_line(aes(x = x, y = estimate, color = Diagnosis), linewidth = 0.9) +
    facet_grid(Region ~ Sex + cond_level) +
    scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(title = paste0(cell_label, ": 5-way interaction"),
         subtitle = paste("Region \u00D7 Sex grid, Diagnosis overlaid,",
                          if (is_age_cond) "at \u00B11 SD Islet" else "at \u00B11 SD Age"),
         x = x_lab, y = paste0("% ", cell_label, " cells")) +
    theme_pub + theme(legend.position = "bottom")
  
  if (!is_age_x) {
    x_breaks <- seq(2.5, 5.5, by = 0.5)
    p <- p + scale_x_continuous(breaks = x_breaks, labels = round(exp(x_breaks)))
  }
  p
}


# #############################################################################
#
#  MASTER EXECUTION
#
# #############################################################################

generate_all_proportion_plots <- function(output_dir = "Proportion/Graphs/All_Comparisons") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING ALL PROPORTION POST-HOC COMPARISON PLOTS\n")
  cat("=============================================================================\n\n")
  
  # --- Load data ---
  cat("Loading data...\n")
  quad_data <- readRDS("Data/quad_data.rds") %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
      Region    = factor(Region, levels = region_levels),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor), ImageID = factor(ImageID)
    )
  cat("  quad_data:", nrow(quad_data), "rows,", n_distinct(quad_data$Donor), "donors\n")
  
  # --- Load all models ---
  cat("Loading models...\n")
  ins_powered <- readRDS("Proportion/Models/ins_powered_model.rds")
  ins_mid <- readRDS("Proportion/Models/ins_mid_model.rds")
  ins_full <- readRDS("Proportion/Models/ins_full_model.rds")
  
  glu_powered <- readRDS("Proportion/Models/glu_powered_model.rds")
  glu_mid <- readRDS("Proportion/Models/glu_mid_model.rds")
  glu_full <- readRDS("Proportion/Models/glu_full_model.rds")
  
  soma_powered <- readRDS("Proportion/Models/soma_powered_model.rds")
  soma_mid <- readRDS("Proportion/Models/soma_mid_model.rds")
  soma_full <- readRDS("Proportion/Models/soma_full_model.rds")
  
  pp_powered <- readRDS("Proportion/Models/pp_powered_model.rds")
  pp_mid <- readRDS("Proportion/Models/pp_mid_model.rds")
  pp_full <- readRDS("Proportion/Models/pp_full_model.rds")
  
  cell_specs <- list(
    beta  = list(powered = ins_powered,  mid = ins_mid,  full = ins_full,
                 var = "Percent.ins",  label = "Beta Cell",  y_max = 100),
    alpha = list(powered = glu_powered,  mid = glu_mid,  full = glu_full,
                 var = "Percent.glu",  label = "Alpha Cell", y_max = 100),
    delta = list(powered = soma_powered, mid = soma_mid, full = soma_full,
                 var = "Percent.soma", label = "Delta Cell", y_max = 60),
    pp    = list(powered = pp_powered,   mid = pp_mid,   full = pp_full,
                 var = "Percent.PP",   label = "PP Cell",    y_max = 90)
  )
  
  # Mid model supplementary functions (3-way DRS, DSAge, DSIslet, RSAge, RSIslet, 4-way)
  # are defined in plot_celltype_supplementary_figures_no_stars.R
  cat("Sourcing supplementary figure functions...\n")
  source("Proportion/R/plot_celltype_supplementary_figures_no_stars.R")
  
  for (ct in names(cell_specs)) {
    spec <- cell_specs[[ct]]
    v <- spec$var; lab <- spec$label; ym <- spec$y_max
    cat("\n=== ", toupper(lab), " ===\n")
    
    dir_pow  <- file.path(output_dir, ct, "Powered")
    dir_mid  <- file.path(output_dir, ct, "Mid")
    dir_full <- file.path(output_dir, ct, "Full")
    
    # ====================== POWERED ======================
    cat("  --- Powered (21 plots) ---\n")
    m <- spec$powered
    
    save_plot(plot_main_diagnosis(quad_data, m, v, lab), "01_main_diagnosis", dir_pow)
    save_plot(plot_main_region(quad_data, m, v, lab), "02_main_region", dir_pow)
    save_plot(plot_main_sex(quad_data, m, v, lab), "03_main_sex", dir_pow)
    save_plot(plot_age_overall(quad_data, m, v, lab), "04_age_overall", dir_pow)
    save_plot(plot_islet_overall(quad_data, m, v, lab), "05_islet_overall", dir_pow)
    
    save_plot(plot_diagnosis_by_region(quad_data, m, v, lab), "06_diagnosis_by_region", dir_pow, 8)
    save_plot(plot_region_by_diagnosis(quad_data, m, v, lab), "07_region_by_diagnosis", dir_pow, 9, 5)
    save_plot(plot_diagnosis_by_sex(quad_data, m, v, lab), "08_diagnosis_by_sex", dir_pow)
    save_plot(plot_sex_by_diagnosis(quad_data, m, v, lab), "09_sex_by_diagnosis", dir_pow)
    save_plot(plot_region_by_sex(quad_data, m, v, lab), "10_region_by_sex", dir_pow, 9, 5)
    save_plot(plot_sex_by_region(quad_data, m, v, lab), "11_sex_by_region", dir_pow, 8)
    
    save_plot(plot_cont_by_cat(quad_data, m, v, lab, "Age_c", "Diagnosis", colors_diag, "Diagnosis"),
              "12_age_by_diagnosis", dir_pow)
    save_plot(plot_cont_by_cat(quad_data, m, v, lab, "Age_c", "Region", colors_region, "Region"),
              "13_age_by_region", dir_pow)
    save_plot(plot_cont_by_cat(quad_data, m, v, lab, "Age_c", "Sex", colors_sex, "Sex"),
              "14_age_by_sex", dir_pow)
    save_plot(plot_cont_by_cat(quad_data, m, v, lab, "log_Islet.Cells_c", "Diagnosis", colors_diag, "Diagnosis"),
              "15_islet_by_diagnosis", dir_pow)
    save_plot(plot_cont_by_cat(quad_data, m, v, lab, "log_Islet.Cells_c", "Region", colors_region, "Region"),
              "16_islet_by_region", dir_pow)
    save_plot(plot_cont_by_cat(quad_data, m, v, lab, "log_Islet.Cells_c", "Sex", colors_sex, "Sex"),
              "17_islet_by_sex", dir_pow)
    
    save_plot(plot_3way_facet_region(quad_data, m, v, lab, "log_Islet.Cells_c", ym),
              "18_3way_islet_facet_region", dir_pow, 14, 5)
    save_plot(plot_3way_facet_diagnosis(quad_data, m, v, lab, "log_Islet.Cells_c", ym),
              "19_3way_islet_facet_diagnosis", dir_pow, 12, 6)
    save_plot(plot_3way_facet_region(quad_data, m, v, lab, "Age_c", ym),
              "20_3way_age_facet_region", dir_pow, 14, 5)
    save_plot(plot_3way_facet_diagnosis(quad_data, m, v, lab, "Age_c", ym),
              "21_3way_age_facet_diagnosis", dir_pow, 12, 6)
    
    # ====================== MID (supplementary extras) ======================
    cat("  --- Mid (12 supplementary plots) ---\n")
    m2 <- spec$mid
    
    # Source supplementary functions if not already loaded
    # These use the exact patterns from plot_celltype_supplementary_figures_no_stars.R
    save_plot(plot_diag_region_sex_A(quad_data, m2, v, lab, ym),
              "22_diag_region_sex_byRegion", dir_mid, 14, 6)
    save_plot(plot_diag_region_sex_B(quad_data, m2, v, lab, ym),
              "23_diag_region_sex_byDiag", dir_mid, 12, 6)
    save_plot(plot_diag_sex_cont_A(quad_data, m2, v, lab, "Age_c", ym),
              "24_diag_sex_age_bySex", dir_mid, 12, 5)
    save_plot(plot_diag_sex_cont_B(quad_data, m2, v, lab, "Age_c", ym),
              "25_diag_sex_age_byDiag", dir_mid, 12, 5)
    save_plot(plot_diag_sex_cont_A(quad_data, m2, v, lab, "log_Islet.Cells_c", ym),
              "26_diag_sex_islet_bySex", dir_mid, 12, 5)
    save_plot(plot_diag_sex_cont_B(quad_data, m2, v, lab, "log_Islet.Cells_c", ym),
              "27_diag_sex_islet_byDiag", dir_mid, 12, 5)
    save_plot(plot_region_sex_cont_A(quad_data, m2, v, lab, "Age_c", ym),
              "28_region_sex_age_bySex", dir_mid, 12, 6)
    save_plot(plot_region_sex_cont_B(quad_data, m2, v, lab, "Age_c", ym),
              "29_region_sex_age_byRegion", dir_mid, 14, 5)
    save_plot(plot_region_sex_cont_A(quad_data, m2, v, lab, "log_Islet.Cells_c", ym),
              "30_region_sex_islet_bySex", dir_mid, 12, 6)
    save_plot(plot_region_sex_cont_B(quad_data, m2, v, lab, "log_Islet.Cells_c", ym),
              "31_region_sex_islet_byRegion", dir_mid, 14, 5)
    save_plot(plot_4way(quad_data, m2, v, lab, "Age_c", ym),
              "32_4way_age", dir_mid, 10, 10)
    save_plot(plot_4way(quad_data, m2, v, lab, "log_Islet.Cells_c", ym),
              "33_4way_islet", dir_mid, 10, 10)
    
    # ====================== FULL (proof of concept) ======================
    cat("  --- Full (2 proof-of-concept plots) ---\n")
    m3 <- spec$full
    save_plot(plot_5way(quad_data, m3, v, lab, "Age_c", "log_Islet.Cells_c", ym),
              "34_5way_age_x_islet_cond", dir_full, 16, 10)
    save_plot(plot_5way(quad_data, m3, v, lab, "log_Islet.Cells_c", "Age_c", ym),
              "35_5way_islet_x_age_cond", dir_full, 16, 10)
    
    rm(m, m2, m3); gc()
  }
  
  cat("\n=============================================================================\n")
  cat("COMPLETE! 35 plots × 4 cell types = 140 plots\n")
  cat("Saved to", output_dir, "\n")
  cat("=============================================================================\n")
}


# --- Run ---
if (interactive()) {
  cat("\nTo generate all proportion comparison plots, run:\n")
  cat("  source('plot_celltype_supplementary_figures_no_stars.R')  # for mid model functions\n")
  cat("  generate_all_proportion_plots()\n\n")
  cat("Or call individual functions, e.g.:\n")
  cat("  plot_main_diagnosis(quad_data, ins_powered, 'Percent.ins', 'Beta Cell')\n")
}