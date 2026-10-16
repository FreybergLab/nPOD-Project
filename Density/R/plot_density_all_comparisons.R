# =============================================================================
# plot_density_all_comparisons.R
# All possible post-hoc comparison plots for Density models
# =============================================================================
#
# Generates every meaningful post-hoc comparison plot:
#
# DIAGNOSIS MODEL (density_diag):
#   Formula: islet_density ~ (Diagnosis + Region + Sex + Age_c)^2 +
#            Diagnosis:Region:Age_c + (1|Donor)
#   Family: Gamma(log)
#
#   Main effects (3):
#     1.  Diagnosis
#     2.  Region
#     3.  Sex
#
#   Two-way categorical (6 views):
#     4.  Diagnosis | Region  (Diagnosis within each Region)
#     5.  Region | Diagnosis  (Region within each Diagnosis)
#     6.  Diagnosis | Sex     (Diagnosis within each Sex)
#     7.  Sex | Diagnosis     (Sex within each Diagnosis)
#     8.  Region | Sex        (Region within each Sex)
#     9.  Sex | Region        (Sex within each Region)
#
#   Age slopes (4):
#    10.  Age overall
#    11.  Age × Diagnosis
#    12.  Age × Region
#    13.  Age × Sex
#
#   Three-way Diagnosis × Region × Age (2):
#    14.  Facet by Region,    colored by Diagnosis
#    15.  Facet by Diagnosis, colored by Region
#
# T1D MODEL (density_t1d):
#   Formula: islet_density ~ (Disease.Duration_c + Region + Sex + age_at_onset_c)^2 +
#            Disease.Duration_c:Region:age_at_onset_c + (1|Donor)
#   Family: Gamma(log)
#
#   Main effects (2):
#    16.  Region  (T1D)
#    17.  Sex     (T1D)
#
#   Two-way categorical (2):
#    18.  Region | Sex  (T1D)
#    19.  Sex | Region  (T1D)
#
#   DD spline (1):
#    20.  DD spline (marginalized, with CI ribbon)
#
#   DD slopes (3):
#    21.  DD overall
#    22.  DD × Region
#    23.  DD × Sex
#
#   AO slopes (3):
#    24.  AO overall
#    25.  AO × Region
#    26.  AO × Sex
#
#   DD × AO interaction (2):
#    27.  DD slopes at low/high AO
#    28.  AO slopes at low/high DD
#
#   Three-way DD × Region × AO (2):
#    29.  DD by Region, faceted at AO ± 1 SD
#    30.  AO by Region, faceted at DD ± 1 SD
#
# Total: 30 individual plots
#
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(patchwork)
library(splines)

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

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

p_to_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
}

find_ci <- function(df) {
  nms <- names(df)
  lower <- upper <- NULL
  for (lc in c("asymp.LCL", "lower.CL", "lower.HPD", "lower.CI"))
    if (lc %in% nms) { lower <- lc; break }
  for (uc in c("asymp.UCL", "upper.CL", "upper.HPD", "upper.CI"))
    if (uc %in% nms) { upper <- uc; break }
  list(lower = lower, upper = upper)
}

# Marginalized predictions on response scale for Gamma(log) models
get_marginalized_preds <- function(model, data, focal_grid, avg_vars = NULL) {
  model_vars <- setdiff(all.vars(formula(model)[-2]), c("Donor", "ImageID"))
  focal_vars <- names(focal_grid)
  if (is.null(avg_vars)) {
    avg_vars <- setdiff(intersect(model_vars, c("Diagnosis", "Region", "Sex")), focal_vars)
  }
  avg_levels <- list()
  for (v in avg_vars) avg_levels[[v]] <- if (is.factor(data[[v]])) levels(data[[v]]) else na.omit(unique(data[[v]]))
  avg_grid <- if (length(avg_levels) > 0) {
    expand.grid(avg_levels, stringsAsFactors = FALSE)
  } else {
    data.frame(.d = 1)
  }
  all_preds <- all_se <- list()
  for (i in seq_len(nrow(avg_grid))) {
    fg <- focal_grid
    if (length(avg_vars) > 0)
      for (v in avg_vars) fg[[v]] <- avg_grid[[v]][i]
    for (v in model_vars)
      if (!v %in% names(fg))
        fg[[v]] <- if (is.numeric(data[[v]])) 0 else levels(data[[v]])[1]
    for (v in names(fg))
      if (is.factor(data[[v]])) fg[[v]] <- factor(fg[[v]], levels = levels(data[[v]]))
    fg$Donor <- NA
    p <- predict(model, newdata = fg, type = "link", re.form = NA,
                 se.fit = TRUE, allow.new.levels = TRUE)
    all_preds[[i]] <- p$fit; all_se[[i]] <- p$se.fit
  }
  ap <- colMeans(do.call(rbind, all_preds))
  as <- sqrt(colMeans(do.call(rbind, all_se)^2))
  focal_grid$estimate <- exp(ap)
  focal_grid$lower    <- exp(ap - 1.96 * as)
  focal_grid$upper    <- exp(ap + 1.96 * as)
  focal_grid
}

# Save a single plot in 3 formats
save_plot <- function(plt, filename, output_dir, width = 7, height = 5) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("pdf", "png", "tiff")) {
    ggsave(file.path(output_dir, paste0(filename, ".", ext)), plt,
           width = width, height = height, dpi = 300,
           bg = if (ext != "pdf") "white" else NULL)
  }
  cat("  Saved:", filename, "\n")
}


# #############################################################################
#
#  SECTION 1:  DIAGNOSIS MODEL — density_diag
#
# #############################################################################


# =============================================================================
# 1. MAIN EFFECT: Diagnosis
# =============================================================================

plot_diag_main_diagnosis <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Diagnosis) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Diagnosis, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Diagnosis), adjust = "none", reverse = TRUE))
  stars <- p_to_stars(prs$p.value[1])

  ggplot() +
    geom_violin(data = donor_means, aes(x = Diagnosis, y = md, fill = Diagnosis),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means, aes(x = Diagnosis, y = md, color = Diagnosis),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    { if (nchar(stars) > 0) annotate("text", x = 1.5, y = max(donor_means$md) * 0.95,
                                      label = stars, size = 8) } +
    labs(x = NULL, y = "Islet density (islets/mm\u00B2)",
         title = "Diagnosis") +
    theme_pub
}

# =============================================================================
# 2. MAIN EFFECT: Region
# =============================================================================

plot_diag_main_region <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Region) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Region, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Region), adjust = "tukey"))

  sig_labels <- prs %>% mutate(label = paste(contrast, p_to_stars(p.value))) %>%
    filter(p.value < 0.05)
  sig_text <- paste(sig_labels$label, collapse = "\n")

  ggplot() +
    geom_violin(data = donor_means, aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means, aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    { if (nchar(sig_text) > 0) annotate("text", x = 2, y = max(donor_means$md) * 0.95,
                                         label = sig_text, size = 4, lineheight = 0.9) } +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Region") +
    theme_pub
}

# =============================================================================
# 3. MAIN EFFECT: Sex
# =============================================================================

plot_diag_main_sex <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Sex, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Sex), adjust = "none"))
  stars <- p_to_stars(prs$p.value[1])

  ggplot() +
    geom_violin(data = donor_means, aes(x = Sex, y = md, fill = Sex),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means, aes(x = Sex, y = md, color = Sex),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_sex, guide = "none") +
    scale_color_manual(values = colors_sex, guide = "none") +
    { if (nchar(stars) > 0) annotate("text", x = 1.5, y = max(donor_means$md) * 0.95,
                                      label = stars, size = 8) } +
    labs(x = NULL, y = "Islet density (islets/mm\u00B2)",
         title = "Sex") +
    theme_pub
}


# =============================================================================
# 4. Diagnosis | Region  (Diagnosis contrast within each Region)
# =============================================================================

plot_diag_diagnosis_by_region <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Diagnosis | Region, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Diagnosis | Region), adjust = "none", reverse = TRUE))
  prs$stars <- p_to_stars(prs$p.value)
  y_max <- max(donor_means$md, na.rm = TRUE)

  ggplot() +
    geom_violin(data = donor_means,
                aes(x = Region, y = md, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = donor_means,
                aes(x = Region, y = md, color = Diagnosis, group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) +
    scale_color_manual(values = colors_diag) +
    annotate("text", x = 1:3, y = y_max * c(0.98, 0.93, 0.88),
             label = prs$stars, size = 8) +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Diagnosis \u00D7 Region",
         subtitle = "Diagnosis contrast within each Region",
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# 5. Region | Diagnosis  (Region contrasts within each Diagnosis)
# =============================================================================

plot_diag_region_by_diagnosis <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Region) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Region | Diagnosis, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Region | Diagnosis), adjust = "tukey"))

  sig_labels <- prs %>%
    mutate(label = paste(contrast, p_to_stars(p.value))) %>%
    filter(p.value < 0.05) %>%
    group_by(Diagnosis) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  # Fill in all Diagnosis levels
  all_diag <- data.frame(Diagnosis = factor(c("ND", "T1D"), levels = levels(data$Diagnosis)))
  sig_labels <- left_join(all_diag, sig_labels, by = "Diagnosis") %>%
    mutate(label = ifelse(is.na(label), "", label))

  y_max <- max(donor_means$md, na.rm = TRUE) * 1.05

  ggplot() +
    geom_violin(data = donor_means,
                aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means,
                aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ Diagnosis) +
    geom_text(data = sig_labels,
              aes(x = 2, y = y_max * 0.95, label = label),
              size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Region \u00D7 Diagnosis",
         subtitle = "Region contrasts within each Diagnosis") +
    theme_pub
}

# =============================================================================
# 6. Diagnosis | Sex  (Diagnosis contrast within each Sex)
# =============================================================================

plot_diag_diagnosis_by_sex <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Diagnosis | Sex, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Diagnosis | Sex), adjust = "none", reverse = TRUE))
  prs$stars <- p_to_stars(prs$p.value)
  y_max <- max(donor_means$md, na.rm = TRUE)

  ggplot() +
    geom_violin(data = donor_means,
                aes(x = Sex, y = md, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = donor_means,
                aes(x = Sex, y = md, color = Diagnosis, group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) +
    scale_color_manual(values = colors_diag) +
    annotate("text", x = 1:2, y = y_max * c(0.98, 0.93),
             label = prs$stars, size = 8) +
    labs(x = NULL, y = "Islet density (islets/mm\u00B2)",
         title = "Diagnosis \u00D7 Sex",
         subtitle = "Diagnosis contrast within each Sex",
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# 7. Sex | Diagnosis  (Sex contrast within each Diagnosis)
# =============================================================================

plot_diag_sex_by_diagnosis <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Diagnosis, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Sex | Diagnosis, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Sex | Diagnosis), adjust = "none"))
  prs$stars <- p_to_stars(prs$p.value)

  y_max <- max(donor_means$md, na.rm = TRUE)

  ggplot() +
    geom_violin(data = donor_means,
                aes(x = Diagnosis, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = donor_means,
                aes(x = Diagnosis, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) +
    scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:2, y = y_max * c(0.98, 0.90),
             label = prs$stars, size = 8) +
    labs(x = NULL, y = "Islet density (islets/mm\u00B2)",
         title = "Sex \u00D7 Diagnosis",
         subtitle = "Sex contrast within each Diagnosis",
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# 8. Region | Sex  (Region contrasts within each Sex)
# =============================================================================

plot_diag_region_by_sex <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Region, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Region | Sex, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Region | Sex), adjust = "tukey"))

  sig_labels <- prs %>%
    mutate(label = paste(contrast, p_to_stars(p.value))) %>%
    filter(p.value < 0.05) %>%
    group_by(Sex) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_sex <- data.frame(Sex = factor(c("Female", "Male"), levels = levels(data$Sex)))
  sig_labels <- left_join(all_sex, sig_labels, by = "Sex") %>%
    mutate(label = ifelse(is.na(label), "", label))

  y_max <- max(donor_means$md, na.rm = TRUE) * 1.05

  ggplot() +
    geom_violin(data = donor_means,
                aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means,
                aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ Sex) +
    geom_text(data = sig_labels,
              aes(x = 2, y = y_max * 0.95, label = label),
              size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Region \u00D7 Sex",
         subtitle = "Region contrasts within each Sex") +
    theme_pub
}

# =============================================================================
# 9. Sex | Region  (Sex contrast within each Region)
# =============================================================================

plot_diag_sex_by_region <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Region, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Sex | Region, type = "response"))
  prs <- as.data.frame(pairs(emmeans(model, ~ Sex | Region), adjust = "none"))
  prs$stars <- p_to_stars(prs$p.value)
  y_max <- max(donor_means$md, na.rm = TRUE)

  ggplot() +
    geom_violin(data = donor_means,
                aes(x = Region, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = donor_means,
                aes(x = Region, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) +
    scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:3, y = y_max * c(0.98, 0.93, 0.88),
             label = prs$stars, size = 8) +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Sex \u00D7 Region",
         subtitle = "Sex contrast within each Region",
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}


# =============================================================================
# 10. Age — overall marginalized slope
# =============================================================================

plot_diag_age_overall <- function(data, model) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- data.frame(Age_c = seq(ar[1] - ma, ar[2] - ma, length.out = 100)) %>%
    mutate(Age = Age_c + ma)
  preds <- get_marginalized_preds(model, data, fg,
                                  avg_vars = c("Diagnosis", "Region", "Sex"))
  dm <- data %>% group_by(Donor, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trend <- emtrends(model, ~ 1, var = "Age_c", data = data)
  p_val <- as.data.frame(test(trend))$p.value[1]
  stars <- p_to_stars(p_val)

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper),
                fill = "grey60", alpha = 0.3) +
    geom_line(data = preds, aes(x = Age, y = estimate), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md), alpha = 0.5, size = 1.8) +
    labs(x = "Age (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Age (overall)") +
    theme_pub

  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = ar[1] + 2, y = max(dm$md) * 0.95,
                           label = stars, size = 8, hjust = 0)
  plt
}

# =============================================================================
# 11. Age × Diagnosis
# =============================================================================

plot_diag_age_by_diagnosis <- function(data, model) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(
    Age_c = seq(ar[1] - ma, ar[2] - ma, length.out = 100),
    Diagnosis = factor(c("ND", "T1D"), levels = levels(data$Diagnosis)),
    stringsAsFactors = FALSE
  ) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = c("Region", "Sex"))
  dm <- data %>% group_by(Donor, Diagnosis, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trends <- emtrends(model, pairwise ~ Diagnosis, var = "Age_c", data = data)
  p_val <- as.data.frame(trends$contrasts)$p.value[1]
  stars <- p_to_stars(p_val)

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Diagnosis),
                alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Diagnosis), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Diagnosis), alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
    labs(x = "Age (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Diagnosis \u00D7 Age") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = ar[1] + 2, y = max(dm$md) * 0.95,
                           label = stars, size = 8, hjust = 0)
  plt
}

# =============================================================================
# 12. Age × Region
# =============================================================================

plot_diag_age_by_region <- function(data, model) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(
    Age_c = seq(ar[1] - ma, ar[2] - ma, length.out = 100),
    Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
    stringsAsFactors = FALSE
  ) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = c("Diagnosis", "Sex"))
  dm <- data %>% group_by(Donor, Region, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trends <- emtrends(model, pairwise ~ Region, var = "Age_c", data = data)
  pairs_df <- as.data.frame(trends$contrasts)
  sig_pairs <- pairs_df %>%
    filter(p.value < 0.05) %>%
    mutate(t = sprintf("%s %s", contrast, p_to_stars(p.value)))
  sig_text <- paste(sig_pairs$t, collapse = "\n")

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Region),
                alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Region), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Region), alpha = 0.4, size = 1.5) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    labs(x = "Age (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Region \u00D7 Age") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(sig_text) > 0)
    plt <- plt + annotate("text", x = ar[1] + 2, y = max(dm$md) * 0.95,
                           label = sig_text, size = 4, hjust = 0, lineheight = 0.9)
  plt
}

# =============================================================================
# 13. Age × Sex
# =============================================================================

plot_diag_age_by_sex <- function(data, model) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(
    Age_c = seq(ar[1] - ma, ar[2] - ma, length.out = 100),
    Sex = factor(c("Female", "Male"), levels = levels(data$Sex)),
    stringsAsFactors = FALSE
  ) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = c("Diagnosis", "Region"))
  dm <- data %>% group_by(Donor, Sex, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trends <- emtrends(model, pairwise ~ Sex, var = "Age_c", data = data)
  p_val <- as.data.frame(trends$contrasts)$p.value[1]
  stars <- p_to_stars(p_val)

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Sex),
                alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Sex), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Sex), alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors_sex) + scale_fill_manual(values = colors_sex) +
    labs(x = "Age (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Sex \u00D7 Age") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = ar[1] + 2, y = max(dm$md) * 0.95,
                           label = stars, size = 8, hjust = 0)
  plt
}


# =============================================================================
# 14. Diagnosis × Region × Age — faceted by Region, colored by Diagnosis
# =============================================================================

plot_diag_3way_facet_region <- function(data, model) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(
    Age_c = seq(ar[1] - ma, ar[2] - ma, length.out = 100),
    Diagnosis = factor(c("ND", "T1D"), levels = levels(data$Diagnosis)),
    Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
    stringsAsFactors = FALSE
  ) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = "Sex")
  dm <- data %>% group_by(Donor, Diagnosis, Region, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  # Diagnosis contrasts on age slopes within each Region
  pairs_df <- as.data.frame(emtrends(model, pairwise ~ Diagnosis | Region,
                                      var = "Age_c", data = data)$contrasts)
  sig_labels <- pairs_df %>%
    mutate(label = p_to_stars(p.value)) %>%
    filter(nchar(label) > 0)
  # Match to Region factor
  all_reg <- data.frame(Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)))
  sig_labels <- left_join(all_reg, sig_labels, by = "Region") %>%
    mutate(label = ifelse(is.na(label), "", label))

  y_max <- max(c(dm$md, preds$upper), na.rm = TRUE) * 1.1

  ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Diagnosis),
                alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Diagnosis), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Diagnosis), alpha = 0.4, size = 1.5) +
    facet_wrap(~ Region, ncol = 3) +
    geom_text(data = sig_labels,
              aes(x = ar[1] + 2, y = y_max * 0.95, label = label),
              size = 6, hjust = 0) +
    scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(x = "Age (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Diagnosis \u00D7 Region \u00D7 Age",
         subtitle = "Diagnosis contrasts on age slopes within each Region") +
    theme_pub + theme(legend.position = "bottom")
}

# =============================================================================
# 15. Diagnosis × Region × Age — faceted by Diagnosis, colored by Region
# =============================================================================

plot_diag_3way_facet_diagnosis <- function(data, model) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(
    Age_c = seq(ar[1] - ma, ar[2] - ma, length.out = 100),
    Diagnosis = factor(c("ND", "T1D"), levels = levels(data$Diagnosis)),
    Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
    stringsAsFactors = FALSE
  ) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = "Sex")
  dm <- data %>% group_by(Donor, Diagnosis, Region, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  # Region contrasts on age slopes within each Diagnosis
  pairs_df <- as.data.frame(emtrends(model, pairwise ~ Region | Diagnosis,
                                      var = "Age_c", data = data)$contrasts)
  sig_labels <- pairs_df %>%
    filter(p.value < 0.05) %>%
    mutate(label = sprintf("%s %s", contrast, p_to_stars(p.value))) %>%
    group_by(Diagnosis) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_diag <- data.frame(Diagnosis = factor(c("ND", "T1D"), levels = levels(data$Diagnosis)))
  sig_labels <- left_join(all_diag, sig_labels, by = "Diagnosis") %>%
    mutate(label = ifelse(is.na(label), "", label))

  y_max <- max(c(dm$md, preds$upper), na.rm = TRUE) * 1.1

  ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Region),
                alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Region), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Region), alpha = 0.4, size = 1.5) +
    facet_wrap(~ Diagnosis, ncol = 2) +
    geom_text(data = sig_labels,
              aes(x = ar[1] + 2, y = y_max * 0.95, label = label),
              size = 4, hjust = 0, vjust = 1, lineheight = 0.9) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(x = "Age (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Diagnosis \u00D7 Region \u00D7 Age",
         subtitle = "Region contrasts on age slopes within each Diagnosis") +
    theme_pub + theme(legend.position = "bottom")
}


# #############################################################################
#
#  SECTION 2:  T1D MODEL — density_t1d
#
# #############################################################################


# =============================================================================
# 16. Region (T1D)
# =============================================================================

plot_t1d_main_region <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Region) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Region, type = "response", data = data))
  prs <- as.data.frame(pairs(emmeans(model, ~ Region, data = data), adjust = "tukey"))

  sig_labels <- prs %>% mutate(label = paste(contrast, p_to_stars(p.value))) %>%
    filter(p.value < 0.05)
  sig_text <- paste(sig_labels$label, collapse = "\n")

  ggplot() +
    geom_violin(data = donor_means, aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means, aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    { if (nchar(sig_text) > 0) annotate("text", x = 2, y = max(donor_means$md) * 0.95,
                                         label = sig_text, size = 4, lineheight = 0.9) } +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Region (T1D)") +
    theme_pub
}

# =============================================================================
# 17. Sex (T1D)
# =============================================================================

plot_t1d_main_sex <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Sex, type = "response", data = data))
  prs <- as.data.frame(pairs(emmeans(model, ~ Sex, data = data), adjust = "none"))
  stars <- p_to_stars(prs$p.value[1])

  ggplot() +
    geom_violin(data = donor_means, aes(x = Sex, y = md, fill = Sex),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means, aes(x = Sex, y = md, color = Sex),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_sex, guide = "none") +
    scale_color_manual(values = colors_sex, guide = "none") +
    { if (nchar(stars) > 0) annotate("text", x = 1.5, y = max(donor_means$md) * 0.95,
                                      label = stars, size = 8) } +
    labs(x = NULL, y = "Islet density (islets/mm\u00B2)",
         title = "Sex (T1D)") +
    theme_pub
}

# =============================================================================
# 18. Region | Sex  (T1D)
# =============================================================================

plot_t1d_region_by_sex <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Region, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Region | Sex, type = "response", data = data))
  prs <- as.data.frame(pairs(emmeans(model, ~ Region | Sex, data = data), adjust = "tukey"))

  sig_labels <- prs %>%
    mutate(label = paste(contrast, p_to_stars(p.value))) %>%
    filter(p.value < 0.05) %>%
    group_by(Sex) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_sex <- data.frame(Sex = factor(c("Female", "Male"), levels = levels(data$Sex)))
  sig_labels <- left_join(all_sex, sig_labels, by = "Sex") %>%
    mutate(label = ifelse(is.na(label), "", label))

  y_max <- max(donor_means$md, na.rm = TRUE) * 1.05

  ggplot() +
    geom_violin(data = donor_means, aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means, aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ Sex) +
    geom_text(data = sig_labels,
              aes(x = 2, y = y_max * 0.95, label = label),
              size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Region \u00D7 Sex (T1D)",
         subtitle = "Region contrasts within each Sex") +
    theme_pub
}

# =============================================================================
# 19. Sex | Region  (T1D)
# =============================================================================

plot_t1d_sex_by_region <- function(data, model) {
  donor_means <- data %>%
    group_by(Donor, Region, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  emm <- as.data.frame(emmeans(model, ~ Sex | Region, type = "response", data = data))
  prs <- as.data.frame(pairs(emmeans(model, ~ Sex | Region, data = data), adjust = "none"))
  prs$stars <- p_to_stars(prs$p.value)
  y_max <- max(donor_means$md, na.rm = TRUE)

  ggplot() +
    geom_violin(data = donor_means,
                aes(x = Region, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = donor_means,
                aes(x = Region, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) +
    scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:3, y = y_max * c(0.98, 0.93, 0.88),
             label = prs$stars, size = 8) +
    labs(x = "Pancreatic region", y = "Islet density (islets/mm\u00B2)",
         title = "Sex \u00D7 Region (T1D)",
         subtitle = "Sex contrast within each Region",
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.90, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}


# =============================================================================
# 20. DD spline (marginalized, with CI ribbon)
# =============================================================================

plot_t1d_dd_spline <- function(data, model_spline) {
  data <- data %>% mutate(Region = factor(Region), Sex = factor(Sex))
  donor_means <- data %>%
    group_by(Donor, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  mean_dur <- mean(data$Disease.Duration, na.rm = TRUE)
  dur_raw <- seq(0.1, max(data$Disease.Duration, na.rm = TRUE), length.out = 100)
  dur_c <- dur_raw - mean_dur

  reg_lvls <- levels(data$Region); sex_lvls <- levels(data$Sex)
  avg_levels <- expand.grid(Region = reg_lvls, Sex = sex_lvls, stringsAsFactors = FALSE)
  all_preds <- all_se <- list()
  for (i in seq_len(nrow(avg_levels))) {
    nd <- data.frame(Disease.Duration_c = dur_c, Age_c = 0,
                     Region = factor(avg_levels$Region[i], levels = reg_lvls),
                     Sex = factor(avg_levels$Sex[i], levels = sex_lvls),
                     Donor = NA)
    p <- predict(model_spline, newdata = nd, type = "link", re.form = NA,
                 se.fit = TRUE, allow.new.levels = TRUE)
    all_preds[[i]] <- p$fit; all_se[[i]] <- p$se.fit
  }
  avg_p <- colMeans(do.call(rbind, all_preds))
  avg_s <- sqrt(colMeans(do.call(rbind, all_se)^2))
  pred_df <- data.frame(x = dur_raw, est = exp(avg_p),
                         lo = exp(avg_p - 1.96 * avg_s),
                         hi = exp(avg_p + 1.96 * avg_s))

  ggplot() +
    geom_ribbon(data = pred_df, aes(x = x, ymin = lo, ymax = hi),
                fill = col_t1d, alpha = 0.2) +
    geom_line(data = pred_df, aes(x = x, y = est), linewidth = 0.8, color = "black") +
    geom_point(data = donor_means, aes(x = Disease.Duration, y = md),
               color = col_t1d, alpha = 0.6, size = 2) +
    labs(x = "Disease duration (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Disease duration spline (T1D)") +
    theme_pub
}

# =============================================================================
# 21. DD overall (linear slope)
# =============================================================================

plot_t1d_dd_overall <- function(data, model) {
  dr <- range(data$Disease.Duration, na.rm = TRUE)
  md <- mean(data$Disease.Duration, na.rm = TRUE)
  fg <- data.frame(Disease.Duration_c = seq(dr[1] - md, dr[2] - md, length.out = 100)) %>%
    mutate(Disease.Duration = Disease.Duration_c + md)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = c("Region", "Sex"))

  donor_means <- data %>% group_by(Donor, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trend <- emtrends(model, ~ 1, var = "Disease.Duration_c", data = data)
  p_val <- as.data.frame(test(trend))$p.value[1]
  stars <- p_to_stars(p_val)

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Disease.Duration, ymin = lower, ymax = upper),
                fill = col_t1d, alpha = 0.2) +
    geom_line(data = preds, aes(x = Disease.Duration, y = estimate),
              linewidth = 0.9, color = "black") +
    geom_point(data = donor_means, aes(x = Disease.Duration, y = md),
               color = col_t1d, alpha = 0.6, size = 2) +
    labs(x = "Disease duration (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Disease duration (T1D)") +
    theme_pub

  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = dr[1] + 1, y = max(donor_means$md) * 0.95,
                           label = stars, size = 8, hjust = 0)
  plt
}

# =============================================================================
# 22. DD × Region (T1D)
# =============================================================================

plot_t1d_dd_by_region <- function(data, model) {
  dr <- range(data$Disease.Duration, na.rm = TRUE)
  md <- mean(data$Disease.Duration, na.rm = TRUE)
  fg <- expand.grid(
    Disease.Duration_c = seq(dr[1] - md, dr[2] - md, length.out = 100),
    Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
    stringsAsFactors = FALSE
  ) %>% mutate(Disease.Duration = Disease.Duration_c + md)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = "Sex")
  dm <- data %>% group_by(Donor, Region, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trends <- emtrends(model, pairwise ~ Region, var = "Disease.Duration_c", data = data)
  pairs_df <- as.data.frame(trends$contrasts)
  sig_pairs <- pairs_df %>%
    filter(p.value < 0.05) %>%
    mutate(t = sprintf("%s %s", contrast, p_to_stars(p.value)))
  sig_text <- paste(sig_pairs$t, collapse = "\n")

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Disease.Duration, ymin = lower, ymax = upper,
                                   fill = Region), alpha = 0.2) +
    geom_line(data = preds, aes(x = Disease.Duration, y = estimate, color = Region),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = Disease.Duration, y = md, color = Region),
               alpha = 0.4, size = 1.5) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    labs(x = "Disease duration (years)", y = "Islet density (islets/mm\u00B2)",
         title = "DD \u00D7 Region (T1D)") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(sig_text) > 0)
    plt <- plt + annotate("text", x = dr[1] + 1, y = max(dm$md) * 0.95,
                           label = sig_text, size = 4, hjust = 0, lineheight = 0.9)
  plt
}

# =============================================================================
# 23. DD × Sex (T1D)
# =============================================================================

plot_t1d_dd_by_sex <- function(data, model) {
  dr <- range(data$Disease.Duration, na.rm = TRUE)
  md <- mean(data$Disease.Duration, na.rm = TRUE)
  fg <- expand.grid(
    Disease.Duration_c = seq(dr[1] - md, dr[2] - md, length.out = 100),
    Sex = factor(c("Female", "Male"), levels = levels(data$Sex)),
    stringsAsFactors = FALSE
  ) %>% mutate(Disease.Duration = Disease.Duration_c + md)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = "Region")
  dm <- data %>% group_by(Donor, Sex, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trends <- emtrends(model, pairwise ~ Sex, var = "Disease.Duration_c", data = data)
  p_val <- as.data.frame(trends$contrasts)$p.value[1]
  stars <- p_to_stars(p_val)

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Disease.Duration, ymin = lower, ymax = upper,
                                   fill = Sex), alpha = 0.2) +
    geom_line(data = preds, aes(x = Disease.Duration, y = estimate, color = Sex),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = Disease.Duration, y = md, color = Sex),
               alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors_sex) + scale_fill_manual(values = colors_sex) +
    labs(x = "Disease duration (years)", y = "Islet density (islets/mm\u00B2)",
         title = "DD \u00D7 Sex (T1D)") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = dr[1] + 1, y = max(dm$md) * 0.95,
                           label = stars, size = 8, hjust = 0)
  plt
}

# =============================================================================
# 24. AO overall (T1D)
# =============================================================================

plot_t1d_ao_overall <- function(data, model) {
  ar <- range(data$age_at_onset, na.rm = TRUE)
  mao <- mean(data$age_at_onset, na.rm = TRUE)
  fg <- data.frame(age_at_onset_c = seq(ar[1] - mao, ar[2] - mao, length.out = 100)) %>%
    mutate(age_at_onset = age_at_onset_c + mao)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = c("Region", "Sex"))

  dm <- data %>% group_by(Donor, age_at_onset, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop") %>%
    mutate(duration_bin = cut(Disease.Duration, breaks = c(0, 5, 15, Inf),
                              labels = c("0\u20135 yr", "5\u201315 yr", ">15 yr")))

  trend <- emtrends(model, ~ 1, var = "age_at_onset_c", data = data)
  p_val <- as.data.frame(test(trend))$p.value[1]
  stars <- p_to_stars(p_val)

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = age_at_onset, ymin = lower, ymax = upper),
                fill = "grey60", alpha = 0.3) +
    geom_line(data = preds, aes(x = age_at_onset, y = estimate),
              linewidth = 0.9, color = "black") +
    geom_point(data = dm, aes(x = age_at_onset, y = md, color = duration_bin),
               alpha = 0.7, size = 2) +
    scale_color_manual(values = colors_duration, name = "Disease\nduration") +
    labs(x = "Age at onset (years)", y = "Islet density (islets/mm\u00B2)",
         title = "Age at onset (T1D)") +
    theme_pub +
    theme(legend.position = c(0.85, 0.80),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = mean(ar), y = max(dm$md) * 0.95,
                           label = stars, size = 8)
  plt
}

# =============================================================================
# 25. AO × Region (T1D)
# =============================================================================

plot_t1d_ao_by_region <- function(data, model) {
  ar <- range(data$age_at_onset, na.rm = TRUE)
  mao <- mean(data$age_at_onset, na.rm = TRUE)
  fg <- expand.grid(
    age_at_onset_c = seq(ar[1] - mao, ar[2] - mao, length.out = 100),
    Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
    stringsAsFactors = FALSE
  ) %>% mutate(age_at_onset = age_at_onset_c + mao)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = "Sex")
  dm <- data %>% group_by(Donor, Region, age_at_onset) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trends <- emtrends(model, pairwise ~ Region, var = "age_at_onset_c", data = data)
  pairs_df <- as.data.frame(trends$contrasts)
  sig_pairs <- pairs_df %>%
    filter(p.value < 0.05) %>%
    mutate(t = sprintf("%s %s", contrast, p_to_stars(p.value)))
  sig_text <- paste(sig_pairs$t, collapse = "\n")

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = age_at_onset, ymin = lower, ymax = upper,
                                   fill = Region), alpha = 0.2) +
    geom_line(data = preds, aes(x = age_at_onset, y = estimate, color = Region),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = age_at_onset, y = md, color = Region),
               alpha = 0.4, size = 1.5) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    labs(x = "Age at onset (years)", y = "Islet density (islets/mm\u00B2)",
         title = "AO \u00D7 Region (T1D)") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(sig_text) > 0)
    plt <- plt + annotate("text", x = ar[1] + 1, y = max(dm$md) * 0.95,
                           label = sig_text, size = 4, hjust = 0, lineheight = 0.9)
  plt
}

# =============================================================================
# 26. AO × Sex (T1D)
# =============================================================================

plot_t1d_ao_by_sex <- function(data, model) {
  ar <- range(data$age_at_onset, na.rm = TRUE)
  mao <- mean(data$age_at_onset, na.rm = TRUE)
  fg <- expand.grid(
    age_at_onset_c = seq(ar[1] - mao, ar[2] - mao, length.out = 100),
    Sex = factor(c("Female", "Male"), levels = levels(data$Sex)),
    stringsAsFactors = FALSE
  ) %>% mutate(age_at_onset = age_at_onset_c + mao)
  preds <- get_marginalized_preds(model, data, fg, avg_vars = "Region")
  dm <- data %>% group_by(Donor, Sex, age_at_onset) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  trends <- emtrends(model, pairwise ~ Sex, var = "age_at_onset_c", data = data)
  p_val <- as.data.frame(trends$contrasts)$p.value[1]
  stars <- p_to_stars(p_val)

  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = age_at_onset, ymin = lower, ymax = upper,
                                   fill = Sex), alpha = 0.2) +
    geom_line(data = preds, aes(x = age_at_onset, y = estimate, color = Sex),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = age_at_onset, y = md, color = Sex),
               alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors_sex) + scale_fill_manual(values = colors_sex) +
    labs(x = "Age at onset (years)", y = "Islet density (islets/mm\u00B2)",
         title = "AO \u00D7 Sex (T1D)") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))

  if (nchar(stars) > 0)
    plt <- plt + annotate("text", x = ar[1] + 1, y = max(dm$md) * 0.95,
                           label = stars, size = 8, hjust = 0)
  plt
}

# =============================================================================
# 27. DD slopes at low/high AO (T1D)
# =============================================================================

plot_t1d_dd_at_ao_levels <- function(data, model) {
  dr <- range(data$Disease.Duration, na.rm = TRUE)
  md <- mean(data$Disease.Duration, na.rm = TRUE)
  ao_sd <- sd(data$age_at_onset_c, na.rm = TRUE)

  preds_list <- list()
  for (ao_lev in c(-1, 1)) {
    fg <- data.frame(
      Disease.Duration_c = seq(dr[1] - md, dr[2] - md, length.out = 100),
      age_at_onset_c = ao_lev * ao_sd
    ) %>% mutate(Disease.Duration = Disease.Duration_c + md,
                  ao_label = ifelse(ao_lev == -1, "AO = mean \u2212 1 SD", "AO = mean + 1 SD"))
    preds_list[[length(preds_list) + 1]] <-
      get_marginalized_preds(model, data, fg, avg_vars = c("Region", "Sex"))
  }
  preds <- bind_rows(preds_list)

  dm <- data %>% group_by(Donor, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  ggplot() +
    geom_ribbon(data = preds, aes(x = Disease.Duration, ymin = lower, ymax = upper,
                                   fill = ao_label), alpha = 0.2) +
    geom_line(data = preds, aes(x = Disease.Duration, y = estimate, color = ao_label),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = Disease.Duration, y = md),
               color = col_t1d, alpha = 0.4, size = 1.5) +
    labs(x = "Disease duration (years)", y = "Islet density (islets/mm\u00B2)",
         title = "DD slopes at low/high age at onset (T1D)",
         color = NULL, fill = NULL) +
    theme_pub +
    theme(legend.position = c(0.75, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# 28. AO slopes at low/high DD (T1D)
# =============================================================================

plot_t1d_ao_at_dd_levels <- function(data, model) {
  ar <- range(data$age_at_onset, na.rm = TRUE)
  mao <- mean(data$age_at_onset, na.rm = TRUE)
  dd_sd <- sd(data$Disease.Duration_c, na.rm = TRUE)

  preds_list <- list()
  for (dd_lev in c(-1, 1)) {
    fg <- data.frame(
      age_at_onset_c = seq(ar[1] - mao, ar[2] - mao, length.out = 100),
      Disease.Duration_c = dd_lev * dd_sd
    ) %>% mutate(age_at_onset = age_at_onset_c + mao,
                  dd_label = ifelse(dd_lev == -1, "DD = mean \u2212 1 SD", "DD = mean + 1 SD"))
    preds_list[[length(preds_list) + 1]] <-
      get_marginalized_preds(model, data, fg, avg_vars = c("Region", "Sex"))
  }
  preds <- bind_rows(preds_list)

  dm <- data %>% group_by(Donor, age_at_onset, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop") %>%
    mutate(duration_bin = cut(Disease.Duration, breaks = c(0, 5, 15, Inf),
                              labels = c("0\u20135 yr", "5\u201315 yr", ">15 yr")))

  ggplot() +
    geom_ribbon(data = preds, aes(x = age_at_onset, ymin = lower, ymax = upper,
                                   fill = dd_label), alpha = 0.2) +
    geom_line(data = preds, aes(x = age_at_onset, y = estimate, color = dd_label),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = age_at_onset, y = md),
               color = col_t1d, alpha = 0.4, size = 1.5) +
    labs(x = "Age at onset (years)", y = "Islet density (islets/mm\u00B2)",
         title = "AO slopes at low/high disease duration (T1D)",
         color = NULL, fill = NULL) +
    theme_pub +
    theme(legend.position = c(0.75, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}


# =============================================================================
# 29. DD × Region × AO — DD by Region, faceted at AO ± 1 SD
# =============================================================================

plot_t1d_3way_dd_region_ao <- function(data, model) {
  dr <- range(data$Disease.Duration, na.rm = TRUE)
  md <- mean(data$Disease.Duration, na.rm = TRUE)
  ao_sd <- sd(data$age_at_onset_c, na.rm = TRUE)

  preds_list <- list()
  for (ao_lev in c(-1, 1)) {
    fg <- expand.grid(
      Disease.Duration_c = seq(dr[1] - md, dr[2] - md, length.out = 100),
      Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
      stringsAsFactors = FALSE
    ) %>% mutate(
      Disease.Duration = Disease.Duration_c + md,
      age_at_onset_c = ao_lev * ao_sd,
      ao_label = factor(ifelse(ao_lev == -1, "AO = mean \u2212 1 SD", "AO = mean + 1 SD"),
                         levels = c("AO = mean \u2212 1 SD", "AO = mean + 1 SD"))
    )
    preds_list[[length(preds_list) + 1]] <-
      get_marginalized_preds(model, data, fg, avg_vars = "Sex")
  }
  preds <- bind_rows(preds_list)
  dm <- data %>% group_by(Donor, Region, Disease.Duration) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  y_max <- max(c(dm$md, preds$upper), na.rm = TRUE) * 1.1

  ggplot() +
    geom_ribbon(data = preds, aes(x = Disease.Duration, ymin = lower, ymax = upper,
                                   fill = Region), alpha = 0.2) +
    geom_line(data = preds, aes(x = Disease.Duration, y = estimate, color = Region),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = Disease.Duration, y = md, color = Region),
               alpha = 0.3, size = 1.3) +
    facet_wrap(~ ao_label) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(x = "Disease duration (years)", y = "Islet density (islets/mm\u00B2)",
         title = "DD \u00D7 Region \u00D7 AO (T1D)",
         subtitle = "DD slopes by Region at low/high age at onset") +
    theme_pub + theme(legend.position = "bottom")
}

# =============================================================================
# 30. DD × Region × AO — AO by Region, faceted at DD ± 1 SD
# =============================================================================

plot_t1d_3way_ao_region_dd <- function(data, model) {
  ar <- range(data$age_at_onset, na.rm = TRUE)
  mao <- mean(data$age_at_onset, na.rm = TRUE)
  dd_sd <- sd(data$Disease.Duration_c, na.rm = TRUE)

  preds_list <- list()
  for (dd_lev in c(-1, 1)) {
    fg <- expand.grid(
      age_at_onset_c = seq(ar[1] - mao, ar[2] - mao, length.out = 100),
      Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
      stringsAsFactors = FALSE
    ) %>% mutate(
      age_at_onset = age_at_onset_c + mao,
      Disease.Duration_c = dd_lev * dd_sd,
      dd_label = factor(ifelse(dd_lev == -1, "DD = mean \u2212 1 SD", "DD = mean + 1 SD"),
                         levels = c("DD = mean \u2212 1 SD", "DD = mean + 1 SD"))
    )
    preds_list[[length(preds_list) + 1]] <-
      get_marginalized_preds(model, data, fg, avg_vars = "Sex")
  }
  preds <- bind_rows(preds_list)
  dm <- data %>% group_by(Donor, Region, age_at_onset) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")

  y_max <- max(c(dm$md, preds$upper), na.rm = TRUE) * 1.1

  ggplot() +
    geom_ribbon(data = preds, aes(x = age_at_onset, ymin = lower, ymax = upper,
                                   fill = Region), alpha = 0.2) +
    geom_line(data = preds, aes(x = age_at_onset, y = estimate, color = Region),
              linewidth = 0.9) +
    geom_point(data = dm, aes(x = age_at_onset, y = md, color = Region),
               alpha = 0.3, size = 1.3) +
    facet_wrap(~ dd_label) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    scale_y_continuous(limits = c(0, y_max)) +
    labs(x = "Age at onset (years)", y = "Islet density (islets/mm\u00B2)",
         title = "AO \u00D7 Region \u00D7 DD (T1D)",
         subtitle = "AO slopes by Region at low/high disease duration") +
    theme_pub + theme(legend.position = "bottom")
}


# #############################################################################
#
#  MAIN EXECUTION
#
# #############################################################################

generate_all_density_plots <- function(output_dir = "Density/Graphs/All_Comparisons") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL DENSITY POST-HOC COMPARISON PLOTS\n")
  cat("=============================================================================\n\n")

  # --- Load data ---
  cat("Loading data and models...\n")
  density_data <- readRDS("Data/density_data.rds") %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
      Region    = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor)
    )
  density_t1d_data <- readRDS("Data/density_data_t1d.rds") %>%
    mutate(
      Region = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex    = factor(Sex, levels = c("Female", "Male")),
      Donor  = factor(Donor)
    ) %>% droplevels()

  density_diag   <- readRDS("Density/Models/density_diag.rds")
  density_t1d    <- readRDS("Density/Models/density_t1d.rds")
  density_spline <- readRDS("Density/Models/density_spline.rds")

  cat("  density_data:", nrow(density_data), "rows,", n_distinct(density_data$Donor), "donors\n")
  cat("  density_t1d: ", nrow(density_t1d_data), "rows,", n_distinct(density_t1d_data$Donor), "donors\n\n")

  dir_diag <- file.path(output_dir, "Diagnosis")
  dir_t1d  <- file.path(output_dir, "T1D")

  # =========================================================================
  # DIAGNOSIS MODEL
  # =========================================================================
  cat("--- DIAGNOSIS MODEL (15 plots) ---\n")

  save_plot(plot_diag_main_diagnosis(density_data, density_diag),
            "01_main_diagnosis", dir_diag)
  save_plot(plot_diag_main_region(density_data, density_diag),
            "02_main_region", dir_diag)
  save_plot(plot_diag_main_sex(density_data, density_diag),
            "03_main_sex", dir_diag)

  save_plot(plot_diag_diagnosis_by_region(density_data, density_diag),
            "04_diagnosis_by_region", dir_diag, width = 8)
  save_plot(plot_diag_region_by_diagnosis(density_data, density_diag),
            "05_region_by_diagnosis", dir_diag, width = 9, height = 5)
  save_plot(plot_diag_diagnosis_by_sex(density_data, density_diag),
            "06_diagnosis_by_sex", dir_diag)
  save_plot(plot_diag_sex_by_diagnosis(density_data, density_diag),
            "07_sex_by_diagnosis", dir_diag)
  save_plot(plot_diag_region_by_sex(density_data, density_diag),
            "08_region_by_sex", dir_diag, width = 9, height = 5)
  save_plot(plot_diag_sex_by_region(density_data, density_diag),
            "09_sex_by_region", dir_diag, width = 8)

  save_plot(plot_diag_age_overall(density_data, density_diag),
            "10_age_overall", dir_diag)
  save_plot(plot_diag_age_by_diagnosis(density_data, density_diag),
            "11_age_by_diagnosis", dir_diag)
  save_plot(plot_diag_age_by_region(density_data, density_diag),
            "12_age_by_region", dir_diag)
  save_plot(plot_diag_age_by_sex(density_data, density_diag),
            "13_age_by_sex", dir_diag)

  save_plot(plot_diag_3way_facet_region(density_data, density_diag),
            "14_3way_diag_region_age_facet_region", dir_diag, width = 12, height = 5)
  save_plot(plot_diag_3way_facet_diagnosis(density_data, density_diag),
            "15_3way_diag_region_age_facet_diagnosis", dir_diag, width = 10, height = 5)

  # =========================================================================
  # T1D MODEL
  # =========================================================================
  cat("\n--- T1D MODEL (15 plots) ---\n")

  save_plot(plot_t1d_main_region(density_t1d_data, density_t1d),
            "16_t1d_main_region", dir_t1d)
  save_plot(plot_t1d_main_sex(density_t1d_data, density_t1d),
            "17_t1d_main_sex", dir_t1d)
  save_plot(plot_t1d_region_by_sex(density_t1d_data, density_t1d),
            "18_t1d_region_by_sex", dir_t1d, width = 9, height = 5)
  save_plot(plot_t1d_sex_by_region(density_t1d_data, density_t1d),
            "19_t1d_sex_by_region", dir_t1d, width = 8)

  save_plot(plot_t1d_dd_spline(density_t1d_data, density_spline),
            "20_t1d_dd_spline", dir_t1d)
  save_plot(plot_t1d_dd_overall(density_t1d_data, density_t1d),
            "21_t1d_dd_overall", dir_t1d)
  save_plot(plot_t1d_dd_by_region(density_t1d_data, density_t1d),
            "22_t1d_dd_by_region", dir_t1d)
  save_plot(plot_t1d_dd_by_sex(density_t1d_data, density_t1d),
            "23_t1d_dd_by_sex", dir_t1d)

  save_plot(plot_t1d_ao_overall(density_t1d_data, density_t1d),
            "24_t1d_ao_overall", dir_t1d)
  save_plot(plot_t1d_ao_by_region(density_t1d_data, density_t1d),
            "25_t1d_ao_by_region", dir_t1d)
  save_plot(plot_t1d_ao_by_sex(density_t1d_data, density_t1d),
            "26_t1d_ao_by_sex", dir_t1d)

  save_plot(plot_t1d_dd_at_ao_levels(density_t1d_data, density_t1d),
            "27_t1d_dd_at_ao_levels", dir_t1d)
  save_plot(plot_t1d_ao_at_dd_levels(density_t1d_data, density_t1d),
            "28_t1d_ao_at_dd_levels", dir_t1d)

  save_plot(plot_t1d_3way_dd_region_ao(density_t1d_data, density_t1d),
            "29_t1d_3way_dd_region_ao", dir_t1d, width = 10, height = 5)
  save_plot(plot_t1d_3way_ao_region_dd(density_t1d_data, density_t1d),
            "30_t1d_3way_ao_region_dd", dir_t1d, width = 10, height = 5)

  # =========================================================================
  # CLEANUP
  # =========================================================================
  rm(density_diag, density_t1d, density_spline); gc()

  cat("\n=============================================================================\n")
  cat("COMPLETE! 30 plots saved to", output_dir, "\n")
  cat("  Diagnosis/ (15 plots)\n")
  cat("  T1D/       (15 plots)\n")
  cat("=============================================================================\n")
}


# --- Run if interactive ---
if (interactive()) {
  cat("\nTo generate all density comparison plots, run:\n")
  cat("  generate_all_density_plots()\n\n")
  cat("Or call individual plot functions, e.g.:\n")
  cat("  plot_diag_main_diagnosis(data, model)\n")
}
