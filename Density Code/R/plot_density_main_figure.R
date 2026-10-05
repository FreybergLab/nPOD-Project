# =============================================================================
# plot_density_main_figure.R
# PI-approved style: clean stars, no annotation boxes, legends inside panels
# Only changes from original Density.Rmd: master palette + marginalized CI on B
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(patchwork)
library(splines)

source("colors_master.R")

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

p_to_stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

# ── Panel A: ND vs T1D ───────────────────────────────────────────────────

plot_panel_a <- function(Density, model_dx) {
  donor_means <- Density %>%
    group_by(Donor, Diagnosis) %>%
    summarise(mean_density = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  emm_dx <- as.data.frame(emmeans(model_dx, ~ Diagnosis, type = "response"))
  pairs_dx <- as.data.frame(pairs(emmeans(model_dx, ~ Diagnosis, type = "response")))
  stars <- p_to_stars(pairs_dx$p.value[1])
  
  ggplot() +
    geom_violin(data = donor_means, aes(x = Diagnosis, y = mean_density, fill = Diagnosis),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = donor_means, aes(x = Diagnosis, y = mean_density, color = Diagnosis),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    annotate("text", x = 1.5, y = max(donor_means$mean_density) * 0.95,
             label = stars, size = 8) +
    labs(x = NULL, y = "Islet density\n(islets/mm\u00B2)", tag = "A",
         title = "ND vs T1D") +
    theme_pub
}

# ── Panel B: Disease duration (T1D, marginalized spline + CI) ────────────

plot_panel_b <- function(Density_T1D, model_spline) {
  # Ensure factors
  Density_T1D <- Density_T1D %>%
    mutate(Region = factor(Region), Sex = factor(Sex))
  
  donor_means_t1d <- Density_T1D %>%
    group_by(Donor, Disease.Duration) %>%
    summarise(mean_density = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  # Build prediction grid on centered scale
  mean_dur <- mean(Density_T1D$Disease.Duration, na.rm = TRUE)
  dur_raw_seq <- seq(0.1, max(Density_T1D$Disease.Duration, na.rm = TRUE), length.out = 100)
  dur_c_seq <- dur_raw_seq - mean_dur
  
  # Marginalize over Region and Sex
  reg_lvls <- levels(Density_T1D$Region)
  sex_lvls <- levels(Density_T1D$Sex)
  avg_levels <- expand.grid(Region = reg_lvls, Sex = sex_lvls, stringsAsFactors = FALSE)
  all_preds <- all_se <- list()
  for (i in seq_len(nrow(avg_levels))) {
    nd <- data.frame(Disease.Duration_c = dur_c_seq, Age_c = 0,
                     Region = factor(avg_levels$Region[i], levels = reg_lvls),
                     Sex = factor(avg_levels$Sex[i], levels = sex_lvls),
                     Donor = NA)
    p <- predict(model_spline, newdata = nd, type = "link", re.form = NA,
                 se.fit = TRUE, allow.new.levels = TRUE)
    all_preds[[i]] <- p$fit; all_se[[i]] <- p$se.fit
  }
  avg_p <- colMeans(do.call(rbind, all_preds))
  avg_s <- sqrt(colMeans(do.call(rbind, all_se)^2))
  pred_df <- data.frame(x = dur_raw_seq, est = exp(avg_p),
                        lo = exp(avg_p - 1.96 * avg_s),
                        hi = exp(avg_p + 1.96 * avg_s))
  
  ggplot() +
    geom_ribbon(data = pred_df, aes(x = x, ymin = lo, ymax = hi),
                fill = col_t1d, alpha = 0.2) +
    geom_line(data = pred_df, aes(x = x, y = est),
              linewidth = 0.8, color = "black") +
    geom_point(data = donor_means_t1d,
               aes(x = Disease.Duration, y = mean_density),
               color = col_t1d, alpha = 0.6, size = 2) +
    labs(x = "Disease duration (years)", y = "Islet density\n(islets/mm\u00B2)",
         tag = "B", title = "Density vs disease duration (T1D)") +
    theme_pub
}

# ── Panel C: Age at onset (T1D) ──────────────────────────────────────────

plot_panel_c <- function(Density_T1D, model_t1d) {
  donor_means_t1d <- Density_T1D %>%
    group_by(Donor, age_at_onset, Disease.Duration) %>%
    summarise(mean_density = mean(islet_density, na.rm = TRUE), .groups = "drop") %>%
    mutate(duration_bin = cut(Disease.Duration, breaks = c(0, 5, 15, Inf),
                              labels = c("0\u20135 yr", "5\u201315 yr", ">15 yr")))
  
  # Test marginalized AO slope (averaged over Region, Sex, at DD = 0)
  ao_trend <- emtrends(model_t1d, ~ 1, var = "age_at_onset_c", data = Density_T1D)
  ao_test  <- as.data.frame(test(ao_trend))
  stars    <- p_to_stars(ao_test$p.value[1])
  
  plt <- ggplot(donor_means_t1d, aes(x = age_at_onset, y = mean_density, color = duration_bin)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "black", linewidth = 0.7, alpha = 0.15) +
    scale_color_manual(values = colors_duration, name = "Disease\nduration") +
    labs(x = "Age at onset (years)", y = "Islet density\n(islets/mm\u00B2)",
         tag = "C", title = "Density vs age at T1D onset") +
    theme_pub +
    theme(legend.position = c(0.85, 0.80),
          legend.background = element_rect(fill = "white", color = NA))
  
  if (nchar(stars) > 0) {
    plt <- plt + annotate("text", x = mean(range(donor_means_t1d$age_at_onset, na.rm = TRUE)),
                          y = max(donor_means_t1d$mean_density, na.rm = TRUE) * 0.95,
                          label = stars, size = 8)
  }
  plt
}

# ── Panel D: Region × Diagnosis ──────────────────────────────────────────

plot_panel_d <- function(Density, model_region) {
  emm_region <- as.data.frame(emmeans(model_region, ~ Diagnosis | Region, type = "response"))
  donor_means_region <- Density %>%
    group_by(Donor, Diagnosis, Region) %>%
    summarise(mean_density = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  pairs_df <- as.data.frame(pairs(emmeans(model_region, ~ Diagnosis | Region, type = "response")))
  stars <- sapply(pairs_df$p.value, p_to_stars)
  y_max <- max(donor_means_region$mean_density, na.rm = TRUE)
  
  ggplot() +
    geom_violin(data = donor_means_region,
                aes(x = Region, y = mean_density, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = donor_means_region,
                aes(x = Region, y = mean_density, color = Diagnosis, group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) +
    scale_color_manual(values = colors_diag) +
    annotate("text", x = 1, y = y_max * 0.98, label = stars[1], size = 8) +
    annotate("text", x = 2, y = y_max * 0.93, label = stars[2], size = 8) +
    annotate("text", x = 3, y = y_max * 0.88, label = stars[3], size = 8) +
    labs(x = "Pancreatic region", y = "Islet density\n(islets/mm\u00B2)",
         tag = "D", title = "Regional heterogeneity of density loss",
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.95, 0.60),
          legend.background = element_rect(fill = "white", color = NA))
}

# ── Panel E: Sex × Diagnosis ────────────────────────────────────────────

plot_panel_e <- function(Density, model_sex) {
  emm_sex <- as.data.frame(emmeans(model_sex, ~ Sex | Diagnosis, type = "response"))
  donor_means_sex <- Density %>%
    group_by(Donor, Diagnosis, Sex) %>%
    summarise(mean_density = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  pairs_df <- as.data.frame(pairs(emmeans(model_sex, ~ Sex | Diagnosis, type = "response")))
  t1d_stars <- p_to_stars(pairs_df$p.value[pairs_df$Diagnosis == "T1D"])
  
  ggplot() +
    geom_violin(data = donor_means_sex,
                aes(x = Diagnosis, y = mean_density, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = donor_means_sex,
                aes(x = Diagnosis, y = mean_density, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.5, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) +
    scale_color_manual(values = colors_sex) +
    annotate("text", x = 2, y = max(donor_means_sex$mean_density) * 0.90,
             label = t1d_stars, size = 8) +
    labs(x = NULL, y = "Islet density\n(islets/mm\u00B2)",
         tag = "E", title = "Sex \u00D7 diagnosis interaction",
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.95, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# GENERATE
# =============================================================================

generate_main_figure <- function(Density, Density_T1D,
                                  model_full, model_spline, model_t1d,
                                  save = TRUE, output_dir = "Density/Graphs") {
  cat("Generating main figure...\n")
  p_a <- plot_panel_a(Density, model_full);          cat("  A\n")
  p_b <- plot_panel_b(Density_T1D, model_spline);    cat("  B\n")
  p_c <- plot_panel_c(Density_T1D, model_t1d);     cat("  C\n")
  p_d <- plot_panel_d(Density, model_full);           cat("  D\n")
  p_e <- plot_panel_e(Density, model_full);           cat("  E\n")
  
  fig <- (p_a | p_b | p_c) / (p_d | p_e) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
  
  if (save) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(output_dir, "density_main_figure.pdf"), fig,
           width = 15, height = 10, dpi = 300)
    ggsave(file.path(output_dir, "density_main_figure.png"), fig,
           width = 15, height = 10, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "density_main_figure.tiff"), fig,
           width = 15, height = 10, dpi = 300, bg = "white")
    cat("Saved to", output_dir, "\n")
  }
  
  list(fig = fig, panels = list(a = p_a, b = p_b, c = p_c, d = p_d, e = p_e))
}
