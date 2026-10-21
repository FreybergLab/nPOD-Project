# =============================================================================
# plot_marker_all_comparisons.R
# Every post-hoc comparison for the six marker rate models.
#
# Style follows plot_density_all_comparisons.R: master palette, donor-level
# violin + jitter, clean stars, save_plot() to pdf/png/tiff at 300 dpi.
#
# Output: Marker/Graphs/All_Comparisons/<Model>/NN_<comparison>.{pdf,png,tiff}
#
# Donor-level points are each donor's own rate, 100 * sum(numerator) /
# sum(denominator) over that donor's islets — not the mean of per-islet ratios,
# which is undefined whenever a donor has an islet with no cells of that type.
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)

source("colors_master.R")
source("Marker/R/marker_common.R")   # registry, data loading, model resolution
source("Marker/R/marker_effects.R")  # verified marginal effects

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

p_to_stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                                    p < 0.05 ~ "*", TRUE ~ "")

# The Region, Sex and interaction panels still derive their stars from the
# emmeans path that produced the wrong Diagnosis marginals (4.65-5.07% for
# every outcome). Those contrasts have not been checked against the data, so
# their stars are suppressed; the donor-level violins and points beneath them
# are computed directly and are correct. The Diagnosis panel uses the verified
# marginal_diagnosis() and keeps its stars.
#
# Set TRUE only after rebuilding rate_means() and rate_contrasts() on
# predict(type = "response") and re-running validate_marker_rates().
SHOW_UNVERIFIED_STARS <- FALSE

save_plot <- function(plt, filename, output_dir, width = 7, height = 5) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("pdf", "png", "tiff")) {
    ggsave(file.path(output_dir, paste0(filename, ".", ext)), plt,
           width = width, height = height, dpi = 300,
           bg = if (ext != "pdf") "white" else NULL)
  }
  cat("  Saved:", filename, "\n")
}

# =============================================================================
# HELPERS
# =============================================================================

has_zi <- function(model) length(fixef(model)$zi) > 0

# Rates per 100 cells, with the offset pinned so emmeans returns a rate.
emm_rate <- function(model, specs, data, per = PER, at = list()) {
  if (is.null(at$Age_c)) at$Age_c <- 0
  out <- tryCatch(
    emmeans(model, specs, data = data, offset = log(per),
            component = "response", at = at),
    error = function(e) NULL)
  if (is.null(out)) {
    out <- emmeans(model, specs, data = data, offset = log(per), at = at)
  }
  out
}

rate_ratio_p <- function(emm, method = "revpairwise", adjust = "none") {
  if (!SHOW_UNVERIFIED_STARS) return(NULL)   # see the note at the top
  ct <- tryCatch(as.data.frame(summary(contrast(
    regrid(emm, transform = "log"), method, adjust = adjust,
    type = "response"))), error = function(e) NULL)
  if (is.null(ct)) return(NULL)
  ct
}

# Each donor's own rate (%), optionally split by grouping variables.
donor_rates <- function(data, def, extra = character(0)) {
  data %>%
    group_by(across(all_of(c("Donor", "Diagnosis", extra)))) %>%
    summarise(num = sum(.data[[def$num]]), den = sum(.data[[def$den]]),
              .groups = "drop") %>%
    filter(den > 0) %>%
    mutate(rate = PER * num / den)
}

y_lab <- function(def) sprintf("%s+ rate\n(%% of %s cells)", def$marker, def$cell)
ttl   <- function(def, what) sprintf("%s in %s cells — %s", def$marker, def$cell, what)

# =============================================================================
# 1. MAIN EFFECT: Diagnosis
# =============================================================================

plot_rate_main_diagnosis <- function(data, model, def) {
  dr <- donor_rates(data, def)
  # Shared with the main figure so the two never disagree.
  eff <- marginal_diagnosis(model)   # verified; see marker_effects.R
  stars <- if (!is.null(eff)) p_to_stars(eff$p) else ""

  ggplot() +
    geom_violin(data = dr, aes(x = Diagnosis, y = rate, fill = Diagnosis),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dr, aes(x = Diagnosis, y = rate, color = Diagnosis),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    { if (nchar(stars) > 0) annotate("text", x = 1.5, y = max(dr$rate) * 0.95,
                                      label = stars, size = 8) } +
    labs(x = NULL, y = y_lab(def), title = ttl(def, "Diagnosis")) +
    theme_pub
}

# =============================================================================
# 2. MAIN EFFECT: Region
# =============================================================================

plot_rate_main_region <- function(data, model, def) {
  dr <- donor_rates(data, def, "Region")
  ct <- rate_ratio_p(emm_rate(model, ~ Region, data), "pairwise", "tukey")
  sig_text <- if (!is.null(ct)) {
    s <- ct %>% mutate(label = paste(contrast, p_to_stars(p.value))) %>%
      filter(p.value < 0.05)
    paste(s$label, collapse = "\n")
  } else ""

  ggplot() +
    geom_violin(data = dr, aes(x = Region, y = rate, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dr, aes(x = Region, y = rate, color = Region),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    { if (nchar(sig_text) > 0) annotate("text", x = 2, y = max(dr$rate) * 0.95,
                                         label = sig_text, size = 4,
                                         lineheight = 0.9) } +
    labs(x = "Pancreatic region", y = y_lab(def), title = ttl(def, "Region")) +
    theme_pub
}

# =============================================================================
# 3. MAIN EFFECT: Sex
# =============================================================================

plot_rate_main_sex <- function(data, model, def) {
  dr <- donor_rates(data, def, "Sex")
  ct <- rate_ratio_p(emm_rate(model, ~ Sex, data))
  stars <- if (!is.null(ct)) p_to_stars(ct$p.value[1]) else ""

  ggplot() +
    geom_violin(data = dr, aes(x = Sex, y = rate, fill = Sex),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dr, aes(x = Sex, y = rate, color = Sex),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_sex, guide = "none") +
    scale_color_manual(values = colors_sex, guide = "none") +
    { if (nchar(stars) > 0) annotate("text", x = 1.5, y = max(dr$rate) * 0.95,
                                      label = stars, size = 8) } +
    labs(x = NULL, y = y_lab(def), title = ttl(def, "Sex")) +
    theme_pub
}

# =============================================================================
# 4. DIAGNOSIS WITHIN REGION
# =============================================================================

plot_rate_diagnosis_by_region <- function(data, model, def) {
  dr <- donor_rates(data, def, "Region")
  ct <- rate_ratio_p(emm_rate(model, ~ Diagnosis | Region, data))
  y_max <- max(dr$rate)
  ann <- if (!is.null(ct)) {
    tibble(Region = factor(ct$Region, levels = levels(dr$Region)),
           lab = p_to_stars(ct$p.value), y = y_max * 0.95)
  } else NULL

  p <- ggplot() +
    geom_violin(data = dr, aes(x = Region, y = rate, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8),
                width = 0.7) +
    geom_jitter(data = dr, aes(x = Region, y = rate, color = Diagnosis,
                               group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12,
                                                dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) +
    scale_color_manual(values = colors_diag) +
    labs(x = "Pancreatic region", y = y_lab(def),
         title = ttl(def, "Diagnosis within Region"),
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.92, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  if (!is.null(ann)) p <- p + geom_text(data = ann, aes(x = Region, y = y,
                                                        label = lab),
                                        size = 8, inherit.aes = FALSE)
  p
}

# =============================================================================
# 5. REGION WITHIN DIAGNOSIS
# =============================================================================

plot_rate_region_by_diagnosis <- function(data, model, def) {
  dr <- donor_rates(data, def, "Region")
  ct <- rate_ratio_p(emm_rate(model, ~ Region | Diagnosis, data),
                     "pairwise", "tukey")
  lab_df <- if (!is.null(ct)) {
    ct %>% mutate(label = paste(contrast, p_to_stars(p.value))) %>%
      filter(p.value < 0.05) %>%
      group_by(Diagnosis) %>%
      summarise(lab = paste(label, collapse = "\n"), .groups = "drop") %>%
      mutate(Diagnosis = factor(Diagnosis, levels = levels(dr$Diagnosis)))
  } else NULL

  p <- ggplot(dr, aes(x = Region, y = rate, fill = Region)) +
    geom_violin(alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(aes(color = Region), width = 0.15, size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    facet_wrap(~ Diagnosis) +
    labs(x = "Pancreatic region", y = y_lab(def),
         title = ttl(def, "Region within Diagnosis")) +
    theme_pub
  if (!is.null(lab_df) && nrow(lab_df) > 0) {
    # x = 2 is the middle level of the discrete Region scale; lab_df carries
    # Diagnosis so the label lands in the right facet.
    p <- p + geom_text(data = lab_df,
                       aes(x = 2, y = max(dr$rate) * 0.95, label = lab),
                       size = 3.5, lineheight = 0.9, inherit.aes = FALSE)
  }
  p
}

# =============================================================================
# 6. DIAGNOSIS WITHIN SEX
# =============================================================================

plot_rate_diagnosis_by_sex <- function(data, model, def) {
  dr <- donor_rates(data, def, "Sex")
  ct <- rate_ratio_p(emm_rate(model, ~ Diagnosis | Sex, data))
  ann <- if (!is.null(ct)) {
    tibble(Sex = factor(ct$Sex, levels = levels(dr$Sex)),
           lab = p_to_stars(ct$p.value), y = max(dr$rate) * 0.95)
  } else NULL

  p <- ggplot() +
    geom_violin(data = dr, aes(x = Sex, y = rate, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8),
                width = 0.7) +
    geom_jitter(data = dr, aes(x = Sex, y = rate, color = Diagnosis,
                               group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12,
                                                dodge.width = 0.8),
                size = 1.4, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) +
    scale_color_manual(values = colors_diag) +
    labs(x = NULL, y = y_lab(def),
         title = ttl(def, "Diagnosis within Sex"),
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.92, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  if (!is.null(ann)) p <- p + geom_text(data = ann, aes(x = Sex, y = y,
                                                        label = lab),
                                        size = 8, inherit.aes = FALSE)
  p
}

# =============================================================================
# 7. SEX WITHIN DIAGNOSIS
# =============================================================================

plot_rate_sex_by_diagnosis <- function(data, model, def) {
  dr <- donor_rates(data, def, "Sex")
  ct <- rate_ratio_p(emm_rate(model, ~ Sex | Diagnosis, data))
  ann <- if (!is.null(ct)) {
    tibble(Diagnosis = factor(ct$Diagnosis, levels = levels(dr$Diagnosis)),
           lab = p_to_stars(ct$p.value), y = max(dr$rate) * 0.95)
  } else NULL

  p <- ggplot() +
    geom_violin(data = dr, aes(x = Diagnosis, y = rate, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8),
                width = 0.7) +
    geom_jitter(data = dr, aes(x = Diagnosis, y = rate, color = Sex,
                               group = Sex),
                position = position_jitterdodge(jitter.width = 0.12,
                                                dodge.width = 0.8),
                size = 1.4, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) +
    scale_color_manual(values = colors_sex) +
    labs(x = NULL, y = y_lab(def), title = ttl(def, "Sex within Diagnosis"),
         fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.92, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  if (!is.null(ann)) p <- p + geom_text(data = ann, aes(x = Diagnosis, y = y,
                                                        label = lab),
                                        size = 8, inherit.aes = FALSE)
  p
}

# =============================================================================
# 8. REGION WITHIN SEX
# =============================================================================

plot_rate_region_by_sex <- function(data, model, def) {
  dr <- donor_rates(data, def, c("Region", "Sex"))
  ggplot(dr, aes(x = Region, y = rate, fill = Region)) +
    geom_violin(alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(aes(color = Region), width = 0.15, size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    facet_wrap(~ Sex) +
    labs(x = "Pancreatic region", y = y_lab(def),
         title = ttl(def, "Region within Sex")) +
    theme_pub
}

# =============================================================================
# 9. SEX WITHIN REGION
# =============================================================================

plot_rate_sex_by_region <- function(data, model, def) {
  dr <- donor_rates(data, def, c("Region", "Sex"))
  ggplot(dr, aes(x = Sex, y = rate, fill = Sex)) +
    geom_violin(alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(aes(color = Sex), width = 0.15, size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_sex, guide = "none") +
    scale_color_manual(values = colors_sex, guide = "none") +
    facet_wrap(~ Region) +
    labs(x = NULL, y = y_lab(def), title = ttl(def, "Sex within Region")) +
    theme_pub
}

# =============================================================================
# 10. AGE, OVERALL
# =============================================================================

plot_rate_age_overall <- function(data, model, def) {
  dr <- data %>%
    group_by(Donor, Diagnosis, Age) %>%
    summarise(num = sum(.data[[def$num]]), den = sum(.data[[def$den]]),
              .groups = "drop") %>%
    filter(den > 0) %>% mutate(rate = PER * num / den)

  tr <- tryCatch(as.data.frame(summary(
    emtrends(model, ~ 1, var = "Age_c", data = data), infer = c(TRUE, TRUE))),
    error = function(e) NULL)
  sub <- if (!is.null(tr)) {
    sl <- grep("trend$", names(tr), value = TRUE)[1]
    sprintf("%.2f%% change in rate per year (p = %.4f)",
            100 * (exp(tr[[sl]][1]) - 1), tr$p.value[1])
  } else NULL

  ggplot(dr, aes(x = Age, y = rate)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7,
                alpha = 0.15) +
    geom_point(aes(color = Diagnosis), size = 2, alpha = 0.7) +
    scale_color_manual(values = colors_diag) +
    labs(x = "Donor age (years)", y = y_lab(def),
         title = ttl(def, "Age"), subtitle = sub, color = "Diagnosis") +
    theme_pub +
    theme(legend.position = c(0.9, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# 11. AGE WITHIN DIAGNOSIS
# =============================================================================

plot_rate_age_by_diagnosis <- function(data, model, def) {
  dr <- data %>%
    group_by(Donor, Diagnosis, Age) %>%
    summarise(num = sum(.data[[def$num]]), den = sum(.data[[def$den]]),
              .groups = "drop") %>%
    filter(den > 0) %>% mutate(rate = PER * num / den)

  ggplot(dr, aes(x = Age, y = rate, color = Diagnosis, fill = Diagnosis)) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = colors_diag) +
    scale_fill_manual(values = colors_diag) +
    facet_wrap(~ Diagnosis) +
    labs(x = "Donor age (years)", y = y_lab(def),
         title = ttl(def, "Age within Diagnosis")) +
    theme_pub + theme(legend.position = "none")
}

# =============================================================================
# 12. AGE WITHIN REGION
# =============================================================================

plot_rate_age_by_region <- function(data, model, def) {
  dr <- data %>%
    group_by(Donor, Diagnosis, Region, Age) %>%
    summarise(num = sum(.data[[def$num]]), den = sum(.data[[def$den]]),
              .groups = "drop") %>%
    filter(den > 0) %>% mutate(rate = PER * num / den)

  ggplot(dr, aes(x = Age, y = rate)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6,
                alpha = 0.15) +
    geom_point(aes(color = Diagnosis), size = 1.6, alpha = 0.7) +
    scale_color_manual(values = colors_diag) +
    facet_wrap(~ Region) +
    labs(x = "Donor age (years)", y = y_lab(def),
         title = ttl(def, "Age within Region"), color = "Diagnosis") +
    theme_pub + theme(legend.position = "bottom")
}

# =============================================================================
# 13. ZERO-INFLATION
# =============================================================================

plot_rate_zero_inflation <- function(data, model, def) {
  if (!has_zi(model)) return(NULL)
  zb <- fixef(model)$zi
  p_nd  <- plogis(zb[["(Intercept)"]])
  p_t1d <- if ("DiagnosisT1D" %in% names(zb))
    plogis(zb[["(Intercept)"]] + zb[["DiagnosisT1D"]]) else p_nd

  zi <- tibble(Diagnosis = factor(c("ND", "T1D"), levels = c("ND", "T1D")),
               p = c(p_nd, p_t1d))
  obs <- data %>%
    group_by(Donor, Diagnosis) %>%
    summarise(p_zero = mean(.data[[def$num]] == 0), .groups = "drop")

  ggplot() +
    geom_col(data = zi, aes(x = Diagnosis, y = 100 * p, fill = Diagnosis),
             alpha = 0.6, width = 0.6) +
    geom_jitter(data = obs, aes(x = Diagnosis, y = 100 * p_zero),
                width = 0.12, size = 1.4, alpha = 0.5, color = "gray25") +
    scale_fill_manual(values = colors_diag, guide = "none") +
    labs(x = NULL, y = "Islets (%)",
         title = ttl(def, "Zero-inflation"),
         subtitle = paste("Bars: fitted P(structural zero).",
                          "Points: each donor's observed % of islets with no",
                          "positive cells (includes sampling zeros)")) +
    theme_pub
}

# =============================================================================
# GENERATE ALL
# =============================================================================

generate_all_marker_plots <- function(
    output_dir = "Marker/Graphs/All_Comparisons",
    model_dir  = "Marker_Rate_Models",
    stain_handling = "double_only") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL MARKER RATE POST-HOC COMPARISON PLOTS\n")
  cat("=============================================================================\n\n")

  if (!SHOW_UNVERIFIED_STARS) {
    cat("NOTE: stars suppressed on the Region, Sex and interaction panels —\n")
    cat("      those contrasts come from the unverified emmeans path. The\n")
    cat("      Diagnosis panel keeps its stars (marginal_diagnosis).\n\n")
  }
  cat("Loading data and models...\n")
  raw <- load_marker_data(stain_handling)
  for (tag in names(raw)) {
    cat(sprintf("  %-5s %7d islets, %2d donors\n", tag, nrow(raw[[tag]]),
                dplyr::n_distinct(raw[[tag]]$Donor)))
  }
  cat("\n")

  for (def in marker_defs) {
    full <- raw[[sub("_.*", "", def$key)]]
    sub  <- marker_subset(full, def)

    r <- load_rate_model(def, model_dir)
    if (is.null(r)) {
      cat("--- SKIPPING", def$label, "- no cached model ---\n"); next
    }
    model <- r$model

    dir_m <- file.path(output_dir, def$label)
    cat(sprintf("--- %s: %s (%s), %d islets, %d donors ---\n", def$label,
                r$family, r$how, nrow(sub), dplyr::n_distinct(sub$Donor)))

    save_plot(plot_rate_main_diagnosis(sub, model, def),
              "01_main_diagnosis", dir_m)
    save_plot(plot_rate_main_region(sub, model, def),
              "02_main_region", dir_m)
    save_plot(plot_rate_main_sex(sub, model, def),
              "03_main_sex", dir_m)
    save_plot(plot_rate_diagnosis_by_region(sub, model, def),
              "04_diagnosis_by_region", dir_m, width = 8)
    save_plot(plot_rate_region_by_diagnosis(sub, model, def),
              "05_region_by_diagnosis", dir_m, width = 9)
    save_plot(plot_rate_diagnosis_by_sex(sub, model, def),
              "06_diagnosis_by_sex", dir_m, width = 8)
    save_plot(plot_rate_sex_by_diagnosis(sub, model, def),
              "07_sex_by_diagnosis", dir_m, width = 8)
    save_plot(plot_rate_region_by_sex(sub, model, def),
              "08_region_by_sex", dir_m, width = 9)
    save_plot(plot_rate_sex_by_region(sub, model, def),
              "09_sex_by_region", dir_m, width = 10)
    save_plot(plot_rate_age_overall(sub, model, def),
              "10_age_overall", dir_m)
    save_plot(plot_rate_age_by_diagnosis(sub, model, def),
              "11_age_by_diagnosis", dir_m, width = 9)
    save_plot(plot_rate_age_by_region(sub, model, def),
              "12_age_by_region", dir_m, width = 10)

    zi_plot <- plot_rate_zero_inflation(sub, model, def)
    if (!is.null(zi_plot)) {
      save_plot(zi_plot, "13_zero_inflation", dir_m)
    } else {
      cat("  (no zero-inflation component — plot 13 skipped)\n")
    }
    cat("\n")
  }

  cat("=============================================================================\n")
  cat("DONE — plots written under", output_dir, "\n")
  cat("=============================================================================\n")
}
