# =============================================================================
# Generate Comprehensive Excel Results for Count Models (glmmTMB)
# =============================================================================
#
# Models (10 total):
#   Diagnosis: count_ins, count_glu, count_soma, count_pp (nbinom1),
#              count_islet (nbinom2, no log_Islet.Cells_c)
#   T1D (AO):  count_ins_t1d, count_glu_t1d, count_soma_t1d,
#              count_pp_t1d, count_islet_t1d
#
# Note: DD + Age and DD + AO models are algebraically equivalent.
# Only the AO parameterization is used; DD slopes are extracted from it.
#
# All glmmTMB with log link, (1|Donor) + (1|Donor:ImageID)
# Data: quad_data (diagnosis) or quad_t1d (T1D)
#
# Output: 10 Excel files in Count/Results/
#
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(writexl)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

fmt <- function(x, digits = 4) ifelse(is.na(x), NA, round(x, digits))

add_stars <- function(p) {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
}

find_ci <- function(df) {
  nms <- names(df)
  lower <- NULL; upper <- NULL
  for (lc in c("lower.CL", "asymp.LCL", "lower.HPD", "lower.CI"))
    if (lc %in% nms) { lower <- lc; break }
  for (uc in c("upper.CL", "asymp.UCL", "upper.HPD", "upper.CI"))
    if (uc %in% nms) { upper <- uc; break }
  list(lower = lower, upper = upper)
}

# All count models use log link, so emmeans default = log scale
extract_emm_response <- function(em_summary_df) {
  ci <- find_ci(em_summary_df)
  em_summary_df$mean_count <- exp(em_summary_df$emmean)
  em_summary_df$lower.CI   <- exp(em_summary_df[[ci$lower]])
  em_summary_df$upper.CI   <- exp(em_summary_df[[ci$upper]])
  em_summary_df
}

extract_contrast_results <- function(contr_summary_df) {
  contr_summary_df$log_estimate <- contr_summary_df$estimate
  contr_summary_df$pct_diff     <- (exp(contr_summary_df$estimate) - 1) * 100
  contr_summary_df
}

extract_fe <- function(model) {
  fe <- summary(model)$coefficients$cond
  est <- fe[, "Estimate"]
  se  <- fe[, "Std. Error"]
  ci_lo <- est - 1.96 * se
  ci_hi <- est + 1.96 * se
  
  tibble(
    term = rownames(fe), estimate = fmt(est, 4),
    std.error = fmt(se, 4),
    `l-95% CI` = fmt(ci_lo, 4),
    `u-95% CI` = fmt(ci_hi, 4),
    z.value = fmt(fe[, "z value"], 2), p.value = fmt(fe[, "Pr(>|z|)"], 4),
    pct_change = fmt((exp(est) - 1) * 100, 1),
    sig = add_stars(fe[, "Pr(>|z|)"])
  )
}

extract_re <- function(model) {
  vc <- VarCorr(model)$cond
  re_df <- tibble(group = names(vc),
                  variance = fmt(sapply(vc, function(x) x[1]), 4),
                  std_dev = fmt(sqrt(sapply(vc, function(x) x[1])), 4))
  bind_rows(re_df, tibble(group = "Residual",
                          variance = fmt(sigma(model)^2, 4),
                          std_dev = fmt(sigma(model), 4)))
}

format_trend_summary <- function(em_obj, group_cols) {
  em_df <- as.data.frame(summary(em_obj))
  trend_col <- grep("\\.trend$", names(em_df), value = TRUE)
  if (length(trend_col) == 0) trend_col <- grep("trend", names(em_df), value = TRUE, ignore.case = TRUE)
  if (length(trend_col) == 0) { num_cols <- setdiff(names(em_df)[sapply(em_df, is.numeric)], c("SE","df")); trend_col <- num_cols[1] }
  ci <- find_ci(em_df)
  if (is.null(ci$lower) || is.null(ci$upper)) {
    ci_df <- tryCatch(as.data.frame(confint(em_obj)), error = function(e) NULL)
    if (!is.null(ci_df)) { ci2 <- find_ci(ci_df); if (!is.null(ci2$lower)) {
      em_df$lower.CI <- ci_df[[ci2$lower]]; em_df$upper.CI <- ci_df[[ci2$upper]]
      ci <- list(lower = "lower.CI", upper = "upper.CI") } }
  }
  if (is.null(ci$lower) || is.null(ci$upper)) {
    em_df$lower.CI <- em_df[[trend_col]] - 1.96 * em_df$SE
    em_df$upper.CI <- em_df[[trend_col]] + 1.96 * em_df$SE
    ci <- list(lower = "lower.CI", upper = "upper.CI")
  }
  result <- tibble(slope = fmt(em_df[[trend_col]], 4), SE = fmt(em_df$SE, 4),
                   `lower.CI` = fmt(em_df[[ci$lower]], 4), `upper.CI` = fmt(em_df[[ci$upper]], 4))
  for (gc in group_cols) result[[gc]] <- em_df[[gc]]
  result %>% select(all_of(group_cols), everything())
}

format_trend_contrasts <- function(contr_obj) {
  contr_df <- as.data.frame(summary(contr_obj))
  ci <- find_ci(contr_df)
  if (is.null(ci$lower) || is.null(ci$upper)) {
    ci_df <- tryCatch(as.data.frame(confint(contr_obj)), error = function(e) NULL)
    if (!is.null(ci_df)) { ci2 <- find_ci(ci_df); if (!is.null(ci2$lower)) {
      contr_df$lower.CI <- ci_df[[ci2$lower]]; contr_df$upper.CI <- ci_df[[ci2$upper]]
      ci <- list(lower = "lower.CI", upper = "upper.CI") } }
  }
  if (is.null(ci$lower) || is.null(ci$upper)) {
    contr_df$lower.CI <- contr_df$estimate - 1.96 * contr_df$SE
    contr_df$upper.CI <- contr_df$estimate + 1.96 * contr_df$SE
    ci <- list(lower = "lower.CI", upper = "upper.CI")
  }
  tibble(contrast = as.character(contr_df$contrast),
         estimate = fmt(contr_df$estimate, 4), SE = fmt(contr_df$SE, 4),
         `lower.CI` = fmt(contr_df[[ci$lower]], 4), `upper.CI` = fmt(contr_df[[ci$upper]], 4),
         p.value = fmt(contr_df$p.value, 4), sig = add_stars(contr_df$p.value))
}

# Convenience: emmeans with explicit data
emm <- function(model, specs, data, ...) emmeans(model, specs, data = data, ...)
emt <- function(model, specs, var, data, ...) emtrends(model, specs, var = var, data = data, ...)

# Standard block: marginal means + pairwise contrasts
get_means_contrasts <- function(model, spec, data, adjust = "tukey", reverse = FALSE) {
  em <- emm(model, spec, data)
  em_df <- extract_emm_response(as.data.frame(summary(em)))
  ci <- find_ci(em_df)
  fac_cols <- setdiff(names(em_df), c("emmean", "SE", "df", ci$lower, ci$upper,
                                       "mean_count", "lower.CI", "upper.CI"))
  means <- em_df[c(fac_cols, "mean_count", "lower.CI", "upper.CI")] %>%
    mutate(mean_count = fmt(mean_count, 2), lower.CI = fmt(lower.CI, 2), upper.CI = fmt(upper.CI, 2))

  contr <- pairs(em, adjust = adjust, reverse = reverse)
  contr_df <- extract_contrast_results(as.data.frame(summary(contr)))
  contrasts <- tibble(contrast = as.character(contr_df$contrast),
                      log_estimate = fmt(contr_df$log_estimate, 4), SE = fmt(contr_df$SE, 4),
                      p.value = fmt(contr_df$p.value, 4), pct_diff = fmt(contr_df$pct_diff, 1),
                      sig = add_stars(contr_df$p.value))
  # Add conditioning variable if present
  contr_summ <- as.data.frame(summary(contr))
  for (fc in fac_cols) {
    if (fc %in% names(contr_summ)) contrasts[[fc]] <- contr_summ[[fc]]
  }
  list(means = means, contrasts = contrasts, em = em)
}


# =============================================================================
# 1. DIAGNOSIS MODEL
# =============================================================================
#
# Cell counts: (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
#              Diagnosis:Region:log_Islet.Cells_c + Diagnosis:Region:Age_c
# Islet count: (Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c
#
# Sheets:
#   1.  Model Summary
#   2.  Diagnosis
#   3.  Region
#   4.  Sex
#   5.  Diagnosis x Region
#   6.  Diagnosis x Sex
#   7.  Region x Sex
#   8.  Age Slopes
#   9.  Islet Size Slopes (cell count models only)
#  10.  Diag x Region x Age (3-way)
#  11.  Diag x Region x Islet (3-way, cell count models only)
#  12.  Random Effects
# =============================================================================

generate_count_diag_excel <- function(model, data, cell_label, response_col,
                                      has_islet = TRUE,
                                      output_dir = "Count/Results") {

  cat("\n=============================================================================\n")
  cat("GENERATING EXCEL: COUNT —", toupper(cell_label), "DIAGNOSIS\n")
  cat("=============================================================================\n\n")

  fam <- family(model)$family
  formula_text <- if (has_islet) {
    paste0(response_col, " ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +\n",
           "    Diagnosis:Region:log_Islet.Cells_c + Diagnosis:Region:Age_c + (1|Donor) + (1|Donor:ImageID)")
  } else {
    paste0(response_col, " ~ (Diagnosis + Region + Sex + Age_c)^2 +\n",
           "    Diagnosis:Region:Age_c + (1|Donor) + (1|Donor:ImageID)")
  }

  # SHEET 1
  cat("Creating Sheet 1: Model Summary...\n")
  fe_df <- extract_fe(model)
  model_notes <- tibble(Note = c(
    paste("CELL COUNT ANALYSIS:", toupper(cell_label), "— DIAGNOSIS MODEL"),
    "", "Generalized linear mixed model (glmmTMB)",
    paste("Response:", response_col), paste("Family:", fam, "(log link)"),
    paste("Formula:", formula_text), "",
    "Reference levels: ND (Diagnosis), Head (Region), Female (Sex)",
    "Coefficients on log scale. Exponentiate for rate ratios.",
    "% change = (exp(estimate) - 1) * 100", "",
    "Multiplicity: Tukey for Region pairwise, Holm for simple effects, None for single contrasts.",
    paste("N observations:", nrow(data)), paste("N donors:", n_distinct(data$Donor)), ""))
  model_summary <- bind_rows(model_notes, fe_df %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEETS 2-4: Main effects
  cat("Creating Sheets 2-4: Main effects...\n")
  diag_res    <- get_means_contrasts(model, ~ Diagnosis, data, "none", reverse = TRUE)
  region_res  <- get_means_contrasts(model, ~ Region, data, "tukey")
  sex_res     <- get_means_contrasts(model, ~ Sex, data, "none")

  make_main_sheet <- function(title, note, means, contrasts, contr_label) {
    bind_rows(
      tibble(Note = title), tibble(Note = note), tibble(Note = NA),
      tibble(Note = "MARGINAL MEANS (predicted count)"),
      means %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = contr_label),
      contrasts %>% mutate(Note = NA) %>% select(Note, everything()))
  }

  diagnosis_sheet <- make_main_sheet("DIAGNOSIS EFFECT",
    paste("Main effect of Diagnosis on", cell_label, "count, averaged over Region and Sex."),
    diag_res$means, diag_res$contrasts, "CONTRAST (T1D vs ND, no adjustment)")
  region_sheet <- make_main_sheet("REGION EFFECT",
    paste("Main effect of Region on", cell_label, "count, averaged over Diagnosis and Sex."),
    region_res$means, region_res$contrasts, "PAIRWISE CONTRASTS (Tukey-adjusted)")
  sex_sheet <- make_main_sheet("SEX EFFECT",
    paste("Main effect of Sex on", cell_label, "count, averaged over Diagnosis and Region."),
    sex_res$means, sex_res$contrasts, "CONTRAST (Female vs Male, no adjustment)")

  # SHEET 5: Diagnosis x Region
  cat("Creating Sheet 5: Diagnosis x Region...\n")
  dr_means <- get_means_contrasts(model, ~ Diagnosis * Region, data, "none")$means
  dr_diag  <- get_means_contrasts(model, ~ Diagnosis | Region, data, "holm", reverse = TRUE)
  dr_reg   <- get_means_contrasts(model, ~ Region | Diagnosis, data, "tukey")

  diag_region_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x REGION INTERACTION"),
    tibble(Note = "Marginal means and contrasts averaged over Sex."), tibble(Note = NA),
    tibble(Note = "CELL MEANS (predicted count)"),
    dr_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT WITHIN EACH REGION (T1D vs ND, Holm-adjusted)"),
    dr_diag$contrasts %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS WITHIN EACH DIAGNOSIS (Tukey-adjusted)"),
    dr_reg$contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 6: Diagnosis x Sex
  cat("Creating Sheet 6: Diagnosis x Sex...\n")
  ds_means <- get_means_contrasts(model, ~ Diagnosis * Sex, data, "none")$means
  ds_diag  <- get_means_contrasts(model, ~ Diagnosis | Sex, data, "holm", reverse = TRUE)
  ds_sex   <- get_means_contrasts(model, ~ Sex | Diagnosis, data, "holm")

  diag_sex_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x SEX INTERACTION"),
    tibble(Note = "Marginal means and contrasts averaged over Region."), tibble(Note = NA),
    tibble(Note = "CELL MEANS (predicted count)"),
    ds_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT WITHIN EACH SEX (T1D vs ND, Holm-adjusted)"),
    ds_diag$contrasts %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SEX EFFECT WITHIN EACH DIAGNOSIS (Female vs Male, Holm-adjusted)"),
    ds_sex$contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 7: Region x Sex
  cat("Creating Sheet 7: Region x Sex...\n")
  rs_means <- get_means_contrasts(model, ~ Region * Sex, data, "none")$means
  rs_sex   <- get_means_contrasts(model, ~ Sex | Region, data, "holm")
  rs_reg   <- get_means_contrasts(model, ~ Region | Sex, data, "tukey")

  region_sex_sheet <- bind_rows(
    tibble(Note = "REGION x SEX INTERACTION"),
    tibble(Note = "Marginal means and contrasts averaged over Diagnosis."), tibble(Note = NA),
    tibble(Note = "CELL MEANS (predicted count)"),
    rs_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SEX EFFECT WITHIN EACH REGION (Female vs Male, Holm-adjusted)"),
    rs_sex$contrasts %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS WITHIN EACH SEX (Tukey-adjusted)"),
    rs_reg$contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 8: Age Slopes
  cat("Creating Sheet 8: Age Slopes...\n")
  age_overall <- format_trend_summary(emt(model, ~ 1, "Age_c", data), character(0))
  age_diag    <- format_trend_summary(emt(model, ~ Diagnosis, "Age_c", data), "Diagnosis")
  age_diag_c  <- format_trend_contrasts(pairs(emt(model, ~ Diagnosis, "Age_c", data), adjust = "none"))
  age_region  <- format_trend_summary(emt(model, ~ Region, "Age_c", data), "Region")
  age_region_c <- format_trend_contrasts(pairs(emt(model, ~ Region, "Age_c", data), adjust = "tukey"))
  age_sex     <- format_trend_summary(emt(model, ~ Sex, "Age_c", data), "Sex")
  age_sex_c   <- format_trend_contrasts(pairs(emt(model, ~ Sex, "Age_c", data), adjust = "none"))

  age_sheet <- bind_rows(
    tibble(Note = "AGE SLOPES"),
    tibble(Note = "Slopes on log scale. % change per year = (exp(slope) - 1) * 100."),
    tibble(Note = NA),
    tibble(Note = "OVERALL AGE SLOPE"),
    age_overall %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "AGE SLOPES BY DIAGNOSIS"),
    age_diag %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "Diagnosis contrast on age slopes:"),
    age_diag_c %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "AGE SLOPES BY REGION"),
    age_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "Region contrasts on age slopes (Tukey):"),
    age_region_c %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "AGE SLOPES BY SEX"),
    age_sex %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "Sex contrast on age slopes:"),
    age_sex_c %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 9: Islet Size Slopes (cell count models only)
  islet_sheet <- NULL
  if (has_islet) {
    cat("Creating Sheet 9: Islet Size Slopes...\n")
    islet_overall <- format_trend_summary(emt(model, ~ 1, "log_Islet.Cells_c", data), character(0))
    islet_diag    <- format_trend_summary(emt(model, ~ Diagnosis, "log_Islet.Cells_c", data), "Diagnosis")
    islet_diag_c  <- format_trend_contrasts(pairs(emt(model, ~ Diagnosis, "log_Islet.Cells_c", data), adjust = "none"))
    islet_region  <- format_trend_summary(emt(model, ~ Region, "log_Islet.Cells_c", data), "Region")
    islet_region_c <- format_trend_contrasts(pairs(emt(model, ~ Region, "log_Islet.Cells_c", data), adjust = "tukey"))
    islet_sex     <- format_trend_summary(emt(model, ~ Sex, "log_Islet.Cells_c", data), "Sex")
    islet_sex_c   <- format_trend_contrasts(pairs(emt(model, ~ Sex, "log_Islet.Cells_c", data), adjust = "none"))

    islet_sheet <- bind_rows(
      tibble(Note = "ISLET SIZE SLOPES"),
      tibble(Note = "Change in log(count) per 1 unit increase in log(Islet.Cells)."),
      tibble(Note = "1 unit on log scale ≈ 2.72-fold increase in cell count."),
      tibble(Note = NA),
      tibble(Note = "OVERALL ISLET SIZE SLOPE"),
      islet_overall %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "ISLET SLOPES BY DIAGNOSIS"),
      islet_diag %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "Diagnosis contrast on islet slopes:"),
      islet_diag_c %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "ISLET SLOPES BY REGION"),
      islet_region %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "Region contrasts on islet slopes (Tukey):"),
      islet_region_c %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "ISLET SLOPES BY SEX"),
      islet_sex %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "Sex contrast on islet slopes:"),
      islet_sex_c %>% mutate(Note = NA) %>% select(Note, everything()))
  }

  # SHEET 10: Diag x Region x Age (3-way)
  cat("Creating Sheet 10: Diag x Region x Age...\n")
  age_dr     <- format_trend_summary(emt(model, ~ Diagnosis * Region, "Age_c", data), c("Diagnosis", "Region"))
  age_dr_d   <- emt(model, ~ Diagnosis | Region, "Age_c", data)
  age_dr_d_c <- format_trend_contrasts(pairs(age_dr_d, adjust = "holm", reverse = TRUE))
  age_dr_d_c$Region <- as.data.frame(summary(pairs(age_dr_d, adjust = "holm", reverse = TRUE)))$Region
  age_dr_r   <- emt(model, ~ Region | Diagnosis, "Age_c", data)
  age_dr_r_c <- format_trend_contrasts(pairs(age_dr_r, adjust = "tukey"))
  age_dr_r_c$Diagnosis <- as.data.frame(summary(pairs(age_dr_r, adjust = "tukey")))$Diagnosis

  diag_region_age_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x REGION x AGE SLOPES (3-way)"),
    tibble(Note = "Age slopes per Diagnosis × Region cell, marginalized over Sex."),
    tibble(Note = NA), tibble(Note = "AGE SLOPES PER CELL"),
    age_dr %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT ON AGE SLOPES WITHIN EACH REGION (Holm-adjusted)"),
    age_dr_d_c %>% mutate(Note = NA) %>% select(Note, Region, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS ON AGE SLOPES WITHIN EACH DIAGNOSIS (Tukey-adjusted)"),
    age_dr_r_c %>% mutate(Note = NA) %>% select(Note, Diagnosis, everything()))

  # SHEET 11: Diag x Region x Islet (3-way, cell count models only)
  diag_region_islet_sheet <- NULL
  if (has_islet) {
    cat("Creating Sheet 11: Diag x Region x Islet...\n")
    islet_dr     <- format_trend_summary(emt(model, ~ Diagnosis * Region, "log_Islet.Cells_c", data), c("Diagnosis", "Region"))
    islet_dr_d   <- emt(model, ~ Diagnosis | Region, "log_Islet.Cells_c", data)
    islet_dr_d_c <- format_trend_contrasts(pairs(islet_dr_d, adjust = "holm", reverse = TRUE))
    islet_dr_d_c$Region <- as.data.frame(summary(pairs(islet_dr_d, adjust = "holm", reverse = TRUE)))$Region
    islet_dr_r   <- emt(model, ~ Region | Diagnosis, "log_Islet.Cells_c", data)
    islet_dr_r_c <- format_trend_contrasts(pairs(islet_dr_r, adjust = "tukey"))
    islet_dr_r_c$Diagnosis <- as.data.frame(summary(pairs(islet_dr_r, adjust = "tukey")))$Diagnosis

    diag_region_islet_sheet <- bind_rows(
      tibble(Note = "DIAGNOSIS x REGION x ISLET SIZE SLOPES (3-way)"),
      tibble(Note = "Islet size slopes per Diagnosis × Region cell, marginalized over Sex."),
      tibble(Note = NA), tibble(Note = "ISLET SLOPES PER CELL"),
      islet_dr %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA),
      tibble(Note = "DIAGNOSIS EFFECT ON ISLET SLOPES WITHIN EACH REGION (Holm-adjusted)"),
      islet_dr_d_c %>% mutate(Note = NA) %>% select(Note, Region, everything()),
      tibble(Note = NA),
      tibble(Note = "REGION CONTRASTS ON ISLET SLOPES WITHIN EACH DIAGNOSIS (Tukey-adjusted)"),
      islet_dr_r_c %>% mutate(Note = NA) %>% select(Note, Diagnosis, everything()))
  }

  # SHEET 12: Random Effects
  cat("Creating Random Effects sheet...\n")
  re_sheet <- bind_rows(
    tibble(Note = "RANDOM EFFECTS VARIANCE COMPONENTS"), tibble(Note = NA),
    extract_re(model) %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "SAMPLE SIZES"),
    tibble(n_obs = nrow(data), n_donors = n_distinct(data$Donor),
           n_images = n_distinct(data$ImageID)) %>% mutate(Note = NA) %>% select(Note, everything()))

  # ASSEMBLE
  cat("Assembling workbook...\n")
  sheets <- list(
    "Model Summary"       = as.data.frame(model_summary),
    "Diagnosis"           = as.data.frame(diagnosis_sheet),
    "Region"              = as.data.frame(region_sheet),
    "Sex"                 = as.data.frame(sex_sheet),
    "Diagnosis x Region"  = as.data.frame(diag_region_sheet),
    "Diagnosis x Sex"     = as.data.frame(diag_sex_sheet),
    "Region x Sex"        = as.data.frame(region_sex_sheet),
    "Age Slopes"          = as.data.frame(age_sheet)
  )
  if (!is.null(islet_sheet)) sheets[["Islet Size Slopes"]] <- as.data.frame(islet_sheet)
  sheets[["Diag x Region x Age"]] <- as.data.frame(diag_region_age_sheet)
  if (!is.null(diag_region_islet_sheet)) sheets[["Diag x Region x Islet"]] <- as.data.frame(diag_region_islet_sheet)
  sheets[["Random Effects"]] <- as.data.frame(re_sheet)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0(cell_label, "_Count_Diagnosis_Analysis.xlsx"))
  write_xlsx(sheets, output_file)
  cat("Saved to:", output_file, "\n")
  return(invisible(sheets))
}


# =============================================================================
# 2. DD / AO MODELS (T1D only)
# =============================================================================
#
# Cell counts: (DD + Region + Sex + Age/AO + log_Islet.Cells_c)^2 +
#              DD:Region:log_Islet.Cells_c + DD:Region:Age/AO
# Islet count: (DD + Region + Sex + Age/AO)^2 + DD:Region:Age/AO
#
# Sheets:
#   1.  Model Summary
#   2.  Region
#   3.  Sex
#   4.  Region x Sex
#   5.  DD Slopes
#   6.  Age / AO Slopes
#   7.  Islet Size Slopes (cell count models only)
#   8.  DD x Region x Age/AO (3-way)
#   9.  DD x Region x Islet (3-way, cell count models only)
#  10.  Random Effects
# =============================================================================

generate_count_t1d_excel <- function(model, data, cell_label, response_col,
                                     has_islet = TRUE, model_type = "dd",
                                     output_dir = "Count/Results") {

  dd_var    <- "Disease.Duration_c"
  age_var   <- if (model_type == "dd") "Age_c" else "age_at_onset_c"
  age_label <- if (model_type == "dd") "Age" else "Age at Onset"
  type_label <- if (model_type == "dd") "Disease Duration" else "Age at Onset"
  type_tag   <- toupper(model_type)

  cat("\n=============================================================================\n")
  cat("GENERATING EXCEL: COUNT —", toupper(cell_label), type_tag, "(T1D)\n")
  cat("=============================================================================\n\n")

  fam <- family(model)$family
  formula_text <- if (has_islet) {
    paste0(response_col, " ~ (", dd_var, " + Region + Sex + ", age_var, " + log_Islet.Cells_c)^2 +\n",
           "    ", dd_var, ":Region:log_Islet.Cells_c + ", dd_var, ":Region:", age_var,
           " + (1|Donor) + (1|Donor:ImageID)")
  } else {
    paste0(response_col, " ~ (", dd_var, " + Region + Sex + ", age_var, ")^2 +\n",
           "    ", dd_var, ":Region:", age_var, " + (1|Donor) + (1|Donor:ImageID)")
  }

  # SHEET 1
  cat("Creating Sheet 1: Model Summary...\n")
  fe_df <- extract_fe(model)
  model_notes <- tibble(Note = c(
    paste("CELL COUNT ANALYSIS:", toupper(cell_label), "—", toupper(type_label), "(T1D only)"),
    "", "Generalized linear mixed model (glmmTMB)",
    paste("Response:", response_col), paste("Family:", fam, "(log link)"),
    paste("Formula:", formula_text), "",
    "Reference levels: Head (Region), Female (Sex)",
    "Coefficients on log scale. % change = (exp(estimate) - 1) * 100.",
    paste("N observations:", nrow(data)), paste("N donors:", n_distinct(data$Donor)), ""))
  model_summary <- bind_rows(model_notes, fe_df %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEETS 2-3: Region, Sex
  cat("Creating Sheets 2-3: Region and Sex...\n")
  region_res <- get_means_contrasts(model, ~ Region, data, "tukey")
  sex_res    <- get_means_contrasts(model, ~ Sex, data, "none")

  region_sheet <- bind_rows(
    tibble(Note = "REGION EFFECT (T1D only)"),
    tibble(Note = "Averaged over Sex."), tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (predicted count)"),
    region_res$means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "PAIRWISE CONTRASTS (Tukey-adjusted)"),
    region_res$contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  sex_sheet <- bind_rows(
    tibble(Note = "SEX EFFECT (T1D only)"),
    tibble(Note = "Averaged over Region."), tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (predicted count)"),
    sex_res$means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "CONTRAST (Female vs Male)"),
    sex_res$contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 4: Region x Sex
  cat("Creating Sheet 4: Region x Sex...\n")
  rs_means <- get_means_contrasts(model, ~ Region * Sex, data, "none")$means
  rs_sex   <- get_means_contrasts(model, ~ Sex | Region, data, "holm")
  rs_reg   <- get_means_contrasts(model, ~ Region | Sex, data, "tukey")

  region_sex_sheet <- bind_rows(
    tibble(Note = "REGION x SEX INTERACTION"), tibble(Note = NA),
    tibble(Note = "CELL MEANS (predicted count)"),
    rs_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SEX EFFECT WITHIN EACH REGION (Holm-adjusted)"),
    rs_sex$contrasts %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS WITHIN EACH SEX (Tukey-adjusted)"),
    rs_reg$contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 5: DD Slopes
  cat("Creating Sheet 5: DD Slopes...\n")
  dd_overall  <- format_trend_summary(emt(model, ~ 1, dd_var, data), character(0))
  dd_region   <- format_trend_summary(emt(model, ~ Region, dd_var, data), "Region")
  dd_region_c <- format_trend_contrasts(pairs(emt(model, ~ Region, dd_var, data), adjust = "tukey"))
  dd_sex      <- format_trend_summary(emt(model, ~ Sex, dd_var, data), "Sex")
  dd_sex_c    <- format_trend_contrasts(pairs(emt(model, ~ Sex, dd_var, data), adjust = "none"))

  dd_sheet <- bind_rows(
    tibble(Note = "DISEASE DURATION SLOPES (T1D only)"),
    tibble(Note = "Slopes on log scale. % change per year = (exp(slope) - 1) * 100."),
    tibble(Note = NA),
    tibble(Note = "OVERALL DD SLOPE"),
    dd_overall %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "DD SLOPES BY REGION"),
    dd_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "Region contrasts on DD slopes (Tukey):"),
    dd_region_c %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "DD SLOPES BY SEX"),
    dd_sex %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "Sex contrast on DD slopes:"),
    dd_sex_c %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 6: Age / AO Slopes
  cat("Creating Sheet 6:", age_label, "Slopes...\n")
  age_overall  <- format_trend_summary(emt(model, ~ 1, age_var, data), character(0))
  age_region   <- format_trend_summary(emt(model, ~ Region, age_var, data), "Region")
  age_region_c <- format_trend_contrasts(pairs(emt(model, ~ Region, age_var, data), adjust = "tukey"))
  age_sex      <- format_trend_summary(emt(model, ~ Sex, age_var, data), "Sex")
  age_sex_c    <- format_trend_contrasts(pairs(emt(model, ~ Sex, age_var, data), adjust = "none"))

  age_sheet <- bind_rows(
    tibble(Note = paste(toupper(age_label), "SLOPES (T1D only)")),
    tibble(Note = "Slopes on log scale."), tibble(Note = NA),
    tibble(Note = paste("OVERALL", toupper(age_label), "SLOPE")),
    age_overall %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = paste(toupper(age_label), "SLOPES BY REGION")),
    age_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = paste("Region contrasts on", age_label, "slopes (Tukey):")),
    age_region_c %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = paste(toupper(age_label), "SLOPES BY SEX")),
    age_sex %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = paste("Sex contrast on", age_label, "slopes:")),
    age_sex_c %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 7: Islet Size Slopes (cell count models only)
  islet_sheet <- NULL
  if (has_islet) {
    cat("Creating Sheet 7: Islet Size Slopes...\n")
    isl_overall  <- format_trend_summary(emt(model, ~ 1, "log_Islet.Cells_c", data), character(0))
    isl_region   <- format_trend_summary(emt(model, ~ Region, "log_Islet.Cells_c", data), "Region")
    isl_region_c <- format_trend_contrasts(pairs(emt(model, ~ Region, "log_Islet.Cells_c", data), adjust = "tukey"))
    isl_sex      <- format_trend_summary(emt(model, ~ Sex, "log_Islet.Cells_c", data), "Sex")
    isl_sex_c    <- format_trend_contrasts(pairs(emt(model, ~ Sex, "log_Islet.Cells_c", data), adjust = "none"))

    islet_sheet <- bind_rows(
      tibble(Note = "ISLET SIZE SLOPES (T1D only)"),
      tibble(Note = "1 unit on log scale ≈ 2.72-fold increase in cell count."),
      tibble(Note = NA),
      tibble(Note = "OVERALL ISLET SIZE SLOPE"),
      isl_overall %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "ISLET SLOPES BY REGION"),
      isl_region %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "Region contrasts on islet slopes (Tukey):"),
      isl_region_c %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "ISLET SLOPES BY SEX"),
      isl_sex %>% mutate(Note = NA) %>% select(Note, everything()),
      tibble(Note = NA), tibble(Note = "Sex contrast on islet slopes:"),
      isl_sex_c %>% mutate(Note = NA) %>% select(Note, everything()))
  }

  # SHEET 8: DD x Region x Age/AO (3-way)
  cat("Creating Sheet 8: DD x Region x", age_label, "...\n")
  # DD slopes at low and high age
  dd_at_lo <- format_trend_summary(emt(model, ~ Region, dd_var, data, at = setNames(list(-1), age_var)), "Region") %>%
    mutate(level = paste(age_label, "= mean - 1"))
  dd_at_hi <- format_trend_summary(emt(model, ~ Region, dd_var, data, at = setNames(list(1), age_var)), "Region") %>%
    mutate(level = paste(age_label, "= mean + 1"))
  dd_at_age <- bind_rows(dd_at_lo, dd_at_hi) %>% select(Region, level, everything())

  three_way_age_sheet <- bind_rows(
    tibble(Note = paste("DD × REGION ×", toupper(age_label), "(3-way)")),
    tibble(Note = paste("DD slopes evaluated at", age_label, "= mean ± 1, by Region.")),
    tibble(Note = NA),
    dd_at_age %>% mutate(Note = NA) %>% select(Note, everything()))

  # SHEET 9: DD x Region x Islet (3-way, cell count models only)
  three_way_islet_sheet <- NULL
  if (has_islet) {
    cat("Creating Sheet 9: DD x Region x Islet...\n")
    dd_at_isl_lo <- format_trend_summary(
      emt(model, ~ Region, dd_var, data, at = list(log_Islet.Cells_c = -1)), "Region") %>%
      mutate(level = "Islet = mean - 1")
    dd_at_isl_hi <- format_trend_summary(
      emt(model, ~ Region, dd_var, data, at = list(log_Islet.Cells_c = 1)), "Region") %>%
      mutate(level = "Islet = mean + 1")
    dd_at_islet <- bind_rows(dd_at_isl_lo, dd_at_isl_hi) %>% select(Region, level, everything())

    three_way_islet_sheet <- bind_rows(
      tibble(Note = "DD × REGION × ISLET SIZE (3-way)"),
      tibble(Note = "DD slopes evaluated at log_Islet.Cells_c = mean ± 1, by Region."),
      tibble(Note = NA),
      dd_at_islet %>% mutate(Note = NA) %>% select(Note, everything()))
  }

  # SHEET 10: Random Effects
  cat("Creating Random Effects sheet...\n")
  re_sheet <- bind_rows(
    tibble(Note = "RANDOM EFFECTS VARIANCE COMPONENTS"), tibble(Note = NA),
    extract_re(model) %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "SAMPLE SIZES"),
    tibble(n_obs = nrow(data), n_donors = n_distinct(data$Donor)) %>%
      mutate(Note = NA) %>% select(Note, everything()))

  # ASSEMBLE
  cat("Assembling workbook...\n")
  age_sheet_name <- paste(age_label, "Slopes")

  sheets <- list(
    "Model Summary" = as.data.frame(model_summary),
    "Region"        = as.data.frame(region_sheet),
    "Sex"           = as.data.frame(sex_sheet),
    "Region x Sex"  = as.data.frame(region_sex_sheet)
  )
  sheets[["DD Slopes"]] <- as.data.frame(dd_sheet)
  sheets[[age_sheet_name]] <- as.data.frame(age_sheet)
  if (!is.null(islet_sheet)) sheets[["Islet Size Slopes"]] <- as.data.frame(islet_sheet)
  three_way_age_name <- substr(paste("DD x Reg x", age_label), 1, 31)
  sheets[[three_way_age_name]] <- as.data.frame(three_way_age_sheet)
  if (!is.null(three_way_islet_sheet)) sheets[["DD x Reg x Islet"]] <- as.data.frame(three_way_islet_sheet)
  sheets[["Random Effects"]] <- as.data.frame(re_sheet)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0(cell_label, "_Count_", type_tag, "_Analysis.xlsx"))
  write_xlsx(sheets, output_file)
  cat("Saved to:", output_file, "\n")
  return(invisible(sheets))
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

generate_all_count_excel <- function(output_dir = "Count/Results") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL COUNT EXCEL RESULTS\n")
  cat("=============================================================================\n\n")

  # Load data
  cat("Loading data...\n")
  quad_data <- readRDS("Data/quad_data.rds")
  quad_data <- quad_data %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
      Region    = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor),
      ImageID   = factor(ImageID)
    )
  quad_t1d <- readRDS("Data/quad_data_t1d.rds")
  quad_t1d <- quad_t1d %>%
    mutate(
      Region  = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex     = factor(Sex, levels = c("Female", "Male")),
      Donor   = factor(Donor),
      ImageID = factor(ImageID)
    ) %>% droplevels()

  cat("quad_data:", nrow(quad_data), "rows,", n_distinct(quad_data$Donor), "donors\n")
  cat("quad_t1d: ", nrow(quad_t1d), "rows,", n_distinct(quad_t1d$Donor), "donors\n\n")

  # Cell type definitions
  cell_defs <- list(
    list(label = "Insulin",       response = "Count.ins",  file_prefix = "count_ins",    has_islet = TRUE),
    list(label = "Glucagon",      response = "Count.glu",  file_prefix = "count_glu",    has_islet = TRUE),
    list(label = "Somatostatin",  response = "Count.soma", file_prefix = "count_soma",   has_islet = TRUE),
    list(label = "PP",            response = "Count.PP",   file_prefix = "count_pp",     has_islet = TRUE),
    list(label = "Islet_Cells",   response = "Islet.Cells",file_prefix = "count_islet",  has_islet = FALSE)
  )

  # --- Diagnosis ---
  cat("--- DIAGNOSIS MODELS ---\n")
  for (cd in cell_defs) {
    model_path <- paste0("Count/Models/", cd$file_prefix, ".rds")
    cat("Loading", model_path, "...\n")
    model <- readRDS(model_path)
    generate_count_diag_excel(model, quad_data, cd$label, cd$response,
                              has_islet = cd$has_islet, output_dir = output_dir)
    rm(model); gc()
  }

  # --- T1D (AO parameterization — DD slopes extracted from same model) ---
  cat("\n--- T1D MODELS (AO parameterization) ---\n")
  for (cd in cell_defs) {
    model_path <- paste0("Count/Models/", cd$file_prefix, "_t1d.rds")
    cat("Loading", model_path, "...\n")
    model <- readRDS(model_path)
    generate_count_t1d_excel(model, quad_t1d, cd$label, cd$response,
                             has_islet = cd$has_islet, model_type = "ao", output_dir = output_dir)
    rm(model); gc()
  }

  cat("\n=============================================================================\n")
  cat("COMPLETE! 10 Excel files created in", output_dir, ":\n\n")
  cat("DIAGNOSIS (5):\n")
  for (cd in cell_defs) cat("  ", cd$label, "_Count_Diagnosis_Analysis.xlsx\n")
  cat("\nT1D — AO (5):\n")
  for (cd in cell_defs) cat("  ", cd$label, "_Count_AO_Analysis.xlsx\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all count Excel files, run:\n")
  cat("  generate_all_count_excel()\n")
}
