# =============================================================================
# Generate Comprehensive Excel Results for Density Models (glmmTMB)
# =============================================================================
#
# Models:
#   1. density_diag  — Diagnosis (ND vs T1D), all donors
#   2. density_t1d   — T1D-only model (DD + AO parameterization)
#      (replaces separate DD and AO models — they are algebraically equivalent
#       reparameterizations since Age = age_at_onset + Disease.Duration)
#
# All models: Gamma(log link), response = islet_density, random = (1|Donor)
#
# Output: 2 Excel files in Density/Results/
#
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(writexl)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

fmt <- function(x, digits = 4) {
  ifelse(is.na(x), NA, round(x, digits))
}

add_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ ""
  )
}

# All density models are Gamma(log), so emmeans default = log scale
# Exponentiate for response-scale means; contrasts are log-ratios
extract_emm_response <- function(em_summary_df) {
  ci <- find_ci(em_summary_df)
  em_summary_df$mean_density <- exp(em_summary_df$emmean)
  em_summary_df$lower.CI     <- exp(em_summary_df[[ci$lower]])
  em_summary_df$upper.CI     <- exp(em_summary_df[[ci$upper]])
  em_summary_df
}

extract_contrast_results <- function(contr_summary_df) {
  contr_summary_df$log_estimate <- contr_summary_df$estimate
  contr_summary_df$pct_diff     <- (exp(contr_summary_df$estimate) - 1) * 100
  contr_summary_df
}

extract_fe <- function(model) {
  fe <- summary(model)$coefficients$cond
  ci <- tryCatch({
    confint(model, method = "Wald", parm = "beta_")
  }, error = function(e) {
    est <- fe[, "Estimate"]
    se  <- fe[, "Std. Error"]
    cbind(est - 1.96 * se, est + 1.96 * se)
  })
  
  tibble(
    term      = rownames(fe),
    estimate  = fmt(fe[, "Estimate"], 4),
    std.error = fmt(fe[, "Std. Error"], 4),
    `l-95% CI` = fmt(ci[rownames(fe), 1], 4),
    `u-95% CI` = fmt(ci[rownames(fe), 2], 4),
    z.value   = fmt(fe[, "z value"], 2),
    p.value   = fmt(fe[, "Pr(>|z|)"], 4),
    pct_change = fmt((exp(fe[, "Estimate"]) - 1) * 100, 1),
    sig       = add_stars(fe[, "Pr(>|z|)"])
  )
}

extract_re <- function(model, family_name = "gaussian") {
  vc <- VarCorr(model)$cond
  re_df <- tibble(
    group    = names(vc),
    variance = fmt(sapply(vc, function(x) x[1]), 4),
    std_dev  = fmt(sqrt(sapply(vc, function(x) x[1])), 4)
  )
  
  # For Gamma and NB families, sigma() is a dispersion parameter, not a
  # residual SD.  Label accordingly so Excel readers aren't misled.
  if (grepl("gamma", family_name, ignore.case = TRUE)) {
    disp <- sigma(model)
    re_df <- bind_rows(re_df, tibble(
      group    = "Gamma dispersion (1/shape)",
      variance = fmt(disp, 4),
      std_dev  = NA
    ))
  } else {
    re_df <- bind_rows(re_df, tibble(
      group    = "Residual",
      variance = fmt(sigma(model)^2, 4),
      std_dev  = fmt(sigma(model), 4)
    ))
  }
  re_df
}

# Helper: find CI columns in an emmeans/emtrends data frame
find_ci <- function(df) {
  nms <- names(df)
  lower <- NULL; upper <- NULL
  
  lower_candidates <- c("lower.CL", "asymp.LCL", "lower.HPD", "lower.CI")
  upper_candidates <- c("upper.CL", "asymp.UCL", "upper.HPD", "upper.CI")
  
  for (lc in lower_candidates) {
    if (lc %in% nms) { lower <- lc; break }
  }
  for (uc in upper_candidates) {
    if (uc %in% nms) { upper <- uc; break }
  }
  
  list(lower = lower, upper = upper)
}

# Helper: find the trend column in an emtrends summary
find_trend_col <- function(df) {
  nms <- names(df)
  trend_col <- grep("\\.trend$", nms, value = TRUE)
  if (length(trend_col) == 0) trend_col <- grep("trend", nms, value = TRUE, ignore.case = TRUE)
  if (length(trend_col) == 0) return(NULL)
  trend_col[1]
}

# Helper: build a contrast tibble from emtrends pairs
format_trend_contrasts <- function(contr_obj) {
  contr_df <- as.data.frame(summary(contr_obj))
  ci <- find_ci(contr_df)
  
  if (is.null(ci$lower) || is.null(ci$upper)) {
    ci_df <- tryCatch(as.data.frame(confint(contr_obj)), error = function(e) NULL)
    if (!is.null(ci_df)) {
      ci <- find_ci(ci_df)
      if (!is.null(ci$lower)) {
        contr_df$lower.CI <- ci_df[[ci$lower]]
        contr_df$upper.CI <- ci_df[[ci$upper]]
        ci <- list(lower = "lower.CI", upper = "upper.CI")
      }
    }
  }
  
  if (is.null(ci$lower) || is.null(ci$upper)) {
    contr_df$lower.CI <- contr_df$estimate - 1.96 * contr_df$SE
    contr_df$upper.CI <- contr_df$estimate + 1.96 * contr_df$SE
    ci <- list(lower = "lower.CI", upper = "upper.CI")
  }
  
  tibble(
    contrast     = as.character(contr_df$contrast),
    estimate     = fmt(contr_df$estimate, 4),
    SE           = fmt(contr_df$SE, 4),
    `lower.CI`   = fmt(contr_df[[ci$lower]], 4),
    `upper.CI`   = fmt(contr_df[[ci$upper]], 4),
    p.value      = fmt(contr_df$p.value, 4),
    sig          = add_stars(contr_df$p.value)
  )
}

# Helper: format emtrends summary with p-values for slope != 0
format_trend_summary <- function(em_obj, group_cols) {
  em_df <- as.data.frame(summary(em_obj))
  
  # Get test of slope != 0
  test_df <- tryCatch(as.data.frame(test(em_obj)), error = function(e) NULL)
  
  trend_col <- find_trend_col(em_df)
  if (is.null(trend_col)) {
    num_cols <- names(em_df)[sapply(em_df, is.numeric)]
    num_cols <- setdiff(num_cols, c("SE", "df"))
    trend_col <- num_cols[1]
  }
  
  ci <- find_ci(em_df)
  
  if (is.null(ci$lower) || is.null(ci$upper)) {
    ci_df <- tryCatch(as.data.frame(confint(em_obj)), error = function(e) NULL)
    if (!is.null(ci_df)) {
      ci2 <- find_ci(ci_df)
      if (!is.null(ci2$lower)) {
        em_df$lower.CI <- ci_df[[ci2$lower]]
        em_df$upper.CI <- ci_df[[ci2$upper]]
        ci <- list(lower = "lower.CI", upper = "upper.CI")
      }
    }
  }
  
  if (is.null(ci$lower) || is.null(ci$upper)) {
    em_df$lower.CI <- em_df[[trend_col]] - 1.96 * em_df$SE
    em_df$upper.CI <- em_df[[trend_col]] + 1.96 * em_df$SE
    ci <- list(lower = "lower.CI", upper = "upper.CI")
  }
  
  result <- tibble(
    slope      = fmt(em_df[[trend_col]], 4),
    SE         = fmt(em_df$SE, 4),
    `lower.CI` = fmt(em_df[[ci$lower]], 4),
    `upper.CI` = fmt(em_df[[ci$upper]], 4)
  )
  
  # Add p-value for slope != 0 if available
  if (!is.null(test_df) && "p.value" %in% names(test_df)) {
    result$p.value <- fmt(test_df$p.value, 4)
    result$sig     <- add_stars(test_df$p.value)
  } else if ("z.ratio" %in% names(em_df)) {
    pv <- 2 * pnorm(abs(em_df$z.ratio), lower.tail = FALSE)
    result$p.value <- fmt(pv, 4)
    result$sig     <- add_stars(pv)
  }
  
  for (gc in group_cols) {
    result[[gc]] <- em_df[[gc]]
  }
  result <- result %>% select(all_of(group_cols), everything())
  result
}


# =============================================================================
# 1. DIAGNOSIS MODEL (unchanged)
# =============================================================================
#
# islet_density ~ (Diagnosis + Region + Sex + Age_c)^2 +
#                 Diagnosis:Region:Age_c + (1|Donor)
#
# Sheets:
#   1. Model Summary
#   2. Diagnosis
#   3. Region
#   4. Sex
#   5. Diagnosis x Region
#   6. Diagnosis x Sex
#   7. Region x Sex
#   8. Age Slopes (overall, by Diagnosis, Region, Sex)
#   9. Diagnosis x Region x Age Slopes (3-way)
#  10. Random Effects & Sample Sizes
# =============================================================================

generate_density_diag_excel <- function(model, data, output_dir = "Density/Results") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING EXCEL: DENSITY — DIAGNOSIS MODEL\n")
  cat("=============================================================================\n\n")
  
  # -------------------------------------------------------------------------
  # SHEET 1: MODEL SUMMARY
  # -------------------------------------------------------------------------
  cat("Creating Sheet 1: Model Summary...\n")
  
  fe_df <- extract_fe(model)
  
  model_notes <- tibble(
    Note = c(
      "ISLET DENSITY ANALYSIS — DIAGNOSIS MODEL",
      "",
      "Generalized linear mixed model (glmmTMB)",
      "Response: islet_density (islets per mm²)",
      "Family: Gamma(log link)",
      "Formula: islet_density ~ (Diagnosis + Region + Sex + Age_c)^2 +",
      "         Diagnosis:Region:Age_c + (1|Donor)",
      "",
      "Reference levels: ND (Diagnosis), Head (Region), Female (Sex)",
      "Coefficients on log scale. Exponentiate for multiplicative effects.",
      "% change = (exp(estimate) - 1) * 100",
      "",
      "Multiplicity adjustments:",
      "  Tukey: all-pairwise comparisons within a factor (Region 3-way)",
      "  Holm: same contrast repeated across levels of a conditioning factor",
      "  None: single planned contrast (T1D vs ND main effect)",
      "",
      paste("N observations:", nrow(data)),
      paste("N donors:", n_distinct(data$Donor)),
      ""
    )
  )
  
  model_summary <- bind_rows(
    model_notes,
    fe_df %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 2: DIAGNOSIS
  # -------------------------------------------------------------------------
  cat("Creating Sheet 2: Diagnosis...\n")
  
  em_diag <- emmeans(model, ~ Diagnosis)
  em_diag_df <- extract_emm_response(as.data.frame(summary(em_diag)))
  
  diag_means <- tibble(
    Diagnosis = em_diag_df$Diagnosis,
    mean_density = fmt(em_diag_df$mean_density, 2),
    `lower.CI` = fmt(em_diag_df$lower.CI, 2),
    `upper.CI` = fmt(em_diag_df$upper.CI, 2)
  )
  
  contr_diag <- pairs(em_diag, adjust = "none", reverse = TRUE)
  contr_diag_df <- extract_contrast_results(as.data.frame(summary(contr_diag)))
  
  diag_contrast <- tibble(
    contrast = as.character(contr_diag_df$contrast),
    log_estimate = fmt(contr_diag_df$log_estimate, 4),
    SE = fmt(contr_diag_df$SE, 4),
    p.value = fmt(contr_diag_df$p.value, 4),
    pct_diff = fmt(contr_diag_df$pct_diff, 1),
    sig = add_stars(contr_diag_df$p.value)
  )
  
  diagnosis_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS EFFECT"),
    tibble(Note = "Main effect of Diagnosis on islet density, averaged over Region and Sex."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (islets/mm²)"),
    diag_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "CONTRAST (T1D vs ND, no adjustment)"),
    diag_contrast %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 3: REGION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: Region...\n")
  
  em_region <- emmeans(model, ~ Region)
  em_region_df <- extract_emm_response(as.data.frame(summary(em_region)))
  
  region_means <- tibble(
    Region = em_region_df$Region,
    mean_density = fmt(em_region_df$mean_density, 2),
    `lower.CI` = fmt(em_region_df$lower.CI, 2),
    `upper.CI` = fmt(em_region_df$upper.CI, 2)
  )
  
  contr_region <- pairs(em_region, adjust = "tukey")
  contr_region_df <- extract_contrast_results(as.data.frame(summary(contr_region)))
  
  region_contrasts <- tibble(
    contrast = as.character(contr_region_df$contrast),
    log_estimate = fmt(contr_region_df$log_estimate, 4),
    SE = fmt(contr_region_df$SE, 4),
    p.value = fmt(contr_region_df$p.value, 4),
    pct_diff = fmt(contr_region_df$pct_diff, 1),
    sig = add_stars(contr_region_df$p.value)
  )
  
  region_sheet <- bind_rows(
    tibble(Note = "REGION EFFECT"),
    tibble(Note = "Main effect of pancreatic Region on islet density, averaged over Diagnosis and Sex."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (islets/mm²)"),
    region_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "PAIRWISE CONTRASTS (Tukey-adjusted)"),
    region_contrasts %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 4: SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 4: Sex...\n")
  
  em_sex <- emmeans(model, ~ Sex)
  em_sex_df <- extract_emm_response(as.data.frame(summary(em_sex)))
  
  sex_means <- tibble(
    Sex = em_sex_df$Sex,
    mean_density = fmt(em_sex_df$mean_density, 2),
    `lower.CI` = fmt(em_sex_df$lower.CI, 2),
    `upper.CI` = fmt(em_sex_df$upper.CI, 2)
  )
  
  contr_sex <- pairs(em_sex, adjust = "none")
  contr_sex_df <- extract_contrast_results(as.data.frame(summary(contr_sex)))
  
  sex_contrast <- tibble(
    contrast = as.character(contr_sex_df$contrast),
    log_estimate = fmt(contr_sex_df$log_estimate, 4),
    SE = fmt(contr_sex_df$SE, 4),
    p.value = fmt(contr_sex_df$p.value, 4),
    pct_diff = fmt(contr_sex_df$pct_diff, 1),
    sig = add_stars(contr_sex_df$p.value)
  )
  
  sex_sheet <- bind_rows(
    tibble(Note = "SEX EFFECT"),
    tibble(Note = "Main effect of Sex on islet density, averaged over Diagnosis and Region."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (islets/mm²)"),
    sex_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "CONTRAST (Female vs Male, no adjustment)"),
    sex_contrast %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 5: DIAGNOSIS x REGION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 5: Diagnosis x Region...\n")
  
  em_dr <- emmeans(model, ~ Diagnosis * Region)
  em_dr_df <- extract_emm_response(as.data.frame(summary(em_dr)))
  
  dr_means <- tibble(
    Diagnosis = em_dr_df$Diagnosis,
    Region = em_dr_df$Region,
    mean_density = fmt(em_dr_df$mean_density, 2),
    `lower.CI` = fmt(em_dr_df$lower.CI, 2),
    `upper.CI` = fmt(em_dr_df$upper.CI, 2)
  )
  
  em_diag_by_region <- emmeans(model, ~ Diagnosis | Region)
  contr_diag_by_region <- pairs(em_diag_by_region, adjust = "holm", reverse = TRUE)
  contr_dbr_df <- extract_contrast_results(as.data.frame(summary(contr_diag_by_region)))
  
  diag_within_region <- tibble(
    Region = contr_dbr_df$Region,
    contrast = as.character(contr_dbr_df$contrast),
    log_estimate = fmt(contr_dbr_df$log_estimate, 4),
    SE = fmt(contr_dbr_df$SE, 4),
    p.value = fmt(contr_dbr_df$p.value, 4),
    pct_diff = fmt(contr_dbr_df$pct_diff, 1),
    sig = add_stars(contr_dbr_df$p.value)
  )
  
  em_region_by_diag <- emmeans(model, ~ Region | Diagnosis)
  contr_region_by_diag <- pairs(em_region_by_diag, adjust = "tukey")
  contr_rbd_df <- extract_contrast_results(as.data.frame(summary(contr_region_by_diag)))
  
  region_within_diag <- tibble(
    Diagnosis = contr_rbd_df$Diagnosis,
    contrast = as.character(contr_rbd_df$contrast),
    log_estimate = fmt(contr_rbd_df$log_estimate, 4),
    SE = fmt(contr_rbd_df$SE, 4),
    p.value = fmt(contr_rbd_df$p.value, 4),
    pct_diff = fmt(contr_rbd_df$pct_diff, 1),
    sig = add_stars(contr_rbd_df$p.value)
  )
  
  diag_region_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x REGION INTERACTION"),
    tibble(Note = "Marginal means and contrasts averaged over Sex."),
    tibble(Note = NA),
    tibble(Note = "CELL MEANS (islets/mm²)"),
    dr_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT WITHIN EACH REGION (T1D vs ND, Holm-adjusted)"),
    diag_within_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS WITHIN EACH DIAGNOSIS (Tukey-adjusted)"),
    region_within_diag %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 6: DIAGNOSIS x SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 6: Diagnosis x Sex...\n")
  
  em_ds <- emmeans(model, ~ Diagnosis * Sex)
  em_ds_df <- extract_emm_response(as.data.frame(summary(em_ds)))
  
  ds_means <- tibble(
    Diagnosis = em_ds_df$Diagnosis,
    Sex = em_ds_df$Sex,
    mean_density = fmt(em_ds_df$mean_density, 2),
    `lower.CI` = fmt(em_ds_df$lower.CI, 2),
    `upper.CI` = fmt(em_ds_df$upper.CI, 2)
  )
  
  em_diag_by_sex <- emmeans(model, ~ Diagnosis | Sex)
  contr_diag_by_sex <- pairs(em_diag_by_sex, adjust = "holm", reverse = TRUE)
  contr_dbs_df <- extract_contrast_results(as.data.frame(summary(contr_diag_by_sex)))
  
  diag_within_sex <- tibble(
    Sex = contr_dbs_df$Sex,
    contrast = as.character(contr_dbs_df$contrast),
    log_estimate = fmt(contr_dbs_df$log_estimate, 4),
    SE = fmt(contr_dbs_df$SE, 4),
    p.value = fmt(contr_dbs_df$p.value, 4),
    pct_diff = fmt(contr_dbs_df$pct_diff, 1),
    sig = add_stars(contr_dbs_df$p.value)
  )
  
  em_sex_by_diag <- emmeans(model, ~ Sex | Diagnosis)
  contr_sex_by_diag <- pairs(em_sex_by_diag, adjust = "holm")
  contr_sbd_df <- extract_contrast_results(as.data.frame(summary(contr_sex_by_diag)))
  
  sex_within_diag <- tibble(
    Diagnosis = contr_sbd_df$Diagnosis,
    contrast = as.character(contr_sbd_df$contrast),
    log_estimate = fmt(contr_sbd_df$log_estimate, 4),
    SE = fmt(contr_sbd_df$SE, 4),
    p.value = fmt(contr_sbd_df$p.value, 4),
    pct_diff = fmt(contr_sbd_df$pct_diff, 1),
    sig = add_stars(contr_sbd_df$p.value)
  )
  
  diag_sex_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x SEX INTERACTION"),
    tibble(Note = "Marginal means and contrasts averaged over Region."),
    tibble(Note = NA),
    tibble(Note = "CELL MEANS (islets/mm²)"),
    ds_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT WITHIN EACH SEX (T1D vs ND, Holm-adjusted)"),
    diag_within_sex %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SEX EFFECT WITHIN EACH DIAGNOSIS (Female vs Male, Holm-adjusted)"),
    sex_within_diag %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 7: REGION x SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 7: Region x Sex...\n")
  
  em_rs <- emmeans(model, ~ Region * Sex)
  em_rs_df <- extract_emm_response(as.data.frame(summary(em_rs)))
  
  rs_means <- tibble(
    Region = em_rs_df$Region,
    Sex = em_rs_df$Sex,
    mean_density = fmt(em_rs_df$mean_density, 2),
    `lower.CI` = fmt(em_rs_df$lower.CI, 2),
    `upper.CI` = fmt(em_rs_df$upper.CI, 2)
  )
  
  em_sex_by_region <- emmeans(model, ~ Sex | Region)
  contr_sex_by_region <- pairs(em_sex_by_region, adjust = "holm")
  contr_sbr_df <- extract_contrast_results(as.data.frame(summary(contr_sex_by_region)))
  
  sex_within_region <- tibble(
    Region = contr_sbr_df$Region,
    contrast = as.character(contr_sbr_df$contrast),
    log_estimate = fmt(contr_sbr_df$log_estimate, 4),
    SE = fmt(contr_sbr_df$SE, 4),
    p.value = fmt(contr_sbr_df$p.value, 4),
    pct_diff = fmt(contr_sbr_df$pct_diff, 1),
    sig = add_stars(contr_sbr_df$p.value)
  )
  
  em_region_by_sex <- emmeans(model, ~ Region | Sex)
  contr_region_by_sex <- pairs(em_region_by_sex, adjust = "tukey")
  contr_rbs_df <- extract_contrast_results(as.data.frame(summary(contr_region_by_sex)))
  
  region_within_sex <- tibble(
    Sex = contr_rbs_df$Sex,
    contrast = as.character(contr_rbs_df$contrast),
    log_estimate = fmt(contr_rbs_df$log_estimate, 4),
    SE = fmt(contr_rbs_df$SE, 4),
    p.value = fmt(contr_rbs_df$p.value, 4),
    pct_diff = fmt(contr_rbs_df$pct_diff, 1),
    sig = add_stars(contr_rbs_df$p.value)
  )
  
  region_sex_sheet <- bind_rows(
    tibble(Note = "REGION x SEX INTERACTION"),
    tibble(Note = "Marginal means and contrasts averaged over Diagnosis."),
    tibble(Note = NA),
    tibble(Note = "CELL MEANS (islets/mm²)"),
    rs_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SEX EFFECT WITHIN EACH REGION (Female vs Male, Holm-adjusted)"),
    sex_within_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS WITHIN EACH SEX (Tukey-adjusted)"),
    region_within_sex %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 8: AGE SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 8: Age Slopes...\n")
  
  em_age_overall <- emtrends(model, ~ 1, var = "Age_c")
  age_overall_df <- format_trend_summary(em_age_overall, character(0))
  
  em_age_diag <- emtrends(model, ~ Diagnosis, var = "Age_c")
  age_diag_df <- format_trend_summary(em_age_diag, "Diagnosis")
  contr_age_diag <- pairs(em_age_diag, adjust = "none")
  age_diag_contr <- format_trend_contrasts(contr_age_diag)
  
  em_age_region <- emtrends(model, ~ Region, var = "Age_c")
  age_region_df <- format_trend_summary(em_age_region, "Region")
  contr_age_region <- pairs(em_age_region, adjust = "tukey")
  age_region_contr <- format_trend_contrasts(contr_age_region)
  
  em_age_sex <- emtrends(model, ~ Sex, var = "Age_c")
  age_sex_df <- format_trend_summary(em_age_sex, "Sex")
  contr_age_sex <- pairs(em_age_sex, adjust = "none")
  age_sex_contr <- format_trend_contrasts(contr_age_sex)
  
  age_sheet <- bind_rows(
    tibble(Note = "AGE SLOPES"),
    tibble(Note = "Slopes of Age_c on log(islet_density) from emtrends."),
    tibble(Note = "On the log scale: slope = change in log(density) per 1-year increase in Age (centered)."),
    tibble(Note = "Positive slope => density increases with age; negative => decreases."),
    tibble(Note = "Approximate % change per year = (exp(slope) - 1) * 100."),
    tibble(Note = NA),
    tibble(Note = "OVERALL AGE SLOPE (marginalized over Diagnosis, Region, Sex)"),
    age_overall_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "AGE SLOPES BY DIAGNOSIS"),
    age_diag_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Diagnosis contrast on age slopes (no adjustment):"),
    age_diag_contr %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "AGE SLOPES BY REGION"),
    age_region_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Region pairwise contrasts on age slopes (Tukey-adjusted):"),
    age_region_contr %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "AGE SLOPES BY SEX"),
    age_sex_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Sex contrast on age slopes (no adjustment):"),
    age_sex_contr %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 9: DIAGNOSIS x REGION x AGE (3-way)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 9: Diagnosis x Region x Age Slopes...\n")
  
  em_age_dr <- emtrends(model, ~ Diagnosis * Region, var = "Age_c")
  age_dr_df <- format_trend_summary(em_age_dr, c("Diagnosis", "Region"))
  
  em_age_diag_by_region <- emtrends(model, ~ Diagnosis | Region, var = "Age_c")
  contr_age_diag_by_region <- pairs(em_age_diag_by_region, adjust = "holm", reverse = TRUE)
  age_diag_within_region <- format_trend_contrasts(contr_age_diag_by_region)
  age_diag_within_region$Region <- as.data.frame(summary(contr_age_diag_by_region))$Region
  age_diag_within_region <- age_diag_within_region %>% select(Region, everything())
  
  em_age_region_by_diag <- emtrends(model, ~ Region | Diagnosis, var = "Age_c")
  contr_age_region_by_diag <- pairs(em_age_region_by_diag, adjust = "tukey")
  age_region_within_diag <- format_trend_contrasts(contr_age_region_by_diag)
  age_region_within_diag$Diagnosis <- as.data.frame(summary(contr_age_region_by_diag))$Diagnosis
  age_region_within_diag <- age_region_within_diag %>% select(Diagnosis, everything())
  
  diag_region_age_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x REGION x AGE SLOPES (3-way interaction)"),
    tibble(Note = "Age slopes for each Diagnosis x Region cell, marginalized over Sex."),
    tibble(Note = NA),
    tibble(Note = "AGE SLOPES PER DIAGNOSIS x REGION CELL"),
    age_dr_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT ON AGE SLOPES WITHIN EACH REGION (T1D vs ND, Holm-adjusted)"),
    tibble(Note = "Tests whether the age-density relationship differs between ND and T1D in each Region."),
    age_diag_within_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS ON AGE SLOPES WITHIN EACH DIAGNOSIS (Tukey-adjusted)"),
    tibble(Note = "Tests whether the age-density relationship differs across Regions within each Diagnosis."),
    age_region_within_diag %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 10: RANDOM EFFECTS & SAMPLE SIZES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 10: Random Effects & Sample Sizes...\n")
  
  re_df <- extract_re(model, family_name = "gamma")
  
  sample_sizes <- tibble(
    n_obs    = nrow(data),
    n_donors = n_distinct(data$Donor),
    n_images = if ("ImageID" %in% names(data)) n_distinct(data$ImageID) else NA
  )
  
  re_sheet <- bind_rows(
    tibble(Note = "RANDOM EFFECTS VARIANCE COMPONENTS"),
    tibble(Note = NA),
    re_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SAMPLE SIZES"),
    sample_sizes %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # ASSEMBLE AND SAVE
  # =========================================================================
  cat("\nAssembling workbook...\n")
  
  sheets <- list(
    "Model Summary"       = as.data.frame(model_summary),
    "Diagnosis"           = as.data.frame(diagnosis_sheet),
    "Region"              = as.data.frame(region_sheet),
    "Sex"                 = as.data.frame(sex_sheet),
    "Diagnosis x Region"  = as.data.frame(diag_region_sheet),
    "Diagnosis x Sex"     = as.data.frame(diag_sex_sheet),
    "Region x Sex"        = as.data.frame(region_sex_sheet),
    "Age Slopes"          = as.data.frame(age_sheet),
    "Diag x Region x Age" = as.data.frame(diag_region_age_sheet),
    "Random Effects"      = as.data.frame(re_sheet)
  )
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, "Islet_Density_Diagnosis_Analysis.xlsx")
  write_xlsx(sheets, output_file)
  cat("Saved to:", output_file, "\n")
  
  return(invisible(sheets))
}


# =============================================================================
# 2. T1D MODEL — AGE AT ONSET PARAMETERIZATION
# =============================================================================
#
# islet_density ~ (Disease.Duration_c + Region + Sex + age_at_onset_c)^2 +
#                 Disease.Duration_c:Region:age_at_onset_c + (1|Donor)
#
# This single model replaces the separate DD and AO models.  Because
# Age = age_at_onset + Disease.Duration, the DD + Age model and the
# DD + AO model span the same column space and produce identical fits.
# The AO parameterization is preferred because both coefficients have
# clean clinical interpretations:
#   - DD effect: longer disease holding age at onset constant
#   - AO effect: earlier vs later onset holding duration constant
#
# Two-way interactions from (...)^2:
#   DD:Region, DD:Sex, DD:AO, Region:Sex, Region:AO, Sex:AO
# Three-way interaction:
#   DD:Region:AO
#
# Sheets:
#   1. Model Summary (fixed effects)
#   2. Region
#   3. Sex
#   4. Region x Sex
#   5. DD Slopes (overall, by Region, by Sex + contrasts)
#   6. AO Slopes (overall, by Region, by Sex + contrasts)
#   7. DD x AO Interaction (DD slopes at representative AO values & vice versa)
#   8. DD x Region x AO (3-way decomposition)
#   9. Random Effects & Sample Sizes
# =============================================================================

generate_density_t1d_excel <- function(model, data, output_dir = "Density/Results") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING EXCEL: DENSITY — T1D MODEL (Age at Onset parameterization)\n")
  cat("=============================================================================\n\n")
  
  # --- Compute representative values for continuous × continuous probing ---
  ao_sd <- sd(data$age_at_onset_c, na.rm = TRUE)
  dd_sd <- sd(data$Disease.Duration_c, na.rm = TRUE)
  ao_at <- c(-ao_sd, 0, ao_sd)
  dd_at <- c(-dd_sd, 0, dd_sd)
  
  # Uncentered means for readable labels
  ao_mean_raw <- mean(data$age_at_onset, na.rm = TRUE)
  dd_mean_raw <- mean(data$Disease.Duration, na.rm = TRUE)
  
  ao_labels <- c(
    paste0("Low (-1 SD, ~", round(ao_mean_raw - ao_sd, 1), " yrs)"),
    paste0("Mean (~", round(ao_mean_raw, 1), " yrs)"),
    paste0("High (+1 SD, ~", round(ao_mean_raw + ao_sd, 1), " yrs)")
  )
  dd_labels <- c(
    paste0("Short (-1 SD, ~", round(dd_mean_raw - dd_sd, 1), " yrs)"),
    paste0("Mean (~", round(dd_mean_raw, 1), " yrs)"),
    paste0("Long (+1 SD, ~", round(dd_mean_raw + dd_sd, 1), " yrs)")
  )
  
  # -------------------------------------------------------------------------
  # SHEET 1: MODEL SUMMARY
  # -------------------------------------------------------------------------
  cat("Creating Sheet 1: Model Summary...\n")
  
  fe_df <- extract_fe(model)
  
  model_notes <- tibble(
    Note = c(
      "ISLET DENSITY ANALYSIS — T1D MODEL (Age at Onset parameterization)",
      "",
      "Generalized linear mixed model (glmmTMB)",
      "Response: islet_density (islets per mm²)",
      "Family: Gamma(log link)",
      "Formula: islet_density ~ (Disease.Duration_c + Region + Sex + age_at_onset_c)^2 +",
      "         Disease.Duration_c:Region:age_at_onset_c + (1|Donor)",
      "",
      "Two-way interactions: DD:Region, DD:Sex, DD:AO, Region:Sex, Region:AO, Sex:AO",
      "Three-way interaction: DD:Region:AO",
      "",
      "Note: Age = age_at_onset + Disease.Duration, so this model is algebraically",
      "equivalent to the DD + Age_c parameterization.  The AO version is preferred",
      "because both slopes have direct clinical meaning:",
      "  DD  = effect of longer duration, holding age at onset constant",
      "  AO  = effect of later onset, holding duration constant",
      "",
      "Reference levels: Head (Region), Female (Sex)",
      "Coefficients on log scale. Exponentiate for multiplicative effects.",
      "% change = (exp(estimate) - 1) * 100",
      "",
      "Multiplicity adjustments:",
      "  Tukey: all-pairwise comparisons within Region (3 levels)",
      "  Holm: same contrast repeated across levels of a conditioning factor",
      "  None: single planned contrast (2-level factor or overall slope)",
      "",
      paste("N observations:", nrow(data)),
      paste("N donors:", n_distinct(data$Donor)),
      paste("Age at onset — mean:", round(ao_mean_raw, 1),
            "yrs, SD:", round(ao_sd, 1), "yrs"),
      paste("Disease duration — mean:", round(dd_mean_raw, 1),
            "yrs, SD:", round(dd_sd, 1), "yrs"),
      ""
    )
  )
  
  model_summary <- bind_rows(
    model_notes,
    fe_df %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 2: REGION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 2: Region...\n")
  
  em_region <- emmeans(model, ~ Region)
  em_region_df <- extract_emm_response(as.data.frame(summary(em_region)))
  
  region_means <- tibble(
    Region = em_region_df$Region,
    mean_density = fmt(em_region_df$mean_density, 2),
    `lower.CI` = fmt(em_region_df$lower.CI, 2),
    `upper.CI` = fmt(em_region_df$upper.CI, 2)
  )
  
  contr_region <- pairs(em_region, adjust = "tukey")
  contr_region_df <- extract_contrast_results(as.data.frame(summary(contr_region)))
  
  region_contrasts <- tibble(
    contrast = as.character(contr_region_df$contrast),
    log_estimate = fmt(contr_region_df$log_estimate, 4),
    SE = fmt(contr_region_df$SE, 4),
    p.value = fmt(contr_region_df$p.value, 4),
    pct_diff = fmt(contr_region_df$pct_diff, 1),
    sig = add_stars(contr_region_df$p.value)
  )
  
  region_sheet <- bind_rows(
    tibble(Note = "REGION EFFECT (T1D only)"),
    tibble(Note = "Main effect of pancreatic Region on islet density, averaged over Sex."),
    tibble(Note = "Evaluated at centered means of DD and AO (both = 0)."),
    tibble(Note = "Region differences may vary at other DD/AO values — see interaction sheets."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (islets/mm²)"),
    region_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "PAIRWISE CONTRASTS (Tukey-adjusted)"),
    region_contrasts %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 3: SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: Sex...\n")
  
  em_sex <- emmeans(model, ~ Sex)
  em_sex_df <- extract_emm_response(as.data.frame(summary(em_sex)))
  
  sex_means <- tibble(
    Sex = em_sex_df$Sex,
    mean_density = fmt(em_sex_df$mean_density, 2),
    `lower.CI` = fmt(em_sex_df$lower.CI, 2),
    `upper.CI` = fmt(em_sex_df$upper.CI, 2)
  )
  
  contr_sex <- pairs(em_sex, adjust = "none")
  contr_sex_df <- extract_contrast_results(as.data.frame(summary(contr_sex)))
  
  sex_contrast <- tibble(
    contrast = as.character(contr_sex_df$contrast),
    log_estimate = fmt(contr_sex_df$log_estimate, 4),
    SE = fmt(contr_sex_df$SE, 4),
    p.value = fmt(contr_sex_df$p.value, 4),
    pct_diff = fmt(contr_sex_df$pct_diff, 1),
    sig = add_stars(contr_sex_df$p.value)
  )
  
  sex_sheet <- bind_rows(
    tibble(Note = "SEX EFFECT (T1D only)"),
    tibble(Note = "Main effect of Sex on islet density, averaged over Region."),
    tibble(Note = "Evaluated at centered means of DD and AO."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (islets/mm²)"),
    sex_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "CONTRAST (Female vs Male, no adjustment)"),
    sex_contrast %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 4: REGION x SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 4: Region x Sex...\n")
  
  em_rs <- emmeans(model, ~ Region * Sex)
  em_rs_df <- extract_emm_response(as.data.frame(summary(em_rs)))
  
  rs_means <- tibble(
    Region = em_rs_df$Region,
    Sex = em_rs_df$Sex,
    mean_density = fmt(em_rs_df$mean_density, 2),
    `lower.CI` = fmt(em_rs_df$lower.CI, 2),
    `upper.CI` = fmt(em_rs_df$upper.CI, 2)
  )
  
  em_sex_by_region <- emmeans(model, ~ Sex | Region)
  contr_sex_by_region <- pairs(em_sex_by_region, adjust = "holm")
  contr_sbr_df <- extract_contrast_results(as.data.frame(summary(contr_sex_by_region)))
  
  sex_within_region <- tibble(
    Region = contr_sbr_df$Region,
    contrast = as.character(contr_sbr_df$contrast),
    log_estimate = fmt(contr_sbr_df$log_estimate, 4),
    SE = fmt(contr_sbr_df$SE, 4),
    p.value = fmt(contr_sbr_df$p.value, 4),
    pct_diff = fmt(contr_sbr_df$pct_diff, 1),
    sig = add_stars(contr_sbr_df$p.value)
  )
  
  em_region_by_sex <- emmeans(model, ~ Region | Sex)
  contr_region_by_sex <- pairs(em_region_by_sex, adjust = "tukey")
  contr_rbs_df <- extract_contrast_results(as.data.frame(summary(contr_region_by_sex)))
  
  region_within_sex <- tibble(
    Sex = contr_rbs_df$Sex,
    contrast = as.character(contr_rbs_df$contrast),
    log_estimate = fmt(contr_rbs_df$log_estimate, 4),
    SE = fmt(contr_rbs_df$SE, 4),
    p.value = fmt(contr_rbs_df$p.value, 4),
    pct_diff = fmt(contr_rbs_df$pct_diff, 1),
    sig = add_stars(contr_rbs_df$p.value)
  )
  
  region_sex_sheet <- bind_rows(
    tibble(Note = "REGION x SEX INTERACTION (T1D only)"),
    tibble(Note = "Marginal means and contrasts evaluated at centered means of DD and AO."),
    tibble(Note = NA),
    tibble(Note = "CELL MEANS (islets/mm²)"),
    rs_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SEX EFFECT WITHIN EACH REGION (Female vs Male, Holm-adjusted)"),
    sex_within_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS WITHIN EACH SEX (Tukey-adjusted)"),
    region_within_sex %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 5: DISEASE DURATION SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 5: DD Slopes...\n")
  
  # Overall DD slope (marginalized over Region, Sex; at AO = 0)
  em_dd_overall <- emtrends(model, ~ 1, var = "Disease.Duration_c")
  dd_overall_df <- format_trend_summary(em_dd_overall, character(0))
  
  # DD slopes by Region (from DD:Region interaction; marginalized over Sex, at AO = 0)
  em_dd_by_region <- emtrends(model, ~ Region, var = "Disease.Duration_c")
  dd_by_region_df <- format_trend_summary(em_dd_by_region, "Region")
  contr_dd_region <- pairs(em_dd_by_region, adjust = "tukey")
  dd_region_contr <- format_trend_contrasts(contr_dd_region)
  
  # DD slopes by Sex (from DD:Sex interaction; marginalized over Region, at AO = 0)
  em_dd_by_sex <- emtrends(model, ~ Sex, var = "Disease.Duration_c")
  dd_by_sex_df <- format_trend_summary(em_dd_by_sex, "Sex")
  contr_dd_sex <- pairs(em_dd_by_sex, adjust = "none")
  dd_sex_contr <- format_trend_contrasts(contr_dd_sex)
  
  dd_slopes_sheet <- bind_rows(
    tibble(Note = "DISEASE DURATION SLOPES (T1D only)"),
    tibble(Note = "Slopes of Disease.Duration_c on log(islet_density) from emtrends."),
    tibble(Note = "Slope = change in log(density) per 1-year increase in DD (centered)."),
    tibble(Note = "Approximate % change per year = (exp(slope) - 1) * 100."),
    tibble(Note = ""),
    tibble(Note = "The model includes DD:Region and DD:Sex interactions, so the DD"),
    tibble(Note = "slope varies by both Region and Sex.  All slopes evaluated at"),
    tibble(Note = "AO = 0 (sample mean).  See DD x AO sheet for how slopes change with AO."),
    tibble(Note = NA),
    tibble(Note = "OVERALL DD SLOPE (marginalized over Region and Sex)"),
    dd_overall_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DD SLOPES BY REGION (marginalized over Sex, Tukey-adjusted contrasts)"),
    tibble(Note = "How the DD-density relationship differs across pancreatic regions."),
    dd_by_region_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Region pairwise contrasts on DD slopes:"),
    dd_region_contr %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DD SLOPES BY SEX (marginalized over Region)"),
    dd_by_sex_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Sex contrast on DD slopes (no adjustment):"),
    dd_sex_contr %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 6: AGE AT ONSET SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 6: AO Slopes...\n")
  
  # Overall AO slope (marginalized over Region, Sex; at DD = 0)
  em_ao_overall <- emtrends(model, ~ 1, var = "age_at_onset_c")
  ao_overall_df <- format_trend_summary(em_ao_overall, character(0))
  
  # AO slopes by Region (from AO:Region interaction)
  em_ao_by_region <- emtrends(model, ~ Region, var = "age_at_onset_c")
  ao_by_region_df <- format_trend_summary(em_ao_by_region, "Region")
  contr_ao_region <- pairs(em_ao_by_region, adjust = "tukey")
  ao_region_contr <- format_trend_contrasts(contr_ao_region)
  
  # AO slopes by Sex (from AO:Sex interaction)
  em_ao_by_sex <- emtrends(model, ~ Sex, var = "age_at_onset_c")
  ao_by_sex_df <- format_trend_summary(em_ao_by_sex, "Sex")
  contr_ao_sex <- pairs(em_ao_by_sex, adjust = "none")
  ao_sex_contr <- format_trend_contrasts(contr_ao_sex)
  
  ao_slopes_sheet <- bind_rows(
    tibble(Note = "AGE AT ONSET SLOPES (T1D only)"),
    tibble(Note = "Slopes of age_at_onset_c on log(islet_density) from emtrends."),
    tibble(Note = "Slope = change in log(density) per 1-year increase in AO (centered)."),
    tibble(Note = "Positive slope => later onset = higher density."),
    tibble(Note = "Approximate % change per year = (exp(slope) - 1) * 100."),
    tibble(Note = ""),
    tibble(Note = "The model includes AO:Region and AO:Sex interactions, so the AO"),
    tibble(Note = "slope varies by both Region and Sex.  All slopes evaluated at"),
    tibble(Note = "DD = 0 (sample mean).  See DD x AO sheet for how slopes change with DD."),
    tibble(Note = NA),
    tibble(Note = "OVERALL AO SLOPE (marginalized over Region and Sex)"),
    ao_overall_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "AO SLOPES BY REGION (marginalized over Sex, Tukey-adjusted contrasts)"),
    tibble(Note = "How the AO-density relationship differs across pancreatic regions."),
    ao_by_region_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Region pairwise contrasts on AO slopes:"),
    ao_region_contr %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "AO SLOPES BY SEX (marginalized over Region)"),
    ao_by_sex_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Sex contrast on AO slopes (no adjustment):"),
    ao_sex_contr %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 7: DD x AO INTERACTION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 7: DD x AO Interaction...\n")
  
  # DD slopes at low / mean / high AO
  em_dd_at_ao <- emtrends(model, ~ age_at_onset_c, var = "Disease.Duration_c",
                           at = list(age_at_onset_c = ao_at))
  dd_at_ao_df <- format_trend_summary(em_dd_at_ao, "age_at_onset_c")
  dd_at_ao_df$AO_level <- ao_labels
  dd_at_ao_df <- dd_at_ao_df %>%
    select(AO_level, age_at_onset_c, everything())
  
  contr_dd_across_ao <- pairs(em_dd_at_ao, adjust = "tukey")
  dd_across_ao_contr <- format_trend_contrasts(contr_dd_across_ao)
  
  # AO slopes at low / mean / high DD
  em_ao_at_dd <- emtrends(model, ~ Disease.Duration_c, var = "age_at_onset_c",
                           at = list(Disease.Duration_c = dd_at))
  ao_at_dd_df <- format_trend_summary(em_ao_at_dd, "Disease.Duration_c")
  ao_at_dd_df$DD_level <- dd_labels
  ao_at_dd_df <- ao_at_dd_df %>%
    select(DD_level, Disease.Duration_c, everything())
  
  contr_ao_across_dd <- pairs(em_ao_at_dd, adjust = "tukey")
  ao_across_dd_contr <- format_trend_contrasts(contr_ao_across_dd)
  
  dd_ao_sheet <- bind_rows(
    tibble(Note = "DD x AO INTERACTION (T1D only)"),
    tibble(Note = "The DD:age_at_onset_c interaction means the DD slope changes with AO and vice versa."),
    tibble(Note = "Below we probe the interaction by evaluating slopes at -1 SD, mean, and +1 SD"),
    tibble(Note = "of the other variable (marginalized over Region and Sex)."),
    tibble(Note = NA),
    tibble(Note = "DD SLOPES AT REPRESENTATIVE AGE-AT-ONSET VALUES"),
    tibble(Note = "How does the effect of disease duration on density change with age at onset?"),
    dd_at_ao_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Pairwise contrasts on DD slopes across AO levels (Tukey-adjusted):"),
    dd_across_ao_contr %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "AO SLOPES AT REPRESENTATIVE DISEASE DURATION VALUES"),
    tibble(Note = "How does the effect of age at onset on density change with disease duration?"),
    ao_at_dd_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "Pairwise contrasts on AO slopes across DD levels (Tukey-adjusted):"),
    ao_across_dd_contr %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 8: DD x REGION x AO (3-way)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 8: DD x Region x AO (3-way)...\n")
  
  # DD slopes for each Region at low / mean / high AO
  em_dd_region_ao <- emtrends(model, ~ Region | age_at_onset_c,
                               var = "Disease.Duration_c",
                               at = list(age_at_onset_c = ao_at))
  dd_region_ao_df <- as.data.frame(summary(em_dd_region_ao))
  
  # Find the trend column
  trend_col <- find_trend_col(dd_region_ao_df)
  ci <- find_ci(dd_region_ao_df)
  if (is.null(ci$lower) || is.null(ci$upper)) {
    dd_region_ao_df$lower.CI <- dd_region_ao_df[[trend_col]] - 1.96 * dd_region_ao_df$SE
    dd_region_ao_df$upper.CI <- dd_region_ao_df[[trend_col]] + 1.96 * dd_region_ao_df$SE
    ci <- list(lower = "lower.CI", upper = "upper.CI")
  }
  
  # Get test p-values
  test_rao <- tryCatch(as.data.frame(test(em_dd_region_ao)), error = function(e) NULL)
  
  threeway_slopes <- tibble(
    Region = dd_region_ao_df$Region,
    age_at_onset_c = round(dd_region_ao_df$age_at_onset_c, 2),
    AO_level = ao_labels[match(round(dd_region_ao_df$age_at_onset_c, 2),
                                round(ao_at, 2))],
    DD_slope = fmt(dd_region_ao_df[[trend_col]], 4),
    SE       = fmt(dd_region_ao_df$SE, 4),
    lower.CI = fmt(dd_region_ao_df[[ci$lower]], 4),
    upper.CI = fmt(dd_region_ao_df[[ci$upper]], 4)
  )
  if (!is.null(test_rao) && "p.value" %in% names(test_rao)) {
    threeway_slopes$p.value <- fmt(test_rao$p.value, 4)
    threeway_slopes$sig     <- add_stars(test_rao$p.value)
  }
  
  # Region contrasts on DD slopes at each AO level (Tukey)
  contr_region_at_ao <- pairs(em_dd_region_ao, adjust = "tukey")
  contr_rao_df <- as.data.frame(summary(contr_region_at_ao))
  ci_rao <- find_ci(contr_rao_df)
  if (is.null(ci_rao$lower) || is.null(ci_rao$upper)) {
    contr_rao_df$lower.CI <- contr_rao_df$estimate - 1.96 * contr_rao_df$SE
    contr_rao_df$upper.CI <- contr_rao_df$estimate + 1.96 * contr_rao_df$SE
    ci_rao <- list(lower = "lower.CI", upper = "upper.CI")
  }
  
  region_contr_at_ao <- tibble(
    age_at_onset_c = round(contr_rao_df$age_at_onset_c, 2),
    AO_level = ao_labels[match(round(contr_rao_df$age_at_onset_c, 2),
                                round(ao_at, 2))],
    contrast = as.character(contr_rao_df$contrast),
    estimate = fmt(contr_rao_df$estimate, 4),
    SE       = fmt(contr_rao_df$SE, 4),
    lower.CI = fmt(contr_rao_df[[ci_rao$lower]], 4),
    upper.CI = fmt(contr_rao_df[[ci_rao$upper]], 4),
    p.value  = fmt(contr_rao_df$p.value, 4),
    sig      = add_stars(contr_rao_df$p.value)
  )
  
  threeway_sheet <- bind_rows(
    tibble(Note = "DD x REGION x AO — THREE-WAY INTERACTION (T1D only)"),
    tibble(Note = "The Disease.Duration_c:Region:age_at_onset_c three-way means"),
    tibble(Note = "the DD slope differs by Region AND those Region differences change with AO."),
    tibble(Note = ""),
    tibble(Note = "We probe this by showing DD slopes for each Region at -1 SD, mean,"),
    tibble(Note = "and +1 SD of age at onset, marginalized over Sex."),
    tibble(Note = NA),
    tibble(Note = "DD SLOPES PER REGION x AO LEVEL"),
    tibble(Note = "Each row = DD slope for one Region at one AO value."),
    tibble(Note = "p-value tests slope != 0."),
    threeway_slopes %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS ON DD SLOPES AT EACH AO LEVEL (Tukey-adjusted)"),
    tibble(Note = "Do Region differences in the DD-density relationship depend on age at onset?"),
    region_contr_at_ao %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # -------------------------------------------------------------------------
  # SHEET 9: RANDOM EFFECTS & SAMPLE SIZES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 9: Random Effects & Sample Sizes...\n")
  
  re_df <- extract_re(model, family_name = "gamma")
  
  sample_sizes <- tibble(
    n_obs    = nrow(data),
    n_donors = n_distinct(data$Donor)
  )
  
  re_sheet <- bind_rows(
    tibble(Note = "RANDOM EFFECTS VARIANCE COMPONENTS"),
    tibble(Note = "Donor variance is on the log scale (log-link model)."),
    tibble(Note = "Gamma dispersion = 1/shape; controls the coefficient of variation"),
    tibble(Note = "of the response (CV ≈ sqrt(dispersion))."),
    tibble(Note = NA),
    re_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SAMPLE SIZES"),
    sample_sizes %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # ASSEMBLE AND SAVE
  # =========================================================================
  cat("\nAssembling workbook...\n")
  
  sheets <- list(
    "Model Summary"     = as.data.frame(model_summary),
    "Region"            = as.data.frame(region_sheet),
    "Sex"               = as.data.frame(sex_sheet),
    "Region x Sex"      = as.data.frame(region_sex_sheet),
    "DD Slopes"         = as.data.frame(dd_slopes_sheet),
    "AO Slopes"         = as.data.frame(ao_slopes_sheet),
    "DD x AO"           = as.data.frame(dd_ao_sheet),
    "DD x Region x AO"  = as.data.frame(threeway_sheet),
    "Random Effects"    = as.data.frame(re_sheet)
  )
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, "Islet_Density_T1D_Analysis.xlsx")
  write_xlsx(sheets, output_file)
  cat("Saved to:", output_file, "\n")
  
  return(invisible(sheets))
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

generate_all_density_excel <- function(output_dir = "Density/Results") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING ALL DENSITY EXCEL RESULTS\n")
  cat("=============================================================================\n\n")
  
  # Load data
  cat("Loading data...\n")
  density_dataset <- readRDS("Data/density_data.rds")
  density_dataset <- density_dataset %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
      Region    = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor)
    )
  # Load the pre-built T1D dataset (re-centered _c variables on T1D subset)
  # rather than filtering the full dataset, which keeps full-sample centering
  density_t1d <- readRDS("Data/density_data_t1d.rds") %>%
    mutate(
      Region = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex    = factor(Sex, levels = c("Female", "Male")),
      Donor  = factor(Donor)
    ) %>% droplevels()
  
  # Load models
  cat("Loading models...\n")
  density_diag  <- readRDS("Density/Models/density_diag.rds")
  density_t1d_model <- readRDS("Density/Models/density_t1d.rds")
  
  # Generate Excel files
  generate_density_diag_excel(density_diag, density_dataset, output_dir)
  generate_density_t1d_excel(density_t1d_model, density_t1d, output_dir)
  
  cat("\n=============================================================================\n")
  cat("COMPLETE! Two Excel files created in", output_dir, ":\n")
  cat("  1. Islet_Density_Diagnosis_Analysis.xlsx\n")
  cat("  2. Islet_Density_T1D_Analysis.xlsx\n")
  cat("=============================================================================\n")
}

# Run if executed directly
if (interactive()) {
  cat("\nTo generate all density Excel files, run:\n")
  cat("  generate_all_density_excel()\n\n")
  cat("Or individually:\n")
  cat("  density_diag <- readRDS('Density/Models/density_diag.rds')\n")
  cat("  density_dataset <- readRDS('Data/density_data.rds')\n")
  cat("  generate_density_diag_excel(density_diag, density_dataset)\n\n")
  cat("  density_t1d_model <- readRDS('Density/Models/density_t1d.rds')\n")
  cat("  density_t1d_data  <- readRDS('Data/density_data_t1d.rds')   # use pre-centered T1D data\n")
  cat("  generate_density_t1d_excel(density_t1d_model, density_t1d_data)\n")
}
