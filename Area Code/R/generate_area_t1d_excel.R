# =============================================================================
# Generate Comprehensive Excel Results for Cell Area T1D Models (glmmTMB)
# AO parameterization — DD and AO slopes extracted from single model
# =============================================================================
#
# For each cell type (Beta, Alpha, Delta, PP):
#   3 subset models (Singlet, SEO, Islet) + 1 combined model
#   = 4 Excel files total
#
# All models are T1D only, AO parameterization.
# Note: DD + Age and DD + AO models are algebraically equivalent.
# The AO model is preferred because both DD and AO slopes have direct
# clinical interpretations.  DD slopes are extracted from this model directly.
#
# Mixed families: Gamma(log), Gaussian(raw), Gaussian(log-transformed)
# All emmeans extraction unified on log scale via regrid for Gaussian raw.
#
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(writexl)

# =============================================================================
# HELPER FUNCTIONS (matching generate_area_excel_v3.R)
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

detect_model_type <- function(model) {
  fam <- family(model)$family
  if (fam == "Gamma") return("gamma_log")
  resp <- as.character(formula(model)[[2]])
  if (grepl("^log_", resp)) return("gaussian_log")
  return("gaussian_raw")
}

get_emm_log_scale <- function(model, specs, data = NULL, ...) {
  mt <- detect_model_type(model)
  if (!is.null(data)) {
    em <- emmeans(model, specs, data = data, ...)
  } else {
    em <- emmeans(model, specs, ...)
  }
  if (mt == "gaussian_raw") em <- regrid(em, transform = "log")
  em
}

extract_emm_response <- function(em_summary_df) {
  ci <- find_ci(em_summary_df)
  em_summary_df$mean_area <- exp(em_summary_df$emmean)
  em_summary_df$lower.CI  <- exp(em_summary_df[[ci$lower]])
  em_summary_df$upper.CI  <- exp(em_summary_df[[ci$upper]])
  em_summary_df
}

extract_contrast_results <- function(contr_summary_df) {
  contr_summary_df$log_estimate <- contr_summary_df$estimate
  contr_summary_df$pct_diff     <- (exp(contr_summary_df$estimate) - 1) * 100
  contr_summary_df
}

extract_fe <- function(model, subset_name) {
  mt <- detect_model_type(model)
  fe <- summary(model)$coefficients$cond
  ci <- tryCatch(confint(model, method = "Wald", parm = "beta_"), error = function(e) {
    est <- fe[, "Estimate"]; se <- fe[, "Std. Error"]
    cbind(est - 1.96 * se, est + 1.96 * se)
  })
  tibble(
    subset = subset_name, family = mt, term = rownames(fe),
    estimate = fmt(fe[, "Estimate"], 4), std.error = fmt(fe[, "Std. Error"], 4),
    `l-95% CI` = fmt(ci[rownames(fe), 1], 4), `u-95% CI` = fmt(ci[rownames(fe), 2], 4),
    z.value = fmt(fe[, "z value"], 2), p.value = fmt(fe[, "Pr(>|z|)"], 4),
    sig = add_stars(fe[, "Pr(>|z|)"])
  )
}

extract_re <- function(model, subset_name) {
  vc <- VarCorr(model)$cond
  re_df <- tibble(subset = subset_name, family = detect_model_type(model),
                  group = names(vc),
                  variance = fmt(sapply(vc, function(x) x[1]), 4),
                  std_dev = fmt(sqrt(sapply(vc, function(x) x[1])), 4))
  bind_rows(re_df, tibble(subset = subset_name, family = detect_model_type(model),
                          group = "Residual",
                          variance = fmt(sigma(model)^2, 4),
                          std_dev = fmt(sigma(model), 4)))
}

# emtrends on log scale (for slope extraction)
# emtrends — slopes are always returned on the model's native scale
# For Gamma(log): slopes are on the log scale (change in log(area) per unit)
# For Gaussian(log-transformed): slopes are on the log scale
# For Gaussian(raw): slopes are on the raw scale (change in µm² per unit)
# NOTE: regrid does NOT work on emtrends objects, so we do not attempt it
get_emtrends_log_scale <- function(model, specs, var, data = NULL, ...) {
  if (!is.null(data)) {
    em <- emtrends(model, specs, var = var, data = data, ...)
  } else {
    em <- emtrends(model, specs, var = var, ...)
  }
  em
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


# =============================================================================
# MAIN FUNCTION (parameterized for DD or AO)
# =============================================================================

generate_area_t1d_excel <- function(
    singlet_model, multicell_model, islet_model, combined_model,
    singlet_data, multicell_data, islet_data, combined_data,
    cell_type, cell_label, area_col,
    model_type = "dd",
    output_dir = "Area/Results/Total_Area"
) {

  dd_var    <- "Disease.Duration_c"
  age_var   <- if (model_type == "dd") "Age_c" else "age_at_onset_c"
  age_label <- if (model_type == "dd") "Age" else "Age at Onset"
  type_label <- if (model_type == "dd") "Disease Duration" else "Age at Onset"
  type_tag   <- toupper(model_type)

  cat("\n=============================================================================\n")
  cat("GENERATING EXCEL:", toupper(cell_label), "CELL AREA —", type_tag, "\n")
  cat("=============================================================================\n\n")

  models <- list(Singlet = singlet_model, `SEO` = multicell_model, Islet = islet_model)
  datasets <- list(Singlet = singlet_data, `SEO` = multicell_data, Islet = islet_data)
  subset_names <- c("Singlet", "SEO", "Islet")
  model_types <- sapply(models, detect_model_type)
  combined_type <- detect_model_type(combined_model)

  cat("Model families:\n")
  for (nm in names(model_types)) cat("  ", nm, "->", model_types[[nm]], "\n")
  cat("  Combined ->", combined_type, "\n\n")

  # =========================================================================
  # SHEET 1: MODEL SUMMARY
  # =========================================================================
  cat("Creating Sheet 1: Model Summary...\n")

  all_fe <- bind_rows(extract_fe(singlet_model, "Singlet"),
                      extract_fe(multicell_model, "SEO"),
                      extract_fe(islet_model, "Islet"))

  family_desc <- paste(sapply(names(model_types), function(nm) paste0(nm, ": ", model_types[[nm]])), collapse = "; ")

  model_notes <- tibble(Note = c(
    paste("CELL AREA ANALYSIS:", toupper(cell_label), "CELLS —", toupper(type_label), "(T1D only)"),
    "", "Generalized linear mixed models (glmmTMB)",
    paste("Response:", area_col, "in µm² (or log-transformed where noted)"),
    paste("Families:", family_desc),
    paste0("Subset models: area ~ (", dd_var, " + Region + Sex + ", age_var, " [+ log_Islet.Cells_c])^2 +"),
    paste0("               ", dd_var, ":Region:", age_var, " [+ ", dd_var, ":Region:log_Islet.Cells_c] + (1|Donor) + (1|Donor:ImageID)"),
    "", "Reference levels: Head (Region), Female (Sex), Singlets (object_type)",
    "Contrasts reported as log-ratios and % change for all families.",
    "Multiplicity: Tukey for Region pairwise, Holm for simple effects, None for single contrasts.", ""))

  model_summary <- bind_rows(model_notes,
    all_fe %>% mutate(across(everything(), as.character)) %>%
      mutate(Note = NA_character_) %>% select(Note, everything()))

  # =========================================================================
  # SHEET 2: OBJECT TYPE (combined model)
  # =========================================================================
  cat("Creating Sheet 2: Object Type...\n")

  em_obj <- get_emm_log_scale(combined_model, ~ object_type, data = combined_data)
  em_obj_df <- extract_emm_response(as.data.frame(summary(em_obj)))
  obj_means <- tibble(object_type = em_obj_df$object_type,
                      mean_area = fmt(em_obj_df$mean_area, 1),
                      `lower.CI` = fmt(em_obj_df$lower.CI, 1),
                      `upper.CI` = fmt(em_obj_df$upper.CI, 1))
  contr_obj <- pairs(em_obj, adjust = "tukey")
  contr_obj_df <- extract_contrast_results(as.data.frame(summary(contr_obj)))
  obj_contrasts <- tibble(contrast = as.character(contr_obj_df$contrast),
                          log_estimate = fmt(contr_obj_df$log_estimate, 4),
                          SE = fmt(contr_obj_df$SE, 4),
                          p.value = fmt(contr_obj_df$p.value, 4),
                          pct_diff = fmt(contr_obj_df$pct_diff, 1),
                          sig = add_stars(contr_obj_df$p.value))

  object_type_sheet <- bind_rows(
    tibble(Note = "OBJECT TYPE COMPARISON (from combined model)"),
    tibble(Note = "Cell area across Singlet, SEO, and Islet, averaged over Region and Sex."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (µm²)"),
    obj_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "PAIRWISE CONTRASTS (Tukey-adjusted)"),
    obj_contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # SHEET 3: OBJECT TYPE × DD SLOPES (combined model)
  # =========================================================================
  cat("Creating Sheet 3: Object Type x DD Slopes...\n")

  em_dd_by_obj <- get_emtrends_log_scale(combined_model, ~ object_type, var = dd_var, data = combined_data)
  dd_by_obj_df <- format_trend_summary(em_dd_by_obj, "object_type")
  contr_dd_obj <- pairs(em_dd_by_obj, adjust = "tukey")
  dd_obj_contr <- format_trend_contrasts(contr_dd_obj)

  obj_dd_sheet <- bind_rows(
    tibble(Note = "OBJECT TYPE × DISEASE DURATION SLOPES (combined model)"),
    tibble(Note = "How the Disease Duration effect on cell area varies by object type."),
    tibble(Note = NA),
    tibble(Note = "DISEASE DURATION SLOPES BY OBJECT TYPE"),
    dd_by_obj_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "PAIRWISE CONTRASTS ON DD SLOPES (Tukey-adjusted)"),
    dd_obj_contr %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # SHEET 4: OBJECT TYPE × AO SLOPES (combined model)
  # =========================================================================
  cat("Creating Sheet 4: Object Type x AO Slopes...\n")

  em_ao_by_obj <- get_emtrends_log_scale(combined_model, ~ object_type, var = age_var, data = combined_data)
  ao_by_obj_df <- format_trend_summary(em_ao_by_obj, "object_type")
  contr_ao_obj <- pairs(em_ao_by_obj, adjust = "tukey")
  ao_obj_contr <- format_trend_contrasts(contr_ao_obj)

  obj_ao_sheet <- bind_rows(
    tibble(Note = paste0("OBJECT TYPE × ", toupper(age_label), " SLOPES (combined model)")),
    tibble(Note = paste0("How the ", age_label, " effect on cell area varies by object type.")),
    tibble(Note = NA),
    tibble(Note = paste0(toupper(age_label), " SLOPES BY OBJECT TYPE")),
    ao_by_obj_df %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = paste0("PAIRWISE CONTRASTS ON ", toupper(age_label), " SLOPES (Tukey-adjusted)")),
    ao_obj_contr %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # SHEETS 5-6: REGION AND SEX (from subset models)
  # =========================================================================
  cat("Creating Sheets 4-5: Region and Sex...\n")

  get_factor_results <- function(model, subset_name, spec, data, adjust = "tukey") {
    em <- get_emm_log_scale(model, spec, data = data)
    em_df <- extract_emm_response(as.data.frame(summary(em)))
    ci <- find_ci(em_df)
    fac_col <- setdiff(names(em_df), c("emmean", "SE", "df", ci$lower, ci$upper, "mean_area", "lower.CI", "upper.CI"))
    means <- tibble(subset = subset_name, family = detect_model_type(model)) %>%
      bind_cols(em_df[fac_col]) %>%
      mutate(mean_area = fmt(em_df$mean_area, 1), `lower.CI` = fmt(em_df$lower.CI, 1), `upper.CI` = fmt(em_df$upper.CI, 1))
    contr <- pairs(em, adjust = adjust, reverse = if (adjust == "none") TRUE else FALSE)
    contr_df <- extract_contrast_results(as.data.frame(summary(contr)))
    contrasts <- tibble(subset = subset_name, family = detect_model_type(model),
                        contrast = as.character(contr_df$contrast),
                        log_estimate = fmt(contr_df$log_estimate, 4), SE = fmt(contr_df$SE, 4),
                        p.value = fmt(contr_df$p.value, 4), pct_diff = fmt(contr_df$pct_diff, 1),
                        sig = add_stars(contr_df$p.value))
    list(means = means, contrasts = contrasts)
  }

  region_results <- map(subset_names, ~get_factor_results(models[[.x]], .x, ~ Region, datasets[[.x]], "tukey"))
  region_means <- bind_rows(map(region_results, "means"))
  region_contrasts <- bind_rows(map(region_results, "contrasts"))

  region_sheet <- bind_rows(
    tibble(Note = "REGION EFFECT (from subset models, T1D only)"),
    tibble(Note = "Main effect of Region, averaged over Sex."), tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (µm²)"),
    region_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "PAIRWISE CONTRASTS (Tukey-adjusted)"),
    region_contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  sex_results <- map(subset_names, ~get_factor_results(models[[.x]], .x, ~ Sex, datasets[[.x]], "none"))
  sex_means <- bind_rows(map(sex_results, "means"))
  sex_contrasts <- bind_rows(map(sex_results, "contrasts"))

  sex_sheet <- bind_rows(
    tibble(Note = "SEX EFFECT (from subset models, T1D only)"),
    tibble(Note = "Main effect of Sex, averaged over Region."), tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (µm²)"),
    sex_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "CONTRAST (Female vs Male)"),
    sex_contrasts %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # SHEET 6: DD SLOPES (from subset models)
  # =========================================================================
  cat("Creating Sheet 6: DD Slopes...\n")

  get_dd_slopes <- function(model, subset_name, data) {
    mt <- detect_model_type(model)
    em_overall <- get_emtrends_log_scale(model, ~ 1, var = dd_var, data = data)
    overall_df <- format_trend_summary(em_overall, character(0)) %>% mutate(subset = subset_name, family = mt, .before = 1)
    em_by_region <- get_emtrends_log_scale(model, ~ Region, var = dd_var, data = data)
    by_region_df <- format_trend_summary(em_by_region, "Region") %>% mutate(subset = subset_name, family = mt, .before = 1)
    contr_region <- pairs(em_by_region, adjust = "tukey")
    region_contr_df <- format_trend_contrasts(contr_region) %>% mutate(subset = subset_name, family = mt, .before = 1)
    list(overall = overall_df, by_region = by_region_df, region_contr = region_contr_df)
  }

  dd_slope_results <- map(subset_names, ~get_dd_slopes(models[[.x]], .x, datasets[[.x]]))

  dd_slopes_sheet <- bind_rows(
    tibble(Note = "DISEASE DURATION SLOPES (from subset models, T1D only)"),
    tibble(Note = "For gamma_log and gaussian_log: slopes are on the log scale."),
    tibble(Note = "  Approximate % change per year = (exp(slope) - 1) * 100."),
    tibble(Note = "For gaussian_raw: slopes are on the raw µm² scale (additive change per year)."),
    tibble(Note = "Check the 'family' column to determine the scale for each subset."),
    tibble(Note = "Because of DD × Region interactions, slopes vary by Region."),
    tibble(Note = NA),
    tibble(Note = "OVERALL DD SLOPES (marginalized over Region and Sex)"),
    bind_rows(map(dd_slope_results, "overall")) %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DD SLOPES BY REGION"),
    bind_rows(map(dd_slope_results, "by_region")) %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS ON DD SLOPES (Tukey-adjusted)"),
    bind_rows(map(dd_slope_results, "region_contr")) %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # SHEET 7: AGE/AO SLOPES (from subset models)
  # =========================================================================
  cat("Creating Sheet 7:", age_label, "Slopes...\n")

  get_age_slopes <- function(model, subset_name, data) {
    mt <- detect_model_type(model)
    em_overall <- get_emtrends_log_scale(model, ~ 1, var = age_var, data = data)
    overall_df <- format_trend_summary(em_overall, character(0)) %>% mutate(subset = subset_name, family = mt, .before = 1)
    em_by_region <- get_emtrends_log_scale(model, ~ Region, var = age_var, data = data)
    by_region_df <- format_trend_summary(em_by_region, "Region") %>% mutate(subset = subset_name, family = mt, .before = 1)
    contr_region <- pairs(em_by_region, adjust = "tukey")
    region_contr_df <- format_trend_contrasts(contr_region) %>% mutate(subset = subset_name, family = mt, .before = 1)
    list(overall = overall_df, by_region = by_region_df, region_contr = region_contr_df)
  }

  age_slope_results <- map(subset_names, ~get_age_slopes(models[[.x]], .x, datasets[[.x]]))

  age_slopes_sheet <- bind_rows(
    tibble(Note = paste(toupper(age_label), "SLOPES (from subset models, T1D only)")),
    tibble(Note = "For gamma_log and gaussian_log: slopes are on the log scale."),
    tibble(Note = "For gaussian_raw: slopes are on the raw µm² scale."),
    tibble(Note = "Check the 'family' column to determine the scale for each subset."),
    tibble(Note = NA),
    tibble(Note = paste("OVERALL", toupper(age_label), "SLOPES")),
    bind_rows(map(age_slope_results, "overall")) %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = paste(toupper(age_label), "SLOPES BY REGION")),
    bind_rows(map(age_slope_results, "by_region")) %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = paste("REGION CONTRASTS ON", toupper(age_label), "SLOPES (Tukey-adjusted)")),
    bind_rows(map(age_slope_results, "region_contr")) %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # SHEET 8: ISLET SIZE SLOPES
  # =========================================================================
  cat("Creating Sheet 8: Islet Size Slopes...\n")

  get_islet_slopes <- function(model, subset_name, data) {
    mt <- detect_model_type(model)
    fe <- summary(model)$coefficients$cond
    islet_row <- which(rownames(fe) == "log_Islet.Cells_c")
    if (length(islet_row) == 0) return(NULL)
    est <- fe[islet_row, "Estimate"]
    pct <- if (mt == "gaussian_raw") {
      emm_mean <- tryCatch({ emm <- emmeans(model, ~1, data = data); summary(emm)$emmean[1] },
                           error = function(e) fe["(Intercept)", "Estimate"])
      (est / emm_mean) * 100
    } else { (exp(est) - 1) * 100 }
    tibble(subset = subset_name, family = mt, predictor = "log_Islet.Cells_c",
           estimate = fmt(est, 4), SE = fmt(fe[islet_row, "Std. Error"], 4),
           p.value = fmt(fe[islet_row, "Pr(>|z|)"], 4),
           pct_per_log_unit = fmt(pct, 2), sig = add_stars(fe[islet_row, "Pr(>|z|)"]))
  }

  islet_slopes <- bind_rows(
    get_islet_slopes(multicell_model, "SEO", multicell_data),
    get_islet_slopes(islet_model, "Islet", islet_data),
    get_islet_slopes(combined_model, "Combined", combined_data)
  ) %>% compact()

  islet_sheet <- if (nrow(islet_slopes) > 0) {
    bind_rows(
      tibble(Note = "ISLET SIZE EFFECT"),
      tibble(Note = "Change in cell area per 1 unit increase in log(Islet.Cells)."),
      tibble(Note = "1 unit on log scale ≈ 2.72-fold increase in cell count."),
      tibble(Note = NA),
      islet_slopes %>% mutate(Note = NA) %>% select(Note, everything()))
  } else { tibble(Note = "Islet size slopes not available.") }

  # =========================================================================
  # SHEET 9: DD × REGION × AGE/AO (3-way from subset models)
  # =========================================================================
  cat("Creating Sheet 9: DD x Region x", age_label, "...\n")

  get_3way_slopes <- function(model, subset_name, data) {
    mt <- detect_model_type(model)
    em_lo <- tryCatch({
      get_emtrends_log_scale(model, ~ Region, var = dd_var, data = data, at = setNames(list(-1), age_var))
    }, error = function(e) NULL)
    em_hi <- tryCatch({
      get_emtrends_log_scale(model, ~ Region, var = dd_var, data = data, at = setNames(list(1), age_var))
    }, error = function(e) NULL)
    if (is.null(em_lo) || is.null(em_hi)) return(NULL)

    lo_df <- format_trend_summary(em_lo, "Region") %>%
      mutate(subset = subset_name, family = mt, level = paste0(age_label, " = mean - 1"), .before = 1)
    hi_df <- format_trend_summary(em_hi, "Region") %>%
      mutate(subset = subset_name, family = mt, level = paste0(age_label, " = mean + 1"), .before = 1)
    bind_rows(lo_df, hi_df)
  }

  three_way_results <- map(subset_names, ~get_3way_slopes(models[[.x]], .x, datasets[[.x]])) %>% compact()
  three_way_df <- if (length(three_way_results) > 0) bind_rows(three_way_results) else NULL

  three_way_sheet <- if (!is.null(three_way_df)) {
    bind_rows(
      tibble(Note = paste("DD × REGION ×", toupper(age_label), "(3-way interaction)")),
      tibble(Note = paste("DD slopes evaluated at", age_label, "= mean ± 1, by Region, from each subset model.")),
      tibble(Note = NA),
      three_way_df %>% mutate(Note = NA) %>% select(Note, everything()))
  } else { tibble(Note = "3-way interaction slopes not available.") }

  # =========================================================================
  # SHEET 10: COMBINED MODEL SUMMARY
  # =========================================================================
  cat("Creating Sheet 10: Combined Model Summary...\n")

  fe_combined <- summary(combined_model)$coefficients$cond
  ci_combined <- tryCatch(confint(combined_model, method = "Wald", parm = "beta_"), error = function(e) {
    est <- fe_combined[, "Estimate"]; se <- fe_combined[, "Std. Error"]
    cbind(est - 1.96 * se, est + 1.96 * se)
  })
  combined_fe <- tibble(
    term = rownames(fe_combined), estimate = fmt(fe_combined[, "Estimate"], 4),
    std.error = fmt(fe_combined[, "Std. Error"], 4),
    `l-95% CI` = fmt(ci_combined[rownames(fe_combined), 1], 4),
    `u-95% CI` = fmt(ci_combined[rownames(fe_combined), 2], 4),
    z.value = fmt(fe_combined[, "z value"], 2),
    p.value = fmt(fe_combined[, "Pr(>|z|)"], 4),
    pct_change = fmt((exp(fe_combined[, "Estimate"]) - 1) * 100, 1),
    sig = add_stars(fe_combined[, "Pr(>|z|)"]))

  combined_sheet <- bind_rows(
    tibble(Note = "COMBINED MODEL (all object types, T1D only)"),
    tibble(Note = paste0("Model: area ~ (", dd_var, " + Region + Sex + ", age_var, ")^2 +")),
    tibble(Note = paste0("       ", dd_var, ":Region:", age_var, " + log_Islet.Cells_c +")),
    tibble(Note = paste0("       object_type * ", dd_var, " + object_type * Region + (1|Donor) + (1|Donor:ImageID)")),
    tibble(Note = paste("Family:", combined_type)),
    tibble(Note = paste("N =", nrow(combined_data))),
    tibble(Note = NA),
    combined_fe %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # SHEET 11: RANDOM EFFECTS & SAMPLE SIZES
  # =========================================================================
  cat("Creating Sheet 11: Random Effects & Sample Sizes...\n")

  random_effects <- bind_rows(
    extract_re(singlet_model, "Singlet"), extract_re(multicell_model, "SEO"),
    extract_re(islet_model, "Islet"), extract_re(combined_model, "Combined"))

  sample_sizes <- tibble(
    subset = c("Singlet", "SEO", "Islet", "Combined"),
    n_obs = c(nrow(singlet_data), nrow(multicell_data), nrow(islet_data), nrow(combined_data)),
    n_donors = c(length(unique(singlet_data$Donor)), length(unique(multicell_data$Donor)),
                 length(unique(islet_data$Donor)), length(unique(combined_data$Donor))))

  re_sheet <- bind_rows(
    tibble(Note = "RANDOM EFFECTS VARIANCE COMPONENTS"), tibble(Note = NA),
    random_effects %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA), tibble(Note = "SAMPLE SIZES"),
    sample_sizes %>% mutate(Note = NA) %>% select(Note, everything()))

  # =========================================================================
  # ASSEMBLE AND SAVE
  # =========================================================================
  cat("\nAssembling workbook...\n")

  sheet_data <- list(
    as.data.frame(model_summary),
    as.data.frame(object_type_sheet),
    as.data.frame(obj_dd_sheet),
    as.data.frame(obj_ao_sheet),
    as.data.frame(region_sheet),
    as.data.frame(sex_sheet),
    as.data.frame(dd_slopes_sheet),
    as.data.frame(age_slopes_sheet),
    as.data.frame(islet_sheet),
    as.data.frame(three_way_sheet),
    as.data.frame(combined_sheet),
    as.data.frame(re_sheet)
  )
  sheet_names <- c(
    "Model Summary",
    "Object Type",
    "Object Type x DD",
    paste("Object Type x", age_label),
    "Region",
    "Sex",
    "DD Slopes",
    paste(age_label, "Slopes"),
    "Islet Size Slopes",
    paste("DD x Reg x", age_label),
    "Combined Model",
    "Random Effects"
  )
  # Truncate to Excel 31-char limit
  names(sheet_data) <- substr(sheet_names, 1, 31)
  sheets <- sheet_data

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0(cell_label, "_Cell_Area_", type_tag, "_Analysis.xlsx"))
  write_xlsx(sheets, output_file)
  cat("\nSaved to:", output_file, "\n")
  return(invisible(sheets))
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

generate_all_area_t1d_excel <- function(output_dir = "Area/Results/Total_Area") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL CELL AREA T1D EXCEL RESULTS (AO parameterization)\n")
  cat("=============================================================================\n\n")

  # --- Helper: same split_area as models (filters area > 0 AND count > 0) ---
  area_to_count <- c(ins_area = "Count.ins", glu_area = "Count.glu",
                     soma_area = "Count.soma", pp_area = "Count.PP")

  split_area <- function(data, cell_col) {
    log_col   <- paste0("log_", cell_col)
    count_col <- area_to_count[[cell_col]]
    data %>%
      filter(!is.na(!!sym(cell_col)) & !!sym(cell_col) > 0 &
             !!sym(count_col) > 0) %>%
      mutate(!!log_col := log(!!sym(cell_col)))
  }

  set_factors <- function(df) {
    df %>% mutate(
      Region    = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor),
      ImageID   = factor(ImageID)) %>%
    droplevels()
  }

  add_object_type <- function(df) {
    df %>% mutate(
      object_type = factor(case_when(
        Islet.Cells == 1  ~ "Singlets",
        Islet.Cells >= 2 & Islet.Cells <= 14 ~ "SEOs",
        Islet.Cells >= 15 ~ "Islets"),
        levels = c("Singlets", "SEOs", "Islets")))
  }

  # --- Load subset-specific T1D datasets (matching model training data) ---
  cat("Loading subset-specific T1D data...\n")
  area_singlets_t1d <- set_factors(readRDS("Data/area_data_singlets_t1d.rds"))
  area_endobs_t1d   <- set_factors(readRDS("Data/area_data_endobs_t1d.rds"))
  area_islets_t1d   <- set_factors(readRDS("Data/area_data_islets_t1d.rds"))

  beta_combined_t1d  <- add_object_type(set_factors(readRDS("Data/area_data_beta_combined_t1d.rds")))
  alpha_combined_t1d <- add_object_type(set_factors(readRDS("Data/area_data_alpha_combined_t1d.rds")))
  delta_combined_t1d <- add_object_type(set_factors(readRDS("Data/area_data_delta_combined_t1d.rds")))
  pp_combined_t1d    <- add_object_type(set_factors(readRDS("Data/area_data_pp_combined_t1d.rds")))
  cat("  Data loaded.\n\n")

  combined_data_map <- list(
    ins  = beta_combined_t1d,
    glu  = alpha_combined_t1d,
    soma = delta_combined_t1d,
    pp   = pp_combined_t1d
  )

  # Cell type definitions — model files use _t1d suffix
  cell_defs <- list(
    list(type = "ins",  label = "β",  area_col = "ins_area",
         prefix_1 = "area_beta_1cell", prefix_2 = "area_beta_2to14",
         prefix_15 = "area_beta_15plus", prefix_comb = "beta_combined"),
    list(type = "glu",  label = "α", area_col = "glu_area",
         prefix_1 = "area_alpha_1cell", prefix_2 = "area_alpha_2to14",
         prefix_15 = "area_alpha_15plus", prefix_comb = "alpha_combined"),
    list(type = "soma", label = "δ", area_col = "soma_area",
         prefix_1 = "area_delta_1cell", prefix_2 = "area_delta_2to14",
         prefix_15 = "area_delta_15plus", prefix_comb = "delta_combined"),
    list(type = "pp",   label = "PP",    area_col = "pp_area",
         prefix_1 = "area_pp_1cell", prefix_2 = "area_pp_2to14",
         prefix_15 = "area_pp_15plus", prefix_comb = "pp_combined")
  )

  cat("--- AGE AT ONSET MODELS (T1D) ---\n")
  for (cd in cell_defs) {
    cat("\nProcessing", cd$label, "cell...\n")

    # Load models — all use _t1d suffix
    m1  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_1, "_t1d.rds"))
    m2  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_2, "_t1d.rds"))
    m15 <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_15, "_t1d.rds"))
    mc  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_comb, "_t1d.rds"))

    # Prepare data subsets from correct source datasets (matching models)
    d1  <- split_area(area_singlets_t1d, cd$area_col)
    d2  <- split_area(area_endobs_t1d, cd$area_col)
    d15 <- split_area(area_islets_t1d, cd$area_col)
    dc  <- split_area(combined_data_map[[cd$type]], cd$area_col)

    generate_area_t1d_excel(
      singlet_model = m1, multicell_model = m2,
      islet_model = m15, combined_model = mc,
      singlet_data = d1, multicell_data = d2,
      islet_data = d15, combined_data = dc,
      cell_type = cd$type, cell_label = cd$label,
      area_col = cd$area_col,
      model_type = "ao", output_dir = output_dir
    )
    rm(m1, m2, m15, mc); gc()
  }

  cat("\n=============================================================================\n")
  cat("COMPLETE! Four Excel files created in", output_dir, ":\n\n")
  for (cd in cell_defs) cat("  ", cd$label, "_Cell_Area_AO_Analysis.xlsx\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all area T1D Excel files, run:\n")
  cat("  generate_all_area_t1d_excel()\n")
}
