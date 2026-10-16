# =============================================================================
# Generate Comprehensive Excel Results for Cell Area Models (glmmTMB)
# Version 3: Handles mixed distributional families (Gamma, Gaussian raw/log)
# All emmeans extraction is unified on the log scale via regrid for Gaussian raw
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

get_ci_cols <- function(df) {
  lower <- if ("lower.CL" %in% names(df)) "lower.CL" else if ("asymp.LCL" %in% names(df)) "asymp.LCL" else NULL
  upper <- if ("upper.CL" %in% names(df)) "upper.CL" else if ("asymp.UCL" %in% names(df)) "asymp.UCL" else NULL
  list(lower = lower, upper = upper)
}

# ---------------------------------------------------------------------------
# Model type detection
# ---------------------------------------------------------------------------
detect_model_type <- function(model) {
  fam <- family(model)$family
  if (fam == "Gamma") return("gamma_log")
  # Gaussian — check response variable name
  resp <- as.character(formula(model)[[2]])
  if (grepl("^log_", resp)) return("gaussian_log")
  return("gaussian_raw")
}

# ---------------------------------------------------------------------------
# get_emm_log_scale: returns emmGrid on log scale for ANY model type
#   - Gamma(log link): emmeans without type="response" → already log scale
#   - Gaussian(log-transformed): emmeans → already log scale (response IS log)
#   - Gaussian(raw): emmeans → raw scale, regrid to log via delta method
# ---------------------------------------------------------------------------
get_emm_log_scale <- function(model, specs, data = NULL, ...) {
  mt <- detect_model_type(model)
  if (!is.null(data)) {
    em <- emmeans(model, specs, data = data, ...)
  } else {
    em <- emmeans(model, specs, ...)
  }
  if (mt == "gaussian_raw") {
    em <- regrid(em, transform = "log")
  }
  em
}

# ---------------------------------------------------------------------------
# extract_emm_response: converts log-scale emmGrid summary to response scale
# All emmGrids arrive on log scale (via get_emm_log_scale), so always exponentiate
# ---------------------------------------------------------------------------
extract_emm_response <- function(em_summary_df) {
  ci_cols <- get_ci_cols(em_summary_df)
  
  # All inputs are on log scale — exponentiate
  em_summary_df$mean_area <- exp(em_summary_df$emmean)
  em_summary_df$lower.CI  <- exp(em_summary_df[[ci_cols$lower]])
  em_summary_df$upper.CI  <- exp(em_summary_df[[ci_cols$upper]])
  
  em_summary_df
}

# ---------------------------------------------------------------------------
# extract_contrast_results: all contrasts are on log scale (log-differences)
# Compute pct_diff = (exp(estimate) - 1) * 100
# ---------------------------------------------------------------------------
extract_contrast_results <- function(contr_summary_df) {
  contr_summary_df$log_estimate <- contr_summary_df$estimate
  contr_summary_df$pct_diff     <- (exp(contr_summary_df$estimate) - 1) * 100
  contr_summary_df
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

generate_cell_type_excel <- function(
    singlet_model, 
    multicell_model, 
    islet_model,
    combined_model,
    singlet_data,
    multicell_data,
    islet_data,
    combined_data,
    cell_type,
    cell_label,
    area_col,
    output_dir = "Area/Results/Total_Area"
) {
  
  cat("\n")
  cat("=============================================================================\n")
  cat("GENERATING EXCEL FOR:", toupper(cell_label), "CELL AREA\n")
  cat("=============================================================================\n\n")
  
  # Store models and data in lists
  models <- list(Singlet = singlet_model, `SEO` = multicell_model, Islet = islet_model)
  datasets <- list(Singlet = singlet_data, `SEO` = multicell_data, Islet = islet_data)
  subset_names <- c("Singlet", "SEO", "Islet")
  
  # Detect families for each model
  model_types <- sapply(models, detect_model_type)
  combined_type <- detect_model_type(combined_model)
  
  cat("Model families:\n")
  for (nm in names(model_types)) cat("  ", nm, "->", model_types[[nm]], "\n")
  cat("  Combined ->", combined_type, "\n\n")
  
  # =========================================================================
  # SHEET 1: MODEL SUMMARY
  # =========================================================================
  cat("Creating Sheet 1: Model Summary...\n")
  
  extract_fe <- function(model, subset_name) {
    mt <- detect_model_type(model)
    fe <- summary(model)$coefficients$cond
    ci <- tryCatch({
      confint(model, method = "Wald", parm = "beta_")
    }, error = function(e) {
      est <- fe[, "Estimate"]
      se <- fe[, "Std. Error"]
      cbind(est - 1.96 * se, est + 1.96 * se)
    })
    
    tibble(
      subset = subset_name,
      family = mt,
      term = rownames(fe),
      estimate = fmt(fe[, "Estimate"], 4),
      std.error = fmt(fe[, "Std. Error"], 4),
      `l-95% CI` = fmt(ci[rownames(fe), 1], 4),
      `u-95% CI` = fmt(ci[rownames(fe), 2], 4),
      z.value = fmt(fe[, "z value"], 2),
      p.value = fmt(fe[, "Pr(>|z|)"], 4),
      sig = add_stars(fe[, "Pr(>|z|)"])
    )
  }
  
  all_fe <- bind_rows(
    extract_fe(singlet_model, "Singlet"),
    extract_fe(multicell_model, "SEO"),
    extract_fe(islet_model, "Islet")
  )
  
  family_desc <- paste(
    sapply(names(model_types), function(nm) paste0(nm, ": ", model_types[[nm]])),
    collapse = "; "
  )
  
  model_notes <- tibble(
    Note = c(
      paste("CELL AREA ANALYSIS:", toupper(cell_label), "CELLS"),
      "",
      "Generalized linear mixed models (glmmTMB)",
      paste("Response:", area_col, "in µm² (or log-transformed where noted)"),
      paste("Families:", family_desc),
      "Subset models: area ~ (Diagnosis + Region + Sex + Age_c + [log_Islet.Cells_c])^2 +",
      "               Diagnosis:Region:[log_Islet.Cells_c] + Diagnosis:Region:Age_c + (1|Donor) + (1|Donor:ImageID)",
      "",
      "Reference levels: ND (Diagnosis), Body (Region), Female (Sex), Singlets (object_type)",
      "For Gamma(log link) and Gaussian(log-transformed): coefficients are on the log scale.",
      "  Exponentiate for multiplicative effects (ratios). % change = (exp(estimate) - 1) * 100",
      "For Gaussian(raw): coefficients are on the µm² scale (additive effects).",
      "Marginal means are reported on the response scale (µm²) for all families.",
      "Contrasts are reported as log-ratios and % change for all families",
      "  (Gaussian raw uses delta-method regridding to the log scale).",
      "",
      "Multiplicity adjustments:",
      "  Tukey: all-pairwise comparisons within a factor (Region 3-way, Object Type 3-way)",
      "  Holm: same contrast repeated across levels of a conditioning factor",
      "  None: single planned contrast (T1D vs ND main effect)",
      ""
    )
  )
  
  all_fe_char <- all_fe %>%
    mutate(across(everything(), as.character)) %>%
    mutate(Note = NA_character_) %>%
    select(Note, everything())
  
  model_summary <- bind_rows(
    model_notes %>% mutate(
      subset = NA_character_, family = NA_character_, term = NA_character_, 
      estimate = NA_character_,
      std.error = NA_character_, `l-95% CI` = NA_character_, `u-95% CI` = NA_character_,
      z.value = NA_character_, p.value = NA_character_, sig = NA_character_
    ),
    all_fe_char
  )
  
  # =========================================================================
  # SHEET 2: OBJECT TYPE COMPARISON (from combined model)
  # =========================================================================
  cat("Creating Sheet 2: Object Type Comparison...\n")
  
  em_obj <- get_emm_log_scale(combined_model, ~ object_type, data = combined_data)
  em_obj_df <- extract_emm_response(as.data.frame(summary(em_obj)))
  
  obj_means <- tibble(
    object_type = em_obj_df$object_type,
    mean_area = fmt(em_obj_df$mean_area, 1),
    `lower.CI` = fmt(em_obj_df$lower.CI, 1),
    `upper.CI` = fmt(em_obj_df$upper.CI, 1)
  )
  
  contr_obj <- pairs(em_obj, adjust = "tukey")
  contr_obj_df <- extract_contrast_results(as.data.frame(summary(contr_obj)))
  
  obj_contrasts <- tibble(
    contrast = as.character(contr_obj_df$contrast),
    log_estimate = fmt(contr_obj_df$log_estimate, 4),
    SE = fmt(contr_obj_df$SE, 4),
    p.value = fmt(contr_obj_df$p.value, 4),
    pct_diff = fmt(contr_obj_df$pct_diff, 1),
    sig = add_stars(contr_obj_df$p.value)
  )
  
  object_type_sheet <- bind_rows(
    tibble(Note = "OBJECT TYPE COMPARISON"),
    tibble(Note = "Comparison of cell area across Singlet, SEO, and Islet object types."),
    tibble(Note = "From combined model with all three subsets."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (µm²) - averaged over Diagnosis, Region, Sex"),
    obj_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "PAIRWISE CONTRASTS (Tukey-adjusted)"),
    obj_contrasts %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 3: OBJECT TYPE x DIAGNOSIS
  # =========================================================================
  cat("Creating Sheet 3: Object Type x Diagnosis...\n")
  
  em_obj_diag <- get_emm_log_scale(combined_model, ~ object_type * Diagnosis, data = combined_data)
  em_obj_diag_df <- extract_emm_response(as.data.frame(summary(em_obj_diag)))
  
  obj_diag_means <- tibble(
    object_type = em_obj_diag_df$object_type,
    Diagnosis = em_obj_diag_df$Diagnosis,
    mean_area = fmt(em_obj_diag_df$mean_area, 1),
    `lower.CI` = fmt(em_obj_diag_df$lower.CI, 1),
    `upper.CI` = fmt(em_obj_diag_df$upper.CI, 1)
  )
  
  em_diag_by_obj <- get_emm_log_scale(combined_model, ~ Diagnosis | object_type, data = combined_data)
  contr_diag_by_obj <- pairs(em_diag_by_obj, adjust = "holm", reverse = TRUE)
  contr_diag_by_obj_df <- extract_contrast_results(as.data.frame(summary(contr_diag_by_obj)))
  
  diag_within_obj <- tibble(
    object_type = contr_diag_by_obj_df$object_type,
    contrast = as.character(contr_diag_by_obj_df$contrast),
    log_estimate = fmt(contr_diag_by_obj_df$log_estimate, 4),
    SE = fmt(contr_diag_by_obj_df$SE, 4),
    p.value = fmt(contr_diag_by_obj_df$p.value, 4),
    pct_diff = fmt(contr_diag_by_obj_df$pct_diff, 1),
    sig = add_stars(contr_diag_by_obj_df$p.value)
  )
  
  em_obj_by_diag <- get_emm_log_scale(combined_model, ~ object_type | Diagnosis, data = combined_data)
  contr_obj_by_diag <- pairs(em_obj_by_diag, adjust = "tukey")
  contr_obj_by_diag_df <- extract_contrast_results(as.data.frame(summary(contr_obj_by_diag)))
  
  obj_within_diag <- tibble(
    Diagnosis = contr_obj_by_diag_df$Diagnosis,
    contrast = as.character(contr_obj_by_diag_df$contrast),
    log_estimate = fmt(contr_obj_by_diag_df$log_estimate, 4),
    SE = fmt(contr_obj_by_diag_df$SE, 4),
    p.value = fmt(contr_obj_by_diag_df$p.value, 4),
    pct_diff = fmt(contr_obj_by_diag_df$pct_diff, 1),
    sig = add_stars(contr_obj_by_diag_df$p.value)
  )
  
  obj_diag_sheet <- bind_rows(
    tibble(Note = "OBJECT TYPE x DIAGNOSIS INTERACTION"),
    tibble(Note = "Tests whether the T1D effect differs across object types."),
    tibble(Note = NA),
    tibble(Note = "CELL MEANS (µm²)"),
    obj_diag_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "T1D EFFECT WITHIN EACH OBJECT TYPE (ND vs T1D, Holm-adjusted)"),
    diag_within_obj %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "OBJECT TYPE CONTRASTS WITHIN EACH DIAGNOSIS (Tukey-adjusted)"),
    obj_within_diag %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 4: DIAGNOSIS x REGION SUMMARY
  # =========================================================================
  cat("Creating Sheet 4: Diagnosis x Region Summary...\n")
  
  get_diag_region_summary <- function(model, subset_name, data) {
    em <- get_emm_log_scale(model, ~ Diagnosis * Region, data = data)
    em_df <- extract_emm_response(as.data.frame(summary(em)))
    
    tibble(
      subset = subset_name,
      family = detect_model_type(model),
      Diagnosis = em_df$Diagnosis,
      Region = em_df$Region,
      mean_area = fmt(em_df$mean_area, 1),
      `lower.CI` = fmt(em_df$lower.CI, 1),
      `upper.CI` = fmt(em_df$upper.CI, 1)
    )
  }
  
  diag_region_summary <- bind_rows(
    get_diag_region_summary(singlet_model, "Singlet", singlet_data),
    get_diag_region_summary(multicell_model, "SEO", multicell_data),
    get_diag_region_summary(islet_model, "Islet", islet_data)
  )
  
  diag_region_wide <- diag_region_summary %>%
    mutate(cell = paste(Diagnosis, Region, sep = "_")) %>%
    select(subset, family, cell, mean_area) %>%
    pivot_wider(names_from = cell, values_from = mean_area)
  
  diag_region_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x REGION CELL MEANS SUMMARY"),
    tibble(Note = "Mean cell area (µm²) for each Diagnosis x Region combination, by object type subset."),
    tibble(Note = NA),
    tibble(Note = "DETAILED VIEW (long format)"),
    diag_region_summary %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "COMPARISON VIEW (wide format)"),
    diag_region_wide %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 5: DIAGNOSIS (from subset models)
  # =========================================================================
  cat("Creating Sheet 5: Diagnosis...\n")
  
  get_diagnosis_results <- function(model, subset_name, data) {
    em <- get_emm_log_scale(model, ~ Diagnosis, data = data)
    em_df <- extract_emm_response(as.data.frame(summary(em)))
    
    means <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      Diagnosis = em_df$Diagnosis,
      mean_area = fmt(em_df$mean_area, 1),
      `lower.CI` = fmt(em_df$lower.CI, 1),
      `upper.CI` = fmt(em_df$upper.CI, 1)
    )
    
    contr <- pairs(em, adjust = "none", reverse = TRUE)
    contr_df <- extract_contrast_results(as.data.frame(summary(contr)))
    
    contrasts <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      contrast = as.character(contr_df$contrast),
      log_estimate = fmt(contr_df$log_estimate, 4),
      SE = fmt(contr_df$SE, 4),
      p.value = fmt(contr_df$p.value, 4),
      pct_diff = fmt(contr_df$pct_diff, 1),
      sig = add_stars(contr_df$p.value)
    )
    
    list(means = means, contrasts = contrasts)
  }
  
  diag_results <- map(subset_names, ~get_diagnosis_results(models[[.x]], .x, datasets[[.x]]))
  diag_means <- bind_rows(map(diag_results, "means"))
  diag_contrasts <- bind_rows(map(diag_results, "contrasts"))
  
  diagnosis_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS EFFECT (from subset-specific models)"),
    tibble(Note = "Main effect of Diagnosis on cell area, averaged over Region and Sex."),
    tibble(Note = "Gaussian(raw) contrasts use delta-method log-scale regridding for comparable ratios."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (µm²)"),
    diag_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "CONTRASTS (ND vs T1D)"),
    diag_contrasts %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 6: REGION
  # =========================================================================
  cat("Creating Sheet 6: Region...\n")
  
  get_region_results <- function(model, subset_name, data) {
    em <- get_emm_log_scale(model, ~ Region, data = data)
    em_df <- extract_emm_response(as.data.frame(summary(em)))
    
    means <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      Region = em_df$Region,
      mean_area = fmt(em_df$mean_area, 1),
      `lower.CI` = fmt(em_df$lower.CI, 1),
      `upper.CI` = fmt(em_df$upper.CI, 1)
    )
    
    contr <- pairs(em, adjust = "tukey")
    contr_df <- extract_contrast_results(as.data.frame(summary(contr)))
    
    contrasts <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      contrast = as.character(contr_df$contrast),
      log_estimate = fmt(contr_df$log_estimate, 4),
      SE = fmt(contr_df$SE, 4),
      p.value = fmt(contr_df$p.value, 4),
      pct_diff = fmt(contr_df$pct_diff, 1),
      sig = add_stars(contr_df$p.value)
    )
    
    list(means = means, contrasts = contrasts)
  }
  
  region_results <- map(subset_names, ~get_region_results(models[[.x]], .x, datasets[[.x]]))
  region_means <- bind_rows(map(region_results, "means"))
  region_contrasts <- bind_rows(map(region_results, "contrasts"))
  
  region_sheet <- bind_rows(
    tibble(Note = "REGION EFFECT"),
    tibble(Note = "Main effect of pancreatic Region on cell area, averaged over Diagnosis and Sex."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (µm²)"),
    region_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "PAIRWISE CONTRASTS (Tukey-adjusted)"),
    region_contrasts %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 7: SEX
  # =========================================================================
  cat("Creating Sheet 7: Sex...\n")
  
  get_sex_results <- function(model, subset_name, data) {
    em <- get_emm_log_scale(model, ~ Sex, data = data)
    em_df <- extract_emm_response(as.data.frame(summary(em)))
    
    means <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      Sex = em_df$Sex,
      mean_area = fmt(em_df$mean_area, 1),
      `lower.CI` = fmt(em_df$lower.CI, 1),
      `upper.CI` = fmt(em_df$upper.CI, 1)
    )
    
    contr <- pairs(em, adjust = "none")
    contr_df <- extract_contrast_results(as.data.frame(summary(contr)))
    
    contrasts <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      contrast = as.character(contr_df$contrast),
      log_estimate = fmt(contr_df$log_estimate, 4),
      SE = fmt(contr_df$SE, 4),
      p.value = fmt(contr_df$p.value, 4),
      pct_diff = fmt(contr_df$pct_diff, 1),
      sig = add_stars(contr_df$p.value)
    )
    
    list(means = means, contrasts = contrasts)
  }
  
  sex_results <- map(subset_names, ~get_sex_results(models[[.x]], .x, datasets[[.x]]))
  sex_means <- bind_rows(map(sex_results, "means"))
  sex_contrasts <- bind_rows(map(sex_results, "contrasts"))
  
  sex_sheet <- bind_rows(
    tibble(Note = "SEX EFFECT"),
    tibble(Note = "Main effect of Sex on cell area, averaged over Diagnosis and Region."),
    tibble(Note = NA),
    tibble(Note = "MARGINAL MEANS (µm²)"),
    sex_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "CONTRAST (Female vs Male)"),
    sex_contrasts %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 8: DIAGNOSIS x REGION (detailed)
  # =========================================================================
  cat("Creating Sheet 8: Diagnosis x Region (detailed)...\n")
  
  get_diag_region_results <- function(model, subset_name, data) {
    em <- get_emm_log_scale(model, ~ Diagnosis * Region, data = data)
    em_df <- extract_emm_response(as.data.frame(summary(em)))
    
    means <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      Diagnosis = em_df$Diagnosis,
      Region = em_df$Region,
      mean_area = fmt(em_df$mean_area, 1),
      `lower.CI` = fmt(em_df$lower.CI, 1),
      `upper.CI` = fmt(em_df$upper.CI, 1)
    )
    
    em_by_region <- get_emm_log_scale(model, ~ Diagnosis | Region, data = data)
    contr_diag <- pairs(em_by_region, adjust = "holm", reverse = TRUE)
    contr_diag_df <- extract_contrast_results(as.data.frame(summary(contr_diag)))
    
    diag_within_region <- tibble(
      subset = subset_name,
      Region = contr_diag_df$Region,
      contrast = as.character(contr_diag_df$contrast),
      log_estimate = fmt(contr_diag_df$log_estimate, 4),
      SE = fmt(contr_diag_df$SE, 4),
      p.value = fmt(contr_diag_df$p.value, 4),
      pct_diff = fmt(contr_diag_df$pct_diff, 1),
      sig = add_stars(contr_diag_df$p.value)
    )
    
    em_by_diag <- get_emm_log_scale(model, ~ Region | Diagnosis, data = data)
    contr_region <- pairs(em_by_diag, adjust = "tukey")
    contr_region_df <- extract_contrast_results(as.data.frame(summary(contr_region)))
    region_within_diag <- tibble(
      subset = subset_name,
      Diagnosis = contr_region_df$Diagnosis,
      contrast = as.character(contr_region_df$contrast),
      log_estimate = fmt(contr_region_df$log_estimate, 4),
      SE = fmt(contr_region_df$SE, 4),
      p.value = fmt(contr_region_df$p.value, 4),
      pct_diff = fmt(contr_region_df$pct_diff, 1),
      sig = add_stars(contr_region_df$p.value)
    )
    
    list(means = means, diag_within_region = diag_within_region, region_within_diag = region_within_diag)
  }
  
  diag_region_results <- map(subset_names, ~get_diag_region_results(models[[.x]], .x, datasets[[.x]]))
  diag_region_means <- bind_rows(map(diag_region_results, "means"))
  diag_within_region <- bind_rows(map(diag_region_results, "diag_within_region"))
  region_within_diag <- bind_rows(map(diag_region_results, "region_within_diag"))
  
  diag_region_detail_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x REGION INTERACTION (detailed)"),
    tibble(Note = "Marginal means and contrasts averaged over Sex."),
    tibble(Note = NA),
    tibble(Note = "CELL MEANS (µm²)"),
    diag_region_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT WITHIN EACH REGION (ND vs T1D, Holm-adjusted)"),
    diag_within_region %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "REGION CONTRASTS WITHIN EACH DIAGNOSIS (Tukey-adjusted)"),
    region_within_diag %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 9: DIAGNOSIS x SEX
  # =========================================================================
  cat("Creating Sheet 9: Diagnosis x Sex...\n")
  
  get_diag_sex_results <- function(model, subset_name, data) {
    em <- get_emm_log_scale(model, ~ Diagnosis * Sex, data = data)
    em_df <- extract_emm_response(as.data.frame(summary(em)))
    
    means <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      Diagnosis = em_df$Diagnosis,
      Sex = em_df$Sex,
      mean_area = fmt(em_df$mean_area, 1),
      `lower.CI` = fmt(em_df$lower.CI, 1),
      `upper.CI` = fmt(em_df$upper.CI, 1)
    )
    
    em_by_sex <- get_emm_log_scale(model, ~ Diagnosis | Sex, data = data)
    contr_diag <- pairs(em_by_sex, adjust = "holm", reverse = TRUE)
    contr_diag_df <- extract_contrast_results(as.data.frame(summary(contr_diag)))
    
    diag_within_sex <- tibble(
      subset = subset_name,
      Sex = contr_diag_df$Sex,
      contrast = as.character(contr_diag_df$contrast),
      log_estimate = fmt(contr_diag_df$log_estimate, 4),
      SE = fmt(contr_diag_df$SE, 4),
      p.value = fmt(contr_diag_df$p.value, 4),
      pct_diff = fmt(contr_diag_df$pct_diff, 1),
      sig = add_stars(contr_diag_df$p.value)
    )
    
    em_by_diag <- get_emm_log_scale(model, ~ Sex | Diagnosis, data = data)
    contr_sex <- pairs(em_by_diag, adjust = "holm")
    contr_sex_df <- extract_contrast_results(as.data.frame(summary(contr_sex)))
    
    sex_within_diag <- tibble(
      subset = subset_name,
      Diagnosis = contr_sex_df$Diagnosis,
      contrast = as.character(contr_sex_df$contrast),
      log_estimate = fmt(contr_sex_df$log_estimate, 4),
      SE = fmt(contr_sex_df$SE, 4),
      p.value = fmt(contr_sex_df$p.value, 4),
      pct_diff = fmt(contr_sex_df$pct_diff, 1),
      sig = add_stars(contr_sex_df$p.value)
    )
    
    list(means = means, diag_within_sex = diag_within_sex, sex_within_diag = sex_within_diag)
  }
  
  diag_sex_results <- map(subset_names, ~get_diag_sex_results(models[[.x]], .x, datasets[[.x]]))
  diag_sex_means <- bind_rows(map(diag_sex_results, "means"))
  diag_within_sex <- bind_rows(map(diag_sex_results, "diag_within_sex"))
  sex_within_diag <- bind_rows(map(diag_sex_results, "sex_within_diag"))
  
  diag_sex_sheet <- bind_rows(
    tibble(Note = "DIAGNOSIS x SEX INTERACTION"),
    tibble(Note = "Marginal means and contrasts averaged over Region."),
    tibble(Note = NA),
    tibble(Note = "CELL MEANS (µm²)"),
    diag_sex_means %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "DIAGNOSIS EFFECT WITHIN EACH SEX (ND vs T1D, Holm-adjusted)"),
    diag_within_sex %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SEX EFFECT WITHIN EACH DIAGNOSIS (Female vs Male, Holm-adjusted)"),
    sex_within_diag %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 10: ISLET SIZE SLOPES
  # =========================================================================
  cat("Creating Sheet 10: Islet Size Slopes...\n")
  
  get_islet_slopes <- function(model, subset_name, data) {
    mt <- detect_model_type(model)
    fe <- summary(model)$coefficients$cond
    islet_row <- which(rownames(fe) == "log_Islet.Cells_c")
    
    if (length(islet_row) == 0) return(NULL)
    
    est <- fe[islet_row, "Estimate"]
    
    # For log-scale models (gamma, gaussian_log), exponentiate for % change
    # For gaussian_raw, divide by marginal mean to get approximate % change
    pct <- if (mt == "gaussian_raw") {
      emm_mean <- tryCatch({
        emm <- emmeans(model, ~1, data = data)
        summary(emm)$emmean[1]
      }, error = function(e) fe["(Intercept)", "Estimate"])
      (est / emm_mean) * 100
    } else {
      (exp(est) - 1) * 100
    }
    
    tibble(
      subset = subset_name,
      family = mt,
      predictor = "log_Islet.Cells_c",
      estimate = fmt(est, 4),
      SE = fmt(fe[islet_row, "Std. Error"], 4),
      p.value = fmt(fe[islet_row, "Pr(>|z|)"], 4),
      pct_per_log_unit = fmt(pct, 2),
      sig = add_stars(fe[islet_row, "Pr(>|z|)"])
    )
  }
  
  islet_slopes <- bind_rows(
    get_islet_slopes(multicell_model, "SEO", multicell_data),
    get_islet_slopes(islet_model, "Islet", islet_data),
    get_islet_slopes(combined_model, "Combined", combined_data)
  ) %>% compact()
  
  if (nrow(islet_slopes) > 0) {
    islet_sheet <- bind_rows(
      tibble(Note = "ISLET SIZE EFFECT"),
      tibble(Note = "Change in cell area per 1 unit increase in log(Islet.Cells)."),
      tibble(Note = "For log-scale models: coefficient is on log scale. pct_per_log_unit = (exp(estimate)-1)*100."),
      tibble(Note = "For Gaussian(raw): coefficient is in µm² per log-unit. pct_per_log_unit = (estimate/marginal_mean)*100."),
      tibble(Note = "1 unit on log scale ≈ 2.72-fold increase in cell count."),
      tibble(Note = "SEO and Islet rows are from subset-specific models."),
      tibble(Note = "Combined row is from the pooled model (singlets fixed at log_Islet.Cells_c = 0)."),
      tibble(Note = NA),
      islet_slopes %>% mutate(Note = NA) %>% select(Note, everything())
    )
  } else {
    islet_sheet <- tibble(Note = "Islet size slopes not available.")
  }
  
  # =========================================================================
  # SHEET 11: COMBINED MODEL SUMMARY
  # =========================================================================
  cat("Creating Sheet 11: Combined Model Summary...\n")
  
  fe_combined <- summary(combined_model)$coefficients$cond
  ci_combined <- tryCatch({
    confint(combined_model, method = "Wald", parm = "beta_")
  }, error = function(e) {
    est <- fe_combined[, "Estimate"]
    se <- fe_combined[, "Std. Error"]
    cbind(est - 1.96 * se, est + 1.96 * se)
  })
  
  combined_fe <- tibble(
    term = rownames(fe_combined),
    estimate = fmt(fe_combined[, "Estimate"], 4),
    std.error = fmt(fe_combined[, "Std. Error"], 4),
    `l-95% CI` = fmt(ci_combined[rownames(fe_combined), 1], 4),
    `u-95% CI` = fmt(ci_combined[rownames(fe_combined), 2], 4),
    z.value = fmt(fe_combined[, "z value"], 2),
    p.value = fmt(fe_combined[, "Pr(>|z|)"], 4),
    pct_change = fmt((exp(fe_combined[, "Estimate"]) - 1) * 100, 1),
    sig = add_stars(fe_combined[, "Pr(>|z|)"])
  )
  
  combined_sheet <- bind_rows(
    tibble(Note = "COMBINED MODEL (all object types)"),
    tibble(Note = "Model: area ~ (Diagnosis + Region + Sex + Age_c)^2 +"),
    tibble(Note = "       Diagnosis:Region:Age_c + log_Islet.Cells_c +"),
    tibble(Note = "       object_type * Diagnosis + object_type * Region + (1|Donor) + (1|Donor:ImageID)"),
    tibble(Note = "Family: Gamma(log link)"),
    tibble(Note = "log_Islet.Cells_c included as main-effect-only adjustment (no interactions, to avoid"),
    tibble(Note = "collinearity with object_type). Singlets set to 0 (centered grand mean)."),
    tibble(Note = ""),
    tibble(Note = paste("N =", nrow(combined_data))),
    tibble(Note = NA),
    combined_fe %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # SHEET 12: RANDOM EFFECTS & SAMPLE SIZES
  # =========================================================================
  cat("Creating Sheet 12: Random Effects & Sample Sizes...\n")
  
  extract_re <- function(model, subset_name) {
    vc <- VarCorr(model)$cond
    
    re_df <- tibble(
      subset = subset_name,
      family = detect_model_type(model),
      group = names(vc),
      variance = fmt(sapply(vc, function(x) x[1]), 4),
      std_dev = fmt(sqrt(sapply(vc, function(x) x[1])), 4)
    )
    
    bind_rows(re_df, tibble(
      subset = subset_name,
      family = detect_model_type(model),
      group = "Residual",
      variance = fmt(sigma(model)^2, 4),
      std_dev = fmt(sigma(model), 4)
    ))
  }
  
  random_effects <- bind_rows(
    extract_re(singlet_model, "Singlet"),
    extract_re(multicell_model, "SEO"),
    extract_re(islet_model, "Islet"),
    extract_re(combined_model, "Combined")
  )
  
  sample_sizes <- tibble(
    subset = c("Singlet", "SEO", "Islet", "Combined"),
    n_obs = c(nrow(singlet_data), nrow(multicell_data), nrow(islet_data), nrow(combined_data)),
    n_donors = c(
      length(unique(singlet_data$Donor)),
      length(unique(multicell_data$Donor)),
      length(unique(islet_data$Donor)),
      length(unique(combined_data$Donor))
    )
  )
  
  re_sheet <- bind_rows(
    tibble(Note = "RANDOM EFFECTS VARIANCE COMPONENTS"),
    tibble(Note = NA),
    random_effects %>% mutate(Note = NA) %>% select(Note, everything()),
    tibble(Note = NA),
    tibble(Note = "SAMPLE SIZES"),
    sample_sizes %>% mutate(Note = NA) %>% select(Note, everything())
  )
  
  # =========================================================================
  # ASSEMBLE AND SAVE
  # =========================================================================
  cat("\nAssembling workbook...\n")
  
  sheets <- list(
    "Model Summary" = as.data.frame(model_summary),
    "Object Type" = as.data.frame(object_type_sheet),
    "Object Type x Diagnosis" = as.data.frame(obj_diag_sheet),
    "Diag x Region Summary" = as.data.frame(diag_region_sheet),
    "Diagnosis" = as.data.frame(diagnosis_sheet),
    "Region" = as.data.frame(region_sheet),
    "Sex" = as.data.frame(sex_sheet),
    "Diagnosis x Region" = as.data.frame(diag_region_detail_sheet),
    "Diagnosis x Sex" = as.data.frame(diag_sex_sheet),
    "Islet Size Slopes" = as.data.frame(islet_sheet),
    "Combined Model" = as.data.frame(combined_sheet),
    "Random Effects" = as.data.frame(re_sheet)
  )
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0(cell_label, "_Cell_Area_Analysis.xlsx"))
  
  write_xlsx(sheets, output_file)
  
  cat("\nSaved to:", output_file, "\n")
  
  return(invisible(sheets))
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

generate_all_area_excel <- function(output_dir = "Area/Results/Total_Area") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING ALL CELL AREA EXCEL RESULTS (DIAGNOSIS)\n")
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
      Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
      Region    = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor),
      ImageID   = factor(ImageID))
  }
  
  add_object_type <- function(df) {
    df %>% mutate(
      object_type = factor(case_when(
        Islet.Cells == 1  ~ "Singlets",
        Islet.Cells >= 2 & Islet.Cells <= 14 ~ "SEOs",
        Islet.Cells >= 15 ~ "Islets"),
        levels = c("Singlets", "SEOs", "Islets")))
  }
  
  # --- Load subset-specific datasets (matching model training data) ---
  cat("Loading subset-specific data...\n")
  area_singlets <- set_factors(readRDS("Data/area_data_singlets.rds"))
  area_endobs   <- set_factors(readRDS("Data/area_data_endobs.rds"))
  area_islets   <- set_factors(readRDS("Data/area_data_islets.rds"))
  
  beta_combined  <- add_object_type(set_factors(readRDS("Data/area_data_beta_combined.rds")))
  alpha_combined <- add_object_type(set_factors(readRDS("Data/area_data_alpha_combined.rds")))
  delta_combined <- add_object_type(set_factors(readRDS("Data/area_data_delta_combined.rds")))
  pp_combined    <- add_object_type(set_factors(readRDS("Data/area_data_pp_combined.rds")))
  cat("  Data loaded.\n\n")
  
  combined_data_map <- list(
    ins  = beta_combined,
    glu  = alpha_combined,
    soma = delta_combined,
    pp   = pp_combined
  )
  
  # Cell type definitions
  cell_defs <- list(
    list(type = "ins",  label = "β",  area_col = "ins_area",
         prefix_1 = "area_beta_1cell", prefix_2 = "area_beta_2to14",
         prefix_15 = "area_beta_15plus", prefix_comb = "beta_combined_model"),
    list(type = "glu",  label = "α", area_col = "glu_area",
         prefix_1 = "area_alpha_1cell", prefix_2 = "area_alpha_2to14",
         prefix_15 = "area_alpha_15plus", prefix_comb = "alpha_combined_model"),
    list(type = "soma", label = "δ", area_col = "soma_area",
         prefix_1 = "area_delta_1cell", prefix_2 = "area_delta_2to14",
         prefix_15 = "area_delta_15plus", prefix_comb = "delta_combined_model"),
    list(type = "pp",   label = "PP",    area_col = "pp_area",
         prefix_1 = "area_pp_1cell", prefix_2 = "area_pp_2to14",
         prefix_15 = "area_pp_15plus", prefix_comb = "pp_combined_model")
  )
  
  for (cd in cell_defs) {
    cat("\nProcessing", cd$label, "cell...\n")
    
    # Load models
    m1  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_1, ".rds"))
    m2  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_2, ".rds"))
    m15 <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_15, ".rds"))
    mc  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_comb, ".rds"))
    
    # Prepare data subsets from correct source datasets (matching models)
    d1  <- split_area(area_singlets, cd$area_col)
    d2  <- split_area(area_endobs, cd$area_col)
    d15 <- split_area(area_islets, cd$area_col)
    dc  <- split_area(combined_data_map[[cd$type]], cd$area_col)
    
    generate_cell_type_excel(
      singlet_model = m1, multicell_model = m2,
      islet_model = m15, combined_model = mc,
      singlet_data = d1, multicell_data = d2,
      islet_data = d15, combined_data = dc,
      cell_type = cd$type, cell_label = cd$label,
      area_col = cd$area_col, output_dir = output_dir
    )
    rm(m1, m2, m15, mc); gc()
  }
  
  cat("\n=============================================================================\n")
  cat("COMPLETE! Four Excel files created in", output_dir, "\n")
  cat("=============================================================================\n")
}

# Usage
if (interactive()) {
  cat("\nTo generate all Excel files, run:\n")
  cat("  generate_all_area_excel()\n")
}
