# =============================================================================
# Generate Comprehensive Excel Results for Proportion T1D Models (ordbetareg)
# Age at Onset parameterization (DD + AO slopes from single model)
# =============================================================================
#
# Models (4 total): ins_ao, glu_ao, soma_ao, pp_ao
#
# All models: ordbetareg (Bayesian ordered beta regression)
#   Formula: Percent.X ~ (Disease.Duration_c + Region + Sex + age_at_onset_c + log_Islet.Cells_c)^2 +
#            Disease.Duration_c:Region:log_Islet.Cells_c +
#            Disease.Duration_c:Region:age_at_onset_c + (1|Donor) + (1|Donor:ImageID)
#
# Note: The DD + Age_c parameterization is algebraically equivalent (same fit).
# The AO model is preferred because both DD and AO slopes have direct clinical
# interpretations.  DD slopes are extracted from this model directly.
#
# Data: quad_t1d (T1D only)
#
# Output: 4 Excel files in Proportion/Results/
#
# =============================================================================

library(tidyverse)
library(brms)
library(coda)
library(openxlsx)

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

compute_contrast <- function(draws1, draws2, label) {
  diff <- draws1 - draws2
  summ <- summarize_posterior(diff)
  summ$contrast <- label
  summ$prob_greater_0 <- mean(diff > 0)
  summ$prob_less_0 <- mean(diff < 0)
  return(summ)
}

fmt <- function(x, digits = 4) round(x, digits)


# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================
#
# Handles both DD and AO models via the `model_type` parameter.
#
# model_type = "dd" → continuous predictors: Disease.Duration_c, Age_c, log_Islet.Cells_c
#   3-way interactions: DD:Region:Age_c, DD:Region:log_Islet.Cells_c
#
# model_type = "ao" → continuous predictors: Disease.Duration_c, age_at_onset_c, log_Islet.Cells_c
#   3-way interactions: DD:Region:age_at_onset_c, DD:Region:log_Islet.Cells_c
#
# Sheets:
#   1. Model Summary
#   2. Region
#   3. Sex
#   4. Region × Sex
#   5. Disease Duration Slopes
#   6. Age / Age at Onset Slopes
#   7. Islet Size Slopes
#   8. DD × Region × Age/AO (3-way)
#   9. DD × Region × Islet (3-way)
#  10. Supplementary Pairwise (Region × Sex cells)
# =============================================================================

generate_proportion_t1d_excel <- function(
    model, data, cell_type, cell_label, model_type = "dd",
    output_dir = "Proportion/Results"
) {

  # Determine the "second" continuous predictor name
  age_var <- if (model_type == "dd") "Age_c" else "age_at_onset_c"
  age_label <- if (model_type == "dd") "Age" else "Age at Onset"
  type_label <- if (model_type == "dd") "Disease Duration" else "Age at Onset"
  type_tag   <- toupper(model_type)

  cat("\n")
  cat("=============================================================================\n")
  cat("GENERATING", type_tag, "ANALYSIS EXCEL:", toupper(cell_label), "\n")
  cat("=============================================================================\n\n")

  # -------------------------------------------------------------------------
  # MODEL SUMMARY
  # -------------------------------------------------------------------------
  cat("Creating Sheet 1: Model Summary...\n")

  model_summary_raw <- summary(model)$fixed
  model_summary <- data.frame(
    term = rownames(model_summary_raw),
    estimate = model_summary_raw$Estimate,
    Est.Error = model_summary_raw$Est.Error,
    l.95.CI = model_summary_raw$`l-95% CI`,
    u.95.CI = model_summary_raw$`u-95% CI`,
    Rhat = model_summary_raw$Rhat,
    Bulk_ESS = model_summary_raw$Bulk_ESS,
    Tail_ESS = model_summary_raw$Tail_ESS
  )

  # -------------------------------------------------------------------------
  # CATEGORICAL GRID (Region × Sex)
  # All continuous predictors at their centered mean (0)
  # -------------------------------------------------------------------------
  cat("Setting up prediction grid and getting posterior draws...\n")

  cat_grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Disease.Duration_c = 0,
    log_Islet.Cells_c = 0
  )
  cat_grid[[age_var]] <- 0

  post_fitted <- get_posterior_fitted(model, cat_grid)
  col_names <- paste(cat_grid$Region, cat_grid$Sex, sep = "_")
  colnames(post_fitted) <- col_names

  # -------------------------------------------------------------------------
  # SHEET 2: REGION MAIN EFFECT
  # -------------------------------------------------------------------------
  cat("Creating Sheet 2: Region...\n")

  head_idx <- which(cat_grid$Region == "Head")
  body_idx <- which(cat_grid$Region == "Body")
  tail_idx <- which(cat_grid$Region == "Tail")

  head_draws <- rowMeans(post_fitted[, head_idx])
  body_draws <- rowMeans(post_fitted[, body_idx])
  tail_draws <- rowMeans(post_fitted[, tail_idx])

  region_marginal <- bind_rows(
    summarize_posterior(head_draws) %>% mutate(Region = "Head", .before = 1),
    summarize_posterior(body_draws) %>% mutate(Region = "Body", .before = 1),
    summarize_posterior(tail_draws) %>% mutate(Region = "Tail", .before = 1)
  ) %>% rename(response = estimate)

  region_contrasts <- bind_rows(
    compute_contrast(head_draws, body_draws, "Head - Body"),
    compute_contrast(head_draws, tail_draws, "Head - Tail"),
    compute_contrast(body_draws, tail_draws, "Body - Tail")
  )

  # -------------------------------------------------------------------------
  # SHEET 3: SEX MAIN EFFECT
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: Sex...\n")

  female_idx <- which(cat_grid$Sex == "Female")
  male_idx   <- which(cat_grid$Sex == "Male")

  female_draws <- rowMeans(post_fitted[, female_idx])
  male_draws   <- rowMeans(post_fitted[, male_idx])

  sex_marginal <- bind_rows(
    summarize_posterior(female_draws) %>% mutate(Sex = "Female", .before = 1),
    summarize_posterior(male_draws) %>% mutate(Sex = "Male", .before = 1)
  ) %>% rename(response = estimate)

  sex_contrast <- compute_contrast(female_draws, male_draws, "Female - Male")

  # -------------------------------------------------------------------------
  # SHEET 4: REGION × SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 4: Region x Sex...\n")

  region_sex_marginal <- list()
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      idx <- which(cat_grid$Region == reg & cat_grid$Sex == s)
      draws <- rowMeans(post_fitted[, idx, drop = FALSE])
      summ <- summarize_posterior(draws)
      summ$Region <- reg
      summ$Sex <- s
      region_sex_marginal[[paste(reg, s, sep = "_")]] <- summ
    }
  }
  region_sex_marginal_df <- bind_rows(region_sex_marginal) %>%
    select(Region, Sex, estimate, lower.HPD, upper.HPD) %>%
    rename(response = estimate)

  # Sex within each Region
  sex_within_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    f_idx <- which(cat_grid$Region == reg & cat_grid$Sex == "Female")
    m_idx <- which(cat_grid$Region == reg & cat_grid$Sex == "Male")
    f_draws <- rowMeans(post_fitted[, f_idx, drop = FALSE])
    m_draws <- rowMeans(post_fitted[, m_idx, drop = FALSE])
    contrast <- compute_contrast(f_draws, m_draws, "Female - Male")
    contrast$Region <- reg
    sex_within_region[[reg]] <- contrast
  }
  sex_within_region_df <- bind_rows(sex_within_region) %>%
    select(contrast, Region, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # Region within each Sex
  region_within_sex <- list()
  for (s in c("Female", "Male")) {
    h_idx <- which(cat_grid$Region == "Head" & cat_grid$Sex == s)
    b_idx <- which(cat_grid$Region == "Body" & cat_grid$Sex == s)
    t_idx <- which(cat_grid$Region == "Tail" & cat_grid$Sex == s)

    h_draws <- rowMeans(post_fitted[, h_idx, drop = FALSE])
    b_draws <- rowMeans(post_fitted[, b_idx, drop = FALSE])
    t_draws <- rowMeans(post_fitted[, t_idx, drop = FALSE])

    hb <- compute_contrast(h_draws, b_draws, "Head - Body"); hb$Sex <- s
    ht <- compute_contrast(h_draws, t_draws, "Head - Tail"); ht$Sex <- s
    bt <- compute_contrast(b_draws, t_draws, "Body - Tail"); bt$Sex <- s

    region_within_sex[[paste0(s, "_HB")]] <- hb
    region_within_sex[[paste0(s, "_HT")]] <- ht
    region_within_sex[[paste0(s, "_BT")]] <- bt
  }
  region_within_sex_df <- bind_rows(region_within_sex) %>%
    select(contrast, Sex, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # CONTINUOUS PREDICTOR SLOPES SETUP
  # Finite-difference approach: predict at ±1, slope = (high - low) / 2
  # -------------------------------------------------------------------------
  cat("Setting up continuous predictor grids for slopes...\n")

  # --- Disease Duration slopes ---
  dd_grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Disease.Duration_c = c(-1, 1),
    log_Islet.Cells_c = 0
  )
  dd_grid[[age_var]] <- 0

  # --- Age / Age at Onset slopes ---
  age_grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Disease.Duration_c = 0,
    log_Islet.Cells_c = 0
  )
  age_grid_low  <- age_grid; age_grid_low[[age_var]]  <- -1
  age_grid_high <- age_grid; age_grid_high[[age_var]] <-  1
  age_grid_full <- bind_rows(age_grid_low, age_grid_high)

  # --- Islet Size slopes ---
  islet_grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Disease.Duration_c = 0,
    log_Islet.Cells_c = c(-1, 1)
  )
  islet_grid[[age_var]] <- 0

  cat("Getting posterior predictions for slopes...\n")
  post_dd    <- get_posterior_fitted(model, dd_grid)
  post_age   <- get_posterior_fitted(model, age_grid_full)
  post_islet <- get_posterior_fitted(model, islet_grid)

  # Compute cell-level slopes for each Region × Sex combination
  dd_slopes_list    <- list()
  age_slopes_list   <- list()
  islet_slopes_list <- list()

  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      key <- paste(reg, s, sep = "_")

      # DD slopes
      low_idx  <- which(dd_grid$Region == reg & dd_grid$Sex == s & dd_grid$Disease.Duration_c == -1)
      high_idx <- which(dd_grid$Region == reg & dd_grid$Sex == s & dd_grid$Disease.Duration_c == 1)
      dd_slopes_list[[key]] <- (post_dd[, high_idx] - post_dd[, low_idx]) / 2

      # Age / AO slopes
      n_cells <- nrow(age_grid)  # 6 cells
      low_idx  <- which(age_grid_low$Region == reg  & age_grid_low$Sex == s)
      high_idx <- which(age_grid_high$Region == reg & age_grid_high$Sex == s) + n_cells
      age_slopes_list[[key]] <- (post_age[, high_idx] - post_age[, low_idx]) / 2

      # Islet slopes
      low_idx  <- which(islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == -1)
      high_idx <- which(islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == 1)
      islet_slopes_list[[key]] <- (post_islet[, high_idx] - post_islet[, low_idx]) / 2
    }
  }

  # -------------------------------------------------------------------------
  # Generic slope summarizer
  # Given a slopes_list (keyed by "Region_Sex"), marginalizes over groups
  # -------------------------------------------------------------------------
  summarize_slopes <- function(slopes_list, group_var, group_levels, marginal_over) {
    results_list <- list()
    draws_list <- list()
    for (g in group_levels) {
      subset_keys <- names(slopes_list)
      if (group_var == "Region") {
        subset_keys <- subset_keys[grepl(paste0("^", g, "_"), subset_keys)]
      } else if (group_var == "Sex") {
        subset_keys <- subset_keys[grepl(paste0("_", g, "$"), subset_keys)]
      } else {
        # overall — use all keys
        subset_keys <- names(slopes_list)
      }
      slope_draws <- rowMeans(do.call(cbind, slopes_list[subset_keys]))
      summ <- summarize_posterior(slope_draws)
      summ[[group_var]] <- g
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      results_list[[g]] <- summ
      draws_list[[g]] <- slope_draws
    }
    list(summary = bind_rows(results_list), draws = draws_list)
  }

  # -------------------------------------------------------------------------
  # SHEET 5: DISEASE DURATION SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 5: Disease Duration Slopes...\n")

  # Overall
  all_dd_slopes <- rowMeans(do.call(cbind, dd_slopes_list))
  overall_dd <- summarize_posterior(all_dd_slopes)
  overall_dd$predictor <- "Disease.Duration_c"
  overall_dd$prob_positive <- mean(all_dd_slopes > 0)
  overall_dd$prob_negative <- mean(all_dd_slopes < 0)

  # By Region
  dd_by_region <- summarize_slopes(dd_slopes_list, "Region", c("Head", "Body", "Tail"), "Sex")
  dd_by_region_df <- dd_by_region$summary %>%
    select(Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  dd_region_contrasts_df <- bind_rows(
    compute_contrast(dd_by_region$draws[["Head"]], dd_by_region$draws[["Body"]], "Head - Body"),
    compute_contrast(dd_by_region$draws[["Head"]], dd_by_region$draws[["Tail"]], "Head - Tail"),
    compute_contrast(dd_by_region$draws[["Body"]], dd_by_region$draws[["Tail"]], "Body - Tail")
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # By Sex
  dd_by_sex <- summarize_slopes(dd_slopes_list, "Sex", c("Female", "Male"), "Region")
  dd_by_sex_df <- dd_by_sex$summary %>%
    select(Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  dd_sex_contrast_df <- compute_contrast(
    dd_by_sex$draws[["Female"]], dd_by_sex$draws[["Male"]], "Female - Male"
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 6: AGE / AGE AT ONSET SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 6:", age_label, "Slopes...\n")

  # Overall
  all_age_slopes <- rowMeans(do.call(cbind, age_slopes_list))
  overall_age <- summarize_posterior(all_age_slopes)
  overall_age$predictor <- age_var
  overall_age$prob_positive <- mean(all_age_slopes > 0)
  overall_age$prob_negative <- mean(all_age_slopes < 0)

  # By Region
  age_by_region <- summarize_slopes(age_slopes_list, "Region", c("Head", "Body", "Tail"), "Sex")
  age_by_region_df <- age_by_region$summary %>%
    select(Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  age_region_contrasts_df <- bind_rows(
    compute_contrast(age_by_region$draws[["Head"]], age_by_region$draws[["Body"]], "Head - Body"),
    compute_contrast(age_by_region$draws[["Head"]], age_by_region$draws[["Tail"]], "Head - Tail"),
    compute_contrast(age_by_region$draws[["Body"]], age_by_region$draws[["Tail"]], "Body - Tail")
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # By Sex
  age_by_sex <- summarize_slopes(age_slopes_list, "Sex", c("Female", "Male"), "Region")
  age_by_sex_df <- age_by_sex$summary %>%
    select(Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  age_sex_contrast_df <- compute_contrast(
    age_by_sex$draws[["Female"]], age_by_sex$draws[["Male"]], "Female - Male"
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 7: ISLET SIZE SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 7: Islet Size Slopes...\n")

  # Overall
  all_islet_slopes <- rowMeans(do.call(cbind, islet_slopes_list))
  overall_islet <- summarize_posterior(all_islet_slopes)
  overall_islet$predictor <- "log_Islet.Cells_c"
  overall_islet$prob_positive <- mean(all_islet_slopes > 0)
  overall_islet$prob_negative <- mean(all_islet_slopes < 0)

  # By Region
  islet_by_region <- summarize_slopes(islet_slopes_list, "Region", c("Head", "Body", "Tail"), "Sex")
  islet_by_region_df <- islet_by_region$summary %>%
    select(Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  islet_region_contrasts_df <- bind_rows(
    compute_contrast(islet_by_region$draws[["Head"]], islet_by_region$draws[["Body"]], "Head - Body"),
    compute_contrast(islet_by_region$draws[["Head"]], islet_by_region$draws[["Tail"]], "Head - Tail"),
    compute_contrast(islet_by_region$draws[["Body"]], islet_by_region$draws[["Tail"]], "Body - Tail")
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # By Sex
  islet_by_sex <- summarize_slopes(islet_slopes_list, "Sex", c("Female", "Male"), "Region")
  islet_by_sex_df <- islet_by_sex$summary %>%
    select(Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  islet_sex_contrast_df <- compute_contrast(
    islet_by_sex$draws[["Female"]], islet_by_sex$draws[["Male"]], "Female - Male"
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 8: DD × REGION × AGE/AO (3-way)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 8: DD x Region x", age_label, "(3-way)...\n")

  # DD slopes for each Region, evaluated at ±1 of the age variable
  # This captures the 3-way: how the DD slope changes with Age/AO across Regions
  dd_age_grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Disease.Duration_c = c(-1, 1),
    log_Islet.Cells_c = 0
  )

  # Low age/AO
  dd_age_grid_low <- dd_age_grid
  dd_age_grid_low[[age_var]] <- -1
  # High age/AO
  dd_age_grid_high <- dd_age_grid
  dd_age_grid_high[[age_var]] <- 1

  post_dd_age_low  <- get_posterior_fitted(model, dd_age_grid_low)
  post_dd_age_high <- get_posterior_fitted(model, dd_age_grid_high)

  # DD slopes at low vs high age/AO for each Region (marginalized over Sex)
  dd_slope_at_age <- list()
  dd_slope_at_age_draws <- list()
  for (age_level in c("low", "high")) {
    post_mat <- if (age_level == "low") post_dd_age_low else post_dd_age_high
    grid_ref <- dd_age_grid  # same structure for indexing
    for (reg in c("Head", "Body", "Tail")) {
      slopes_for_sex <- list()
      for (s in c("Female", "Male")) {
        low_idx  <- which(grid_ref$Region == reg & grid_ref$Sex == s & grid_ref$Disease.Duration_c == -1)
        high_idx <- which(grid_ref$Region == reg & grid_ref$Sex == s & grid_ref$Disease.Duration_c == 1)
        slopes_for_sex[[s]] <- (post_mat[, high_idx] - post_mat[, low_idx]) / 2
      }
      slope_draws <- rowMeans(do.call(cbind, slopes_for_sex))
      key <- paste(reg, age_level, sep = "_")
      summ <- summarize_posterior(slope_draws)
      summ$Region <- reg
      summ$age_level <- paste0(age_label, " = mean ", ifelse(age_level == "low", "- 1", "+ 1"))
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      dd_slope_at_age[[key]] <- summ
      dd_slope_at_age_draws[[key]] <- slope_draws
    }
  }
  dd_slope_at_age_df <- bind_rows(dd_slope_at_age) %>%
    select(Region, age_level, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

  # Difference in DD slope: high age - low age, within each Region
  dd_age_interaction_by_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff_draws <- dd_slope_at_age_draws[[paste(reg, "high", sep = "_")]] -
      dd_slope_at_age_draws[[paste(reg, "low", sep = "_")]]
    summ <- summarize_posterior(diff_draws)
    summ$Region <- reg
    summ$contrast <- paste0("DD slope at high ", age_label, " - low ", age_label)
    summ$prob_greater_0 <- mean(diff_draws > 0)
    summ$prob_less_0 <- mean(diff_draws < 0)
    dd_age_interaction_by_region[[reg]] <- summ
  }
  dd_age_interaction_df <- bind_rows(dd_age_interaction_by_region) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # Region contrasts on DD slopes (at mean age), from Sheet 5
  # Already computed in dd_by_region above — reuse

  # -------------------------------------------------------------------------
  # SHEET 9: DD × REGION × ISLET (3-way)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 9: DD x Region x Islet (3-way)...\n")

  dd_islet_grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Disease.Duration_c = c(-1, 1)
  )
  dd_islet_grid[[age_var]] <- 0

  # Low islet
  dd_islet_grid_low <- dd_islet_grid
  dd_islet_grid_low$log_Islet.Cells_c <- -1
  # High islet
  dd_islet_grid_high <- dd_islet_grid
  dd_islet_grid_high$log_Islet.Cells_c <- 1

  post_dd_islet_low  <- get_posterior_fitted(model, dd_islet_grid_low)
  post_dd_islet_high <- get_posterior_fitted(model, dd_islet_grid_high)

  dd_slope_at_islet <- list()
  dd_slope_at_islet_draws <- list()
  for (islet_level in c("low", "high")) {
    post_mat <- if (islet_level == "low") post_dd_islet_low else post_dd_islet_high
    grid_ref <- dd_islet_grid
    for (reg in c("Head", "Body", "Tail")) {
      slopes_for_sex <- list()
      for (s in c("Female", "Male")) {
        low_idx  <- which(grid_ref$Region == reg & grid_ref$Sex == s & grid_ref$Disease.Duration_c == -1)
        high_idx <- which(grid_ref$Region == reg & grid_ref$Sex == s & grid_ref$Disease.Duration_c == 1)
        slopes_for_sex[[s]] <- (post_mat[, high_idx] - post_mat[, low_idx]) / 2
      }
      slope_draws <- rowMeans(do.call(cbind, slopes_for_sex))
      key <- paste(reg, islet_level, sep = "_")
      summ <- summarize_posterior(slope_draws)
      summ$Region <- reg
      summ$islet_level <- paste0("Islet size = mean ", ifelse(islet_level == "low", "- 1", "+ 1"))
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      dd_slope_at_islet[[key]] <- summ
      dd_slope_at_islet_draws[[key]] <- slope_draws
    }
  }
  dd_slope_at_islet_df <- bind_rows(dd_slope_at_islet) %>%
    select(Region, islet_level, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

  # Difference: high islet - low islet, within each Region
  dd_islet_interaction_by_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff_draws <- dd_slope_at_islet_draws[[paste(reg, "high", sep = "_")]] -
      dd_slope_at_islet_draws[[paste(reg, "low", sep = "_")]]
    summ <- summarize_posterior(diff_draws)
    summ$Region <- reg
    summ$contrast <- "DD slope at high Islet - low Islet"
    summ$prob_greater_0 <- mean(diff_draws > 0)
    summ$prob_less_0 <- mean(diff_draws < 0)
    dd_islet_interaction_by_region[[reg]] <- summ
  }
  dd_islet_interaction_df <- bind_rows(dd_islet_interaction_by_region) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 10: SUPPLEMENTARY PAIRWISE (Region × Sex cells)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 10: Supplementary Pairwise...\n")

  cell_combos <- expand.grid(
    r1 = c("Head", "Body", "Tail"),
    s1 = c("Female", "Male"),
    r2 = c("Head", "Body", "Tail"),
    s2 = c("Female", "Male"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(c1 = paste(r1, s1), c2 = paste(r2, s2)) %>%
    filter(c1 < c2) %>%
    select(r1, s1, r2, s2)

  full_pairwise <- list()
  for (i in 1:nrow(cell_combos)) {
    row <- cell_combos[i, ]
    idx1 <- which(cat_grid$Region == row$r1 & cat_grid$Sex == row$s1)
    idx2 <- which(cat_grid$Region == row$r2 & cat_grid$Sex == row$s2)
    draws1 <- rowMeans(post_fitted[, idx1, drop = FALSE])
    draws2 <- rowMeans(post_fitted[, idx2, drop = FALSE])
    label <- paste0(row$r1, " ", row$s1, " - ", row$r2, " ", row$s2)
    full_pairwise[[label]] <- compute_contrast(draws1, draws2, label)
  }
  full_pairwise_df <- bind_rows(full_pairwise) %>%
    select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # CREATE EXCEL WORKBOOK
  # -------------------------------------------------------------------------
  cat("\nCreating Excel workbook...\n")

  wb <- createWorkbook()

  # Sheet 1: Model Summary
  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", data.frame(
    Note = c(
      paste(toupper(type_label), "ANALYSIS:", cell_label, "Content (T1D only)"),
      "Bayesian ordered beta regression (ordbetareg)",
      paste0("Response: Percent.", cell_type, " (0-1 scale)"),
      paste0("Model: (Disease.Duration_c + Region + Sex + ", age_var, " + log_Islet.Cells_c)^2 +"),
      paste0("       Disease.Duration_c:Region:log_Islet.Cells_c + Disease.Duration_c:Region:", age_var),
      "       + (1|Donor) + (1|Donor:ImageID)",
      "",
      "All estimates are posterior medians with 95% Highest Posterior Density (HPD) intervals.",
      "prob_greater_0 / prob_less_0 are posterior probabilities that the contrast is positive/negative.",
      "Slopes computed via finite differences at ±1 of the centered predictor.",
      "")
  ))
  writeData(wb, "Model Summary", model_summary, startRow = 13)

  # Sheet 2: Region
  addWorksheet(wb, "Region")
  writeData(wb, "Region", data.frame(
    Note = paste("Main effect of Region on predicted percent", tolower(cell_label), "content, averaged over Sex.")))
  writeData(wb, "Region", region_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region", data.frame(x = ""), startRow = 7)
  writeData(wb, "Region", data.frame(Section = "Pairwise contrasts:"), startRow = 8)
  writeData(wb, "Region", region_contrasts %>%
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 9)

  # Sheet 3: Sex
  addWorksheet(wb, "Sex")
  writeData(wb, "Sex", data.frame(
    Note = paste("Main effect of Sex on predicted percent", tolower(cell_label), "content, averaged over Region.")))
  writeData(wb, "Sex", sex_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Sex", data.frame(x = ""), startRow = 6)
  writeData(wb, "Sex", data.frame(Section = "Contrast:"), startRow = 7)
  writeData(wb, "Sex", sex_contrast %>%
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 8)

  # Sheet 4: Region x Sex
  addWorksheet(wb, "Region x Sex")
  writeData(wb, "Region x Sex", data.frame(
    Note = "Region x Sex interaction. Marginal means and contrasts."))
  writeData(wb, "Region x Sex", region_sex_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region x Sex", data.frame(x = ""), startRow = 10)
  writeData(wb, "Region x Sex", data.frame(Section = "Sex effect (Female - Male) within each Region:"), startRow = 11)
  writeData(wb, "Region x Sex", sex_within_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 12)
  writeData(wb, "Region x Sex", data.frame(x = ""), startRow = 16)
  writeData(wb, "Region x Sex", data.frame(Section = "Region contrasts within each Sex:"), startRow = 17)
  writeData(wb, "Region x Sex", region_within_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 18)

  # Sheet 5: Disease Duration Slopes
  addWorksheet(wb, "DD Slopes")
  writeData(wb, "DD Slopes", data.frame(
    Note = paste("Disease Duration slopes: change in predicted percent", tolower(cell_label),
                 "per 1 year increase in Disease Duration (centered).")))
  writeData(wb, "DD Slopes", data.frame(Section = "Overall DD slope (averaged over Region and Sex):"), startRow = 3)
  writeData(wb, "DD Slopes", overall_dd %>%
              select(predictor, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 4)
  writeData(wb, "DD Slopes", data.frame(x = ""), startRow = 6)
  writeData(wb, "DD Slopes", data.frame(Section = "DD slopes by Region:"), startRow = 7)
  writeData(wb, "DD Slopes", dd_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 8)
  writeData(wb, "DD Slopes", data.frame(x = ""), startRow = 12)
  writeData(wb, "DD Slopes", data.frame(Section = "Region contrasts (DD slopes):"), startRow = 13)
  writeData(wb, "DD Slopes", dd_region_contrasts_df %>% mutate(across(where(is.numeric), fmt)), startRow = 14)
  writeData(wb, "DD Slopes", data.frame(x = ""), startRow = 18)
  writeData(wb, "DD Slopes", data.frame(Section = "DD slopes by Sex:"), startRow = 19)
  writeData(wb, "DD Slopes", dd_by_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 20)
  writeData(wb, "DD Slopes", data.frame(x = ""), startRow = 23)
  writeData(wb, "DD Slopes", data.frame(Section = "Sex contrast (DD slopes):"), startRow = 24)
  writeData(wb, "DD Slopes", dd_sex_contrast_df %>% mutate(across(where(is.numeric), fmt)), startRow = 25)

  # Sheet 6: Age / Age at Onset Slopes
  age_sheet_name <- paste0(age_label, " Slopes")
  addWorksheet(wb, age_sheet_name)
  writeData(wb, age_sheet_name, data.frame(
    Note = paste(age_label, "slopes: change in predicted percent", tolower(cell_label),
                 "per 1 year increase in", age_label, "(centered).")))
  writeData(wb, age_sheet_name, data.frame(Section = paste0("Overall ", age_label, " slope (averaged over Region and Sex):")), startRow = 3)
  writeData(wb, age_sheet_name, overall_age %>%
              select(predictor, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 4)
  writeData(wb, age_sheet_name, data.frame(x = ""), startRow = 6)
  writeData(wb, age_sheet_name, data.frame(Section = paste(age_label, "slopes by Region:")), startRow = 7)
  writeData(wb, age_sheet_name, age_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 8)
  writeData(wb, age_sheet_name, data.frame(x = ""), startRow = 12)
  writeData(wb, age_sheet_name, data.frame(Section = paste("Region contrasts (", age_label, "slopes):")), startRow = 13)
  writeData(wb, age_sheet_name, age_region_contrasts_df %>% mutate(across(where(is.numeric), fmt)), startRow = 14)
  writeData(wb, age_sheet_name, data.frame(x = ""), startRow = 18)
  writeData(wb, age_sheet_name, data.frame(Section = paste(age_label, "slopes by Sex:")), startRow = 19)
  writeData(wb, age_sheet_name, age_by_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 20)
  writeData(wb, age_sheet_name, data.frame(x = ""), startRow = 23)
  writeData(wb, age_sheet_name, data.frame(Section = paste("Sex contrast (", age_label, "slopes):")), startRow = 24)
  writeData(wb, age_sheet_name, age_sex_contrast_df %>% mutate(across(where(is.numeric), fmt)), startRow = 25)

  # Sheet 7: Islet Size Slopes
  addWorksheet(wb, "Islet Size Slopes")
  writeData(wb, "Islet Size Slopes", data.frame(
    Note = paste("Islet size slopes: change in predicted percent", tolower(cell_label),
                 "per 1 unit increase in log(Islet Cells). 1 unit on log scale ≈ 2.72-fold increase in cell count.")))
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Overall Islet slope (averaged over Region and Sex):"), startRow = 3)
  writeData(wb, "Islet Size Slopes", overall_islet %>%
              select(predictor, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 4)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 6)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Islet slopes by Region:"), startRow = 7)
  writeData(wb, "Islet Size Slopes", islet_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 8)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 12)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Region contrasts (Islet slopes):"), startRow = 13)
  writeData(wb, "Islet Size Slopes", islet_region_contrasts_df %>% mutate(across(where(is.numeric), fmt)), startRow = 14)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 18)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Islet slopes by Sex:"), startRow = 19)
  writeData(wb, "Islet Size Slopes", islet_by_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 20)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 23)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Sex contrast (Islet slopes):"), startRow = 24)
  writeData(wb, "Islet Size Slopes", islet_sex_contrast_df %>% mutate(across(where(is.numeric), fmt)), startRow = 25)

  # Sheet 8: DD x Region x Age/AO (3-way)
  three_way_age_name <- paste0("DD x Region x ", age_label)
  # Truncate sheet name to 31 chars (Excel limit)
  three_way_age_name <- substr(three_way_age_name, 1, 31)
  addWorksheet(wb, three_way_age_name)
  writeData(wb, three_way_age_name, data.frame(
    Note = c(
      paste0("3-WAY INTERACTION: Disease Duration x Region x ", age_label),
      paste0("How the DD-", tolower(cell_label), " relationship changes with ", age_label, " across Regions."),
      paste0("DD slopes evaluated at ", age_label, " = mean ± 1 (centered), marginalized over Sex.")
    )))
  writeData(wb, three_way_age_name, data.frame(x = ""), startRow = 5)
  writeData(wb, three_way_age_name, data.frame(Section = paste0("DD slopes at low vs high ", age_label, ":")), startRow = 6)
  writeData(wb, three_way_age_name, dd_slope_at_age_df %>% mutate(across(where(is.numeric), fmt)), startRow = 7)
  writeData(wb, three_way_age_name, data.frame(x = ""), startRow = 14)
  writeData(wb, three_way_age_name, data.frame(
    Section = paste0("Interaction test: difference in DD slope (high ", age_label, " - low ", age_label, ") within each Region:")), startRow = 15)
  writeData(wb, three_way_age_name, dd_age_interaction_df %>% mutate(across(where(is.numeric), fmt)), startRow = 16)

  # Sheet 9: DD x Region x Islet (3-way)
  addWorksheet(wb, "DD x Region x Islet")
  writeData(wb, "DD x Region x Islet", data.frame(
    Note = c(
      "3-WAY INTERACTION: Disease Duration x Region x Islet Size",
      paste0("How the DD-", tolower(cell_label), " relationship changes with islet size across Regions."),
      "DD slopes evaluated at log_Islet.Cells_c = mean ± 1, marginalized over Sex."
    )))
  writeData(wb, "DD x Region x Islet", data.frame(x = ""), startRow = 5)
  writeData(wb, "DD x Region x Islet", data.frame(Section = "DD slopes at low vs high Islet size:"), startRow = 6)
  writeData(wb, "DD x Region x Islet", dd_slope_at_islet_df %>% mutate(across(where(is.numeric), fmt)), startRow = 7)
  writeData(wb, "DD x Region x Islet", data.frame(x = ""), startRow = 14)
  writeData(wb, "DD x Region x Islet", data.frame(
    Section = "Interaction test: difference in DD slope (high Islet - low Islet) within each Region:"), startRow = 15)
  writeData(wb, "DD x Region x Islet", dd_islet_interaction_df %>% mutate(across(where(is.numeric), fmt)), startRow = 16)

  # Sheet 10: Full Pairwise
  addWorksheet(wb, "Pairwise Region x Sex")
  writeData(wb, "Pairwise Region x Sex", data.frame(
    Note = c("All pairwise comparisons between the 6 Region x Sex cells.",
             "15 total comparisons. Posterior probabilities do not require multiple comparison corrections.",
             "")))
  writeData(wb, "Pairwise Region x Sex",
            full_pairwise_df %>% mutate(across(where(is.numeric), fmt)), startRow = 5)

  # Save
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0(
    cell_label, "_Content_", type_tag, "_Analysis.xlsx"))
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Saved to:", output_file, "\n")

  return(invisible(wb))
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

generate_all_proportion_t1d_excel <- function(output_dir = "Proportion/Results") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL PROPORTION T1D EXCEL RESULTS (AO parameterization)\n")
  cat("=============================================================================\n\n")

  # Load data
  cat("Loading data...\n")
  quad_t1d <- readRDS("Data/quad_data_t1d.rds")
  quad_t1d <- quad_t1d %>%
    mutate(
      Region  = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex     = factor(Sex, levels = c("Female", "Male")),
      Donor   = factor(Donor),
      ImageID = factor(ImageID)
    ) %>%
    droplevels()

  cat("quad_t1d:", nrow(quad_t1d), "rows,", n_distinct(quad_t1d$Donor), "donors\n\n")

  # Define cell types
  cell_types <- list(
    list(type = "ins",  label = "Insulin"),
    list(type = "glu",  label = "Glucagon"),
    list(type = "soma", label = "Somatostatin"),
    list(type = "pp",   label = "PP")
  )

  # Generate AO Excel files (DD slopes are extracted from the AO model directly)
  cat("--- AGE AT ONSET MODELS ---\n")
  for (ct in cell_types) {
    model_path <- paste0("Proportion/Models/", ct$type, "_ao_model.rds")
    cat("Loading", model_path, "...\n")
    model <- readRDS(model_path)
    generate_proportion_t1d_excel(
      model = model, data = quad_t1d,
      cell_type = ct$type, cell_label = ct$label,
      model_type = "ao", output_dir = output_dir
    )
    rm(model); gc()
  }

  # Summary
  cat("\n=============================================================================\n")
  cat("COMPLETE! Four Excel files created in", output_dir, ":\n\n")
  cat("  1. Insulin_Content_AO_Analysis.xlsx\n")
  cat("  2. Glucagon_Content_AO_Analysis.xlsx\n")
  cat("  3. Somatostatin_Content_AO_Analysis.xlsx\n")
  cat("  4. PP_Content_AO_Analysis.xlsx\n")
  cat("=============================================================================\n")
}

# Run
if (interactive()) {
  cat("\nTo generate all T1D proportion Excel files, run:\n")
  cat("  generate_all_proportion_t1d_excel()\n\n")
  cat("Or individually:\n")
  cat("  model <- readRDS('Proportion/Models/ins_ao_model.rds')\n")
  cat("  data  <- readRDS('Data/quad_data_t1d.rds')\n")
  cat("  generate_proportion_t1d_excel(model, data, 'ins', 'Insulin', model_type = 'ao')\n")
}
