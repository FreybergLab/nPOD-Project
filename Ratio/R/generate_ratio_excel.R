# =============================================================================
# Generate Comprehensive Excel Results for Ratio Models (ordbetareg)
# =============================================================================
#
# Models:
#   Diagnosis: ratio_powered (primary), ratio_mid (supplementary)
#   T1D only:  ratio_dd, ratio_ao
#
# All models: ordbetareg (Bayesian ordered beta regression)
# Response: beta_alpha_ratio = beta / (alpha + beta)
# Data: cutoff_data, Stain != "Single", Donor != "6473"
#
# Output: 4 Excel files in Results/Ratio/
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

# Generic slope summarizer for T1D models
summarize_slopes <- function(slopes_list, group_var, group_levels) {
  results_list <- list()
  draws_list <- list()
  for (g in group_levels) {
    subset_keys <- names(slopes_list)
    if (group_var == "Region") {
      subset_keys <- subset_keys[grepl(paste0("^", g, "_"), subset_keys)]
    } else if (group_var == "Sex") {
      subset_keys <- subset_keys[grepl(paste0("_", g, "$"), subset_keys)]
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


# =============================================================================
# 1. DIAGNOSIS MODEL — PRIMARY ANALYSIS (ratio_powered)
# =============================================================================
#
# beta_alpha_ratio ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
#   Diagnosis:Region:log_Islet.Cells_c + Diagnosis:Region:Age_c +
#   (1|Donor) + (1|Donor:ImageID)
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
#   9.  Islet Size Slopes
#  10.  Diag x Region x Age (3-way)
#  11.  Diag x Region x Islet (3-way)
#  12.  Pairwise Diag x Sex
#  13.  Pairwise Diag x Region
# =============================================================================

generate_ratio_diag_excel <- function(model, data, output_dir = "Ratio/Results") {

  cat("\n=============================================================================\n")
  cat("GENERATING PRIMARY ANALYSIS EXCEL: BETA:ALPHA RATIO — DIAGNOSIS\n")
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
  # CATEGORICAL GRID (Diagnosis × Region × Sex)
  # -------------------------------------------------------------------------
  cat("Setting up prediction grid...\n")

  cat_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Age_c = 0,
    log_Islet.Cells_c = 0
  )

  post_fitted <- get_posterior_fitted(model, cat_grid)
  col_names <- paste(cat_grid$Diagnosis, cat_grid$Region, cat_grid$Sex, sep = "_")
  colnames(post_fitted) <- col_names

  # -------------------------------------------------------------------------
  # SHEET 2: DIAGNOSIS
  # -------------------------------------------------------------------------
  cat("Creating Sheet 2: Diagnosis...\n")

  nd_idx  <- which(cat_grid$Diagnosis == "ND")
  t1d_idx <- which(cat_grid$Diagnosis == "T1D")
  nd_draws  <- rowMeans(post_fitted[, nd_idx])
  t1d_draws <- rowMeans(post_fitted[, t1d_idx])

  diagnosis_marginal <- bind_rows(
    summarize_posterior(nd_draws) %>% mutate(Diagnosis = "ND", .before = 1),
    summarize_posterior(t1d_draws) %>% mutate(Diagnosis = "T1D", .before = 1)
  ) %>% rename(response = estimate)

  diagnosis_contrast <- compute_contrast(nd_draws, t1d_draws, "ND - T1D")

  # -------------------------------------------------------------------------
  # SHEET 3: REGION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: Region...\n")

  head_draws <- rowMeans(post_fitted[, which(cat_grid$Region == "Head")])
  body_draws <- rowMeans(post_fitted[, which(cat_grid$Region == "Body")])
  tail_draws <- rowMeans(post_fitted[, which(cat_grid$Region == "Tail")])

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
  # SHEET 4: SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 4: Sex...\n")

  female_draws <- rowMeans(post_fitted[, which(cat_grid$Sex == "Female")])
  male_draws   <- rowMeans(post_fitted[, which(cat_grid$Sex == "Male")])

  sex_marginal <- bind_rows(
    summarize_posterior(female_draws) %>% mutate(Sex = "Female", .before = 1),
    summarize_posterior(male_draws) %>% mutate(Sex = "Male", .before = 1)
  ) %>% rename(response = estimate)

  sex_contrast <- compute_contrast(female_draws, male_draws, "Female - Male")

  # -------------------------------------------------------------------------
  # SHEET 5: DIAGNOSIS x REGION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 5: Diagnosis x Region...\n")

  diag_region_marginal <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == reg)
      draws <- rowMeans(post_fitted[, idx, drop = FALSE])
      summ <- summarize_posterior(draws)
      summ$Diagnosis <- diag; summ$Region <- reg
      diag_region_marginal[[paste(diag, reg, sep = "_")]] <- summ
    }
  }
  diag_region_marginal_df <- bind_rows(diag_region_marginal) %>%
    select(Diagnosis, Region, estimate, lower.HPD, upper.HPD) %>%
    rename(response = estimate)

  # Diagnosis within each Region
  diag_within_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    nd_d <- rowMeans(post_fitted[, which(cat_grid$Diagnosis == "ND" & cat_grid$Region == reg), drop = FALSE])
    t1d_d <- rowMeans(post_fitted[, which(cat_grid$Diagnosis == "T1D" & cat_grid$Region == reg), drop = FALSE])
    contrast <- compute_contrast(nd_d, t1d_d, "ND - T1D")
    contrast$Region <- reg
    diag_within_region[[reg]] <- contrast
  }
  diag_within_region_df <- bind_rows(diag_within_region) %>%
    select(contrast, Region, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # Region within each Diagnosis
  region_within_diag <- list()
  for (d in c("ND", "T1D")) {
    h <- rowMeans(post_fitted[, which(cat_grid$Region == "Head" & cat_grid$Diagnosis == d), drop = FALSE])
    b <- rowMeans(post_fitted[, which(cat_grid$Region == "Body" & cat_grid$Diagnosis == d), drop = FALSE])
    t <- rowMeans(post_fitted[, which(cat_grid$Region == "Tail" & cat_grid$Diagnosis == d), drop = FALSE])
    hb <- compute_contrast(h, b, "Head - Body"); hb$Diagnosis <- d
    ht <- compute_contrast(h, t, "Head - Tail"); ht$Diagnosis <- d
    bt <- compute_contrast(b, t, "Body - Tail"); bt$Diagnosis <- d
    region_within_diag[[paste0(d, "_HB")]] <- hb
    region_within_diag[[paste0(d, "_HT")]] <- ht
    region_within_diag[[paste0(d, "_BT")]] <- bt
  }
  region_within_diag_df <- bind_rows(region_within_diag) %>%
    select(contrast, Diagnosis, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 6: DIAGNOSIS x SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 6: Diagnosis x Sex...\n")

  diag_sex_marginal <- list()
  for (diag in c("ND", "T1D")) {
    for (s in c("Female", "Male")) {
      idx <- which(cat_grid$Diagnosis == diag & cat_grid$Sex == s)
      draws <- rowMeans(post_fitted[, idx, drop = FALSE])
      summ <- summarize_posterior(draws)
      summ$Diagnosis <- diag; summ$Sex <- s
      diag_sex_marginal[[paste(diag, s, sep = "_")]] <- summ
    }
  }
  diag_sex_marginal_df <- bind_rows(diag_sex_marginal) %>%
    select(Diagnosis, Sex, estimate, lower.HPD, upper.HPD) %>% rename(response = estimate)

  diag_within_sex <- list()
  for (s in c("Female", "Male")) {
    nd_d <- rowMeans(post_fitted[, which(cat_grid$Diagnosis == "ND" & cat_grid$Sex == s), drop = FALSE])
    t1d_d <- rowMeans(post_fitted[, which(cat_grid$Diagnosis == "T1D" & cat_grid$Sex == s), drop = FALSE])
    contrast <- compute_contrast(nd_d, t1d_d, "ND - T1D"); contrast$Sex <- s
    diag_within_sex[[s]] <- contrast
  }
  diag_within_sex_df <- bind_rows(diag_within_sex) %>%
    select(contrast, Sex, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  sex_within_diag <- list()
  for (d in c("ND", "T1D")) {
    f_d <- rowMeans(post_fitted[, which(cat_grid$Sex == "Female" & cat_grid$Diagnosis == d), drop = FALSE])
    m_d <- rowMeans(post_fitted[, which(cat_grid$Sex == "Male" & cat_grid$Diagnosis == d), drop = FALSE])
    contrast <- compute_contrast(f_d, m_d, "Female - Male"); contrast$Diagnosis <- d
    sex_within_diag[[d]] <- contrast
  }
  sex_within_diag_df <- bind_rows(sex_within_diag) %>%
    select(contrast, Diagnosis, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 7: REGION x SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 7: Region x Sex...\n")

  region_sex_marginal <- list()
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      idx <- which(cat_grid$Region == reg & cat_grid$Sex == s)
      draws <- rowMeans(post_fitted[, idx, drop = FALSE])
      summ <- summarize_posterior(draws)
      summ$Region <- reg; summ$Sex <- s
      region_sex_marginal[[paste(reg, s, sep = "_")]] <- summ
    }
  }
  region_sex_marginal_df <- bind_rows(region_sex_marginal) %>%
    select(Region, Sex, estimate, lower.HPD, upper.HPD) %>% rename(response = estimate)

  sex_within_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    f_d <- rowMeans(post_fitted[, which(cat_grid$Region == reg & cat_grid$Sex == "Female"), drop = FALSE])
    m_d <- rowMeans(post_fitted[, which(cat_grid$Region == reg & cat_grid$Sex == "Male"), drop = FALSE])
    contrast <- compute_contrast(f_d, m_d, "Female - Male"); contrast$Region <- reg
    sex_within_region[[reg]] <- contrast
  }
  sex_within_region_df <- bind_rows(sex_within_region) %>%
    select(contrast, Region, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  region_within_sex <- list()
  for (s in c("Female", "Male")) {
    h <- rowMeans(post_fitted[, which(cat_grid$Region == "Head" & cat_grid$Sex == s), drop = FALSE])
    b <- rowMeans(post_fitted[, which(cat_grid$Region == "Body" & cat_grid$Sex == s), drop = FALSE])
    t <- rowMeans(post_fitted[, which(cat_grid$Region == "Tail" & cat_grid$Sex == s), drop = FALSE])
    hb <- compute_contrast(h, b, "Head - Body"); hb$Sex <- s
    ht <- compute_contrast(h, t, "Head - Tail"); ht$Sex <- s
    bt <- compute_contrast(b, t, "Body - Tail"); bt$Sex <- s
    region_within_sex[[paste0(s, "_HB")]] <- hb
    region_within_sex[[paste0(s, "_HT")]] <- ht
    region_within_sex[[paste0(s, "_BT")]] <- bt
  }
  region_within_sex_df <- bind_rows(region_within_sex) %>%
    select(contrast, Sex, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # CONTINUOUS SLOPES SETUP
  # -------------------------------------------------------------------------
  cat("Setting up continuous predictor grids...\n")

  age_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), Age_c = c(-1, 1), log_Islet.Cells_c = 0
  )
  islet_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), Age_c = 0, log_Islet.Cells_c = c(-1, 1)
  )

  post_age   <- get_posterior_fitted(model, age_grid)
  post_islet <- get_posterior_fitted(model, islet_grid)

  age_slopes_list   <- list()
  islet_slopes_list <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        key <- paste(diag, reg, s, sep = "_")
        lo <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & age_grid$Sex == s & age_grid$Age_c == -1)
        hi <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & age_grid$Sex == s & age_grid$Age_c == 1)
        age_slopes_list[[key]] <- (post_age[, hi] - post_age[, lo]) / 2
        lo <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == -1)
        hi <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == 1)
        islet_slopes_list[[key]] <- (post_islet[, hi] - post_islet[, lo]) / 2
      }
    }
  }

  # -------------------------------------------------------------------------
  # Generic slope decomposition (reused for Age and Islet)
  # -------------------------------------------------------------------------
  decompose_slopes <- function(slopes_list) {
    # Overall
    all_slopes <- rowMeans(do.call(cbind, slopes_list))
    overall <- summarize_posterior(all_slopes)
    overall$prob_positive <- mean(all_slopes > 0)
    overall$prob_negative <- mean(all_slopes < 0)

    # By Diagnosis
    by_diag <- list(); by_diag_draws <- list()
    for (diag in c("ND", "T1D")) {
      keys <- grep(paste0("^", diag, "_"), names(slopes_list), value = TRUE)
      draws <- rowMeans(do.call(cbind, slopes_list[keys]))
      summ <- summarize_posterior(draws); summ$Diagnosis <- diag
      summ$prob_positive <- mean(draws > 0); summ$prob_negative <- mean(draws < 0)
      by_diag[[diag]] <- summ; by_diag_draws[[diag]] <- draws
    }
    by_diag_df <- bind_rows(by_diag) %>% select(Diagnosis, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
    diag_diff <- compute_contrast(by_diag_draws[["T1D"]], by_diag_draws[["ND"]], "T1D - ND")

    # By Region
    by_region_draws <- list()
    by_region <- list()
    for (reg in c("Head", "Body", "Tail")) {
      keys <- grep(paste0("_", reg, "_"), names(slopes_list), value = TRUE)
      draws <- rowMeans(do.call(cbind, slopes_list[keys]))
      summ <- summarize_posterior(draws); summ$Region <- reg
      summ$prob_positive <- mean(draws > 0); summ$prob_negative <- mean(draws < 0)
      by_region[[reg]] <- summ; by_region_draws[[reg]] <- draws
    }
    by_region_df <- bind_rows(by_region) %>% select(Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
    region_contr <- bind_rows(
      compute_contrast(by_region_draws[["Head"]], by_region_draws[["Body"]], "Head - Body"),
      compute_contrast(by_region_draws[["Head"]], by_region_draws[["Tail"]], "Head - Tail"),
      compute_contrast(by_region_draws[["Body"]], by_region_draws[["Tail"]], "Body - Tail")
    ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

    # By Sex
    by_sex_draws <- list(); by_sex <- list()
    for (s in c("Female", "Male")) {
      keys <- grep(paste0("_", s, "$"), names(slopes_list), value = TRUE)
      draws <- rowMeans(do.call(cbind, slopes_list[keys]))
      summ <- summarize_posterior(draws); summ$Sex <- s
      summ$prob_positive <- mean(draws > 0); summ$prob_negative <- mean(draws < 0)
      by_sex[[s]] <- summ; by_sex_draws[[s]] <- draws
    }
    by_sex_df <- bind_rows(by_sex) %>% select(Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
    sex_contr <- compute_contrast(by_sex_draws[["Female"]], by_sex_draws[["Male"]], "Female - Male") %>%
      select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

    list(overall = overall, by_diag = by_diag_df, diag_diff = diag_diff,
         by_region = by_region_df, region_contr = region_contr,
         by_sex = by_sex_df, sex_contr = sex_contr,
         by_diag_draws = by_diag_draws, by_region_draws = by_region_draws)
  }

  # -------------------------------------------------------------------------
  # SHEET 8: AGE SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 8: Age Slopes...\n")
  age_decomp <- decompose_slopes(age_slopes_list)

  # -------------------------------------------------------------------------
  # SHEET 9: ISLET SIZE SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 9: Islet Size Slopes...\n")
  islet_decomp <- decompose_slopes(islet_slopes_list)

  # -------------------------------------------------------------------------
  # SHEET 10: DIAG x REGION x AGE (3-way)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 10: Diag x Region x Age (3-way)...\n")

  age_diag_region <- list(); age_diag_region_draws <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      keys <- paste(diag, reg, c("Female", "Male"), sep = "_")
      draws <- rowMeans(do.call(cbind, age_slopes_list[keys]))
      summ <- summarize_posterior(draws)
      summ$Diagnosis <- diag; summ$Region <- reg
      summ$prob_positive <- mean(draws > 0); summ$prob_negative <- mean(draws < 0)
      key <- paste(diag, reg, sep = "_")
      age_diag_region[[key]] <- summ; age_diag_region_draws[[key]] <- draws
    }
  }
  age_diag_region_df <- bind_rows(age_diag_region) %>%
    select(Diagnosis, Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

  age_slope_diff_by_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff <- age_diag_region_draws[[paste("T1D", reg, sep = "_")]] -
      age_diag_region_draws[[paste("ND", reg, sep = "_")]]
    summ <- summarize_posterior(diff); summ$Region <- reg; summ$contrast <- "T1D - ND"
    summ$prob_positive <- mean(diff > 0); summ$prob_negative <- mean(diff < 0)
    age_slope_diff_by_region[[reg]] <- summ
  }
  age_slope_diff_by_region_df <- bind_rows(age_slope_diff_by_region) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

  # -------------------------------------------------------------------------
  # SHEET 11: DIAG x REGION x ISLET (3-way)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 11: Diag x Region x Islet (3-way)...\n")

  islet_diag_region <- list(); islet_diag_region_draws <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      keys <- paste(diag, reg, c("Female", "Male"), sep = "_")
      draws <- rowMeans(do.call(cbind, islet_slopes_list[keys]))
      summ <- summarize_posterior(draws)
      summ$Diagnosis <- diag; summ$Region <- reg
      summ$prob_positive <- mean(draws > 0); summ$prob_negative <- mean(draws < 0)
      key <- paste(diag, reg, sep = "_")
      islet_diag_region[[key]] <- summ; islet_diag_region_draws[[key]] <- draws
    }
  }
  islet_diag_region_df <- bind_rows(islet_diag_region) %>%
    select(Diagnosis, Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

  islet_slope_diff_by_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff <- islet_diag_region_draws[[paste("T1D", reg, sep = "_")]] -
      islet_diag_region_draws[[paste("ND", reg, sep = "_")]]
    summ <- summarize_posterior(diff); summ$Region <- reg; summ$contrast <- "T1D - ND"
    summ$prob_positive <- mean(diff > 0); summ$prob_negative <- mean(diff < 0)
    islet_slope_diff_by_region[[reg]] <- summ
  }
  islet_slope_diff_by_region_df <- bind_rows(islet_slope_diff_by_region) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

  # Region contrasts on islet slopes within each Diagnosis
  islet_region_within_diag <- list()
  for (d in c("ND", "T1D")) {
    hd <- islet_diag_region_draws[[paste(d, "Head", sep = "_")]]
    bd <- islet_diag_region_draws[[paste(d, "Body", sep = "_")]]
    td <- islet_diag_region_draws[[paste(d, "Tail", sep = "_")]]
    hb <- compute_contrast(hd, bd, "Head - Body"); hb$Diagnosis <- d
    ht <- compute_contrast(hd, td, "Head - Tail"); ht$Diagnosis <- d
    bt <- compute_contrast(bd, td, "Body - Tail"); bt$Diagnosis <- d
    islet_region_within_diag[[paste0(d, "_HB")]] <- hb
    islet_region_within_diag[[paste0(d, "_HT")]] <- ht
    islet_region_within_diag[[paste0(d, "_BT")]] <- bt
  }
  islet_region_within_diag_df <- bind_rows(islet_region_within_diag) %>%
    select(Diagnosis, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEETS 12-13: SUPPLEMENTARY PAIRWISE
  # -------------------------------------------------------------------------
  cat("Creating Sheets 12-13: Supplementary Pairwise...\n")

  # Diagnosis x Sex pairwise
  pairwise_diag_sex <- list()
  combos <- list(c("ND","Female","T1D","Female"), c("ND","Female","ND","Male"),
                 c("ND","Female","T1D","Male"), c("T1D","Female","ND","Male"),
                 c("T1D","Female","T1D","Male"), c("ND","Male","T1D","Male"))
  for (combo in combos) {
    idx1 <- which(cat_grid$Diagnosis == combo[1] & cat_grid$Sex == combo[2])
    idx2 <- which(cat_grid$Diagnosis == combo[3] & cat_grid$Sex == combo[4])
    d1 <- rowMeans(post_fitted[, idx1, drop = FALSE])
    d2 <- rowMeans(post_fitted[, idx2, drop = FALSE])
    label <- paste0(combo[1], " ", combo[2], " - ", combo[3], " ", combo[4])
    pairwise_diag_sex[[label]] <- compute_contrast(d1, d2, label)
  }
  pairwise_diag_sex_df <- bind_rows(pairwise_diag_sex) %>%
    select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # Diagnosis x Region pairwise
  diag_region_combos <- expand.grid(d1 = c("ND","T1D"), r1 = c("Head","Body","Tail"),
                                    d2 = c("ND","T1D"), r2 = c("Head","Body","Tail"),
                                    stringsAsFactors = FALSE) %>%
    mutate(c1 = paste(d1, r1), c2 = paste(d2, r2)) %>% filter(c1 < c2) %>% select(d1, r1, d2, r2)

  pairwise_diag_region <- list()
  for (i in 1:nrow(diag_region_combos)) {
    row <- diag_region_combos[i, ]
    idx1 <- which(cat_grid$Diagnosis == row$d1 & cat_grid$Region == row$r1)
    idx2 <- which(cat_grid$Diagnosis == row$d2 & cat_grid$Region == row$r2)
    d1 <- rowMeans(post_fitted[, idx1, drop = FALSE])
    d2 <- rowMeans(post_fitted[, idx2, drop = FALSE])
    label <- paste0(row$d1, " ", row$r1, " - ", row$d2, " ", row$r2)
    pairwise_diag_region[[label]] <- compute_contrast(d1, d2, label)
  }
  pairwise_diag_region_df <- bind_rows(pairwise_diag_region) %>%
    select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # CREATE WORKBOOK
  # -------------------------------------------------------------------------
  cat("\nCreating Excel workbook...\n")
  wb <- createWorkbook()

  # Helper for slope sheets
  write_slope_sheet <- function(wb, sheet_name, predictor_label, decomp) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data.frame(
      Note = paste(predictor_label, "slopes: change in predicted beta:alpha ratio per 1 unit increase.")))
    writeData(wb, sheet_name, data.frame(Section = paste0("Overall ", predictor_label, " slope:")), startRow = 3)
    writeData(wb, sheet_name, decomp$overall %>%
                select(estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
                mutate(across(where(is.numeric), fmt)), startRow = 4)
    writeData(wb, sheet_name, data.frame(Section = paste(predictor_label, "slopes by Diagnosis:")), startRow = 7)
    writeData(wb, sheet_name, decomp$by_diag %>% mutate(across(where(is.numeric), fmt)), startRow = 8)
    writeData(wb, sheet_name, data.frame(Section = "Diagnosis contrast (T1D - ND):"), startRow = 11)
    writeData(wb, sheet_name, decomp$diag_diff %>%
                select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
                mutate(across(where(is.numeric), fmt)), startRow = 12)
    writeData(wb, sheet_name, data.frame(Section = paste(predictor_label, "slopes by Region:")), startRow = 15)
    writeData(wb, sheet_name, decomp$by_region %>% mutate(across(where(is.numeric), fmt)), startRow = 16)
    writeData(wb, sheet_name, data.frame(Section = "Region contrasts:"), startRow = 20)
    writeData(wb, sheet_name, decomp$region_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 21)
    writeData(wb, sheet_name, data.frame(Section = paste(predictor_label, "slopes by Sex:")), startRow = 25)
    writeData(wb, sheet_name, decomp$by_sex %>% mutate(across(where(is.numeric), fmt)), startRow = 26)
    writeData(wb, sheet_name, data.frame(Section = "Sex contrast:"), startRow = 29)
    writeData(wb, sheet_name, decomp$sex_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 30)
  }

  # Sheet 1
  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", data.frame(
    Note = c("PRIMARY ANALYSIS: Beta:Alpha Ratio — Diagnosis",
             "Bayesian ordered beta regression (ordbetareg)",
             "Response: beta_alpha_ratio = beta / (alpha + beta), 0-1 scale",
             "Model: (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +",
             "       Diagnosis:Region:Age_c + Diagnosis:Region:log_Islet.Cells_c + (1|Donor) + (1|Donor:ImageID)",
             "All estimates are posterior medians with 95% HPD intervals.",
             "prob_greater_0 / prob_less_0 = posterior probability contrast is positive/negative.",
             "")))
  writeData(wb, "Model Summary", model_summary, startRow = 10)

  # Sheet 2: Diagnosis
  addWorksheet(wb, "Diagnosis")
  writeData(wb, "Diagnosis", data.frame(Note = "Main effect of Diagnosis on beta:alpha ratio, averaged over Region and Sex."))
  writeData(wb, "Diagnosis", diagnosis_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diagnosis", data.frame(Section = "Contrast:"), startRow = 6)
  writeData(wb, "Diagnosis", diagnosis_contrast %>%
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 7)

  # Sheet 3: Region
  addWorksheet(wb, "Region")
  writeData(wb, "Region", data.frame(Note = "Main effect of Region, averaged over Diagnosis and Sex."))
  writeData(wb, "Region", region_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region", data.frame(Section = "Pairwise contrasts:"), startRow = 7)
  writeData(wb, "Region", region_contrasts %>%
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 8)

  # Sheet 4: Sex
  addWorksheet(wb, "Sex")
  writeData(wb, "Sex", data.frame(Note = "Main effect of Sex, averaged over Diagnosis and Region."))
  writeData(wb, "Sex", sex_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Sex", data.frame(Section = "Contrast:"), startRow = 6)
  writeData(wb, "Sex", sex_contrast %>%
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 7)

  # Sheet 5: Diagnosis x Region
  addWorksheet(wb, "Diagnosis x Region")
  writeData(wb, "Diagnosis x Region", data.frame(Note = "Diagnosis x Region interaction, averaged over Sex."))
  writeData(wb, "Diagnosis x Region", diag_region_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diagnosis x Region", data.frame(Section = "Diagnosis effect (ND - T1D) within each Region:"), startRow = 10)
  writeData(wb, "Diagnosis x Region", diag_within_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 11)
  writeData(wb, "Diagnosis x Region", data.frame(Section = "Region contrasts within each Diagnosis:"), startRow = 15)
  writeData(wb, "Diagnosis x Region", region_within_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 16)

  # Sheet 6: Diagnosis x Sex
  addWorksheet(wb, "Diagnosis x Sex")
  writeData(wb, "Diagnosis x Sex", data.frame(Note = "Diagnosis x Sex interaction, averaged over Region."))
  writeData(wb, "Diagnosis x Sex", diag_sex_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diagnosis x Sex", data.frame(Section = "Diagnosis effect (ND - T1D) within each Sex:"), startRow = 8)
  writeData(wb, "Diagnosis x Sex", diag_within_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 9)
  writeData(wb, "Diagnosis x Sex", data.frame(Section = "Sex effect (Female - Male) within each Diagnosis:"), startRow = 12)
  writeData(wb, "Diagnosis x Sex", sex_within_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 13)

  # Sheet 7: Region x Sex
  addWorksheet(wb, "Region x Sex")
  writeData(wb, "Region x Sex", data.frame(Note = "Region x Sex interaction, averaged over Diagnosis."))
  writeData(wb, "Region x Sex", region_sex_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region x Sex", data.frame(Section = "Sex effect (Female - Male) within each Region:"), startRow = 10)
  writeData(wb, "Region x Sex", sex_within_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 11)
  writeData(wb, "Region x Sex", data.frame(Section = "Region contrasts within each Sex:"), startRow = 15)
  writeData(wb, "Region x Sex", region_within_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 16)

  # Sheets 8-9: Slopes
  write_slope_sheet(wb, "Age Slopes", "Age", age_decomp)
  write_slope_sheet(wb, "Islet Size Slopes", "Islet size", islet_decomp)

  # Sheet 10: Diag x Region x Age
  addWorksheet(wb, "Diag x Region x Age")
  writeData(wb, "Diag x Region x Age", data.frame(
    Note = "3-WAY: Diagnosis x Region x Age. Age slopes per Diagnosis x Region cell (averaged over Sex)."))
  writeData(wb, "Diag x Region x Age", age_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diag x Region x Age", data.frame(Section = "Diagnosis effect (T1D - ND) on Age slopes within each Region:"), startRow = 10)
  writeData(wb, "Diag x Region x Age", age_slope_diff_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 11)

  # Sheet 11: Diag x Region x Islet
  addWorksheet(wb, "Diag x Region x Islet")
  writeData(wb, "Diag x Region x Islet", data.frame(
    Note = "3-WAY: Diagnosis x Region x Islet Size. Islet slopes per Diagnosis x Region cell (averaged over Sex)."))
  writeData(wb, "Diag x Region x Islet", islet_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diag x Region x Islet", data.frame(Section = "Diagnosis effect (T1D - ND) on Islet slopes within each Region:"), startRow = 10)
  writeData(wb, "Diag x Region x Islet", islet_slope_diff_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 11)
  writeData(wb, "Diag x Region x Islet", data.frame(Section = "Region contrasts on Islet slopes within each Diagnosis:"), startRow = 15)
  writeData(wb, "Diag x Region x Islet", islet_region_within_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 16)

  # Sheet 12: Pairwise Diagnosis x Sex
  addWorksheet(wb, "Pairwise Diag x Sex")
  writeData(wb, "Pairwise Diag x Sex", data.frame(
    Note = c("All pairwise comparisons between the 4 Diagnosis x Sex cells.", "")))
  writeData(wb, "Pairwise Diag x Sex", pairwise_diag_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 4)

  # Sheet 13: Pairwise Diagnosis x Region
  addWorksheet(wb, "Pairwise Diag x Region")
  writeData(wb, "Pairwise Diag x Region", data.frame(
    Note = c("All pairwise comparisons between the 6 Diagnosis x Region cells (15 total).", "")))
  writeData(wb, "Pairwise Diag x Region", pairwise_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 4)

  # Save
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, "Ratio_Primary_Analysis.xlsx")
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Saved to:", output_file, "\n")
  return(invisible(wb))
}


# =============================================================================
# 2. T1D MODELS — DD & AO (ratio_dd, ratio_ao)
# =============================================================================
#
# Identical structure to proportion T1D models but with beta_alpha_ratio outcome.
# Reuses the parameterized DD/AO approach.
#
# DD: beta_alpha_ratio ~ (Disease.Duration_c + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
#     Disease.Duration_c:Region:log_Islet.Cells_c + Disease.Duration_c:Region:Age_c +
#     (1|Donor) + (1|Donor:ImageID)
#
# AO: same but age_at_onset_c replaces Age_c
# =============================================================================

generate_ratio_t1d_excel <- function(model, data, model_type = "dd",
                                     output_dir = "Ratio/Results") {

  age_var   <- if (model_type == "dd") "Age_c" else "age_at_onset_c"
  age_label <- if (model_type == "dd") "Age" else "Age at Onset"
  type_label <- if (model_type == "dd") "Disease Duration" else "Age at Onset"
  type_tag   <- toupper(model_type)

  cat("\n=============================================================================\n")
  cat("GENERATING", type_tag, "ANALYSIS EXCEL: BETA:ALPHA RATIO (T1D)\n")
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
  # -------------------------------------------------------------------------
  cat("Setting up prediction grid...\n")
  cat_grid <- expand.grid(
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Disease.Duration_c = 0, log_Islet.Cells_c = 0
  )
  cat_grid[[age_var]] <- 0

  post_fitted <- get_posterior_fitted(model, cat_grid)
  col_names <- paste(cat_grid$Region, cat_grid$Sex, sep = "_")
  colnames(post_fitted) <- col_names

  # -------------------------------------------------------------------------
  # SHEET 2: REGION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 2: Region...\n")
  head_draws <- rowMeans(post_fitted[, which(cat_grid$Region == "Head")])
  body_draws <- rowMeans(post_fitted[, which(cat_grid$Region == "Body")])
  tail_draws <- rowMeans(post_fitted[, which(cat_grid$Region == "Tail")])

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
  # SHEET 3: SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: Sex...\n")
  female_draws <- rowMeans(post_fitted[, which(cat_grid$Sex == "Female")])
  male_draws   <- rowMeans(post_fitted[, which(cat_grid$Sex == "Male")])

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
      summ <- summarize_posterior(draws); summ$Region <- reg; summ$Sex <- s
      region_sex_marginal[[paste(reg, s, sep = "_")]] <- summ
    }
  }
  region_sex_marginal_df <- bind_rows(region_sex_marginal) %>%
    select(Region, Sex, estimate, lower.HPD, upper.HPD) %>% rename(response = estimate)

  sex_within_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    f <- rowMeans(post_fitted[, which(cat_grid$Region == reg & cat_grid$Sex == "Female"), drop = FALSE])
    m <- rowMeans(post_fitted[, which(cat_grid$Region == reg & cat_grid$Sex == "Male"), drop = FALSE])
    c <- compute_contrast(f, m, "Female - Male"); c$Region <- reg
    sex_within_region[[reg]] <- c
  }
  sex_within_region_df <- bind_rows(sex_within_region) %>%
    select(contrast, Region, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  region_within_sex <- list()
  for (s in c("Female", "Male")) {
    h <- rowMeans(post_fitted[, which(cat_grid$Region == "Head" & cat_grid$Sex == s), drop = FALSE])
    b <- rowMeans(post_fitted[, which(cat_grid$Region == "Body" & cat_grid$Sex == s), drop = FALSE])
    t <- rowMeans(post_fitted[, which(cat_grid$Region == "Tail" & cat_grid$Sex == s), drop = FALSE])
    hb <- compute_contrast(h, b, "Head - Body"); hb$Sex <- s
    ht <- compute_contrast(h, t, "Head - Tail"); ht$Sex <- s
    bt <- compute_contrast(b, t, "Body - Tail"); bt$Sex <- s
    region_within_sex[[paste0(s, "_HB")]] <- hb
    region_within_sex[[paste0(s, "_HT")]] <- ht
    region_within_sex[[paste0(s, "_BT")]] <- bt
  }
  region_within_sex_df <- bind_rows(region_within_sex) %>%
    select(contrast, Sex, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # CONTINUOUS SLOPES
  # -------------------------------------------------------------------------
  cat("Computing slopes...\n")

  # DD slopes
  dd_grid <- expand.grid(Region = c("Head","Body","Tail"), Sex = c("Female","Male"),
                         Disease.Duration_c = c(-1, 1), log_Islet.Cells_c = 0)
  dd_grid[[age_var]] <- 0

  # Age/AO slopes
  age_grid <- expand.grid(Region = c("Head","Body","Tail"), Sex = c("Female","Male"),
                          Disease.Duration_c = 0, log_Islet.Cells_c = 0)
  age_grid_low <- age_grid; age_grid_low[[age_var]] <- -1
  age_grid_high <- age_grid; age_grid_high[[age_var]] <- 1
  age_grid_full <- bind_rows(age_grid_low, age_grid_high)

  # Islet slopes
  islet_grid <- expand.grid(Region = c("Head","Body","Tail"), Sex = c("Female","Male"),
                            Disease.Duration_c = 0, log_Islet.Cells_c = c(-1, 1))
  islet_grid[[age_var]] <- 0

  post_dd    <- get_posterior_fitted(model, dd_grid)
  post_age   <- get_posterior_fitted(model, age_grid_full)
  post_islet <- get_posterior_fitted(model, islet_grid)

  dd_slopes_list <- list(); age_slopes_list <- list(); islet_slopes_list <- list()
  n_cells <- nrow(age_grid)
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      key <- paste(reg, s, sep = "_")
      lo <- which(dd_grid$Region == reg & dd_grid$Sex == s & dd_grid$Disease.Duration_c == -1)
      hi <- which(dd_grid$Region == reg & dd_grid$Sex == s & dd_grid$Disease.Duration_c == 1)
      dd_slopes_list[[key]] <- (post_dd[, hi] - post_dd[, lo]) / 2

      lo <- which(age_grid_low$Region == reg & age_grid_low$Sex == s)
      hi <- which(age_grid_high$Region == reg & age_grid_high$Sex == s) + n_cells
      age_slopes_list[[key]] <- (post_age[, hi] - post_age[, lo]) / 2

      lo <- which(islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == -1)
      hi <- which(islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == 1)
      islet_slopes_list[[key]] <- (post_islet[, hi] - post_islet[, lo]) / 2
    }
  }

  # Summarize each predictor
  summarize_one_slope <- function(slopes_list, predictor_name) {
    overall <- rowMeans(do.call(cbind, slopes_list))
    overall_summ <- summarize_posterior(overall)
    overall_summ$predictor <- predictor_name
    overall_summ$prob_positive <- mean(overall > 0)
    overall_summ$prob_negative <- mean(overall < 0)

    by_region <- summarize_slopes(slopes_list, "Region", c("Head", "Body", "Tail"))
    region_contr <- bind_rows(
      compute_contrast(by_region$draws[["Head"]], by_region$draws[["Body"]], "Head - Body"),
      compute_contrast(by_region$draws[["Head"]], by_region$draws[["Tail"]], "Head - Tail"),
      compute_contrast(by_region$draws[["Body"]], by_region$draws[["Tail"]], "Body - Tail")
    ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

    by_sex <- summarize_slopes(slopes_list, "Sex", c("Female", "Male"))
    sex_contr <- compute_contrast(by_sex$draws[["Female"]], by_sex$draws[["Male"]], "Female - Male") %>%
      select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

    list(overall = overall_summ,
         by_region = by_region$summary %>% select(Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative),
         region_contr = region_contr,
         by_sex = by_sex$summary %>% select(Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative),
         sex_contr = sex_contr)
  }

  dd_summ    <- summarize_one_slope(dd_slopes_list, "Disease.Duration_c")
  age_summ   <- summarize_one_slope(age_slopes_list, age_var)
  islet_summ <- summarize_one_slope(islet_slopes_list, "log_Islet.Cells_c")

  # -------------------------------------------------------------------------
  # 3-WAY: DD × Region × Age/AO
  # -------------------------------------------------------------------------
  cat("Computing 3-way interactions...\n")

  compute_3way <- function(slopes_list, model, dd_grid_template, modifier_var, modifier_label) {
    grid_lo <- dd_grid_template; grid_lo[[modifier_var]] <- -1
    grid_hi <- dd_grid_template; grid_hi[[modifier_var]] <- 1
    post_lo <- get_posterior_fitted(model, grid_lo)
    post_hi <- get_posterior_fitted(model, grid_hi)

    slope_at <- list(); slope_at_draws <- list()
    for (level in c("low", "high")) {
      pm <- if (level == "low") post_lo else post_hi
      for (reg in c("Head", "Body", "Tail")) {
        sex_slopes <- list()
        for (s in c("Female", "Male")) {
          lo <- which(dd_grid_template$Region == reg & dd_grid_template$Sex == s & dd_grid_template$Disease.Duration_c == -1)
          hi <- which(dd_grid_template$Region == reg & dd_grid_template$Sex == s & dd_grid_template$Disease.Duration_c == 1)
          sex_slopes[[s]] <- (pm[, hi] - pm[, lo]) / 2
        }
        draws <- rowMeans(do.call(cbind, sex_slopes))
        key <- paste(reg, level, sep = "_")
        summ <- summarize_posterior(draws)
        summ$Region <- reg
        summ$level <- paste0(modifier_label, " = mean ", ifelse(level == "low", "- 1", "+ 1"))
        summ$prob_positive <- mean(draws > 0); summ$prob_negative <- mean(draws < 0)
        slope_at[[key]] <- summ; slope_at_draws[[key]] <- draws
      }
    }
    slopes_df <- bind_rows(slope_at) %>%
      select(Region, level, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

    interaction <- list()
    for (reg in c("Head", "Body", "Tail")) {
      diff <- slope_at_draws[[paste(reg, "high", sep = "_")]] -
        slope_at_draws[[paste(reg, "low", sep = "_")]]
      summ <- summarize_posterior(diff); summ$Region <- reg
      summ$contrast <- paste0("DD slope at high ", modifier_label, " - low ", modifier_label)
      summ$prob_greater_0 <- mean(diff > 0); summ$prob_less_0 <- mean(diff < 0)
      interaction[[reg]] <- summ
    }
    interaction_df <- bind_rows(interaction) %>%
      select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

    list(slopes = slopes_df, interaction = interaction_df)
  }

  dd_grid_base <- expand.grid(Region = c("Head","Body","Tail"), Sex = c("Female","Male"),
                              Disease.Duration_c = c(-1, 1), log_Islet.Cells_c = 0)
  dd_grid_base[[age_var]] <- 0

  three_way_age <- compute_3way(dd_slopes_list, model, dd_grid_base, age_var, age_label)
  three_way_islet <- compute_3way(dd_slopes_list, model, dd_grid_base, "log_Islet.Cells_c", "Islet size")

  # -------------------------------------------------------------------------
  # CREATE WORKBOOK
  # -------------------------------------------------------------------------
  cat("\nCreating Excel workbook...\n")
  wb <- createWorkbook()

  # Helper for slope sheets
  write_t1d_slope_sheet <- function(wb, sheet_name, pred_label, summ) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data.frame(
      Note = paste(pred_label, "slopes: change in predicted ratio per 1 unit increase (T1D only).")))
    writeData(wb, sheet_name, data.frame(Section = paste0("Overall ", pred_label, " slope:")), startRow = 3)
    writeData(wb, sheet_name, summ$overall %>%
                select(predictor, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
                mutate(across(where(is.numeric), fmt)), startRow = 4)
    writeData(wb, sheet_name, data.frame(Section = paste(pred_label, "slopes by Region:")), startRow = 7)
    writeData(wb, sheet_name, summ$by_region %>% mutate(across(where(is.numeric), fmt)), startRow = 8)
    writeData(wb, sheet_name, data.frame(Section = "Region contrasts:"), startRow = 12)
    writeData(wb, sheet_name, summ$region_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 13)
    writeData(wb, sheet_name, data.frame(Section = paste(pred_label, "slopes by Sex:")), startRow = 17)
    writeData(wb, sheet_name, summ$by_sex %>% mutate(across(where(is.numeric), fmt)), startRow = 18)
    writeData(wb, sheet_name, data.frame(Section = "Sex contrast:"), startRow = 21)
    writeData(wb, sheet_name, summ$sex_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 22)
  }

  # Sheet 1: Model Summary
  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", data.frame(
    Note = c(paste(toupper(type_label), "ANALYSIS: Beta:Alpha Ratio (T1D only)"),
             "Bayesian ordered beta regression (ordbetareg)",
             "Response: beta_alpha_ratio = beta / (alpha + beta), 0-1 scale",
             paste0("Model: (Disease.Duration_c + Region + Sex + ", age_var, " + log_Islet.Cells_c)^2 +"),
             paste0("       Disease.Duration_c:Region:log_Islet.Cells_c + Disease.Duration_c:Region:", age_var),
             "       + (1|Donor) + (1|Donor:ImageID)",
             "All estimates are posterior medians with 95% HPD intervals.",
             "")))
  writeData(wb, "Model Summary", model_summary, startRow = 10)

  # Sheet 2: Region
  addWorksheet(wb, "Region")
  writeData(wb, "Region", data.frame(Note = "Main effect of Region, averaged over Sex."))
  writeData(wb, "Region", region_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region", data.frame(Section = "Pairwise contrasts:"), startRow = 7)
  writeData(wb, "Region", region_contrasts %>%
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 8)

  # Sheet 3: Sex
  addWorksheet(wb, "Sex")
  writeData(wb, "Sex", data.frame(Note = "Main effect of Sex, averaged over Region."))
  writeData(wb, "Sex", sex_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Sex", data.frame(Section = "Contrast:"), startRow = 6)
  writeData(wb, "Sex", sex_contrast %>%
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 7)

  # Sheet 4: Region x Sex
  addWorksheet(wb, "Region x Sex")
  writeData(wb, "Region x Sex", data.frame(Note = "Region x Sex interaction."))
  writeData(wb, "Region x Sex", region_sex_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region x Sex", data.frame(Section = "Sex effect (Female - Male) within each Region:"), startRow = 10)
  writeData(wb, "Region x Sex", sex_within_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 11)
  writeData(wb, "Region x Sex", data.frame(Section = "Region contrasts within each Sex:"), startRow = 15)
  writeData(wb, "Region x Sex", region_within_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 16)

  # Sheets 5-7: Slopes
  write_t1d_slope_sheet(wb, "DD Slopes", "Disease Duration", dd_summ)
  age_sheet_name <- paste0(age_label, " Slopes")
  write_t1d_slope_sheet(wb, age_sheet_name, age_label, age_summ)
  write_t1d_slope_sheet(wb, "Islet Size Slopes", "Islet size", islet_summ)

  # Sheet 8: DD x Region x Age/AO
  twa_name <- substr(paste0("DD x Region x ", age_label), 1, 31)
  addWorksheet(wb, twa_name)
  writeData(wb, twa_name, data.frame(
    Note = c(paste0("3-WAY: Disease Duration x Region x ", age_label),
             paste0("DD slopes evaluated at ", age_label, " = mean ± 1, marginalized over Sex."))))
  writeData(wb, twa_name, data.frame(Section = paste0("DD slopes at low vs high ", age_label, ":")), startRow = 4)
  writeData(wb, twa_name, three_way_age$slopes %>% mutate(across(where(is.numeric), fmt)), startRow = 5)
  writeData(wb, twa_name, data.frame(Section = "Interaction test (high - low) within each Region:"), startRow = 12)
  writeData(wb, twa_name, three_way_age$interaction %>% mutate(across(where(is.numeric), fmt)), startRow = 13)

  # Sheet 9: DD x Region x Islet
  addWorksheet(wb, "DD x Region x Islet")
  writeData(wb, "DD x Region x Islet", data.frame(
    Note = c("3-WAY: Disease Duration x Region x Islet Size",
             "DD slopes evaluated at log_Islet.Cells_c = mean ± 1, marginalized over Sex.")))
  writeData(wb, "DD x Region x Islet", data.frame(Section = "DD slopes at low vs high Islet size:"), startRow = 4)
  writeData(wb, "DD x Region x Islet", three_way_islet$slopes %>% mutate(across(where(is.numeric), fmt)), startRow = 5)
  writeData(wb, "DD x Region x Islet", data.frame(Section = "Interaction test (high - low) within each Region:"), startRow = 12)
  writeData(wb, "DD x Region x Islet", three_way_islet$interaction %>% mutate(across(where(is.numeric), fmt)), startRow = 13)

  # Save
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0("Ratio_", type_tag, "_Analysis.xlsx"))
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Saved to:", output_file, "\n")
  return(invisible(wb))
}


# =============================================================================
# 3. DIAGNOSIS MODEL — SUPPLEMENTARY EXPLORATORY (ratio_mid)
# =============================================================================
#
# beta_alpha_ratio ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
#   3-way (Diag:Region:Sex, Diag:Region:Age, etc.) +
#   4-way (Diag:Region:Sex:Age, Diag:Region:Sex:log_Islet.Cells) +
#   (1|Donor) + (1|Donor:ImageID)
#
# Sheets:
#   1.  Model Summary
#   2.  Donor Counts
#   3.  Diag x Region x Sex Cell Means
#   4.  Diagnosis within Region x Sex
#   5.  Sex within Diagnosis x Region
#   6.  Region within Diagnosis x Sex
#   7.  3-Way Diag x Sex Age Slopes
#   8.  3-Way Diag x Sex Islet Slopes
#   9.  3-Way Region x Sex Age Slopes
#  10.  3-Way Region x Sex Islet Slopes
#  11.  4-Way Age Slopes
#  12.  4-Way Islet Slopes
#  13.  Full Pairwise All 12 Cells
# =============================================================================

generate_ratio_supplementary_excel <- function(model, data, output_dir = "Ratio/Results") {

  cat("\n=============================================================================\n")
  cat("GENERATING SUPPLEMENTARY EXPLORATORY EXCEL: BETA:ALPHA RATIO\n")
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
  # DONOR COUNTS
  # -------------------------------------------------------------------------
  cat("Creating Sheet 2: Donor Counts...\n")
  donor_counts <- data %>%
    distinct(Donor, Diagnosis, Region, Sex) %>%
    count(Diagnosis, Region, Sex) %>%
    pivot_wider(names_from = Region, values_from = n)

  # -------------------------------------------------------------------------
  # PREDICTION GRIDS
  # -------------------------------------------------------------------------
  cat("Setting up prediction grids...\n")

  cat_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), Age_c = 0, log_Islet.Cells_c = 0
  )
  post_fitted <- get_posterior_fitted(model, cat_grid)
  col_names <- paste(cat_grid$Diagnosis, cat_grid$Region, cat_grid$Sex, sep = "_")
  colnames(post_fitted) <- col_names

  age_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), Age_c = c(-1, 1), log_Islet.Cells_c = 0
  )
  islet_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"), Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"), Age_c = 0, log_Islet.Cells_c = c(-1, 1)
  )
  post_age   <- get_posterior_fitted(model, age_grid)
  post_islet <- get_posterior_fitted(model, islet_grid)

  # -------------------------------------------------------------------------
  # SHEET 3: 3-WAY CELL MEANS
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: 3-Way Cell Means...\n")

  cell_means_3way <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == reg & cat_grid$Sex == s)
        draws <- post_fitted[, idx]
        summ <- summarize_posterior(draws)
        summ$Diagnosis <- diag; summ$Region <- reg; summ$Sex <- s
        cell_means_3way[[paste(diag, reg, s, sep = "_")]] <- summ
      }
    }
  }
  cell_means_3way_df <- bind_rows(cell_means_3way) %>%
    select(Diagnosis, Region, Sex, estimate, lower.HPD, upper.HPD) %>%
    rename(response = estimate)

  # -------------------------------------------------------------------------
  # SHEET 4: DIAGNOSIS WITHIN REGION x SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 4: Diagnosis within Region x Sex...\n")

  diag_within_region_sex <- list()
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      nd_d <- post_fitted[, which(cat_grid$Diagnosis == "ND" & cat_grid$Region == reg & cat_grid$Sex == s)]
      t1d_d <- post_fitted[, which(cat_grid$Diagnosis == "T1D" & cat_grid$Region == reg & cat_grid$Sex == s)]
      contrast <- compute_contrast(nd_d, t1d_d, "ND - T1D")
      contrast$Region <- reg; contrast$Sex <- s
      diag_within_region_sex[[paste(reg, s, sep = "_")]] <- contrast
    }
  }
  diag_within_region_sex_df <- bind_rows(diag_within_region_sex) %>%
    select(Region, Sex, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 5: SEX WITHIN DIAGNOSIS x REGION
  # -------------------------------------------------------------------------
  cat("Creating Sheet 5: Sex within Diagnosis x Region...\n")

  sex_within_diag_region <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      f_d <- post_fitted[, which(cat_grid$Diagnosis == diag & cat_grid$Region == reg & cat_grid$Sex == "Female")]
      m_d <- post_fitted[, which(cat_grid$Diagnosis == diag & cat_grid$Region == reg & cat_grid$Sex == "Male")]
      contrast <- compute_contrast(f_d, m_d, "Female - Male")
      contrast$Diagnosis <- diag; contrast$Region <- reg
      sex_within_diag_region[[paste(diag, reg, sep = "_")]] <- contrast
    }
  }
  sex_within_diag_region_df <- bind_rows(sex_within_diag_region) %>%
    select(Diagnosis, Region, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # SHEET 6: REGION WITHIN DIAGNOSIS x SEX
  # -------------------------------------------------------------------------
  cat("Creating Sheet 6: Region within Diagnosis x Sex...\n")

  region_within_diag_sex <- list()
  for (diag in c("ND", "T1D")) {
    for (s in c("Female", "Male")) {
      hd <- post_fitted[, which(cat_grid$Diagnosis == diag & cat_grid$Region == "Head" & cat_grid$Sex == s)]
      bd <- post_fitted[, which(cat_grid$Diagnosis == diag & cat_grid$Region == "Body" & cat_grid$Sex == s)]
      td <- post_fitted[, which(cat_grid$Diagnosis == diag & cat_grid$Region == "Tail" & cat_grid$Sex == s)]
      hb <- compute_contrast(hd, bd, "Head - Body"); hb$Diagnosis <- diag; hb$Sex <- s
      ht <- compute_contrast(hd, td, "Head - Tail"); ht$Diagnosis <- diag; ht$Sex <- s
      bt <- compute_contrast(bd, td, "Body - Tail"); bt$Diagnosis <- diag; bt$Sex <- s
      region_within_diag_sex[[paste(diag, s, "HB", sep = "_")]] <- hb
      region_within_diag_sex[[paste(diag, s, "HT", sep = "_")]] <- ht
      region_within_diag_sex[[paste(diag, s, "BT", sep = "_")]] <- bt
    }
  }
  region_within_diag_sex_df <- bind_rows(region_within_diag_sex) %>%
    select(Diagnosis, Sex, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # COMPUTE ALL 12-CELL SLOPES
  # -------------------------------------------------------------------------
  cat("Computing 4-way slopes...\n")

  age_slopes_list <- list(); islet_slopes_list <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        key <- paste(diag, reg, s, sep = "_")
        lo <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & age_grid$Sex == s & age_grid$Age_c == -1)
        hi <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & age_grid$Sex == s & age_grid$Age_c == 1)
        age_slopes_list[[key]] <- (post_age[, hi] - post_age[, lo]) / 2
        lo <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == -1)
        hi <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == 1)
        islet_slopes_list[[key]] <- (post_islet[, hi] - post_islet[, lo]) / 2
      }
    }
  }

  # -------------------------------------------------------------------------
  # 3-WAY MARGINALIZED SLOPES
  # -------------------------------------------------------------------------
  cat("Computing 3-way marginalized slopes...\n")

  # Generic 3-way slope builder: given slopes_list, compute slopes for
  # each combo of var1 x var2 (marginalized over the third factor)
  compute_3way_slopes <- function(slopes_list, var1, var1_levels, var2, var2_levels, margin_var, margin_levels) {
    slopes <- list(); draws_list <- list()
    for (v1 in var1_levels) {
      for (v2 in var2_levels) {
        margin_slopes <- list()
        for (mv in margin_levels) {
          # Build key: "Diagnosis_Region_Sex" format
          parts <- list(Diagnosis = NULL, Region = NULL, Sex = NULL)
          parts[[var1]] <- v1; parts[[var2]] <- v2; parts[[margin_var]] <- mv
          k <- paste(parts$Diagnosis, parts$Region, parts$Sex, sep = "_")
          margin_slopes[[mv]] <- slopes_list[[k]]
        }
        slope_draws <- rowMeans(do.call(cbind, margin_slopes))
        summ <- summarize_posterior(slope_draws)
        summ[[var1]] <- v1; summ[[var2]] <- v2
        summ$prob_positive <- mean(slope_draws > 0); summ$prob_negative <- mean(slope_draws < 0)
        key <- paste(v1, v2, sep = "_")
        slopes[[key]] <- list(summary = summ, draws = slope_draws)
      }
    }
    list(slopes = slopes,
         summary_df = bind_rows(lapply(slopes, `[[`, "summary")) %>%
           select(all_of(c(var1, var2)), estimate, lower.HPD, upper.HPD, prob_positive, prob_negative))
  }

  # Compute contrasts for a 3-way slope set
  compute_3way_contrasts <- function(slopes_3way, var1, var1_levels, var2, var2_levels) {
    # var1 contrasts within each var2
    v1_within_v2 <- list()
    if (length(var1_levels) == 2) {
      for (v2 in var2_levels) {
        k1 <- paste(var1_levels[1], v2, sep = "_")
        k2 <- paste(var1_levels[2], v2, sep = "_")
        # For Diagnosis: T1D - ND; for Sex: Female - Male
        label <- if (var1 == "Diagnosis") "T1D - ND" else paste(var1_levels[1], "-", var1_levels[2])
        diff <- slopes_3way$slopes[[k2]]$draws - slopes_3way$slopes[[k1]]$draws
        if (var1 == "Diagnosis") diff <- slopes_3way$slopes[[paste("T1D", v2, sep = "_")]]$draws -
            slopes_3way$slopes[[paste("ND", v2, sep = "_")]]$draws
        summ <- summarize_posterior(diff)
        summ[[var2]] <- v2; summ$contrast <- label
        summ$prob_greater_0 <- mean(diff > 0); summ$prob_less_0 <- mean(diff < 0)
        v1_within_v2[[v2]] <- summ
      }
    }

    # var2 contrasts within each var1
    v2_within_v1 <- list()
    if (length(var2_levels) == 2) {
      for (v1 in var1_levels) {
        k1 <- paste(v1, var2_levels[1], sep = "_")
        k2 <- paste(v1, var2_levels[2], sep = "_")
        label <- paste(var2_levels[1], "-", var2_levels[2])
        diff <- slopes_3way$slopes[[k1]]$draws - slopes_3way$slopes[[k2]]$draws
        summ <- summarize_posterior(diff)
        summ[[var1]] <- v1; summ$contrast <- label
        summ$prob_greater_0 <- mean(diff > 0); summ$prob_less_0 <- mean(diff < 0)
        v2_within_v1[[v1]] <- summ
      }
    } else if (length(var2_levels) == 3) {
      # Region pairwise
      for (v1 in var1_levels) {
        for (pair in list(c("Head","Body"), c("Head","Tail"), c("Body","Tail"))) {
          k1 <- paste(v1, pair[1], sep = "_")
          k2 <- paste(v1, pair[2], sep = "_")
          diff <- slopes_3way$slopes[[k1]]$draws - slopes_3way$slopes[[k2]]$draws
          label <- paste(pair[1], "-", pair[2])
          summ <- summarize_posterior(diff)
          summ[[var1]] <- v1; summ$contrast <- label
          summ$prob_greater_0 <- mean(diff > 0); summ$prob_less_0 <- mean(diff < 0)
          v2_within_v1[[paste(v1, label, sep = "_")]] <- summ
        }
      }
    }

    list(
      v1_within_v2 = if (length(v1_within_v2) > 0) bind_rows(v1_within_v2) %>%
        select(all_of(var2), contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) else NULL,
      v2_within_v1 = if (length(v2_within_v1) > 0) bind_rows(v2_within_v1) %>%
        select(all_of(var1), contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) else NULL
    )
  }

  # Diag x Sex (marginalized over Region) — Age
  age_ds <- compute_3way_slopes(age_slopes_list, "Diagnosis", c("ND","T1D"), "Sex", c("Female","Male"), "Region", c("Head","Body","Tail"))
  age_ds_contr <- compute_3way_contrasts(age_ds, "Diagnosis", c("ND","T1D"), "Sex", c("Female","Male"))

  # Diag x Sex (marginalized over Region) — Islet
  islet_ds <- compute_3way_slopes(islet_slopes_list, "Diagnosis", c("ND","T1D"), "Sex", c("Female","Male"), "Region", c("Head","Body","Tail"))
  islet_ds_contr <- compute_3way_contrasts(islet_ds, "Diagnosis", c("ND","T1D"), "Sex", c("Female","Male"))

  # Region x Sex (marginalized over Diagnosis) — Age
  age_rs <- compute_3way_slopes(age_slopes_list, "Region", c("Head","Body","Tail"), "Sex", c("Female","Male"), "Diagnosis", c("ND","T1D"))
  age_rs_contr <- compute_3way_contrasts(age_rs, "Region", c("Head","Body","Tail"), "Sex", c("Female","Male"))

  # Region x Sex (marginalized over Diagnosis) — Islet
  islet_rs <- compute_3way_slopes(islet_slopes_list, "Region", c("Head","Body","Tail"), "Sex", c("Female","Male"), "Diagnosis", c("ND","T1D"))
  islet_rs_contr <- compute_3way_contrasts(islet_rs, "Region", c("Head","Body","Tail"), "Sex", c("Female","Male"))

  # -------------------------------------------------------------------------
  # 4-WAY SLOPES
  # -------------------------------------------------------------------------
  cat("Creating 4-Way slopes...\n")

  build_4way <- function(slopes_list) {
    slopes_df <- list(); slope_draws_all <- list()
    for (diag in c("ND", "T1D")) {
      for (reg in c("Head", "Body", "Tail")) {
        for (s in c("Female", "Male")) {
          key <- paste(diag, reg, s, sep = "_")
          draws <- slopes_list[[key]]
          summ <- summarize_posterior(draws)
          summ$Diagnosis <- diag; summ$Region <- reg; summ$Sex <- s
          summ$prob_positive <- mean(draws > 0); summ$prob_negative <- mean(draws < 0)
          slopes_df[[key]] <- summ
          slope_draws_all[[key]] <- draws
        }
      }
    }
    slopes_4way_df <- bind_rows(slopes_df) %>%
      select(Diagnosis, Region, Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)

    # Diagnosis effect within Region x Sex
    diag_contr <- list()
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        diff <- slope_draws_all[[paste("T1D", reg, s, sep = "_")]] -
          slope_draws_all[[paste("ND", reg, s, sep = "_")]]
        summ <- summarize_posterior(diff)
        summ$Region <- reg; summ$Sex <- s; summ$contrast <- "T1D - ND"
        summ$prob_positive <- mean(diff > 0)
        diag_contr[[paste(reg, s, sep = "_")]] <- summ
      }
    }
    diag_contr_df <- bind_rows(diag_contr) %>%
      select(Region, Sex, contrast, estimate, lower.HPD, upper.HPD, prob_positive)

    # Sex effect within Diagnosis x Region
    sex_contr <- list()
    for (diag in c("ND", "T1D")) {
      for (reg in c("Head", "Body", "Tail")) {
        diff <- slope_draws_all[[paste(diag, reg, "Female", sep = "_")]] -
          slope_draws_all[[paste(diag, reg, "Male", sep = "_")]]
        summ <- summarize_posterior(diff)
        summ$Diagnosis <- diag; summ$Region <- reg; summ$contrast <- "Female - Male"
        summ$prob_positive <- mean(diff > 0)
        sex_contr[[paste(diag, reg, sep = "_")]] <- summ
      }
    }
    sex_contr_df <- bind_rows(sex_contr) %>%
      select(Diagnosis, Region, contrast, estimate, lower.HPD, upper.HPD, prob_positive)

    list(slopes = slopes_4way_df, diag_contr = diag_contr_df, sex_contr = sex_contr_df)
  }

  age_4way <- build_4way(age_slopes_list)
  islet_4way <- build_4way(islet_slopes_list)

  # -------------------------------------------------------------------------
  # FULL PAIRWISE (ALL 12 CELLS)
  # -------------------------------------------------------------------------
  cat("Creating Full Pairwise (66 comparisons)...\n")

  all_cells <- expand.grid(Diagnosis = c("ND","T1D"), Region = c("Head","Body","Tail"),
                           Sex = c("Female","Male"), stringsAsFactors = FALSE) %>%
    mutate(cell_name = paste(Diagnosis, Region, Sex))

  full_pairwise <- list()
  for (i in 1:(nrow(all_cells) - 1)) {
    for (j in (i + 1):nrow(all_cells)) {
      c1 <- all_cells[i, ]; c2 <- all_cells[j, ]
      idx1 <- which(cat_grid$Diagnosis == c1$Diagnosis & cat_grid$Region == c1$Region & cat_grid$Sex == c1$Sex)
      idx2 <- which(cat_grid$Diagnosis == c2$Diagnosis & cat_grid$Region == c2$Region & cat_grid$Sex == c2$Sex)
      label <- paste0(c1$cell_name, " - ", c2$cell_name)
      full_pairwise[[label]] <- compute_contrast(post_fitted[, idx1], post_fitted[, idx2], label)
    }
  }
  full_pairwise_df <- bind_rows(full_pairwise) %>%
    select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)

  # -------------------------------------------------------------------------
  # CREATE WORKBOOK
  # -------------------------------------------------------------------------
  cat("\nCreating Supplementary Excel workbook...\n")
  wb <- createWorkbook()

  # Sheet 1: Model Summary
  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", data.frame(
    Note = c("SUPPLEMENTARY EXPLORATORY ANALYSIS: Beta:Alpha Ratio",
             "", "CAUTION: This analysis is EXPLORATORY due to limited sample sizes in some cells.",
             "Results should be interpreted with caution, particularly for cells with <20 donors.",
             "", "Bayesian ordered beta regression (ordbetareg)",
             "Response: beta_alpha_ratio = beta / (alpha + beta), 0-1 scale",
             "Model includes 3-way and 4-way interactions",
             "", "All estimates are posterior medians with 95% HPD intervals.",
             "prob_greater_0 / prob_less_0 = posterior probability contrast is positive/negative.",
             "")))
  writeData(wb, "Model Summary", model_summary, startRow = 14)

  # Sheet 2: Donor Counts
  addWorksheet(wb, "Donor Counts")
  writeData(wb, "Donor Counts", data.frame(
    Note = c("Donor counts by cell. Cells with <20 donors have limited power.", "")))
  writeData(wb, "Donor Counts", donor_counts, startRow = 4)

  # Sheet 3: 3-Way Cell Means
  addWorksheet(wb, "Diag x Reg x Sex Means")
  writeData(wb, "Diag x Reg x Sex Means", data.frame(
    Note = "EXPLORATORY: 3-way cell means for Diagnosis x Region x Sex."))
  writeData(wb, "Diag x Reg x Sex Means", cell_means_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)

  # Sheet 4
  addWorksheet(wb, "Diag within Reg x Sex")
  writeData(wb, "Diag within Reg x Sex", data.frame(
    Note = "EXPLORATORY: Diagnosis effect (ND - T1D) within each Region x Sex cell."))
  writeData(wb, "Diag within Reg x Sex", diag_within_region_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)

  # Sheet 5
  addWorksheet(wb, "Sex within Diag x Reg")
  writeData(wb, "Sex within Diag x Reg", data.frame(
    Note = "EXPLORATORY: Sex effect (Female - Male) within each Diagnosis x Region cell."))
  writeData(wb, "Sex within Diag x Reg", sex_within_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)

  # Sheet 6
  addWorksheet(wb, "Region within Diag x Sex")
  writeData(wb, "Region within Diag x Sex", data.frame(
    Note = "EXPLORATORY: Region contrasts within each Diagnosis x Sex cell."))
  writeData(wb, "Region within Diag x Sex", region_within_diag_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)

  # Helper for 3-way slope sheets
  write_3way_slope_sheet <- function(wb, sheet_name, note, slopes_df, v1_contr, v1_label, v2_contr, v2_label) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data.frame(Note = note))
    writeData(wb, sheet_name, slopes_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
    r <- nrow(slopes_df) + 4
    if (!is.null(v1_contr)) {
      writeData(wb, sheet_name, data.frame(Section = v1_label), startRow = r); r <- r + 1
      writeData(wb, sheet_name, v1_contr %>% mutate(across(where(is.numeric), fmt)), startRow = r)
      r <- r + nrow(v1_contr) + 1
    }
    if (!is.null(v2_contr)) {
      writeData(wb, sheet_name, data.frame(Section = v2_label), startRow = r); r <- r + 1
      writeData(wb, sheet_name, v2_contr %>% mutate(across(where(is.numeric), fmt)), startRow = r)
    }
  }

  # Sheet 7: 3-Way Diag x Sex Age Slopes
  write_3way_slope_sheet(wb, "3Way DiagxSex Age Slopes",
    "EXPLORATORY: Diagnosis x Sex age slopes (marginalized over Region).",
    age_ds$summary_df,
    age_ds_contr$v1_within_v2, "Diagnosis effect (T1D - ND) on age slopes within each Sex:",
    age_ds_contr$v2_within_v1, "Sex effect (Female - Male) on age slopes within each Diagnosis:")

  # Sheet 8: 3-Way Diag x Sex Islet Slopes
  write_3way_slope_sheet(wb, "3Way DiagxSex Islet Slopes",
    "EXPLORATORY: Diagnosis x Sex islet size slopes (marginalized over Region).",
    islet_ds$summary_df,
    islet_ds_contr$v1_within_v2, "Diagnosis effect (T1D - ND) on islet slopes within each Sex:",
    islet_ds_contr$v2_within_v1, "Sex effect (Female - Male) on islet slopes within each Diagnosis:")

  # Sheet 9: 3-Way Region x Sex Age Slopes
  write_3way_slope_sheet(wb, "3Way RegxSex Age Slopes",
    "EXPLORATORY: Region x Sex age slopes (marginalized over Diagnosis).",
    age_rs$summary_df,
    age_rs_contr$v2_within_v1, "Sex effect (Female - Male) on age slopes within each Region:",
    age_rs_contr$v1_within_v2, "Region pairwise contrasts on age slopes within each Sex:")

  # Sheet 10: 3-Way Region x Sex Islet Slopes
  write_3way_slope_sheet(wb, "3Way RegxSex Islet Slopes",
    "EXPLORATORY: Region x Sex islet size slopes (marginalized over Diagnosis).",
    islet_rs$summary_df,
    islet_rs_contr$v2_within_v1, "Sex effect (Female - Male) on islet slopes within each Region:",
    islet_rs_contr$v1_within_v2, "Region pairwise contrasts on islet slopes within each Sex:")

  # Sheet 11: 4-Way Age Slopes
  addWorksheet(wb, "4-Way Age Slopes")
  writeData(wb, "4-Way Age Slopes", data.frame(
    Note = "EXPLORATORY: Age slopes for each Diagnosis x Region x Sex cell."))
  writeData(wb, "4-Way Age Slopes", age_4way$slopes %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "4-Way Age Slopes", data.frame(Section = "Diagnosis effect (T1D - ND) on Age slopes within each Region x Sex:"), startRow = 16)
  writeData(wb, "4-Way Age Slopes", age_4way$diag_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 17)
  writeData(wb, "4-Way Age Slopes", data.frame(Section = "Sex effect (Female - Male) on Age slopes within each Diagnosis x Region:"), startRow = 24)
  writeData(wb, "4-Way Age Slopes", age_4way$sex_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 25)

  # Sheet 12: 4-Way Islet Slopes
  addWorksheet(wb, "4-Way Islet Slopes")
  writeData(wb, "4-Way Islet Slopes", data.frame(
    Note = "EXPLORATORY: Islet size slopes for each Diagnosis x Region x Sex cell."))
  writeData(wb, "4-Way Islet Slopes", islet_4way$slopes %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "4-Way Islet Slopes", data.frame(Section = "Diagnosis effect (T1D - ND) on Islet slopes within each Region x Sex:"), startRow = 16)
  writeData(wb, "4-Way Islet Slopes", islet_4way$diag_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 17)
  writeData(wb, "4-Way Islet Slopes", data.frame(Section = "Sex effect (Female - Male) on Islet slopes within each Diagnosis x Region:"), startRow = 24)
  writeData(wb, "4-Way Islet Slopes", islet_4way$sex_contr %>% mutate(across(where(is.numeric), fmt)), startRow = 25)

  # Sheet 13: Full Pairwise
  addWorksheet(wb, "Full Pairwise All Cells")
  writeData(wb, "Full Pairwise All Cells", data.frame(
    Note = c("EXPLORATORY: All pairwise comparisons between the 12 Diagnosis x Region x Sex cells.",
             "66 total comparisons. Posterior probabilities do not require multiple comparison corrections.",
             "")))
  writeData(wb, "Full Pairwise All Cells", full_pairwise_df %>% mutate(across(where(is.numeric), fmt)), startRow = 5)

  # Save
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, "Ratio_Supplementary_Exploratory.xlsx")
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Saved to:", output_file, "\n")
  return(invisible(wb))
}


# =============================================================================
# MAIN EXECUTION
# =============================================================================

generate_all_ratio_excel <- function(output_dir = "Ratio/Results") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL RATIO EXCEL RESULTS\n")
  cat("=============================================================================\n\n")

  # Load data
  cat("Loading data...\n")
  ratio_data     <- readRDS("Data/ratio_data.rds")
  ratio_data_t1d <- readRDS("Data/ratio_data_t1d.rds")
  
  ratio_data <- ratio_data %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
      Region    = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor),
      ImageID   = factor(ImageID)
    )
  ratio_data_t1d <- ratio_data_t1d %>%
    mutate(
      Region  = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex     = factor(Sex, levels = c("Female", "Male")),
      Donor   = factor(Donor),
      ImageID = factor(ImageID)
    ) %>% droplevels()
  
  cat("ratio_data:", nrow(ratio_data), "rows,", n_distinct(ratio_data$Donor), "donors\n")
  cat("ratio_data_t1d:", nrow(ratio_data_t1d), "rows,", n_distinct(ratio_data_t1d$Donor), "donors\n\n")
  
  # --- Diagnosis (primary) ---
  cat("Loading ratio_powered model...\n")
  ratio_powered <- readRDS("Ratio/Models/ratio_powered_model.rds")
  generate_ratio_diag_excel(ratio_powered, ratio_data, output_dir)
  rm(ratio_powered); gc()
  
  # --- Diagnosis (supplementary) ---
  cat("\nLoading ratio_mid model...\n")
  ratio_mid <- readRDS("Ratio/Models/ratio_mid_model.rds")
  generate_ratio_supplementary_excel(ratio_mid, ratio_data, output_dir)
  rm(ratio_mid); gc()
  
  # --- T1D (AO parameterization only — DD is algebraically equivalent) ---
  cat("\nLoading ratio_ao model...\n")
  ratio_ao <- readRDS("Ratio/Models/ratio_ao_model.rds")
  generate_ratio_t1d_excel(ratio_ao, ratio_data_t1d, model_type = "ao", output_dir = output_dir)
  rm(ratio_ao); gc()
  
  cat("\n=============================================================================\n")
  cat("COMPLETE! Three Excel files created in", output_dir, ":\n")
  cat("  1. Ratio_Primary_Analysis.xlsx              (Diagnosis — powered model)\n")
  cat("  2. Ratio_Supplementary_Exploratory.xlsx     (Diagnosis — mid model)\n")
  cat("  3. Ratio_AO_Analysis.xlsx                   (T1D — AO parameterization)\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all ratio Excel files, run:\n")
  cat("  generate_all_ratio_excel()\n")
}
