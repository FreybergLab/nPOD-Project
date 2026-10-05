# =============================================================================
# Generate Comprehensive Excel Results for ALL Cell Types
# =============================================================================
# 
# This script generates Excel workbooks for all four cell types:
#   - Insulin (ins)
#   - Glucagon (glu)
#   - Somatostatin (soma)
#   - Pancreatic Polypeptide (pp)
#
# For each cell type, TWO Excel files are created:
#   1. Primary Analysis (*_powered model) - Confirmatory results
#   2. Supplementary Analysis (*_mid model) - Exploratory 3-way/4-way interactions
#
# Total output: 8 Excel files
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
# MAIN ANALYSIS FUNCTION (for *_powered models)
# =============================================================================

generate_primary_excel <- function(model, data, cell_type, cell_label, output_dir = "Proportion/Results") {
  
  cat("\n")
  cat("=============================================================================\n")
  cat("GENERATING PRIMARY ANALYSIS EXCEL:", toupper(cell_label), "\n
")
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
  # CATEGORICAL GRID SETUP
  # -------------------------------------------------------------------------
  cat("Setting up prediction grid and getting posterior draws...\n")
  
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
  # SHEET 2: DIAGNOSIS MAIN EFFECT
  # -------------------------------------------------------------------------
  cat("Creating Sheet 2: Diagnosis...\n")
  
  nd_idx <- which(cat_grid$Diagnosis == "ND")
  t1d_idx <- which(cat_grid$Diagnosis == "T1D")
  
  nd_draws <- rowMeans(post_fitted[, nd_idx])
  t1d_draws <- rowMeans(post_fitted[, t1d_idx])
  
  diagnosis_marginal <- bind_rows(
    summarize_posterior(nd_draws) %>% mutate(Diagnosis = "ND", .before = 1),
    summarize_posterior(t1d_draws) %>% mutate(Diagnosis = "T1D", .before = 1)
  ) %>% rename(response = estimate)
  
  diagnosis_contrast <- compute_contrast(nd_draws, t1d_draws, "ND - T1D")
  
  # -------------------------------------------------------------------------
  # SHEET 3: REGION MAIN EFFECT
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: Region...\n")
  
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
  # SHEET 4: SEX MAIN EFFECT
  # -------------------------------------------------------------------------
  cat("Creating Sheet 4: Sex...\n")
  
  female_idx <- which(cat_grid$Sex == "Female")
  male_idx <- which(cat_grid$Sex == "Male")
  
  female_draws <- rowMeans(post_fitted[, female_idx])
  male_draws <- rowMeans(post_fitted[, male_idx])
  
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
      summ$Diagnosis <- diag
      summ$Region <- reg
      diag_region_marginal[[paste(diag, reg, sep = "_")]] <- summ
    }
  }
  diag_region_marginal_df <- bind_rows(diag_region_marginal) %>%
    select(Diagnosis, Region, estimate, lower.HPD, upper.HPD) %>%
    rename(response = estimate)
  
  # Diagnosis effect within each Region
  diag_within_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    nd_idx <- which(cat_grid$Diagnosis == "ND" & cat_grid$Region == reg)
    t1d_idx <- which(cat_grid$Diagnosis == "T1D" & cat_grid$Region == reg)
    nd_draws <- rowMeans(post_fitted[, nd_idx, drop = FALSE])
    t1d_draws <- rowMeans(post_fitted[, t1d_idx, drop = FALSE])
    contrast <- compute_contrast(nd_draws, t1d_draws, "ND - T1D")
    contrast$Region <- reg
    diag_within_region[[reg]] <- contrast
  }
  diag_within_region_df <- bind_rows(diag_within_region) %>%
    select(contrast, Region, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Region effects within each Diagnosis
  region_within_diag <- list()
  for (d in c("ND", "T1D")) {
    head_idx <- which(cat_grid$Region == "Head" & cat_grid$Diagnosis == d)
    body_idx <- which(cat_grid$Region == "Body" & cat_grid$Diagnosis == d)
    tail_idx <- which(cat_grid$Region == "Tail" & cat_grid$Diagnosis == d)
    
    head_draws <- rowMeans(post_fitted[, head_idx, drop = FALSE])
    body_draws <- rowMeans(post_fitted[, body_idx, drop = FALSE])
    tail_draws <- rowMeans(post_fitted[, tail_idx, drop = FALSE])
    
    hb <- compute_contrast(head_draws, body_draws, "Head - Body")
    hb$Diagnosis <- d
    ht <- compute_contrast(head_draws, tail_draws, "Head - Tail")
    ht$Diagnosis <- d
    bt <- compute_contrast(body_draws, tail_draws, "Body - Tail")
    bt$Diagnosis <- d
    
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
      summ$Diagnosis <- diag
      summ$Sex <- s
      diag_sex_marginal[[paste(diag, s, sep = "_")]] <- summ
    }
  }
  diag_sex_marginal_df <- bind_rows(diag_sex_marginal) %>%
    select(Diagnosis, Sex, estimate, lower.HPD, upper.HPD) %>%
    rename(response = estimate)
  
  # Diagnosis effect within each Sex
  diag_within_sex <- list()
  for (s in c("Female", "Male")) {
    nd_idx <- which(cat_grid$Diagnosis == "ND" & cat_grid$Sex == s)
    t1d_idx <- which(cat_grid$Diagnosis == "T1D" & cat_grid$Sex == s)
    nd_draws <- rowMeans(post_fitted[, nd_idx, drop = FALSE])
    t1d_draws <- rowMeans(post_fitted[, t1d_idx, drop = FALSE])
    contrast <- compute_contrast(nd_draws, t1d_draws, "ND - T1D")
    contrast$Sex <- s
    diag_within_sex[[s]] <- contrast
  }
  diag_within_sex_df <- bind_rows(diag_within_sex) %>%
    select(contrast, Sex, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Sex effect within each Diagnosis
  sex_within_diag <- list()
  for (d in c("ND", "T1D")) {
    f_idx <- which(cat_grid$Sex == "Female" & cat_grid$Diagnosis == d)
    m_idx <- which(cat_grid$Sex == "Male" & cat_grid$Diagnosis == d)
    f_draws <- rowMeans(post_fitted[, f_idx, drop = FALSE])
    m_draws <- rowMeans(post_fitted[, m_idx, drop = FALSE])
    contrast <- compute_contrast(f_draws, m_draws, "Female - Male")
    contrast$Diagnosis <- d
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
      summ$Region <- reg
      summ$Sex <- s
      region_sex_marginal[[paste(reg, s, sep = "_")]] <- summ
    }
  }
  region_sex_marginal_df <- bind_rows(region_sex_marginal) %>%
    select(Region, Sex, estimate, lower.HPD, upper.HPD) %>%
    rename(response = estimate)
  
  # Sex effect within each Region
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
  
  # Region effects within each Sex
  region_within_sex <- list()
  for (s in c("Female", "Male")) {
    head_idx <- which(cat_grid$Region == "Head" & cat_grid$Sex == s)
    body_idx <- which(cat_grid$Region == "Body" & cat_grid$Sex == s)
    tail_idx <- which(cat_grid$Region == "Tail" & cat_grid$Sex == s)
    
    head_draws <- rowMeans(post_fitted[, head_idx, drop = FALSE])
    body_draws <- rowMeans(post_fitted[, body_idx, drop = FALSE])
    tail_draws <- rowMeans(post_fitted[, tail_idx, drop = FALSE])
    
    hb <- compute_contrast(head_draws, body_draws, "Head - Body")
    hb$Sex <- s
    ht <- compute_contrast(head_draws, tail_draws, "Head - Tail")
    ht$Sex <- s
    bt <- compute_contrast(body_draws, tail_draws, "Body - Tail")
    bt$Sex <- s
    
    region_within_sex[[paste0(s, "_HB")]] <- hb
    region_within_sex[[paste0(s, "_HT")]] <- ht
    region_within_sex[[paste0(s, "_BT")]] <- bt
  }
  region_within_sex_df <- bind_rows(region_within_sex) %>%
    select(contrast, Sex, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # CONTINUOUS PREDICTOR SLOPES SETUP
  # -------------------------------------------------------------------------
  cat("Setting up continuous predictor grids...\n")
  
  age_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Age_c = c(-1, 1),
    log_Islet.Cells_c = 0
  )
  
  islet_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Age_c = 0,
    log_Islet.Cells_c = c(-1, 1)
  )
  
  cat("Getting posterior predictions for slopes...\n")
  post_age <- get_posterior_fitted(model, age_grid)
  post_islet <- get_posterior_fitted(model, islet_grid)
  
  # Store all cell-level slopes
  age_slopes_list <- list()
  islet_slopes_list <- list()
  
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        # Age slopes
        low_idx <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & 
                           age_grid$Sex == s & age_grid$Age_c == -1)
        high_idx <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & 
                            age_grid$Sex == s & age_grid$Age_c == 1)
        age_slopes_list[[paste(diag, reg, s, sep = "_")]] <- (post_age[, high_idx] - post_age[, low_idx]) / 2
        
        # Islet slopes
        low_idx <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & 
                           islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == -1)
        high_idx <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & 
                            islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == 1)
        islet_slopes_list[[paste(diag, reg, s, sep = "_")]] <- (post_islet[, high_idx] - post_islet[, low_idx]) / 2
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # SHEET 8: AGE SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 8: Age Slopes...\n")
  
  # Overall Age slope
  all_age_slopes <- do.call(cbind, age_slopes_list)
  overall_age <- rowMeans(all_age_slopes)
  overall_age_summ <- summarize_posterior(overall_age)
  overall_age_summ$predictor <- "Age_c"
  overall_age_summ$prob_positive <- mean(overall_age > 0)
  overall_age_summ$prob_negative <- mean(overall_age < 0)
  
  # Age slopes by Diagnosis
  age_by_diag <- list()
  for (diag in c("ND", "T1D")) {
    slopes_subset <- list()
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        slopes_subset[[paste(reg, s, sep = "_")]] <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
    }
    slope_draws <- rowMeans(do.call(cbind, slopes_subset))
    summ <- summarize_posterior(slope_draws)
    summ$Diagnosis <- diag
    summ$prob_positive <- mean(slope_draws > 0)
    summ$prob_negative <- mean(slope_draws < 0)
    age_by_diag[[diag]] <- list(summary = summ, draws = slope_draws)
  }
  
  age_by_diag_df <- bind_rows(
    age_by_diag[["ND"]]$summary,
    age_by_diag[["T1D"]]$summary
  ) %>% select(Diagnosis, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # Age slope difference by Diagnosis
  age_diag_diff <- age_by_diag[["T1D"]]$draws - age_by_diag[["ND"]]$draws
  age_diag_diff_summ <- summarize_posterior(age_diag_diff)
  age_diag_diff_summ$contrast <- "T1D - ND"
  age_diag_diff_summ$prob_positive <- mean(age_diag_diff > 0)
  age_diag_diff_summ$prob_negative <- mean(age_diag_diff < 0)
  
  # Age slopes by Region
  age_by_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    slopes_subset <- list()
    for (diag in c("ND", "T1D")) {
      for (s in c("Female", "Male")) {
        slopes_subset[[paste(diag, s, sep = "_")]] <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
    }
    slope_draws <- rowMeans(do.call(cbind, slopes_subset))
    summ <- summarize_posterior(slope_draws)
    summ$Region <- reg
    summ$prob_positive <- mean(slope_draws > 0)
    summ$prob_negative <- mean(slope_draws < 0)
    age_by_region[[reg]] <- summ
  }
  age_by_region_df <- bind_rows(age_by_region) %>%
    select(Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # Age region contrasts
  age_region_contrasts <- bind_rows(
    {d <- age_by_region[["Head"]]$estimate - age_by_region[["Body"]]$estimate; 
    tibble(contrast = "Head - Body", estimate = d)},  # Placeholder - need draws
    {d <- age_by_region[["Head"]]$estimate - age_by_region[["Tail"]]$estimate;
    tibble(contrast = "Head - Tail", estimate = d)},
    {d <- age_by_region[["Body"]]$estimate - age_by_region[["Tail"]]$estimate;
    tibble(contrast = "Body - Tail", estimate = d)}
  )
  
  # Recompute with proper draws for region contrasts
  age_by_region_draws <- list()
  for (reg in c("Head", "Body", "Tail")) {
    slopes_subset <- list()
    for (diag in c("ND", "T1D")) {
      for (s in c("Female", "Male")) {
        slopes_subset[[paste(diag, s, sep = "_")]] <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
    }
    age_by_region_draws[[reg]] <- rowMeans(do.call(cbind, slopes_subset))
  }
  
  age_region_contrasts_df <- bind_rows(
    compute_contrast(age_by_region_draws[["Head"]], age_by_region_draws[["Body"]], "Head - Body"),
    compute_contrast(age_by_region_draws[["Head"]], age_by_region_draws[["Tail"]], "Head - Tail"),
    compute_contrast(age_by_region_draws[["Body"]], age_by_region_draws[["Tail"]], "Body - Tail")
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Age slopes by Sex
  age_by_sex <- list()
  for (s in c("Female", "Male")) {
    slopes_subset <- list()
    for (diag in c("ND", "T1D")) {
      for (reg in c("Head", "Body", "Tail")) {
        slopes_subset[[paste(diag, reg, sep = "_")]] <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
    }
    slope_draws <- rowMeans(do.call(cbind, slopes_subset))
    summ <- summarize_posterior(slope_draws)
    summ$Sex <- s
    summ$prob_positive <- mean(slope_draws > 0)
    summ$prob_negative <- mean(slope_draws < 0)
    age_by_sex[[s]] <- list(summary = summ, draws = slope_draws)
  }
  age_by_sex_df <- bind_rows(
    age_by_sex[["Female"]]$summary,
    age_by_sex[["Male"]]$summary
  ) %>% select(Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  age_sex_contrast_df <- compute_contrast(
    age_by_sex[["Female"]]$draws, 
    age_by_sex[["Male"]]$draws, 
    "Female - Male"
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # SHEET 9: ISLET SIZE SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 9: Islet Size Slopes...\n")
  
  # Overall Islet slope
  all_islet_slopes <- do.call(cbind, islet_slopes_list)
  overall_islet <- rowMeans(all_islet_slopes)
  overall_islet_summ <- summarize_posterior(overall_islet)
  overall_islet_summ$predictor <- "log_Islet.Cells_c"
  overall_islet_summ$prob_positive <- mean(overall_islet > 0)
  overall_islet_summ$prob_negative <- mean(overall_islet < 0)
  
  # Islet slopes by Diagnosis
  islet_by_diag <- list()
  for (diag in c("ND", "T1D")) {
    slopes_subset <- list()
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        slopes_subset[[paste(reg, s, sep = "_")]] <- islet_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
    }
    slope_draws <- rowMeans(do.call(cbind, slopes_subset))
    summ <- summarize_posterior(slope_draws)
    summ$Diagnosis <- diag
    summ$prob_positive <- mean(slope_draws > 0)
    summ$prob_negative <- mean(slope_draws < 0)
    islet_by_diag[[diag]] <- list(summary = summ, draws = slope_draws)
  }
  
  islet_by_diag_df <- bind_rows(
    islet_by_diag[["ND"]]$summary,
    islet_by_diag[["T1D"]]$summary
  ) %>% select(Diagnosis, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  islet_diag_diff <- islet_by_diag[["T1D"]]$draws - islet_by_diag[["ND"]]$draws
  islet_diag_diff_summ <- summarize_posterior(islet_diag_diff)
  islet_diag_diff_summ$contrast <- "T1D - ND"
  islet_diag_diff_summ$prob_positive <- mean(islet_diag_diff > 0)
  islet_diag_diff_summ$prob_negative <- mean(islet_diag_diff < 0)
  
  # Islet slopes by Region
  islet_by_region <- list()
  islet_by_region_draws <- list()
  for (reg in c("Head", "Body", "Tail")) {
    slopes_subset <- list()
    for (diag in c("ND", "T1D")) {
      for (s in c("Female", "Male")) {
        slopes_subset[[paste(diag, s, sep = "_")]] <- islet_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
    }
    slope_draws <- rowMeans(do.call(cbind, slopes_subset))
    summ <- summarize_posterior(slope_draws)
    summ$Region <- reg
    summ$prob_positive <- mean(slope_draws > 0)
    summ$prob_negative <- mean(slope_draws < 0)
    islet_by_region[[reg]] <- summ
    islet_by_region_draws[[reg]] <- slope_draws
  }
  islet_by_region_df <- bind_rows(islet_by_region) %>%
    select(Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  islet_region_contrasts_df <- bind_rows(
    compute_contrast(islet_by_region_draws[["Head"]], islet_by_region_draws[["Body"]], "Head - Body"),
    compute_contrast(islet_by_region_draws[["Head"]], islet_by_region_draws[["Tail"]], "Head - Tail"),
    compute_contrast(islet_by_region_draws[["Body"]], islet_by_region_draws[["Tail"]], "Body - Tail")
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Islet slopes by Sex
  islet_by_sex <- list()
  for (s in c("Female", "Male")) {
    slopes_subset <- list()
    for (diag in c("ND", "T1D")) {
      for (reg in c("Head", "Body", "Tail")) {
        slopes_subset[[paste(diag, reg, sep = "_")]] <- islet_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
    }
    slope_draws <- rowMeans(do.call(cbind, slopes_subset))
    summ <- summarize_posterior(slope_draws)
    summ$Sex <- s
    summ$prob_positive <- mean(slope_draws > 0)
    summ$prob_negative <- mean(slope_draws < 0)
    islet_by_sex[[s]] <- list(summary = summ, draws = slope_draws)
  }
  islet_by_sex_df <- bind_rows(
    islet_by_sex[["Female"]]$summary,
    islet_by_sex[["Male"]]$summary
  ) %>% select(Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  islet_sex_contrast_df <- compute_contrast(
    islet_by_sex[["Female"]]$draws, 
    islet_by_sex[["Male"]]$draws, 
    "Female - Male"
  ) %>% select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # SHEET 10: 3-WAY INTERACTION - Diag x Region x Age
  # -------------------------------------------------------------------------
  cat("Creating Sheet 10: Diag x Region x Age (3-way)...\n")
  
  age_diag_region <- list()
  age_diag_region_draws <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      slopes_subset <- list()
      for (s in c("Female", "Male")) {
        slopes_subset[[s]] <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
      slope_draws <- rowMeans(do.call(cbind, slopes_subset))
      summ <- summarize_posterior(slope_draws)
      summ$Diagnosis <- diag
      summ$Region <- reg
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      key <- paste(diag, reg, sep = "_")
      age_diag_region[[key]] <- summ
      age_diag_region_draws[[key]] <- slope_draws
    }
  }
  age_diag_region_df <- bind_rows(age_diag_region) %>%
    select(Diagnosis, Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # Age slope difference (T1D - ND) by Region
  age_slope_diff_by_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff_draws <- age_diag_region_draws[[paste("T1D", reg, sep = "_")]] - 
      age_diag_region_draws[[paste("ND", reg, sep = "_")]]
    summ <- summarize_posterior(diff_draws)
    summ$Region <- reg
    summ$contrast <- "T1D - ND"
    summ$prob_positive <- mean(diff_draws > 0)
    summ$prob_negative <- mean(diff_draws < 0)
    age_slope_diff_by_region[[reg]] <- summ
  }
  age_slope_diff_by_region_df <- bind_rows(age_slope_diff_by_region) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # -------------------------------------------------------------------------
  # SHEET 11: 3-WAY INTERACTION - Diag x Region x Islet
  # -------------------------------------------------------------------------
  cat("Creating Sheet 11: Diag x Region x Islet (3-way)...\n")
  
  islet_diag_region <- list()
  islet_diag_region_draws <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      slopes_subset <- list()
      for (s in c("Female", "Male")) {
        slopes_subset[[s]] <- islet_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
      slope_draws <- rowMeans(do.call(cbind, slopes_subset))
      summ <- summarize_posterior(slope_draws)
      summ$Diagnosis <- diag
      summ$Region <- reg
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      key <- paste(diag, reg, sep = "_")
      islet_diag_region[[key]] <- summ
      islet_diag_region_draws[[key]] <- slope_draws
    }
  }
  islet_diag_region_df <- bind_rows(islet_diag_region) %>%
    select(Diagnosis, Region, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # Islet slope difference (T1D - ND) by Region
  islet_slope_diff_by_region <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff_draws <- islet_diag_region_draws[[paste("T1D", reg, sep = "_")]] - 
      islet_diag_region_draws[[paste("ND", reg, sep = "_")]]
    summ <- summarize_posterior(diff_draws)
    summ$Region <- reg
    summ$contrast <- "T1D - ND"
    summ$prob_positive <- mean(diff_draws > 0)
    summ$prob_negative <- mean(diff_draws < 0)
    islet_slope_diff_by_region[[reg]] <- summ
  }
  islet_slope_diff_by_region_df <- bind_rows(islet_slope_diff_by_region) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # Region contrasts for Islet slopes within each Diagnosis
  islet_region_within_diag <- list()
  for (d in c("ND", "T1D")) {
    head_draws <- islet_diag_region_draws[[paste(d, "Head", sep = "_")]]
    body_draws <- islet_diag_region_draws[[paste(d, "Body", sep = "_")]]
    tail_draws <- islet_diag_region_draws[[paste(d, "Tail", sep = "_")]]
    
    hb <- compute_contrast(head_draws, body_draws, "Head - Body")
    hb$Diagnosis <- d
    ht <- compute_contrast(head_draws, tail_draws, "Head - Tail")
    ht$Diagnosis <- d
    bt <- compute_contrast(body_draws, tail_draws, "Body - Tail")
    bt$Diagnosis <- d
    
    islet_region_within_diag[[paste0(d, "_HB")]] <- hb
    islet_region_within_diag[[paste0(d, "_HT")]] <- ht
    islet_region_within_diag[[paste0(d, "_BT")]] <- bt
  }
  islet_region_within_diag_df <- bind_rows(islet_region_within_diag) %>%
    select(Diagnosis, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # SHEETS 12-13: SUPPLEMENTARY PAIRWISE COMPARISONS
  # -------------------------------------------------------------------------
  cat("Creating Sheets 12-13: Supplementary Pairwise Comparisons...\n")
  
  # Pairwise: Diagnosis x Sex
  pairwise_diag_sex <- list()
  combos <- list(
    c("ND", "Female", "T1D", "Female"),
    c("ND", "Female", "ND", "Male"),
    c("ND", "Female", "T1D", "Male"),
    c("T1D", "Female", "ND", "Male"),
    c("T1D", "Female", "T1D", "Male"),
    c("ND", "Male", "T1D", "Male")
  )
  
  for (combo in combos) {
    idx1 <- which(cat_grid$Diagnosis == combo[1] & cat_grid$Sex == combo[2])
    idx2 <- which(cat_grid$Diagnosis == combo[3] & cat_grid$Sex == combo[4])
    draws1 <- rowMeans(post_fitted[, idx1, drop = FALSE])
    draws2 <- rowMeans(post_fitted[, idx2, drop = FALSE])
    label <- paste0(combo[1], " ", combo[2], " - ", combo[3], " ", combo[4])
    contrast <- compute_contrast(draws1, draws2, label)
    pairwise_diag_sex[[label]] <- contrast
  }
  pairwise_diag_sex_df <- bind_rows(pairwise_diag_sex) %>%
    select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Pairwise: Diagnosis x Region
  pairwise_diag_region <- list()
  diag_region_combos <- expand.grid(
    d1 = c("ND", "T1D"),
    r1 = c("Head", "Body", "Tail"),
    d2 = c("ND", "T1D"),
    r2 = c("Head", "Body", "Tail"),
    stringsAsFactors = FALSE
  )
  diag_region_combos <- diag_region_combos %>%
    mutate(combo1 = paste(d1, r1), combo2 = paste(d2, r2)) %>%
    filter(combo1 < combo2) %>%
    select(d1, r1, d2, r2)
  
  for (i in 1:nrow(diag_region_combos)) {
    row <- diag_region_combos[i, ]
    idx1 <- which(cat_grid$Diagnosis == row$d1 & cat_grid$Region == row$r1)
    idx2 <- which(cat_grid$Diagnosis == row$d2 & cat_grid$Region == row$r2)
    draws1 <- rowMeans(post_fitted[, idx1, drop = FALSE])
    draws2 <- rowMeans(post_fitted[, idx2, drop = FALSE])
    label <- paste0(row$d1, " ", row$r1, " - ", row$d2, " ", row$r2)
    contrast <- compute_contrast(draws1, draws2, label)
    pairwise_diag_region[[label]] <- contrast
  }
  pairwise_diag_region_df <- bind_rows(pairwise_diag_region) %>%
    select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # CREATE EXCEL WORKBOOK
  # -------------------------------------------------------------------------
  cat("\nCreating Primary Analysis Excel workbook...\n")
  
  wb <- createWorkbook()
  
  # Sheet 1: Model Summary
  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", data.frame(
    Note = c(paste("PRIMARY ANALYSIS:", cell_label, "Content"),
             "Bayesian ordered beta regression (ordbetareg)",
             paste("Response: Percent", cell_label, "Content (0-1 scale)"),
             "Model: (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 + Diagnosis:Region:Age_c + Diagnosis:Region:log_Islet.Cells_c",
             "All estimates are posterior medians with 95% Highest Posterior Density (HPD) intervals.",
             "prob_greater_0 and prob_less_0 are posterior probabilities that the contrast is positive/negative.",
             "")
  ))
  writeData(wb, "Model Summary", model_summary, startRow = 9)
  
  # Sheet 2: Diagnosis
  addWorksheet(wb, "Diagnosis")
  writeData(wb, "Diagnosis", data.frame(
    Note = paste("Main effect of Diagnosis on predicted percent", tolower(cell_label), "content, averaged over Region and Sex.")))
  writeData(wb, "Diagnosis", diagnosis_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diagnosis", data.frame(x = ""), startRow = 6)
  writeData(wb, "Diagnosis", data.frame(Section = "Contrast:"), startRow = 7)
  writeData(wb, "Diagnosis", diagnosis_contrast %>% 
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 8)
  
  # Sheet 3: Region
  addWorksheet(wb, "Region")
  writeData(wb, "Region", data.frame(
    Note = paste("Main effect of pancreatic Region on predicted percent", tolower(cell_label), "content, averaged over Diagnosis and Sex.")))
  writeData(wb, "Region", region_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region", data.frame(x = ""), startRow = 7)
  writeData(wb, "Region", data.frame(Section = "Pairwise contrasts:"), startRow = 8)
  writeData(wb, "Region", region_contrasts %>% 
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 9)
  
  # Sheet 4: Sex
  addWorksheet(wb, "Sex")
  writeData(wb, "Sex", data.frame(
    Note = paste("Main effect of Sex on predicted percent", tolower(cell_label), "content, averaged over Diagnosis and Region.")))
  writeData(wb, "Sex", sex_marginal %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Sex", data.frame(x = ""), startRow = 6)
  writeData(wb, "Sex", data.frame(Section = "Contrast:"), startRow = 7)
  writeData(wb, "Sex", sex_contrast %>% 
              select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 8)
  
  # Sheet 5: Diagnosis x Region
  addWorksheet(wb, "Diagnosis x Region")
  writeData(wb, "Diagnosis x Region", data.frame(
    Note = "Diagnosis x Region interaction. Marginal means and contrasts averaged over Sex."))
  writeData(wb, "Diagnosis x Region", diag_region_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diagnosis x Region", data.frame(x = ""), startRow = 10)
  writeData(wb, "Diagnosis x Region", data.frame(Section = "Diagnosis effect (ND - T1D) within each Region:"), startRow = 11)
  writeData(wb, "Diagnosis x Region", diag_within_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 12)
  writeData(wb, "Diagnosis x Region", data.frame(x = ""), startRow = 16)
  writeData(wb, "Diagnosis x Region", data.frame(Section = "Region contrasts within each Diagnosis:"), startRow = 17)
  writeData(wb, "Diagnosis x Region", region_within_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 18)
  
  # Sheet 6: Diagnosis x Sex
  addWorksheet(wb, "Diagnosis x Sex")
  writeData(wb, "Diagnosis x Sex", data.frame(
    Note = "Diagnosis x Sex interaction. Marginal means and contrasts averaged over Region."))
  writeData(wb, "Diagnosis x Sex", diag_sex_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diagnosis x Sex", data.frame(x = ""), startRow = 8)
  writeData(wb, "Diagnosis x Sex", data.frame(Section = "Diagnosis effect (ND - T1D) within each Sex:"), startRow = 9)
  writeData(wb, "Diagnosis x Sex", diag_within_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 10)
  writeData(wb, "Diagnosis x Sex", data.frame(x = ""), startRow = 13)
  writeData(wb, "Diagnosis x Sex", data.frame(Section = "Sex effect (Female - Male) within each Diagnosis:"), startRow = 14)
  writeData(wb, "Diagnosis x Sex", sex_within_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 15)
  
  # Sheet 7: Region x Sex
  addWorksheet(wb, "Region x Sex")
  writeData(wb, "Region x Sex", data.frame(
    Note = "Region x Sex interaction. Marginal means and contrasts averaged over Diagnosis."))
  writeData(wb, "Region x Sex", region_sex_marginal_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Region x Sex", data.frame(x = ""), startRow = 10)
  writeData(wb, "Region x Sex", data.frame(Section = "Sex effect (Female - Male) within each Region:"), startRow = 11)
  writeData(wb, "Region x Sex", sex_within_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 12)
  writeData(wb, "Region x Sex", data.frame(x = ""), startRow = 16)
  writeData(wb, "Region x Sex", data.frame(Section = "Region contrasts within each Sex:"), startRow = 17)
  writeData(wb, "Region x Sex", region_within_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 18)
  
  # Sheet 8: Age Slopes
  addWorksheet(wb, "Age Slopes")
  writeData(wb, "Age Slopes", data.frame(
    Note = paste("Age slopes: change in predicted percent", tolower(cell_label), "per 1 year increase in Age.")))
  writeData(wb, "Age Slopes", data.frame(Section = "Overall Age slope (averaged over all factors):"), startRow = 3)
  writeData(wb, "Age Slopes", overall_age_summ %>% 
              select(predictor, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 4)
  writeData(wb, "Age Slopes", data.frame(x = ""), startRow = 6)
  writeData(wb, "Age Slopes", data.frame(Section = "Age slopes by Diagnosis:"), startRow = 7)
  writeData(wb, "Age Slopes", age_by_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 8)
  writeData(wb, "Age Slopes", data.frame(x = ""), startRow = 11)
  writeData(wb, "Age Slopes", data.frame(Section = "Age slope difference (T1D - ND):"), startRow = 12)
  writeData(wb, "Age Slopes", age_diag_diff_summ %>% 
              select(contrast, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 13)
  writeData(wb, "Age Slopes", data.frame(x = ""), startRow = 15)
  writeData(wb, "Age Slopes", data.frame(Section = "Age slopes by Region:"), startRow = 16)
  writeData(wb, "Age Slopes", age_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 17)
  writeData(wb, "Age Slopes", data.frame(x = ""), startRow = 21)
  writeData(wb, "Age Slopes", data.frame(Section = "Region contrasts (Age slopes):"), startRow = 22)
  writeData(wb, "Age Slopes", age_region_contrasts_df %>% mutate(across(where(is.numeric), fmt)), startRow = 23)
  writeData(wb, "Age Slopes", data.frame(x = ""), startRow = 27)
  writeData(wb, "Age Slopes", data.frame(Section = "Age slopes by Sex:"), startRow = 28)
  writeData(wb, "Age Slopes", age_by_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 29)
  writeData(wb, "Age Slopes", data.frame(x = ""), startRow = 32)
  writeData(wb, "Age Slopes", data.frame(Section = "Sex contrast (Age slopes):"), startRow = 33)
  writeData(wb, "Age Slopes", age_sex_contrast_df %>% mutate(across(where(is.numeric), fmt)), startRow = 34)
  
  # Sheet 9: Islet Size Slopes
  addWorksheet(wb, "Islet Size Slopes")
  writeData(wb, "Islet Size Slopes", data.frame(
    Note = paste("Islet size slopes: change in predicted percent", tolower(cell_label), "per 1 unit increase in log(Islet Cells). Note: 1 unit on log scale ≈ 2.72-fold increase in cell count.")))
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Overall Islet size slope (averaged over all factors):"), startRow = 3)
  writeData(wb, "Islet Size Slopes", overall_islet_summ %>% 
              select(predictor, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 4)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 6)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Islet size slopes by Diagnosis:"), startRow = 7)
  writeData(wb, "Islet Size Slopes", islet_by_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 8)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 11)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Islet slope difference (T1D - ND):"), startRow = 12)
  writeData(wb, "Islet Size Slopes", islet_diag_diff_summ %>% 
              select(contrast, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
              mutate(across(where(is.numeric), fmt)), startRow = 13)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 15)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Islet size slopes by Region:"), startRow = 16)
  writeData(wb, "Islet Size Slopes", islet_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 17)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 21)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Region contrasts (Islet slopes):"), startRow = 22)
  writeData(wb, "Islet Size Slopes", islet_region_contrasts_df %>% mutate(across(where(is.numeric), fmt)), startRow = 23)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 27)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Islet size slopes by Sex:"), startRow = 28)
  writeData(wb, "Islet Size Slopes", islet_by_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 29)
  writeData(wb, "Islet Size Slopes", data.frame(x = ""), startRow = 32)
  writeData(wb, "Islet Size Slopes", data.frame(Section = "Sex contrast (Islet slopes):"), startRow = 33)
  writeData(wb, "Islet Size Slopes", islet_sex_contrast_df %>% mutate(across(where(is.numeric), fmt)), startRow = 34)
  
  # Sheet 10: Diag x Region x Age
  addWorksheet(wb, "Diag x Region x Age")
  writeData(wb, "Diag x Region x Age", data.frame(
    Note = "KEY 3-WAY INTERACTION: Diagnosis x Region x Age. Age slopes for each Diagnosis x Region cell (averaged over Sex)."))
  writeData(wb, "Diag x Region x Age", age_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diag x Region x Age", data.frame(x = ""), startRow = 10)
  writeData(wb, "Diag x Region x Age", data.frame(Section = "Diagnosis effect (T1D - ND) on Age slopes within each Region:"), startRow = 11)
  writeData(wb, "Diag x Region x Age", age_slope_diff_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 12)
  
  # Sheet 11: Diag x Region x Islet
  addWorksheet(wb, "Diag x Region x Islet")
  writeData(wb, "Diag x Region x Islet", data.frame(
    Note = "KEY 3-WAY INTERACTION: Diagnosis x Region x Islet Size. Islet slopes for each Diagnosis x Region cell (averaged over Sex)."))
  writeData(wb, "Diag x Region x Islet", islet_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "Diag x Region x Islet", data.frame(x = ""), startRow = 10)
  writeData(wb, "Diag x Region x Islet", data.frame(Section = "Diagnosis effect (T1D - ND) on Islet slopes within each Region:"), startRow = 11)
  writeData(wb, "Diag x Region x Islet", islet_slope_diff_by_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 12)
  writeData(wb, "Diag x Region x Islet", data.frame(x = ""), startRow = 16)
  writeData(wb, "Diag x Region x Islet", data.frame(Section = "Region contrasts for Islet slopes within each Diagnosis:"), startRow = 17)
  writeData(wb, "Diag x Region x Islet", islet_region_within_diag_df %>% mutate(across(where(is.numeric), fmt)), startRow = 18)
  
  # Sheet 12: Supp Pairwise Diag x Sex
  addWorksheet(wb, "Supp Pairwise Diag x Sex")
  writeData(wb, "Supp Pairwise Diag x Sex", data.frame(
    Note = "Supplementary: All pairwise comparisons of Diagnosis x Sex cells (averaged over Region)."))
  writeData(wb, "Supp Pairwise Diag x Sex", pairwise_diag_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  
  # Sheet 13: Supp Pairwise Diag x Region
  addWorksheet(wb, "Supp Pairwise Diag x Reg")
  writeData(wb, "Supp Pairwise Diag x Reg", data.frame(
    Note = "Supplementary: All pairwise comparisons of Diagnosis x Region cells (averaged over Sex)."))
  writeData(wb, "Supp Pairwise Diag x Reg", pairwise_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  
  # Save workbook
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0(cell_label, "_Content_Primary_Analysis.xlsx"))
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Primary analysis saved to:", output_file, "\n")
  
  return(invisible(wb))
}

# =============================================================================
# SUPPLEMENTARY ANALYSIS FUNCTION (for *_mid models)
# =============================================================================

generate_supplementary_excel <- function(model, data, cell_type, cell_label, output_dir = "Proportion/Results") {
  
  cat("\n")
  cat("=============================================================================\n")
  cat("GENERATING SUPPLEMENTARY ANALYSIS EXCEL:", toupper(cell_label), "\n")
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
  # CATEGORICAL GRID SETUP
  # -------------------------------------------------------------------------
  cat("Setting up prediction grids...\n")
  
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
  
  # Age and islet grids
  age_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Age_c = c(-1, 1),
    log_Islet.Cells_c = 0
  )
  
  islet_grid <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    Age_c = 0,
    log_Islet.Cells_c = c(-1, 1)
  )
  
  post_age <- get_posterior_fitted(model, age_grid)
  post_islet <- get_posterior_fitted(model, islet_grid)
  
  # -------------------------------------------------------------------------
  # SHEET 3: 3-WAY CELL MEANS (Diagnosis x Region x Sex)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 3: 3-Way Cell Means...\n")
  
  cell_means_3way <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == reg & cat_grid$Sex == s)
        draws <- post_fitted[, idx]
        summ <- summarize_posterior(draws)
        summ$Diagnosis <- diag
        summ$Region <- reg
        summ$Sex <- s
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
      nd_idx <- which(cat_grid$Diagnosis == "ND" & cat_grid$Region == reg & cat_grid$Sex == s)
      t1d_idx <- which(cat_grid$Diagnosis == "T1D" & cat_grid$Region == reg & cat_grid$Sex == s)
      nd_draws <- post_fitted[, nd_idx]
      t1d_draws <- post_fitted[, t1d_idx]
      contrast <- compute_contrast(nd_draws, t1d_draws, "ND - T1D")
      contrast$Region <- reg
      contrast$Sex <- s
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
      f_idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == reg & cat_grid$Sex == "Female")
      m_idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == reg & cat_grid$Sex == "Male")
      f_draws <- post_fitted[, f_idx]
      m_draws <- post_fitted[, m_idx]
      contrast <- compute_contrast(f_draws, m_draws, "Female - Male")
      contrast$Diagnosis <- diag
      contrast$Region <- reg
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
      head_idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == "Head" & cat_grid$Sex == s)
      body_idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == "Body" & cat_grid$Sex == s)
      tail_idx <- which(cat_grid$Diagnosis == diag & cat_grid$Region == "Tail" & cat_grid$Sex == s)
      
      head_draws <- post_fitted[, head_idx]
      body_draws <- post_fitted[, body_idx]
      tail_draws <- post_fitted[, tail_idx]
      
      hb <- compute_contrast(head_draws, body_draws, "Head - Body")
      hb$Diagnosis <- diag
      hb$Sex <- s
      ht <- compute_contrast(head_draws, tail_draws, "Head - Tail")
      ht$Diagnosis <- diag
      ht$Sex <- s
      bt <- compute_contrast(body_draws, tail_draws, "Body - Tail")
      bt$Diagnosis <- diag
      bt$Sex <- s
      
      region_within_diag_sex[[paste(diag, s, "HB", sep = "_")]] <- hb
      region_within_diag_sex[[paste(diag, s, "HT", sep = "_")]] <- ht
      region_within_diag_sex[[paste(diag, s, "BT", sep = "_")]] <- bt
    }
  }
  region_within_diag_sex_df <- bind_rows(region_within_diag_sex) %>%
    select(Diagnosis, Sex, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # COMPUTE ALL CELL-LEVEL SLOPES
  # -------------------------------------------------------------------------
  cat("Computing 4-way slopes...\n")
  
  age_slopes_list <- list()
  islet_slopes_list <- list()
  
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        # Age slopes
        low_idx <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & 
                           age_grid$Sex == s & age_grid$Age_c == -1)
        high_idx <- which(age_grid$Diagnosis == diag & age_grid$Region == reg & 
                            age_grid$Sex == s & age_grid$Age_c == 1)
        age_slopes_list[[paste(diag, reg, s, sep = "_")]] <- (post_age[, high_idx] - post_age[, low_idx]) / 2
        
        # Islet slopes
        low_idx <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & 
                           islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == -1)
        high_idx <- which(islet_grid$Diagnosis == diag & islet_grid$Region == reg & 
                            islet_grid$Sex == s & islet_grid$log_Islet.Cells_c == 1)
        islet_slopes_list[[paste(diag, reg, s, sep = "_")]] <- (post_islet[, high_idx] - post_islet[, low_idx]) / 2
      }
    }
  }
  
  # =============================================================================
  # ADDITIONAL 3-WAY MARGINALIZED SLOPE SHEETS
  # =============================================================================
  #
  # Insert this code into generate_supplementary_excel() AFTER the 4-way slope
  # computation block (after islet_slopes_list is fully populated, ~line 1136)
  # and BEFORE the "SHEET 7: 4-WAY AGE SLOPES" section (~line 1139).
  #
  # This adds four new sheets to each Supplementary Excel workbook:
  #   - "3-Way Diag x Sex Age Slopes"    (Panel B of 3-way figure)
  #   - "3-Way Diag x Sex Islet Slopes"  (Panel C of 3-way figure)
  #   - "3-Way Reg x Sex Age Slopes"     (Panel D of 3-way figure)
  #   - "3-Way Reg x Sex Islet Slopes"   (Panel E of 3-way figure)
  #
  # Each sheet provides:
  #   1. Marginalized slopes for each cell
  #   2. Relevant pairwise contrasts
  # =============================================================================
  
  # -------------------------------------------------------------------------
  # 3-WAY MARGINALIZED SLOPES
  # -------------------------------------------------------------------------
  cat("Computing 3-way marginalized slopes...\n")
  
  # =========================================================================
  # 3-Way: Diagnosis x Sex (marginalized over Region) - AGE SLOPES
  # For Panel B: Diagnosis x Sex x Age
  # =========================================================================
  
  age_diag_sex_slopes <- list()
  for (diag in c("ND", "T1D")) {
    for (s in c("Female", "Male")) {
      # Average age slope draws across the 3 regions
      region_slopes <- list()
      for (reg in c("Head", "Body", "Tail")) {
        region_slopes[[reg]] <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
      slope_draws <- rowMeans(do.call(cbind, region_slopes))
      summ <- summarize_posterior(slope_draws)
      summ$Diagnosis <- diag
      summ$Sex <- s
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      age_diag_sex_slopes[[paste(diag, s, sep = "_")]] <- list(summary = summ, draws = slope_draws)
    }
  }
  
  age_diag_sex_slopes_df <- bind_rows(lapply(age_diag_sex_slopes, `[[`, "summary")) %>%
    select(Diagnosis, Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
    rename(Age_slope = estimate, Age_lower.HPD = lower.HPD, Age_upper.HPD = upper.HPD)
  
  # Contrasts: Diagnosis effect (T1D - ND) within each Sex
  age_diag_within_sex_3way <- list()
  for (s in c("Female", "Male")) {
    diff <- age_diag_sex_slopes[[paste("T1D", s, sep = "_")]]$draws - 
      age_diag_sex_slopes[[paste("ND", s, sep = "_")]]$draws
    summ <- summarize_posterior(diff)
    summ$Sex <- s
    summ$contrast <- "T1D - ND"
    summ$prob_greater_0 <- mean(diff > 0)
    summ$prob_less_0 <- mean(diff < 0)
    age_diag_within_sex_3way[[s]] <- summ
  }
  age_diag_within_sex_3way_df <- bind_rows(age_diag_within_sex_3way) %>%
    select(Sex, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Contrasts: Sex effect (Female - Male) within each Diagnosis
  age_sex_within_diag_3way <- list()
  for (diag in c("ND", "T1D")) {
    diff <- age_diag_sex_slopes[[paste(diag, "Female", sep = "_")]]$draws - 
      age_diag_sex_slopes[[paste(diag, "Male", sep = "_")]]$draws
    summ <- summarize_posterior(diff)
    summ$Diagnosis <- diag
    summ$contrast <- "Female - Male"
    summ$prob_greater_0 <- mean(diff > 0)
    summ$prob_less_0 <- mean(diff < 0)
    age_sex_within_diag_3way[[diag]] <- summ
  }
  age_sex_within_diag_3way_df <- bind_rows(age_sex_within_diag_3way) %>%
    select(Diagnosis, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # =========================================================================
  # 3-Way: Diagnosis x Sex (marginalized over Region) - ISLET SLOPES
  # For Panel C: Diagnosis x Sex x Islet Size
  # =========================================================================
  
  islet_diag_sex_slopes <- list()
  for (diag in c("ND", "T1D")) {
    for (s in c("Female", "Male")) {
      region_slopes <- list()
      for (reg in c("Head", "Body", "Tail")) {
        region_slopes[[reg]] <- islet_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
      slope_draws <- rowMeans(do.call(cbind, region_slopes))
      summ <- summarize_posterior(slope_draws)
      summ$Diagnosis <- diag
      summ$Sex <- s
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      islet_diag_sex_slopes[[paste(diag, s, sep = "_")]] <- list(summary = summ, draws = slope_draws)
    }
  }
  
  islet_diag_sex_slopes_df <- bind_rows(lapply(islet_diag_sex_slopes, `[[`, "summary")) %>%
    select(Diagnosis, Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
    rename(Islet_slope = estimate, Islet_lower.HPD = lower.HPD, Islet_upper.HPD = upper.HPD)
  
  # Contrasts: Diagnosis effect (T1D - ND) within each Sex
  islet_diag_within_sex_3way <- list()
  for (s in c("Female", "Male")) {
    diff <- islet_diag_sex_slopes[[paste("T1D", s, sep = "_")]]$draws - 
      islet_diag_sex_slopes[[paste("ND", s, sep = "_")]]$draws
    summ <- summarize_posterior(diff)
    summ$Sex <- s
    summ$contrast <- "T1D - ND"
    summ$prob_greater_0 <- mean(diff > 0)
    summ$prob_less_0 <- mean(diff < 0)
    islet_diag_within_sex_3way[[s]] <- summ
  }
  islet_diag_within_sex_3way_df <- bind_rows(islet_diag_within_sex_3way) %>%
    select(Sex, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Contrasts: Sex effect (Female - Male) within each Diagnosis
  islet_sex_within_diag_3way <- list()
  for (diag in c("ND", "T1D")) {
    diff <- islet_diag_sex_slopes[[paste(diag, "Female", sep = "_")]]$draws - 
      islet_diag_sex_slopes[[paste(diag, "Male", sep = "_")]]$draws
    summ <- summarize_posterior(diff)
    summ$Diagnosis <- diag
    summ$contrast <- "Female - Male"
    summ$prob_greater_0 <- mean(diff > 0)
    summ$prob_less_0 <- mean(diff < 0)
    islet_sex_within_diag_3way[[diag]] <- summ
  }
  islet_sex_within_diag_3way_df <- bind_rows(islet_sex_within_diag_3way) %>%
    select(Diagnosis, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # =========================================================================
  # 3-Way: Region x Sex (marginalized over Diagnosis) - AGE SLOPES
  # For Panel D: Region x Sex x Age
  # =========================================================================
  
  age_reg_sex_slopes <- list()
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      diag_slopes <- list()
      for (diag in c("ND", "T1D")) {
        diag_slopes[[diag]] <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
      slope_draws <- rowMeans(do.call(cbind, diag_slopes))
      summ <- summarize_posterior(slope_draws)
      summ$Region <- reg
      summ$Sex <- s
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      age_reg_sex_slopes[[paste(reg, s, sep = "_")]] <- list(summary = summ, draws = slope_draws)
    }
  }
  
  age_reg_sex_slopes_df <- bind_rows(lapply(age_reg_sex_slopes, `[[`, "summary")) %>%
    select(Region, Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
    rename(Age_slope = estimate, Age_lower.HPD = lower.HPD, Age_upper.HPD = upper.HPD)
  
  # Contrasts: Region pairwise within each Sex
  age_region_within_sex_3way <- list()
  for (s in c("Female", "Male")) {
    for (pair in list(c("Head", "Body"), c("Head", "Tail"), c("Body", "Tail"))) {
      diff <- age_reg_sex_slopes[[paste(pair[1], s, sep = "_")]]$draws - 
        age_reg_sex_slopes[[paste(pair[2], s, sep = "_")]]$draws
      label <- paste(pair[1], "-", pair[2])
      summ <- summarize_posterior(diff)
      summ$Sex <- s
      summ$contrast <- label
      summ$prob_greater_0 <- mean(diff > 0)
      summ$prob_less_0 <- mean(diff < 0)
      age_region_within_sex_3way[[paste(s, label, sep = "_")]] <- summ
    }
  }
  age_region_within_sex_3way_df <- bind_rows(age_region_within_sex_3way) %>%
    select(Sex, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Contrasts: Sex effect (Female - Male) within each Region
  age_sex_within_reg_3way <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff <- age_reg_sex_slopes[[paste(reg, "Female", sep = "_")]]$draws - 
      age_reg_sex_slopes[[paste(reg, "Male", sep = "_")]]$draws
    summ <- summarize_posterior(diff)
    summ$Region <- reg
    summ$contrast <- "Female - Male"
    summ$prob_greater_0 <- mean(diff > 0)
    summ$prob_less_0 <- mean(diff < 0)
    age_sex_within_reg_3way[[reg]] <- summ
  }
  age_sex_within_reg_3way_df <- bind_rows(age_sex_within_reg_3way) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # =========================================================================
  # 3-Way: Region x Sex (marginalized over Diagnosis) - ISLET SLOPES
  # For Panel E: Region x Sex x Islet Size
  # =========================================================================
  
  islet_reg_sex_slopes <- list()
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      diag_slopes <- list()
      for (diag in c("ND", "T1D")) {
        diag_slopes[[diag]] <- islet_slopes_list[[paste(diag, reg, s, sep = "_")]]
      }
      slope_draws <- rowMeans(do.call(cbind, diag_slopes))
      summ <- summarize_posterior(slope_draws)
      summ$Region <- reg
      summ$Sex <- s
      summ$prob_positive <- mean(slope_draws > 0)
      summ$prob_negative <- mean(slope_draws < 0)
      islet_reg_sex_slopes[[paste(reg, s, sep = "_")]] <- list(summary = summ, draws = slope_draws)
    }
  }
  
  islet_reg_sex_slopes_df <- bind_rows(lapply(islet_reg_sex_slopes, `[[`, "summary")) %>%
    select(Region, Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative) %>%
    rename(Islet_slope = estimate, Islet_lower.HPD = lower.HPD, Islet_upper.HPD = upper.HPD)
  
  # Contrasts: Region pairwise within each Sex
  islet_region_within_sex_3way <- list()
  for (s in c("Female", "Male")) {
    for (pair in list(c("Head", "Body"), c("Head", "Tail"), c("Body", "Tail"))) {
      diff <- islet_reg_sex_slopes[[paste(pair[1], s, sep = "_")]]$draws - 
        islet_reg_sex_slopes[[paste(pair[2], s, sep = "_")]]$draws
      label <- paste(pair[1], "-", pair[2])
      summ <- summarize_posterior(diff)
      summ$Sex <- s
      summ$contrast <- label
      summ$prob_greater_0 <- mean(diff > 0)
      summ$prob_less_0 <- mean(diff < 0)
      islet_region_within_sex_3way[[paste(s, label, sep = "_")]] <- summ
    }
  }
  islet_region_within_sex_3way_df <- bind_rows(islet_region_within_sex_3way) %>%
    select(Sex, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # Contrasts: Sex effect (Female - Male) within each Region
  islet_sex_within_reg_3way <- list()
  for (reg in c("Head", "Body", "Tail")) {
    diff <- islet_reg_sex_slopes[[paste(reg, "Female", sep = "_")]]$draws - 
      islet_reg_sex_slopes[[paste(reg, "Male", sep = "_")]]$draws
    summ <- summarize_posterior(diff)
    summ$Region <- reg
    summ$contrast <- "Female - Male"
    summ$prob_greater_0 <- mean(diff > 0)
    summ$prob_less_0 <- mean(diff < 0)
    islet_sex_within_reg_3way[[reg]] <- summ
  }
  islet_sex_within_reg_3way_df <- bind_rows(islet_sex_within_reg_3way) %>%
    select(Region, contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # SHEET 7: 4-WAY AGE SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 7: 4-Way Age Slopes...\n")
  
  age_4way <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        slope_draws <- age_slopes_list[[paste(diag, reg, s, sep = "_")]]
        summ <- summarize_posterior(slope_draws)
        summ$Diagnosis <- diag
        summ$Region <- reg
        summ$Sex <- s
        summ$prob_positive <- mean(slope_draws > 0)
        summ$prob_negative <- mean(slope_draws < 0)
        age_4way[[paste(diag, reg, s, sep = "_")]] <- summ
      }
    }
  }
  age_4way_df <- bind_rows(age_4way) %>%
    select(Diagnosis, Region, Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # Age slope contrasts
  # Diagnosis effect on age slopes within Region x Sex
  age_slope_diag_within_reg_sex <- list()
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      nd_draws <- age_slopes_list[[paste("ND", reg, s, sep = "_")]]
      t1d_draws <- age_slopes_list[[paste("T1D", reg, s, sep = "_")]]
      diff <- t1d_draws - nd_draws
      summ <- summarize_posterior(diff)
      summ$Region <- reg
      summ$Sex <- s
      summ$contrast <- "T1D - ND"
      summ$prob_positive <- mean(diff > 0)
      age_slope_diag_within_reg_sex[[paste(reg, s, sep = "_")]] <- summ
    }
  }
  age_slope_diag_within_reg_sex_df <- bind_rows(age_slope_diag_within_reg_sex) %>%
    select(Region, Sex, contrast, estimate, lower.HPD, upper.HPD, prob_positive)
  
  # Sex effect on age slopes within Diagnosis x Region
  age_slope_sex_within_diag_reg <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      f_draws <- age_slopes_list[[paste(diag, reg, "Female", sep = "_")]]
      m_draws <- age_slopes_list[[paste(diag, reg, "Male", sep = "_")]]
      diff <- f_draws - m_draws
      summ <- summarize_posterior(diff)
      summ$Diagnosis <- diag
      summ$Region <- reg
      summ$contrast <- "Female - Male"
      summ$prob_positive <- mean(diff > 0)
      age_slope_sex_within_diag_reg[[paste(diag, reg, sep = "_")]] <- summ
    }
  }
  age_slope_sex_within_diag_reg_df <- bind_rows(age_slope_sex_within_diag_reg) %>%
    select(Diagnosis, Region, contrast, estimate, lower.HPD, upper.HPD, prob_positive)
  
  # -------------------------------------------------------------------------
  # SHEET 8: 4-WAY ISLET SLOPES
  # -------------------------------------------------------------------------
  cat("Creating Sheet 8: 4-Way Islet Slopes...\n")
  
  islet_4way <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      for (s in c("Female", "Male")) {
        slope_draws <- islet_slopes_list[[paste(diag, reg, s, sep = "_")]]
        summ <- summarize_posterior(slope_draws)
        summ$Diagnosis <- diag
        summ$Region <- reg
        summ$Sex <- s
        summ$prob_positive <- mean(slope_draws > 0)
        summ$prob_negative <- mean(slope_draws < 0)
        islet_4way[[paste(diag, reg, s, sep = "_")]] <- summ
      }
    }
  }
  islet_4way_df <- bind_rows(islet_4way) %>%
    select(Diagnosis, Region, Sex, estimate, lower.HPD, upper.HPD, prob_positive, prob_negative)
  
  # Islet slope contrasts
  islet_slope_diag_within_reg_sex <- list()
  for (reg in c("Head", "Body", "Tail")) {
    for (s in c("Female", "Male")) {
      nd_draws <- islet_slopes_list[[paste("ND", reg, s, sep = "_")]]
      t1d_draws <- islet_slopes_list[[paste("T1D", reg, s, sep = "_")]]
      diff <- t1d_draws - nd_draws
      summ <- summarize_posterior(diff)
      summ$Region <- reg
      summ$Sex <- s
      summ$contrast <- "T1D - ND"
      summ$prob_positive <- mean(diff > 0)
      islet_slope_diag_within_reg_sex[[paste(reg, s, sep = "_")]] <- summ
    }
  }
  islet_slope_diag_within_reg_sex_df <- bind_rows(islet_slope_diag_within_reg_sex) %>%
    select(Region, Sex, contrast, estimate, lower.HPD, upper.HPD, prob_positive)
  
  islet_slope_sex_within_diag_reg <- list()
  for (diag in c("ND", "T1D")) {
    for (reg in c("Head", "Body", "Tail")) {
      f_draws <- islet_slopes_list[[paste(diag, reg, "Female", sep = "_")]]
      m_draws <- islet_slopes_list[[paste(diag, reg, "Male", sep = "_")]]
      diff <- f_draws - m_draws
      summ <- summarize_posterior(diff)
      summ$Diagnosis <- diag
      summ$Region <- reg
      summ$contrast <- "Female - Male"
      summ$prob_positive <- mean(diff > 0)
      islet_slope_sex_within_diag_reg[[paste(diag, reg, sep = "_")]] <- summ
    }
  }
  islet_slope_sex_within_diag_reg_df <- bind_rows(islet_slope_sex_within_diag_reg) %>%
    select(Diagnosis, Region, contrast, estimate, lower.HPD, upper.HPD, prob_positive)
  
  # -------------------------------------------------------------------------
  # SHEET 9: FULL PAIRWISE (ALL 12 CELLS)
  # -------------------------------------------------------------------------
  cat("Creating Sheet 9: Full Pairwise All Cells...\n")
  
  all_cells <- expand.grid(
    Diagnosis = c("ND", "T1D"),
    Region = c("Head", "Body", "Tail"),
    Sex = c("Female", "Male"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(cell_name = paste(Diagnosis, Region, Sex))
  
  full_pairwise <- list()
  for (i in 1:(nrow(all_cells) - 1)) {
    for (j in (i + 1):nrow(all_cells)) {
      cell1 <- all_cells[i, ]
      cell2 <- all_cells[j, ]
      
      idx1 <- which(cat_grid$Diagnosis == cell1$Diagnosis & 
                      cat_grid$Region == cell1$Region & 
                      cat_grid$Sex == cell1$Sex)
      idx2 <- which(cat_grid$Diagnosis == cell2$Diagnosis & 
                      cat_grid$Region == cell2$Region & 
                      cat_grid$Sex == cell2$Sex)
      
      draws1 <- post_fitted[, idx1]
      draws2 <- post_fitted[, idx2]
      
      label <- paste0(cell1$cell_name, " - ", cell2$cell_name)
      contrast <- compute_contrast(draws1, draws2, label)
      full_pairwise[[label]] <- contrast
    }
  }
  full_pairwise_df <- bind_rows(full_pairwise) %>%
    select(contrast, estimate, lower.HPD, upper.HPD, prob_greater_0, prob_less_0)
  
  # -------------------------------------------------------------------------
  # CREATE EXCEL WORKBOOK
  # -------------------------------------------------------------------------
  cat("\nCreating Supplementary Analysis Excel workbook...\n")
  
  wb <- createWorkbook()
  
  # Sheet 1: Model Summary
  addWorksheet(wb, "Model Summary")
  writeData(wb, "Model Summary", data.frame(
    Note = c(paste("SUPPLEMENTARY EXPLORATORY ANALYSIS:", cell_label, "Content"),
             "",
             "CAUTION: This analysis is EXPLORATORY due to limited sample sizes in some cells.",
             "Results should be interpreted with caution, particularly for cells with <20 donors.",
             "",
             "Bayesian ordered beta regression (ordbetareg)",
             paste("Response: Percent", cell_label, "Content (0-1 scale)"),
             "Model includes 3-way (Diagnosis:Region:Sex) and 4-way (Diagnosis:Region:Sex:Age/Islet) interactions",
             "",
             "All estimates are posterior medians with 95% HPD intervals.",
             "prob_greater_0 and prob_less_0 are posterior probabilities.",
             "")
  ))
  writeData(wb, "Model Summary", model_summary, startRow = 14)
  
  # Sheet 2: Donor Counts
  addWorksheet(wb, "Donor Counts")
  writeData(wb, "Donor Counts", data.frame(
    Note = c("IMPORTANT: Donor counts by cell. Cells with <20 donors (highlighted) have limited power.",
             "Interpret results for these cells with caution.",
             "")
  ))
  writeData(wb, "Donor Counts", donor_counts, startRow = 5)
  
  # Sheet 3: 3-Way Cell Means
  addWorksheet(wb, "Diag x Reg x Sex Means")
  writeData(wb, "Diag x Reg x Sex Means", data.frame(
    Note = paste("EXPLORATORY: 3-way cell means for Diagnosis x Region x Sex from", cell_type, "mid model.")))
  writeData(wb, "Diag x Reg x Sex Means", cell_means_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  
  # Sheet 4: Diagnosis within Region x Sex
  addWorksheet(wb, "Diag within Reg x Sex")
  writeData(wb, "Diag within Reg x Sex", data.frame(
    Note = "EXPLORATORY: Diagnosis effect (ND - T1D) within each Region x Sex cell."))
  writeData(wb, "Diag within Reg x Sex", diag_within_region_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  
  # Sheet 5: Sex within Diagnosis x Region
  addWorksheet(wb, "Sex within Diag x Reg")
  writeData(wb, "Sex within Diag x Reg", data.frame(
    Note = "EXPLORATORY: Sex effect (Female - Male) within each Diagnosis x Region cell."))
  writeData(wb, "Sex within Diag x Reg", sex_within_diag_region_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  
  # Sheet 6: Region within Diagnosis x Sex
  addWorksheet(wb, "Region within Diag x Sex")
  writeData(wb, "Region within Diag x Sex", data.frame(
    Note = "EXPLORATORY: Region contrasts within each Diagnosis x Sex cell."))
  writeData(wb, "Region within Diag x Sex", region_within_diag_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  
  # Sheet 10: 3-Way Diag x Sex Age Slopes (Panel B)
  addWorksheet(wb, "3-Way Diag x Sex Age Slopes")
  writeData(wb, "3-Way Diag x Sex Age Slopes", data.frame(
    Note = "EXPLORATORY: Diagnosis x Sex age slopes (marginalized over Region). For 3-way figure Panel B."))
  writeData(wb, "3-Way Diag x Sex Age Slopes", 
            age_diag_sex_slopes_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "3-Way Diag x Sex Age Slopes", data.frame(x = ""), startRow = 8)
  writeData(wb, "3-Way Diag x Sex Age Slopes", data.frame(
    Section = "Diagnosis effect (T1D - ND) on age slopes within each Sex:"), startRow = 9)
  writeData(wb, "3-Way Diag x Sex Age Slopes", 
            age_diag_within_sex_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 10)
  writeData(wb, "3-Way Diag x Sex Age Slopes", data.frame(x = ""), startRow = 13)
  writeData(wb, "3-Way Diag x Sex Age Slopes", data.frame(
    Section = "Sex effect (Female - Male) on age slopes within each Diagnosis:"), startRow = 14)
  writeData(wb, "3-Way Diag x Sex Age Slopes", 
            age_sex_within_diag_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 15)
  
  # Sheet 11: 3-Way Diag x Sex Islet Slopes (Panel C)
  addWorksheet(wb, "3-Way Diag x Sex Islet Slopes")
  writeData(wb, "3-Way Diag x Sex Islet Slopes", data.frame(
    Note = "EXPLORATORY: Diagnosis x Sex islet size slopes (marginalized over Region). For 3-way figure Panel C."))
  writeData(wb, "3-Way Diag x Sex Islet Slopes", 
            islet_diag_sex_slopes_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "3-Way Diag x Sex Islet Slopes", data.frame(x = ""), startRow = 8)
  writeData(wb, "3-Way Diag x Sex Islet Slopes", data.frame(
    Section = "Diagnosis effect (T1D - ND) on islet slopes within each Sex:"), startRow = 9)
  writeData(wb, "3-Way Diag x Sex Islet Slopes", 
            islet_diag_within_sex_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 10)
  writeData(wb, "3-Way Diag x Sex Islet Slopes", data.frame(x = ""), startRow = 13)
  writeData(wb, "3-Way Diag x Sex Islet Slopes", data.frame(
    Section = "Sex effect (Female - Male) on islet slopes within each Diagnosis:"), startRow = 14)
  writeData(wb, "3-Way Diag x Sex Islet Slopes", 
            islet_sex_within_diag_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 15)
  
  # Sheet 12: 3-Way Region x Sex Age Slopes (Panel D)
  addWorksheet(wb, "3-Way Reg x Sex Age Slopes")
  writeData(wb, "3-Way Reg x Sex Age Slopes", data.frame(
    Note = "EXPLORATORY: Region x Sex age slopes (marginalized over Diagnosis). For 3-way figure Panel D."))
  writeData(wb, "3-Way Reg x Sex Age Slopes", 
            age_reg_sex_slopes_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "3-Way Reg x Sex Age Slopes", data.frame(x = ""), startRow = 10)
  writeData(wb, "3-Way Reg x Sex Age Slopes", data.frame(
    Section = "Region pairwise contrasts on age slopes within each Sex:"), startRow = 11)
  writeData(wb, "3-Way Reg x Sex Age Slopes", 
            age_region_within_sex_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 12)
  writeData(wb, "3-Way Reg x Sex Age Slopes", data.frame(x = ""), startRow = 19)
  writeData(wb, "3-Way Reg x Sex Age Slopes", data.frame(
    Section = "Sex effect (Female - Male) on age slopes within each Region:"), startRow = 20)
  writeData(wb, "3-Way Reg x Sex Age Slopes", 
            age_sex_within_reg_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 21)
  
  # Sheet 13: 3-Way Region x Sex Islet Slopes (Panel E)
  addWorksheet(wb, "3-Way Reg x Sex Islet Slopes")
  writeData(wb, "3-Way Reg x Sex Islet Slopes", data.frame(
    Note = "EXPLORATORY: Region x Sex islet size slopes (marginalized over Diagnosis). For 3-way figure Panel E."))
  writeData(wb, "3-Way Reg x Sex Islet Slopes", 
            islet_reg_sex_slopes_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "3-Way Reg x Sex Islet Slopes", data.frame(x = ""), startRow = 10)
  writeData(wb, "3-Way Reg x Sex Islet Slopes", data.frame(
    Section = "Region pairwise contrasts on islet slopes within each Sex:"), startRow = 11)
  writeData(wb, "3-Way Reg x Sex Islet Slopes", 
            islet_region_within_sex_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 12)
  writeData(wb, "3-Way Reg x Sex Islet Slopes", data.frame(x = ""), startRow = 19)
  writeData(wb, "3-Way Reg x Sex Islet Slopes", data.frame(
    Section = "Sex effect (Female - Male) on islet slopes within each Region:"), startRow = 20)
  writeData(wb, "3-Way Reg x Sex Islet Slopes", 
            islet_sex_within_reg_3way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 21)
  
  # Sheet 7: 4-Way Age Slopes
  addWorksheet(wb, "4-Way Age Slopes")
  writeData(wb, "4-Way Age Slopes", data.frame(
    Note = "EXPLORATORY: Age slopes for each Diagnosis x Region x Sex cell."))
  writeData(wb, "4-Way Age Slopes", age_4way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "4-Way Age Slopes", data.frame(x = ""), startRow = 16)
  writeData(wb, "4-Way Age Slopes", data.frame(Section = "Diagnosis effect (T1D - ND) on Age slopes within each Region x Sex:"), startRow = 17)
  writeData(wb, "4-Way Age Slopes", age_slope_diag_within_reg_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 18)
  writeData(wb, "4-Way Age Slopes", data.frame(x = ""), startRow = 25)
  writeData(wb, "4-Way Age Slopes", data.frame(Section = "Sex effect (Female - Male) on Age slopes within each Diagnosis x Region:"), startRow = 26)
  writeData(wb, "4-Way Age Slopes", age_slope_sex_within_diag_reg_df %>% mutate(across(where(is.numeric), fmt)), startRow = 27)
  
  # Sheet 8: 4-Way Islet Slopes
  addWorksheet(wb, "4-Way Islet Slopes")
  writeData(wb, "4-Way Islet Slopes", data.frame(
    Note = "EXPLORATORY: Islet size slopes for each Diagnosis x Region x Sex cell."))
  writeData(wb, "4-Way Islet Slopes", islet_4way_df %>% mutate(across(where(is.numeric), fmt)), startRow = 3)
  writeData(wb, "4-Way Islet Slopes", data.frame(x = ""), startRow = 16)
  writeData(wb, "4-Way Islet Slopes", data.frame(Section = "Diagnosis effect (T1D - ND) on Islet slopes within each Region x Sex:"), startRow = 17)
  writeData(wb, "4-Way Islet Slopes", islet_slope_diag_within_reg_sex_df %>% mutate(across(where(is.numeric), fmt)), startRow = 18)
  writeData(wb, "4-Way Islet Slopes", data.frame(x = ""), startRow = 25)
  writeData(wb, "4-Way Islet Slopes", data.frame(Section = "Sex effect (Female - Male) on Islet slopes within each Diagnosis x Region:"), startRow = 26)
  writeData(wb, "4-Way Islet Slopes", islet_slope_sex_within_diag_reg_df %>% mutate(across(where(is.numeric), fmt)), startRow = 27)
  
  # Sheet 9: Full Pairwise All Cells
  addWorksheet(wb, "Full Pairwise All Cells")
  writeData(wb, "Full Pairwise All Cells", data.frame(
    Note = c("EXPLORATORY: All pairwise comparisons between the 12 Diagnosis x Region x Sex cells.",
             "Note: 66 total comparisons. Posterior probabilities do not require multiple comparison corrections.",
             "")))
  writeData(wb, "Full Pairwise All Cells", full_pairwise_df %>% mutate(across(where(is.numeric), fmt)), startRow = 5)
  
  # Save workbook
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0(cell_label, "_Content_Supplementary_Exploratory.xlsx"))
  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Supplementary analysis saved to:", output_file, "\n")
  
  return(invisible(wb))
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

generate_diag_excel <- function(output_dir = "Proportion/Results") {
  
  cat("\n")
  cat("=============================================================================\n")
  cat("GENERATING ALL PROPORTION DIAGNOSIS EXCEL RESULTS\n")
  cat("=============================================================================\n")
  
  # Load data
  cat("\nLoading data...\n")
  quad_data <- readRDS("Data/quad_data.rds")
  quad_data <- quad_data %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
      Region    = factor(Region, levels = c("Head", "Body", "Tail")),
      Sex       = factor(Sex, levels = c("Female", "Male")),
      Donor     = factor(Donor),
      ImageID   = factor(ImageID)
    )
  
  cat("quad_data:", nrow(quad_data), "rows,", n_distinct(quad_data$Donor), "donors\n")
  
  # Define cell types
  cell_types <- list(
    list(type = "ins",  label = "Insulin"),
    list(type = "glu",  label = "Glucagon"),
    list(type = "soma", label = "Somatostatin"),
    list(type = "pp",   label = "PP")
  )
  
  # Generate Excel files for each cell type
  for (ct in cell_types) {
    
    # Primary analysis (powered model)
    powered_path <- paste0("Proportion/Models/", ct$type, "_powered_model.rds")
    cat("\nLoading", powered_path, "...\n")
    model_powered <- readRDS(powered_path)
    generate_primary_excel(
      model = model_powered,
      data = quad_data,
      cell_type = ct$type,
      cell_label = ct$label,
      output_dir = output_dir
    )
    rm(model_powered); gc()
    
    # Supplementary analysis (mid model)
    mid_path <- paste0("Proportion/Models/", ct$type, "_mid_model.rds")
    cat("\nLoading", mid_path, "...\n")
    model_mid <- readRDS(mid_path)
    generate_supplementary_excel(
      model = model_mid,
      data = quad_data,
      cell_type = ct$type,
      cell_label = ct$label,
      output_dir = output_dir
    )
    rm(model_mid); gc()
  }
  
  # Summary
  cat("\n")
  cat("=============================================================================\n")
  cat("COMPLETE!\n")
  cat("=============================================================================\n")
  cat("\nEight Excel files created in", output_dir, ":\n\n")
  cat("PRIMARY ANALYSES (*_powered models):\n")
  cat("  1. Insulin_Content_Primary_Analysis.xlsx\n")
  cat("  2. Glucagon_Content_Primary_Analysis.xlsx\n")
  cat("  3. Somatostatin_Content_Primary_Analysis.xlsx\n")
  cat("  4. PP_Content_Primary_Analysis.xlsx\n")
  cat("\nSUPPLEMENTARY ANALYSES (*_mid models):\n")
  cat("  5. Insulin_Content_Supplementary_Exploratory.xlsx\n")
  cat("  6. Glucagon_Content_Supplementary_Exploratory.xlsx\n")
  cat("  7. Somatostatin_Content_Supplementary_Exploratory.xlsx\n")
  cat("  8. PP_Content_Supplementary_Exploratory.xlsx\n")
  cat("\n=============================================================================\n")
}

# Run if executed directly
if (interactive()) {
  cat("\nTo generate all proportion diagnosis Excel files, run:\n")
  cat("  generate_diag_excel()\n\n")
  cat("Or individually:\n")
  cat("  model <- readRDS('Proportion/Models/ins_powered_model.rds')\n")
  cat("  data  <- readRDS('Data/quad_data.rds')\n")
  cat("  generate_primary_excel(model, data, 'ins', 'Insulin')\n")
}
