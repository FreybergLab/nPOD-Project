# =============================================================================
# fit_count_adj_models.R
# Fits count-adjusted total area models for SEO, Islet, and Combined subsets.
# Both Diagnosis (ND vs T1D) and T1D-only (age-at-onset / disease duration).
#
# Singlet models are NOT refitted — count is identically 1 for singlets, so
# count-adjustment is a no-op. The existing Area/Models/Total_Area/area_*_1cell*.rds
# files are reused downstream.
#
# Output: Area/Models/Total_Area_AdjCount/   (does NOT touch existing model dirs)
# =============================================================================

library(glmmTMB)
library(dplyr)
library(rlang)
library(readr)

area_to_count <- list(ins_area = "Count.ins", glu_area = "Count.glu",
                      soma_area = "Count.soma", pp_area = "Count.PP")

# ---- Data prep helper -------------------------------------------------------
# Filters cell-type-positive rows; computes log(area) and log_count_c.
# For combined data (mixing singlets + multi-cell), uses hybrid centering:
#   singlets pinned at 0, multi-cell centered on multi-cell mean
#   (mirroring how log_Islet.Cells_c is centered in the existing combined data).
# For subset data (SEO or Islet only), centers within subset.
split_total_area <- function(data, cell_col, is_combined = FALSE) {
  count_col    <- area_to_count[[cell_col]]
  log_area_col <- paste0("log_", cell_col)

  d <- data %>%
    filter(!is.na(.data[[cell_col]]), .data[[cell_col]] > 0,
           .data[[count_col]] > 0) %>%
    mutate(!!log_area_col := log(.data[[cell_col]]),
           log_count       = log(.data[[count_col]]))

  if (is_combined) {
    multi_mean <- mean(d$log_count[d$Islet.Cells > 1])
    d <- d %>% mutate(
      log_count_c = if_else(Islet.Cells == 1, 0, log_count - multi_mean)
    )
  } else {
    d <- d %>% mutate(log_count_c = log_count - mean(log_count))
  }

  d
}

# ---- Formulas (clean — log_count_c appears exactly once each) ---------------

# Diagnosis SEO: log_count_c as main effect only
ta_seo_formula <- function(log_out) {
  as.formula(paste0(log_out, " ~ (Diagnosis + Region + Sex + Age_c)^2 +
    Diagnosis:Region:Age_c + log_count_c +
    (1|Donor) + (1|Donor:ImageID)"))
}

# Diagnosis Islet: log_count_c with full two-way interactions and a three-way
ta_full_formula <- function(log_out) {
  as.formula(paste0(log_out, " ~ (Diagnosis + Region + Sex + Age_c + log_count_c)^2 +
    Diagnosis:Region:log_count_c + Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID)"))
}

# Diagnosis Combined: parallel to existing combined formula but with log_count_c
# replacing log_Islet.Cells_c (which was collinear).
ta_combined_formula <- function(log_out) {
  as.formula(paste0(log_out, " ~ (Diagnosis + Region + Sex + Age_c)^2 +
    Diagnosis:Region:Age_c + log_count_c +
    object_type * Diagnosis + object_type * Region +
    (1|Donor) + (1|Donor:ImageID)"))
}

# T1D versions
ta_t1d_seo_formula <- function(log_out) {
  as.formula(paste0(log_out, " ~ (Disease.Duration_c + Region + Sex + age_at_onset_c)^2 +
    Disease.Duration_c:Region:age_at_onset_c + log_count_c +
    (1|Donor) + (1|Donor:ImageID)"))
}

ta_t1d_full_formula <- function(log_out) {
  as.formula(paste0(log_out, " ~ (Disease.Duration_c + Region + Sex + age_at_onset_c + log_count_c)^2 +
    Disease.Duration_c:Region:log_count_c + Disease.Duration_c:Region:age_at_onset_c +
    (1|Donor) + (1|Donor:ImageID)"))
}

ta_t1d_combined_formula <- function(log_out) {
  as.formula(paste0(log_out, " ~ (Disease.Duration_c + Region + Sex + age_at_onset_c)^2 +
    Disease.Duration_c:Region:age_at_onset_c + log_count_c +
    object_type * Disease.Duration_c + object_type * Region +
    (1|Donor) + (1|Donor:ImageID)"))
}

# ---- Factor / object_type helpers (match existing code) ---------------------
set_factors <- function(df) {
  df %>% mutate(
    Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
    Region    = factor(Region, levels = c("Head", "Body", "Tail")),
    Sex       = factor(Sex, levels = c("Female", "Male")),
    Donor     = factor(Donor),
    ImageID   = factor(ImageID))
}
set_factors_t1d <- function(df) set_factors(df) %>% droplevels()

add_object_type <- function(df) {
  df %>% mutate(
    object_type = factor(case_when(
      Islet.Cells == 1 ~ "Singlets",
      Islet.Cells >= 2 & Islet.Cells <= 14 ~ "SEOs",
      Islet.Cells >= 15 ~ "Islets"),
      levels = c("Singlets", "SEOs", "Islets")))
}

# ---- Main fit driver --------------------------------------------------------
fit_all_count_adj_models <- function(output_dir = "Area/Models/Total_Area_AdjCount") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n=============================================================================\n")
  cat("FITTING COUNT-ADJUSTED TOTAL AREA MODELS\n")
  cat("=============================================================================\n\n")
  cat("Output directory:", output_dir, "\n\n")

  cat("Loading data...\n")
  area_endobs     <- set_factors(readRDS("Data/area_data_endobs.rds"))
  area_islets     <- set_factors(readRDS("Data/area_data_islets.rds"))
  area_endobs_t1d <- set_factors_t1d(readRDS("Data/area_data_endobs_t1d.rds"))
  area_islets_t1d <- set_factors_t1d(readRDS("Data/area_data_islets_t1d.rds"))

  beta_combined  <- add_object_type(set_factors(readRDS("Data/area_data_beta_combined.rds")))
  alpha_combined <- add_object_type(set_factors(readRDS("Data/area_data_alpha_combined.rds")))
  delta_combined <- add_object_type(set_factors(readRDS("Data/area_data_delta_combined.rds")))
  pp_combined    <- add_object_type(set_factors(readRDS("Data/area_data_pp_combined.rds")))

  beta_combined_t1d  <- add_object_type(set_factors_t1d(readRDS("Data/area_data_beta_combined_t1d.rds")))
  alpha_combined_t1d <- add_object_type(set_factors_t1d(readRDS("Data/area_data_alpha_combined_t1d.rds")))
  delta_combined_t1d <- add_object_type(set_factors_t1d(readRDS("Data/area_data_delta_combined_t1d.rds")))
  pp_combined_t1d    <- add_object_type(set_factors_t1d(readRDS("Data/area_data_pp_combined_t1d.rds")))
  cat("  Data loaded.\n\n")

  combined_data <- list(
    ins  = list(diag = beta_combined,  t1d = beta_combined_t1d),
    glu  = list(diag = alpha_combined, t1d = alpha_combined_t1d),
    soma = list(diag = delta_combined, t1d = delta_combined_t1d),
    pp   = list(diag = pp_combined,    t1d = pp_combined_t1d)
  )

  cell_specs <- list(
    list(type = "ins",  short = "beta",  area = "ins_area"),
    list(type = "glu",  short = "alpha", area = "glu_area"),
    list(type = "soma", short = "delta", area = "soma_area"),
    list(type = "pp",   short = "pp",    area = "pp_area")
  )

  # Bumped iteration limits — prevents the convergence warning we saw on alpha islet
  fit_ctrl <- glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4))

  fit_one <- function(formula_fn, data, log_area_col, save_path) {
    m <- glmmTMB(formula_fn(log_area_col), data = data,
                 family = gaussian(), control = fit_ctrl)
    saveRDS(m, save_path)
    invisible(NULL)
  }

  for (cs in cell_specs) {
    cat("===", toupper(cs$short), "===\n")
    log_area_col <- paste0("log_", cs$area)

    # ---- Diagnosis ----
    cat("  Diag SEO...\n")
    fit_one(ta_seo_formula,
            split_total_area(area_endobs, cs$area, is_combined = FALSE),
            log_area_col,
            file.path(output_dir, paste0("ta_", cs$short, "_2to14.rds")))
    gc()

    cat("  Diag Islet...\n")
    fit_one(ta_full_formula,
            split_total_area(area_islets, cs$area, is_combined = FALSE),
            log_area_col,
            file.path(output_dir, paste0("ta_", cs$short, "_15plus.rds")))
    gc()

    cat("  Diag Combined...\n")
    fit_one(ta_combined_formula,
            split_total_area(combined_data[[cs$type]]$diag, cs$area, is_combined = TRUE),
            log_area_col,
            file.path(output_dir, paste0("ta_", cs$short, "_combined.rds")))
    gc()

    # ---- T1D ----
    cat("  T1D SEO...\n")
    fit_one(ta_t1d_seo_formula,
            split_total_area(area_endobs_t1d, cs$area, is_combined = FALSE),
            log_area_col,
            file.path(output_dir, paste0("ta_", cs$short, "_2to14_t1d.rds")))
    gc()

    cat("  T1D Islet...\n")
    fit_one(ta_t1d_full_formula,
            split_total_area(area_islets_t1d, cs$area, is_combined = FALSE),
            log_area_col,
            file.path(output_dir, paste0("ta_", cs$short, "_15plus_t1d.rds")))
    gc()

    cat("  T1D Combined...\n")
    fit_one(ta_t1d_combined_formula,
            split_total_area(combined_data[[cs$type]]$t1d, cs$area, is_combined = TRUE),
            log_area_col,
            file.path(output_dir, paste0("ta_", cs$short, "_combined_t1d.rds")))
    gc()

    cat("\n")
  }

  cat("=============================================================================\n")
  cat("MODEL FITTING COMPLETE.\n")
  cat("Files saved to:", output_dir, "\n")
  cat("  24 models total (4 cell types × 3 subsets × 2 conditions [Diag, T1D])\n")
  cat("  Singlets are NOT here — reuse Area/Models/Total_Area/area_*_1cell*.rds\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo fit all count-adjusted models:\n")
  cat("  source('Area/R/fit_count_adj_models.R')\n")
  cat("  fit_all_count_adj_models()\n")
}
