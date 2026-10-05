# =============================================================================
# plot_count_adj_comparisons.R
# Post-hoc comparison plots for count-adjusted total area models.
#
# Singlets reuse existing Area/Models/Total_Area/ singlet models.
# SEO, Islet, Combined come from Area/Models/Total_Area_AdjCount/.
#
# IMPORTANT: has_islet is set to FALSE for SEO/Islet plots here. The new
# count-adjusted models do NOT have log_Islet.Cells_c as a term, so the
# islet-size-stratified plots from the existing plotting functions cannot
# be generated. (The new continuous predictor is log_count_c; if you want
# log_count_c stratified plots, that would require extending the plotting
# functions — see comments at bottom of this script.)
#
# Output: Area/Graphs/Total_Area_AdjCount/  (separate from existing Graphs)
# =============================================================================

source("Area/R/plot_area_all_comparisons.R")

area_to_count <- c(ins_area = "Count.ins", glu_area = "Count.glu",
                   soma_area = "Count.soma", pp_area = "Count.PP")

# ---- Data prep helpers ------------------------------------------------------
split_total_area_subset <- function(data, cell_col) {
  count_col <- area_to_count[[cell_col]]
  log_col   <- paste0("log_", cell_col)
  data %>%
    filter(!is.na(.data[[cell_col]]), .data[[cell_col]] > 0,
           .data[[count_col]] > 0) %>%
    mutate(!!log_col := log(.data[[cell_col]]),
           log_count   = log(.data[[count_col]]),
           log_count_c = log_count - mean(log_count))
}

split_total_area_combined <- function(data, cell_col) {
  count_col <- area_to_count[[cell_col]]
  log_col   <- paste0("log_", cell_col)
  d <- data %>%
    filter(!is.na(.data[[cell_col]]), .data[[cell_col]] > 0,
           .data[[count_col]] > 0) %>%
    mutate(!!log_col := log(.data[[cell_col]]),
           log_count   = log(.data[[count_col]]))
  multi_mean <- mean(d$log_count[d$Islet.Cells > 1])
  d %>% mutate(log_count_c = if_else(Islet.Cells == 1, 0, log_count - multi_mean))
}

split_area_singlet <- function(data, cell_col) {
  count_col <- area_to_count[[cell_col]]
  log_col   <- paste0("log_", cell_col)
  data %>%
    filter(!is.na(.data[[cell_col]]), .data[[cell_col]] > 0,
           .data[[count_col]] > 0) %>%
    mutate(!!log_col := log(.data[[cell_col]]),
           log_count_c = 0)
}

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

# ---- Main plotting driver ---------------------------------------------------
generate_all_count_adj_plots <- function(
    new_model_dir     = "Area/Models/Total_Area_AdjCount",
    singlet_model_dir = "Area/Models/Total_Area",
    output_dir        = "Area/Graphs/Total_Area_AdjCount") {

  cat("\n=============================================================================\n")
  cat("GENERATING COUNT-ADJUSTED TOTAL AREA PLOTS\n")
  cat("=============================================================================\n\n")

  cat("Loading data...\n")
  area_singlets     <- set_factors(readRDS("Data/area_data_singlets.rds"))
  area_endobs       <- set_factors(readRDS("Data/area_data_endobs.rds"))
  area_islets       <- set_factors(readRDS("Data/area_data_islets.rds"))
  area_singlets_t1d <- set_factors_t1d(readRDS("Data/area_data_singlets_t1d.rds"))
  area_endobs_t1d   <- set_factors_t1d(readRDS("Data/area_data_endobs_t1d.rds"))
  area_islets_t1d   <- set_factors_t1d(readRDS("Data/area_data_islets_t1d.rds"))

  beta_combined  <- add_object_type(set_factors(readRDS("Data/area_data_beta_combined.rds")))
  alpha_combined <- add_object_type(set_factors(readRDS("Data/area_data_alpha_combined.rds")))
  delta_combined <- add_object_type(set_factors(readRDS("Data/area_data_delta_combined.rds")))
  pp_combined    <- add_object_type(set_factors(readRDS("Data/area_data_pp_combined.rds")))
  beta_combined_t1d  <- add_object_type(set_factors_t1d(readRDS("Data/area_data_beta_combined_t1d.rds")))
  alpha_combined_t1d <- add_object_type(set_factors_t1d(readRDS("Data/area_data_alpha_combined_t1d.rds")))
  delta_combined_t1d <- add_object_type(set_factors_t1d(readRDS("Data/area_data_delta_combined_t1d.rds")))
  pp_combined_t1d    <- add_object_type(set_factors_t1d(readRDS("Data/area_data_pp_combined_t1d.rds")))
  cat("  Data loaded.\n\n")

  cell_types <- list(
    beta  = list(area_col = "ins_area",  cell_col = "ins_area",  label = "β cell area (count-adj)"),
    alpha = list(area_col = "glu_area",  cell_col = "glu_area",  label = "α cell area (count-adj)"),
    delta = list(area_col = "soma_area", cell_col = "soma_area", label = "δ cell area (count-adj)"),
    pp    = list(area_col = "pp_area",   cell_col = "pp_area",   label = "PP cell area (count-adj)")
  )

  combined_data <- list(
    beta  = list(diag = beta_combined,  t1d = beta_combined_t1d),
    alpha = list(diag = alpha_combined, t1d = alpha_combined_t1d),
    delta = list(diag = delta_combined, t1d = delta_combined_t1d),
    pp    = list(diag = pp_combined,    t1d = pp_combined_t1d)
  )

  for (ct in names(cell_types)) {
    spec <- cell_types[[ct]]; ac <- spec$area_col; cc <- spec$cell_col; lab <- spec$label
    cat("===", toupper(lab), "===\n")

    # ---- Diagnosis: Singlets (reuse existing total-area singlet model) ----
    cat("  Diag Singlets...\n")
    m <- readRDS(file.path(singlet_model_dir, paste0("area_", ct, "_1cell.rds")))
    d <- split_area_singlet(area_singlets, cc)
    generate_area_diag_plots(m, d, paste(lab, "(Singlets)"), ac,
                              has_islet = FALSE,
                              file.path(output_dir, ct, "Diagnosis", "1cell"))
    rm(m, d); gc()

    # ---- Diagnosis: SEO (count-adjusted) ----
    cat("  Diag SEO...\n")
    m <- readRDS(file.path(new_model_dir, paste0("ta_", ct, "_2to14.rds")))
    d <- split_total_area_subset(area_endobs, cc)
    generate_area_diag_plots(m, d, paste(lab, "(SEOs)"), ac,
                              has_islet = FALSE,
                              file.path(output_dir, ct, "Diagnosis", "2to14"))
    rm(m, d); gc()

    # ---- Diagnosis: Islet (count-adjusted) ----
    cat("  Diag Islet...\n")
    m <- readRDS(file.path(new_model_dir, paste0("ta_", ct, "_15plus.rds")))
    d <- split_total_area_subset(area_islets, cc)
    generate_area_diag_plots(m, d, paste(lab, "(Islets)"), ac,
                              has_islet = FALSE,
                              file.path(output_dir, ct, "Diagnosis", "15plus"))
    rm(m, d); gc()

    # ---- Diagnosis: Combined ----
    cat("  Diag Combined...\n")
    m <- readRDS(file.path(new_model_dir, paste0("ta_", ct, "_combined.rds")))
    d <- split_total_area_combined(combined_data[[ct]]$diag, cc)
    generate_area_combined_diag_plots(m, d, paste(lab, "(combined)"), ac,
                                       file.path(output_dir, ct, "Diagnosis", "combined"))
    rm(m, d); gc()

    # ---- T1D: Singlets ----
    cat("  T1D Singlets...\n")
    m <- readRDS(file.path(singlet_model_dir, paste0("area_", ct, "_1cell_t1d.rds")))
    d <- split_area_singlet(area_singlets_t1d, cc)
    generate_area_t1d_plots(m, d, paste(lab, "(Singlets)"), ac,
                             has_islet = FALSE,
                             file.path(output_dir, ct, "T1D", "1cell"))
    rm(m, d); gc()

    # ---- T1D: SEO ----
    cat("  T1D SEO...\n")
    m <- readRDS(file.path(new_model_dir, paste0("ta_", ct, "_2to14_t1d.rds")))
    d <- split_total_area_subset(area_endobs_t1d, cc)
    generate_area_t1d_plots(m, d, paste(lab, "(SEOs)"), ac,
                             has_islet = FALSE,
                             file.path(output_dir, ct, "T1D", "2to14"))
    rm(m, d); gc()

    # ---- T1D: Islet ----
    cat("  T1D Islet...\n")
    m <- readRDS(file.path(new_model_dir, paste0("ta_", ct, "_15plus_t1d.rds")))
    d <- split_total_area_subset(area_islets_t1d, cc)
    generate_area_t1d_plots(m, d, paste(lab, "(Islets)"), ac,
                             has_islet = FALSE,
                             file.path(output_dir, ct, "T1D", "15plus"))
    rm(m, d); gc()

    # ---- T1D: Combined ----
    cat("  T1D Combined...\n")
    m <- readRDS(file.path(new_model_dir, paste0("ta_", ct, "_combined_t1d.rds")))
    d <- split_total_area_combined(combined_data[[ct]]$t1d, cc)
    generate_area_combined_t1d_plots(m, d, paste(lab, "(combined)"), ac,
                                      file.path(output_dir, ct, "T1D", "combined"))
    rm(m, d); gc()
  }

  cat("\n=============================================================================\n")
  cat("PLOTTING COMPLETE.\n")
  cat("Output directory:", output_dir, "\n")
  cat("NOTE: has_islet = FALSE for all subsets — log_Islet.Cells_c plots are skipped\n")
  cat("      (the new models use log_count_c, which the plotting functions don't\n")
  cat("      currently visualize). See comments at bottom of this script for how\n")
  cat("      to extend.\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate count-adjusted plots:\n")
  cat("  source('Area/R/plot_count_adj_comparisons.R')\n")
  cat("  generate_all_count_adj_plots()\n")
}

# =============================================================================
# Extending to plot log_count_c interactions (optional, future work):
# The existing generate_area_diag_plots() / _t1d_plots() functions take a
# `has_islet` flag that switches in plots of log_Islet.Cells_c continuous
# slopes and its interactions. To get equivalent plots for log_count_c, you
# would need to either:
#   (a) Add a `has_count` flag and parallel plotting blocks in those functions,
#       referencing "log_count_c" instead of "log_Islet.Cells_c"; or
#   (b) Write a new function generate_area_count_slope_plots() that reads
#       emtrends(model, ~ Diagnosis | Region, var = "log_count_c") and plots
#       the result.
# Neither is hard, but both are out of scope for this drop-in. The Excel
# Model Summary sheet still reports the log_count_c coefficient.
# =============================================================================
