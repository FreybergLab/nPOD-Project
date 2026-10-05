# =============================================================================
# plot_pc_area_all_comparisons.R
# All post-hoc comparison plots for Per-Cell Area models
# Sources plot_area_all_comparisons.R for all plotting functions
#
# Singlet plots are SKIPPED — per-cell ≡ total area for 1-cell objects.
# See Area/Graphs for existing singlet plots.
# =============================================================================

source("Area/R/plot_area_all_comparisons.R")

# Override y-axis label for per-cell area
y_lab_area <- expression(paste("Per-cell area (", mu, "m"^2, "/cell)"))

generate_all_pc_area_plots <- function(model_dir = "Area/Models/Per_Cell_Area",
                                        output_dir = "Area/Graphs/Per_Pell_Area") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL PER-CELL AREA POST-HOC COMPARISON PLOTS\n")
  cat("=============================================================================\n\n")

  area_to_count <- c(ins_area = "Count.ins", glu_area = "Count.glu",
                     soma_area = "Count.soma", pp_area = "Count.PP")

  # --- Data helpers ---
  split_per_cell <- function(data, cell_col) {
    count_col <- area_to_count[[cell_col]]
    pc_col <- paste0(cell_col, "_pc")
    log_col <- paste0("log_", cell_col, "_pc")
    data %>%
      filter(!is.na(!!sym(cell_col)) & !!sym(cell_col) > 0 &
             !!sym(count_col) > 0) %>%
      mutate(!!pc_col := !!sym(cell_col) / !!sym(count_col),
             !!log_col := log(!!sym(cell_col) / !!sym(count_col)))
  }

  add_object_type <- function(df) {
    df %>% mutate(object_type = factor(case_when(
      Islet.Cells == 1 ~ "Singlets",
      Islet.Cells >= 2 & Islet.Cells <= 14 ~ "SEOs",
      Islet.Cells >= 15 ~ "Islets"),
      levels = c("Singlets", "SEOs", "Islets")))
  }

  set_factors <- function(df) {
    df %>% mutate(Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
                  Region = factor(Region, levels = c("Head", "Body", "Tail")),
                  Sex = factor(Sex, levels = c("Female", "Male")),
                  Donor = factor(Donor), ImageID = factor(ImageID))
  }
  set_factors_t1d <- function(df) set_factors(df) %>% droplevels()

  # --- Load base datasets ---
  cat("Loading data...\n")
  area_endobs       <- set_factors(readRDS("Data/area_data_endobs.rds"))
  area_islets       <- set_factors(readRDS("Data/area_data_islets.rds"))
  area_endobs_t1d   <- set_factors_t1d(readRDS("Data/area_data_endobs_t1d.rds"))
  area_islets_t1d   <- set_factors_t1d(readRDS("Data/area_data_islets_t1d.rds"))
  beta_combined     <- add_object_type(set_factors(readRDS("Data/area_data_beta_combined.rds")))
  alpha_combined    <- add_object_type(set_factors(readRDS("Data/area_data_alpha_combined.rds")))
  delta_combined    <- add_object_type(set_factors(readRDS("Data/area_data_delta_combined.rds")))
  pp_combined       <- add_object_type(set_factors(readRDS("Data/area_data_pp_combined.rds")))
  beta_combined_t1d  <- add_object_type(set_factors_t1d(readRDS("Data/area_data_beta_combined_t1d.rds")))
  alpha_combined_t1d <- add_object_type(set_factors_t1d(readRDS("Data/area_data_alpha_combined_t1d.rds")))
  delta_combined_t1d <- add_object_type(set_factors_t1d(readRDS("Data/area_data_delta_combined_t1d.rds")))
  pp_combined_t1d    <- add_object_type(set_factors_t1d(readRDS("Data/area_data_pp_combined_t1d.rds")))
  cat("  Data loaded.\n\n")

  cell_types <- list(
    beta  = list(area_col = "ins_area",  cell_col = "ins_area",  label = "β per-cell area"),
    alpha = list(area_col = "glu_area",  cell_col = "glu_area",  label = "α per-cell area"),
    delta = list(area_col = "soma_area", cell_col = "soma_area", label = "δ per-cell area"),
    pp    = list(area_col = "pp_area",   cell_col = "pp_area",   label = "PP per-cell area")
  )

  # Per-cell area column names (used as area_col in plotting functions)
  pc_area_cols <- list(
    beta = "ins_area_pc", alpha = "glu_area_pc",
    delta = "soma_area_pc", pp = "pp_area_pc"
  )

  combined_data <- list(
    beta = list(diag = beta_combined, t1d = beta_combined_t1d),
    alpha = list(diag = alpha_combined, t1d = alpha_combined_t1d),
    delta = list(diag = delta_combined, t1d = delta_combined_t1d),
    pp = list(diag = pp_combined, t1d = pp_combined_t1d)
  )

  for (ct in names(cell_types)) {
    spec <- cell_types[[ct]]; cc <- spec$cell_col; lab <- spec$label
    pc_col <- pc_area_cols[[ct]]
    cat("===", toupper(lab), "===\n")
    cat("  Singlets SKIPPED (per-cell ≡ total area; see Area/Graphs)\n")

    # --- DIAGNOSIS: Small EOs ---
    cat("  Diag 2to14...\n")
    d2 <- split_per_cell(area_endobs, cc)
    m2 <- readRDS(file.path(model_dir, paste0("pc_", ct, "_2to14.rds")))
    generate_area_diag_plots(m2, d2, paste(lab, "(SEOs)"), pc_col, TRUE,
                              file.path(output_dir, ct, "Diagnosis", "2to14"))
    rm(m2); gc()

    # --- DIAGNOSIS: Islets ---
    cat("  Diag 15plus...\n")
    d15 <- split_per_cell(area_islets, cc)
    m15 <- readRDS(file.path(model_dir, paste0("pc_", ct, "_15plus.rds")))
    generate_area_diag_plots(m15, d15, paste(lab, "(Islets)"), pc_col, TRUE,
                              file.path(output_dir, ct, "Diagnosis", "15plus"))
    rm(m15); gc()

    # --- DIAGNOSIS: Combined ---
    cat("  Diag combined...\n")
    mc <- readRDS(file.path(model_dir, paste0("pc_", ct, "_combined.rds")))
    dc <- split_per_cell(combined_data[[ct]]$diag, cc)
    generate_area_combined_diag_plots(mc, dc, paste(lab, "(combined)"), pc_col,
                                       file.path(output_dir, ct, "Diagnosis", "combined"))
    rm(mc); gc()

    # --- T1D: Small EOs ---
    cat("  T1D 2to14...\n")
    d2t <- split_per_cell(area_endobs_t1d, cc)
    m2t <- readRDS(file.path(model_dir, paste0("pc_", ct, "_2to14_t1d.rds")))
    generate_area_t1d_plots(m2t, d2t, paste(lab, "(SEOs)"), pc_col, TRUE,
                             file.path(output_dir, ct, "T1D", "2to14"))
    rm(m2t); gc()

    # --- T1D: Islets ---
    cat("  T1D 15plus...\n")
    d15t <- split_per_cell(area_islets_t1d, cc)
    m15t <- readRDS(file.path(model_dir, paste0("pc_", ct, "_15plus_t1d.rds")))
    generate_area_t1d_plots(m15t, d15t, paste(lab, "(Islets)"), pc_col, TRUE,
                             file.path(output_dir, ct, "T1D", "15plus"))
    rm(m15t); gc()

    # --- T1D: Combined ---
    cat("  T1D combined...\n")
    mct <- readRDS(file.path(model_dir, paste0("pc_", ct, "_combined_t1d.rds")))
    dct <- split_per_cell(combined_data[[ct]]$t1d, cc)
    generate_area_combined_t1d_plots(mct, dct, paste(lab, "(combined)"), pc_col,
                                      file.path(output_dir, ct, "T1D", "combined"))
    rm(mct); gc()
  }

  cat("\n=============================================================================\n")
  cat("COMPLETE! Per-cell area plots saved to", output_dir, "\n")
  cat("Singlet plots not generated (identical to total area — see Area/Graphs)\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all per-cell area plots:\n")
  cat("  generate_all_pc_area_plots()\n")
}
