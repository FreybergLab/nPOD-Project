# =============================================================================
# generate_count_adj_excel.R
# Excel generation for count-adjusted total area models — Diagnosis AND T1D.
#
# Reuses generate_cell_type_excel() (Diagnosis) and generate_area_t1d_excel()
# (T1D) workhorses, pointed at the new count-adjusted models.
#
# Singlets reuse existing singlet models from Area/Models/Total_Area/ — for a
# singlet, count = 1, log_count_c = 0, so count adjustment is a no-op.
#
# Output: Area/Results/Total_Area_AdjCount/
#   Diagnosis:  {label}_AdjCount_Cell_Area_Analysis.xlsx
#   T1D:        {label}_AdjCount_Cell_Area_AO_Analysis.xlsx
# (Distinct filenames from existing files — safe to run alongside.)
# =============================================================================

source("Area/R/generate_area_excel_diag.R")
source("Area/R/generate_area_t1d_excel.R")

area_to_count <- c(ins_area = "Count.ins", glu_area = "Count.glu",
                   soma_area = "Count.soma", pp_area = "Count.PP")

# ---- Data prep helpers (mirror the fit script) ------------------------------
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

# Singlet model doesn't reference log_count_c, but we add it as 0 so the data
# frame has a consistent shape with the other subsets (defensive).
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

cell_defs <- list(
  list(type = "ins",  short = "beta",  label = "β_AdjCount",  orig_area = "ins_area"),
  list(type = "glu",  short = "alpha", label = "α_AdjCount",  orig_area = "glu_area"),
  list(type = "soma", short = "delta", label = "δ_AdjCount",  orig_area = "soma_area"),
  list(type = "pp",   short = "pp",    label = "PP_AdjCount", orig_area = "pp_area")
)

# ---- DIAGNOSIS wrapper ------------------------------------------------------
generate_all_count_adj_excel_diag <- function(
    output_dir       = "Area/Results/Total_Area_AdjCount/",
    new_model_dir    = "Area/Models/Total_Area_AdjCount/",
    singlet_model_dir = "Area/Models/Total_Area/") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n=============================================================================\n")
  cat("GENERATING COUNT-ADJUSTED EXCEL — DIAGNOSIS\n")
  cat("=============================================================================\n\n")

  area_singlets <- set_factors(readRDS("Data/area_data_singlets.rds"))
  area_endobs   <- set_factors(readRDS("Data/area_data_endobs.rds"))
  area_islets   <- set_factors(readRDS("Data/area_data_islets.rds"))

  combined_data_map <- list(
    ins  = add_object_type(set_factors(readRDS("Data/area_data_beta_combined.rds"))),
    glu  = add_object_type(set_factors(readRDS("Data/area_data_alpha_combined.rds"))),
    soma = add_object_type(set_factors(readRDS("Data/area_data_delta_combined.rds"))),
    pp   = add_object_type(set_factors(readRDS("Data/area_data_pp_combined.rds")))
  )

  for (cd in cell_defs) {
    cat("\nProcessing", cd$label, "...\n")

    m1  <- readRDS(file.path(singlet_model_dir, paste0("area_", cd$short, "_1cell.rds")))
    m2  <- readRDS(file.path(new_model_dir,     paste0("ta_",   cd$short, "_2to14.rds")))
    m15 <- readRDS(file.path(new_model_dir,     paste0("ta_",   cd$short, "_15plus.rds")))
    mc  <- readRDS(file.path(new_model_dir,     paste0("ta_",   cd$short, "_combined.rds")))

    d1  <- split_area_singlet(area_singlets, cd$orig_area)
    d2  <- split_total_area_subset(area_endobs, cd$orig_area)
    d15 <- split_total_area_subset(area_islets, cd$orig_area)
    dc  <- split_total_area_combined(combined_data_map[[cd$type]], cd$orig_area)

    generate_cell_type_excel(
      singlet_model = m1, multicell_model = m2,
      islet_model   = m15, combined_model = mc,
      singlet_data  = d1, multicell_data  = d2,
      islet_data    = d15, combined_data  = dc,
      cell_type     = cd$type, cell_label = cd$label,
      area_col      = cd$orig_area,
      output_dir    = output_dir
    )
    rm(m1, m2, m15, mc); gc()
  }

  cat("\nDIAGNOSIS COMPLETE. Files in", output_dir, ":\n")
  for (cd in cell_defs) cat("  ", cd$label, "_Cell_Area_Analysis.xlsx\n")
}

# ---- T1D wrapper ------------------------------------------------------------
generate_all_count_adj_excel_t1d <- function(
    output_dir       = "Area/Results/Total_Area_AdjCount/",
    new_model_dir    = "Area/Models/Total_Area_AdjCount/",
    singlet_model_dir = "Area/Models/Total_Area/") {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n=============================================================================\n")
  cat("GENERATING COUNT-ADJUSTED EXCEL — T1D (AGE AT ONSET)\n")
  cat("=============================================================================\n\n")

  area_singlets_t1d <- set_factors_t1d(readRDS("Data/area_data_singlets_t1d.rds"))
  area_endobs_t1d   <- set_factors_t1d(readRDS("Data/area_data_endobs_t1d.rds"))
  area_islets_t1d   <- set_factors_t1d(readRDS("Data/area_data_islets_t1d.rds"))

  combined_data_map <- list(
    ins  = add_object_type(set_factors_t1d(readRDS("Data/area_data_beta_combined_t1d.rds"))),
    glu  = add_object_type(set_factors_t1d(readRDS("Data/area_data_alpha_combined_t1d.rds"))),
    soma = add_object_type(set_factors_t1d(readRDS("Data/area_data_delta_combined_t1d.rds"))),
    pp   = add_object_type(set_factors_t1d(readRDS("Data/area_data_pp_combined_t1d.rds")))
  )

  for (cd in cell_defs) {
    cat("\nProcessing", cd$label, "...\n")

    m1  <- readRDS(file.path(singlet_model_dir, paste0("area_", cd$short, "_1cell_t1d.rds")))
    m2  <- readRDS(file.path(new_model_dir,     paste0("ta_",   cd$short, "_2to14_t1d.rds")))
    m15 <- readRDS(file.path(new_model_dir,     paste0("ta_",   cd$short, "_15plus_t1d.rds")))
    mc  <- readRDS(file.path(new_model_dir,     paste0("ta_",   cd$short, "_combined_t1d.rds")))

    d1  <- split_area_singlet(area_singlets_t1d, cd$orig_area)
    d2  <- split_total_area_subset(area_endobs_t1d, cd$orig_area)
    d15 <- split_total_area_subset(area_islets_t1d, cd$orig_area)
    dc  <- split_total_area_combined(combined_data_map[[cd$type]], cd$orig_area)

    generate_area_t1d_excel(
      singlet_model = m1, multicell_model = m2,
      islet_model   = m15, combined_model = mc,
      singlet_data  = d1, multicell_data  = d2,
      islet_data    = d15, combined_data  = dc,
      cell_type     = cd$type, cell_label = cd$label,
      area_col      = cd$orig_area,
      model_type    = "ao",
      output_dir    = output_dir
    )
    rm(m1, m2, m15, mc); gc()
  }

  cat("\nT1D COMPLETE. Files in", output_dir, ":\n")
  for (cd in cell_defs) cat("  ", cd$label, "_Cell_Area_AO_Analysis.xlsx\n")
}

# ---- Convenience: run both --------------------------------------------------
generate_all_count_adj_excel <- function(
    output_dir = "Area/Results/Total_Area_AdjCount/") {
  generate_all_count_adj_excel_diag(output_dir = output_dir)
  generate_all_count_adj_excel_t1d(output_dir = output_dir)
  cat("\n=============================================================================\n")
  cat("ALL COUNT-ADJUSTED EXCEL FILES COMPLETE.\n")
  cat("Output directory:", output_dir, "\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all count-adjusted Excel files:\n")
  cat("  source('Area/R/generate_count_adj_excel.R')\n")
  cat("  generate_all_count_adj_excel()\n")
}
