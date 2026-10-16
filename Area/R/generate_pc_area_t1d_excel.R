# =============================================================================
# Generate Excel Results for Per-Cell Area T1D Models — AO parameterization
# Reuses generate_area_t1d_excel.R main function with per-cell models/data
#
# Singlet models are reused from Area/Models (per-cell ≡ total for singlets).
# SEO and Islet models come from PerCellArea/Models.
# =============================================================================

source("Area/R/generate_area_t1d_excel.R")

generate_all_pc_area_t1d_excel <- function(output_dir = "Area/Results/Per_Cell_Area") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING ALL PER-CELL AREA T1D EXCEL RESULTS (AO parameterization)\n")
  cat("=============================================================================\n\n")
  
  # --- Helper: same split_per_cell as models ---
  area_to_count <- c(ins_area = "Count.ins", glu_area = "Count.glu",
                     soma_area = "Count.soma", pp_area = "Count.PP")
  
  split_per_cell <- function(data, cell_col) {
    count_col <- area_to_count[[cell_col]]
    pc_col    <- paste0(cell_col, "_pc")
    log_col   <- paste0("log_", cell_col, "_pc")
    data %>%
      filter(!is.na(!!sym(cell_col)) & !!sym(cell_col) > 0 &
             !!sym(count_col) > 0) %>%
      mutate(!!pc_col  := !!sym(cell_col) / !!sym(count_col),
             !!log_col := log(!!sym(cell_col) / !!sym(count_col)))
  }
  
  # For singlets, per-cell ≡ total area
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
  area_singlets_t1d  <- set_factors(readRDS("Data/area_data_singlets_t1d.rds"))
  area_endobs_t1d    <- set_factors(readRDS("Data/area_data_endobs_t1d.rds"))
  area_islets_t1d    <- set_factors(readRDS("Data/area_data_islets_t1d.rds"))
  
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
  
  cell_defs <- list(
    list(type = "ins",  label = "β_PC",  area_col = "ins_area_pc",
         prefix_1 = "area_beta_1cell",
         prefix_2 = "pc_beta_2to14",
         prefix_15 = "pc_beta_15plus",
         prefix_comb = "pc_beta_combined",
         orig_area = "ins_area"),
    list(type = "glu",  label = "α_PC", area_col = "glu_area_pc",
         prefix_1 = "area_alpha_1cell",
         prefix_2 = "pc_alpha_2to14",
         prefix_15 = "pc_alpha_15plus",
         prefix_comb = "pc_alpha_combined",
         orig_area = "glu_area"),
    list(type = "soma", label = "δ_PC", area_col = "soma_area_pc",
         prefix_1 = "area_delta_1cell",
         prefix_2 = "pc_delta_2to14",
         prefix_15 = "pc_delta_15plus",
         prefix_comb = "pc_delta_combined",
         orig_area = "soma_area"),
    list(type = "pp",   label = "PP_PC",    area_col = "pp_area_pc",
         prefix_1 = "area_pp_1cell",
         prefix_2 = "pc_pp_2to14",
         prefix_15 = "pc_pp_15plus",
         prefix_comb = "pc_pp_combined",
         orig_area = "pp_area")
  )
  
  cat("--- AGE AT ONSET MODELS (T1D, Per-Cell Area) ---\n")
  for (cd in cell_defs) {
    cat("\nProcessing", cd$label, "...\n")
    
    # Singlet from Area/Models (reuse), rest from PerCellArea/Models
    m1  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_1, "_t1d.rds"))
    m2  <- readRDS(paste0("Area/Models/Per_Cell_Area/", cd$prefix_2, "_t1d.rds"))
    m15 <- readRDS(paste0("Area/Models/Per_Cell_Area/", cd$prefix_15, "_t1d.rds"))
    mc  <- readRDS(paste0("Area/Models/Per_Cell_Area/", cd$prefix_comb, "_t1d.rds"))
    
    # Singlet data uses original area col from singlets dataset
    d1  <- split_area(area_singlets_t1d, cd$orig_area)
    # SEO and Islet from subset-specific datasets
    d2  <- split_per_cell(area_endobs_t1d, cd$orig_area)
    d15 <- split_per_cell(area_islets_t1d, cd$orig_area)
    # Combined from cell-type-specific combined dataset
    dc  <- split_per_cell(combined_data_map[[cd$type]], cd$orig_area)
    
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
  cat("COMPLETE! Four per-cell area T1D Excel files created in", output_dir, ":\n\n")
  for (cd in cell_defs) cat("  ", cd$label, "_Cell_Area_AO_Analysis.xlsx\n")
  cat("\nNote: Singlet results reuse total area models (per-cell ≡ total for 1-cell objects)\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all per-cell area T1D Excel files, run:\n")
  cat("  generate_all_pc_area_t1d_excel()\n")
}
