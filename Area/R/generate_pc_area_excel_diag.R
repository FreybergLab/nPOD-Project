# =============================================================================
# Generate Excel Results for Per-Cell Area Models — Diagnosis
# Reuses generate_area_excel_diag.R main function with per-cell models/data
#
# Singlet models are reused from Area/Models (per-cell ≡ total for singlets).
# SEO and Islet models come from PerCellArea/Models.
# =============================================================================

source("Area/R/generate_area_excel_diag.R")

generate_all_pc_area_excel <- function(output_dir = "Area/Results/Per_Cell_Area/") {
  
  cat("\n=============================================================================\n")
  cat("GENERATING ALL PER-CELL AREA EXCEL RESULTS (DIAGNOSIS)\n")
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
  
  # For singlets, per-cell ≡ total area; use same split_area as total area models
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
  area_singlets  <- set_factors(readRDS("Data/area_data_singlets.rds"))
  area_endobs    <- set_factors(readRDS("Data/area_data_endobs.rds"))
  area_islets    <- set_factors(readRDS("Data/area_data_islets.rds"))
  
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
  # Singlet models: reuse from Area/Models (per-cell ≡ total area for singlets)
  # SEO and Islet models: from PerCellArea/Models
  cell_defs <- list(
    list(type = "ins",  label = "β_PC",  area_col = "ins_area_pc",
         prefix_1 = "area_beta_1cell",    # Area/Models (reuse)
         prefix_2 = "pc_beta_2to14",      # PerCellArea/Models
         prefix_15 = "pc_beta_15plus",    # PerCellArea/Models
         prefix_comb = "pc_beta_combined", # PerCellArea/Models
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
  
  for (cd in cell_defs) {
    cat("\nProcessing", cd$label, "...\n")
    
    # Singlet model from Area/Models (reuse — per-cell ≡ total for singlets)
    m1  <- readRDS(paste0("Area/Models/Total_Area/", cd$prefix_1, ".rds"))
    # SEO and Islet from PerCellArea/Models
    m2  <- readRDS(paste0("Area/Models/Per_Cell_Area/", cd$prefix_2, ".rds"))
    m15 <- readRDS(paste0("Area/Models/Per_Cell_Area/", cd$prefix_15, ".rds"))
    mc  <- readRDS(paste0("Area/Models/Per_Cell_Area/", cd$prefix_comb, ".rds"))
    
    # Singlet data — uses original area col (identical to per-cell for singlets)
    d1  <- split_area(area_singlets, cd$orig_area)
    # SEO and Islet data — per-cell area from subset-specific datasets
    d2  <- split_per_cell(area_endobs, cd$orig_area)
    d15 <- split_per_cell(area_islets, cd$orig_area)
    # Combined data
    dc  <- split_per_cell(combined_data_map[[cd$type]], cd$orig_area)
    
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
  cat("COMPLETE! Four per-cell area Excel files created in", output_dir, "\n")
  cat("Note: Singlet results reuse total area models (per-cell ≡ total for 1-cell objects)\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all per-cell area Diagnosis Excel files, run:\n")
  cat("  generate_all_pc_area_excel()\n")
}
