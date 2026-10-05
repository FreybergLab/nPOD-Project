# =============================================================================
# finish_pp_area.R
#
# Recovery script for the audit pipeline. plot_area_all_comparisons.R failed
# right before PP T1D 15plus due to Windows's C-runtime stdio handle limit
# (Cairo/PDF font handles accumulate across 400+ plots). All other cell types
# already finished and their files are saved in the quarantine.
#
# This script:
#   1. Generates ONLY the two missing PP plot batches (T1D 15plus + T1D combined)
#   2. Marks plot_area_all_comparisons.R::generate_all_area_plots as complete
#      in the audit state file, so re-sourcing the main audit skips it.
#
# Usage (in a fresh R session — close & reopen RStudio first):
#   setwd("F:/Work/nPod_final_v4")
#   source("finish_pp_area.R")
#   # ...then to continue with the rest of the audit:
#   source("audit_pipeline_outputs.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmmTMB)
  library(brms)
})

LIVE       <- getwd()
QUARANTINE <- file.path(LIVE, "_audit_quarantine")
STATE_FILE <- file.path(QUARANTINE, "_audit_state.txt")

if (!dir.exists(QUARANTINE)) stop("Quarantine not found at ", QUARANTINE)

# Work inside the quarantine so all relative paths resolve correctly
old_wd <- setwd(QUARANTINE)
on.exit(setwd(old_wd), add = TRUE)

# Source the area plot script — gives us generate_area_t1d_plots and
# generate_area_combined_t1d_plots at top level.
script_path <- file.path(LIVE, "Area", "R", "plot_area_all_comparisons.R")
if (!file.exists(script_path)) stop("Not found: ", script_path)
source(script_path)

# Replicate the inner helpers (defined inside generate_all_area_plots, so not
# accessible after sourcing — small enough to duplicate cleanly).
split_area <- function(data, cell_col) {
  area_to_count <- c(ins_area = "Count.ins", glu_area = "Count.glu",
                     soma_area = "Count.soma", pp_area = "Count.PP")
  log_col   <- paste0("log_", cell_col)
  count_col <- area_to_count[[cell_col]]
  data %>%
    filter(!is.na(.data[[cell_col]]) & .data[[cell_col]] > 0 & .data[[count_col]] > 0) %>%
    mutate(!!log_col := log(.data[[cell_col]]))
}

set_factors <- function(df) {
  df %>% mutate(
    Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
    Region    = factor(Region,    levels = c("Head", "Body", "Tail")),
    Sex       = factor(Sex,       levels = c("Female", "Male")),
    Donor     = factor(Donor),
    ImageID   = factor(ImageID))
}
set_factors_t1d <- function(df) set_factors(df) %>% droplevels()

add_object_type <- function(df) {
  df %>% mutate(
    object_type = factor(case_when(
      Islet.Cells == 1                       ~ "Singlets",
      Islet.Cells >= 2 & Islet.Cells <= 14   ~ "SEOs",
      Islet.Cells >= 15                      ~ "Islets"),
      levels = c("Singlets", "SEOs", "Islets")))
}

cat("Loading PP T1D data...\n")
area_islets_t1d <- set_factors_t1d(readRDS("Data/area_data_islets_t1d.rds"))
pp_combined_t1d <- add_object_type(set_factors_t1d(readRDS("Data/area_data_pp_combined_t1d.rds")))

OUT_DIR   <- "Area/Graphs/total_area/pp/T1D"
MODEL_DIR <- "Area/Models/Total_Area"

# ----- PP T1D 15plus -----
cat("\n=== PP T1D 15plus ===\n")
m <- readRDS(file.path(MODEL_DIR, "area_pp_15plus_t1d.rds"))
dat <- split_area(area_islets_t1d, "pp_area")
generate_area_t1d_plots(m, dat,
                        cell_label = "PP cell area (Islets)",
                        area_col   = "pp_area",
                        has_islet  = TRUE,
                        output_dir = file.path(OUT_DIR, "15plus"))
rm(m, dat); closeAllConnections(); graphics.off(); gc(verbose = FALSE)

# ----- PP T1D combined -----
cat("\n=== PP T1D combined ===\n")
m <- readRDS(file.path(MODEL_DIR, "pp_combined_t1d.rds"))
generate_area_combined_t1d_plots(m, split_area(pp_combined_t1d, "pp_area"),
                                 cell_label = "PP cell area (combined)",
                                 area_col   = "pp_area",
                                 output_dir = file.path(OUT_DIR, "combined"))
rm(m); closeAllConnections(); graphics.off(); gc(verbose = FALSE)

# ----- Mark the area entry as complete in the state file -----
entry_id <- "plot_area_all_comparisons.R::generate_all_area_plots"
existing <- if (file.exists(STATE_FILE)) readLines(STATE_FILE) else character(0)
if (!entry_id %in% existing) {
  cat(entry_id, "\n", file = STATE_FILE, sep = "", append = TRUE)
  cat("\nMarked ", entry_id, " as complete in ", STATE_FILE, "\n", sep = "")
} else {
  cat("\n", entry_id, " was already marked complete.\n", sep = "")
}

cat("\nDone. Now re-source audit_pipeline_outputs.R to finish the remaining entries.\n")
