# =============================================================================
# marker_common.R
# Shared registry, data loading, model resolution and a preflight check.
# Sourced by the other Marker/R scripts — source this one first if you are
# working interactively.
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)

# Where rate_models.Rmd cached its fits. Change here if you ran that file from
# a different working directory.
MARKER_MODEL_DIR       <- "Marker_Rate_Models"   # rate_models.Rmd
MARKER_COUNT_MODEL_DIR <- "Marker_Models_v6"     # marker_models_v6.Rmd
PER <- 100   # report rates per 100 cells (= percent)

# The two model families answer different questions and are reported side by
# side, not as alternatives:
#
#   RATE  (offset(log(denominator)))  "of the cells that could be positive,
#                                      what fraction are?"   <- primary
#   COUNT (log_Islet.Cells_c covariate) "how many positive cells per islet?"
#                                                            <- abundance
#
# Reading the count models as rates is what produced the four reversed
# conclusions: Ki-67+ beta cells per islet fell 28x while the per-beta-cell
# rate was unchanged (p = 0.88). Both numbers are true; only one is a rate.
#
# The twelve T1D duration count models are deliberately NOT reported: every
# adequately powered duration term was null, the one borderline term flipped
# its 0.05 verdict in 11 of 21 leave-one-donor-out folds, and the insulin and
# acinar variants fit 25 fixed effects to 639-1,238 non-zero islets.

marker_defs <- list(
  list(key = "ki67_ins", marker = "Ki-67", cell = "beta",   label = "Ki67_Beta",
       num = "insulin_cells",  den = "Count.ins"),
  list(key = "ki67_glu", marker = "Ki-67", cell = "alpha",  label = "Ki67_Alpha",
       num = "glucagon_cells", den = "Count.glu"),
  list(key = "ki67_aci", marker = "Ki-67", cell = "acinar", label = "Ki67_Acinar",
       num = "acinar_cells",   den = "Count.Acinar"),
  list(key = "cd3_ins",  marker = "CD3",   cell = "beta",   label = "CD3_Beta",
       num = "insulin_cells",  den = "Count.ins"),
  list(key = "cd3_glu",  marker = "CD3",   cell = "alpha",  label = "CD3_Alpha",
       num = "glucagon_cells", den = "Count.glu"),
  list(key = "cd3_aci",  marker = "CD3",   cell = "acinar", label = "CD3_Acinar",
       num = "acinar_cells",   den = "Count.Acinar")
)

marker_set_factors <- function(df) {
  df %>% mutate(
    Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
    Region    = factor(Region, levels = c("Head", "Body", "Tail")),
    Sex       = factor(Sex, levels = c("Female", "Male")),
    Donor     = factor(Donor),
    ImageID   = factor(ImageID))
}

# Full marker datasets, Single-stain removed to match ratio_data in the count
# pipeline (a Single run detects insulin only, so Count.glu is structurally 0).
load_marker_data <- function(stain_handling = "double_only") {
  out <- list()
  for (tag in c("ki67", "cd3")) {
    p <- sprintf("Data/%s_data.rds", tag)
    if (!file.exists(p)) stop("missing ", p, " — run Create_marker_datasets_v1.Rmd")
    d <- marker_set_factors(readRDS(p))
    if (identical(stain_handling, "double_only")) {
      d <- droplevels(d[!is.na(d$Stain) & d$Stain == "Double", , drop = FALSE])
    }
    out[[tag]] <- d
  }
  out
}

# Islets with a usable denominator, which is the data each model was fitted to.
marker_subset <- function(full, def) {
  sub <- droplevels(full[!is.na(full[[def$den]]) & full[[def$den]] > 0, ,
                         drop = FALSE])
  attr(sub, "full_data") <- full
  sub
}

# -----------------------------------------------------------------------------
# Model resolution
# -----------------------------------------------------------------------------
# rate_models.Rmd caches EVERY candidate family during selection, so the
# directory contains e.g. cd3_aci_rate_nb1.rds, _nb2.rds and _zinb1.rds side by
# side. Picking by filename order would load whichever came first rather than
# the family that won, so the selected family is read from rate_results.rds.
# If that file is absent, fall back to choosing by AIC rather than by order.
# mode = "rate"  -> <key>_rate_<fam>.rds, selected family in rate_results.rds
# mode = "count" -> <key>_<fam>.rds,       selected family in selection_v6.rds
load_marker_model <- function(def, mode = c("rate", "count"), model_dir = NULL) {
  mode <- match.arg(mode)
  if (is.null(model_dir)) {
    model_dir <- if (mode == "rate") MARKER_MODEL_DIR else MARKER_COUNT_MODEL_DIR
  }
  id  <- if (mode == "rate") paste0(def$key, "_rate") else def$key
  sel <- file.path(model_dir,
                   if (mode == "rate") "rate_results.rds" else "selection_v6.rds")

  if (file.exists(sel)) {
    res <- tryCatch(readRDS(sel), error = function(e) NULL)
    fam <- if (!is.null(res$families) && id %in% names(res$families))
      res$families[[id]] else NULL
    if (!is.null(fam) && !is.na(fam)) {
      p <- file.path(model_dir, sprintf("%s_%s.rds", id, fam))
      if (file.exists(p)) {
        return(list(model = readRDS(p), family = fam, how = "selected"))
      }
    }
  }

  cands <- list()
  for (f in c("nb1", "nb2", "zinb1", "zinb2")) {
    p <- file.path(model_dir, sprintf("%s_%s.rds", id, f))
    if (file.exists(p)) cands[[f]] <- tryCatch(readRDS(p), error = function(e) NULL)
  }
  cands <- cands[!vapply(cands, is.null, logical(1))]
  if (!length(cands)) return(NULL)
  aics <- vapply(cands, function(m) tryCatch(AIC(m), error = function(e) NA_real_), 0)
  if (all(is.na(aics))) return(NULL)
  best <- names(cands)[which.min(aics)]
  list(model = cands[[best]], family = best, how = "AIC fallback")
}

# Back-compatible alias.
load_rate_model <- function(def, model_dir = MARKER_MODEL_DIR) {
  load_marker_model(def, "rate", model_dir)
}

# -----------------------------------------------------------------------------
# Diagnosis effect — DELIBERATELY NOT DEFINED HERE
# -----------------------------------------------------------------------------
# marker_diagnosis_effect() used to live here, built on
# emmeans -> regrid(transform = "log") -> contrast(type = "response").
# That chain returned marginal rates of 4.65-5.07% for all six outcomes,
# including 4.70% for Ki-67 in alpha cells where the observed rate is 0.62%,
# and ratios that were just the quotient of the two wrong means. It has been
# removed rather than fixed so nothing can call it by accident.
#
# Use marginal_diagnosis() from Marker/R/marker_effects.R. That computes the
# same estimand as an exact linear contrast on fixef() and vcov(), and
# validate_marker_rates() checks every result against the observed rates.

# -----------------------------------------------------------------------------
# Preflight: run this first if anything produces no output
# -----------------------------------------------------------------------------
marker_preflight <- function(model_dir = MARKER_MODEL_DIR,
                             count_model_dir = MARKER_COUNT_MODEL_DIR) {
  cat("\n=============================================================================\n")
  cat("MARKER PREFLIGHT\n")
  cat("=============================================================================\n\n")
  cat("working directory:", getwd(), "\n")
  cat("  (must be the PROJECT ROOT — the folder containing Data/ and Marker/)\n\n")

  ok <- TRUE
  chk <- function(path, what) {
    e <- file.exists(path)
    cat(sprintf("  %-44s %s\n", path, if (e) "found" else paste("MISSING —", what)))
    if (!e) ok <<- FALSE
    e
  }

  cat("Required inputs:\n")
  chk("colors_master.R", "needed by every plot script")
  chk("Data/ki67_data.rds", "run Create_marker_datasets_v1.Rmd")
  chk("Data/cd3_data.rds",  "run Create_marker_datasets_v1.Rmd")
  chk("Data/ratio_data.rds",   "run Create_count_datasets_v1.Rmd")
  chk("Data/density_data.rds", "run Create_count_datasets_v1.Rmd")

  cat("\nModel directory:", model_dir, "\n")
  if (!dir.exists(model_dir)) {
    cat("  MISSING — run rate_models.Rmd, or set MARKER_MODEL_DIR to where it\n")
    cat("  wrote its cache (it uses MODEL_DIR <- \"Marker_Rate_Models\" relative\n")
    cat("  to wherever that file was knitted).\n")
    ok <- FALSE
  } else {
    n <- length(list.files(model_dir, pattern = "\\.rds$"))
    cat(sprintf("  %d .rds files present\n", n))
    cat(sprintf("  rate_results.rds: %s\n",
                if (file.exists(file.path(model_dir, "rate_results.rds")))
                  "found (selected families will be used)" else
                  "MISSING (will fall back to choosing by AIC)"))
    cat("\n  Model resolution:\n")
    for (def in marker_defs) {
      r <- load_rate_model(def, model_dir)
      if (is.null(r)) {
        cat(sprintf("    %-12s NOT FOUND\n", def$label)); ok <- FALSE
      } else {
        cat(sprintf("    %-12s %-6s (%s)  AIC %10.1f\n", def$label, r$family,
                    r$how, AIC(r$model)))
      }
    }
  }

  cat("\nCount model directory (abundance results):", count_model_dir, "\n")
  if (!dir.exists(count_model_dir)) {
    cat("  not present — the count workbooks will be skipped. These are the\n")
    cat("  marker_models_v6.Rmd fits; optional, since the rate models are the\n")
    cat("  primary results.\n")
  } else {
    for (def in marker_defs) {
      r <- load_marker_model(def, "count", count_model_dir)
      if (is.null(r)) {
        cat(sprintf("    %-12s NOT FOUND\n", def$label))
      } else {
        cat(sprintf("    %-12s %-6s (%s)  AIC %10.1f\n", def$label, r$family,
                    r$how, AIC(r$model)))
      }
    }
  }

  cat("\nData check:\n")
  dat <- tryCatch(load_marker_data(), error = function(e) {
    cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (!is.null(dat)) {
    for (tag in names(dat)) {
      d <- dat[[tag]]
      cat(sprintf("  %-5s %7d islets, %2d donors (ND %d / T1D %d)\n", tag,
                  nrow(d), dplyr::n_distinct(d$Donor),
                  dplyr::n_distinct(d$Donor[d$Diagnosis == "ND"]),
                  dplyr::n_distinct(d$Donor[d$Diagnosis == "T1D"])))
    }
  } else ok <- FALSE

  cat("\n")
  if (ok) {
    cat("ALL CHECKS PASSED — safe to run:\n")
    cat("  source(\"Marker/R/generate_marker_excel.R\");  generate_all_marker_excel()\n")
    cat("  source(\"Marker/R/plot_marker_main_figure.R\"); generate_marker_main_figure()\n")
    cat("  source(\"Marker/R/plot_marker_all_comparisons.R\"); generate_all_marker_plots()\n")
  } else {
    cat("SOMETHING IS MISSING — fix the items above before running the rest.\n")
  }
  cat("=============================================================================\n")
  invisible(ok)
}
