# sensitivity_composition_HTC.R
#
# Ordered-beta sensitivity refits for the islet-composition models. Packages,
# sampler settings and formulas are those of proportion_diagnosis_HTC.R (the
# *_powered fits) and proportion_t1d_HTC.R (the *_ao fits); each variant
# changes exactly one thing.
#
# Before the fits, once (installs any missing packages, checks the data, fits
# nothing):
#   Rscript sensitivity_composition_HTC.R check
# One fit per call, so variants can go out as a job array:
#   Rscript sensitivity_composition_HTC.R <outcome> <variant> [model]
#     outcome: ins | glu | soma | pp
#     variant: slopes | no_atypical | no_hpap | age_ge5 | no_cap | min10 | min20
#     model:   diag (default; compare with Models/<outcome>_powered_model.rds)
#              ao   (T1D only; compare with Models/<outcome>_ao_model.rds)
# Once fits exist:
#   Rscript sensitivity_composition_HTC.R summarize
#     -> Sensitivity/composition_sensitivity.csv
#
# Needs Data/no_cutoff_data.rds (the size variants start before the size
# filter) and Data/quad_data.rds. Every fit first rebuilds the published quad
# set and stops if its row count differs from Data/quad_data.rds; summarize
# checks that the published fits reproduce two Table S2 values.

# ============================================================
# PACKAGES
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install into the personal library. R uses it only if the directory exists,
# and a batch job cannot offer to create it, so install.packages() would
# otherwise fall back to the module's read-only system library.
user_lib <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep)[[1]][1]
if (!is.na(user_lib) && nzchar(user_lib)) {
  user_lib <- path.expand(user_lib)
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_lib, .libPaths()))
  cat("installing any missing packages into", user_lib, "\n")
}

for (pkg in c("tidyverse", "brms", "ordbetareg")) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) install.packages(pkg)
}
if (!require("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}
if (!dir.exists(cmdstanr::cmdstan_path())) cmdstanr::install_cmdstan()

library(ordbetareg)
library(brms)
library(cmdstanr)
library(tidyverse)

# ============================================================
# SPECIFICATION
# ============================================================

OUTCOMES <- c(ins = "Percent.ins", glu = "Percent.glu", soma = "Percent.soma", pp = "Percent.PP")
VARIANTS <- c("slopes", "no_atypical", "no_hpap", "age_ge5", "no_cap", "min10", "min20")
ATYPICAL <- c("6405", "6456", "6466")   # pregnant at death; 6 days postpartum; autopsy, severe autolysis

# Verbatim from the published HTC scripts.
RHS <- c(
  diag = paste("(Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +",
               "Diagnosis:Region:log_Islet.Cells_c + Diagnosis:Region:Age_c"),
  ao   = paste("(Disease.Duration_c + Region + Sex + age_at_onset_c + log_Islet.Cells_c)^2 +",
               "Disease.Duration_c:Region:log_Islet.Cells_c + Disease.Duration_c:Region:age_at_onset_c")
)
RE <- c(base   = "(1|Donor) + (1|Donor:ImageID)",
        slopes = "(1 + log_Islet.Cells_c|Donor) + (1|Donor:ImageID)")
PUBLISHED <- c(diag = "Models/%s_powered_model.rds", ao = "Models/%s_ao_model.rds")

# Table S2 values the summary must reproduce from the published fits: every
# quantity it reports, from the workbooks Table S2 is built from
# (Results/<cell> Proportions/*_Content_Primary_Analysis.xlsx for the diagnosis
# models, *_Content_AO_Analysis.xlsx for the T1D models). T1D - ND is the
# negative of the workbooks' ND - T1D.
S2_CHECK <- tribble(
  ~outcome, ~model, ~quantity,         ~level,     ~s2,
  "ins",    "diag", "T1D - ND",         "T1D - ND", -0.5814,
  "ins",    "diag", "islet_size_slope", "ND",       +0.0037,
  "ins",    "diag", "islet_size_slope", "T1D",      +0.0081,
  "ins",    "ao",   "duration_slope",   "overall",  -0.0174,
  "ins",    "ao",   "islet_size_slope", "overall",  +0.0104,
  "glu",    "diag", "T1D - ND",         "T1D - ND", +0.4613,
  "glu",    "diag", "islet_size_slope", "ND",       +0.0542,
  "glu",    "diag", "islet_size_slope", "T1D",      +0.0312,
  "glu",    "ao",   "duration_slope",   "overall",  +0.0261,
  "glu",    "ao",   "islet_size_slope", "overall",  +0.0294,
  "soma",   "diag", "T1D - ND",         "T1D - ND", +0.0037,
  "soma",   "diag", "islet_size_slope", "ND",       -0.0375,
  "soma",   "diag", "islet_size_slope", "T1D",      -0.0155,
  "soma",   "ao",   "duration_slope",   "overall",  -0.0061,
  "soma",   "ao",   "islet_size_slope", "overall",  -0.0118,
  "pp",     "diag", "T1D - ND",         "T1D - ND", +0.0051,
  "pp",     "diag", "islet_size_slope", "ND",       -0.0154,
  "pp",     "diag", "islet_size_slope", "T1D",      -0.0228,
  "pp",     "ao",   "duration_slope",   "overall",  +0.0034,
  "pp",     "ao",   "islet_size_slope", "overall",  -0.0203 
)
# The published fits are the same draws summarised the same way, so they must
# match to the precision Table S2 reports: half a unit in the 4th decimal
# (0.00005), plus a little for floating point.
S2_TOL <- 6e-5

# ============================================================
# DATA AND FIT
# ============================================================

build_data <- function(variant, model) {
  quad_all <- readRDS("Data/no_cutoff_data.rds") %>%
    filter(Stain == "Quad") %>%
    mutate(Donor = as.character(Donor))

  n_rebuilt <- sum(quad_all$Islet.Cells >= 15 & quad_all$log_Islet.Cells <= 5.5)
  n_pub     <- nrow(readRDS("Data/quad_data.rds"))
  if (n_rebuilt != n_pub) {
    stop(sprintf("Rebuilt quad set has %d rows; Data/quad_data.rds has %d.", n_rebuilt, n_pub))
  }

  min_cells <- switch(variant, min10 = 10, min20 = 20, 15)
  max_log   <- if (variant == "no_cap") Inf else 5.5
  d <- quad_all %>% filter(Islet.Cells >= min_cells, log_Islet.Cells <= max_log)

  d <- switch(variant,
              no_atypical = d %>% filter(!Donor %in% ATYPICAL),
              no_hpap     = d %>% filter(!str_detect(Donor, "^HPAP")),
              age_ge5     = d %>% filter(Age >= 5),
              d)
  if (model == "ao") d <- d %>% filter(Diagnosis == "T1D")

  # Re-centred on the analysis set, as in Create_count_datasets_v1.Rmd.
  d %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
           Region    = factor(Region, levels = c("Head", "Body", "Tail")),
           Sex       = factor(Sex, levels = c("Female", "Male")),
           Donor     = factor(Donor),
           ImageID   = factor(ImageID),
           Age_c              = scale(Age, center = TRUE, scale = FALSE)[, 1],
           log_Islet.Cells_c  = scale(log_Islet.Cells, center = TRUE, scale = FALSE)[, 1],
           Disease.Duration_c = scale(Disease.Duration, center = TRUE, scale = FALSE)[, 1],
           age_at_onset_c     = scale(Age - Disease.Duration, center = TRUE, scale = FALSE)[, 1]) %>%
    droplevels()
}

fit_one <- function(outcome, variant, model) {
  d <- build_data(variant, model)
  f <- as.formula(paste(OUTCOMES[[outcome]], "~", RHS[[model]], "+",
                        RE[[if (variant == "slopes") "slopes" else "base"]]))
  cat(sprintf("%s | %s | %s: %d islets, %d donors\n",
              outcome, model, variant, nrow(d), n_distinct(d$Donor)))
  ordbetareg(f, data = d,
             cores = 4, chains = 4, threads = threading(4),
             iter = 2000, warmup = 1000,
             backend = "cmdstanr", seed = 123,
             file = file.path("Models_sens", paste(outcome, model, variant, sep = "_")))
}

# ============================================================
# SUMMARY
# ============================================================

# Posterior median, 95% HPD and probability of direction for the quantities
# the paper reports, plus the sampler diagnostics it does not yet report.
#
# Same method as the published tables (Proportion/R/generate_diag_excel_results.R
# and generate_proportion_t1d_excel.R), so a published fit reproduces Table S2:
# expected proportions without random effects, other continuous predictors at
# the published data's averages, in every Region x Sex cell; a slope is
# (fit at +1 - fit at -1) / 2; cells are averaged with equal weight.
#
# Grid values are offsets from those averages. `shift` moves them into this
# fit's own centred units just before predicting (see shift_for() below).
summarise_fit <- function(fit, model, shift = NULL) {
  cells  <- expand_grid(Region = c("Head", "Body", "Tail"), Sex = c("Female", "Male"))
  fitted_draws <- function(grid) {
    for (v in names(shift)) if (v %in% names(grid)) grid[[v]] <- grid[[v]] + shift[[v]]
    fitted(fit, newdata = grid, re_formula = NA, summary = FALSE)
  }
  # Slope over var: the grid's cells at var = -1 and +1, in the same cell order.
  slope_draws <- function(base, var, keep = rep(TRUE, nrow(base))) {
    g <- bind_rows(mutate(base, !!var := -1), mutate(base, !!var := 1))
    e <- fitted_draws(g)
    lo <- which(g[[var]] == -1 & rep(keep, 2))
    hi <- which(g[[var]] ==  1 & rep(keep, 2))
    rowMeans((e[, hi, drop = FALSE] - e[, lo, drop = FALSE]) / 2)
  }
  if (model == "diag") {
    base <- expand_grid(Diagnosis = c("ND", "T1D"), cells, Age_c = 0, log_Islet.Cells_c = 0)
    e    <- fitted_draws(base)
    quantities <- list(
      list("T1D - ND", "T1D - ND",
           rowMeans(e[, base$Diagnosis == "T1D", drop = FALSE]) - rowMeans(e[, base$Diagnosis == "ND", drop = FALSE])),
      list("islet_size_slope", "ND",  slope_draws(base, "log_Islet.Cells_c", base$Diagnosis == "ND")),
      list("islet_size_slope", "T1D", slope_draws(base, "log_Islet.Cells_c", base$Diagnosis == "T1D")))
  } else {
    base <- expand_grid(cells, Disease.Duration_c = 0, age_at_onset_c = 0, log_Islet.Cells_c = 0)
    quantities <- list(
      list("duration_slope",   "overall", slope_draws(base, "Disease.Duration_c")),
      list("islet_size_slope", "overall", slope_draws(base, "log_Islet.Cells_c")))
  }
  np <- nuts_params(fit)
  map_dfr(quantities, function(q) {
    hpd <- coda::HPDinterval(coda::as.mcmc(q[[3]]), prob = 0.95)
    tibble(quantity  = q[[1]],
           level     = q[[2]],
           estimate  = median(q[[3]]),
           lower.HPD = hpd[1, "lower"],
           upper.HPD = hpd[1, "upper"],
           pd        = max(mean(q[[3]] > 0), mean(q[[3]] < 0)))
  }) %>%
    mutate(n_obs         = nobs(fit),
           max_rhat      = max(rhat(fit), na.rm = TRUE),
           min_ess_ratio = min(neff_ratio(fit), na.rm = TRUE),
           divergences   = sum(np$Value[np$Parameter == "divergent__"]))
}

# ============================================================
# MAIN
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

check_data_files <- function() {
  missing_data <- Filter(Negate(file.exists), c("Data/no_cutoff_data.rds", "Data/quad_data.rds"))
  if (length(missing_data) > 0) {
    stop("not found under ", getwd(), ": ", paste(missing_data, collapse = ", "), call. = FALSE)
  }
}

if (length(args) >= 1 && args[1] == "check") {
  # Everything a fit needs except the fit: packages (loaded above), data files,
  # and the rebuild of the published quad set with its row-count check.
  check_data_files()
  d <- build_data("slopes", "diag")
  cat(sprintf("ready: packages loaded, data found, quad set rebuilt (%d islets, %d donors)\n",
              nrow(d), n_distinct(d$Donor)))
} else if (length(args) >= 1 && args[1] == "summarize") {
  dir.create("Sensitivity", showWarnings = FALSE)

  sens <- tibble(path = Sys.glob("Models_sens/*.rds")) %>%
    mutate(key = str_remove(basename(path), "\\.rds$")) %>%
    separate(key, into = c("outcome", "model", "variant"), sep = "_", extra = "merge")
  published <- expand_grid(outcome = names(OUTCOMES), model = names(PUBLISHED)) %>%
    mutate(path = sprintf(PUBLISHED[model], outcome), variant = "published")
  if (!all(file.exists(published$path))) {
    warning("published fits not found, so not summarized: ",
            paste(published$path[!file.exists(published$path)], collapse = ", "))
  }
  published <- published %>% filter(file.exists(path))
  if (nrow(published) + nrow(sens) == 0) {
    stop("no fits found in Models/ or Models_sens/ under ", getwd(),
         ". Run Rscript from the folder that holds them.")
  }

  # Every fit is measured at the raw values where Table S2 measures the
  # published fits: the diagnosis models at their own 0, and the T1D models at
  # the T1D averages. So a sensitivity row differs from the published row only
  # through its data, not through a moved reference point. The published T1D
  # models were centred on ND + T1D data, so their own 0 is not the T1D average
  # (for duration it is ~1.4 years); each refit is centred on its own subset.
  raw_of <- c(Age_c = "Age", log_Islet.Cells_c = "log_Islet.Cells",
              Disease.Duration_c = "Disease.Duration", age_at_onset_c = "age_at_onset")
  with_onset <- function(d) mutate(d, age_at_onset = Age - Disease.Duration)
  # The raw value where each of a fit's centred predictors is 0, from the fit
  # itself: mean(raw) - mean(centred) over the rows it was fitted on. The
  # centred columns in Data/quad_data.rds are not relied on, because they need
  # not be the ones a published fit used.
  zero_of <- function(fit, d) {
    shared <- intersect(names(fit$data), names(d))
    d <- d[complete.cases(d[shared]), , drop = FALSE]   # rows brms would keep
    if (nrow(d) != nrow(fit$data)) {
      stop("a fit has ", nrow(fit$data), " rows but its data set has ", nrow(d))
    }
    d <- with_onset(d)
    sapply(names(raw_of), function(v) {
      if (!v %in% names(fit$data)) return(NA_real_)
      mean(d[[raw_of[[v]]]], na.rm = TRUE) - mean(fit$data[[v]], na.rm = TRUE)
    })
  }
  pub_sets  <- list(diag = readRDS("Data/quad_data.rds"))
  pub_sets$ao <- filter(pub_sets$diag, Diagnosis == "T1D")
  t1d_average <- with_onset(pub_sets$ao) %>%
    { setNames(sapply(raw_of, function(r) mean(.[[r]], na.rm = TRUE)), names(raw_of)) }
  # Reference points, taken from each published fit (summarized before its
  # refits): a diagnosis model's own 0, or the T1D averages for a T1D model.
  reference <- list()
  shift_for <- function(fit, outcome, model, variant) {
    d    <- if (variant == "published") pub_sets[[model]] else build_data(variant, model)
    zero <- zero_of(fit, d)
    key  <- paste(outcome, model)
    if (variant == "published") reference[[key]] <<- if (model == "diag") zero else t1d_average
    if (is.null(reference[[key]])) stop("no published ", key, " fit to take the reference point from")
    s <- reference[[key]] - zero
    s[is.na(s)] <- 0   # predictors this model does not use
    s
  }

  out <- pmap_dfr(bind_rows(published, sens), function(path, outcome, model, variant) {
    cat("summarizing", path, "\n")
    fit <- readRDS(path)
    res <- summarise_fit(fit, model, shift_for(fit, outcome, model, variant)) %>%
      mutate(outcome = outcome, model = model, variant = variant, .before = 1)
    rm(fit); gc()
    res
  }) %>%
    left_join(S2_CHECK, by = c("outcome", "model", "quantity", "level")) %>%
    mutate(reproduces_S2 = if_else(variant == "published" & !is.na(s2),
                                   abs(estimate - s2) <= S2_TOL, NA))

  if (any(out$reproduces_S2 %in% FALSE)) {
    warning("A published fit does not reproduce Table S2; the summary method differs from yours.")
  }
  write_csv(out, "Sensitivity/composition_sensitivity.csv")
  print(out %>%
          dplyr::select(outcome, model, variant, quantity, level, estimate,
                        lower.HPD, upper.HPD, pd, max_rhat, divergences, reproduces_S2),
        n = Inf, width = Inf)
} else {
  model <- if (length(args) >= 3) args[3] else "diag"
  if (length(args) < 2 || !args[1] %in% names(OUTCOMES) ||
      !args[2] %in% VARIANTS || !model %in% names(RHS)) {
    stop("bad arguments: '", paste(args, collapse = " "), "'\n",
         "Usage: Rscript sensitivity_composition_HTC.R <outcome> <variant> [model]\n",
         "  outcome: ", paste(names(OUTCOMES), collapse = " | "), "\n",
         "  variant: ", paste(VARIANTS, collapse = " | "), "\n",
         "  model:   ", paste(names(RHS), collapse = " | "), " (default diag)\n",
         "  or:      Rscript sensitivity_composition_HTC.R check | summarize",
         call. = FALSE)
  }
  check_data_files()
  dir.create("Models_sens", showWarnings = FALSE)
  fit_one(args[1], args[2], model)
  cat("done\n")
}
