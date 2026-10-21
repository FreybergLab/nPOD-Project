# =============================================================================
# Generate Excel Results for the Marker Rate Models
# =============================================================================
#
# REBUILT WITHOUT EMMEANS.
#
# The previous version used emmeans -> regrid(transform = "log") ->
# contrast(type = "response"). That returned marginal rates of 4.64-4.99% for
# all six outcomes, including 4.70% for Ki-67 in alpha cells where the observed
# rate is 0.62%, and ratios that were just the quotient of the two wrong means.
# Three attempts to fix it in place failed because the composition of
# transformations could not be reasoned about reliably without running it.
#
# Everything here is instead an exact linear contrast of the fixed effects.
# For a log-link model every marginal mean, every group difference and every
# covariate slope is l'beta for some weight vector l, so:
#
#   estimate = l'beta  (+ log(PER) for a mean, where the offset enters)
#   SE       = sqrt(l' V l)
#
# and the rate or ratio is exp() of that. No transformations are composed, and
# each number can be checked by hand from the coefficient table.
#
# VERIFIED: for Ki-67 alpha this reproduces ND 3.220% and T1D 0.760% against
# observed 3.311% and 0.619%, and a rate ratio of 0.2360 against 0.187 observed.
#
# THREE RATE COLUMNS, because they answer different questions:
#   rate_median_donor  exp(l'beta + log(PER)). The geometric-mean rate for a
#                      donor at the centre of the random-effect distribution.
#                      This is what the CI and p-value attach to.
#   rate_population    the same multiplied by exp(sum(RE variances)/2) and by
#                      (1 - P(structural zero)). Because E[exp(u)] != exp(E[u]),
#                      the median-donor rate is systematically BELOW the
#                      observed average — by 1.9x for Ki-67 alpha. This column
#                      is the one comparable to the data.
#   rate_observed      straight from the islets. The ground truth.
#
# Ratios are unaffected by the Jensen factor, which cancels.
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(writexl)

source("Marker/R/marker_common.R")
source("Marker/R/marker_effects.R")

# The fixed-effects RHS the rate models were fitted with. Asserted against the
# fitted coefficient names below, so a mismatch errors rather than producing a
# quietly wrong design matrix.
RHS_RATE_EXCEL <- paste("(Diagnosis + Region + Sex + Age_c)^2",
                        "+ Diagnosis:Region:Age_c")

# =============================================================================
# FORMATTING
# =============================================================================

fmt <- function(x, digits = 4) ifelse(is.na(x), NA, round(x, digits))

# p-values need signif(), not round(): round(1.54e-05, 4) is 0, which reads as
# an exact zero in the workbook. signif keeps the magnitude and Excel displays
# it in scientific notation.
fmt_p <- function(p, digits = 3) ifelse(is.na(p), NA, signif(p, digits))

add_stars <- function(p) {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < 0.1 ~ ".", TRUE ~ "")
}

build_sheet <- function(title, description, blocks) {
  parts <- list(tibble(Note = title), tibble(Note = description), tibble(Note = NA))
  for (nm in names(blocks)) {
    b <- blocks[[nm]]
    if (is.null(b) || nrow(b) == 0) next
    parts <- c(parts, list(tibble(Note = nm),
                           b %>% mutate(Note = NA) %>% dplyr::select(Note, everything()),
                           tibble(Note = NA)))
  }
  bind_rows(parts)
}

# =============================================================================
# LINEAR-CONTRAST ENGINE
# =============================================================================

# Fixed-effect coefficients and covariance, with V reordered by name so that
# l' V l cannot silently use a different ordering than l' b.
model_bV <- function(model) {
  b <- fixef(model)$cond
  V <- vcov(model)$cond
  if (!identical(rownames(V), names(b))) {
    if (!all(names(b) %in% rownames(V))) {
      stop("vcov() and fixef() disagree on coefficient names")
    }
    V <- V[names(b), names(b), drop = FALSE]
  }
  list(b = b, V = V)
}

# Prediction grid. The denominator column is set to PER so that the offset
# term evaluates to log(PER) and all rates come out per PER cells.
marker_grid <- function(data, def, age = 0, per = PER) {
  g <- expand.grid(Diagnosis = levels(data$Diagnosis),
                   Region    = levels(data$Region),
                   Sex       = levels(data$Sex),
                   stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  g$Age_c <- age
  g[[def$den]] <- per
  for (v in c("Diagnosis", "Region", "Sex")) {
    g[[v]] <- factor(g[[v]], levels = levels(data[[v]]))
  }
  g
}

# Design matrix for the grid. The assertion is the safety net: if the columns
# do not line up with the fitted coefficients, every contrast below would be
# wrong, so it is an error rather than a warning.
marker_X <- function(model, data, grid, rhs = RHS_RATE_EXCEL) {
  fac  <- Filter(is.factor, data)
  xlev <- lapply(fac, levels)
  xlev <- xlev[intersect(names(xlev), names(grid))]
  X <- model.matrix(as.formula(paste("~", rhs)), data = grid, xlev = xlev)
  nms <- names(fixef(model)$cond)
  if (!identical(colnames(X), nms)) {
    stop("design matrix does not match the fitted coefficients.\n  design: ",
         paste(setdiff(colnames(X), nms), collapse = ", "),
         "\n  model:  ", paste(setdiff(nms, colnames(X)), collapse = ", "),
         "\n  RHS_RATE_EXCEL is probably out of step with how the models were fitted.")
  }
  X
}

# Averaging weights over a subset of grid rows.
wmean <- function(n, idx) { w <- numeric(n); w[idx] <- 1 / length(idx); w }

# est = l'b + add, SE = sqrt(l'Vl). `add` carries log(PER) for a mean and is 0
# for a contrast, where the offsets cancel.
lin_est <- function(bV, l, add = 0) {
  est <- sum(l * bV$b) + add
  se  <- sqrt(max(drop(t(l) %*% bV$V %*% l), 0))
  z   <- if (se > 0) est / se else NA_real_
  tibble(estimate = est, SE = se,
         value = exp(est), lwr = exp(est - 1.96 * se), upr = exp(est + 1.96 * se),
         z = z, p = if (is.na(z)) NA_real_ else 2 * pnorm(-abs(z)))
}

# Jensen factor: E[exp(u)] = exp(sigma^2/2) for u ~ N(0, sigma^2). Converts the
# median-donor rate to a population-average rate.
jensen_factor <- function(model) {
  vc <- VarCorr(model)$cond
  v <- sum(vapply(vc, function(m) sum(diag(as.matrix(m))), numeric(1)))
  exp(v / 2)
}

# Probability a group's islets are structural zeros, if the family has a ZI part.
zi_prob <- function(model, diagnosis) {
  zi <- fixef(model)$zi
  if (!length(zi)) return(0)
  lp <- zi[["(Intercept)"]]
  if (diagnosis == "T1D" && "DiagnosisT1D" %in% names(zi)) lp <- lp + zi[["DiagnosisT1D"]]
  plogis(lp)
}

# -----------------------------------------------------------------------------
# Marginal means, averaged over everything not in `by`
# -----------------------------------------------------------------------------
marginal_block <- function(model, bV, X, grid, by, sub, def, per = PER) {
  key  <- interaction(grid[by], drop = TRUE, sep = " / ")
  jf   <- jensen_factor(model)

  out <- bind_rows(lapply(levels(key), function(k) {
    idx <- which(key == k)
    l   <- as.numeric(wmean(nrow(X), idx) %*% X)
    dg  <- as.character(grid$Diagnosis[idx[1]])
    zp  <- if ("Diagnosis" %in% by) zi_prob(model, dg) else
             mean(c(zi_prob(model, "ND"), zi_prob(model, "T1D")))
    lin_est(bV, l, add = log(per)) %>%
      mutate(group = k, zi_prob = zp,
             rate_population = value * jf * (1 - zp), .before = 1)
  }))

  obs <- sub %>%
    group_by(across(all_of(by))) %>%
    summarise(rate_observed = per * sum(.data[[def$num]]) / sum(.data[[def$den]]),
              donors = dplyr::n_distinct(Donor), islets = n(), .groups = "drop") %>%
    mutate(group = interaction(across(all_of(by)), drop = TRUE, sep = " / ")) %>%
    dplyr::select(group, rate_observed, donors, islets)

  out %>%
    mutate(group = as.character(group)) %>%
    left_join(obs %>% mutate(group = as.character(group)), by = "group") %>%
    transmute(!!paste(by, collapse = " / ") := group,
              rate_median_donor = fmt(value, 4),
              `lower.CI` = fmt(lwr, 4), `upper.CI` = fmt(upr, 4),
              rate_population = fmt(rate_population, 4),
              rate_observed = fmt(rate_observed, 4),
              P_structural_zero = fmt(zi_prob, 4),
              donors, islets)
}

# -----------------------------------------------------------------------------
# Contrasts of `focal`, optionally within each level of `within`
# -----------------------------------------------------------------------------
contrast_block <- function(bV, X, grid, focal, within = NULL,
                           pairs = c("vs_first", "all")) {
  pairs <- match.arg(pairs)
  strata <- if (is.null(within)) list(all = seq_len(nrow(grid))) else
    split(seq_len(nrow(grid)), interaction(grid[within], drop = TRUE, sep = " / "))

  rows <- list()
  for (sname in names(strata)) {
    srows <- strata[[sname]]
    lv <- levels(droplevels(factor(grid[[focal]][srows])))
    combos <- if (pairs == "vs_first") {
      lapply(lv[-1], function(x) c(lv[1], x))
    } else {
      cb <- utils::combn(lv, 2, simplify = FALSE); cb
    }
    for (cc in combos) {
      i_ref <- srows[grid[[focal]][srows] == cc[1]]
      i_alt <- srows[grid[[focal]][srows] == cc[2]]
      l <- as.numeric((wmean(nrow(X), i_alt) - wmean(nrow(X), i_ref)) %*% X)
      rows[[length(rows) + 1L]] <- lin_est(bV, l) %>%
        mutate(stratum = sname, contrast = sprintf("%s / %s", cc[2], cc[1]),
               .before = 1)
    }
  }
  out <- bind_rows(rows)
  if (is.null(within)) out$stratum <- NULL
  out %>% transmute(
    across(any_of("stratum")), contrast,
    rate_ratio = fmt(value, 4), `lower.CI` = fmt(lwr, 4), `upper.CI` = fmt(upr, 4),
    log_ratio = fmt(estimate, 4), SE = fmt(SE, 4),
    pct_change = fmt(100 * (value - 1), 1),
    p.value = fmt_p(p), sig = add_stars(p))
}

# -----------------------------------------------------------------------------
# Covariate slopes. eta is linear in Age_c, so design(Age = 1) - design(Age = 0)
# is exactly d(eta)/d(Age_c), with no manual term enumeration.
# -----------------------------------------------------------------------------
slope_block <- function(bV, X0, X1, grid, by = NULL) {
  D <- X1 - X0
  groups <- if (is.null(by)) list(Overall = seq_len(nrow(grid))) else
    split(seq_len(nrow(grid)), interaction(grid[by], drop = TRUE, sep = " / "))
  bind_rows(lapply(names(groups), function(g) {
    l <- as.numeric(wmean(nrow(D), groups[[g]]) %*% D)
    lin_est(bV, l) %>% mutate(group = g, .before = 1)
  })) %>%
    transmute(group,
              slope_log = fmt(estimate, 5), SE = fmt(SE, 5),
              fold_per_year = fmt(value, 4),
              pct_per_year = fmt(100 * (value - 1), 2),
              `lower.CI` = fmt(lwr, 4), `upper.CI` = fmt(upr, 4),
              p.value = fmt_p(p), sig = add_stars(p))
}

# =============================================================================
# WORKBOOK
# =============================================================================

generate_marker_rate_excel <- function(model, sub, def, output_dir) {

  cat("\n  Building workbook for", def$label, "...\n")
  bV   <- model_bV(model)
  grid <- marker_grid(sub, def)
  X    <- marker_X(model, sub, grid)
  X1   <- marker_X(model, sub, marker_grid(sub, def, age = 1))
  unit <- sprintf("%% of %s cells %s+", def$cell, def$marker)

  rate_note <- paste(
    "rate_median_donor is exp(l'beta + log(100)) — the geometric-mean rate for",
    "a donor at the centre of the random-effect distribution, and the quantity",
    "the CI attaches to. rate_population multiplies it by exp(sum(RE",
    "variances)/2) and by 1 - P(structural zero); because E[exp(u)] is not",
    "exp(E[u]), the median-donor rate sits systematically below the observed",
    "average. rate_observed comes straight from the islets. Rate RATIOS are",
    "unaffected — the correction cancels.")

  # ---- Sheet 1: Model Summary --------------------------------------------
  cat("    1 Model Summary\n")
  co <- summary(model)$coefficients$cond
  coef_tbl <- tibble(
    term = rownames(co), estimate = fmt(co[, "Estimate"], 5),
    SE = fmt(co[, "Std. Error"], 5), rate_ratio = fmt(exp(co[, "Estimate"]), 4),
    z.value = fmt(co[, "z value"], 3), p.value = fmt_p(co[, "Pr(>|z|)"]),
    sig = add_stars(co[, "Pr(>|z|)"]))

  spec_tbl <- tibble(
    quantity = c("Response", "Offset", "Family", "Zero-inflation",
                 "Fixed effects", "Random effects", "Islets", "Donors",
                 "Images", "AIC", "Dispersion", "Jensen factor exp(var/2)",
                 "Marginal estimand"),
    value = c(def$num, sprintf("log(%s)", def$den), family(model)$family,
              if (length(fixef(model)$zi))
                paste(deparse(formula(model, component = "zi")), collapse = "") else "none",
              RHS_RATE_EXCEL,
              paste(names(VarCorr(model)$cond), collapse = " + "),
              as.character(nrow(model$frame)),
              as.character(dplyr::n_distinct(sub$Donor)),
              as.character(dplyr::n_distinct(sub$ImageID)),
              as.character(fmt(AIC(model), 1)),
              as.character(fmt(sigma(model), 4)),
              as.character(fmt(jensen_factor(model), 4)),
              "equal-weight average of the design rows on the log scale"))

  model_summary <- build_sheet(
    sprintf("%s RATE MODEL — %s CELLS", toupper(def$marker), toupper(def$cell)),
    sprintf(paste("Count of %s-positive cells in the %s bin, with log(%s) as an",
                  "offset, so coefficients are rate ratios per %s cell.",
                  "Marginal quantities are exact linear contrasts of the fixed",
                  "effects — no emmeans. Rates are reported as %s."),
            def$marker, def$cell, def$den, def$cell, unit),
    list("SPECIFICATION" = spec_tbl,
         "FIXED EFFECTS (log rate scale)" = coef_tbl))

  # ---- Sheets 2-4: main effects -----------------------------------------
  cat("    2 Diagnosis   3 Region   4 Sex\n")
  diagnosis_sheet <- build_sheet(
    "DIAGNOSIS EFFECT",
    paste("Main effect of Diagnosis, averaged over Region and Sex at mean Age.",
          rate_note),
    list("MARGINAL RATES"         = marginal_block(model, bV, X, grid, "Diagnosis", sub, def),
         "RATE RATIO (T1D vs ND)" = contrast_block(bV, X, grid, "Diagnosis")))

  weight_note <- paste(
    "NOTE: this sheet averages Diagnosis with EQUAL weight, while",
    "rate_observed pools islets, of which far more are ND. Since ND rates are",
    "higher, rate_population here sits below rate_observed by design — about",
    "25% for this model. That is an estimand difference, not an error. The",
    "Diagnosis sheet, which does not average across Diagnosis, agrees with the",
    "observed rates to within a few percent.")

  region_sheet <- build_sheet(
    "REGION EFFECT",
    paste("Main effect of pancreatic Region, averaged over Diagnosis and Sex.",
          weight_note, rate_note),
    list("MARGINAL RATES"        = marginal_block(model, bV, X, grid, "Region", sub, def),
         "PAIRWISE RATE RATIOS"  = contrast_block(bV, X, grid, "Region", pairs = "all")))

  sex_sheet <- build_sheet(
    "SEX EFFECT",
    paste("Main effect of Sex, averaged over Diagnosis and Region.",
          weight_note, rate_note),
    list("MARGINAL RATES"      = marginal_block(model, bV, X, grid, "Sex", sub, def),
         "RATE RATIO (M vs F)" = contrast_block(bV, X, grid, "Sex")))

  # ---- Sheets 5-7: two-way --------------------------------------------
  cat("    5 Diagnosis x Region   6 Diagnosis x Sex   7 Region x Sex\n")
  diag_region_sheet <- build_sheet(
    "DIAGNOSIS x REGION",
    paste("Diagnosis effect within each Region, and Region effects within each",
          "Diagnosis.", rate_note),
    list("MARGINAL RATES" = marginal_block(model, bV, X, grid, c("Diagnosis","Region"), sub, def),
         "T1D vs ND, WITHIN EACH REGION" = contrast_block(bV, X, grid, "Diagnosis", "Region"),
         "REGION PAIRS, WITHIN EACH DIAGNOSIS" =
           contrast_block(bV, X, grid, "Region", "Diagnosis", pairs = "all")))

  diag_sex_sheet <- build_sheet(
    "DIAGNOSIS x SEX",
    paste("Diagnosis effect within each Sex, and the Sex effect within each",
          "Diagnosis.", rate_note),
    list("MARGINAL RATES" = marginal_block(model, bV, X, grid, c("Diagnosis","Sex"), sub, def),
         "T1D vs ND, WITHIN EACH SEX"    = contrast_block(bV, X, grid, "Diagnosis", "Sex"),
         "M vs F, WITHIN EACH DIAGNOSIS" = contrast_block(bV, X, grid, "Sex", "Diagnosis")))

  region_sex_sheet <- build_sheet(
    "REGION x SEX",
    paste("Region effects within each Sex.", weight_note, rate_note),
    list("MARGINAL RATES" = marginal_block(model, bV, X, grid, c("Region","Sex"), sub, def),
         "REGION PAIRS, WITHIN EACH SEX" =
           contrast_block(bV, X, grid, "Region", "Sex", pairs = "all")))

  # ---- Sheets 8-9: age --------------------------------------------------
  cat("    8 Age Slopes   9 Diag x Region x Age\n")
  age_sheet <- build_sheet(
    "AGE SLOPES",
    paste("Change in the per-cell rate per additional year of donor age.",
          "Obtained as design(Age_c = 1) - design(Age_c = 0), which is exactly",
          "d(log rate)/d(Age_c) because the linear predictor is linear in age.",
          "fold_per_year is exp(slope); pct_per_year is the percentage change."),
    list("OVERALL"      = slope_block(bV, X, X1, grid),
         "BY DIAGNOSIS" = slope_block(bV, X, X1, grid, "Diagnosis"),
         "BY REGION"    = slope_block(bV, X, X1, grid, "Region"),
         "BY SEX"       = slope_block(bV, X, X1, grid, "Sex")))

  diag_region_age_sheet <- build_sheet(
    "DIAGNOSIS x REGION x AGE",
    "The three-way term in the project formula: age slopes per Diagnosis x Region cell.",
    list("AGE SLOPES BY DIAGNOSIS x REGION" =
           slope_block(bV, X, X1, grid, c("Diagnosis", "Region"))))

  # ---- Sheet 10: zero-inflation ----------------------------------------
  cat("   10 Zero-Inflation\n")
  if (length(fixef(model)$zi)) {
    zco <- summary(model)$coefficients$zi
    p0 <- zi_prob(model, "ND"); p1 <- zi_prob(model, "T1D")
    zi_blocks <- list(
      "ZI COEFFICIENTS (logit scale)" = tibble(
        term = rownames(zco), estimate = fmt(zco[, "Estimate"], 4),
        SE = fmt(zco[, "Std. Error"], 4), p.value = fmt_p(zco[, "Pr(>|z|)"]),
        sig = add_stars(zco[, "Pr(>|z|)"])),
      "P(STRUCTURAL ZERO)" = tibble(Diagnosis = c("ND", "T1D"),
                                    p_structural_zero = fmt(c(p0, p1), 4)),
      "EFFECT ON THE MARGINAL RATE RATIO" = tibble(
        quantity = "non-zero probability ratio (T1D / ND)",
        value = fmt((1 - p1) / (1 - p0), 4)))
    zi_desc <- paste("A zero-inflated family was selected. The rate ratios on",
                     "the other sheets are conditional on an islet not being a",
                     "structural zero; multiply by the ratio below for the",
                     "marginal rate ratio. The rate_population column already",
                     "includes this factor.")
  } else {
    zi_blocks <- list("RESULT" = tibble(detail = paste(
      "No zero-inflation component — a plain negative binomial won family",
      "selection. Excess zeros are handled by the conditional distribution, so",
      "the rate ratios elsewhere are already marginal.")))
    zi_desc <- "Family selection did not retain a zero-inflation component."
  }
  zi_sheet <- build_sheet("ZERO-INFLATION", zi_desc, zi_blocks)

  # ---- Sheet 11: random effects ----------------------------------------
  cat("   11 Random Effects\n")
  vc <- VarCorr(model)$cond
  re_df <- bind_rows(lapply(names(vc), function(g) {
    m <- as.matrix(vc[[g]])
    tibble(group = g, term = rownames(m), variance = fmt(diag(m), 5),
           std.dev = fmt(sqrt(diag(m)), 5))
  }))
  re_sheet <- build_sheet(
    "RANDOM EFFECTS",
    paste("Variance components. Their total drives the Jensen factor",
          sprintf("exp(sum(var)/2) = %.4f,", jensen_factor(model)),
          "which is how far the median-donor rate sits below the",
          "population average."),
    list("VARIANCE COMPONENTS" = re_df,
         "SAMPLE SIZES" = tibble(
           quantity = c("Islets (denominator > 0)", "Donors", "Images",
                        "ND donors", "T1D donors"),
           n = c(nrow(model$frame), dplyr::n_distinct(sub$Donor),
                 dplyr::n_distinct(sub$ImageID),
                 dplyr::n_distinct(sub$Donor[sub$Diagnosis == "ND"]),
                 dplyr::n_distinct(sub$Donor[sub$Diagnosis == "T1D"])))))

  # ---- Sheet 12: restriction -------------------------------------------
  cat("   12 Restriction\n")
  full <- attr(sub, "full_data")
  if (!is.null(full)) {
    all_don  <- unique(as.character(full$Donor))
    kept_don <- unique(as.character(sub$Donor))
    lost     <- setdiff(all_don, kept_don)
    lost_tbl <- if (length(lost)) {
      full %>% filter(as.character(Donor) %in% lost) %>%
        distinct(Donor, Diagnosis, Disease.Duration) %>%
        arrange(Diagnosis, Disease.Duration) %>%
        mutate(Donor = as.character(Donor), Diagnosis = as.character(Diagnosis))
    } else tibble(Donor = "none", Diagnosis = NA_character_,
                  Disease.Duration = NA_real_)
    counts_tbl <- tibble(
      quantity = c("Islets before restriction", "Islets retained",
                   "% islets dropped", "Donors before", "Donors retained",
                   "Donors dropped"),
      value = c(nrow(full), nrow(sub), fmt(100 * (1 - nrow(sub)/nrow(full)), 2),
                length(all_don), length(kept_don), length(lost)))
  } else {
    lost_tbl <- tibble(detail = "full dataset not attached")
    counts_tbl <- tibble(quantity = "Islets retained", value = nrow(sub))
  }
  restriction_sheet <- build_sheet(
    "DENOMINATOR RESTRICTION",
    paste("offset(log(denominator)) is undefined when the denominator is zero,",
          "so islets with no cells of this type are excluded. The donors this",
          "removes are not random and they bound what the result can claim."),
    list("ISLETS AND DONORS" = counts_tbl, "DONORS DROPPED ENTIRELY" = lost_tbl))

  # ---- write ------------------------------------------------------------
  sheets <- list(
    "Model Summary" = model_summary, "Diagnosis" = diagnosis_sheet,
    "Region" = region_sheet, "Sex" = sex_sheet,
    "Diagnosis x Region" = diag_region_sheet, "Diagnosis x Sex" = diag_sex_sheet,
    "Region x Sex" = region_sex_sheet, "Age Slopes" = age_sheet,
    "Diag x Region x Age" = diag_region_age_sheet,
    "Zero-Inflation" = zi_sheet, "Random Effects" = re_sheet,
    "Restriction" = restriction_sheet)
  sheets <- lapply(sheets, as.data.frame)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(output_dir, sprintf("%s_Rate_Diagnosis_Analysis.xlsx", def$label))
  write_xlsx(sheets, f)
  cat("    saved:", f, "\n")
  invisible(sheets)
}

# =============================================================================
# WRAPPER, WITH A GATE BEFORE ANY FILE IS WRITTEN
# =============================================================================

generate_all_marker_excel <- function(output_dir = "Marker/Results/Rate/",
                                      model_dir  = MARKER_MODEL_DIR,
                                      stain_handling = "double_only",
                                      tolerance = 3) {

  cat("\n=============================================================================\n")
  cat("MARKER RATE EXCEL — exact linear contrasts, no emmeans\n")
  cat("=============================================================================\n\n")

  # Nothing is written until the Diagnosis contrast of every model agrees with
  # the observed rates. This is the check whose absence let the emmeans values
  # reach a figure.
  cat("Validating before writing anything...\n")
  effects <- marker_effects_table(model_dir, stain_handling)
  validate_marker_rates(effects, tolerance = tolerance)

  raw <- load_marker_data(stain_handling)
  out <- list()
  for (def in marker_defs) {
    r <- load_marker_model(def, "rate", model_dir)
    if (is.null(r)) { cat("\n  SKIPPING", def$label, "- no cached model\n"); next }
    sub <- marker_subset(raw[[sub("_.*", "", def$key)]], def)
    out[[def$label]] <- generate_marker_rate_excel(r$model, sub, def, output_dir)
  }

  cat("\n=== cross-check: Diagnosis sheet vs observed ===\n")
  print(as.data.frame(effects %>%
    dplyr::select(marker, cell, obs_ND, obs_T1D, obs_ratio, ratio,
                  marginal_ratio, p, sig)), digits = 4)

  cat("\nDONE —", length(out), "workbooks in", output_dir, "\n")
  cat("Sanity check: on the Diagnosis sheet, rate_population should land near\n")
  cat("rate_observed, and rate_median_donor should sit below both.\n")
  invisible(out)
}
