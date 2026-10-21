# =============================================================================
# marker_effects.R
# Marginal Diagnosis effects for the rate models, plus a validation gate.
#
# WHY THIS EXISTS
# The emmeans + regrid + contrast chain in generate_marker_excel.R returned
# marginal "rates" of 4.65-5.07% for all six outcomes, including 4.70% for
# Ki-67 in alpha cells where the observed rate is 0.62%, and rate ratios that
# were simply the quotient of those two wrong means (0.94 where the truth is
# about 0.24, 1.09 where the truth is about 9.2). That path is not used here.
#
# Everything below is computed from fixef() and vcov() as an exact linear
# contrast, and every model estimate is printed beside the observed rate so a
# scale error cannot pass unnoticed. validate_marker_rates() turns that
# comparison into a hard stop.
# =============================================================================

library(tidyverse)
library(glmmTMB)

# Coefficients are matched on the SET of their components, so the order R uses
# for interaction names does not have to be guessed.
find_term <- function(nms, parts) {
  hit <- nms[vapply(strsplit(nms, ":", fixed = TRUE),
                    function(x) setequal(x, parts), logical(1))]
  if (length(hit) == 1) hit else NA_character_
}

# Marginal T1D-vs-ND log rate ratio, averaged with equal weights over the three
# regions and two sexes.
#
# NOTE ON THE ESTIMAND. The average is taken on the LOG scale, so the result is
# the geometric mean of the six region-by-sex rate ratios. That is the standard
# convention for a log-link model and matches what emmeans does by default. It
# is NOT the same quantity as the Proportion/Area workbooks, which average
# posterior fitted values on the response scale (arithmetic mean). Both are
# legitimate definitions of "averaged over region and sex"; they differ because
# a mean of logs is not the log of a mean. Worth one sentence in methods so the
# two domains are not assumed identical.
#
# With Age_c = 0 and the formula
#   (Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c
# the log rate difference in region r, sex s is
#   b[DiagnosisT1D] + b[DiagnosisT1D:r] + b[DiagnosisT1D:s]
# (Age terms vanish at Age_c = 0). Averaging over r and s gives
#   b[DiagnosisT1D] + (1/3)(b[:RegionBody] + b[:RegionTail]) + (1/2) b[:SexMale]
# with Head and Female as references. SE uses the same weights against vcov(),
# so the interval and p-value are exact rather than approximated.
# strict = TRUE: every term the equal-weight average needs must be present, and
#   a missing one is an error — for the primary models a missing term would mean
#   the formula was not what we think it was.
# strict = FALSE: missing interaction terms are allowed and simply carry zero
#   weight, which is the CORRECT average for a model that does not contain them
#   (the Diagnosis effect is then constant across region and sex). Used by the
#   sensitivity analyses, where a restricted subset may force a reduced formula.
#   n_terms in the output records how many were found, so the reader can see it.
marginal_diagnosis <- function(model, verbose = FALSE, strict = TRUE) {
  b <- fixef(model)$cond
  V <- vcov(model)$cond
  nms <- names(b)

  # t(w) %*% V %*% w ignores names, so V must be in the same order as b. Reorder
  # by name rather than trusting that it already is; a silent mismatch would
  # give a plausible-looking but wrong standard error.
  if (!identical(rownames(V), nms)) {
    if (!all(nms %in% rownames(V))) {
      stop("vcov() and fixef() disagree on coefficient names — cannot build ",
           "the contrast safely")
    }
    V <- V[nms, nms, drop = FALSE]
  }

  main <- find_term(nms, "DiagnosisT1D")
  if (is.na(main)) stop("no DiagnosisT1D coefficient in this model")

  # Every term the average needs. A missing one is NOT skipped quietly: leaving
  # it out of w changes the estimand, so it has to be an error.
  wanted <- list(list(parts = "DiagnosisT1D",                      wt = 1),
                 list(parts = c("DiagnosisT1D", "RegionBody"),     wt = 1/3),
                 list(parts = c("DiagnosisT1D", "RegionTail"),     wt = 1/3),
                 list(parts = c("DiagnosisT1D", "SexMale"),        wt = 1/2))

  w <- setNames(numeric(length(b)), nms)
  found <- character(0); missing <- character(0)
  for (x in wanted) {
    t <- find_term(nms, x$parts)
    if (is.na(t)) {
      missing <- c(missing, paste(x$parts, collapse = ":"))
    } else {
      w[t] <- x$wt
      found <- c(found, t)
    }
  }
  if (length(missing) > 0 && strict) {
    stop("cannot form the marginal average — these terms were not found in the ",
         "model: ", paste(missing, collapse = ", "),
         "\n  Averaging over region and sex requires all of them. Check that ",
         "the model\n  used the project formula ",
         "(Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c.",
         "\n  If this is a reduced sensitivity model, pass strict = FALSE.",
         call. = FALSE)
  }
  if (length(missing) > 0 && verbose) {
    cat("    absent (zero weight):", paste(missing, collapse = ", "), "\n")
  }
  if (verbose) cat("    contrast terms:", paste(found, collapse = " + "), "\n")

  est <- sum(w * b)
  se  <- sqrt(drop(t(w) %*% V %*% w))
  z   <- est / se

  # Zero-inflation scales the marginal rate by the ratio of non-zero
  # probabilities. Its uncertainty is NOT propagated, so quote the conditional
  # interval with the conditional ratio and report the ZI shift separately.
  zi <- fixef(model)$zi
  p_nd <- NA_real_; p_t1d <- NA_real_; zi_scale <- NA_real_
  if (length(zi) && "DiagnosisT1D" %in% names(zi)) {
    p_nd     <- plogis(zi[["(Intercept)"]])
    p_t1d    <- plogis(zi[["(Intercept)"]] + zi[["DiagnosisT1D"]])
    zi_scale <- (1 - p_t1d) / (1 - p_nd)
  }

  tibble(log_ratio = est, SE = se,
         ratio = exp(est),
         lwr = exp(est - 1.96 * se), upr = exp(est + 1.96 * se),
         pct_change = 100 * (exp(est) - 1),
         z = z, p = 2 * pnorm(-abs(z)),
         zi_ND = p_nd, zi_T1D = p_t1d, zi_scale = zi_scale,
         marginal_ratio = exp(est) * ifelse(is.na(zi_scale), 1, zi_scale),
         n_terms = length(found))
}

# Observed donor-level rate: marker-positive cells summed over a donor's
# islets, divided by the summed count of cells of that type. This is what the
# points in the figures show. It needs no back-transformation, so unlike a
# model-predicted mean it cannot be wrong on a scale — which is what makes it
# usable as a check.
observed_rates <- function(sub, def) {
  sub %>%
    group_by(Diagnosis) %>%
    summarise(rate   = 100 * sum(.data[[def$num]]) / sum(.data[[def$den]]),
              donors = dplyr::n_distinct(Donor),
              islets = n(), .groups = "drop")
}

marker_stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",
                           ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))

# All six models, model estimate and observed rate side by side.
marker_effects_table <- function(model_dir = MARKER_MODEL_DIR,
                                 stain_handling = "double_only") {
  raw <- load_marker_data(stain_handling)
  bind_rows(lapply(marker_defs, function(def) {
    r <- load_marker_model(def, "rate", model_dir)
    if (is.null(r)) { warning("no model for ", def$label); return(NULL) }
    sub <- marker_subset(raw[[sub("_.*", "", def$key)]], def)
    obs <- observed_rates(sub, def)
    marginal_diagnosis(r$model, verbose = TRUE) %>%
      mutate(label = def$label, marker = def$marker, cell = def$cell,
             family = r$family,
             obs_ND  = obs$rate[obs$Diagnosis == "ND"],
             obs_T1D = obs$rate[obs$Diagnosis == "T1D"],
             obs_ratio = obs_T1D / obs_ND,
             n_ND  = obs$donors[obs$Diagnosis == "ND"],
             n_T1D = obs$donors[obs$Diagnosis == "T1D"],
             sig = marker_stars(p), .before = 1)
  }))
}

# -----------------------------------------------------------------------------
# VALIDATION GATE
# -----------------------------------------------------------------------------
# Purpose: catch a scale or transformation bug of the kind the emmeans path
# produced, where the model-derived ratio bore no relation to the data (0.94
# where the truth was 0.24, 1.09 where it was 9.2 — disagreements of 5.2x and
# 0.17x on effects the model itself called significant).
#
# Two things it must NOT fire on, both learned the hard way:
#
#   1. A null effect whose point estimates straddle 1. CD3 in alpha cells gives
#      observed 0.83 against model 1.11 — only 33% apart, but on opposite sides
#      of 1. The sign of a null estimate is undefined, so a direction test on it
#      fires at random. Direction is therefore only checked when the model
#      effect is significant.
#
#   2. A genuine crude-versus-adjusted divergence. Ki-67 in acinar cells gives
#      observed 0.46 against model 1.43. That is not a bug: the crude ratio
#      weights donors by how many islets they contribute, while the model
#      shrinks donors toward a common mean, and with a donor intercept SD of
#      1.26 on the log scale those two weightings can legitimately part company.
#      It means the estimate is UNRELIABLE, not that it is wrongly computed, so
#      it warns rather than stops.
#
# A hard stop is reserved for the combination that would actually put a wrong
# number in a manuscript: a large disagreement on an effect the model reports
# as significant.
validate_marker_rates <- function(tab, tolerance = 3, alpha = 0.05) {
  cat("\n=== VALIDATION: model vs observed rate ratio ===\n")
  cat("Stop   = disagreement >", tolerance, "x on a significant effect.\n")
  cat("Warn   = disagreement >", tolerance, "x on a null effect (unreliable,\n")
  cat("         report the observed rates instead).\n")
  cat("Direction is only checked where the model effect is significant.\n\n")

  chk <- tab %>%
    mutate(
      significant       = !is.na(p) & p < alpha,
      fold_disagreement = marginal_ratio / obs_ratio,
      big_gap           = !is.finite(fold_disagreement) |
                            fold_disagreement > tolerance |
                            fold_disagreement < 1 / tolerance,
      direction_clash   = significant &
                            sign(log(marginal_ratio)) != sign(log(obs_ratio)),
      verdict = case_when(
        !is.finite(fold_disagreement)  ~ "STOP — cannot check",
        direction_clash                ~ "STOP — direction clash on a significant effect",
        big_gap &  significant         ~ "STOP — scale error on a significant effect",
        big_gap & !significant         ~ "WARN — unreliable, do not report",
        TRUE                           ~ "ok"))

  print(as.data.frame(chk %>% dplyr::select(
    marker, cell, obs_ND, obs_T1D, obs_ratio, ratio, marginal_ratio,
    fold_disagreement, p, verdict)), digits = 3)

  warn <- chk %>% filter(grepl("^WARN", verdict))
  if (nrow(warn) > 0) {
    cat("\nUNRELIABLE (", nrow(warn), "): ",
        paste(sprintf("%s %s", warn$marker, warn$cell), collapse = ", "),
        "\n  Crude and adjusted estimates disagree on a non-significant effect.",
        "\n  Treat as uninformative rather than null, and quote the observed",
        "\n  rates if these appear anywhere.\n", sep = "")
  }

  bad <- chk %>% filter(grepl("^STOP", verdict))
  if (nrow(bad) > 0) {
    cat("\n")
    stop(nrow(bad), " model(s) FAILED: ",
         paste(sprintf("%s %s", bad$marker, bad$cell), collapse = ", "),
         "\n  A significant effect disagrees with the observed rates. This is the",
         "\n  signature of a scale or transformation bug. Do not use these",
         "\n  numbers. The observed columns are the ground truth.", call. = FALSE)
  }

  cat("\n", sum(chk$verdict == "ok"), " of ", nrow(chk),
      " models passed cleanly.\n", sep = "")
  invisible(chk)
}
