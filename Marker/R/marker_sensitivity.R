# =============================================================================
# marker_sensitivity.R
#
# Two sensitivity analyses for the marker rate models, plus a function that
# saves the descriptive values quoted in the text so they are reproducible
# rather than existing only in a chat log.
#
#   age_ge5     Excluding donors under 5 years of age. Addresses the concern
#               that very young ND donors drive the ND proliferation rate,
#               since Ki-67 positivity falls steeply with age.
#   early_t1d   T1D donors with disease duration <= 1.5 years versus all ND
#               donors. Asks whether the effect is present near onset rather
#               than only in long-standing disease.
#
# Each refits the primary model on the restricted subset with the SAME formula,
# family, zero-inflation structure and random effects, so the only thing that
# changes is the data. The Diagnosis contrast is then recomputed with the same
# marginal_diagnosis() used for the primary results.
#
# A restricted subset can be too thin for 17 fixed effects. Rather than report
# a degenerate fit, the code steps down a formula ladder and records which
# formula was used, so a reduced fit is visible instead of silent.
# =============================================================================

library(tidyverse)
library(glmmTMB)

source("Marker/R/marker_common.R")
source("Marker/R/marker_effects.R")

SENS_TMB_CTRL <- glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1000))
SENS_RE <- "(1 | Donor) + (1 | Donor:ImageID)"

# Formula ladder, most complete first.
SENS_RHS <- c(
  full    = "(Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c",
  reduced = paste("Diagnosis + Region + Sex + Age_c",
                  "+ Diagnosis:Region + Diagnosis:Sex + Diagnosis:Age_c"),
  minimal = "Diagnosis + Region + Sex + Age_c")

SENS_DEFS <- list(
  list(id = "age_ge5",
       label = "Excluding donors < 5 years of age",
       apply = function(d) d[!is.na(d$Age) & d$Age >= 5, , drop = FALSE]),
  list(id = "early_t1d",
       label = "T1D duration <= 1.5 years vs ND",
       apply = function(d) d[d$Diagnosis == "ND" |
                              (!is.na(d$Disease.Duration) &
                                 d$Disease.Duration <= 1.5), , drop = FALSE])
)

# -----------------------------------------------------------------------------

converged_ok <- function(m) {
  if (is.null(m)) return(FALSE)
  tryCatch(isTRUE(m$sdr$pdHess) &&
             identical(as.integer(m$fit$convergence), 0L) &&
             all(is.finite(sqrt(diag(vcov(m, full = TRUE))))),
           error = function(e) FALSE)
}

# A fit can converge and still be useless: a huge estimate with a huge SE is a
# term the data cannot identify. Treat that as a failure so the ladder steps on.
degenerate <- function(m, est_max = 15, se_max = 50) {
  co <- tryCatch(summary(m)$coefficients$cond, error = function(e) NULL)
  if (is.null(co)) return(TRUE)
  any(!is.finite(co[, 1]) | !is.finite(co[, 2]) |
        abs(co[, 1]) > est_max | co[, 2] > se_max)
}

# Reconstruct the primary model's family rather than relying on a cached
# attribute, so this works on models fitted before fit_spec was added.
primary_family <- function(model) {
  fam <- family(model)$family
  has_zi <- length(fixef(model)$zi) > 0
  list(family = switch(fam, nbinom1 = nbinom1(), nbinom2 = nbinom2(),
                       stop("unexpected family: ", fam)),
       name = paste0(if (has_zi) "zi" else "", fam),
       ziformula = if (has_zi) ~ Diagnosis else ~ 0)
}

fit_restricted <- function(def, sub, fam) {
  for (nm in names(SENS_RHS)) {
    f <- as.formula(sprintf("%s ~ %s + offset(log(%s)) + %s",
                            def$num, SENS_RHS[[nm]], def$den, SENS_RE))
    environment(f) <- globalenv()
    m <- tryCatch(glmmTMB(f, data = sub, family = fam$family,
                          ziformula = fam$ziformula, control = SENS_TMB_CTRL),
                  error = function(e) NULL)
    if (converged_ok(m) && !degenerate(m)) {
      return(list(model = m, formula = nm))
    }
  }
  list(model = NULL, formula = NA_character_)
}

# -----------------------------------------------------------------------------

run_marker_sensitivity <- function(
    panels    = c("ki67_ins", "ki67_glu", "cd3_ins", "cd3_glu"),
    model_dir = MARKER_MODEL_DIR,
    out_dir   = "Marker/Results",
    stain_handling = "double_only") {

  cat("\n=============================================================================\n")
  cat("MARKER SENSITIVITY ANALYSES\n")
  cat("=============================================================================\n")

  raw <- load_marker_data(stain_handling)
  rows <- list()

  for (def in marker_defs) {
    if (!def$key %in% panels) next
    r <- load_marker_model(def, "rate", model_dir)
    if (is.null(r)) { cat("\n[skip]", def$label, "- no cached model\n"); next }

    full <- marker_subset(raw[[sub("_.*", "", def$key)]], def)
    fam  <- primary_family(r$model)
    base <- marginal_diagnosis(r$model, strict = TRUE)

    cat(sprintf("\n--- %s (%s) ---\n", def$label, fam$name))
    cat(sprintf("  primary: RR %.3f (%.3f-%.3f), p = %.3g, %d donors\n",
                base$ratio, base$lwr, base$upr, base$p,
                dplyr::n_distinct(full$Donor)))

    rows[[length(rows) + 1L]] <- tibble(
      panel = def$label, marker = def$marker, cell = def$cell,
      analysis = "Primary", formula_used = "full",
      n_ND = dplyr::n_distinct(full$Donor[full$Diagnosis == "ND"]),
      n_T1D = dplyr::n_distinct(full$Donor[full$Diagnosis == "T1D"]),
      islets = nrow(full),
      obs_ND = 100 * sum(full[[def$num]][full$Diagnosis == "ND"]) /
        sum(full[[def$den]][full$Diagnosis == "ND"]),
      obs_T1D = 100 * sum(full[[def$num]][full$Diagnosis == "T1D"]) /
        sum(full[[def$den]][full$Diagnosis == "T1D"]),
      ratio = base$ratio, lwr = base$lwr, upr = base$upr, p = base$p,
      sig = marker_stars(base$p), note = NA_character_)

    for (sd in SENS_DEFS) {
      subr <- droplevels(sd$apply(full))
      ndn  <- dplyr::n_distinct(subr$Donor[subr$Diagnosis == "ND"])
      t1dn <- dplyr::n_distinct(subr$Donor[subr$Diagnosis == "T1D"])
      dropped <- setdiff(as.character(full$Donor), as.character(subr$Donor))

      note <- NA_character_
      # A restricted subset can lose a whole factor level, which silently
      # changes the design. Record it rather than letting it pass.
      lost <- c(if (dplyr::n_distinct(subr$Region) < 3) "a Region level",
                if (dplyr::n_distinct(subr$Sex) < 2) "a Sex level",
                if (dplyr::n_distinct(subr$Diagnosis) < 2) "a Diagnosis level")
      if (length(lost)) note <- paste("lost", paste(lost, collapse = " and "))

      if (ndn < 3 || t1dn < 3) {
        cat(sprintf("  %-34s SKIPPED - only %d ND / %d T1D donors\n",
                    sd$label, ndn, t1dn))
        rows[[length(rows) + 1L]] <- tibble(
          panel = def$label, marker = def$marker, cell = def$cell,
          analysis = sd$label, formula_used = NA_character_,
          n_ND = ndn, n_T1D = t1dn, islets = nrow(subr),
          obs_ND = NA_real_, obs_T1D = NA_real_, ratio = NA_real_,
          lwr = NA_real_, upr = NA_real_, p = NA_real_, sig = "",
          note = "too few donors to fit")
        next
      }

      fit <- fit_restricted(def, subr, fam)
      if (is.null(fit$model)) {
        cat(sprintf("  %-34s FAILED to fit at any formula\n", sd$label))
        rows[[length(rows) + 1L]] <- tibble(
          panel = def$label, marker = def$marker, cell = def$cell,
          analysis = sd$label, formula_used = NA_character_,
          n_ND = ndn, n_T1D = t1dn, islets = nrow(subr),
          obs_ND = NA_real_, obs_T1D = NA_real_, ratio = NA_real_,
          lwr = NA_real_, upr = NA_real_, p = NA_real_, sig = "",
          note = paste(na.omit(c("did not converge", note)), collapse = "; "))
        next
      }

      # strict = FALSE: a reduced formula legitimately lacks some interaction
      # terms, which then carry zero weight in the average.
      e <- marginal_diagnosis(fit$model, strict = FALSE)
      oND <- 100 * sum(subr[[def$num]][subr$Diagnosis == "ND"]) /
        sum(subr[[def$den]][subr$Diagnosis == "ND"])
      oT1 <- 100 * sum(subr[[def$num]][subr$Diagnosis == "T1D"]) /
        sum(subr[[def$den]][subr$Diagnosis == "T1D"])

      cat(sprintf("  %-34s RR %.3f (%.3f-%.3f), p = %.3g  [%s formula, %d ND / %d T1D]\n",
                  sd$label, e$ratio, e$lwr, e$upr, e$p, fit$formula, ndn, t1dn))
      if (length(dropped)) {
        cat(sprintf("      donors dropped (%d): %s\n", length(dropped),
                    paste(head(dropped, 10), collapse = ", ")))
      }

      rows[[length(rows) + 1L]] <- tibble(
        panel = def$label, marker = def$marker, cell = def$cell,
        analysis = sd$label, formula_used = fit$formula,
        n_ND = ndn, n_T1D = t1dn, islets = nrow(subr),
        obs_ND = oND, obs_T1D = oT1,
        ratio = e$ratio, lwr = e$lwr, upr = e$upr, p = e$p,
        sig = marker_stars(e$p), note = note)
    }
  }

  out <- bind_rows(rows)

  # Does each sensitivity reach the same conclusion as its primary?
  verdict <- out %>%
    group_by(panel) %>%
    mutate(prim_ratio = ratio[analysis == "Primary"][1],
           prim_sig   = p[analysis == "Primary"][1] < 0.05) %>%
    ungroup() %>%
    filter(analysis != "Primary") %>%
    mutate(same_direction = is.na(ratio) | is.na(prim_ratio) |
             sign(log(ratio)) == sign(log(prim_ratio)),
           same_verdict = is.na(p) | ((p < 0.05) == prim_sig),
           holds = case_when(
             is.na(ratio)                      ~ "not estimable",
             same_direction & same_verdict      ~ "holds",
             same_direction & !same_verdict     ~ "same direction, significance changes",
             TRUE                               ~ "DIRECTION CHANGES"))

  cat("\n=== summary ===\n")
  print(as.data.frame(out %>% dplyr::select(
    marker, cell, analysis, formula_used, n_ND, n_T1D,
    obs_ND, obs_T1D, ratio, lwr, upr, p, sig)), digits = 3)

  cat("\n=== does each sensitivity reach the same conclusion? ===\n")
  print(as.data.frame(verdict %>% dplyr::select(
    marker, cell, analysis, prim_ratio, ratio, p, holds)), digits = 3)
  cat("\nA loss of significance on a much smaller subset is expected and is not\n")
  cat("evidence against the primary result; a change of DIRECTION would be.\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(out, file.path(out_dir, "marker_sensitivity.csv"))
  cat("\nSaved", file.path(out_dir, "marker_sensitivity.csv"), "\n")
  invisible(out)
}

# =============================================================================
# DESCRIPTIVE VALUES QUOTED IN THE TEXT AND SUPPLEMENTARY TABLE
# =============================================================================
# These were previously computed ad hoc. Saving them makes every number in the
# table traceable to a file.
#
#   cells per islet, by type         -> the "alpha-cells per islet" fold changes
#   marker+ cells per 100 islet cells -> "CD3+ cells per islet cell"
#   relative enrichment              -> "beta-cell co-localization per CD3+ cell"
#
# Relative enrichment = (share of marker+ cells in a bin) / (share of cells of
# that type). 1.0 means the marker is distributed in proportion to what is
# available, so the bin carries no information beyond availability. Computed
# over ALL retained islets, which is the basis the manuscript quotes.
# =============================================================================

marker_descriptives <- function(out_dir = "Marker/Results",
                                stain_handling = "double_only") {
  raw <- load_marker_data(stain_handling)

  comp <- bind_rows(lapply(names(raw), function(tag) {
    raw[[tag]] %>%
      group_by(Diagnosis) %>%
      summarise(islets = n(), donors = dplyr::n_distinct(Donor),
                beta_per_islet   = mean(Count.ins),
                alpha_per_islet  = mean(Count.glu),
                acinar_per_islet = mean(Count.Acinar),
                total_endocrine  = mean(Islet.Cells), .groups = "drop") %>%
      mutate(marker = tag, .before = 1)
  }))

  folds <- comp %>%
    dplyr::select(marker, Diagnosis, beta_per_islet, alpha_per_islet,
                  acinar_per_islet, total_endocrine) %>%
    pivot_longer(-c(marker, Diagnosis), names_to = "measure") %>%
    pivot_wider(names_from = Diagnosis, values_from = value) %>%
    mutate(fold_T1D_over_ND = T1D / ND)

  enrich <- bind_rows(lapply(names(raw), function(tag) {
    d <- raw[[tag]]
    d <- d %>% mutate(av = Count.ins + Count.glu + Count.Acinar)
    bind_rows(lapply(list(c("beta","insulin_cells","Count.ins"),
                          c("alpha","glucagon_cells","Count.glu"),
                          c("acinar","acinar_cells","Count.Acinar")), function(x) {
      cell <- x[1]; num <- x[2]; den <- x[3]
      g <- d %>% group_by(Diagnosis) %>%
        summarise(share_of_marker = sum(.data[[num]]) / sum(total_cells),
                  share_of_cells  = sum(.data[[den]]) / sum(av),
                  marker_per_100_cells = 100 * sum(total_cells) / sum(av),
                  .groups = "drop") %>%
        mutate(enrichment = share_of_marker / share_of_cells)
      tibble(marker = tag, cell = cell,
             marker_per_100_cells_ND  = g$marker_per_100_cells[g$Diagnosis == "ND"],
             marker_per_100_cells_T1D = g$marker_per_100_cells[g$Diagnosis == "T1D"],
             share_of_marker_ND  = g$share_of_marker[g$Diagnosis == "ND"],
             share_of_marker_T1D = g$share_of_marker[g$Diagnosis == "T1D"],
             enrichment_ND  = g$enrichment[g$Diagnosis == "ND"],
             enrichment_T1D = g$enrichment[g$Diagnosis == "T1D"],
             enrichment_fold = enrichment_T1D / enrichment_ND)
    }))
  }))

  cat("\n=== cells per islet, and fold changes ===\n")
  print(as.data.frame(folds), digits = 4)
  cat("\n=== marker density and relative enrichment ===\n")
  print(as.data.frame(enrich), digits = 4)
  cat("\nmarker_per_100_cells is the 'CD3+ cells per islet cell' row;\n")
  cat("enrichment_fold for cell = beta is the 'per CD3+ cell' row.\n")

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(folds,  file.path(out_dir, "marker_descriptives_folds.csv"))
  write_csv(enrich, file.path(out_dir, "marker_descriptives_enrichment.csv"))
  cat("\nSaved marker_descriptives_folds.csv and marker_descriptives_enrichment.csv\n")
  invisible(list(folds = folds, enrichment = enrich))
}
