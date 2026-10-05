# =============================================================================
# plot_area_all_comparisons.R
# All possible post-hoc comparison plots for Cell Area models (glmmTMB)
# =============================================================================
#
# All models: Gaussian on log(area) — emmeans returns log scale, exponentiate
#
# Per cell type (beta, alpha, delta, pp):
#   DIAGNOSIS tier (4 models):
#     Singlets  (reduced: no log_IC_c) → 15 plots
#     Small EOs (seo: log_IC_c main effect) → 15 plots
#     Islets    (full: has log_IC_c)   → 21 plots
#     Combined  (object_type interactions) → ~25 plots
#   T1D tier (4 models, AO parameterization):
#     Singlets  (reduced) → 12 plots
#     Small EOs (seo: log_IC_c main effect) → 12 plots
#     Islets    (full)    → 19 plots
#     Combined  (object_type) → ~23 plots
#
# Total: ~488 plots across 32 models
#
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(patchwork)

source("colors_master.R")

# =============================================================================
# THEME
# =============================================================================

theme_pub <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.margin = margin(8, 12, 8, 8),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )

p_to_stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")
emm_fn <- function(model, specs, data, ...) emmeans(model, specs, data = data, ...)
emt_fn <- function(model, specs, var, data, ...) emtrends(model, specs, var = var, data = data, ...)

y_lab_area <- expression(paste("Cell area (", mu, "m"^2, ")"))

# =============================================================================
# CORE: Exponentiate emmeans (Gaussian on log-area)
# =============================================================================

find_ci <- function(df) {
  lower <- if ("asymp.LCL" %in% names(df)) "asymp.LCL" else "lower.CL"
  upper <- if ("asymp.UCL" %in% names(df)) "asymp.UCL" else "upper.CL"
  list(lower = lower, upper = upper)
}

exp_emm <- function(model, specs, data = NULL, at = NULL) {
  args <- list(object = model, specs = specs)
  if (!is.null(data)) args$data <- data
  if (!is.null(at)) args$at <- at
  e <- as.data.frame(summary(do.call(emmeans, args)))
  ci <- find_ci(e)
  e$response <- exp(e$emmean)
  e[[ci$lower]] <- exp(e[[ci$lower]])
  e[[ci$upper]] <- exp(e[[ci$upper]])
  e$lo <- e[[ci$lower]]; e$hi <- e[[ci$upper]]
  e
}

save_plot <- function(plt, filename, output_dir, width = 7, height = 5) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("tiff"))
    ggsave(file.path(output_dir, paste0(filename, ".", ext)), plt,
           width = width, height = height, dpi = 300,
           bg = if (ext != "pdf") "white" else NULL)
  cat("    Saved:", filename, "\n")
}


# #############################################################################
#
#  DIAGNOSIS SUBSET PLOTS (singlets, endobs, islets)
#
# #############################################################################

generate_area_diag_plots <- function(model, data, cell_label, area_col,
                                      has_islet, output_dir) {
  cat("    Generating diagnosis plots for", cell_label, "...\n")

  # Donor geometric means (matches model's log-scale estimates)
  dm_cat <- function(...) {
    data %>% group_by(Donor, ...) %>%
      summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)), .groups = "drop")
  }

  # --- Main effects ---
  for (fac in c("Diagnosis", "Region", "Sex")) {
    dm <- dm_cat(!!sym(fac))
    e <- exp_emm(model, as.formula(paste("~", fac)), data)
    cols <- if (fac == "Diagnosis") colors_diag else if (fac == "Region") colors_region else colors_sex
    if (fac == "Region") {
      prs <- as.data.frame(pairs(emmeans(model, ~ Region, data = data), adjust = "tukey"))
      sig <- prs %>% filter(p.value < 0.05) %>% mutate(label = paste(contrast, p_to_stars(p.value)))
      sig_text <- paste(sig$label, collapse = "\n")
    } else {
      prs <- as.data.frame(pairs(emmeans(model, as.formula(paste("~", fac)), data = data),
                                  adjust = "none", reverse = (fac == "Diagnosis")))
      sig_text <- p_to_stars(prs$p.value[1])
    }

    plt <- ggplot() +
      geom_violin(data = dm, aes(x = .data[[fac]], y = md, fill = .data[[fac]]),
                  alpha = 0.3, color = NA, width = 0.7) +
      geom_jitter(data = dm, aes(x = .data[[fac]], y = md, color = .data[[fac]]),
                  width = 0.15, size = 1.5, alpha = 0.6) +
      scale_fill_manual(values = cols, guide = "none") +
      scale_color_manual(values = cols, guide = "none") +
      labs(x = if (fac == "Region") "Pancreatic region" else NULL, y = y_lab_area,
           title = paste(cell_label, "\u2014", fac)) + theme_pub
    if (nchar(sig_text) > 0)
      plt <- plt + annotate("text", x = if (fac == "Region") 2 else 1.5,
                             y = max(dm$md) * 0.95, label = sig_text,
                             size = if (fac == "Region") 4 else 8, lineheight = 0.9)
    save_plot(plt, paste0("main_", tolower(fac)), output_dir)
  }

  # --- Two-way categorical (6 views) ---
  two_way <- list(
    list(fac = "Diagnosis", cond = "Region"), list(fac = "Region", cond = "Diagnosis"),
    list(fac = "Diagnosis", cond = "Sex"),    list(fac = "Sex", cond = "Diagnosis"),
    list(fac = "Region", cond = "Sex"),       list(fac = "Sex", cond = "Region")
  )
  for (tw in two_way) {
    fac <- tw$fac; cond <- tw$cond
    dm <- dm_cat(!!sym(fac), !!sym(cond))
    e <- exp_emm(model, as.formula(paste("~", fac, "|", cond)), data)
    cols_f <- if (fac == "Diagnosis") colors_diag else if (fac == "Region") colors_region else colors_sex
    n_fac <- length(unique(data[[fac]]))

    if (n_fac <= 2) {
      prs <- as.data.frame(pairs(emmeans(model, as.formula(paste("~", fac, "|", cond)), data = data),
                                  adjust = "none", reverse = (fac == "Diagnosis")))
      n_c <- length(unique(data[[cond]]))
      plt <- ggplot() +
        geom_violin(data = dm, aes(x = .data[[cond]], y = md, fill = .data[[fac]]),
                    alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
        geom_jitter(data = dm, aes(x = .data[[cond]], y = md, color = .data[[fac]], group = .data[[fac]]),
                    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.2, alpha = 0.5) +
        scale_fill_manual(values = cols_f) + scale_color_manual(values = cols_f) +
        annotate("text", x = seq_len(n_c), y = max(dm$md) * seq(0.98, by = -0.05, length.out = n_c),
                 label = sapply(prs$p.value, p_to_stars), size = 8) +
        labs(x = if (cond == "Region") "Pancreatic region" else NULL, y = y_lab_area,
             title = paste(cell_label, "\u2014", fac, "\u00D7", cond),
             fill = fac, color = fac) +
        theme_pub + theme(legend.position = c(0.90, 0.85), legend.background = element_rect(fill = "white", color = NA))
    } else {
      prs <- as.data.frame(pairs(emmeans(model, as.formula(paste("~", fac, "|", cond)), data = data), adjust = "tukey"))
      sig_labels <- prs %>% filter(p.value < 0.05) %>%
        mutate(label = paste(contrast, p_to_stars(p.value))) %>%
        group_by(!!sym(cond)) %>% summarise(label = paste(label, collapse = "\n"), .groups = "drop")
      all_c <- tibble(!!cond := unique(data[[cond]]))
      sig_labels <- left_join(all_c, sig_labels, by = cond) %>% mutate(label = replace_na(label, ""))

      plt <- ggplot() +
        geom_violin(data = dm, aes(x = .data[[fac]], y = md, fill = .data[[fac]]),
                    alpha = 0.3, color = NA, width = 0.7) +
        geom_jitter(data = dm, aes(x = .data[[fac]], y = md, color = .data[[fac]]),
                    width = 0.15, size = 1.2, alpha = 0.5) +
        facet_wrap(as.formula(paste("~", cond))) +
        geom_text(data = sig_labels, aes(x = 2, y = max(dm$md) * 0.95, label = label),
                  size = 3.5, lineheight = 0.9) +
        scale_fill_manual(values = cols_f, guide = "none") + scale_color_manual(values = cols_f, guide = "none") +
        labs(x = if (fac == "Region") "Pancreatic region" else NULL, y = y_lab_area,
             title = paste(cell_label, "\u2014", fac, "\u00D7", cond)) + theme_pub
    }
    save_plot(plt, paste0(tolower(fac), "_by_", tolower(cond)), output_dir,
              if (fac == "Region" || cond == "Region") 9 else 7, 5)
  }

  # --- Continuous slopes (emmeans-based) ---
  cont_specs <- list(list(var = "Age_c", raw = "Age", lab = "Age (years)", tag = "age"))
  if (has_islet) cont_specs[[2]] <- list(var = "log_Islet.Cells_c", raw = "log_Islet.Cells",
                                          lab = "Islet size (log cells)", tag = "islet")
  for (cs in cont_specs) {
    cv <- cs$var; rv <- cs$raw; xl <- cs$lab; tag <- cs$tag
    ar <- range(data[[rv]], na.rm = TRUE); mc <- mean(data[[rv]], na.rm = TRUE)
    cv_seq <- seq(ar[1] - mc, ar[2] - mc, length.out = 50)
    at_list <- setNames(list(cv_seq), cv)
    cat("      ", tag, "slopes...\n")

    # Overall
    e_ov <- exp_emm(model, as.formula(paste("~", cv)), data, at = at_list)
    e_ov$x <- e_ov[[cv]] + mc
    dm <- data %>% group_by(Donor) %>% summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)),
                                                   x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
    trend <- as.data.frame(test(emt_fn(model, ~ 1, cv, data))); stars <- p_to_stars(trend$p.value[1])
    plt <- ggplot() +
      geom_ribbon(data = e_ov, aes(x = x, ymin = lo, ymax = hi), fill = "grey60", alpha = 0.3) +
      geom_line(data = e_ov, aes(x = x, y = response), linewidth = 0.9) +
      geom_point(data = dm, aes(x = x, y = md), alpha = 0.5, size = 1.8) +
      labs(x = xl, y = y_lab_area, title = paste(cell_label, "\u2014", tag, "(overall)")) + theme_pub
    if (nchar(stars) > 0) plt <- plt + annotate("text", x = ar[1] + 1, y = max(dm$md) * 0.95,
                                                  label = stars, size = 8, hjust = 0)
    save_plot(plt, paste0(tag, "_overall"), output_dir)

    # By categorical
    for (bi in list(list(v = "Diagnosis", cols = colors_diag),
                     list(v = "Region", cols = colors_region),
                     list(v = "Sex", cols = colors_sex))) {
      bv <- bi$v
      e_by <- exp_emm(model, as.formula(paste("~", cv, "|", bv)), data, at = at_list)
      e_by$x <- e_by[[cv]] + mc
      dm2 <- data %>% group_by(Donor, !!sym(bv)) %>%
        summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)),
                  x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
      plt2 <- ggplot() +
        geom_ribbon(data = e_by, aes(x = x, ymin = lo, ymax = hi, fill = .data[[bv]]), alpha = 0.2) +
        geom_line(data = e_by, aes(x = x, y = response, color = .data[[bv]]), linewidth = 0.9) +
        geom_point(data = dm2, aes(x = x, y = md, color = .data[[bv]]), alpha = 0.5, size = 1.8) +
        scale_color_manual(values = bi$cols) + scale_fill_manual(values = bi$cols) +
        labs(x = xl, y = y_lab_area, title = paste0(cell_label, " \u2014 ", bv, " \u00D7 ", tag)) +
        theme_pub + theme(legend.position = c(0.85, 0.85), legend.background = element_rect(fill = "white", color = NA))
      save_plot(plt2, paste0(tag, "_by_", tolower(bv)), output_dir)
    }

    # 3-way: facet Region + Diagnosis
    e_3 <- exp_emm(model, as.formula(paste("~", cv, "| Diagnosis * Region")), data, at = at_list)
    e_3$x <- e_3[[cv]] + mc
    dm3 <- data %>% group_by(Donor, Diagnosis, Region) %>%
      summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)),
                x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
    ym <- max(c(dm3$md, e_3$hi), na.rm = TRUE) * 1.1

    plt_a <- ggplot() +
      geom_ribbon(data = e_3, aes(x = x, ymin = lo, ymax = hi, fill = Diagnosis), alpha = 0.2) +
      geom_line(data = e_3, aes(x = x, y = response, color = Diagnosis), linewidth = 0.9) +
      geom_point(data = dm3, aes(x = x, y = md, color = Diagnosis), alpha = 0.4, size = 1.5) +
      facet_wrap(~ Region, ncol = 3) +
      scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
      scale_y_continuous(limits = c(0, ym)) +
      labs(x = xl, y = y_lab_area,
           title = paste0(cell_label, ": Diag \u00D7 Region \u00D7 ", tag),
           subtitle = "Facet by Region") +
      theme_pub + theme(legend.position = "bottom")
    save_plot(plt_a, paste0("3way_", tag, "_facet_region"), output_dir, 12, 5)

    plt_b <- ggplot() +
      geom_ribbon(data = e_3, aes(x = x, ymin = lo, ymax = hi, fill = Region), alpha = 0.2) +
      geom_line(data = e_3, aes(x = x, y = response, color = Region), linewidth = 0.9) +
      geom_point(data = dm3, aes(x = x, y = md, color = Region), alpha = 0.4, size = 1.5) +
      facet_wrap(~ Diagnosis, ncol = 2) +
      scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
      scale_y_continuous(limits = c(0, ym)) +
      labs(x = xl, y = y_lab_area,
           title = paste0(cell_label, ": Diag \u00D7 Region \u00D7 ", tag),
           subtitle = "Facet by Diagnosis") +
      theme_pub + theme(legend.position = "bottom")
    save_plot(plt_b, paste0("3way_", tag, "_facet_diagnosis"), output_dir, 10, 5)
  }
}


# #############################################################################
#
#  DIAGNOSIS COMBINED PLOTS (object_type interactions)
#
# #############################################################################

generate_area_combined_diag_plots <- function(model, data, cell_label, area_col, output_dir) {
  cat("    Generating combined diagnosis plots for", cell_label, "...\n")

  # First: all standard plots (has_islet = FALSE since log_IC_c is additive, not in ^2)
  # The combined model has log_IC_c as a covariate but NOT in the interaction expansion
  generate_area_diag_plots(model, data, cell_label, area_col, has_islet = FALSE, output_dir)

  # Additional: object_type effects
  dm_ot <- data %>% group_by(Donor, object_type) %>%
    summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)), .groups = "drop")
  e_ot <- exp_emm(model, ~ object_type, data)
  prs_ot <- as.data.frame(pairs(emmeans(model, ~ object_type, data = data), adjust = "tukey"))
  sig <- prs_ot %>% filter(p.value < 0.05) %>% mutate(label = paste(contrast, p_to_stars(p.value)))
  sig_text <- paste(sig$label, collapse = "\n")

  plt <- ggplot() +
    geom_violin(data = dm_ot, aes(x = object_type, y = md, fill = object_type),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm_ot, aes(x = object_type, y = md, color = object_type),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    scale_color_brewer(palette = "Set2", guide = "none") +
    { if (nchar(sig_text) > 0) annotate("text", x = 2, y = max(dm_ot$md) * 0.95,
                                         label = sig_text, size = 4, lineheight = 0.9) } +
    labs(x = "Object type", y = y_lab_area,
         title = paste(cell_label, "\u2014 Object type")) + theme_pub
  save_plot(plt, "combined_object_type", output_dir)

  # object_type × Diagnosis
  dm_od <- data %>% group_by(Donor, object_type, Diagnosis) %>%
    summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)), .groups = "drop")
  e_od <- exp_emm(model, ~ Diagnosis | object_type, data)
  prs_od <- as.data.frame(pairs(emmeans(model, ~ Diagnosis | object_type, data = data),
                                 adjust = "none", reverse = TRUE))
  plt <- ggplot() +
    geom_violin(data = dm_od, aes(x = object_type, y = md, fill = Diagnosis),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm_od, aes(x = object_type, y = md, color = Diagnosis, group = Diagnosis),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_diag) + scale_color_manual(values = colors_diag) +
    annotate("text", x = 1:3, y = max(dm_od$md) * c(0.98, 0.93, 0.88),
             label = sapply(prs_od$p.value, p_to_stars), size = 8) +
    labs(x = "Object type", y = y_lab_area,
         title = paste(cell_label, "\u2014 Diagnosis \u00D7 Object type"),
         fill = "Diagnosis", color = "Diagnosis") +
    theme_pub + theme(legend.position = c(0.90, 0.85), legend.background = element_rect(fill = "white", color = NA))
  save_plot(plt, "combined_diagnosis_by_objecttype", output_dir, 8)

  # object_type × Region
  dm_or <- data %>% group_by(Donor, object_type, Region) %>%
    summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)), .groups = "drop")
  e_or <- exp_emm(model, ~ Region | object_type, data)

  plt <- ggplot() +
    geom_violin(data = dm_or, aes(x = Region, y = md, fill = Region),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm_or, aes(x = Region, y = md, color = Region),
                width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ object_type, ncol = 3) +
    scale_fill_manual(values = colors_region, guide = "none") +
    scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = y_lab_area,
         title = paste(cell_label, "\u2014 Region \u00D7 Object type")) + theme_pub
  save_plot(plt, "combined_region_by_objecttype", output_dir, 12, 5)
}


# #############################################################################
#
#  T1D SUBSET PLOTS (AO parameterization)
#
# #############################################################################

generate_area_t1d_plots <- function(model, data, cell_label, area_col,
                                     has_islet, output_dir) {
  cat("    Generating T1D plots for", cell_label, "...\n")

  dm_cat <- function(...) {
    data %>% group_by(Donor, ...) %>%
      summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)), .groups = "drop")
  }

  # --- Categorical: Region, Sex, Region×Sex both ways ---
  for (fac in c("Region", "Sex")) {
    dm <- dm_cat(!!sym(fac))
    e <- exp_emm(model, as.formula(paste("~", fac)), data)
    cols <- if (fac == "Region") colors_region else colors_sex
    if (fac == "Region") {
      prs <- as.data.frame(pairs(emmeans(model, ~ Region, data = data), adjust = "tukey"))
      sig <- prs %>% filter(p.value < 0.05) %>% mutate(label = paste(contrast, p_to_stars(p.value)))
      sig_text <- paste(sig$label, collapse = "\n")
    } else {
      prs <- as.data.frame(pairs(emmeans(model, ~ Sex, data = data), adjust = "none"))
      sig_text <- p_to_stars(prs$p.value[1])
    }
    plt <- ggplot() +
      geom_violin(data = dm, aes(x = .data[[fac]], y = md, fill = .data[[fac]]),
                  alpha = 0.3, color = NA, width = 0.7) +
      geom_jitter(data = dm, aes(x = .data[[fac]], y = md, color = .data[[fac]]),
                  width = 0.15, size = 1.5, alpha = 0.6) +
      scale_fill_manual(values = cols, guide = "none") + scale_color_manual(values = cols, guide = "none") +
      labs(x = if (fac == "Region") "Pancreatic region" else NULL, y = y_lab_area,
           title = paste(cell_label, "\u2014", fac, "(T1D)")) + theme_pub
    if (nchar(sig_text) > 0)
      plt <- plt + annotate("text", x = if (fac == "Region") 2 else 1.5, y = max(dm$md) * 0.95,
                             label = sig_text, size = if (fac == "Region") 4 else 8, lineheight = 0.9)
    save_plot(plt, paste0("main_", tolower(fac)), output_dir)
  }

  # Region|Sex and Sex|Region
  dm_rs <- dm_cat(Region, Sex)
  e_rs <- exp_emm(model, ~ Region | Sex, data)
  prs_rs <- as.data.frame(pairs(emmeans(model, ~ Region | Sex, data = data), adjust = "tukey"))
  sig_rs <- prs_rs %>% filter(p.value < 0.05) %>%
    mutate(label = paste(contrast, p_to_stars(p.value))) %>%
    group_by(Sex) %>% summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_sex <- tibble(Sex = factor(c("Female", "Male"), levels = levels(data$Sex)))
  sig_rs <- left_join(all_sex, sig_rs, by = "Sex") %>% mutate(label = replace_na(label, ""))

  plt <- ggplot() +
    geom_violin(data = dm_rs, aes(x = Region, y = md, fill = Region), alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm_rs, aes(x = Region, y = md, color = Region), width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ Sex) +
    geom_text(data = sig_rs, aes(x = 2, y = max(dm_rs$md) * 0.95, label = label), size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = colors_region, guide = "none") + scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = y_lab_area, title = paste(cell_label, "\u2014 Region \u00D7 Sex (T1D)")) + theme_pub
  save_plot(plt, "region_by_sex", output_dir, 9, 5)

  e_sr <- exp_emm(model, ~ Sex | Region, data)
  prs_sr <- as.data.frame(pairs(emmeans(model, ~ Sex | Region, data = data), adjust = "none"))
  plt <- ggplot() +
    geom_violin(data = dm_rs, aes(x = Region, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm_rs, aes(x = Region, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) + scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:3, y = max(dm_rs$md) * c(0.98, 0.93, 0.88),
             label = sapply(prs_sr$p.value, p_to_stars), size = 8) +
    labs(x = "Pancreatic region", y = y_lab_area, title = paste(cell_label, "\u2014 Sex \u00D7 Region (T1D)"),
         fill = "Sex", color = "Sex") +
    theme_pub + theme(legend.position = c(0.90, 0.85), legend.background = element_rect(fill = "white", color = NA))
  save_plot(plt, "sex_by_region", output_dir, 8)

  # --- Continuous slopes ---
  cont_specs <- list(
    list(var = "Disease.Duration_c", raw = "Disease.Duration", lab = "Disease duration (years)", tag = "dd"),
    list(var = "age_at_onset_c", raw = "age_at_onset", lab = "Age at onset (years)", tag = "ao")
  )
  if (has_islet) cont_specs[[3]] <- list(var = "log_Islet.Cells_c", raw = "log_Islet.Cells",
                                          lab = "Islet size (log cells)", tag = "islet")
  for (cs in cont_specs) {
    cv <- cs$var; rv <- cs$raw; xl <- cs$lab; tag <- cs$tag
    ar <- range(data[[rv]], na.rm = TRUE); mc <- mean(data[[rv]], na.rm = TRUE)
    cv_seq <- seq(ar[1] - mc, ar[2] - mc, length.out = 50)
    at_list <- setNames(list(cv_seq), cv)
    cat("      ", tag, "slopes...\n")

    e_ov <- exp_emm(model, as.formula(paste("~", cv)), data, at = at_list)
    e_ov$x <- e_ov[[cv]] + mc
    dm <- data %>% group_by(Donor) %>% summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)),
                                                   x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
    trend <- as.data.frame(test(emt_fn(model, ~ 1, cv, data))); stars <- p_to_stars(trend$p.value[1])
    plt <- ggplot() +
      geom_ribbon(data = e_ov, aes(x = x, ymin = lo, ymax = hi), fill = col_t1d, alpha = 0.2) +
      geom_line(data = e_ov, aes(x = x, y = response), linewidth = 0.9, color = "black") +
      geom_point(data = dm, aes(x = x, y = md), color = col_t1d, alpha = 0.6, size = 2) +
      labs(x = xl, y = y_lab_area, title = paste(cell_label, "\u2014", tag, "overall (T1D)")) + theme_pub
    if (nchar(stars) > 0) plt <- plt + annotate("text", x = ar[1] + 1, y = max(dm$md) * 0.95,
                                                  label = stars, size = 8, hjust = 0)
    save_plot(plt, paste0(tag, "_overall"), output_dir)

    for (bi in list(list(v = "Region", cols = colors_region), list(v = "Sex", cols = colors_sex))) {
      bv <- bi$v
      e_by <- exp_emm(model, as.formula(paste("~", cv, "|", bv)), data, at = at_list)
      e_by$x <- e_by[[cv]] + mc
      dm2 <- data %>% group_by(Donor, !!sym(bv)) %>%
        summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)),
                  x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
      plt2 <- ggplot() +
        geom_ribbon(data = e_by, aes(x = x, ymin = lo, ymax = hi, fill = .data[[bv]]), alpha = 0.2) +
        geom_line(data = e_by, aes(x = x, y = response, color = .data[[bv]]), linewidth = 0.9) +
        geom_point(data = dm2, aes(x = x, y = md, color = .data[[bv]]), alpha = 0.5, size = 1.8) +
        scale_color_manual(values = bi$cols) + scale_fill_manual(values = bi$cols) +
        labs(x = xl, y = y_lab_area, title = paste0(cell_label, " \u2014 ", tag, " \u00D7 ", bv, " (T1D)")) +
        theme_pub + theme(legend.position = c(0.85, 0.85), legend.background = element_rect(fill = "white", color = NA))
      save_plot(plt2, paste0(tag, "_by_", tolower(bv)), output_dir)
    }
  }

  # --- 3-way: DD × Region × AO ---
  cat("      3-way interactions...\n")
  for (focal in list(list(var = "Disease.Duration_c", raw = "Disease.Duration",
                           lab = "Disease duration (years)", cond = "age_at_onset_c", tag = "dd_region_ao"),
                      list(var = "age_at_onset_c", raw = "age_at_onset",
                           lab = "Age at onset (years)", cond = "Disease.Duration_c", tag = "ao_region_dd"))) {
    ar <- range(data[[focal$raw]], na.rm = TRUE); mc <- mean(data[[focal$raw]], na.rm = TRUE)
    cond_sd <- sd(data[[focal$cond]], na.rm = TRUE)
    cv_seq <- seq(ar[1] - mc, ar[2] - mc, length.out = 50)
    preds_list <- list()
    for (lev in c(-1, 1)) {
      at_list <- setNames(list(cv_seq, lev * cond_sd), c(focal$var, focal$cond))
      e3 <- exp_emm(model, as.formula(paste("~", focal$var, "| Region")), data, at = at_list)
      e3$x <- e3[[focal$var]] + mc
      cond_lab <- sub("_c$", "", focal$cond)
      e3$cond_level <- factor(ifelse(lev == -1, paste0(cond_lab, " = mean \u2212 1 SD"),
                                                  paste0(cond_lab, " = mean + 1 SD")))
      preds_list[[length(preds_list) + 1]] <- e3
    }
    pa <- bind_rows(preds_list); ym <- max(pa$hi, na.rm = TRUE) * 1.1
    plt <- ggplot(pa) +
      geom_ribbon(aes(x = x, ymin = lo, ymax = hi, fill = Region), alpha = 0.2) +
      geom_line(aes(x = x, y = response, color = Region), linewidth = 0.9) +
      facet_wrap(~ cond_level) +
      scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
      scale_y_continuous(limits = c(0, ym)) +
      labs(x = focal$lab, y = y_lab_area, title = paste0(cell_label, ": ", focal$tag, " (T1D)")) +
      theme_pub + theme(legend.position = "bottom")
    save_plot(plt, paste0("3way_", focal$tag), output_dir, 10, 5)
  }

  # --- 3-way with Islet (if has_islet) ---
  if (has_islet) {
    for (focal in list(list(var = "Disease.Duration_c", raw = "Disease.Duration",
                             lab = "Disease duration (years)", tag = "dd_region_islet"),
                        list(var = "age_at_onset_c", raw = "age_at_onset",
                             lab = "Age at onset (years)", tag = "ao_region_islet"))) {
      ar <- range(data[[focal$raw]], na.rm = TRUE); mc <- mean(data[[focal$raw]], na.rm = TRUE)
      isd <- sd(data$log_Islet.Cells_c, na.rm = TRUE)
      cv_seq <- seq(ar[1] - mc, ar[2] - mc, length.out = 50)
      preds_list <- list()
      for (lev in c(-1, 1)) {
        at_list <- setNames(list(cv_seq, lev * isd), c(focal$var, "log_Islet.Cells_c"))
        e3 <- exp_emm(model, as.formula(paste("~", focal$var, "| Region")), data, at = at_list)
        e3$x <- e3[[focal$var]] + mc
        e3$cond_level <- factor(ifelse(lev == -1, "Islet = mean \u2212 1 SD", "Islet = mean + 1 SD"))
        preds_list[[length(preds_list) + 1]] <- e3
      }
      pa <- bind_rows(preds_list); ym <- max(pa$hi, na.rm = TRUE) * 1.1
      plt <- ggplot(pa) +
        geom_ribbon(aes(x = x, ymin = lo, ymax = hi, fill = Region), alpha = 0.2) +
        geom_line(aes(x = x, y = response, color = Region), linewidth = 0.9) +
        facet_wrap(~ cond_level) +
        scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
        scale_y_continuous(limits = c(0, ym)) +
        labs(x = focal$lab, y = y_lab_area, title = paste0(cell_label, ": ", focal$tag, " (T1D)")) +
        theme_pub + theme(legend.position = "bottom")
      save_plot(plt, paste0("3way_", focal$tag), output_dir, 10, 5)
    }
  }
}


# T1D combined — reuses t1d subset + adds object_type plots
generate_area_combined_t1d_plots <- function(model, data, cell_label, area_col, output_dir) {
  cat("    Generating combined T1D plots for", cell_label, "...\n")
  generate_area_t1d_plots(model, data, cell_label, area_col, has_islet = FALSE, output_dir)

  # object_type main
  dm_ot <- data %>% group_by(Donor, object_type) %>%
    summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)), .groups = "drop")
  e_ot <- exp_emm(model, ~ object_type, data)
  prs_ot <- as.data.frame(pairs(emmeans(model, ~ object_type, data = data), adjust = "tukey"))
  sig <- prs_ot %>% filter(p.value < 0.05) %>% mutate(label = paste(contrast, p_to_stars(p.value)))
  plt <- ggplot() +
    geom_violin(data = dm_ot, aes(x = object_type, y = md, fill = object_type),
                alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm_ot, aes(x = object_type, y = md, color = object_type),
                width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_brewer(palette = "Set2", guide = "none") + scale_color_brewer(palette = "Set2", guide = "none") +
    { if (nrow(sig) > 0) annotate("text", x = 2, y = max(dm_ot$md) * 0.95,
                                   label = paste(sig$label, collapse = "\n"), size = 4, lineheight = 0.9) } +
    labs(x = "Object type", y = y_lab_area, title = paste(cell_label, "\u2014 Object type (T1D)")) + theme_pub
  save_plot(plt, "combined_object_type", output_dir)

  # object_type × DD slope (faceted by object type)
  for (cs in list(list(var = "Disease.Duration_c", raw = "Disease.Duration",
                        lab = "Disease duration (years)", tag = "dd"),
                   list(var = "age_at_onset_c", raw = "age_at_onset",
                        lab = "Age at onset (years)", tag = "ao"))) {
    cv <- cs$var; rv <- cs$raw; xl <- cs$lab; tag <- cs$tag
    ar <- range(data[[rv]], na.rm = TRUE); mc <- mean(data[[rv]], na.rm = TRUE)
    cv_seq <- seq(ar[1] - mc, ar[2] - mc, length.out = 50)
    at_list <- setNames(list(cv_seq), cv)

    e_by_ot <- exp_emm(model, as.formula(paste("~", cv, "| object_type")), data, at = at_list)
    e_by_ot$x <- e_by_ot[[cv]] + mc
    dm_ot_cont <- data %>% group_by(Donor, object_type) %>%
      summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)),
                x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
    trends <- as.data.frame(test(emt_fn(model, ~ object_type, cv, data)))
    trend_labels <- tibble(
      object_type = factor(trends$object_type, levels = levels(data$object_type)),
      label = sapply(trends$p.value, p_to_stars))

    plt <- ggplot() +
      geom_ribbon(data = e_by_ot, aes(x = x, ymin = lo, ymax = hi),
                  fill = col_t1d, alpha = 0.2) +
      geom_line(data = e_by_ot, aes(x = x, y = response), linewidth = 0.9, color = "black") +
      geom_point(data = dm_ot_cont, aes(x = x, y = md), color = col_t1d, alpha = 0.6, size = 1.8) +
      facet_wrap(~ object_type, ncol = 3) +
      geom_text(data = trend_labels, aes(x = ar[1] + diff(ar) * 0.1,
                y = max(dm_ot_cont$md) * 0.95, label = label),
                size = 6, hjust = 0) +
      labs(x = xl, y = y_lab_area,
           title = paste(cell_label, "\u2014", tag, "\u00D7 Object type (T1D)")) + theme_pub
    save_plot(plt, paste0("combined_", tag, "_by_objecttype"), output_dir, 12, 5)
  }

  # object_type × Region
  dm_or <- data %>% group_by(Donor, object_type, Region) %>%
    summarise(md = exp(mean(log(.data[[area_col]]), na.rm = TRUE)), .groups = "drop")
  e_or <- exp_emm(model, ~ Region | object_type, data)
  plt <- ggplot() +
    geom_violin(data = dm_or, aes(x = Region, y = md, fill = Region), alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm_or, aes(x = Region, y = md, color = Region), width = 0.15, size = 1.2, alpha = 0.5) +
    facet_wrap(~ object_type, ncol = 3) +
    scale_fill_manual(values = colors_region, guide = "none") + scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = y_lab_area,
         title = paste(cell_label, "\u2014 Region \u00D7 Object type (T1D)")) + theme_pub
  save_plot(plt, "combined_region_by_objecttype", output_dir, 12, 5)
}


# #############################################################################
#
#  MASTER EXECUTION
#
# #############################################################################

generate_all_area_plots <- function(model_dir = "Area/Models/Total_Area",
                                     output_dir = "Area/Graphs/total_area") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL AREA POST-HOC COMPARISON PLOTS\n")
  cat("=============================================================================\n\n")

  # --- Data loading helper ---
  split_area <- function(data, cell_col) {
    area_to_count <- c(ins_area = "Count.ins", glu_area = "Count.glu",
                       soma_area = "Count.soma", pp_area = "Count.PP")
    log_col <- paste0("log_", cell_col)
    count_col <- area_to_count[[cell_col]]
    data %>% filter(!is.na(!!sym(cell_col)) & !!sym(cell_col) > 0 & !!sym(count_col) > 0) %>%
      mutate(!!log_col := log(!!sym(cell_col)))
  }

  add_object_type <- function(df) {
    df %>% mutate(
      object_type = factor(case_when(
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

  # Load base datasets
  cat("Loading data...\n")
  area_singlets    <- set_factors(readRDS("Data/area_data_singlets.rds"))
  area_endobs      <- set_factors(readRDS("Data/area_data_endobs.rds"))
  area_islets      <- set_factors(readRDS("Data/area_data_islets.rds"))
  area_singlets_t1d <- set_factors_t1d(readRDS("Data/area_data_singlets_t1d.rds"))
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
    beta  = list(area_col = "ins_area",  cell_col = "ins_area",  label = "β cell area"),
    alpha = list(area_col = "glu_area",  cell_col = "glu_area",  label = "α cell area"),
    delta = list(area_col = "soma_area", cell_col = "soma_area", label = "δ cell area"),
    pp    = list(area_col = "pp_area",   cell_col = "pp_area",   label = "PP cell area")
  )

  combined_data <- list(
    beta = list(diag = beta_combined, t1d = beta_combined_t1d),
    alpha = list(diag = alpha_combined, t1d = alpha_combined_t1d),
    delta = list(diag = delta_combined, t1d = delta_combined_t1d),
    pp = list(diag = pp_combined, t1d = pp_combined_t1d)
  )

  combined_model_names <- list(
    beta = "beta_combined_model", alpha = "alpha_combined_model",
    delta = "delta_combined_model", pp = "pp_combined_model"
  )
  combined_t1d_names <- list(
    beta = "beta_combined_t1d", alpha = "alpha_combined_t1d",
    delta = "delta_combined_t1d", pp = "pp_combined_t1d"
  )

  for (ct in names(cell_types)) {
    spec <- cell_types[[ct]]; ac <- spec$area_col; cc <- spec$cell_col; lab <- spec$label
    cat("===", toupper(lab), "===\n")

    # --- DIAGNOSIS: 3 subsets ---
    for (sz in list(list(tag = "1cell", display = "Singlets", data = split_area(area_singlets, cc),
                          model_file = paste0("area_", ct, "_1cell.rds"), has_islet = FALSE),
                     list(tag = "2to14", display = "SEOs", data = split_area(area_endobs, cc),
                          model_file = paste0("area_", ct, "_2to14.rds"), has_islet = TRUE),
                     list(tag = "15plus", display = "Islets", data = split_area(area_islets, cc),
                          model_file = paste0("area_", ct, "_15plus.rds"), has_islet = TRUE))) {
      cat("  Diag", sz$tag, "...\n")
      m <- readRDS(file.path(model_dir, sz$model_file))
      generate_area_diag_plots(m, sz$data, paste(lab, paste0("(", sz$display, ")")), ac, sz$has_islet,
                                file.path(output_dir, ct, "Diagnosis", sz$tag))
      rm(m); gc()
    }

    # --- DIAGNOSIS: Combined ---
    cat("  Diag combined...\n")
    m <- readRDS(file.path(model_dir, paste0(combined_model_names[[ct]], ".rds")))
    generate_area_combined_diag_plots(m, split_area(combined_data[[ct]]$diag, cc),
                                       paste(lab, "(combined)"), ac,
                                       file.path(output_dir, ct, "Diagnosis", "combined"))
    rm(m); gc()

    # --- T1D: 3 subsets ---
    for (sz in list(list(tag = "1cell", display = "Singlets", data = split_area(area_singlets_t1d, cc),
                          model_file = paste0("area_", ct, "_1cell_t1d.rds"), has_islet = FALSE),
                     list(tag = "2to14", display = "SEOs", data = split_area(area_endobs_t1d, cc),
                          model_file = paste0("area_", ct, "_2to14_t1d.rds"), has_islet = TRUE),
                     list(tag = "15plus", display = "Islets", data = split_area(area_islets_t1d, cc),
                          model_file = paste0("area_", ct, "_15plus_t1d.rds"), has_islet = TRUE))) {
      cat("  T1D", sz$tag, "...\n")
      m <- readRDS(file.path(model_dir, sz$model_file))
      generate_area_t1d_plots(m, sz$data, paste(lab, paste0("(", sz$display, ")")), ac, sz$has_islet,
                               file.path(output_dir, ct, "T1D", sz$tag))
      rm(m); gc()
    }

    # --- T1D: Combined ---
    cat("  T1D combined...\n")
    m <- readRDS(file.path(model_dir, paste0(combined_t1d_names[[ct]], ".rds")))
    generate_area_combined_t1d_plots(m, split_area(combined_data[[ct]]$t1d, cc),
                                      paste(lab, "(combined)"), ac,
                                      file.path(output_dir, ct, "T1D", "combined"))
    rm(m); gc()
  }

  cat("\n=============================================================================\n")
  cat("COMPLETE! ~488 plots saved to", output_dir, "\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all area plots:\n")
  cat("  generate_all_area_plots()\n")
}
