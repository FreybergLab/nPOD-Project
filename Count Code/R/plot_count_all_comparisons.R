# =============================================================================
# plot_count_all_comparisons.R
# All possible post-hoc comparison plots for Count models (glmmTMB)
# =============================================================================
#
# DIAGNOSIS MODELS (5):
#   Cell counts (ins, glu, soma, pp): NB1, log link
#     (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
#     Diagnosis:Region:log_Islet.Cells_c + Diagnosis:Region:Age_c
#   Islet.Cells: NB2, log link (NO log_Islet.Cells_c)
#     (Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c
#
#   Per cell count model (21 plots):
#     Main effects (3), two-way categorical (6), Age slopes (4),
#     Islet slopes (4), 3-way Age (2), 3-way Islet (2)
#   Per Islet.Cells model (15 plots):
#     Main effects (3), two-way categorical (6), Age slopes (4), 3-way Age (2)
#
# T1D MODELS (5, AO parameterization — DD slopes extracted):
#   Cell counts: (DD_c + Region + Sex + AO_c + log_IC_c)^2 + DD:Region:logIC + DD:Region:AO
#   Islet.Cells: (DD_c + Region + Sex + AO_c)^2 + DD:Region:AO
#
#   Per cell count T1D (19 plots):
#     Categorical (4), DD slopes (3), AO slopes (3), Islet slopes (3),
#     3-way DD×Region×AO (2), 3-way DD×Region×Islet (2), 3-way AO×Region×Islet (2)
#   Per Islet.Cells T1D (12 plots):
#     Categorical (4), DD slopes (3), AO slopes (3), 3-way DD×Region×AO (2)
#
# Total: 4×21 + 15 + 4×19 + 12 = 187 plots
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
    axis.title       = element_text(size = 11, face = "bold"),
    axis.text        = element_text(size = 9, color = "black"),
    axis.line        = element_line(linewidth = 0.4),
    axis.ticks       = element_line(linewidth = 0.3),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 10, color = "gray40"),
    plot.tag         = element_text(size = 18, face = "bold"),
    plot.margin      = margin(8, 12, 8, 8),
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 10),
    strip.text       = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )

p_to_stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

# =============================================================================
# HELPERS — emmeans on log link, exponentiate for response scale
# =============================================================================

emm <- function(model, specs, data, ...) emmeans(model, specs, data = data, ...)
emt <- function(model, specs, var, data, ...) emtrends(model, specs, var = var, data = data, ...)

# Marginalized predictions (works for any log-link glmmTMB)
get_marg_preds <- function(model, data, focal_grid, avg_vars = NULL) {
  model_vars <- setdiff(all.vars(formula(model)[-2]), c("Donor", "ImageID"))
  focal_vars <- names(focal_grid)
  if (is.null(avg_vars))
    avg_vars <- setdiff(intersect(model_vars, c("Diagnosis", "Region", "Sex")), focal_vars)
  avg_levels <- list()
  for (v in avg_vars) avg_levels[[v]] <- unique(data[[v]])
  avg_grid <- if (length(avg_levels) > 0) expand.grid(avg_levels, stringsAsFactors = FALSE) else data.frame(.d = 1)
  all_preds <- all_se <- list()
  for (i in seq_len(nrow(avg_grid))) {
    fg <- focal_grid
    if (length(avg_vars) > 0) for (v in avg_vars) fg[[v]] <- avg_grid[[v]][i]
    for (v in model_vars) if (!v %in% names(fg)) fg[[v]] <- if (is.numeric(data[[v]])) 0 else levels(data[[v]])[1]
    for (v in names(fg)) if (is.factor(data[[v]])) fg[[v]] <- factor(fg[[v]], levels = levels(data[[v]]))
    fg$Donor <- NA
    p <- predict(model, newdata = fg, type = "link", re.form = NA, se.fit = TRUE, allow.new.levels = TRUE)
    all_preds[[i]] <- p$fit; all_se[[i]] <- p$se.fit
  }
  ap <- colMeans(do.call(rbind, all_preds)); as <- sqrt(colMeans(do.call(rbind, all_se)^2))
  focal_grid$estimate <- exp(ap); focal_grid$lower <- exp(ap - 1.96 * as); focal_grid$upper <- exp(ap + 1.96 * as)
  focal_grid
}

save_plot <- function(plt, filename, output_dir, width = 7, height = 5) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("pdf", "png", "tiff"))
    ggsave(file.path(output_dir, paste0(filename, ".", ext)), plt,
           width = width, height = height, dpi = 300,
           bg = if (ext != "pdf") "white" else NULL)
  cat("    Saved:", filename, "\n")
}

fmt_sig <- function(pairs_df) {
  pairs_df %>% filter(p.value < 0.05) %>%
    mutate(label = paste(contrast, p_to_stars(p.value)))
}


# #############################################################################
#
#  DIAGNOSIS MODEL PLOTS
#
# #############################################################################

generate_diag_plots <- function(model, data, cell_label, response_col,
                                 has_islet, output_dir) {
  y_lab <- paste(cell_label, "count")

  # Donor means helper
  dm_cat <- function(...) {
    data %>% group_by(Donor, ...) %>%
      summarise(md = mean(.data[[response_col]], na.rm = TRUE), .groups = "drop")
  }

  # --- Main effects ---
  for (fac in c("Diagnosis", "Region", "Sex")) {
    dm <- dm_cat(!!sym(fac))
    e <- as.data.frame(emmeans(model, as.formula(paste("~", fac)), type = "response", data = data))
    if (fac == "Region") {
      prs <- as.data.frame(pairs(emmeans(model, ~ Region, data = data), adjust = "tukey"))
      sig <- fmt_sig(prs); sig_text <- paste(sig$label, collapse = "\n")
    } else {
      prs <- as.data.frame(pairs(emmeans(model, as.formula(paste("~", fac)), data = data),
                                  adjust = "none", reverse = (fac == "Diagnosis")))
      sig_text <- p_to_stars(prs$p.value[1])
    }
    cols <- if (fac == "Diagnosis") colors_diag else if (fac == "Region") colors_region else colors_sex

    plt <- ggplot() +
      geom_violin(data = dm, aes(x = .data[[fac]], y = md, fill = .data[[fac]]),
                  alpha = 0.3, color = NA, width = 0.7) +
      geom_jitter(data = dm, aes(x = .data[[fac]], y = md, color = .data[[fac]]),
                  width = 0.15, size = 1.5, alpha = 0.6) +
      geom_pointrange(data = e, aes(x = .data[[fac]], y = response, ymin = asymp.LCL, ymax = asymp.UCL),
                      size = 0.5, linewidth = 0.8, color = "black") +
      scale_fill_manual(values = cols, guide = "none") +
      scale_color_manual(values = cols, guide = "none") +
      labs(x = if (fac == "Region") "Pancreatic region" else NULL, y = y_lab,
           title = paste(cell_label, "\u2014", fac)) + theme_pub

    if (nchar(sig_text) > 0) {
      xpos <- if (fac == "Region") 2 else 1.5
      plt <- plt + annotate("text", x = xpos, y = max(dm$md) * 0.95,
                             label = sig_text, size = if (fac == "Region") 4 else 8,
                             lineheight = 0.9)
    }
    save_plot(plt, paste0("main_", tolower(fac)), output_dir)
  }

  # --- Two-way categorical (all 6 views) ---
  two_way_specs <- list(
    list(fac = "Diagnosis", cond = "Region", rev = TRUE, adj_cond = "none", adj_within = "tukey"),
    list(fac = "Region", cond = "Diagnosis", rev = FALSE, adj_cond = "tukey", adj_within = "none"),
    list(fac = "Diagnosis", cond = "Sex", rev = TRUE, adj_cond = "none", adj_within = "none"),
    list(fac = "Sex", cond = "Diagnosis", rev = FALSE, adj_cond = "none", adj_within = "none"),
    list(fac = "Region", cond = "Sex", rev = FALSE, adj_cond = "tukey", adj_within = "none"),
    list(fac = "Sex", cond = "Region", rev = FALSE, adj_cond = "none", adj_within = "none")
  )

  for (spec in two_way_specs) {
    fac <- spec$fac; cond <- spec$cond
    dm <- dm_cat(!!sym(fac), !!sym(cond))
    e <- as.data.frame(emmeans(model, as.formula(paste("~", fac, "|", cond)),
                                type = "response", data = data))
    cols_fac <- if (fac == "Diagnosis") colors_diag else if (fac == "Region") colors_region else colors_sex
    n_cond <- length(unique(data[[cond]]))

    n_fac_levels <- length(unique(data[[fac]]))

    if (n_fac_levels <= 2) {
      # Dodged view
      prs <- as.data.frame(pairs(emmeans(model, as.formula(paste("~", fac, "|", cond)),
                                          data = data), adjust = spec$adj_cond, reverse = spec$rev))
      stars <- sapply(prs$p.value, p_to_stars)
      y_max <- max(dm$md, na.rm = TRUE)

      plt <- ggplot() +
        geom_violin(data = dm, aes(x = .data[[cond]], y = md, fill = .data[[fac]]),
                    alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
        geom_jitter(data = dm, aes(x = .data[[cond]], y = md, color = .data[[fac]], group = .data[[fac]]),
                    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                    size = 1.2, alpha = 0.5) +
        geom_pointrange(data = e, aes(x = .data[[cond]], y = response, ymin = asymp.LCL, ymax = asymp.UCL,
                                       group = .data[[fac]]),
                        position = position_dodge(0.8), size = 0.4, linewidth = 0.7, color = "black") +
        scale_fill_manual(values = cols_fac) + scale_color_manual(values = cols_fac) +
        annotate("text", x = seq_len(n_cond),
                 y = y_max * seq(0.98, by = -0.05, length.out = n_cond),
                 label = stars, size = 8) +
        labs(x = if (cond == "Region") "Pancreatic region" else NULL, y = y_lab,
             title = paste(cell_label, "\u2014", fac, "\u00D7", cond),
             subtitle = paste(fac, "contrast within each", cond),
             fill = fac, color = fac) +
        theme_pub + theme(legend.position = c(0.90, 0.85),
                           legend.background = element_rect(fill = "white", color = NA))
    } else {
      # Faceted view (Region within Diagnosis or Sex)
      prs <- as.data.frame(pairs(emmeans(model, as.formula(paste("~", fac, "|", cond)),
                                          data = data), adjust = spec$adj_cond))
      sig_labels <- prs %>% filter(p.value < 0.05) %>%
        mutate(label = paste(contrast, p_to_stars(p.value))) %>%
        group_by(!!sym(cond)) %>% summarise(label = paste(label, collapse = "\n"), .groups = "drop")
      all_cond <- tibble(!!cond := unique(data[[cond]]))
      sig_labels <- left_join(all_cond, sig_labels, by = cond) %>% mutate(label = replace_na(label, ""))
      y_max <- max(dm$md, na.rm = TRUE) * 1.05

      plt <- ggplot() +
        geom_violin(data = dm, aes(x = .data[[fac]], y = md, fill = .data[[fac]]),
                    alpha = 0.3, color = NA, width = 0.7) +
        geom_jitter(data = dm, aes(x = .data[[fac]], y = md, color = .data[[fac]]),
                    width = 0.15, size = 1.2, alpha = 0.5) +
        geom_pointrange(data = e, aes(x = .data[[fac]], y = response, ymin = asymp.LCL, ymax = asymp.UCL),
                        size = 0.4, linewidth = 0.7, color = "black") +
        facet_wrap(as.formula(paste("~", cond))) +
        geom_text(data = sig_labels, aes(x = 2, y = y_max * 0.95, label = label),
                  size = 3.5, lineheight = 0.9) +
        scale_fill_manual(values = cols_fac, guide = "none") +
        scale_color_manual(values = cols_fac, guide = "none") +
        labs(x = if (fac == "Region") "Pancreatic region" else NULL, y = y_lab,
             title = paste(cell_label, "\u2014", fac, "\u00D7", cond),
             subtitle = paste(fac, "contrasts within each", cond)) + theme_pub
    }
    w <- if (n_cond == 3 || fac == "Region") 9 else 7
    save_plot(plt, paste0(tolower(fac), "_by_", tolower(cond)), output_dir, w, 5)
  }

  # --- Age slopes ---
  cont_specs <- list(list(var = "Age_c", raw = "Age", lab = "Age (years)"))
  if (has_islet) cont_specs[[2]] <- list(var = "log_Islet.Cells_c", raw = "log_Islet.Cells", lab = "Islet size (log cells)")

  for (cs in cont_specs) {
    cv <- cs$var; rv <- cs$raw; xl <- cs$lab
    cont_name <- if (cv == "Age_c") "age" else "islet"
    ar <- range(data[[rv]], na.rm = TRUE); mc <- mean(data[[rv]], na.rm = TRUE)

    # Overall
    fg <- data.frame(x_c = seq(ar[1] - mc, ar[2] - mc, length.out = 100)); names(fg) <- cv
    fg[[rv]] <- fg[[cv]] + mc
    preds <- get_marg_preds(model, data, fg, avg_vars = c("Diagnosis", "Region", "Sex"))
    preds[[rv]] <- preds[[cv]] + mc
    dm <- data %>% group_by(Donor) %>%
      summarise(md = mean(.data[[response_col]], na.rm = TRUE),
                x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
    trend <- as.data.frame(test(emt(model, ~ 1, cv, data)))
    stars <- p_to_stars(trend$p.value[1])

    plt <- ggplot() +
      geom_ribbon(data = preds, aes(x = .data[[rv]], ymin = lower, ymax = upper),
                  fill = "grey60", alpha = 0.3) +
      geom_line(data = preds, aes(x = .data[[rv]], y = estimate), linewidth = 0.9) +
      geom_point(data = dm, aes(x = x, y = md), alpha = 0.5, size = 1.8) +
      labs(x = xl, y = y_lab, title = paste(cell_label, "\u2014", if (cv == "Age_c") "Age" else "Islet size", "(overall)")) +
      theme_pub
    if (nchar(stars) > 0) plt <- plt + annotate("text", x = ar[1] + 1, y = max(dm$md) * 0.95,
                                                  label = stars, size = 8, hjust = 0)
    save_plot(plt, paste0(cont_name, "_overall"), output_dir)

    # By each categorical
    for (by_info in list(list(v = "Diagnosis", cols = colors_diag),
                          list(v = "Region", cols = colors_region),
                          list(v = "Sex", cols = colors_sex))) {
      bv <- by_info$v
      fg2 <- expand.grid(x_c = seq(ar[1] - mc, ar[2] - mc, length.out = 100),
                           fac = levels(data[[bv]]), stringsAsFactors = FALSE)
      names(fg2) <- c(cv, bv); fg2[[rv]] <- fg2[[cv]] + mc
      avg <- setdiff(c("Diagnosis", "Region", "Sex"), bv)
      preds2 <- get_marg_preds(model, data, fg2, avg_vars = avg)
      preds2[[rv]] <- preds2[[cv]] + mc

      dm2 <- data %>% group_by(Donor, !!sym(bv)) %>%
        summarise(md = mean(.data[[response_col]], na.rm = TRUE),
                  x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")

      plt2 <- ggplot() +
        geom_ribbon(data = preds2, aes(x = .data[[rv]], ymin = lower, ymax = upper,
                                        fill = .data[[bv]]), alpha = 0.2) +
        geom_line(data = preds2, aes(x = .data[[rv]], y = estimate, color = .data[[bv]]),
                  linewidth = 0.9) +
        geom_point(data = dm2, aes(x = x, y = md, color = .data[[bv]]),
                   alpha = 0.5, size = 1.8) +
        scale_color_manual(values = by_info$cols) + scale_fill_manual(values = by_info$cols) +
        labs(x = xl, y = y_lab,
             title = paste0(cell_label, " \u2014 ", bv, " \u00D7 ", if (cv == "Age_c") "Age" else "Islet")) +
        theme_pub +
        theme(legend.position = c(0.85, 0.85),
              legend.background = element_rect(fill = "white", color = NA))
      save_plot(plt2, paste0(cont_name, "_by_", tolower(bv)), output_dir)
    }

    # 3-way: Diagnosis × Region × continuous (2 views)
    fg3 <- expand.grid(x_c = seq(ar[1] - mc, ar[2] - mc, length.out = 80),
                        Diagnosis = factor(c("ND", "T1D"), levels = levels(data$Diagnosis)),
                        Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
                        stringsAsFactors = FALSE)
    names(fg3)[1] <- cv; fg3[[rv]] <- fg3[[cv]] + mc
    preds3 <- get_marg_preds(model, data, fg3, avg_vars = "Sex")
    preds3[[rv]] <- preds3[[cv]] + mc
    dm3 <- data %>% group_by(Donor, Diagnosis, Region) %>%
      summarise(md = mean(.data[[response_col]], na.rm = TRUE),
                x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
    y_max <- max(c(dm3$md, preds3$upper), na.rm = TRUE) * 1.1

    # View A: facet Region, color Diagnosis
    plt_a <- ggplot() +
      geom_ribbon(data = preds3, aes(x = .data[[rv]], ymin = lower, ymax = upper, fill = Diagnosis), alpha = 0.2) +
      geom_line(data = preds3, aes(x = .data[[rv]], y = estimate, color = Diagnosis), linewidth = 0.9) +
      geom_point(data = dm3, aes(x = x, y = md, color = Diagnosis), alpha = 0.4, size = 1.5) +
      facet_wrap(~ Region, ncol = 3) +
      scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
      scale_y_continuous(limits = c(0, y_max)) +
      labs(x = xl, y = y_lab,
           title = paste0(cell_label, ": Diag \u00D7 Region \u00D7 ", if (cv == "Age_c") "Age" else "Islet"),
           subtitle = "Facet by Region, Diagnosis overlaid") +
      theme_pub + theme(legend.position = "bottom")
    save_plot(plt_a, paste0("3way_", cont_name, "_facet_region"), output_dir, 12, 5)

    # View B: facet Diagnosis, color Region
    plt_b <- ggplot() +
      geom_ribbon(data = preds3, aes(x = .data[[rv]], ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
      geom_line(data = preds3, aes(x = .data[[rv]], y = estimate, color = Region), linewidth = 0.9) +
      geom_point(data = dm3, aes(x = x, y = md, color = Region), alpha = 0.4, size = 1.5) +
      facet_wrap(~ Diagnosis, ncol = 2) +
      scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
      scale_y_continuous(limits = c(0, y_max)) +
      labs(x = xl, y = y_lab,
           title = paste0(cell_label, ": Diag \u00D7 Region \u00D7 ", if (cv == "Age_c") "Age" else "Islet"),
           subtitle = "Facet by Diagnosis, Region overlaid") +
      theme_pub + theme(legend.position = "bottom")
    save_plot(plt_b, paste0("3way_", cont_name, "_facet_diagnosis"), output_dir, 10, 5)
  }
}


# #############################################################################
#
#  T1D MODEL PLOTS (AO parameterization — DD slopes extracted)
#
# #############################################################################

generate_t1d_plots <- function(model, data, cell_label, response_col,
                                has_islet, output_dir) {
  y_lab <- paste(cell_label, "count")
  dm_cat <- function(...) {
    data %>% group_by(Donor, ...) %>%
      summarise(md = mean(.data[[response_col]], na.rm = TRUE), .groups = "drop")
  }

  # --- Categorical: Region, Sex, Region×Sex ---
  for (fac in c("Region", "Sex")) {
    dm <- dm_cat(!!sym(fac))
    e <- as.data.frame(emmeans(model, as.formula(paste("~", fac)), type = "response", data = data))
    cols <- if (fac == "Region") colors_region else colors_sex
    if (fac == "Region") {
      prs <- as.data.frame(pairs(emm(model, ~ Region, data), adjust = "tukey"))
      sig <- fmt_sig(prs); sig_text <- paste(sig$label, collapse = "\n")
    } else {
      prs <- as.data.frame(pairs(emm(model, ~ Sex, data), adjust = "none"))
      sig_text <- p_to_stars(prs$p.value[1])
    }
    plt <- ggplot() +
      geom_violin(data = dm, aes(x = .data[[fac]], y = md, fill = .data[[fac]]),
                  alpha = 0.3, color = NA, width = 0.7) +
      geom_jitter(data = dm, aes(x = .data[[fac]], y = md, color = .data[[fac]]),
                  width = 0.15, size = 1.5, alpha = 0.6) +
      geom_pointrange(data = e, aes(x = .data[[fac]], y = response, ymin = asymp.LCL, ymax = asymp.UCL),
                      size = 0.5, linewidth = 0.8, color = "black") +
      scale_fill_manual(values = cols, guide = "none") + scale_color_manual(values = cols, guide = "none") +
      labs(x = if (fac == "Region") "Pancreatic region" else NULL, y = y_lab,
           title = paste(cell_label, "\u2014", fac, "(T1D)")) + theme_pub
    if (nchar(sig_text) > 0)
      plt <- plt + annotate("text", x = if (fac == "Region") 2 else 1.5,
                             y = max(dm$md) * 0.95, label = sig_text,
                             size = if (fac == "Region") 4 else 8, lineheight = 0.9)
    save_plot(plt, paste0("main_", tolower(fac)), output_dir)
  }

  # Region | Sex (faceted)
  dm_rs <- dm_cat(Region, Sex)
  e_rs <- as.data.frame(emmeans(model, ~ Region | Sex, type = "response", data = data))
  prs_rs <- as.data.frame(pairs(emm(model, ~ Region | Sex, data), adjust = "tukey"))
  sig_rs <- prs_rs %>% filter(p.value < 0.05) %>%
    mutate(label = paste(contrast, p_to_stars(p.value))) %>%
    group_by(Sex) %>% summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  all_sex <- tibble(Sex = factor(c("Female", "Male"), levels = levels(data$Sex)))
  sig_rs <- left_join(all_sex, sig_rs, by = "Sex") %>% mutate(label = replace_na(label, ""))

  plt <- ggplot() +
    geom_violin(data = dm_rs, aes(x = Region, y = md, fill = Region), alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(data = dm_rs, aes(x = Region, y = md, color = Region), width = 0.15, size = 1.2, alpha = 0.5) +
    geom_pointrange(data = e_rs, aes(x = Region, y = response, ymin = asymp.LCL, ymax = asymp.UCL),
                    size = 0.4, linewidth = 0.7, color = "black") +
    facet_wrap(~ Sex) +
    geom_text(data = sig_rs, aes(x = 2, y = max(dm_rs$md) * 0.95, label = label), size = 3.5, lineheight = 0.9) +
    scale_fill_manual(values = colors_region, guide = "none") + scale_color_manual(values = colors_region, guide = "none") +
    labs(x = "Pancreatic region", y = y_lab, title = paste(cell_label, "\u2014 Region \u00D7 Sex (T1D)")) + theme_pub
  save_plot(plt, "region_by_sex", output_dir, 9, 5)

  # Sex | Region (dodged)
  e_sr <- as.data.frame(emmeans(model, ~ Sex | Region, type = "response", data = data))
  prs_sr <- as.data.frame(pairs(emm(model, ~ Sex | Region, data), adjust = "none"))
  plt <- ggplot() +
    geom_violin(data = dm_rs, aes(x = Region, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm_rs, aes(x = Region, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8), size = 1.2, alpha = 0.5) +
    geom_pointrange(data = e_sr, aes(x = Region, y = response, ymin = asymp.LCL, ymax = asymp.UCL, group = Sex),
                    position = position_dodge(0.8), size = 0.4, linewidth = 0.7, color = "black") +
    scale_fill_manual(values = colors_sex) + scale_color_manual(values = colors_sex) +
    annotate("text", x = 1:3, y = max(dm_rs$md) * c(0.98, 0.93, 0.88),
             label = sapply(prs_sr$p.value, p_to_stars), size = 8) +
    labs(x = "Pancreatic region", y = y_lab, title = paste(cell_label, "\u2014 Sex \u00D7 Region (T1D)"),
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

    # Overall
    fg <- data.frame(x_c = seq(ar[1] - mc, ar[2] - mc, length.out = 100)); names(fg) <- cv
    preds <- get_marg_preds(model, data, fg, avg_vars = c("Region", "Sex"))
    preds$x <- preds[[cv]] + mc
    dm <- data %>% group_by(Donor) %>%
      summarise(md = mean(.data[[response_col]], na.rm = TRUE),
                x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")
    trend <- as.data.frame(test(emt(model, ~ 1, cv, data)))
    stars <- p_to_stars(trend$p.value[1])

    plt <- ggplot() +
      geom_ribbon(data = preds, aes(x = x, ymin = lower, ymax = upper), fill = col_t1d, alpha = 0.2) +
      geom_line(data = preds, aes(x = x, y = estimate), linewidth = 0.9, color = "black") +
      geom_point(data = dm, aes(x = x, y = md), color = col_t1d, alpha = 0.6, size = 2) +
      labs(x = xl, y = y_lab, title = paste(cell_label, "\u2014", tag, "overall (T1D)")) + theme_pub
    if (nchar(stars) > 0) plt <- plt + annotate("text", x = ar[1] + 1, y = max(dm$md) * 0.95,
                                                  label = stars, size = 8, hjust = 0)
    save_plot(plt, paste0(tag, "_overall"), output_dir)

    # By Region, by Sex
    for (by_info in list(list(v = "Region", cols = colors_region),
                          list(v = "Sex", cols = colors_sex))) {
      bv <- by_info$v
      fg2 <- expand.grid(x_c = seq(ar[1] - mc, ar[2] - mc, length.out = 100),
                           fac = levels(data[[bv]]), stringsAsFactors = FALSE)
      names(fg2) <- c(cv, bv)
      preds2 <- get_marg_preds(model, data, fg2, avg_vars = setdiff(c("Region", "Sex"), bv))
      preds2$x <- preds2[[cv]] + mc
      dm2 <- data %>% group_by(Donor, !!sym(bv)) %>%
        summarise(md = mean(.data[[response_col]], na.rm = TRUE),
                  x = mean(.data[[rv]], na.rm = TRUE), .groups = "drop")

      plt2 <- ggplot() +
        geom_ribbon(data = preds2, aes(x = x, ymin = lower, ymax = upper, fill = .data[[bv]]), alpha = 0.2) +
        geom_line(data = preds2, aes(x = x, y = estimate, color = .data[[bv]]), linewidth = 0.9) +
        geom_point(data = dm2, aes(x = x, y = md, color = .data[[bv]]), alpha = 0.5, size = 1.8) +
        scale_color_manual(values = by_info$cols) + scale_fill_manual(values = by_info$cols) +
        labs(x = xl, y = y_lab, title = paste0(cell_label, " \u2014 ", tag, " \u00D7 ", bv, " (T1D)")) +
        theme_pub + theme(legend.position = c(0.85, 0.85),
                           legend.background = element_rect(fill = "white", color = NA))
      save_plot(plt2, paste0(tag, "_by_", tolower(bv)), output_dir)
    }
  }

  # --- 3-way: DD × Region × AO (always present) ---
  for (focal in list(list(var = "Disease.Duration_c", raw = "Disease.Duration",
                           lab = "Disease duration (years)", cond = "age_at_onset_c", tag = "dd_region_ao"),
                      list(var = "age_at_onset_c", raw = "age_at_onset",
                           lab = "Age at onset (years)", cond = "Disease.Duration_c", tag = "ao_region_dd"))) {
    ar <- range(data[[focal$raw]], na.rm = TRUE); mc <- mean(data[[focal$raw]], na.rm = TRUE)
    cond_sd <- sd(data[[focal$cond]], na.rm = TRUE)

    preds_list <- list()
    for (lev in c(-1, 1)) {
      fg <- expand.grid(x_c = seq(ar[1] - mc, ar[2] - mc, length.out = 80),
                          Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
                          stringsAsFactors = FALSE)
      names(fg)[1] <- focal$var
      avg <- setdiff(c("Region", "Sex"), "Region")
      fg[[focal$cond]] <- lev * cond_sd
      preds <- get_marg_preds(model, data, fg, avg_vars = "Sex")
      preds$x <- preds[[focal$var]] + mc
      cond_lab <- sub("_c$", "", focal$cond)
      preds$cond_level <- factor(ifelse(lev == -1,
                                         paste0(cond_lab, " = mean \u2212 1 SD"),
                                         paste0(cond_lab, " = mean + 1 SD")))
      preds_list[[length(preds_list) + 1]] <- preds
    }
    preds_all <- bind_rows(preds_list)
    y_max <- max(preds_all$upper, na.rm = TRUE) * 1.1

    plt <- ggplot(preds_all) +
      geom_ribbon(aes(x = x, ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
      geom_line(aes(x = x, y = estimate, color = Region), linewidth = 0.9) +
      facet_wrap(~ cond_level) +
      scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
      scale_y_continuous(limits = c(0, y_max)) +
      labs(x = focal$lab, y = y_lab,
           title = paste0(cell_label, ": ", toupper(focal$tag), " (T1D)")) +
      theme_pub + theme(legend.position = "bottom")
    save_plot(plt, paste0("3way_", focal$tag), output_dir, 10, 5)
  }

  # --- 3-way: DD × Region × Islet + AO × Region × Islet (if has_islet) ---
  if (has_islet) {
    for (focal in list(list(var = "Disease.Duration_c", raw = "Disease.Duration",
                             lab = "Disease duration (years)", tag = "dd_region_islet"),
                        list(var = "age_at_onset_c", raw = "age_at_onset",
                             lab = "Age at onset (years)", tag = "ao_region_islet"))) {
      ar <- range(data[[focal$raw]], na.rm = TRUE); mc <- mean(data[[focal$raw]], na.rm = TRUE)
      islet_sd <- sd(data$log_Islet.Cells_c, na.rm = TRUE)

      preds_list <- list()
      for (lev in c(-1, 1)) {
        fg <- expand.grid(x_c = seq(ar[1] - mc, ar[2] - mc, length.out = 80),
                            Region = factor(c("Head", "Body", "Tail"), levels = levels(data$Region)),
                            stringsAsFactors = FALSE)
        names(fg)[1] <- focal$var
        fg$log_Islet.Cells_c <- lev * islet_sd
        preds <- get_marg_preds(model, data, fg, avg_vars = "Sex")
        preds$x <- preds[[focal$var]] + mc
        preds$cond_level <- factor(ifelse(lev == -1, "Islet = mean \u2212 1 SD", "Islet = mean + 1 SD"))
        preds_list[[length(preds_list) + 1]] <- preds
      }
      preds_all <- bind_rows(preds_list)
      y_max <- max(preds_all$upper, na.rm = TRUE) * 1.1

      plt <- ggplot(preds_all) +
        geom_ribbon(aes(x = x, ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
        geom_line(aes(x = x, y = estimate, color = Region), linewidth = 0.9) +
        facet_wrap(~ cond_level) +
        scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
        scale_y_continuous(limits = c(0, y_max)) +
        labs(x = focal$lab, y = y_lab,
             title = paste0(cell_label, ": ", toupper(focal$tag), " (T1D)")) +
        theme_pub + theme(legend.position = "bottom")
      save_plot(plt, paste0("3way_", focal$tag), output_dir, 10, 5)
    }
  }
}


# #############################################################################
#
#  MASTER EXECUTION
#
# #############################################################################

generate_all_count_plots <- function(model_dir = "Count/Models",
                                      output_dir = "Count/Graphs/All_Comparisons") {

  cat("\n=============================================================================\n")
  cat("GENERATING ALL COUNT POST-HOC COMPARISON PLOTS\n")
  cat("=============================================================================\n\n")

  cat("Loading data...\n")
  quad_data <- readRDS("Data/quad_data.rds") %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
           Region = factor(Region, levels = c("Head", "Body", "Tail")),
           Sex = factor(Sex, levels = c("Female", "Male")),
           Donor = factor(Donor), ImageID = factor(ImageID))
  quad_t1d <- readRDS("Data/quad_data_t1d.rds") %>%
    mutate(Region = factor(Region, levels = c("Head", "Body", "Tail")),
           Sex = factor(Sex, levels = c("Female", "Male")),
           Donor = factor(Donor), ImageID = factor(ImageID)) %>% droplevels()

  cat("  quad_data:", nrow(quad_data), "rows,", n_distinct(quad_data$Donor), "donors\n")
  cat("  quad_t1d: ", nrow(quad_t1d), "rows,", n_distinct(quad_t1d$Donor), "donors\n\n")

  cell_defs <- list(
    list(label = "Insulin",      resp = "Count.ins",   prefix = "count_ins",   has_islet = TRUE),
    list(label = "Glucagon",     resp = "Count.glu",   prefix = "count_glu",   has_islet = TRUE),
    list(label = "Somatostatin", resp = "Count.soma",  prefix = "count_soma",  has_islet = TRUE),
    list(label = "PP",           resp = "Count.PP",    prefix = "count_pp",    has_islet = TRUE),
    list(label = "Islet Cells",  resp = "Islet.Cells", prefix = "count_islet", has_islet = FALSE)
  )

  for (cd in cell_defs) {
    cat("===", toupper(cd$label), "===\n")

    # Diagnosis model
    cat("  --- Diagnosis ---\n")
    m <- readRDS(file.path(model_dir, paste0(cd$prefix, ".rds")))
    generate_diag_plots(m, quad_data, cd$label, cd$resp, cd$has_islet,
                         file.path(output_dir, cd$prefix, "Diagnosis"))
    rm(m); gc()

    # T1D model (AO parameterization)
    cat("  --- T1D ---\n")
    m <- readRDS(file.path(model_dir, paste0(cd$prefix, "_t1d.rds")))
    generate_t1d_plots(m, quad_t1d, cd$label, cd$resp, cd$has_islet,
                        file.path(output_dir, cd$prefix, "T1D"))
    rm(m); gc()
  }

  cat("\n=============================================================================\n")
  cat("COMPLETE! ~187 plots saved to", output_dir, "\n")
  cat("=============================================================================\n")
}

if (interactive()) {
  cat("\nTo generate all count plots:\n")
  cat("  generate_all_count_plots()\n")
}
