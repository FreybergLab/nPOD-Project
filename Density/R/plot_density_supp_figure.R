# =============================================================================
# plot_density_supp_figure.R
# PI-approved minimalist style: simple p-value text, no boxes
# Marginalized predictions + CI ribbons on all continuous panels
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(emmeans)
library(patchwork)

source("colors_master.R")

theme_pub <- theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.margin = margin(8, 12, 8, 8),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_blank()
  )

# --- Marginalized predictions helper ---
get_marginalized_predictions <- function(model, data, focal_grid, avg_vars = NULL) {
  model_vars <- setdiff(all.vars(formula(model)[-2]), c("Donor", "ImageID"))
  focal_vars <- names(focal_grid)
  if (is.null(avg_vars)) avg_vars <- setdiff(intersect(model_vars, c("Diagnosis","Region","Sex")), focal_vars)
  avg_levels <- list(); for (v in avg_vars) avg_levels[[v]] <- if (is.factor(data[[v]])) levels(data[[v]]) else na.omit(unique(data[[v]]))
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

p_to_stars <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

# ── Panel A: Diagnosis × Age ─────────────────────────────────────────────

plot_supp_a <- function(data, density_diag) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(Age_c = seq(ar[1]-ma, ar[2]-ma, length.out = 100),
                    Diagnosis = factor(c("ND","T1D"), levels = levels(data$Diagnosis)),
                    stringsAsFactors = FALSE) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_predictions(density_diag, data, fg, avg_vars = c("Region","Sex"))
  dm <- data %>% group_by(Donor, Diagnosis, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  # Only annotate if significant
  trends <- emtrends(density_diag, pairwise ~ Diagnosis, var = "Age_c")
  p_val <- as.data.frame(trends$contrasts)$p.value[1]
  stars <- p_to_stars(p_val)
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Diagnosis), alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Diagnosis), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Diagnosis), alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors_diag) + scale_fill_manual(values = colors_diag) +
    labs(x = "Age (years)", y = "Islet density\n(islets/mm\u00B2)", tag = "A",
         title = "Diagnosis \u00D7 Age") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  
  if (nchar(stars) > 0) {
    plt <- plt + annotate("text", x = ar[1] + 2, y = max(dm$md) * 0.95,
                          label = stars, size = 8, hjust = 0)
  }
  plt
}

# ── Panel B: Region × Age ────────────────────────────────────────────────

plot_supp_b <- function(data, density_diag) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(Age_c = seq(ar[1]-ma, ar[2]-ma, length.out = 100),
                    Region = factor(c("Head","Body","Tail"), levels = levels(data$Region)),
                    stringsAsFactors = FALSE) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_predictions(density_diag, data, fg, avg_vars = c("Diagnosis","Sex"))
  dm <- data %>% group_by(Donor, Region, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  trends <- emtrends(density_diag, pairwise ~ Region, var = "Age_c")
  pairs_df <- as.data.frame(trends$contrasts)
  sig_pairs <- pairs_df %>%
    filter(p.value < 0.05) %>%
    mutate(t = sprintf("%s %s", contrast, p_to_stars(p.value)))
  sig_text <- paste(sig_pairs$t, collapse = "\n")
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Region), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Region), alpha = 0.4, size = 1.5) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    labs(x = "Age (years)", y = "Islet density\n(islets/mm\u00B2)", tag = "B",
         title = "Region \u00D7 Age") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  
  if (nchar(sig_text) > 0) {
    plt <- plt + annotate("text", x = ar[1] + 2, y = max(dm$md) * 0.95,
                          label = sig_text, size = 4, hjust = 0, lineheight = 0.9)
  }
  plt
}

# ── Panel C: Diagnosis × Region × Age (facet by Diagnosis) ──────────────

plot_supp_c <- function(data, density_diag) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(Age_c = seq(ar[1]-ma, ar[2]-ma, length.out = 100),
                    Diagnosis = factor(c("ND","T1D"), levels = levels(data$Diagnosis)),
                    Region = factor(c("Head","Body","Tail"), levels = levels(data$Region)),
                    stringsAsFactors = FALSE) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_predictions(density_diag, data, fg, avg_vars = "Sex")
  dm <- data %>% group_by(Donor, Diagnosis, Region, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  # Region contrasts within each Diagnosis - only show significant ones as stars
  pairs_df <- as.data.frame(emtrends(density_diag, pairwise ~ Region | Diagnosis, var = "Age_c")$contrasts)
  
  sig_labels <- pairs_df %>%
    filter(p.value < 0.05) %>%
    mutate(label = sprintf("%s %s", contrast, p_to_stars(p.value))) %>%
    group_by(Diagnosis) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")
  
  # Fill in empty labels for non-significant Diagnosis levels
  all_diag <- data.frame(Diagnosis = factor(c("ND","T1D"), levels = levels(data$Diagnosis)))
  sig_labels <- left_join(all_diag, sig_labels, by = "Diagnosis") %>%
    mutate(label = ifelse(is.na(label), "", label))
  
  ym <- max(c(dm$md, preds$upper), na.rm = TRUE) * 1.1
  
  ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Region), alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Region), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Region), alpha = 0.4, size = 1.5) +
    facet_wrap(~ Diagnosis, ncol = 2) +
    geom_text(data = sig_labels, aes(x = ar[1] + 2, y = ym * 0.95, label = label),
              size = 4, hjust = 0, vjust = 1, lineheight = 0.9) +
    scale_color_manual(values = colors_region) + scale_fill_manual(values = colors_region) +
    scale_y_continuous(limits = c(0, ym)) +
    labs(x = "Age (years)", y = "Islet density\n(islets/mm\u00B2)", tag = "C",
         title = "Diagnosis \u00D7 Region \u00D7 Age") +
    theme_pub +
    theme(legend.position = "bottom")
}

# ── Panel D: Sex × Age ───────────────────────────────────────────────────

plot_supp_d <- function(data, density_diag) {
  ar <- range(data$Age, na.rm = TRUE); ma <- mean(data$Age, na.rm = TRUE)
  fg <- expand.grid(Age_c = seq(ar[1]-ma, ar[2]-ma, length.out = 100),
                    Sex = factor(c("Female","Male"), levels = levels(data$Sex)),
                    stringsAsFactors = FALSE) %>% mutate(Age = Age_c + ma)
  preds <- get_marginalized_predictions(density_diag, data, fg, avg_vars = c("Diagnosis","Region"))
  dm <- data %>% group_by(Donor, Sex, Age) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  trends <- emtrends(density_diag, pairwise ~ Sex, var = "Age_c")
  p_val <- as.data.frame(trends$contrasts)$p.value[1]
  stars <- p_to_stars(p_val)
  
  plt <- ggplot() +
    geom_ribbon(data = preds, aes(x = Age, ymin = lower, ymax = upper, fill = Sex), alpha = 0.2) +
    geom_line(data = preds, aes(x = Age, y = estimate, color = Sex), linewidth = 0.9) +
    geom_point(data = dm, aes(x = Age, y = md, color = Sex), alpha = 0.5, size = 1.8) +
    scale_color_manual(values = colors_sex) + scale_fill_manual(values = colors_sex) +
    labs(x = "Age (years)", y = "Islet density\n(islets/mm\u00B2)", tag = "D",
         title = "Sex \u00D7 Age") +
    theme_pub +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
  
  if (nchar(stars) > 0) {
    plt <- plt + annotate("text", x = ar[1] + 2, y = max(dm$md) * 0.95,
                          label = stars, size = 8, hjust = 0)
  }
  plt
}

# ── Panel E: Region × Sex (categorical) ──────────────────────────────────

plot_supp_e <- function(data, density_diag) {
  emm <- as.data.frame(emmeans(density_diag, ~ Sex | Region, type = "response"))
  dm <- data %>% group_by(Donor, Region, Sex) %>%
    summarise(md = mean(islet_density, na.rm = TRUE), .groups = "drop")
  
  ggplot() +
    geom_violin(data = dm, aes(x = Region, y = md, fill = Sex),
                alpha = 0.25, color = NA, position = position_dodge(0.8), width = 0.7) +
    geom_jitter(data = dm, aes(x = Region, y = md, color = Sex, group = Sex),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
                size = 1.2, alpha = 0.5) +
    scale_fill_manual(values = colors_sex) + scale_color_manual(values = colors_sex) +
    labs(x = "Pancreatic region", y = "Islet density\n(islets/mm\u00B2)", tag = "E",
         title = "Region \u00D7 Sex", fill = "Sex", color = "Sex") +
    theme_pub +
    theme(legend.position = c(0.95, 0.85),
          legend.background = element_rect(fill = "white", color = NA))
}

# =============================================================================
# GENERATE
# =============================================================================

generate_supp_figure <- function(Density, density_diag, save = TRUE, output_dir = "Density/Graphs") {
  cat("Generating supplementary figure...\n")
  pa <- plot_supp_a(Density, density_diag); cat("  A\n")
  pb <- plot_supp_b(Density, density_diag); cat("  B\n")
  pc <- plot_supp_c(Density, density_diag); cat("  C\n")
  pd <- plot_supp_d(Density, density_diag); cat("  D\n")
  pe <- plot_supp_e(Density, density_diag); cat("  E\n")
  
  fig <- (pa | pb) / pc / (pd | pe) +
    plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
  
  if (save) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(output_dir, "density_supp_figure.pdf"), fig,
           width = 14, height = 14, dpi = 300)
    ggsave(file.path(output_dir, "density_supp_figure.png"), fig,
           width = 14, height = 14, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, "density_supp_figure.tiff"), fig,
           width = 14, height = 14, dpi = 300, bg = "white")
    cat("Saved to", output_dir, "\n")
  }
  
  list(fig = fig, panels = list(a = pa, b = pb, c = pc, d = pd, e = pe))
}
