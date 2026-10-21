# plot_size_threshold_figure.R — Supplementary Figure S4: the 15-cell islet threshold.
#
#   (A) size distribution of all endocrine objects, split at the 15-cell cutoff
#   (B) size distribution of the included islets (>= 15 cells)
#   (C) how much one reclassified cell changes an object's composition (100/n %)
#   (D) numbers of objects excluded and included by the cutoff
#
# Usage, from the repository root:
#   source("Threshold/R/plot_size_threshold_figure.R")
#   s4 <- generate_size_threshold_figure(readRDS("Data/no_cutoff_data.rds"))
#
# Counts are of unique endocrine objects (Donor, ImageID, Object.Id), so a file
# exported twice cannot count its objects twice. Objects larger than 250 cells
# are counted in (A), (B) and (D) but lie beyond the plotted x-axis.
# Requires: dplyr, ggplot2 (>= 3.5), patchwork, scales.
library(dplyr)
library(ggplot2)
library(patchwork)

S4_CUTOFF    <- 15
S4_BIN_WIDTH <- 5
s4_col_excluded <- "#D55E00"
s4_col_included <- "#0072B2"
s4_groups <- c("< 15 cells (excluded)", "\u2265 15 cells (included)")
s4_two_lines <- function(x) sub(" (", "\n(", x, fixed = TRUE)

s4_theme <- theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title    = element_text(size = 12, face = "bold"),
    axis.text     = element_text(size = 10, color = "black"),
    axis.line     = element_line(linewidth = 0.5, color = "black"),
    axis.ticks    = element_line(linewidth = 0.4, color = "black"),
    plot.title    = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0),
    plot.tag      = element_text(size = 14, face = "bold"),
    plot.margin   = margin(10, 15, 10, 10)
  )

s4_group <- function(cells) {
  factor(ifelse(cells < S4_CUTOFF, s4_groups[1], s4_groups[2]), levels = s4_groups)
}

# Objects per 5-cell bin (1-4, 5-9, 10-14, 15-19, ...), so no bin straddles the cutoff.
size_bins <- function(objects) {
  objects %>%
    mutate(bin_start = floor(Islet.Cells / S4_BIN_WIDTH) * S4_BIN_WIDTH) %>%
    count(bin_start, name = "objects") %>%
    mutate(bin_end = bin_start + S4_BIN_WIDTH - 1,
           group   = s4_group(bin_start)) %>%
    select(bin_start, bin_end, objects, group)
}

s4_size_axes <- list(
  scale_x_continuous(breaks = seq(0, 250, 25), expand = expansion(mult = c(0, 0.02))),
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.08))),
  coord_cartesian(xlim = c(0, 250))
)

plot_size_distribution <- function(bins) {
  ggplot(bins, aes(x = bin_start + S4_BIN_WIDTH / 2, y = objects, fill = group)) +
    geom_col(width = S4_BIN_WIDTH, color = "white", linewidth = 0.2) +
    geom_vline(xintercept = S4_CUTOFF, linetype = "dashed", linewidth = 0.7, color = "black") +
    annotate("text", x = S4_CUTOFF + 2, y = Inf, label = "Cutoff = 15 cells",
             hjust = 0, vjust = 1.5, size = 3.5, fontface = "italic") +
    scale_fill_manual(values = setNames(c(s4_col_excluded, s4_col_included), s4_groups),
                      labels = s4_two_lines, name = NULL) +
    s4_size_axes +
    labs(x = "Islet cell count", y = "Frequency",
         title = "Distribution of endocrine object size", tag = "A") +
    s4_theme +
    theme(legend.position = "inside", legend.position.inside = c(0.75, 0.75),
          legend.background = element_rect(fill = "white", color = NA),
          legend.text = element_text(size = 9))
}

plot_included_sizes <- function(bins) {
  ggplot(filter(bins, bin_start >= S4_CUTOFF), aes(x = bin_start + S4_BIN_WIDTH / 2, y = objects)) +
    geom_col(width = S4_BIN_WIDTH, fill = s4_col_included, color = "white", linewidth = 0.2) +
    s4_size_axes +
    labs(x = "Islet cell count", y = "Frequency",
         title = "Included islets (\u2265 15 cells)", tag = "B") +
    s4_theme
}

plot_reclassification_sensitivity <- function() {
  sens <- tibble(islet_size = c(5, 10, 15, 20, 30, 50, 75, 100, 150),
                 pct_change = 100 / islet_size)
  # Values beside the points they belong to; the 10-cell label goes on the left
  # so it stays clear of the cutoff line.
  labs_df <- filter(sens, islet_size <= S4_CUTOFF) %>%
    mutate(label   = paste0(signif(pct_change, 2), "%"),
           left    = islet_size == 10,
           label_x = islet_size + ifelse(left, -2, 2),
           hjust   = ifelse(left, 1, 0))
  ggplot(sens, aes(x = islet_size, y = pct_change)) +
    geom_line(linewidth = 0.8, color = "grey30") +
    geom_point(aes(color = islet_size < S4_CUTOFF), size = 3) +
    geom_vline(xintercept = S4_CUTOFF, linetype = "dashed", linewidth = 0.7, color = "black") +
    geom_text(data = labs_df,
              aes(x = label_x, label = label, hjust = hjust, color = islet_size < S4_CUTOFF),
              size = 3.2, fontface = "bold") +
    scale_color_manual(values = c("TRUE" = s4_col_excluded, "FALSE" = s4_col_included), guide = "none") +
    scale_x_continuous(breaks = sens$islet_size) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x = "Islet cell count", y = "Composition change\nfrom 1 cell reclassification",
         title = "Proportional sensitivity by islet size", tag = "C") +
    s4_theme
}

plot_cutoff_summary <- function(objects) {
  summary_df <- objects %>%
    count(group = s4_group(Islet.Cells), name = "objects") %>%
    mutate(pct = round(100 * objects / sum(objects), 1))
  ggplot(summary_df, aes(x = group, y = objects, fill = group)) +
    geom_col(width = 0.6, color = "white", linewidth = 0.3) +
    geom_text(aes(label = paste0(scales::comma(objects), "\n(", pct, "%)")),
              vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = setNames(c(s4_col_excluded, s4_col_included), s4_groups), guide = "none") +
    scale_x_discrete(labels = s4_two_lines) +
    scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.2))) +
    labs(x = NULL, y = "Number of endocrine objects",
         title = "Objects excluded by 15-cell cutoff", tag = "D") +
    s4_theme
}

generate_size_threshold_figure <- function(objects, save = TRUE, output_dir = "Threshold/Graphs") {
  objects <- distinct(objects, Donor, ImageID, Object.Id, .keep_all = TRUE)
  bins    <- size_bins(objects)
  panels  <- list(a = plot_size_distribution(bins),
                  b = plot_included_sizes(bins),
                  c = plot_reclassification_sensitivity(),
                  d = plot_cutoff_summary(objects))
  fig <- (panels$a | panels$b) / (panels$c | panels$d)

  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(output_dir, "Supplementary_Figure_S4.pdf"), fig, width = 12, height = 9)
    ggsave(file.path(output_dir, "Supplementary_Figure_S4.png"), fig, width = 12, height = 9, dpi = 300)
  }

  n_ex <- sum(objects$Islet.Cells < S4_CUTOFF)
  cat(sprintf("Endocrine objects: %s | excluded (< %d cells): %s (%.1f%%) | included: %s (%.1f%%)\n",
              scales::comma(nrow(objects)), S4_CUTOFF,
              scales::comma(n_ex), 100 * n_ex / nrow(objects),
              scales::comma(nrow(objects) - n_ex), 100 * (1 - n_ex / nrow(objects))))
  invisible(list(fig = fig, panels = panels))
}
