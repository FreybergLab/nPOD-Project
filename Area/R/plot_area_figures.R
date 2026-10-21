# plot_area_figures.R — the published area figures, drawn from the data that
# total_area_analysis.Rmd and per_cell_area_analysis.Rmd save for each cell
# type (figure_data.rds), with their final titles, axis labels and layout:
#
#   Figure 8   total area, ND vs T1D
#   Figure S17 per-cell area, ND vs T1D
#   Figure S18 total area vs disease duration (T1D)
#   Figure S19 per-cell area vs disease duration (T1D)
#
# Each figure has a row per cell type (A-D: beta, alpha, delta, PP) and a column
# per object size. Because the data are read from disk, the cell types can be
# run in separate R sessions; this script needs only the four small files.
#
# Usage, from the repository root:
#   source("Area/R/plot_area_figures.R")
#   generate_area_figures("Area/Graphs/total_area", measure = "total")       # Fig 8, S18
#   generate_area_figures("Area/Graphs/Per_Cell_Area", measure = "per_cell") # S17, S19
# Requires: ggplot2, patchwork; theme_pub and the colours come from
# plot_area_all_comparisons.R and colors_master.R.
library(ggplot2)
library(patchwork)
source("Area/R/plot_area_all_comparisons.R")

area_cell_names <- c(beta = "\u03b2-cells", alpha = "\u03b1-cells", delta = "\u03b4-cells", pp = "PP-cells")
area_size_names <- c("1cell" = "Single", "2to14" = "SEOs", "15plus" = "Islets")
area_y_lab      <- "Cell area (\u00b5m\u00b2)"

# ND vs T1D: donor geometric-mean areas, as violins with points.
area_diag_panel <- function(d, title) {
  p <- ggplot(d$points, aes(x = Diagnosis, y = md)) +
    geom_violin(aes(fill = Diagnosis), alpha = 0.3, color = NA, width = 0.7) +
    geom_jitter(aes(color = Diagnosis), width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors_diag, guide = "none") +
    scale_color_manual(values = colors_diag, guide = "none") +
    labs(x = NULL, y = area_y_lab, title = title) + theme_pub
  if (nchar(d$stars) > 0) {
    p <- p + annotate("text", x = 1.5, y = Inf, vjust = 1.3, label = d$stars, size = 8)
  }
  p
}

# T1D: donor geometric-mean areas against disease duration, with the model fit.
area_t1d_panel <- function(d, title) {
  p <- ggplot() +
    geom_ribbon(data = d$curve, aes(x = x, ymin = lo, ymax = hi), fill = col_t1d, alpha = 0.2) +
    geom_line(data = d$curve, aes(x = x, y = response), linewidth = 0.9, color = "black") +
    geom_point(data = d$points, aes(x = x, y = md), color = col_t1d, alpha = 0.6, size = 2) +
    labs(x = "Disease duration (years)", y = area_y_lab, title = title) + theme_pub
  if (nchar(d$stars) > 0) {
    p <- p + annotate("text", x = mean(range(d$curve$x)), y = Inf, vjust = 1.3,
                      label = d$stars, size = 8)
  }
  p
}

# One figure: a row per cell type, tagged A-D, and a column per object size.
assemble_area_figure <- function(fig_data, stage, sizes, title_suffix) {
  draw <- if (stage == "diag") area_diag_panel else area_t1d_panel
  panels <- lapply(setNames(nm = names(area_cell_names)), function(ct) {
    lapply(setNames(nm = sizes), function(sz) {
      draw(fig_data[[ct]][[stage]][[sz]],
           paste0(area_cell_names[[ct]], ": ", area_size_names[[sz]], title_suffix))
    })
  })
  rows <- lapply(panels, function(row) wrap_elements(full = wrap_plots(row, nrow = 1)))
  figure <- wrap_plots(rows, ncol = 1) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 16, face = "bold"))
  list(figure = figure, panels = panels)
}

generate_area_figures <- function(data_dir, measure = c("total", "per_cell"),
                                  save = TRUE, output_dir = "Area/Graphs/Figures") {
  measure <- match.arg(measure)
  spec <- if (measure == "total") {
    list(sizes = c("1cell", "2to14", "15plus"), suffix = "",
         files = c(diag = "Figure_8", t1d = "Supplementary_Figure_S18"))
  } else {
    list(sizes = c("2to14", "15plus"), suffix = " (Per-cell)",
         files = c(diag = "Supplementary_Figure_S17", t1d = "Supplementary_Figure_S19"))
  }
  files <- file.path(data_dir, names(area_cell_names), "figure_data.rds")
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Run these cell types first (no figure data yet): ", paste(missing, collapse = ", "))
  fig_data <- setNames(lapply(files, readRDS), names(area_cell_names))

  set.seed(2026)  # the jittered points land in the same place every run
  out <- list(diag = assemble_area_figure(fig_data, "diag", spec$sizes, spec$suffix),
              t1d  = assemble_area_figure(fig_data, "t1d",  spec$sizes, spec$suffix))
  if (save) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    width <- 3.6 * length(spec$sizes)
    for (k in names(out)) {
      f <- file.path(output_dir, spec$files[[k]])
      # cairo_pdf draws the Greek letters and the superscript; the default pdf() cannot.
      ggsave(paste0(f, ".pdf"), out[[k]]$figure, width = width, height = 13, device = cairo_pdf)
      ggsave(paste0(f, ".tiff"), out[[k]]$figure, width = width, height = 13, dpi = 300,
             compression = "lzw", bg = "white")
      cat("  Saved:", f, "(.pdf and .tiff)\n")
    }
  }
  invisible(list(figures = lapply(out, function(o) o$figure),
                 panels  = lapply(out, function(o) o$panels),
                 data    = fig_data))
}
