# sensitivity_local.R
#
# Data checks and islet-density sensitivity refits for the nPOD manuscript.
# Local run (glmmTMB + emmeans, no Stan). Run from the project root, i.e. the
# directory holding Data/ for all_npod_models_v4.Rmd, after
# Create_count_datasets_v1.Rmd.
#
# Part 1 fits nothing: it prints checks [1]-[5].
# Part 2 refits the published density model, stops unless it reproduces
# Table S2 (T1D - ND = -0.3382 on the log scale), then refits it with exactly
# one change per variant. It then tests density by object-size threshold and
# compares the alpha fraction and composition between Double and Quad runs.
# Part 3 compares cells per islet with the beta-loss counterfactual.
# Part 4 compares segmented single-cell areas between T1D and ND.
# Part 5 tests T1D cells per islet against the counterfactual directly.
#
# Outputs (all in Sensitivity/): checks_*.csv, density_sensitivity.csv,
# density_multicellular_objects.csv, density_by_threshold.csv,
# alpha_double_vs_quad*.csv, composition_by_diagnosis.csv,
# counterfactual_by_donor.csv, cells_per_islet_ratio.csv,
# singlet_cell_area.csv, cells_per_islet_T1D_vs_counterfactual.csv,
# counterfactual_delta_pp_tests.csv

library(tidyverse)
library(glmmTMB)
library(emmeans)

dir.create("Sensitivity", showWarnings = FALSE)

# ImageID is a factor in the count data and numeric in density_data; compare
# them through numeric, as Part C4 does.
img_key <- function(x) as.character(as.numeric(as.character(x)))

no_cutoff <- readRDS("Data/no_cutoff_data.rds") %>%
  mutate(Donor = as.character(Donor), ImageID = img_key(ImageID))
cutoff <- readRDS("Data/cutoff_data.rds") %>%
  mutate(Donor = as.character(Donor), ImageID = img_key(ImageID))
quad <- readRDS("Data/quad_data.rds") %>%
  mutate(Donor = as.character(Donor), ImageID = img_key(ImageID))
density <- readRDS("Data/density_data.rds") %>%
  mutate(Donor = as.character(Donor), ImageID = img_key(ImageID))


# =============================================================================
# Part 1: checks
# =============================================================================

# [1] Classifier runs per image. Object.Id restarts in every run and Part A
# never collapses runs, so an image exported from k runs can contribute the
# same physical islet up to k times to cutoff_data, to the density islet counts
# (C3 counts every cutoff_data row per ImageID) and to ratio_data.
runs <- no_cutoff %>% distinct(Diagnosis, ImageID, source_file, Stain)
runs_per_image <- runs %>%
  group_by(Diagnosis, ImageID) %>%
  summarise(n_runs = n(),
            stains = paste(sort(unique(Stain)), collapse = "+"),
            .groups = "drop")
multi <- runs_per_image$ImageID[runs_per_image$n_runs > 1]

cat("\n[1] Classifier runs per image\n")
print(runs_per_image %>% count(Diagnosis, n_runs, stains), n = Inf)
cat(sprintf("Images with >1 run: %d of %d. Rows from them: %d of %d (no_cutoff), %d of %d (cutoff)\n",
            length(multi), nrow(runs_per_image),
            sum(no_cutoff$ImageID %in% multi), nrow(no_cutoff),
            sum(cutoff$ImageID %in% multi), nrow(cutoff)))
write_csv(runs_per_image, "Sensitivity/checks_runs_per_image.csv")

# [2] Rows per physical object, from the bounding boxes that
# Create_count_datasets_v1 saves. 1.00 means no duplication across runs.
if (file.exists("Data/phenotyper_geometry.rds")) {
  geom <- readRDS("Data/phenotyper_geometry.rds") %>%
    transmute(source_file = as.character(source_file),
              Object.Id   = as.integer(Object.Id),
              islet_uid   = paste(ImageID, XMin, XMax, YMin, YMax, sep = "_")) %>%
    # A file exported twice repeats its geometry rows; keep one per object so
    # the join cannot multiply data rows.
    distinct(source_file, Object.Id, .keep_all = TRUE)
  rows_per_object <- function(d) {
    d %>%
      mutate(source_file = as.character(source_file),
             Object.Id   = as.integer(Object.Id)) %>%
      left_join(geom, by = c("source_file", "Object.Id")) %>%
      summarise(rows             = n(),
                no_geometry      = sum(is.na(islet_uid)),
                physical_objects = n_distinct(islet_uid[!is.na(islet_uid)])) %>%
      mutate(rows_per_object = (rows - no_geometry) / physical_objects)
  }
  ratio_rows <- cutoff %>%
    filter(Stain != "Single", Donor != "6473", Count.ins > 0 | Count.glu > 0)
  cat("\n[2] Rows per physical object (1.00 = no duplication)\n")
  print(bind_rows(no_cutoff = rows_per_object(no_cutoff),
                  cutoff    = rows_per_object(cutoff),
                  ratio     = rows_per_object(ratio_rows),
                  quad      = rows_per_object(quad),
                  .id = "dataset"))
} else {
  cat("\n[2] Data/phenotyper_geometry.rds not found - run Create_count_datasets_v1.Rmd\n")
}

# [3] Density inputs: the classifier run(s) behind each density image, and the
# images Part C drops because no object passed the size filter. On an
# insulin-only run Islet.Cells counts only insulin+ cells, so an
# insulin-deficient islet can never reach 15 cells there.
dens_chk <- density %>%
  left_join(runs_per_image %>% dplyr::select(ImageID, stains), by = "ImageID")
cat(sprintf("\n[3] Density rows: %d; distinct images: %d\n",
            nrow(dens_chk), n_distinct(dens_chk$ImageID)))
print(dens_chk %>%
        count(Diagnosis, Sex, stains) %>%
        pivot_wider(names_from = stains, values_from = n, values_fill = 0))
cat("Images with endocrine objects but no 15-244-cell islet (dropped, not entered as 0):\n")
print(runs_per_image %>% filter(!ImageID %in% density$ImageID) %>% count(Diagnosis, stains))

# [4] Part B4 joins area to counts on (Donor, ImageID, Object.Id). A key shared
# by more than one run can pair an area row with another object's counts.
shared_keys <- no_cutoff %>% count(Donor, ImageID, Object.Id) %>% filter(n > 1)
cat(sprintf("\n[4] (Donor, ImageID, Object.Id) keys shared by >1 run: %d of %d\n",
            nrow(shared_keys), nrow(distinct(no_cutoff, Donor, ImageID, Object.Id))))

# [5] Quad-run objects by size class and diagnosis: R2's islet numbers by
# diagnosis, what the 244-cell cap (log_Islet.Cells <= 5.5) removes, and how
# large those objects are. Also shows whether T1D objects shift below 15 cells.
area_col <- grep("^Region\\.Area", names(no_cutoff), value = TRUE)[1]
size_classes <- no_cutoff %>%
  filter(Stain == "Quad") %>%
  mutate(size_class = case_when(Islet.Cells == 1       ~ "1 single cell",
                                Islet.Cells <= 14      ~ "2 SEO (2-14)",
                                log_Islet.Cells <= 5.5 ~ "3 islet (15-244)",
                                TRUE                   ~ "4 above cap (>244)"),
         eq_diameter_um = if (is.na(area_col)) NA_real_ else 2 * sqrt(.data[[area_col]] / pi)) %>%
  group_by(Diagnosis, size_class) %>%
  summarise(objects            = n(),
            donors             = n_distinct(Donor),
            endocrine_cells    = sum(Islet.Cells),
            median_diameter_um = median(eq_diameter_um, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(Diagnosis) %>%
  mutate(pct_objects         = 100 * objects / sum(objects),
         pct_endocrine_cells = 100 * endocrine_cells / sum(endocrine_cells)) %>%
  ungroup()
cat("\n[5] Quad-run objects by size class\n")
print(size_classes, width = Inf)
write_csv(size_classes, "Sensitivity/checks_size_classes.csv")


# =============================================================================
# Part 2: density refits
# =============================================================================

# Published model (all_npod_models_v4.Rmd, chunk density-diag), verbatim.
FE_DENSITY <- "(Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c"
S2_T1D_ND  <- -0.3382                    # Table S2, Fig 2A, log scale
ATYPICAL   <- c("6405", "6456", "6466")  # pregnant at death; 6 days postpartum; autopsy, severe autolysis

# Islets per image from a single classifier run, preferring Quad > Double >
# Single as the marker pipeline does. A physical-stain lookup would be better.
one_run <- cutoff %>%
  mutate(rank = coalesce(match(Stain, c("Quad", "Double", "Single")), 4L)) %>%
  group_by(ImageID) %>%
  filter(source_file == source_file[which.min(rank)]) %>%
  summarise(islet_count_1run = n(), stain_1run = first(Stain), .groups = "drop")

dens <- density %>% left_join(one_run, by = "ImageID")

prep <- function(d) {
  d %>%
    mutate(Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
           Region    = factor(Region, levels = c("Head", "Body", "Tail")),
           Sex       = factor(Sex, levels = c("Female", "Male")),
           Donor     = factor(Donor),
           Age_c     = scale(Age, center = TRUE, scale = FALSE)[, 1]) %>%
    droplevels()
}

# One change per variant relative to the published fit.
variants <- list(
  baseline     = function(d) d,
  one_run      = function(d) d %>% mutate(islet_density = islet_count_1run / full_tissue_area),
  glu_stained  = function(d) d %>% filter(stain_1run %in% c("Double", "Quad")),
  age_ge5      = function(d) d %>% filter(Age >= 5),
  age_ge13     = function(d) d %>% filter(Age >= 13),          # youngest T1D male is 13.2
  age_log      = function(d) d %>% mutate(Age = log(Age + 1)), # Age_c becomes centred log(age + 1)
  no_hpap      = function(d) d %>% filter(!str_detect(Donor, "^HPAP")),
  no_atypical  = function(d) d %>% filter(!Donor %in% ATYPICAL),
  region_slope = function(d) d
)
RE_OPTIONS <- list(region_slope = c("(1 + Region | Donor)",
                                    "(1 | Donor) + (1 | Donor:Region)"))

fit_density <- function(d, re_options) {
  for (re in re_options) {
    m <- tryCatch(glmmTMB(as.formula(paste("islet_density ~", FE_DENSITY, "+", re)),
                          data = d, family = Gamma(link = "log")),
                  error = function(e) NULL)
    if (!is.null(m) && isTRUE(m$sdr$pdHess)) {
      attr(m, "re") <- re
      return(m)
    }
  }
  stop("no random-effect structure converged: ", paste(re_options, collapse = " | "))
}

# Contrasts as in Table S2: log scale, averaged over the other factors at mean
# Age_c, Holm-adjusted across levels of the conditioning factor. `d` is passed
# to emmeans because the model was fitted inside a function, where emmeans
# cannot always re-evaluate the call's data argument.
density_contrasts <- function(m, d) {
  holm <- function(con) summary(con, infer = TRUE, by = NULL, adjust = "holm")
  bind_rows(
    summary(pairs(emmeans(m, ~ Diagnosis, data = d), reverse = TRUE), infer = TRUE) %>%
      as.data.frame() %>% mutate(within = "all"),
    holm(pairs(emmeans(m, ~ Diagnosis | Sex, data = d), reverse = TRUE)) %>%
      as.data.frame() %>% mutate(within = paste("Sex =", Sex)),
    holm(pairs(emmeans(m, ~ Sex | Diagnosis, data = d))) %>%
      as.data.frame() %>% filter(Diagnosis == "T1D") %>% mutate(within = "Diagnosis = T1D"),
    holm(pairs(emmeans(m, ~ Diagnosis | Region, data = d), reverse = TRUE)) %>%
      as.data.frame() %>% mutate(within = paste("Region =", Region))
  ) %>%
    transmute(contrast = as.character(contrast), within,
              log_estimate = estimate, SE, p.value,
              pct_diff = 100 * (exp(estimate) - 1))
}

cat("\n=== Density refits ===\n")
density_results <- imap_dfr(variants, function(f, name) {
  d   <- dens %>% f() %>% prep()
  m   <- fit_density(d, RE_OPTIONS[[name]] %||% "(1 | Donor)")
  con <- density_contrasts(m, d)
  if (name == "baseline") {
    got <- con$log_estimate[con$within == "all"]
    if (abs(got - S2_T1D_ND) > 0.002) {
      stop(sprintf("Baseline gives T1D - ND = %.4f; Table S2 has %.4f. Resolve this before reading any variant.",
                   got, S2_T1D_ND))
    }
    cat(sprintf("Baseline reproduces Table S2 (T1D - ND = %.4f)\n", got))
  }
  t1d_sex <- d %>% filter(Diagnosis == "T1D") %>% distinct(Donor, Sex) %>% count(Sex)
  cat(sprintf("  %-12s %4d images, %3d donors, RE %s\n",
              name, nrow(d), n_distinct(d$Donor), attr(m, "re")))
  con %>%
    mutate(variant = name, re = attr(m, "re"),
           n_images = nrow(d), n_donors = n_distinct(d$Donor),
           T1D_female_donors = sum(t1d_sex$n[t1d_sex$Sex == "Female"]),
           T1D_male_donors   = sum(t1d_sex$n[t1d_sex$Sex == "Male"]),
           .before = 1)
})
write_csv(density_results, "Sensitivity/density_sensitivity.csv")

# The density claims in the paper, per variant: overall T1D - ND, T1D - ND
# within each sex, and Female - Male within T1D. as_tibble() because
# density_results is a data.frame, and print.data.frame passes n = Inf on to
# print.default, where it partially matches na.print and errors.
density_results %>%
  as_tibble() %>%
  filter(within %in% c("all", "Sex = Female", "Sex = Male", "Diagnosis = T1D")) %>%
  transmute(variant, within, contrast, pct = round(pct_diff, 1), p = signif(p.value, 2)) %>%
  print(n = Inf)

# Rows sharing a source file and Object.Id: an export saved twice.
duplicate_objects <- no_cutoff %>% add_count(source_file, Object.Id) %>% filter(n > 1) %>% count(Donor, Diagnosis, Stain)
write_csv(duplicate_objects, "Sensitivity/checks_duplicate_objects.csv")
print(duplicate_objects)

eo <- no_cutoff %>%
  filter(Islet.Cells >= 2) %>%
  count(ImageID, name = "eo_count") %>%
  inner_join(dens %>% filter(stain_1run %in% c("Double", "Quad")), by = "ImageID") %>%
  mutate(islet_density = eo_count / full_tissue_area) %>%
  prep()
m_eo <- fit_density(eo, "(1 | Donor)")
density_multicellular <- density_contrasts(m_eo, eo) %>% as_tibble() %>% filter(within %in% c("all", "Sex = Female", "Sex = Male"))
write_csv(density_multicellular, "Sensitivity/density_multicellular_objects.csv")
print(density_multicellular)

by_threshold <- map_dfr(c(2, 5, 10, 15, 25), function(k) {
  d <- no_cutoff %>%
    filter(Islet.Cells >= k, log_Islet.Cells <= 5.5) %>%
    count(ImageID, name = "n_obj") %>%
    inner_join(dens %>% filter(stain_1run %in% c("Double", "Quad")), by = "ImageID") %>%
    mutate(islet_density = n_obj / full_tissue_area) %>%
    prep()
  density_contrasts(fit_density(d, "(1 | Donor)"), d) %>%
    as_tibble() %>%
    filter(within %in% c("all", "Sex = Female", "Sex = Male")) %>%
    mutate(min_cells = k, images = nrow(d), .before = 1)
})
write_csv(by_threshold, "Sensitivity/density_by_threshold.csv")
by_threshold %>%
  transmute(min_cells, within, pct = round(pct_diff, 1), p = signif(p.value, 2), images) %>%
  print(n = Inf)

alpha_by_stain <- no_cutoff %>%
  filter(Stain %in% c("Double", "Quad"),
         Count.ins + Count.glu >= 15, Count.ins + Count.glu <= 244) %>%
  group_by(Diagnosis, Donor, Region, Stain) %>%
  summarise(alpha_frac = sum(Count.glu) / sum(Count.ins + Count.glu), .groups = "drop") %>%
  pivot_wider(names_from = Stain, values_from = alpha_frac) %>%
  filter(!is.na(Double), !is.na(Quad)) %>%
  group_by(Diagnosis, Donor) %>%
  summarise(Double = mean(Double), Quad = mean(Quad), .groups = "drop")

write_csv(alpha_by_stain, "Sensitivity/alpha_double_vs_quad_by_donor.csv")
alpha_double_vs_quad <- alpha_by_stain %>%
  group_by(Diagnosis) %>%
  summarise(donors        = n(),
            median_diff   = median(Quad - Double),
            p             = wilcox.test(Quad, Double, paired = TRUE)$p.value,
            double_median = median(Double),
            quad_median   = median(Quad))
write_csv(alpha_double_vs_quad, "Sensitivity/alpha_double_vs_quad.csv")
print(alpha_double_vs_quad)

# Each donor's cell-weighted composition, then the mean and SD across donors.
# Each donor's four fractions sum to 1, so the four means do too.
composition_capped <- no_cutoff %>%
  filter(Stain == "Quad", Islet.Cells >= 15, log_Islet.Cells <= 5.5) %>%
  group_by(Diagnosis, Donor) %>%
  summarise(across(c(Count.ins, Count.glu, Count.soma, Count.PP), sum), .groups = "drop") %>%
  mutate(total = Count.ins + Count.glu + Count.soma + Count.PP,
         across(c(Count.ins, Count.glu, Count.soma, Count.PP), ~ .x / total)) %>%
  group_by(Diagnosis) %>%
  summarise(donors = n(),
            across(c(Count.ins, Count.glu, Count.soma, Count.PP), list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}"))

composition_uncapped <- no_cutoff %>%
  filter(Stain == "Quad", Islet.Cells >= 15) %>%
  group_by(Diagnosis, Donor) %>%
  summarise(across(c(Count.ins, Count.glu, Count.soma, Count.PP), sum), .groups = "drop") %>%
  mutate(total = Count.ins + Count.glu + Count.soma + Count.PP,
         across(c(Count.ins, Count.glu, Count.soma, Count.PP), ~ .x / total)) %>%
  group_by(Diagnosis) %>%
  summarise(donors = n(),
            across(c(Count.ins, Count.glu, Count.soma, Count.PP), list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}"))

composition_by_diagnosis <- bind_rows(capped_244 = composition_capped, uncapped = composition_uncapped,
                                      .id = "window")
write_csv(composition_by_diagnosis, "Sensitivity/composition_by_diagnosis.csv")
print(composition_by_diagnosis, width = Inf)

#  =============================================================================
# Part 3: cells per islet against the beta-loss counterfactual
# =============================================================================

# Quad objects of every size. (Part 1's `quad` is the islet-only quad_data.)
quad_all <- no_cutoff %>% filter(Stain == "Quad")

# ND objects with every beta cell removed: T1D if beta cells were lost and
# nothing replaced them. Same >= 15-cell rule applied to every group.
nd_minus_beta <- quad_all %>%
  filter(Diagnosis == "ND") %>%
  mutate(Count.ins = 0, Islet.Cells = Count.glu + Count.soma + Count.PP) %>%
  filter(Islet.Cells > 0)

counterfactual <- bind_rows(quad_all %>% mutate(group = as.character(Diagnosis)),
                            nd_minus_beta %>% mutate(group = "ND minus beta")) %>%
  group_by(group, Donor) %>%
  summarise(cells_per_islet     = mean(Islet.Cells[Islet.Cells >= 15]),
            pct_cells_in_islets = 100 * sum(Islet.Cells[Islet.Cells >= 15]) / sum(Islet.Cells),
            delta_pct = 100 * sum(Count.soma[Islet.Cells >= 15]) / sum(Islet.Cells[Islet.Cells >= 15]),
            pp_pct    = 100 * sum(Count.PP[Islet.Cells >= 15])   / sum(Islet.Cells[Islet.Cells >= 15]),
            .groups = "drop")
write_csv(counterfactual, "Sensitivity/counterfactual_by_donor.csv")

counterfactual_summary <- counterfactual %>%
  group_by(group) %>%
  summarise(donors = n(), across(cells_per_islet:pp_pct, ~ median(.x, na.rm = TRUE)))
cat("\n=== ND, T1D and ND with beta cells removed (donor medians) ===\n")
print(counterfactual_summary)

# Reference value for the fitted ratio below (a ratio of donor medians, not a
# fitted contrast).
cf_ratio <- with(counterfactual_summary,
                 cells_per_islet[group == "ND minus beta"] / cells_per_islet[group == "ND"])

# Published total-cell model (all_npod_models_v4.Rmd, chunk count-diag), verbatim.
FE_COUNT        <- "(Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c"
S2_COUNT_T1D_ND <- -0.049   # Fig 7A, log scale (53.4 vs 50.9 cells, -4.8%)

fit_count <- function(d) {
  glmmTMB(as.formula(paste("Islet.Cells ~", FE_COUNT, "+ (1 | Donor) + (1 | Donor:ImageID)")),
          data = d, family = nbinom2)
}
count_ratio <- function(m, d) {
  summary(pairs(emmeans(m, ~ Diagnosis, data = d), reverse = TRUE),
          type = "response", infer = TRUE) %>%
    as.data.frame()
}
islets <- function(max_log) {
  quad_all %>%
    filter(Islet.Cells >= 15, log_Islet.Cells <= max_log) %>%
    prep() %>%
    mutate(ImageID = factor(ImageID))
}

# Published 15-244 window first; stop unless it reproduces Fig 7A.
capped   <- islets(5.5)
r_capped <- count_ratio(fit_count(capped), capped)
if (abs(log(r_capped$ratio) - S2_COUNT_T1D_ND) > 0.002) {
  stop(sprintf("Capped refit gives log(T1D/ND) = %.4f; Fig 7A has %.4f. Resolve before reading the uncapped fit.",
               log(r_capped$ratio), S2_COUNT_T1D_ND))
}
cat(sprintf("Capped refit reproduces Fig 7A (log T1D/ND = %.4f)\n", log(r_capped$ratio)))

uncapped   <- islets(Inf)
r_uncapped <- count_ratio(fit_count(uncapped), uncapped)

count_results <- bind_rows(capped = r_capped, uncapped = r_uncapped, .id = "window") %>%
  mutate(counterfactual_ratio    = cf_ratio,
         excludes_counterfactual = asymp.LCL > cf_ratio)
write_csv(count_results, "Sensitivity/cells_per_islet_ratio.csv")
cat("\n=== T1D / ND cells per islet (published NB2 model) vs the beta-loss counterfactual ===\n")
print(count_results %>%
        dplyr::select(window, ratio, asymp.LCL, asymp.UCL, p.value,
                      counterfactual_ratio, excludes_counterfactual))
# =============================================================================
# Part 4: segmented single-cell area, T1D vs ND
# =============================================================================

# Single-cell quad objects with no acinar cells in the object, so the object's
# average cell area is that one endocrine cell's segmented area.
singlets <- no_cutoff %>%
  filter(Stain == "Quad", Islet.Cells == 1, Count.Acinar == 0,
         !is.na(Avg.Cell.Area..um.2.), Avg.Cell.Area..um.2. > 0) %>%
  mutate(cell_type = case_when(Count.ins  == 1 ~ "beta", Count.glu == 1 ~ "alpha",
                               Count.soma == 1 ~ "delta", Count.PP  == 1 ~ "PP"),
         log_cell_area = log(Avg.Cell.Area..um.2.),
         Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
         Region    = factor(Region, levels = c("Head", "Body", "Tail")),
         Sex       = factor(Sex, levels = c("Female", "Male")),
         Donor     = factor(Donor),
         ImageID   = factor(ImageID))

# Published single-cell stained-area changes (Results), for comparison.
# (Published beta/alpha singlets also included Double slides; these are Quad only.)
published_stained <- c(beta = -41.2, alpha = -31.7, delta = -8.3, PP = 131.6)

singlet_size <- map_dfr(names(published_stained), function(ct) {
  d <- singlets %>%
    filter(cell_type == ct) %>%
    mutate(Age_c = scale(Age, center = TRUE, scale = FALSE)[, 1]) %>%
    droplevels()
  m <- glmmTMB(log_cell_area ~ (Diagnosis + Region + Sex + Age_c)^2 + Diagnosis:Region:Age_c +
                 (1 | Donor) + (1 | Donor:ImageID),
               data = d, family = gaussian())
  summary(pairs(emmeans(m, ~ Diagnosis, data = d), reverse = TRUE), infer = TRUE) %>%
    as.data.frame() %>%
    transmute(cell_type = ct, singlets = nrow(d),
              segmented_T1D_vs_ND_pct = 100 * (exp(estimate) - 1),
              lower = 100 * (exp(asymp.LCL) - 1), upper = 100 * (exp(asymp.UCL) - 1),
              p.value, published_stained_pct = published_stained[[ct]])
})
write_csv(singlet_size, "Sensitivity/singlet_cell_area.csv")
print(singlet_size)


# =============================================================================
# Part 5: T1D cells per islet against the counterfactual, directly
# =============================================================================

t1d_vs_cf <- bind_rows(nd_minus_beta, quad_all %>% filter(Diagnosis == "T1D")) %>%
  filter(Islet.Cells >= 15) %>%
  prep() %>%
  mutate(ImageID = factor(ImageID))
# "T1D / ND" in this ratio is T1D / (ND minus beta).
t1d_vs_counterfactual <- count_ratio(fit_count(t1d_vs_cf), t1d_vs_cf)
write_csv(as.data.frame(t1d_vs_counterfactual), "Sensitivity/cells_per_islet_T1D_vs_counterfactual.csv")
print(t1d_vs_counterfactual)

cf <- counterfactual %>% filter(group %in% c("T1D", "ND minus beta"))
counterfactual_tests <- tibble(
  cell       = c("delta", "PP"),
  wilcoxon_p = c(wilcox.test(delta_pct ~ group, data = cf)$p.value,
                 wilcox.test(pp_pct ~ group, data = cf)$p.value))
write_csv(counterfactual_tests, "Sensitivity/counterfactual_delta_pp_tests.csv")
print(counterfactual_tests)
