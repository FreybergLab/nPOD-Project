################################################################################
#  Proportion Models — T1D Disease Duration & Age at Onset (ordbetareg)        #
#  Run on HTC cluster                                                          #
#  Outcomes: Percent.ins, Percent.glu, Percent.soma, Percent.PP                #
#  Data: quad_data, Diagnosis == "T1D"                                         #
#  Saves to: Models/                                                           #
################################################################################

# ============================================================
# PACKAGE INSTALLATION
# ============================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!require("tidyverse", quietly = TRUE))   install.packages("tidyverse")
if (!require("brms", quietly = TRUE))        install.packages("brms")
if (!require("ordbetareg", quietly = TRUE))  install.packages("ordbetareg")
if (!require("cmdstanr", quietly = TRUE))    install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

if (!dir.exists(cmdstanr::cmdstan_path())) cmdstanr::install_cmdstan()

# ============================================================
# LOAD LIBRARIES
# ============================================================

library(ordbetareg)
library(brms)
library(cmdstanr)
library(tidyverse)

# ============================================================
# DATA PREPARATION
# ============================================================

if (!dir.exists("Models")) dir.create("Models")

quad_data <- readRDS("Data/quad_data.rds")

quad_t1d <- quad_data %>%
  filter(Diagnosis == "T1D") %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
    Region    = factor(Region, levels = c("Head", "Body", "Tail")),
    Sex       = factor(Sex, levels = c("Female", "Male")),
    Donor     = factor(Donor),
    ImageID   = factor(ImageID)
  ) %>%
  droplevels()

rm(quad_data); gc()
cat("quad_t1d: ", nrow(quad_t1d), "rows,", n_distinct(quad_t1d$Donor), "donors\n")

# ============================================================
# DISEASE DURATION MODELS
# ============================================================

# --- Beta Cell (DD) ---
cat("Starting ins_dd model...\n")
ins_dd <- ordbetareg(
  Percent.ins ~ (Disease.Duration_c + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ins_dd_model"
)
cat("ins_dd complete!\n")
rm(ins_dd); gc()

# --- Alpha Cell (DD) ---
cat("Starting glu_dd model...\n")
glu_dd <- ordbetareg(
  Percent.glu ~ (Disease.Duration_c + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/glu_dd_model"
)
cat("glu_dd complete!\n")
rm(glu_dd); gc()

# --- Delta Cell (DD) ---
cat("Starting soma_dd model...\n")
soma_dd <- ordbetareg(
  Percent.soma ~ (Disease.Duration_c + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/soma_dd_model"
)
cat("soma_dd complete!\n")
rm(soma_dd); gc()

# --- PP Cell (DD) ---
cat("Starting pp_dd model...\n")
pp_dd <- ordbetareg(
  Percent.PP ~ (Disease.Duration_c + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/pp_dd_model"
)
cat("pp_dd complete!\n")
rm(pp_dd); gc()

# ============================================================
# AGE AT ONSET MODELS
# ============================================================

# --- Beta Cell (AO) ---
cat("Starting ins_ao model...\n")
ins_ao <- ordbetareg(
  Percent.ins ~ (Disease.Duration_c + Region + Sex + age_at_onset_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:age_at_onset_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ins_ao_model"
)
cat("ins_ao complete!\n")
rm(ins_ao); gc()

# --- Alpha Cell (AO) ---
cat("Starting glu_ao model...\n")
glu_ao <- ordbetareg(
  Percent.glu ~ (Disease.Duration_c + Region + Sex + age_at_onset_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:age_at_onset_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/glu_ao_model"
)
cat("glu_ao complete!\n")
rm(glu_ao); gc()

# --- Delta Cell (AO) ---
cat("Starting soma_ao model...\n")
soma_ao <- ordbetareg(
  Percent.soma ~ (Disease.Duration_c + Region + Sex + age_at_onset_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:age_at_onset_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/soma_ao_model"
)
cat("soma_ao complete!\n")
rm(soma_ao); gc()

# --- PP Cell (AO) ---
cat("Starting pp_ao model...\n")
pp_ao <- ordbetareg(
  Percent.PP ~ (Disease.Duration_c + Region + Sex + age_at_onset_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:age_at_onset_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/pp_ao_model"
)
cat("pp_ao complete!\n")
rm(pp_ao); gc()

cat("\nAll T1D proportion models complete! (8 models)\n")
cat("Saved to Models/:\n")
list.files("Models/", pattern = "(ins|glu|soma|pp)_(dd|ao)") %>% cat(sep = "\n")
