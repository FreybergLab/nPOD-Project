################################################################################
#  Ratio Models — T1D Disease Duration & Age at Onset (ordbetareg)             #
#  Run on HTC cluster                                                          #
#  Outcome: beta_alpha_ratio = beta/(alpha+beta)                               #
#  Data: full_data, Stain != "Single", Diagnosis == "T1D"                      #
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

full_data <- readRDS("Data/cutoff_data.rds")

ratio_data_t1d <- full_data %>%
  filter(Stain != "Single") %>%
  filter(Donor != "6473") %>%
  filter(Diagnosis == "T1D") %>%
  filter(Count.ins > 0 | Count.glu > 0) %>%
  mutate(beta_alpha_ratio = Count.ins / (Count.ins + Count.glu)) %>%
  droplevels()

rm(full_data); gc()
cat("ratio_data_t1d: ", nrow(ratio_data_t1d), "rows,", n_distinct(ratio_data_t1d$Donor), "donors\n")

# ============================================================
# DISEASE DURATION MODEL
# ============================================================

cat("Starting ratio_dd model...\n")
ratio_dd <- ordbetareg(
  beta_alpha_ratio ~ (Disease.Duration_c + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = ratio_data_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ratio_dd_model"
)
cat("ratio_dd complete!\n")
rm(ratio_dd); gc()

# ============================================================
# AGE AT ONSET MODEL
# ============================================================

cat("Starting ratio_ao model...\n")
ratio_ao <- ordbetareg(
  beta_alpha_ratio ~ (Disease.Duration_c + Region + Sex + age_at_onset_c + log_Islet.Cells_c)^2 +
    Disease.Duration_c:Region:log_Islet.Cells_c +
    Disease.Duration_c:Region:age_at_onset_c +
    (1|Donor) + (1|Donor:ImageID),
  data = ratio_data_t1d,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ratio_ao_model"
)
cat("ratio_ao complete!\n")
rm(ratio_ao); gc()

cat("\nAll T1D ratio models complete! (2 models)\n")
cat("Saved to Models/:\n")
list.files("Models/", pattern = "ratio_(dd|ao)") %>% cat(sep = "\n")
