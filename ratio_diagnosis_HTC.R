################################################################################
#  Ratio Models — Diagnosis (ordbetareg)                                       #
#  Run on HTC cluster                                                          #
#  Outcome: beta_alpha_ratio = beta/(alpha+beta)                               #
#  Data: full_data, Stain != "Single"                                          #
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

ratio_data <- full_data %>%
  filter(Stain != "Single") %>%
  filter(Donor != "6473") %>%
  filter(Count.ins > 0 | Count.glu > 0) %>%
  mutate(beta_alpha_ratio = Count.ins / (Count.ins + Count.glu))

rm(full_data); gc()
cat("ratio_data: ", nrow(ratio_data), "rows,", n_distinct(ratio_data$Donor), "donors\n")

# ============================================================
# MODELS
# ============================================================

# --- Ratio Powered ---
cat("Starting ratio_powered model...\n")
ratio_powered <- ordbetareg(
  beta_alpha_ratio ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = ratio_data,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ratio_powered_model"
)
cat("ratio_powered complete!\n")
rm(ratio_powered); gc()

# --- Ratio Mid ---
cat("Starting ratio_mid model...\n")
ratio_mid <- ordbetareg(
  beta_alpha_ratio ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:Sex +
    Diagnosis:Region:Age_c +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Sex:Age_c +
    Diagnosis:Sex:log_Islet.Cells_c +
    Region:Sex:Age_c +
    Region:Sex:log_Islet.Cells_c +
    Diagnosis:Region:Sex:Age_c +
    Diagnosis:Region:Sex:log_Islet.Cells_c +
    (1|Donor) + (1|Donor:ImageID),
  data = ratio_data,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ratio_mid_model"
)
cat("ratio_mid complete!\n")
rm(ratio_mid); gc()

# --- Ratio Full ---
cat("Starting ratio_full model...\n")
ratio_full <- ordbetareg(
  beta_alpha_ratio ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^5 +
    (1|Donor) + (1|Donor:ImageID),
  data = ratio_data,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ratio_full_model"
)
cat("ratio_full complete!\n")
rm(ratio_full); gc()

cat("All ratio Diagnosis models complete!\n")
cat("Saved to Models/:\n")
list.files("Models/", pattern = "ratio_") %>% cat(sep = "\n")
