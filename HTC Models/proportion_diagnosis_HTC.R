################################################################################
#  Proportion Models — Diagnosis (ordbetareg)                                  #
#  Run on HTC cluster                                                          #
#  Outcomes: Percent.ins, Percent.glu, Percent.soma, Percent.PP                #
#  Data: quad_data                                                             #
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

quad_data <- read_csv("Data/quad_data.csv")

quad_data <- quad_data %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
    Region    = factor(Region, levels = c("Head", "Body", "Tail")),
    Sex       = factor(Sex, levels = c("Female", "Male")),
    Donor     = factor(Donor),
    ImageID   = factor(ImageID)
  )

cat("quad_data: ", nrow(quad_data), "rows,", n_distinct(quad_data$Donor), "donors\n")

# ============================================================
# MODELS
# ============================================================

# --- Beta Cell Proportion ---
cat("Starting ins_powered model...\n")
ins_powered <- ordbetareg(
  Percent.ins ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/ins_powered_model"
)
cat("ins_powered complete!\n")
rm(ins_powered); gc()

# --- Alpha Cell Proportion ---
cat("Starting glu_powered model...\n")
glu_powered <- ordbetareg(
  Percent.glu ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/glu_powered_model"
)
cat("glu_powered complete!\n")
rm(glu_powered); gc()

# --- Delta Cell Proportion ---
cat("Starting soma_powered model...\n")
soma_powered <- ordbetareg(
  Percent.soma ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/soma_powered_model"
)
cat("soma_powered complete!\n")
rm(soma_powered); gc()

# --- PP Cell Proportion ---
cat("Starting pp_powered model...\n")
pp_powered <- ordbetareg(
  Percent.PP ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4, chains = 4, threads = threading(4),
  iter = 2000, warmup = 1000,
  backend = "cmdstanr", seed = 123,
  file = "Models/pp_powered_model"
)
cat("pp_powered complete!\n")
rm(pp_powered); gc()

cat("All proportion Diagnosis models complete!\n")
cat("Saved to Models/:\n")
list.files("Models/", pattern = "(ins|glu|soma|pp)_powered") %>% cat(sep = "\n")
