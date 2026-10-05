# ============================================================
# PACKAGE INSTALLATION
# ============================================================

# Set a CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install packages if not already installed
if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!require("brms", quietly = TRUE)) {
  install.packages("brms")
}

if (!require("ordbetareg", quietly = TRUE)) {
  install.packages("ordbetareg")
}

if (!require("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}

# cmdstanr requires CmdStan to be installed separately
# This checks if it's installed and installs if needed
if (!dir.exists(cmdstanr::cmdstan_path())) {
  cmdstanr::install_cmdstan()
}

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

quad_data <- read_csv("/ihome/hpark/het29/nPOD/quad_data.csv")

quad_data$log_Islet.Cells <- as.numeric(quad_data$log_Islet.Cells)

# Set factor levels explicitly
quad_data <- quad_data %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("ND", "T1D")),
    Region = factor(Region, levels = c("Head", "Body", "Tail")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    Donor = factor(Donor),
    ImageID = factor(ImageID)
  ) %>%
  mutate(
    Age_c = scale(Age, center = TRUE, scale = FALSE)[,1],
    log_Islet.Cells_c = scale(log_Islet.Cells, center = TRUE, scale = FALSE)[,1]
  )

# Decile calculations
log_islet_values <- quantile(
  quad_data$log_Islet.Cells,
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
  na.rm = TRUE
)

quad_data$log_islet_values <- cut(
  quad_data$log_Islet.Cells,
  breaks = c(-Inf, log_islet_values, Inf),
  labels = c("<10th", "10th-20th", "20th-30th", "30th-40th", "40th-50th", 
             "50th-60th", "60th-70th", "70th-80th", "80th-90th", ">90th"),
  include.lowest = TRUE
)

decile_lookup <- data.frame(
  log_Islet.Cells = log_islet_values,
  decile_label = c("10th", "20th", "30th", "40th", "50th", "60th", "70th", "80th", "90th")
)

# ============================================================
# MODELS - Each with UNIQUE file name for caching
# ============================================================

cat("Starting ins_powered model...\n")
ins_powered <- ordbetareg(
  Percent.ins ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "ins_powered_model"
)
cat("ins_powered complete!\n")

cat("Starting ins_mid model...\n")
ins_mid <- ordbetareg(
  Percent.ins ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    # 3-way interactions (excluding any with both Age_c and log_Islet.Cells_c)
    Diagnosis:Region:Sex +
    Diagnosis:Region:Age_c +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Sex:Age_c +
    Diagnosis:Sex:log_Islet.Cells_c +
    Region:Sex:Age_c +
    Region:Sex:log_Islet.Cells_c +
    # 4-way interactions
    Diagnosis:Region:Sex:Age_c +
    Diagnosis:Region:Sex:log_Islet.Cells_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "ins_mid_model"
)
cat("ins_mid complete!\n")

cat("Starting ins_full model...\n")
ins_full <- ordbetareg(
  Percent.ins ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^5 +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "ins_full_model"
)
cat("ins_full complete!\n")

cat("Starting glu_powered model...\n")
glu_powered <- ordbetareg(
  Percent.glu ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "glu_powered_model"
)
cat("glu_powered complete!\n")

cat("Starting glu_mid model...\n")
glu_mid <- ordbetareg(
  Percent.glu ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    # 3-way interactions (excluding any with both Age_c and log_Islet.Cells_c)
    Diagnosis:Region:Sex +
    Diagnosis:Region:Age_c +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Sex:Age_c +
    Diagnosis:Sex:log_Islet.Cells_c +
    Region:Sex:Age_c +
    Region:Sex:log_Islet.Cells_c +
    # 4-way interactions
    Diagnosis:Region:Sex:Age_c +
    Diagnosis:Region:Sex:log_Islet.Cells_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "glu_mid_model"
)
cat("glu_mid complete!\n")

cat("Starting glu_full model...\n")
glu_full <- ordbetareg(
  Percent.glu ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^5 +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "glu_full_model"
)
cat("glu_full complete!\n")

cat("Starting soma_powered model...\n")
soma_powered <- ordbetareg(
  Percent.soma ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "soma_powered_model"
)
cat("soma_powered complete!\n")

cat("Starting soma_mid model...\n")
soma_mid <- ordbetareg(
  Percent.soma ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    # 3-way interactions (excluding any with both Age_c and log_Islet.Cells_c)
    Diagnosis:Region:Sex +
    Diagnosis:Region:Age_c +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Sex:Age_c +
    Diagnosis:Sex:log_Islet.Cells_c +
    Region:Sex:Age_c +
    Region:Sex:log_Islet.Cells_c +
    # 4-way interactions
    Diagnosis:Region:Sex:Age_c +
    Diagnosis:Region:Sex:log_Islet.Cells_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "soma_mid_model"
)
cat("soma_mid complete!\n")

cat("Starting soma_full model...\n")
soma_full <- ordbetareg(
  Percent.soma ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^5 +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "soma_full_model"
)
cat("soma_full complete!\n")

cat("Starting pp_powered model...\n")
pp_powered <- ordbetareg(
  Percent.PP ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Region:Age_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "pp_powered_model"
)
cat("pp_powered complete!\n")

cat("Starting pp_mid model...\n")
pp_mid <- ordbetareg(
  Percent.PP ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^2 +
    # 3-way interactions (excluding any with both Age_c and log_Islet.Cells_c)
    Diagnosis:Region:Sex +
    Diagnosis:Region:Age_c +
    Diagnosis:Region:log_Islet.Cells_c +
    Diagnosis:Sex:Age_c +
    Diagnosis:Sex:log_Islet.Cells_c +
    Region:Sex:Age_c +
    Region:Sex:log_Islet.Cells_c +
    # 4-way interactions
    Diagnosis:Region:Sex:Age_c +
    Diagnosis:Region:Sex:log_Islet.Cells_c +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "pp_mid_model"
)
cat("pp_mid complete!\n")

cat("Starting pp_full model...\n")
pp_full <- ordbetareg(
  Percent.PP ~ (Diagnosis + Region + Sex + Age_c + log_Islet.Cells_c)^5 +
    (1|Donor) + (1|Donor:ImageID),
  data = quad_data,
  cores = 4,
  chains = 4,
  threads = threading(4),
  iter = 2000,
  warmup = 1000,
  backend = "cmdstanr",
  seed = 123,
  file = "pp_full_model"
)
cat("pp_full complete!\n")

# Save all models together in one RData file as backup
save(ins_powered, ins_mid, ins_full, 
     glu_powered, glu_mid, glu_full, 
     soma_powered, soma_mid, soma_full, 
     pp_powered, pp_mid, pp_full,
     quad_data, decile_lookup,
     file = "all_ordbeta_models.RData")

cat("All models complete and saved!\n")