# nPOD-Project

R code for the analysis reported in:

> Ward, Banks-Tibbs, Thorpe, *et al.* Artificial intelligence–guided analysis of human islet sociology reveals compositional remodeling in type 1 diabetes.

The code builds analysis datasets from HALO image-analysis exports of pancreatic sections from nPOD and HPAP organ donors. It then fits the statistical models: mixed models with glmmTMB locally, and Bayesian ordered beta regression with brms/cmdstanr on a computing cluster. Finally, it produces the figures, Table 1, and Supplementary Table 2. The raw image-analysis exports are not included (see [Data](#data)).

## Repository structure

```
nPOD-Project/
├── Create_count_datasets_v1.Rmd    Step 1: build analysis datasets from the HALO exports
├── Create_marker_datasets_v1.Rmd   Step 2: build Ki-67 / CD3 / CD45 datasets
├── all_npod_models_v4.Rmd          Step 3: density, count and area models (glmmTMB)
├── HTC Models/                     Step 4: composition and β/(α+β) models (brms, cluster)
├── Marker/                         Step 5: Ki-67 and CD3 rate models
├── Density/  Proportion/  Count/  Ratio/  Area/
│                                   Step 6: result workbooks and figures for each analysis
├── summary_dot_map_v3.Rmd          Supplementary Figure 23 (results summary)
├── Demo_table.Rmd                  Table 1
├── colors_master.R                 Shared plotting palette
├── Pseudotime/                     [TO ADD] UMAP and pseudotime analysis
├── Sensitivity/                    [TO ADD] sensitivity analyses
└── Results/                        Model-output workbooks behind Supplementary Table 2
```

Each analysis folder (`Density/`, `Proportion/`, `Count/`, `Ratio/`, `Area/`) has the same layout: an orchestrating `.Rmd`, an `R/` folder with the functions it sources, and `Models/`, `Results/` and `Graphs/` folders that the scripts fill.

## Data

**Raw inputs (not included; too large for GitHub).** These are the HALO per-object and summary exports for three analyses: the phenotyper (cell counts per object), area quantification, and the Ki-67, CD3 and CD45 marker analyses. The tissue comes from nPOD and HPAP. 

Data is available at this link: DOI:10.5281/zenodo.23423201

Set the raw-data locations before running:

| Variable | File | Points to |
|---|---|---|
| `count_base_path` | `Create_count_datasets_v1.Rmd` | Phenotyper exports (default: `corrected nPOD Data/` in the repository root) |
| `area_base_path` | `Create_count_datasets_v1.Rmd` | Area-quantification exports |
| Ki-67, CD3 and CD45 paths | `Create_marker_datasets_v1.Rmd` (top of file) | Marker exports |

**Metadata inputs.** Place these in the repository root. **[TO ADD to the repository, if nPOD data-use terms allow.]**

- `01282025_Full.nPOD.csv`: donor clinical metadata
- `remove_list.xlsx`: objects excluded after manual image QC
- `No_stain_data_stain_added.xlsx`: stain assignments for images whose exports lacked them
- `Supplemental_Table_1.xlsx`: Supplemental Table 1, read by `Demo_table.Rmd`

**Generated data.** Steps 1–2 write everything else to `Data/` as `.rds`, with `.csv` copies. Files ending in `_t1d` contain T1D donors only.

| File | Contents |
|---|---|
| `no_cutoff_data` | All endocrine objects |
| `cutoff_data` | Objects passing the islet size filter |
| `quad_data` | Islets from quadruple-stained sections (composition and count analyses) |
| `ratio_data` | Islets used for the β/(α+β) analyses |
| `density_data` | Islet density per image |
| `area_data_*` | Hormone-positive area by cell type and object class |
| `phenotyper_geometry` | Object bounding boxes, used to match marker and phenotyper objects |
| `ki67_data`, `cd3_data`, `cd45_data` | Marker analyses |

Endocrine-object size classes and the islet thresholds are defined in `Create_count_datasets_v1.Rmd`, in the section that applies the cutoffs.

## Requirements

**R.** Local analyses: R 4.4.2. Cluster models: R 4.5.0.

**Packages:**
- **Data handling and plotting:** tidyverse (dplyr, readr, stringr, ggplot2 3.5.1), janitor, readxl, openxlsx, writexl, patchwork, cowplot, ggtext, rlang, knitr, officer, flextable
- **Mixed models:** glmmTMB 1.1.9, emmeans 1.10.4, car
- **Bayesian models:** brms 2.22.0, ordbetareg 0.8, cmdstanr 0.9.0 with CmdStan 2.37.0, coda

## Running the analysis

All scripts assume the working directory is the repository root. The analysis files in `Density/`, `Proportion/`, `Count/`, `Ratio/`, `Area/` and `Marker/` set this themselves when knitted. Run the root-level files from the root, for example by opening the folder as an RStudio project.

**Step 0: create the output folders.** `all_npod_models_v4.Rmd` and `Create_count_datasets_v1.Rmd` do not create their output folders, so create them first.

```r
for (d in c("Data", "Density/Models", "Count/Models",
            "Area/Models/Total_Area", "Area/Models/Per_Cell_Area",
            "Area/Models/Total_Area_AdjCount",
            "Proportion/Models", "Ratio/Models")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
```

**Step 1: build the analysis datasets.** Knit `Create_count_datasets_v1.Rmd`. It reads the phenotyper and area exports plus the metadata files, applies image QC and the size filters, and writes the datasets above to `Data/`.

**Step 2: build the marker datasets.** Knit `Create_marker_datasets_v1.Rmd`. It matches each Ki-67, CD3 and CD45 marker object to its phenotyper islet by bounding box, then writes `Data/{ki67,cd3,cd45}_data.{rds,csv}`.

**Step 3: fit the local models.** Knit `all_npod_models_v4.Rmd`. It fits the islet density models (Gamma), the cell-count models (negative binomial) and the hormone-positive area models (Gaussian on log area), saving them to `Density/Models/`, `Count/Models/` and `Area/Models/`.

**Step 4: fit the Bayesian models on a cluster.** Copy `Data/quad_data.csv`, `Data/quad_data.rds` and `Data/cutoff_data.rds` into a `Data/` folder on the cluster, then run the scripts in `HTC Models/`. Every model is a cmdstanr fit with 4 chains of 2,000 iterations (1,000 warm-up), 4 threads per chain and seed 123. The published runs used 16 cores and 64 GB per job, and each fit can take many hours. Copy the resulting `.rds` files back as listed:

| Script | Models | Copy to |
|---|---|---|
| `proportion_diagnosis_HTC.R` | `{ins,glu,soma,pp}_powered_model` | `Proportion/Models/` |
| `ord_beta_reg_HTC.R` | `{ins,glu,soma,pp}_{mid,full}_model` | `Proportion/Models/` |
| `proportion_t1d_HTC.R` | `{ins,glu,soma,pp}_{ao,dd}_model` | `Proportion/Models/` |
| `ratio_diagnosis_HTC.R` | `ratio_{powered,mid,full}_model` | `Ratio/Models/` |
| `ratio_t1d_HTC.R` | `ratio_{ao,dd}_model` | `Ratio/Models/` |

All five scripts read from `Data/`. `ord_beta_reg_HTC.R` writes its model files to the working directory; the others write to `Models/`.

**Step 5: run the Ki-67 and CD3 analyses** (`Marker/`).

| File | What it does |
|---|---|
| `rate_models.Rmd` | Fits the per-cell rate models (negative binomial with a log cell-count offset); saves them and `rate_results.rds` to `Marker_Rate_Models/` |
| `rate_units.Rmd` | Tabulates rate-model effects (`rate_models_effects.csv`, `rate_models_units.csv`) |
| `survivor_selection_test.Rmd` | Tests whether survivor selection could explain the T1D composition (`selection_burden_per_donor.csv`, `selection_density_inequality.csv`) |
| `Marker_analysis.Rmd` and `R/` | **[TO ADD]** Result workbooks and Supplementary Figure 16 |

**Step 6: generate the result workbooks and figures.** Knit each file below. They write workbooks to `<folder>/Results/` and figures to `<folder>/Graphs/`.

| File | Produces |
|---|---|
| `Density/Density_analysis.Rmd` | Figure 2; Supplementary Figure 5 |
| `Proportion/Proportion_analysis.Rmd` | Figures 3–6 and 7E; Supplementary Figures 6–9 and 11–15 |
| `Count/Count_analysis.Rmd` | Figure 7A–D |
| `Ratio/Ratio_analysis.Rmd` | Supplementary Figure 10 |
| `Area/total_area_analysis.Rmd` | Figure 8; Supplementary Figure 18 |
| `Area/per_cell_area_analysis.Rmd` | Supplementary Figures 17 and 19 |
| `summary_dot_map_v3.Rmd` | Supplementary Figure 23 |
| `Demo_table.Rmd` | Table 1 |

`Area/R/fit_count_adj_models.R` fits the count-adjusted area models, whose outputs are in `Results/Total Area AdjCount/`.

**Step 7: pseudotime.** **[TO ADD]** UMAP and pseudotime code for Figure 9 and Supplementary Figures 20–22.

**Step 8: sensitivity analyses.** **[TO ADD]**

## Figures not produced by code in this repository

- **Figure 1** is a schematic.
- **Supplementary Figure 1** shows annotated images.
- **Supplementary Figures 2–3** (classifier validation) were made in GraphPad Prism 10.4.0 from manual versus AI cell counts.
- **Supplementary Figure 4** (endocrine-object size distribution and thresholds).

## Results folder

`Results/` contains the model-output workbooks that Supplementary Table 2 is built from, organized by analysis: density, cell-type proportions and counts, the β/(α+β) ratio, total and per-cell area, and the Ki-67 and CD3 rates.

## Citation, license and contact

- **Citation:** https://www.biorxiv.org/content/10.64898/2026.06.05.730081v1
- **Contact:** freyberg@pitt.edu
