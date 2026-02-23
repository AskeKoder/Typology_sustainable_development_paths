# Sustainable Development Paths Article (R analysis)

This repository contains the code and data used for the work: "Turning resources into wellbeing: a global typology of sustainable development paths".
The workflow builds a multi-source panel dataset (1990-2020), imputes missing values,
performs clustering experiments, and runs well-being, environmental, and colonial
analyses with figures and summary outputs.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18750965.svg)](https://doi.org/10.5281/zenodo.18750965)

## Contents

- `1_DataPreparation.R` - loads and harmonizes raw sources (SPI + external indicators)
  and writes a scaled/merged dataset (`extendedDataScaled.csv`).
- `2.0_ImputeData.R` - explores missingness and imputes missing values (MICE + some
  interpolation), producing `2_ImputedData.csv`.
- `2.1_ImputationComparison.R` - compares imputation variants.
- `3_ClusterExperimentation_ver2.R` - runs clustering experiments on imputed data.
- `3.1_Cluster_PostProcessing.R` - post-processes/renames clusters and joins HDI.
- `4.*.R` - analysis scripts for well-being, environmental, colonial, and decoupling tests.

Data folders:
- `Data_WellBeing/`
- `Data_Environmental/`
- `Data_Colonial/`
- `Figures/` (generated outputs)
- `old/` (legacy scripts/outputs)

## System requirements

Software:
- R (latest available on Windows as of January 28, 2026)
- RStudio 2026.01.0 Build 392 (Posit Software)
- R packages (list versions if available):
  - `dplyr`
  - `readxl`
  - `tidyr`
  - `mice`
  - `ggplot2`
  - `NbClust`
  - `partitionComparison`
  - `rnaturalearth`
  - `countrycode`
  - `zoo`
  - other packages as needed by the scripts

Operating systems:
- Windows (Microsoft) or macOS (Apple) or Linux (open source)

Hardware:
- No non-standard hardware required.

## Installation guide

1. Install R (and RStudio if desired).
2. Install required packages:
   ```r
   install.packages(c(
     "dplyr","readxl","tidyr","mice","ggplot2","NbClust",
     "partitionComparison","rnaturalearth","countrycode","zoo"
   ))
   ```
3. Open the repo root as the working directory in R/RStudio.

Typical install time on a "normal" desktop computer:
- few minutes

## Instructions for adding new data (development, environmntal or colonial)

1. Place your input data in the appropriate `Data_*` folders and match the
   expected file names/columns referenced in `1_DataPreparation.R`.
2. Run:
   ```r
   source("1_DataPreparation.R")
   source("2.0_ImputeData.R")
   source("3_ClusterExperimentation_ver2.R")
   ```
3. Run any analysis script in `4.*.R` relevant to your study.

## Reproducibility (optional but recommended)

To reproduce all quantitative results in the manuscript:
- Run all the R scripts in the provided script order, with parameters, and any fixed random seeds.

## Outputs

- `extendedDataScaled.csv` - merged/scaled panel dataset.
- `2_ImputedData.csv` - imputed dataset used for clustering and analysis.
- `3_ClusterResults.RDS`, `3_clusterVariations.RDS`, `4_RankedClusters.RDS` -
  clustering and ranking artifacts.
- `Figures/` - plots, maps, and exported visuals.
