# R code for 'Contrasted biodiversity responses to drying in European rivers under climate change'

This repository contains R code associated with a manuscript in review.

## Abstract

Climate change is reshaping river flow regimes, expanding the spatial and temporal extent of river drying worldwide. Yet we lack quantitative predictions of how these hydrological shifts will impact freshwater biodiversity. Here we couple novel hydrological modeling with large-scale standardized biodiversity sampling to forecast drying-induced changes in local taxonomic diversity for multiple organism groups across six European river networks. We sampled aquatic macroinvertebrate, diatom, bacteria, and fungal communities up to six times over one year across 20-26 reaches in each of the Genal (Spain), Butižnica (Croatia), Velička (Czechia), Bükkösdi-víz (Hungary), Albarine (France), and Lepsämänjoki (Finland) river networks. These river networks represent a gradient of climatic conditions, biogeographic regions, and drying regimes. Sensitivity to drying duration, frequency, onset timing, and hydrological network connectivity differed among organism groups but remained consistent within groups over a 2700-km latitudinal gradient. Long-term hydrology better explained local biodiversity than short-term flow conditions, with biogeographical context modulating impacts. Under a high-emissions scenario, changes in drying patterns are projected to result in a decline of macroinvertebrate and diatom richness by an average of 14% and 11%, respectively, from 1991-2020 to 2071-2100, while fungal richness may increase by 35% — reshaping ecosystem structure and function. These findings underscore the urgency of integrating flow protection into climate adaptation planning to sustain resilient river ecosystems.

## Analysis structure and underlying data

This analysis relies as much as possible on [good enough practices in scientific computing](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510), which users are encouraged to read.

**Structure**: the overall project directory is structured with the following sub-directories:\
data/ (raw data, read-only, not to be altered)\
results/ (results of the analysis, mostly reproduceable through code execution. However, also includes manually modified results) R/ (code written for the project)\
\_targets.R bio_impacts_river_drying_europe.Rproj

All scripts rely on this structure.

**R Workflow**: this project is setup with a [targets workflow](https://docs.ropensci.org/targets/), ensuring reproducibility. In the `targets` philosophy, every action is a function, and every R object resulting from a workflow step is a "target" with dependencies. Intermediate targets/objects are stored in a `_targets` directory.

**Dependency management**: the R library of this project is managed by [renv](https://rstudio.github.io/renv/articles/renv.html). This makes sure that the exact same package versions are used when recreating the project. When calling `renv::restore()`, all required packages will be installed with their specific version. Please note that this project was built with R version 4.4 on a Windows 10 operating system.

**Syntax**: this analysis relies on the [data.table](https://rdatatable.gitlab.io/data.table/) syntax, which provides a high-performance version of data.frame. It is concise, faster, and more memory efficient than conventional data.frames and the tidyverse syntax.

## Getting started

### [Download and unzip the repository](https://anonymous.4open.science/api/repo/biodiversity_impacts_of_river_drying_in_europe-33DD/zip)

### Repository structure

- [**data/data_xxxxx**]() - supporting biological data
- [**R/**]() — core of the analysis
  - [*functions.R*]() - all custom functions used in the data formatting and analysis.
  - [*packages.R*]() - all packages used in the workflow.
- [*.Rprofile*]() — used to activate renv for new R sessions launched in the project.
- [*bio_impacts_river_drying_europe.Rproj*]() — R project file.
- [*LICENSE*]() - terms of use, modification and sharing for this software.
- [*README.md*]() — README for Github (this file)
- [*\_targets.R*]() — configuration script for targets workflow, this specific file name is required by the targets package. Contains the targets “plan”, the high-level catalog of all the steps in the workflow (see the corresponding chapter in the targets user manual). This plan defines the order of functions to use, their inputs and outputs (usually, targets), and the relationship among targets and steps.
- [*renv.lock*]() — renv lockfile, describing the state of the project’s library (installed packages and their version).
- [*report_wp2_multiorganism.qmd*]() — report file in quarto format summarizing data and information to support manuscript writing

## Running the analysis

Provided that your were given the necessary data, the entire analysis can simply be re-run with the following code found in

```{r rmake, eval = FALSE}
source('_targets.R')
tar_make()
```

`tar_make()` is the central function of the targets approach. It runs all the steps of the workflow in the correct order, skipping any work that is already up to date. Because of how targets tracks global functions and objects as dependencies of targets, the use of `tar_make()` is needed to run the analysis pipeline in a clean reproducible environment. If all targets are up to date in the caching directory, then nothing will be run.

## Inspecting results

If you were provided intermediate targets (i.e., a `_targets/` directory; or once you have re-run the analysis), you can load individual targets in the environment with the following commands (even if the targets are not up to date due to e.g. a change in source path).

```{r loadtarg, eval = FALSE}
tar_load(summary_multiorganism_richness) #Load target in memory (R environment) with original target name as variable name 
summary_results_list <- tar_read(summary_multiorganism_richness) #Load target in memory with new variable name
```

## Workflow (in `_targets.R`)

#### **1. Project Initialization and Configuration**

- Load required R libraries
- Set root directory and working directory based on project structure
- Source custom R scripts:
  - `R/packages.R` (package management)
  - `R/functions.R` (core utility functions)
  - `R/SpaTemp_function_M_edit.R` (Spatial-temporal connectivity functions)
  - `bin/03_diversity_metrics.R` (Biodiversity calculation functions)
- Define directory paths
- Configure parallel processing
- Define metadata

#### **2. Data Ingestion Pipeline**

- Define data paths
- Download external data
- Read raw data (reach attributes, hydrological modeling, environmental data, site metadata, biological data)

#### **3. River Network Preprocessing**

- Clean networks
- Compute Strahler stream order
- Reassign reach IDs to match hydrological modeling data
- Create formatted spatial datasets

#### **4. Hydrological Metric Calculation**

- Historical Period Metrics
- Spatial-Temporal Connectivity (STcon)
- Flow Distance Metrics (Fdist)
- Compile and Summarize - Combine hydrological statistics and connectivity metrics into a single dataset (`hydrocon_sites_compiled`)
- Compute annual averages of hydrological metrics during the sampling year for each site
- Compute mean annual hydrological metrics for each reach in the network over the historical period (1960–2021)
- Compute mean annual hydrological metrics at sampling sites for future periods (1990–2100) across all scenarios (SSP126, SSP245, SSP370, SSP585) and GCMs

#### **5. Environmental and Biological Data Processing**

- Environmental Data: summarize environmental variables at sampling sites, averaged across campaigns - Incorporate drainage area data for Genal (Spain) sites
- Taxonomic Diversity Metrics: Compute local alpha diversity for each combination of: Sampling site, sampling campaign, DRN (Drying River Network), Organism group (macroinvertebrates, diatoms, fungi, bacteria), then merge ecological, environmental, and hydrological data
- Format data for Spatial Stream Network (SSN) Models

#### **6. Data Exploration and Visualization**

#### **7. Statistical Modeling: Sites × Date Analysis**

**7.1. Model Setup**

**7.2. Model Selection Workflow**

- **Step 1: Covariance Structure Screening:**
- Run initial SSN models with a single hydrological variable (`DurD365past`) for each organism group
- Test all 144 covariance structure combinations
- Select top 5 covariance structures per organism based on AIC
- **Step 2: Full Model Fitting**
  - For each organism group and diversity metric (`richness`, `invsimpson`): fit models with:
    - Base formula: `log10(basin_area_km2) + log10(basin_area_km2):country`
    - Partition formula: `~ as.factor(campaign)`
    - Random effects: `~ country`
    - All 92 hydrological variables (one at a time)
    - Top 5 covariance structures from Step 1
    - Gaussian family, ML estimation
  - For fungi and bacteria: Test both linear and parabolic relationships
- **Step 3: Model Selection**
- Select best covariance structure for each model
- Generate performance tables (AIC, R², p-values, coefficients)
- Decompose model variance (fixed effects, random effects, spatial covariance)
- **Step 4: Post-Hoc Analysis**
  - Compute estimated marginal trends (emtrends)
  - Generate predictions for best models
  - Plot: Observed vs. predicted values, Predictor vs. fitted values, Variance decomposition, and Estimated marginal means and trends

**7.3. Hypothesis Testing**

**7.5. Results Aggregation**

#### **8. Statistical Modeling: Summarized (Site-Averaged) Analysis**

- Model Setup
- Model Fitting by Organism Group
- Model Evaluation and Diagnostics
- Future Projections

#### **9. Final Outputs and Deliverables**

- Summary Tables
- Figures
- Hypothesis Testing Results

### **Lexicon**

| **Term** | **Definition** |
|----|----|
| **DRN** | Drying River Network. A river network characterized by intermittent flow (periods of flowing, pooling, and drying). The study focuses on 6 DRNs across Europe. |
| **miv, ept, och, dia, fun, bac** | Macroinvertebrates, EPT (Ephemeroptera, Plecoptera, Trichoptera), OCH (Odonata, Coleoptera, Hephemeroptera), diatoms, fungi, and bacteria |
| **sedi and biof** | Sediment and biofilm sampling substrates. Samples collected from riverbed sediments. |
| **STcon** | Spatial-Temporal Connectivity. A metric quantifying the degree to which river reaches are hydrologically connected over time. Can be **directed** (upstream only) or **undirected** (full network). |
| **Fdist** | Flow distance. The distance (in river network distance) to the nearest wet (flowing) reach. Can be **directed** (upstream only) or **undirected** (any direction). |
| **SSN** | Spatial Stream Network. A statistical framework for modeling data collected from river networks, accounting for spatial dependencies (upstream/downstream and Euclidean distance). |
| **tar_target** | A target in the `targets` R package pipeline. Represents a step in the workflow that produces an output (e.g., a dataset, model, or plot). Targets are only re-run if their dependencies change. |
| **hydromod** | Hydrological modeling data. Output from hydrological models simulating flow state and discharge in river reaches. |
| **reaches** | River segments. Discrete units of a river network, each with unique attributes (ID, length, slope, etc.). |
| **sspX** | Shared Socioeconomic Pathway. Future scenarios used in climate projections (e.g., SSP126, SSP245, SSP370, SSP585). Higher numbers indicate more severe climate change. |
| **GCM** | General Circulation Model. Climate models used to project future hydrological conditions. |
| **pseudo-R²** | Pseudo coefficient of determination. A measure of model fit for mixed-effects models, representing variance explained by fixed effects. |
| **emtrends** | Estimated Marginal Trends. The estimated relationship between a predictor and response variable, averaged over other predictors in the model. |
| **emmeans** | Estimated Marginal Means. The predicted response value for a given predictor level, averaged over other predictors. |
