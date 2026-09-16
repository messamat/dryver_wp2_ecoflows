# R code for 'Contrasted biodiversity responses to drying in European rivers under climate change'

This repository contains R code associated with a manuscript in review.

## Abstract
Climate change is reshaping river flow regimes, expanding the spatial and
temporal extent of river drying worldwide. Yet we lack quantitative predictions 
of how these hydrological shifts will impact freshwater biodiversity. 
Here we couple novel hydrological modeling with large-scale standardized 
biodiversity sampling to forecast drying-induced changes in local taxonomic
diversity for multiple organism groups across six European river networks. 
We sampled aquatic macroinvertebrate, diatom, bacteria, and fungal communities
up to six times over one year across 20-26 reaches in each of the Genal (Spain),
Butižnica (Croatia), Velička (Czechia), Bükkösdi-víz (Hungary), Albarine (France), 
and Lepsämänjoki (Finland) river networks. These river networks represent a gradient
of climatic conditions, biogeographic regions, and drying regimes. Sensitivity to
drying duration, frequency, onset timing, and hydrological network connectivity 
differed among organism groups but remained consistent within groups over a 2700-km
latitudinal gradient. Long-term hydrology better explained local biodiversity than
short-term flow conditions, with biogeographical context modulating impacts. Under 
a high-emissions scenario, changes in drying patterns are projected to result in a
decline of macroinvertebrate and diatom richness by an average of 14% and 11%,
respectively, from 1991-2020 to 2071-2100, while fungal richness may increase by 35%
— reshaping ecosystem structure and function. These findings underscore the urgency of 
integrating flow protection into climate adaptation planning to sustain resilient river ecosystems.


## Analysis structure and underlying data

This analysis relies as much as possible on [good enough practices in scientific computing](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510), which users are encouraged to read.

**Structure**: the overall project directory is structured with the following sub-directories:  
data/ (raw data, read-only, not to be altered)  
results/ (results of the analysis, mostly reproduceable through code execution. However, also includes manually modified results)
R/ (code written for the project)  
_targets.R
bio_impacts_river_drying_europe.Rproj


All scripts rely on this structure.

**R Workflow**: this project is setup with a [targets workflow](https://docs.ropensci.org/targets/), ensuring reproducibility.
In the `targets` philosophy, every action is a function, and every R object resulting from a workflow step is a "target" with dependencies.
Intermediate targets/objects are stored in a `_targets` directory. 

**Dependency management**: the R library of this project is managed by [renv](https://rstudio.github.io/renv/articles/renv.html).
This makes sure that the exact same package versions are used when recreating the project.
When calling `renv::restore()`, all required packages will be installed with their specific version. 
Please note that this project was built with R version 4.4 on a Windows 10 operating system.

**Syntax**: this analysis relies on the [data.table](https://rdatatable.gitlab.io/data.table/) syntax, which provides a high-performance version of data.frame. It is concise, faster, and more memory efficient than conventional data.frames and the tidyverse syntax.

## Getting started
### [Download and unzip the repository](https://anonymous.4open.science/api/repo/biodiversity_impacts_of_river_drying_in_europe-33DD/zip)

### Repository structure
- [**data/data_annika**]() - supporting biological data
- [**R/**]() — core of the analysis
  - [*functions.R*]() - all custom functions used in the data formatting and analysis. 
  - [*packages.R*]() - all packages used in the workflow.
- [*.Rprofile*]() — used to activate renv for new R sessions launched in the project.
- [*bio_impacts_river_drying_europe.Rproj*]() — R project file.
- [*LICENSE*]() - terms of use, modification and sharing for this software.
- [*README.md*]() — README for Github (this file)
- [*\_targets.R*]() — configuration script for targets workflow,  this specific file name is required by the targets package. Contains the targets “plan”, the high-level catalog of all the steps in the workflow (see the corresponding chapter in the targets user manual). This plan defines the order of functions to use, their inputs and outputs (usually, targets), and the relationship among targets and steps.
- [*renv.lock*]() — renv lockfile, describing the state of the project’s library (installed packages and their version).
- [*report_wp2_multiorganism.qmd*]()  — report file in quarto format summarizing data and information to support manuscript writing

## Running the analysis
Provided that your were given the necessary data, the entire analysis can simply be re-run with the following code found in 
```{r rmake, eval = FALSE}
source('_targets.R')
tar_make()
```
`tar_make()` is the central function of the targets approach. It runs all the steps of the workflow in the correct order, skipping any work that is already up to date. Because of how targets tracks global functions and objects as dependencies of targets, the use of `tar_make()`  is needed to run the analysis pipeline in a clean reproducible environment. If all targets are up to date in the caching directory, then nothing will be run.

## Inspecting results
If you were provided intermediate targets (i.e., a `_targets/` directory; or once you have re-run the analysis), you can load individual targets in the environment with the following commands (even if the targets are not up to date due to e.g. a change in source path). 
``` {r loadtarg, eval = FALSE}
tar_load(summary_multiorganism_richness) #Load target in memory (R environment) with original target name as variable name 
summary_results_list <- tar_read(summary_multiorganism_richness) #Load target in memory with new variable name
```
