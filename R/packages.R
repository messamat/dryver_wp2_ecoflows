pkg_vec <- list(
  'ade4',
  'adespatial',
  'assertthat',
  'dbscan',
  'dplyr',
  'forecast',
  'future.apply',
  'ggcorrplot',
  'ggnewscale',
  'ggnetwork',
  'ggplot2',
  'gridExtra',
  'Hmisc',
  'igraph',
  'lmerTest', 
  'lme4', 
  'lubridate',
  'lwgeom',
  'magrittr',
  'memuse',
  'MuMIn', 
  'ncdf4',
  'phangorn',
  'plotly',
  'predictmeans',
  'qs',
  'readxl',
  'reshape2',
  'rlang',
  'rprojroot',
  'sf',
  'sfnetworks',
  'skimr',
  'sjPlot',
  'SSN2',
  'SSNbler',
  'stringr',
  'tarchetypes',
  'targets',
  'terra',
  'tibble',
  'tidygraph',
  'tidyterra',
  'tidytext',
  'vegan',
  'data.table')
lapply(pkg_vec, function(p) {
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, character.only = TRUE, warn.conflicts = FALSE, quietly = FALSE)
    ))
})
