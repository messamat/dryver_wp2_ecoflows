pkg_vec <- list(
  'ade4',
  'adespatial',
  'assertthat',
  'dbscan',
  'dplyr',
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
  'qs',
  'readxl',
  'reshape2',
  'rlang',
  'rprojroot',
  'sf',
  'sfnetworks',
  'skimr',
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
  'mapview',
  'data.table')
lapply(pkg_vec, function(p) {
  suppressWarnings(
    suppressPackageStartupMessages(
      library(p, character.only = TRUE, warn.conflicts = FALSE, quietly = FALSE)
    ))
})
