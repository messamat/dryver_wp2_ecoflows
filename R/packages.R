pkg_vec <- list(
  'ggnewscale',
  'ggplot2',
  'gridExtra',
  'lmerTest', 
  'lme4', 
  'lubridate',
  'magrittr',
  'MuMIn', 
  'ncdf4',
  'qs',
  'readxl',
  'reshape2',
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
  'vegan',
  'mapview',
  'data.table')
lapply(pkg_vec, function(p) {
  suppressWarnings(
    suppressPackageStartupMessages(
      require(p, character.only = TRUE, warn.conflicts = FALSE, quietly = FALSE)
    ))
})



