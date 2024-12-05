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
  'skimr',
  'stringr',
  'tarchetypes',
  'targets',
  'terra',
  'tibble',
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



