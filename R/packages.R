pkg_vec <- list('ggplot2',
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
             'vegan',
             'mapview',
             'data.table')
lapply(pkg_vec, function(p) {
  suppressPackageStartupMessages(
    require(p, character.only = TRUE, warn.conflicts = FALSE, quietly = FALSE)
  )
})



