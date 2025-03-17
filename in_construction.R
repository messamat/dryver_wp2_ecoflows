library(rprojroot)
rootdir <- rprojroot::find_root(has_dir("R"))
setwd(rootdir)

source("R/packages.R")
source("R/functions.R")

hydromod_present_dir <- file.path("data", "wp1", "Results_present_period_final", "data")
bio_dir <- file.path("data", "wp2", "01_WP2 final data")
# datdir <- file.path('data', 'data_annika')
resdir <- "results"

tar_option_set(format = "qs")

metacols <- c(
  "campaign", "site", "running_id", "country", "date",
  "summa_sample", "sample.id", "sample_id", "organism"
)
drn_dt <- data.table(
  country = c("Croatia", "Czech", "Finland", "France", "Hungary", "Spain"),
  catchment = c("Butiznica", "Velicka", "Lepsamaanjoki", "Albarine", "Bukkosdi", "Genal")
)

hydro_combi <- expand.grid(
  in_country = drn_dt$country,
  in_varname = c("isflowing", "qsim"),
  stringsAsFactors = FALSE
)

#Problem: 'cat' or 'reach_id' in the shapefile, associated with lines across 
#confluences. Makes it impossible to correctly use network hydrology data.

in_country <- 'Croatia'



#Identify real pseudonodes (that are reflected in hydrological models):
#   - Temporarily remove cats that are not represented in the shapefile
#   - Identify pairs of cats that are linked by a pseudonode
#Prepare a network without splits between confluences  by removing all other pseudonodes

#Intersect with original network to check length of overlap for each 'cat'

#Then, iteratively:
#Starting upstream in first order streams and the shortest lines
#Assign the ID that is a first-order stream in reach.par of which this line represents the largest portion in the network
#Deal with all first-order streams this way, making sure that any stream that was dealt with before cannot be reassigned.

check <- readRDS("C:\\DRYvER_wp2\\dryver_wp2_ecoflows\\data\\lysandre\\dryver-main\\Data\\01_albarine_hydro_area_constr.rds")
plot(check$river_graph)
gp <- ggplot(check$river_graph, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "black") +
  geom_nodes(color = "black", size = 8) +
  geom_nodetext(aes(label = name), color='red',
                fontface = "bold") +
  theme_blank()
  ggplotly(gp)
  
check$river_graph
