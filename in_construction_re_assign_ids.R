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


#Prepare a network without splits between confluences  by removing all pseudonodes


#Then, iteratively:
#Starting upstream in first order streams with IDs that are only on one line, and among them, the shortest lines
#Assign the ID that is a first-order stream in reach.par of which this line represents the largest portion in the network
#Deal with all first-order streams this way, making sure that any stream that was dealt with before cannot be reassigned.

#If a reach ID disappeared because of this process. Grab Q of downstream reach id, and assign "dry" in any of the two are dry.

#Parameters

in_country <- 'Croatia'
rivnet_path <- tar_read(network_ssnready_gpkg_list)[[in_country]]
strahler_dt <- tar_read(network_strahler)[[in_country]] 
in_reaches_dt <- tar_read(reaches_dt)[country==in_country,] %>%
  .[, .(ID, to_reach, length)] %>%
  setnames(c('ID_hydromod', 'to_reach_hydromod', 'length_hydromod'))

#Read network and join with hydromod data
rivnet <- st_read(rivnet_path) %>%
  merge(strahler_dt, by='UID') %>%
  merge(in_reaches_dt, by.x='cat', by.y='ID_hydromod')

#Remove pseudonodes
rivnet_smooth <- as_sfnetwork(rivnet) %>%
  activate("edges") %>%
  convert(to_spatial_smooth) %>%
  st_as_sf() 
rivnet_smooth$UID_smooth <- seq_len(nrow(rivnet_smooth))

# write_sf(rivnet_smooth[, c('geometry','strahler', 'UID_smooth')],
#          file.path(resdir, 'scratch.gdb'), layer='test_6')

#Intersect with original network to check length of overlap for each 'cat'
rivnet_smooth_inters <- st_intersection(rivnet, 
                                        rivnet_smooth[, c('geometry', 'UID_smooth')]) %>%
  .[st_geometry_type(.) != 'POINT',] %>%
  st_cast('MULTILINESTRING') %>%
  st_cast('LINESTRING')

# Remove pseudo-nodes among segments sections of the same cat -------------------
rivnet_smooth_inters_nopseudo <- as_sfnetwork(rivnet_smooth_inters) %>%
  activate(edges) %>%
  convert(to_spatial_smooth, 
          require_equal = c("cat", "UID_smooth"),
          summarise_attributes = "first") %>%
  st_as_sf()

# write_sf(rivnet_smooth_inters_nopseudo[
#   , c('geometry', 'cat', 'UID', 'strahler', 'UID_smooth')],
#          file.path(resdir, 'scratch.gdb'), layer='test_inters3')

#Compute UID length
rivnet_smooth_inters_nopseudo$length_uid <- as.numeric(st_length(rivnet_smooth_inters_nopseudo))

#Compute other lengths
rivnet_inters_dt <- as.data.table(rivnet_smooth_inters_nopseudo) 
rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]
rivnet_inters_dt[, length_smooth := sum(length_uid), by=UID_smooth]

#Compute the percentage length of that cat ("hydromod" reach) that is represented by
#that segment section (identified by UID)
rivnet_inters_dt[, length_uid_catper := length_uid/length_cat]
#Check whether this is the largest segment section for that cat 
rivnet_inters_dt[, length_uid_catmax := (length_uid == max(length_uid)), by=cat]

#Computer number of full segments that a cat overlaps
rivnet_inters_dt[, n_seg_overlap := length(unique(UID_smooth)), by=cat]

#For first order segments where the cat is only represented by that segment,
#keep that cat for this section of the segment
#making sure that length_uid_catmax is TRUE is to deal with cats that are split
#between non-contiguous segment sections (see cat==2938 for Croatia)
rivnet_inters_dt[n_seg_overlap == 1 & strahler == 1 & length_uid_catmax, 
                 cat_smooth := cat]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For full segments that only have one cat associated with them, confirm cat
rivnet_inters_dt[length_uid == length_smooth & !(cat %in% cats_assigned), 
                 cat_smooth := cat]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For first order segments, for cats for which at least 50% of the length is 
#in that segment, assign the cat to that segment section
rivnet_inters_dt[strahler == 1 & length_uid_catper > 0.5, ]

############TO BE CONTINUED

##########################################################
########################################################

#Re-assign cat from cat_smooth
rivnet_inters_dt[!is.na(cat_smooth), cat := cat_smooth]

#Re-compute total cat length
rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]

#Then, iteratively:
#Starting upstream in first order streams with IDs that are only on one line, and among them, the shortest lines
#Assign the ID that is a first-order stream in reach.par of which this line represents the largest portion in the network
#Deal with all first-order streams this way, making sure that any stream that was dealt with before cannot be reassigned.


rivnet_inters_dt[cat==2960,]
rivnet_inters_dt[UID_smooth==569,]





ggplotly(
ggplot(rivnet_smooth_inters) +
  geom_sf(aes(color=cat))
)


#Only consider hydromod_reaches that are represented in the shapefile
reaches_sub <- in_reaches_dt[ID %in% rivnet$cat,]













#Preformat basic network
sfnet_ini <- rivnet %>%
  as_sfnetwork %>%
  activate("edges") %>%
  filter(!edge_is_multiple()) %>% #Keep shortest of edges that connect the same pair of nodes
  filter(!edge_is_loop()) #Remove obvious loops: edges that start and end at the same node

#------------------ Split lines at intersections -----------------------------
#Get confluence nodes (nodes of third degree: with at least 3 intersecting edges)
#2nd degree nodes are pseudonodes and 1st degree nodes are dangling
splitting_nodes <- activate(sfnet_ini, nodes) %>%
  mutate(degree = igraph::degree(.)) %>%
  filter(degree >= 3) %>%
  st_as_sf("nodes")

