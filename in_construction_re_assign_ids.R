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

#Create data.table to work with
rivnet_inters_dt <- as.data.table(rivnet_smooth_inters_nopseudo) 

#Back up cat before edit
rivnet_inters_dt[, cat_copy := cat]

#Compute other lengths
rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]
rivnet_inters_dt[, length_smooth := sum(length_uid), by=UID_smooth]

#Compute the percentage length of that cat ("hydromod" reach) that is represented by
#that segment section (identified by UID)
rivnet_inters_dt[, length_uid_catper := length_uid/length_cat]
#Check whether this is the largest segment section for that cat 
rivnet_inters_dt[, length_uid_catmax := (length_uid == max(length_uid)), by=cat]
#Compute the percentage length of that full segment that is represented by
#that segment section (identified by UID)
rivnet_inters_dt[, length_uid_smoothper := length_uid/length_smooth]

#Computer number of full segments that a cat overlaps
rivnet_inters_dt[, n_seg_overlap := length(unique(UID_smooth)), by=cat]

#For first order segments where the cat is only represented by that segment,
#keep that cat for this section of the segment
#making sure that length_uid_catmax is TRUE is to deal with cats that are split
#between non-contiguous segment sections (see cat==2938 for Croatia)
rivnet_inters_dt[n_seg_overlap == 1 & strahler == 1, 
                 cat_smooth := cat]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For full segments that only have one cat associated with them, confirm cat
rivnet_inters_dt[length_uid == length_smooth & !(cat %in% cats_assigned), 
                 cat_smooth := cat]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For cats that are already assigned, remove their "cat" from other segments
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#Re-compute length_uid_smoothper, excluding segment sections whose initial
#cat was assigned elsewhere
rivnet_inters_dt[!is.na(cat), length_uid_smoothper := length_uid/sum(length_uid),
                 by=UID_smooth]

#For those segments where only one cat remains, 
#assign this cat to the entire segment
cat_smooth_toassign <-  rivnet_inters_dt[
  is.na(cat_smooth) & length_uid_smoothper == 1 & length_uid/length_smooth > (1/3),
  .(UID_smooth, cat)] %>%
  setnames('cat', 'cat_smooth')
rivnet_inters_dt[is.na(cat_smooth), 
                 cat_smooth := cat_smooth_toassign[.SD, on='UID_smooth', 
                                                   cat_smooth]]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For cats that are already assigned, remove their "cat" from other segments
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#For segment section with cat==NAs in first order streams, assign cat_smooth of
#upstream segment
cat_smooth_toassign <- rivnet_inters_dt[!is.na(cat_smooth) & (strahler == 1) & 
                   (to %in% rivnet_inters_dt[is.na(cat) & strahler == 1, unique(from)]), 
                 .(to, UID_smooth, cat_smooth)] %>%
  setnames('to', 'from')
rivnet_inters_dt[is.na(cat) & strahler == 1,
                 cat_smooth := cat_smooth_toassign[
                   .SD, on=c('from', 'UID_smooth'), cat_smooth]]

#Re-compute length_uid_smoothper, excluding segment sections whose initial
#cat was assigned elsewhere
rivnet_inters_dt[!is.na(cat), length_uid_smoothper := length_uid/sum(length_uid),
                 by=UID_smooth]

#For those segments where only one cat remains, 
#assign this cat to the entire segment
cat_smooth_toassign <-  rivnet_inters_dt[
  is.na(cat_smooth) & length_uid_smoothper == 1 & length_uid/length_smooth > (1/3),
  .(UID_smooth, cat)] %>%
  setnames('cat', 'cat_smooth')
rivnet_inters_dt[is.na(cat_smooth), 
                 cat_smooth := cat_smooth_toassign[.SD, on='UID_smooth', 
                                                   cat_smooth]]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For cats that are already assigned, remove their "cat" from other segments
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#Computer number of full segments that a cat overlaps
rivnet_inters_dt[, n_seg_overlap := length(unique(UID_smooth)), by=cat]

#For second order segments where the cat is only represented by that segment,
#keep that cat for this section of the segment
rivnet_inters_dt[n_seg_overlap == 1 & strahler <= 2 & is.na(cat_smooth), 
                 cat_smooth := cat]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For cats that are already assigned, remove their "cat" from other segments
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#For segment section with cat==NAs in first order streams, assign cat_smooth of
#upstream segment
cat_smooth_toassign <- rivnet_inters_dt[!is.na(cat_smooth) & (strahler == 1) & 
                                          (to %in% rivnet_inters_dt[is.na(cat) & strahler == 1, unique(from)]), 
                                        .(to, UID_smooth, cat_smooth)] %>%
  setnames('to', 'from')
rivnet_inters_dt[is.na(cat) & strahler == 1,
                 cat_smooth := cat_smooth_toassign[
                   .SD, on=c('from', 'UID_smooth'), cat_smooth]]

#Re-compute length_uid_smoothper, excluding segment sections whose initial
#cat was assigned elsewhere
rivnet_inters_dt[!is.na(cat), length_uid_smoothper := length_uid/sum(length_uid),
                 by=UID_smooth]

#For those segments where only one cat remains, 
#assign this cat to the entire segment
cat_smooth_toassign <-  rivnet_inters_dt[
  is.na(cat_smooth) & length_uid_smoothper == 1 & length_uid/length_smooth > (1/3),
  .(UID_smooth, cat)] %>%
  setnames('cat', 'cat_smooth')
rivnet_inters_dt[is.na(cat_smooth), 
                 cat_smooth := cat_smooth_toassign[.SD, on='UID_smooth', 
                                                   cat_smooth]]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For cats that are already assigned, remove their "cat" from other segments
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#For sections that represent over 75% of the total length of that cat, assign cat_smooth
rivnet_inters_dt[length_uid_catper>0.7 & is.na(cat_smooth), cat_smooth := cat]


#For cats that are already assigned, remove their "cat" from other segments
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#Re-compute length_uid_smoothper, excluding segment sections whose initial
#cat was assigned elsewhere
rivnet_inters_dt[!is.na(cat), length_uid_smoothper := length_uid/sum(length_uid),
                 by=UID_smooth]

#For those segments where only one cat remains, 
#assign this cat to the entire segment
cat_smooth_toassign <-  rivnet_inters_dt[
  is.na(cat_smooth) & length_uid_smoothper == 1 & length_uid/length_smooth > (1/3),
  .(UID_smooth, cat)] %>%
  setnames('cat', 'cat_smooth')
rivnet_inters_dt[is.na(cat_smooth), 
                 cat_smooth := cat_smooth_toassign[.SD, on='UID_smooth', 
                                                   cat_smooth]]

cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]

#For cats that are already assigned, remove their "cat" from other segments
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#For segment sections in strahler==1 & n_seg_overlap > 2, cat := NA
rivnet_inters_dt[strahler == 1 & n_seg_overlap > 2 & is.na(cat_smooth), 
                 cat := NA] 

#if only representative of cat, assign cat_smooth
rivnet_inters_dt[!is.na(cat_smooth), cat := cat_smooth]
rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]
rivnet_inters_dt[, length_uid_catper := length_uid/length_cat, by=cat]
rivnet_inters_dt[length_uid_catper > 0.7 & is.na(cat_smooth), cat_smooth := cat]

#For cats that are already assigned, remove their "cat" from other segments
cats_assigned <- rivnet_inters_dt[!is.na(cat_smooth), unique(cat_smooth)]
rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_smooth), cat := NA]

#Re-compute length_uid_smoothper, excluding segment sections whose initial
#cat was assigned elsewhere
rivnet_inters_dt[!is.na(cat_smooth), 
                 length_uid_smoothper := length_uid/sum(length_uid),
                 by=UID_smooth]
cat_smooth_toassign <-  rivnet_inters_dt[
  is.na(cat_smooth) & length_uid_smoothper == 1,
  .(UID_smooth, cat)] %>%
  setnames('cat', 'cat_smooth')
rivnet_inters_dt[is.na(cat_smooth), 
                 cat_smooth := cat_smooth_toassign[.SD, on='UID_smooth', 
                                                   cat_smooth]]

#For the remaining segments with all NAs, 
#assign the cat of the section representing the highest percent of the segment
rivnet_inters_dt[is.na(cat_smooth), NAlength := sum(length_uid),
                 by=UID_smooth]
rivnet_inters_dt[(NAlength==length_smooth), 
                 cat_smooth := fifelse(length_uid == max(length_uid), cat_copy, NA), 
                 by=UID_smooth]

#For all remaining sections, remove their "cat" 
rivnet_inters_dt[is.na(cat_smooth), cat := NA]

for (i in 1:2) {
  #For segment section with cat==NAs and cat_smooth==NA, 
  #assign cat_smooth of upstream segment of same strahler order
  cat_smooth_toassign <- rivnet_inters_dt[
    !is.na(cat_smooth) &
      (from %in% rivnet_inters_dt[is.na(cat) & is.na(cat_smooth), unique(to)]), 
    .(from, UID_smooth, cat_smooth, strahler)] %>%
    setnames('from', 'to')
  rivnet_inters_dt[is.na(cat) & is.na(cat_smooth),
                   cat_smooth := cat_smooth_toassign[
                     .SD, on=c('to', 'UID_smooth'), cat_smooth]]
  
  #For is.na(cat) & is.na(cat_smooth), 
  #assign downstream cat_smooth of same UID_smooth
  cat_smooth_toassign <- rivnet_inters_dt[
    !is.na(cat_smooth) &
      (to %in% rivnet_inters_dt[is.na(cat) & is.na(cat_smooth), unique(from)]), 
    .(to, UID_smooth, cat_smooth, strahler)] %>%
    setnames('to', 'from')
  rivnet_inters_dt[is.na(cat) & is.na(cat_smooth),
                   cat_smooth := cat_smooth_toassign[
                     .SD, on=c('from', 'UID_smooth'), cat_smooth]]
}

#Check
write_sf(merge(rivnet_smooth_inters_nopseudo, 
               rivnet_inters_dt[, .(cat_smooth, UID)],
               by='UID')[, c('UID', 'cat', 'strahler', 'UID_smooth', 'cat_smooth')],
         file.path(resdir, 'scratch.gdb'), layer='test_cat16')
