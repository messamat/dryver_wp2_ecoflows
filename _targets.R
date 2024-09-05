library(rprojroot)
rootdir <- rprojroot::find_root(has_dir('R'))
setwd(rootdir)

source('R/packages.R')
source("R/functions.R")

datdir <- file.path('data', 'data_annika')
resdir <- 'results'

tar_option_set(format = "qs")

#------------------------- DEFINE PATHS ----------------------------------------
#-------- WP2 data -------------------------------------------------------------
bio_data_paths <- list(
  dia_sedi = file.path(datdir, "diatoms_dna_sediment_removed_pools_and_zero_rows_and_columns.csv"),
  dia_biof = file.path(datdir, "diatoms_dna_biofilm_removed_pools_and_zero_rows_and_columns.csv"),
  fun_sedi = file.path(datdir, "fungi_dna_sediment_removed_pools_zero_rows_and_columns.csv"),
  fun_biof = file.path(datdir, "fungi_dna_Biof_removed_pools_zero_rows_and_columns.csv"),
  bac_sedi = file.path(datdir, "Sediment_bacteria_removed_pools_and_zero_species.csv"),
  bac_biof = file.path(datdir, "Biofilm_bacteria_removed_pools_and_zero_species.csv"),
  miv = file.path(datdir, "DRYvER_MIV_data_genus_reduced_taxa.csv"),
  miv_nopools = file.path(datdir, "DRYvER_MIV_data_genus_reduced_taxa_removed_pools.csv"),
  miv_nopools_flying = file.path(datdir, "DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_flying.csv"),
  miv_nopools_nonflying = file.path(datdir, "DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_nonflying.csv")
)

metacols <- c('Campaign', 'Site', 'running_id', 'Country', 'Date',
               'summa_sample', 'Sample.ID', 'Sample_ID', 'organism')

#--------------------------  Define targets plan -------------------------------
list(
  #------------------------------- Read files ----------------------------------
  #Read in metadata
  tar_target(
    bio_dt,
    read_biodt(path_list=bio_data_paths)
  ),

  tar_target(
   sprich,
   lapply(bio_dt, function(dt) {
     calc_sprich(in_biodt=dt, in_metacols=metacols)
   }) %>% rbindlist
  )
)
  


# 
#   
#   
# )
# 
# # UPDATE 1: to update bacteria, open old diversity data
# divdata <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity.csv")
# 
# 
# names(bac_biof_S_mean)
# 
# names(divdata)
# divdata$V1 <- NULL
# names(divdata)
# 
# divupdate <- rbind(divdata, bac_biof_S_mean)
# 
# # reorder columns in bacteria sediment diversity dataframes
# names(bac_sed_S_mean)
# bac_sed_S_mean2 <- bac_sed_S_mean[,c(3,1,2)]
# names(bac_sed_S_mean2)
# 
# divupdate <- rbind(divupdate, bac_sed_S_mean2)
# 
# #save as csv
# write.csv(divupdate, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity_updated.csv")
# 
# 
# 
# # UPDATE 2, 15.4.2024: update flying and non-flying macroinvertebrates, with no pools data
# 
# divdata <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity_updated.csv")
# 
# names(flymac_S_mean)
# names(nonflymac_S_mean)
# names(divdata)
# divdata$V1 <- NULL
# names(divdata)
# 
# 
# divupdate2 <- rbind(divdata, flymac_S_mean, nonflymac_S_mean)
# 
# 
# #save as csv
# write.csv(divupdate2, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity_updated2.csv", row.names = F)