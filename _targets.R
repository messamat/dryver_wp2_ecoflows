library(rprojroot)
rootdir <- rprojroot::find_root(has_dir('dryver_wp2_ecoflows'))
setwd(rootdir)

source('dryver_wp2_ecoflows/R/packages.R')
source("dryver_wp2_ecoflows/R/functions.R")

datdir <- 'WP2 - Predicting biodiversity changes in DRNs'
resdir <- 'results'

tar_option_set(format = "qs")

#------------------------- DEFINE PATHS ----------------------------------------
#------------ WP1 data ---------------------------------------------------------


#------------ WP2 data ---------------------------------------------------------
findatdir <- file.path(datdir, '01_WP2 final data')
data_description <- file.path(findatdir, '_Data description information', 
                              'Environmental variables_DRYvER.xlsx')

envi_rawdat <- file.path(findatdir, 
                         '_Environmental data', 
                         'ENV_all_1622023')

fish_rawdat <- list.files(file.path(findatdir,
                                    '_Fish_data')) %>%
  grep('[.]xlsx$', ., value=T) 

minv_rawdat <- list.files(file.path(findatdir,
                                    '_Macroinvertebrate data')) %>%
  grep('[.](xlsx|csv)$', ., value=T) 

bact_rawdat <- file.path(findatdir, 
                         '_Microbes data', 
                         'Metabarcoding_Bacteria_bact02', 
                         'dryver-bact02-clust_motus_samples.rds')

euka_rawdat <- file.path(findatdir, 
                         '_Microbes data', 
                         'Metabarcoding_Eukaryota_euka02', 
                         'dryver-euka02-clust_motus_samples.rds')


#------------------ IMPORT DATA ------------------------------------------------
check <- load("C:\\Users\\mamessager\\OneDrive - INRAE SharePoint SE\\DRyVER\\WP2 - Predicting biodiversity changes in DRNs\\upscaled_diversities\\future_drn\\Albarine\\Albarine_scenario_ssp370_model_ipsl-cm6a-lr_Member_10.Rdata")
