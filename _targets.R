#Small changes to implement before next run
  #L22 functions: capitalize S in spain
  #L56 and L158 functions: catch metacols regardless of capitalization
  #Standardize country names across all input datasets
  #Standardize metacols across all input datasets


library(rprojroot)
rootdir <- rprojroot::find_root(has_dir('R'))
setwd(rootdir)

source('R/packages.R')
source("R/functions.R")

datdir <- file.path('data', 'data_annika')
resdir <- 'results'

tar_option_set(format = "qs")

metacols <- c('Campaign', 'Site', 'running_id', 'Country', 'Date',
              'summa_sample', 'Sample.ID', 'Sample_ID', 'organism')
countries_list <- c('Croatia', 'Czech Republic', 'Finland',
                    'France', 'Hungary', 'Spain')

#--------------------------  Define targets plan -------------------------------
list(
  #------------------------------- Define paths ----------------------------------
  tar_target(
    env_data_path, 
    file.path(datdir, "ENV_all_NAs_as_blank.csv"),
    format='file'
  ),
  
  tar_target(
    metadata_edna_path,
    file.path('data', 'metadata_DRYvER_eDNA_updated.xlsx'),
    format='file'
  ),
  
  tar_target(
    bio_data_paths,
    list(
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
    #, format='file'
  ),
  
  tar_target(
    interm90_data_paths,
    {path_list <- sapply(countries_list,
                         function(country) {
                           file.path(datdir, 'Datasets', 'Intermittence_Data',
                                     paste0(country, '_Local_Interm_90_d.csv'))
                         },
                         USE.NAMES = TRUE)
    path_list['France'] <- gsub('90_d', '90_d_corrected', path_list['France'])
    return(path_list)
    }
  ),
  
  #------------------------------- Read in data ----------------------------------
  tar_target(
    env_dt,
    fread(env_data_path)
  ),
  
  tar_target(
    interm90_dt,
    lapply(interm90_data_paths, fread) %>% 
      rbindlist 
  ),
  
  tar_target(
    metadata_edna,
    read.xlsx(metadata_edna_path, 
              sheetName = 'metadataDNA') %>%
      as.data.table
  ),
  
  tar_target(
    bio_dt,
    read_biodt(path_list = bio_data_paths,
               in_metadata_edna = metadata_edna)
  ),
  
  tar_target(
    sprich,
    lapply(bio_dt, function(dt) {
      calc_sprich(in_biodt=dt, 
                  in_metacols=metacols)
    }) %>% rbindlist
  )
  ,
  
  tar_target(
    alphadat_merged,
    merge_alphadat(in_env_dt = env_dt,
                   in_interm90_dt = interm90_dt,
                   in_sprich = sprich)
  )
  ,
  
  tar_target(
    alpha_cor,
    alphadat_merged[, list(meanS_totdur90_cor = stats::cor(mean_S, totdur90),
                           meanS_discharge_cor = stats::cor(mean_S, discharge)
    ), by=Country]
  )
  ,
  
  tar_target(
    alpha_cor_plots,
    plot_alpha_cor(alphadat_merged,
                   out_dir = file.path(resdir, 'Null_models'))
  )
  ,
  
  # fungi_biof_path <- file.path(datdir, 'Datasets', 'fungi_dna_Biof_removed_zero_rows_and_columns.csv')
  # fungi_biof <- fread(fungi_biof_path)
  # fungi_biof[Country == 'Czech Republic', Country := 'Czech']
  # countries_list <- unique(fungi_biof$Country)
  # fungi_biof[1:10, 1:10]
  
  tar_target(
    null_models,
    lapply(
      names(bio_dt)[!(names(bio_dt) %in% c('bac_sedi', 'bac_biof'))],
      function(organism_type) {
        print(organism_type)
        organism_dt <- bio_dt[[organism_type]]

        if (organism_type %in% c("bac_biof_nopools", "bac_sedi_nopools")) {
          in_nsimul = 99; in_method = 'greedyqswap'; in_thin = 100
        } else {
          in_nsimul = 999; in_method = 'quasiswap'; in_thin = 1
        }

        organism_dt[, compute_null_model_inner(
          in_dt = .SD,
          in_metacols = metacols,
          min_siteN = 2,
          nsimul = in_nsimul,
          method = in_method,
          thin = in_thin)
          , by = Country] %>%
          .[, organism := organism_type]
      }
    )  %>% rbindlist
  )
)

