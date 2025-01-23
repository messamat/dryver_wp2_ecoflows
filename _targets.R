#Small changes to implement before next run
#L22 functions: capitalize S in spain
#L56 and L158 functions: catch metacols regardless of capitalization
#Standardize country names across all input datasets
#Standardize metacols across all input datasets
#site names in "data/wp1/Results_present_period_final/data/Genal/Genal_sampling_sites_ReachIDs.csv" are incorrect
#update readxl
#Need to know how sites were snapped to network "C:\DRYvER_wp2\WP2 - Predicting biodiversity changes in DRNs\Coordinates\Shapefiles with sites moved to the river network\Croatia_near_coords.shp"
#3s flow acc for Europe: https://data.hydrosheds.org/file/hydrosheds-v1-acc/eu_acc_3s.zip

#Make sure that biological data are standardized by area to get densities
#Get ancillary catchment data


library(rprojroot)
rootdir <- rprojroot::find_root(has_dir('R'))
setwd(rootdir)

source('R/packages.R')
source("R/functions.R")

hydromod_present_dir <- file.path('data', 'wp1', 'Results_present_period_final', 'data')
bio_dir <- file.path('data', 'wp2', '01_WP2 final data')
#datdir <- file.path('data', 'data_annika')
resdir <- 'results'

tar_option_set(format = "qs")

metacols <- c('campaign', 'site', 'running_id', 'country', 'date',
              'summa_sample', 'sample.id', 'sample_id', 'organism')
drn_dt <- data.table(
  country = c("Croatia", "Czech", "Finland", "France",  "Hungary", "Spain"),
  catchment = c("Butiznica", "Velicka", "Lepsamaanjoki", "Albarine", "Bukkosdi", "Genal")
)

hydro_combi <- expand.grid(
  in_country = drn_dt$country,
  in_varname =  c('isflowing', 'qsim'),
  stringsAsFactors = FALSE)

#--------------------------  Define targets plan -------------------------------
list(
  #------------------------------- Define paths ----------------------------------
  #
  #Path to local environmental data
  tar_target(
    env_data_path_annika,
    file.path(rootdir, 'data', 'data_annika', "ENV_all_NAs_as_blank.csv"),
    format='file'
  ),
  
  tar_target(
    env_data_path_common,
    file.path(bio_dir,"_Environmental data", "ENV_all_1622023.csv"),
    format='file'
  ),
  
  #Path to modeled hydrological data
  tar_target(
    hydromod_paths_dt,
    define_hydromod_paths(in_hydromod_dir = hydromod_present_dir)
  ),
  
  #Path to metadata accompanying eDNA data
  tar_target(
    metadata_edna_path,
    file.path(bio_dir, '_Microbes data', 'metadata DRYvER eDNA updated.xlsx'),
    format='file'
  ),
  
  #Paths to pre-processed biological sampling data
  tar_target(
    bio_data_paths,
    list(
      dia_sedi = "diatoms_dna_sediment_removed_pools_and_zero_rows_and_columns.csv",
      dia_biof = "diatoms_dna_biofilm_removed_pools_and_zero_rows_and_columns.csv",
      fun_sedi = "fungi_dna_sediment_removed_pools_zero_rows_and_columns.csv",
      fun_biof = "fungi_dna_Biof_removed_pools_zero_rows_and_columns.csv",
      bac_sedi = "Sediment_bacteria_removed_pools_and_zero_species.csv",
      bac_biof = "Biofilm_bacteria_removed_pools_and_zero_species.csv",
      miv = "DRYvER_MIV_data_genus_reduced_taxa.csv",
      miv_nopools = "DRYvER_MIV_data_genus_reduced_taxa_removed_pools.csv",
      miv_nopools_flying = "DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_flying.csv",
      miv_nopools_nonflying = "DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_nonflying.csv"
    ) %>%
      lapply(function(path) file.path(rootdir, 'data', 'data_annika', path))
    #, format='file'
  )
  ,
  
  #------------------------------- Download data -------------------------------
  #Download amber river barriers dataset
  tar_target(
    amber_path,
    {
      amber_dir_path <- download_unzip(
        url =  "https://figshare.com/ndownloader/articles/12629051/versions/5",
        out_dir = file.path('data', 'amber'), 
        out_zip=NULL)
      unzip(file.path(amber_dir_path, 'Fig1_AMBER_BARRIER_ATLAS_V1.zip'),
            exdir = file.path('data', 'amber'))
      return(file.path(amber_dir_path, 'AMBER_BARRIER_ATLAS_V1.csv'))
    }
  ),
  
  #------------------------------- Read in data ----------------------------------
  #Read reach data
  tar_target(
    reaches_dt,
    lapply(hydromod_paths_dt$country, function(in_country) {
      hydromod_paths_dt[country==in_country, 
                        fread(reaches_path, sep="\t", skip=1, 
                              header=T, drop=c(17,18)) %>%
                          .[-c(1,2,3),] %>%
                          .[ID != '# end of reach.par',] %>%
                          .[, lapply(.SD, as.numeric)] %>%
                          setnames(gsub('[-]', '_', names(.)))
      ] %>%
        .[, country := in_country]
    }) %>% rbindlist(use.names=T, fill=T) 
    #Reach that flows to 9999 is outlet in Czech, Finland, France, and Hungary
    #Reach that flows to 0 is outlet in Spain and Croatia
  ),
  
  #Read local environmental data
  tar_target(
    env_dt,
    read_envdt(in_env_data_path_annika = env_data_path_annika,
               in_env_data_path_common = env_data_path_common)
  ),
  
  #Correct river site names 
  tar_target(
    sites_dt,
    hydromod_paths_dt[, format_site_dt(in_path = sites_reachids, 
                                       in_country = country),
                      by = country] 
  ),
  
  #Subset river network shapefiles to only keep sections within which there
  #are sampling site-reaches (based on sub-catchment file given by countries)
  tar_target(
    network_sub_gpkg_list,
    subset_network(in_hydromod_paths_dt = hydromod_paths_dt,
                   out_dir = file.path('results', 'gis'),
                   overwrite = T)
  )
  ,
  
  #Clean networks
  tar_target(
    network_clean_gpkg_list,
    lapply(names(network_sub_gpkg_list), function(in_country) {
      #Set size of simplifying radius to remove loops. See function
      clustering_dist_list <- list(Croatia=50, Czech=60, Finland=50, 
                                   France=50, Hungary=50, Spain=40)
      
      out_net_path <- clean_network(
        rivnet_path = network_sub_gpkg_list[[in_country]],
        idcol = 'cat',
        node_clustering_dist = clustering_dist_list[[in_country]],
        min_segment_length = 20,
        outdir = file.path(resdir, 'gis'),
        save_gpkg = TRUE, 
        return_path = TRUE)
      
      #Manual corrections
      if (in_country == 'Czech') {
        clean_net <- st_read(out_net_path) %>%
          .[!(.[['cat']] %in% c(40112, 40106, 40478)),]
        st_write(clean_net, out_net_path, append=F)
      } else if (in_country == 'Croatia') {
        clean_net <- st_read(out_net_path) %>%
          .[!(.[['UID']] %in% c(1102, 2988)),]
        st_write(clean_net, out_net_path, append=F)
      }
      
      return(out_net_path)
    }) %>% setNames(names(network_sub_gpkg_list)) 
  )
  ,
  
  tar_target(
    network_directed_gpkg_list,
    lapply(names(network_clean_gpkg_list), function(in_country) {
      #Set size of simplifying radius to remove loops. See function
      outlet_uid_list <- list(Croatia = 458, Czech = 4, Finland = 682,
                              France = 1, Hungary = 5, Spain = 86)
      
      out_net_path <- direct_network(
        rivnet_path = network_clean_gpkg_list[[in_country]],
        idcol = 'UID',
        outletid = outlet_uid_list[[in_country]],
        outdir = file.path(resdir, 'gis'), 
        save_gpkg = TRUE) 
      
      return(out_net_path)
    }) %>% setNames(names(network_clean_gpkg_list))
  )
  ,
  
  
  tar_target(
    network_ssnready_gpkg_list,
    lapply(names(network_directed_gpkg_list), function(in_country) {
      fix_complex_confluences(
        rivnet_path = network_directed_gpkg_list[[in_country]], 
        max_node_shift = 5,
        out_path= file.path(resdir, 'gis', 
                            paste0(tolower(in_country)
                                   , '_river_network_ssnready',
                                   format(Sys.time(), "%Y%m%d"), '.gpkg'
                            ))
      )
    }) %>% setNames(names(network_directed_gpkg_list))
  )
  ,
  
  #Compute strahler stream order
  tar_target(
    network_strahler,
    lapply(names(network_ssnready_gpkg_list), function(in_country) {
      assign_strahler_order(
        in_rivnet = network_ssnready_gpkg_list[[in_country]], 
        idcol = 'UID')
    }) %>% setNames(names(network_ssnready_gpkg_list))
  )
  ,
  
  #Re-assign correct IDs to match with hydrological data
  tar_target(
    network_reided,
    lapply(names(network_ssnready_gpkg_list), function(in_country) {
      reassign_netids(rivnet_path = network_ssnready_gpkg_list[[in_country]], 
                      strahler_dt = network_strahler[[in_country]], 
                      in_reaches_hydromod_dt = reaches_dt[country==in_country,], 
                      outdir = file.path(resdir, 'gis')
      )
    }) %>% setNames(names(network_ssnready_gpkg_list))
  )
  ,

  #Create shapefile of sampling site-reaches
  tar_target(
    site_reaches_gpkg_list,
    create_sites_gpkg(in_hydromod_paths_dt = hydromod_paths_dt,
                      in_sites_dt = sites_dt,
                      out_dir = file.path('results', 'gis'),
                      geom= 'reaches',
                      overwrite = T)
  )
  ,
  
  tar_target(
    site_points_gpkg_list,
    create_sites_gpkg(in_hydromod_paths_dt = hydromod_paths_dt,
                      in_sites_dt = sites_dt,
                      out_dir = file.path('results', 'gis'),
                      geom = 'points',
                      overwrite = T)
  )
  ,
  
  tar_target(
    site_snapped_gpkg_list,
    lapply(names(site_points_gpkg_list), function(in_country) {
      snap_river_sites(in_sites_path = site_points_gpkg_list[[in_country]], 
                       in_network_path = network_ssnready_gpkg_list[[in_country]],
                       custom_proj = F,
                       overwrite = T)
    }) %>% setNames(names(site_points_gpkg_list))
  ),
  
  #Subset amber river barrier dataset to only keep barriers on DRNs
  tar_target(
    barrier_points_gpkg_list,
    subset_amber(amber_path = amber_path, 
                 in_hydromod_paths_dt = hydromod_paths_dt, 
                 out_dir = file.path(resdir, 'gis'),
                 overwrite = T) 
  )
  ,
  
  tar_target(
    barrier_snapped_gpkg_list,
    lapply(names(barrier_points_gpkg_list), function(in_country) {
      snap_barrier_sites(in_sites_path = barrier_points_gpkg_list[[in_country]], 
                         in_network_path = network_ssnready_gpkg_list[[in_country]],
                         out_snapped_sites_path=NULL, 
                         custom_proj = F,
                         overwrite = T)
    }) %>% setNames(names(site_points_gpkg_list))
  )
  # ,
  # 
  # #Read hydrological modeling data for flow intermittence and discharge
  # tarchetypes::tar_map(
  #   values = hydro_combi,
  #   tar_target(
  #     hydromod_dt,
  #     get_drn_hydromod(hydromod_path = hydromod_paths_dt[country==in_country,
  #                                                        all_sims_path],
  #                      varname = in_varname,
  #                      selected_sims = 1:20)
  #   ),
  # 
  #   tar_target(
  #     hydrostats,
  #     compute_hydrostats_drn(
  #       in_network_path = network_ssnready_gpkg_list[[in_country]],
  #       in_sites_dt = sites_dt[country == in_country,],
  #       varname = in_varname,
  #       in_hydromod_drn = hydromod_dt)
  #   )
  # )
  # ,
  # 
  # #Read metadata accompanying eDNA data
  # tar_target(
  #   metadata_edna,
  #   read_xlsx(metadata_edna_path,
  #             sheet = 'metadataDNA') %>%
  #     as.data.table %>%
  #     setnames(tolower(names(.)))
  # ),
  # 
  # #Read pre-processed biological sampling data
  # tar_target(
  #   bio_dt,
  #   read_biodt(path_list = bio_data_paths,
  #              in_metadata_edna = metadata_edna)
  # ),
  # 
  # #Compute local species richness
  # tar_target(
  #   sprich,
  #   lapply(bio_dt, function(dt) {
  #     calc_sprich(in_biodt=dt,
  #                 in_metacols=metacols)
  #   }) %>%
  #     rbindlist %>%
  #     .[, running_id := paste0(site, '_', campaign)] %>%
  #     merge(env_dt, by=c('site', 'campaign', 'running_id'))
  # )
  #,
  
  
  # #Read intermittence data (90 days)
  # tar_target(
  #   interm90_dt,
  #   lapply(interm90_data_paths, fread) %>%
  #     rbindlist
  # ),
  #
  #
  
  #
  # #Merge species richness with local environmental and intermittence data
  # tar_target(
  #   alphadat_merged,
  #   merge_alphadat(in_env_dt = env_dt,
  #                  in_interm90_dt = interm90_dt,
  #                  in_sprich = sprich)
  # )
  # ,
  #
  #
  # tar_target(
  #   alpha_cor,
  #   alphadat_merged[, list(meanS_totdur90_cor = stats::cor(mean_S, TotDur90),
  #                          meanS_discharge_cor = stats::cor(mean_S, discharge)
  #   ), by=Country]
  # )
  # ,
  #
  # tar_target(
  #   alpha_cor_plots_wrap,
  #   plot_alpha_cor(alphadat_merged,
  #                  out_dir = file.path(resdir, 'Null_models'),
  #                  facet_wrap = TRUE)
  # )
  # ,
  #
  # tar_target(
  #   alpha_cor_plots_all,
  #   plot_alpha_cor(alphadat_merged,
  #                  out_dir = file.path(resdir, 'Null_models'),
  #                  facet_wrap = FALSE)
  # )
  # ,
  #
  # tar_target(
  #   lmer_S_int,
  #   compute_lmer_mods(
  #     in_dt = unique(alphadat_merged, by=c('site', 'organism')),
  #     in_yvar = 'mean_S')
  # ),
  #
  # # fungi_biof_path <- file.path(datdir, 'Datasets', 'fungi_dna_Biof_removed_zero_rows_and_columns.csv')
  # # fungi_biof <- fread(fungi_biof_path)
  # # fungi_biof[Country == 'Czech Republic', Country := 'Czech']
  # # countries_list <- unique(fungi_biof$Country)
  # # fungi_biof[1:10, 1:10]
  #
  # tar_target(
  #   null_models,
  #   lapply(
  #     names(bio_dt)[!(names(bio_dt) %in% c('bac_sedi', 'bac_biof'))],
  #     function(organism_type) {
  #       print(organism_type)
  #       organism_dt <- bio_dt[[organism_type]]
  #
  #       if (organism_type %in% c("bac_biof_nopools", "bac_sedi_nopools")) {
  #         in_nsimul = 99; in_method = 'greedyqswap'; in_thin = 100
  #       } else {
  #         in_nsimul = 999; in_method = 'quasiswap'; in_thin = 1
  #       }
  #
  #       organism_dt[, compute_null_model_inner(
  #         in_dt = .SD,
  #         in_metacols = metacols,
  #         min_siteN = 2,
  #         nsimul = in_nsimul,
  #         method = in_method,
  #         thin = in_thin)
  #         , by = Country] %>%
  #         .[, organism := organism_type]
  #     }
  #   )  %>% rbindlist
  # ),
  #
  # tar_target(
  #   env_null_models_dt,
  #   merge_env_null_models(in_null_models = null_models,
  #                         in_env = env_dt,
  #                         in_int = interm90_dt)
  # ),
  #
  # tar_target(
  #   p_z_by_stream_type,
  #   for (in_organism in unique(env_null_models_dt$organism)) {
  #     print(in_organism)
  #     plot_z_by_stream_type(in_env_null_models_dt = env_null_models_dt[organism==in_organism,],
  #                           outdir = resdir)
  #   }
  # ),
  #
  # tar_map(
  #   values = list(env_var_mapped = c('TotDur90', 'TotLeng90', 'discharge')),
  #   tar_target(
  #     p_by_env,
  #     plot_z_by_env(in_env_null_models_dt = env_null_models_dt,
  #                   env_var = env_var_mapped,
  #                   outdir = resdir),
  #
  #   )
  # ),
  #
  # tar_target(
  #   lmer_z_int,
  #   compute_lmer_mods(in_dt = env_null_models_dt,
  #                     in_yvar = 'z')
  # ),
  #
  # tar_target(
  #   p_z_jitter,
  #   plot_z_jitter_by_organism(in_env_null_models_dt = env_null_models_dt,
  #                             outdir = resdir)
  # )
)

######################## UNUSED TARGETS ########################################
#Paths to intermittence data (90 days)
# tar_target(
#   interm90_data_paths,
#   {path_list <- sapply(countries_list,
#                        function(country) {
#                          file.path(datdir, 'Datasets', 'Intermittence_Data',
#                                    paste0(country, '_Local_Interm_90_d.csv'))
#                        },
#                        USE.NAMES = TRUE)
#   path_list['France'] <- gsub('90_d', '90_d_corrected', path_list['France'])
#   return(path_list)
#   }
# ),

