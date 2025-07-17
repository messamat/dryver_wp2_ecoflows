#Small changes to implement before next run
#L22 functions: capitalize S in spain
#L56 and L158 functions: catch metacols regardless of capitalization
#Impute GEN04, campaign 1 date
#Correct avg velocity 25 m/s

library(rprojroot)
rootdir <- rprojroot::find_root(has_dir('R'))
setwd(rootdir)

source('R/packages.R')
source("R/functions.R")
source("R/SpaTemp_function_Mathis_edit.R")
# if (!file.exists("bin/03_diversity_metrics.R")) {
#   download.file(url = "https://github.com/LysandreJ/dryver/blob/main/Script/03_diversity_metrics.R",
#                 destfil = file.path("bin", "03_diversity_metrics.R")
#   )
# }
source("bin/03_diversity_metrics.R") #From 

hydromod_present_dir <- file.path('data', 'wp1', 'Results_present_period_final')
wp1_data_gouv_dir <- file.path('data', 'wp1', 'data_gouv') #Official WP1 data published later, with fuller attributes (upstream drainage area)

bio_dir <- file.path('data', 'wp2', '01_WP2 final data')
#datdir <- file.path('data', 'data_annika')
resdir <- 'results'

#Set more efficient data storage
tar_option_set(format = "qs")

metacols <- c('campaign', 'site', 'running_id', 'country', 'date',
              'summa_sample', 'sample.id', 'sample_id', 'organism')

drn_dt <- data.table(
  country = c("Croatia", "Czech", "Finland", "France",  "Hungary", "Spain"),
  catchment = c("Butiznica", "Velicka", "Lepsamaanjoki", "Albarine", "Bukkosdi", "Genal"),
  color = c("#8c510a", "#bf812d", "#01665e", "#80cdc1", "#8073ac", "#543005")
)

hydro_combi <- expand.grid(
  in_country = drn_dt$country,
  in_varname =  c('isflowing', 'qsim'),
  stringsAsFactors = FALSE)

if (!interactive()) {
  perf_ratio <- 0.7 #Set how much you want to push your computer (% of cores and RAM)
  nthreads <- round(parallel::detectCores(logical=F)*perf_ratio)
  future::plan("future::multisession", workers=nthreads)
  total_ram <- memuse::Sys.meminfo()$totalram@size*(10^9) #In GiB #ADJUST BASED ON PLATFORM
  options(future.globals.maxSize = perf_ratio*total_ram)
  # tar_option_set(controller = crew_controller_local(workers = nthreads)) #Set up parallel computing in targets
}

#--------------------------  Define targets plan -------------------------------
#-------------- Preformatting targets ------------------------------------------
preformatting_targets <- list(
  ##############################################################################
  ### DEFINE PATHS #############################################################
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
    file.path(rootdir, 'data', 'data_annika', 'metadata_DRYvER_eDNA_updated.csv'),
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
      bac_sedi = "bacteria_dna_sediment_removed_pools_and_zero_species.csv",
      bac_biof = "bacteria_dna_biofilm_removed_pools_and_zero_species.csv",
      miv = "DRYvER_MIV_data_genus_reduced_taxa.csv",
      miv_nopools = "DRYvER_MIV_data_genus_reduced_taxa_removed_pools.csv",
      miv_nopools_flying = "DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_flying.csv",
      miv_nopools_nonflying = "DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_nonflying.csv"
    ) %>%
      lapply(function(path) file.path(rootdir, 'data', 'data_annika', path))
    #, format='file'
  )
  ,
  
  #Path to DEM for Genal catchment in Spain (basin area is missing)
  tar_target(
    dem_genal_rediam_path,
    file.path(rootdir, 'data', 'dem_genal', 'rediam', 'MDT_2010_11_AND.tif'),
    format='file'
  ),
  
  
  ##############################################################################
  ### DOWNLOAD DATA ############################################################
  #Download HydroSHEDS flow direction at 90 m for Europe
  tar_target(
    flowdir_hydrosheds90m_path,
    {
      out_tif <- file.path('data', 'hydrosheds', 'eu_dir_3s', 'eu_dir_3s.tif')
      if (!file.exists(file.path('data', 'hydrosheds', 'eu_dir_3s', 'eu_dir_3s.tif'))) {
        hs_dir_path <- download_unzip(
          url =  "https://data.hydrosheds.org/file/hydrosheds-v1-dir/eu_dir_3s.zip",
          out_dir = file.path('data', 'hydrosheds'), 
          out_zip=NULL)
        unzip(file.path(hs_dir_path, 'eu_dir_3s.zip'),
              exdir = file.path('data', 'hydrosheds'))
      }
      
      return(out_tif)
    }, format = 'file'
  )
  ,
  
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
  
  ##############################################################################
  ### READ DATA ################################################################
  #Read new reach data
  tar_target(
    reaches_attri,
    lapply(drn_dt$catchment, function(in_catchment) {
      file.path(wp1_data_gouv_dir, 'gis_data',
                paste0('reaches_', in_catchment, '.csv')) %>%
        fread
    }) %>% setNames(drn_dt$country)
  )
  ,
  
  
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
  
  #Read and correct river site names 
  tar_target(
    sites_dt,
    hydromod_paths_dt[, format_site_dt(in_path = sites_reachids, 
                                       in_env_dt = env_dt,
                                       in_country = country),
                      by = country] 
  ),
  
  #Read metadata accompanying eDNA data
  tar_target(
    metadata_edna,
    fread(metadata_edna_path,
          colClasses = c(rep('character', 5), 'integer', 
                         rep('character',4))) %>%
      setnames(tolower(names(.))) %>%
      setnames(c('matchingenvid', 'sampling_date'),
               c('running_id', 'date')) %>%
      .[, date := as.Date(date, "%d.%m.%Y")]
  ),
  
  #Read pre-processed biological sampling data
  tar_target(
    bio_dt,
    read_biodt(path_list = bio_data_paths,
               in_metadata_edna = metadata_edna,
               include_bacteria = T) 
  )
  ,
  
  ##############################################################################
  ### PRE-FORMAT HYDROGRAPHIC DATA #############################################
  
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
          .[!(.[['UID']] %in% c(1102, 2988, 457, 458, 462, 464)),] %>%
          #Change startpoint of disconnected line where there is a site 
          #to match downstream line (startpoint, because the line dir is reversed)
          manual_clean_croatia
        
        st_write(clean_net, out_net_path, append=F)
      }
      return(out_net_path)
    }) %>% setNames(names(network_sub_gpkg_list)) 
  )
  ,
  
  #Direct network
  tar_target(
    network_directed_gpkg_list,
    lapply(names(network_clean_gpkg_list), function(in_country) {
      outlet_uid_list <- list(Croatia = 463, Czech = 4, Finland = 682,
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
  
  #Fix complex confluences
  tar_target(
    network_nocomplexconf_gpkg_list,
    lapply(names(network_directed_gpkg_list), function(in_country) {
      fix_complex_confluences(
        rivnet_path = network_directed_gpkg_list[[in_country]], 
        max_node_shift = 5,
        out_path= file.path(resdir, 'gis', 
                            paste0(tolower(in_country)
                                   , '_river_network_nocomplexconf',
                                   format(Sys.time(), "%Y%m%d"), '.gpkg'
                            ))
      )
    }) %>% setNames(names(network_directed_gpkg_list))
  )
  ,
  
  #Compute strahler stream order
  tar_target(
    network_strahler,
    lapply(names(network_nocomplexconf_gpkg_list), function(in_country) {
      assign_strahler_order(
        in_rivnet = network_nocomplexconf_gpkg_list[[in_country]], 
        idcol = 'UID')
    }) %>% setNames(names(network_nocomplexconf_gpkg_list))
  )
  ,
  
  #Re-assign correct IDs to match with hydrological data
  tar_target(
    network_ssnready_gpkg_list,
    lapply(names(network_nocomplexconf_gpkg_list), function(in_country) {
      out_net_path <- reassign_netids(
        rivnet_path = network_nocomplexconf_gpkg_list[[in_country]], 
        strahler_dt = network_strahler[[in_country]], 
        in_reaches_hydromod_dt = reaches_dt[country==in_country,], 
        in_reaches_attri_dt = reaches_attri[[in_country]],
        outdir = file.path(resdir, 'gis'),
        country = in_country
      )
      return(out_net_path)
    }) %>% setNames(names(network_nocomplexconf_gpkg_list))
  )
  ,
  
  #Copy gpkg to shapefiles for sharing (and removing extraneous UIDs)
  tar_target(
    network_ssnready_shp_list,
    lapply(network_ssnready_gpkg_list, function(path) {
      old_cols <- c('UID', 'strahler', 'length_uid', 'cat_cor', 'from', 'to',
                    'to_reach_shpcor', 'to_reach_hydromod', 'upstream_area_net',
                    'hydromod_shpcor_match', 'geom')
      new_cols <- c('UID', 'strahler', 'length_m', 'cat', 'from', 'to',
                    'to_cat_shp', 'to_cat_mod', 'upstream_area_net',
                    'mod_match', 'geom')
      lyr <- st_read(path)[,old_cols]
      names(lyr) <- new_cols
      lyr$UID <- seq_along(lyr$UID)
      lyr$strahler <- as.integer(lyr$strahler)
      lyr$cat <- as.integer(lyr$cat)
      lyr$to_cat_shp <- as.integer(lyr$to_cat_shp)
      lyr$to_cat_mod <- as.integer(lyr$to_cat_mod)
      write_sf(lyr, paste0(tools::file_path_sans_ext(path), '.shp'))
    })
  )
  ,
  
  #Create shapefile of sampling sites (reach lines)
  tar_target(
    site_reaches_gpkg_list,
    create_sites_gpkg(in_hydromod_paths_dt = hydromod_paths_dt,
                      in_sites_dt = sites_dt,
                      out_dir = file.path('results', 'gis'),
                      geom = 'reaches',
                      in_network_path_list = network_ssnready_shp_list,
                      in_network_idcol = 'cat',
                      overwrite = T)
  )
  ,
  
  #Create shapefile of sampling sites (points)
  tar_target(
    site_points_gpkg_list,
    create_sites_gpkg(in_hydromod_paths_dt = hydromod_paths_dt,
                      in_sites_dt = sites_dt,
                      out_dir = file.path('results', 'gis'),
                      geom = 'points',
                      overwrite = T)
  )
  ,
  
  #Snap sites to corresponding reach in network
  #Note that BUT12 is located several 100 m from 
  #the corresponding reach in the corrected network
  tar_target(
    site_snapped_gpkg_list,
    lapply(names(site_points_gpkg_list), function(in_country) {
      snap_river_sites(in_sites_path = site_points_gpkg_list[[in_country]], 
                       in_network_path = network_ssnready_shp_list[[in_country]],
                       custom_proj = F,
                       in_sites_unique_id = 'site',
                       in_network_unique_id = 'UID',
                       in_sites_idcol_tomatch = 'reach_id',
                       in_network_idcol_tomatch = 'cat',
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
  
  #Snap barriers to river network
  tar_target(
    barrier_snapped_gpkg_list,
    lapply(names(barrier_points_gpkg_list), function(in_country) {
      snap_barrier_sites(in_sites_path = barrier_points_gpkg_list[[in_country]], 
                         in_network_path = network_ssnready_shp_list[[in_country]],
                         out_snapped_sites_path=NULL, 
                         in_sites_idcol = 'GUID',
                         attri_to_join = c('cat', 'UID'),
                         custom_proj = F,
                         overwrite = T)
    }) %>% setNames(names(site_points_gpkg_list))
  )
  ,
  
  tar_target(
    genal_sites_upa_dt,
    get_genal_drainage_area(in_flowdir_path = flowdir_hydrosheds90m_path, 
                            outdir = file.path(resdir, 'gis')
    )
  )
)  

################################################################################
### COMPUTE METRICS ############################################################

mapped_hydrotargets <- tarchetypes::tar_map(
  values = hydro_combi,
  #Read hydrological modeling data for flow intermittence and discharge
  tar_target(
    hydromod_hist_dt,
    get_drn_hydromod(hydromod_path = hydromod_paths_dt[country==in_country,
                                                       sel_sim_path],
                     varname = in_varname)
    #selected_sims = 1:20)
  ),
  
  #Compute hydrological statistics for a given DRN for all dates
  tar_target(
    hydrostats_sites_hist,
    compute_hydrostats_drn(
      in_network_path = network_ssnready_shp_list[[in_country]],
      in_sites_dt = sites_dt[country == in_country,],
      varname = in_varname,
      in_hydromod_drn = hydromod_hist_dt,
      in_network_idcol = 'cat')
  ),
  
  #Subset statistics to keep only those for sampling site-dates combinations
  tar_target(
    hydrostats_sites_tsub,
    subset_hydrostats(hydrostats_sites_hist, 
                      in_bio_dt = bio_dt)
  )
)

combined_hydrotargets <- list(
  tar_combine(
    hydromod_comb_hist,
    mapped_hydrotargets[['hydromod_hist_dt']],
    command = list(!!!.x)
  ),
  tar_combine(
    hydrostats_sites_tsub_comb,
    mapped_hydrotargets[['hydrostats_sites_tsub']],
    command = list(!!!.x)
  )
)

analysis_targets <- list(
  #Prepare data for STcon
  tar_target(
    preformatted_data_STcon,
    lapply(names(network_ssnready_shp_list), function(in_country) {
      #print(in_country)
      prepare_data_for_STcon(
        in_hydromod_drn = hydromod_comb_hist[[paste0(
          "hydromod_hist_dt_", in_country, '_isflowing')]], 
        in_net_shp_path = network_ssnready_shp_list[[in_country]]
      )
    }) %>% setNames(names(network_ssnready_shp_list))
  )
  ,
  
  #Compute directed STcon 
  tar_target(
    STcon_directed_list,
    {
      out_STcon_list <- future_lapply(
        names(preformatted_data_STcon), function(in_country) {
          print(in_country)
          
          unique_sampling_dates <- lapply(bio_dt, function(org_dt) {
            org_dt[country==in_country, .(date)]
          }) %>% rbindlist %>% unique
          in_dates <- unique_sampling_dates
          
          window_size_list <- c(10, 30, 90, 365) #30, 60, 90, 180, 365
          lapply(window_size_list, function(in_window) {
            print(in_window)
            if (in_window == 365) {
              in_output <- 'all'
            } else {
              in_output <- 'STcon'
            }
            
            compute_STcon_rolling(
              in_preformatted_data = preformatted_data_STcon[[in_country]],
              ref = FALSE,
              #in_nsim = hydromod_paths_dt[country == in_country,]$best_sim,
              in_dates = in_dates,
              window = in_window,
              output = in_output,
              direction = 'directed',
              routing_mode = 'in',
              weighting = TRUE,
              rounding_factor = 1)
          }) %>% setNames(paste0('STcon_m', window_size_list)) #60, 90, 180,
        })
      setNames(out_STcon_list, names(preformatted_data_STcon))
    }
  )
  ,
  
  #Compute undirected STcon 
  tar_target(
    STcon_undirected_list, 
    {
      out_STcon_list <- future_lapply(
        names(preformatted_data_STcon), function(in_country) {
          print(in_country)
          
          unique_sampling_dates <- lapply(bio_dt, function(org_dt) {
            org_dt[country==in_country, .(date)]
          }) %>% rbindlist %>% unique
          in_dates <- unique_sampling_dates
          
          window_size_list <- c(10, 30, 90, 365) #, 30, 90, 120, 180, 365) 
          lapply(window_size_list, function(in_window) { 
            print(in_window)
            if (in_window == 365) {
              in_output <- 'all'
            } else {
              in_output <- 'STcon'
            }
            
            compute_STcon_rolling(
              in_preformatted_data = preformatted_data_STcon[[in_country]], 
              ref = FALSE,
              #in_nsim = hydromod_paths_dt[country == in_country,]$best_sim, 
              in_dates = in_dates, 
              window = in_window, 
              output = in_output,
              direction = 'undirected',
              routing_mode = 'all',
              weighting = TRUE,
              rounding_factor = 1)
          }) %>% setNames(paste0('STcon_m', window_size_list)) #60, 90, 180, 
        })
      setNames(out_STcon_list, names(preformatted_data_STcon))
    }
  )
  ,
  
  #Format directed STcon for subsequent use in models
  tar_target(
    STcon_directed_formatted,
    future_lapply(names(STcon_directed_list), function(in_country) {
      postprocess_STcon(in_STcon = STcon_directed_list[[in_country]],
                        in_net_shp_path = network_ssnready_shp_list[[in_country]],
                        standardize_STcon = FALSE, in_STcon_ref = NULL)
    }) %>% setNames(names(STcon_directed_list))
  ),
  
  #Format undirected STcon for subsequent use in models
  tar_target(
    STcon_undirected_formatted,
    future_lapply(names(STcon_undirected_list), function(in_country) {
      postprocess_STcon(in_STcon = STcon_undirected_list[[in_country]],
                        in_net_shp_path = network_ssnready_shp_list[[in_country]],
                        standardize_STcon = FALSE, in_STcon_ref = NULL) 
    }) %>% setNames(names(STcon_undirected_list))
  ),
  
  #Compute distance to nearest wet site for each reach and time step
  tar_target(
    Fdist_undirected,
    future_lapply(names(preformatted_data_STcon), function(in_country) {
      compute_Fdist(
        sites_status_matrix = preformatted_data_STcon[[in_country]]$sites_status_matrix,
        network_structure = preformatted_data_STcon[[in_country]]$network_structure, 
        routing_mode = 'all', 
        raw_dist_matrix = preformatted_data_STcon[[in_country]]$river_dist_mat, 
        in_net_shp_path = network_ssnready_shp_list[[in_country]]
      ) %>%
        compute_Fdist_rolling( #Compute average and max distance within rolling windows
          in_sites_dt=as.data.table(vect(site_snapped_gpkg_list[[in_country]]))) %>%
        melt(id.vars=c('date', 'UID'), #Convert to long format
             measure.vars = grep('Fdist_', names(.), value=T),
             value.name = 'fdist_value')
    }) %>% setNames(names(preformatted_data_STcon))
  ),
  
  #Compute distance to nearest wet site for each reach and time step
  #To avoid issues with infinites when reach has no wet upstream segments, 
  #convert infinite to twice the maximum distance to wet upstream segment otherwise
  tar_target(
    Fdist_directed,
    future_lapply(names(preformatted_data_STcon), function(in_country) {
      compute_Fdist(
        sites_status_matrix = preformatted_data_STcon[[in_country]]$sites_status_matrix,
        network_structure = preformatted_data_STcon[[in_country]]$network_structure, 
        routing_mode = 'in', 
        raw_dist_matrix = preformatted_data_STcon[[in_country]]$river_dist_mat, 
        in_net_shp_path = network_ssnready_shp_list[[in_country]]
        ) %>%
        .[is.infinite(Fdist), #Replace infinite Fdist with twice max otherwise
          Fdist := .[!is.infinite(Fdist), 2*max(Fdist)]] %>%
        compute_Fdist_rolling( #Compute average and max distance within rolling windows
          in_sites_dt=as.data.table(vect(site_snapped_gpkg_list[[in_country]]))) %>%
        melt(id.vars=c('date', 'UID'), #Convert to long format
             measure.vars = grep('Fdist_', names(.), value=T),
             value.name = 'fdist_value') 
    }) %>% setNames(names(preformatted_data_STcon))
  ),
  
  #Compile hydrostats and connectivity data together as dt
  tar_target(
    hydrocon_sites_compiled,
    lapply(drn_dt$country, function(in_country) {
      compile_hydrocon_sites_country(hydrostats_sites_tsub_comb,
                                     STcon_directed_formatted,
                                     STcon_undirected_formatted,
                                     Fdist_directed,
                                     Fdist_undirected,
                                     site_snapped_gpkg_list, 
                                     in_country)
    }) %>% rbindlist 
  )
  ,
  
  tar_target(
    hydrocon_sites_summarized,
    summarize_sites_hydrocon(
      in_hydrocon_compiled = hydrocon_sites_compiled
      #date_range = hydrocon_sites_compiled[, c(min(date), min(date)+years(1))])
    )
  )
  ,
  
  #Compute summary hydrostats for each reach in network 
  tar_target(
    hydrostats_net_hist,
    lapply(names(network_ssnready_shp_list), function(in_country) {
      summarize_network_hydrostats(
        in_country = in_country,
        in_hydromod = hydromod_comb_hist,
        in_all_date_range = c(as.Date('1960-10-01', '%Y-%m-%d'),
                              as.Date('2021-10-01', '%Y-%m-%d')),
        in_samp_date_range = hydrocon_sites_compiled[, range(date)] #c(min(date), min(date)+years(1))]
      )
    }) %>% rbindlist
  )
  ,
  
  
  tar_target(
    hydrostats_net_proj,
    lapply(list.files(file.path(wp1_data_gouv_dir, 'projections'), 
                      pattern='flowstate.*2015[-]2100',
                      full.names = TRUE),
           function(in_path) {
             print(in_path)
             summarize_drn_hydroproj_stats(hydroproj_path=in_path)
           }) %>% 
      rbindlist %>%
      .[catchment == "Lepsamanjoki", catchment := "Lepsamaanjoki"] %>%
      merge(drn_dt[, .(country, catchment)], by='catchment')
  )
  ,
  
  tar_target(
    env_summarized,
    summarize_env(in_env_dt = env_dt)
  )
  ,
  
  #Compute local taxonomic diversity
  tar_target(
    spdiv_local,
    lapply(names(bio_dt), function(in_org) {
      print(in_org)
      bio_dt[[in_org]][, calc_spdiv(in_biodt = .SD,
                                    in_metacols = metacols,
                                    level = 'local'),
                       by=country]
    }) %>%
      rbindlist
  )
  ,

  #Compute regional taxonomic diversity
  tar_target(
    spdiv_drn,
    lapply(names(bio_dt), function(in_org) {
      print(in_org)
      bio_dt[[in_org]][, calc_spdiv(in_biodt = .SD,
                                    in_metacols = metacols,
                                    level = 'regional'),
                       by=country]
    }) %>%
      rbindlist
  )
  ,
  
  #Merge ecological, environmental and hydrological data
  tar_target(
    allvars_merged,
    merge_allvars_sites(in_spdiv_local = spdiv_local, 
                        in_spdiv_drn = spdiv_drn,
                        in_hydrocon_compiled = hydrocon_sites_compiled,
                        in_hydrocon_summarized = hydrocon_sites_summarized,
                        in_env_dt = env_dt,
                        in_env_summarized = env_summarized,
                        in_genal_upa = genal_sites_upa_dt)
  )
  ,
  
  #Prepare prediction points and data for SSN
  tar_target(
    ssn_pred_pts,
    create_ssn_preds(in_network_path = network_ssnready_shp_list,
                     in_hydrostats_net_hist = hydrostats_net_hist,
                     in_hydrostats_net_proj = hydrostats_net_proj)
  )
  ,
  
  ##############################################################################
  # ANALYZE DATA
  ##############################################################################
  
  #Visualize percentage of flowing sites by DRN over time
  tar_target(
    relF_bydrn_plot,
    compare_drn_hydro(in_hydrocon_compiled = hydrocon_sites_compiled, 
                      in_sites_dt = sites_dt) 
  )
  ,
  
  tar_target(
    biof_vs_sedi_plots,
    plot_edna_biof_vs_sedi(in_allvars_merged=allvars_merged) 
  )
  ,

  #For sites x dates: Create matrices of correlations between predictors and responses, and among predictors
  tar_target(
    cor_matrices_list,
    compute_cor_matrix(allvars_merged)
  )
  ,
  
  #For data averaged by site: Create matrices of correlations between predictors and responses, and among predictors
  tar_target(
    cor_matrices_list_summarized,
    compute_cor_matrix_summarized(allvars_merged)
  )
  ,
  
  #Create matrices of correlations between predictors and responses, only for non-perennial sites
  tar_target(
    cor_matrices_list_ires, {
      allvars_merged_ires <- copy(allvars_merged)
      allvars_merged_ires$dt <- allvars_merged$dt[stream_type=='TR',]
      
      return(compute_cor_matrix(allvars_merged_ires))
    }
  )
  ,
  
  #Create correlation heatmaps across all variables
  tar_target(
    cor_heatmaps,
    plot_cor_heatmaps(in_cor_matrices = cor_matrices_list, 
                      p_threshold = 0.05)
  )
  ,
  
  tar_target(
    cor_heatmaps_summarized,
    plot_cor_heatmaps(in_cor_matrices = cor_matrices_list_summarized, 
                      p_threshold = 0.05)
  )
  ,
  
  #Check whether richness is related to habitat volumne (for miv and biofilm)
  #-> No, good
  tar_target(
    corplots_div_habvol,
    check_cor_div_habvol(in_allvars_merged = allvars_merged)
  )
  ,
  
  #All organisms: max 12
  tar_target(
    organism_list,
    c('miv_nopools', 'miv_nopools_flying', 'miv_nopools_nonflying', 
      'fun', 'dia', 'bac')
  )
  ,
  
  #Define all types of spatial covariance structure to test: 144 unique combinations
  tar_target(
    ssn_covtypes,
    expand.grid(
      c("none", "linear", "spherical", "exponential", "mariah", "epa"),
      c("none", "linear", "spherical", "exponential", "mariah", "epa"),
      c("none", "spherical", "exponential", "gaussian")) %>%
      as.data.table %>%
      setnames(c('down', 'up', 'euc')) %>%
      .[, label := paste(down, up, euc, sep='_')]
  ),
  
  ##############################################################################
  # MODEL SITES x DATE
  ##############################################################################
  #Ordinate local environmental variables to use axes in regression models
  tar_target(
    local_env_pca,
    ordinate_local_env(in_allvars_dt = allvars_merged$dt)
  ),

  #Create Spatial Stream Network (SSN) object
  # tar_target(
  #   ssn_eu,
  #   create_ssn_europe(in_network_path = network_ssnready_shp_list,
  #                     in_sites_path = site_snapped_gpkg_list,
  #                     in_allvars_dt = allvars_merged$dt,
  #                     in_local_env_pca = local_env_pca,
  #                     in_barriers_path = barrier_snapped_gpkg_list,
  #                     in_hydrostats_net_hist = hydrostats_net_hist,
  #                     in_pred_pts =  NULL,
  #                     out_dir = file.path(resdir, 'ssn'),
  #                     out_ssn_name = 'ssn_eu',
  #                     overwrite = T)
  # )
  # ,
  # 
  # #Define all hydrological variables: 16 (max ~73)
  # tar_target(
  #   hydro_vars_forssn,
  #   {
  #     hydrovar_grid <- expand.grid(
  #       c('DurD', 'FreD'), #'PDurD', 'FreD', 'PFreD', 'uQ90', 'oQ10', 'maxPQ', 'PmeanQ'
  #       paste0(c(30, 60, 90, 365, 3650), 'past')
  #     )
  #     stcon_grid <- expand.grid(paste0('STcon_m', c(30, 60, 90, 180, 365)),
  #                               c('_directed', '_undirected'))
  #     fdist_grid <- expand.grid(paste0('Fdist_mean_', c(30, 60, 90, 180, 365, 3650), 'past'),
  #                               c('_directed', '_undirected'))
  #     hydro_regex_list <- c(
  #       paste0(hydrovar_grid$Var1,'.*', hydrovar_grid$Var2),
  #       paste0(stcon_grid$Var1, stcon_grid$Var2),
  #       paste0(fdist_grid$Var1, fdist_grid$Var2)
  #     )
  # 
  #     lapply(hydro_regex_list, function(var_str) {
  #       grep(paste0('^', var_str),
  #            names(ssn_eu[[1]]$ssn$obs),
  #            value=T)
  #     }) %>% unlist
  #   }
  # )
  # ,
  # 
  # # Run a first SSN with a single hydrological variable for each organism
  # # to determine the top spatial covariance types
  # tar_target(
  #   ssn_richness_covtype,
  #   model_ssn_hydrowindow(
  #     in_ssn = ssn_eu,
  #     organism = organism_list,
  #     formula_root = '~ log10(basin_area_km2) + log10(basin_area_km2):country',
  #     hydro_var = 'DurD365past',
  #     response_var = 'richness',
  #     ssn_covtypes = ssn_covtypes
  #   ),
  #   pattern = map(organism_list),
  #   iteration = "list"
  # )
  # ,
  # 
  # #Select the top 5 covariance structures for each organism based on AIC and
  # #structure the models to run
  # tar_target(
  #   ssn_richness_models_to_run,
  #   {
  #     #Get top 5 covtypes by organism
  #     selected_covtypes <- lapply(ssn_richness_covtype, `[[`, "ssn_glance") %>%
  #       rbindlist %>%
  #       .[, .SD[order(AIC)][1:5, .(covtypes)], by = organism]
  # 
  #     #Convert to named list by organism
  #     covtypes_by_organism <- selected_covtypes[, .(covtypes = list(covtypes)), by = organism]
  # 
  #     #Build a list of model setups (one per organism × hydro_var)
  #     model_setups <- CJ(
  #       organism = covtypes_by_organism$organism,
  #       hydro_var = hydro_vars_forssn
  #     )
  # 
  #     #Attach covtypes to each setup
  #     model_setup_list <- list()
  #     for (org in unique(covtypes_by_organism$organism)) {
  #       for (hv in hydro_vars_forssn) {
  #         model_setup_list[[length(model_setup_list) + 1]] <- list(
  #           organism = org,
  #           hydro_var = hv,
  #           covtypes = covtypes_by_organism[organism==org,]$covtypes[[1]]
  #         )
  #       }
  #     }
  # 
  #     return(model_setup_list)
  #   }
  # )
  # ,
  # 
  # #Run SSN for each chosen variable and time window
  # tar_target(
  #   ssn_richness_hydrowindow,
  #   future_lapply(
  #     ssn_richness_models_to_run,
  #     function(model_setup) {
  #       print(model_setup)
  #       model_ssn_hydrowindow(
  #         in_ssn = ssn_eu,
  #         organism = model_setup$organism,
  #         formula_root = '~ log10(basin_area_km2) + log10(basin_area_km2):country',
  #         hydro_var = model_setup$hydro_var,
  #         response_var = 'richness',
  #         ssn_covtypes = ssn_covtypes[label %in% model_setup$covtypes, ]
  #       )
  #     })
  # )
  # ,
  # 
  # tar_target(
  #   ssn_covtype_selected,
  #   select_ssn_covariance(in_ssnmodels=ssn_richness_hydrowindow)
  # )
  # ,
  # 
  # tar_target(
  #   ssn_richness_hydrowindow_formatted,
  #   {
  #     ssn_model_names <- do.call(rbind, ssn_richness_models_to_run)[,1:2] %>%
  #       as.data.table
  #     ssnmodels <- cbind(ssn_model_names, ssn_richness_hydrowindow)
  # 
  #     out_list <- lapply(organism_list, function(in_organism) {
  #       print(in_organism)
  #       format_ssn_hydrowindow(in_ssnmodels = ssnmodels,
  #                              in_organism = in_organism,
  #                              in_covtype_selected = ssn_covtype_selected)
  #     })
  #     names(out_list) <- organism_list
  # 
  #     return(out_list)
  #   }
  # )
  # ,
  ##############################################################################
  # MODEL SITES SUMMARIZED
  ##############################################################################
  #Ordinate local environmental variables to use axes in regression models
  tar_target(
    local_env_pca_summarized,
    ordinate_local_env(in_allvars_dt = allvars_merged$dt_summarized)
  )
  ,
  
  #Create Spatial Stream Network (SSN) object
  tar_target(
    ssn_eu_summarized,
    create_ssn_europe(in_network_path = network_ssnready_shp_list,
                      in_sites_path = site_snapped_gpkg_list,
                      in_allvars_dt = allvars_merged$dt_summarized,
                      in_local_env_pca = local_env_pca_summarized,
                      in_barriers_path = barrier_snapped_gpkg_list,
                      in_hydrostats_net_hist = hydrostats_net_hist,
                      in_pred_pts = ssn_pred_pts,
                      out_dir = file.path(resdir, 'ssn'),
                      out_ssn_name = 'ssn_eu_summarized',
                      overwrite = T)
  )
  ,
  
  tar_target(
    ssn_summarized_maps,
    map_ssn_summarized(in_ssn_summarized = ssn_eu_summarized,
                      in_allvars_merged = allvars_merged,
                      selected_organism_list = organism_list)
  )
  ,
  
  tar_target(
    ssn_mods_miv_yr,
    model_miv_yr(in_ssn_eu_summarized = ssn_eu_summarized,
                 in_allvars_merged = allvars_merged,
                 in_cor_matrices = cor_matrices_list_summarized, 
                 ssn_covtypes = ssn_covtypes)
  )
)


#Run for other diversity indices
  
  # combined_ssntargets <- list(
  #   tar_combine(
  #     ssnmodels_combined,
  #     mapped_ssntargets[["ssn_richness_hydrowindow"]],
  #     command = list(!!!.x)
  #   )
  # )
  

list(preformatting_targets, mapped_hydrotargets, 
     combined_hydrotargets, analysis_targets) %>%
  unlist(recursive = FALSE)



######################### OLD TARGETS ##########################################
#   tar_target(
#     alpha_cor_plots_wrap,
#     plot_alpha_cor(alphadat_merged,
#                    out_dir = file.path(resdir, 'Null_models'),
#                    facet_wrap = TRUE)
#   )
#   ,
# 
# tar_target(
#   alpha_cor_plots_all,
#   plot_alpha_cor(alphadat_merged,
#                  out_dir = file.path(resdir, 'Null_models'),
#                  facet_wrap = FALSE)
# )

# tar_target(
#   alpha_cor,
#   alphadat_merged[, list(meanS_totdur90_cor = stats::cor(mean_S, TotDur90),
#                          meanS_discharge_cor = stats::cor(mean_S, discharge)
#   ), by=Country]
# )
# ,
#

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


######################## UNUSED TARGETS ########################################
#Plot spearman's correlation between hydro metrics by time window 
#and site-specific richness and t-minus1 betadiv
# tar_target(
#   corplots_div_hydrowindow, {
#     hydrovar_list <- c('DurD', 'PDurD', 'FreD', 'PFreD', 
#                        'uQ90', 'oQ10', 'maxPQ', 'PmeanQ',
#                        'STcon.*_directed', 'STcon.*_undirected',
#                        'Fdist.*_directed', 'Fdist.*_undirected')
#     lapply(hydrovar_list, function(in_var_substr) {
#       plot_cor_hydrowindow(in_cor_dt = cor_matrices_list$div_bydrn, 
#                            temporal_var_substr = in_var_substr, 
#                            response_var_list = c('richness', 'invsimpson','Jtm1'),
#                            colors_list = drn_dt$color,
#                            save_plot=T,
#                            out_dir = resdir)
#     }) %>% setNames(hydrovar_list)
#     
#     lapply(hydrovar_list, function(in_var_substr) {
#       plot_cor_hydrowindow(in_cor_dt = cor_matrices_list_ires$div_bydrn, 
#                            temporal_var_substr = in_var_substr, 
#                            response_var_list = c('richness','invsimpson','Jtm1'),
#                            colors_list = drn_dt$color,
#                            save_plot=T,
#                            plot_name_suffix = '_IRES',
#                            out_dir = resdir)
#     }) %>% setNames(paste0(hydrovar_list, '_IRES'))
#   }
# )
# ,
# 
#
# tar_target(
#   STcon_rolling_ref_list,
#   {
#     out_STcon_list <- future_lapply(
#       names(preformatted_data_STcon), function(in_country) {
#         print(in_country)
# 
#         unique_sampling_dates <- lapply(bio_dt, function(org_dt) {
#           org_dt[country==in_country, .(date)]
#         }) %>% rbindlist %>% unique
#         in_dates <- unique_sampling_dates
# 
#         window_size_list <-  c(10, 365) #30, 60, 90, 180,
#         lapply(window_size_list, function(in_window) {
#           print(in_window)
# 
#           if (in_window == 365) {
#             in_output <- 'all'
#           } else {
#             in_output <- 'STcon'
#           }
# 
#           compute_STcon_rolling(in_preformatted_data = preformatted_data_STcon[[in_country]],
#                                 ref = TRUE,
#                                 in_nsim = NULL,
#                                 in_dates = in_dates,
#                                 window = in_window,
#                                 output = in_output,
#                                 direction = 'directed',
#                                 routing_mode = 'in',
#                                 weighting = FALSE)
#         }) %>% setNames(paste0('STcon_m', window_size_list))
#       })
#     setNames(out_STcon_list, names(preformatted_data_STcon))
#   }
# )
# ,

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

