library(data.table)
library(dplyr)
library(igraph)
library(magrittr)
library(MASS)
library(mvtnorm)
library(OCNet)
library(sf)
library(SSN2)
library(SSNbler)
library(targets)

in_network_path = tar_read(network_ssnready_shp_list)
in_hydrostats_net_hist = tar_read(hydrostats_net_hist)
lsn_path = "troubleshoot_ssn//lsn"
in_sites_path = tar_read(site_snapped_gpkg_list)
in_allvars_dt = tar_read(allvars_merged)$dt
out_ssn_path = 'troubleshoot_ssn//reprex_ssn.ssn'

###############################################################################
create_ssn_europe(in_network_path = network_ssnready_shp_list,
                  in_sites_path = site_snapped_gpkg_list,
                  in_allvars_dt = allvars_merged$dt,
                  in_local_env_pca = local_env_pca,
                  in_barriers_path = barrier_snapped_gpkg_list,
                  in_hydrostats_net_hist = hydrostats_net_hist,
                  in_pred_pts = NULL,
                  out_dir = file.path(resdir, 'ssn'),
                  out_ssn_name = 'ssn_eu',
                  overwrite = T)


in_miv_ssn <- SSN2::ssn_import("results/ssn/ssn_eu_miv_nopools.ssn", overwrite=TRUE)



net_eu <- SSN2::ssn_get_data(in_miv_ssn) 
net_eu <- net_eu[net_eu$country %in% c('France', 'Spain'),
                 c('rid', 'site', 'country', 'campaign', 'richness', 
                   'basin_area_km2', 'netgeom', 'qsim', 'afv_qsqrt', 'geometry')]


###############################################################################

#Read input network
net_eu <- lapply(c('France', 'Spain'), function(in_country) {
  #print(in_country)
  net_proj <- in_network_path[[in_country]] %>%
    st_cast("LINESTRING") %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') %>%
    st_transform(3035)
  
  net_hydro <- merge(net_proj, in_hydrostats_net_hist[country==in_country,],
                     by.x = 'cat', by.y = 'reach_id')
  
  return(net_hydro)
}) %>% do.call(rbind, .) 


edges_lsn <- SSNbler::lines_to_lsn(
  streams = net_eu, 
  lsn_path = lsn_path, 
  check_topology = TRUE, 
  snap_tolerance = 0.1, 
  topo_tolerance = 20, 
  overwrite = overwrite)


edges_lsn <- SSNbler::updist_edges(edges = edges_lsn, 
                                   save_local = TRUE, 
                                   lsn_path = lsn_path,
                                   calc_length = TRUE)


sites_eu <- lapply(c('France', 'Spain'), function(in_country) {
  st_read(in_sites_path[[in_country]]) %>%
    st_transform(3035)  
}) %>% 
  do.call(rbind, .) %>%
  rename(country=country_sub)

sites_lsn <- SSNbler::sites_to_lsn(sites = sites_eu, 
                                   edges = edges_lsn, 
                                   lsn_path = lsn_path,
                                   file_name = 'sites',
                                   snap_tolerance = 5,
                                   save_local = TRUE, 
                                   overwrite = overwrite)



sites_lsn_attri <- merge(sites_lsn,
                         in_allvars_dt[, c('country', 'site', 'upstream_area_net',
                                           'campaign','organism','richness', 'basin_area_km2'),
                                       with=F], 
                         by=c('country', 'site', 'upstream_area_net'),
                         all.y=F) 

sites_list <- list(sites = sites_lsn_attri)


#  Calculate upstream distance -----
edges_lsn <- updist_edges(
  edges =  edges_lsn,
  save_local = TRUE,
  lsn_path = lsn_path,
  calc_length = TRUE
)

sites_list_lsn <- updist_sites(
  sites = sites_list,
  edges = edges_lsn,
  length_col = "Length",
  save_local = TRUE,
  lsn_path = lsn_path
)
# ---------------------- Compute AFVs ---------------------- #
edges_lsn$qsim_avg_sqrt <- sqrt(edges_lsn$qsim_avg)

edges_lsn <- afv_edges(
  edges = edges_lsn,
  infl_col = "qsim_avg_sqrt",
  segpi_col = "pi_qsqrt",
  afv_col = "afv_qsqrt",
  lsn_path = lsn_path
)

sites_list_lsn <- afv_sites(
  sites = sites_list_lsn,
  edges = edges_lsn,
  afv_col = "afv_qsqrt",
  save_local = TRUE,
  lsn_path = lsn_path
)
# ---------------------- Assemble SSN ---------------------- #
out_ssn <- ssn_assemble(
  edges = edges_lsn,
  lsn_path = lsn_path,
  obs_sites = sites_list_lsn$sites[sites_list_lsn$sites$organism == 'miv_nopools',],
  ssn_path = out_ssn_path,
  import = TRUE,
  check = TRUE,
  afv_col = "afv_qsqrt",
  overwrite = overwrite
)