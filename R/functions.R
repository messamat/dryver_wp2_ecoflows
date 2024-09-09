#-------------- workflow functions ---------------------------------------------
#path_list = tar_read(bio_data_paths)

#------ read_biodt -------------------------------------------------------------
read_biodt <- function(path_list, in_metadata_edna) {
  #Read and name all data tables
  dt_list <- mapply(function(in_path, in_name) {
    fread(in_path) %>% 
      .[, organism := in_name]},
    path_list, names(path_list))
  
  #Add Campaign and Site to bacteria data, then remove pool sites
  dt_list$bac_sedi[, c('Site', 'Campaign') := tstrsplit(V1, '_')] %>%
    setnames('V1', 'running_id')
  dt_list$bac_biof[, c('Site', 'Campaign') := tstrsplit(V1, '_')] %>%
    setnames('V1', 'running_id')
  
  #Remove bacteria in pools
  dt_list$bac_sedi_nopools <- merge(dt_list$bac_sedi, 
                                    in_metadata_edna[Sample_type=='sediment', 
                                                     .(matchingEnvID, Habitat)],
                                    by.x='running_id', by.y='matchingEnvID') %>%
    .[Habitat!='pool',] %>%
    .[, Habitat := NULL]
  
  dt_list$bac_biof_nopools <- merge(dt_list$bac_biof,
                                    in_metadata_edna[Sample_type=='biofilm',
                                                     .(matchingEnvID, Habitat)],
                                    by.x='running_id', by.y='matchingEnvID') %>%
    .[Habitat!='pool',] %>%
    .[, Habitat := NULL]
  
  return(dt_list)
}

#------ calc_sprich ------------------------------------------------------------
calc_sprich <- function(in_biodt, in_metacols) {
  metacols_sub <- names(in_biodt)[names(in_biodt) %in% in_metacols]
  biodt_melt <- melt(in_biodt, id.vars = metacols_sub)
  biodt_sprich <- biodt_melt[, list(S=sum(value>0)), 
                             by=.(Site, Campaign, organism)] %>%
    .[, mean_S := mean(S), by=Site]
  return(biodt_sprich)
}

#------ format_envinterm -------------------------------------------------------
# in_env_dt <- tar_read(env_dt)
# in_interm90_dt <- tar_read(interm90_dt)
# in_sprich = tar_read(sprich)

merge_alphadat <- function(in_env_dt, in_interm90_dt) {

  interm90_mean <- in_interm90_dt[
    , list(totdur90 = mean(TotDur, na.rm=T),
           TotLeng90 = mean(TotLeng, na.rm=T)),
    by=Sites] %>%
    setnames('Sites', 'site')

  env_mean <- in_env_dt[
    , list(discharge = mean(discharge_l_s, na.rm=T),
           moss = mean(moss_cover, na.rm=T),
           particlesize = mean(particle_size, na.rm=T)),
    by=.(DRN, site, stream_type)]
  
  #data.table::setnames(in_sprich, 'Site', 'site') #not working for some reason 
  #SET_STRING_ELT() can only be applied to a 'character vector', not a 'char'
  in_sprich[, `:=`(site = Site,
                   Site = NULL)]
  
  merged_dat <- merge(in_sprich, interm90_mean, by= 'site', all.x=T) %>%
    merge(env_mean, by='site', all.x=T) %>%
    .[site != 'BUK52',] %>% #Missing intermittence indicators according to Annika?
    .[DRN == 'Czech Republic', DRN := 'Czech'] %>%
    setnames('DRN', 'Country')
  
  return(merged_dat)
}
