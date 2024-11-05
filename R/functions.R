#-------------- workflow functions ---------------------------------------------
# path_list = tar_read(bio_data_paths)
# in_metadata_edna <- tar_read(metadata_edna)

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
  
  #Standardize country names
  country_standard <- data.table(
    original = c("CRO", "FRA", "SPA", "CZ",  "HUN", "FIN"),
    new = c("Croatia", 'France', "spain", "Czech", "Hungary", 'Finland')
  )
  in_metadata_edna <- merge(in_metadata_edna, country_standard,
                            by.x='Country', by.y='original') %>%
    .[, Country := new] %>%
    .[, `:=`(new = NULL)]
                                 
  #Remove bacteria in pools
  dt_list$bac_sedi <- merge(dt_list$bac_sedi, 
                            in_metadata_edna[Sample_type=='sediment', 
                                             .(matchingEnvID, Habitat, Country)],
                            by.x='running_id', by.y='matchingEnvID')
  dt_list$bac_biof <- merge(dt_list$bac_biof,
                            in_metadata_edna[Sample_type=='biofilm',
                                             .(matchingEnvID, Habitat, Country)],
                            by.x='running_id', by.y='matchingEnvID')
  
  dt_list$bac_biof_nopools <- dt_list$bac_biof %>%
    .[Habitat!='pool',] %>%
    .[, Habitat := NULL]
  dt_list$bac_biof[, Habitat := NULL]
  
  
  dt_list$bac_sedi_nopools <- dt_list$bac_sedi %>%
    .[Habitat!='pool',] %>%
    .[, Habitat := NULL]
  dt_list$bac_sedi[, Habitat := NULL]
  
  return(dt_list)
}

#------ calc_sprich ------------------------------------------------------------
calc_sprich <- function(in_biodt, in_metacols) {
  #Get metadata columns (all except species data)
  metacols_sub <- names(in_biodt)[names(in_biodt) %in% in_metacols]
  biodt_melt <- melt(in_biodt, id.vars = metacols_sub)
  #Compute species richness
  biodt_sprich <- biodt_melt[, list(S=sum(value>0)), 
                             by=.(Site, Campaign, organism)] %>%
    .[, mean_S := mean(S), by=Site]
  return(biodt_sprich)
}

#------ format_envinterm -------------------------------------------------------
# in_env_dt <- tar_read(env_dt)
# in_interm90_dt <- tar_read(interm90_dt)
# in_sprich = tar_read(sprich)

merge_alphadat <- function(in_env_dt, in_interm90_dt, in_sprich) {
  #Compute mean 90-day drying duration and event length
  interm90_mean <- in_interm90_dt[
    , list(totdur90 = mean(TotDur, na.rm=T),
           totleng90 = mean(TotLeng, na.rm=T)),
    by=Sites] %>%
    setnames('Sites', 'site')

  #Average environmental variables
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

#------- plot_alpha_cor --------------------------------------------------------
# tar_load(alphadat_merged)

plot_alpha_cor_inner <- function(in_alphadat_merged_organism, x_var) {
  #Compute simple linear regression
  in_alphadat_merged_organism[
    , lm_pval_ltype := fifelse(
      coef(summary(lm(mean_S~get(x_var))))[2,4] < 0.05,
      'solid', 'dashed'
    ),
    by=Country] 
  
  alpha_plots <- ggplot(in_alphadat_merged_organism, aes(x=get(x_var), y=mean_S)) + 
    geom_point(size = 2) + 
    geom_smooth(aes(linetype=lm_pval_ltype), method='lm', linewidth = 0.5, se = F) +
    scale_linetype_identity() +
    labs(x=x_var) +
    facet_wrap(~Country) +
    theme_classic()
  
  return(alpha_plots)
}  
  
plot_alpha_cor <- function(in_alphadat_merged, out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  #Plot relationship between each organism alpha div and totdur90 for each country
  lapply(unique(in_alphadat_merged$organism), function(org_sel) {
    ggsave(
      filename = file.path(out_dir, paste0(org_sel, '_mean_S_vs_totdur90_lm_sig.png')),
      plot = plot_alpha_cor_inner(in_alphadat_merged[organism == org_sel,],
                                  x_var = 'totdur90') + scale_x_log10(),
      width=10, height=10
    )
  })
  
  #Plot relationship between each organism alpha div and discharge for each country
  lapply(unique(in_alphadat_merged$organism), function(org_sel) {
    ggsave(
      filename = file.path(out_dir, paste0(org_sel, '_mean_S_vs_discharge_lm_sig.png')),
      plot = plot_alpha_cor_inner(in_alphadat_merged[organism == org_sel,],
                                  x_var = 'discharge') + scale_x_log10(),
      width=10, height=10
    )
  })
}

#------ compute_null_model_inner -----------------------------------------------
compute_null_model_inner <- function(in_dt, 
                                     in_metacols,
                                     min_siteN, 
                                     method,
                                     thin,
                                     nsimul=999,
                                     in_dist='jac') {
  
  #Remove sampling site/dates with 0 abundance across all species
  sp_cols <- names(in_dt)[!(names(in_dt) %in% in_metacols)]
  nonnull_sites <- in_dt[,rowSums(.SD), .SDcols=sp_cols]>0
  in_dt_sub <- in_dt[nonnull_sites,] 
  #Remove rate sites
  nonrare_sites <- in_dt_sub[, .N, by=Site][N >= min_siteN, Site]
  in_dt_sub <- in_dt_sub[Site %in% nonrare_sites,]
  site_factorID <- as.factor(in_dt_sub$Site) #Get factor for each site
  
  #Only keep species that occur in at least one of the sites, and sites with at least one species
  nonnull_cols <- in_dt_sub[, lapply(.SD, sum)>0, .SDcols=sp_cols] %>%
    names(.)[.]
  com <- in_dt_sub[, nonnull_cols, with=F]
  
  #Compute null models and compare to simulations
  foo <-function(x, groups, ...) {
    diag(meandist(vegdist(x, in_dist, binary =TRUE), grouping = groups))
  }
  oecosimu_out <- oecosimu(comm = com, nestfun = foo, method = method, 
                           thin = thin, nsimul = nsimul, groups = site_factorID)

  
  #Format null model outputs
  oecosimu_dt <- as.data.table(oecosimu_out$oecosimu)[
    , .(statistic, z, means, pval)] %>%
    .[, `:=`(Site = names(oecosimu_out$statistic),
             dist = in_dist
             )]
  
  return(oecosimu_dt)
}
