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


#------------------------ merge_env_mod ----------------------------------------
merge_env_null_models <- function(in_null_models, in_env, in_int) {
  res <- in_null_models[, Country := str_to_title(Country)] %>%
    .[Country == 'Czech', Country := 'Czech Republic']
  
  env <- setDT(in_env) %>%
    data.table::setnames(c('site', 'DRN'),
                         c('Site', 'Country')) %>%
    .[Country == 'Czech', Country := 'Czech Republic']
  
  int <- setDT(in_int) %>%
    data.table::setnames('Sites', 'Site')
  
  int90 <- int[, list(TotDur90 = mean(TotDur, na.rm=T),
                      TotLeng90 = mean(TotLeng, na.rm=T)), by=.(Site, Country)]
  env_mean <- env[, list(discharge = mean(discharge_l_s, na.rm=T),
                         moss = mean(moss_cover, na.rm=T),
                         particle_size = mean(particle_size, na.rm=T)
  )
  , by=.(Site, Country, stream_type)]
  
  env_mods_dt <- merge(res, int90, by=c("Site", "Country"), all.x=T, sort=F) %>%
    merge(env_mean, by=c("Site", "Country"), all.x=T, sort=F) %>%
    .[Site != 'BUK52', ] %>% #Intermittence indicators are missing here 
    .[, Country := as.factor(Country)] %>%
    .[, significance := fifelse(pval <= 0.05, "sig", "nonsig")]
  
  return(env_mods_dt)
}

#-------------------------- plot_z_by_stream_type ------------------------------
plot_z_by_stream_type <- function(in_env_null_models_dt, outdir) {
  plots <- list()
  for(i in levels(in_env_null_models_dt$Country)){
    d <- subset(in_env_null_models_dt, Country == i)
    plots[[paste0(i)]] <- ggplot(d, aes(x=z, y=Site, color=stream_type)) + 
      scale_colour_manual(values = c("steelblue","orange")) +
      geom_boxplot() + 
      coord_flip() +
      ggtitle(paste(i)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  pdf(file.path(outdir, "Null_models", 
                paste0(unique(in_env_null_models_dt$organism), "_SES_vs_streamtype.pdf")), 
      height=10, width=10)
  do.call('grid.arrange', c(plots))
  dev.off()
}

#------------------------- plot_z_by_env ----------------------------------------
# in_env_null_models_dt <- tar_read(env_null_models_dt)
# env_var <- 'TotDur90'
# outdir <- resdir 

plot_z_by_env <- function(in_env_null_models_dt, env_var, outdir) {
  #### scatterplot, all countries in the same plot, separately for each organism group ####
  for (in_organism in unique(in_env_null_models_dt$organism)) {
    print(in_organism)
    ss <- in_env_null_models_dt[organism == in_organism] %>%
      .[, pval_lm := coef(summary(lm(z~get(env_var), data=.SD)))[2,4], 
        by=Country]
    ss2 <- ss[significance == "sig",]
    
    p1 <- ggplot() + 
      geom_point(aes(get(env_var), z, colour=Country), data = ss, size = 1, alpha=0.2) +
      geom_point(aes(get(env_var), z, colour=Country), data = ss2, size = 1, alpha=0.5) +
      geom_smooth(data=ss[pval_lm <= 0.05,], 
                  aes(get(env_var), z, colour=Country), 
                  method = "lm", linewidth = 0.75, se = F) + 
      geom_smooth(data=ss[pval_lm > 0.05,], 
                  aes(get(env_var), z, colour=Country), 
                  method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
      geom_smooth(aes(get(env_var), z), data=ss, 
                  colour="black", method = "lm", linewidth = 1.1, se = F) + # , linetype="dashed"
      theme_classic() +
      ggtitle(in_organism) + 
      xlab(env_var) +
      ylab("z") +
      scale_color_manual(values = c("Croatia" = "#ef476f",
                                    "Czech Republic" = "#f78c6b", 
                                    "Finland" = "#ffd166", 
                                    "France" = "#06d6a0", 
                                    "Hungary" = "#118ab2",
                                    "Spain" = "#073b4c")) +
      theme(legend.position = "none")
    
    out_path <- file.path(outdir, "Null_models", 
                          paste0("All_", in_organism, "_SES_vs_",
                                 env_var, "_lm2.pdf"))
    pdf(out_path, height=3, width=4)
    print(p1)
    dev.off()
  }
}

#------ compute_lmer_mods ------------------------------------------------------
compute_lmer_mods <- function(in_env_null_models_dt) {
  lmer_int <- in_env_null_models_dt[
    ,  list(
      TotDur90_full = list(lmer(z ~ TotDur90 + (1|Country), data=.SD)),
      TotDur90_null = list(lmer(z ~ (1|Country), data=.SD)),
      TotDur90_ML = list(anova(lmer(z ~ TotDur90 + (1|Country), data=.SD),
                               lmer(z ~ (1|Country), data=.SD))),
      TotLeng90_full = list(lmer(z ~ TotLeng90 + (1|Country), data=.SD)),
      TotLeng90_null = list(lmer(z ~ (1|Country), data=.SD)),
      TotLeng90_ML = list(anova(lmer(z ~ TotLeng90 + (1|Country), data=.SD),
                                lmer(z ~ (1|Country), data=.SD)))
    )
    , by=organism] 
  return(lmer_int)
}

#------ plot_z_jitter_by_organism ----------------------------------------------
plot_z_jitter_by_organism_inner <- function(in_dt) {
  jitter_p <- ggplot(in_dt) + 
    scale_y_continuous() +
    geom_jitter(aes(Country, z, colour=significance), 
                shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    labs(y = "z") +
    ggtitle(unique(in_dt$organism)) +
    scale_color_manual(values = c("sig" = "mediumblue",
                                  "nonsig" = "lightslateblue")) +
    labs(colour = "Departure from zero") +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1),
          legend.position = 'non')
  return(jitter_p)
}

plot_z_jitter_by_organism <- function(in_env_null_models_dt, outdir) {
  plots_jitter <- list()
  for (in_organism in unique(in_env_null_models_dt$organism)) {
    print(in_organism)
    plots_jitter[[in_organism]] <- plot_z_jitter_by_organism_inner(
      in_dt = in_env_null_models_dt[organism==in_organism,])
  }
  
  pdf(file.path(outdir, "Null_models", "Jitterplots_SES_significance.pdf"),
      height=10, width=6)
  do.call("grid.arrange", c(plots_jitter, ncol=2))
  dev.off()
}


