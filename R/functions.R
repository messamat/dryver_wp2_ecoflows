#-------------- workflow functions ---------------------------------------------
read_biodt <- function(path_list) {
  #Read and name all data tables
  dt_list <- mapply(function(in_path, in_name) {
    fread(in_path) %>% 
      .[, organism := in_name]},
    path_list, names(path_list))
  
  #Add Campaign and Site to bacteria data
  dt_list$bac_sedi[, c('Site', 'Campaign') := tstrsplit(V1, '_')]
  dt_list$bac_biof[, c('Site', 'Campaign') := tstrsplit(V1, '_')]
  
  return(dt_list)
}

calc_sprich <- function(in_biodt, in_metacols) {
  metacols_sub <- names(in_biodt)[names(in_biodt) %in% in_metacols]
  biodt_melt <- melt(in_biodt, id.vars = metacols_sub)
  biodt_sprich <- biodt_melt[, list(S=sum(value>0)), 
                             by=.(Site, Campaign, organism)] %>%
    .[, mean_S := mean(S), by=Site]
  return(biodt_sprich)
}


# 
# 
# com <- dia_sed[,c(3,8:1301)]
# meta <- dia_sed[,1:7]
# # first column as row names
# com <- data.frame(com[,-1], row.names=com$running_id)
# 
# # diversity indices for running_id
# dia_sed_S <- data.frame(specnumber(com))
# dia_sed_S <- cbind(dia_sed_S, meta)
# # mean diversity for each site
# dia_sed_S_mean<- dia_sed_S %>%
#   group_by(Site) %>%
#   summarise(mean_S=mean(specnumber.com., na.rm = T))
# dia_sed_S_mean$Organism <- "Sediment_diatoms"
# 
# rm(dia_sed, dia_sed_S)
# rm(com2)