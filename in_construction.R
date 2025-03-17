library(data.table)
library(magrittr)
library(MASS) # For boxcox

# Function to fill NAs hierarchically
fill_nas_hierarchical <- function(dt, cols_to_fill, site_col, country_col) {
  
  # Create a copy to avoid modifying the original data.table in place
  dt_copy <- copy(dt)
  
  # 1. Fill NAs with site averages
  dt_copy[,
          (cols_to_fill) := lapply(.SD, function(x) {
            fifelse(is.na(x), mean(x, na.rm = TRUE), x)
          }),
          .SDcols = cols_to_fill,
          by = site_col
  ]
  
  # 2. Fill remaining NAs with country averages
  dt_copy[,
          (cols_to_fill) := lapply(.SD, function(x) {
            fifelse(is.na(x), mean(x, na.rm = TRUE), x)
          }),
          .SDcols = cols_to_fill,
          by = country_col
  ]
  
  # 3. Fill any remaining NAs with overall averages
  overall_means <- dt_copy[, lapply(.SD, mean, na.rm = TRUE),
                           .SDcols = cols_to_fill]
  for (col in cols_to_fill) {
    dt_copy[, (col) := fifelse(is.na(get(col)), overall_means[[col]], get(col))]
  }
  
  return(dt_copy)
}

# Function for Box-Cox Transformation, Z-standardization, and PCA
trans_pca <- function(in_dt, in_cols_to_ordinate, num_pca_axes = 4) {
  dt_copy <- copy(in_dt)
  
  # 1. Determine optimal lambda for Box-Cox
  dt_copy[, (in_cols_to_ordinate) := lapply(.SD, function(x) {
    min_positive <- min(x[x > 0], na.rm = TRUE)
    trans_shift <- if (is.finite(min_positive)) min_positive * 0.01 else 1
    return(x + trans_shift)
  }), .SDcols = in_cols_to_ordinate]
  
  bc_lambdas <- dt_copy[, lapply(.SD, forecast::BoxCox.lambda),
                        .SDcols = in_cols_to_ordinate] %>%
    round(1)
  
  # 2. Apply Box-Cox and Z-standardization
  for (in_col in in_cols_to_ordinate) {
    dt_copy[, (in_col) := base::scale(
      forecast::BoxCox(get(in_col), 
                       bc_lambdas[[in_col]])
    )]
  } 
  
  # 3. PCA (using .SD for transformed columns only, and ensuring non-NA data)
  pca_out <- in_dt[, stats::prcomp(.SD, center=FALSE, scale=FALSE), 
                   .SDcols = in_cols_to_ordinate] 
  pca_out_dt <- as.data.table(pca_out$x) %>%
    setnames(paste0('env_', names(.)))  %>%
    .[, .SD, .SDcols = seq(1, num_pca_axes)]
  
  return(list(
    pca = pca_out,
    pca_dt = pca_out_dt,
    trans_dt = dt_copy
  ))
}

trans_pca_wrapper <- function(in_dt, in_cols_to_ordinate, id_cols, 
                              group_cols = NULL, num_pca_axes = 4) {
  if (is.null(group_cols)) {
    pca_out <- trans_pca(in_dt = in_dt, 
                         in_cols_to_ordinate = in_cols_to_ordinate, 
                         num_pca_axes = num_pca_axes)
    out_dt <- cbind(in_dt[, id_cols, with=F], pca_out$pca_dt)
    out_plot <- ordiplot(pca_out$pca, 
             choices = c(1, 2), 
             type="text", 
             display='species')
    
    return(list(
      pca = pca_out$pca,
      trans_dt = pca_out$trans_dt,
      dt = out_dt,
      plot = out_plot
    ))
  } else {
    out_dt <- cbind(
      in_dt[, trans_pca(in_dt = .SD, 
                   in_cols_to_ordinate = in_cols_to_ordinate, 
                   num_pca_axes = num_pca_axes)$pca_dt
       , by=group_cols],
      in_dt[, .SD, .SDcols = id_cols, by=group_cols][, id_cols, with=F]
    )
    
    return(list(
      pca = NULL,
      trans_dt = NULL, #Could be easily done
      dt = out_dt,
      plot = NULL
    ))
  }
}
  

#------------ Analysis ---------------------------------------------------
tar_load(allvars_merged)
dt <- allvars_merged$dt[, n_nas := rowSums(is.na(.SD))] %>%
  setorderv(c('country', 'site', 'campaign', 'n_nas')) %>%
  .[!duplicated(paste0(site, campaign)), ] %>%
  .[!(site == 'GEN04' & campaign == 1),]
cols_to_transform <- setdiff(allvars_merged$cols$env, 
                             c(allvars_merged$cols$exclude_cols,
                               allvars_merged$cols$group_cols))
by_group = NULL
num_pca_axes = 2

in_allvars_merged <- tar_read(allvars_merged)

reduce_env <- function(in_allvars_merged) {
  #1. Compute PCA for miv_nopools ----------------------------------------------
  env_cols_miv <- c('avg_velocity_macroinvertebrates', 'embeddedness',
                    'bedrock', 'particle_size', 'oxygen_sat', 'filamentous_algae',
                    'incrusted_algae', 'macrophyte_cover', 'leaf_litter_cover',
                    'moss_cover', 'wood_cover', 'riparian_cover_in_the_riparian_area',
                    'shade', 'hydromorphological_alteration', 'm2_biofilm',
                    'conductivity_micros_cm', 'oxygen_mg_l', 'ph')
  
  #Convert all columns to numeric (rather than integer)
  in_allvars_merged$dt[, (env_cols_miv) := lapply(.SD, as.numeric), 
                    .SDcols = env_cols_miv]
  
  #Fill NAs hierarchically. First by site, then by country, then overall
  miv_nopools_dt <- fill_nas_hierarchical(
    dt = in_allvars_merged$dt[organism == 'miv_nopools'], 
    cols_to_fill = env_cols_miv, 
    site_col = 'site', 
    country_col = 'country')
  
  #Check distributions by country
  dt_miv_envmelt <- melt(in_allvars_merged$dt[organism == 'miv',], 
                         id.vars = c('country', 'site', 'campaign'),
                         measure.vars = env_cols_miv)
  ggplot(dt_miv_envmelt[variable=='conductivity_micros_cm',],
         aes(x=country, y=value)) +
    geom_jitter() + 
    facet_wrap(~variable, scales='free_y') +
    scale_y_log10()
  
  ggplot(dt_miv_envmelt, aes(x=value)) +
    geom_density() + 
    facet_wrap(~variable, scales='free')
  
  out_list_miv <- trans_pca_wrapper(in_dt = miv_nopools_dt, 
                    in_cols_to_ordinate = env_cols_miv, 
                    id_cols = c('site', 'date', 'country'), 
                    group_cols = NULL, 
                    num_pca_axes = 4)
  
  out_list_miv_country <- trans_pca_wrapper(in_dt = miv_nopools_dt, 
                                in_cols_to_ordinate = env_cols_miv, 
                                id_cols = c('site', 'date'), 
                                group_cols = 'country', 
                                num_pca_axes = 4)
  
  #1. Compute PCA for microbes -------------------------------------------------
  
  return(out_list)
}
