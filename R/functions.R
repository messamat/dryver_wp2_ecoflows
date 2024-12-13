#-------------- utility functions ----------------------------------------------
#------ hsaHydroYearSeasons ----------------------------------------------------
# Compute a data.frame containing the hydrological year, and 
# optionnally the days of the hydrological year and the season which
# each day belongs to. (developped by Ivan Horner)

# dates:        a date/time vector (will be coerced to date using 
#               as.Dates())
# month:        month defining the start of a hydrological year
# day_of_the_year: should the days of the (hydrological) year be
#               computed?
# seasons:      list defining the seasons in month number (1:12)
#               (or days if 'seasons_by_day' is TRUE).
# seasons_by_day: are seasons defined by days (of the year)?
# minimal:      should only the minimal desired results be returned
#               or intermediate results as well?
hsaHydroYearSeasons <- function(dates, month = 9, day_of_the_year = TRUE, 
                                seasons = list("SummerFall" = 5:10, "WinterSpring" = c(11:12, 1:4)),
                                seasons_by_day = FALSE, minimal = FALSE) {
  dates <- as.Date(dates)
  
  # get hydrological year
  m <- as.numeric(format(dates, format = "%m"))
  hy <- y <- as.numeric(format(dates, format = "%Y"))
  m_prevy <- !m%in%c(month:12)
  hy[m_prevy] <- hy[m_prevy] - 1
  
  # get days of the year
  if (day_of_the_year) {
    j <- as.numeric(format(dates, format = "%j"))
    start_hy <- as.Date(paste0(y, "-", month, "-1"))
    start_y <- as.Date(paste0(y, "-1-1"))
    j_hy <- as.numeric(format(dates - start_hy + start_y, format = "%j"))
  } else {
    j_hy <- NA
  }
  
  # get seasons
  if (!is.null(seasons) && is.list(seasons)) {
    # check season formatting
    if (!all(vapply(seasons, is.atomic, logical(1L)))) {
      stop("Wrong formatting for 'seasons': it should be a named list containing only atomic vectors (e.g. not lists).")
    }
    s_names <- names(seasons)
    s_names <- factor(s_names, levels = s_names)
    # build a table with character string and corresponding factors for matching with the whole time series
    s_table <- data.frame(i = unname(unlist(seasons)), season = rep(s_names, unlist(lapply(seasons, length))))
    if (seasons_by_day && day_of_the_year){
      s <- s_table[match(j, s_table[, 1]), 2]
    } else {
      if (seasons_by_day) stop("If 'seasons_by_day' is TRUE, you must also have 'day_of_the_year' set to TRUE.")
      
      s <- s_table[match(m, s_table[, 1]), 2]
    }
    if (minimal) return(data.frame(hy = hy, j_hy = j_hy, s = s))
    data.frame(dates = dates, m = m, y = y, hy = hy, j = j, j_hy = j_hy, s = s)
  } else {
    if (minimal) return(data.frame(hy = hy, j_hy = j_hy))
    data.frame(dates = dates, m = m, y = y, hy = hy, j = j, j_hy = j_hy)
  }
}

#------ hsaValidHydroYear ------------------------------------------------------
# Given a hydrological year vector and a matrix with the same number
# of rows (and any number of column), this function returns a logical
# vector of the same length as the hydrological year vector that
# indicates whether or not the time steps are part of a valid
# hydrological year. Validity of a hydrological year is assessed 
# according the two following rules:
#  - is the year complete (i.e. is there at least 'n' (365) days)?
#  - is there less than 'na.th' (proportion) missing values in 
#    the year?
# hy:           a hydrological year vector (will be coerced to factor)
# x:            a matrix (or a vector) with as many row (elements) as the length 
#               of 'hy' used to look for missing values
# n:            the minimum length of a complete hydrological year
#               this is also use to compute the proportion of missing 
#               values
# na.th:        the minimal tolerated proportion of missing values
hsaValidHydroYear <- function(hy, x = NULL, n = 365, na.th = 0.005) {
  # only full year
  rle_res <- rle(hy)
  unique_hy<- rle_res$values
  valid_hy_1 <- rep(rle_res$lengths >= n, rle_res$lengths)
  
  if (!is.null(x))  {
    # if x is provided, only year where the percentage of missing value is below na.th
    if (is.null(dim(x))) x <- matrix(x, length(x), 1)
    na <- apply(apply(x, 2, is.na), 1, any)
    valid_hy_2 <- rep(as.vector(tapply(na, hy, sum) / n <= na.th), rle_res$lengths)
    valid_hy_1 & valid_hy_2
  } else {
    valid_hy_1
  }
}


#------ get_nc_var_present -----------------------------------------------------
get_nc_var_present <- function(nc, varname, reachID, dates, selected_sims=1:20) {
  nc_data <- ncvar_get(nc, varname)
  
  dates_format <- hsaHydroYearSeasons(
    dates, month = 10, day_of_the_year = TRUE, # compute hydrological years, and other time periods (seasons, month, civil year, day of the year)
    seasons = list("SummerFall" = 5:10, "WinterSpring" = c(11:12, 1:4)),
    seasons_by_day = FALSE, minimal = FALSE) %>%
    setDT %>%
    .[, complete_year := hsaValidHydroYear(y, n = 365, na.th = 0.005)] %>%
    setnames(c('dates', 'j', 'm', 's', 'y'), 
             c('date', 'doy', 'month', 'season', 'year'))
  
  get_nc_var_inner <- function(in_nc_data) {
    as.data.table(in_nc_data) %>%
      setnames(as.character(reachID)) %>%
      .[, date := dates] %>%
      data.table::melt(id.vars = 'date', 
                       variable.name = 'reach_id',
                       value.name = varname,
                       variable.factor = FALSE) %>%
      .[, reach_id := as.integer(reach_id)]
  }
  
  #Check if there are multiple layers/sims in the netcdf
  nsims <- dim(nc_data)[3]
  if (!is.na(nsims)) {
    sims <- seq(1, nsims)
    out_dt <- lapply(sims[sims %in% selected_sims], function(in_sim) {
      #print(in_sim)
      get_nc_var_inner(nc_data[,,in_sim]) %>%
        .[, nsim := in_sim]
    }) %>%
      rbindlist
  } else {
    out_dt <- get_nc_var_inner(nc_data)
  }
  
  return(list(
    data_all = out_dt,
    dates_format = dates_format
  ))
}

#------ zero_lomf -----------------
#' Last \[non-zero\] Observation Moved Forward (lomf)
#'
#' Finds the index, for each row, of the previous row with a non-zero value
#'
#' @param x Numeric vector.
#' @param first (logical) Whether to consider first value as a non-zero value
#'   whose index is moved forward even if it is zero. This prevents having NAs
#'   in the results and somewhat assumes that, for a time series, the day prior
#'   to the first value is non-zero.
#'
#' @return Numeric vector of the indices of the previous non-zero for each
#'   element of the input vector.
#'
#' @examples
#' test1 <- c(1,1,1,0,0,0,0,1,1)
#' zero_lomf(test1)
#' test2 <- c(0,0,0,0,0,1,1,0,1)
#' zero_lomf(test2, first=FALSE)
#' zero_lomf(test2, first=TRUE)
#'
#' @export
#' 
zero_lomf <- function(x, first=TRUE) {
  if (length(x) > 0) {
    non.zero.idx <- which(x != 0)
    if (x[1]==0) {
      if(first==T) {
        non.zero.idx=c(1,non.zero.idx)
        #Repeat index of previous row with non-zero as many times gap until next non-zero values
        rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
      } else {
        non.zero.idx=c(0,non.zero.idx)
        #Repeat index of previous row with non-zero as many times gap until next non-zero values
        rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))[-1]
      } 
    } else {
      rep.int(non.zero.idx, diff(c(non.zero.idx, length(x) + 1)))
    }
  }
}
#------ create_sitepoints_raw --------------------------------------------------
create_sitepoints_raw <- function(in_dt, lon_col, lat_col, out_points_path,
                                  columns_to_include=NULL) {
  #Create point feature class from formatted site data
  sitesp <- terra::vect(in_dt,
                        geom = c(lon_col, lat_col),
                        crs = "+proj=longlat +datum=WGS84")
  if (is.null(columns_to_include)) {
    columns_to_include <- names(sitesp)
  }
  terra::writeVector(sitesp[, columns_to_include], 
                     out_points_path, overwrite=TRUE)
  return(out_points_path)
}
#------ compute_ecdf_lookup ----------------------------------------------------
#Batch-Compute left-continuous ecdf values (if nine 0s and one 1, then 0s are given 0 and 1 is given 0.9)
#and compare
#https://stats.stackexchange.com/questions/585291/is-there-an-equivalent-to-an-ecdf-with-a-sign
#https://math.stackexchange.com/questions/1807120/why-arent-cdfs-left-continuous/1807136#1807136
#Change the ecdf to be left-continuous.
compute_ecdf_lookup <- function(in_dt, ecdf_column, grouping_columns, 
                                na.rm=TRUE) {
  # Filter out NAs in the specified ECDF column if na.rm is TRUE
  #get frequency of each value of the ecdf column for each multicategory
  #of grouping columns (e.g., number days with a given discharge value for each
  #reach and day of year across the record). 
  dt_freq <- in_dt[if (na.rm) !is.na(get(ecdf_column)) else TRUE,
                   list(freq = .N), by = c(ecdf_column, grouping_columns)]
  
  # Sort data
  setorderv(dt_freq, ecdf_column)
  
  # Compute relative cumulative frequencies for the ECDF
  dt_freq[, paste0('P', ecdf_column) := (cumsum(freq) - freq) / sum(freq),
          by = grouping_columns]
  
  return(dt_freq[, c(paste0('P', ecdf_column), ecdf_column, grouping_columns), 
                 with=F])
}

#------ compute_ecdf_multimerge ------------------------------------------------
compute_ecdf_multimerge <- function(in_dt, ecdf_columns, grouping_columns, 
                                    keep_column, na.rm=TRUE) {
  ecdf_lookup_list <- lapply(ecdf_columns, function(in_ecdf_column) {
    compute_ecdf_lookup(in_dt = in_dt,
                        ecdf_column = in_ecdf_column, 
                        grouping_columns = grouping_columns, 
                        na.rm=na.rm) %>%
      merge(in_dt, ., by=c(in_ecdf_column, grouping_columns), all.x=T) %>%
      .[, c(paste0('P', in_ecdf_column), grouping_columns, keep_column), with=F]
  })
  
  ecdf_out <- Reduce(
    function(...) {
      merge(..., by = c(grouping_columns, keep_column), all = TRUE, sort = FALSE)
    }, c(list(in_dt), ecdf_lookup_list))
  
  return(ecdf_out)
}

#------ compute_hydrostats_intermittence ---------------------------------------
compute_hydrostats_intermittence <- function(in_hydromod_dt,
                                             in_sites_dt) {
  # -- Compute network-wide statistics ---------------------------------------
  #minimum 7-day and 30-day average over previous year
  #average in previous 10, 30, 45, 60, 90, 120, 180, 365, 365*5, 365*10 
  
  total_reach_length <- unique(in_hydromod_dt, by='reach_id')[, sum(reach_length)]
  
  #RelFlow: Proportion of network length with flowing conditions (opposite of RelInt)
  relF_dt <- in_hydromod_dt[isflowing == 1,
                            list(relF = sum(reach_length)/total_reach_length)
                            , by=.(date, nsim)] %>%
    setorder(date)
  
  #ggplot(relF_dt, aes(x=date, y=relF, color=nsim)) + geom_line()
  
  #Compute min 7-day and 30-day relF
  relF_dt[, `:=`(
    relF7mean = frollmean(x=relF, n=7, align='center', na.rm=T),
    relF30mean = frollmean(x=relF, n=30, align='center', na.rm=T)
  ), by = nsim] %>%
    .[, `:=`(
      relF7mean_yrmin = frollapply(relF7mean, n=365, FUN=min, align='right'),
      relF30mean_yrmin = frollapply(relF30mean, n=365, FUN=min, align='right')
    ), by = nsim]
  
  #Compute previous mean over many windows
  meanstep <- c(10, 30, 60, 90, 120, 180, 365, 365*5, 365*10)
  relF_dt[, paste0("relF", meanstep, "past") := frollmean(relF, meanstep, na.rm=T), 
          by = nsim]
  
  #PatchC: Patchiness of steady and intermittent flow conditions?
  #proportion of model-derived reach length with changing flow conditions 
  #compared to downstream reaches
  #-> not sure how to compute it. Hard to access JAM source files
  
  #Size of flowing connected patch
  
  # -- Compute statistics for specific reaches -------------------------------
  # DurD: DryDuration
  # PDurD: DryDuration_relative_to_longterm
  # FreD: Drying frequency - absolute or relative number of drying events per time interval
  # PFreD: FreD_relative_to_longterm
  
  hydromod_dt_sites <- in_hydromod_dt[reach_id %in% unique(in_sites_dt$reach_id),] %>%
    setorderv(c('nsim', 'reach_id', 'date'))
  
  #Compute duration of no-flow periods and time since last no-flow period #PrdD: prior days to last dry/pool/flowing event
  hydromod_dt_sites[, `:=`(noflow_period = rleid(isflowing==0), #Compute no-flow periods
                           last_noflow = zero_lomf(isflowing==0, first=FALSE) #Compute row index of last no flow day
  ), by=.(nsim, reach_id)] %>%
    .[isflowing == 1, noflow_period := NA] %>% 
    .[last_noflow ==0, last_noflow := NA] %>%
    .[!is.na(noflow_period), noflow_period_dur := .N,  #Compute duration of each no-flow period
      by = .(noflow_period, nsim, reach_id)] %>%
    .[, last_noflowdate := .SD[last_noflow, date], by = .(nsim, reach_id)] %>% #Convert to data
    .[, PrdD := difftime(date, last_noflowdate, units='days'), #Compute duration of each no-flow period
      by = .(nsim, reach_id)] %>%
    .[, last_noflow := NULL]
  
  # ggplot(hydromod_dt_sites[nsim==1,], aes(x=date, y=time_to_lastzero,color=nsim)) +
  #   geom_line() +
  #   facet_wrap(~id)
  
  #Compute monthly statistics --------------------------------------------------
  # qstats_absolute <- hydromod_dt_sites[,
  #                             list(
  #                               DurD = sum(isflowing==0),
  #                               FreD = length(na.omit(unique(noflow_period, na.rm=T))) #faster than uniqueN
  #                             ),
  #                             by = .(month, hy, reach_id, nsim)
  # ]
  # 
  # monthly_qstats_relative <- compute_ecdf_multimerge(
  #   in_dt = qstats_absolute,
  #   ecdf_columns = c('DurD', 'FreD'),
  #   grouping_columns = c('nsim', 'reach_id', 'month'), 
  #   keep_column = 'hy',
  #   na.rm=TRUE)
  
  #Compute moving-window statistics --------------------------------------------
  rollingstep_short <- c(10, 30, 60, 90, 120, 180)
  rollingstep_long <- c(365, 365*5, 365*10)
  rollingstep <- c(rollingstep_short, rollingstep_long)
  hydromod_dt_sites[, paste0("DurD", rollingstep, "past") :=  
                      frollapply(isflowing, n=rollingstep, 
                                 FUN=function(x) sum(x==0), 
                                 align='right'), by = .(reach_id, nsim)
  ]
  hydromod_dt_sites[, paste0("FreD", rollingstep, "past") := 
                      frollapply(noflow_period, n=rollingstep, 
                                 FUN = function(x) { #45% faster than using na.omit; 75% faster than uniqueN
                                   ux <- unique(x[!is.na(x)])
                                   length(ux)
                                 },
                                 align='right')  
                    , by = .(reach_id, nsim)
  ]
  
  #Compute ecdf value by day of year for moving windows
  rolling_column_names_short <- expand.grid(c('DurD', 'FreD'), 
                                            rollingstep_short) %>%
    setDT %>% 
    .[, paste0(Var1, Var2, 'past')]
  
  hydromod_dt_sites <- compute_ecdf_multimerge(
    in_dt = hydromod_dt_sites,
    ecdf_columns = rolling_column_names_short,
    grouping_columns = c('nsim', 'reach_id', 'doy'), 
    keep_column = 'date',
    na.rm=TRUE)
  
  #Compute ecdf value across entire records for moving windows > 365 days
  rolling_column_names_long <- expand.grid(c('DurD', 'FreD'), 
                                           rollingstep_long) %>% 
    setDT %>% 
    .[, paste0(Var1, Var2, 'past')]
  
  hydromod_dt_sites <- compute_ecdf_multimerge(
    in_dt = hydromod_dt_sites,
    ecdf_columns = rolling_column_names_long,
    grouping_columns = c('nsim', 'reach_id'), 
    keep_column = 'date',
    na.rm=TRUE)
  
  hydromod_dt_sites <- merge(in_sites_dt[, .(id, reach_id)], 
                             hydromod_dt_sites, 
                             by = 'reach_id', all.x = T, all.y = F,
                             allow.cartesian = T) %>%
    setorderv(c('nsim', 'reach_id', 'date'))
  
  return(list(
    network = relF_dt,
    sites = hydromod_dt_sites
  ))
}

#------ compute_hydrostats_q ---------------------------------------------------
compute_hydrostats_q <- function(in_hydromod_dt = hydromod_dt_sites) {
  
  setorderv(in_hydromod_dt, c('reach_id','date'))
  
  #Compute mean Q in the past X days for each day
  meanstep <- c(10, 30, 60, 90, 120, 180, 365, 365*5, 365*10)
  in_hydromod_dt[, paste0("meanQ", meanstep, "past") := 
                   frollmean(qsim, meanstep, na.rm=T), 
                 by = reach_id]
  
  #Compute daily flow percentile relative to the entire FDC
  in_hydromod_dt <- merge(in_hydromod_dt,
                          compute_ecdf_lookup(in_dt = in_hydromod_dt, 
                                              ecdf_column = 'qsim', 
                                              grouping_columns = 'reach_id', 
                                              na.rm=TRUE),
                          by=c('reach_id', 'qsim')
  )
  # ggplot(compute_ecdf_lookup(in_dt = in_hydromod_dt,
  #                            ecdf_column = 'qsim',
  #                            grouping_columns = 'reach_id',
  #                            na.rm=TRUE), 
  #        aes(x=qsim, y=Pqsim, color=as.factor(reach_id), group=reach_id)) +
  #   geom_line() +
  #   scale_x_log10()
  
  setorderv(in_hydromod_dt, c('reach_id','date'))
  
  #Proportion of days below overall 10th percentile in the past c(10, 30, 60, 90, 120, 180, 365) days
  rollingstep <- c(10, 30, 60, 90, 120, 180, 365)
  in_hydromod_dt[, paste0("uQ90_", rollingstep, "past") := 
                   frollapply(Pqsim, n = rollingstep, 
                              FUN = function(x) { #45% faster than using na.omit; 75% faster than uniqueN
                                length(x[x<=0.1])/length(x)
                              },
                              align='right'), 
                 by = reach_id]
  
  #Proportion of days above overall 90th percentile in the past c(10, 30, 60, 90, 120, 180) days
  in_hydromod_dt[, paste0("uQ10_", rollingstep, "past") := 
                   frollapply(Pqsim, n = rollingstep, 
                              FUN = function(x) { #45% faster than using na.omit; 75% faster than uniqueN
                                length(x[x>0.9])/length(x)
                              },
                              align='right'), 
                 by = reach_id]
  
  # ggplot(in_hydromod_dt[reach_id==1033201,], aes(x=date)) +
  #   geom_line(aes(y=Pqsim), color='blue') +
  #   geom_line(aes(y=uQ10_30past, group=reach_id), color='orange') +
  #   geom_line(aes(y=uQ10_365past, group=reach_id), color='red') 
  
  #Maximum flow percentile (based on long-term record) in the previous 3, 10, 30, 60, 90, 120, 180 days
  in_hydromod_dt[, paste0("maxPQ_", c(3, rollingstep), "past") := 
                   frollapply(Pqsim, n = c(3, rollingstep), FUN=max, align='right'), 
                 by = reach_id]
  
  #Compute percentiles of mean Q within past time windows
  #For 10 (to 180) days, for example, compute an FDC of the mean Q in the previous 
  #10 (to 180) days for each day of the year. Then compute for each day the probability
  #on the FDC (i.e., comparing to other years for this DOY)
  in_hydromod_dt <- compute_ecdf_multimerge(
    in_dt = in_hydromod_dt,
    ecdf_columns = paste0("meanQ", c(10, 30, 60, 90, 120, 180), "past"),
    grouping_columns = c('reach_id', 'doy'), 
    keep_column = 'date',
    na.rm=TRUE)
  
  #For 365 (to 3650) days, compute an FDC for the entire record of the mean Q
  #in the previous 365 (to 3650) days
  in_hydromod_dt <- compute_ecdf_multimerge(
    in_dt = in_hydromod_dt,
    ecdf_columns = paste0("meanQ", c(365, 365*5, 365*10), "past"),
    grouping_columns = c('reach_id'), 
    keep_column = 'date',
    na.rm=TRUE)
  
  setorderv(in_hydromod_dt, c('reach_id','date'))
  
  return(in_hydromod_dt)
}

#-------------- workflow functions ---------------------------------------------
# path_list = tar_read(bio_data_paths)
# in_metadata_edna <- tar_read(metadata_edna)

#------ define_hydromod_paths --------------------------------------------------
#in_hydromod_dir <- hydromod_present_dir

#List data paths for hydrological data
define_hydromod_paths <- function(in_hydromod_dir) {
  hydro_drn_paths_dt <- data.table(
    country = c("Croatia", "Czech", "Finland", "France",  "Hungary", "Spain"),
    catchment = c("Butiznica", "Velicka", "Lepsamaanjoki", "Albarine", "Bukkosdi", "Genal"), 
    all_sims_filename = c(
      "Butiznica_2022-12-15_option0.nc", #_run8_final
      "Velicka_2023-02-01_option0.nc",
      "Lepsamaanjoki_2022-12-16_option0.nc",
      'Albarine_2022-12-16_option0.nc',
      "Bukkosdi_2022-12-16_option0.nc", #_run8_final
      "Genal_2023-01-18_option0.nc"
    ),
    sel_sim_filename = c(
      "Butiznica_2022-12-15_option0_run8_final.nc",
      "Velicka_2023-02-01_option0_run9_final.nc",
      "Lepsamaanjoki_2022-12-16_option0_run20_final.nc",
      "Albarine_2022-12-16_option0_run3_final.nc",
      "Bukkosdi_2022-12-16_option0_run8_final.nc", 
      "Genal_2023-01-18_option0_run15_final.nc"
    )
  ) %>%
    .[, `:=`(all_sims_path = file.path(in_hydromod_dir, catchment, 
                                       "Results_present_period", all_sims_filename),
             sel_sim_path = file.path(in_hydromod_dir, catchment, 
                                      "Results_present_period", sel_sim_filename),
             catchment_path = file.path(in_hydromod_dir, catchment,
                                        "watershed_small_catchment.shp"),
             network_path = file.path(in_hydromod_dir, catchment,
                                      "river_network.shp"),
             sites_reachids = file.path(in_hydromod_dir, catchment,
                                        paste0(catchment, 
                                               '_sampling_sites_ReachIDs.csv')))]
}

#------ get_drn_hydromod -------------------------------------------------------
get_drn_hydromod <- function(hydromod_path, varname, selected_sims) {
  nc <- nc_open(hydromod_path) # open netcdf file
  reachID <- ncvar_get(nc, "reachID") # get list of reaches IDs
  dates <- ncvar_get(nc, "date") # get dates of simulation period
  dates <- as.Date(dates, origin="1950-01-01") # convert dates into R date format
  
  out_dt <- get_nc_var_present(nc = nc, varname = varname, # 0=dry, 1=flowing
                               reachID = reachID, dates = dates,
                               selected_sims = selected_sims) 
  return(out_dt)
}

#------ read_envdt -------------------------------------------------------------
# in_env_data_path_annika <- tar_read(env_data_path_annika)
# in_env_data_path_common <- tar_read(env_data_path_common)

read_envdt <- function(in_env_data_path_annika, 
                       in_env_data_path_common) {
  env_dt_annika <- fread(in_env_data_path_annika)
  env_dt_common <- fread(in_env_data_path_common)
  
  #Analyze missing data:
  #Finland, campaign #1, LEP21, not sampled
  #France, campaign #1, AL06, not sampled for macroinvertebrates
  #       campaign #3, RA01, impossible to access due to floods
  #       campaign #4, CA02, substrate, depth and velocity miv not recorded 
  #       campaign #5, lots of sites could not be sampled due to dangerously high flows
  #                     AL03, substrate, depth and velocity miv not recorded 
  #       campaign #6, data sheet lost for AL02
  #                     JO01 and RO01 were frozen
  #                     AL04, substrate, depth and velocity miv not recorded 
  #Hungary, campaign #4, BUK10 too little water to measure depth, velocity
  #Spain, campaign #1, GEN-04 was not sampled
  #       campaign #6, GEN-21 was not sampled
  
  #No oxygen saturation measure for Finland
  #No basin area for Spain
  #No upstream geology for Hungary
  #No embeddedness, substrate type, oxygen sat, etc. for Dry sites
  
  #Merge Annika's env data and that from the final data directory from teams
  #Several typoes/issues have been corrected in the final data but the rounding
  #has been badly performed, so keep Annika's for most columns, aside for those
  #where there have been significant upates to the data
  compare_env_dts <- merge(
    env_dt_annika, 
    env_dt_common[, names(env_dt_annika), with=F],
    by='running_id', suffixes=c('_annika', '_common')) %>%
    .[, sort(names(.)), with=F] %>%
    data.table::melt(id.var='running_id') %>%
    .[, value := trimws(value, which='both')] %>%
    .[, orig_dt := grepl('.*_annika$', variable)] %>%
    .[, variable := gsub('(_annika)|(_common)', '', variable)] %>%
    data.table::dcast(formula=running_id+variable~orig_dt) %>%
    setnames(c('FALSE', 'TRUE'), c('common', 'annika'))
  
  check_diff <- compare_env_dts[
    common!=annika & 
      (is.na(as.numeric(common)) |
         (!is.na(as.numeric(common)) & 
            (abs(as.numeric(common)-as.numeric(annika))/as.numeric(annika))>0.1)),] 
  
  #Use updated data from the following columns because errors have been corrected
  cols_common <- c('Average_wetted_width_m',
                   'Bankfull_at_max__wetted_width_m',
                   'Bankfull_at_min__wetted_width_m',
                   'conductivity_micros_cm',
                   'Max_wetted_width_m',
                   'Min_wetted_width_m',
                   'native_riparian_species_richness',
                   'oxygen_mg_l',
                   'ph',
                   'state_of_flow',
                   'temperature_C')
  
  env_dt_merged <- merge(env_dt_annika[, -cols_common, with=F],
                         env_dt_common[, c('running_id', cols_common), with = F]) %>%
    setnames(tolower(names(.))) %>% #Convert all columns to lower case
    .[!(is.na(date) & is.na(state_of_flow)),] %>% #Remove records that were not sampled at all, but keep those that were simply dry
    setnames(c('bankfull_at_max__wetted_width_m', 'bankfull_at_min__wetted_width_m'),
             c('bankfull_at_max_wetted_width_m', 'bankfull_at_min_wetted_width_m'))
  
  #Substitute date for those that are dry
  #env_dt_merged[is.na(date),]
  env_dt_merged[running_id == 'AL02_6', date := as.Date('20/01/2022')]
  env_dt_merged[running_id %in% c('GEN04_6', 'GEN10_6', 'GEN11_6', 'GEN13_6'),
                date := as.Date('10/02/2022')]
  
  #Trim white spaces in character columns
  env_dt_merged[, names(.SD) := lapply(.SD, str_trim), .SDcols=is.character]
  
  #Check for empty strings
  col_w_empty_strings <- env_dt_merged[
    ,sapply(.SD, function(x) any(x=='', na.rm=T))] %>%
    names(.)[.]
  #if_ip_number_and_size_2_axes_+_depth_of_the_pools is 
  #the only column with empty strings
  
  #Correct issues
  env_dt_merged[running_id %in% c('AL02_4', 'AL03_4'), #Correct typo
                state_of_flow := 'D']
  env_dt_merged[running_id %in% c('JO01_6', 'RO01_6'), #Replace NA with frozen
                state_of_flow := 'FROZEN']
  #Remove un-analyzable records
  env_dt_merged <- env_dt_merged[state_of_flow != 'FROZEN',]
  
  env_dt_merged[running_id == 'BUK10_4', avg_depth_macroinvertebrates := 1] #too small to measure
  env_dt_merged[state_of_flow == 'IP' & is.na(avg_velocity_macroinvertebrates), #assign 0 velocity to pools
                avg_velocity_macroinvertebrates := 0]
  
  #In a dozen site in Hungary, "could not measure velocity" -- assign very low value just below existing minimum
  env_dt_merged[
    drn=='Hungary' & is.na(avg_velocity_macroinvertebrates) & state_of_flow == 'F',
    avg_velocity_macroinvertebrates := 0.001]
  
  #Replace missing values for geomorphological measures with average across campaigns
  invisible(lapply(c('embeddedness', 'maximum_depth_cm', "filamentous_algae", 
                     "incrusted_algae", "macrophyte_cover", "leaf_litter_cover",
                     "moss_cover","wood_cover"), 
                   function(col) {
                     env_dt_merged[state_of_flow %in% c('F'), 
                                   (col) := nafill(get(col), type = "const", 
                                                   fill = mean(get(col), na.rm = TRUE)),
                                   by=site]
                   }))
  
  #Other geomorphological attributes are NA for dry periods
  env_dt_merged[state_of_flow == 'D', `:=`(
    maximum_depth_cm = 0,
    discharge_l_s = 0,
    average_wetted_width_m = 0,
    max_wetted_width_m = 0,
    min_wetted_width_m = 0,
    avg_depth_macroinvertebrates = 0,
    avg_velocity_macroinvertebrates = 0
  )]
  
  env_dt_merged[state_of_flow == 'D', `:=`(
    conductivity_micros_cm = NA,
    oxygen_mg_l = NA,
    oxygen_sat = NA,
    ph = NA,
    temperature_c = NA
  )]
  
  #Fill missing data for average depth and velocity for macroinvertebrates based on discharge
  #using the average ratio between log(discharge) and these measures in other dates
  avg_ratio_sampling_hydraulic_dis <- env_dt_merged[
    drn=='France' & avg_depth_macroinvertebrates>0 & discharge_l_s>1, 
    list(
      site_id = site,
      ratio_depth = avg_depth_macroinvertebrates/log10(discharge_l_s),
      ratio_velo = avg_velocity_macroinvertebrates/log10(discharge_l_s)
    ), by=site
  ] %>%
    .[, list(mean_ratio_d = mean(ratio_depth),
             mean_ratio_v = mean(ratio_velo)),
      by=site_id]
  
  env_dt_merged[
    state_of_flow == 'F' & is.na(avg_depth_macroinvertebrates) & 
      discharge_l_s>1,
    avg_depth_macroinvertebrates := log10(discharge_l_s)*
      avg_ratio_sampling_hydraulic_dis[avg_ratio_sampling_hydraulic_dis$site_id == site, 'mean_ratio_d'],
    by=site]
  
  env_dt_merged[
    state_of_flow == 'F' & is.na(avg_velocity_macroinvertebrates) & 
      discharge_l_s>1,
    avg_velocity_macroinvertebrates := log10(discharge_l_s)*
      avg_ratio_sampling_hydraulic_dis[avg_ratio_sampling_hydraulic_dis$site_id == site, 'mean_ratio_v'],
    by=site]
  
  #Compute simple discharge for the site
  env_dt_merged[running_id == 'AL07_1',
                discharge_l_s := avg_velocity_macroinvertebrates*
                  (average_wetted_width_m*avg_depth_macroinvertebrates/100)]
  
  #Complete altitude for GEN04 based on topographic map from Spanish IGN
  env_dt_merged[site == 'GEN04', altitude_m := 610]
  
  #Fill average wetted width with average ratio between max and min for that site
  mean_width_ratio <- env_dt_merged[site == 'GEN09', mean(
    (average_wetted_width_m - min_wetted_width_m)/
      (max_wetted_width_m - min_wetted_width_m),
    na.rm=T)]
  env_dt_merged[running_id == 'GEN09_3',
                average_wetted_width_m := mean_width_ratio*
                  (max_wetted_width_m - min_wetted_width_m)]
  
  #Fill conductivity data for the few Hungarian sites for campaign #5
  #by computing their conductivity compared to the mean across sites
  avg_conductivity_ratio <- env_dt_merged[
    drn=='Hungary' & campaign != 5  & state_of_flow == 'F',
    list(
      site_id = site,
      ratio=conductivity_micros_cm/mean(conductivity_micros_cm, na.rm=T),
      mean_conduct = mean(conductivity_micros_cm, na.rm=T)
    ), by=campaign] %>%
    .[, mean(ratio), by=site_id]
  
  avg_conduct_5 <- env_dt_merged[
    drn=='Hungary' & campaign == 5  & state_of_flow == 'F', 
    mean(conductivity_micros_cm, na.rm=T)]
  
  env_dt_merged[running_id %in% c('BUK09_5', 'BUK10_5','BUK50_5','BUK52_5'),
                conductivity_micros_cm := avg_conduct_5*
                  avg_conductivity_ratio[avg_conductivity_ratio$site_id == site, 'V1'],
                by=site]
  
  #Correct riparian area based on observation of satellite imagery
  env_dt_merged[running_id == 'BUK36_1', riparian_cover_in_the_riparian_area := 90]
  
  #Use corrected data for some records
  env_dt_merged[running_id == 'BUK42_3', bankfull_at_max_wetted_width_m := 4.7] #from original data sheet uploaded on Teams
  
  #Check that max is larger than min
  check <- env_dt_merged[max_wetted_width_m < min_wetted_width_m,]
  env_dt_merged[running_id=='BUK50_4', max_wetted_width_m := min_wetted_width_m]
  env_dt_merged[running_id=='BUK50_4', min_wetted_width_m := 2]
  substitute_width <- env_dt_merged[average_wetted_width_m < min_wetted_width_m, min_wetted_width_m]
  env_dt_merged[average_wetted_width_m < min_wetted_width_m,
                `:=`(min_wetted_width_m = average_wetted_width_m,
                     average_wetted_width_m =  substitute_width)]
  env_dt_merged[running_id=='BUK30_2', 
                max_wetted_width_m := average_wetted_width_m]
  env_dt_merged[running_id=='BUK30_2', 
                average_wetted_width_m := 3.05]
  
  #Inspect data
  #skim(env_dt_merged)
  #skim(env_dt_merged[state_of_flow == 'F',])
  
  return(env_dt_merged)
}
#------ read_biodt -------------------------------------------------------------
# path_list = tar_read(bio_data_paths)
# in_metadata_edna = tar_read(metadata_edna)

read_biodt <- function(path_list, in_metadata_edna) {
  #Read and name all data tables
  dt_list <- mapply(function(in_path, in_name) {
    fread(in_path) %>% 
      .[, organism := in_name] %>%
      setnames(tolower(names(.)))},
    path_list, names(path_list))
  
  #Add Campaign and Site to bacteria data, then remove pool sites
  dt_list$bac_sedi[, c('site', 'campaign') := tstrsplit(v1, '_')] %>%
    setnames('v1', 'running_id') %>%
    .[, campaign := as.integer(campaign)]
  dt_list$bac_biof[, c('site', 'campaign') := tstrsplit(v1, '_')] %>%
    setnames('v1', 'running_id') %>%
    .[, campaign := as.integer(campaign)]
  
  #Standardize country names
  country_standard <- data.table(
    original = c("CRO", "FRA", "SPA", "CZ",  "HUN", "FIN"),
    new = c("Croatia", 'France', "spain", "Czech", "Hungary", 'Finland')
  )
  in_metadata_edna <- merge(in_metadata_edna, country_standard,
                            by.x='country', by.y='original') %>%
    .[, country := new] %>%
    .[, `:=`(new = NULL)]
  
  #Remove bacteria in pools
  dt_list$bac_sedi <- merge(dt_list$bac_sedi, 
                            in_metadata_edna[sample_type=='sediment', 
                                             .(matchingenvid, habitat, country)],
                            by.x='running_id', by.y='matchingenvid')
  dt_list$bac_biof <- merge(dt_list$bac_biof,
                            in_metadata_edna[sample_type=='biofilm',
                                             .(matchingenvid, habitat, country)],
                            by.x='running_id', by.y='matchingenvid')
  
  dt_list$bac_biof_nopools <- dt_list$bac_biof %>%
    .[habitat!='pool',] %>%
    .[, habitat := NULL]
  dt_list$bac_biof[, habitat := NULL]
  
  
  dt_list$bac_sedi_nopools <- dt_list$bac_sedi %>%
    .[habitat!='pool',] %>%
    .[, habitat := NULL]
  dt_list$bac_sedi[, habitat := NULL]
  
  return(dt_list)
}

#------ calc_sprich ------------------------------------------------------------
# in_biodt <- tar_read(bio_dt)[[1]]
# in_metacols <- metacols

calc_sprich <- function(in_biodt, in_metacols) {
  #Get metadata columns (all except species data)
  metacols_sub <- names(in_biodt)[names(in_biodt) %in% in_metacols]
  biodt_melt <- melt(in_biodt, id.vars = metacols_sub)
  #Compute species richness
  biodt_sprich <- biodt_melt[, list(richness = sum(value > 0)), 
                             by=.(site, campaign, organism)] %>%
    .[, mean_richness := mean(richness), by=site]
  return(biodt_sprich)
}

#------ sprich_plot ------------------------------------------------------------
# in_sprich <- tar_read(sprich)
# in_envdt <- tar_read(env_dt)

plot_sprich <- function(in_sprich, in_envdt) {
  
  sprich_hydroobs <- merge(
    in_sprich[, list(mean_drn_richness = mean(richness)),
              by=.(campaign, drn, organism)],
    in_envdt[, list(per_flowing = 100*.SD[state_of_flow=='F', .N]/.N),
             by=.(campaign, drn)],
    by=c('campaign', 'drn')) %>%
    .[, mean_drn_richness_relative := mean_drn_richness/max(mean_drn_richness),
      by=.(drn, organism)]
  
  ggplot(sprich_hydroobs[
    !(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=campaign, y=mean_drn_richness_relative)) +
    geom_line(aes(group=drn, color=drn), size=1.2) +
    new_scale_color() +
    geom_point(aes(color=per_flowing), size=2) +
    scale_color_distiller(palette='Spectral', direction=1) +
    facet_wrap(~organism, scales='free_y') +
    theme_classic()
  
  ggplot(sprich_hydroobs[
    !(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=per_flowing, y=mean_drn_richness_relative, color=drn)) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    facet_wrap(~organism, scales='free_y')
  
  overall_lines <- ggplot(
    in_sprich[!(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=campaign, y=richness)) +
    #geom_point() +
    geom_line(aes(group=site, color=site)) +
    geom_smooth(aes(group=drn)) +
    facet_wrap(drn~organism, scales = 'free_y') +
    theme(legend.position='none')
  
  ggplot(
    in_sprich[!(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=campaign, y=richness)) +
    #geom_point() +
    geom_boxplot(aes(group=site, color=state_of_flow)) +
    facet_wrap(drn~organism, scales = 'free_y') 
  
  ggplot(
    in_sprich[(organism %in% c('miv_nopools')),],
    aes(x=campaign, y=richness)) +
    #geom_point() +
    geom_line(aes(group=site, color=site)) +
    geom_smooth(aes(group=drn)) +
    facet_wrap(drn~organism, scales = 'free_y') +
    theme(legend.position='none')
}
#------ subset_network -----------------------------------------------------------
# in_hydromod_paths_dt <- tar_read(hydromod_paths_dt)
# out_dir <- file.path('results', 'gis')
# overwrite = FALSE

subset_network <- function(in_hydromod_paths_dt, out_dir, overwrite=FALSE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  in_hydromod_paths_dt[
    , network_sub_path := file.path(out_dir, 
                                    paste0(tolower(country), 
                                           '_river_network_sub.shp'))]
  
  in_hydromod_paths_dt[, {
    if (!file.exists(network_sub_path) | overwrite) {
      terra::vect(network_path) %>%
        .[relate(terra::vect(catchment_path), .,"contains")[1,],] %>%
        terra::writeVector(filename = network_sub_path)
    }
  }, by=country]
  
  return(in_hydromod_paths_dt[, stats::setNames(network_sub_path, country)])
}

#------ format_sites_dt ----------------------------------------------------------
# in_country <- 'Finland'
# in_path <- tar_read(hydromod_paths_dt)[country == in_country, sites_reachids]
# check <- format_site_dt(in_path, in_country)
# print(check)

format_site_dt <- function(in_path, in_country) {
  sites_dt <- fread(in_path) %>%
    setnames(tolower(names(.))) %>%
    .[, reach_id := as.integer(reach_id)]

  if (in_country == 'Croatia') {
    sites_dt[, id := sub('BUT', '', id) %>%
               str_pad(width=2, side='left', pad = 0) %>%
               paste0('BUT', .)]
  } 
  
  if (in_country == 'Czech') {
    sites_dt[, id := sub('S', '', id) %>%
               str_pad(width=2, side='left', pad = 0) %>%
               paste0('VEL', .)]
  } 
  
  if (in_country == 'France') {
    #Update reach_ids based on documentation in WP1 
    #"data\wp1\Results_present_period_final\intermittence_indicators\results until sampling sites of each DRN (small catchment)\Albarine\readme.txt"
    sites_dt[id=='VI01', reach_id := 2438600]
    sites_dt[id=='AL01', reach_id := 2455400]
  }
  
  if (in_country == 'Finland') {
      sites_dt[, id := sub('FIN', '', id) %>%
                 str_pad(width=2, side='left', pad = 0) %>%
                 paste0('LEP', .)]
  }
  
  if (in_country == 'Hungary') {
    sites_dt <- merge(sites_dt,
                      data.table(
                        reach_id_old = c(659601, 674601, 676601, 665201, 666401, 652001),
                        reach_id_new = c(659600, 674600, 676600, 665200, 666400, 652000)
                      ), by.x='reach_id', by.y='reach_id_old', all.x=T) %>%
      .[!is.na(reach_id_new), reach_id := reach_id_new] %>%
      .[, reach_id_new := NULL]
  }
  
  if (in_country == 'Spain') {
    sites_dt[, id := sub('spa', 'GEN', id)]
  }
  
  return(sites_dt[, .(id, lat, lon, reach_id)])
}

#------ create_sites_shp --------------------------------------------------------
# in_hydromod_paths_dt = tar_read(hydromod_paths_dt)
# in_sites_dt = tar_read(sites_dt)
# out_dir = file.path('results', 'gis')
# geom = 'points'
# overwrite = TRUE

create_sites_shp <- function(in_hydromod_paths_dt,
                             in_sites_dt,
                             out_dir, 
                             geom,
                             overwrite = FALSE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  setnames(setDT(in_sites_dt), 'country', 'country_sub')
  
  in_hydromod_paths_dt[
    , sites_shp_path := file.path(out_dir, 
                                  paste0(tolower(country), 
                                         '_site_', geom, '.shp'))]
  
  if (geom == 'reaches') {
    #Create site reaches
    in_hydromod_paths_dt[, {
      if (!file.exists(sites_shp_path) | overwrite) {
        terra::vect(network_path) %>%
          aggregate(by='cat') %>%
          merge(in_sites_dt[country_sub==country,],
                by.x='cat', by.y='reach_id', all.x=F) %>%
          terra::writeVector(filename = sites_shp_path, 
                             overwrite = T)
      }
    }, by=country]
  }

  if (geom == 'points') {
    #Create site points
    in_hydromod_paths_dt[, {
      if (!file.exists(sites_shp_path) | overwrite) {
        create_sitepoints_raw(in_dt = in_sites_dt[country_sub==country,], 
                              lon_col = 'lon', lat_col = 'lat',
                              out_points_path = sites_shp_path) 
      }
    }, by=country]
  }

  return(in_hydromod_paths_dt[, stats::setNames(sites_shp_path, country)])
}


#------ compute_hydrostats_drn -------------------------------------------------
# in_drn <- 'Czech'
# varname <- 'qsim'#'isflowing'
# in_sites_dt <- tar_read(sites_dt)[country == in_drn,]
# in_network_path <- tar_read(network_sub_shp_list)[[in_drn]]
# in_hydromod_drn <- tar_read(hydromod_dt_Czech_qsim)

compute_hydrostats_drn <- function(in_network_path,
                                   in_sites_dt,
                                   varname,
                                   in_hydromod_drn) {
  setDT(in_hydromod_drn$data_all)
  setDT(in_hydromod_drn$dates_format) 
  
  #-------------------- Compute intermittence statistics -----------------------
  if (varname == 'isflowing') {
    #Import network shapefiles and get their length
    network_v <- terra::vect(in_network_path)
    network_v$reach_length <- terra::perim(network_v)
    reach_dt <- as.data.table(network_v[, c('cat', 'reach_length')]) %>%
      setnames('cat', 'reach_id') %>%
      .[, list(reach_length = sum(reach_length)), by=reach_id]
    remove(network_v)
    
    #Merge hydro data and reach length for computing network-wide statistics
    #qdat[, unique(reach_id)[!(unique(reach_id) %in% reach_dt$reach_id)]] 
    #Two IDs are not in the shapefile for Finland? maybe a hydrological unit not associated
    intermod_dt <- merge(in_hydromod_drn$data_all, reach_dt, 
                         by='reach_id', all.x=F, all.y=F) %>%
      merge(in_hydromod_drn$dates_format[, .(date, month, hy, doy)],
            by='date', all.x=T, all.y=F)
    
    q_stats <- compute_hydrostats_intermittence(
      in_hydromod_dt = intermod_dt,
      in_sites_dt = in_sites_dt)
  }
  
  #-------------------- Compute discharge statistics ---------------------------
  if (varname == 'qsim') {
    #Keep only hydrological data for sampled reaches
    hydromod_dt_sites <- in_hydromod_drn$data_all[
      reach_id %in% unique(in_sites_dt$reach_id),] %>%
      setorderv(c('reach_id', 'date')) %>%
      merge(in_hydromod_drn$dates_format[, .(date, month, hy, doy)],
            by='date', all.x=T, all.y=F)
    
    q_stats <- compute_hydrostats_q(in_hydromod_dt = hydromod_dt_sites) %>%
      merge(in_sites_dt[, .(id, reach_id)], .,
            by = 'reach_id', all.x = T, all.y = F,
            allow.cartesian = T) %>%
      setorderv(c('reach_id', 'date'))
  }
  
  return(q_stats)
}

#------ format_envinterm -------------------------------------------------------
# in_env_dt <- tar_read(env_dt)
# in_interm90_dt <- tar_read(interm90_dt)
# in_sprich = tar_read(sprich)

merge_alphadat <- function(in_env_dt, in_interm90_dt, in_sprich) {
  #Compute mean 90-day drying duration and event length
  interm90_mean <- in_interm90_dt[
    , list(TotDur90 = mean(TotDur, na.rm=T),
           TotLeng90 = mean(TotLeng, na.rm=T)),
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



#------ plot_alpha_cor --------------------------------------------------------
# tar_load(alphadat_merged)

plot_alpha_cor_inner <- function(in_alphadat_merged_organism, x_var, facet_wrap=F) {
  #Compute simple linear regression
  in_alphadat_merged_organism[
    , lm_pval_ltype := fifelse(
      coef(summary(lm(mean_S~get(x_var))))[2,4] < 0.05,
      'solid', 'dashed'
    ),
    by=Country] 
  
  if (facet_wrap) {
    alpha_plots <- ggplot(in_alphadat_merged_organism, aes(x=get(x_var), y=mean_S)) + 
      geom_point(size = 2) + 
      geom_smooth(aes(linetype=lm_pval_ltype), method='lm', linewidth = 0.5, se = F) +
      scale_linetype_identity() +
      labs(x=x_var) +
      facet_wrap(~Country) +
      theme_classic()
  } else {
    alpha_plots <- ggplot(in_alphadat_merged_organism, aes(get(x_var), mean_S)) + 
      geom_point(aes(colour=Country), size = 1, alpha=0.5) +
      geom_smooth(aes(linetype=lm_pval_ltype, colour=Country), method='lm', linewidth = 0.5, se = F) +
      geom_smooth(colour="black", method = "lm", linewidth = 1.1, se = F) + 
      scale_linetype_identity() +   
      scale_color_manual(values = c("Croatia" = "#ef476f",
                                    "Czech Republic" = "#f78c6b", 
                                    "Finland" = "#ffd166", 
                                    "France" = "#06d6a0", 
                                    "Hungary" = "#118ab2",
                                    "Spain" = "#073b4c")) +
      labs(x=x_var) +
      theme_classic()
    #+  theme(legend.position = "none")
  }
  
  return(alpha_plots)
}  

plot_alpha_cor <- function(in_alphadat_merged, out_dir, facet_wrap=F) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  organism_list <- unique(in_alphadat_merged$organism)
  
  #Plot relationship between each organism alpha div and totdur90 for each country
  plotlist_totdur90 <- lapply(organism_list, function(org_sel) {
    out_p <- plot_alpha_cor_inner(in_alphadat_merged[organism == org_sel,],
                                  x_var = 'TotDur90',
                                  facet_wrap = facet_wrap) +
      ggtitle(org_sel)
    
    out_suffix <- paste0(ifelse(facet_wrap==F, '_all', ''), 
                         '_mean_S_vs_totdur90_lm_sig.png') 
    ggsave(
      filename = file.path(out_dir, paste0(org_sel, out_suffix)),
      plot = out_p, 
      width=10, height=10)
    
    return(out_p)
  })
  names(plotlist_totdur90) <- organism_list
  
  #Plot relationship between each organism alpha div and discharge for each country
  plotlist_discharge <- lapply(organism_list, function(org_sel) {
    out_p <- plot_alpha_cor_inner(in_alphadat_merged[organism == org_sel,],
                                  x_var = 'discharge',
                                  facet_wrap = facet_wrap) + 
      scale_x_log10() +
      ggtitle(org_sel)
    
    out_suffix <- paste0(ifelse(facet_wrap==F, '_all', ''), 
                         '_mean_S_vs_discharge_lm_sig.png') 
    ggsave(
      filename = file.path(out_dir, 
                           paste0(org_sel, '_mean_S_vs_discharge_lm_sig.png')),
      plot = out_p,
      width=10, height=10
    )
  })
  names(plotlist_discharge) <- organism_list
  
  return(list(
    totdur90 = plotlist_totdur90,
    discharge = plotlist_discharge
  ))
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

#------ merge_env_mod ----------------------------------------
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

#------ plot_z_by_stream_type ------------------------------
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

#------ plot_z_by_env ----------------------------------------
# in_env_null_models_dt <- tar_read(env_null_models_dt)
# env_var <- 'TotDur90'
# outdir <- resdir 

plot_z_by_env <- function(in_env_null_models_dt, env_var, outdir) {
  #### scatterplot, all countries in the same plot, separately for each organism group ####
  organism_list <- unique(in_env_null_models_dt$organism)
  
  plotlist <- lapply(organism_list, function(in_organism) {
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
    
    return(p1)
  })
  names(plotlist) <- organism_list
  
  return(plotlist)
}

#------ compute_lmer_mods ------------------------------------------------------
compute_lmer_mods <- function(in_dt, in_yvar) {
  lmer_int <- in_dt[
    ,  list(
      #TotDur90 models
      TotDur90_full = list(
        lmer(as.formula(paste(in_yvar, 
                              "~ TotDur90 + (1|Country)")), data=.SD)),
      TotDur90_null = list(lmer(as.formula(paste(in_yvar, 
                                                 "~ (1|Country)")), data=.SD)),
      TotDur90_ML = list(anova(lmer(as.formula(paste(in_yvar, 
                                                     "~ TotDur90 + (1|Country)")), data=.SD),
                               lmer(as.formula(paste(in_yvar, 
                                                     "~ (1|Country)")), data=.SD))),
      #TotLeng90 models
      TotLeng90_full = list(lmer(as.formula(paste(in_yvar, 
                                                  "~ TotLeng90 + (1|Country)")), data=.SD)),
      TotLeng90_null = list(lmer(as.formula(paste(in_yvar, 
                                                  "~ (1|Country)")), data=.SD)),
      TotLeng90_ML = list(anova(lmer(as.formula(paste(in_yvar, 
                                                      "~ TotLeng90 + (1|Country)")), data=.SD),
                                lmer(as.formula(paste(in_yvar, 
                                                      "~ (1|Country)")), data=.SD)))
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


