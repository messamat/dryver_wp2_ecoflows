#-------------- utility functions ----------------------------------------------
#------ download_unzip ---------------------------------------
download_unzip <- function(url, out_dir, download_mode='wb', out_zip=NULL) {
  # if (!dir.exists(out_dir)) {
  #   dir.create(out_dir)
  # }
  
  if (is.null(out_zip)) {
    out_zip <- paste0(out_dir, '.zip')
  }
  
  download.file(url = url,
                mode = download_mode,
                destfile = out_zip
  )
  
  unzip(out_zip, exdir = out_dir)
  
  return(out_dir)
}

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
get_nc_var_present <- function(nc, varname, reachID, dates, selected_sims=NULL) {
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
#------ mergDTlist -----------------------------------------
#Function to merge a list of data.tables (dt), adding a suffix to all columns
#by dt based on the name of the input data.table
mergeDTlist <- function(dt_list, by = NULL, all = TRUE, sort = FALSE,
                        set_suffix=TRUE) {
  
  if (set_suffix) {
    dt_list <-  Map(function(dt_name, dt) {
      dt_copy <- copy(dt)
      cols_to_rename <- names(dt)[!(names(dt) %in% by)]
      setnames(dt_copy,
               old=cols_to_rename,
               new=paste(cols_to_rename, dt_name, sep='_'))
      return(dt_copy)
    },
    names(dt_list), dt_list
    )
  }
  
  Reduce(
    function(...) {
      merge(..., by = by, all = all, sort = sort)
    }, dt_list)
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
                                             in_sites_dt,
                                             scale = 'all') {
  # -- Compute network-wide statistics ---------------------------------------
  if (scale %in% c('all', 'drn')) {
    #minimum 7-day and 30-day average over previous year
    #average in previous 10, 30, 45, 60, 90, 120, 180, 365, 365*5, 365*10 
    
    total_reach_length <- unique(in_hydromod_dt, by='reach_id')[, sum(reach_length)]
    
    #RelFlow: Proportion of network length with flowing conditions (opposite of RelInt)
    relF_dt <- in_hydromod_dt[isflowing == 1,
                              list(relF = sum(reach_length)/total_reach_length)
                              , by=.(date)] %>%
      setorder(date)
    
    #ggplot(relF_dt, aes(x=date, y=relF)) + geom_line()
    
    #Compute min 7-day and 30-day relF
    relF_dt[, `:=`(
      relF7mean = frollmean(x=relF, n=7, align='center', na.rm=T),
      relF30mean = frollmean(x=relF, n=30, align='center', na.rm=T)
    )] %>%
      .[, `:=`(
        relF7mean_yrmin = frollapply(relF7mean, n=365, FUN=min, align='right'),
        relF30mean_yrmin = frollapply(relF30mean, n=365, FUN=min, align='right')
      )]
    
    #Compute previous mean over many windows
    meanstep <- c(10, 30, 60, 90, 120, 180, 365, 365*5, 365*10)
    relF_dt[, paste0("relF", meanstep, "past") := frollmean(relF, meanstep, na.rm=T)]
    
    #PatchC: Patchiness of steady and intermittent flow conditions?
    #proportion of model-derived reach length with changing flow conditions 
    #compared to downstream reaches
    #-> not sure how to compute it. Hard to access JAM source files
    
    #Size of flowing connected patch
  }
  
  # -- Compute statistics for specific reaches -------------------------------
  # DurD: DryDuration
  # PDurD: DryDuration_relative_to_longterm
  # FreD: Drying frequency - absolute or relative number of drying events per time interval
  # PFreD: FreD_relative_to_longterm
  if (scale %in% c('all', 'site')) {
    hydromod_dt_sites <- in_hydromod_dt[reach_id %in% unique(in_sites_dt$reach_id),] %>%
      setorderv(c('reach_id', 'date'))
    
    #Compute duration of no-flow periods and time since last no-flow period #PrdD: prior days to last dry/pool/flowing event
    hydromod_dt_sites[, `:=`(noflow_period = rleid(isflowing==0), #Compute no-flow periods
                             last_noflow = zero_lomf(isflowing==0, first=FALSE) #Compute row index of last no flow day
    ), by=.(reach_id)] %>%
      .[isflowing == 1, noflow_period := NA] %>% 
      .[last_noflow ==0, last_noflow := NA] %>%
      .[!is.na(noflow_period), noflow_period_dur := .N,  #Compute duration of each no-flow period
        by = .(noflow_period, reach_id)] %>%
      .[, last_noflowdate := .SD[last_noflow, date], by = .(reach_id)] %>% #Convert to data
      .[, PrdD := difftime(date, last_noflowdate, units='days'), #Compute duration of each no-flow period
        by = .(reach_id)] %>%
      .[, last_noflow := NULL]
    
    # ggplot(hydromod_dt_sites, aes(x=date, y=time_to_lastzero)) +
    #   geom_line() +
    #   facet_wrap(~id)
    
    #Compute monthly statistics --------------------------------------------------
    # qstats_absolute <- hydromod_dt_sites[,
    #                             list(
    #                               DurD = sum(isflowing==0),
    #                               FreD = length(na.omit(unique(noflow_period, na.rm=T))) #faster than uniqueN
    #                             ),
    #                             by = .(month, hy, reach_id)
    # ]
    # 
    # monthly_qstats_relative <- compute_ecdf_multimerge(
    #   in_dt = qstats_absolute,
    #   ecdf_columns = c('DurD', 'FreD'),
    #   grouping_columns = c('reach_id', 'month'), 
    #   keep_column = 'hy',
    #   na.rm=TRUE)
    
    #Compute moving-window DurD and FreD statistics --------------------------------------------
    rollingstep_short <- c(10, 30, 60, 90, 120, 180)
    rollingstep_long <- c(365, 365*5, 365*10)
    rollingstep <- c(rollingstep_short, rollingstep_long)
    hydromod_dt_sites[, paste0("DurD", rollingstep, "past") :=  
                        frollapply(isflowing, n=rollingstep, 
                                   FUN=function(x) sum(x==0), 
                                   align='right'), by = .(reach_id)
    ]
    hydromod_dt_sites[, paste0("FreD", rollingstep, "past") := 
                        frollapply(noflow_period, n=rollingstep, 
                                   FUN = function(x) { #45% faster than using na.omit; 75% faster than uniqueN
                                     ux <- unique(x[!is.na(x)])
                                     length(ux)
                                   },
                                   align='right')  
                      , by = .(reach_id)
    ]
    
    #Compute ecdf value by day of year for moving windows
    rolling_column_names_short <- expand.grid(c('DurD', 'FreD'), 
                                              rollingstep_short) %>%
      setDT %>% 
      .[, paste0(Var1, Var2, 'past')]
    
    hydromod_dt_sites <- compute_ecdf_multimerge(
      in_dt = hydromod_dt_sites,
      ecdf_columns = rolling_column_names_short,
      grouping_columns = c('reach_id', 'doy'), 
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
      grouping_columns = c('reach_id'), 
      keep_column = 'date',
      na.rm=TRUE)
    
    hydromod_dt_sites <- merge(in_sites_dt[, .(site, reach_id)], 
                               hydromod_dt_sites, 
                               by = 'reach_id', all.x = T, all.y = F,
                               allow.cartesian = T) %>%
      setorderv(c('reach_id', 'date'))
  }
  
  # Output -----------------------------
  out_list <- list()
  if (scale %in% c('all', 'site')) {
    out_list <- append(out_list, list(site = hydromod_dt_sites))
  }
  if (scale %in% c('all', 'drn')) {
    out_list <- append(out_list, list(drn = relF_dt))
  }
  
  return(out_list)
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
  
  #Proportion of days over overall 90th percentile in the past c(10, 30, 60, 90, 120, 180) days
  in_hydromod_dt[, paste0("oQ10_", rollingstep, "past") := 
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

#------ dist_proj  -------------------------------------------------
# Define standard two-point equidistance projection for a given bounding box
#https://gis.stackexchange.com/questions/313721/automatically-get-an-adequate-projection-system-based-on-a-bounding-box
## distance projection (tpeqd - two-point equidistant) with projection parameters 
## derived from feature extent
dist_proj <- function(x) {
  bb <- sf::st_bbox(x)
  paste0("+proj=tpeqd +lat_1=",
         bb[2],
         " +lon_1=", 
         bb[1],
         " +lat_2=",
         bb[4], 
         " +lon_2=", 
         bb[3],
         " +x_0=0",
         " +y_0=0",
         " +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
}

#------ snap_points_inner ------------------------------------------------------
snap_points_inner <- function(in_pts,
                              in_target,
                              sites_idcol,
                              attri_to_join=NULL
) {
  #Snap points (fastest custom way in R, it seems):
  #first computing a line between site and snapping place on nearest segment
  sitesnap_l <- terra::nearest(in_pts, in_target, centroids = F, lines = T)
  values(sitesnap_l) <- values(in_pts)
  sitesnap_l$snap_dist_m <- perim(sitesnap_l)
  
  #convert the line to a point (the line's end point)
  sitesnap_p <- terra::as.points(sitesnap_l) %>%
    .[duplicated(values(.)[, sites_idcol]),]
  
  #Join attributes of nearest line to that point
  if (!is.null(attri_to_join)) {
    if ('all' %in% attri_to_join) { 
      sitesnap_p[, names(in_target)] <- terra::nearby(
        sitesnap_p, in_target, k=1, centroids=FALSE)[,'k1'] %>% #Could grab the nth nearest or place a distance limit
        as.data.frame(in_target)[.,] 
    } else {
      sitesnap_p[, attri_to_join] <- terra::nearby(
        sitesnap_p, in_target, k=1, centroids=FALSE)[,'k1'] %>%
        as.data.frame(in_target)[., attri_to_join] 
    }
  }
  
  return(sitesnap_p)
}

#------ split_sp_line----------------------------------------------------------
# Original author: Miguel Porto
# #From https://github.com/miguel-porto/fix-streams
# splits a Lines object (or coordinate matrix) in n segments of length length
#(starting in the begining) plus the remaining segment (what is left)

split_sp_line <- function(line, n, length, debug = F) {
  # splits a Lines object (or coordinate matrix) in n segments of length length (starting in the begining) plus the remaining segment (what is left)
  if (debug) plot(line)
  coo <- sp::coordinates(line)
  
  if (inherits(coo, "list")) {
    if (length(coo) > 1) {
      stop("Multiple lines not allowed")
    } else {
      coo <- coo[[1]]
      if (!inherits(coo, "matrix")) {
        if (!inherits(coo, "list")) stop("Invalid line object")
        if (length(coo) > 1) stop("Multiple lines not allowed")
        coo <- coo[[1]]
      }
    }
  } else {
    if (!inherits(coo, "matrix")) stop("Invalid line object")
  }
  
  pieces <- list()
  accum <- 0
  i <- 1
  remainder <- 0
  newcoords <- matrix(nc = 2, nr = 0)
  
  repeat {
    newcoords <- rbind(newcoords, coo[i, ])
    v <- c(coo[i + 1, ] - coo[i, ])
    hyp <- sqrt(sum(v^2))
    accum <- accum + hyp
    if (accum > length && length(pieces) < n) { # cut in this segment
      oript <- coo[i, ]
      # 			accum=hyp
      repeat {
        newpt <- oript + v / hyp * (length - remainder)
        newcoords <- rbind(newcoords, newpt)
        pieces <- c(pieces, list(newcoords))
        newcoords <- matrix(newpt, nr = 1)
        remainder <- 0
        if (debug) points(newpt[1], newpt[2], pch = 19)
        accum <- accum - length
        if (accum < length || length(pieces) >= n) break
        oript <- newpt
      }
      remainder <- accum
    } else {
      remainder <- accum
    }
    
    i <- i + 1
    if (i >= dim(coo)[1]) {
      newcoords <- rbind(newcoords, coo[dim(coo)[1], ])
      pieces <- c(pieces, list(newcoords))
      break
    }
  }
  
  if (debug) {
    for (i in 1:length(pieces)) {
      pieces[[i]][, 1] <- pieces[[i]][, 1] - 20 * i
      lines(pieces[[i]], col = "red")
      points(pieces[[i]][1, 1], pieces[[i]][1, 2])
    }
  }
  
  return(pieces)
}

#------ fix_confluences_inner ---------------------------------------------------------
# Fixes complex confluences in stream networks
# Original author: Miguel Porto (only cosmetic changes were performed)
# From https://github.com/miguel-porto/fix-streams/blob/master/fix-streams.r
# All nodes which have >2 streams flowing to it are corrected. The outermost streams' end vertices
# are adjusted by "step" meters along the downgoing stream. New nodes are created at suitable places,
# and existing lines suitably split.
#*** Requires a SpatialLinesDataFrame with the proper FROM_NODE and TO_NODE fields. 
# The network is assumed to be correct in all other aspects, there is no error checking.
###### USAGE EXAMPLE
# rios=readOGR("streams_Pt.shp","streams_Pt")
# correctedshp=fix.streams(rios,step=10)
# writeOGR(correctedshp,"streams_corrected.shp","streams_corrected","ESRI Shapefile")

fix_confluences_inner <- function(shp, from = "FROM_NODE", to = "TO_NODE", 
                                  step = 10, fields_to_keep=NULL) {
  # step is the desired length (in map units) by which the river sinks are adjusted (separated) downstream.
  pieces <- list()
  probrivers <- list()
  removeindexes <- integer(0)
  CRS <- shp@proj4string
  
  # find multiple confluence nodes
  nv <- table(shp@data[, to])
  mc <- as.numeric(names(nv[nv > 2])) # these are the nodes with >2 rivers flowing to	them
  maxnode <- max(c(shp@data[, to], 
                   shp@data[, from])) # max node ID, for creating new nodes
  
  if (length(mc) > 0) {
    cat(length(mc), "complex confluences found.\n")
    cat("Cutting lines and tweaking vertices...\n")
    flush.console()
    
    for (i in mc) { # for each problematic node
      # msrc=shp[shp@data[,to] %in% i,]	# get source rivers flowing to it
      privers <- which(shp@data[, to] == i) # get source rivers flowing to it (problematic rivers)
      sinkriver <- which(shp@data[, from] == i) # get sink river
      msrc <- shp[privers, ]
      mto <- shp[sinkriver, ]
      delta <- sum(shp@data[, to] %in% i) - 2 # how many problematic rivers flow to there? leave only two, the others correct
      newnodes <- c(mto@data[, from], (maxnode + 1):(maxnode + delta), mto@data[, to]) # the IDs of the nodes that will be created (the first remains the same for the two "good" rivers)
      
      coo <- sp::coordinates(mto)[[1]][[1]] # coordinates of the sink river
      
      # order the rivers by their angle, so that the outermost rivers are first adjusted (alternating the side)
      v1 <- coo[1, ] - coo[2, ]
      v2 <- matrix(nc = 2, nr = length(msrc))
      for (j in 1:length(msrc)) { #For each source river
        tmp <- sp::coordinates(msrc[j, ])[[1]][[1]] #Get coordinates of all of its nodes
        v2[j, ] <- tmp[dim(tmp)[1] - 1, ] - tmp[dim(tmp)[1], ] #Compute diff in lon and lat between the second to last and last nodes
      }
      ang <- atan2(v2[, 2], v2[, 1]) - atan2(v1[2], v1[1]) #Compute angle between source rivers and sink river
      ang <- (ang + pi) %% (2 * pi) - pi
      # ang[ang > 0] <- ang[ang > 0] %% pi
      # ang[ang < 0] <- -((-ang[ang < 0]) %% pi)
      names(ang) <- privers
      ang <- ang[order(abs(ang), decreasing = T)]
      
      if (sum(ang > 0) > sum(ang < 0)) {
        angneg <- c(as.numeric(names(ang[ang < 0])), rep(NA, sum(ang > 0) - sum(ang < 0)))
        angpos <- as.numeric(names(ang[ang >= 0]))
      } else {
        angpos <- c(as.numeric(names(ang[ang >= 0])), rep(NA, sum(ang < 0) - sum(ang > 0)))
        angneg <- as.numeric(names(ang[ang < 0]))
      }
      angs <- matrix(c(angpos, angneg), nc = 2)
      privers <- na.omit(as.vector(t(angs)))[1:delta]
      
      # split sink river in delta pieces plus the remainder
      tmplines <- split_sp_line(line = coo, 
                                n = delta, 
                                length = step,
                                debug = F)
      
      if (length(tmplines) <= delta) stop("You must decrease step: river ", sinkriver)
      
      for (j in 1:length(tmplines)) {
        # cut sink river into pieces, as many as necessary
        rid <- runif(1, 10^6, 10^7) # random ID for the piece
        # create a new piece with the j'th (step+1) vertices of the sink river
        piece <- sp::SpatialLines(list(sp::Lines(list(sp::Line(tmplines[[j]])), rid)), 
                                  proj4string = CRS)
        newdata <- mto@data
        newdata[1, ] <- NA
        newdata[1, from] <- newnodes[j]
        newdata[1, to] <- newnodes[j + 1]
        rownames(newdata) <- rid
        
        #Re-assign kept data 
        if (!is.null(fields_to_keep)) {
          newdata[[fields_to_keep]] <- mto[[fields_to_keep]]
        }
        
        pieces <- c(pieces, 
                    list(sp::SpatialLinesDataFrame(piece, newdata, match = F))) # save pieces for later use
        
        if (j > 1) {
          # now change coords of problematic rivers
          pri <- privers[j - 1] # pick the j'th problematic river
          tmp <- sp::coordinates(shp[pri, ])[[1]][[1]]
          tmp[dim(tmp)[1], ] <- tmplines[[j]][1, ] # change the coordinate of the last vertex of problematic river
          tmp1 <- sp::SpatialLines(list(sp::Lines(list(sp::Line(tmp)), 
                                                  shp[pri, ]@lines[[1]]@ID)), 
                                   proj4string = CRS) # keep same ID (original will be removed)
          tmp1 <- sp::SpatialLinesDataFrame(tmp1, shp[pri, ]@data)
          tmp1@data[1, to] <- newnodes[j]
          probrivers <- c(probrivers, list(tmp1)) # collect new rivers to replace old
        }
      }
      
      removeindexes <- c(removeindexes, c(privers[1:delta], sinkriver))
      maxnode <- maxnode + delta
    }
    
    cat("Now reassembling shape...\n")
    flush.console()
    newlines <- pieces[[1]]
    
    for (j in 2:length(pieces)) {
      newlines <- rbind(newlines, pieces[[j]])
    }
    
    for (j in 1:length(probrivers)) {
      newlines <- rbind(newlines, probrivers[[j]])
    }
    
    # remove all problematic + sink rivers
    newshp <- shp[-removeindexes, ]
    newshp <- rbind(newshp, newlines)
  } else {
    print("No complex confluences found. Returning input network unchanged.")
    newshp <- shp
  }
  
  return(newshp)
}

#------ compute_cor_matrix_inner -----------------------------------------------
# Function to compute and flatten correlations
compute_cor_matrix_inner <- function(in_dt, group_vars = NULL,
                                     x_cols, y_cols = NULL,
                                     correlation_type = "spearman",
                                     exclude_diagonal = TRUE) {
  
  # If y_cols is NULL, compute correlation within x_cols
  if (is.null(y_cols)) {
    y_cols <- x_cols
    single_group <- TRUE
  } else {
    single_group <- FALSE
  }
  
  # Compute correlation within each group (or overall if group_vars is NULL)
  cor_results <- in_dt[, {
    if (single_group) {
      cor_result <- Hmisc::rcorr(as.matrix(.SD[, x_cols, with=FALSE]), type = correlation_type)
    } else {
      cor_result <- Hmisc::rcorr(as.matrix(.SD[, x_cols, with=FALSE]),
                                 as.matrix(.SD[, y_cols, with=FALSE]),
                                 type = correlation_type)
    }
    
    # Flatten within the data.table operation, using consistent variable names
    combinations <- CJ(variable1 = rownames(cor_result$r),
                       variable2 = colnames(cor_result$r))
    list(
      variable1 = combinations$variable1,
      variable2 = combinations$variable2,
      correlation = cor_result$r[cbind(combinations$variable1, combinations$variable2)],
      p_value = cor_result$P[cbind(combinations$variable1, combinations$variable2)]
    )
    
  }, by = group_vars]
  
  if (exclude_diagonal) {
    cor_results <- cor_results[variable1 != variable2,]
  }
  
  return(cor_results)  # Return the result directly
}
#------ fill_nas_hierarchical --------------------------------------------------
# Function to fill NAs hierarchically (by site -> by country -> overall)
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
#------ trans_pca --------------------------------------------------------------
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
  pca_out <- dt_copy[, stats::prcomp(.SD, center=FALSE, scale=FALSE), 
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
#------ trans_pca_wrapper ------------------------------------------------------
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
    best_sim = c(8, 9, 20, 3, 8, 15),
    all_sims_filename = c(
      "Butiznica_2022-12-15_option0.nc", #_run8_final
      "Velicka_2023-02-01_option0.nc",
      "Lepsamaanjoki_2022-12-16_option0.nc",
      'Albarine_2022-12-16_option0.nc',
      "Bukkosdi_2022-12-16_option0.nc", #_run8_final
      "Genal_2023-01-18_option0.nc"
    )
  ) %>%
    .[, `:=`(all_sims_path = file.path(in_hydromod_dir, 'data', catchment, 
                                       "Results_present_period", all_sims_filename),
             sel_sim_path = file.path(in_hydromod_dir, 
                                      'flow_intermittence_discharge_baseflow_1960_2022',
                                      paste0(catchment, 
                                             '_flow_intermittence_1960_oct2022_corr.nc')),
             catchment_path = file.path(in_hydromod_dir, 'data',catchment,
                                        "watershed_small_catchment.shp"),
             network_path = file.path(in_hydromod_dir, 'data', catchment,
                                      "river_network.shp"),
             reaches_path = file.path(in_hydromod_dir,'data', catchment, "reach.par"),
             sites_reachids = file.path(in_hydromod_dir, 'data', catchment,
                                        paste0(catchment, 
                                               '_sampling_sites_ReachIDs.csv')))]
  
  return(hydro_drn_paths_dt)
}

#------ get_genal_drainage_area -----------------------------------
#in_flowdir_path <- tar_read(flowdir_hydrosheds90m_path)
# outdir <- file.path(resdir, 'gis')
# in_sites_snapped <- tar_read(site_snapped_gpkg_list)

get_genal_drainage_area <- function(in_flowdir_path, outdir) {
  out_upa_path <- file.path(outdir, 'genal_hydrosheds90m_upa.tif')
  
  #Try computing upstream drainage area based on Modelo Digital de Elevaciones (MDE) de Andalucia
  # if (!file.exists(out_upa_path)) {
  #   dem <- terra::rast(in_dem_path) %>%
  #     project(y="epsg:25830", method = 'bilinear')
  #   flowdir <- terra::terrain(dem, v = 'flowdir', neighbors = 8)
  #   cell_area <- terra::cellSize(flowdir, unit = "km")
  #   upstream_area <- terra::flowAccumulation(flowdir, weight = cell_area, 
  #                                            filename = out_upa_path)
  # } else {
  #   upstream_area <- terra::rast(out_upa_path)
  # }
  #But too raw, would need to be hydrologically conditioned, etc.
  
  #Use HydroSHEDS 90m flow direction instead
  if (!file.exists(out_upa_path)) {
    flowdir <- terra::rast(in_flowdir_path) %>%
      crop(ext(-6, -4, 36, 37))
    
    cell_area <- terra::cellSize(flowdir, unit = "km", transform=TRUE)
    upstream_area <- terra::flowAccumulation(flowdir, weight = cell_area, 
                                             filename = out_upa_path)
  }
  
  #Extract area for each site by visual match
  genal_sites_upa_dt <- data.table(
    site = paste0('GEN',str_pad(seq(1, 26), width=2, pad='0')),
    basin_area_km2 = c(5.57, 5.80, 60.15, 0.75, 7.00, 14.38, 73.71, 88.83,
                       2.17, 13.70, 8.04, 142.69, 7.54, 5.91, 157.33, 3.48,
                       17.03, 37.57, 8.75, 13.75, 2.89, 250.77, 22.40, 277.11,
                       320.47, 336.78
    )
  )
  
  return(genal_sites_upa_dt)
}

#------ get_drn_hydromod -------------------------------------------------------
get_drn_hydromod <- function(hydromod_path, varname, selected_sims=NULL) {
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
  #No basin area for Spain - solved
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
  
  #Make sure velocity is OK
  env_dt_merged[avg_velocity_macroinvertebrates > 5, avg_velocity_macroinvertebrates := NA]
  
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
  env_dt_merged[running_id == 'BUK36_1', 
                riparian_cover_in_the_riparian_area := 90]
  
  #Use corrected data for some records
  env_dt_merged[running_id == 'BUK42_3',
                bankfull_at_max_wetted_width_m := 4.7] #from original data sheet uploaded on Teams
  
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
  
  #Convert dates from character to Date format
  env_dt_merged[, date :=  as.Date(date, format='%d.%m.%Y')]
  
  return(env_dt_merged)
}
#------ read_biodt -------------------------------------------------------------
# path_list = tar_read(bio_data_paths)
# in_metadata_edna = tar_read(metadata_edna)

read_biodt <- function(path_list, in_metadata_edna, include_bacteria=T) {
  #Read and name all data tables
  dt_list <- mapply(function(in_path, in_name) {
    print(in_path)
    out_dt <- fread(in_path) %>% 
      .[, organism := in_name] %>%
      setnames(tolower(names(.))) 
    
    if ('country' %in% names(out_dt)) {
      out_dt[country == 'Czech Republic', country := 'Czech']
    }
    
    if (is.character(out_dt$date)) {
      out_dt[, date := as.Date(date, '%d.%m.%Y')]
    }
    
    if (all(c('site', 'campaign') %in% names(out_dt))) {
      out_dt[, running_id := paste0(site, '_', campaign)]
    }
    
    return(out_dt)
  }, path_list, names(path_list)
  )
  
  #Remove spaces in running id
  in_metadata_edna[, running_id := gsub(' ', '', running_id)]

  #Fill NAs in dates with eDNA metadata if possible
  dt_list$dia_sedi[
    is.na(date),
    date := in_metadata_edna[sample_type=='sediment', .(running_id, date)][
      .SD, on='running_id', x.date]]
  dt_list$fun_sedi[
    is.na(date),
    date := in_metadata_edna[sample_type=='sediment', .(running_id, date)][
      .SD, on='running_id', x.date]] 
  dt_list$fun_biof[
    is.na(date),
    date := in_metadata_edna[sample_type=='biofilm', .(running_id, date)][
      .SD, on='running_id', x.date]] 
  
  #Replace dates in dia_biof, which seem erroneous
  dt_list$dia_biof[
    ,
    date := in_metadata_edna[sample_type=='biofilm', .(running_id, date)][
      .SD, on='running_id', x.date]] 

  #Add Campaign and Site to bacteria data, then remove pool sites
  dt_list$bac_sedi[, c('site', 'campaign') := tstrsplit(v1, '_')] %>%
    setnames('v1', 'running_id') %>%
    .[, campaign := as.integer(campaign)] 
  dt_list$bac_biof[, c('site', 'campaign') := tstrsplit(v1, '_')] %>%
    setnames('v1', 'running_id') %>%
    .[, campaign := as.integer(campaign)]
  
  # #Standardize country names
  # country_standard <- data.table(
  #   original = c("CRO", "FRA", "SPA", "CZ",  "HUN", "FIN"),
  #   new = c("Croatia", 'France', "spain", "Czech", "Hungary", 'Finland')
  # )
  # in_metadata_edna <- merge(in_metadata_edna, country_standard,
  #                           by.x='country', by.y='original') %>%
  #   .[, country := new] %>%
  #   .[, `:=`(new = NULL)]
  
  #Remove bacteria in pools
  dt_list$bac_sedi <- merge(dt_list$bac_sedi, 
                            in_metadata_edna[sample_type=='sediment', 
                                             .(running_id, date, habitat, country)],
                            by='running_id') %>%
    .[country == 'Czech Republic', country := 'Czech']
  
  dt_list$bac_biof <- merge(dt_list$bac_biof,
                            in_metadata_edna[sample_type=='biofilm',
                                             .(running_id, date, habitat, country)],
                            by='running_id') %>%
    .[country == 'Czech Republic', country := 'Czech']
  
  dt_list$bac_biof_nopools <- dt_list$bac_biof %>%
    .[habitat!='pool',] %>%
    .[, `:=`(habitat = NULL,
             organism = 'bac_biof_nopools')]
  dt_list$bac_biof[, habitat := NULL]
  
  
  dt_list$bac_sedi_nopools <- dt_list$bac_sedi %>%
    .[habitat!='pool',] %>%
    .[, `:=`(habitat = NULL,
             organism = 'bac_sedi_nopools')]
  dt_list$bac_sedi[, habitat := NULL]
  
  #Correct typo in date for fungi  
  dt_list$fun_sedi[date == as.Date("2012-02-25"), 
                   date := as.Date("2021-02-25")]

  dt_list$fun_biof[date == as.Date("2012-02-25"), 
                   date := as.Date("2021-02-25")]
  
  if (!include_bacteria) {
    dt_list <- dt_list[!grepl('bac_', names(dt_list))]
  }
  
  return(dt_list)
}

#------ calc_spdiv -------------------------------------------------------------
# in_country <- 'Croatia'
# in_biodt <- tar_read(bio_dt)[['dia_biof']][country == in_country,]
# in_metacols <- metacols
# level='local'

#calc_spdiv(in_biodt, in_metacols, level = 'local')

#Compute nestedness and turnover based on temporal beta diversity between t and t-1
#https://www.rdocumentation.org/packages/adespatial/versions/0.3-24/topics/beta.div.comp
comp_richrepl_inner <- function(dt, spcols, beta_div_coef, quant) {
  sub_dt <- copy(dt)
  beta_div_format <- beta_div_coef
  if (quant) {
    beta_div_format <- ifelse(beta_div_coef == 'J', 'R', 'BC')
  }
  #J: Jaccard, S: SOrense, R: Ruzicka, O: Bray-Curtis
  
  repl_rich_list <- adespatial::beta.div.comp(
    as.matrix(sub_dt[, spcols, with=F]),
    coef = beta_div_coef, quant = quant)
  
  out_dt <- sub_dt[2:.N, `:=`(
    repl_tm1 ={as.matrix(repl_rich_list$repl) %>% #Get replacement compared to t-1; off-diagonal t2 to t1, t3 to t2, etc. 
        .[row(.)==(col(.)+1)]},
    rich_tm1 = {as.matrix(repl_rich_list$rich) %>% 
        .[row(.)==(col(.)+1)]},
    tm1 = {as.matrix(repl_rich_list$D) %>% 
        .[row(.)==(col(.)+1)]}
  )] %>%
    .[, c('repl_tm1', 'rich_tm1', 'tm1'), with=F] %>%
    cbind(data.table(t(repl_rich_list$part))) %>%
    setnames(paste0(beta_div_format , names(.)))
  
  out_dt[, campaign :=  sub_dt$campaign]
}

calc_spdiv <- function(in_biodt, in_metacols, level = 'local') {
  #Get metadata columns (all except species data)
  metacols_sub <- names(in_biodt)[names(in_biodt) %in% in_metacols]
  spcols <- names(in_biodt)[!(names(in_biodt) %in% in_metacols)]
  
  #Set order
  setorderv(in_biodt, c('campaign', 'site'))
  
  #Compute and partition taxonomic gamma diversity 
  #between alpha and beta for the overall period, for individual sites
  #and for the metacom
  spxp <- as.matrix(in_biodt[,-metacols_sub, with=F])
  bio_structure <- data.frame(space = as.factor(in_biodt$site),
                              time = as.factor(in_biodt$campaign))
  z = rowSums(spxp > 0)
  decomp <- try(HierAnodiv(spxp = spxp[z>0,], 
                           structure = bio_structure[z>0,], 
                           phy = NULL, weight = NULL, check = T, q = 1))
  
  if (level == 'local') {
    #Compute species richness
    biodt_div <-in_biodt[, list(richness = rowSums(.SD > 0)
    ), by=.(site, campaign, date, organism), .SDcols = spcols] %>%
      .[, mean_richness := mean(richness), by=site]
    
    #Keep only sites x campaigns with species
    biodt_copy <- merge(in_biodt, 
                        biodt_div[richness>0, .(site, campaign)],
                        by=c('site', 'campaign'),
                        all.x = F
    )
    
    #Fourth-root transform data
    biodt_copy[, (spcols) := lapply(.SD, function(x) x^(1/4)), 
               .SDcols = spcols]
    
    # ggplot(biodt_melt, aes(x=(value^(1/4)))) +
    #   geom_histogram() +
    #   scale_x_continuous() +
    #   facet_wrap(~country)
    # 
    # ggplot(biodt_melt, aes(x=log10(value+1))) +
    #   geom_histogram() +
    #   scale_x_continuous() +
    #   facet_wrap(~country)
    
    #Compute invsimpson Index and inverse Simpson (alpha diversity) by site x time step
    sha_simp <- biodt_copy[, list(
      campaign = campaign,
      site = site,
      shannon = vegan::diversity(as.matrix(.SD), index = "shannon"),
      invsimpson = vegan::diversity(as.matrix(.SD), index = "invsimpson")
    ), .SDcols = spcols]
    
    #Compute nestedness and turnover based on temporal beta diversity
    multicampaign_sites <- biodt_copy[, .N, by=site][N>1, site]
    J_richrepl <- biodt_copy[site %in% multicampaign_sites, 
                             comp_richrepl_inner(dt = .SD, spcols = spcols,
                                                 beta_div_coef = 'J', quant=F), 
                             by=site]
    R_richrepl <- biodt_copy[site %in% multicampaign_sites,
                             comp_richrepl_inner(dt = .SD, spcols = spcols,
                                                 beta_div_coef = 'J', quant=T),
                             by=site]
    
    #Format local diversity decomposition
    localdiv_decomp <- as.data.table(decomp$tab_by_site) %>%
      setnames(c('nsite'), 'ncampaigns')
    localdiv_decomp[, site := unique(biodt_copy$site)]
    
    out_dt <- mergeDTlist(
      list(biodt_div, sha_simp, J_richrepl, R_richrepl), 
      by=c('site', 'campaign'), all=T, sort = T, set_suffix=F) %>%
      merge(localdiv_decomp, by='site', all.x=T)
    
  } else if (level == 'regional') {
    nsites <- nrow(decomp$tab_by_site) #Number of sites
    nsteps <- max(decomp$tab_by_site[,'nsite']) #Number of sampling campaigns
    
    out_dt <- t(decomp[[1]]) %>%
      data.table %>%
      setnames(c('Gamma', 'Beta1', 'Beta2in1', 'mAlpha'),
               paste0(c('gamma', 'beta_s', 'beta_t_s', 'malpha'), '_drn')
      ) %>% 
      .[, `:=`(organism = in_biodt[1, .(organism)][[1]],
               beta_s_drn_std = (beta_s_drn-1)/(nsites-1),
               beta_t_s_drn_std = (beta_t_s_drn-1)/(nsteps-1) 
               )] 
  }
  
  return(out_dt)
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
    scale_color_manual(values=c("#8c510a", "#bf812d", "#01665e", "#80cdc1", "#8073ac", "#543005")) +
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
    scale_color_manual(values=c("#8c510a", "#bf812d", "#01665e", "#80cdc1", "#8073ac", "#543005")) +
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
    , network_sub_path := file.path(
      out_dir, 
      paste0(tolower(country), '_river_network_sub_',
             format(Sys.time(), "%Y%m%d"),'.gpkg'))]
  
  in_hydromod_paths_dt[, {
    if (!file.exists(network_sub_path) | overwrite) {
      terra::vect(network_path) %>%
        .[relate(terra::vect(catchment_path), .,"intersects")[1,],] %>% #Contains removes some segments
        terra::writeVector(filename = network_sub_path,
                           overwrite = T)
    }
  }, by=country]
  
  return(in_hydromod_paths_dt[, stats::setNames(network_sub_path, country)])
}

#------ clean_network ------------------------------
# in_country <- 'Croatia'
# rivnet_path <- tar_read(network_sub_gpkg_list)[[in_country]]
# idcol <- 'cat'
# node_clustering_dist = 50
# min_segment_length = 20
# outdir = file.path(resdir, 'gis')
# save_gpkg = TRUE

clean_network <- function(rivnet_path, idcol, 
                          node_clustering_dist,
                          min_segment_length = 20,
                          outdir=NULL, save_gpkg=FALSE, 
                          return_path=FALSE) {
  #Read input network
  rivnet <- st_read(rivnet_path) %>%
    st_cast("LINESTRING") %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') 
  
  #Preformat basic network
  sfnet_ini <- rivnet %>%
    as_sfnetwork %>%
    activate("edges") %>%
    filter(!edge_is_multiple()) %>% #Keep shortest of edges that connect the same pair of nodes
    filter(!edge_is_loop()) #Remove obvious loops: edges that start and end at the same node
  
  #------------------ Split lines at intersections -----------------------------
  #Get confluence nodes (nodes of third degree: with at least 3 intersecting edges)
  #2nd degree nodes are pseudonodes and 1st degree nodes are dangling
  splitting_nodes <- activate(sfnet_ini, nodes) %>%
    mutate(degree = igraph::degree(.)) %>%
    filter(degree >= 3) %>%
    st_as_sf("nodes")
  
  #Visualize interactively
  # ggplotly(
  #   ggplot(rivnet) +
  #     geom_sf() +
  #     geom_sf(data=splitting_nodes, color='red')
  # )
  #Write split nodes to double check
  #st_write(splitting_nodes, 'split_nodes.gpkg')
  
  #Simplify network by first fully dissolving and then splitting at confluences
  rivnet_agg <- sf::st_union(rivnet) %>%
    sf::st_line_merge(.) %>%
    lwgeom::st_split(., splitting_nodes)  %>%
    sf::st_collection_extract("LINESTRING")
  
  #------------------ Remove loops by clustering nearby nodes--------------------
  sfnet <- as_sfnetwork(rivnet_agg)
  
  node_coords <- sfnet %>%
    activate("nodes") %>%
    st_coordinates()
  
  # Cluster the nodes with the DBSCAN spatial clustering algorithm.
  # We set eps = 40 such that:
  # Nodes within a distance of 40 m from each other will be in the same cluster.
  # We set minPts = 1 such that:
  # A node is assigned a cluster even if it is the only member of that cluster.
  clusters = dbscan::dbscan(node_coords, 
                            eps = node_clustering_dist, minPts = 1)$cluster
  
  # Add the cluster information to the nodes of the network.
  clustered = sfnet %>%
    activate("nodes") %>%
    mutate(cls = clusters)
  
  #contracts groups of nodes based on cluster number as a grouping variable. 
  #The geometry of each contracted node is the centroid of the original 
  #group members’ geometries. Moreover, the geometries of the edges that start
  #or end at a contracted node are updated such that their boundaries match
  #the new node geometries
  sfnet_clustered <- tidygraph::convert( #
    clustered,
    sfnetworks::to_spatial_contracted,
    cls,
    simplify = TRUE
  ) %>%
    tidygraph::convert(sfnetworks::to_spatial_smooth) #Remove pseudo nodes
  
  #------------------ Trim spits (segments with dangle point under minimum length ) ------------
  #Convert back to sf
  rivnet_clustered <- st_as_sf(sfnet_clustered, "edges")
  
  #Re-aggregate and split the network
  rivnet_clustered_agg <- sf::st_union(rivnet_clustered) %>%
    sf::st_line_merge(.) %>%
    lwgeom::st_split(., splitting_nodes)  %>%
    sf::st_collection_extract("LINESTRING") %>%
    st_sf(geometry=.)
  
  #Compute segment length
  rivnet_clustered_agg$length <- as.numeric(st_length(rivnet_clustered_agg))
  
  #Identify dangle points
  agg_endpts <- c(lwgeom::st_startpoint(rivnet_clustered_agg), 
                  lwgeom::st_endpoint(rivnet_clustered_agg)) 
  
  danglepts <- agg_endpts[!(duplicated(agg_endpts) | 
                              duplicated(agg_endpts, fromLast=TRUE)),] %>%
    st_sf %>%
    mutate(dangle = 'dangle')
  
  #Remove spits under min length
  rivnet_clustered_aggsub <- st_join(x = rivnet_clustered_agg,
                                     y = danglepts, 
                                     join = st_intersects,
                                     left = TRUE,
                                     all.x = TRUE) %>%
    mutate(dangle = replace_na(dangle, "connected")) %>%
    dplyr::filter(!((dangle == 'dangle') &
                      (length < min_segment_length))) %>%
    select(-c(length, dangle))
  
  #------------------ Re-assign original IDs to all segments ---------------
  #Split back into component linestrings
  rivnet_endpts <- c(lwgeom::st_startpoint(rivnet), 
                     lwgeom::st_endpoint(rivnet))
  rivnet_clustered_resplit <- lwgeom::st_split(rivnet_clustered_aggsub, 
                                               rivnet_endpts)  %>%
    sf::st_collection_extract("LINESTRING") %>%
    dplyr::distinct(.) #Remove duplicate lines
  #.[-unlist(st_equals(., retain_unique=TRUE)),] 
  
  rivnet_clustered_resplit$UID <- seq_along(rivnet_clustered_resplit$geometry)
  
  #Join those that have not moved
  rivnet_clustered_joinini <- st_join(
    rivnet_clustered_resplit,
    rivnet,
    join = st_equals,
    suffix = c(".ini", ".clustered"),
    left = TRUE,
    largest = FALSE
  )
  
  #Convert those without match to points every 10 meters
  rivnet_nomatch_pts <- rivnet_clustered_joinini %>%
    .[is.na(rivnet_clustered_joinini[[idcol]]),] %>%
    st_line_sample(density = 1/10) %>%
    st_sf(geometry = ., crs = st_crs(rivnet_clustered_joinini))
  rivnet_nomatch_pts$UID <- rivnet_clustered_joinini[
    is.na(rivnet_clustered_joinini[[idcol]]),][['UID']]
  rivnet_nomatch_pts <- st_cast(rivnet_nomatch_pts, 'POINT')
  
  #For each point, get id of nearest line in initial river network
  nearest_segix <- rivnet_nomatch_pts %>%
    st_nearest_feature(., rivnet)
  rivnet_nomatch_pts$nearest_id <- rivnet[nearest_segix,][[idcol]]
  
  #Compute distance to that nearest line
  rivnet_nomatch_pts$nearest_dist <- st_distance(
    rivnet_nomatch_pts, rivnet[nearest_segix,], 
    by_element=TRUE)
  
  #Keep line for which the sum of the inverse distance to the points is greatest              
  rivnet_nomatch_selid <- as.data.table(rivnet_nomatch_pts) %>%
    .[, list(inverse_dist_sum = sum(1/nearest_dist)), 
      by=.(UID, nearest_id)] %>%
    .[, list(nearest_id=.SD[which.max(inverse_dist_sum), 
                            nearest_id]), 
      by=UID] 
  
  #Fill NAs with those
  rivnet_clustered_joinall <- merge(rivnet_clustered_joinini,
                                    rivnet_nomatch_selid, 
                                    by='UID', all.x=T)
  rivnet_clustered_joinall[[idcol]] <- dplyr::coalesce(
    rivnet_clustered_joinall[[idcol]],
    rivnet_clustered_joinall$nearest_id)
  
  #------------------ Write out results ------------------------------------------
  out_net <- rivnet_clustered_joinall[
    , c('UID', names(rivnet)[names(rivnet) != 'geometry'])]
  
  out_path <- file.path(outdir,
                        paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                               '_clean',
                               format(Sys.time(), "%Y%m%d"),
                               '.gpkg')
  )
  
  if (save_gpkg) {
    st_write(out_net, out_path, append=F)
  }
  
  return(ifelse(return_path, out_path, out_net))
}




#------ manual_clean_croatia ---------------------------------------------------
manual_clean_croatia <- function(in_net) {
  line_to_edit <- in_net[in_net$UID==651,]$geom 
  in_net[in_net$UID==651,]$geom <- st_sfc(st_linestring(
    rbind(
      c(X=596968.793, Y=4899716.907),
      st_coordinates(line_to_edit)[-1, c('X', 'Y')]
    )),
    crs = st_crs(in_net)) %>%
    st_reverse()
  
  in_net <- st_snap(in_net, in_net, tolerance = 0.005)
  
  #Split connected line at intersection
  line_to_split <- in_net[in_net$UID==622,]
  sink_parts <- st_collection_extract(
    st_split(line_to_split$geom, 
             st_sfc(st_point(c(X=596968.793, Y=4899716.907)), 
                    crs = st_crs(in_net))
    ),"LINESTRING")
  in_net[in_net$UID==622,]$geom <- sink_parts[1,]
  line_to_split$UID <- max(in_net$UID) + 1
  line_to_split$geom <- sink_parts[2,]
  out_net <- rbind(in_net, line_to_split)
  
  return(out_net)
}
#------ direct_network -----------------------------------------
# Define helper functions
get_endpoint <- function(line) st_coordinates(line)[nrow(st_coordinates(line)), ]
get_startpoint <- function(line) st_coordinates(line)[1, ]

direct_network_inner <- function(segment, in_network, idcol, visited = NULL) {
  # #Reverse upstream segments recursively
  # visited <- NULL
  # segment <- rivnet[rivnet[[idcol]] == 22,] 
  # in_network = rivnet
  
  #print(segment[[idcol]])
  visited <- c(visited, segment[[idcol]])
  
  # Find connected segments 
  inseg_startpoint <- get_startpoint(segment$geometry)
  
  connected <- in_network[!(in_network[[idcol]] %in% visited),] %>%
    .[as.vector(st_intersects(., lwgeom::st_startpoint(segment$geometry), 
                              sparse=F)),]
  
  # ggplotly(ggplot(in_network[!(in_network[[idcol]] %in% visited),]) +
  #            geom_sf() +
  #            geom_sf(data=segment, color='red') +
  #            geom_sf(data=lwgeom::st_startpoint(segment$geometry)))
  
  # Reverse and recurse
  for (i in seq_len(nrow(connected))) {
    seg_id <- connected[[idcol]][i]
    seg_geom <- connected$geometry[i]
    # Reverse if not aligned
    if (!(all.equal(get_endpoint(seg_geom), inseg_startpoint) == TRUE)) {
      in_network[in_network[[idcol]]==seg_id,]$geometry <- st_reverse(seg_geom)
    }
    # Recursively process upstream
    in_network <- direct_network_inner(
      segment = in_network[in_network[[idcol]]==seg_id,], 
      in_network, idcol, visited)
  }
  
  return(in_network)
}


# outlet_uid_list <- list(Croatia = 458,
#                         Czech = 4,
#                         Finland = 682,
#                         France = 1,
#                         Hungary = 5,
#                         Spain = 86
# )
# 
# in_country <- 'Croatia'
# rivnet_path <- tar_read(network_clean_gpkg_list)[[in_country]]
# idcol <- 'UID'
# outletid <- outlet_uid_list[[in_country]]
# outdir = file.path(resdir, 'gis')
# save_gpkg = TRUE

direct_network <- function(rivnet_path, idcol,
                           outletid, outdir=NULL, 
                           save_gpkg=FALSE) {
  #Read input network
  rivnet <- st_read(rivnet_path) %>%
    st_cast("LINESTRING") %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') 
  
  #Interactively plot network
  # (ggplot(rivnet) +
  #   geom_sf(size = 0.1,
  #           arrow = arrow(angle = 30,
  #                         length = unit(0.075, 'inches'),
  #                         ends = "last",
  #                         type = "closed")))
  
  #Identify dangle points
  agg_endpts <- c(lwgeom::st_startpoint(rivnet), 
                  lwgeom::st_endpoint(rivnet)) 
  
  danglepts <- agg_endpts[!(duplicated(agg_endpts) | 
                              duplicated(agg_endpts, fromLast=TRUE)),] %>%
    st_sf %>%
    mutate(dangle = 'dangle')
  
  #Make sure outlet segment is in the right direction
  outlet_endpt <- lwgeom::st_endpoint(
    rivnet[rivnet[[idcol]] == outletid,])
  
  if (!(outlet_endpt %in% danglepts$geometry)) {
    rivnet[rivnet[[idcol]] == outletid, 'geometry'] <- st_reverse(
      rivnet[rivnet[[idcol]] == outletid, 'geometry'])
  }
  
  outlet_seg <- rivnet[rivnet[[idcol]] == outletid,]
  
  # Apply to the entire network
  out_net <- direct_network_inner(segment = outlet_seg, 
                                  in_network = rivnet, 
                                  idcol = idcol,
                                  visited = NULL)
  
  #------------------ Write out results ------------------------------------------
  out_path <- file.path(outdir,
                        paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                               '_directed',
                               format(Sys.time(), "%Y%m%d"),
                               '.gpkg')
  )
  
  if (save_gpkg) {
    st_write(out_net, out_path, append=F)
  }
  
  return(out_path)
}

#------ fix_complex_confluences ------------------------------------------------
# in_country <- 'Croatia'
# rivnet_path = tar_read(network_directed_gpkg_list)[[in_country]]
# outdir = file.path(resdir, 'gis')
# max_node_shift = 5

#Network must be directed & topologically correct aside from the complex confluences
fix_complex_confluences <- function(rivnet_path, max_node_shift = 5,
                                    outdir=NULL, out_path=NULL) {
  # Read input network
  rivnet <- st_read(rivnet_path) %>%
    st_cast("LINESTRING") %>%
    st_set_geometry("geometry")
  
  #Compute from-to fields
  net_fromto <- as_sfnetwork(rivnet) %>%
    activate(edges) %>%
    as.data.table %>%
    .[, .(from, to, UID)]
  
  rivnet_fromto <- merge(rivnet, net_fromto, by='UID', all.x=T)
  
  #Run fix streamlines from Miguel Porto
  rivnet_fixed <- fix_confluences_inner(shp = as_Spatial(rivnet_fromto), 
                                        from = "from",
                                        to = "to", 
                                        step = max_node_shift,
                                        fields_to_keep = 'cat') %>%
    .[, !(names(.) %in% c('from', 'to'))] %>%
    st_as_sf
  
  rivnet_fixed[is.na(rivnet_fixed$UID),]$UID <- max(rivnet_fixed$UID, na.rm=T) + 
    seq_len(sum(is.na(rivnet_fixed$UID)))
  
  #------------------ Write out results ------------------------------------------
  if (is.null(out_path)) {
    out_path <- file.path(outdir,
                          paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                                 '_conflufixed',
                                 format(Sys.time(), "%Y%m%d"),
                                 '.gpkg')
    )
  }
  
  st_write(rivnet_fixed, out_path, append=F)
  
  return(out_path)
}

#------ assign_strahler_order --------------------------------------------------
# in_country <- 'Croatia'
# in_rivnet = network_ssnready_gpkg_list[[in_country]]
# idcol = 'UID'

assign_strahler_order <- function(in_rivnet, idcol, verbose = F) {
  if (is.character(in_rivnet)) {
    rivnet <- st_read(in_rivnet)
    
    #Compute from-to fields
    net_fromto <- as_sfnetwork(rivnet) %>%
      activate(edges) %>%
      as.data.table %>%
      .[, c('from', 'to', idcol), with=F] %>%
      merge(.[, list(nsource = .N), by=to], 
            by.x='from', by.y='to', all.x=T) %>%
      .[is.na(nsource), nsource := 0]
    
    #Join to spatial network
    rivnet_fromto <- merge(rivnet, net_fromto, by=idcol)
    
    #Convert to data.table for speed and syntax
    rivnet_fromto_dt <- as.data.table(rivnet_fromto)
    
  } else if (is.data.table(in_rivnet)) {
    rivnet_fromto_dt <- merge(in_rivnet,
                              in_rivnet[, list(nsource = .N), by=to], 
                              by.x='from', by.y='to', all.x=T) %>%
      .[is.na(nsource), nsource := 0]
  }
  
  #Assign strahler 1 to lines with no source line
  rivnet_fromto_dt[nsource == 0, strahler := 1]
  
  # Compute Strahler order iteratively
  #while (sum(is.na(rivnet_fromto_dt$strahler)) > 417) {
  while (any(is.na(rivnet_fromto_dt$strahler))) {
    if (verbose) { print(sum(is.na(rivnet_fromto_dt$strahler)))}
    # Identify lines whose sources' Strahler orders are all assigned
    rivnet_fromto_dt <-  merge(rivnet_fromto_dt, rivnet_fromto_dt[, .(to, strahler)], 
                               by.x='from', by.y='to', suffixes = c('_down', '_up'),
                               all.x=T)
    
    #Extend strahler order downstream for consecutive sections on the same segment 
    #(only one source section)
    rivnet_fromto_dt[is.na(strahler_down) & !is.na(strahler_up) & nsource == 1,
                     strahler_down := strahler_up]
    
    #Flag as eligible those segments for which all upstream segments have
    #a strahler order
    rivnet_fromto_dt[is.na(strahler_down) & nsource >= 2 & !is.na(strahler_up),
                     eligible := ((.N>=2) & .N==nsource), by=idcol]
    
    #Identify max upstream strahler and number of upstream sections with that strahler order
    rivnet_fromto_dt[eligible & !is.na(eligible), 
                     max_strahler_u := max(strahler_up),
                     by=idcol]
    rivnet_fromto_dt[eligible & !is.na(eligible),
                     n_max_strahler_u := .SD[strahler_up==max_strahler_u, .N],
                     by=idcol]
    
    rivnet_fromto_dt[eligible & !is.na(eligible),
                     strahler_down := fifelse(
                       n_max_strahler_u >= 2,
                       max_strahler_u + 1,
                       max_strahler_u),
                     by = idcol]
    
    setnames(rivnet_fromto_dt, 'strahler_down', 'strahler')
    
    rivnet_fromto_dt <- rivnet_fromto_dt[
      !duplicated(get(idcol)), 
      -c('strahler_up', 'eligible', 'max_strahler_u', 'n_max_strahler_u'), 
      with=F]
  }
  
  return(rivnet_fromto_dt[
    , c(idcol, 'from', 'to', 'strahler', 'nsource'), with=F])
}

#------ reassign_netids --------------------------------------------------
#Problem: 'cat' or 'reach_id' in the shapefile, associated with lines across 
#confluences. Makes it impossible to correctly use network hydrology data.

#There are three spatial units in this code:
#Segment sections: uniquely identified by UID, these are the individual lines
#                 contained in the shapefile generated by WP1 after hydrological
#                 modeling.
#Segments: uniquely identifed by UID_fullseg, representing lines extending between 
#           two confluences. They often contain multiple sections.
#Cat: uniquely identified by cat (or reach_id), these are unique reaches as modeled
#     by the hydrological model, but not well represented/assigned to segment
#     sections  in the shapefile. Some may be omitted by the shapefile, 
#     many span multiple segment sections and even multiple segments across 
#     confluences and Strahler stream order.

#The goal of this function is to re-assign a topologically/hydrologically
#logical cat to each segment section to be able to join the hydrological model 
#outputs to the shapefile.

#Helper functions
#For those segments where only one cat remains, assign cat to the entire segment
#at least a given proportion of the length of the full segment
assign_singlecat_to_seg <- function(net_dt, length_threshold = 1/3) {
  # Identify segments with only one remaining section
  #with a cat (length_uid_fullsegper == 1) and where the actual length of that
  #section is at least a given proportion of the length of the full segment
  #including all sections.
  cat_toassign <- net_dt[
    is.na(cat_cor) & (length_uid_fullsegper == 1) &
      (length_uid / length_fullseg > length_threshold),
    .(UID_fullseg, cat)
  ] %>%
    setnames('cat', 'cat_cor')
  
  #Assign those cat_cor to all sections within those segments
  net_dt[is.na(cat_cor), 
         cat_cor := cat_toassign[.SD, on = 'UID_fullseg', cat_cor]]
}

#For segment sections with cat==NAs and cat_cor==NA, 
#get cat_cor from upstream section in the same segment
get_nearby_cat <- function(in_dt, source_direction, 
                           strahler_list = seq(1,15)) {
  if (source_direction == 'upstream') {
    source_col <- 'to' #Column of section to get cat_cor from
    sink_col <- 'from' #Column of section to assign cat_cor to
  } else if (source_direction == 'downstream') {
    source_col <- 'from' #Column of section to get cat_cor from
    sink_col <- 'to' #Column of section to assign cat_cor to
  }
  
  # Identify sections to fill
  sections_tofill <- in_dt[
    , which(is.na(cat) & is.na(cat_cor) & (strahler %in% strahler_list))]
  
  # Prepare cat_cor values to assign
  cat_cor_toassign <- in_dt[
    !is.na(cat_cor) & 
      (strahler %in% strahler_list) & 
      (get(source_col) %in% in_dt[sections_tofill, unique(get(sink_col))]), 
    c(source_col, 'UID_fullseg', 'cat_cor'), with=F] %>%
    setnames(source_col, sink_col)
  
  # Assign cat_cor to sections_tofill using a join
  in_dt[sections_tofill,
        cat_cor := cat_cor_toassign[
          .SD, on=c(sink_col, 'UID_fullseg'), cat_cor]]
}

#Remove pseudo-nodes among segments sections provided equal attributes
#Re-compute from-to and length_uid
remove_pseudonodes <- function(in_net, equal_cols = FALSE, 
                               summarise_attributes ='first') {
  out_net <- as_sfnetwork(in_net) %>%
    activate(edges) %>%
    convert(to_spatial_smooth, 
            require_equal = equal_cols,
            summarise_attributes = summarise_attributes) %>%
    st_as_sf() %>%
    mutate(length_uid = as.numeric(st_length(.)))
  out_net
}

#Parameters
# in_country <- 'Croatia'
# rivnet_path <- tar_read(network_nocomplexconf_gpkg_list)[[in_country]]
# strahler_dt <- tar_read(network_strahler)[[in_country]]
# in_reaches_hydromod_dt <- tar_read(reaches_dt)[country==in_country,]

reassign_netids <- function(rivnet_path, strahler_dt, 
                            in_reaches_hydromod_dt, outdir,
                            country = NULL, in_ext='gpkg') {
  #---------- Prepare data -------------------------------------------------------
  #Read network and join with hydromod data
  rivnet <- st_read(rivnet_path) %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') %>%
    merge(strahler_dt, by='UID') 
  
  #Remove pseudonodes to define full segments between confluences
  rivnet_fullseg <- as_sfnetwork(rivnet) %>%
    activate("edges") %>%
    convert(to_spatial_smooth) %>%
    st_as_sf() 
  
  #Give these full segments between confluences a separate ID: UID_fullseg
  rivnet_fullseg$UID_fullseg <- seq_len(nrow(rivnet_fullseg))
  
  #Intersect with original network to check length of overlap for each 'cat'
  rivnet_fullseg_inters <- st_intersection(
    rivnet, rivnet_fullseg[, c('geometry', 'UID_fullseg')]) %>%
    .[st_geometry_type(.) != 'POINT',] %>%
    st_cast('MULTILINESTRING') %>%
    st_cast('LINESTRING')
  
  #Remove pseudo-nodes among segments sections of the same cat
  #and compute length of segment sections (UID)
  rivnet_fullseg_inters_nopseudo <- remove_pseudonodes(
    in_net = rivnet_fullseg_inters, 
    equal_cols = c("cat", "UID_fullseg"), 
    summarise_attributes='first')
  
  #Create data.table to work with
  rivnet_inters_dt <- as.data.table(rivnet_fullseg_inters_nopseudo) 
  
  #Back up cat before editing it
  rivnet_inters_dt[, cat_copy := cat]
  
  #Compute other lengths
  rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]
  rivnet_inters_dt[, length_fullseg := sum(length_uid), by=UID_fullseg]
  
  #Compute the percentage length of that cat ("hydromod" reach) that is 
  #represented by that segment section (identified by UID)
  rivnet_inters_dt[, length_uid_catper := length_uid/length_cat]
  #Check whether this is the largest segment section for that cat 
  rivnet_inters_dt[, length_uid_catmax := (length_uid == max(length_uid)), by=cat]
  #Compute the percentage length of that full segment that is represented by
  #that segment section (identified by UID)
  rivnet_inters_dt[, length_uid_fullsegper := length_uid/length_fullseg]
  
  #Computer number of full segments that a cat overlaps
  rivnet_inters_dt[, n_seg_overlap := length(unique(UID_fullseg)), by=cat]
  
  #Format network data from hydrological model
  reaches_hydromod_format <- in_reaches_hydromod_dt[
    , .(ID, to_reach, length, strahler)] %>%
    setnames(c('ID_hydromod', 'to_reach_hydromod', 
               'length_hydromod', 'strahler_hydromod')) %>%
    .[, `:=`(from = ID_hydromod, to = to_reach_hydromod)] %>%
    merge( #Compute actual strahler order
      assign_strahler_order(in_rivnet = ., idcol = 'ID_hydromod', verbose = F),
      by = 'ID_hydromod'
    ) %>%
    setnames(c('strahler', 'nsource'),
             c('strahler_hydromod_recalc', 'nsource_hydromod'))
  
  #---------- Initial assignment of correct cat -------------------------------------------------
  #For first order segments where the cat is only represented by that segment,
  #keep that cat for this section of the segment (see cat==2938 for Croatia)
  rivnet_inters_dt[n_seg_overlap == 1 & strahler == 1, 
                   cat_cor := cat]
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For full segments that only have one cat associated with them, confirm cat
  rivnet_inters_dt[length_uid == length_fullseg, 
                   cat_cor := cat]
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For segment sections with cat==NAs and cat_cor==NA in first order streams, 
  #get cat_cor from upstream segment in same segment
  get_nearby_cat(rivnet_inters_dt, 
                 source_direction = 'upstream', 
                 strahler_list = 1) 
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Computer number of full segments that a cat overlaps
  rivnet_inters_dt[, n_seg_overlap := length(unique(UID_fullseg)), by=cat]
  
  #For first and second order segments where the cat is now only represented 
  #by that segment, keep that cat for this section of the segment
  rivnet_inters_dt[n_seg_overlap == 1 & strahler <= 2 & is.na(cat_cor), 
                   cat_cor := cat]
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For segment section with cat==NAs in first order streams, assign cat_cor of
  #upstream segment
  get_nearby_cat(rivnet_inters_dt, 
                 source_direction = 'upstream', 
                 strahler_list = 1) 
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For sections that represent over 75% of the total length of that cat, assign cat_cor
  rivnet_inters_dt[length_uid_catper>0.7 & is.na(cat_cor), cat_cor := cat]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  #For those segments where only one cat remains, assign cat to the entire segment
  assign_singlecat_to_seg(rivnet_inters_dt)
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  
  #For cats that are already assigned, remove their "cat" from other segments
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #For segment sections in strahler==1 & n_seg_overlap > 2, cat := NA 
  rivnet_inters_dt[strahler == 1 & n_seg_overlap > 2 & is.na(cat_cor), 
                   cat := NA] 
  
  #if main representative of cat, assign cat_cor
  rivnet_inters_dt[!is.na(cat_cor), cat := cat_cor]
  rivnet_inters_dt[, length_cat := sum(length_uid), by=cat]
  rivnet_inters_dt[, length_uid_catper := length_uid/length_cat, by=cat]
  rivnet_inters_dt[length_uid_catper > 0.7 & is.na(cat_cor), cat_cor := cat]
  
  #For cats that are already assigned, remove their "cat" from other segments
  cats_assigned <- rivnet_inters_dt[!is.na(cat_cor), unique(cat_cor)]
  rivnet_inters_dt[cat %in% cats_assigned & is.na(cat_cor), cat := NA]
  
  #Re-compute length_uid_fullsegper, excluding segment sections whose initial
  #cat was assigned elsewhere
  rivnet_inters_dt[!is.na(cat_cor), 
                   length_uid_fullsegper := length_uid/sum(length_uid),
                   by=UID_fullseg]
  
  ##For those segments where only one cat remains, assign cat to the entire segment
  #regardless of the proportion of the segment it represents
  cat_cor_toassign <-  rivnet_inters_dt[
    is.na(cat_cor) & length_uid_fullsegper == 1,
    .(UID_fullseg, cat)] %>%
    setnames('cat', 'cat_cor')
  rivnet_inters_dt[is.na(cat_cor), 
                   cat_cor := cat_cor_toassign[.SD, on='UID_fullseg', 
                                               cat_cor]]
  
  #For the remaining segments where all cats are NAs (have been assigned elsewhere), 
  #assign the cat of the section representing the highest percent of the segment
  #These are usually near loops that have been removed 
  rivnet_inters_dt[is.na(cat_cor), NAlength := sum(length_uid),
                   by=UID_fullseg]
  rivnet_inters_dt[(NAlength==length_fullseg), 
                   cat_cor := fifelse(length_uid == max(length_uid), cat_copy, NA), 
                   by=UID_fullseg]
  
  #For all remaining sections, remove their "cat" 
  rivnet_inters_dt[is.na(cat_cor), cat := NA]
  
  for (i in 1:3) {
    #For segment section with cat==NAs and cat_cor==NA, 
    #assign cat_cor of upstream segment of same strahler order
    get_nearby_cat(rivnet_inters_dt, source_direction = 'upstream')
    
    #For is.na(cat) & is.na(cat_cor), 
    #assign downstream cat_cor of same UID_fullseg
    get_nearby_cat(rivnet_inters_dt, source_direction = 'downstream') 
  }
  
  #---------- Merge back with spatial data and remove pseudo nodes ---------------
  #Merge processed attributes with pre-formatted spatial network
  rivnet_catcor <- merge(
    rivnet_fullseg_inters_nopseudo, 
    rivnet_inters_dt[, .(UID, cat_cor, cat_copy, length_uid)],
    by=c('UID', 'length_uid')) %>%
    .[, c('UID', 'cat_cor', 'cat_copy', 'strahler', 
          'nsource', 'UID_fullseg', 'geometry')]
  
  #Remove pseudo-nodes among segments sections of the same cat_cor
  #re-compute length uid
  rivnet_catcor_smooth <- remove_pseudonodes(
    in_net = rivnet_catcor, 
    equal_cols = c("cat_cor", "UID_fullseg"))
  
  #st_write( rivnet_catcor_smooth[,c('from', 'to', 'UID', 'cat_cor', 'cat_copy', 'geometry')], 'results/rivnet_catcor_smooth.shp')
  
  #---------- Adjust results based on topology from hydrological model network ---
  #Merge with hydrological model network topology data
  rivnet_catcor_hydromod <- merge(
    rivnet_catcor_smooth,
    reaches_hydromod_format[
      , .(ID_hydromod, to_reach_hydromod, length_hydromod, 
          strahler_hydromod_recalc, nsource_hydromod)],
    by.x = 'cat_cor', by.y = 'ID_hydromod',
    all.x = T
  ) %>%
    as.data.table
  
  if (rivnet_catcor_hydromod[is.na(cat_cor), .N] == 0) {
    #Remove cats that are not supposed to be downstream of any other cat
    #when there are sections of the correct stream order on that segment
    rivnet_catcor_hydromod[, n_cats_fullseg := length(unique(cat_cor)),
                           by=UID_fullseg]
    
    rivnet_catcor_hydromod[n_cats_fullseg > 1 & 
                             nsource_hydromod == 0 & nsource > 0,
                           `:=`(cat = NA, cat_cor  = NA)]
    
    #assign cat_cor of upstream segment of same strahler order
    for (i in 1:3) {
      get_nearby_cat(in_dt=rivnet_catcor_hydromod, source_direction = 'upstream')
    }
    get_nearby_cat(in_dt=rivnet_catcor_hydromod, source_direction = 'downstream')
    
    rivnet_catcor_hydromod <- st_as_sf(rivnet_catcor_hydromod) %>%
      remove_pseudonodes(equal_cols = c("cat_cor", "UID_fullseg")) %>%
      as.data.table
  }
  
  #---------- Check results and correct based on hydrological model topology --------
  #Compare downstream segment cat
  to_reach_shpcor <-  rivnet_catcor_hydromod[
    , .(from, cat_cor)] %>%
    setnames(c('to', 'to_reach_shpcor'))
  
  rivnet_catcor_hydromod <- merge(rivnet_catcor_hydromod, 
                                  to_reach_shpcor, by='to', all.x=T) %>%
    .[, hydromod_shpcor_match := (to_reach_hydromod == to_reach_shpcor)]
  
  #If two upstream sections both should point to the same reach but do not
  #correct if downstream cat is in network topology but missing in shapefile
  #or if the assigned cat_cor is duplicated
  dupli_catcor <- rivnet_catcor_hydromod[, .N, by=cat_cor][N>=2,]$cat_cor
  missing_catcor <- reaches_hydromod_format[
    !(ID_hydromod %in% unique(rivnet_catcor_hydromod$cat_cor)),]$ID_hydromod
  
  #Remove first-order sections with no upstream sections whose downstream cat_cor
  #and strahler order does not match network topology from hydrological model
  rivnet_catcor_hydromod <- rivnet_catcor_hydromod[
    !(hydromod_shpcor_match==FALSE & nsource == 0 & strahler_hydromod_recalc>1),]
  
  #Identify pairs of reaches that are both upstream of the wrong reach
  double_downstream_mismatch_correct <- rivnet_catcor_hydromod[
    (to_reach_hydromod != to_reach_shpcor),
    list(n1=.N, to, to_reach_shpcor), by=to_reach_hydromod] %>%
    .[n1>=2, list(n2=.N, to_reach_hydromod, to_reach_shpcor), by=to] %>%
    .[n2>=2 & ((to_reach_hydromod %in% missing_catcor) |
                 (to_reach_shpcor %in% dupli_catcor)),] %>%
    unique %>%
    setnames(c('to_reach_shpcor', 'to'), c('cat_cor', 'from'))
  
  rivnet_catcor_hydromod[from %in% double_downstream_mismatch_correct$from &
                           cat_cor %in% double_downstream_mismatch_correct$cat_cor
                         , cat_cor := double_downstream_mismatch_correct[.SD, on=c('cat_cor', 'from'), 
                                                                         to_reach_hydromod]]
  
  #---------- Implement a few manual corrections  ---------------------------------
  if (country == 'Croatia') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      filter(!(UID %in% c(489, 90))) %>% #cat_cor 2082, and 2480, respectively) %>%
      mutate(cat_cor = case_match(
        UID,
        103 ~ 1776,  #UID 103 (catcor 1832) -> catcor 1776
        116 ~ 1832, #UID 116 (catcor 1836) -> catcor 1832
        31 ~ 1788, #UID 31 (catcor 1792) -> catcor 1788
        1045 ~ 2898, #UID 1045 (cat_cor 2908) -> catcor 2898
        .default = cat_cor
      )
      )
  } else if (country == 'Czech') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      mutate(cat_cor = case_match(
        UID,
        298 ~ 40232, #UID 298 (catcor 40126) -> catcor 40232
        158 ~ 40258, #UID 158 (40260) -> catcor 40258
        .default = cat_cor
      )
      )
  } else if (country == 'Finland') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      mutate(cat_cor = case_match(
        UID,
        538 ~ 1049400, #UID 538 (catcor 1051800) -> catcor 1049400
        392 ~ 1069001, #UID 392 (catcor 1086800) -> catcor 1069001
        397 ~ 1069001, #UID 397 (catcor 1086600) -> catcor 1069001 
        .default = cat_cor
      )
      )
  } else if (country == 'France') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      filter(!(UID %in% c(529, 630, 642, 991))) %>% #2489000, 2475800, 2476200, 2475800
      mutate(cat_cor = case_match(
        UID,
        205 ~ 2422600, #UID 205  (catcor 2497800) -> catcor 2422600
        693 ~ 2444000, #UID 693 (catcor 2465400) -> catcor 2444000  
        805 ~ 2467600, #UID 805 (catcor 2467800) -> catcor 2467600
        500 ~ 2485400, #UID 500 (catcor 2485600) -> catcor 2485400
        .default = cat_cor
      )
      )
  } else if (country == 'Hungary') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod %>%
      mutate(cat_cor = case_match(
        UID,
        51 ~ 650800,  #UID 51 (catcor 651000) -> catcor 650800
        331 ~ 652001,  #UID 331 (catcor 652200) -> catcor 652001
        379 ~ 673000,  #UID 379 (catcor 673400) -> catcor 673000
        .default = cat_cor
      )
      )
  } else if (country == 'Spain') {
    rivnet_catcor_manual <- rivnet_catcor_hydromod 
    rivnet_catcor_manual[rivnet_catcor_manual$UID == 25, 'cat_cor'] <- 5426 #UID 25 (catcor 5428) -> catcor 5426
  }
  
  
  #------ Final processing -------------------------------------------------------
  #Remove pseudo-nodes among segments sections of the same cat_cor
  #re-compute length uid, and re-assign from-to
  out_rivnet<- remove_pseudonodes(
    in_net = st_as_sf(rivnet_catcor_manual), 
    equal_cols = c("cat_cor", "UID_fullseg"))
  
  #re-assign downstream segment cat
  to_reach_shpcor_dt <-  as.data.table(out_rivnet)[
    , .(from, cat_cor)] %>%
    setnames(c('to', 'to_reach_shpcor'))
  
  out_rivnet <- merge(
    as.data.table(out_rivnet)[, -c('to_reach_shpcor', 'to_reach_hydromod'), with=F], 
    to_reach_shpcor_dt, by='to', all.x=T) %>%
    merge(reaches_hydromod_format[, .(ID_hydromod, to_reach_hydromod)],
          by.x='cat_cor', by.y='ID_hydromod', all.x=T) %>%
    .[, hydromod_shpcor_match := (to_reach_hydromod == to_reach_shpcor)]
  
  #---------- Write out results ------------------------------------------------------
  out_path <- file.path(outdir,
                        paste0(tools::file_path_sans_ext(basename(rivnet_path)),
                               '_reided',
                               format(Sys.time(), "%Y%m%d"),
                               '.', in_ext)
  )
  
  #Export results to gpkg
  write_sf(st_as_sf(out_rivnet)[
    , c('UID', 'strahler','length_uid', 'UID_fullseg', 'cat_cor', 'from', 'to',
        'to_reach_shpcor', 'to_reach_hydromod', 'hydromod_shpcor_match')],
    out_path)
  
  return(out_path)
}
#------ compute_hydrostats_drn -------------------------------------------------
# in_drn <- 'Czech'
# varname <-  'qsim' #'isflowing' #qsim
# in_sites_dt <- tar_read(sites_dt)[country == in_drn,]
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[in_drn]]
# in_hydromod_drn <- tar_read_raw((paste0('hydromod_dt_', in_drn, '_', varname)))
# in_network_idcol = 'cat'

compute_hydrostats_drn <- function(in_network_path,
                                   in_sites_dt,
                                   varname,
                                   in_hydromod_drn,
                                   in_network_idcol = 'cat') {
  setDT(in_hydromod_drn$data_all)
  setDT(in_hydromod_drn$dates_format) 
  
  #-------------------- Compute intermittence statistics -----------------------
  if (varname == 'isflowing') {
    #Import network shapefiles and get their length
    network_v <- terra::vect(in_network_path)
    network_v$reach_length <- terra::perim(network_v)
    reach_dt <- as.data.table(network_v[, c(in_network_idcol, 'reach_length')]) %>%
      setnames(in_network_idcol, 'reach_id') %>%
      .[, list(reach_length = sum(reach_length)), by=reach_id]
    remove(network_v)
    
    #Merge hydro data and reach length for computing network-wide statistics
    #qdat[, unique(reach_id)[!(unique(reach_id) %in% reach_dt$reach_id)]] 
    #Two IDs are not in the shapefile for Finland? maybe a hydrological unit not associated
    intermod_dt <- merge(in_hydromod_drn$data_all, reach_dt, 
                         by='reach_id', all.x=F, all.y=F) %>%
      merge(in_hydromod_drn$dates_format[, .(date, month, hy, doy)],
            by='date', all.x=T, all.y=F)
    
    if ('nsim' %in% names(intermod_dt)) {
      q_stats <- list()
      #Compute hydrological statistics for entire DRN
      q_stats$drn <- intermod_dt[, 
                                 compute_hydrostats_intermittence(
                                   in_hydromod_dt = .SD,
                                   in_sites_dt = in_sites_dt,
                                   scale = 'drn')$drn, 
                                 by = nsim] 
      
      #Compute hydrological statistics for individual sites
      q_stats$site <- intermod_dt[, 
                                  compute_hydrostats_intermittence(
                                    in_hydromod_dt = .SD,
                                    in_sites_dt = in_sites_dt,
                                    scale = 'site')$site
                                  , by = nsim]
    } else {
      q_stats <-compute_hydrostats_intermittence(
        in_hydromod_dt = intermod_dt,
        in_sites_dt = in_sites_dt)
    }
    
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
      merge(in_sites_dt[, .(site, reach_id)], .,
            by = 'reach_id', all.x = T, all.y = F,
            allow.cartesian = T) %>%
      setorderv(c('reach_id', 'date'))
  }
  
  return(q_stats)
}

#------ subset_hydrostats ------------------------------------------------------
# in_country <- in_drn <- 'Spain'
# hydrostats <- tar_read_raw(paste0("hydrostats_", in_country, '_qsim'))
# in_bio_dt <- tar_read(bio_dt)

subset_hydrostats <- function(hydrostats, in_country, in_bio_dt) {
  unique_sampling_dates <- lapply(in_bio_dt, function(org_dt) {
    org_dt[country==in_country, .(date)]
  }) %>% rbindlist %>% unique %>% .$date
  
  min_date <- min(unique_sampling_dates, na.rm=T)
  max_date <- max(unique_sampling_dates, na.rm=T)
  
  if ('drn' %in% names(hydrostats)) {
    hydrostats[['site']] <- hydrostats[['site']][date>=min_date & date<=max_date,]
    hydrostats[['drn']] <- hydrostats[['drn']][date>=min_date & date<=max_date,]
  } else {
    hydrostats <- hydrostats[date>=min_date & date<=max_date,]
  }
  return(hydrostats)
}

#------ format_sites_dt ----------------------------------------------------------
# in_country <- 'Spain'
# in_path <- tar_read(hydromod_paths_dt)[country == in_country, sites_reachids]
# #check <- format_site_dt(in_path, in_country)
# in_env_dt <- tar_read(env_dt)

format_site_dt <- function(in_path, in_env_dt, in_country) {
  sites_dt <- fread(in_path) %>%
    setnames(tolower(names(.))) %>%
    .[, reach_id := as.integer(reach_id)]
  
  if (in_country == 'Croatia') {
    sites_dt[, id := sub('BUT', '', id) %>%
               str_pad(width=2, side='left', pad = 0) %>%
               paste0('BUT', .)] %>%
      .[id=='BUT19', reach_id := 2868] #Noticed after correcting the network topology
  } 
  
  if (in_country == 'Czech') {
    sites_dt[, id := sub('S', '', id) %>%
               str_pad(width=2, side='left', pad = 0) %>%
               paste0('VEL', .)]
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
  
  #Fill missing sites info from WP1 with info from field data
  #For Hungary and Spain
  sites_dt_filled <- in_env_dt[drn==in_country, .(site, latitude, longitude)] %>%
    setnames(c('latitude', 'longitude'), c('latitude_field', 'longitude_field')) %>%
    unique %>%
    merge(sites_dt, by.x='site', by.y='id', all.x=T) %>%
    .[is.na(lat), `:=`(lat = latitude_field, lon = longitude_field)] %>%
    .[, `:=`(latitude_field = NULL, longitude_field = NULL)]
  
  sites_dt_filled[, reach_id := fcase(
    site=='BUK01', 657801L,
    site=='BUK36', 665000L,
    site=='GEN05', 5780L,
    site=='GEN11', 5594L,
    site=='GEN12', 5588L,
    site=='GEN17', 5882L,
    site=='GEN23', 5732L,
    site=='GEN26', 5392L,
    site=='VI01', 2438600L, #based on documentation in WP1 
    site=='AL01', 2455400L, #based on documentation in WP1 
    default=reach_id
  )]
  
  return(sites_dt_filled[, .(site, lat, lon, reach_id)])
}

#------ create_sites_gpkg --------------------------------------------------------
# in_hydromod_paths_dt = tar_read(hydromod_paths_dt)
# in_sites_dt = tar_read(sites_dt)
# out_dir = file.path('results', 'gis')
# geom = 'reaches'
# overwrite = TRUE
# in_network_path_list = tar_read(network_ssnready_gpkg_list)

create_sites_gpkg <- function(in_hydromod_paths_dt,
                              in_sites_dt,
                              out_dir, 
                              geom,
                              in_network_path_list = NULL,
                              in_network_idcol = 'cat',
                              overwrite = FALSE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  setnames(setDT(in_sites_dt), 'country', 'country_sub')
  
  in_hydromod_paths_dt[
    , sites_gpkg_path := file.path(
      out_dir, 
      paste0(tolower(country), '_site_', geom, 
             format(Sys.time(), "%Y%m%d"), '.gpkg'))]
  
  if (geom == 'reaches' & is.character(in_network_path_list)) {
    #Create site reaches
    in_hydromod_paths_dt[, {
      if (!file.exists(sites_gpkg_path) | overwrite) {
        terra::vect(in_network_path_list[[country]]) %>%
          aggregate(by=in_network_idcol) %>%
          merge(in_sites_dt[country_sub==country,],
                by.x=in_network_idcol, by.y='reach_id', all.x=F) %>%
          terra::writeVector(filename = sites_gpkg_path, 
                             overwrite = T)
      }
    }, by=country]
  }
  
  if (geom == 'points') {
    #Create site points
    in_hydromod_paths_dt[, {
      if (!file.exists(sites_gpkg_path) | overwrite) {
        create_sitepoints_raw(in_dt = in_sites_dt[country_sub==country,], 
                              lon_col = 'lon', lat_col = 'lat',
                              out_points_path = sites_gpkg_path) 
      }
    }, by=country]
  }
  
  return(in_hydromod_paths_dt[, stats::setNames(sites_gpkg_path, country)])
}


#------ snap_sites -------------------------------------------------------------
# drn <- 'Croatia'
# in_sites_path <- tar_read(site_points_gpkg_list)[[drn]]
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[drn]]
# out_snapped_sites_path = NULL
# overwrite = T
# custom_proj = F
# in_sites_unique_id = 'site'
# in_network_unique_id = 'UID'
# in_sites_idcol_tomatch = 'reach_id'
# in_network_idcol_tomatch = 'cat'
# proj_back = F

snap_river_sites <- function(in_sites_path, 
                             in_network_path,
                             out_snapped_sites_path=NULL, 
                             custom_proj = F,
                             proj_back = F,
                             in_sites_unique_id = 'site',
                             in_network_unique_id = 'UID',
                             in_sites_idcol_tomatch = 'reach_id',
                             in_network_idcol_tomatch = 'cat',
                             overwrite = F) {
  
  if (is.null(out_snapped_sites_path)) {
    out_snapped_sites_path <- sub(
      '[.](?=(shp|gpkg)$)', 
      paste0('_snap', format(Sys.time(), "%Y%m%d"), '.'),
      in_sites_path, perl=T)
  }
  
  if (!file.exists(out_snapped_sites_path) | overwrite) {
    if (length(in_sites_path) > 1) {
      sitesp <- do.call(rbind, lapply(in_sites_path, vect)) 
    } else {
      sitesp <- terra::vect(in_sites_path)
    }
    
    #Read network
    target <- terra::vect(in_network_path)
    
    #Project sites 
    # Global datasets tend to be in geographic coordinates. The unit of these 
    # coordinates are decimal degrees, whose west-east length decreases
    # with increasing latitude. Therefore, for identifying the nearest line,
    # which is based on distance calculation, the point dataset needs to be 
    # projected. However, no single projection is valid for the entire planet. 
    # Consequently, for each basin, the sites are projected using a custom 
    # projection which minimizes distortions for distance calculations within the
    # network bounding box.
    if (custom_proj) {
      if (nrow(target) > 1 & ((xmax(target) != xmin(target)) | 
                              (ymax(target) != ymin(target)))
      ){
        target_proj <- terra::project(target, 
                                      dist_proj(target))
      } else {
        #if only one target object, project to UTM
        target_proj <- terra::project(
          target,
          paste0('+proj=utm +zone=', 
                 floor((xmin(target) + 180) / 6) + 1,
                 ' +datum=WGS84 +units=m +no_defs +ellps=WGS84')
        )
      } 
    } else {
      target_proj <- target
    } 
    remove(target)
    
    #Project sites to custom projection
    sitesp_proj <- terra::project(sitesp, crs(target_proj))
    
    #Snap each site to the nearest point on the reach that it is paired with
    sites_unique_id_list <- as.data.frame(sitesp_proj)[[in_sites_unique_id]]
    sitesnap_p <- lapply(
      sites_unique_id_list, 
      function(in_pt_id) {
        pt <- unique(sitesp_proj[sites_unique_id_list == in_pt_id,])
        tar <- target_proj[
          target_proj[[in_network_idcol_tomatch]] == pt[[in_sites_idcol_tomatch]][[1]],]
        
        if (!is.empty(tar) && nrow(pt) > 0) {
          out_p <- snap_points_inner(
            in_pts = pt,
            in_target = tar,
            sites_idcol = in_sites_unique_id,
            attri_to_join = c(in_network_idcol_tomatch, in_network_unique_id)
          )
        } 
        return(out_p)
      }) %>%
      vect(.)
    
    if (proj_back) {
      #Reproject points to original proj
      sitesnap_p <- terra::project(sitesnap_p, sitesp)
    }
    
    #Write it out
    terra::writeVector(sitesnap_p,
                       out_snapped_sites_path,
                       overwrite=overwrite)
  }
  
  return(out_snapped_sites_path) #Path to layer containing site points with attribute data
}

#------ subset_amber -----------------------------------------------------------
# amber_path <- tar_read(amber_path)
# in_hydromod_paths_dt <- tar_read(hydromod_paths_dt)
# out_dir <- file.path(resdir, 'gis')
# overwrite = T

subset_amber <- function(amber_path, in_hydromod_paths_dt, out_dir,
                         overwrite = T) {
  amber_dt <- fread(amber_path)
  
  country_list <- in_hydromod_paths_dt$country
  
  amber_countries<- amber_dt[
    grepl(paste(toupper(country_list), collapse='|'), 
          Country),]
  
  amber_pts_path <- file.path(out_dir,
                              'amber_sub_pts.gpkg')
  
  create_sitepoints_raw(in_dt = amber_countries, 
                        lon_col = 'Longitude_WGS84', lat_col = 'Latitude_WGS84',
                        out_points_path = amber_pts_path) 
  
  amber_vec <- vect(amber_pts_path)
  
  amber_country_list <- lapply(country_list, function(in_country) {
    amber_country_path <- file.path(out_dir, 
                                    paste0(tolower(in_country), 
                                           '_amber_pts.gpkg'))
    if (!file.exists(amber_country_path) | overwrite) {
      catchment_vec <- vect(in_hydromod_paths_dt[country==in_country, catchment_path])
      
      amber_proj <- terra::project(amber_vec, crs(catchment_vec))
      
      amber_proj %>%
        .[ relate(catchment_vec, .,"intersects")[1,],] %>%
        terra::writeVector(filename = amber_country_path,
                           overwrite = T)
      
      return(amber_country_path)
    }
  }) %>% setNames(country_list)
  
  return(amber_country_list)
}



#------ snap_barriers ----------------------------------------------------------
# drn <- 'France'
# in_sites_path <- tar_read(barrier_points_gpkg_list)[[drn]]
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[drn]]
# out_snapped_sites_path = NULL
# overwrite = T
# custom_proj = F

snap_barrier_sites <- function(in_sites_path, 
                               in_network_path,
                               in_sites_idcol,
                               attri_to_join,
                               out_snapped_sites_path=NULL, 
                               custom_proj = T,
                               overwrite = F) {
  
  if (is.null(out_snapped_sites_path)) {
    out_snapped_sites_path <- sub(
      '[.](?=(shp|gpkg)$)', 
      paste0('_snap', format(Sys.time(), "%Y%m%d"), '.'),
      in_sites_path, perl=T)
  }
  
  if (!file.exists(out_snapped_sites_path) | overwrite) {
    #Iterate over every basin where there is a site
    if (length(in_sites_path) > 1) {
      sitesp <- do.call(rbind, lapply(in_sites_path, vect)) 
    } else {
      sitesp <- terra::vect(in_sites_path)
    }
    
    #Read network
    target <- terra::vect(in_network_path)
    
    #Project sites to custom projection
    sitesp_proj <- terra::project(sitesp, crs(target))
    
    sitesnap_p <- snap_points_inner(
      in_pts = sitesp_proj,
      in_target =  target,
      sites_idcol = in_sites_idcol,
      attri_to_join = attri_to_join
    )
    
    #Reproject points to original proj
    sitesnap_p <- terra::project(sitesnap_p, crs(sitesp))
    
    #Write it out
    terra::writeVector(sitesnap_p,
                       out_snapped_sites_path,
                       overwrite=overwrite)
  }
  
  return(out_snapped_sites_path) #Path to layer containing site points with attribute data
}

#------ prepare_data_for_STcon ---------------------------------------------------
# in_country <- in_drn <- 'Croatia'
# in_hydromod_drn <- tar_read(hydromod_comb)[[paste0(
#   "hydromod_dt_", in_country, '_isflowing')]]
# in_net_shp_path <- tar_read(network_ssnready_shp_list)[[in_drn]]

prepare_data_for_STcon <- function(in_hydromod_drn, in_net_shp_path) {
  net <- in_net_shp_path
  
  #Build the adjacency list to built the graph
  net_dt <- net  %>% 
    arrange(from) %>% 
    select(UID, cat, from, to, length_m, to_cat_shp) %>%
    as.data.table
  
  #Identify the outlet (NA in to_cat_shp)
  outlet_to <- net_dt[is.na(net_dt$to_cat_shp),]$to 
  
  #Create the graph from the DRN_adj_list (from - to)
  net_graph <- igraph::graph_from_edgelist(
    as.matrix(net_dt[, .(from, to)]), 
    directed = TRUE) 
  
  # Incorporate a new vertex corresponding the outlet
  # We name the cat 11111 and it will be our last vertex and where the outlet will be directed.  
  nodes_cat <- rbind(net_dt, 
                     data.table(UID = 11111, from = outlet_to), fill=T) %>%
    setorder(from)
  
  # We assign vertex attributes according to the cat_shape, ordered following the "from" value from the cat_shape, 
  # which is not exactly "from anymore" it is an ID of the vertex
  V(net_graph)$UID <- as.character(nodes_cat$UID)
  
  #Compute edge length for distance matrices
  E(net_graph)$weight <- as.integer(net_dt$length_m)
  
  #Compute distance matrix (in integer, lighter to handle -- one-meter difference does not matter)
  river_dist_mat <- igraph::distances(net_graph)
  
  #Format intermittence data
  setDT(in_hydromod_drn$data_all) %>%
    setnames('reach_id', 'cat')
  
  # We create the "End_point" site that will correspond to the 111111 in the flow_intermittence dataset
  # this point will be added with the same frequency of any other reach
  end_point_interm <- in_hydromod_drn$data_all %>%
    .[cat == .[1, cat],] %>% #select whatever reach
    .[, `:=`(UID = 11111,  from = outlet_to, value = 1)] %>% #format
    .[, .(date, UID, isflowing, from)] 
  
  # Merge shapefile with flow intermittence data
  # Merge the Endpoint site and "pivot_wide" the table to obtain the TRUE intermittence table, 
  # where each row corresponds to a day (dates as factors) and columns to all nodes of the network.
  nsims <- isTRUE('nsim' %in% names(in_hydromod_drn$data_all))
  
  cast_formula <- if (nsims) {as.formula('date + nsim ~ from')} else {'date ~ from'}
  sites_status_matrix <- merge(net_dt, in_hydromod_drn$data_all,
                               by='cat', all.x=T) %>%
    rbind(end_point_interm, fill=T) %>%
    .[, isflowing := as.integer(isflowing)] %>%
    setorder(from) %>%
    dcast(formula = cast_formula, value.var = 'isflowing') %>%
    .[!is.na(date),] %>%
    setnafill(type='const', fill=1L) #FILL NAs in original hydrological data
  
  # We built the matrix of the network structure for the STcon, which is the "base" on which connectivity will be assessed. 
  network_structure <- as.matrix(as_adjacency_matrix(net_graph)) 
  
  #Create "reference river" datasets
  interm_ref <- sites_status_matrix %>%
    replace(.==0, 1L) 
  if (nsims) {
    interm_ref <- filter(interm_ref, nsim==1)
  }
  
  return(list(
    sites_status_matrix = sites_status_matrix,
    network_structure = network_structure,
    river_dist_mat = river_dist_mat,
    interm_ref = interm_ref
  ))
}


#------ compute_STcon_rolling ---------------------------------------------------
# in_drn <- 'Croatia'
# in_preformatted_data <- tar_read(preformatted_data_STcon)[[in_drn]]
# in_bio_dt <- tar_read(bio_dt)
# in_nsim <- NULL#tar_read(hydromod_paths_dt)[country == in_drn,]$best_sim
# 
# unique_sampling_dates <- lapply(in_bio_dt, function(org_dt) {
#   org_dt[country==in_drn, .(date)]
# }) %>% rbindlist %>% unique
# in_dates <- unique_sampling_dates
# window <- 10
# direction <- 'undirected'
# weighting <- TRUE
# routing_mode <- 'in'
# output <- 'all'
# verbose <- F
# ref = F

compute_STcon_rolling <- function(in_preformatted_data, ref = F, in_nsim = NULL, 
                                  in_dates, window, output,
                                  direction, routing_mode, weighting, rounding_factor,
                                  verbose = F,
                                  ...) {
  
  #Make sure that date has been modeled
  in_dates <- in_dates[date <= max(in_preformatted_data$sites_status_matrix$date)]
  
  if (isTRUE(in_nsim == 'all')) { #Wrap conditional statement to avoid error in in_nsim == NULL
    in_nsim <- unique(in_preformatted_data$sites_status_matrix$nsim)
  }
  
  #Select reference data or 
  #actual data with intermittence (for which a simulation is selected)
  interm_dt_sel <- if (ref) {in_preformatted_data$interm_ref
  } else {
    if (!is.null(in_nsim)) {
      in_preformatted_data$sites_status_matrix[nsim %in% in_nsim,]
    } else {
      in_preformatted_data$sites_status_matrix
    }
  }
  
  #Compute STcon for every date and previous window of time across simulations
  compute_stcon_date_wrapper <- function(end_date, nsim_sel=NULL) {
    if (verbose) print(end_date)
    
    interm_sub <- interm_dt_sel[
      (date %in% seq(from=end_date-(window-1), to=end_date, by='day')),] %>%
      setorder(date)
    
    if (weighting) {
      value_s_link = value_t_link = 0.1
      value_no_s_link = value_no_t_link = 1
    } else {
      value_no_s_link = value_no_t_link = 0L
      value_s_link = value_t_link = 1L
    }
    
    STcon_list <- compute_stcon(
      sites_status_matrix = as.matrix(interm_sub[, !c('date', 'nsim'), with = FALSE]),
      network_structure = in_preformatted_data$network_structure, 
      direction = direction, 
      routing_mode = routing_mode,
      weighting = weighting,
      dist_matrix = in_preformatted_data$river_dist_mat, # Weighting pairs
      indirect_dispersal = FALSE,
      standardize_neighbors = FALSE,
      value_s_link = value_s_link,
      value_t_link = value_t_link, # Values to links
      value_no_s_link = value_no_s_link,
      value_no_t_link = value_no_t_link, # Values to links
      convert_to_integer = T,
      rounding_factor = rounding_factor, #Because the summed numbers are so big with distances in meters
      output = output,
      verbose = verbose
    ) %>%
      append(
        list(IDs = names(interm_dt_sel[, !c('date', 'nsim'), with = FALSE]))
      )
    
    if (!is.null(nsim_sel)) {
      STcon_list <- append(STcon_list, list(nsim = nsim_sel))
    }
    
    return(STcon_list)
  }
  
  #If multiple nsims to run
  if (!is.null(in_nsim)) {
    STcon_datelist <-  lapply(in_nsim, function(nsim_sel) {
      lapply(in_dates$date, function(end_date) {
        compute_stcon_date_wrapper(end_date, nsim_sel)
      }) %>% setNames(in_dates$date)
    }) %>% setNames(in_nsim)
  } else {
    #otherwise, simply run function for every focus date
    STcon_datelist <- lapply(in_dates$date, function(end_date) { 
      compute_stcon_date_wrapper(end_date)
    }) %>% setNames(in_dates$date)
  }
  
  return(STcon_datelist)
}


#------ postprocess_STcon ------------------------------------------------------
# tar_load(STcon_undirected_list)
# tar_load(STcon_directed_list)
# tar_load(STcon_rolling_ref_list)
# tar_load(preformatted_data_STcon)
# 
# in_country <- 'Croatia'
# in_STcon <- STcon_undirected_list[[in_country]]
# #in_STcon_ref <- STcon_directed_ref_list[[in_country]]
# in_preformatted_data_STcon <- preformatted_data_STcon[[in_country]]
# window_name <- 'STcon_m10'
# date <- '2021-02-25'
# standardize_STcon = FALSE 
# in_STcon_ref = NULL
# in_net_shp_path <- tar_read(network_ssnready_shp_list)[[in_country]]

postprocess_STcon <- function(in_STcon, in_net_shp_path,
                              standardize_STcon = FALSE, in_STcon_ref = NULL) {
  
  #Get original network data
  net_dt <- in_net_shp_path %>%
    as.data.table %>%
    setorder(from)
  
  #Identify the outlet (NA in to_cat_shp)
  outlet_from <- net_dt[is.na(net_dt$to_cat_shp),]$to 
  
  #Compile STcon data in long format by window size, hydrological simulation, date, and ID
  STcon_dt <- lapply(names(in_STcon), function(window_name) {
    lapply(names(in_STcon[[window_name]]), function(date) {
      if (standardize_STcon & !is.null(in_STcon_ref)) {
        #Standardize by STcon in a network without intermittence
        out_STcon <- (in_STcon[[window_name]][[date]]$STcon/
                        in_STcon_ref[[window_name]][[date]]$STcon)
      } else {
        out_STcon <- in_STcon[[window_name]][[date]]$STcon
      }
      
      out_dt <- data.table(variable = window_name,
                           date = as.Date(date),
                           from = in_STcon[[window_name]][[date]]$IDs,
                           stcon_value = out_STcon
      ) %>% .[from != outlet_from,] #Remove outlet
    }) %>% rbindlist
  }) %>% rbindlist %>%
    merge(net_dt[, list(from=as.character(from), UID)], by='from') %>% #Replace from IDs with UID
    .[, from := NULL]
  
  #Compile STcon matrices by window size and date 
  #Prepare UIDs to assign to col/rownames (removing outlet)
  UIDs_order <- as.integer(in_STcon[[1]][[1]][[1]]$IDs) %>%
    setdiff(outlet_from) %>%
    match(net_dt$from) 
  UIDs_to_assign <- net_dt[UIDs_order, UID]
  
  STcon_mat <- lapply(names(in_STcon), function(window_name) {
    lapply(names(in_STcon[[window_name]]), function(date) {
      if (!is.null(in_STcon[[window_name]][[date]]$STconmat)) {
        if (standardize_STcon & !is.null(in_STcon_ref)) {
          out_STconmat <- (in_STcon[[window_name]][[date]]$STconmat/
                             in_STcon_ref[[window_name]][[date]]$STconmat)
        } else {
          out_STconmat <- in_STcon[[window_name]][[date]]$STconmat
        }
        
        #Remove outlet
        out_STconmat <- out_STconmat[rownames(out_STconmat) != 
                                       as.character(outlet_from), 
                                     colnames(out_STconmat) != 
                                       as.character(outlet_from)]
        #Assign UIDs
        colnames(out_STconmat) <- rownames(out_STconmat) <- UIDs_to_assign
        
        data.table(variable = window_name,
                   date = as.Date(date),
                   STconmat = list(out_STconmat)
        )
      }
    }) %>% rbindlist
  }) %>% rbindlist
  
  return(list(
    STcon_dt = STcon_dt,
    STcon_mat = STcon_mat
  ))
}

#------ plot_STcon ------------------------------------------------------------
# tar_load(STcon_directed_formatted)
# in_country <- 'France'
# in_STcon_list <- STcon_directed_formatted[['France']]
# in_date <-in_STcon_list$STcon_dt[1000, date]
# in_net_shp_path = tar_read(network_ssnready_shp_list)[[in_country]]

plot_STcon <- function(in_STcon_list, in_date, in_window=10, 
                       in_net_shp_path, reverse_weighted_stcon = TRUE) {
  
  stcon_sel <- in_STcon_list$STcon_dt[
    (date == in_date) & (variable == paste0('STcon_m', in_window)),]
  
  
  net_v_stcon <- in_net_shp_path %>%
    merge(stcon_sel, by='UID')
  
  net_v_stcon$inverse_stcon <-  1 - (
    (net_v_stcon$stcon_value-min(net_v_stcon$stcon_value))/
      (max(net_v_stcon$stcon_value)-min(net_v_stcon$stcon_value))
  )
  
  out_p <- ggplot(data=net_v_stcon)+
    scale_color_distiller(palette='Spectral') +
    theme_classic()
  
  if (reverse_weighted_stcon) {
    out_p <- out_p +
      geom_sf(aes(color=inverse_stcon)) 
  } else {
    out_p <- out_p + geom_sf(aes(color=stcon_value)) 
  }
  
  return(out_p)
}


#------ compute_Fdist ----------------------------------------------------
#' Calculate distance to nearest active site
#'
#' This function calculates the distance to the nearest active site based on 
#' spatiotemporal graph networks. 
#'
#' The initial use for this function is to quantify the distance to the nearest
#' perennial site, but it does not necessarily need to quantify water presence.
#' Each cell value must represent a feature defining connectivity "on" or "off" 
#' and that can be transmitted to build meaningful links in a spatiotemporal 
#' graph.
#'
#' Requires: assertthat, data.table, igraph, magrittr, purrr
#'
#' @param sites_status_matrix A matrix representing the status of the site (wet/dry, 
#' active/inactive). 
#' The dataset should have columns for each site and rows for each monitored day.
#' Warning: no other columns should be included (e.g., date, other IDs)
#' @param network_structure A square matrix representing the basic connections among 
#' sites (adjacency matrix): for a given site (row), each adjacent connected 
#' site (column) is given a value of 1, all others 0. Must have the same number
#' of rows and columns as there are columns in sites_status_matrix.
#' @param routing_mode The direction for graph connectivity when directed, 
#' can be "in" (routing from upstream if directed), "out" (routing from 
#' downstream if directed), or "all". See ?igraph or ?igraph::closeness
#' for a better understanding. 
#' @param dist_matrix A distance matrix representing the distances between sites.
#' Can be any type of distance (euclidean, environmental, topographic, ...) 
#' between pairs of sites
#' 

# in_country <- 'Spain'
# in_preformatted_data = tar_read(preformatted_data_STcon)[[in_country]]
# in_net_shp_path <- tar_read(network_ssnready_shp_list)[[in_country]]
# sites_status_matrix = in_preformatted_data$sites_status_matrix
# network_structure = in_preformatted_data$network_structure
# raw_dist_matrix <- in_preformatted_data$river_dist_mat
# routing_mode = 'in'

compute_Fdist <- function(sites_status_matrix, network_structure, 
                                routing_mode, raw_dist_matrix, in_net_shp_path) {
  
  #Get a matrix of adjacency relationship among segments that includes routing mode
  #(1 if connected, Inf if not connected -- depending on routing mode)
  routed_adjacency <- graph_from_adjacency_matrix(network_structure, 
                                   mode = 'directed', 
                                   diag = FALSE) %>%
    distances(mode = routing_mode, 
              algorithm = "unweighted")
  routed_adjacency[!is.infinite(routed_adjacency)] <- 1
  
  #To turn on or off the "routing mode" parameter
  #When routing mode is "in", Inf values for sites downstream (inverse when mode is "out")
  #value of 0 for site to itself, and actual distance value for sites upstream (downstream if mode is out)
  directed_dist_matrix <- routed_adjacency*raw_dist_matrix
  
  #For each time step and reach, compute nearest wet site
  #Convert self-distance to Inf to avoid taking it in account with column minimums
  sites_status_matrix[sites_status_matrix==0] <- Inf 
  
  dist_to_nearest_wet <- sites_status_matrix[,{
    pair_dist_status <- directed_dist_matrix*as.numeric(as.matrix(.SD)) #Multiple distance by site status (Inf is dry, 1 is wet)
    pair_dist_status[is.na(pair_dist_status)] <- Inf
    as.list(Rfast::colMins(pair_dist_status, value=T))
  }
  ,  by=date] %>% 
    setnames(names(sites_status_matrix))
  
  #Merge to UID for subsequent use
  dist_to_nearest_wet_melt <- melt(dist_to_nearest_wet,
                                   id.vars = 'date', 
                                   variable.name = 'ID',
                                   variable.factor = FALSE) %>%
    .[, ID := as.integer(ID)] %>%
    setnames('value', 'Fdist')
  
  #Get original network data
  net_dt <- in_net_shp_path %>%
    as.data.table %>%
    setorder(from)
  
  #Identify the outlet (NA in to_cat_shp)
  outlet_from <- net_dt[is.na(net_dt$to_cat_shp),]$to 
  
  IDs_sel <- unique(dist_to_nearest_wet_melt$ID) %>%
    setdiff(outlet_from) 
  UIDs_to_assign <- net_dt[match(IDs_sel, net_dt$from), list(ID=IDs_sel, UID)]
  
  dist_to_nearest_wet_UID <- merge(dist_to_nearest_wet_melt, 
                                   UIDs_to_assign, 
                                   by='ID')

  return(dist_to_nearest_wet_UID)
}

#------ compute_Fdist_rolling ----------------------------------------------------
# in_country <- 'Czech'
# in_Fdist_dt <- tar_read(dist_to_wet_directed)[[in_country]]
# in_sites_dt <- as.data.table(vect(tar_read(site_snapped_gpkg_list)[[in_country]]))
# setnames(setDT(in_Fdist_dt), 'value', 'Fdist')
# 
# check <- compute_Fdist_rolling(in_Fdist_dt, in_sites_dt)

compute_Fdist_rolling <- function(in_Fdist_dt, in_sites_dt) {
  rollingstep_short <-  c(10, 30, 60, 90, 120, 180)
  rollingstep_long <- c(365, 365*5, 365*10)
  rollingstep <- c(rollingstep_short, rollingstep_long)
  
  #Compute mean Fdist within rolling window, only for sites
  Fdist_sites_rolling <- in_Fdist_dt[UID %in% unique(in_sites_dt$UID),] %>% 
    .[,
      paste0("Fdist_mean_", rollingstep, "past") :=
        frollapply(Fdist, n=rollingstep, 
                   FUN=mean, 
                   align='right'), by = .(UID)
    ]
  
  #Compute max Fdist within rolling window, only for sites
  Fdist_sites_rolling[,
                       paste0("Fdist_max_", rollingstep, "past") :=
                         frollapply(Fdist, n=rollingstep, 
                                    FUN=max, 
                                    align='right'), by = .(UID)
  ]
  
  return(Fdist_sites_rolling)
}


#------ compile_hydrocon_country -----------------------------------------------
# in_hydrostats_sub_comb <- tar_read(hydrostats_sub_comb)
# in_STcon_directed <- tar_read(STcon_directed_formatted)
# in_STcon_undirected <- tar_read(STcon_undirected_formatted)
# in_Fdist_directed <- tar_read(Fdist_network_directed)
# in_Fdist_undirected <- tar_read(Fdist_network_undirected)
# in_site_snapped_gpkg_list <- tar_read(site_snapped_gpkg_list)
# in_country <- 'Spain'

compile_hydrocon_country <- function(in_hydrostats_sub_comb, 
                                     in_STcon_directed,
                                     in_STcon_undirected, 
                                     in_Fdist_directed,
                                     in_Fdist_undirected,
                                     in_site_snapped_gpkg_list,
                                     in_country) {
  
  in_sites_dt <- as.data.table(vect(in_site_snapped_gpkg_list[[in_country]]))
  
  hydrostats_isflowing_site <- in_hydrostats_sub_comb[[
    paste0('hydrostats_sub_', in_country, '_isflowing')]][['site']] %>%
    setDT
  
  hydrostats_isflowing_drn <- in_hydrostats_sub_comb[[
    paste0('hydrostats_sub_', in_country, '_isflowing')]][['drn']] %>%
    setDT
  
  hydrostats_qsim <- in_hydrostats_sub_comb[[
    paste0('hydrostats_sub_', in_country, '_qsim')]] %>%
    setDT
  
  cols_to_keep_site <- c(setdiff(names(hydrostats_isflowing_site),
                                 names(hydrostats_qsim)),
                         c('date', 'site'))
  
  cols_to_keep_drn <- c(setdiff(names(hydrostats_isflowing_drn),
                                names(hydrostats_qsim)),
                        c('date'))
  
  #Format STcon data
  in_STcon_directed[[in_country]]$STcon_dt[
    , variable := paste0(variable, '_directed')]
  in_STcon_undirected[[in_country]]$STcon_dt[
    , variable := paste0(variable, '_undirected')]
  
  STcon_cast <- rbind(in_STcon_directed[[in_country]]$STcon_dt,
                      in_STcon_undirected[[in_country]]$STcon_dt) %>%
    dcast(date+UID~variable, value.var = 'stcon_value')
  #With nsim
  # STcon_cast <- in_STcon_directed[[in_country]]$STcon_dt %>%
  #   dcast(date+UID+nsim~variable, value.var = 'stcon_value') %>%
  #   .[, nsim := as.integer(nsim)]
  
  #Format Fdist data
  in_Fdist_directed[[in_country]][, variable := paste0(variable, '_directed')]
  in_Fdist_undirected[[in_country]][, variable := paste0(variable, '_undirected')]
  Fdist_cast <- rbind(in_Fdist_directed[[in_country]],
                      in_Fdist_undirected[[in_country]]) %>%
    dcast(date+UID~variable, value.var = 'fdist_value')

  #Merge all statistics
  out_dt <- merge(hydrostats_isflowing_site[, cols_to_keep_site, with=F],
                  hydrostats_isflowing_drn[, cols_to_keep_drn, with=F],
                  by='date') %>%
    merge(hydrostats_qsim, by=c('date', 'site')) %>%
    merge(in_sites_dt[country_sub == in_country,
                      .(site, UID)], ., by='site') %>%
    merge(STcon_cast, by=c('date', 'UID'), all.x=T) %>%
    merge(Fdist_cast, by=c('date', 'UID'), all.x=T)
    
  # merge(STcon_cast, by=c('date', 'UID', 'nsim'), all.x=T) #with nsim
  
  return(out_dt)  
}

#------ merge_allvars_sites ----------------------------------------------------
# in_country <- 'Spain'
# in_spdiv_local <- tar_read(spdiv_local)
# in_spdiv_drn <- tar_read(spdiv_drn)
# in_hydrocon_compiled <- tar_read(hydrocon_compiled)
# in_ssn_eu <- tar_read(ssn_eu)
# in_env_dt <- tar_read(env_dt)
# in_genal_upa = tar_read(genal_sites_upa_dt)

merge_allvars_sites <- function(in_spdiv_local, in_spdiv_drn,
                                in_hydrocon_compiled, in_env_dt,
                                in_genal_upa) {
  
  #Fill basin area NAs in environmental data for Genal basin in Spain
  #https://stackoverflow.com/questions/72940045/replace-na-in-a-table-with-values-in-column-with-another-table-by-conditions-in
  in_env_dt[is.na(basin_area_km2),
                    basin_area_km2 :=  in_genal_upa[
                      .SD, on='site', x.basin_area_km2]]
  
  #Compute state of flow in previous time step (could be inserted much before in the workflow later)
  in_env_dt[, state_of_flow_tm1 := lag(state_of_flow, n=1L), by=site]
  #Assume that the state of flow prior to first time step is the same as during the first time step
  in_env_dt[is.na(state_of_flow_tm1), state_of_flow_tm1 := state_of_flow]
  
  #Merge diversity data
  setDT(in_spdiv_drn)
  spdiv <- merge(in_spdiv_local, in_spdiv_drn,
                 by=c('country', 'organism')) %>%
    .[!(site=='GEN04' & campaign=='1'),]
  #Remove GEN04_1 from all organisms, no local environmental data, sampled only for eDNA. Too unsure.

  #Create "organism_class" column for labeling/coloring/merging
  spdiv[, organism_class := gsub('_[a-z]+', '', organism)]
  
  
  #List column names by originating dt
  dtcols <- list(
    div = setdiff(names(spdiv), 
                  c(names(in_hydrocon_compiled), names(in_env_dt))),
    hydro_con = setdiff(names(in_hydrocon_compiled), 
                        c(names(spdiv), names(in_env_dt))),
    env = setdiff(names(in_env_dt),
                  c(names(spdiv), names(in_hydrocon_compiled))),
    group_cols = c("running_id", "site", "date", "campaign", "organism", 
                   "organism_class", "country", "UID", "stream_type", 
                   "state_of_flow", "state_of_flow_tm1"),
    exclude_cols = c("ncampaigns", "name", "isflowing", "reach_length",
                     "noflow_period", "noflow_period_dur", "last_noflowdate", "drn",
                     "if_ip_number_and_size_2_axes_+_depth_of_the_pools",
                     "latitude", "longitude", "reach_id", "hy", "month",
                     'min_wetted_width', 'left_river_bank_slope', 'right_river_bank_slope',
                     'qsim')
  )

  #Compute average metric between sediment and biofilm for eDNA
  avg_spdiv_edna <- spdiv[
    organism_class %in% c('dia', 'fun', 'bac'), 
    lapply(.SD, function(x) mean(x, na.rm=T)), 
    by = eval(names(spdiv)[names(spdiv) %in% 
                             setdiff(dtcols$group_cols, 'organism')]),
    .SDcols = setdiff(dtcols$div, 
                      c(dtcols$group_cols, dtcols$exclude_cols))
  ] %>%
    .[, organism := organism_class]
  
  spdiv <- rbind(spdiv, avg_spdiv_edna, use.names=T, fill=T)
  
  #Merge diversity metrics with hydro_con
  spdiv_hydro_con <- merge(spdiv, in_hydrocon_compiled,
                           by=c('date', 'site'), all.x=T) 
  
  #Merge environmental variables
  setDT(in_env_dt)
  env_subcols <-  c('site', 'campaign', 
                    setdiff(names(in_env_dt), names(spdiv_hydro_con)))
  all_vars_merged <- merge(spdiv_hydro_con, in_env_dt[, env_subcols, with=F],
                           by=c('site', 'campaign'), all.x=T)
  
  return(list(
    dt = all_vars_merged,
    cols = dtcols)
  )
}


#------ plot_edna_biof_vs_sedi -------------------------------------------------
#in_allvars_merged <- tar_read(allvars_merged)$dt

plot_edna_biof_vs_sedi <- function(in_allvars_merged) {
  allvars_edna <- setDT(in_allvars_merged$dt)[organism_class %in% c('dia', 'fun', 'bac'),] 
  allvars_edna[, edna_source := gsub('^[a-z]+_', '', organism)]
  
  #Compare richness by source of edna (biofilm vs sediment), organism and country
  rich_cast <- dcast(allvars_edna, 
                    country+site+date+organism_class~edna_source, 
                    value.var = 'richness')
  p_richness <- ggplot(rich_cast, aes(x=biof, y=sedi)) + 
    geom_point() + 
    geom_abline() + 
    coord_fixed() +
    facet_grid(organism_class~country)
  
  #Compare gamma by source of edna (biofilm vs sediment), organism and country
  site_gamma_cast <- dcast(allvars_edna, 
                    country+site+date+organism_class~edna_source, 
                    value.var = 'Gamma')
  p_site_gamma <- ggplot(site_gamma_cast, aes(x=biof, y=sedi)) + 
    geom_point() + 
    geom_abline() + 
    coord_fixed() +
    facet_grid(organism_class~country)
  
  return(list(
    richness=p_richness,
    site_gamma=p_site_gamma
  ))
  
  
}


#------ ordinate_local_env ----------------------------------------------------
# autoplot(local_env_pca$miv_nopools$pca, data = local_env_pca$miv_nopools$trans_dt, 
#          loadings = TRUE, loadings.colour = 'blue',
#          loadings.label = TRUE, loadings.label.size = 5) + theme_bw()
# in_allvars_merged <- tar_read(allvars_merged)
ordinate_local_env <- function(in_allvars_merged) {
  allvars_merged_copy <- copy(in_allvars_merged)
  
  #1. Compute PCA for miv_nopools ----------------------------------------------
  env_cols_miv <- c('avg_velocity_macroinvertebrates', 'embeddedness',
                    'bedrock', 'particle_size', 'oxygen_sat', 'filamentous_algae',
                    'incrusted_algae', 'macrophyte_cover', 'leaf_litter_cover',
                    'moss_cover', 'wood_cover', 'riparian_cover_in_the_riparian_area',
                    'shade', 'hydromorphological_alteration', 'm2_biofilm',
                    'conductivity_micros_cm', 'ph', 'temperature_c')
  #'oxygen_mg_l'
  
  #Convert all columns to numeric (rather than integer)
  allvars_merged_copy$dt[, (env_cols_miv) := lapply(.SD, as.numeric), 
                       .SDcols = env_cols_miv] 
  
  #Fill NAs hierarchically. First by site, then by country, then overall
  miv_nopools_dt <- fill_nas_hierarchical(
    dt = allvars_merged_copy$dt[organism == 'miv_nopools'], 
    cols_to_fill = env_cols_miv, 
    site_col = 'site', 
    country_col = 'country')
  
  #Check distributions by country
  dt_miv_envmelt <- melt(allvars_merged_copy$dt[organism == 'miv',], 
                         id.vars = c('country', 'site', 'campaign'),
                         measure.vars = env_cols_miv)
  # ggplot(dt_miv_envmelt, #[variable=='conductivity_micros_cm',],
  #        aes(x=country, y=value, color=country)) +
  #   geom_jitter() + 
  #   facet_wrap(~variable, scales='free_y') +
  #   scale_y_sqrt()
  # 
  # ggplot(dt_miv_envmelt, aes(x=value)) +
  #   geom_density() + 
  #   facet_wrap(~variable, scales='free')
  
  out_list_miv <- trans_pca_wrapper(in_dt = miv_nopools_dt, 
                                    in_cols_to_ordinate = env_cols_miv, 
                                    id_cols = c('site', 'date', 'country'), 
                                    group_cols = NULL, 
                                    num_pca_axes = 4)
  
  # out_list_miv_country <- trans_pca_wrapper(in_dt = miv_nopools_dt, 
  #                                           in_cols_to_ordinate = env_cols_miv, 
  #                                           id_cols = c('site', 'date'), 
  #                                           group_cols = 'country', 
  #                                           num_pca_axes = 4)
  
  #2. Compute PCA for eDNA data ------------------------------------------------
  edna_orglist <- c('dia_sedi', 'dia_biof', 'dia', 
                    'fun_sedi', 'fun_biof',  'fun',
                    'bac_sedi_nopools', 'bac_biof_nopools', 'bac')
  
  env_cols_edna <- c('filamentous_algae', 'incrusted_algae', 
                     'macrophyte_cover', 'leaf_litter_cover','moss_cover', 
                     'wood_cover', 'riparian_cover_in_the_riparian_area',
                     'shade', 'm2_biofilm', 'conductivity_micros_cm',
                     'oxygen_sat', 'ph','temperature_c')
  
  #Convert all columns to numeric (rather than integer)
  allvars_merged_copy$dt[, (env_cols_edna) := lapply(.SD, as.numeric), 
                       .SDcols = env_cols_edna]
  
  #Check distributions by country
  dt_edna_envmelt <- unique(allvars_merged_copy$dt[organism %in% edna_orglist,],
                            by=c('site', 'date')) %>%
    melt(id.vars = c('country', 'site', 'campaign', 'organism'),
         measure.vars = env_cols_edna)
  # ggplot(dt_edna_envmelt, #[variable=='conductivity_micros_cm',],
  #        aes(x=country, y=value, color=country)) +
  #   geom_jitter() + 
  #   facet_grid(organism~variable, scales='free_y') +
  #   scale_y_log10()
  # 
  # ggplot(dt_edna_envmelt, aes(x=value)) +
  #   geom_density() + 
  #   facet_wrap(~variable, scales='free')
  
  #Fill NAs hierarchically. First by site, then by country, then overall
  #this gives the average conditions when there is flow to dry samples. 
  #will need to think about how to use this
  edna_dt <- unique(allvars_merged_copy$dt[organism %in% edna_orglist,],
                    by=c('site', 'date')) %>%
    .[, c('country', 'site', 'date', env_cols_edna), with=F] %>%
    fill_nas_hierarchical(cols_to_fill = env_cols_edna, 
                          site_col = 'site', 
                          country_col = 'country') 
  
  out_list_edna <- trans_pca_wrapper(in_dt = edna_dt, 
                                     in_cols_to_ordinate = env_cols_edna, 
                                     id_cols = c('site', 'date', 'country'), 
                                     group_cols = NULL, 
                                     num_pca_axes = 4)
  
  #Merge dts
  out_dt_miv <- merge(
    allvars_merged_copy$dt[organism == 'miv_nopools',
                         .(site, date, country, organism)],
    out_list_miv$dt)
  
  out_dt_edna <- merge(
    allvars_merged_copy$dt[organism %in% edna_orglist,
                         .(site, date, country, organism)],
    out_list_edna$dt)
  

  out_dt <- rbind(out_dt_miv, out_dt_edna) %>%
    .[, organism_class := gsub('_[a-z]+', '', organism)] %>%
    data.table::unique(by=c('country', 'site', 'date', 'organism_class'))

  return(list(
    miv_nopools = out_list_miv,
    edna = out_list_edna,
    dt_all = out_dt
  ))
}
#------ create_ssn_europe ------------------------------------------------------
# in_network_path = tar_read(network_ssnready_shp_list)
# in_sites_path = tar_read(site_snapped_gpkg_list)
# in_barriers_path = tar_read(barrier_snapped_gpkg_list)
# in_allvars_merged = tar_read(allvars_merged)
# in_local_env_pca = ordinate_local_env(in_allvars_merged)
# in_hydromod = tar_read(hydromod_comb)
# out_dir = 'results/ssn'
# out_ssn_name = 'all_drns'
# overwrite=T

create_ssn_europe <- function(in_network_path,
                              in_sites_path,
                              in_allvars_merged,
                              in_local_env_pca,
                              in_barriers_path,
                              in_hydromod,
                              out_dir,
                              out_ssn_name,
                              overwrite = T) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  lsn_path <- file.path(out_dir,
                        paste0(out_ssn_name, '_lsn')
  )
  
  #Build landscape network (lsn) -----------------------------------------------
  #Read input network
  net_eu <- lapply(names(in_network_path), function(in_country) {
    #print(in_country)
    net_proj <- in_network_path[[in_country]] %>%
      st_cast("LINESTRING") %>%
      #Make sure that the geometry column is equally named regardless 
      #of file format (see https://github.com/r-spatial/sf/issues/719)
      st_set_geometry('geometry') %>%
      st_transform(3035)
    
    hydromod_country <- in_hydromod[[
      paste0('hydromod_dt_', in_country, '_qsim')]]
    
    net_hydro <- merge(net_proj,
                       hydromod_country$data_all[
                         date < as.Date('2021-10-01', '%Y-%m-%d'),  #Link q data
                         list(mean_qsim = mean(qsim, na.rm=T)), 
                         by=reach_id],
                       by.x = 'cat', by.y = 'reach_id')
    return(net_hydro)
  }) %>% do.call(rbind, .)
  
  edges_lsn <- SSNbler::lines_to_lsn(
    streams = net_eu,
    lsn_path = lsn_path,
    check_topology = TRUE,
    snap_tolerance = 0.1,
    topo_tolerance = 20,
    overwrite = overwrite
  )
  
  #Incorporate sites into the landscape network --------------------------------
  sites_eu <- lapply(names(in_sites_path), function(in_country) {
    st_read(in_sites_path[[in_country]]) %>%
      st_transform(3035)  
  }) %>% 
    do.call(rbind, .) %>%
    rename(country=country_sub)
  
  sites_lsn <- SSNbler::sites_to_lsn(
    sites = sites_eu,
    edges =  edges_lsn,
    lsn_path = lsn_path,
    file_name = "sites",
    snap_tolerance = 5,
    save_local = TRUE,
    overwrite = TRUE
  )
  
  setDT(in_allvars_merged$dt)[country == 'Czech Republic', country := 'Czech']
  sites_lsn_attri <- merge(sites_lsn,
                           in_allvars_merged$dt, 
                          by=c('country', 'site')) %>%
    merge(in_local_env_pca$dt_all,
          by=c('site', 'date', 'country', 'organism_class'))
  
  sites_list <- list(sites = sites_lsn_attri)
  
  #Incorporate barriers into the landscape network
  #Only keep barriers over 2 m and under 100 m snap from network
  barriers_eu_sub <- lapply(names(in_barriers_path), function(in_country) {
    st_read(in_barriers_path[[in_country]]) %>%
      st_transform(3035) %>%
      filter((!is.na(Height) & Height > 2) & snap_dist_m < 100)
  }) %>% do.call(rbind, .)
  
  if (nrow(barriers_eu_sub) > 0) {
    barriers_lsn <- sites_to_lsn(
      sites = barriers_eu_sub,
      edges =  edges_lsn,
      lsn_path = lsn_path,
      file_name = "barriers",
      snap_tolerance = 5,
      save_local = TRUE,
      overwrite = TRUE
    )
    
    sites_list$barriers <- barriers_lsn
  }
  
  #Calculate upstream distance
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
  
  #Compute segment Proportional Influence (PI) and Additive Function Values (AFVs)
  if (min(net_eu$mean_qsim) > 0) {
    edges_lsn$mean_qsim_sqrt <- sqrt(edges_lsn$mean_qsim)
    
    edges_lsn <- afv_edges(
      edges = edges_lsn,
      infl_col = "mean_qsim_sqrt",
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
    
  } else {
    stop("Trying to use mean discharge to compute Additive Function Values (AFVs),
         but there are 0s in the discharge column.")
  }
  
  
  #Assemble an SSN for each organism 
  #(so that it only includes data for the  corresponding sites and dates)
  out_ssn_list <- lapply(unique(sites_list_lsn$sites$organism), function(in_org) {
    out_ssn_path <- paste0(out_dir, '_', out_ssn_name, '_', in_org, '.ssn')
    
    out_ssn <- ssn_assemble(
      edges = edges_lsn,
      lsn_path = lsn_path,
      obs_sites = sites_list_lsn$sites[sites_list_lsn$sites$organism == in_org,],
      ssn_path = out_ssn_path,
      import = TRUE,
      check = TRUE,
      afv_col = "afv_qsqrt",
      overwrite = TRUE
    )
    
    return(list(
      path = out_ssn_path,
      ssn = out_ssn
    ))
  }) %>% setNames(unique(sites_list_lsn$sites$organism))
  
  return(out_ssn_list)
}

#------ compute_cor_matrix -----------------------------------------------------
#in_allvars_merged <- tar_read(allvars_merged)
compute_cor_matrix <- function(in_allvars_merged) {
  dt <- in_allvars_merged$dt
  cols_by_origin <- in_allvars_merged$cols
  exclude_cols <- cols_by_origin$exclude_cols
  group_cols <- cols_by_origin$group_cols
  
  #Pre-formatting
  dt[is.na(noflow_period_dur), noflow_period_dur := 0]
  dt[, avg_bank_slope := rowMeans(.SD),
     .SDcols=c('left_river_bank_slope', 'right_river_bank_slope')]
  dt[, PrdD := as.numeric(PrdD)]
  
  #Define columns to group by and columns to correlate
  div_cols <- setdiff(cols_by_origin$div, c(exclude_cols, group_cols))
  hydrocon_cols <- setdiff(cols_by_origin$hydro_con, c(exclude_cols, group_cols))
  env_cols <- c(setdiff(cols_by_origin$env, c(exclude_cols, group_cols)),
                'avg_bank_slope')
  
  # --- Calculate Correlations ---
  # 1. Overall Correlation (Hydro + Env)
  cor_hydroenv <- compute_cor_matrix_inner(
    dt,
    x_cols = c(hydrocon_cols, env_cols),
    exclude_diagonal = FALSE) 
  
  # 2. By Organism (Div x (Hydro + Env))
  cor_div <- compute_cor_matrix_inner(
    dt,
    group_vars = "organism",
    x_cols = div_cols,
    y_cols = c(hydrocon_cols, env_cols),
    exclude_diagonal = FALSE) 
  
  # 3. By Country (Hydro + Env)
  cor_hydroenv_bydrn <- compute_cor_matrix_inner(
    dt,
    group_vars = "country",
    x_cols = c(hydrocon_cols, env_cols),
    exclude_diagonal = FALSE)  
  
  # 4. By Organism and Country (Div x (Hydro + Env))
  cor_div_bydrn <- compute_cor_matrix_inner(
    dt,
    group_vars = c("organism", "country"),
    x_cols = div_cols,
    y_cols = c(hydrocon_cols, env_cols),
    exclude_diagonal = FALSE) 
  
  return(list(hydroenv = cor_hydroenv,
              div = cor_div,
              hydroenv_bydrn = cor_hydroenv_bydrn,
              div_bydrn = cor_div_bydrn,
              div_cols = div_cols,
              hydrocon_cols = hydrocon_cols,
              env_col = env_cols
  ))
}

#------ plot_cor_heatmaps -------------------------------------------------------
in_cor_matrices <- tar_read(cor_matrices_list)
in_allvars_merged <- tar_read(allvars_merged)
p_threshold <- 0.05

create_correlation_heatmap <- function(cor_matrix, p_matrix, title,
                                       p_threshold = 1,
                                       is_square = TRUE, hc_order = TRUE) { # Add is_square argument
  
  # Create the heatmap
  heatmap <- ggcorrplot(cor_matrix,
                        p.mat = p_matrix,
                        hc.order = (is_square && hc_order),
                        hc.method = 'average',
                        lab = TRUE, lab_size = 3,
                        digits = 1, insig = 'blank', sig.level = p_threshold,
                        outline.color = "white",
                        type = if (is_square) "full" else "upper",  # Control square/triangle
                        show.diag = if (is_square) TRUE else FALSE) + # Control diagonal
    scale_fill_distiller(
      name = str_wrap("Correlation coefficient Spearman's rho", 20),
      palette = 'RdBu',
      breaks = c(-0.7, -0.5, 0, 0.5, 0.7),
      limits = c(-1, 1)) +
    ggtitle(title)
  
  return(heatmap)
}

plot_cor_heatmaps <- function(in_cor_matrices, in_allvars_merged,
                              p_threshold = 1) {
  org_list <- unique(in_cor_matrices$div$organism) %>%
    setdiff('miv')
  country_list <- unique(in_cor_matrices$hydroenv_bydrn$country)
  
  # # --- 1. Hydroenv Correlation ---
  # hydroenv_sub <- in_cor_matrices$hydroenv[p_value <= p_threshold & 
  #                                            correlation >= cor_threshold,]
  # 
  # Extract the correlation and p-value matrices
  cor_matrix_hydroenv <- dcast(in_cor_matrices$hydroenv, 
                               variable1 ~ variable2, 
                               value.var = "correlation")
  rnames <- cor_matrix_hydroenv$variable1
  cor_matrix_hydroenv <- as.matrix(cor_matrix_hydroenv[, -1])
  rownames(cor_matrix_hydroenv) <- rnames
  cor_matrix_hydroenv[is.na(cor_matrix_hydroenv)] <- 0
  
  p_matrix_hydroenv <- dcast(in_cor_matrices$hydroenv, 
                             variable1 ~ variable2, 
                             value.var = "p_value")
  p_matrix_hydroenv <- as.matrix(p_matrix_hydroenv[, -1])
  rownames(p_matrix_hydroenv) <- rnames
  p_matrix_hydroenv[is.na(p_matrix_hydroenv)] <- 1L
  
  # Create and print the heatmap
  hydroenv_heatmap <- create_correlation_heatmap(
    cor_matrix = cor_matrix_hydroenv, 
    p_matrix = p_matrix_hydroenv,
    title = "Hydro-Env Correlations",
    p_threshold = p_threshold,
    is_square = TRUE)
  #ggplotly(hydroenv_heatmap)
  
  # --- 2. By Organism (Div x (Hydro + Env)) ---
  div_heatmaps <- lapply(org_list, function(org) {
    div_sub <- in_cor_matrices$div[organism == org,] %>%
      .[variable1 %in% in_allvars_merged$cols$div,]
    
    # Convert to matrices
    cor_matrix_org <- dcast(div_sub, 
                            variable1 ~ variable2, 
                            value.var = "correlation")
    rnames_org <- cor_matrix_org$variable1
    cor_matrix_org <- as.matrix(cor_matrix_org[, -1])
    rownames(cor_matrix_org) <- rnames_org
    cor_matrix_org[is.na(cor_matrix_org)] <- 0L
    
    p_matrix_org <- dcast(div_sub,
                          variable1 ~ variable2, 
                          value.var = "p_value")
    p_matrix_org <- as.matrix(p_matrix_org[, -1])
    rownames(p_matrix_org) <- rnames_org
    p_matrix_org[is.na(p_matrix_org)] <- 1L
    
    heatmap <- create_correlation_heatmap(
      cor_matrix = cor_matrix_org, 
      p_matrix = p_matrix_org,
      title = paste("Div vs. Hydro/Env - Organism:", org),
      p_threshold = p_threshold,
      is_square = FALSE)
  }) %>% setNames(org_list)
  
  # --- 3. Hydroenv Correlation by Country (Square) ---
  hydroenv_country_heatmaps <- lapply(country_list, function(in_country) {
    hydroenv_sub <- in_cor_matrices$hydroenv_bydrn[country == in_country]
    
    cor_matrix_hydroenv <- dcast(hydroenv_sub, variable1 ~ variable2, 
                                 value.var = "correlation")
    rnames <- cor_matrix_hydroenv$variable1
    cor_matrix_hydroenv <- as.matrix(cor_matrix_hydroenv[, -1])
    rownames(cor_matrix_hydroenv) <- rnames
    cor_matrix_hydroenv[is.na(cor_matrix_hydroenv)] <- 0
    
    p_matrix_hydroenv <- dcast(hydroenv_sub, variable1 ~ variable2, 
                               value.var = "p_value")
    p_matrix_hydroenv <- as.matrix(p_matrix_hydroenv[, -1])
    rownames(p_matrix_hydroenv) <- rnames
    p_matrix_hydroenv[is.na(p_matrix_hydroenv)] <- 1
    
    hydroenv_heatmap <- create_correlation_heatmap(
      cor_matrix = cor_matrix_hydroenv,
      p_matrix = p_matrix_hydroenv,
      title = paste("Hydro-Env Correlations - Country:", in_country),
      p_threshold = p_threshold,
      is_square = TRUE)
  }) %>% setNames(country_list)
  
  # --- 4. Div x Hydroenv Correlation by Country (Non-square) ---
  div_country_heatmaps <- lapply(country_list, function(in_country) {
    lapply(org_list, function(org) {
      div_hydro_sub <- in_cor_matrices$div_bydrn[
        country == in_country & organism == org,]
      
      cor_matrix_div_hydro <- dcast(div_hydro_sub, variable1 ~ variable2, 
                                    value.var = "correlation")
      rnames <- cor_matrix_div_hydro$variable1
      cor_matrix_div_hydro <- as.matrix(cor_matrix_div_hydro[, -1])
      rownames(cor_matrix_div_hydro) <- rnames
      cor_matrix_div_hydro[is.na(cor_matrix_div_hydro)] <- 0L
      
      p_matrix_div_hydro <- dcast(div_hydro_sub, variable1 ~ variable2, 
                                  value.var = "p_value")
      p_matrix_div_hydro <- as.matrix(p_matrix_div_hydro[, -1])
      rownames(p_matrix_div_hydro) <- rnames
      p_matrix_div_hydro [is.na(p_matrix_div_hydro )] <- 1L
      
      heatmap_div_hydro <- create_correlation_heatmap(
        cor_matrix = cor_matrix_div_hydro,
        p_matrix = p_matrix_div_hydro,
        title = paste("Div x Hydro/Env - Country:", in_country,
                      "- Organism:", org),
        p_threshold = p_threshold,
        is_square = FALSE)
      
    }) %>% setNames(org_list)
  }) %>% setNames(country_list)
  
  return(list(
    hydroenv = hydroenv_heatmap,
    div = div_heatmaps,
    hydroenv_country = hydroenv_country_heatmaps,
    div_country = div_country_heatmaps
  ))
}

#------ compare_drn_hydro ------------------------------------------------------
# in_sites_dt <- tar_read(sites_dt)
# in_hydrocon_compiled <- tar_read(hydrocon_compiled)

compare_drn_hydro <- function(in_hydrocon_compiled, in_sites_dt) {
  hydrocon_compiled_country <- merge(in_hydrocon_compiled, 
                                     in_sites_dt[, .(site, country)], by='site')
  sampling_date_lims <- hydrocon_compiled_country[, list(
    max_date = max(date, na.rm=T),
    min_date = min(date, na.rm=T)), by=country]
  
  plot_relF90past <- ggplot(
    hydrocon_compiled_country[!duplicated(paste(country ,date)),]) +
    geom_area(aes(x=date, y=relF90past), fill='#2b8cbe') +
    geom_rect(data=sampling_date_lims, 
              aes(ymin=0, ymax=1, 
                  xmin=max_date, xmax=max(sampling_date_lims$max_date)),
              fill = 'grey') +
    geom_rect(data=sampling_date_lims, 
              aes(ymin=0, ymax=1, 
                  xmin=min(sampling_date_lims$min_date), xmax=min_date),
              fill = 'grey') +
    coord_cartesian(expand = FALSE) +
    facet_grid(~country) +
    theme_classic() +
    theme(panel.background = element_rect(fill='#feb24c'))
  
  # ggplot(hydrocon_compiled_country,
  #        aes(x=date, y=DurD90past, color = site)) +
  #   geom_line() +
  #   geom_smooth(color='black', method='gam') +
  #   coord_cartesian(expand = FALSE) +
  #   facet_grid(~country) +
  #   theme_classic() +
  #   theme(legend.position = 'none')
  # 
  # ggplot(hydrocon_compiled_country,
  #        aes(x=date, y=PDurD90past, color = site)) +
  #   geom_line() +
  #   geom_smooth(color='black', method='gam') +
  #   coord_cartesian(expand = FALSE) +
  #   facet_grid(~country) +
  #   theme_classic() +
  #   theme(legend.position = 'none')
  
  return(plot_relF90past)
}

#------ check_diversity_dependence_on_habitat_volume ---------------------------
#n_allvars_merged <- tar_read(allvars_merged)

check_cor_div_habvol <- function(in_allvars_merged) {
  #names(in_allvars_merged$dt)
  in_allvars_merged$dt[, avg_vol_miv := avg_depth_macroinvertebrates*average_wetted_width_m]
  in_allvars_merged$dt[, mean_avg_vol_miv := mean(avg_vol_miv), by=site]
  in_allvars_merged$dt[, mean_velo_miv := mean(avg_velocity_macroinvertebrates), by=site]
  
  p_miv_vol <- ggplot(in_allvars_merged$dt[organism=='miv',], 
                      aes(x=avg_vol_miv/mean_avg_vol_miv, 
                          y=richness/mean_richness)) +
    geom_point() +
    geom_smooth(aes(color=site), method='lm', se=F) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  p_miv_vel <- ggplot(in_allvars_merged$dt[organism=='miv_nopools' & 
                                             avg_velocity_macroinvertebrates<5,], 
                      aes(x=avg_velocity_macroinvertebrates, 
                          y=richness/mean_richness)) +
    geom_point() +
    geom_smooth(aes(color=site), method='lm', se=F) +
    geom_smooth(color='black', method='lm', linewidth=1.2) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  p_dia_biof_vol <- ggplot(in_allvars_merged$dt[organism=='dia_biof',], 
                           aes(x=volume_biofilm, 
                               y=richness/mean_richness)) +
    geom_point() +
    geom_smooth(aes(color=site), method='lm', se=F) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  p_dia_biof_m2 <- ggplot(in_allvars_merged$dt[organism=='dia_biof',], 
                          aes(x=m2_biofilm, 
                              y=richness/mean_richness)) +
    geom_point() +
    geom_smooth(aes(color=site), method='lm', se=F) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  p_fun_biof_m2 <- ggplot(in_allvars_merged$dt[organism=='fun_biof',], 
                          aes(x=m2_biofilm, 
                              y=richness/mean_richness)) +
    geom_point() +
    geom_smooth(aes(color=site), method='lm', se=F) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  return(list(
    miv_vol= p_miv_vol,
    miv_vel = p_miv_vel,
    dia_biof_vol = p_dia_biof_vol,
    dia_biof_m2 = p_dia_biof_m2,
    fun_biof_m2 = p_fun_biof_m2
  ))
}


#------ plot_cor_hydrowindow  --------------------------------------------------
#For hydrological variables check by time window

# hydrovar_list <- c('DurD', 'PDurD', 'FreD', 'PFreD',
#                'uQ90', 'oQ10', 'maxPQ', 'PmeanQ',
#                'STcon.*_directed', 'STcon.*_undirected')
# var_substr <- 'STcon.*_directed'
# in_cor_dt <- tar_read(cor_matrices_list)$div_bydrn

plot_cor_hydrowindow <-  function(in_cor_dt, temporal_var_substr, response_var_list,
                                  colors_list, save_plot=T, plot_name_suffix="", out_dir) {
  sub_dt <- in_cor_dt[grep(paste0('^', temporal_var_substr), variable2),] %>%
    .[organism != 'miv',] %>%
    .[(variable1 %in% response_var_list),] 
  
  sub_dt[, window_d := str_extract(variable2,
                                   '([0-9]+(?=past))|((?<=m)[0-9]+)')] %>%
    .[, window_d := factor(window_d, levels=sort(unique(as.integer(window_d))))]
  
  sub_dt[,`:=`(mean_cor = mean(`correlation`),
               sd_cor = sd(`correlation`, na.rm=T)
  ), by=c('window_d', 'organism', 'variable1')]
  
  out_p <- ggplot(sub_dt, aes(x=window_d, y=correlation)) +
    geom_hline(yintercept = 0, color='grey') +
    geom_point(aes(y=mean_cor), size=3, color='darkgrey', alpha=1/3) +
    geom_segment(aes(xend=window_d, y=mean_cor-sd_cor, yend=mean_cor+sd_cor),
                 color='darkgrey') +
    geom_line(aes(group=country), color='darkgrey', linewidth=1) +
    geom_line(data=sub_dt[p_value < 0.05,], 
              aes(color=country, group=country), linewidth=1) +
    scale_color_manual(values=colors_list) +
    facet_grid(organism~variable1) + #, scales='free_y') +
    theme_bw() +
    ggtitle(paste0(temporal_var_substr, plot_name_suffix))
  
  if (save_plot) {
    out_path <- file.path(out_dir, paste0('plot_cor_hydrowindow_', temporal_var_substr,
                                          plot_name_suffix, '.png'))
    ggsave(out_path, plot = out_p, 
           height = 15, width = 7.5, units='in', dpi = 300)
  }
  
  return(out_p)
}

#------ plot_scatter_lm --------------------------------------------------
# hydrovar_list <- c('DurD', 'PDurD', 'FreD', 'PFreD', 
#                    'uQ90', 'oQ10', 'maxPQ', 'PmeanQ',
#                    'STcon.*_directed', 'STcon.*_undirected')
# temporal_var_substr <- 'DurD'
# in_allvars_merged <- tar_read(allvars_merged)
# response_var = 'richness' #c('richness', 'invsimpson','Jtm1'),
# colors_list = drn_dt$color
# plot=T
# out_dir = resdir

plot_scatter_lm <-  function(in_allvars_merged, temporal_var_substr, response_var,
                          colors_list, save_plot=T, plot_name_suffix="", out_dir) {
  
  dt <- in_allvars_merged$dt
  x_cols <- grep(paste0('^', temporal_var_substr), names(dt), value=T) %>%
    .[seq(1, length(.), 2)]
  
  dt_melt <- melt(dt, 
                  id.vars=c(in_allvars_merged$cols$group_cols, response_var), 
                  measure.vars=x_cols)

  p_lm <- ggplot(dt_melt, 
                 aes_string(x='value', y=response_var, color='country')) +
    geom_point(alpha=0.5) +
    geom_smooth(method='lm', se=F) +
    scale_color_manual(values=colors_list) +
    scale_y_sqrt() +
    facet_grid(organism~variable, scales='free') +
    theme_bw()

  if (save_plot) {
    out_path <- file.path(out_dir, paste0('plot_lm_', temporal_var_substr, '_',
                                          response_var, plot_name_suffix, '.png'))
    ggsave(out_path, plot = p_lm, 
           height = 15, width = 15, units='in', dpi = 300)
  }
  
  return(p_lm)
}
  

#------ quick_ssn ------
# in_ssn_eu <- tar_read(ssn_eu)
# 
# in_formula = 'richness ~ log10(basin_area_km2) + log10(basin_area_km2):country + DurD60past + DurD60past:country'
# in_ssn <- in_ssn_eu$miv_nopools$ssn
# tar_load(ssn_covtypes)

quick_ssn <- function(in_ssn, in_formula, ssn_covtypes, estmethod = "ml") {
  SSN2::ssn_create_distmat(in_ssn)

  #summary(lm(formula = as.formula(in_formula), data=in_ssn$obs))

  ssn_list <- mapply(function(down_type, up_type, euc_type) {
    print(paste(down_type, up_type, euc_type))
    out_ssn <- ssn_lm(
      formula = as.formula(in_formula),
      ssn.object = in_ssn,
      taildown_type = down_type,
      tailup_type = up_type,
      euclid_type = euc_type,
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = ~ country,
      estmethod = estmethod
    )
    return(out_ssn)
  },
  down_type = ssn_covtypes$down,
  up_type = ssn_covtypes$up,
  euc_type = ssn_covtypes$euc,
  SIMPLIFY = FALSE) %>%
    setNames(ssn_covtypes$label)
  
  return(ssn_list)
}
 
#------ model_ssn_hydrowindow --------------------------------------------------
# hydrovar_list <- c('DurD', 'PDurD', 'FreD', 'PFreD',
#                'uQ90', 'oQ10', 'maxPQ', 'PmeanQ',
#                'STcon.*_directed', 'STcon.*_undirected')
# var_substr <- 'STcon.*_directed'

in_ssn = tar_read(ssn_eu)
organism = 'miv_nopools_flying'
formula_root = '~ log10(basin_area_km2) + log10(basin_area_km2):country'
hydro_var = 'DurD365past'
response_var = 'richness'
tar_load(ssn_covtypes)

model_ssn_hydrowindow <- function(in_ssn, organism, formula_root, 
                                  hydro_var, response_var, ssn_covtypes) {
  
  # hydro_var <- grep(paste0('^', hydro_var_str), 
  #                   names(in_ssn[[organism]]$ssn$obs), 
  #                   value=T)
  
  full_formula <- paste0(response_var, formula_root,' + ', hydro_var, ' + ',
                         hydro_var, ':country')
  
  
  ssn_list <- quick_ssn(in_ssn = in_ssn[[organism]]$ssn, 
                        in_formula = as.formula(full_formula),
                        ssn_covtypes = ssn_covtypes, #ssn_covtypes[sample(144, 10)],
                        estmethod = "ml")
  
  ssn_glance <- purrr::map(ssn_list, SSN2::glance,
                           .id=names(ssn_list)) %>%
    rbindlist(id='covtypes') %>%
    setorder(AIC) %>%
    as.data.table %>%
    .[, `:=`(
      organism = organism,
      response_var = response_var,
      hydro_var = hydro_var
    )]
  
  return(list(
    ssn_list=ssn_list,
    ssn_glance=ssn_glance
  ))
}

#------ compare standard hydro metrics with flow-duration curve-based metrics ---
# vars_list <- c('DurD', 'FreD')
# var_substr <- vars_list[[1]]
# in_cor_dt <- tar_read(cor_matrices_list)$div_bydrn
# color_list = drn_dt$color

plot_hydro_comparison <- function(in_cor_dt, color_list) {
  sub_dt_compare <- in_cor_dt[grep(var_substr, variable2),] %>%
    .[organism %in% org_list,] %>%
    .[(variable1 %in% c('invsimpson','Jtm1')),] 
  
  sub_dt_compare[, window_d := str_extract(variable2, '[0-9]+')] %>%
    .[, window_d := factor(window_d, levels=sort(unique(as.integer(window_d))))]
  sub_dt_compare[, root_var := gsub('[0-9]+past', '', variable2)]
  compare_cast <- dcast(sub_dt_compare[p_value < 0.1,], 
                        variable1+organism+country+window_d~root_var, 
                        value.var = 'correlation')
  ggplot(compare_cast, aes(x=abs(DurD), y=abs(PDurD))) +
    geom_point(aes(color=country)) +
    geom_abline() +
    facet_grid(organism~variable1) +
    coord_fixed() +
    theme_bw()
}  

#------ format_ssn_hydrowindow -------------------------------------------------
format_ssn_hydrowindow <- function(in_ssnmodels_combined) {
  in_ssnmodels_combined <- tar_read(ssnmodels_combined)
# 
#   lapply(
#     in_ssnmodels_combined, function(x) {
#       setnames(setDT(x[['ssn_glance']]), 'down_euc_types', 'covtypes')})
  
  ssn_hydrowindow_perf_allvars <- lapply(
    in_ssnmodels_combined, function(x) {
      merge(
        x[['ssn_glance']],
        data.table(covtypes=names(x[['ssn_list']]),
                   mod=x[['ssn_list']]),
        by='covtypes'
      )
  }) %>% 
    do.call(rbind, .) 
  
  #Identify covariance structure with the lowest AICc across all time windows 
  #for each hydrological variable
  ssn_hydrowindow_perf_allvars[, hydro_var_root := gsub(
    "(_*[0-9]+past)|(_*m[0-9]+)", "", hydro_var)] %>%
    .[, delta_AICc := (AICc - min(AICc)), by=hydro_var] %>%
    .[, delta_AICc_covtypemean := mean(delta_AICc), 
      by=.(covtypes, hydro_var_root)] %>%
    .[, window_d := str_extract(hydro_var, '([0-9]+(?=past))|((?<=m)[0-9]+)')] %>%
    .[, window_d := factor(window_d, levels=sort(unique(as.integer(window_d))))]
  
  #Plot model delta_AICc (epa stands for Epanechnikov kernel model)
  #Ordered globally
  ggplot(ssn_hydrowindow_perf_allvars, 
         aes(x = tidytext::reorder_within(
           as.factor(covtypes), 
           delta_AICc_covtypemean,
           hydro_var_root),
           y = delta_AICc, 
           fill = hydro_var_root)) +
    geom_boxplot() +
    scale_y_sqrt() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~hydro_var_root, scale='free')
  
  #Get best model table for all hydrowindows with embedded models
  ssn_hydrowindow_perf_sub <- ssn_hydrowindow_perf_allvars[, {
    min_aic_covtype <- .SD[which.min(delta_AICc_covtypemean), covtypes]
    .SD[covtypes == min_aic_covtype]
  }, by = hydro_var_root]
  
  #Get variance decomposition estimated by model
  ssn_hydrowindow_varcomp <- ssn_hydrowindow_perf_sub[
    , SSN2::varcomp(mod[[1]]), 
    by=.(hydro_var, covtypes, hydro_var_root, window_d)] 
  
  #Plot variance decomposition
  plot_ssn_hydrowindow_varcomp <- ggplot(ssn_hydrowindow_varcomp,
                                         aes(x=window_d, y=proportion, fill=varcomp)) +
    geom_bar(stat = 'identity', position='stack') +
    facet_wrap(~hydro_var_root, scales='free') +
    coord_cartesian(expand = FALSE) +
    theme_bw()
  
  #Plot best model for each variable for single window
  ssn_hydrowindow_best <- ssn_hydrowindow_perf_sub[
    , .SD[which.min(AICc),], 
    by = hydro_var_root]
  
  #Get model predictions
  mod_preds <- lapply(seq(nrow(ssn_hydrowindow_best)), function(i) {
    # Subset the data.table to get the current row
    aug_data <- SSN2::augment(ssn_hydrowindow_best[i, mod[[1]]], 
                              drop=FALSE) %>%
      cbind(ssn_hydrowindow_best[i, .(hydro_var, covtypes, hydro_var_root, window_d)])
    return(aug_data)
  }) %>% rbindlist
  
  ggplot(mod_preds, aes(x=.fitted, y=richness, color=country)) +
    geom_abline() +
    geom_point() +
    geom_smooth(method='lm') +
    facet_grid(country~hydro_var) +
    theme_bw()

  plot(ssn_hydrowindow_perf_sub[1, mod[[1]]])
  ssn_hydrowindow_best[1, SSN2::tidy(mod[[1]])]
  ssn_hydrowindow_best[1, SSN2::summary(mod[[1]])]
  
  
  ssn_hydrowindow_perf_sub[1, sjPlot::plot_model(mod[[1]], type = "pred")]
}

#------ model_miv_t ------------------------------------------------------------
#Model for macroinvertebrates for individual sampling dates

# Regularization: Techniques like LASSO (L1 regularization), Ridge Regression 
# (L2 regularization), and Elastic Net (a combination) can help prevent overfitting 
# by penalizing model complexity. glmnet is a great package for this, and it's easily integrated with caret.

# in_allvars_merged <- tar_read(allvars_merged)
# in_ssn_eu <- tar_read(ssn_eu)
# in_cor_matrices <- tar_read(cor_matrices_list)
# library(glmulti)

model_miv_t <- function(in_ssn_eu, in_allvars_merged, 
                        in_cor_matrices, ssn_covtypes) {
  
  #Subset SSN and create distance matrices-----
  ssn_miv <- in_ssn_eu$miv_nopools$ssn
  SSN2::ssn_create_distmat(ssn_miv)
  
  allvars_dt <- as.data.table(ssn_miv$obs) %>%
    setorderv(c('country', 'site', 'campaign')) 

  hydro_candidates <- setdiff(
    in_allvars_merged$cols$hydro_con,
    c('doy', 'month', 'hy', 'reach_id', 'qsim', 'UID',
      'isflowing', 'reach_length', 'noflow_period',
      'noflow_period_dur', 'last_noflowdate', 'PrdD')
  )

  #Scale data-----
  allvars_dt[
    , (hydro_candidates) := lapply(.SD,
                                   function(x) base::scale(x, center=T, scale=T)),
    .SDcols = hydro_candidates]

  #Basic function to train an lm or lmer model-----
  basic_train_mod <- function(in_dt, mod_formula, mod_name, 
                              id_cols = c('site', 'date'), color_var='country') {
    dt_copy <- copy(in_dt)
    
    if (grepl('[|]', mod_formula)) {
      out_mod <- lmer(formula= as.formula(mod_formula),
                      data=in_dt)
    } else {
      out_mod <- lm(formula= as.formula(mod_formula),
                      data=in_dt)
    }

    resp_col <- trimws(gsub('[~].*', '', mod_formula))
    
    predm_col <- paste0('mod_', mod_name, '_predm')
    predc_col <- paste0('mod_', mod_name, '_predc')
    resid_col <- paste0('mod_', mod_name, '_resid')
    dt_copy [, `:=`(predm = predict(out_mod, re.form = NA),
                 predc = predict(out_mod, re.form = NULL),
                 resid = resid(out_mod)
                 )]
    
    comp_p <- ggplot(dt_copy ,
                    aes_string(x='predc', y=resp_col, color=color_var)) +
      geom_abline() +
      geom_point() +
      geom_smooth(se=F, method='lm') +
      coord_fixed() +
      labs(x='Predicted', y='Observed') +
      theme_bw()
      
    print(summary(out_mod))
    print(comp_p)
    print(paste('AIC:', AIC(out_mod)))
    print(MuMIn::r.squaredGLMM((out_mod)))
    
    #str_extract_all()
    
    return(list(
      mod = out_mod,
      summ = summary(out_mod),
      dt = dt_copy,
      plot = comp_p,
      AIC = AIC(out_mod),
      mod_name = mod_name
    ))
  }
  
  compute_resid_corr <- function(in_mod, in_group_vars = c("organism", "country"),
                         in_x_cols) {

    resid_corr <- compute_cor_matrix_inner(
      in_mod$dt,
      group_vars = in_group_vars,
      y_cols = 'resid',
      x_cols = in_x_cols,
      exclude_diagonal = FALSE) %>%
      .[variable1 == 'resid' & !is.na(correlation) & !is.na(p_value),]
    
    resid_abscor_geomean <- resid_corr[
      , exp(mean(log(abs(correlation)))),
      by='variable2']
    
    return(list(
      all=resid_corr,
      abs_geomean = resid_abscor_geomean
    ))
  }

  #1.1. Model invsimpson diversity without space x all countries ------------------
  #Examine correlations with invsimpson
  topcors_overall_invsimpson <- in_cor_matrices$div[
    variable1 == 'invsimpson' & organism == 'miv_nopools' &
      variable2 %in% hydro_candidates & !is.na(correlation),] %>%
    .[, abs_cor := abs(correlation)] %>%
    setorder(-abs_cor)

  topcors_bydrn_invsimpson <- in_cor_matrices$div_bydrn[
    variable1 == 'invsimpson' & organism == 'miv_nopools' &
      variable2 %in% hydro_candidates & !is.na(correlation),] %>%
    .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
             var_label = paste(variable2, round(correlation, 2))
    ), by=country] %>%
    dcast(cor_order~country, value.var = 'var_label', fill=NA)

  #Check n-1 temporal autocorrelation
  acf_invsimpson_dt <- allvars_dt[, list(
    ac1 = cor(invsimpson,
              lag(invsimpson),
              use='pairwise.complete.obs'),
    n = .N-1)
    , by=c('site', 'country')] %>% .[n > 2,]

  acf_invsimpson_p <- ggplot(acf_invsimpson_dt,
                             aes(x=country, y=ac1,color=country)) +
    geom_hline(yintercept=0) +
    geom_boxplot(alpha=1/2) +
    geom_text(aes(label=as.numeric(gsub('^[A-Z]+', '', site))), 
              position = position_jitter(seed = 1, width=0.2, height=0)
              ) +
    theme_classic()

  # To test:
  # DurD3650past, FreD3650past, DurD365past, DurD180past, FreD365past, FreD180past,
  # PDurD180past, PFreD180past, relF90past, maxPQ_30past, oQ10_10past, maxPQ_365past,
  # oQ10_90past, relF1825past, PmeanQ10past, relF60past, relF7mean

  invsimpson_miv_nopools_modlist <- list()

  invsimpson_miv_nopools_modlist[['null']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)',
    mod_name = 'null')

  # ggplot(allvars_dt, aes(x=env_PC1, y=invsimpson)) +
  #   geom_point(aes(color=country)) +
  #   geom_smooth(color='black', method='lm') +
  #   geom_smooth(aes(color=country), se=F, method='lm')

  invsimpson_miv_nopools_modlist[['env1_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country) + env_PC1',
    mod_name = 'env1')

  invsimpson_miv_nopools_modlist[['env2_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country) + env_PC1 + env_PC2',
    mod_name = 'env2')

  invsimpson_miv_nopools_modlist[['all1_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country) + env_PC1 + state_of_flow_tm1',
    mod_name = 'all1')
  
  invsimpson_miv_nopools_modlist[['all2_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country) + env_PC1 + stream_type',
    mod_name = 'all2')
  resid_check_2 <- compute_resid_corr(
    in_mod = invsimpson_miv_nopools_modlist[['all2_mod']],
    in_group_vars = c("country"),
    in_x_cols = c(hydro_candidates, 
                  in_cor_matrices$env_cols))$abs_geomean
  
  ggplot(invsimpson_miv_nopools_modlist[['all2_mod']]$dt,
         aes(x=meanQ1825past, y=resid, color=country)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    scale_x_log10() +
    theme_bw() +
    facet_wrap(~country)
  
  invsimpson_miv_nopools_modlist[['all4_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country) + meanQ1825past + env_PC1 + stream_type',
    mod_name = 'all4')
  
  invsimpson_miv_nopools_modlist[['all5_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country) + meanQ1825past + stream_type',
    mod_name = 'all5')
  
  invsimpson_miv_nopools_modlist[['all6_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1 + meanQ1825past||country) + stream_type',
    mod_name = 'all6')
  
  invsimpson_miv_nopools_modlist[['all6_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + stream_type',
    mod_name = 'all6')
  
  # The fact that Model 1 has a substantially lower AIC suggests that its approach
  # of fitting distinct, fixed slopes for each country provides a significantly
  # better fit to the data that outweighs its potentially higher number of fixed 
  # parameters compared to Model 2. This may happen for two reasons:
  #   The actual slopes for each country are quite different and don't conform well 
  #   to being random draws from a single distribution (as assumed by Model 2's random slopes).
  #   The number of groups (countries = 6) is small. Estimating a variance for
  #   random slopes based on only 6 groups can be unstable or less informative 
  #   than directly estimating the 6 slopes as fixed effects. Model 1 doesn't try
  #   to generalize how slopes vary; it just estimates each one.
  
  invsimpson_miv_nopools_modlist[['all7_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + stream_type + env_PC2',
    mod_name = 'all7')
  
  resid_check_6 <- compute_resid_corr(
    in_mod = invsimpson_miv_nopools_modlist[['all6_mod']],
    in_group_vars = c("organism", "country"),
    in_x_cols = c(hydro_candidates, 
                  in_cor_matrices$env_cols))$abs_geomean

  ggplot(invsimpson_miv_nopools_modlist[['all6_mod']]$dt,
         aes(x=relF60past, y=resid, color=country)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    theme_bw() +
    facet_wrap(~country)
  
  ggplot(invsimpson_miv_nopools_modlist[['all6_mod']]$dt,
         aes(x=relF60past, y=invsimpson, color=country)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    theme_bw() +
    facet_wrap(~country)
  
  ggplot(invsimpson_miv_nopools_modlist[['all6_mod']]$dt,
         aes(x=STcon_m90_undirected, y=resid, color=country)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    scale_x_log10() +
    theme_bw() +
    facet_wrap(~country)
  
  #There are opposite patterns between countries, but these patterns are visible
  #in both relF and STcon: residuals increase with proportion of network flowing 
  #for Croatia, France, and Hungary; decreases in Czech, Finland and Spain
  #These patterns are globally consistant across time windows
  
  invsimpson_miv_nopools_modlist[['all8_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + stream_type + relF60past',
    mod_name = 'all8')
  
  invsimpson_miv_nopools_modlist[['all9_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1 + relF60past ||country) +  meanQ1825past + meanQ1825past:country + stream_type',
    mod_name = 'all9')
  
  invsimpson_miv_nopools_modlist[['all10_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + stream_type + relF60past:country',
    mod_name = 'all10')
  car::vif(invsimpson_miv_nopools_modlist[['all10_mod']]$mod)
  
  invsimpson_miv_nopools_modlist[['all10_mod']]$plot +
    facet_wrap(~country)
  
  resid_check_10 <- compute_resid_corr(
    in_mod = invsimpson_miv_nopools_modlist[['all10_mod']],
    in_group_vars = c("organism", "country"),
    in_x_cols = c(hydro_candidates, 
                  in_cor_matrices$env_cols))$abs_geomean
  
  ggplot(invsimpson_miv_nopools_modlist[['all10_mod']]$dt,
         aes(x=oQ10_365past, y=invsimpson)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    theme_bw() 
  
  ggplot(invsimpson_miv_nopools_modlist[['all10_mod']]$dt,
         aes(x=oQ10_365past, y=invsimpson, color=country)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    theme_bw() 
  
  ggplot(invsimpson_miv_nopools_modlist[['all10_mod']]$dt,
         aes(x=DurD1825past, y=resid)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    theme_bw() 
  
  ggplot(invsimpson_miv_nopools_modlist[['all10_mod']]$dt,
         aes(x=FreD3650past, y=resid, color=country)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    theme_bw() 
  
  ggplot(invsimpson_miv_nopools_modlist[['all10_mod']]$dt,
         aes(x=FreD3650past, y=invsimpson, color=country)) +
    geom_point() +
    geom_smooth(method='lm', span=1, se=F) +
    theme_bw() 
  
  invsimpson_miv_nopools_modlist[['all11_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + relF60past + relF60past:country + DurD3650past',
    mod_name = 'all11')

  invsimpson_miv_nopools_modlist[['all12_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + relF60past + relF60past:country + stream_type + DurD1825past',
    mod_name = 'all12')
  
  invsimpson_miv_nopools_modlist[['all13_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + relF60past + relF60past:country + stream_type + FreD1825past',
    mod_name = 'all13')
  
  invsimpson_miv_nopools_modlist[['all15_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + relF60past+ relF60past:country + stream_type +  PmeanQ60past',
    mod_name = 'all15')
 
  invsimpson_miv_nopools_modlist[['all16_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + relF60past+ relF60past:country + stream_type +  PDurD365past',
    mod_name = 'all16')
  
  #Selected non-spatial model
  invsimpson_miv_nopools_modlist[['all14_mod']] <- basic_train_mod(
    in_dt=allvars_dt,
    mod_formula='invsimpson ~ (1|country)  + meanQ1825past + meanQ1825past:country + relF60past + relF60past:country + stream_type + FreD3650past',
    mod_name = 'all14')
  sjPlot::plot_model(invsimpson_miv_nopools_modlist[['all14_mod']]$mod, type='diag')
  sjPlot::plot_model(invsimpson_miv_nopools_modlist[['all14_mod']]$mod, type='re')
  sjPlot::plot_model(invsimpson_miv_nopools_modlist[['all14_mod']]$mod, 
                     type = "pred",  
                     terms = c("meanQ1825past", "country")) +
    coord_cartesian(ylim=c(0, 100)) +
    theme_bw()
  
  sjPlot::plot_model(invsimpson_miv_nopools_modlist[['all14_mod']]$mod, 
                     type = "pred",  
                     terms = c("relF60past", "country")) +
    coord_cartesian(ylim=c(0, 50)) +
    theme_bw()
  
  
  #1.2. Model invsimpson diversity with space x all countries-------------------
  #Examine Torgegram ------
  #describes how the semivariance (i.e. halved average squared difference)
  #between observations or residuals from the model change
  #with hydrologic or Euclidean distances
  tg_miv_invsimpson_null <- Torgegram(
    formula = invsimpson ~ country,
    ssn.object = ssn_miv,
    type = c("flowcon", "flowuncon", "euclid"),
    bins = 15,
    cutoff = 40000,
    partition_factor = ~ as.factor(campaign)
  ) %>%
    rbindlist(idcol = 'dist_type')

  ggplot(tg_miv_invsimpson_null, aes(x=dist, y=gamma, color=dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(span=1, method='lm', se=F) +
    scale_x_log10() +
    facet_wrap(~dist_type, scales='free')

  tg_miv_invsimpson_mod14 <- Torgegram(
    formula = invsimpson ~ country  + meanQ1825past:country + relF60past:country + stream_type + FreD3650past,
    ssn.object = ssn_miv,
    type = c("flowcon", "flowuncon", "euclid"),
    bins = 10,
    cutoff = 40000,
    partition_factor = ~ as.factor(campaign)
  ) %>%
    rbindlist(idcol = 'dist_type')

  ggplot(tg_miv_invsimpson_mod14, aes(x=dist, y=gamma, color=dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(span=1, method='lm', se=F) +
    scale_x_log10() +
    facet_wrap(~dist_type, scales='free')

  #Train SSN models -----------
  #-- Null model ----------------------------------------------------------------
  ssn_mod_null <- ssn_lm(
    formula = invsimpson ~ 1,
    ssn.object = ssn_miv,
    taildown_type = "spherical",
    euclid_type = "gaussian",
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = ~ country
  )
  summary(ssn_mod_null)
  SSN2::varcomp(ssn_mod_null)
  
  #-- Model selected based on LME before --------------------------------------
  ssn_mod14 <- mapply(function(down_type, euc_type) {
    print(paste(down_type, euc_type))
    out_ssn <- ssn_lm(
      formula = invsimpson ~ meanQ1825past:country + relF60past:country + stream_type + FreD3650past,
      ssn.object = ssn_miv,
      taildown_type = down_type,
      euclid_type = euc_type,
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = ~ country
    )
    return(out_ssn)
  },
  down_type = ssn_covtypes$down,
  euc_type = ssn_covtypes$euc,
  SIMPLIFY = FALSE) %>%
    setNames(ssn_covtypes$label)

  mod14_types_glance <- purrr::map(ssn_mod14, SSN2::glance,
                                      .id=names(ssn_mod14)) %>%
    rbindlist(id='down_euc_types')%>%
    setorder(AIC)


  summary(ssn_mod14[['mariah_gaussian']])
  SSN2::varcomp(ssn_mod14[['mariah_gaussian']])
  SSN2::tidy(ssn_mod14[['mariah_gaussian']], conf.int = TRUE)

  #preds <- predict(ssn_mod20[['linear_spherical']], ssn_miv$obs)

  #-- Simpler model ------------------------------------------------------------
  ssn_modsimp <- mapply(function(down_type, euc_type) {
    print(paste(down_type, euc_type))
    out_ssn <- ssn_lm(
      formula = invsimpson ~ env_PC2 + DurD3650past + DurD3650past:country + relF10past,
      ssn.object = ssn_miv,
      taildown_type = down_type,
      euclid_type = euc_type,
      additive = "afv_qsqrt",
      partition_factor = ~ country
    )
    return(out_ssn)
  },
  down_type = ssn_covtypes$down,
  euc_type = ssn_covtypes$euc,
  SIMPLIFY = FALSE) %>%
    setNames(ssn_covtypes$label)

  modsimp_types_glance <- purrr::map(ssn_modsimp, SSN2::glance,
                                   .id=names(ssn_modsimp)) %>%
    rbindlist(id='down_euc_types') %>%
    setorder(AIC)

  summary(ssn_modsimp[['linear_gaussian']])
  SSN2::varcomp(ssn_modsimp[['linear_gaussian ']])
  SSN2::tidy(ssn_modsimp[['linear_gaussian']], conf.int = TRUE)


  #2.1. Model Tmj1 diversity without space x all countries ------------------
  #Examine correlations with invsimpson
  topcors_overall_Jtm1 <- in_cor_matrices$div[
    variable1 == 'Jtm1' & organism == 'miv_nopools' &
      variable2 %in% hydro_candidates & !is.na(correlation),] %>%
    .[, abs_cor := abs(correlation)] %>%
    setorder(-abs_cor)

  topcors_bydrn_Jtm1 <- in_cor_matrices$div_bydrn[
    variable1 == 'invsimpson' & organism == 'miv_nopools' &
      variable2 %in% hydro_candidates & !is.na(correlation),] %>%
    .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
             var_label = paste(variable2, round(correlation, 2))
    ), by=country] %>%
    dcast(cor_order~country, value.var = 'var_label', fill=NA)

  invsimpson_miv_nopools_modlist <- list()




  #Check distributions of residuals
  if (length(in_mod$coefficients)>1) {
    nsp_diag <- gg_diagnose(in_mod)
  } else {
    nsp_diag <- NULL
  }
  MAE = list(
    lm = mae(in_mod$model$ddt_to_bdtopo_ddratio_ceind,
             in_mod$fitted.values)
  )

  MAPE = list(
    lm = mape(in_mod$model$ddt_to_bdtopo_ddratio_ceind,
              in_mod$fitted.values)
  )

  vars_to_try <- sample(hydro_candidates, 5)
  glmulti_result <- glmulti(invsimpson ~ .^2,
                            maxsize = 2,
                            data = allvars_dt[, c('invsimpson',
                                                  vars_to_try ,
                                                  paste0('env_PC', seq(1,4))),
                                              with=F],
                            crit = "aic", level = 2, fitfunction = "glm",
                            method='h')




}

################################################################################
################################################################################

#   alpha_env_melt <- melt(
#     alphadat_env_dt,
#     id.vars = idcols,
#     measure.vars = predcols
#   )
#   
#   alpha_env_cor_drn_type <- alpha_env_melt[, list(
#     cor = .SD[!is.na(value),
#               cor(richness, value, method='spearman')]),
#     , by=c('drn', 'variable', 'organism', 'stream_type', 'nsim')] %>%
#     .[,  list(mean_spr = mean(cor, na.rm=T),
#               sd_spr = sd(cor, na.rm=T)),
#       by=c('drn', 'variable', 'organism', 'stream_type')]
#   
#   #Compute correlation between all predictor variables and species richness
#   #by organism and drn, then average across DRNs
#   alpha_env_cor_drn <- alpha_env_melt[, list(
#     cor = .SD[!is.na(value),
#               cor(richness, value, method='spearman')]),
#     , by=c('drn', 'variable', 'organism', 'nsim')] %>%
#     .[,  list(mean_spr = mean(cor, na.rm=T), #Average across RF sims
#               sd_spr = sd(cor, na.rm=T)),
#       by=c('drn', 'variable', 'organism')]
#   
#   alpha_env_cor_avg <- alpha_env_cor_drn[
#     , list(mean_spr = mean(mean_spr, na.rm=T)),
#     by = c('variable', 'organism')]
#   
#   ggplot(alpha_env_cor_avg[
#     abs(mean_spr) > 0.2 &
#       !(organism %in% c('miv_nopools', 'miv_nopools_flying',
#                         'miv_nopools_nonflying', 'bac_biof')),],
#     aes(x=tidytext::reorder_within(variable, mean_spr, within=organism),
#         y=mean_spr)) +
#     geom_bar(aes(fill = mean_spr), stat='identity') +
#     scale_fill_distiller(palette='Spectral', direction=1,
#                          breaks = seq(-1, 1, 0.1)) +
#     tidytext::scale_x_reordered(name = 'Candidate predictor variable') +
#     scale_y_discrete(name = "Mean Spearman's correlation across DRNs", expand=c(0,0)) +
#     facet_wrap(~organism, scales='free', nrow = 1) +
#     coord_flip() +
#     theme_bw()
#   
#   
# }
# 
# env_dd_dep_cormat <- dcast(env_dd_dep_cor, NOM+INSEE_DEP~description,
#                            value.var='cor')
# mat_names <-  env_dd_dep_cormat$NOM
# env_dd_dep_cormat <- as.matrix(env_dd_dep_cormat[, -c('NOM', 'INSEE_DEP')])
# row.names(env_dd_dep_cormat) <- mat_names
# 
# #Heatmap of correlation between drainage density ratio and variables----------
# colnames(env_dd_dep_cormat) <- 
#   gsub("gw", "gw",
#        gsub("surface water", "sw",
#             gsub("water withdrawals", "ww", colnames(env_dd_dep_cormat))))
# 
# max(env_dd_dep_cormat, na.rm=T)
# min(env_dd_dep_cormat, na.rm=T)
# 
# env_dd_dep_cormat_avg_morecl <- env_dd_dep_cormat[env_dd_dep_hclust_avg$order,]
# #as.data.table(env_dd_dendo_avg_morecl[[1]])[order(ID),order(gclass)],]
# 
# class_colors_avg_morecl <- classcol[as.data.table(
#   env_dd_dendo_avg_morecl[[1]])[order(gclass, ID),gclass]]
# 
# env_ddratio_corheatmap_avg_morecl <- 
#   ggcorrplot(corr=env_dd_dep_cormat_avg_morecl, #method = "circle", 
#              #hc.order = TRUE, hc.method = 'average',
#              lab=T, lab_size =3,
#              digits=1, insig='blank',
#              outline.color = "white") +
#   scale_fill_distiller(
#     name=str_wrap("Correlation coefficient Spearman's rho", 20),
#     palette='RdBu', 
#     limits=c(-0.8, 0.81), 
#     breaks=c(-0.7, -0.5, 0, 0.5, 0.7))  +
#   coord_flip() +
#   theme(axis.text.y = ggtext::element_markdown(
#     colour = class_colors_avg_morecl))
# 
# #Plot heatmap of env-dd correlations, clustered
# env_dd_dep_cormat_ward_morecl <- env_dd_dep_cormat[env_dd_dep_hclust_ward$order,]
# 
# class_colors_ward_morecl <- classcol[as.data.table(
#   env_dd_dendo_ward_morecl[[1]])[order(gclass, ID),gclass]]
# 
# 
# env_ddratio_corheatmap_ward_morecl <- 
#   ggcorrplot(env_dd_dep_cormat_ward_morecl, #method = "circle", 
#              #hc.order = TRUE, hc.method = 'average',
#              lab=T, lab_size =3,
#              digits=1, insig='blank',
#              outline.color = "white") +
#   scale_fill_distiller(
#     name=str_wrap("Correlation coefficient Spearman's rho", 20),
#     palette='RdBu', 
#     limits=c(-0.81, 0.81), 
#     breaks=c(-0.7, -0.5, 0, 0.5, 0.7))  +
#   coord_flip() +
#   theme(axis.text.y = element_text(
#     colour = class_colors_ward_morecl))

#------ tabulate_cor_matrix ----------------------------------------------------
#------ plot_spdiv_local -------------------------------------------------------
#------ plot_spdiv_drn ---------------------------------------------------------
#------ train_lme_aspatial -----------------------------------------------------
#------ analyze_residuals_autocor ----------------------------------------------
#------ train_lme_ssn ----------------------------------------------------------


#------ model_SSN --------------------------------------------------------------

# in_ssn <- tar_read(ssn_list)[[in_country]]$ssn
# SSN2::ssn_create_distmat(in_ssn, overwrite=TRUE)



# merge_
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

###################### NOT USED #####################################################
#------ create_ssn_drn ---------------------------------------------------------
# in_country <- 'Croatia'
# in_network_path = tar_read(network_ssnready_gpkg_list)[[in_country]]
# in_sites_path = tar_read(site_snapped_gpkg_list)[[in_country]]
# in_barriers_path = tar_read(barrier_snapped_gpkg_list)[[in_country]]
# out_dir = 'results/ssn'
# out_ssn_name = paste0(in_country, '_drn')
# overwrite=T

create_ssn_drn <- function(in_network_path,
                           in_sites_path,
                           in_barriers_path,
                           in_hydromod,
                           out_dir,
                           out_ssn_name,
                           overwrite = T) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  lsn_path <- file.path(out_dir,
                        tools::file_path_sans_ext(
                          sub('[.](?=(shp|gpkg)$)', '_lsn.',
                              basename(in_network_path), perl=T)
                        ))
  
  #Read input network
  net <- st_read(in_network_path) %>%
    st_cast("LINESTRING") %>%
    #Make sure that the geometry column is equally named regardless 
    #of file format (see https://github.com/r-spatial/sf/issues/719)
    st_set_geometry('geometry') %>%
    merge(in_hydromod$data_all[date < as.Date('2021-10-01', '%Y-%m-%d'),  #Link q data
                               list(mean_qsim = mean(qsim, na.rm=T)), 
                               by=reach_id],
          by.x = 'cat', by.y = 'reach_id') 
  
  #Build landscape network
  edges_lsn <- SSNbler::lines_to_lsn(
    streams = net,
    lsn_path = lsn_path,
    check_topology = TRUE,
    snap_tolerance = 0.05,
    topo_tolerance = 20,
    overwrite = overwrite
  )
  
  #Incorporate sites into the landscape network
  sites_lsn <- SSNbler::sites_to_lsn(
    sites = st_read(in_sites_path),
    edges =  edges_lsn,
    lsn_path = lsn_path,
    file_name = "sites",
    snap_tolerance = 5,
    save_local = TRUE,
    overwrite = TRUE
  )
  
  sites_list <- list(sites = sites_lsn)
  
  #Incorporate barriers into the landscape network
  #Only keep barriers over 2 m and under 100 m snap from network
  barriers_sf_sub <- read_sf(in_barriers_path) %>%
    filter((!is.na(Height) & Height > 2) & snap_dist_m < 100)
  
  if (nrow(barriers_sf_sub) > 0) {
    barriers_lsn <- sites_to_lsn(
      sites = barriers_sf_sub,
      edges =  edges_lsn,
      lsn_path = lsn_path,
      file_name = "barriers",
      snap_tolerance = 5,
      save_local = TRUE,
      overwrite = TRUE
    )
    
    sites_list$barriers <- barriers_lsn
  }
  
  #Calculate upstream distance
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
  
  #Compute segment Proportional Influence (PI) and Additive Function Values (AFVs)
  if (min(net$mean_qsim) > 0) {
    edges_lsn$mean_qsim_sqrt <- sqrt(edges_lsn$mean_qsim)
    
    edges_lsn <- afv_edges(
      edges = edges_lsn,
      infl_col = "mean_qsim_sqrt",
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
    
  } else {
    stop("Trying to use mean discharge to compute Additive Function Values (AFVs),
         but there are 0s in the discharge column.")
  }
  
  
  #Assemble the SSN
  out_ssn_path <- paste0(out_dir, out_ssn_name)
  
  out_ssn <- ssn_assemble(
    edges = edges_lsn,
    lsn_path = lsn_path,
    obs_sites = sites_list_lsn$sites,
    ssn_path = out_ssn_path,
    import = TRUE,
    check = TRUE,
    afv_col = "afv_qsqrt",
    overwrite = TRUE
  )
  
  # ggplot() +
  #   geom_sf(
  #     data = out_ssn$edges,
  #     color = "medium blue",
  #     aes(linewidth = mean_qsim_sqrt)
  #   ) +
  #   scale_linewidth(range = c(0.1, 2.5)) +
  #   geom_sf(
  #     data = out_ssn$obs,
  #     size = 1.7,
  #     aes(color = id)
  #   ) +
  #   coord_sf(datum = st_crs(out_ssn$edges)) +
  #   labs(color = "ID", linewidth = "sqrt(Q average)") +
  #   theme(
  #     legend.text = element_text(size = 6),
  #     legend.title = element_text(size = 8)
  #   )
  
  return(list(
    path = out_ssn_path,
    ssn = out_ssn)
  )
}



###################### OLD #####################################################
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


