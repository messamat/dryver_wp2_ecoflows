################################################################################
#---------------------------------- utility functions ----------------------------------------------
#------ download_unzip ---------------------------------------
#' Download and unzip a file from a URL to a specified directory.
#'
#' @param url URL of the file to download.
#' @param out_dir Output directory where the file will be unzipped.
#' @param download_mode Mode for download.file (default 'wb' for binary).
#' @param out_zip Optional: Path for the downloaded zip file.
#' @return Output directory.
download_unzip <- function(url, out_dir, download_mode='wb', out_zip=NULL) {
  # if (!dir.exists(out_dir)) {
  #   dir.create(out_dir)
  # }
  # download and unzip a file from a URL to a directory
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
#'
#' @param dates Date/time vector (will be coerced to date using 
#               as.Dates())
#' @param month Month defining start of hydrological year (default 9).
#' @param day_of_the_year Logical: Compute day of hydrological year.
#' @param seasons Named list of months for seasons (1:12; or days if '
#' seasons_by_day' is TRUE).
#' @param seasons_by_day Logical: Are seasons defined by days.
#' @param minimal Logical: Return minimal result (default FALSE).
#' @return Data.frame with hydrological year, days, and seasons.
hsaHydroYearSeasons <- function(dates, month = 9, day_of_the_year = TRUE, 
                                seasons = list("SummerFall" = 5:10, "WinterSpring" = c(11:12, 1:4)),
                                seasons_by_day = FALSE, minimal = FALSE) {
  dates <- as.Date(dates)
  
  # get hydrological year based on month threshold
  m <- as.numeric(format(dates, format = "%m"))
  hy <- y <- as.numeric(format(dates, format = "%Y"))
  m_prevy <- !m%in%c(month:12)
  hy[m_prevy] <- hy[m_prevy] - 1
  
  # get days of the year
  if (day_of_the_year) {
    j <- as.numeric(format(dates, format = "%j"))
    start_hy <- as.Date(paste0(y, "-", month, "-1"))
    start_y <- as.Date(paste0(y, "-1-1"))
    #Calculate hydrological day of year for each date
    j_hy <- as.numeric(format(dates - start_hy + start_y, format = "%j"))
  } else {
    j_hy <- NA
  }
  
  # get seasons
  if (!is.null(seasons) && is.list(seasons)) {
    # check season formatting (must be atomic vectors, not lists)
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
#' Assess validity of hydrological years based on completeness and missing values.
#'
#' @param hy Vector. Hydrological year assignments (will be coerced to factor)
#' @param x Matrix or vector with same length as hy, used to check missing values.
#' @param n Integer. Minimum length of a complete hydrological year.
#' @param na.th Numeric. Maximum tolerated proportion (threshold) of missing values.
#' @return Logical vector indicating valid time steps.
#'  
#' @details Given a hydrological year vector and a matrix with the same number
#' of rows (and any number of column), this function returns a logical
#' vector of the same length as the hydrological year vector that
#' indicates whether or not the time steps are part of a valid
#' hydrological year. Validity of a hydrological year is assessed 
#' according the two following rules:
#'  - is the year complete (i.e. is there at least 'n' (365) days)?
#'  - is there less than 'na.th' (proportion) missing values in 
#'    the year?
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
#' Read variable from NetCDF and format with hydrological year and season info.
#'
#' @param nc NetCDF object.
#' @param varname Character. Variable name to extract.
#' @param reachID Vector of reach IDs.
#' @param dates Vector of dates.
#' @param selected_sims Optional: which simulation layers to select.
#' @return List: formatted data and date info.
#' 
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
  
  # Internal helper to format NetCDF variable array into long-table.
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
  #Otherwie process single layer
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
#' Merge a list of data.tables, adding suffixes by source table name.
#'
#' @param dt_list List of data.tables.
#' @param by Columns to merge by.
#' @param all Logical: Merge all rows.
#' @param sort Logical: Sort merged table.
#' @param set_suffix Logical: Add suffix by name.
#' @return Merged data.table.
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
#' Create spatial point features from site data.
#'
#' @param in_dt Input data.table of sites.
#' @param lon_col Longitude column name.
#' @param lat_col Latitude column name.
#' @param out_points_path Path for output points file.
#' @param columns_to_include Optional: columns to include in output.
#' @return Path to output points file.
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
#' Compute left-continuous ECDF values by group.
#'
#' @param in_dt Input data.table.
#' @param ecdf_column Column to compute ECDF for.
#' @param grouping_columns Columns to group by.
#' @param na.rm Remove NAs in ECDF column.
#' @return Data.table with ECDF values by group.
#' 
#' @details If nine 0s and one 1, then 0s are given 0 and 1 is given 0.9
#' @source #https://stats.stackexchange.com/questions/585291/is-there-an-equivalent-to-an-ecdf-with-a-sign
#' https://math.stackexchange.com/questions/1807120/why-arent-cdfs-left-continuous/1807136#1807136
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
#' Merge multiple ECDF lookups into one table.
#'
#' @param in_dt Input data.table.
#' @param ecdf_columns List of columns to compute ECDF for.
#' @param grouping_columns Columns to group by.
#' @param keep_column Column to keep in output.
#' @param na.rm Logical: remove NAs.
#' @return Merged data.table with ECDF columns.
#' 
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

#------ identify_drywet6mo   ------------------------------------------------
identify_drywet6mo <- function(in_dt, flow_col='isflowing', jday_col = 'doy') {
  
  #Add 3 months before and after record to avoid having NAs on the edges
  fill_dt <- data.table(rep(NA,91)) %>% setnames(flow_col)
  
  in_dt <- rbind(fill_dt,
                 rbind(in_dt, fill_dt, 
                       fill=T),
                 fill=T)
  
  #Computer total number of no-flow days in the 6-month period centered around
  #each day in the record
  in_dt <- in_dt[, noflow_6morollsum := frollsum(
    isflowing==0,
    n=183, align='center', hasNA=T, na.rm=T)] %>%
    .[ 92:(.N-91), ]
  
  #Find the Julian day with the most number of no-flow days on interannual average
  driest_6mocenter <- in_dt[get(jday_col) < 366 #can lead to artefacts if drying was infrequent but occurred during that time
                            , as.Date(
                              unique(get(jday_col))[
                                which.max(.SD[, mean(noflow_6morollsum, na.rm=T),by=jday_col]$V1)],
                              origin=as.Date("1970-01-01"))
  ]
  
  #Identify the 6-month period centered on that julian day as the dry period
  driest_6moperiod <- data.table(
    dry_6mo = T,
    jday = as.numeric(format(
      seq(driest_6mocenter-91, driest_6mocenter+91, by='day'), '%j'))
  )
  
  #Identify the other julian days as the wet period
  in_dt <- merge(in_dt, driest_6moperiod, 
                 by.x=jday_col, by.y='jday', all.x=T) %>%
    .[is.na(dry_6mo), dry_6mo := F]
  
  in_dt[, noflow_6morollsum := NULL]
  
  return(in_dt[order(date),])
}

#------ compute_hydrostats_intermittence ---------------------------------------
#' Compute intermittence statistics for hydrological data.
#'
#' @param in_hydromod_dt Hydrological time series data.table.
#' @param in_sites_dt Sites data.table.
#' @param scale Scale for computation: 'all', 'drn', or 'site'.
#' @return List of statistics by scale.
compute_hydrostats_intermittence <- function(in_hydromod_dt,
                                             in_sites_dt,
                                             scale = 'all') {
  rollingstep_yr <- c(10, 30)
  
  # -- Compute network-wide statistics ---------------------------------------
  if (scale %in% c('all', 'drn')) {
    #minimum 7-day and 30-day average over previous year
    #average in previous 10, 30, 45, 60, 90, 120, 180, 365, 365*5, 365*10 
    
    #Total network length for normalization
    total_reach_length <- unique(in_hydromod_dt, by='reach_id')[, sum(reach_length)]
    
    #RelFlow: Proportion of network length with flowing conditions (opposite of RelInt)
    relF_dt <- in_hydromod_dt[isflowing == 1,
                              list(relF = sum(reach_length)/total_reach_length)
                              , by=.(date, hy)] %>%
      setorder(date)
    
    #ggplot(relF_dt, aes(x=date, y=relF)) + geom_line()
    
    #Compute min 7-day and 30-day relF
    relF_dt[, `:=`(
      relF7mean = frollmean(x=relF, n=7, align='center', na.rm=T),
      relF30mean = frollmean(x=relF, n=30, align='center', na.rm=T)
    )] %>%
      .[, `:=`(
        relF7mean_yrmin = min(relF7mean, na.rm=T),
        relF30mean_yrmin = min(relF30mean, na.rm=T)
      ), by=hy] 
    
    
    #Compute previous mean over many windows
    meanstep <- c(10, 30, 60, 90, 120, 180, 365, 365*5, 365*10)
    relF_dt[, paste0("relF", meanstep, "past") := frollmean(relF, meanstep, na.rm=T)]
    
    
    relF_dt_yr <- relF_dt[!duplicated(hy), .(relF7mean_yrmin, hy)] %>%
      .[, paste0("relF7mean_yrmin_cv", rollingstep_yr, "yrpast") :=  
          frollapply(relF7mean_yrmin, n=rollingstep_yr, 
                     FUN=function(x) sd(x, na.rm=T)/mean(x, na.rm=T), 
                     align='right')
      ] %>% 
      .[, relF7mean_yrmin := NULL]
    
    relF_dt <- merge(relF_dt, relF_dt_yr, by='hy')
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
    
    # hydromod_dt_sites[, paste0("MaxConD", rollingstep, "past") :=
    #                     nafill(
    #                       frollapply(noflow_period_dur, n=rollingstep,
    #                                  FUN = function(x) max(x, na.rm=T),
    #                                  align='right')
    #                       , fill = 0
    #                     ),
    #                   by = .(reach_id)
    # ]
    
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
    
    
    #Compute annual timing and interannual variability
    hydromod_dt_yr <- hydromod_dt_sites[, list(
      DurD_yr = sum(isflowing==0)/.N,
      FreD_yr = length(unique(noflow_period[!is.na(noflow_period)])),
      FstDrE = .SD[isflowing==0, min(doy, na.rm=T)],
      meanConD_yr = .SD[!duplicated(noflow_period), mean(noflow_period_dur, na.rm=T)]             
    ), by=.(reach_id, year = as.integer(format(date, '%Y')))] %>%
      .[is.infinite(FstDrE), FstDrE := NA]
    
    
    #CV of annual number of no-flow days
    hydromod_dt_yr[
      , paste0("DurD_CV", rollingstep_yr, "yrpast") :=  
        frollapply(DurD_yr, n=rollingstep_yr, 
                   FUN=function(x) fifelse(mean(x, na.rm=T)==0,
                                           0,
                                           sd(x, na.rm=T)/mean(x, na.rm=T)
                   ), 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #CV of annual number of no-flow events
    hydromod_dt_yr[
      , paste0("FreD_CV", rollingstep_yr, "yrpast") :=  
        frollapply(FreD_yr, n=rollingstep_yr, 
                   FUN=function(x) fifelse(mean(x, na.rm=T)==0,
                                           0,
                                           sd(x, na.rm=T)/mean(x, na.rm=T)
                   ), 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #CV of average annual event duration
    hydromod_dt_yr[
      , paste0("meanConD_CV", rollingstep_yr, "yrpast") :=  
        frollapply(meanConD_yr, n=rollingstep_yr, 
                   FUN=function(x) fifelse(mean(x, na.rm=T)==0,
                                           0,
                                           sd(x, na.rm=T)/mean(x, na.rm=T)
                   ), 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #SD No-flow start date - Julian day of the first drying event
    hydromod_dt_yr[
      , paste0("FstDrE_SD", rollingstep_yr, "yrpast") :=  
        frollapply(FstDrE, n=rollingstep_yr, 
                   FUN=function(x) sd(x, na.rm=T)
                   , 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #SD6: seasonal predictability of no-flow events (Gallart et al. 2012)
    #Identify contiguous six months with the most zero-flow days for computing Sd6
    hydromod_dt_sites <- hydromod_dt_sites[, identify_drywet6mo(in_dt=.SD), 
                                           by=.(reach_id)]
    
    base_dt <- hydromod_dt_sites[, ym := format(date, '%Y%m')] %>%
      .[, any(!is.na(noflow_period)), 
        by=.(year = as.integer(format(date, '%Y')), ym, dry_6mo, reach_id)] %>%
      .[, .(n_months = sum(V1, na.rm=TRUE)), by=.(year, dry_6mo, reach_id)] %>%
      dcast(year+reach_id~dry_6mo, value.var = 'n_months')
    
    rolling_sd6 <- function(dt, k) {
      colname <- paste0("sd6_", k, "yrpast")
      
      dt[, (colname) :=
           frollapply(year, n = k, align = "right",
                      FUN = function(yrs) {
                        subdt <- .SD[year %in% yrs]   # only this reach_id’s rows
                        val <- 1 - subdt[, (mean(`FALSE`, na.rm=TRUE) / mean(`TRUE`, na.rm=TRUE))]
                        ifelse(is.na(val) | mean(`TRUE`, na.rm=TRUE) == 0, 1, val)
                      }),
         by = reach_id] 
    }
    
    sd6_dt <- copy(base_dt)
    sd6_dt <- rolling_sd6(sd6_dt, 10)
    sd6_dt <- rolling_sd6(sd6_dt, 30)
    
    #Merge everything
    hydromod_dt_sites <- merge(in_sites_dt[, .(site, reach_id)], 
                               hydromod_dt_sites, 
                               by = 'reach_id', all.x = T, all.y = F,
                               allow.cartesian = T) %>%
      .[, year := as.integer(format(date, '%Y'))] %>%
      merge(hydromod_dt_yr, by=c('reach_id', 'year')) %>%
      merge(sd6_dt, by=c('reach_id', 'year')) %>%
      .[, `:=`(ym=NULL, dry_6mo=NULL, `FALSE`=NULL, `TRUE`=NULL)] %>%
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
#' Compute hydrological statistics for discharge.
#'
#' @param in_hydromod_dt Hydrological time series data.table.
#' @return Data.table with discharge statistics.
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
#' Create a two-point equidistant projection string for a bounding box
#'
#' Defines a standard two-point equidistant (tpeqd) projection based on the
#' bounding box of a spatial object. Useful for distance calculations that
#' require a locally suitable projection.
#'
#' @param x An `sf` object (or anything accepted by `sf::st_bbox`).
#'
#' @return A PROJ.4 projection string for a two-point equidistant projection
#'   with WGS84 ellipsoid and units in meters.
#'
#' @references
#' Based on: \url{https://gis.stackexchange.com/questions/313721/automatically-get-an-adequate-projection-system-based-on-a-bounding-box}
#'
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
#' Snap points to the nearest target geometry
#'
#' Snaps input points to the nearest target features (e.g., lines or polygons),
#' creating snapped points with optional attribute joins.
#'
#' @param in_pts A `SpatVector` of points (from `terra`).
#' @param in_target A `SpatVector` target geometry (e.g., lines).
#' @param sites_idcol Name of the column containing point IDs (used to filter duplicates).
#' @param attri_to_join Optional. Either a character vector of column names from
#'   `in_target` to join, or `"all"` to join all attributes. Default: NULL 
#'   (i.e. no attribute is joined)
#'
#' @return A `SpatVector` of snapped points with optional joined attributes.
#'
#' @details
#' - Snapping is performed using `terra::nearest` with line output, then
#'   converted to endpoints.
#' - If `attri_to_join` is provided, attributes of the nearest target are joined.
#'
#' @note
#' Uses a heuristic: filtering duplicates with `duplicated(values(.)[, sites_idcol])`.
#' This may drop unintended points if IDs are not unique.
#' 
snap_points_inner <- function(in_pts,
                              in_target,
                              sites_idcol,
                              attri_to_join=NULL
) {
  #Snap points (fastest custom way in R, it seems):
  #Compute line between site and snapping place on nearest segment
  sitesnap_l <- terra::nearest(in_pts, in_target, centroids = F, lines = T)
  values(sitesnap_l) <- values(in_pts)
  sitesnap_l$snap_dist_m <- perim(sitesnap_l)
  
  #convert the line to a point (the line's end point)
  sitesnap_p <- terra::as.points(sitesnap_l) %>%
    .[duplicated(values(.)[, sites_idcol]),]
  
  #Join attributes of nearest line to that point (optional)
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
#' Split a spatial line into equal-length segments
#'
#' Splits a `Lines` object (or coordinate matrix) into `n` segments of
#' a specified length, plus a remainder segment if the total length does not
#' divide evenly.
#'
#' @param line A `sp::Lines` object or a coordinate matrix.
#' @param n Number of equal-length segments to create.
#' @param length Desired length of each segment (in the same units as `line`).
#' @param debug Logical. If `TRUE`, plots the line and split points.
#'
#' @return A list of coordinate matrices, each representing a segment.
#'
#' @note
#' - Multiple lines are not supported; the function stops if multiple lines are detected.
#' - The name `length` shadows the base R function `length()`. This is not an error,
#'   but may cause confusion.
#' - Original author: Miguel Porto, from https://github.com/miguel-porto/fix-streams
split_sp_line <- function(line, n, length, debug = F) {
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
#' Fix complex confluences in stream networks
#'
#' Adjusts complex stream confluences (nodes with >2 incoming streams) in a
#' `SpatialLinesDataFrame`. Outermost streams' end vertices are shifted 
#' downstream by a given step distance, new nodes are created, 
#' and affected lines are split.
#'
#' @param shp A `SpatialLinesDataFrame` with `FROM_NODE` and `TO_NODE` fields.
#' @param from Column name of upstream node IDs (default `"FROM_NODE"`).
#' @param to Column name of downstream node IDs (default `"TO_NODE"`).
#' @param step Numeric. Distance (in map units) by which streams are adjusted downstream.
#' @param fields_to_keep Optional. Column names to keep from the original attributes.
#'
#' @return A corrected `SpatialLinesDataFrame` with confluences fixed.
#'
#' @details
#' - The function assumes the network is otherwise topologically correct.
#' - No input validation is performed; invalid inputs may cause errors.
#'
#' @references
#' Based on: \url{https://github.com/miguel-porto/fix-streams}
#'
#' @examples
#' \dontrun{
#' rios <- rgdal::readOGR("streams_Pt.shp", "streams_Pt")
#' corrected <- fix_confluences_inner(rios, step=10)
#' rgdal::writeOGR(corrected, "streams_corrected.shp", "streams_corrected", "ESRI Shapefile")
#' }
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
#' Compute and flatten correlation matrices
#'
#' Computes pairwise correlations (Spearman or Pearson) between sets of columns,
#' optionally grouped by specified variables, and returns a long-format data.table.
#'
#' @param in_dt A `data.table` containing the data.
#' @param group_vars Optional. Character vector of grouping column names.
#' @param x_cols Character vector of column names for the first set of variables.
#' @param y_cols Optional. Character vector of column names for the second set of
#'   variables. If `NULL`, correlations are computed within `x_cols`.
#' @param correlation_type Character. Type of correlation: `"spearman"` (default)
#'   or `"pearson"`.
#' @param exclude_diagonal Logical. If `TRUE` (default), removes self-correlations
#'   (variable1 == variable2).
#'
#' @return A `data.table` with columns:
#' - `variable1`, `variable2`: Names of the correlated variables.
#' - `correlation`: Correlation coefficient.
#' - `p_value`: Associated p-value.
#'
#' @note
#' Relies on `Hmisc::rcorr`, which may return `NA` values for correlations with
#' insufficient data.
#' 
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
#' Fill NA values hierarchically: by site, then by country, then overall mean.
#'
#' This function fills missing values in specified columns, first using the site average,
#' then the country average, and finally the overall average, ensuring no NA remains.
#'
#' @param dt data.table. The data to fill.
#' @param cols_to_fill Character vector. Names of columns to fill.
#' @param site_col Character. Site column name.
#' @param country_col Character. Country column name.
#' @return data.table. Copy of input with NAs filled hierarchically.
#' 
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
#'
#' @param in_dt data.table. Data to transform.
#' @param in_cols_to_ordinate Character vector. Columns to include in ordination.
#' @param num_pca_axes Integer. Number of principal axes to keep (default 4).
#' @return List. Contains PCA object, scores as data.table, and transformed data.
trans_pca <- function(in_dt, in_cols_to_ordinate, num_pca_axes = 4) {
  dt_copy <- copy(in_dt)
  
  # 1. Determine optimal lambda for Box-Cox
  dt_copy[, (in_cols_to_ordinate) := lapply(.SD, function(x) {
    # Ensure all values > 0 for Box-Cox
    min_positive <- min(x[x > 0], na.rm = TRUE)
    trans_shift <- if (is.finite(min_positive)) min_positive * 0.01 else 1
    return(x + trans_shift)
  }), .SDcols = in_cols_to_ordinate]
  
  # Find optimal lambda for each column
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
    .[, .SD, .SDcols = seq(1, num_pca_axes)] #keep only selected # of axes
  
  return(list(
    pca = pca_out,
    pca_dt = pca_out_dt,
    trans_dt = dt_copy
  ))
}
#------ trans_pca_wrapper ------------------------------------------------------
#' Wrapper for PCA transformation by group or overall.
#'
#' Performs PCA using trans_pca, optionally by group, and returns scores merged with IDs.
#'
#' @param in_dt data.table. Data to transform.
#' @param in_cols_to_ordinate Character vector. Columns to include in ordination.
#' @param id_cols Character vector. Columns to use as IDs.
#' @param group_cols Character vector or NULL. If provided, PCA is run within each group.
#' @param num_pca_axes Integer. Number of principal components to keep.
#' @return List. Includes PCA object, transformed data, scores, and plot (if available).
#'
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
#------ get_hydro_var_root -----------------------------------------------------
get_hydro_var_root <- function(dt, in_place=TRUE) {
  if (!in_place) {
    dt <- copy(dt)
  }
  dt[, hydro_var_root := gsub(
    "(_*[0-9]+(yr)*past)|(_*m[0-9]+)", "", hydro_var)] %>%
    .[, window_d := str_extract(hydro_var,
                                '([0-9]+(?=yrpast))|([0-9]+(?=past))|((?<=m)[0-9]+)')] %>%
    .[, window_d := factor(window_d, levels=sort(unique(as.integer(window_d))))]
  
  if (!in_place) {
    return(dt)
  }
}

#------ get_full_hydrolabel ----------------------------------------------------
get_full_hydrolabel <- function(in_hydro_vars_dt,
                                in_hydro_var) {
  h_unit <- ifelse(grepl(pattern='.*yrpast$', in_hydro_var), 'y', 'd')
  
  return(in_hydro_vars_dt[hydro_var == in_hydro_var,
                          paste(hydro_label, ': past',
                                window_d, h_unit)])
}

#------ get_ssn_emmeans -------------------------------------------------------
get_ssn_emmeans <- function(in_mod, 
                            in_pred_var, interaction_var,
                            in_pred_var_label=NULL,
                            in_emm_dt=NULL, plot=T, label_pred_var=T) {
  
  if (is.null(in_emm_dt)) {
    in_ssn <- in_mod$ssn.object
    newdat <- list(
      seq(min(in_ssn$obs[[in_pred_var]], na.rm = TRUE),
          max(in_ssn$obs[[in_pred_var]], na.rm = TRUE),
          length.out = 50)
    )
    names(newdat) <- in_pred_var
    
    response_var <- all.vars(in_mod$formula)[[1]]
    
    #Marginal means plot
    emm_dt <- emmeans(in_mod, 
                      data=in_mod$ssn.object$obs, 
                      specs= as.formula(paste0('~', in_pred_var, '|', interaction_var)),
                      at = newdat,
                      type = "response") %>%
      data.frame %>%
      setDT %>%
      merge(as.data.table(in_ssn$obs)[, list(max_pred_var = max(get(in_pred_var))),
                                      by=interaction_var], by=interaction_var) %>%
      .[get(in_pred_var) < max_pred_var] %>%
      .[, `:=`(response_var=response_var,
               pred_var_name = in_pred_var)] %>%
      setnames(in_pred_var, 'pred_var')
    
  } else if (is.data.table(in_emm_dt) && 
             all(c('pred_var', interaction_var, 'emmean')  %in% names(in_emm_dt))
             ){
    emm_dt <- in_emm_dt
    response_var <- emm_dt[1, response_var]
  }
  
  if (plot) {
    emm_plot <- emm_dt %>%
      ggplot(aes(x=pred_var, fill=get(interaction_var))) +
      geom_line(aes(y=emmean, color=get(interaction_var)), linewidth=2) +
      geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=0.1) +
      scale_y_continuous(name=paste('Estimated marginal mean:', 
                                    Hmisc::capitalize(response_var))) +
      scale_color_discrete(name=Hmisc::capitalize(interaction_var)) + 
      scale_fill_discrete(name=Hmisc::capitalize(interaction_var)) + 
      theme_classic()
    
    if (label_pred_var && !is.null(in_pred_var_label)) {
      emm_plot <- emm_plot + scale_x_continuous(name=in_pred_var_label)
    }
  } else {
    emm_plot <- NULL
  }
  
  return(list(
    dt=emm_dt,
    plot=emm_plot
  ))
}

#------ get_ssn_emtrends ------------------------------------------------------
get_ssn_emtrends <- function(in_mod, in_pred_var, 
                              in_pred_var_label, interaction_var,
                              in_emtrends_dt=NULL, plot=T) {
  if (is.null(in_emtrends_dt)) {
    response_var <- all.vars(in_mod$formula)[[1]]
    
    emtrends_dt <- emtrends(object = in_mod, 
                            specs = as.formula(paste('~', interaction_var)), 
                            var = in_pred_var,
                            data=in_mod$ssn.object$obs) %>% 
      as.data.frame %>%
      setDT %>%
      .[, `:=`(response_var=response_var,
               pred_var_name = in_pred_var)] %>%
      setnames(paste0(in_pred_var, '.trend') , 'trend')
    
  }  else if (is.data.table(in_emm_dt)) {
    emtrends_dt <- in_emtrends_dt
    response_var <- emtrends_dt[1, response_var]
  }

  if (plot) {
    emtrends_plot <- ggplot(emtrends_dt, aes(x = get(interaction_var), 
                                             y = trend)) +
      geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL,
                          color=get(interaction_var))) +
      geom_hline(yintercept=0, linetype=2) +
      theme_minimal() +
      labs(y = paste("Estimated slope of", response_var),
           x = Hmisc::capitalize(interaction_var)) +
      coord_flip()
  } else {
    emtrends_plot <- NULL
  }
  
  return(list(
    dt = emtrends_dt,
    plot = emtrends_plot
  ))
}

################################################################################
#---------------------------------- workflow functions ---------------------------------------------
# path_list = tar_read(bio_data_paths)
# in_metadata_edna <- tar_read(metadata_edna)

#------ define_hydromod_paths --------------------------------------------------
#in_hydromod_dir <- hydromod_present_dir

#' Define hydrological model input and output data paths for the 6 DRNs
#' 
#' Creates a table of file paths related to hydrological
#' simulations for several European catchments. It standardizes the paths
#' to key files: NetCDF outputs, shapefiles, parameter files, and reach IDs.
#'
#' @param in_hydromod_dir Character. Path to the root hydrological model directory.
#'
#' @return A `data.table` with country, catchment, best simulation number,
#'   and standardized file paths for each catchment.
#'   
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

#' Compute upstream drainage area for the Genal basin
#'
#' Uses HydroSHEDS 90m flow direction data to compute upstream drainage
#' areas for the Genal River (southern Spain). Optionally reuses a
#' pre-computed raster if it exists on disk.
#'
#' @param in_flowdir_path Character. Path to HydroSHEDS flow direction raster.
#' @param outdir Character. Output directory for saving computed drainage area.
#'
#' @return A `data.table` with site codes (GEN01–GEN26) and estimated
#'   basin drainage areas in km².
#'
#' @note The DEM-based method is commented out, as it requires
#'   hydrological conditioning. HydroSHEDS is used instead.
#' @note Site areas are currently based on **visual matches** — this may
#'   introduce subjectivity. Consider automating extraction via snapping.
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
#' Extract hydrological model data from NetCDF
#'
#' Reads a NetCDF hydrological model file and extracts a variable (e.g., 
#' flow intermittence) across reaches and simulation dates.
#'
#' @param hydromod_path Character. Path to the NetCDF file.
#' @param varname Character. Variable name to extract from the NetCDF file.
#' @param selected_sims Optional. Integer vector of simulation indices to 
#'   subset (default = `NULL`, all simulations used).
#'
#' @return A `data.table` with reach IDs, dates, and values for the requested variable.
#'
#' @details
#' - Dates are converted from days since 1950-01-01 (NetCDF convention).
#' - Relies on an external helper function `get_nc_var_present()` 
#'   
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

#' Read and clean environmental datasets
#'
#' Merges Annika’s version of environmental data with the common 
#' project dataset, applies multiple corrections, fills missing 
#' values, and harmonizes formats. This prepares the environmental 
#' dataset for downstream biodiversity and hydrology analyses.
#'
#' @param in_env_data_path_annika Character. Path to Annika’s dataset.
#' @param in_env_data_path_common Character. Path to the common dataset.
#'
#' @return A cleaned `data.table` of environmental data with harmonized variables.
#'
#' @details
#' - Handles numerous data quality issues: missing values, typos, dry sites, 
#'   impossible measurements, etc.
#' - Applies imputations for missing geomorphological measures and 
#'   macroinvertebrate hydraulic attributes.
#' - Standardizes date formats and variable naming conventions.
#'
#' @note
#' - Some imputations (e.g., Hungary velocities, conductivity ratios) 
#'   are heuristic and could bias results.
#' - Uses campaign/site-specific manual corrections (e.g., GEN04 altitude).
#'
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
  
  # --- merging & checking differences between Annika’s and common datasets ---
  
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
  
  # --- merge into a cleaned dataset ---
  env_dt_merged <- merge(env_dt_annika[, -cols_common, with=F],
                         env_dt_common[, c('running_id', cols_common), with = F]) %>%
    setnames(tolower(names(.))) %>% #Convert all columns to lower case
    .[!(is.na(date) & is.na(state_of_flow)),] %>% #Remove records that were not sampled at all, but keep those that were simply dry
    setnames(c('bankfull_at_max__wetted_width_m', 'bankfull_at_min__wetted_width_m'),
             c('bankfull_at_max_wetted_width_m', 'bankfull_at_min_wetted_width_m'))
  
  # --- corrections for missing data ---
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

#' Read and clean biodiversity datasets
#'
#' Imports multiple biodiversity datasets (diatoms, fungi, macroinvertebrates, 
#' bacteria) and harmonizes formats, dates, site/campaign IDs, and metadata.
#'
#' @param path_list Named list of file paths. Names should correspond to organism types.
#' @param in_metadata_edna Data.table of eDNA metadata (sample type, site info).
#' @param include_bacteria Logical. If `FALSE`, bacterial datasets are dropped. Default = TRUE.
#'
#' @return A named list of cleaned `data.table`s, one per organism group.
#'
#' @details
#' - Fills missing sampling dates from eDNA metadata.
#' - Splits running IDs into site + campaign when needed.
#' - Removes bacteria pool samples (keeps separate `_nopools` tables).
#' - Corrects a few known typos in dates (e.g., fungi 2012 → 2021).
#'
read_biodt <- function(path_list, in_metadata_edna, 
                       in_miv_full_dt, include_bacteria=T) {
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
  
  #Separate OCH from EPT--------------------------------------------
  #Prepare reduced data for joining with full/raw data that contains family/taxagroup info
  miv_melt <- melt(dt_list$miv_nopools, 
                   id.vars=c('organism', 'running_id', "campaign"
                             ,"site", "country", "date", "sample.id" ),
                   variable.name = "taxon_name",
                   value.name = 'density'
                   )
  #Check that can safely clean
  gsub('[.]', '', gsub('(gen[.])|(sp[.])', '', names(dt_list$miv_nopools))) %>% 
    duplicated %>% sum
  #Clean name for joining
  miv_melt[, taxon_name_format := gsub('[.]|([.]sp)|\\s', '',
                                       gsub('(sp)[.]', '', taxon_name))
           ]
  
  #Clean raw data
  taxo_info <- setDT(in_miv_full_dt) %>%
    .[, taxon_name := tolower(`Genus / Higher taxonomic group`)] %>%
    .[, taxon_name_format := gsub('[.]|([.]sp)|[//]|\\s', '',
                           gsub('(sp|lv|ad)[.]', '', taxon_name))
    ] %>%
    .[!duplicated(taxon_name_format), 
      .(taxon_name_format, Taxagroup, Family, Subfamily)] %>%
    rbind(.[taxon_name_format=='limnephilus',] %>% #Add a missing family-level identified taxon
            .[, taxon_name_format := 'limnephilinaegen'])
  
  #Merge them
  miv_full <- merge(miv_melt, taxo_info, by='taxon_name_format', all.x=T)

  #check that everything was matched
  check <- miv_full[is.na(Taxagroup),]
  check[!duplicated(taxon_name_format), .(taxon_name_format, taxon_name)]

  #Divide OCH and EPT
  no_ept_taxon_list <- miv_full[
    !(Taxagroup %in% c('Ephemeroptera', 'Plecoptera', 'Trichoptera')),
    as.character(unique(taxon_name))]
  no_och_taxon_list <- miv_full[
    !(Taxagroup %in% c('Odonata', 'Coleoptera', 'Heteroptera')),
    as.character(unique(taxon_name))]
  
  #Subset original table's columns based on selected taxa
  dt_list$miv_nopools_ept <- dt_list$miv_nopools[
    ,.SD, .SDcols = !no_ept_taxon_list] %>% 
    .[, organism := 'miv_nopools_ept']
  dt_list$miv_nopools_och <- dt_list$miv_nopools[
    ,.SD, .SDcols = !no_och_taxon_list]%>% 
    .[, organism := 'miv_nopools_och']
  
  #Fill NAs in dates with eDNA metadata if possible - remove pools--------------
  dt_list$dia_sedi[
    is.na(date),
    date := in_metadata_edna[sample_type=='sediment' & habitat != 'pool', 
                             .(running_id, date)][
                               .SD, on='running_id', x.date]]
  dt_list$fun_sedi[
    is.na(date),
    date := in_metadata_edna[sample_type=='sediment' & habitat != 'pool',
                             .(running_id, date)][
                               .SD, on='running_id', x.date]] 
  dt_list$fun_biof[
    is.na(date),
    date := in_metadata_edna[sample_type=='biofilm' & habitat != 'pool',
                             .(running_id, date)][
                               .SD, on='running_id', x.date]] 
  
  #Replace dates in dia_biof, which seem erroneous
  dt_list$dia_biof[
    ,
    date := in_metadata_edna[sample_type=='biofilm' & habitat != 'pool', 
                             .(running_id, date)][
                               .SD, on='running_id', x.date]
  ] 
  
  #Add Campaign and Site to bacteria data
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
  
  #Remove pools
  for(org_medium in c('dia_sedi', 'dia_biof', 'fun_sedi', 
                      'fun_biof', 'bac_sedi', 'bac_biof') 
  ){
    print(paste('Removing pools from', org_medium))
    nsplit <- str_split(org_medium, '_')[[1]]
    org <- nsplit[1]
    medium <- ifelse(nsplit[2]=='sedi', 'sediment', 'biofilm')
    new_org_medium <- paste0(org_medium, '_nopools')
    
    if (org=='bac') {
      metacols_toget <- c('running_id', 'date', 'habitat', 'country')
    } else {
      metacols_toget <- c('running_id', 'habitat')
    }
    
    dt_list[[org_medium]] <- merge(dt_list[[org_medium]], 
                                   in_metadata_edna[sample_type==medium & 
                                                      habitat != 'pool',
                                                    metacols_toget, with=F],
                                   by='running_id') 
    if (org=='bac') {
      dt_list[[org_medium]][country == 'Czech Republic', 
                            country := 'Czech']
    }
    
    if (org=='fun') {
      print('Correcting fungi dates')
      dt_list[[org_medium]][date == as.Date("2012-02-25"), 
                            date := as.Date("2021-02-25")]
    }
    
    dt_list[[new_org_medium]] <- dt_list[[org_medium]] %>%
      .[habitat!='pool',] %>%
      .[, `:=`(habitat = NULL,
               organism = new_org_medium)]
    dt_list[[org_medium]][, habitat := NULL]
  } %>% invisible
  
  if (!include_bacteria) {
    dt_list <- dt_list[!grepl('bac_', names(dt_list))]
  }
  
  return(dt_list)
}

#------ calc_spdiv -------------------------------------------------------------
# in_country <- 'Croatia'
# in_biodt <- tar_read(bio_dt)[['dia_sedi']][country == in_country,]
# in_metacols <- metacols
# level='local'

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

#' Calculate 'species' (taxonomic) diversity metrics
#'
#' Computes diversity metrics (richness, exponential Shannon, inverse Simpson) 
#' and partitions beta-diversity into turnover and nestedness at either the local or 
#' regional level.
#'
#' @param in_biodt A biodiversity `data.table` (species × samples).
#' @param in_metacols Character vector of metadata column names.
#' @param level Character. `"local"` (default) or `"regional"`.
#'
#' @return A `data.table` with richness, alpha/beta/gamma diversity, 
#'   and turnover/nestedness metrics depending on `level`.
#'
#' @details
#' - Local level: computes site-level richness, Shannon, Simpson, and 
#'   temporal beta-diversity.
#' - Regional level: uses `HierAnodiv` to partition gamma diversity 
#'   into spatial and temporal beta components.
#'
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
                        biodt_div[mean_richness>0, .(site, campaign)],
                        by=c('site', 'campaign'),
                        all.x = F
    )
    
    #Fourth-root transform data - DO NOT TRANSFORM DATA IN THE END
    # biodt_copy[, (spcols) := lapply(.SD, function(x) x^(1/4)), 
    #            .SDcols = spcols]
    #Example of publication that 4th-root transforms:
    #Garcıa-Roger et al. (2011). Aquatic Sciences 73, 567–579. doi:10.1007/S00027-011-0218-3
    #Elias et al. (2015). Marine and Freshwater Research, 2015, 66, 469–480 http://dx.doi.org/10.1071/MF13312
    
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
      expshannon = exp(vegan::diversity(as.matrix(.SD), index = "shannon")),
      invsimpson = vegan::diversity(as.matrix(.SD), index = "invsimpson")
    ), .SDcols = spcols] %>%
      .[is.infinite(invsimpson), invsimpson := 0]  %>%
      .[, `:=`(
        mean_expshannon = mean(expshannon),
        mean_invsimpson = mean(invsimpson)
      ), by=site]
    
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

#' Plot species richness against hydrological observations
#'
#' Produces ggplot-based visualizations comparing species richness 
#' with hydrological intermittence (percent flowing).
#'
#' @param in_sprich A `data.table` with species richness values.
#' @param in_envdt A `data.table` with environmental observations 
#'   including `state_of_flow`.
#'
#' @return Several ggplot objects (printed directly).
#'
#' @note
#' - Excludes macroinvertebrate subsets (flying/nonflying pools).
#' - Relies on consistent naming across `in_sprich` and `in_envdt`.
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

#' Subset river networks by DRN
#'
#' Creates country-level subsets of hydromod river networks, clipped to catchment
#' boundaries, and writes them to GeoPackage files.
#'
#' @param in_hydromod_paths_dt `data.table` with columns:
#'   - `country`: Country name
#'   - `network_path`: Path to full river network file
#'   - `catchment_path`: Path to catchment file
#' @param out_dir Directory where subset GeoPackages will be written.
#' @param overwrite Logical; whether to overwrite existing files.
#'
#' @return A named character vector of file paths to the subset network GeoPackages,
#'   named by `country`.
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

#' Clean river network
#'
#' Cleans a river network to ensure topological consistency, remove spurious
#' segments, and cluster nearby nodes. Steps include: splitting at confluences,
#' clustering nodes within a distance threshold, removing short dangling
#' segments, and re-assigning original IDs.
#'
#' @param rivnet_path Path to input river network file (GeoPackage).
#' @param idcol Column name for segment IDs.
#' @param node_clustering_dist Numeric; distance threshold for clustering
#'   nearby nodes (in CRS units, e.g. meters).
#' @param min_segment_length Minimum length (in CRS units) for dangling
#'   segments to be retained.
#' @param outdir Optional output directory.
#' @param save_gpkg Logical; whether to save cleaned network as GeoPackage.
#' @param return_path Logical; if TRUE, return file path instead of sf object.
#'
#' @return Cleaned river network as an `sf` object or file path (if `return_path=TRUE`).
#' @export
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
  # Remove duplicate/multiple edges and loops
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
#' Manual cleaning for Croatia river network
#'
#' Fixes geometrical issues with specific segments (hard-coded IDs) by replacing
#' geometries and splitting problem lines.
#'
#' @param rivnet `sf` object of the river network.
#'
#' @return A manually corrected `sf` object.
#' @export
manual_clean_croatia <- function(in_net) {
  # NOTE: Hard-coded UID values here are dataset-specific.
  
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
  
  # Combine corrected piece back with the rest
  out_net <- rbind(in_net, line_to_split)
  
  return(out_net)
}
#------ direct_network -----------------------------------------
# Define helper functions
get_endpoint <- function(line) st_coordinates(line)[nrow(st_coordinates(line)), ]
get_startpoint <- function(line) st_coordinates(line)[1, ]


#' Recursively ensures river segments are oriented downstream, starting from a
#' specified segment (typically the outlet). If a segment is reversed relative
#' to its downstream neighbor, its geometry is reversed.
#'
#' @param net `sf` object of river segments.
#' @param segid Integer; ID of the current segment to check.
#' @param check_ids Integer vector of segment IDs already checked (used to
#'   avoid infinite recursion).
#'
#' @return `sf` object with geometries reversed as needed.
#' @keywords internal
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



#' Direct entire river network
#'
#' Orients all river segments in a network to flow downstream, beginning from a
#' specified outlet segment. The function assumes the network has a single outlet
#' (a dangling segment at the river mouth).
#'
#' @param net `sf` object of river segments.
#' @param outlet_id Integer; ID of the outlet segment (`UID`).
#'
#' @return `sf` object with all segments oriented downstream.
#' @export

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

#' @title Fix complex confluences in a river network
#' @description This function takes a directed river network shapefile and 
#' fixes complex confluences by shifting nodes to resolve topological errors.
#' It relies on the `fix_confluences_inner` function, a helper function.
#' 
#' @param rivnet_path A character string specifying the path to the input river 
#'     network shapefile. The network must be directed and topologically correct 
#'     except for the complex confluences.
#' @param max_node_shift A numeric value specifying the maximum distance to 
#'     shift nodes in order to fix confluences. Defaults to 5.
#' @param outdir A character string specifying the output directory for the 
#'     fixed network file. If not specified, the output file will be written to 
#'     the same directory as the input.
#' @param out_path A character string specifying the full path for the output file.
#'     If provided, `outdir` is ignored.
#' @return A character string with the file path to the newly created fixed network shapefile.
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

#' @title Assign Strahler stream order
#' @description This function computes the Strahler stream order for a river 
#'    network, either from a shapefile path or a data.table. It uses an iterative 
#'    approach to propagate stream orders downstream.
#' @param in_rivnet A character string specifying the path to the river network 
#'     shapefile or a data.table object representing the network.
#' @param idcol A character string specifying the name of the unique identifier 
#'     column for each stream segment.
#' @param verbose A logical value. If `TRUE`, prints the number of stream sections 
#'     with unassigned Strahler orders at each iteration. Defaults to `FALSE`.
#' @return A data.table containing the stream network with computed Strahler 
#'     orders and other topological information.
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
# in_reaches_attri_dt <- tar_read(reaches_attri)[[in_country]]
# country <- in_country

#' @title Re-assign topologically logical identifiers
#' @description This function re-assigns a topologically/hydrologically logical 'cat' or 'reach_id' to each segment section in a river network.
#' The goal is to correct a common problem where these identifiers are incorrectly assigned across confluences in the source data.
#' The function relies on comparing the spatial network with an independent hydrological model network to make corrections.
#' @param rivnet_path A character string specifying the path to the input river network shapefile.
#' @param strahler_dt A data.table containing the Strahler order for each segment section, typically from the `assign_strahler_order` function.
#' @param in_reaches_hydromod_dt A data.table with topological data from the hydrological model.
#' @param in_reaches_attri_dt A data.table with additional attributes for the reaches from the hydrological model.
#' @param outdir A character string specifying the output directory.
#' @param country A character string specifying the country. This is used for a series of hard-coded manual corrections for specific countries.
#' @param in_ext A character string specifying the file extension for the output file (e.g., 'gpkg').
#' @return A character string with the file path to the newly created and corrected network shapefile.
reassign_netids <- function(rivnet_path, strahler_dt, 
                            in_reaches_hydromod_dt,
                            in_reaches_attri_dt,
                            outdir, country = NULL, in_ext='gpkg') {
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
  
  #Rename reaches attribute columns
  setDT(in_reaches_attri_dt) %>%
    setnames(c('to-reach', 'length', 'slope', 'upstream_area'),
             c('to_reach_datagouv', 'length_datagouv', 'slope_net', 'upstream_area_net')
    )
  
  #Merge all data
  out_rivnet <- merge(
    as.data.table(out_rivnet)[, -c('to_reach_shpcor', 'to_reach_hydromod'), with=F], 
    to_reach_shpcor_dt, by='to', all.x=T) %>%
    merge(reaches_hydromod_format[, .(ID_hydromod, to_reach_hydromod)],
          by.x='cat_cor', by.y='ID_hydromod', all.x=T) %>%
    .[, hydromod_shpcor_match := (to_reach_hydromod == to_reach_shpcor)] %>%
    merge(in_reaches_attri_dt, by.x='cat_cor', by.y='ID', all.x=T)
  
  #Check for discrepancies between old and new reach data... few
  check <- out_rivnet[length_hydromod != length_datagouv,]
  
  
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
        'to_reach_shpcor', 'to_reach_hydromod', 'upstream_area_net', 
        'hydromod_shpcor_match')],
    out_path)
  
  return(out_path)
}
#------ compute_hydrostats_drn -------------------------------------------------
# in_drn <- 'Croatia'
# varname <-  'isflowing' #qsim
# in_sites_dt <- tar_read(sites_dt)[country == in_drn,]
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[in_drn]]
# in_hydromod_drn <- tar_read_raw((paste0('hydromod_dt_', in_drn, '_', varname)))
# in_network_idcol = 'cat_cor'

#' @title Compute hydrological statistics for a DRN
#' @description This function computes hydrological statistics of intermittence 
#' and discharge for a given DRN.
#' It can compute statistics at the drying river network (DRN) scale or at the 
#' scale of individual sites.
#' @param in_network_path A character string specifying the path to the input network shapefile.
#' @param in_sites_dt A data.table containing site information, including the `reach_id` for each site.
#' @param varname A character string specifying the variable to compute statistics for. 
#'     Can be 'isflowing' for intermittence or 'qsim' for discharge.
#' @param in_hydromod_drn A list containing hydrological model data, 
#'     typically including a data.table of time series data (`data_all`) and a
#'     data.table with date formats (`dates_format`).
#' @param in_network_idcol A character string specifying the name of the column 
#'     in the network that serves as a unique identifier for reaches. Defaults to 'cat'.
#' @return A list or data.table containing the computed hydrological statistics.
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
      #Compute hydrological statistics at the DRN scale (Relative flowing length)
      q_stats$drn <- intermod_dt[, 
                                 compute_hydrostats_intermittence(
                                   in_hydromod_dt = .SD,
                                   in_sites_dt = in_sites_dt,
                                   scale = 'drn')$drn, 
                                 by = nsim] 
      
      #Compute hydrological statistics at the scale of individual sites
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
# in_country <- in_drn <- 'Czech'
# hydrostats <- tar_read_raw(paste0("hydrostats_", in_country, '_qsim'))
# in_bio_dt <- tar_read(bio_dt)

#' @title Subset hydrological statistics by date
#' @description This function filters a list of hydrological statistics (`hydrostats`)
#'      to only include data that falls within the date range of biological sampling.
#' @param hydrostats A list or data.table containing hydrological statistics. 
#'     It is expected to have a 'date' column.
#' @param in_bio_dt A list of data.tables, where each data.table contains 
#'     biological data and a 'date' column.
#' @return The `hydrostats` object, filtered to the date range of the biological data.
subset_hydrostats <- function(hydrostats, in_bio_dt) {
  unique_sampling_dates <- lapply(in_bio_dt, function(org_dt) {
    org_dt[, .(date)]
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

#' @title Format site data table
#' @description This function reads a site data table, performs country-specific
#'    corrections to site IDs and reach IDs, and fills in missing coordinates using 
#'    a separate environmental data table.
#' @param in_path A character string specifying the path to the raw site data file.
#' @param in_env_dt A data.table containing environmental data with site information.
#' @param in_country A character string specifying the country. Used to apply specific formatting rules.
#' @return A data.table with cleaned and formatted site information, including site ID, coordinates, and reach ID.
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

#' @title Create sites GeoPackage
#' @description This function creates a GeoPackage (`.gpkg`) file containing
#'     either site points or river reaches associated with specific countries.
#' @param in_hydromod_paths_dt A data.table with a list of paths to hydrological 
#'      model input/output data for each country.
#' @param in_sites_dt A data.table containing site information.
#' @param out_dir A character string specifying the output directory for the GeoPackage files.
#' @param geom A character string specifying the geometry type to create, either 'reaches' or 'points'.
#' @param in_network_path_list A list of character strings specifying the paths 
#'     to network shapefiles for each country. Required if `geom` is 'reaches'.
#' @param in_network_idcol A character string specifying the ID column for network reaches. Defaults to 'cat'.
#' @param overwrite A logical value. If `TRUE`, existing files will be overwritten. Defaults to `FALSE`.
#' @return A named vector of character strings, where names are countries and 
#'     values are the paths to the created GeoPackage files.
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

#' @title Snap river sites to a river network
#' @description This function takes a vector of site points and snaps them to 
#'     the nearest point on a corresponding river network.
#'     It handles projection differences and joins relevant attributes from 
#'     the network to the snapped points.
#' @param in_sites_path A character string or vector of character strings 
#'     specifying the path(s) to the input site points shapefile(s).
#' @param in_network_path A character string specifying the path to the input network shapefile.
#' @param out_snapped_sites_path A character string specifying the output path 
#'    for the snapped sites file. If `NULL`, a default path is created.
#' @param custom_proj A logical value. If `TRUE`, a custom projection is created 
#'     for distance calculations. Defaults to `FALSE`.
#' @param proj_back A logical value. If `TRUE`, the snapped points are re-projected 
#'     back to their original projection. Defaults to `FALSE`.
#' @param in_sites_unique_id A character string specifying the unique ID column 
#'     for the sites. Defaults to 'site'.
#' @param in_network_unique_id A character string specifying the unique ID column 
#'     for the network segments. Defaults to 'UID'.
#' @param in_sites_idcol_tomatch A character string specifying the ID column in 
#'     the sites data to match with the network. Defaults to 'reach_id'.
#' @param in_network_idcol_tomatch A character string specifying the ID column 
#'     in the network to match with the sites. Defaults to 'cat'.
#' @param overwrite A logical value. If `TRUE`, existing files will be overwritten. Defaults to `FALSE`.
#' @return A character string with the path to the newly created snapped sites file.
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
            attri_to_join = c(in_network_idcol_tomatch, in_network_unique_id,
                              'upstream_area_net')
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

#' @title Subset AMBER database and create country-specific files
#' @description This function reads the full AMBER database, filters it for 
#'     specified countries, and creates a GeoPackage (`.gpkg`) file of the barrier 
#'     locations for each country. The points are subsetted based on intersection
#'      with country-specific catchment boundaries.
#' @param amber_path A character string specifying the path to the full AMBER data file.
#' @param in_hydromod_paths_dt A data.table containing paths to hydrological 
#'     model data for each country, including catchment boundaries.
#' @param out_dir A character string specifying the output directory for the GeoPackage files.
#' @param overwrite A logical value. If `TRUE`, existing files will be overwritten. Defaults to `TRUE`.
#' @return A named list of character strings, where names are countries and 
#'      values are the paths to the created AMBER points GeoPackage files.
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

#' @title Snap barrier sites to a river network
#' @description This function takes a vector of barrier points and snaps them to 
#'     the nearest point on a corresponding river network.
#' @param in_sites_path A character string or vector of character strings 
#'     specifying the path(s) to the input barrier points file(s).
#' @param in_network_path A character string specifying the path to the input network shapefile.
#' @param in_sites_idcol A character string specifying the unique ID column for the sites.
#' @param attri_to_join A character vector of attribute names from the network 
#'     to be joined to the snapped points.
#' @param out_snapped_sites_path A character string specifying the output path 
#'     for the snapped sites file. If `NULL`, a default path is created.
#' @param custom_proj A logical value. If `TRUE`, a custom projection is created 
#'     for distance calculations. Defaults to `TRUE`.
#' @param overwrite A logical value. If `TRUE`, existing files will be overwritten. Defaults to `FALSE`.
#' @return A character string with the path to the newly created snapped sites file.
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
# in_hydromod_drn <- tar_read(hydromod_comb_hist)[[paste0(
#   "hydromod_dt_", in_country, '_isflowing')]]
# in_net_shp_path <- tar_read(network_ssnready_shp_list)[[in_drn]]

#' @title Prepare data for STcon analysis
#' @description This function takes hydrological and network data, builds a 
#'     directed graph of the river network, calculates a river distance matrix, 
#'     and formats flow intermittence data into a matrix suitable for STcon analysis.
#' @param in_hydromod_drn A list containing a data.table with hydrological data.
#' @param in_net_shp_path A character string specifying the path to the network shapefile.
#' @return A list containing the sites status matrix, network structure matrix, 
#'     river distance matrix, and a reference intermittence data.table.
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
  # Name the cat 11111 and it will be our last vertex and where the outlet will be directed.  
  nodes_cat <- rbind(net_dt, 
                     data.table(UID = 11111, from = outlet_to), fill=T) %>%
    setorder(from)
  
  # Assign vertex attributes according to the cat_shape, ordered following
  # the "from" value from the cat_shape, 
  # which is not exactly "from anymore" it is an ID of the vertex
  V(net_graph)$UID <- as.character(nodes_cat$UID)
  
  #Compute edge length for distance matrices
  E(net_graph)$weight <- as.integer(net_dt$length_m)
  
  #Compute distance matrix (in integer, lighter to handle -- one-meter difference does not matter)
  river_dist_mat <- igraph::distances(net_graph)
  
  #Format intermittence data
  setDT(in_hydromod_drn$data_all) %>%
    setnames('reach_id', 'cat')
  
  # Create the "End_point" site that will correspond to the 111111 in the flow_intermittence dataset
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

#' @title Compute rolling spatio-temporal connectivity
#' @description This function computes spatio-temporal connectivity (STcon) for 
#'     a given time window, iterating over a set of dates and hydrological simulations.
#' @param in_preformatted_data A list of preformatted data for STcon analysis, 
#'     typically the output of `prepare_data_for_STcon`.
#' @param ref A logical value. If `TRUE`, a reference dataset without intermittence 
#'     is used. Defaults to `FALSE`.
#' @param in_nsim A character or numeric vector of simulation numbers to be included.
#'      Use 'all' for all simulations, or `NULL` for a single simulation.
#' @param in_dates A data.table or data.frame with a 'date' column representing 
#'     the end dates of the rolling windows.
#' @param window An integer specifying the size of the rolling window in days.
#' @param output A character string specifying the desired output from the `compute_stcon` function.
#' @param direction A character string specifying the direction of flow ('undirected', 'in', 'out').
#' @param routing_mode A character string specifying the routing mode.
#' @param weighting A logical value. If `TRUE`, weighting is applied.
#' @param rounding_factor An integer for rounding the results.
#' @param verbose A logical value. If `TRUE`, prints progress messages. Defaults to `FALSE`.
#' @param ... Additional arguments to be passed to the `compute_stcon` function.
#' @return A nested list of STcon results, organized by simulation and date.

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

#' @title Post-process STcon results
#' @description This function takes the raw output of STcon computation, 
#'     standardizes the values, and formats the results into a long data table 
#'     and a list of matrices for easier use.
#' @param in_STcon A list of STcon results, typically the output of `compute_STcon_rolling`.
#' @param in_net_shp_path A character string specifying the path to the network shapefile.
#' @param standardize_STcon A logical value. If `TRUE`, STcon values are 
#'    standardized by a reference dataset. Defaults to `FALSE`.
#' @param in_STcon_ref A list of reference STcon results for standardization.
#'    Required if `standardize_STcon` is `TRUE`.
#' @return A list containing a long-format data table of STcon values (`STcon_dt`) 
#'    and a list of STcon matrices (`STcon_mat`).
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

#------ plot_STcon [NOT USED] ------------------------------------------------------------
# tar_load(STcon_directed_formatted)
# in_country <- 'France'
# in_STcon_list <- STcon_directed_formatted[['France']]
# in_date <-in_STcon_list$STcon_dt[1000, date]
# in_net_shp_path = tar_read(network_ssnready_shp_list)[[in_country]]

#' @title Plot STcon values on a river network
#' @description This function takes a formatted list of STcon results, 
#'     a specific date, and a window size, and generates a `ggplot` object to 
#'     visualize the STcon values on a river network.
#' @param in_STcon_list A list containing a data.table of STcon values, typically the output of `postprocess_STcon`.
#' @param in_date A Date object specifying the date to plot.
#' @param in_window An integer specifying the window size to plot. Defaults to `10`.
#' @param in_net_shp_path A character string specifying the path to the network shapefile.
#' @param reverse_weighted_stcon A logical value. If `TRUE`, 
#'     the STcon values are inverted for plotting (e.g., lower STcon corresponds to a darker color). 
#'     Defaults to `TRUE`.
#' @return A `ggplot` object representing the plot.
plot_STcon <- function(in_STcon_list, in_date, in_window=10, 
                       in_net_shp_path, reverse_weighted_stcon = TRUE) {
  
  stcon_sel <- in_STcon_list$STcon_dt[
    (date == in_date) & (variable == paste0('STcon_m', in_window)),]
  
  
  net_v_stcon <- in_net_shp_path %>%
    sf::st_read %>%
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
# in_country <- 'Spain'
# in_preformatted_data = tar_read(preformatted_data_STcon)[[in_country]]
# in_net_shp_path <- tar_read(network_ssnready_shp_list)[[in_country]]
# sites_status_matrix = in_preformatted_data$sites_status_matrix
# network_structure = in_preformatted_data$network_structure
# raw_dist_matrix <- in_preformatted_data$river_dist_mat
# routing_mode = 'in'

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
compute_Fdist <- function(sites_status_matrix, 
                          network_structure, 
                          routing_mode, 
                          raw_dist_matrix, 
                          in_net_shp_path) {
  
  #Get a matrix of adjacency relationship among segments that includes routing mode
  #(1 if connected, Inf if not connected -- depending on routing mode)
  routed_adjacency <- graph_from_adjacency_matrix(network_structure, 
                                                  mode = 'directed', 
                                                  diag = FALSE) %>%
    distances(mode = routing_mode, 
              algorithm = "unweighted")
  # set all reachable nodes to 1 for matrix multiplication later
  routed_adjacency[!is.infinite(routed_adjacency)] <- 1
  
  #Multiply the raw distance matrix by the routed adjacency to get directed distances
  #To turn on or off the "routing mode" parameter
  #When routing mode is "in", Inf values for sites downstream (inverse when mode is "out")
  #value of 0 for site to itself, and actual distance value for sites upstream (downstream if mode is out)
  directed_dist_matrix <- routed_adjacency*raw_dist_matrix
  
  #Convert self-distance to Inf to avoid taking it in account with column minimums
  sites_status_matrix[sites_status_matrix==0] <- Inf 
  
  #For each time step and reach, compute nearest wet site
  dist_to_nearest_wet <- sites_status_matrix[,{
    pair_dist_status <- directed_dist_matrix*as.numeric(as.matrix(.SD)) #Multiple distance by site status (Inf is dry, 1 is wet)
    pair_dist_status[is.na(pair_dist_status)] <- Inf # replace NA with Inf
    as.list(Rfast::colMins(pair_dist_status, value=T))  # get the minimum distance to a wet site for each column (site)
  }
  ,  by=date] %>% 
    setnames(names(sites_status_matrix))  # rename columns to match original status matrix
  
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
  
  # select valid IDs and map them to UIDs
  IDs_sel <- unique(dist_to_nearest_wet_melt$ID) %>%
    setdiff(outlet_from) 
  UIDs_to_assign <- net_dt[match(IDs_sel, net_dt$from), list(ID=IDs_sel, UID)]
  
  # merge the calculated distances with UID information
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

#' @title Compute rolling Fdist
#' @description This function computes rolling mean and maximum values of Fdist 
#'     (distance to nearest wet site) for a set of sites.
#' @param in_Fdist_dt A data.table containing Fdist values, with 'date' and 'UID' columns.
#' @param in_sites_dt A data.table of site information, used to filter for relevant UIDs.
#' @return A data.table with new columns for the rolling mean and max Fdist values.
compute_Fdist_rolling <- function(in_Fdist_dt, in_sites_dt) {
  
  # define rolling window sizes in days
  rollingstep_short <-  c(10, 30, 60, 90, 120, 180)
  rollingstep_long <- c(365, 365*5, 365*10)
  rollingstep <- c(rollingstep_short, rollingstep_long)
  
  #Compute mean Fdist within rolling window, only for sites
  # filtering data to only include sites present in `in_sites_dt`
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


#------ compile_hydrocon_sites_country -----------------------------------------------
# in_hydrostats_sub_comb <- tar_read(hydrostats_sites_tsub_comb)
# in_STcon_directed <- tar_read(STcon_directed_formatted)
# in_STcon_undirected <- tar_read(STcon_undirected_formatted)
# in_Fdist_directed <- tar_read(Fdist_directed)
# in_Fdist_undirected <- tar_read(Fdist_undirected)
# in_site_snapped_gpkg_list <- tar_read(site_snapped_gpkg_list)
# in_country <- 'Spain'

#' @title Compile hydrological and connectivity data for a country
#' @description This function merges various hydrological and connectivity metrics 
#'     (hydrological stats, STcon, and Fdist) for a specific country's sites.
#' @param in_hydrostats_sub_comb A list containing hydrological statistics for sites and river network segments.
#' @param in_STcon_directed A list of directed STcon data, formatted as an output of `postprocess_STcon`.
#' @param in_STcon_undirected A list of undirected STcon data.
#' @param in_Fdist_directed A list of directed Fdist data.
#' @param in_Fdist_undirected A list of undirected Fdist data.
#' @param in_site_snapped_gpkg_list A list of paths to snapped site geopackage files.
#' @param in_country A character string specifying the country to compile data for.
#' @return A data.table containing all merged hydrological and connectivity metrics.
compile_hydrocon_sites_country <- function(in_hydrostats_sub_comb, 
                                           in_STcon_directed,
                                           in_STcon_undirected, 
                                           in_Fdist_directed,
                                           in_Fdist_undirected,
                                           in_site_snapped_gpkg_list,
                                           in_country) {
  # load and format site data
  in_sites_dt <- as.data.table(vect(in_site_snapped_gpkg_list[[in_country]]))
  
  # load and format hydrological statistics for sites and DRNs
  hydrostats_isflowing_site <- in_hydrostats_sub_comb[[
    paste0('hydrostats_sites_tsub_', in_country, '_isflowing')]][['site']] %>%
    setDT
  
  hydrostats_isflowing_drn <- in_hydrostats_sub_comb[[
    paste0('hydrostats_sites_tsub_', in_country, '_isflowing')]][['drn']] %>%
    setDT
  
  hydrostats_qsim <- in_hydrostats_sub_comb[[
    paste0('hydrostats_sites_tsub_', in_country, '_qsim')]] %>%
    setDT
  
  # define columns to keep for merging
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
                      .(site, UID, upstream_area_net)], ., by='site') %>%
    merge(STcon_cast, by=c('date', 'UID'), all.x=T) %>%
    merge(Fdist_cast, by=c('date', 'UID'), all.x=T)
  
  # merge(STcon_cast, by=c('date', 'UID', 'nsim'), all.x=T) #with nsim
  
  return(out_dt)  
}

#------ summarize_sites_hydrocon --------------------------------------------
# in_hydrocon_compiled <- tar_read(hydrocon_compiled)


#' @title Summarize hydrological and connectivity data 
#' @description This function computes summary statistics (e.g., mean, max) 
#'     for various hydrological and connectivity metrics for each site, 
#'     across all dates.
#' @param in_hydrocon_compiled A data.table containing compiled hydrological and 
#'    connectivity data, typically the output of `compile_hydrocon_sites_country`.
#' @return A data.table with a single row per site, containing summarized statistics.
#' 
summarize_sites_hydrocon <- function(in_hydrocon_compiled
                                     #, date_range
) {
  hydrocon_summmarized <- in_hydrocon_compiled[
    #(date >= min(date_range)) &   (date <= max(date_range))
    , list( #`:=`
      DurD_samp = sum(isflowing==0)/.N, #Total number of no-flow days 
      DurD3650past = .SD[.N, DurD3650past]/3650,
      PDurD365past = .SD[.N, PDurD365past],
      FreD_samp = length(.SD[!is.na(noflow_period), unique(noflow_period)])*.N/365, #Total number of no-flow periods 
      FreD3650past = .SD[.N, FreD3650past]/10,
      PFreD365past = .SD[.N, PFreD365past],
      DurD_max_samp = nafill(as.integer(max(noflow_period_dur, na.rm=T)), fill=0), #Maximum duration of no-flow period
      DurD_avg_samp = nafill(.SD[!duplicated(noflow_period), mean(noflow_period_dur, na.rm=T)], fill=0), #AVerage duration of no-flow period (same as DurD/FreD)
      RelF_avg_samp = mean(relF),
      RelF_min_samp = min(relF),
      relF3650past = .SD[.N, relF3650past],
      qsim_avg_samp = mean(qsim),
      Pqsim_avg_samp = mean(Pqsim),
      PmeanQ10past_max_samp = max(PmeanQ10past),
      STcon_m10_directed_avg_samp = mean(STcon_m10_directed, na.rm=T),
      STcon_m10_undirected_avg_samp = mean(STcon_m10_undirected, na.rm=T),
      Fdist_mean_10past_directed_avg_samp = mean(Fdist_mean_10past_directed),
      Fdist_mean_10past_undirected_avg_samp = mean(Fdist_mean_10past_undirected)
      #PDurD_samp =
      #PDurdmax_samp = 
    ), by=.(site, UID)]
  
  return(hydrocon_summmarized)
}

#Summarized----
#Annual DurCV
#Annual FreCV
#Mean event duration 
#CV of event duration
#Average No-flow start date - Julian day of the first drying event
#SD No-flow start date - Julian day of the first drying event
#SD6 (Gallart et al. 2012)
#CV RelF
#CV mean distance

#------ summarize_network_hydrostats -------------------------------------------
# in_country = 'Croatia'
# in_hydromod = tar_read(hydromod_comb_hist)
# in_all_date_range = c(as.Date('1960-10-01', '%Y-%m-%d'),
#                       as.Date('2021-10-01', '%Y-%m-%d'))
# in_samp_date_range = tar_read(hydrocon_compiled)[, range(date)]

#' @title Summarize network hydrological statistics
#' @description This function takes hydrological model output and summarizes 
#' basic statistics for a river network across a defined date range. 
#' @param in_hydromod A list containing historical hydrological data for different countries.
#' @param in_all_date_range A vector of two dates defining the full historical period.
#' @param in_samp_date_range A vector of two dates defining the specific sampling period.
#' @param in_country A character string specifying the country for which to summarize data.
#' @return A data.table containing summarized hydrological statistics for each river reach.
summarize_network_hydrostats <- function(
    in_hydromod,
    in_all_date_range,
    in_samp_date_range,
    in_country) {
  
  #Extract discharge and isflowing data for the specified country
  qsim_country <- in_hydromod[[
    paste0('hydromod_hist_dt_', in_country, '_qsim')]]
  isflowing_country <- in_hydromod[[
    paste0('hydromod_hist_dt_', in_country, '_isflowing')]]
  
  #Link q data - keep only full hydrological years, 
  #Exclude 2022 because includes period after sampling
  
  #Calculate average discharge over the full historical period
  qsim_all_dt <- qsim_country$data_all[
    (date >= min(in_all_date_range)) &
      (date < max(in_all_date_range)),  
    list(qsim_avg = mean(qsim, na.rm=T)), 
    by=reach_id] 
  
  #Calculate average discharge over the sampling period
  qsim_samp_dt <- qsim_country$data_all[
    (date >= min(in_samp_date_range)) &
      (date <= max(in_samp_date_range)),  
    list(qsim_avg_samp = mean(qsim, na.rm=T)), 
    by=reach_id] 
  
  #Calculate flow duration (proportion of dry days) over the sampling period
  isflowing_samp_dt <- isflowing_country$data_all[
    (date >= min(in_samp_date_range)) &
      (date <= max(in_samp_date_range)),  
    list(DurD_samp = sum(isflowing==0)/.N), 
    by=reach_id] 
  
  #Merge all summary data tables and add country information
  hydrostats_net <- mergeDTlist(
    list(qsim_all_dt, qsim_samp_dt, isflowing_samp_dt),
    by='reach_id', set_suffix=F) %>%
    .[, country := in_country]
  
  return(hydrostats_net)
}

#------ summarize_env ---------------------------------------------------------
# in_env_dt <- tar_read(env_dt)

#' @title Summarize environmental data
#' @description This function takes raw environmental data and computes the mean 
#'     value for a set of variables, grouped by drainage basin, site, and stream type. 
#'     It also generates a boxplot for visualization.
#' @param in_env_dt A data.table containing raw environmental data with a 'state_of_flow' column.
#' @return A list containing two summarized data tables 
#'     (all flows, and only flowing sites) and a ggplot object of the boxplot.
summarize_env <- function(in_env_dt) {
  # dynamic_vars <- c(
  #   'avg_depth_macroinvertebrates',
  #   'avg_velocity_macroinvertebrates',
  #   'embeddedness',
  #   'bedrock',
  #   'boulders',
  #   'cobbles',
  #   'gravel',
  #   'sand',
  #   'particle size',
  #   'maximum depth',
  #   'oxygen_sat',
  #   'discharge_l_s',
  #   'filamentous algae',
  #   'macrophyte cover',
  #   'leaf_litter_cover',
  #   'moss_cover',
  #   'wood_cover',
  #   'average_wetted_width',
  #   'conductivity_micros_cm',
  #   'max_wetted_width_m',
  #   'min_wetted_width_m',
  #   'oxygen_mg_l',
  #   'ph',
  #   'temperature_c')
  
  # define grouping and exclusion columns
  group_cols <- c('drn', 'site', 'stream_type')
  exclude_cols <- c('if_ip_number_and_size_2_axes_+_depth_of_the_pools',
                    'campaign', 'date', 'state_of_flow')
  
  # identify data columns by excluding grouping and exclusion columns
  dat_cols <- setdiff(names(in_env_dt), 
                      c(group_cols, exclude_cols, 'running_id'))
  
  #str(in_env_dt[, dat_cols, with=F])
  
  # calculate mean for each data column, excluding dry sites
  env_summarized <- in_env_dt[state_of_flow != 'D', 
                              lapply(.SD, function(x) mean(x, na.rm=T)),
                              .SDcols = dat_cols,
                              by=group_cols]
  
  # calculate mean for each data column, only for 'F' (flowing) sites
  env_summarized_nopools <- in_env_dt[state_of_flow == 'F',
                                      lapply(.SD, function(x) mean(x, na.rm=T)),
                                      .SDcols = dat_cols,
                                      by=group_cols]
  
  # create a boxplot to visualize the summarized environmental data
  country_env_plot <- ggplot(melt(env_summarized, id.vars=group_cols),
                             aes(x=drn, y=value, fill=drn)) +
    geom_boxplot() +
    facet_wrap(~variable, scales='free')
  
  return(list(
    dt_all = env_summarized,
    dt_nopools = env_summarized_nopools,
    plot = country_env_plot
  ))
}


#------ summarize_drn_hydroproj_stats ------------------------------------------------
# ** 4 simulated variables:
#   * Discharge [m3/s]: Spatially distributed discharges at the reach level simulated at daily time step (outputs of the hydrological model JAMS-J2000)
# * Baseflow (groundwater contribution to the discharge) [m3/s]: Spatially distributed baseflows at the reach level simulated at daily time step (outputs of the hydrological model JAMS-J2000)
# * State of flow (binary variable: 0=dry, 1=flowing): Spatially distributed state of flow at the reach level simulated at daily time step (outputs of the flow intermittence random forest model)
# * Other hydroclimatic variables (temperature [°C], precipitation [mm], rainfall [mm], snowfall [mm], potential evapotranspiration [mm], actual evapotranspiration [mm], vegetation interception [mm], snow water equivalent [mm], saturation of the soil layer [0-1], saturation of the grounwater layer [0-1]): Spatially aggregated variables at the catchment scale simulated at daily time step (outputs of the hydrological model JAMS-J2000)
# 
# ** Characteristics of the projection simulations:
#   * 3 SSP scenarios: SSP1-2.6, SSP3-7.0, SSP5-8.5
# * Global Climate Models: gfdl-esm4, ipsl-cm6a-lr, mpi-esm1-2-hr, mri-esm2-0, ukesm1-0-ll
# * The analogue downscaling method produced 20-members ensembles for each combination of GCMs and SSP scenarios for the 6 studied catchments. As the uncertainty related to the downscaling method is rather negligeable, only the 10th member at the center of the distribution is given in this dataset (see Mimeau et al. 2024).
# * Reference period: 1985-2014 (data for the reference period are only given for SSP3-7.0)
# * Projection period: 2015-2100
# * JAMS-J2000 spatially distributed hydrological model: j2k_Guadiaro
# * Flow intermittence model: random_forest_Genal_2023-01-18_option0_15.RData
# 
# ** Name of the ncdf files:
#   {1}_{2}_{3}_{4}_{5}_{6}_{7}.nc
# 
# {1}: name of the DRN
# {2}: variable
# {3}: rprojection (in distinction with the reconstruction simulations)
# {4}: GCM
# {5}: SSP
# {6}: period of the dataset (reference period 1985-2014 or projection period 2015-2100)
# {7}: spatial aggregation (distributed or aggregated)

# hydroproj_path <- file.path(
#   wp1_data_gouv_dir,
#   "projections",
#   "Albarine_flowstate_projection_gfdl-esm4_ssp585_2015-2100_spatially-distributed.nc")
# varname = 'flowstate'

#' @title Summarize drainage basin hydrological projection stats
#' @description This function processes a NetCDF file containing hydrological 
#' projection data. It extracts metadata, reads flow state data, and computes a 
#' yearly summary of flow duration.
#' @param hydroproj_path A character string specifying the path to the NetCDF file.
#' @return A data.table containing summarized flow state statistics by year
#'  for each reach, along with associated metadata.
summarize_drn_hydroproj_stats <- function(hydroproj_path) {
  
  #Decompose name
  metadata_dt <- str_split(basename(hydroproj_path), '_')[[1]] %>%
    setNames(c('catchment', 'varname', 'time_period', 'gcm', 
               'scenario', 'date_range', 'file_end')) %>%
    as.list %>%
    data.frame %>%
    setDT %>%
    .[, path := hydroproj_path]
  
  nc <- nc_open(hydroproj_path) # open netcdf file
  reachID <- ncvar_get(nc, "reachID") # get list of reaches IDs
  dates <- ncvar_get(nc, "date") # get dates of simulation period
  dates <- as.Date(dates, origin="1950-01-01") # convert dates into R date format
  
  # get simulated variable from the NetCDF file
  hydro_dt <- get_nc_var_present(nc = nc, varname = metadata_dt$varname, # 0=dry, 1=flowing
                                 reachID = reachID, dates = dates,
                                 selected_sims = NULL) 
  setnames(hydro_dt$data_all, 'flowstate', 'isflowing', skip_absent = TRUE)
  metadata_dt[varname=='flowstate', varname:='isflowing']
  
  # qsim_samp_dt <- qsim_country$data_all[
  #   (date >= min(in_samp_date_range)) &
  #     (date <= max(in_samp_date_range)),  
  #   list(qsim_avg_samp = mean(qsim, na.rm=T)), 
  #   by=reach_id] 
  
  # compute yearly drying duration stats if the variable is 'isflowing'
  if (metadata_dt$varname == 'isflowing') {
    stats_dt <- hydro_dt$data_all[(date > as.Date('2014-12-31')) &
                                    (date < as.Date('2100-01-01')), 
                                  list(DurD_yr = sum(isflowing==0)/.N), 
                                  by=.(reach_id, format(date, '%Y'))] %>%
      setnames('format', 'year')
  }
  
  return(cbind(stats_dt, metadata_dt[, .(catchment, gcm, scenario)])
  )
}

#------ merge_allvars_sites ----------------------------------------------------
# in_country <- 'Spain'
# in_spdiv_local <- tar_read(spdiv_local)
# in_spdiv_drn <- NULL
# in_hydrocon_compiled <- tar_read(hydrocon_sites_compiled)
# in_hydrocon_summarized <- tar_read(hydrocon_sites_summarized)
# in_env_dt <- tar_read(env_dt)
# in_env_summarized <- tar_read(env_summarized)
# in_genal_upa = tar_read(genal_sites_upa_dt)

#' @title Merge all variables for sites
#' @description This function takes multiple data tables (biodiversity, hydrological connectivity, environmental variables) and merges them into two comprehensive data tables: one at the individual campaign level and one summarized by site.
#' @param in_spdiv_local A data table with biodiversity statistics.
#' @param in_spdiv_drn A data table with DRN biodiversity statistics.
#' @param in_hydrocon_compiled A data table with compiled hydrological and connectivity data by date and site.
#' @param in_hydrocon_summarized A data table with summarized hydrological and connectivity data by site.
#' @param in_env_dt A data table with sampling-time environmental data.
#' @param in_env_summarized A list containing summarized environmental data tables.
#' @param in_genal_upa A data table with upstream area data for a specific site in Genal.
#' @return A list containing the merged data tables and a list of column names.
merge_allvars_sites <- function(in_spdiv_local, in_spdiv_drn=NULL,
                                in_hydrocon_compiled, in_hydrocon_summarized,
                                in_env_dt, in_env_summarized,
                                in_genal_upa) {
  
  #Fill basin area NAs in environmental data for Genal basin in Spain
  #https://stackoverflow.com/questions/72940045/replace-na-in-a-table-with-values-in-column-with-another-table-by-conditions-in
  in_env_dt[is.na(basin_area_km2),
            basin_area_km2 :=  in_genal_upa[
              .SD, on='site', x.basin_area_km2]]
  
  in_env_summarized$dt_nopools[is.na(basin_area_km2),
                               basin_area_km2 :=  in_genal_upa[
                                 .SD, on='site', x.basin_area_km2]]
  
  #Compute state of flow in previous time step (could be inserted much before in the workflow later)
  in_env_dt[, state_of_flow_tm1 := lag(state_of_flow, n=1L), by=site]
  #Assume that the state of flow prior to first time step is the same as during the first time step
  in_env_dt[is.na(state_of_flow_tm1), state_of_flow_tm1 := state_of_flow]
  
  #Remove sites x organism combinations that were always dry
  dry_only_sites <- in_env_dt[, any(state_of_flow=='F'), by=site][V1==F, site]
  
  #Merge diversity data
  if (!is.null(in_spdiv_drn)) {
    setDT(in_spdiv_drn)
    spdiv <- merge(in_spdiv_local, in_spdiv_drn,
                   by=c('country', 'organism')) 
  } else {
    spdiv <- in_spdiv_local
  }
  
  spdiv <- spdiv[!(site=='GEN04' & campaign=='1'),]
  #Remove GEN04_1 from all organisms, no local environmental data, sampled only for eDNA. Too unsure.
  
  #Create "organism_class" column for labeling/coloring/merging
  spdiv[, organism_class := gsub('_[a-z]+', '', organism)]
  
  
  #List column names by originating dt
  group_cols = c("running_id", "site", "date", "campaign", "organism", 
                 "organism_class", "country", "UID", "upstream_area_net")
  exclude_cols = c("ncampaigns", "name", "isflowing", "reach_length",
                   "noflow_period", "noflow_period_dur", "last_noflowdate", "drn",
                   "if_ip_number_and_size_2_axes_+_depth_of_the_pools",
                   "latitude", "longitude", "reach_id", "hy", "month",
                   'min_wetted_width', 'left_river_bank_slope', 'right_river_bank_slope',
                   'qsim')
  
  dtcols <- list(
    div = setdiff(names(spdiv), 
                  c(group_cols, exclude_cols,
                    names(in_hydrocon_compiled), names(in_env_dt))),
    div_summarized = c("mean_richness", "mean_expshannon", "mean_invsimpson",
                       "JBDtotal", "JRepl", "JRichDif", "JRepl/BDtotal", "JRichDif/BDtotal",
                       "RBDtotal", "RRepl", "RRichDif", "RRepl/BDtotal", "RRichDif/BDtotal",
                       "Gamma", "Beta", "mAlpha"), 
    hydro_con = setdiff(names(in_hydrocon_compiled), 
                        c(group_cols, exclude_cols,
                          names(spdiv), names(in_env_dt))),
    hydro_con_summarized = setdiff(names(in_hydrocon_summarized), 
                                   c(group_cols, exclude_cols,
                                     names(spdiv), names(in_env_dt))), 
    env = setdiff(names(in_env_dt),
                  c(group_cols, exclude_cols,
                    names(spdiv), names(in_hydrocon_compiled))),
    env_summarized = setdiff(c(names(in_env_summarized$dt_nopools), 
                               'upstream_area_net'),
                             c(group_cols, exclude_cols)),
    group_cols =  group_cols,
    exclude_cols = exclude_cols
  )
  
  #Delineate numeric environmental columns
  dtcols <- c(dtcols, 
              env_num = list(intersect(dtcols$env,
                                       names(in_env_dt)[
                                         sapply(in_env_dt, class) 
                                         %in% c('numeric', 'integer')])
              ),
              env_summarized_num = list(intersect(
                dtcols$env_summarized,
                names(in_env_summarized$dt_nopools)[
                  sapply(in_env_summarized$dt_nopools, class) 
                  %in% c('numeric', 'integer')])
              )
  )
  
  #Compute average metric between sediment and biofilm for eDNA
  # avg_spdiv_edna <- spdiv[
  #   organism_class %in% c('dia', 'fun', 'bac'), 
  #   lapply(.SD, function(x) mean(x, na.rm=T)), 
  #   by = eval(names(spdiv)[names(spdiv) %in% 
  #                            setdiff(dtcols$group_cols, 'organism')]),
  #   .SDcols = setdiff(dtcols$div, 
  #                     c(dtcols$group_cols, dtcols$exclude_cols))
  # ] %>%
  #   .[, organism := organism_class]
  
  # bind the averaged eDNA data with the original biodiversity data
  # spdiv <- rbind(spdiv, avg_spdiv_edna, use.names=T, fill=T)
  
  
  #------------ Create data.table for individual site and campaigns ------------
  #Merge diversity metrics with hydro_con
  spdiv_hydro_con <- merge(spdiv, in_hydrocon_compiled,
                           by=c('date', 'site'), all.x=T) 
  
  #Merge environmental variables
  setDT(in_env_dt)
  all_vars_merged <- merge(spdiv_hydro_con, 
                           in_env_dt[, c(dtcols$env, 'site', 'campaign'), with=F],
                           by=c('site', 'campaign'), all.x=T) %>%
    .[!(site %in% dry_only_sites), ] %>% #Remove sites that never flowed 
    .[, country := factor(
      country,
      levels = c("Finland", "France",  "Hungary", "Czech", "Croatia", "Spain" ))]
  
  
  #------------ Create data.table summarized by site across campaigns ---------
  # summarize biodiversity data across all dates
  spdiv_summarized <- spdiv[!duplicated(paste(site, organism)), 
                            intersect(names(spdiv), 
                                      c(group_cols, dtcols$div_summarized)), 
                            with=F] %>%
    .[, `:=`(mean_richness = as.integer(round(mean_richness)))] #Round mean richness to an integer (for binomial modelling)
  
  all_vars_summarized <- merge(
    spdiv_summarized,
    in_hydrocon_summarized[, c(dtcols$hydro_con_summarized, 'site'),  with=F], 
    by='site', all.x=T) %>%
    merge(in_hydrocon_compiled[!duplicated(site), .(site, upstream_area_net)],
          by='site') %>%
    merge(in_env_summarized$dt_nopools[, c(dtcols$env_summarized, 'site'), with=F],
          by='site', all.x=T) %>%
    .[, c('campaign', 'date') := NULL] %>%
    .[!(site %in% dry_only_sites), ] %>% #Remove sites that never flowed 
    .[, country := factor(
      country,
      levels = c("Finland", "France",  "Hungary", "Czech", "Croatia", "Spain" ))]
  
  return(list(
    dt = all_vars_merged,
    dt_summarized = all_vars_summarized,
    cols = dtcols)
  )
}

#------ plot_edna_biof_vs_sedi -------------------------------------------------
# in_allvars_merged <- tar_read(allvars_merged)

#' @title Plot eDNA biofilm vs sediment
#' @description This function creates two plots comparing richness and Gamma 
#'     diversity between eDNA samples from biofilm and sediment, faceted by organism 
#'     class and country.
#' @param in_allvars_merged A list containing the merged data tables, specifically `in_allvars_merged$dt`.
#' @return A list containing two ggplot objects: `richness` and `site_gamma`.
plot_edna_biof_vs_sedi <- function(in_allvars_merged) {
  allvars_edna <- setDT(in_allvars_merged$dt)[organism_class %in% c('dia', 'fun', 'bac'),] 
  allvars_edna[, edna_source := gsub('^[a-z]+_', '', organism)]
  
  #Compare diversity metric by source of edna (biofilm vs sediment), organism and country
  alpha_index_list <- c('richness', 'expshannon', 'invsimpson')
  alpha_cast <- 
    melt(allvars_edna, 
         id.vars=c('country', 'site', 'date', 'organism', 'organism_class', 'edna_source'),
         measure.vars=alpha_index_list,
         variable.name='diversity_index') %>%
    dcast(country+site+date+organism_class+diversity_index~edna_source, 
          value.var = 'value')
  
  p_alpha_list <- lapply(alpha_index_list, function(index) {
    ggplot(alpha_cast[diversity_index==index,], aes(x=biof, y=sedi)) + 
      geom_point() + 
      geom_abline() + 
      ggtitle(label = index) +
      facet_wrap(organism_class~country, scales='free', nrow=3)
  }) %>% setNames(alpha_index_list)
  
  #Compare gamma by source of edna (biofilm vs sediment), organism and country
  malpha_index_list <- paste0('mean_', alpha_index_list)
  malpha_cast <- 
    melt(allvars_edna, 
         id.vars=c('country', 'site', 'date', 'organism', 'organism_class', 'edna_source'),
         measure.vars=malpha_index_list,
         variable.name='diversity_index') %>%
    dcast(country+site+date+organism_class+diversity_index~edna_source, 
          value.var = 'value')
  
  p_malpha_list <- lapply(malpha_index_list, function(index) {
    ggplot(malpha_cast[diversity_index==index,], aes(x=biof, y=sedi)) + 
      geom_point() + 
      geom_abline() + 
      ggtitle(label = index) +
      facet_wrap(organism_class~country, scales='free', nrow=3)
  }) %>% setNames(malpha_index_list)
  
  return(list(
    alpha=p_alpha_list,
    mean_alpha=p_malpha_list
  ))
}

#------ compute_cor_matrix -----------------------------------------------------
# in_allvars_merged <- tar_read(allvars_merged)

#' @title Compute grouped correlation matrices
#' @description This function computes correlation matrices for different sets of 
#'     variables (hydrology, environment, diversity) and different groupings 
#'     (overall, by organism, by country).
#' @param in_allvars_merged A list containing the data table (`dt`) and a list 
#'     of column names (`cols`) categorized by their origin.
#' @return A list containing five correlation matrices for different variable 
#'     combinations and groupings, plus the original column list.
compute_cor_matrix <- function(in_allvars_merged) {
  dt <- in_allvars_merged$dt
  cols_by_origin <- in_allvars_merged$cols
  
  #Pre-formatting: fill NA values and convert data types
  dt[is.na(noflow_period_dur), noflow_period_dur := 0]
  dt[, PrdD := as.numeric(PrdD)]
  
  # --- Calculate Correlations for each site and date ---
  # 1. hydro: correlations among hydrological variables
  cor_hydro <- compute_cor_matrix_inner(
    dt,
    x_cols = c(cols_by_origin$hydro_con),
    exclude_diagonal = FALSE) 
  
  # 2. Overall Correlation (Hydro + Env)
  cor_hydroenv <- compute_cor_matrix_inner(
    dt,
    x_cols = c(cols_by_origin$hydro_con, 
               cols_by_origin$env_num),
    exclude_diagonal = FALSE) 
  
  # 3. By Organism (Div x (Hydro + Env))
  cor_div <- compute_cor_matrix_inner(
    dt,
    group_vars = "organism",
    x_cols = cols_by_origin$div,
    y_cols = c(cols_by_origin$hydro_con, 
               cols_by_origin$env_num),
    exclude_diagonal = FALSE) 
  
  # 4. By Country (Hydro + Env)
  cor_hydroenv_bydrn <- compute_cor_matrix_inner(
    dt,
    group_vars = "country",
    x_cols = c(cols_by_origin$hydro_con, 
               cols_by_origin$env_num),
    exclude_diagonal = FALSE)  
  
  # 5. By Organism and Country (Div x (Hydro + Env))
  cor_div_bydrn <- compute_cor_matrix_inner(
    dt,
    group_vars = c("organism", "country"),
    x_cols = cols_by_origin$div,
    y_cols = c(cols_by_origin$hydro_con, 
               cols_by_origin$env_num),
    exclude_diagonal = FALSE) 
  
  return(list(hydro = cor_hydro,
              hydroenv = cor_hydroenv,
              div = cor_div,
              hydroenv_bydrn = cor_hydroenv_bydrn,
              div_bydrn = cor_div_bydrn,
              
              cols_by_origin = cols_by_origin
  ))
}

#------ compute_cor_matrix_summarized ------------------------------------------
# in_allvars_merged <- tar_read(allvars_merged)

#' @title Compute summarized correlation matrices
#' @description This function computes correlation matrices using site-summarized 
#'      data (across all sampling dates), analyzing relationships between 
#'      hydrological, environmental, and diversity metrics at a broader scale.
#' @param in_allvars_merged A list containing the summarized data table (`dt_summarized`) and column list (`cols`).
#' @return A list of correlation matrices for summarized data.
compute_cor_matrix_summarized <- function(in_allvars_merged) {
  dt <- in_allvars_merged$dt_summarized
  cols_by_origin <- in_allvars_merged$cols
  
  # --- Calculate Correlations for each site summarized ---
  # 1. Overall Correlation (Hydro)
  cor_hydro <- compute_cor_matrix_inner(
    dt,
    x_cols = c(cols_by_origin$hydro_con_summarized),
    exclude_diagonal = FALSE) 
  
  # 1. Overall Correlation (Hydro + Env)
  cor_hydroenv <- compute_cor_matrix_inner(
    dt,
    x_cols = c(cols_by_origin$hydro_con_summarized, 
               cols_by_origin$env_summarized_num),
    exclude_diagonal = FALSE) 
  
  # # 2. By Organism (Div x (Hydro + Env))
  cor_div <- compute_cor_matrix_inner(
    dt,
    group_vars = "organism",
    x_cols = cols_by_origin$div_summarized,
    y_cols = c(cols_by_origin$hydro_con_summarized, 
               cols_by_origin$env_summarized_num),
    exclude_diagonal = FALSE)
  
  # # 3. By Country (Hydro + Env)
  cor_hydroenv_bydrn <- compute_cor_matrix_inner(
    dt,
    group_vars = "country",
    x_cols = c(cols_by_origin$hydro_con_summarized, 
               cols_by_origin$env_summarized_num),
    exclude_diagonal = FALSE)
  
  # 4. By Organism and Country (Div x (Hydro + Env))
  cor_div_bydrn <- compute_cor_matrix_inner(
    dt,
    group_vars = c("organism", "country"),
    x_cols = cols_by_origin$div_summarized,
    y_cols = c(cols_by_origin$hydro_con_summarized, 
               cols_by_origin$env_summarized_num),
    exclude_diagonal = FALSE) 
  
  return(list(
    hydro = cor_hydro,
    hydroenv = cor_hydroenv,
    div = cor_div,
    hydroenv_bydrn = cor_hydroenv_bydrn,
    div_bydrn = cor_div_bydrn,
    cols_by_origin = cols_by_origin
  ))
}

#------ plot_cor_heatmaps -------------------------------------------------------
# in_cor_matrices <- tar_read(cor_matrices_list)
# p_threshold <- 0.05

#' @title Plot correlation heatmaps
#' @description This function generates a set of heatmaps to visualize the 
#'     correlation matrices computed previously. It can display correlations for
#'     different variable sets and across various groupings like organism and country.
#' @param in_cor_matrices A list of correlation matrices, typically the output 
#'    of `compute_cor_matrix` or `compute_cor_matrix_summarized`.
#' @param p_threshold A numeric value to filter correlations based on their p-value.
plot_cor_heatmaps <- function(in_cor_matrices,
                              p_threshold = 1) {
  
  create_correlation_heatmap <- function(cor_matrix, p_matrix, title,
                                         p_threshold = 1,
                                         is_square = TRUE, hc_order = TRUE) { # Add is_square argument
    
    # create a helper function to generate a single heatmap
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
  
  # get lists of unique organisms and countries for iteration
  org_list <- unique(in_cor_matrices$div$organism) %>%
    setdiff(c('miv','bac_biof', 'bac_sedi')) #Use no_pools instead
  country_list <- unique(in_cor_matrices$hydroenv_bydrn$country)
  
  # # --- 1. Hydro Correlation ---
  # Extract the correlation and p-value matrices for hydro
  cor_matrix_hydro <- dcast(in_cor_matrices$hydro, 
                            variable1 ~ variable2, 
                            value.var = "correlation")
  rnames <- cor_matrix_hydro$variable1
  cor_matrix_hydro <- as.matrix(cor_matrix_hydro[, -1])
  rownames(cor_matrix_hydro) <- rnames
  cor_matrix_hydro[is.na(cor_matrix_hydro)] <- 0
  
  p_matrix_hydro <- dcast(in_cor_matrices$hydro, 
                          variable1 ~ variable2, 
                          value.var = "p_value")
  p_matrix_hydro <- as.matrix(p_matrix_hydro[, -1])
  rownames(p_matrix_hydro) <- rnames
  p_matrix_hydro[is.na(p_matrix_hydro)] <- 1L
  
  # Create and print the heatmap
  hydro_heatmap <- create_correlation_heatmap(
    cor_matrix = cor_matrix_hydro, 
    p_matrix = p_matrix_hydro,
    title = "Hydro Correlations",
    p_threshold = p_threshold,
    is_square = TRUE)
  
  # # --- 2. Hydroenv Correlation ---
  # hydroenv_sub <- in_cor_matrices$hydroenv[p_value <= p_threshold & 
  #                                            correlation >= cor_threshold,]
  # 
  # Extract the correlation and p-value matrices for hydro env
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
  
  # --- 3. By Organism (Div x (Hydro + Env)) ---
  div_heatmaps <- lapply(org_list, function(org) {
    div_sub <- in_cor_matrices$div[organism == org,] %>%
      .[variable1 %in% in_cor_matrices$cols_by_origin$div,]
    
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
  
  # --- 4. Hydroenv Correlation by Country (Square) ---
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
  
  # --- 5. Div x Hydroenv Correlation by Country (Non-square) ---
  if (!is.null(in_cor_matrices$div_bydrn)) {
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
  } else {
    div_country_heatmaps <- NULL
  }
  
  
  return(list(
    hydro = hydro_heatmap,
    hydroenv = hydroenv_heatmap,
    div = div_heatmaps,
    hydroenv_country = hydroenv_country_heatmaps,
    div_country = div_country_heatmaps
  ))
}

#------ ordinate_local_env ----------------------------------------------------
# autoplot(local_env_pca$miv_nopools$pca, data = local_env_pca$miv_nopools$trans_dt, 
#          loadings = TRUE, loadings.colour = 'blue',
#          loadings.label = TRUE, loadings.label.size = 5) + theme_bw()

# in_allvars_dt <- tar_read(allvars_merged)$dt
# by_date=F

#' @title Ordinate local environmental data
#' @description This function performs Principal Component Analysis (PCA) on a 
#'      dataset of local environmental variables to create a reduced set of 
#'      dimensions (PCA axes) that capture the main sources of environmental variation.
#' @param in_allvars_dt A data table containing environmental data, which may be summarized by site.
#' @return A list containing the PCA results for macroinvertebrates and eDNA groups, 
#'     as well as a combined data table with the new PCA axes.
ordinate_local_env <- function(in_allvars_dt) {
  #1. Compute PCA for miv_nopools ----------------------------------------------
  # define environmental columns for macroinvertebrates
  env_cols_miv <- c('avg_velocity_macroinvertebrates', 'embeddedness',
                    'bedrock', 'particle_size', 'oxygen_sat', 'filamentous_algae',
                    'incrusted_algae', 'macrophyte_cover', 'leaf_litter_cover',
                    'moss_cover', 'wood_cover', 'riparian_cover_in_the_riparian_area',
                    'shade', 'hydromorphological_alteration', 'm2_biofilm',
                    'conductivity_micros_cm', 'ph', 'temperature_c')
  
  #'oxygen_mg_l' missing for entire country
  
  id_cols <- c('site', 'country')
  if ('date' %in% names(in_allvars_dt)) {
    id_cols <- c(id_cols, 'date')
  }
  
  #Convert all columns to numeric (rather than integer)
  in_allvars_dt[, (env_cols_miv) := lapply(.SD, as.numeric), 
                .SDcols = env_cols_miv] 
  
  #Fill NAs hierarchically. First by site, then by country, then overall
  miv_nopools_dt <- fill_nas_hierarchical(
    dt = in_allvars_dt[organism == 'miv_nopools'], 
    cols_to_fill = env_cols_miv, 
    site_col = 'site', 
    country_col = 'country')
  
  #Check distributions by country
  # dt_miv_envmelt <- melt(in_allvars_dt[organism == 'miv',], 
  #                        id.vars = c('country', 'site', 'campaign'),
  #                        measure.vars = env_cols_miv)
  # ggplot(dt_miv_envmelt, #[variable=='conductivity_micros_cm',],
  #        aes(x=country, y=value, color=country)) +
  #   geom_jitter() + 
  #   facet_wrap(~variable, scales='free_y') +
  #   scale_y_sqrt()
  # 
  # ggplot(dt_miv_envmelt, aes(x=value)) +
  #   geom_density() + 
  #   facet_wrap(~variable, scales='free')
  
  # perform PCA using a wrapper function
  out_list_miv <- trans_pca_wrapper(in_dt = miv_nopools_dt, 
                                    in_cols_to_ordinate = env_cols_miv, 
                                    id_cols = id_cols, 
                                    group_cols = NULL, 
                                    num_pca_axes = 4)
  
  # out_list_miv_country <- trans_pca_wrapper(in_dt = miv_nopools_dt, 
  #                                           in_cols_to_ordinate = env_cols_miv, 
  #                                           id_cols = c('site', 'date'), 
  #                                           group_cols = 'country', 
  #                                           num_pca_axes = 4)
  
  #2. Compute PCA for eDNA data ------------------------------------------------
  edna_orglist <- c('dia_sedi_nopools', 'dia_biof_nopools', 
                    'fun_sedi_nopools', 'fun_biof_nopools',
                    'bac_sedi_nopools', 'bac_biof_nopools')
  
  env_cols_edna <- c('filamentous_algae', 'incrusted_algae', 
                     'macrophyte_cover', 'leaf_litter_cover','moss_cover', 
                     'wood_cover', 'riparian_cover_in_the_riparian_area',
                     'shade', 'm2_biofilm', 'conductivity_micros_cm',
                     'oxygen_sat', 'ph','temperature_c')
  
  #Convert all columns to numeric (rather than integer)
  in_allvars_dt[, (env_cols_edna) := lapply(.SD, as.numeric), 
                .SDcols = env_cols_edna]
  
  #Check distributions by country
  # dt_edna_envmelt <- unique(in_allvars_dt[organism %in% edna_orglist,],
  #                           by=c('site', 'date')) %>%
  #   melt(id.vars = c('country', 'site', 'campaign', 'organism'),
  #        measure.vars = env_cols_edna)
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
  edna_dt <- unique(in_allvars_dt[organism %in% edna_orglist,],
                    by=setdiff(id_cols, 'country')) %>%
    .[, c(id_cols, env_cols_edna), with=F] %>%
    fill_nas_hierarchical(cols_to_fill = env_cols_edna, 
                          site_col = 'site', 
                          country_col = 'country') 
  
  out_list_edna <- trans_pca_wrapper(in_dt = edna_dt, 
                                     in_cols_to_ordinate = env_cols_edna, 
                                     id_cols = id_cols, 
                                     group_cols = NULL, 
                                     num_pca_axes = 4)
  
  #Merge dts
  out_dt_miv <- merge(
    in_allvars_dt[organism == 'miv_nopools', c(id_cols, 'organism'), with=F],
    out_list_miv$dt,
    by=id_cols)
  
  out_dt_edna <- merge(
    in_allvars_dt[organism %in% edna_orglist, c(id_cols, 'organism'), with=F],
    out_list_edna$dt,
    by=id_cols)
  
  
  out_dt <- rbind(out_dt_miv, out_dt_edna) %>%
    .[, organism_class := gsub('_[a-z]+', '', organism)] %>%
    .[, organism := NULL] %>%
    unique(by=c(id_cols, 'organism_class')) 
  
  return(list(
    miv_nopools = out_list_miv,
    edna = out_list_edna,
    dt_all = out_dt
  ))
}

#------ create_ssn_preds -------------------------------------------------------
# in_network_path = tar_read(network_ssnready_shp_list)
# in_hydrostats_net_hist = tar_read(hydrostats_net_hist)
# in_hydrostats_net_proj = tar_read(hydrostats_net_proj)

#' @title Create Spatial Stream Network prediction points
#' @description This function prepares spatial data for SSN modeling by creating 
#'     prediction points from a river network and associating them with historical 
#'     and projected hydrological data.
#' @param in_network_path A list of file paths to the river network shapefiles, separated by country.
#' @param in_hydrostats_net_hist A data frame with historical hydrological statistics for the river network.
#' @param in_hydrostats_net_proj A data frame with projected hydrological statistics for the river network.
#' @return A list containing two spatial objects: one for historical prediction 
#'      points and one for projected prediction points.
create_ssn_preds <- function(in_network_path,
                             in_hydrostats_net_hist,
                             in_hydrostats_net_proj) {
  
  # 1. Process the river network data
  # load and combine all country-specific river network shapefiles
  net_proj <- lapply(names(in_network_path), function(in_country) {
    #print(in_country)
    out_sf <- in_network_path[[in_country]] %>%
      st_cast("LINESTRING") %>%
      #Make sure that the geometry column is equally named regardless 
      #of file format (see https://github.com/r-spatial/sf/issues/719)
      st_set_geometry('geometry') %>%
      st_transform(3035) # transform to a common European projection (EPSG:3035)
    out_sf$country <- in_country
    return(out_sf)
  }) %>% do.call(rbind, .)
  
  # create prediction points at the centroid of each river segment
  net_centroids <- net_proj
  net_centroids$geometry <- st_line_sample(net_proj, sample=0.5) %>%
    st_cast('POINT')
  
  # 2. Merge with historical data
  net_predpts_hist <- merge(net_centroids, 
                            in_hydrostats_net_hist,
                            by.x = c('cat', 'country'), 
                            by.y = c('reach_id', 'country'),
                            all.x=T
  ) %>%
    rename(c('basin_area_km2'='upstream_area_net'))
  
  # 3. Merge with projected data
  in_hydrostats_net_proj[, proj_id := paste(gcm, scenario, year, sep='_')]
  proj_cast <- dcast(in_hydrostats_net_proj,
                     reach_id + country ~ proj_id,
                     value.var = 'DurD_yr')
  
  net_predpts_proj <- merge(
    net_centroids, proj_cast,
    by.x = c('cat', 'country'), 
    by.y = c('reach_id', 'country'),
    all.x=T) %>%
    rename(c('basin_area_km2'='upstream_area_net'))
  
  return(list(
    hist = net_predpts_hist,
    proj = net_predpts_proj
  ))
}

#------ create_ssn_europe ------------------------------------------------------
# in_network_path = tar_read(network_ssnready_shp_list)
# in_sites_path = tar_read(site_snapped_gpkg_list)
# in_barriers_path = tar_read(barrier_snapped_gpkg_list)
# in_allvars_dt= tar_read(allvars_merged)$dt
# in_local_env_pca = tar_read(local_env_pca)
# in_hydrostats_net_hist = tar_read(hydrostats_net_hist)
# in_pred_pts = NULL
# out_dir = file.path(resdir, 'ssn')
# out_ssn_name = 'ssn_eu'
# overwrite=T

# in_network_path = tar_read(network_ssnready_shp_list)
# in_sites_path = tar_read(site_snapped_gpkg_list)
# in_allvars_dt = tar_read(allvars_merged)$dt_summarized
# in_local_env_pca = tar_read(local_env_pca_summarized)
# in_barriers_path = tar_read(barrier_snapped_gpkg_list)
# in_hydromod = tar_read(hydromod_comb_hist)
# in_pred_pts = tar_read(ssn_pred_pts)
# out_dir = file.path(resdir, 'ssn')
# out_ssn_name = 'ssn_eu_summarized'
# overwrite = T

#' @title Create a European-scale SSN
#' @description This function assembles a full Spatial Stream Network (SSN) for 
#'      Europe by integrating the river network with observation sites, barriers, 
#'      and prediction points. It calculates key upstream metrics to support spatial modeling.
#' @param in_network_path A list of file paths to the river network shapefiles.
#' @param in_sites_path A list of file paths to the site observation points.
#' @param in_allvars_dt A data table containing all observation data, 
#'      including diversity and environmental metrics.
#' @param in_local_env_pca The output from the `ordinate_local_env` function, 
#'      containing PCA axes for environmental data.
#' @param in_barriers_path A list of file paths to dam and barrier locations.
#' @param in_hydrostats_net_hist Historical hydrological statistics for the network.
#' @param in_pred_pts An optional list of historical and projected prediction points from `create_ssn_preds`.
#' @param out_dir The directory to save the output SSN files.
#' @param out_ssn_name The base name for the output SSN file.
#' @param overwrite A logical value indicating whether to overwrite existing files.
#' @return A list of SSN objects, one for each organism.
create_ssn_europe <- function(in_network_path,
                              in_sites_path,
                              in_allvars_dt,
                              in_local_env_pca,
                              in_barriers_path,
                              in_hydrostats_net_hist,
                              in_pred_pts = NULL,
                              out_dir,
                              out_ssn_name,
                              overwrite = T) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  
  id_cols <- c('site', 'country')
  if ('date' %in% names(in_allvars_dt)) {
    id_cols <- c(id_cols, 'date')
  }
  
  lsn_path <- file.path(out_dir,
                        paste0(out_ssn_name, '_lsn')
  )
  
  #1. Build landscape network (lsn) -----------------------------------------------
  #Read input network
  net_eu <- lapply(names(in_network_path), function(in_country) {
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
    overwrite = overwrite
  )
  
  #Incorporate sites into the landscape network --------------------------------
  #2. Incorporate observation sites-----
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
    overwrite = overwrite
  )
  
  setDT(in_allvars_dt)[country == 'Czech Republic', country := 'Czech']
  
  # merge site data with all variables and PCA results
  sites_lsn_attri <- merge(sites_lsn,
                           in_allvars_dt, 
                           by=c('country', 'site', 'upstream_area_net')) %>%
    merge(in_local_env_pca$dt_all,
          by=c(id_cols, 'organism_class'))
  
  sites_list <- list(sites = sites_lsn_attri)
  
  #3. Incorporate prediction "sites" -----
  if (!(is.null(in_pred_pts))) {
    #Add historical prediction sites
    if (!(crs(in_pred_pts$hist) == crs(sites_eu))) {
      in_pred_pts$hist <- st_transform(in_pred_pts$hist, 3035) 
    }
    
    preds_hist_lsn <- SSNbler::sites_to_lsn(
      sites = in_pred_pts$hist,
      edges =  edges_lsn,
      lsn_path = lsn_path,
      file_name = "preds_hist",
      snap_tolerance = 5,
      save_local = TRUE,
      overwrite = overwrite
    )
    
    sites_list$preds_hist <- preds_hist_lsn
    
    #Add future prediction "sites"
    if (!(crs(in_pred_pts$proj) == crs(sites_eu))) {
      in_pred_pts$proj <- st_transform(in_pred_pts$proj, 3035) 
    }
    
    preds_proj_lsn <- SSNbler::sites_to_lsn(
      sites = in_pred_pts$proj,
      edges =  edges_lsn,
      lsn_path = lsn_path,
      file_name = "preds_proj",
      snap_tolerance = 5,
      save_local = TRUE,
      overwrite = overwrite
    )
    
    sites_list$preds_proj <- preds_proj_lsn
  }
  
  #4. Incorporate barriers into the landscape network -----
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
      overwrite = overwrite
    )
    
    sites_list$barriers <- barriers_lsn
  }
  
  # 5. Calculate upstream distance -----
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
  
  # 6. Compute segment Proportional Influence (PI) and Additive Function Values (AFVs) ----
  if (min(net_eu$qsim_avg) > 0) {
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
    
  } else {
    stop("Trying to use mean discharge to compute Additive Function Values (AFVs),
         but there are 0s in the discharge column.")
  }
  
  
  #7. Assemble an SSN for each organism  -----
  #(so that it only includes data for the  corresponding sites and dates)
  out_ssn_list <- lapply(unique(sites_list_lsn$sites$organism), function(in_org) {
    out_ssn_path <- file.path(out_dir, paste0(out_ssn_name, '_', in_org, '.ssn'))
    
    out_ssn <- ssn_assemble(
      edges = edges_lsn,
      lsn_path = lsn_path,
      obs_sites = sites_list_lsn$sites[sites_list_lsn$sites$organism == in_org,],
      preds_list = if (!is.null(in_pred_pts)) {
        sites_list_lsn[c("preds_hist", "preds_proj")]
      } else {NULL},
      ssn_path = out_ssn_path,
      import = TRUE,
      check = TRUE,
      afv_col = "afv_qsqrt",
      overwrite = overwrite
    )
    
    return(list(
      path = out_ssn_path,
      ssn = out_ssn
    ))
  }) %>% setNames(unique(sites_list_lsn$sites$organism))
  
  #Return --------------
  return(out_ssn_list)
}

#------ map_ssn_util -----------------------------------------------------
#Examples in https://cran.r-project.org/web/packages/SSNbler/vignettes/introduction.html
# in_ssn <- tar_read(ssn_eu_summarized)$miv_nopools$ssn

# 
# in_country='Hungary'
# color_col='upstream_area_net'
# linewidth_col='qsim_avg'
# in_edges <- in_ssn$edges[in_ssn$edges$country==in_country,]
# in_pts <- in_ssn$obs[in_ssn$obs$country==in_country,]

#' @title Plot a single SSN map
#' @description This utility function creates a single ggplot map of an SSN, 
#'     allowing for visualization of network edges and points. 
#'     It supports customization of line width, line color, point color, 
#'     and point shape based on specified data columns.
#' @param in_ssn A single SSN object.
#' @param in_edges The edges (river segments) of the SSN to be plotted.
#' @param linewidth_col The name of the column to control river segment line width.
#' @param linecolor_col An optional column for coloring river segments.
#' @param linecolor_lims Optional limits for the line color scale.
#' @param in_pts Optional points (e.g., sites) to be plotted on the network.
#' @param shape_col An optional column for controlling point shape.
#' @param ptcolor_col An optional column for coloring points.
#' @param ptcolor_lims Optional limits for the point color scale.
#' @return A ggplot object representing the SSN map.
map_ssn_util <- function(in_ssn, 
                         in_edges,
                         linewidth_col='qsim_avg',                                
                         linecolor_col=NULL, linecolor_lims=NULL,
                         in_pts=NULL, 
                         shape_col=NULL,
                         ptcolor_col=NULL, ptcolor_lims=NULL
) {
  
  # Compute limits if not provided
  if (is.null(ptcolor_lims) & !is.null(in_pts)) {
    ptcolor_lims <- range(in_pts[[ptcolor_col]])
  }
  
  if (is.null(linecolor_lims) & !is.null(linecolor_col)) {
    linecolor_lims <- range(in_edges[[linecolor_col]])
  }
  
  # Base plot with in_edges
  out_map <- ggplot() +
    geom_sf(data = in_edges,
            aes(
              linewidth = !!sym(linewidth_col)
            ),
            color = 'grey'
    ) + 
    scale_linewidth(
      name = if (linewidth_col=='qsim_avg') {
        expression('Simulated discharge'~m^3 * s^{-1})  
      } else {
        str_to_sentence(gsub('_', ' ', linewidth_col))
      },
      transform = 'sqrt',
      limits = range(in_ssn$edges[[linewidth_col]]),
      range = c(1, 2.5)
    ) +
    ggspatial::annotation_scale(location = "br", style='ticks') +
    theme_classic() +
    theme(axis.text = element_blank())
  
  # add line color if specified
  if (!is.null(linecolor_col)) {
    out_map <- out_map +
      geom_sf(
        data = in_edges,
        aes(
          linewidth = !!sym(linewidth_col),
          color = !!sym(linecolor_col)
        )
      ) 
    
    if (length(setdiff(sign(linecolor_lims), 0)) == 2) {
      out_map <- out_map +
        scale_color_fermenter(
          name = str_to_sentence(gsub('_', ' ', linecolor_col)),
          limits = linecolor_lims,
          n.breaks = 10,
          palette = "Spectral"
        )
    } else {
      out_map <- out_map +
        scale_color_viridis_b(
          name = str_to_sentence(gsub('_', ' ', linecolor_col)),
          limits = linecolor_lims,
          n.breaks=5) 
    }
    
  }
  
  # add points if specified
  if (!is.null(in_pts)) {
    out_map <- out_map +    
      geom_sf(
        data = in_pts,
        aes(color = get(as.character(ptcolor_col)), 
            shape= if (!is.null(shape_col)) {get(shape_col)} else {'16'}
        ),
        size = 3
      )
    
    if (!is.null(linecolor_col)) {
      out_map <- out_map +  new_scale_color()
    }
    out_map <- out_map + 
      scale_color_viridis_b(
        name = str_to_sentence(gsub('_', ' ', ptcolor_col)),
        limits = ptcolor_lims,
        n.breaks=5) 
  }
  
  # handle point shapes
  if (!is.null(shape_col) && (shape_col == 'stream_type')) {
    out_map <- out_map + 
      scale_shape(name='Stream type',
                  labels = c('Perennial', 'Non-perennial'))
  }
  
  if (is.null(shape_col)) {
    out_map <- out_map +
      scale_shape(guide="none")
  }
  
  return(out_map)
}

#------ pad_ssn_map ------------------------------------------------------------
#SHould generalize: useless to use edges and pts, could just be a list of sf
#then lapplied

#' @title Pad SSN map to a square bounding box
#' @description A utility function that calculates a square bounding box for an SSN map. This ensures consistent aspect ratios across multiple faceted plots, even if the geographic extents of the data differ.
#' @param in_edges The river network edges (an `sf` object).
#' @param in_pts The observation points (an `sf` object).
#' @return A list with `xlim` and `ylim` vectors for the new square plot limits.
pad_ssn_map <- function(in_edges=NULL, in_pts=NULL) {
  #Get bbox for edges + obs
  bb_edges <- if (!is.null(in_edges)) {sf::st_bbox(in_edges)}
  bb_obs   <- if (!is.null(in_pts)) {sf::st_bbox(in_pts)}
  
  # combine to get a single bounding box
  bb <- c(
    xmin = min(bb_edges["xmin"], bb_obs["xmin"]),
    ymin = min(bb_edges["ymin"], bb_obs["ymin"]),
    xmax = max(bb_edges["xmax"], bb_obs["xmax"]),
    ymax = max(bb_edges["ymax"], bb_obs["ymax"])
  )
  
  #Expand to square
  x_range <- bb["xmax"] - bb["xmin"]
  y_range <- bb["ymax"] - bb["ymin"]
  max_range <- max(x_range, y_range)
  
  # center the new bounding box
  x_mid <- (bb["xmax"] + bb["xmin"]) / 2
  y_mid <- (bb["ymax"] + bb["ymin"]) / 2
  
  # new limits
  return(list(
    xlim = c(x_mid - max_range / 2, x_mid + max_range / 2),
    ylim = c(y_mid - max_range / 2, y_mid + max_range / 2)
  ))
}

#------ map_ssn_facets-----------------------------------------------------------
# color_col='mean_richness'
# linewidth_col='qsim_avg'
# facet_col <- 'country'
# page_title <- 'Macroinvertebrates - Mean richness'

# in_ssn_summarized <- tar_read(ssn_eu_summarized)
# in_organism = 'miv_nopools'
# in_color_col <- 'upstream_area_net'
# in_ssn = in_ssn_summarized[[in_organism]]$ssn
# color_col = as.character(in_color_col)
# facet_col = 'country'
# linewidth_col = 'qsim_avg'
# page_title = paste(in_organism, in_color_col, sep=' - ')


# in_ssn = in_ssn_summarized[[1]]$ssn
# in_ptcolor_col = 'DurD_samp'
# in_pts = 'obs'
# ptcolor_col = as.character(in_ptcolor_col)
# facet_col = 'country'
# linewidth_col = 'qsim_avg'
# page_title = in_ptcolor_col
# linecolor_col=NULL
# shape_col=NULL

#' @title Create a faceted SSN map
#' @description This function generates a multi-panel map of an SSN, with each 
#'     panel representing a facet category (e.g., country). It ensures a consistent,
#'     square aspect ratio across all panels for visual comparability.
#' @param in_ssn A single SSN object.
#' @param facet_col The column name to be used for faceting.
#' @param in_pts An optional name of the point set to plot (e.g., 'obs').
#' @param linewidth_col The column name for controlling river segment line width.
#' @param linecolor_col An optional column for coloring river segments.
#' @param shape_col An optional column for controlling point shape.
#' @param ptcolor_col An optional column for coloring points.
#' @param page_title An overall title for the plot.
#' @return A patchwork object containing the faceted maps.
map_ssn_facets <- function(in_ssn, 
                           facet_col,
                           in_pts = NULL,
                           linewidth_col='qsim_avg',                                
                           linecolor_col=NULL, 
                           shape_col=NULL,
                           ptcolor_col=NULL,
                           page_title=NULL) {
  
  assert_that((facet_col %in% names(in_ssn$edges)) &
                (facet_col %in% names(in_ssn$obs)),
              msg = paste0(facet_col, ' not in in_ssn edges or obs'))
  
  #Define the different facet values to iterate over
  #Order countries by level of intermittence (RelF)
  if (facet_col == 'country') {
    in_ssn$edges[[facet_col]] <- factor(
      in_ssn$edges[[facet_col]],
      levels = c("Finland", "France",  "Hungary", "Czech", "Croatia", "Spain" ),
      ordered = T)
    facet_vals <- levels(in_ssn$edges[[facet_col]])
  } else {
    facet_vals <- unique(in_ssn$edges[[facet_col]])
  }
  
  #Define global limits for the color scale
  if (is.null(in_pts)) {
    pts <- NULL
    ptcolor_lims <- NULL
  } else if (in_pts == 'obs') {
    pts <- in_ssn$obs 
    ptcolor_lims <- range(pts[[ptcolor_col]], na.rm=T)
  } else {
    pts <- in_pts
    ptcolor_lims <- range(pts[[ptcolor_col]], na.rm=T)
  }
  
  #Define global limits for the color scale
  if (is.null(linecolor_col)) {
    linecolor_lims <- NULL
  } else {
    linecolor_lims <- range(in_ssn$edges[[linecolor_col]], na.rm=T)
  } 
  
  map_list <- lapply(facet_vals, function(facet_i) {
    #Subset edges and observations for given facet
    edges_i <- in_ssn$edges[in_ssn$edges[[facet_col]] == facet_i, ]
    
    pts_i <- pts[pts[[facet_col]] == facet_i, ]
    
    #Define cartographic extent to make sure the map is square to 
    #allow uniform faceting (despite uneven cartographic scale among facets)
    map_lims <- pad_ssn_map(in_edges = edges_i, 
                            in_pts = pts_i)
    
    #Map individual facet
    map_ssn_util(
      in_ssn = in_ssn, 
      in_edges = edges_i,
      linewidth_col = linewidth_col,                                
      linecolor_col = linecolor_col,
      linecolor_lims = linecolor_lims,
      in_pts = pts_i, 
      shape_col = shape_col,
      ptcolor_col = ptcolor_col, 
      ptcolor_lims = ptcolor_lims
    ) +
      ggtitle(facet_i) +
      coord_sf(xlim = map_lims$xlim, 
               ylim = map_lims$ylim) +
      theme(aspect.ratio = 1)
  })
  
  # combine the maps
  map_patchwork <- patchwork::wrap_plots(map_list, ncol=2) +
    plot_layout(guides = 'collect') +
    plot_annotation(page_title)
  
  return(map_patchwork)
}


#------ map_ssn_summarized -----------------------------------------------------
# in_ssn_summarized <- tar_read(ssn_eu_summarized)
# in_allvars_merged <- tar_read(allvars_merged)
# in_organism_dt = tar_read(organism_dt)
# verbose = T

#' @title Map summarized SSN variables
#' @description Generates a series of faceted maps visualizing diversity and 
#'      physical variables across the European SSN. It creates separate plots 
#'      for each variable, faceted by country.
#' @param in_ssn_summarized A list of SSN objects, typically for different organisms.
#' @param in_allvars_merged A list containing a data table of all variables and column names.
#' @param in_organism_dt A data table mapping organism codes to their labels.
#' @param verbose A logical value to indicate whether to print progress messages.
#' @return A named list of ggplot objects, combining maps of diversity and physical variables.
map_ssn_summarized <- function(in_ssn_summarized,
                               in_allvars_merged,
                               in_organism_dt,
                               verbose = T) {
  
  #Map every diversity variable for each organism
  # create a data frame of all organism-diversity variable combinations
  div_map_params <- expand.grid(intersect(names(in_ssn_summarized), 
                                          in_organism_dt$organism),
                                in_allvars_merged$cols$div_summarized,
                                stringsAsFactors = FALSE) %>%
    setDT %>%
    setnames(c('organism', 'divcol')) 
  
  maps_div <- mapply(
    function(in_organism, in_ptcolor_col) {
      if (verbose) {print(paste('Mapping', in_organism, in_ptcolor_col))}
      
      map_ssn_facets(
        in_ssn = in_ssn_summarized[[in_organism]]$ssn, 
        in_pts = 'obs',
        ptcolor_col = in_ptcolor_col, 
        facet_col = 'country',
        linewidth_col = 'qsim_avg',
        shape_col = 'stream_type',
        page_title = in_organism_dt[organism==in_organism, organism_label]
      )
    }
    ,
    in_organism = div_map_params$organism,
    in_ptcolor_col = div_map_params$divcol
  ) %>% setNames(div_map_params[, paste0(organism, '_', divcol)])
  
  
  #Plot every physical variable
  physvars <- c(in_allvars_merged$cols$hydro_con_summarized,
                in_allvars_merged$cols$env_summarized_num,
                'upstream_area_net')
  maps_physvars <- lapply(physvars, function(in_ptcolor_col) {
    if (verbose) {print(paste('Mapping', in_ptcolor_col))}
    
    map_ssn_facets(in_ssn = in_ssn_summarized[[1]]$ssn, 
                   in_pts = 'obs',
                   ptcolor_col = as.character(in_ptcolor_col), 
                   facet_col = 'country',
                   linewidth_col = 'qsim_avg',
                   page_title = in_ptcolor_col
    )
  }) %>% setNames(physvars)
  
  return(c(maps_div, maps_physvars))
}

#------ save_ssn_summarized_maps -----------------------------------------------
# in_ssn_summarized_maps = tar_read(ssn_summarized_maps)
# out_dir = figdir

#' @title Save SSN summary maps
#' @description Saves a pre-defined list of key SSN maps to a specified output
#'  directory as high-resolution PNG files.
#' @param in_ssn_summarized_maps A named list of ggplot objects to be saved.
#' @param out_dir The output directory for the saved files.
#' @return A list of file paths to the saved PNG files.
save_ssn_summarized_maps <- function(in_ssn_summarized_maps,
                                     out_dir) {
  
  pname_list <- c('miv_nopools_mean_richness',
                  'dia_mean_richness',
                  'fun_mean_richness',
                  'bac_mean_richness',
                  'DurD3650past',
                  'FreD3650past',
                  'STcon_m10_directed_avg_samp',
                  'STcon_m10_undirected_avg_samp',
                  'Fdist_mean_10past_directed_avg_samp',
                  'Fdist_mean_10past_undirected_avg_samp'
  )
  
  lapply(pname_list, function(pname) {
    out_filename <- file.path(out_dir, paste0('map_ssn_',  pname, '.png'))
    
    ggsave(
      filename = out_filename,
      plot = in_ssn_summarized_maps[[pname]],
      width = 9,
      height = 9,
      units='in',
      dpi=600
    )
    
    return(out_filename)
  })
}

#------ plot_drn_hydrodiv ------------------------------------------------------
# in_sites_dt <- tar_read(sites_dt)
# in_hydrocon_compiled <- tar_read(hydrocon_sites_compiled)
# in_allvars_dt <- tar_read(allvars_merged)$dt
# in_drn_dt = drn_dt
# out_dir = file.path(resdir, 'figures')

#' @title Plot stream intermittence and alpha diversity
#' @description Creates plots visualizing the relationship between the proportion
#' of flowing network length and alpha diversity over time, faceted by country.
#' @param in_hydrocon_compiled A data.table of compiled hydrological and connectivity data.
#' @param in_sites_dt A data.table with site information.
#' @param in_allvars_dt A data.table of all merged variables, including diversity and environmental data.
#' @param in_organism_dt A data.table mapping organism codes to their labels.
#' @param alpha_var Character string: name of the alpha diversity column to plot (e.g., "richness", "shannon").
#' @param write_plots Logical: whether to save the plots as PNG files.
#' @param out_dir Output directory for saved plots.
#' @return A list of ggplot objects: hydrograph, hydrograph with sampling dates, and alpha diversity plots per organism.
plot_drn_hydrodiv <- function(in_hydrocon_compiled,
                              in_sites_dt,
                              in_allvars_dt,
                              in_organism_dt,
                              alpha_var = "richness",
                              write_plots = FALSE,
                              out_dir) {
  hydrocon_compiled_country <- merge(
    in_hydrocon_compiled, 
    in_sites_dt[, .(site, country)], 
    by = 'site'
  )
  
  sampling_date_lims <- in_allvars_dt[, .(
    max_date = max(date, na.rm = TRUE),
    min_date = min(date, na.rm = TRUE)
  ), by = country]
  
  hydrocon_compiled_country[, country := factor(
    country,
    levels = c("Finland", "France", "Hungary", "Czech", "Croatia", "Spain"),
    ordered = TRUE
  )]
  
  hydrocon_compiled_country[, relD90past := 1 - relF90past]
  relF_melt <- hydrocon_compiled_country[!duplicated(paste(country, date)), ] %>%
    melt(id.vars = c('country', 'date'), measure.vars = c('relF90past', 'relD90past'))
  
  # 1. Plot mean proportion of flowing river network length over time
  plot_relF90past <- ggplot(relF_melt) +
    geom_area(aes(x = date, y = value, fill = variable), alpha = 0.7) +
    scale_fill_manual(
      name = 'Flow status',
      values = c('#2b8cbe', '#feb24c'),
      labels = c('Flowing', 'Non-flowing')
    ) +
    geom_rect(data = sampling_date_lims, 
              aes(ymin = 0, ymax = 1, 
                  xmin = max_date, xmax = max(sampling_date_lims$max_date)),
              fill = 'white', alpha = 0.8) +
    geom_rect(data = sampling_date_lims, 
              aes(ymin = 0, ymax = 1, 
                  xmin = min(sampling_date_lims$min_date), xmax = min_date),
              fill = 'white', alpha = 0.8) +
    coord_cartesian(expand = FALSE) +
    scale_x_date(name = 'Date') +
    scale_y_continuous(
      name = str_wrap('Mean % of network flowing/non-flowing over previous 90 days', 40),
      labels = scales::label_percent()
    ) +
    facet_wrap(~as.factor(country)) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.justification = c(0, 1),
      legend.box.just = "top",
      legend.box.margin = margin(10, 10, 10, 10),
      plot.margin = margin(10, 50, 10, 10)
    )
  
  # 2. Plot sampling campaigns on top
  plot_relF90past_sampling <- plot_relF90past +
    geom_vline(data = in_allvars_dt, aes(xintercept = date, color = campaign), linewidth = 1) +
    scale_color_gradient(name = 'Campaign', low = 'black', high = 'grey') 
  
  # 3. Plot alpha diversity over time
  in_allvars_dt <- in_allvars_dt %>%
    setorderv(c('organism', 'country', 'site', 'date')) %>%
    .[, paste0(alpha_var, "_tstandard") := get(alpha_var) - .SD[1, get(alpha_var)], 
      by = c('organism', 'country', 'site')] %>%
    .[, mean_campaign_date := mean(date, na.rm = TRUE), by = c('organism', 'country', 'campaign')] %>%
    .[, c(paste0(alpha_var, "_tstandard_scaled"), paste0(alpha_var, "_scaled")) := 
        .(scales::rescale(get(paste0(alpha_var, "_tstandard"))),
          scales::rescale(get(alpha_var))),
      by = 'organism'] %>%
    merge(in_organism_dt, by = 'organism')
  
  alpha_plot_list <- lapply(unique(in_allvars_dt$organism_label), function(in_organism) {
    dt_sub <- in_allvars_dt[organism_label == in_organism, ]
    
    campaign_color_scale <- dt_sub[, uniqueN(campaign)] %>%
      seq(0, 1, length.out = .) %>%
      scales::pal_seq_gradient("black", "black")(.)
    
    alpha_range <- range(dt_sub[[alpha_var]])
    
    plot_alpha <- plot_relF90past + 
      ggnewscale::new_scale('fill') +
      geom_boxplot(data = dt_sub,
                   aes(
                     x = mean_campaign_date,
                     y = get(paste0(alpha_var, "_scaled")),
                     color = factor(campaign),
                     fill = stream_type
                   ),
                   outliers = TRUE,
                   width = 20,
                   position = position_dodge(preserve = "single")) +
      scale_colour_manual(name = 'Campaign', values = campaign_color_scale, guide = 'none') +
      scale_fill_manual(name = 'Stream type',
                        values = c('#2b8cbe', '#feb24c'),
                        labels = c('Perennial', 'Non-perennial')) +
      scale_y_continuous(
        name = "Mean % of network flowing/non-flowing (background)",
        labels = scales::label_percent(),
        sec.axis = sec_axis(
          transform = ~ . * (alpha_range[[2]] - alpha_range[[1]]) + alpha_range[[1]],
          name = paste0(alpha_var, " (boxplots)")
        )
      ) +
      ggtitle(label = in_organism)
    
    return(plot_alpha)
  }) %>% setNames(unique(in_allvars_dt$organism_label))
  
  # Write plots
  if (write_plots) {
    ggsave(file.path(out_dir, 'relF90past_drn.png'), 
           plot = plot_relF90past, 
           width = 10, height = 6, 
           units = 'in', dpi = 600)
    ggsave(file.path(out_dir, 'relF90past_drn_sampling.png'), 
           plot = plot_relF90past_sampling, 
           width = 10, height = 6, 
           units = 'in', dpi = 600)
    
    lapply(names(alpha_plot_list), function(in_organism_label) {
      ggsave(file.path(out_dir, 
                       paste0('relF90past_drn_', alpha_var, '_',
                              in_organism_label, '.png')),
             plot = alpha_plot_list[[in_organism_label]],
             width = 10, height = 6, units = 'in', dpi = 600)
    })
  }
  
  return(list(
    hydro = plot_relF90past,
    hydro_sampling = plot_relF90past_sampling,
    alpha_list = alpha_plot_list
  ))
}
#------ check_diversity_dependence_on_habitat_volume ---------------------------
#n_allvars_merged <- tar_read(allvars_merged)

#' @title Check correlation between alpha diversity and habitat volume
#' @description Creates diagnostic scatter plots to examine the relationship 
#'     between alpha diversity and local habitat volume metrics (e.g., depth, width, velocity),
#'     faceted by country.
#' @param in_allvars_merged A list containing a data.table of all merged variables.
#' @param alpha_var Character string: name of the alpha diversity column to use (e.g., "richness", "shannon").
#' @return A list of ggplot objects showing the relationships for different organism groups.
check_cor_div_habvol <- function(in_allvars_merged, alpha_var = "richness") {
  
  dt <- copy(in_allvars_merged$dt)
  
  # Compute scaled habitat metrics
  dt[, avg_vol_miv := avg_depth_macroinvertebrates * average_wetted_width_m]
  dt[, mean_avg_vol_miv := mean(avg_vol_miv), by = site]
  dt[, mean_velo_miv := mean(avg_velocity_macroinvertebrates), by = site]
  dt[, mean_alpha := mean(get(alpha_var)), by = site]  # scaled alpha diversity
  
  # 1. MIV: volume vs alpha
  p_miv_vol <- ggplot(dt[organism == 'miv', ],
                      aes(x = avg_vol_miv / mean_avg_vol_miv,
                          y = get(alpha_var) / mean_alpha)) +
    geom_point() +
    geom_smooth(aes(color = site), method = 'lm', se = FALSE) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  # 2. MIV: velocity vs alpha
  p_miv_vel <- ggplot(dt[organism == 'miv_nopools' & avg_velocity_macroinvertebrates < 5, ],
                      aes(x = avg_velocity_macroinvertebrates,
                          y = get(alpha_var) / mean_alpha)) +
    geom_point() +
    geom_smooth(aes(color = site), method = 'lm', se = FALSE) +
    geom_smooth(color = 'black', method = 'lm', linewidth = 1.2) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  # 3. Diatoms biofilm: volume vs alpha
  p_dia_biof_vol <- ggplot(dt[organism_class == 'dia_biof', ],
                           aes(x = volume_biofilm,
                               y = get(alpha_var) / mean_alpha)) +
    geom_point() +
    geom_smooth(aes(color = site), method = 'lm', se = FALSE) +
    facet_wrap(organism~country) +
    theme(legend.position = 'none')
  
  # 4. Diatoms biofilm: m2 vs alpha
  p_dia_biof_m2 <- ggplot(dt[organism == 'dia_biof', ],
                          aes(x = m2_biofilm,
                              y = get(alpha_var) / mean_alpha)) +
    geom_point() +
    geom_smooth(aes(color = site), method = 'lm', se = FALSE) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  # 5. Fungi biofilm: m2 vs alpha
  p_fun_biof_m2 <- ggplot(dt[organism == 'fun_biof', ],
                          aes(x = m2_biofilm,
                              y = get(alpha_var) / mean_alpha)) +
    geom_point() +
    geom_smooth(aes(color = site), method = 'lm', se = FALSE) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  # 6. BActeria biofilm: m2 vs alpha
  p_bac_biof_m2 <- ggplot(dt[organism == 'bac_biof', ],
                          aes(x = m2_biofilm,
                              y = get(alpha_var) / mean_alpha)) +
    geom_point() +
    geom_smooth(aes(color = site), method = 'lm', se = FALSE) +
    facet_wrap(~country) +
    theme(legend.position = 'none')
  
  return(list(
    miv_vol = p_miv_vol,
    miv_vel = p_miv_vel,
    dia_biof_vol = p_dia_biof_vol,
    dia_biof_m2 = p_dia_biof_m2,
    fun_biof_m2 = p_fun_biof_m2
  ))
}


#------ plot_areadiv_scatter ---------------------------------------------------
# in_dt=tar_read(allvars_merged)$dt_summarized
# in_organism_dt <- tar_read(organism_dt)
# save_plots=T
# out_dir = figdir

#' @title Plot scatter of alpha diversity vs. basin area and discharge
#' @description Creates scatter plots to visualize the relationship between 
#'      alpha diversity and upstream basin area or simulated discharge. 
#'      The plots are faceted by organism and country.
#' @param in_dt A data.table of summarized variables.
#' @param in_organism_dt A data.table mapping organism codes to their labels.
#' @param alpha_var Character string: name of the alpha diversity column to plot (e.g., "mean_richness", "mean_shannon").
#' @param write_plots Logical: whether to save the plots as PNG files.
#' @param out_dir Output directory for saved plots.
#' @return A list of ggplot objects: one for basin area and one for discharge.
plot_areadiv_scatter <- function(in_dt,
                                 in_organism_dt,
                                 alpha_var = "mean_richness",
                                 write_plots = TRUE,
                                 out_dir = NULL) {
  
  # prepare data and order countries for consistent faceting
  in_dt_labels <- merge(in_dt,
                        in_organism_dt,
                        by = 'organism', all.x = FALSE) %>%
    .[, country := factor(
      country, 
      levels = c("Finland", "France",  "Hungary", "Czech", "Croatia", "Spain"),
      ordered = TRUE)] %>%
    .[!(organism %in% c('miv_nopools_flying', 'miv_nopools_noflying'))]
  
  # 1. Plot alpha diversity vs. basin area
  plot_area <- ggplot(in_dt_labels, 
                      aes(x = basin_area_km2, y = get(alpha_var))) +
    geom_point(aes(color = stream_type)) +
    geom_smooth(aes(color = stream_type), method = 'gam', se = FALSE) +
    geom_smooth(method = 'gam', color = 'black', se = FALSE) +
    scale_x_log10(name = expression('Basin area' ~ km^2)) +
    scale_y_continuous(name = alpha_var) +
    scale_color_manual(name = 'Stream type',
                       values = c('#2b8cbe', '#feb24c'),
                       labels = c('Perennial', 'Non-perennial')) +
    facet_grid(organism_label ~ country, scales = 'free') +
    theme_bw()
  
  # 2. Plot alpha diversity vs. simulated discharge
  plot_qsim <- ggplot(in_dt_labels, 
                      aes(x = qsim_avg_samp, y = get(alpha_var))) +
    geom_point(aes(color = stream_type)) +
    geom_smooth(aes(color = stream_type), method = 'lm', se = FALSE) +
    geom_smooth(method = 'gam', color = 'black', se = FALSE) +
    scale_x_log10(name = expression('Mean simulated discharge during sampling period' ~ m^3 ~ s^-1)) +
    scale_y_continuous(name = alpha_var) +
    scale_color_manual(name = 'Stream type',
                       values = c('#2b8cbe', '#feb24c'),
                       labels = c('Perennial', 'Non-perennial')) +
    facet_grid(organism_label ~ country, scales = 'free') +
    theme_bw()
  
  # 3. Save plots if requested
  if (write_plots) {
    ggsave(filename = file.path(out_dir, paste0('scatter_area_', alpha_var, '.png')),
           plot = plot_area,
           width = 12,
           height = 8,
           units = 'in',
           dpi = 600)
    
    ggsave(filename = file.path(out_dir, paste0('scatter_qsim_', alpha_var, '.png')),
           plot = plot_qsim,
           width = 12,
           height = 8,
           units = 'in',
           dpi = 600)
  }
  
  return(list(
    area = plot_area,
    qsim = plot_qsim
  ))
}

#------ plot_cor_hydrowindow  --------------------------------------------------
#For hydrological variables check by time window

# hydrovar_list <- c('DurD', 'PDurD', 'FreD', 'PFreD',
#                'uQ90', 'oQ10', 'maxPQ', 'PmeanQ',
#                'STcon.*_directed', 'STcon.*_undirected')
# var_substr <- 'STcon.*_directed'
# in_cor_dt <- tar_read(cor_matrices_list)$div_bydrn

#' @title Plot correlation by hydrological time window
#' @description Visualizes the Spearman's correlation between biological variables and 
#' hydrological variables across different time windows.
#' @param in_cor_dt A data table containing correlation coefficients.
#' @param temporal_var_substr A string to filter for specific hydrological variables (e.g., 'DurD').
#' @param response_var_list A character vector of biological response variables to plot.
#' @param colors_list A named character vector of colors for each country.
#' @param save_plot A logical value to save the plot as a PNG file.
#' @param plot_name_suffix A string to append to the output filename.
#' @param out_dir The output directory for the saved plot.
#' @return A ggplot object of the correlation plot.
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
# out_dir = figdir
# in_organism_dt <- tar_read(organism_dt)
# plot_name_suffix = ''
# 
# in_allvars_merged=tar_read(allvars_merged)
# temporal_var_substr='DurD'
# response_var=alpha_var_list
# in_organism_dt = tar_read(organism_dt)
# colors_list=drn_dt$color
# write_plots=T
# plot_name_suffix=""
# out_dir=figdir

#' @title Plot scatter with linear model fit
#' @description Creates scatter plots of a biological response variable versus 
#'      hydrological variables with linear model trend lines, 
#'      faceted by organism and hydrological variable.
#' @param in_allvars_merged A list containing a data table of all merged variables.
#' @param in_organism_dt A data table mapping organism codes to their labels.
#' @param temporal_var_substr A string to filter for specific hydrological variables (e.g., 'DurD').
#' @param response_var A list containing strings specifying the biological response variables (e.g., 'richness').
#' @param colors_list A named character vector of colors for each country.
#' @param write_plots A logical value to save the plot as a PNG file.
#' @param plot_name_suffix A string to append to the output filename.
#' @param out_dir The output directory for the saved plot.
#' @return A ggplot object of the scatter plot with linear model fits.
plot_scatter_lm <-  function(in_allvars_merged, 
                             in_organism_dt,
                             temporal_var_substr, 
                             response_var,
                             colors_list, 
                             write_plots=T,
                             plot_name_suffix="", 
                             out_dir) {
  
  dt <- in_allvars_merged$dt
  x_cols <- grep(paste0('^', temporal_var_substr), names(dt), value=T) %>%
    .[seq(1, length(.), 2)]
  
  dt_melt <- melt(dt, 
                  id.vars=intersect(
                    c(in_allvars_merged$cols$group_cols, response_var), 
                    names(dt)
                  ),
                  measure.vars=x_cols) 
  
  p_list <- lapply(response_var, function(in_var) {
    p_lm <- ggplot(dt_melt[organism %in% in_organism_dt$organism,], 
                   aes_string(x='value', y=in_var, color='country')) +
      geom_point(alpha=0.5) +
      geom_smooth(method='lm', se=F) +
      scale_color_manual(values=colors_list) +
      scale_y_sqrt() +
      facet_grid(organism~variable, scales='free') +
      theme_bw()
    
    if (write_plots) {
      out_path <- file.path(out_dir, paste0('plot_lm_', temporal_var_substr, '_',
                                            in_var, plot_name_suffix, '.png'))
      ggsave(out_path, plot = p_lm, 
             height = 12, width = 12, units='in', dpi = 300)
    }
    return(p_lm)
  }) %>% setNames(response_var)
  
  return(p_list)
  
}


#------ quick_ssn ------
# in_ssn_eu <- tar_read(ssn_eu)
# 
# in_formula = 'richness ~ log10(basin_area_km2) + log10(basin_area_km2):country + DurD60past + DurD60past:country'
# in_ssn <- in_ssn_eu$miv_nopools$ssn
# tar_load(ssn_covtypes)

#' @title Quick SSN model fitting
#' @description A wrapper function to fit linear or generalized linear 
#' spatial stream network (SSN) models with different covariance structures.
#' @param in_ssn A spatial stream network (SSN) object.
#' @param in_formula A formula object or string for the model.
#' @param ssn_covtypes A data table specifying the covariance types 
#'      (functional form of the'taildown', 'tailup', 'euclid' components).
#' @param partition_formula A formula object for a partition factor (e.g., a sampling campaign).
#' @param random_formula A formula object for random effects (e.g., country).
#' @param family A string for the `glm` family, or 'Gaussian' for `lm`.
#' @param estmethod A string for the estimation method (e.g., 'ml' for Maximum Likelihood).
#' @return A named list of fitted SSN model objects, or an object indicating failure.
quick_ssn <- function(in_ssn, in_formula, ssn_covtypes,  
                      partition_formula = as.formula("~ as.factor(campaign)"),
                      random_formula = as.formula("~ country"),
                      family = "Gaussian", # arguments are passed to and evaluated by ssn_lm()
                      estmethod = "ml") {
  # 1. Create distance matrices for the SSN object
  SSN2::ssn_create_distmat(in_ssn)
  
  # 2. Fit models for each combination of covariance structures
  ssn_list <- mapply(function(down_type, up_type, euc_type) {
    label <- paste(down_type, up_type, euc_type, sep = "_")
    message("Fitting model: ", label)
    
    result <- tryCatch({
      if (family %in% c("Gaussian", "gaussian", gaussian())) {
        fit <- ssn_lm(
          formula = as.formula(in_formula),
          ssn.object = in_ssn,
          taildown_type = down_type,
          tailup_type = up_type,
          euclid_type = euc_type,
          additive = "afv_qsqrt",
          partition_factor = partition_formula,
          random = random_formula,
          estmethod = estmethod
        )
      } else {
        fit <- ssn_glm(
          formula = as.formula(in_formula),
          family = eval(family),
          ssn.object = in_ssn,
          taildown_type = down_type,
          tailup_type = up_type,
          euclid_type = euc_type,
          additive = "afv_qsqrt",
          partition_factor = partition_formula,
          random = random_formula,
          estmethod = estmethod
        )
      }
      
      fit$fit_status <- "ok"
      return(fit)
      
    }, error = function(e) {
      warning(paste("Model failed for", label, ":", e$message))
      structure(list(
        fit_status = "failed",
        error_message = e$message,
        label = label,
        taildown = down_type,
        tailup = up_type,
        euclid = euc_type
      ), class = "ssn_glm_failed")
    })
    
    return(result)
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

# in_ssn = tar_read(ssn_eu)
# organism = 'miv_nopools'
# formula_root = ' log10(basin_area_km2) + log10(basin_area_km2):country'
# hydro_var = 'DurD365past'
# response_var = 'expshannon'
# tar_load(ssn_covtypes)
# family = "Gaussian"
# estmethod = "ml"
# partition_formula = as.formula("~ as.factor(campaign)")
# random_formula = as.formula("~ country")

# tar_load(ssn_covtypes)
# in_ssn = tar_read(ssn_eu_summarized)
# organism = c('bac_biof_nopools')
# formula_root = ' sqrt(basin_area_km2)'
# hydro_var = 'DurD3650past'
# response_var = 'mean_expshannon'
# ssn_covtypes = ssn_covtypes[1:10,]
# partition_formula = as.formula('~ country')
# random_formula = NULL
# estmethod='ml'
# family= "Gaussian"

#' @title Model SSN with hydrological variables across time windows
#' @description Automates the process of fitting SSN models to test the effect 
#'      of a specific hydrological variable on a biological response variable, 
#'      collecting model summaries for comparison.
#' @param in_ssn A list of SSN objects, typically one for each organism.
#' @param organism A string specifying the organism to model.
#' @param formula_root A string of fixed effects to be included in all models.
#' @param hydro_var A string specifying the hydrological variable to test.
#' @param response_var A string specifying the biological response variable.
#' @param ssn_covtypes A data table of covariance types for model fitting.
#' @param partition_formula A formula object for a partition factor.
#' @param random_formula A formula object for random effects.
#' @param family A string for the `glm` family, or 'Gaussian' for `lm`.
#' @param estmethod A string for the estimation method.
#' @return A list containing the fitted SSN models and a data table of their summary statistics.
model_ssn_hydrowindow <- function(in_ssn, organism, formula_root, 
                                  hydro_var, response_var, ssn_covtypes,
                                  partition_formula = as.formula("~ as.factor(campaign)"),
                                  random_formula = as.formula("~ country"),
                                  family = "Gaussian",
                                  estmethod = "ml",
                                  include_seasonality=F) {
  
  # hydro_var <- grep(paste0('^', hydro_var_str), 
  #                   names(in_ssn[[organism]]$ssn$obs), 
  #                   value=T)
  
  # 1. Construct the full model formula
  if (!is.null(hydro_var)) {
    full_formula <- paste0(response_var, ' ~ ', hydro_var, ' + ',
                           hydro_var, ':country +', formula_root)
  } else {
    full_formula <- paste0(response_var, ' ~ ', formula_root)
  }
  
  if (include_seasonality) {
    full_formula <- paste(full_formula, '+ doy + I(doy^2)')
  }
  
  # 2. Fit SSN models using quick_ssn
  ssn_list <- quick_ssn(in_ssn = in_ssn[[organism]]$ssn, 
                        in_formula = as.formula(full_formula),
                        family = family,
                        ssn_covtypes = ssn_covtypes, #ssn_covtypes[sample(144, 10)],
                        partition_formula = partition_formula,
                        random_formula = random_formula,
                        estmethod = estmethod)
  
  # 3. Collect model summary statistics (glance)
  ssn_glance <- lapply(names(ssn_list), function(label) {
    model <- ssn_list[[label]]
    
    if (inherits(model, "ssn_glm_failed")) {
      data.table(
        covtypes = label,
        fit_status = "failed",
        error_message = model$error_message,
        AICc=NA
      )
    } else {
      out <- SSN2::glance(model) %>%
        setDT
      out[, covtypes := label]
      out[, fit_status := "ok"]
      out[, error_message := NA_character_]
      return(out)
    }
  }) %>% 
    rbindlist(fill = TRUE) %>%
    setorder(fit_status, AICc) %>%
    .[, `:=`(
      organism = organism,
      response_var = response_var,
      hydro_var = hydro_var,
      family = family
    )]
  
  return(list(
    ssn_list=ssn_list,
    ssn_glance=ssn_glance
  ))
}



#------ plot_hydro_comparison --------------------------------------------------
#------ compare standard hydro metrics with flow-duration curve-based metrics ---
# vars_list <- c('DurD', 'FreD')
# var_substr <- vars_list[[1]]
# in_cor_dt <- tar_read(cor_matrices_list)$div_bydrn
# color_list = drn_dt$color

#' @title Compare standard vs. flow-duration curve-based hydro metrics
#' @description This function generates a scatter plot to compare the correlation 
#'      coefficients of two different hydrological metrics against a biological 
#'      response variable. It helps to visualize how closely two different measures 
#'      of flow, such as `DurD` (duration) and `PDurD` (proportion of duration), 
#'      relate to ecological outcomes.
#' @param in_cor_dt A data.table containing correlation matrices, including `correlation` and `p_value`.
#' @param color_list A named vector of colors for plotting each country.
#' @return A ggplot object of the comparison plot.
plot_hydro_comparison <- function(var_substr, in_cor_dt, color_list) {
  
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

#------ select_ssn_covariance -------------------------------------------------
# in_ssnmodels <- tar_read(ssn_div_hydrowindow_invsimpson)

#' @title Select the best-performing SSN covariance structure
#' @description This function takes a list of fitted SSN models and evaluates 
#'      them to identify the best-performing spatial covariance structure based 
#'      on AICc. It calculates and visualizes the distribution of delta-AICc 
#'      (the difference from the best model's AICc) across all models.
#' @param in_ssnmodels A list of fitted SSN model objects, typically from `model_ssn_hydrowindow`.
#' @return A list containing a data table of covariance statistics, a subset of
#'      the best covariance types, and a boxplot.
select_ssn_covariance <- function(in_ssnmodels) {
  # Combine the 'ssn_glance' data tables from all fitted models into a single data table
  ssn_glance_bind <- lapply(in_ssnmodels, `[[`, "ssn_glance") %>%
    rbindlist(fill=T)
  
  # Calculate delta-AICc for each model, grouped by hydrological variable and organism.
  # delta-AICc is the difference between a model's AICc and the minimum AICc in its group.
  ssn_glance_bind[, delta_AICc := (AICc - min(AICc, na.rm=T)),
                  by=.(hydro_var, organism)] %>%
    # Calculate the mean and median delta-AICc for each covariance type and organism
    .[, `:=`(delta_AICc_covtypemean = mean(delta_AICc, na.rm=T), 
             delta_AICc_covtypemedian = median(delta_AICc, na.rm=T),
             n = .N), by=.(organism, covtypes)]
  
  # Create a unique data table with the aggregated covariance statistics
  covtype_stats <- ssn_glance_bind[, .(delta_AICc_covtypemean,
                                       delta_AICc_covtypemedian,
                                       organism,
                                       covtypes,
                                       n)] %>%
    unique
  
  # Create a boxplot to visualize the distribution of delta-AICc for each covariance type.
  covtype_stats_plot <- ggplot(ssn_glance_bind, 
                               aes(x=tidytext::reorder_within(
                                 as.factor(covtypes), 
                                 delta_AICc_covtypemean,
                                 organism), 
                                 y=delta_AICc)) +
    geom_boxplot() +
    coord_flip() + # Flip the coordinates to make the labels readable
    tidytext::scale_x_reordered() +  
    facet_wrap(~organism, scales='free')
  
  # Select the best covariance type for each organism based on the minimum mean delta-AICc
  covtype_selected <- covtype_stats[, 
                                    .SD[which.min(delta_AICc_covtypemean),], 
                                    by=organism]
  
  return(list(
    dt = covtype_stats,
    dt_sub = covtype_selected,
    plot = covtype_stats_plot
    
  ))
}

#------ format_ssn_hydrowindow -------------------------------------------------
# in_organism <- 'miv_nopools_ept'
# in_response_var = 'richness'
# ssn_div_models_to_run <- tar_read(ssn_div_models_to_run)
# ssn_div_hydrowindow <- tar_read(ssn_div_hydrowindow_richness_miv_nopools_ept)
# in_hydro_vars_dt <- tar_read(hydro_vars_dt)
# in_covtype_selected <- tar_read(ssn_covtype_selected_richness_miv_nopools_ept)
# 
# ssn_model_names <- do.call(rbind, ssn_div_models_to_run)[
#   , c("organism", "hydro_var", "response_var")] %>%
#   as.data.table() %>%
#   .[response_var == in_response_var & organism == in_organism, ]
# 
# in_ssnmodels <- cbind(ssn_model_names, ssn_div_hydrowindow)
# names(in_ssnmodels)[ncol(in_ssnmodels)] <- "ssn_div_models"

#' @title Format and plot SSN model results
#' @description This comprehensive function processes the output of SSN models 
#'      to extract key results, including variance decomposition and model predictions. 
#'      It generates several plots to visualize these results, such as stacked bar 
#'      charts of variance components and scatter plots of observed vs. fitted values.
#' @param in_ssnmodels A data table with fitted SSN models for various hydrological windows.
#' @param in_organism A string specifying the organism to analyze.
#' @param in_covtype_selected A list from `select_ssn_covariance` containing the best covariance types.
#' @return A list containing the formatted data table and several ggplot objects.
format_ssn_hydrowindow <- function(in_ssnmodels, 
                                   in_organism, 
                                   in_covtype_selected,
                                   in_hydro_vars_dt) {
  # Get the selected covariance type for the specified organism
  org_covtype <- in_covtype_selected$dt_sub[organism==in_organism,][['covtypes']]
  
  # Combine model glance tables and attach the full model objects---------------
  ssn_hydrowindow_perf_allvars <- lapply(
    in_ssnmodels[organism==in_organism, ssn_div_models],
    function(x) {
      merge(
        x[['ssn_glance']],
        data.table(covtypes=names(x[['ssn_list']]),
                   mod=x[['ssn_list']]),
        by='covtypes'
      )
    }) %>% 
    rbindlist(fill=T) %>%
    # Filter for the selected covariance type and handle 'null' hydro_var
    .[covtypes==org_covtype,] %>%
    .[is.na(hydro_var), hydro_var := 'null']
  
  
  # Remove the large model object to save memory 
  #(the data table `ssn_hydrowindow_perf_allvars` still contains the `mod` column)
  remove(in_ssnmodels)
  
  #Identify covariance structure with the lowest AICc across all time windows 
  #for each hydrological variable ----------------------------------------------
  get_hydro_var_root(ssn_hydrowindow_perf_allvars)
  
  # Define labels for the hydrological variables for plotting
  ssn_hydrowindow_perf_allvars <- merge(
    ssn_hydrowindow_perf_allvars,
    unique(in_hydro_vars_dt, by='hydro_var_root')[, .(hydro_var_root, hydro_label)],
    by='hydro_var_root')
  
  #Get variance decomposition estimated by model
  varcomp_labels <- c('Remaining variance (nugget)', 
                      'Spatially dependent variance - euclidean',
                      'Spatially dependent variance - taildown',
                      'Spatially dependent variance - tailup',
                      'Random intercept - country',
                      'Fixed effects pseudo-R2'
  )
  
  ssn_hydrowindow_varcomp <- ssn_hydrowindow_perf_allvars[
    ,    {
      if (inherits(mod[[1]], c("ssn_lm", "ssn_glm"))) {
        vc <- SSN2::varcomp(mod[[1]])
        as.data.table(vc)   # multi-row result
      } else {
        data.table()        # return empty DT so it skips
      }
    }
    , by=.(response_var, hydro_var, covtypes, hydro_label, window_d)
  ] %>%
    merge(data.table(
      varcomp=c(
        "nugget",
        "euclid_de", "taildown_de", "tailup_de",         
        "1 | country", "Covariates (PR-sq)"),
      varcomp_label = factor(varcomp_labels, 
                             levels=varcomp_labels, 
                             ordered=T),
      varcomp_color = c('lightgrey', '#b2df8a', '#a6cee3',
                        '#1f78b4', '#f1b6da','#feb24c')
    ),
    by='varcomp')
  
  #Plot variance decomposition as a stacked bar chart --------------------------
  varcomp_subdt <- ssn_hydrowindow_varcomp[proportion > 0 & hydro_label != 'Null',]
  varcomp_null <- ssn_hydrowindow_varcomp[proportion > 0 & hydro_label == 'Null',] %>%
    .[, hydro_label := 'Null: only basin area']
  nrow_pag = 2
  ncol_pag = 3
  page_num <- ceiling(varcomp_subdt[, length(unique(hydro_label))]
                      /(nrow_pag*ncol_pag))
  
  plot_varcomp_list <- lapply(seq_len(page_num), function(in_page) {
    p_main <- ggplot(
      varcomp_subdt,
      aes(x = window_d, y = proportion, fill = varcomp_label)
    ) +
      geom_bar(stat = "identity", position = "stack", alpha=0.75) +
      geom_text(
        aes(label = round(100 * proportion)),
        position = position_stack(vjust = 0.5),
        size=3,
        colour = '#555555') +
      scale_fill_manual(
        name = 'Variance components',
        values=ssn_hydrowindow_varcomp[
          proportion > 0,][order(varcomp_label), unique(varcomp_color)]) +
      scale_x_discrete('Temporal window of aggregation') + 
      scale_y_continuous('Percentage of variance', 
                         breaks=c(0, 0.5, 1),
                         labels = scales::label_percent()) +
      facet_wrap_paginate(~hydro_label, 
                          labeller = label_wrap_gen(width=25),
                          nrow=nrow_pag, ncol=ncol_pag, scales = "free_x",
                          page=in_page) +
      coord_cartesian(expand = FALSE) +
      theme_bw() +
      theme(legend.background=element_blank())
    
    p_null <-  ggplot(
      varcomp_null,
      aes(x = window_d, y = proportion, fill = varcomp_label)
    ) +
      geom_bar(stat = "identity", position = "stack", alpha=0.75) +
      geom_text(
        aes(label = round(100 * proportion)),
        position = position_stack(vjust = 0.5),
        size=3,
        colour = '#555555') +
      scale_fill_manual(
        name = 'Variance components',
        values=ssn_hydrowindow_varcomp[
          proportion > 0,][order(varcomp_label), unique(varcomp_color)]) +
      scale_x_discrete('Temporal window of aggregation (days)') + 
      scale_y_continuous('Percentage of variance', 
                         breaks=c(0, 0.5, 1),
                         labels = scales::label_percent()) +
      coord_cartesian(expand = FALSE) +
      theme_bw() +
      theme(legend.position='none',
            axis.title=element_blank(),
            axis.text=element_blank()) +
      facet_wrap(~hydro_label,
                 labeller = label_wrap_gen(width=15)
                 )
    
    (p_main + guide_area()) + p_null  + 
      plot_layout(
        design = "
    AAAAAAAABBB
    AAAAAAAAC##
    ",
        axes = 'collect',
        guides = 'collect')

  })
  plot_varcomp_list[[3]] 

  
  #Get best model for each variable for single window for each hydrological variable
  ssn_hydrowindow_best <- ssn_hydrowindow_perf_allvars[
    , .SD[which.min(AICc),], 
    by = hydro_var_root]
  

  #Plot marginal means for each "best model" -----------------------------------
  emmeans_dt_all <- lapply(ssn_hydrowindow_best$hydro_var, function(in_pred_var) {
    print(in_pred_var)
    in_mod <- ssn_hydrowindow_best[hydro_var == in_pred_var, mod][[1]]
    if (in_pred_var=='null') {
      in_pred_var <- all.vars(in_mod$formula)[2] #Fragile
    }
    
    in_pred_var_label <-  ifelse((in_pred_var %in% in_hydro_vars_dt$hydro_var) 
                                && (in_pred_var != 'null'),
                                get_full_hydrolabel(in_hydro_vars_dt, in_pred_var),
                                in_pred_var)

    emm_dt <- get_ssn_emmeans(in_mod, in_pred_var, pred_var_label, 
                              interaction_var='country', plot=F)$dt
    
    return(emm_dt)
  }) %>% rbindlist(use.names=T, fill=T)
  
  get_ssn_emmeans(in_emm_dt = emmeans_dt_all, interaction_var='country', plot=T)$plot +
    facet_wrap(~pred_var_name, scales='free')
  

  #Plot estimated coefs and confint for each best model  -----------------------
  emtrends_dt_all <- lapply(ssn_hydrowindow_best$hydro_var, function(in_pred_var) {
    print(in_pred_var)
    in_mod <- ssn_hydrowindow_best[hydro_var == in_pred_var, mod][[1]]
    if (in_pred_var=='null') {
      in_pred_var <- all.vars(in_mod$formula)[2] #Fragile
    }
    
    emm_dt <- get_ssn_emtrends(in_mod, in_pred_var, pred_var_label, 
                               interaction_var='country', plot=F)$dt
    
    return(emm_dt)
  }) %>% rbindlist(use.names=T, fill=T) 
  
  
  emtrends_dt_all[, pred_var_label := get_full_hydrolabel(in_hydro_vars_dt, 
                                                          pred_var_name)]
  
  in_pred_var_label <-  ifelse((in_pred_var %in% in_hydro_vars_dt$hydro_var) 
                               && (in_pred_var != 'null'),
                               ,
                               in_pred_var)
  
  get_ssn_emtrends(in_emtrends_dt = emtrends_dt_all, 
                   interaction_var='country', plot=T)$plot +
    facet_wrap(~pred_var_name, scales='free_x')
  



  plot_ssn_emtrends(in_mod, in_pred_var, 
                    in_pred_var_label, interaction_var, plot=F)
  

  ##############################################################################
  ##############################################################################
  
  #Get model predictions from the best models ----------------------------------
  mod_preds <- lapply(seq(nrow(ssn_hydrowindow_best)), function(i) {
    # Subset the data.table to get the current row
    aug_data <- SSN2::augment(ssn_hydrowindow_best[i, mod[[1]]], 
                              drop=FALSE) %>%
      cbind(ssn_hydrowindow_best[i, .(hydro_var, response_var, covtypes,
                                      hydro_var_root, window_d)])
    return(aug_data)
  }) %>% rbindlist
  
  #Plot observation vs predictions for the best models -------------------------
  resp_var <- ssn_hydrowindow_best$response_var[[1]]
  obs_preds_plot <- ggplot(mod_preds, 
                           aes(x=.fitted, y=get(resp_var), color=country)) +
    geom_abline() +
    geom_point() +
    geom_smooth(method='lm') +
    scale_y_continuous(name=resp_var) +
    facet_grid(country~hydro_var) +
    theme_bw()
  
  # Reshape data for plotting predictor vs. predictions
  mod_preds_melt <- mod_preds[, 
                              c(setdiff(ssn_hydrowindow_best$hydro_var, 'null'),
                                'country',  'basin_area_km2', '.fitted'), 
                              with=F] %>%
    melt(id.vars=c('country', 'basin_area_km2', '.fitted'))
  
  #Plot predictor variable vs predictions
  x_preds_plot <- ggplot(mod_preds_melt, 
                         aes(x = value, y = .fitted, color = country)) +
    geom_point() +
    geom_smooth(method='lm') +
    facet_wrap(country~variable, scales='free') +
    theme_bw()
  
  return(list(
    dt = ssn_hydrowindow_perf_allvars[,.SD, .SDcols = !c('mod')],
    plot_varcomp = plot_ssn_hydrowindow_varcomp,
    plot_obs_preds = obs_preds_plot,
    plot_x_preds = x_preds_plot
    # , plot_marginal = marginal_plot
  ))
}

#------ save_ssn_div_hydrowindow_plots ------------------------------------
# in_ssn_div_hydrowindow_formatted = tar_read(ssn_div_hydrowindow_formatted)
# out_dir = figdir

#' @title Save SSN div hydro-window plots
#' @description This function is a utility to save the plots and data tables generated by `format_ssn_hydrowindow` to a specified output directory. It iterates through the results and saves each plot as a PNG and the performance table as a CSV.
#' @param in_ssn_div_hydrowindow_formatted The formatted list of results from `format_ssn_hydrowindow`.
#' @param out_dir The output directory where files will be saved.
#' @return A list of file paths to the saved plots.
save_ssn_div_hydrowindow_plots <- function(
    in_ssn_div_hydrowindow_formatted,
    in_organism,
    in_response_var,
    out_dir) {
  
  # Save the formatted performance data table as a CSV file. 
  #The 'mod' column is excluded to prevent saving large objects.
  fwrite(in_ssn_div_hydrowindow_formatted$dt, 
         file.path(out_dir, 
                   paste0('ssn_div_hydrowindow_perftable_', 
                          in_organism, '_',
                          in_response_var, '_',
                          format(Sys.Date(), '%Y%m%d'),
                          '.csv'))
  )
  
  # Define a list of plot names to save
  plot_names <- list('plot_varcomp'
                     , 'plot_obs_preds'
                     , 'plot_x_preds'
                     # , 'plot_marginal'
                     )
  
  # Loop through each plot and save it as a PNG file
  lapply(plot_names, function(pname) {
    out_path <- file.path(out_dir, 
                          paste0('ssn_div_hydrowindow_', 
                                 in_organism, '_',
                                 in_response_var, '_',
                                 pname, '_',
                                 format(Sys.Date(), '%Y%m%d'),
                                 '.png'))
    message(paste0('Saving ', out_path))
    
    ggsave(
      filename = out_path,
      plot = in_ssn_div_hydrowindow_formatted[[pname]],
      width = 9,
      height = 9,
      units='in',
      dpi=600
    )
    
    return(out_path)
  })
}

#------ model_miv_t ------------------------------------------------------------
# Regularization: Techniques like LASSO (L1 regularization), Ridge Regression 
# (L2 regularization), and Elastic Net (a combination) can help prevent overfitting 
# by penalizing model complexity. glmnet is a great package for this, and it's easily integrated with caret.

# in_allvars_merged <- tar_read(allvars_merged)
# in_ssn_eu <- tar_read(ssn_eu)
# in_cor_matrices <- tar_read(cor_matrices_list)
# library(glmulti)

#' Macroinvertebrate Model for Individual Sampling Dates
#'
#' This function performs a comprehensive analysis of macroinvertebrate diversity
#' using a combination of linear mixed-effects models (LME) and spatial stream
#' network (SSN) models. It handles data preparation, model training, and
#' diagnostics.
#'
#' @param in_ssn_eu A list containing the SSN object for macroinvertebrates.
#'   Expected to have a structure like `in_ssn_eu$miv_nopools$ssn`.
#' @param in_allvars_merged A list containing merged environmental and hydrological
#'   variables. Expected to have a structure like `in_allvars_merged$cols$hydro_con`.
#' @param in_cor_matrices A list of correlation matrices. Expected to contain
#'   `div` and `div_bydrn` for diversity correlations and `env_cols`.
#' @param ssn_covtypes A data frame or list specifying SSN covariance types
#'   (e.g., `down` and `euc` types for `ssn_lm`).
#'
#' @return A list of various model outputs, including LME models, SSN models,
#'   diagnostic plots, AIC values, and model data tables.
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

#------ model_miv_yr -----------------------------------------------------------
# in_allvars_merged <- tar_read(allvars_merged)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T

#' Macroinvertebrate model for annual data
#'
#' This function models macroinvertebrate alpha diversity averaged across all samplig dates. 
#' It performs data preparation, conducts exploratory data analysis (if `interactive` is TRUE),
#' tests various model structures using a streamlined SSN generalized linear
#' modeling approach, and provides a summary of model performance.
#'
#' The function's primary goal is to identify the best-fitting spatial stream
#' network (SSN) model for predicting macroinvertebrate alpha diversity based on
#' hydrological, environmental, and spatial variables.
#'
#' @param in_ssn_eu_summarized A list containing the summarized SSN object for
#'   macroinvertebrates. Expected to have a structure like `in_ssn_eu_summarized$miv_nopools$ssn`.
#' @param in_allvars_merged A list containing merged environmental and hydrological
#'   variables. This is used to define the list of candidate variables.
#' @param in_cor_matrices A list of correlation matrices, used to check for
#'   relationships between variables and alpha diversity.
#' @param ssn_covtypes A data frame or list specifying SSN covariance types.
#'   This is used to systematically test different spatial covariance structures.
#' @param scale_predictors A logical value. If TRUE, all candidate hydrological 
#'   predictor variables are scaled to have a mean of 0 and SD of 1. Defaults to TRUE.
#' @param interactive A logical value. If TRUE, the function generates a series
#'   of exploratory plots and model diagnostics for interactive inspection. Defaults to FALSE.
#'
#' @return A list containing the following:
#'   \itemize{
#'     \item \strong{model_selection_table}: A data table summarizing the performance
#'       (e.g., AIC, VIF) of all tested models.
#'     \item \strong{ssn_mod_fit}: The final fitted SSN model selected for prediction.
#'     \item \strong{ssn_pred_final}: The augmented data frame with predictions from the
#'       prediction model.
#'     \item \strong{torgegram_mod_pred}: A diagnostic plot (ggplot object) showing the
#'       Torgegram for the prediction model's residuals.
#'     \item \strong{ssn_mod_fit_best}: The "absolute best" model fit, refitted with REML.
#'     \item \strong{ssn_pred_best}: Augmented data for the "best" model.
#'   }
model_miv_yr <- function(in_ssn_eu_summarized,
                         in_allvars_merged,
                         in_cor_matrices, 
                         ssn_covtypes,
                         scale_predictors = T,
                         interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_miv <- in_ssn_eu_summarized$miv_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_miv,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  allvars_dt <- as.data.table(ssn_miv$obs) %>%
    setorderv(c('country', 'site')) 
  
  # Define candidate hydrological variables
  hydro_candidates <- in_allvars_merged$cols$hydro_con_summarized
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  if (scale_predictors) {
    allvars_dt[
      , (hydro_candidates) := lapply(.SD,
                                     function(x) base::scale(x, center=T, scale=T)),
      .SDcols = hydro_candidates]
  }
  
  alpha_var <- 'mean_richness'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_miv$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_miv$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_miv$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_miv$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_miv$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    data.table::melt(ssn_miv$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
  }
  
  #Settle on a basic covariance structure to start  with -----
  #Check correlation
  topcors_overall <- in_cor_matrices$div[
    variable1 == alpha_var & organism == 'miv_nopools' 
    # & variable2 %in% hydro_candidates 
    & !is.na(correlation),] %>%
    .[, abs_cor := abs(correlation)] %>%
    setorder(-abs_cor)
  
  topcors_bydrn <- in_cor_matrices$div_bydrn[
    variable1 == alpha_var & organism == 'miv_nopools' &
      variable2 %in% hydro_candidates & !is.na(correlation),] %>%
    .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
             var_label = paste(variable2, round(correlation, 2))
    ), by=country] %>%
    dcast(cor_order~country, value.var = 'var_label', fill=NA)
  
  topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
    variable1 == alpha_var & organism == 'miv_nopools' &
      variable2 %in% hydro_candidates & !is.na(correlation),] %>%
    .[, mean(correlation), by=variable2] %>%
    .[order(V1),]
  
  #Test all possible covariance structures
  if (interactive) {
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('miv_nopools'),
      formula_root = ' sqrt(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_pois_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",  #Gamma(link = "log"),
      organism = c('miv_nopools'),
      formula_root = ' sqrt(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_nbin_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('miv_nopools'),
      formula_root = ' sqrt(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      ssn_norm_ini$ssn_glance,
      ssn_pois_ini$ssn_glance,
      ssn_nbin_ini$ssn_glance
    ))
    
    #CHoose negative binomial based on the AIC scores
    ssn_nbin_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('miv_nopools'),
      formula_root = ' sqrt(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_nbin_cov$ssn_glance)
    
    loocv(ssn_nbin_cov$ssn_list$none_none_none)
    loocv(ssn_nbin_cov$ssn_list$none_epa_none)
    loocv(ssn_nbin_cov$ssn_list$none_linear_none)
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  
  quick_miv_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'nbinomial',
      ssn.object = ssn_miv,
      taildown_type = 'none',
      tailup_type = 'linear',
      euclid_type = 'none',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = in_random,
      estmethod = estmethod
    )
  }
  
  #Then test multiple models -----
  # This is the main model selection loop, testing various combinations of
  # predictors.
  miv_rich_modlist <- list()
  
  #Null model
  miv_rich_modlist [['null']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  miv_rich_modlist [['mod1']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ sqrt(basin_area_km2)')))
  
  miv_rich_modlist [['mod2']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(basin_area_km2)')))
  
  miv_rich_modlist [['mod3']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ basin_area_km2 + country:sqrt(basin_area_km2)')))
  
  miv_rich_modlist [['mod4']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + sqrt(basin_area_km2) + country:sqrt(basin_area_km2)')))
  
  miv_rich_modlist [['mod5']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod6']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod7']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod8']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod9']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_avg_samp + country:DurD_avg_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod10']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + country:FreD_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod11']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod12']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ Fdist_mean_10past_undirected_avg_samp + country:Fdist_mean_10past_undirected_avg_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod13']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ PFreD365past + country:PFreD365past + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod14']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ STcon_m10_undirected_avg_samp + country:STcon_m10_undirected_avg_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod15']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp) + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod16']] <- quick_miv_ssn(as.formula(paste(alpha_var, '~ DurD3650past + DurD_samp + country:DurD_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod17']] <- quick_miv_ssn(as.formula(paste(alpha_var, '~ DurD3650past + country:DurD3650past + DurD_samp + country:DurD_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod18']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + DurD_samp + country:DurD_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod19']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + FreD_samp + country:FreD_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod20']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + PFreD365past + country:PFreD365past + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod21']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ Fdist_mean_10past_undirected_avg_samp + DurD_samp + country:DurD_samp + sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod22']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ STcon_m10_undirected_avg_samp + DurD_samp + country:DurD_samp + sqrt(basin_area_km2)')))
  
  ###########TRY ADDING ENVIRONMENTAL VARIABLES
  miv_rich_modlist[['mod23']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + DurD_samp + country:DurD_samp + sqrt(basin_area_km2) + env_PC1')))
  
  miv_rich_modlist[['mod24']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + DurD_samp + country:DurD_samp + sqrt(basin_area_km2) + env_PC1 + env_PC2')))
  
  miv_rich_modlist[['mod25']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(basin_area_km2) + env_PC1 + env_PC2')))
  
  miv_rich_modlist[['mod26']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(basin_area_km2) + env_PC2')))
  
  miv_rich_modlist[['mod27']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + DurD_samp + country:DurD_samp + sqrt(qsim_avg_samp) + env_PC1 + env_PC2')))
  
  miv_rich_modlist[['mod28']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + DurD_samp + country:DurD_samp + sqrt(qsim_avg_samp) + env_PC1 + env_PC2')))
  
  miv_rich_modlist[['mod29']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(qsim_avg_samp) + env_PC2')))
  
  miv_rich_modlist[['mod30']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~  DurD3650past + DurD_samp + country:DurD_samp + sqrt(qsim_avg_samp)')))
  
  miv_rich_modlist[['mod31']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~  DurD3650past + DurD_samp + country:DurD_samp + sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp)')))
  
  miv_rich_modlist[['mod32']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(basin_area_km2) + country:sqrt(basin_area_km2) + env_PC1 + env_PC2')))
  
  miv_rich_modlist[['mod33']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(basin_area_km2) + country:sqrt(basin_area_km2) + env_PC2')))
  
  miv_rich_modlist[['mod34']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(basin_area_km2) + country:sqrt(basin_area_km2)')))
  
  miv_rich_modlist[['mod35']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp) + env_PC2')))
  
  miv_rich_modlist[['mod36']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp)')))
  
  miv_rich_modlist[['mod37']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + sqrt(basin_area_km2) + country:sqrt(basin_area_km2) + env_PC2')))
  
  miv_rich_modlist[['mod38']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + sqrt(basin_area_km2) + country:sqrt(basin_area_km2) + env_PC2')))
  
  miv_rich_modlist[['mod39']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp) + env_PC2')))
  
  miv_rich_modlist[['mod40']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp)')))
  
  miv_rich_modlist[['mod41']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + stream_type + country:stream_type + sqrt(basin_area_km2) + country:sqrt(basin_area_km2) + env_PC2'))) 
  
  miv_rich_modlist[['mod42']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + sqrt(basin_area_km2) + country:sqrt(basin_area_km2) + env_PC2'))) 
  
  miv_rich_modlist[['mod43']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + sqrt(basin_area_km2) + country:sqrt(basin_area_km2)'))) 
  
  # miv_rich_modlist[['mod44']] <- quick_miv_ssn(mean_richness ~ DurD_samp*Fdist_mean_10past_undirected_avg_samp + country:DurD_samp + sqrt(basin_area_km2) + country:sqrt(basin_area_km2))
  
  #Summary table
  mod_perf_tab <- lapply(miv_rich_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(miv_rich_modlist)]
  
  # Interactive summaries and diagnostics
  if (interactive) {
    summary( miv_rich_modlist[['mod31']])
    varcomp(miv_rich_modlist[['mod31']])
    summary(miv_rich_modlist[['mod27']])
    varcomp(miv_rich_modlist[['mod27']])
    summary(miv_rich_modlist[['mod28']])
    varcomp(miv_rich_modlist[['mod28']])
    summary(miv_rich_modlist[['mod25']])
    varcomp(miv_rich_modlist[['mod25']])
    summary(miv_rich_modlist[['mod34']])
    varcomp(miv_rich_modlist[['mod34']])
    summary(miv_rich_modlist[['mod36']])
    varcomp(miv_rich_modlist[['mod36']])
    varcomp(miv_rich_modlist[['mod33']])
    summary(miv_rich_modlist[['mod35']])
    summary(miv_rich_modlist[['mod37']])
    summary(miv_rich_modlist[['mod43']])
    varcomp(miv_rich_modlist[['mod43']])
  }
  
  #####CHOOSE MOD 43 for "absolute" best model###############"
  #####CHOOSE MOD 36 OR MOD 34 for continuous predictions###########
  
  #------ For "best" model -------------------------------------------------
  # Refits the selected "best" model using the REML method for better variance estimation.
  ssn_miv_best_final <- quick_miv_ssn(
    as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type + 
      sqrt(basin_area_km2) + country:sqrt(basin_area_km2)')),
    estmethod = 'reml'
  )
  ssn_miv_best_preds <- augment(ssn_miv_best_final,
                                drop=F)
  
  if (interactive) {
    summary(ssn_miv_best_final)
    SSN2::varcomp(ssn_miv_best_final)
    plot(ssn_miv_best_final)
  }
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  
  if (interactive) {
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'nbinomial',
      organism = c('miv_nopools'),
      formula_root = 'sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance_qsim <- ssn_mod_predictions_covtypes$ssn_glance
  }
  
  #Test covariance structures again, this time with REML 
  ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
    in_ssn = in_ssn_eu_summarized,
    family = 'nbinomial',
    organism = c('miv_nopools'),
    formula_root = 'sqrt(basin_area_km2) + country:sqrt(basin_area_km2)',
    hydro_var = 'DurD_samp',
    response_var = alpha_var,
    ssn_covtypes = ssn_covtypes,
    partition_formula = as.formula('~ country'),
    random_formula = NULL,
    estmethod='reml'
  )
  
  selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
  #Turns out area is better than the qsim
  
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = ssn_mod_predictions_covtypes$ssn_list$linear_none_none$formula,
    ssn.object = ssn_miv,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', size=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  summary(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
  varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
  
  summary(ssn_mod_predictions_covtypes$ssn_list$linear_none_none)
  varcomp(ssn_mod_predictions_covtypes$ssn_list$linear_none_none)
  
  #Re-fit final model with REML
  ssn_miv_pred_final <- ssn_glm(
    formula =  ssn_mod_predictions_covtypes$ssn_list$linear_none_none$formula,
    family = 'nbinomial',
    ssn.object = ssn_miv,
    taildown_type = 'linear',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_miv_mod_preds <- augment(ssn_miv_pred_final, drop=F)
  
  if (interactive) {
    tidy(miv_rich_modlist[['mod34']])
    tidy(ssn_miv_pred_final)
    summary(ssn_miv_pred_final)
    plot(ssn_miv_best_final)
    SSN2::varcomp(ssn_miv_pred_final)
    
    #Check predictions
    ggplot(data=ssn_miv_mod_preds, aes(x=.fitted, y=alpha_var)) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_miv_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_miv_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot, #Torgegram for prediction model 
    ssn_mod_fit_best = ssn_miv_best_final, #"Best" model fit with REML
    ssn_pred_best = ssn_miv_best_preds #Augmented data for "best" model
  ))
}

#------ diagnose_ssn_mod -----------------------------------------------------
# in_ssn_mods <- tar_read(ssn_mods_miv_yr)
# write_plots=T
# out_dir = figdir

#' Diagnose SSN Model Performance and Create Diagnostic Plots
#'
#' This function takes a list of fitted SSN models and performs diagnostics. It creates
#' plots of predicted vs. observed values for both the "best" and "final" models
#' and optionally saves them to a specified directory. It also returns variance
#' components for both models.
#'
#' @param in_ssn_mods A list of fitted SSN models, typically generated by another
#'   function (e.g., `model_miv_yr`). It is expected to contain `ssn_pred_best`,
#'   `ssn_pred_final`, `ssn_mod_fit_best`, and `ssn_mod_fit` objects.
#' @param write_plots A logical value. If TRUE, the generated plots are saved to
#'   the directory specified by `out_dir`.
#' @param out_dir A character string specifying the output directory for plots
#'   if `write_plots` is TRUE.
#'
#' @return A list containing the variance components for the "best" and "final" models.
#'   \itemize{
#'     \item \strong{varcomp_best}: Variance components for the best-fit model.
#'     \item \strong{varcomp_final}: Variance components for the final prediction model.
#'   }
diagnose_ssn_mod <- function(in_ssn_mods,
                             write_plots=T,
                             response_var_label,
                             plot_path_prefix='miv_predobs_plot_',
                             out_dir) {
  
  alpha_var <- all.vars(in_ssn_mods$ssn_mod_fit_best$formula)[[1]]
  
  #Check predictions for best model
  predobs_plot_best <- setDT(in_ssn_mods$ssn_pred_best) %>%
    .[, country := factor(
      country, 
      levels = c("Finland", "France",  "Hungary", "Czech", "Croatia", "Spain" ),
      ordered=T)
    ] %>% 
    ggplot(aes(x=.fitted, y=get(alpha_var))) +
    geom_point(aes(color=stream_type)) +
    geom_smooth(method = 'lm',
                color = 'black', se = FALSE) +
    geom_abline(linetype='dashed') +
    scale_x_continuous(name=paste('Predicted', response_var_label)) +
    scale_y_continuous(name=paste('Observed', response_var_label)) +
    scale_color_manual(
      name = 'Stream type',
      labels = c('Perennial', 'Non-perennial'),
      values=c('#2b8cbe', '#feb24c')) +
    coord_fixed() +
    facet_wrap(~country) +
    theme_classic() 
  
  #Check predictions for final model
  predobs_plot_final <- setDT(in_ssn_mods$ssn_pred_final)[, country := factor(
    country, 
    levels = c("Finland", "France",  "Hungary", "Czech", "Croatia", "Spain" ),
    ordered=T)
  ] %>% 
    ggplot(aes(x=.fitted, y=mean_richness)) +
    geom_point(aes(color=stream_type)) +
    geom_smooth(method = 'lm',
                color = 'black', se = FALSE) +
    geom_abline(linetype='dashed') +
    scale_x_continuous(name=paste('Predicted', response_var_label)) +
    scale_y_continuous(name=paste('Observed', response_var_label)) +
    scale_color_manual(
      name = 'Stream type',
      labels = c('Perennial', 'Non-perennial'),
      values=c('#2b8cbe', '#feb24c')) +
    coord_fixed() +
    facet_wrap(~country) +
    theme_classic()
  
  if (write_plots) {
    ggsave(
      filename = file.path(
        out_dir, 
        paste0(plot_path_prefix, '_', alpha_var, '_best.png')),
      plot = predobs_plot_best,
      width = 8,
      height = 6,
      units='in',
      dpi=600
    )
    
    ggsave(
      filename = file.path(
        out_dir, 
        paste0(plot_path_prefix, '_', alpha_var, '_final.png')),
      plot = predobs_plot_final,
      width = 8,
      height = 6,
      units='in',
      dpi=600
    )
  }
  
  return(list(
    varcomp_best = SSN2::varcomp(in_ssn_mods$ssn_mod_fit_best),
    varcomp_final = SSN2::varcomp((in_ssn_mods$ssn_mod_fit))
  )
  )
}

#------ plot_ssn2_marginal_effects ---------------------------------------------
# tar_load(ssn_mods_miv_yr)
# write_plots=T
# out_dir = figdir
# in_mod <- ssn_mods_miv_yr$ssn_mod_fit_best

# #Plot the "marginal" effect of each hydrological variable while
# #keeping surface area at global mean and global intercept
#' Plot Marginal Effects of Variables from a Fitted SSN2 Model
#'
#' This function is designed to visualize the marginal effect of a single
#' hydrological variable on the predicted response, while holding a second
#' variable constant at its mean. It handles both main effects and interaction
#' terms with a grouping factor.
#'
#' @param in_mod A fitted SSN2 model object.
#' @param hydro_var A character string specifying the name of the hydrological
#'   covariate to vary on the x-axis.
#' @param fixed_var A character string specifying the name of the covariate to
#'   hold constant at its mean.
#' @param fixed_var_mean The mean value of the fixed variable.
#' @param group_var A character string specifying the name of the grouping
#'   factor (e.g., "country"). Defaults to "country".
#' @param n_points An integer specifying the number of points to use for the
#'   prediction grid. Defaults to 100.
#'
#' @return A ggplot object visualizing the marginal effects.
plot_ssn_summarized_marginal <- function(
    in_mod,                # fitted SSN2 model
    hydro_var,          # the covariate you want to VARY
    fixed_var,          # the covariate you want to FIX at mean
    fixed_var_mean,
    group_var = "country",
    n_points = 100
) {
  
  in_mod_best <- ssn_mods_miv_yr$ssn_mod_fit_best
  in_ssn <- in_mod_best$ssn.object
  
  emm_best_stream_type <- emmeans(in_mod_best, 
                                  data=in_mod_best$ssn.object$obs, 
                                  specs=~stream_type|country,
                                  type = "response") %>%
    data.frame %>%
    ggplot(aes(x=country, y=prob, color=stream_type)) +
    geom_point(size=2) +
    geom_segment(aes(y=asymp.LCL, yend=asymp.UCL), linewidth=2, alpha=0.5) +
    scale_color_manual(
      name = 'Stream type',
      labels = c('Perennial', 'Non-perennial'),
      values=c('#2b8cbe', '#feb24c')) +
    scale_y_continuous(name='Estimated marginal mean: mean richness') +
    # coord_flip() +
    theme_classic() 
  
  
  ggsave(
    filename = file.path(resdir,
                         'figures',
                         'emm_miv_yr_best_stream_type.png'),
    plot = emm_best_stream_type, 
    width=5, height=5)
  
  
  # Define a grid of basin areas
  newdat_area <- data.frame(
    basin_area_km2 = seq(round(100*min(in_ssn$obs$basin_area_km2, na.rm = TRUE)),
                         round(100*max(in_ssn$obs$basin_area_km2, na.rm = TRUE)),
                         length.out = 50)/100
  )
  
  emm_best_areas <- emmeans(in_mod_best, 
                            data=in_ssn$obs,
                            ~ sqrt(basin_area_km2) | country, 
                            at = list(basin_area_km2 = newdat_area$basin_area_km2),
                            type='response') %>%
    data.frame() %>%
    setDT %>%
    ggplot(aes(x=basin_area_km2, fill=country)) +
    geom_line(aes(y=prob, color=country), linewidth=2) +
    geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=0.1) +
    scale_x_continuous(name=expression('Basin area -'~km^2)) +
    scale_y_continuous(name='Estimated marginal mean: mean richness') +
    theme_bw()
  
  ggsave(
    filename = file.path(resdir,
                         'figures',
                         'emm_miv_yr_best_area.png'),
    plot =emm_best_areas, 
    width=5, height=5)
  
  
  #################################################################
  in_mod_predict <- ssn_mods_miv_yr$ssn_mod_fit
  
  # 
  newdat_durd <- data.frame(
    DurD_samp = seq(0, 0.75, 0.01)
  )
  
  emm_predict_durd <- emmeans(in_mod_predict, 
                              data=in_ssn$obs,
                              ~ DurD_samp | country, 
                              at = list(DurD_samp = newdat_durd$DurD_samp),
                              type='response') %>%
    data.frame() %>%
    setDT %>%
    ggplot(aes(x=DurD_samp, fill=country)) +
    geom_line(aes(y=prob, color=country), size=2) +
    geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=0.1) +
    scale_x_continuous(name='Drying duration during sampling period') +
    scale_y_continuous(name='Estimated marginal mean: mean richness') +
    theme_bw()
  
  ggsave(
    filename = file.path(resdir,
                         'figures',
                         'emm_miv_yr_final_durd_samp.png'),
    plot = emm_predict_durd, 
    width=5, height=5)
  
  # Define a grid of basin areas
  newdat_area <- data.frame(
    basin_area_km2 = seq(round(100*min(in_ssn$obs$basin_area_km2, na.rm = TRUE)),
                         round(100*max(in_ssn$obs$basin_area_km2, na.rm = TRUE)),
                         length.out = 50)/100
  )
  
  emm_predict_areas <- emmeans(in_mod_predict, 
                               data=in_ssn$obs,
                               ~ sqrt(basin_area_km2) | country, 
                               at = list(basin_area_km2 = newdat_area$basin_area_km2),
                               type='response') %>%
    data.frame() %>%
    setDT %>%
    ggplot(aes(x=basin_area_km2, fill=country)) +
    geom_line(aes(y=prob, color=country), size=2) +
    geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=0.1) +
    scale_x_continuous(name=expression('Basin area -'~km^2)) +
    scale_y_continuous(name='Estimated marginal mean: mean richness') +
    theme_bw()
  
  ggsave(
    filename = file.path(resdir,
                         'figures',
                         'emm_miv_yr_best_area.png'),
    plot =emm_best_areas, 
    width=5, height=5)
  
  
  # Get fixed effects
  fixef_dt <- SSN2::tidy(in_mod) %>% as.data.table()
  
  # Identify base terms
  intercept <- fixef_dt[term == "(Intercept)", estimate]
  b_hydro   <- fixef_dt[term == hydro_var, estimate]
  b_fixed   <- fixef_dt[term == fixed_var, estimate]
  
  # Find interaction slopes for hydro
  interaction_pattern <- paste0(group_var, ".*:", hydro_var)
  hydro_interactions <- fixef_dt[grepl(interaction_pattern, term)]
  hydro_interactions[, group := sub(paste0(group_var, "(.*):.*"), "\\1", term)]
  hydro_interactions <- hydro_interactions[, .(group = group, b_hydro_group = estimate)]
  
  # Find interaction slopes for fixed covariate
  fixed_pattern <- paste0(fixed_var, ":", group_var)
  fixed_interactions <- fixef_dt[str_starts(term, fixed(fixed_pattern)),]
  fixed_interactions[, group := sub(paste0(fixed_var, ":.*", group_var, "(.*)"), "\\1", term)]
  fixed_interactions <- fixed_interactions[, .(group = group, b_fixed_group = estimate)]
  
  # Get observed data to build range & mean
  mod_obs <- in_mod$ssn.object$obs %>% as.data.table()
  
  # Build new data grid: hydro varies, fixed stays constant
  grid_hydro <- mod_obs[, .(
    hydro = seq(min(get(hydro_var), na.rm = TRUE),
                max(get(hydro_var), na.rm = TRUE),
                length.out = n_points)
  ), by = group_var]
  
  newdata_dt <- copy(grid_hydro)
  setnames(newdata_dt, old = group_var, new = "group")
  newdata_dt[, fixed := fixed_var_mean]
  
  # Merge interactions
  newdata_dt <- merge(newdata_dt, hydro_interactions, by = "group", all.x = TRUE)
  newdata_dt <- merge(newdata_dt, fixed_interactions, by = "group", all.x = TRUE)
  
  # Fill NA with 0 if no interaction
  newdata_dt[is.na(b_hydro_group), b_hydro_group := 0]
  newdata_dt[is.na(b_fixed_group), b_fixed_group := 0]
  
  # Calculate predicted
  newdata_dt[, predicted :=
               intercept +
               (b_hydro + b_hydro_group) * hydro +
               (b_fixed + b_fixed_group) * fixed
  ]
  
  # Plot
  out_plot <- ggplot(newdata_dt, aes(x = hydro, y = predicted, color = group)) +
    geom_line(linewidth = 1) +
    labs(
      x = hydro_var,
      y = "Predicted",
      color = group_var
    ) +
    theme_bw()
  
  return(out_plot)
}

#------ predict_ssn_mod -------------------------------------------------------
# in_ssn_mods = tar_read(ssn_mods_miv_yr)
# proj_years <- c(2025, 2026, 2030, 2035, 2041)


#' Predict future richness from a fitted SSN model
#'
#' This function takes a fitted SSN model and a list of projection years, then
#' generates predictions for historical and future conditions.
#'
#' The function first generates predictions for historical conditions (year 2021)
#' and then loops through the specified projection years to generate future predictions.
#' It handles the renaming of columns to match the model's formula and subsets
#' the data to include only relevant reaches for efficient computation.
#'
#' @param in_ssn_mods A list of fitted SSN models, typically generated by `model_miv_yr`.
#'   It is expected to contain a `ssn_mod_fit` object.
#' @param proj_years A vector of integers specifying the years for which to generate
#'   future predictions. If NULL, it will use all available years from the projection
#'   data. Defaults to NULL.
#'
#' @return A list containing two data tables:
#'   \itemize{
#'     \item \strong{hist}: A data table with historical predictions.
#'     \item \strong{proj}: A data table with future predictions, including
#'       GCM and scenario information.
#'   }
predict_ssn_mod <- function(in_ssn_mods, proj_years = NULL) {
  #Contemporary predictions ------------------------------------------------------
  # For the year 2021.
  
  preds_hist_dt <- augment(in_ssn_mods$ssn_mod_fit 
                           , newdata = 'preds_hist'
                           , drop = FALSE
                           , type.predict = 'response'
                           , interval = 'prediction'
  ) %>% 
    .[, c('rid', 'country', '.fitted', '.lower', '.upper')] %>%
    st_drop_geometry %>%
    as.data.table %>%
    .[, `:=`(year = 2021,
             colname = 'historical')]
  
  #Future predictions  ---------------------------------------------------------
  preds_proj_colnames <- names(in_ssn_mods$ssn_mod_fit$ssn.object$preds$preds_proj)
  
  # Create a data.table to organize projection column names by GCM, scenario, and year.
  durd_proj_col_dt <- data.table(
    colname = setdiff(
      preds_proj_colnames,
      c(names(in_ssn_mods$ssn_mod_fit$ssn.object$obs),
        names(in_ssn_mods$ssn_mod_fit$ssn.object$edges))
    ) 
  ) %>%
    .[, c('gcm', 'scenario', 'year') := tstrsplit(colname, '_')]
  
  # Determine the projection years to use
  if (is.null(proj_years)) {
    proj_years <- as.integer(unique(durd_proj_col_dt$year))
  }
  
  # Loop through each projection year and generate predictions
  preds_proj_dt <- lapply(
    durd_proj_col_dt[(as.integer(year) %in% proj_years), colname],
    function(in_colname) {
      print(paste0('Predicting for ', in_colname))
      
      ssn_mod_fit_copy <- copy(in_ssn_mods$ssn_mod_fit)
      
      #Subset reaches to only keep those associated with sites for COMPUTATION ONLY
      ssn_mod_fit_copy$ssn.object$preds$preds_proj <- 
        filter(ssn_mod_fit_copy$ssn.object$preds$preds_proj,
               rid %in%  unique(ssn_mod_fit_copy$ssn.object$obs$rid))
      
      # Rename the current projection column to match the model's formula
      names(ssn_mod_fit_copy$ssn.object$preds$preds_proj)[
        preds_proj_colnames == in_colname] <- 'DurD_samp'
      
      # Generate predictions using the augmented data.
      preds_proj_pts <- SSN2::augment(ssn_mod_fit_copy  
                                      , newdata = 'preds_proj'
                                      , drop = FALSE
                                      , type.predict = 'response',
                                      , interval = 'prediction'
      ) %>% 
        .[, c('rid', 'country', '.fitted', '.lower', '.upper')] %>%
        st_drop_geometry %>%
        as.data.table %>%
        .[, colname := in_colname]
      
      return(preds_proj_pts)
    }
  ) %>%
    rbindlist %>%
    merge(durd_proj_col_dt,
          by='colname')
  
  #Write out results  ---------------------------------------------------------
  return(list(
    hist = preds_hist_dt,
    proj = preds_proj_dt
  )
  )
}

#------ map_ssn_mod --------------------------------------------------------
# in_ssn_mods <- tar_read(ssn_mods_miv_yr)
# in_ssn <- tar_read(ssn_eu_summarized)
# in_ssn_preds <- tar_read(ssn_preds)
# out_dir = figdir

#' Map SSN model predictions and diagnostic metrics
#'
#' This function takes a fitted SSN model, a spatial stream network (SSN) object,
#' and model predictions to generate and save a variety of diagnostic plots and maps.
#' It computes future statistics across different GCMs and time periods, and
#' creates plots and maps to visualize these results.
#'
#' @param in_ssn The input SSN object.
#' @param in_ssn_mods A list of fitted SSN models, as returned by `model_miv_yr`.
#' @param in_ssn_preds A list of historical and projected predictions, as returned by
#'   `predict_ssn_mod`.
#' @param out_dir A character string specifying the output directory where the
#'   generated maps and plots will be saved.
#'
#' @return The function is primarily used for its side effects (creating plots and
#'   maps) and does not return a value.
map_ssn_mod <- function(in_ssn,
                        in_ssn_mods,
                        in_ssn_preds,
                        out_dir) {
  
  #Compute future statistics by decade
  preds_period_mean <- in_ssn_preds$proj %>%
    .[, decade := floor(as.integer(year)/10)*10] %>%
    .[, list(fitted_mean_decade = mean(.fitted, na.rm=T),
             lower_mean_decade = mean(.lower, na.rm=T),
             upper_mean_decade = mean(.upper, na.rm=T)),
      by=.(rid, country, gcm, scenario, decade, colname)] %>%
    merge(y=in_ssn_preds$hist[, .(rid, country, .fitted)],
          by=c('rid', 'country')) %>%
    merge(
      data.table(
        decade = c(2010,2020, seq(2040, 2090, 10)),
        period = c(rep('2015-2020', 2), rep('2040-2069', 3), rep('2070-2099', 3))
      ),
      by='decade') %>%
    .[, nyears := .N,
      by=c('country', 'rid', 'period',
           'gcm','scenario', 'decade')]
  
  #Mean and SNR across gcms
  proj_stats_gcm <- preds_period_mean[, list(
    fitted_mean_period=weighted.mean(fitted_mean_decade, nyears, na.rm=T)
  ), by=c('country', 'rid', 'period', 'gcm', 'scenario')] %>%
    dcast(rid+country+gcm+scenario~period, value.var = 'fitted_mean_period') %>%
    .[, preds_diff_2040_2069 := `2040-2069` - `2015-2020`] %>%
    .[, preds_diff_2070_2099 := `2070-2099` - `2015-2020`] %>%
    .[, preds_per_2040_2069 := (`2040-2069` - `2015-2020`)/`2015-2020`] %>%
    .[, preds_per_2070_2099 := (`2070-2099` - `2015-2020`)/`2015-2020`] %>%
    melt(id.vars=c('rid','country','gcm','scenario'))
  
  proj_stats_multigcm <- proj_stats_gcm[, list(
    multigcm_avg = mean(value, na.rm=T),
    snr = mean(value, na.rm=T)/diff(range(value, na.rm=T))
  ), by = c('country', 'rid', 'variable', 'scenario')] %>%
    merge(st_drop_geometry(in_ssn_mods$ssn_mod_fit$ssn.object$obs)[
      , c('rid', 'stream_type')]
      , by='rid')
  
  #------  Make a plot of projections across time periods and countries -----
  ggplot(proj_stats_multigcm[variable %in% c('preds_diff_2040_2099',
                                             'preds_diff_2070_2099'),],
         aes(x=variable, y=multigcm_avg, color=scenario)) +
    geom_boxplot() +
    facet_wrap(~country)
  
  # Another plot to show average percent change in richness
  ggplot(proj_stats_multigcm[variable %in% c('preds_per_2040_2069',
                                             'preds_per_2070_2099'),],
         aes(x=factor(variable), y=multigcm_avg, fill=scenario)) +
    geom_hline(yintercept=0) +
    scale_y_continuous(name='Average % change in richness across GCMs compared to 2015-2020',
                       labels=scales::label_percent()) +
    scale_x_discrete(name='Time Period',
                     labels=c('2040-2060', '2070-2099')) +
    geom_boxplot() +
    scale_fill_brewer(palette='YlOrRd') +
    facet_grid(country~stream_type) +
    theme_classic()
  
  #------------- Makes maps ----------------------------------------------------
  # Make map of 2021 (historical) richness -------------
  ssn_preds_hist <- in_ssn_mods$ssn_mod_fit$ssn.object
  response_var <- all.vars(in_ssn_mods$ssn_mod_fit$formula)[[1]]
  ssn_preds_hist$edges <- merge(
    in_ssn_mods$ssn_mod_fit$ssn.object$edges,
    in_ssn_preds$hist,
    by = c('country', 'rid'), all.x=T)
  
  map_hist <- map_ssn_facets(in_ssn=ssn_preds_hist, 
                             facet_col='country',
                             linewidth_col='qsim_avg',                                
                             linecolor_col='.fitted', 
                             page_title='Macroinvertebrate richness - Historical'
  )
  
  ggsave(
    filename = file.path(figdir, 'ssn_pred_hist_map_miv.png'),
    plot =  map_hist,
    width = 9,
    height = 9,
    units='in',
    dpi=600
  )
  
  #SHOULD REPLACE THIS WITH LOOCV - BUT DOESN't EXIST YET
  ssn_preds_hist$obs <- ssn_preds_hist$obs %>%
    merge(in_ssn_preds$hist[, c('rid', '.fitted', 'country'), with=F], 
          by=c('rid', 'country')) %>%
    mutate(pred_error = `.fitted`-get(response_var))
  
  map_error <- map_ssn_facets(in_ssn=ssn_preds_hist, 
                              in_pts='obs',
                              facet_col='country',
                              ptcolor_col = 'pred_error',
                              shape_col = 'stream_type',
                              linewidth_col='qsim_avg',                                
                              page_title='Macroinvertebrate richness - Error (predicted - observed)'
  )
  
  # Makes maps of projections  ----------------------------
  ssn_preds_proj <- in_ssn_mods$ssn_mod_fit$ssn.object
  
  # ggplot(proj_stats_gcm, aes(x=period, y=preds_diff_multigcm_avg, fill=scenario)) +
  #   geom_boxplot() + 
  #   facet_wrap(~country)
  # 
  
  # Creates a map of a specific projection scenario (2070-2099, ssp585)
  ssn_preds_proj$obs <- ssn_preds_proj$obs %>%
    merge(proj_stats_multigcm[variable=='preds_per_2070_2099' & scenario=='ssp585',], 
          by=c('rid', 'country'))
  map_ssn_facets(in_ssn=ssn_preds_proj, 
                 in_pts='obs',
                 facet_col='country',
                 ptcolor_col = 'multigcm_avg',
                 shape_col = 'stream_type.x',
                 linewidth_col='qsim_avg',                                
                 page_title='Macroinvertebrate mean richness - 2070-2099 - SSP5-8.5'
  )
  
  # Creates a map of the average multimodel richness
  #Plot each period and scenario
  ssn_preds_proj$edges <- merge(
    ssn_preds_proj$edges,
    proj_stats_multigcm[variable=='2070-2099' & scenario=='ssp585',],
    by = c('country', 'rid'), all.x=T)
  
  map_diff_avg <- map_ssn_facets(in_ssn=ssn_preds_proj, 
                                 facet_col='country',
                                 linewidth_col='qsim_avg',                                
                                 linecolor_col='multigcm_avg', 
                                 page_title=NULL) 
  
  # Makes maps of specific projection scenario/year ----------------------------
  # ssn_preds_proj <- in_ssn_mods$ssn_mod_fit$ssn.object
  # in_colname <- 'gfdl-esm4_ssp585_2026'
  # ssn_preds_proj$edges <- merge(
  #   in_ssn_mods$ssn_mod_fit$ssn.object$edges,
  #   in_ssn_preds$proj[colname==in_colname,],
  #   by = c('country', 'rid'), all.x=T)
  # 
  # map_proj <- map_ssn_facets(in_ssn=ssn_preds_proj, 
  #                            facet_col='country',
  #                            linewidth_col='qsim_avg',                                
  #                            linecolor_col='.fitted', 
  #                            page_title=in_colname
  # )
  
  # Map the difference for a specific GCM/scenario/year
  # ssn_preds_proj$edges <- merge(
  #   in_ssn_mods$ssn_mod_fit$ssn.object$edges,
  #   in_ssn_preds$proj_stats[decade==2040 & colname=='gfdl-esm4_ssp585_2040',],
  #   by = c('country', 'rid'), all.x=T)
  # 
  # map_diff_indiv <- map_ssn_facets(in_ssn=ssn_preds_proj,
  #                                  facet_col='country',
  #                                  linewidth_col='qsim_avg',
  #                                  linecolor_col='preds_diff',
  #                                  page_title=unique( ssn_preds_proj$edges$colname)
  # )
  
}

#------ plot_ssn_proj ----------------------------------------------------------
# ggplot(preds_decade_mean, aes(x=scenario, y=preds_diff, fill=as.factor(decade))) +
#   geom_boxplot() +
#   facet_wrap(~country)


################################################################################
################################################################################

#------ tabulate_cor_matrix ----------------------------------------------------
#------ plot_spdiv_drn ---------------------------------------------------------
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

