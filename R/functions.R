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

#------ preformat_intermittence_stats ------------------------------------------
preformat_intermittence_stats <- function(dt) {
  dt[!is.na(isflowing),
     `:=`(noflow_period = rleid(isflowing==0), #Compute no-flow periods
          last_noflow = zero_lomf(isflowing==0, first=FALSE) #Compute row index of last no flow day
     ), by=.(reach_id)] %>%
    .[isflowing == 1, noflow_period := NA] %>% 
    .[last_noflow ==0, last_noflow := NA] %>%
    .[!is.na(noflow_period), noflow_period_dur := .N,  #Compute duration of each no-flow period
      by = .(noflow_period, reach_id)] %>%
    .[, last_noflowdate := .SD[last_noflow, date], by = .(reach_id)] %>% #Convert to data
    .[, PrdD := difftime(date, last_noflowdate, units='days'), #Compute time to last no flow date
      by = .(reach_id)] %>%
    .[, last_noflow := NULL]
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

#------ compute_sd6 ------------------------------------------------------------
compute_rolling_sd6 <- function(dt, time_window_yr) {
  #Identify contiguous six months with the most zero-flow days for computing Sd6
  daily_classified <- dt[, identify_drywet6mo(.SD), by = reach_id]
  daily_classified[, ym := format(date, "%Y%m")] # Add month
  
  # Collapse by year-month based on majority vote (>50% of days in dry period)
  monthly_classified <- daily_classified[
    , .(dry_6mo = mean(dry_6mo, na.rm = TRUE) > 0.5),  
    by = .(reach_id, ym)
  ] %>%
    .[, year := as.integer(substr(ym, 1, 4))]
  
  # For each month, check if no-flow occurred
  noflow_by_month <- dt[
    , .(has_noflow = any(!is.na(noflow_period))),
    by = .(reach_id, ym = format(date, "%Y%m"))
  ]
  
  # Join no-flow with dry/wet classification
  monthly_dt <- merge(noflow_by_month, monthly_classified,
                      by = c('reach_id', 'ym'), all.x=T)
  
  # Count number of months with no-flow per year & season
  base_dt <- monthly_dt[
    , .(n_months = sum(has_noflow, na.rm = TRUE)), 
    by = .(year, dry_6mo, reach_id)
  ] %>%
    dcast(year + reach_id ~ dry_6mo, value.var = "n_months", fill = 0)
  
  #Compute SD6 in a rolling window
  rolling_sd6_inner <- function(dt, k) {
    colname <- paste0("sd6_", k, "yrpast")
    
    dt[, (colname) :=
         frollapply(year, n = k, align = "right",
                    FUN = function(yrs) {
                      subdt <- .SD[year %in% yrs]   # only this reach_id’s rows
                      val <- 1 - subdt[, (mean(`FALSE`, na.rm=TRUE) / mean(`TRUE`, na.rm=TRUE))]
                      fifelse(is.na(val) | mean(`TRUE`, na.rm=TRUE) == 0, 1, val)
                    }),
       by = reach_id] 
  }
  
  out_sd6_dt <- rolling_sd6_inner(base_dt, k = time_window_yr)
  out_sd6_dt[, `:=`(ym=NULL, dry_6mo=NULL, `FALSE`=NULL, `TRUE`=NULL)]
  
  return(out_sd6_dt)
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
  
  # -- Compute statistics for specific reaches -------------------------------
  # DurD: DryDuration
  # PDurD: DryDuration_relative_to_longterm
  # FreD: Drying frequency - absolute or relative number of drying events per time interval
  # PFreD: FreD_relative_to_longterm
  if (scale %in% c('all', 'site')) {
    hydromod_dt_sites <- in_hydromod_dt[reach_id %in% unique(in_sites_dt$reach_id),] %>%
      setorderv(c('reach_id', 'date'))
    
    #Compute duration of no-flow periods and time since last no-flow period #PrdD: prior days to last dry/pool/flowing event
    preformat_intermittence_stats(hydromod_dt_sites)
    
    
    # ggplot(hydromod_dt_sites, aes(x=date, y=time_to_lastzero)) +
    #   geom_line() +
    #   facet_wrap(~id)
    
    #Compute moving-window DurD and FreD statistics --------------------------------------------
    rollingstep_short <- c(10, 30, 60, 90, 120, 180)
    rollingstep_long <- c(365, 365*5, 365*10)
    rollingstep <- c(rollingstep_short, rollingstep_long)
    hydromod_dt_sites[, paste0("DurD", rollingstep, "past") :=  
                        frollapply(isflowing, n=rollingstep, 
                                   FUN=function(x) sum(x==0)/length(x[!is.na(x)]), 
                                   align='right'), by = .(reach_id)
    ]

    hydromod_dt_sites[
      , paste0("FreD", rollingstep, "past") :=
        frollapply(
          seq_len(.N), n = rollingstep,
          FUN = function(idx) {
            # get the slice of the rolling window
            i <- idx
            x_noflow <- noflow_period[i]
            x_valid  <- isflowing[i]  # use isflowing to know which days are valid
            
            valid_idx <- !is.na(x_valid)
            if (!any(valid_idx)) return(NA_real_)
            
            # number of unique no-flow events in this window
            length(unique(x_noflow[!is.na(x_noflow)]))/length(valid_idx)
          },
          align = "right"
        ),
      by = reach_id
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
    
    #Compute annual timing and interannual variability--------------------------
    hydromod_dt_yr <- hydromod_dt_sites[, list(
      DurD_yr = sum(isflowing==0)/.SD[!is.na(isflowing), .N],
      FreD_yr = length(unique(noflow_period[!is.na(noflow_period)]))/.SD[!is.na(isflowing), .N],
      FstDrE = .SD[isflowing==0, min(doy, na.rm=T)],
      meanConD_yr = .SD[!duplicated(noflow_period), mean(noflow_period_dur, na.rm=T)]             
    ), by=.(reach_id, year = as.integer(format(date, '%Y')))] %>%
      .[is.infinite(FstDrE), FstDrE := 365] %>% #For sites that don't dry, set it to 365
      .[is.na(meanConD_yr), meanConD_yr := 0]
    
    
    #CV of annual number of no-flow days----------------------------------------
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
    
    #CV of annual number of no-flow events--------------------------------------
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
    
    #CV of average annual event duration----------------------------------------
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
    
    #SD No-flow start date - Julian day of the first drying event---------------
    hydromod_dt_yr[
      , paste0("FstDrE_SD", rollingstep_yr, "yrpast") :=  
        frollapply(FstDrE, n=rollingstep_yr, 
                   FUN=function(x) sd(x, na.rm=T)
                   , 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #SD6: seasonal predictability of no-flow events (Gallart et al. 2012)-------
    sd6_10yr_dt <- compute_rolling_sd6(dt=hydromod_dt_sites, 10)
    sd6_30yr_dt <- compute_rolling_sd6(dt=hydromod_dt_sites, 30)
    
    #Merge everything
    hydromod_dt_sites <- merge(in_sites_dt[, .(site, reach_id)], 
                               hydromod_dt_sites, 
                               by = 'reach_id', all.x = T, all.y = F,
                               allow.cartesian = T) %>%
      .[, year := as.integer(format(date, '%Y'))] %>%
      merge(hydromod_dt_yr, by=c('reach_id', 'year')) %>%
      merge(sd6_10yr_dt, by=c('reach_id', 'year')) %>%
      merge(sd6_30yr_dt, by=c('reach_id', 'year')) %>%
      setorderv(c('reach_id', 'date'))
    
  }
  
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

#------ import_hydromod_gcm ----------------------------------------------------
# hydroproj_path <- file.path(
#   wp1_data_gouv_dir,
#   "projections",
#   "Albarine_flowstate_projection_gfdl-esm4_ssp585_2015-2100_spatially-distributed.nc")
# 
# hydroref_path <- file.path(
#   wp1_data_gouv_dir,
#   "projections",
#   "Albarine_flowstate_projection_gfdl-esm4_ssp370_1985-2014_spatially-distributed.nc")

# in_drn_dt <- drn_dt

import_hydromod_gcm <- function(hydroproj_path, hydroref_path, in_drn_dt) {
  #Decompose name
  metadata_dt <- str_split(basename(hydroproj_path), '_')[[1]] %>%
    setNames(c('catchment', 'varname', 'time_period', 'gcm', 
               'scenario', 'date_range', 'file_end')) %>%
    as.list %>%
    data.frame %>%
    setDT %>%
    .[catchment == "Lepsamanjoki", catchment := "Lepsamaanjoki"] %>%
    merge(in_drn_dt[, .(catchment, country)], by='catchment') %>%
    .[, path := hydroproj_path]
  
  #Import and format reference hydrological data -------------------------------
  nc_ref <- nc_open(hydroref_path) # open netcdf file
  reachID_ref <- ncvar_get(nc_ref, "reachID") # get list of reaches IDs
  dates_ref <- ncvar_get(nc_ref, "date") # get dates of simulation period
  dates_ref <- as.Date(dates_ref, origin="1950-01-01") # convert dates into R date format
  
  #Get simulated variable from the NetCDF file
  hydro_ref_dt <- get_nc_var_present(nc = nc_ref, varname = metadata_dt$varname, # 0=dry, 1=flowing
                                     reachID = reachID_ref, dates = dates_ref,
                                     selected_sims = NULL)$data_all
  
  
  #Import and format projected hydrological data -------------------------------
  nc <- nc_open(hydroproj_path) # open netcdf file
  reachID <- ncvar_get(nc, "reachID") # get list of reaches IDs
  dates <- ncvar_get(nc, "date") # get dates of simulation period
  dates <- as.Date(dates, origin="1950-01-01") # convert dates into R date format
  
  #Get simulated variable from the NetCDF file
  hydro_proj_dt <- get_nc_var_present(nc = nc, varname = metadata_dt$varname, # 0=dry, 1=flowing
                                      reachID = reachID, dates = dates,
                                      selected_sims = NULL)$data_all 
  
  # Bind reference and projected data ------------------------------------------
  hydro_dt <- rbind(hydro_ref_dt[date<hydro_proj_dt[,min(date)],],
                    hydro_proj_dt[date>hydro_ref_dt[,max(date)],])
  remove(list=c('nc_ref', 'nc', 'hydro_ref_dt', 'hydro_proj_dt'))
  
  setnames(hydro_dt, 
           c('flowstate', 'discharge'),
           c('isflowing', 'qsim'), 
           skip_absent = TRUE)
  
  return(list(
    metadata_dt = metadata_dt,
    hydro_dt = hydro_dt
  ))
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
      cor_result <- Hmisc::rcorr(as.matrix(.SD[, x_cols, with=FALSE]), 
                                 type = correlation_type)
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
#------ scale_ssn_predictors ---------------------------------------------------
# Scale predictor data (mean of 0 and SD of 1)
scale_ssn_predictors <- function(in_ssn, in_vars, scale_ssn_preds) {
  #Scale obs and preds
  vars_in_obs <- intersect(in_vars, names(in_ssn$obs))
  
  in_ssn$obs %<>% mutate(across(all_of(vars_in_obs), ~ scale(.x), 
                                .names = "{.col}_z"))
  
  # Compute scaling parameters from obs
  scaling_means <- in_ssn$obs[paste0(vars_in_obs, '_z')] %>%
    st_drop_geometry() %>%
    sapply(function(x) attr(x, "scaled:center"))
  names(scaling_means) <- gsub('_z$', '', names(scaling_means))
  scaling_sds <- in_ssn$obs[paste0(vars_in_obs, '_z')] %>%
    st_drop_geometry() %>%
    sapply(function(x) attr(x, "scaled:scale"))
  names(scaling_sds) <- gsub('_z$', '', names(scaling_sds))
  
  in_ssn$obs %<>% mutate_at(paste0(vars_in_obs, '_z'), as.numeric)
  
  # Apply same scaling to preds
  if (scale_ssn_preds) {
    vars_in_preds <- intersect(vars_in_obs, names(in_ssn$preds$preds_proj))
    # in_ssn$preds$preds_hist <- in_ssn$preds$preds_hist %>%
    #   mutate(across(all_of(vars_in_preds),
    #                 ~ (.x - scaling_means[cur_column()]) / scaling_sds[cur_column()],
    #                 .names = "{.col}_z"))
    in_ssn$preds$preds_proj %<>%
      mutate(across(all_of(vars_in_preds),
                    ~ (.x - scaling_means[cur_column()]) / scaling_sds[cur_column()],
                    .names = "{.col}_z"))
  }
  
  return(in_ssn)
}
#------ get_hydro_var_root -----------------------------------------------------
get_hydro_var_root <- function(dt, in_place=TRUE) {
  if (!in_place) {
    dt <- copy(dt)
  }
  dt[, hydro_var_root := gsub(
    "(_*[0-9]+(yr)*past)|(_*m[0-9]+)|(_scaled)", "", hydro_var)] %>%
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
  # print(in_hydro_var)
  h_unit <- ifelse(grepl(pattern='.*yrpast.*', in_hydro_var), 'y', 'd')
  
  if (grepl('_scaled', in_hydro_var)) {
    # scaled_ext <- 'scaled - '
    in_hydro_var <- gsub('_scaled', '', in_hydro_var)
  } else {
    # scaled_ext <- ''
  }
  
  
  if (in_hydro_var %in% in_hydro_vars_dt$hydro_var) {
    out_label <- in_hydro_vars_dt[hydro_var == in_hydro_var,] %>%
      .[1, paste(hydro_label, ': past', window_d, h_unit)] 
  } else {
    out_label <- in_hydro_var
    
    if (in_hydro_var == 'basin_area_km2') {
      out_label <- 'Basin area - square kilometers'
    }
  }
  
  return(out_label)
}

#------ get_ssn_emmeans -------------------------------------------------------
get_ssn_emmeans <- function(in_mod, 
                            in_pred_var, interaction_var, 
                            in_drn_dt,
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
      merge(in_drn_dt, by='country')  %>%
      .[, country := factor(
        country,
        levels = c("Finland", "France",  "Hungary", "Czechia", "Croatia", "Spain" ),
        ordered=T)
      ] %>%
      setnames(in_pred_var, 'pred_var')
    
    if (!is.null(in_pred_var_label)) {
      emm_dt[, pred_var_label := in_pred_var_label]
    }
    
  } else if (is.data.table(in_emm_dt) && 
             all(c('pred_var', interaction_var, 'emmean')  %in% names(in_emm_dt))
  ){
    emm_dt <- in_emm_dt
    response_var <- emm_dt[1, response_var]
  }
  
  setnames(emm_dt, c('prob', 'rate'), rep('emmean', 2), skip_absent=TRUE)
  
  color_vec <- emm_dt[!duplicated(get(interaction_var)),
                      setNames(color, get(interaction_var))]
  
  if (plot) {
    emm_plot <- emm_dt %>%
      ggplot(aes(x=pred_var, fill=get(interaction_var))) +
      geom_line(aes(y=emmean, color=get(interaction_var)), 
                linewidth=2, alpha=0.7) +
      geom_ribbon(aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=0.1) +
      scale_y_continuous(name=paste('Estimated marginal mean:', 
                                    Hmisc::capitalize(response_var))) +
      scale_color_manual(name=Hmisc::capitalize(interaction_var),
                         values=color_vec ) +
      scale_fill_manual(name=Hmisc::capitalize(interaction_var),
                        values=color_vec) + 
      theme_minimal()
    
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
                             in_pred_var_label,
                             interaction_var,
                             in_drn_dt,
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
               pred_var_name = in_pred_var)]  %>%
      merge(in_drn_dt, by='country') %>%
      .[, country := factor(
        country,
        levels = c("Finland", "France",  "Hungary", "Czechia", "Croatia", "Spain" ),
        ordered=T)
      ] %>%
      setnames(paste0(in_pred_var, '.trend') , 'trend')
    
  }  else if (is.data.table(in_emtrends_dt)) {
    emtrends_dt <- in_emtrends_dt
    response_var <- emtrends_dt[1, response_var]
  }
  
  if (all(c('lcl_capped', 'ucl_capped') %in% names(emtrends_dt))) {
    p_rangemin <- 'lcl_capped'
    p_rangemax <- 'ucl_capped'
  } else {
    p_rangemin <- 'asymp.LCL'
    p_rangemax <- 'asymp.UCL'
  }
  
  color_vec <- emtrends_dt[!duplicated(get(interaction_var)),
                           setNames(color, get(interaction_var))]
  
  if (plot) {
    emtrends_plot <- ggplot(emtrends_dt, aes(x = get(interaction_var), 
                                             y = trend)) +
      geom_pointrange(aes(ymin = get(p_rangemin), ymax = get(p_rangemax),
                          color=get(interaction_var)),
                      alpha=0.5, fatten=6) +
      geom_hline(yintercept=0, linetype=2) +
      scale_color_manual(name=Hmisc::capitalize(interaction_var),
                         values=color_vec ) +
      labs(y = paste("Estimated slope of", response_var),
           x = Hmisc::capitalize(interaction_var)) +
      coord_flip(clip='off') +
      theme_minimal() 
    
  } else {
    emtrends_plot <- NULL
  }
  
  return(list(
    dt = emtrends_dt,
    plot = emtrends_plot
  ))
}

#------ check_resid_corr -------------------------------------------------------
#Check correlation of residuals with other variables
check_resid_corr <- function(in_ssn_mod, in_idcol='site', 
                             in_response_var, in_candidates) {
  augment_dt <- augment(in_ssn_mod, 
                        drop=FALSE, type.residuals="response") %>%
    st_drop_geometry %>%
    as.data.table
  
  augment_melt <- melt(augment_dt, 
                       id.vars=c(in_idcol, in_response_var, 'country', '.resid'),
                       measure.vars = in_candidates) 
  
  resid_corr_plot <- ggplot(augment_melt, 
                            aes(x=value, y=`.resid`, color=country)) +
    geom_point() +
    geom_smooth(method='lm',se=FALSE) +
    facet_wrap(~variable, scales='free_x')
  
  check_avg_corr <- augment_melt[, spearman(`.resid`, value), 
                                 by=.(country, variable)] %>%
    .[, list(
      abs_avg_corr=abs(mean(V1, na.rm=T)),
      avg_abs_corr=mean(abs(V1), na.rm=T)), by=variable] %>%
    setorder(cols=-abs_avg_corr, na.last=TRUE)
  
  print(head(check_avg_corr, n=10L))
  
  return(resid_corr_plot)
}
#------ get_link_function ------------------------------------------------------
get_inverse_link_function <- function(fam) {
  inv_link <- switch(
    tolower(fam),
    
    # Gaussian / Identity link
    "gaussian"   = identity,
    "identity"   = identity,
    
    # Poisson, Negative Binomial, Lognormal, Exponential
    "poisson"    = exp,
    "nbinomial"  = exp,
    "lognormal"  = exp,
    "exponential"= exp,
    
    # Binomial
    "binomial"   = plogis,   # inverse logit
    
    # Gamma (inverse link)
    "gamma"      = function(eta) 1 / eta,
    
    # Fallback
    stop(sprintf("Family '%s' not recognized or no inverse link defined.", fam))
  )
    
 return(inv_link) 
}

#------ format_ssn_glm_equation ------------------------------------------------
# in_mod_fit <- tar_read(ssn_mods_miv_invsimpson_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_miv_richness_yr)$ssn_mod_fit

# in_mod_fit <- tar_read(ssn_mods_miv_richness_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_miv_invsimpson_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_dia_sedi_richness_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_dia_sedi_invsimpson_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_dia_biof_richness_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_fun_sedi_richness_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_fun_sedi_invsimpson_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_fun_biof_richness_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_fun_biof_invsimpson_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_bac_sedi_richness_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_bac_sedi_invsimpson_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_bac_biof_richness_yr)$ssn_mod_fit
# in_mod_fit <- tar_read(ssn_mods_bac_biof_invsimpson_yr)$ssn_mod_fit
# 
# greek=T
# format_ssn_glm_equation(in_mod_fit, greek = TRUE)

format_ssn_glm_equation <- function(in_mod_fit, greek = TRUE) {
  if (!inherits(in_mod_fit, c('lm','glm','ssn_lm', 'ssn_glm', 'lme', 'merMod'))) {
    return(NULL)
  }
  
  # Extract response and predictors
  f <- formula(in_mod_fit)
  y <- as.character(f[[2]])
  x <- attr(terms(in_mod_fit), "term.labels")
  
  # Define coefficient names
  coefs <- if (greek) paste0("β_", seq_along(x)) else paste0("b", seq_along(x))
  intercept <- if (greek) "β₀" else "b0"
  
  # Detect family (for SSN2)
  fam <- in_mod_fit$family

  if (!is.null(fam)) {
    link <- switch(
      fam,
      "Gaussian"   = "identity",
      "gaussian"   = "identity",
      "poisson"    = "log",
      "nbinomial"  = "log",
      "binomial"   = "logit",
      "gamma"      = "inverse",
      "lognormal"  = "log",
      "exponential"= "log",
      "identity"   = "identity",
      NA_character_
    )
  } else {
    link <- 'identity'
  }
  
  if (length(y) > 1) {
    if (y[[1]]=='log10') {
      link <- "log"
    }
    y <- y[[2]]
  }
  
  # --- Pretty formatting for variable terms ---
  pretty_terms <- x
  
  # Replace I(var^2) → var²
  pretty_terms <- gsub("I\\(([^\\^]+)\\^2\\)", "\\1²", pretty_terms)
  # Replace I(var^3) → var³
  pretty_terms <- gsub("I\\(([^\\^]+)\\^3\\)", "\\1³", pretty_terms)
  # Replace interactions x1:x2 → x₁ × x₂
  pretty_terms <- gsub(":", " × ", pretty_terms)
  # Optional cleanup for nicer rendering
  pretty_terms <- gsub("\\s+", "", pretty_terms)
  
  # Build the right-hand side
  rhs <- paste(c(intercept, paste0(coefs, "*", pretty_terms)), collapse = " + ")
  
  # Construct the full equation
  if (is.na(link) || link == "identity") {
    eq <- paste0(y, " = ", rhs)
  } else {
    eq <- paste0(link, "(", y, ") = ", rhs)
  }
  
  return(eq)
}

#------ plot_ssn_obs_pred ------------------------------------------------------
plot_ssn_obs_pred <- function(in_mod_fit, 
                              in_drn_dt,
                              response_var_label) {
  
  response_var <- all.vars(in_mod_fit$formula)[[1]]
  
  pred_final <- setDT(augment(in_mod_fit, drop=F, type.predict = 'response')) %>%
    merge(in_drn_dt, by='country')  %>%
    .[, country := factor(
      country,
      levels = c("Finland", "France",  "Hungary", "Czechia", "Croatia", "Spain" ),
      ordered=T)
    ]
  
  color_vec <- pred_final[!duplicated(country),
                          setNames(color, country)]
  
  #Check predictions for final model
  predobs_plot <- ggplot(pred_final, aes(x=.fitted, y=get(response_var))) +
    geom_point(aes(color=country, shape=stream_type), alpha=0.7) +
    # geom_smooth(method = 'gam',
    #             color = 'black', se = FALSE) +
    geom_abline(linetype='dashed') +
    scale_x_continuous(name=paste('Predicted', response_var_label)) +
    scale_y_continuous(name=paste('Observed', response_var_label)) +
    scale_color_manual(name='Country', values=color_vec ) +
    scale_shape_discrete(name='Stream type', 
                         labels = c('Perennial', 'Non-perennial')) +
    # scale_color_manual(
    #   name = 'Stream type',
    #   labels = c('Perennial', 'Non-perennial'),
    #   values=c('#2b8cbe', '#feb24c')) +
    # facet_wrap(~country) +
    coord_fixed() +
    theme_classic()
  
  return(predobs_plot)
}

#------ plot_formula_emtrends ------------------------------------------------
# in_mod_fit <- tar_read(ssn_mods_miv_richness_yr)$ssn_mod_fit
# in_drn_dt = drn_dt
# in_hydro_vars_dt = tar_read(hydro_vars_dt)
# plot=T
# verbose=T

plot_formula_emtrends <- function(in_mod_fit,
                                  in_drn_dt,
                                  in_hydro_vars_dt,
                                  plot = TRUE,
                                  verbose = TRUE) {
  
  mf <- model.frame(in_mod_fit)  # model frame to check variable types
  term_labels <- attr(terms(in_mod_fit), "term.labels")
  if (length(term_labels) == 0) {
    stop("No predictor terms found in model.")
  }
  
  # --- helper: pick the best candidate continuous predictor from vars_in_term ---
  pick_continuous <- function(vars_in_term, mf) {
    candidates <- intersect(vars_in_term, names(mf))
    if (length(candidates) > 0) {
      numeric_cand <- candidates[vapply(mf[candidates], is.numeric, logical(1))]
      if (length(numeric_cand) > 0) return(numeric_cand[[1]])
      return(candidates[[1]])  # fallback if none numeric
    }
    return(vars_in_term[[1]])  # fallback if none in mf
  }
  
  # --- helper: pick a factor/grouping variable (if any) ---
  pick_factor_spec <- function(vars_in_term, mf) {
    candidate_vars <- intersect(vars_in_term, names(mf))
    facs <- candidate_vars[vapply(mf[candidate_vars], function(x)
      is.factor(x) || is.character(x), logical(1))]
    if (length(facs) > 0) return(facs[[1]])
    return(NULL)
  }
  
  results <- list()
  skipped <- list()
  
  for (term in term_labels) {
    vars_in_term <- all.vars(as.formula(paste("~", term)))
    vars_in_term <- setdiff(vars_in_term, all.vars(formula(in_mod_fit))[1])
    
    if (length(vars_in_term) == 0) {
      skipped[[term]] <- "no variables parsed from term"
      if (verbose) message("Skipping term '", term, "': no vars parsed.")
      next
    }
    
    pred_var <- pick_continuous(vars_in_term, mf)
    spec_factor <- pick_factor_spec(vars_in_term, mf)
    
    # skip pure factors (no slope to estimate)
    if (pred_var %in% names(mf) && is.factor(mf[[pred_var]])) {
      skipped[[term]] <- paste0("predictor '", pred_var, "' is factor; no slope to estimate")
      if (verbose) message("Skipping term '", term, "': predictor is a factor.")
      next
    }
    
    specs_formula <- if (!is.null(spec_factor)) {
      as.formula(paste("~", spec_factor))
    } else {
      as.formula("~ 1")
    }
    
    # --- Try emtrends extraction ---
    this_dt <- tryCatch({
      out <- get_ssn_emtrends(
        in_mod = in_mod_fit,
        in_pred_var = pred_var,
        in_pred_var_label = pred_var,
        interaction_var = if (!is.null(spec_factor)) spec_factor else "",
        in_drn_dt = in_drn_dt,
        in_emtrends_dt = NULL,
        plot = FALSE
      )
      dt <- copy(out$dt)
      dt[, `:=`(source_term = term,
                predictor = pred_var,
                spec_factor = ifelse(is.null(spec_factor), NA_character_, spec_factor))]
      dt
    }, error = function(e) {
      # fallback: direct emtrends call
      try_direct <- try({
        emtr <- emtrends(object = in_mod_fit,
                         specs = specs_formula,
                         var = pred_var,
                         data = in_mod_fit$ssn.object$obs)
        dt2 <- as.data.frame(emtr) %>% setDT()
        setnames(dt2, paste0(pred_var, ".trend"), "trend", skip_absent = TRUE)
        dt2[, `:=`(
          response_var = all.vars(formula(in_mod_fit))[[1]],
          pred_var_name = pred_var,
          source_term = term,
          predictor = pred_var,
          spec_factor = ifelse(is.null(spec_factor), NA_character_, spec_factor)
        )]
        dt2
      }, silent = TRUE)
      
      if (inherits(try_direct, "try-error")) {
        if (verbose) message("Failed emtrends for term '", term, "': ", e$message)
        skipped[[term]] <- paste("emtrends failed:", e$message)
        NULL
      } else {
        try_direct
      }
    })
    
    if (!is.null(this_dt) && nrow(this_dt) > 0) {
      results[[term]] <- this_dt
      if (verbose) message("Succeeded for term '", term, "' -> predictor '", pred_var,
                           if (!is.null(spec_factor)) paste0(" | specs: ", spec_factor) else " | specs: ~1")
    } else {
      if (is.null(skipped[[term]])) skipped[[term]] <- "no rows returned"
      if (verbose) message("No rows returned for term '", term, "'.")
    }
  }
  
  if (length(results) == 0) {
    stop("No emtrends results produced for any term. See `skipped` for reasons.")
  }
  
  emtrends_all <- rbindlist(results, fill = TRUE, use.names = TRUE)
  
  # normalize uncertainty column names
  if (all(c("lcl_capped", "ucl_capped") %in% names(emtrends_all))) {
    p_rangemin <- "lcl_capped"; p_rangemax <- "ucl_capped"
  } else {
    p_rangemin <- "asymp.LCL"; p_rangemax <- "asymp.UCL"
  }
  
  # --- get formatted variable names -----------------
  emtrends_all[, pred_var_label := get_full_hydrolabel(in_hydro_vars_dt, 
                                                       pred_var_name),
               by=.I]
  
  # --- Build summary plot -----
  summary_plot <- NULL
  if (plot) {
    color_by <- if ("spec_factor" %in% names(emtrends_all) && any(!is.na(emtrends_all$spec_factor))) {
      possible_group_cols <- setdiff(
        names(emtrends_all),
        c("trend", "asymp.LCL", "asymp.UCL", p_rangemin, p_rangemax,
          "response_var", "pred_var_name", "source_term", "predictor",
          "spec_factor", "SE", "df")
      )
      
      group_col <- NULL
      for (nc in possible_group_cols) {
        if (nc %in% names(in_drn_dt) &&
            any(emtrends_all[[nc]] %in% unique(in_drn_dt[[nc]]), na.rm = TRUE)) {
          group_col <- nc; break
        }
      }
      if (is.null(group_col) && length(possible_group_cols) > 0) {
        group_col <- possible_group_cols[[1]]
      }
      group_col
    } else {
      "source_term"
    }
    
    # --- build color vector ---
    color_vec <- NULL
    if (color_by %in% names(in_drn_dt) && "color" %in% names(in_drn_dt)) {
      color_vec <- in_drn_dt[!duplicated(get(color_by)),
                             setNames(color, get(color_by))]
    } else if (color_by %in% names(emtrends_all) && is.character(emtrends_all[[color_by]])) {
      vals <- unique(emtrends_all[[color_by]])
      color_vec <- setNames(scales::hue_pal()(length(vals)), vals)
    }
    
    # clean up missing colors
    if (!color_by %in% names(emtrends_all)) color_by <- "source_term"
    emtrends_all[is.na(get(color_by)), (color_by) := "All"]
    
    if (exists("color_vec") && is.character(emtrends_all[[color_by]])) {
      if ("All" %in% emtrends_all[[color_by]]) {
        color_vec[["All"]] <- "black"
      }
      missing_cols <- setdiff(unique(emtrends_all[[color_by]]), names(color_vec))
      if (length(missing_cols) > 0) {
        extra_cols <- scales::hue_pal()(length(missing_cols))
        names(extra_cols) <- missing_cols
        color_vec <- c(color_vec, extra_cols)
      }
      if (anyNA(names(color_vec))) color_vec <- color_vec[!is.na(names(color_vec))]
    }
    
    summary_plot <- ggplot(emtrends_all, aes(x = pred_var_label, y = trend)) +
      geom_pointrange(aes_string(ymin = p_rangemin, ymax = p_rangemax, color = color_by),
                      position = position_dodge(width = 0.6), alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = 2) +
      scale_x_discrete(labels=label_wrap_gen(width = 25)) +
      facet_wrap(~ predictor, scales = "free_y", ncol = 1) +
      coord_flip() +
      labs(
        x = "Predictor (term)",
        y = paste("Estimated slope of", all.vars(formula(in_mod_fit))[[1]]),
        color = Hmisc::capitalize(color_by)
      ) +
      theme_minimal() +
      theme(strip.text = element_blank())
    
    if (!is.null(color_vec)) {
      summary_plot <- summary_plot + scale_color_manual(values = color_vec)
    }
  }
  
  return(list(
    dt = emtrends_all,
    plot = summary_plot,
    skipped = skipped
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
    country = c("Croatia", "Czechia", "Finland", "France",  "Hungary", "Spain"),
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
             c('bankfull_at_max_wetted_width_m', 'bankfull_at_min_wetted_width_m')) %>%
    setnames('drn', 'country') %>%
    .[country=='Czech', country := 'Czechia']
  
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
    country=='Hungary' & is.na(avg_velocity_macroinvertebrates) & state_of_flow == 'F',
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
    country=='France' & avg_depth_macroinvertebrates>0 & discharge_l_s>1, 
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
    country=='Hungary' & campaign != 5  & state_of_flow == 'F',
    list(
      site_id = site,
      ratio=conductivity_micros_cm/mean(conductivity_micros_cm, na.rm=T),
      mean_conduct = mean(conductivity_micros_cm, na.rm=T)
    ), by=campaign] %>%
    .[, mean(ratio), by=site_id]
  
  avg_conduct_5 <- env_dt_merged[
    country=='Hungary' & campaign == 5  & state_of_flow == 'F', 
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
  
  #Correct stream time
  env_dt_merged[site=='BUK23', stream_type := 'TR']
  
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
# in_miv_full_dt = tar_read(miv_full_dt)
# include_bacteria = T

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
#' - Removes pool samples (keeps separate `_nopools` tables).
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
      out_dt[country == 'Czech Republic', country := 'Czechia']
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
  
  #Remove Nematoda and Nematomorpha phyla, and Hydrachnidae family
  dt_list$miv <- dt_list$miv %>%
    select(-c(hydrachnidae.gen..sp., nematoda.gen..sp., nematomorpha.gen..sp.))
  dt_list$miv_nopools <- dt_list$miv_nopools %>%
    select(-c(hydrachnidae.gen..sp., nematoda.gen..sp., nematomorpha.gen..sp.))

  
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
  #   new = c("Croatia", 'France', "spain", "Czechia", "Hungary", 'Finland')
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
    
    
    if (org=='fun') {
      print('Correcting fungi dates')
      dt_list[[org_medium]][date == as.Date("2012-02-25"), 
                            date := as.Date("2021-02-25")]
    }
    
    dt_list[[new_org_medium]] <- dt_list[[org_medium]] %>%
      .[country %in% c('Czech Republic', 'Czech'), country := 'Czechia'] %>%
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
# in_country <- 'Czech'
# in_biodt <- tar_read(bio_dt)[['dia_biof']][country == in_country,]
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
  
  #Compute and partition taxonomic gamma diversity -----------------------------
  #between alpha and beta for the overall period, for individual sites
  #and for the metacom
  spxp <- as.matrix(in_biodt[, ..spcols])
  bio_structure <- data.frame(space = as.factor(in_biodt$site),
                              time = as.factor(in_biodt$campaign))
  
  #Only compute decomp if there is variation in the species matrix
  all_sites <- unique(in_biodt$site)
  
  if (any(rowSums(spxp > 0) > 0) && 
      ncol(spxp[, colSums(spxp) > 0, drop=FALSE]) >= 2) {
    decomp <- try(HierAnodiv(spxp = spxp[rowSums(spxp > 0) > 0,], 
                             structure = bio_structure[rowSums(spxp > 0) > 0,], 
                             phy = NULL, weight = NULL, check = TRUE, q = 1)) 
    localdiv_decomp <- as.data.table(decomp$tab_by_site) %>%
      setnames(c('nsite'), 'ncampaigns')
    localdiv_decomp$site <- unique(as.character(bio_structure[rowSums(spxp > 0) > 0,]$space))
    
    # Add missing sites with NA
    missing_sites <- setdiff(all_sites, localdiv_decomp$site)
    if (length(missing_sites) > 0) {
      na_rows <- data.table(
        site = missing_sites,
        ncampaigns = NA,
        Gamma = 0,
        Beta = 0,
        mAlpha = 0
      )
      localdiv_decomp <- rbindlist(list(localdiv_decomp, na_rows), use.names = TRUE, fill = TRUE)
    }
    
  } else {
    # No variation -> all sites with NA
    localdiv_decomp <- data.table(
      site = all_sites,
      ncampaigns = NA,
      Gamma = 0,
      Beta = 0,
      mAlpha = 0
    )
  }
  
  if (level == 'local') {
    #Compute species richness --------------------------------------------------
    biodt_div <-in_biodt[, list(richness = rowSums(.SD > 0)
    ), by=.(site, campaign, date, organism), .SDcols = spcols] %>%
      .[, mean_richness := mean(richness), by=site]
    
    biodt_copy <- copy(in_biodt)
    
    #Compute exp shannon and inverse Simpson (alpha diversity) by site x time step -----------------
    sha_simp <- biodt_copy[, list(
      campaign = campaign,
      site = site,
      expshannon = ifelse(rowSums(.SD) == 0, 0, exp(vegan::diversity(as.matrix(.SD), index = "shannon"))),
      invsimpson = ifelse(rowSums(.SD) == 0, 0, vegan::diversity(as.matrix(.SD), index = "invsimpson"))
    ), .SDcols = spcols] %>%
      .[is.infinite(invsimpson), invsimpson := 0]  %>%
      .[, `:=`(
        mean_expshannon = mean(expshannon),
        mean_invsimpson = mean(invsimpson)
      ), by=site]
    
    #Compute nestedness and turnover based on temporal beta diversity-----------
    # Identify sites with multiple campaigns
    multicampaign_sites <- biodt_copy[, .N, by = site][N > 1, site]
    
    # Number of species per site
    site_stats <- biodt_copy[site %in% multicampaign_sites,
                             .(n_species = sum(colSums(.SD) > 0)), by = site, .SDcols = spcols]
    
    valid_sites <- site_stats[n_species >= 2, site]
    zero_sites  <- site_stats[n_species < 2, site]
    
    # Compute comp_richrepl_inner for valid sites
    J_list <- lapply(valid_sites, function(s) {
      sub_dt <- biodt_copy[site == s]
      out <- comp_richrepl_inner(sub_dt, spcols = spcols, beta_div_coef = 'J', quant = FALSE)
      out[, site := s]
      out
    })
    R_list <- lapply(valid_sites, function(s) {
      sub_dt <- biodt_copy[site == s]
      out <- comp_richrepl_inner(sub_dt, spcols = spcols, beta_div_coef = 'J', quant = TRUE)
      out[, site := s]
      out
    })
    
    # Combine valid sites
    J_valid <- rbindlist(J_list, use.names = TRUE, fill = TRUE)
    R_valid <- rbindlist(R_list, use.names = TRUE, fill = TRUE)
    
    # Fill zeros for zero-species sites using same column names as valid sites
    if (length(zero_sites) > 0) {
      zero_template_J <- J_valid[0]  # empty DT with correct column names
      J_zero <- rbindlist(lapply(zero_sites, function(s) {
        sub_dt <- biodt_copy[site == s, .(site, campaign)]
        zero_dt <- zero_template_J[rep(1, nrow(sub_dt))]
        zero_dt[, names(zero_dt) := 0]
        zero_dt[, site := sub_dt$site]
        zero_dt[, campaign := sub_dt$campaign]
        zero_dt
      }))
      
      zero_template_R <- R_valid[0]
      R_zero <- rbindlist(lapply(zero_sites, function(s) {
        sub_dt <- biodt_copy[site == s, .(site, campaign)]
        zero_dt <- zero_template_R[rep(1, nrow(sub_dt))]
        zero_dt[, names(zero_dt) := 0]
        zero_dt[, site := sub_dt$site]
        zero_dt[, campaign := sub_dt$campaign]
        zero_dt
      }))
      
      # Combine
      J_richrepl <- rbindlist(list(J_valid, J_zero), use.names = TRUE, fill = TRUE)
      R_richrepl <- rbindlist(list(R_valid, R_zero), use.names = TRUE, fill = TRUE)
    } else {
      J_richrepl <- J_valid
      R_richrepl <- R_valid
    }
    
    # Optional: sort
    setorderv(J_richrepl, "site")
    setorderv(R_richrepl, "site")
    
    
    
    #Merge all metrics ---------------------------------------------------------
    out_dt <- mergeDTlist(
      list(biodt_div, sha_simp, J_richrepl, R_richrepl), 
      by=c('site', 'campaign'), all=T, sort = T, set_suffix=F) %>%
      merge(localdiv_decomp, by='site', all.x=T)
    
    #Remove GEN04_1 from all organisms, no local environmental data, sampled only for eDNA. Too unsure.
    out_dt <- out_dt[!(site == 'GEN04' & campaign == '1'),]
    #Create organism class for coloring/ordering
    out_dt[, organism_class := gsub('_[a-z]+', '', organism)]
    
  } else if (level == 'regional') {
    nsites <- if(!is.null(decomp)) nrow(decomp$tab_by_site) else 1
    nsteps <- if(!is.null(decomp)) max(decomp$tab_by_site[,'nsite']) else 1
    
    if (!is.null(decomp)) {
      out_dt <- t(decomp[[1]]) %>%
        data.table %>%
        setnames(c('Gamma', 'Beta1', 'Beta2in1', 'mAlpha'),
                 paste0(c('gamma', 'beta_s', 'beta_t_s', 'malpha'), '_drn')
        ) %>% 
        .[, `:=`(organism = in_biodt[1, .(organism)][[1]],
                 beta_s_drn_std = (beta_s_drn-1)/(nsites-1),
                 beta_t_s_drn_std = (beta_t_s_drn-1)/(nsteps-1) 
        )] 
    } else {
      #No species -> all zeros
      out_dt <- data.table(
        gamma_drn = 0,
        beta_s_drn = 0,
        beta_t_s_drn = 0,
        malpha_drn = 0,
        organism = in_biodt[1, .(organism)][[1]],
        beta_s_drn_std = 0,
        beta_t_s_drn_std = 0
      )
    }
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
              by=.(campaign, country, organism)],
    in_envdt[, list(per_flowing = 100*.SD[state_of_flow=='F', .N]/.N),
             by=.(campaign, country)],
    by=c('campaign', 'country')) %>%
    .[, mean_drn_richness_relative := mean_drn_richness/max(mean_drn_richness),
      by=.(country, organism)]
  
  ggplot(sprich_hydroobs[
    !(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=campaign, y=mean_drn_richness_relative)) +
    geom_line(aes(group=country, color=country), size=1.2) +
    scale_color_manual(values=c("#8c510a", "#bf812d", "#01665e", "#80cdc1", "#8073ac", "#543005")) +
    new_scale_color() +
    geom_point(aes(color=per_flowing), size=2) +
    scale_color_distiller(palette='Spectral', direction=1) +
    facet_wrap(~organism, scales='free_y') +
    theme_classic()
  
  ggplot(sprich_hydroobs[
    !(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=per_flowing, y=mean_drn_richness_relative, color=country)) +
    geom_point() +
    geom_smooth(method='lm', se=F) +
    scale_color_manual(values=c("#8c510a", "#bf812d", "#01665e", "#80cdc1", "#8073ac", "#543005")) +
    facet_wrap(~organism, scales='free_y')
  
  overall_lines <- ggplot(
    in_sprich[!(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=campaign, y=richness)) +
    #geom_point() +
    geom_line(aes(group=site, color=site)) +
    geom_smooth(aes(group=country)) +
    facet_wrap(country~organism, scales = 'free_y') +
    theme(legend.position='none')
  
  ggplot(
    in_sprich[!(organism %in% c('miv', 'miv_nopools_flying', 'miv_nopools_nonflying')),],
    aes(x=campaign, y=richness)) +
    #geom_point() +
    geom_boxplot(aes(group=site, color=state_of_flow)) +
    facet_wrap(country~organism, scales = 'free_y') 
  
  ggplot(
    in_sprich[(organism %in% c('miv_nopools')),],
    aes(x=campaign, y=richness)) +
    #geom_point() +
    geom_line(aes(group=site, color=site)) +
    geom_smooth(aes(group=country)) +
    facet_wrap(country~organism, scales = 'free_y') +
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
  rivnet_clustered_joinini <- rivnet_clustered_joinini[
    as.numeric(st_length(rivnet_clustered_joinini)) > 1,]
  
  #Convert those without match to points every 10 meters
  rivnet_nomatch_lines <- rivnet_clustered_joinini %>%
    .[is.na(rivnet_clustered_joinini[[idcol]]),] 
  
  if (nrow(rivnet_nomatch_lines) > 0) {
    rivnet_nomatch_pts <- rivnet_nomatch_lines %>%
      st_line_sample(density = 1/5) %>%
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
    
    rivnet_clustered_joinall <- merge(rivnet_clustered_joinini,
                                      rivnet_nomatch_selid, 
                                      by='UID', all.x=T)
  } else {
    rivnet_clustered_joinall <- rivnet_clustered_joinini
  }
  
  #Fill NAs with those
  
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
# in_country <- 'France'
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
  } else if (country == 'Czechia') {
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
#------ impute_hydromod --------------------------------------------------------
# in_country <- 'Hungary'
# varname <- 'isflowing' #'qsim' #
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[in_country]]
# in_hydromod_drn <- tar_read_raw((paste0('hydromod_hist_dt_', in_country, '_', varname)))
# in_network_idcol = 'cat_cor'

impute_hydromod <- function(in_network_path,
                            varname,
                            in_hydromod_drn,
                            in_network_idcol = 'cat_cor') {
  setDT(in_hydromod_drn$data_all)
  setDT(in_hydromod_drn$dates_format) 
  
  #Import network shapefiles and get their length
  network_v <- terra::vect(in_network_path)
  network_v$reach_length <- terra::perim(network_v)
  reach_dt <- as.data.table(network_v) %>%
    setnames(in_network_idcol, 'reach_id', skip_absent=TRUE)
  remove(network_v)
  
  #Merge hydro data and reach length for computing network-wide statistics
  #qdat[, unique(reach_id)[!(unique(reach_id) %in% reach_dt$reach_id)]] 
  #Two IDs are not in the shapefile for Finland? maybe a hydrological unit not associated
  hydromod_sub <- in_hydromod_drn$data_all[reach_id %in% 
                                             unique(reach_dt$reach_id),]
  
  #For those segments that are in the network but that were not modeled,
  # get the intermittence status of the upstream segment with the largest
  # drainage area
  #(i.e.: their ID exists topologically, and they were corrected to be correct
  # but they were not modeled, so they have no discharge or intermittence data,
  unmodeled_reaches_dt <- reach_dt[!(reach_id %in% 
                                       unique(hydromod_sub$reach_id)),]
  
  if (varname == 'isflowing') {
    
    unmodeled_hydromod <- reach_dt[
      to_reach_hydromod %in% unmodeled_reaches_dt$reach_id, 
      list(
        largest_upstream_cat_cor = .SD[which.max(upstream_area_net), reach_id]
      ),
      by=to_reach_hydromod] %>%
      merge(hydromod_sub, by.x='largest_upstream_cat_cor', by.y='reach_id') %>%
      .[, largest_upstream_cat_cor := NULL] %>%
      setnames('to_reach_hydromod', 'reach_id')
    
  } else  if (varname == 'qsim') {
    
    unmodeled_hydromod <- in_hydromod_drn$data_all[
      reach_id %in% reach_dt[to_reach_hydromod %in% unmodeled_reaches_dt$reach_id, reach_id], ] %>%
      merge(reach_dt[, .(reach_id, to_reach_hydromod)], by='reach_id', all.x=T) %>%
      .[to_reach_hydromod %in% unmodeled_reaches_dt$reach_id, 
        .(qsim = sum(get(varname), na.rm=TRUE)),
        by=.(to_reach_hydromod, date)] %>%
      setnames('to_reach_hydromod', 'reach_id')
    
  }
  
  #Bind back to modeled data and format dates
  hydromod_filled <- rbind(hydromod_sub, unmodeled_hydromod) %>%
    merge(in_hydromod_drn$dates_format[, .(date, month, hy, doy)],
          by='date', all.x=T, all.y=F)
  
  return(hydromod_filled)
  # ggplot(hydromod_filled[reach_id %in% c('2424600', '2490400',
  #                                         '2504000', '2504200',
  #                                         '2424400', '2503800') &
  #                          date < as.Date('1970-01-01'),]) +
  #   geom_line(aes(x=date, y=qsim, color=as.factor(reach_id))) +
  #   scale_y_log10()
}
#------ compute_hydrostats_drn -------------------------------------------------
# in_country <- 'Hungary'
# varname <- 'qsim' #'isflowing'
# in_sites_dt <- tar_read(sites_dt)[country == in_country,]
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[in_country]]
# in_hydromod_filled <- tar_read_raw((paste0('hydromod_hist_filled_dt_', in_country, '_', varname)))
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
                                   in_hydromod_filled,
                                   in_network_idcol = 'cat') {
  #-------------------- Compute intermittence statistics -----------------------
  if (varname == 'isflowing') {
    #Import network shapefiles and get their length
    network_v <- terra::vect(in_network_path)
    network_v$reach_length <- terra::perim(network_v)
    reach_dt <- as.data.table(network_v[, c(in_network_idcol, 'reach_length')]) %>%
      setnames(in_network_idcol, 'reach_id', skip_absent=T) %>%
      .[, list(reach_length = sum(reach_length)), by=reach_id]
    remove(network_v)
    
    #Merge hydro data and reach length for computing network-wide statistics
    #qdat[, unique(reach_id)[!(unique(reach_id) %in% reach_dt$reach_id)]] 
    #Two IDs are not in the shapefile for Finland? maybe a hydrological unit not associated
    intermod_dt <- merge(in_hydromod_filled, reach_dt, 
                         by='reach_id', all.x=T, all.y=F)
    
    if ('nsim' %in% names(intermod_dt)) {
      q_stats <- list()
      #Compute hydrological statistics at the DRN scale (Relative flowing length)
      q_stats$country <- intermod_dt[, 
                                     compute_hydrostats_intermittence(
                                       in_hydromod_dt = .SD,
                                       in_sites_dt = in_sites_dt,
                                       scale = 'drn')$country, 
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
    hydromod_dt_sites <- in_hydromod_filled[
      reach_id %in% unique(in_sites_dt$reach_id),] %>%
      setorderv(c('reach_id', 'date'))
    
    q_stats <- compute_hydrostats_q(in_hydromod_dt = hydromod_dt_sites) %>%
      merge(in_sites_dt[, .(site, reach_id)], .,
            by = 'reach_id', all.x = T, all.y = F,
            allow.cartesian = T) %>%
      setorderv(c('reach_id', 'date'))
  }
  
  return(q_stats)
}

#------ subset_hydrostats ------------------------------------------------------
# in_country <- in_drn <- 'Czechia'
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
# in_country <- 'Czechia'
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
  
  #Not needed anymore: corrected the input table.
  # if (in_country == 'Czechia') {
  #   sites_dt[, id := sub('S', '', id) %>%
  #              str_pad(width=2, side='left', pad = 0) %>%
  #              paste0('VEL', .)]
  # } 
  
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
  sites_dt_filled <- in_env_dt[country==in_country, .(site, latitude, longitude)] %>%
    setnames(c('latitude', 'longitude'), c('latitude_field', 'longitude_field')) %>%
    unique %>%
    merge(sites_dt, by.x='site', by.y='id', all.x=T) %>%
    .[is.na(lat), `:=`(lat = latitude_field, lon = longitude_field)] %>%
    .[, `:=`(latitude_field = NULL, longitude_field = NULL)] %>%
    .[site != 'BUT08',] #BUT08 is not represented in river network. No hydrological information
  
  
  sites_dt_filled[, reach_id := fcase(
    site=='BUK01', 657801L,
    site=='BUK36', 665000L,
    site=='GEN05', 5780L,
    site=='GEN11', 5594L,
    site=='GEN12', 5588L,
    site=='GEN17', 5882L,
    site=='GEN23', 5732L,
    site=='GEN26', 5392L,
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
# country <- 'Croatia'
# in_sites_path <- tar_read(site_points_gpkg_list)[[country]]
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[country]]
# out_snapped_sites_path = NULL
# overwrite = T
# custom_proj = F
# in_sites_unique_id = 'site'
# in_network_unique_id = 'UID'
# in_sites_idcol_tomatch = 'reach_id'
# in_network_idcol_tomatch = 'cat_cor'
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
# country <- 'France'
# in_sites_path <- tar_read(barrier_points_gpkg_list)[[country]]
# in_network_path <- tar_read(network_ssnready_gpkg_list)[[country]]
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
# in_country <- in_drn <- 'Czechia'
# in_hydromod_drn <- tar_read(hydromod_comb_hist)[[paste0("hydromod_hist_filled_dt_", in_country, '_isflowing')]]
# in_net_shp_path <- tar_read(network_ssnready_shp_list)[[in_country]]

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
  setDT(in_hydromod_drn) %>%
    setnames('reach_id', 'cat')
  
  # Create the "End_point" site that will correspond to the 111111 in the flow_intermittence dataset
  # this point will be added with the same frequency of any other reach
  end_point_interm <- in_hydromod_drn %>%
    .[cat == .[1, cat],] %>% #select whatever reach
    .[, `:=`(UID = 11111,  from = outlet_to, value = 1)] %>% #format
    .[, .(date, UID, isflowing, from)] 
  
  # Merge shapefile with flow intermittence data
  # Merge the Endpoint site and "pivot_wide" the table to obtain the TRUE intermittence table, 
  # where each row corresponds to a day (dates as factors) and columns to all nodes of the network.
  nsims <- isTRUE('nsim' %in% names(in_hydromod_drn))
  
  cast_formula <- if (nsims) {as.formula('date + nsim ~ from')} else {'date ~ from'}
  sites_status_matrix <- merge(net_dt, in_hydromod_drn,
                               by='cat', all.x=T, 
                               allow.cartesian=TRUE) %>% #To allow UIDs that have the same cat
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
# in_country <- 'Croatia'
# in_preformatted_data <- tar_read(preformatted_data_STcon)[[in_country]]
# in_bio_dt <- tar_read(bio_dt)
# in_nsim <- NULL#tar_read(hydromod_paths_dt)[country == in_country,]$best_sim

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
# in_country <- 'Czechia'
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


#------ compute_connectivity_proj ----------------------------------------------
# hydroref_path <- file.path(
#   wp1_data_gouv_dir,
#   "projections",
#   "Albarine_flowstate_projection_gfdl-esm4_ssp370_1985-2014_spatially-distributed.nc")
# 
# hydroproj_path <- file.path(
#   wp1_data_gouv_dir,
#   "projections",
#   "Albarine_flowstate_projection_gfdl-esm4_ssp585_2015-2100_spatially-distributed.nc")
# 
# tar_load(network_ssnready_shp_list)
# 
# subset_sites = T
# in_sites_dt <- tar_read(sites_dt)
# in_drn_dt <- drn_dt
# min_date = as.Date('1990-01-01')
# STcon_window=10

compute_connectivity_proj <- function(hydroproj_path, hydroref_path, in_drn_dt,
                                  network_ssnready_shp_list,
                                  subset_sites, in_sites_dt, 
                                  STcon_window=365,
                                  min_date = as.Date('1990-01-01')) {
  
  #Format data
  hydromod_formatted <- import_hydromod_gcm(hydroproj_path, 
                                            hydroref_path,
                                            in_drn_dt)
  
  hydro_dt <- list()
  hydro_dt <- hydromod_formatted$hydro_dt
  
  net_sf <- network_ssnready_shp_list[[hydromod_formatted$metadata_dt$country]]
  
  formatted_data <- prepare_data_for_STcon(in_hydromod_drn=hydro_dt, 
                                           in_net_shp_path = net_sf)
  
  #Get December 31st of every year
  in_dates <- hydro_dt[, .(date = date[.N]), 
                       by = format(date, "%Y")] %>%
    .[, .(date)]
  
  STcon_directed <- list()
  STcon_directed$STcon_m365 <- compute_STcon_rolling(
    in_preformatted_data = formatted_data,
    ref = FALSE,
    #in_nsim = hydromod_paths_dt[country == in_country,]$best_sim,
    in_dates = in_dates,
    window = STcon_window,
    output = 'STcon',
    direction = 'directed',
    routing_mode = 'in',
    weighting = FALSE,
    rounding_factor = 1)
  
  # STcon_directed_ref <- list()
  # STcon_directed_ref[[paste0('STcon_m', in_window)]] <- compute_STcon_rolling(
  #   in_preformatted_data = formatted_data,
  #   ref = TRUE,
  #   #in_nsim = hydromod_paths_dt[country == in_country,]$best_sim,
  #   in_dates = in_dates[1,],
  #   window = in_window,
  #   output = in_output,
  #   direction = 'directed',
  #   routing_mode = 'in',
  #   weighting = FALSE,
  #   rounding_factor = 1)
  
  STcon_directed_formatted <- postprocess_STcon(
    in_STcon = STcon_directed,
    in_net_shp_path = net_sf,
    standardize_STcon = FALSE, 
    in_STcon_ref = NULL)
  
  if (subset_sites) {
    STcon_directed_formatted$STcon_dt <- merge(STcon_directed_formatted$STcon_dt,
                                               st_drop_geometry(net_sf[,c('UID', 'cat')]),
                                               by='UID') %>%
      merge(in_sites_dt[, .(reach_id, site)], by.x='cat', by.y='reach_id', all.x=F)
  }
  
  setnames(STcon_directed_formatted$STcon_dt, 'stcon_value', 'STcon_directed_mean_yr')
  STcon_directed_formatted$STcon_dt[, `:=`(year = format(date, '%Y'),
                                           variable = NULL)]
  
  #Compute STcon_undirected
  STcon_undirected <- list()
  STcon_undirected$STcon_m365 <- compute_STcon_rolling(
    in_preformatted_data = formatted_data,
    ref = FALSE,
    #in_nsim = hydromod_paths_dt[country == in_country,]$best_sim,
    in_dates = in_dates,
    window = STcon_window,
    output = 'STcon',
    direction = 'undirected',
    routing_mode = 'in',
    weighting = FALSE,
    rounding_factor = 1)
  
  STcon_undirected_formatted <- postprocess_STcon(
    in_STcon = STcon_undirected,
    in_net_shp_path = net_sf,
    standardize_STcon = FALSE, 
    in_STcon_ref = NULL)
  
  if (subset_sites) {
    STcon_undirected_formatted$STcon_dt <- merge(STcon_undirected_formatted$STcon_dt,
                                               st_drop_geometry(net_sf[,c('UID', 'cat')]),
                                               by='UID') %>%
      merge(in_sites_dt[, .(reach_id, site)], by.x='cat', by.y='reach_id', all.x=F)
  }
  
  setnames(STcon_undirected_formatted$STcon_dt, 'stcon_value', 'STcon_undirected_mean_yr')
  STcon_undirected_formatted$STcon_dt[, `:=`(year = format(date, '%Y'),
                                           variable = NULL)]
  
  #Compute Fdist undirected
  Fdist_undirected_proj <- compute_Fdist(
    sites_status_matrix = formatted_data$sites_status_matrix,
    network_structure = formatted_data$network_structure, 
    routing_mode = 'all', 
    raw_dist_matrix = formatted_data$river_dist_mat, 
    in_net_shp_path = net_sf
  ) %>%
    .[is.infinite(Fdist), Fdist := .[!is.infinite(Fdist), 2*max(Fdist)]] %>% #Replace infinite Fdist with twice max otherwise
    .[UID %in% STcon_directed_formatted$STcon_dt[, unique(UID)],] %>% 
    .[, list(Fdist_undmean_yr = mean(Fdist)), by = .(UID, year=format(date, '%Y'))]
  
  
  out_dt <- merge(STcon_directed_formatted$STcon_dt,
                  STcon_undirected_formatted$STcon_dt,
                  by=c('UID', 'site', 'cat','date', 'year')) %>%
    merge(Fdist_undirected_proj,
          by=c('UID', 'year')) %>%
    cbind(hydromod_formatted$metadata_dt[, .(country, gcm, scenario)])
  
  return(out_dt)
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
# in_hydrocon_compiled <- tar_read(hydrocon_sites_compiled)


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
  hydrocon_summarized <- in_hydrocon_compiled[
    #(date >= min(date_range)) &   (date <= max(date_range))
    , list( #`:=`
      upstream_area_net = .SD[1, upstream_area_net],
      DurD_samp = sum(isflowing==0)/.N, #Total number of no-flow days 
      DurD3650past = .SD[.N, DurD3650past],
      PDurD365past = .SD[.N, PDurD365past],
      FreD_samp = length(.SD[!is.na(noflow_period), unique(noflow_period)])/.N, #Total number of no-flow periods per day
      FreD3650past = .SD[.N, FreD3650past],
      PFreD365past = .SD[.N, PFreD365past],
      DurD_max_samp = nafill(as.integer(max(noflow_period_dur, na.rm=T)), fill=0), #Maximum duration of no-flow period
      DurD_avg_samp = nafill(.SD[!duplicated(noflow_period), mean(noflow_period_dur, na.rm=T)], fill=0), #AVerage duration of no-flow period (same as DurD/FreD)
      DurD_yr = .SD[1, DurD_yr],
      FreD_yr = .SD[1, FreD_yr],
      FstDrE = .SD[1, FstDrE],
      meanConD_yr = .SD[1, meanConD_yr],
      DurD_CV10yrpast = .SD[.N, DurD_CV10yrpast],
      DurD_CV30yrpast = .SD[.N, DurD_CV30yrpast],
      FreD_CV10yrpast = .SD[.N, FreD_CV10yrpast],
      FreD_CV30yrpast = .SD[.N, FreD_CV30yrpast],
      meanConD_CV10yrpast = .SD[.N, meanConD_CV10yrpast],
      meanConD_CV30yrpast = .SD[.N, meanConD_CV30yrpast],
      FstDrE_SD10yrpast = .SD[.N, FstDrE_SD10yrpast],
      FstDrE_SD30yrpast = .SD[.N, FstDrE_SD30yrpast],
      sd6_10yrpast = .SD[.N, sd6_10yrpast],
      sd6_30yrpast = .SD[.N, sd6_30yrpast], 
      RelF_avg_samp = mean(relF),
      RelF_min_samp = min(relF),
      relF3650past = .SD[.N, relF3650past],
      relF7mean_yrmin_cv10yrpast = .SD[.N, relF7mean_yrmin_cv10yrpast],
      relF7mean_yrmin_cv30yrpast = .SD[.N, relF7mean_yrmin_cv30yrpast],
      meanQ3650past = .SD[.N, meanQ3650past],
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
  
  return(hydrocon_summarized)
}

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
    paste0('hydromod_hist_filled_dt_', in_country, '_qsim')]]
  isflowing_country <- in_hydromod[[
    paste0('hydromod_hist_filled_dt_', in_country, '_isflowing')]]
  
  #Link q data - keep only full hydrological years, 
  #Exclude 2022 because includes period after sampling
  
  #Calculate average discharge over the full historical period
  qsim_all_dt <- qsim_country[
    (date >= min(in_all_date_range)) &
      (date < max(in_all_date_range)),  
    list(qsim_avg = mean(qsim, na.rm=T)), 
    by=reach_id] 
  
  #Calculate average discharge over the sampling period
  qsim_samp_dt <- qsim_country[
    (date >= min(in_samp_date_range)) &
      (date <= max(in_samp_date_range)),  
    list(qsim_avg_samp = mean(qsim, na.rm=T)), 
    by=reach_id] 
  
  #Calculate flow duration (proportion of dry days) over the sampling period
  isflowing_samp_dt <- isflowing_country[
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
# in_genal_upa <- tar_read(genal_sites_upa_dt)

#' @title Summarize environmental data
#' @description This function takes raw environmental data and computes the mean 
#'     value for a set of variables, grouped by drainage basin, site, and stream type. 
#'     It also generates a boxplot for visualization.
#' @param in_env_dt A data.table containing raw environmental data with a 'state_of_flow' column.
#' @return A list containing two summarized data tables 
#'     (all flows, and only flowing sites) and a ggplot object of the boxplot.
summarize_env <- function(in_env_dt,
                          in_genal_upa) {
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
  group_cols <- c('country', 'site', 'stream_type')
  exclude_cols <- c('if_ip_number_and_size_2_axes_+_depth_of_the_pools',
                    'campaign', 'date', 'state_of_flow')
  
  # identify data columns by excluding grouping and exclusion columns
  dat_cols <- setdiff(names(in_env_dt), 
                      c(group_cols, exclude_cols, 'running_id'))
  
  #str(in_env_dt[, dat_cols, with=F])
  
  # Fill basin area NAs for Genal basin
  in_env_dt[is.na(basin_area_km2),  
            basin_area_km2 := in_genal_upa[.SD, on="site", x.basin_area_km2]]
  
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
                             aes(x=country, y=value, fill=country)) +
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
# 
# hydroref_path <- file.path(
#   wp1_data_gouv_dir,
#   "projections",
#   "Albarine_flowstate_projection_gfdl-esm4_ssp370_1985-2014_spatially-distributed.nc")

# hydroproj_path <- file.path(
#     wp1_data_gouv_dir,
#     "projections",
#     "Albarine_discharge_projection_gfdl-esm4_ssp585_2015-2100_spatially-distributed.nc")

# subset_sites = T
# in_sites_dt <- tar_read(sites_dt)
# varname = 'flowstate'
# in_drn_dt <- drn_dt
# min_date = as.Date('1990-01-01')
# include_metadata = FALSE

#' @title Summarize drainage basin hydrological projection stats
#' @description This function processes a NetCDF file containing hydrological 
#' projection data. It extracts metadata, reads flow state data, and computes a 
#' yearly summary of flow duration.
#' @param hydroproj_path A character string specifying the path to the NetCDF file.
#' @return A data.table containing summarized flow state statistics by year
#'  for each reach, along with associated metadata.
summarize_drn_hydroproj_stats <- function(hydroproj_path,
                                          hydroref_path,
                                          subset_sites=TRUE,
                                          in_sites_dt=NULL,
                                          in_drn_dt=NULL,
                                          min_date = as.Date('1990-01-01'),
                                          include_metadata = TRUE) {
  

  
  hydromod_preformatted <- import_hydromod_gcm(hydroproj_path, 
                                               hydroref_path, in_drn_dt) 
  
  metadata_dt <- hydromod_preformatted$metadata_dt
  hydro_dt <-   hydromod_preformatted$hydro_dt
  
  # Subset reaches to those where there are sampling sites ---------------------
  if (subset_sites && (!is.null(in_sites_dt))) {
    #Prepare site data for subsetting hydrological data
    sites_reach_id <- in_sites_dt[country %in% metadata_dt$country, reach_id]
    hydro_dt <- hydro_dt[reach_id %in% sites_reach_id,]
  }
  
  #Add 30 years of NAs at the beginning to avoid NAs at the beginning ----------
  hydro_dt <- rbind(
    hydro_dt[, list(date=seq(min(date)-lubridate::years(30), 
                                      min(date)-lubridate::days(1), 
                                      by="day")),
                      by=reach_id],
    hydro_dt,
    fill=T) %>%
    setorderv(c('reach_id', 'date'))

  #Compute yearly discharge statistics -----------------------------------------
  if (metadata_dt$varname == 'discharge') {
    metadata_dt[varname=='discharge', varname:='qsim']
    
    stats_dt_yr <- hydro_dt[, year := year(date)] %>%
      .[(date >= min_date - lubridate::years(30)) 
        & (date < as.Date('2100-01-01')),
        list(meanQ = mean(qsim, na.rm = TRUE)),
        by = .(reach_id, year)] %>%
      .[ order(reach_id, year),
         meanQ3650past := frollmean(meanQ, 10, na.rm = TRUE, align="right"),
         by = reach_id] %>%
      .[year >= year(min_date),]
  }
  
  # compute yearly drying duration stats ---------------------------------------
  if (metadata_dt$varname == 'flowstate') {
    metadata_dt[varname=='flowstate', varname:='isflowing']
    
    stats_dt <- hydro_dt[(date >= min_date-lubridate::years(30)) &
                                    (date < as.Date('2100-01-01')),]
    stats_dt[, `:=`(doy  = as.numeric(format(date,"%j")),
                    year =  as.integer(format(date, '%Y'))
                    )]
    
    #Compute duration of no-flow periods and time since last no-flow period #PrdD: prior days to last dry/pool/flowing event
    preformat_intermittence_stats(stats_dt)
    
    #Compute annual statistics of duration, freaquency, timing and variability----------
    stats_dt_yr <- stats_dt[, list(
      DurD_yr = sum(isflowing==0, na.rm=T)/.SD[!is.na(isflowing), .N],
      FreD_yr = length(unique(noflow_period[!is.na(noflow_period)]))/.SD[!is.na(isflowing), .N],
      # Store unique noflow_period IDs per year (for 10-year rolling union later)
      noflow_periods = list(unique(noflow_period[!is.na(noflow_period)])),
      n_valid = .SD[!is.na(isflowing), .N],
      FstDrE = .SD[isflowing==0, min(doy, na.rm=T)],
      meanConD_yr = .SD[!duplicated(noflow_period), mean(noflow_period_dur, na.rm=T)]             
    ), by=.(reach_id, year)] %>%
      .[is.infinite(FstDrE), FstDrE := 365] %>% #For sites that don't dry, set it to 365
      .[is.na(meanConD_yr), meanConD_yr := 0] %>%
      setorder(reach_id, year)
    
    #--- Compute FreD3650past (unique noflow_periods over the past 10 years) ------
    stats_dt_yr[
      , FreD3650past :=
        # Apply over each row: union of the current and previous 9 years' lists
        sapply(seq_len(.N), function(i) {
          length(unique(unlist(noflow_periods[pmax(1, i-9):i])))/sum(n_valid[pmax(1, i-9):i])
        }),
      by = reach_id]
    
    #Drop the list column to keep the table tidy
    stats_dt_yr[, `:=`(noflow_periods =NULL,
                       n_valid=NULL)]
    
    #--- Compute ECDFs  ------------------------------------------------------
    stats_dt_yr <- compute_ecdf_multimerge(
      in_dt = stats_dt_yr,
      ecdf_columns = c("DurD_yr", "FreD_yr"),
      grouping_columns = "reach_id",
      keep_column = "year",
      na.rm = TRUE
    )

    #CV of annual number of no-flow days----------------------------------------
    rollingstep_yr <- c(10, 30)
    stats_dt_yr[
      , paste0("DurD_CV", rollingstep_yr, "yrpast") :=  
        frollapply(DurD_yr, n=rollingstep_yr, 
                   FUN=function(x) fifelse(mean(x, na.rm=T)==0,
                                           0,
                                           sd(x, na.rm=T)/mean(x, na.rm=T)
                   ), 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #CV of average annual event duration----------------------------------------
    stats_dt_yr[
      , paste0("meanConD_CV", rollingstep_yr, "yrpast") :=  
        frollapply(meanConD_yr, n=rollingstep_yr, 
                   FUN=function(x) fifelse(mean(x, na.rm=T)==0,
                                           0,
                                           sd(x, na.rm=T)/mean(x, na.rm=T)
                   ), 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #SD No-flow start date - Julian day of the first drying event---------------
    stats_dt_yr[
      , paste0("FstDrE_SD", rollingstep_yr, "yrpast") :=  
        frollapply(FstDrE, n=rollingstep_yr, 
                   FUN=function(x) sd(x, na.rm=T)
                   , 
                   align='right')
      , by=.(reach_id)
    ] 
    
    #Format data before export -------------------------------------------------
    setorderv(stats_dt_yr, c('reach_id', 'year'))
    
  }
  
  if (include_metadata) {
    stats_dt_yr <- cbind(stats_dt_yr, metadata_dt[, .(country, gcm, scenario)])
  }
  
  return(stats_dt_yr[year >= year(min_date),])
}


#------ merge_allvars_sites ----------------------------------------------------
# in_country <- 'Spain'
# in_spdiv_local <- tar_read(spdiv_local)
# in_spdiv_drn <- NULL
# in_hydrocon_compiled <- tar_read(hydrocon_sites_compiled)
# in_env_dt <- tar_read(env_dt)
# in_genal_upa = tar_read(genal_sites_upa_dt)

#' @title Merge all variables for sites (campaign-level)
#' @description Merges biodiversity, hydrological, and environmental data
#'              at the individual campaign × site × organism level.
#' @param in_spdiv_local A biodiversity data.table.
#' @param in_hydrocon_compiled Hydrological/connectivity data by date & site.
#' @param in_env_dt Environmental data at sampling times.
#' @param in_genal_upa Upstream area data for the Genal basin.
merge_allvars_sites <- function(in_spdiv_local, 
                                in_spdiv_drn=NULL,
                                in_hydrocon_compiled,
                                in_env_dt,
                                in_genal_upa) {
  
  # Prepare data  ---------------------------------------------------------------
  setDT(in_env_dt)
  
  #Fill basin area NAs in environmental data for Genal basin in Spain
  #https://stackoverflow.com/questions/72940045/replace-na-in-a-table-with-values-in-column-with-another-table-by-conditions-in
  in_env_dt[is.na(basin_area_km2),
            basin_area_km2 := in_genal_upa[.SD, on='site', x.basin_area_km2]]
  
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
  
  # Define columns -------------------------------------------------------------
  group_cols = c("running_id", "site", "date", "campaign", "organism", 
                 "organism_class", "country", "UID", "upstream_area_net")
  exclude_cols = c("ncampaigns", "name", "isflowing", "reach_length",
                   "noflow_period", "noflow_period_dur", "last_noflowdate", "drn",
                   "if_ip_number_and_size_2_axes_+_depth_of_the_pools",
                   "latitude", "longitude", "reach_id", "hy", "month",
                   "min_wetted_width", "left_river_bank_slope", "right_river_bank_slope",
                   "qsim")
  
  dtcols <- list(
    div = setdiff(names(spdiv), 
                  c(group_cols, exclude_cols,
                    names(in_hydrocon_compiled), names(in_env_dt))),
    hydro_con = setdiff(names(in_hydrocon_compiled), 
                        c(group_cols, exclude_cols,
                          names(spdiv), names(in_env_dt))),
    env = setdiff(names(in_env_dt),
                  c(group_cols, exclude_cols,
                    names(spdiv), names(in_hydrocon_compiled))),
    group_cols = group_cols,
    exclude_cols = exclude_cols
  )
  
  # numeric environmental columns
  dtcols <- c(
    dtcols,
    env_num = list(intersect(dtcols$env,
                             names(in_env_dt)[
                               sapply(in_env_dt, class) %in% c("numeric", "integer")]))
  )
  
  # Merge ----------------------------------------------------------------------
  # Merge biodiversity + hydro connectivity
  spdiv_hydro_con <- merge(spdiv, in_hydrocon_compiled,
                           by=c('date', 'site'), all.x=TRUE)
  
  # Merge environment
  all_vars_merged <- merge(
    spdiv_hydro_con, 
    in_env_dt[, c(dtcols$env, 'site', 'campaign'), with=FALSE],
    by=c('site', 'campaign'), all.x=TRUE) %>%
    .[!(site %in% dry_only_sites), ] %>%
    .[, country := factor(
      country,
      levels = c("Finland","France","Hungary","Czechia","Croatia","Spain"))] 
  
  return(list(
    dt = all_vars_merged,
    dry_only_sites = dry_only_sites,
    cols = dtcols
  ))
}

#------ merge_allvars_summarized ----------------------------------------------
# in_spdiv_local <- tar_read(spdiv_local)
# in_hydrocon_summarized <- tar_read(hydrocon_summarized)
# in_env_summarized <- tar_read(env_summarized)
# dry_only_sites <- tar_read(allvars_sites)$dry_only_sites
# dtcols_sites <- tar_read(allvars_sites)$cols

#' @title Merge all variables summarized by site
#' @description Summarizes biodiversity across campaigns and merges with
#'              site-level hydrological & environmental summaries.
#' @param spdiv A biodiversity data.table (already cleaned).
#' @param in_hydrocon_summarized Hydrological/connectivity data summarized by site.
#' @param in_env_summarized A list of summarized environmental data.tables.
#' @param in_genal_upa Upstream area data for the Genal basin.
#' @param dry_only_sites Sites that never flowed (from merge_allvars_sites).
#' @param dtcols_sites Column lists from merge_allvars_sites.
#' @return A list containing:
#'   \item{dt_summarized}{Merged site-level summarized data}
#'   \item{cols}{List of relevant column vectors for summarized merging}
merge_allvars_summarized <- function(in_spdiv_local,
                                     in_hydrocon_summarized,
                                     in_env_summarized,
                                     dry_only_sites,
                                     dtcols_sites) {
  
  #----------------- Column definitions (summarized) -----------------
  dtcols <- list(
    div_summarized = c("mean_richness", "mean_expshannon", "mean_invsimpson",
                       "JBDtotal", "JRepl", "JRichDif", "JRepl/BDtotal", "JRichDif/BDtotal",
                       "RBDtotal", "RRepl", "RRichDif", "RRepl/BDtotal", "RRichDif/BDtotal",
                       "Gamma", "Beta", "mAlpha"),
    hydro_con_summarized = setdiff(names(in_hydrocon_summarized), 
                                   c(dtcols_sites$group_cols, 
                                     dtcols_sites$exclude_cols,
                                     names(in_spdiv_local))),
    env_summarized = setdiff(c(names(in_env_summarized$dt_nopools),
                               "upstream_area_net"),
                             c(dtcols_sites$group_cols, 
                               dtcols_sites$exclude_cols))
  )
  
  dtcols <- c(
    dtcols,
    env_summarized_num = list(intersect(
      dtcols$env_summarized,
      names(in_env_summarized$dt_nopools)[
        sapply(in_env_summarized$dt_nopools, class) %in% c("numeric", "integer")]))
  )
  
  #------------ Summarize biodiversity -----------------
  spdiv_summarized <- in_spdiv_local[
    !duplicated(paste(site, organism)), 
    intersect(c(dtcols_sites$group_cols, dtcols$div_summarized),
              names(in_spdiv_local)), 
    with=FALSE] %>%
    .[, mean_richness := as.integer(round(mean_richness))]
  
  #------------ Site-level merge -----------------
  all_vars_merged <- spdiv_summarized %>%
    merge(in_hydrocon_summarized[, c(dtcols$hydro_con_summarized, 
                                     "upstream_area_net", "site"), with=FALSE],
          by="site", all.x=TRUE) %>%
    merge(in_env_summarized$dt_nopools[, c(dtcols$env_summarized, "site"), with=FALSE],
          by="site", all.x=TRUE) %>%
    .[, c("campaign", "date") := NULL] %>%
    .[!(site %in% dry_only_sites), ] %>%
    .[, country := factor(country,
                          levels = c("Finland", "France", "Hungary", 
                                     "Czechia", "Croatia", "Spain"))]
  
  return(list(
    dt = all_vars_merged,
    cols = dtcols
  ))
}

#------ create_hydro_vars_dt ---------------------------------------------------      
create_hydro_vars_dt <- function(in_hydro_vars_forssn) {
  dt <- data.table(hydro_var=in_hydro_vars_forssn)
  
  labels_dt <- data.table(
    hydro_var_root = c('DurD', 'PDurD', "PrdD",
                       'FreD', "PFreD",
                       'DurD_CV', "FreD_CV", "meanConD_CV", "FstDrE_SD", 
                       "sd6", "FstDrE",
                       "uQ90", "oQ10", "maxPQ", "PmeanQ",
                       "STcon_directed", "STcon_undirected",
                       "Fdist_mean_directed", "Fdist_mean_undirected"),
    hydro_label = c(
      "Proportion of no-flow days",
      "Proportion of no-flow days (percentile)",
      "Time to last no flow date",
      "Number of no-flow periods",
      "Number of no-flow periods (percentile)",
      "CV of the annual proportion of no-flow days",
      "CV of the annual number of no-flow periods",
      "CV of the mean annual drying event duration",
      "SD of the date of first drying",
      "Seasonality of drying (SD6)",
      "Date of first drying",
      "Number of low-flow days (< Q90)",
      "Number of high-flow days (>Q10)",
      "Maximum flow percentile",
      "Mean flow percentile",
      "Spatio-temporal connectivity - upstream",
      "Spatio-temporal connectivity - undirected",
      "Mean distance to the nearest flowing site - upstream",
      "Mean distance to the nearest flowing site - undirected"),
    hydro_class = c(
      rep('Drying duration', 3),
      rep('Drying frequency', 2),
      rep('Drying predictability', 4),
      rep('Drying timing', 2),
      rep('Flow magnitude', 4),
      rep('Connectivity', 4)
    )
  ) %>%
    .[, metric_num:=seq_along(hydro_label), by=hydro_class]
  
  dt <- get_hydro_var_root(dt, in_place=F) %>%
    merge(., labels_dt, by='hydro_var_root') %>%
    .[, `:=`(hydro_label = factor(hydro_label,
                                  levels=labels_dt$hydro_label,
                                  ordered=T),
             hydro_class = factor(hydro_class,
                                  levels=unique(labels_dt$hydro_class),
                                  ordered=T)
    )] %>%
    rbind(data.table(hydro_var='null', hydro_var_root='null',
                     hydro_label='Null'), fill=T) 
  
  return(dt)
}

#------ plot_edna_biof_vs_sedi -------------------------------------------------
# in_allvars_sites <- tar_read(allvars_sites)

#' @title Plot eDNA biofilm vs sediment
#' @description This function creates two plots comparing richness and Gamma 
#'     diversity between eDNA samples from biofilm and sediment, faceted by organism 
#'     class and country.
#' @param in_allvars_sites A list containing the merged data tables, specifically `in_allvars_sites$dt`.
#' @return A list containing two ggplot objects: `richness` and `site_gamma`.
plot_edna_biof_vs_sedi <- function(in_allvars_sites) {
  allvars_edna <- setDT(in_allvars_sites$dt)[organism_class %in% c('dia', 'fun', 'bac'),] 
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
# in_allvars_sites <- tar_read(allvars_sites)

#' @title Compute grouped correlation matrices
#' @description This function computes correlation matrices for different sets of 
#'     variables (hydrology, environment, diversity) and different groupings 
#'     (overall, by organism, by country).
#' @param in_allvars_sites A list containing the data table (`dt`) and a list 
#'     of column names (`cols`) categorized by their origin.
#' @return A list containing five correlation matrices for different variable 
#'     combinations and groupings, plus the original column list.
compute_cor_matrix <- function(in_allvars_sites) {
  dt <- in_allvars_sites$dt
  cols_by_origin <- in_allvars_sites$cols
  
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
# in_allvars_summarized <- tar_read(allvars_summarized)

#' @title Compute summarized correlation matrices
#' @description This function computes correlation matrices using site-summarized 
#'      data (across all sampling dates), analyzing relationships between 
#'      hydrological, environmental, and diversity metrics at a broader scale.
#' @param in_allvars_summarized A list containing the summarized data table (`dt_summarized`) and column list (`cols`).
#' @return A list of correlation matrices for summarized data.
compute_cor_matrix_summarized <- function(in_allvars_summarized) {
  dt <- in_allvars_summarized$dt
  cols_by_origin <- in_allvars_summarized$cols
  
  # --- Calculate Correlations for each site summarized ---
  # 1. Overall Correlation (Hydro)
  cor_hydro <- compute_cor_matrix_inner(
    in_dt = dt,
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

# in_allvars_dt <- tar_read(allvars_sites)$dt
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

#' @title Create Spatial Stream Network prediction points
#' @description This function prepares spatial data for SSN modeling by creating 
#'     prediction points from a river network and associating them with historical 
#'     and projected hydrological data.
#' @param in_network_path A list of file paths to the river network shapefiles, separated by country.
#' @return A list containing two spatial objects: one for historical prediction 
#'      points and one for projected prediction points.
create_ssn_pred_pts <- function(in_network_path) {
  # Process the river network data
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
  
  return(net_centroids)
}

#------ create_ssn_europe ------------------------------------------------------
# in_network_path = tar_read(network_ssnready_shp_list)
# in_sites_path = tar_read(site_snapped_gpkg_list)
# in_barriers_path = tar_read(barrier_snapped_gpkg_list)
# in_allvars_dt= tar_read(allvars_sites)$dt
# in_local_env_pca = tar_read(local_env_pca)
# in_hydrostats_net_hist = tar_read(hydrostats_net_hist)
# in_pred_pts = NULL
# out_dir = file.path(resdir, 'ssn')
# out_ssn_name = 'ssn_eu'
# overwrite=T

# in_network_path = tar_read(network_ssnready_shp_list)
# in_sites_path = tar_read(site_snapped_gpkg_list)
# in_allvars_dt = tar_read(allvars_summarized)$dt
# in_local_env_pca = tar_read(local_env_pca_summarized)
# in_hydrostats_net_hist = tar_read(hydrostats_net_hist)
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
                              in_barriers_path = NULL,
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
    
    net_hydro <- merge(net_proj, 
                       in_hydrostats_net_hist[country==in_country,],
                       by.x = 'cat', by.y = 'reach_id',
                       all.x=T)
    
    return(net_hydro)
  }) %>% do.call(rbind, .)
  
  
  # net_proj[!net_proj$cat %in% in_hydrostats_net_hist[country=='France',]$reach_id,]
  
  # st_write(net_eu, file.path(resdir, 'test.gpkg'))
  
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
  
  setDT(in_allvars_dt)[country == 'Czech Republic', country := 'Czechia']
  
  # merge site data with all variables and PCA results
  sites_lsn_attri <- merge(sites_lsn,
                           in_allvars_dt, 
                           by=c('country', 'site', 'upstream_area_net')) %>%
    merge(in_local_env_pca$dt_all,
          by=c(id_cols, 'organism_class'))
  
  sites_list <- list(sites = sites_lsn_attri)
  
  #3. Incorporate placeholders for prediction "sites"  -------------------------
  #   but without the actual hydrological data (to keep the SSN light)-----
  if (!(is.null(in_pred_pts))) {
    #Add historical prediction sites
    if (!(crs(in_pred_pts) == crs(sites_eu))) {
      in_pred_pts<- st_transform(in_pred_pts, 3035) 
    }
    
    preds_hist_lsn <- SSNbler::sites_to_lsn(
      sites = in_pred_pts,
      edges =  edges_lsn,
      lsn_path = lsn_path,
      file_name = "preds_hist",
      snap_tolerance = 5,
      save_local = TRUE,
      overwrite = overwrite
    )
    
    sites_list$preds_hist <- preds_hist_lsn
    
    #Add future prediction "sites"
    preds_proj_lsn <- SSNbler::sites_to_lsn(
      sites = in_pred_pts,
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
  if (!is.null(in_barriers_path)) {
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
                         ptcolor_col=NULL, ptcolor_lims=NULL,
                         nbreaks=NULL,
                         cover_zero=FALSE 
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
          n.breaks = ifelse(is.null(nbreaks), 10, nbreaks),
          palette = "Spectral"
        )
    } else {
      out_map <- out_map +
        scale_color_viridis_b(
          name = str_to_sentence(gsub('_', ' ', linecolor_col)),
          limits = linecolor_lims,
          n.breaks=ifelse(is.null(nbreaks), 5, nbreaks),) 
    }
    
    if (cover_zero) {
      out_map <- out_map +
        geom_sf(data = in_edges[in_edges[[linecolor_col]]==0,],
                aes(
                  linewidth = !!sym(linewidth_col)
                ),
                color = 'grey'
        )
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
    
    
    if (length(setdiff(sign(ptcolor_lims), 0)) == 2) {
      # Data spans both positive and negative values
      # Center the limits on 0 and make them symmetric
      max_abs <- max(abs(ptcolor_lims))
      out_map <- out_map +
        scale_color_fermenter(
          name = str_to_sentence(gsub('_', ' ', ptcolor_col)),
          limits = c(-max_abs, max_abs),
          n.breaks = ifelse(is.null(nbreaks), 10, nbreaks),
          palette = "Spectral",
          direction = -1  # Reverse the palette (blue for positive, red for negative)
        )
    } else {
      out_map <- out_map + 
        scale_color_viridis_b(
          name = str_to_sentence(gsub('_', ' ', ptcolor_col)),
          limits = ptcolor_lims,
          n.breaks=ifelse(is.null(nbreaks), 5, nbreaks),) 
    }
    
    
    if (!is.null(linecolor_col)) {
      out_map <- out_map +  new_scale_color()
    }

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
                           nbreaks = NULL,
                           page_title=NULL,
                           cover_zero = FALSE) {
  
  assert_that((facet_col %in% names(in_ssn$edges)) &
                (facet_col %in% names(in_ssn$obs)),
              msg = paste0(facet_col, ' not in in_ssn edges or obs'))
  
  #Define the different facet values to iterate over
  #Order countries by level of intermittence (RelF)
  if (facet_col == 'country') {
    in_ssn$edges[[facet_col]] <- factor(
      in_ssn$edges[[facet_col]],
      levels = c("Finland", "France",  "Hungary", "Czechia", "Croatia", "Spain" ),
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
    facet_map <- map_ssn_util(
      in_ssn = in_ssn, 
      in_edges = edges_i,
      linewidth_col = linewidth_col,                                
      linecolor_col = linecolor_col,
      linecolor_lims = linecolor_lims,
      in_pts = pts_i, 
      shape_col = shape_col,
      ptcolor_col = ptcolor_col, 
      ptcolor_lims = ptcolor_lims,
      nbreaks = nbreaks,
      cover_zero = cover_zero
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
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_organism_dt = tar_read(organism_dt)
# verbose = T

#' @title Map summarized SSN variables
#' @description Generates a series of faceted maps visualizing diversity and 
#'      physical variables across the European SSN. It creates separate plots 
#'      for each variable, faceted by country.
#' @param in_ssn_summarized A list of SSN objects, typically for different organisms.
#' @param in_allvars_summarized A list containing a data table of all variables and column names.
#' @param in_organism_dt A data table mapping organism codes to their labels.
#' @param verbose A logical value to indicate whether to print progress messages.
#' @return A named list of ggplot objects, combining maps of diversity and physical variables.
map_ssn_summarized <- function(in_ssn_summarized,
                               in_allvars_summarized,
                               in_organism_dt,
                               verbose = T) {
  
  #Map every diversity variable for each organism
  # create a data frame of all organism-diversity variable combinations
  div_map_params <- expand.grid(intersect(names(in_ssn_summarized), 
                                          in_organism_dt$organism),
                                in_allvars_summarized$cols$div_summarized,
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
  physvars <- c(in_allvars_summarized$cols$hydro_con_summarized,
                in_allvars_summarized$cols$env_summarized_num,
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
# in_allvars_dt <- tar_read(allvars_sites)$dt
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
    levels = c("Finland", "France", "Hungary", "Czechia", "Croatia", "Spain"),
    ordered = TRUE
  )]
  
  hydrocon_compiled_country[, relD10past := 1 - relF10past]
  
  relF_melt <- hydrocon_compiled_country[!duplicated(paste(country, date)), ] %>%
    melt(id.vars = c('country', 'date'), measure.vars = c('relF10past', 'relD10past'))
  
  # 1. Plot mean proportion of flowing river network length over time
  plot_relF10past <- ggplot(relF_melt) +
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
      name = str_wrap('Mean % of network flowing/non-flowing over previous 10 days', 40),
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
  plot_relF10past_sampling <- plot_relF10past +
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
    
    plot_alpha <- plot_relF10past + 
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
    ggsave(file.path(out_dir, 'relF10past_drn.png'), 
           plot = plot_relF10past, 
           width = 10, height = 6, 
           units = 'in', dpi = 600)
    ggsave(file.path(out_dir, 'relF10past_drn_sampling.png'), 
           plot = plot_relF10past_sampling, 
           width = 10, height = 6, 
           units = 'in', dpi = 600)
    
    lapply(names(alpha_plot_list), function(in_organism_label) {
      ggsave(file.path(out_dir, 
                       paste0('relF10past_drn_', alpha_var, '_',
                              in_organism_label, '.png')),
             plot = alpha_plot_list[[in_organism_label]],
             width = 10, height = 6, units = 'in', dpi = 600)
    })
  }
  
  return(list(
    hydro = plot_relF10past,
    hydro_sampling = plot_relF10past_sampling,
    alpha_list = alpha_plot_list
  ))
}
#------ check_diversity_dependence_on_habitat_volume ---------------------------
#n_allvars_sites <- tar_read(allvars_sites)

#' @title Check correlation between alpha diversity and habitat volume
#' @description Creates diagnostic scatter plots to examine the relationship 
#'     between alpha diversity and local habitat volume metrics (e.g., depth, width, velocity),
#'     faceted by country.
#' @param in_allvars_sites A list containing a data.table of all merged variables.
#' @param alpha_var Character string: name of the alpha diversity column to use (e.g., "richness", "shannon").
#' @return A list of ggplot objects showing the relationships for different organism groups.
check_cor_div_habvol <- function(in_allvars_sites, alpha_var = "richness") {
  
  dt <- copy(in_allvars_sites$dt)
  
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
# in_dt=tar_read(allvars_summarized)$dt
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
      levels = c("Finland", "France",  "Hungary", "Czechia", "Croatia", "Spain"),
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
# in_allvars_sites <- tar_read(allvars_sites)
# response_var = 'richness' #c('richness', 'invsimpson','Jtm1'),
# colors_list = drn_dt$color
# plot=T
# out_dir = figdir
# in_organism_dt <- tar_read(organism_dt)
# plot_name_suffix = ''
# 
# in_allvars_sites=tar_read(allvars_sites)
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
#' @param in_allvars_sites A list containing a data table of all merged variables.
#' @param in_organism_dt A data table mapping organism codes to their labels.
#' @param temporal_var_substr A string to filter for specific hydrological variables (e.g., 'DurD').
#' @param response_var A list containing strings specifying the biological response variables (e.g., 'richness').
#' @param colors_list A named character vector of colors for each country.
#' @param write_plots A logical value to save the plot as a PNG file.
#' @param plot_name_suffix A string to append to the output filename.
#' @param out_dir The output directory for the saved plot.
#' @return A ggplot object of the scatter plot with linear model fits.
plot_scatter_lm <-  function(in_allvars_sites, 
                             in_organism_dt,
                             temporal_var_substr, 
                             response_var,
                             colors_list, 
                             write_plots=T,
                             plot_name_suffix="", 
                             out_dir) {
  
  dt <- in_allvars_sites$dt
  x_cols <- grep(paste0('^', temporal_var_substr), names(dt), value=T) %>%
    .[seq(1, length(.), 2)]
  
  dt_melt <- melt(dt, 
                  id.vars=intersect(
                    c(in_allvars_sites$cols$group_cols, response_var), 
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
# organism = 'miv_nopools_ept'
# formula_root = ' log10(basin_area_km2) + log10(basin_area_km2):country'
# hydro_var = 'DurD365past'
# response_var = 'invsimpson'
# tar_load(ssn_covtypes)
# family = "Gaussian"
# estmethod = "ml"
# partition_formula = as.formula("~ as.factor(campaign)")
# random_formula = as.formula("~ country")
# standardize_hydro_var = T
# include_state_of_flow = F
# include_seasonality = F

# tar_load(ssn_covtypes)
# in_ssn = tar_read(ssn_eu_summarized)
# organism = c('bac_biof_nopools')
# formula_root = ' log10(basin_area_km2)'
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
                                  standardize_hydro_var = T,
                                  include_state_of_flow = F,
                                  include_seasonality = F) {
  
  # hydro_var <- grep(paste0('^', hydro_var_str), 
  #                   names(in_ssn[[organism]]$ssn$obs), 
  #                   value=T)
  
  # Standardize hydrological variable
  if ((standardize_hydro_var) & !is.null(hydro_var)) {
    in_ssn[[organism]]$ssn$obs[[paste0(hydro_var, '_scaled')]] <-
      as.numeric(scale(in_ssn[[organism]]$ssn$obs[[hydro_var]]))
    hydro_var <- paste0(hydro_var, '_scaled')
  }
  
  # Construct the full model formula
  if (!is.null(hydro_var)) {
    full_formula <- paste0(response_var, ' ~ ', hydro_var, ' + ',
                           hydro_var, ':country +', formula_root)
  } else {
    full_formula <- paste0(response_var, ' ~ ', formula_root)
  }
  
  if (include_seasonality) {
    full_formula <- paste(full_formula, '+ doy + I(doy^2)')
  }
  
  if (include_state_of_flow) {
    full_formula <- paste(full_formula, '+ state_of_flow')
  }
  
  # Fit SSN models using quick_ssn
  ssn_list <- quick_ssn(in_ssn = in_ssn[[organism]]$ssn, 
                        in_formula = as.formula(full_formula),
                        family = family,
                        ssn_covtypes = ssn_covtypes, #ssn_covtypes[sample(144, 10)],
                        partition_formula = partition_formula,
                        random_formula = random_formula,
                        estmethod = estmethod)
  
  # Collect model summary statistics (glance)
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

#------ prepare_hydrowindow_perf_table -----------------------------------------------------
# in_response_var <- 'invsimpson'
# in_organism <- 'miv_nopools_ept'
# tar_load(ssn_div_models_to_run)
# ssn_div_hydrowindow <- tar_read_raw(paste0('ssn_div_hydrowindow_', 
#                                            in_response_var, '_', in_organism))
# 
# 
# ssn_model_names <- do.call(rbind, ssn_div_models_to_run)[
#   , c("organism", "hydro_var", "response_var")] %>%
#   as.data.table() %>%
#   .[response_var == in_response_var & organism == in_organism, ]
# 
# ssnmodels <- cbind(ssn_model_names, ssn_div_hydrowindow)
# names(ssnmodels)[ncol(ssnmodels)] <- "ssn_div_models"
# 
# in_ssnmodels <- ssnmodels 
# in_covtype_selected <- tar_read_raw(paste0('ssn_covtype_selected_', 
#                                            in_response_var, '_', in_organism))
# in_hydro_vars_dt <- tar_read(hydro_vars_dt)

#' Prepare performance table summarizing hydrological window SSN models
#'
#' Extracts `glance` tables and fitted model objects, keeps only the selected 
#' covariance type, and attaches hydro variable roots and labels.
#'
#' @param in_ssnmodels A `data.table` with columns `organism`, `hydro_var`, `response_var`,
#'   and `ssn_div_models` (list-column of fitted models).
#' @param in_organism Character. Organism identifier (e.g. `"miv_nopools"`).
#' @param in_covtype_selected A list returned by `select_ssn_covariance()` containing
#'   a `dt_sub` element with best covariance type per organism.
#' @param in_hydro_vars_dt A `data.table` with hydro variable metadata including
#'   columns `hydro_var_root` and `hydro_label`.
#'
#' @return A `data.table` with model performance information (AICc, variance components, etc.).
#' @export
prepare_hydrowindow_perf_table <- function(in_ssnmodels, in_organism, 
                                           in_covtype_selected, in_hydro_vars_dt) {
  # selected covariance type for this organism
  org_covtype <- in_covtype_selected$dt_sub[organism == in_organism, covtypes]
  
  # combine glance tables with models
  perf_dt <- lapply(
    in_ssnmodels[organism == in_organism, ssn_div_models],
    function(x) {
      merge(
        x[["ssn_glance"]],
        data.table(covtypes = names(x[["ssn_list"]]), mod = x[["ssn_list"]]),
        by = "covtypes"
      )
    }) %>%
    rbindlist(fill = TRUE) %>%
    .[covtypes == org_covtype, ] %>%
    .[is.na(hydro_var), hydro_var := "null"]
  
  # attach hydro_var_root
  get_hydro_var_root(perf_dt)
  
  # attach hydro labels
  perf_dt <- merge(
    perf_dt,
    unique(in_hydro_vars_dt, by = "hydro_var_root")[, .(hydro_var_root, hydro_label)],
    by = "hydro_var_root"
  )
  
  #Select best models
  best_dt <- perf_dt[, .SD[which.min(AICc)], by = hydro_var_root]
  
  
  return(list(
    all=perf_dt,
    best=best_dt
  ))
}

#------ get_hydrowindow_varcomp -------------------------------------------------------
# hydrowindown_perf_tables <- tar_read(hydrowindow_perf_tables_richness_fun_sedi_nopools)
# perf_dt <- hydrowindown_perf_tables$all
# nrow_pag = 2
# ncol_pag = 3

#' Extract and plot SSN variance decomposition 
#'
#' Computes variance components from selected SSN models, and Creates stacked 
#' barplots of variance components across hydrological windows.
#'
#' @param perf_dt A `data.table` returned by [prepare_hydrowindow_perf_table()].
#'
#' @return A list containing dt: `data.table` of variance components with labels
#'  and colors for plotting, and plots: a list of stacked barcharts
#' @export
get_hydrowindow_varcomp <- function(perf_dt, nrow_pag = 2, ncol_pag = 3) {
  varcomp_labels <- c(
    "Remaining variance (nugget)",
    "Spatially dependent variance - euclidean",
    "Spatially dependent variance - taildown",
    "Spatially dependent variance - tailup",
    "Random intercept - country",
    "Fixed effects pseudo-R2"
  )
  
  vc_dt <- perf_dt[
    , {
      if (inherits(mod[[1]], c("ssn_lm", "ssn_glm"))) { #To deal with models that did not work
        as.data.table(SSN2::varcomp(mod[[1]])) #Extract variance decomposition
      } else data.table()
    },
    by = .(response_var, hydro_var, covtypes, hydro_label, window_d)
  ] %>%
    merge(
      data.table(
        varcomp = c("nugget","euclid_de","taildown_de",
                    "tailup_de","1 | country","Covariates (PR-sq)"),
        varcomp_label = factor(varcomp_labels,
                               levels = varcomp_labels,
                               ordered = TRUE),
        varcomp_color = c("lightgrey","#b2df8a","#a6cee3",
                          "#1f78b4","#f1b6da","#feb24c")
      ),
      by = "varcomp"
    )
  
  #Make plots
  varcomp_subdt <- vc_dt[proportion > 0 & hydro_label != "Null", ]
  varcomp_null <- vc_dt[proportion > 0 & hydro_label == "Null", ][
    , hydro_label := "Null: only basin area"]
  
  #Determine the number of plot pages required
  page_num <- ceiling(varcomp_subdt[, length(unique(hydro_label))] 
                      / (nrow_pag * ncol_pag))
  
  #Plot multiple pages of variance decomposition plots, once facet per variable
  vc_plot_list <- lapply(seq_len(page_num), function(in_page) {
    p_main <- ggplot(varcomp_subdt, 
                     aes(x = window_d, y = proportion, fill = varcomp_label)) +
      geom_bar(stat = "identity", 
               position = "stack",
               alpha = 0.75) +
      geom_text(aes(label = round(100 * proportion)),
                position = position_stack(vjust = 0.5),
                size = 3, colour = "#555555") +
      scale_fill_manual(name = "Variance components",
                        values = vc_dt[proportion > 0,][
                          order(varcomp_label), unique(varcomp_color)]) +
      scale_x_discrete("Temporal window of aggregation") +
      scale_y_continuous("Percentage of variance", 
                         breaks = c(0,0.5,1), 
                         labels = scales::label_percent()) +
      facet_wrap_paginate(~hydro_label, 
                          nrow = nrow_pag, ncol = ncol_pag,
                          labeller = label_wrap_gen(width = 25), 
                          scales = "free_x",
                          page = in_page) +
      coord_cartesian(expand = FALSE) + 
      theme_bw() +
      theme(legend.background = element_blank())
    
    p_null <- ggplot(varcomp_null, 
                     aes(x = window_d, y = proportion, fill = varcomp_label)) +
      geom_bar(stat = "identity", 
               position = "stack",
               alpha = 0.75)  +
      geom_text(aes(label = round(100 * proportion)),
                position = position_stack(vjust = 0.5),
                size = 3, colour = "#555555") +
      scale_fill_manual(name = "Variance components",
                        values = vc_dt[proportion > 0,][
                          order(varcomp_label), unique(varcomp_color)]) +
      facet_wrap(~hydro_label, 
                 labeller = label_wrap_gen(width = 15)) +
      theme_bw() + 
      theme(legend.position = "none", 
            axis.title = element_blank(),
            axis.text = element_blank())
    
    #Merge plot elements
    (p_main + guide_area()) + p_null +
      plot_layout(design = "AAAAAAAABBB\nAAAAAAAAC##", axes = "collect", guides = "collect")
  })
  
  return(list(
    dt=vc_dt,
    plots=vc_plot_list)
  )
}

#------ get_hydrowindow_emmeans --------------------------------------------------------
# hydrowindown_perf_tables <- tar_read(hydrowindow_perf_tables_richness_fun_sedi_nopools)
# best_dt <- hydrowindown_perf_tables$best
# in_hydro_vars_dt <- tar_read(hydro_vars_dt)
# in_drn_dt <- drn_dt

#' Extract and plot marginal means (EMMeans) from best models
#'
#' Runs [get_ssn_emmeans()] for all predictors in the best models.
#'
#' @param best_dt A `data.table` of best models ($best) returned by [prepare_hydrowindow_perf_table()].
#' @param in_hydro_vars_dt Metadata table with hydro variable names and labels.
#'
#' @return A `data.table` and plots of estimated marginal means across predictors.
#' @export
get_hydrowindow_emmeans <- function(best_dt, in_hydro_vars_dt, in_drn_dt) {
  
  emmeans_dt_all <- lapply(best_dt$hydro_var, function(in_pred_var) {
    in_mod <- best_dt[hydro_var == in_pred_var, mod][[1]]
    if (in_pred_var == "null") in_pred_var <- all.vars(in_mod$formula)[2]
    
    pred_var_label <- ifelse(
      in_pred_var == "null",
      in_pred_var,
      get_full_hydrolabel(in_hydro_vars_dt, in_pred_var)
    )
    
    emmeans_row <- get_ssn_emmeans(in_mod=in_mod,
                                   in_pred_var=in_pred_var, 
                                   in_pred_var_label=pred_var_label,
                                   interaction_var = "country", 
                                   in_drn_dt = in_drn_dt,
                                   plot = FALSE)$dt
    
    return(emmeans_row)
  }) %>% rbindlist(use.names = TRUE, fill = TRUE)  %>%
    setnames('pred_var_name', 'hydro_var') %>%
    get_hydro_var_root(in_place=F) %>%
    merge(in_hydro_vars_dt[!duplicated(hydro_var_root), 
                           .(hydro_var_root, hydro_class)], 
          by='hydro_var_root')
  
  
  emmeans_plot <- get_ssn_emmeans(in_emm_dt = emmeans_dt_all, 
                                  interaction_var='country', plot=T)$plot +
    facet_wrap(~ hydro_class + str_wrap(pred_var_label, 30),
               ncol = 4,
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 list(do.call(paste, c(labels, list(sep = "\n"))))
               }, 
               scales='free')
  
  return(list(
    dt=emmeans_dt_all,
    plot=emmeans_plot
  ))
}

#------ get_hydrowindow_emtrends -------------------------------------------------------
# hydrowindown_perf_tables <- tar_read(hydrowindow_perf_tables_richness_bac_biof_nopools)
# best_dt <- hydrowindown_perf_tables$best
# in_hydro_vars_dt <- tar_read(hydro_vars_dt)
# 
#' Extract and plot marginal slopes (EMTrends) from best models
#'
#' Runs [get_ssn_emtrends()] for all predictors in the best models.
#'
#' @param best_dt A `data.table` ($best) returned by [prepare_hydrowindow_perf_table()].
#' @param in_hydro_vars_dt Metadata table with hydro variable names and labels.
#'
#' @return A `data.table` and plots of estimated slopes across predictors.
#' @export
get_hydrowindow_emtrends <- function(best_dt, in_hydro_vars_dt, in_drn_dt) {
  
  emtrends_dt_all <- lapply(best_dt$hydro_var, function(in_pred_var) {
    print(in_pred_var)
    in_mod <- best_dt[hydro_var == in_pred_var, mod][[1]]
    if (in_pred_var == "null") in_pred_var <- all.vars(in_mod$formula)[2]
    
    get_ssn_emtrends(in_mod, in_pred_var, in_pred_var, in_drn_dt,
                     interaction_var = "country", plot = FALSE)$dt
    
  }) %>% rbindlist(use.names = TRUE, fill = TRUE) %>%
    .[, pred_var_label := get_full_hydrolabel(in_hydro_vars_dt, pred_var_name),
      by=.(pred_var_name, country)] %>%
    setnames('pred_var_name', 'hydro_var') %>%
    get_hydro_var_root(in_place=F) %>%
    merge(in_hydro_vars_dt[!duplicated(hydro_var_root), 
                           .(hydro_var_root, hydro_class)], 
          by='hydro_var_root') %>%
    .[!(pred_var_label %in% c('Null', 'Mean flow percentile')),]
  
  
  emtrends_dt_all[, `:=`(
    lcl_capped = fifelse(asymp.LCL < min(trend), min(trend), asymp.LCL),
    ucl_capped = fifelse(asymp.UCL > max(trend), max(trend), asymp.UCL)
  ), by=response_var]
  
  emtrends_plot <- get_ssn_emtrends(in_emtrends_dt = emtrends_dt_all, 
                                    interaction_var='country', plot=T)$plot +
    facet_wrap(~ hydro_class + str_wrap(pred_var_label, 30),
               ncol = 4,
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 list(do.call(paste, c(labels, list(sep = "\n"))))
               })
  
  return(list(
    dt=emtrends_dt_all,
    plot=emtrends_plot
  ))
}

#------ get_hydrowindow_predictions --------------------------------------------
# in_organism = 'miv_nopools_ept'
# in_response_var = 'invsimpson'
# tar_ext <- paste0('_', in_response_var, '_', in_organism)
# best_dt <- tar_read_raw(paste0('hydrowindow_perf_tables', tar_ext))$best

#' Get predictions from best models
#'
#' Uses `SSN2::augment()` to obtain fitted vs observed values and attaches metadata.
#'
#' @param best_dt A `data.table` ($best) returned by [prepare_hydrowindow_perf_table()].
#'
#' @return A `data.table` of observed and fitted values with hydro variable metadata.
#' @export
get_hydrowindow_predictions <- function(best_dt) {
  lapply(seq(nrow(best_dt)), function(i) {
    print(i)
    SSN2::augment(best_dt[i, mod[[1]]], drop = FALSE) %>%
      cbind(best_dt[i, .(hydro_var, response_var, covtypes, hydro_var_root, window_d)])
  }) %>% rbindlist(fill=T)
}

#------ plot_hydrocon_summarized -----------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_drn_dt <- drn_dt
# in_hydro_vars_dt <- tar_read(hydro_vars_dt)
# outdir = figdir

plot_hydrocon_summarized <- function(in_allvars_summarized,
                                     in_drn_dt, in_hydro_vars_dt, 
                                     outdir, write_plot=T) {
  
  metacols_sub <- c(intersect(metacols, names(in_allvars_summarized$dt)), 'stream_type')
  hydrocon_sub <- in_allvars_summarized$dt %>%
    .[!duplicated(site) & stream_type=='TR',
      c(metacols_sub,
        in_allvars_summarized$cols$hydro_con_summarized), with=F] %>%
    .[, Fdist_mean_10past_undirected_avg_samp_log10 := 
        log10(Fdist_mean_10past_undirected_avg_samp+0.1)] 
  
  hydrocon_melt <- hydrocon_sub %>%
    melt(id.vars=metacols_sub) %>%
    merge(in_drn_dt, by='country')  %>%
    .[, country := factor(
      country,
      levels = c("Finland", "France",  "Hungary", "Czechia", "Croatia", "Spain" ),
      ordered=T)
    ] %>%
    .[, variable_name := get_full_hydrolabel(in_hydro_vars_dt, 
                                             in_hydro_var=variable),
      by=.I]
  
  color_vec <-  hydrocon_melt[!duplicated(country),
                          setNames(color, country)]
  
  vars_sub <- c("DurD_samp", "PDurD365past", "FreD_samp", "meanConD_yr", 
    "FstDrE", "DurD_CV30yrpast", "FstDrE_SD30yrpast", 
    "STcon_m10_undirected_avg_samp", 
    "Fdist_mean_10past_undirected_avg_samp_log10")
  
  hydrocon_summarized_plot <- hydrocon_melt[variable %in% vars_sub,] %>%
    ggplot(aes(x=country, y=value, fill=country)) +
    geom_boxplot(alpha=0.5) +
    scale_fill_manual(name='Country', values=color_vec) +
    facet_wrap(~variable_name,
               scales='free_x',
               labeller = label_wrap_gen(width = 22)) +
    coord_flip() +
    theme_classic() +
    theme(legend.position='none')
  
  out_plot_path <- file.path(outdir, 'hydrocon_summarized_boxplot.png')
  if (!file.exists(out_plot_path)) {
    ggsave(out_plot_path,
           hydrocon_summarized_plot,
           width = 6, height = 6, dpi=300)
  }
  
  out_tab_path <- file.path(outdir, 'hydrocon_summarized.csv')
  if (!file.exists(out_tab_path)) {
    vtable::st(
      hydrocon_sub, 
      group = 'country', 
      group.long = FALSE,
      out = "csv",
      file = out_tab_path)
  }
  
  return(list(
    plot = hydrocon_summarized_plot,
    tab_path = out_tab_path
    )
  )
}


#------ plot_hydrowindow_obs_preds ----------------------------------------------
#' Plot observed vs predicted values
#'
#' Scatterplot of observed vs fitted values for each hydro variable and country.
#'
#' @param preds A `data.table` returned by [get_predictions()].
#' @param resp_var Character. Name of the response variable.
#'
#' @return A `ggplot2` object.
#' @export
plot_hydrowindow_obs_preds <- function(preds, resp_var) {
  ggplot(preds, aes(x = .fitted, y = get(resp_var), color = country)) +
    geom_abline() + geom_point() + geom_smooth(method = "lm") +
    scale_y_continuous(name = resp_var) +
    facet_grid(country ~ hydro_var) + theme_bw()
}

#------ plot_hydrowindow_x_preds -----------------------------------------------
#' Plot predictors vs fitted values
#'
#' Scatterplots of predictor variables against fitted values across countries.
#'
#' @param preds A `data.table` returned by [get_predictions()].
#' @param best_dt A `data.table` ($best) returned by [prepare_hydrowindow_perf_table()].
#'
#' @return A `ggplot2` object.
#' @export
plot_hydrowindow_x_preds <- function(preds, best_dt) {
  preds_melt <- preds[, c(setdiff(best_dt$hydro_var,"null"), "country","basin_area_km2",".fitted"), with = FALSE] %>%
    melt(id.vars = c("country","basin_area_km2",".fitted"))
  
  ggplot(preds_melt, aes(x = value, y = .fitted, color = country)) +
    geom_point() + geom_smooth(method = "lm") +
    facet_wrap(country ~ variable, scales = "free") + theme_bw()
}


#------ save_ssn_div_hydrowindow_plots ------------------------------------
# in_organism = 'miv_nopools_ept'
# in_response_var = 'invsimpson'
# tar_ext <- paste0('_', in_response_var, '_', in_organism)
# 
# perf_table <- targets::tar_read_raw(paste0('hydrowindow_perf_tables', tar_ext))
# plot_varcomp = targets::tar_read_raw(paste0('hydrowindow_varcomp', tar_ext))
# plot_obs_preds = targets::tar_read_raw(paste0('hydrowindow_obs_preds_plot', tar_ext))
# plot_x_preds = targets::tar_read_raw(paste0('hydrowindow_x_preds_plot', tar_ext))
# plot_emmeans = targets::tar_read_raw(paste0('hydrowindow_emmeans', tar_ext))
# plot_emtrends = targets::tar_read_raw(paste0('hydrowindow_emtrends', tar_ext))
# out_dir <- figdir

#' Save SSN hydro-window results (plots and tables)
#'
#' Saves performance tables and plots from the modular hydro-window pipeline to a specified directory.  
#'
#' - The performance table is saved as CSV.  
#' - Variance decomposition plots (list of ggplots, one per page) are saved as multiple PNGs.  
#' - Other plots (`obs_preds`, `x_preds`, `emmeans`, `emtrends`) are saved as individual PNGs if present.  
#'
#' @param perf_table A `data.table` of performance results (from [prepare_hydrowindow_perf_table()]), 
#'   with the `mod` column removed.
#' @param plot_varcomp A list of `ggplot2` objects returned by [get_hydrowindow_varcomp()].  
#' @param plot_obs_preds A `ggplot2` object returned by [plot_hydrowindow_obs_preds()].  
#' @param plot_x_preds A `ggplot2` object returned by [plot_hydrowindow_x_preds()].  
#' @param plot_emmeans A `ggplot2` object returned by [get_all_emmeans()] (with `plot=TRUE`).  
#' @param plot_emtrends A `ggplot2` object returned by [get_all_emtrends()] (with `plot=TRUE`).  
#' @param in_organism Character. Organism identifier.  
#' @param in_response_var Character. Response variable identifier.  
#' @param out_dir Character. Directory to save results.  
#'
#' @return A named list of file paths to saved CSV and PNG files.  
#' @export
save_ssn_div_hydrowindow_plots <- function(
    perf_table,
    plot_varcomp = NULL,
    plot_obs_preds = NULL,
    plot_x_preds = NULL,
    plot_emmeans = NULL,
    plot_emtrends = NULL,
    in_organism,
    in_response_var,
    out_dir) {
  
  # Save performance table
  table_path <- file.path(
    out_dir,
    paste0(
      "ssn_div_hydrowindow_perftable_", 
      in_organism, "_", in_response_var, "_",
      format(Sys.Date(), "%Y%m%d"), ".csv"
    )
  )
  fwrite(perf_table$all[, .SD, .SDcols=!c('mod')], table_path)
  
  out_paths <- list(table = table_path)
  
  # Variance decomposition: multi-page
  if (!is.null(plot_varcomp)) {
    vc_paths <- vapply(seq_along(plot_varcomp$plots), function(i) {
      out_path <- file.path(
        out_dir,
        paste0(
          "ssn_div_hydrowindow_", 
          in_organism, "_", in_response_var, "_plot_varcomp_page", i, "_",
          format(Sys.Date(), "%Y%m%d"), 
          ".png"
        )
      )
      message("Saving ", out_path)
      ggsave(out_path, 
             plot = plot_varcomp$plots[[i]], 
             width = 12, height = 9, units = "in", 
             dpi = 300)
      return(out_path)
    }, character(1))
    out_paths$plot_varcomp <- vc_paths
  }
  
  # Single-page plots
  single_plots <- list(
    plot_obs_preds = plot_obs_preds,
    plot_x_preds   = plot_x_preds,
    plot_emmeans   = plot_emmeans$plot,
    plot_emtrends  = plot_emtrends$plot
  )
  
  for (pname in names(single_plots)) {
    if (!is.null(single_plots[[pname]])) {
      out_path <- file.path(
        out_dir,
        paste0(
          "ssn_div_hydrowindow_", 
          in_organism, "_", in_response_var, "_", pname, "_",
          format(Sys.Date(), "%Y%m%d"), 
          ".png"
        )
      )
      message("Saving ", out_path)
      ggsave(out_path, 
             plot = single_plots[[pname]], 
             width = 12, height = 9, units = "in", 
             dpi = 300)
      out_paths[[pname]] <- out_path
    }
  }
  
  return(out_paths)
}

#------ plot_varcomp_multiorganisms --------------------------------------------
# facet: hydrological variable
# x-axis: % variance explained by fixed effects in addition to basin area
# y-axis: organism
# label: dominant time window

# in_hydrowindow_varcomp_multiorg = tar_read(
#   hydrowindow_varcomp_richness_multiorg)
# in_organism_dt <- tar_read(organism_dt)
# in_hydro_vars_dt <- tar_read(hydro_vars_dt)
# write_plot=T
# out_dir=figdir

plot_varcomp_multiorganisms <- function(in_hydrowindow_varcomp_multiorg,
                                        in_hydro_vars_dt,
                                        in_organism_dt,
                                        write_plot=F,
                                        out_dir=NULL
) {
  
  varcomp_best <- in_hydrowindow_varcomp_multiorg %>%
    .[varcomp=="Covariates (PR-sq)", 
      .SD[which.max(proportion),],
      by=.(hydro_label, organism)] %>%
    .[, null_fixedR2 := .SD[hydro_var=='null', proportion], by=organism] %>%
    .[, marginal_R2 := proportion - null_fixedR2] %>%
    .[!(hydro_label %in% c('Null', 'Mean flow percentile')),] %>%
    merge(in_organism_dt, by='organism') %>%
    merge(in_hydro_vars_dt[!duplicated(hydro_label),
                           .(hydro_label, hydro_class, metric_num)], 
          by='hydro_label') %>%
    .[, h_unit := fifelse(grepl(pattern='.*yrpast.*', hydro_var), 'y', 'd')]
  
  
  varcomp_plot <- ggplot(varcomp_best, 
                         aes(x=organism_class, y=marginal_R2, fill=organism_sub)) +
    geom_bar(stat = "identity", position='dodge', alpha=0.6) +
    geom_text(aes(y=Inf, hjust=0.1,label=paste(window_d, h_unit)), 
              position = position_dodge(width = .9), size=3, color="#777777") +
    scale_y_continuous(name='Marginal explained variance (%)',
                       breaks=c(0, 0.05, 0.10),
                       labels=scales::label_percent()) +
    scale_x_discrete(name='Organism') +
    scale_fill_manual(name='Subset',
                      values= c('#666666', '#7570b3', '#66a61e', '#a6761d', '#1b9e77')
    ) + 
    coord_flip(clip='off') +
    facet_wrap(~ hydro_class + str_wrap(hydro_label, 30), ncol = 4, 
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 list(do.call(paste, c(labels, list(sep = "\n"))))
               }) +
    theme_minimal() +
    theme(
      strip.text = element_text(margin = margin(t = 0, b = 0)), # shrink text padding
      panel.spacing = unit(1.5, "lines"),        # space between facets
      legend.box.margin = margin(l = 10, t = 0, r = 0, b = 0)
    ) 
  
  if (write_plot && !is.null(out_dir)) {
    out_path <- file.path(
      out_dir,
      paste0(
        "ssn_hydrowindow_varcomp_multiorganisms", 
        "_", unique(in_hydrowindow_varcomp_multiorg$response_var),  "_",
        format(Sys.Date(), "%Y%m%d"), 
        ".png"
      )
    )
    message("Saving ", out_path)
    ggsave(out_path, 
           plot = varcomp_plot, 
           width = 12, height = 9, units = "in", 
           dpi = 600)
  }
  
  return(list(
    dt = varcomp_best,
    plot = varcomp_plot
  ))
}
#------ plot_emtrends_multiorganisms --------------------------------------------
# -> em slopes across organisms
# facet: hydrological variable
# y-axis: organism
# x-axis: relative slope
# color: countries

# 
# in_hydro_vars_dt <- tar_read(hydro_vars_dt)
# in_organism_dt <- tar_read(organism_dt)
# 
# in_hydrowindow_best_intercept_dt <- tar_read(hydrowindow_best_intercept_dt)
# 
# emtrends_list = list(
#   miv_nopools = tar_read(hydrowindow_emtrends_richness_miv_nopools)$dt,
#   miv_nopools_ept = tar_read(hydrowindow_emtrends_richness_miv_nopools_ept)$dt,
#   miv_nopools_och = tar_read(hydrowindow_emtrends_richness_miv_nopools_och)$dt,
#   dia_biof_nopools = tar_read(hydrowindow_emtrends_richness_dia_biof_nopools)$dt,
#   dia_sedi_nopools = tar_read(hydrowindow_emtrends_richness_dia_sedi_nopools)$dt,
#   fun_biof_nopools = tar_read(hydrowindow_emtrends_richness_fun_biof_nopools)$dt,
#   fun_sedi_nopools = tar_read(hydrowindow_emtrends_richness_fun_sedi_nopools)$dt,
#   bac_biof_nopools = tar_read(hydrowindow_emtrends_richness_bac_biof_nopools)$dt,
#   bac_sedi_nopools = tar_read(hydrowindow_emtrends_richness_bac_sedi_nopools)$dt
# )
# in_drn_dt <- drn_dt

plot_emtrends_multiorganisms <- function(emtrends_list,
                                         in_hydrowindow_best_intercept_dt,
                                         in_hydro_vars_dt,
                                         in_organism_dt,
                                         in_drn_dt,
                                         write_plot = F,
                                         subset_variables = T, 
                                         out_dir = NULL) {
  
  emtrends_all <- rbindlist(emtrends_list, idcol='organism') 
  
  #Get hydrological variable labels and class
  emtrends_all <- emtrends_all %>%
    get_hydro_var_root(in_place=F) %>%
    merge(in_hydro_vars_dt[!duplicated(hydro_var_root), #Format to merge with hydro_vars_dt
                           .(hydro_var_root, hydro_label)],
          by='hydro_var_root') 
  
  #Normalize coefficients by intercept
  emtrends_all <- emtrends_all %>% 
    merge(in_hydrowindow_best_intercept_dt[, .(hydro_var_root, organism, intercept)], 
          by=c('hydro_var_root', 'organism')) %>%
    .[, `:=`(trend_rel = trend/abs(intercept),
             lcl_rel = asymp.LCL/abs(intercept),
             ucl_rel = asymp.UCL/abs(intercept)
    )] 
  
  #Cap ucl and lcl based on means
  emtrends_all[, `:=`(
    lcl_rel_capped = fifelse(lcl_rel < min(trend_rel), min(trend_rel), lcl_rel),
    ucl_rel_capped = fifelse(ucl_rel > max(trend_rel), max(trend_rel), ucl_rel)
  ), by = hydro_label]
  
  #Format country and labels for plotting
  emtrends_all <- emtrends_all  %>%
    .[!(hydro_label %in% c('Null', 'Mean flow percentile')),] %>% #Remove to have 16 facets
    merge(in_organism_dt, by='organism') %>% #Get organism labels and class
    .[, h_unit := fifelse(grepl(pattern='.*yrpast.*', hydro_var), 'y', 'd')] %>% #Get time window unit
    .[, organism_sub := gsub('Sediment', 'Sedi.', organism_sub)] %>%
    .[, organism_sub := gsub('Biofilm', 'Biof.', organism_sub)] %>%
    setorder('organism') 
  
  response_var <- unique(emtrends_all$response_var)
  
  color_vec <- emtrends_all[!duplicated(country),
                            setNames(color, country)]
  
  #Plot trend coefficients
  emtrends_plot <- ggplot(emtrends_all, aes(x = organism_class, 
                                            y = trend_rel,
                                            group = organism_label)) +
    geom_pointrange(aes(ymin = lcl_rel_capped, ymax = ucl_rel_capped, color= country),
                    alpha=0.5, fatten=4, position = position_dodge(0.7)) +
    geom_hline(yintercept=0, linetype=2) +
    geom_text(aes(y=Inf, hjust=0.1,
                  label=paste(organism_sub, '|', window_d, h_unit)), 
              position = position_dodge(width = .7), size=3, color="darkgrey") +
    scale_color_manual(name='Country',
                       values=color_vec) +
    labs(y = paste("Estimated slope of", response_var),
         x = 'Organism') +
    coord_flip(clip='off') +
    facet_wrap(~ hydro_class + str_wrap(hydro_label, 30), 
               ncol = 4, 
               labeller = function (labels) {
                 labels <- lapply(labels, as.character)
                 list(do.call(paste, c(labels, list(sep = "\n"))))
               },
               scales='free_x') +
    theme_minimal() +
    theme(
      strip.text = element_text(margin = margin(t = 0, b = 0)), # shrink text padding
      panel.spacing = unit(2.3, "lines"),        # space between facets
      legend.box.margin = margin(l = 20, t = 0, r = 0, b = 0)
    ) 
  
  
  #Subset variables
  if (subset_variables) {
    emtrends_sub <- emtrends_all[hydro_var_root %in% 
                                   c('DurD', 'FreD', 'Fdist_mean_undirected',
                                     'STcon_undirected', 'DurD_CV', 'FstDrE', 
                                     'sd6', 'oQ10'),]
    emtrends_plot_sub <- emtrends_plot + emtrends_sub
  }

  
  if (write_plot && !is.null(out_dir)) {
    out_path <- file.path(
      out_dir,
      paste0(
        "ssn_hydrowindow_emtrends_multiorganisms", 
        "_", response_var,  "_",
        format(Sys.Date(), "%Y%m%d"), 
        ".png"
      )
    )
    message("Saving ", out_path)
    ggsave(out_path, 
           plot = emtrends_plot, 
           width = 14, height = 11, units = "in", 
           dpi = 600)
    
    if (subset_variables) {
      ggsave(gsub(response_var, paste0(response_var, '_subvars'), out_path), 
             plot = emtrends_plot_sub, 
             width = 11, height = 7, units = "in", 
             dpi = 600)
      
    }
  }
  
  return(list(
    dt=emtrends_all,
    plot=emtrends_plot
  ))
  
}

#------ model_miv_richness_yr -----------------------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = F
# scale_predictors = T

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
#' @param in_allvars_summarized A list containing merged environmental and hydrological
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
#'   
model_miv_richness_yr <- function(in_ssn_eu_summarized,
                                  in_allvars_summarized,
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
  
  #Transform some variables  
  ssn_miv$obs %<>%
    mutate(meanQ3650past_sqrt = sqrt(meanQ3650past))
  # ssn_miv$preds$preds_hist %>%
  #   mutate(meanQ3650past_sqrt = sqrt(meanQ3650past))
  # ssn_miv$preds$preds_proj %>%
  #   mutate(meanQ3650past_sqrt = sqrt(meanQ3650past))
    
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'meanQ3650past_sqrt')
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_miv <- scale_ssn_predictors(in_ssn = ssn_miv,
                                  in_vars = hydro_candidates,
                                  scale_ssn_preds = FALSE) 
  
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
    
    ggplot(ssn_miv$obs, aes(x=basin_area_km2, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10()
    
    ggplot(ssn_miv$obs, aes(x=meanQ3650past, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_miv$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_miv$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_miv$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    data.table::melt(as.data.table(ssn_miv$obs), 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
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
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('miv_nopools'),
      formula_root = ' log10(basin_area_km2)',
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
      formula_root = ' log10(basin_area_km2)',
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
      formula_root = ' log10(basin_area_km2)',
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
    
    ssn_cov_glance[family=='nbinomial' & which.min(AICc),]
    # plot(ssn_norm_ini[[1]][[1]])
    # plot(ssn_nbin_ini[[1]][[1]])
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose negative binomial based on the AIC scores
    ssn_nbin_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('miv_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    SSN2::loocv(ssn_nbin_cov$ssn_list$none_none_none)
    SSN2::loocv(ssn_nbin_cov$ssn_list$none_none_spherical)
    SSN2::loocv(ssn_nbin_cov$ssn_list$none_mariah_none)
    lapply(list(
      ssn_nbin_cov$ssn_list$none_none_none,
      ssn_nbin_cov$ssn_list$none_none_spherical,
      ssn_nbin_cov$ssn_list$none_mariah_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_miv_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'nbinomial',
      ssn.object = ssn_miv,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'spherical',
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
  miv_rich_modlist[['null']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  miv_rich_modlist[['mod1']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod2']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod3']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ basin_area_km2 + country:log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod4']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  miv_rich_modlist[['mod5']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  miv_rich_modlist[['mod6']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod7']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ meanQ3650past + country:sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod8']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country*sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod9']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + country*sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod10']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past + sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod11']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod12']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + country*log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod13']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod14']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod15']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod16']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod17']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_avg_samp + country:DurD_avg_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod18']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_max_samp + country:DurD_max_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod19']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + country:FreD_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod20']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ meanConD_yr + country:meanConD_yr + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod21']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod22']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + FstDrE_SD30yrpast + country:FstDrE_SD30yrpast + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod23']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_avg_samp + country:DurD_avg_samp + FreD_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod24']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sd6_10yrpast + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod25']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ Fdist_mean_10past_undirected_avg_samp + country:Fdist_mean_10past_undirected_avg_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod26']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ STcon_m10_undirected_avg_samp + country:STcon_m10_undirected_avg_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod27']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + PFreD365past + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod28']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + STcon_m10_undirected_avg_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod29']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + STcon_m10_directed_avg_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod30']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp*sd6_10yrpast + country:DurD_samp + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod31']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + FstDrE + FstDrE:country + log10(basin_area_km2)'))) #Overspec
  
  miv_rich_modlist[['mod32']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + RelF_avg_samp + log10(basin_area_km2)'))) 
  
  miv_rich_modlist[['mod33']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + FstDrE + FstDrE:country + log10(basin_area_km2)'))) 
  
  miv_rich_modlist[['mod34']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod35']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + FstDrE + log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod36']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + FreD3650past + FreD3650past:country + sqrt(meanQ3650past)')))
  
  miv_rich_modlist[['mod37']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_CV10yrpast + DurD_CV10yrpast:country + log10(basin_area_km2)'))) 
  
  miv_rich_modlist[['mod38']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast + FstDrE_SD10yrpast:country + log10(basin_area_km2)'))) 
  
  miv_rich_modlist[['mod39']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past*PFreD365past + country:DurD3650past + sqrt(meanQ3650past)'))) 
  
  miv_rich_modlist[['mod40']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + PFreD365past + country:DurD3650past + country:PFreD365past + sqrt(meanQ3650past)'))) 
  
  miv_rich_modlist[['mod41']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ STcon_m10_undirected_avg_samp + country:STcon_m10_undirected_avg_samp + log10(basin_area_km2)'))) 
  
  miv_rich_modlist[['mod42']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + STcon_m10_undirected_avg_samp + country:STcon_m10_undirected_avg_samp + log10(basin_area_km2)'))) 
  
  miv_rich_modlist[['mod43']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + STcon_m10_undirected_avg_samp + country:STcon_m10_undirected_avg_samp + log10(basin_area_km2)'))) 
  #Overspecified
  
  miv_rich_modlist[['mod44']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(STcon_m10_directed_avg_samp) + country:sqrt(STcon_m10_directed_avg_samp)'))) 
  #Overspecified
 
  miv_rich_modlist[['mod45']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + sqrt(STcon_m10_directed_avg_samp) + log10(basin_area_km2)'))) 
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- miv_rich_modlist[['mod36']]
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    #Check residuals, variance decomposition, summary, and predictions
    plot(chosen_mod)
    SSN2::varcomp(chosen_mod)
    summary(chosen_mod)
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'FreD3650past',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'FreD3650past',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # get_ssn_emtrends(in_mod=chosen_mod , 
    #                  in_pred_var = 'sd6_10yrpast',
    #                  interaction_var = 'country', 
    #                  in_drn_dt = drn_dt,
    #                  plot=T)
    # 
    # get_ssn_emmeans(in_mod=chosen_mod , 
    #                 in_pred_var = 'FstDrE',
    #                 interaction_var = 'DurD_samp', 
    #                 in_drn_dt = drn_dt,
    #                 plot=T, label_pred_var=T)
  }
  
  #Add environmental variables and "stream_type"
  miv_rich_modlist[['mod46']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1')))
  
  miv_rich_modlist[['mod47']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1 + env_PC2')))
  
  miv_rich_modlist[['mod48']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC2')))
  
  miv_rich_modlist[['mod49']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1*env_PC2')))

  miv_rich_modlist[['mod50']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod51']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ stream_type + DurD_samp + country:DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  miv_rich_modlist[['mod52']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + log10(basin_area_km2) + country:log10(basin_area_km2)'))) 
  
  miv_rich_modlist[['mod53']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)'))) 
  
  miv_rich_modlist[['mod54']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + stream_type + log10(basin_area_km2) + country:log10(basin_area_km2) + env_PC2'))) 
  
  miv_rich_modlist[['mod55']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + stream_type + log10(basin_area_km2) + country:log10(basin_area_km2)'))) 
  
  # allvars_dt <- as.data.table(ssn_miv$obs) %>%
  #   setorderv(c('country', 'site')) 
  # 
  # allvars_dt[
  #   , (hydro_candidates) := lapply(.SD,
  #                                  function(x) base::scale(x, center=T, scale=T)),
  #   .SDcols = hydro_candidates]
  # 
  # check_rf <- ranger::ranger(mean_richness ~ ., 
  #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
  # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
  # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) + 
  #   geom_point() +
  #   facet_wrap(~country)
  
  #Summary table
  mod_perf_tab <- lapply(miv_rich_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(miv_rich_modlist)]
  
  #####CHOOSE MOD 52 for "absolute" best model###############"
  #####CHOOSE MOD 36 for continuous predictions###########
  
  #------ For "best" model -------------------------------------------------
  # Refits the selected "best" model using the REML method for better variance estimation.
  ssn_miv_best_final <- quick_miv_ssn(
    as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type +
      log10(basin_area_km2) + country:log10(basin_area_km2)')),
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
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ DurD_samp + FreD3650past + FreD3650past:country + sqrt(meanQ3650past)',
    ssn.object = ssn_miv,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'nbinomial',
      organism = c('miv_nopools'),
      formula_root = 'sqrt(meanQ3650past)',
      hydro_var = 'FreD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
    SSN2::varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    SSN2::varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_mariah_gaussian)
    SSN2::varcomp(ssn_mod_predictions_covtypes$ssn_list$none_mariah_gaussian)
  }
  
  #Re-fit final model with REML
  ssn_miv_pred_final <- ssn_glm(
    formula = mean_richness ~ DurD_samp_z + FreD3650past_z + FreD3650past_z:country + meanQ3650past_sqrt_z,
    family = 'nbinomial',
    ssn.object = ssn_miv,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'spherical',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_miv_mod_preds <- augment(ssn_miv_pred_final, drop=F, type.predict = 'response')
  
  if (interactive) {
    tidy(miv_rich_modlist[['mod36']])
    tidy(ssn_miv_pred_final)
    summary(ssn_miv_pred_final)
    plot(ssn_miv_pred_final)
    SSN2::varcomp(ssn_miv_pred_final)
    
    #Check predictions
    ggplot(data=ssn_miv_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
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

#------ model_miv_invsimpson_yr ------------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

model_miv_invsimpson_yr <- function(in_ssn_eu_summarized,
                                    in_allvars_summarized,
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
  
  ssn_miv$obs$meanQ3650past_sqrt <- sqrt(ssn_miv$obs$meanQ3650past)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'meanQ3650past_sqrt')

  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_miv <- scale_ssn_predictors(in_ssn = ssn_miv,
                                  in_vars = hydro_candidates,
                                  scale_ssn_preds = FALSE) 
  
  alpha_var <- 'mean_invsimpson'
  
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
    
    ggplot(ssn_miv$obs, aes(x=basin_area_km2, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_miv$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
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
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('miv_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_gamma_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gamma',
      organism = c('miv_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_lognorm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('miv_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      ssn_norm_ini$ssn_glance,
      ssn_gamma_ini$ssn_glance,
      ssn_lognorm_ini$ssn_glance
    ))
    
    #CHoose lognormal
    ssn_lognorm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('miv_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    ssn_lognorm_cov$ssn_glance
    
    SSN2::loocv(ssn_lognorm_cov$ssn_list$none_none_none)
    SSN2::loocv(ssn_lognorm_cov$ssn_list$none_none_gaussian)
    SSN2::loocv(ssn_lognorm_cov$ssn_list$none_none_spherical)
    SSN2::loocv(ssn_lognorm_cov$ssn_list$epa_none_none)
    lapply(list(
      ssn_lognorm_cov$ssn_list$none_none_none,
      ssn_lognorm_cov$ssn_list$none_none_spherical,
      ssn_lognorm_cov$ssn_list$none_none_gaussian,
      ssn_lognorm_cov$ssn_list$epa_none_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_miv_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_miv,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'spherical',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = in_random,
      estmethod = estmethod
    )
  }
  
  #Then test multiple models -----
  # This is the main model selection loop, testing various combinations of
  # predictors.
  ssn_miv$obs$mean_invsimpson_log10 <- log10(ssn_miv$obs$mean_invsimpson)
  alpha_var <- 'mean_invsimpson_log10'
  miv_simps_modlist <- list()
  
  #Null model
  miv_simps_modlist[['null']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  miv_simps_modlist[['mod1']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  miv_simps_modlist[['mod2']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  miv_simps_modlist[['mod3']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ basin_area_km2 + country:log10(basin_area_km2)')))
  
  miv_simps_modlist[['mod4']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  miv_simps_modlist[['mod5']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  miv_simps_modlist[['mod6']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(meanQ3650past)')))
  
  miv_simps_modlist[['mod7']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ meanQ3650past + country:sqrt(meanQ3650past)')))
  
  miv_simps_modlist[['mod8']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ meanQ3650past + country:log10(meanQ3650past)')))
  
  miv_simps_modlist[['mod9']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + country:FreD_samp + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  miv_simps_modlist[['mod10']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  miv_simps_modlist[['mod11']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  miv_simps_modlist[['mod12']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  miv_simps_modlist[['mod13']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  miv_simps_modlist[['mod14']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  #Add environmental variables and "stream_type"
  miv_simps_modlist[['mod15']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + sqrt(meanQ3650past) + country:sqrt(meanQ3650past) + env_PC1')))
  
  miv_simps_modlist[['mod16']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + sqrt(meanQ3650past) + country:sqrt(meanQ3650past) + env_PC1 + env_PC2')))
  
  miv_simps_modlist[['mod17']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + sqrt(meanQ3650past) + country:sqrt(meanQ3650past) + env_PC2')))
  
  miv_simps_modlist[['mod18']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  miv_simps_modlist[['mod19']] <- quick_miv_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  #Summary table
  mod_perf_tab <- lapply(miv_simps_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(miv_simps_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- miv_simps_modlist[['mod10']]
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    #Check residuals, variance decomposition, summary, and predictions
    plot(chosen_mod)
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=log10(mean_invsimpson))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'FreD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'FreD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
  }
    
    #RF does not do so well
    # check_rf <- ranger::ranger(mean_invsimpson ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_invsimpson', 'country'), with=F])
    # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness)) +
    #   geom_point(aes(color=stream_type)) + 
    #   geom_smooth(method='gam') +
    #   facet_wrap(~country)
    
    #------ For prediction model -------------------------------------------------
    # This section focuses on creating a model suitable for prediction to new reaches,
    # using only variables available for all reaches.
    
    #Check Torgegram
    tg_selected <- SSN2::Torgegram(
      formula = 'mean_invsimpson_log10 ~ PDurD365past + country:PDurD365past + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)',
      ssn.object = ssn_miv,
      type = c("flowcon", "flowuncon", "euclid"),
      partition_factor = as.formula('~ country'),
      bins = 15
    ) %>%
      rbindlist(idcol='dist_type')
    
    tg_plot <- ggplot(tg_selected, 
                      aes(x=dist, y=gamma, color = dist_type)) +
      geom_point(aes(size=np)) +
      geom_smooth(method='lm', linewidth=2, se=F) +
      scale_size_continuous(range=c(0.5, 7)) +
      scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
      theme_classic() +
      facet_wrap(~dist_type)
    
    if (interactive) {
      #Test covariance structures again, this time with REML 
      in_ssn_eu_summarized$miv_nopools$ssn$obs %<>%
        mutate(meanQ3650past_sqrt = sqrt(meanQ3650past),
               mean_invsimpson_log10 = log10(mean_invsimpson))
      
      ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
        in_ssn = in_ssn_eu_summarized,
        family = 'Gaussian',
        organism = c('miv_nopools'),
        formula_root = 'sqrt(meanQ3650past)',
        hydro_var = 'PDurD365past',
        response_var = alpha_var,
        ssn_covtypes = ssn_covtypes,
        partition_formula = as.formula('~ country'),
        random_formula = NULL,
        estmethod='reml'
      )
      
      selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
      
      summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
      varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
      
      summary(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
      varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    }
    
    #Re-fit final model with REML
    ssn_miv_pred_final <- ssn_lm(
      formula = mean_invsimpson_log10 ~ PDurD365past_z + country:PDurD365past_z + 
        meanQ3650past_sqrt_z + country:meanQ3650past_sqrt_z,
      family = 'Gaussian',
      ssn.object = ssn_miv,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'none',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = NULL,
      estmethod = 'reml'
    )
    
    ssn_miv_mod_preds <- augment(ssn_miv_pred_final, drop=F, 
                                 type.predict = 'response')
    
    if (interactive) {
      tidy(ssn_miv_pred_final)
      summary(ssn_miv_pred_final)
      plot(ssn_miv_pred_final)
      SSN2::varcomp(ssn_miv_pred_final)
      
      get_ssn_emmeans(in_mod=ssn_miv_pred_final, 
                      in_pred_var = 'PDurD365past',
                      interaction_var = 'country', 
                      in_drn_dt = drn_dt,
                      plot=T, label_pred_var=T)$plot
      
      
      get_ssn_emtrends(in_mod=ssn_miv_pred_final, 
                       in_pred_var = 'PDurD365past',
                       interaction_var = 'country', 
                       in_drn_dt = drn_dt,
                       plot=T)$plot
      
      #Check predictions
      ggplot(data=ssn_miv_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
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
    torgegram_mod_pred = tg_plot #Torgegram for prediction model
  ))
}

#------ model_ept_richness_yr --------------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

model_ept_richness_yr <- function(in_ssn_eu_summarized,
                                  in_allvars_summarized,
                                  in_cor_matrices, 
                                  ssn_covtypes,
                                  scale_predictors = T,
                                  interactive = F) {
  
  in_organism <- 'miv_nopools_ept'
  alpha_var <- 'mean_richness'
  
  #Subset SSN and create distance matrices-----
  ssn_miv <- in_ssn_eu_summarized$miv_nopools_ept$ssn
  SSN2::ssn_create_distmat(
    ssn_miv,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_miv$obs$basin_area_km2_log10 <- log10(ssn_miv$obs$basin_area_km2)
  ssn_miv$obs$meanQ3650past_log10 <- log10(ssn_miv$obs$meanQ3650past)
  
  # ssn_miv$preds$preds_hist$basin_area_km2_log10 <- log10(ssn_miv$preds$preds_hist$basin_area_km2)
  # ssn_miv$preds$preds_proj$basin_area_km2_log10 <- log10(ssn_miv$preds$preds_proj$basin_area_km2)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'basin_area_km2_log10',
                        'meanQ3650past_log10')
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_miv <- scale_ssn_predictors(in_ssn = ssn_miv,
                                  in_vars = hydro_candidates,
                                  scale_ssn_preds = FALSE) 
  
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
    
    ggplot(ssn_miv$obs, aes(x=meanQ3650past, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_miv$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_miv$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
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
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == in_organism 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == in_organism &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == in_organism &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    #Use negative binomial, like for all miv
    
    #CHoose negative binomial based on the AIC scores
    ssn_nbin_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c(in_organism),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = ~ country,
      estmethod='reml'
    )
    View(ssn_nbin_cov$ssn_glance)
    
    SSN2::loocv(ssn_nbin_cov$ssn_list$none_none_spherical)
    SSN2::loocv(ssn_nbin_cov$ssn_list$linear_none_none)
    SSN2::loocv(ssn_nbin_cov$ssn_list$epa_none_none)
    lapply(list(
      ssn_nbin_cov$ssn_list$none_none_spherical,
      ssn_nbin_cov$ssn_list$linear_none_none,
      ssn_nbin_cov$ssn_list$epa_none_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_ept_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'nbinomial',
      ssn.object = ssn_miv,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'spherical',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = in_random,
      estmethod = estmethod
    )
  }
  
  #Then test multiple models -----
  # This is the main model selection loop, testing various combinations of
  # predictors.
  ept_rich_modlist <- list()
  
  #Null model
  ept_rich_modlist[['null']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Null model with country random effect
  ept_rich_modlist[['mod1']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ 1')), in_random=~country)
  
  #Initial model
  ept_rich_modlist[['mod2']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod3']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ basin_area_km2 + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod4']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  ept_rich_modlist[['mod5']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod6']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod7']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod8']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod9']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod10']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ Fdist_mean_10past_undirected_avg_samp + country:Fdist_mean_10past_undirected_avg_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod11']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ STcon_m10_undirected_avg_samp + country:STcon_m10_undirected_avg_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod12']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod13']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod14']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod15']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod16']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod17']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past +  sd6_30yrpast + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  #sd6_30yrpast; se too large
  
  ept_rich_modlist[['mod18']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + PmeanQ10past_max_samp + country:PmeanQ10past_max_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod19']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + DurD3650past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod20']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past +  FstDrE + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod21']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + FreD3650past + FreD3650past + country:log10(meanQ3650past)')))
  
  ept_rich_modlist[['mod22']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + PmeanQ10past_max_samp + country:PmeanQ10past_max_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  #PmeanQ10past_max: standard errors too large
  
  ept_rich_modlist[['mod23']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + FreD_CV10yrpast + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  #Add environmental variables and "stream_type"
  ept_rich_modlist[['mod24']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1')))
  
  ept_rich_modlist[['mod25']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1 + env_PC2')))
  
  ept_rich_modlist[['mod26']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC2')))
  
  ept_rich_modlist[['mod27']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1*env_PC2')))
  
  ept_rich_modlist[['mod28']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod29']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ stream_type + DurD_samp + country:DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  ept_rich_modlist[['mod30']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + log10(basin_area_km2) + country:log10(basin_area_km2)'))) 
  
  ept_rich_modlist[['mod31']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + log10(meanQ3650past) + country:log10(meanQ3650past)'))) 
  
  ept_rich_modlist[['mod32']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + stream_type + log10(basin_area_km2) + country:log10(basin_area_km2) + env_PC2'))) 
  
  ept_rich_modlist[['mod33']] <- quick_ept_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + stream_type + log10(basin_area_km2) + country:log10(basin_area_km2)'))) 
  
  # check_rf <- ranger::ranger(mean_richness ~ .,
  #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
  # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
  # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) +
  #   geom_point() +
  #   facet_wrap(~country)
  # 
  #Summary table
  mod_perf_tab <- lapply(ept_rich_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(ept_rich_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- ept_rich_modlist[['mod9']]
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    #Check residuals, variance decomposition, summary, and predictions
    plot(chosen_mod)
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)$plot
    
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)$plot
    
    get_ssn_emtrends(chosen_mod, 
                     in_pred_var = 'meanQ3650past', 
                     interaction_var='country',
                     in_drn_dt = drn_dt)$plot
    
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  
  #####CHOOSE MOD 30 for "absolute" best model###############"
  #####CHOOSE MOD 9 for continuous predictions###########
  
  
  #------ For "best" model -------------------------------------------------
  # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_ept_best_final <- quick_ept_ssn(
  #   as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type + 
  #     log10(basin_area_km2) + country:log10(basin_area_km2)')),
  #   estmethod = 'reml'
  # )
  # ssn_ept_best_preds <- augment(ssn_ept_best_final,
  #                               drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_ept_best_final)
  #   SSN2::varcomp(ssn_ept_best_final)
  #   plot(ssn_ept_best_final)
  # }
  # 
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ DurD_samp + country:DurD_samp + log10(meanQ3650past) + country:log10(meanQ3650past)',
    ssn.object = ssn_miv,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'nbinomial',
      organism = c(in_organism),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
  }
  
  #Re-fit final model with REML
  ssn_ept_pred_final <- ssn_glm(
    formula = mean_richness ~ DurD_samp_z + country:DurD_samp_z + 
      meanQ3650past_log10_z + country:meanQ3650past_log10_z,
    family = 'nbinomial',
    ssn.object = ssn_miv,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'spherical',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_ept_mod_preds <- augment(ssn_ept_pred_final, drop=F, 
                               type.predict = 'response')
  
  if (interactive) {
    tidy(ssn_ept_pred_final)
    summary(ssn_ept_pred_final)
    plot(ssn_ept_pred_final)
    SSN2::varcomp(ssn_ept_pred_final)
    
    get_ssn_emmeans(in_mod=ssn_ept_pred_final, 
                    in_pred_var = 'DurD_samp_z',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)$plot
    
    
    get_ssn_emtrends(in_mod=ssn_ept_pred_final, 
                     in_pred_var = 'DurD_samp_z',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)$plot
    
    #Check predictions
    ggplot(data=ssn_ept_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_ept_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_ept_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    # , ssn_mod_fit_best = ssn_ept_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_ept_best_preds #Augmented data for "best" model
  ))
}

#------ model_och_richness_yr --------------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = F
# scale_predictors = T

model_och_richness_yr <- function(in_ssn_eu_summarized,
                                  in_allvars_summarized,
                                  in_cor_matrices, 
                                  ssn_covtypes,
                                  scale_predictors = T,
                                  interactive = F) {
  
  in_organism <- 'miv_nopools_och'
  alpha_var <- 'mean_richness'
  
  #Subset SSN and create distance matrices-----
  ssn_miv <- in_ssn_eu_summarized$miv_nopools_och$ssn
  SSN2::ssn_create_distmat(
    ssn_miv,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_miv$obs$basin_area_km2_log10 <- log10(ssn_miv$obs$basin_area_km2)
  ssn_miv$obs$meanQ3650past_log10 <- log10(ssn_miv$obs$meanQ3650past)
  
  # ssn_miv$preds$preds_hist$basin_area_km2_log10 <- log10(ssn_miv$preds$preds_hist$basin_area_km2)
  # ssn_miv$preds$preds_proj$basin_area_km2_log10 <- log10(ssn_miv$preds$preds_proj$basin_area_km2)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'basin_area_km2_log10',
                        'meanQ3650past_log10')
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_miv <- scale_ssn_predictors(in_ssn = ssn_miv,
                                  in_vars = hydro_candidates,
                                  scale_ssn_preds = FALSE) 
  
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
    
    ggplot(ssn_miv$obs, aes(x=meanQ3650past, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_miv$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_miv$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    ggplot(ssn_miv$obs, aes(x=FstDrE, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm')+
      facet_wrap(~country)
    
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
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == in_organism 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == in_organism &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == in_organism &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    #Use negative binomial, like for all miv
    
    #CHoose negative binomial based on the AIC scores
    ssn_nbin_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c(in_organism),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = ~ country,
      estmethod='reml'
    )
    View(ssn_nbin_cov$ssn_glance)
    
    SSN2::loocv(ssn_nbin_cov$ssn_list$spherical_none_none)
    SSN2::loocv(ssn_nbin_cov$ssn_list$none_none_none)
    lapply(list(
      ssn_nbin_cov$ssn_list$spherical_none_none,
      ssn_nbin_cov$ssn_list$none_none_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_och_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'nbinomial',
      ssn.object = ssn_miv,
      taildown_type = 'none',
      tailup_type = 'none',
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
  och_rich_modlist <- list()
  
  #Null model
  och_rich_modlist[['null']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Null model with country random effect
  och_rich_modlist[['mod1']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ 1')), in_random=~country)
  
  #Initial model
  och_rich_modlist[['mod2']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  och_rich_modlist[['mod3']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ basin_area_km2 + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod4']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  och_rich_modlist[['mod5']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  och_rich_modlist[['mod6']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod7']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod8']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod9']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod10']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ Fdist_mean_10past_undirected_avg_samp + country:Fdist_mean_10past_undirected_avg_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod11']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ STcon_m10_undirected_avg_samp + country:STcon_m10_undirected_avg_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod12']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod13']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod14']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod15']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod16']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod17']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past +  sd6_30yrpast + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  #sd6_30yrpast; se too large
  
  och_rich_modlist[['mod18']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FreD3650past + country:FreD3650past + PmeanQ10past_max_samp + country:PmeanQ10past_max_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod19']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + DurD3650past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod20']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past +  FstDrE + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod21']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + FreD3650past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  och_rich_modlist[['mod22']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + FreD3650past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  #PmeanQ10past_max: standard errors too large
  
  och_rich_modlist[['mod23']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + FreD_CV10yrpast + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  #Add environmental variables and "stream_type"
  och_rich_modlist[['mod24']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1')))
  
  och_rich_modlist[['mod25']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1 + env_PC2')))
  
  och_rich_modlist[['mod26']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC2')))
  
  och_rich_modlist[['mod27']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + env_PC1*env_PC2')))
  
  och_rich_modlist[['mod28']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ stream_type + country:stream_type + DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod29']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ stream_type + DurD_samp + country:DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  och_rich_modlist[['mod30']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + log10(basin_area_km2) + country:log10(basin_area_km2)'))) 
  
  och_rich_modlist[['mod31']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + stream_type + country:stream_type + log10(meanQ3650past) + country:log10(meanQ3650past)'))) 
  
  och_rich_modlist[['mod32']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + stream_type + log10(basin_area_km2) + country:log10(basin_area_km2) + env_PC2'))) 
  
  och_rich_modlist[['mod33']] <- quick_och_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + stream_type + log10(basin_area_km2) + country:log10(basin_area_km2)'))) 
  
  
  allvars_dt <- as.data.table(ssn_miv$obs)
  check_rf <- ranger::ranger(mean_richness ~ .,
                             data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
  allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
  ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) +
    geom_point() +
    geom_abline() +
    facet_wrap(~country)
  # 
  #Summary table
  mod_perf_tab <- lapply(och_rich_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(och_rich_modlist)]
  
  if (interactive) {
    #Check model in additional dochh (applied to most models before continuing)
    chosen_mod <- och_rich_modlist[['mod15']]
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    #Check residuals, variance decomposition, summary, and predictions
    plot(chosen_mod)
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)$plot
    
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)$plot
    
    get_ssn_emtrends(chosen_mod, 
                     in_pred_var = 'meanQ3650past', 
                     interaction_var='country',
                     in_drn_dt = drn_dt)$plot
    
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  
  #####CHOOSE MOD 15 for continuous predictions###########
  
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)',
    ssn.object = ssn_miv,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'nbinomial',
      organism = c(in_organism),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
  }
  
  #Re-fit final model with REML
  ssn_och_pred_final <- ssn_glm(
    formula = mean_richness ~ DurD_samp_z + country:DurD_samp_z + 
      basin_area_km2_log10_z + country:basin_area_km2_log10_z,
    family = 'nbinomial',
    ssn.object = ssn_miv,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_och_mod_preds <- augment(ssn_och_pred_final, drop=F, 
                               type.predict = 'response')
  
  if (interactive) {
    tidy(ssn_och_pred_final)
    summary(ssn_och_pred_final)
    plot(ssn_och_pred_final)
    SSN2::varcomp(ssn_och_pred_final)
    
    get_ssn_emmeans(in_mod=ssn_och_pred_final, 
                    in_pred_var = 'DurD_samp_z',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)$plot
    
    
    get_ssn_emtrends(in_mod=ssn_och_pred_final, 
                     in_pred_var = 'DurD_samp_z',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)$plot
    
    #Check predictions
    ggplot(data=ssn_och_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_och_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_och_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
  ))
}
#------ model_dia_sedi_invsimpson_yr --------------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

model_dia_sedi_invsimpson_yr <- function(in_ssn_eu_summarized,
                                         in_allvars_summarized,
                                         in_cor_matrices, 
                                         ssn_covtypes,
                                         scale_predictors = T,
                                         interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_dia_sedi <- in_ssn_eu_summarized$dia_sedi_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_dia_sedi,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_dia_sedi$obs$meanQ3650past_log10 <- log10(ssn_dia_sedi$obs$meanQ3650past)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'meanQ3650past_log10')
  
  ssn_dia_sedi <- scale_ssn_predictors(in_ssn = ssn_dia_sedi,
                                       in_vars = c(hydro_candidates,
                                                   'meanQ3650past_log10'),
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_invsimpson'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_dia_sedi$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_dia_sedi$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_dia_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_dia_sedi$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_dia_sedi$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    data.table::melt(ssn_dia_sedi$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'dia_sedi_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('dia_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_gamma_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gamma',
      organism = c('dia_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_lognorm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('dia_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      ssn_norm_ini$ssn_glance,
      ssn_gamma_ini$ssn_glance,
      ssn_lognorm_ini$ssn_glance
    ))
    
    norm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='mean_invsimpson',][which.min(AICc),]
    lognorm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='log10(mean_invsimpson)',][which.min(AICc),]
    
    plot(ssn_lognorm_ini[[1]][[1]])
    plot(ssn_norm_ini[[1]][[1]])
    
    # plot(ssn_nbin_ini[[1]][[1]])
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose lognormal
    ssn_lognorm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('dia_sedi_nopools'),
      formula_root = ' sqrt(meanQ3650past)',
      hydro_var = 'DurD3650past',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_lognorm_cov$ssn_glance)
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_dia_sedi_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_dia_sedi,
      taildown_type = 'none',
      tailup_type = 'none',
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
  ssn_dia_sedi$obs$mean_invsimpson_log10 <- log10(ssn_dia_sedi$obs$mean_invsimpson)
  alpha_var <- 'mean_invsimpson_log10'
  dia_sedi_simps_modlist <- list()
  
  #Null model
  dia_sedi_simps_modlist[['null']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  dia_sedi_simps_modlist[['mod1']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  dia_sedi_simps_modlist[['mod2']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ basin_area_km2 + country:log10(basin_area_km2)')))
  
  dia_sedi_simps_modlist[['mod3']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  dia_sedi_simps_modlist[['mod4']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_simps_modlist[['mod5']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast + country:FstDrE_SD10yrpast + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_simps_modlist[['mod6']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast + country:FstDrE_SD10yrpast + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_simps_modlist[['mod7']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FreD_CV30yrpast + country:FreD_CV30yrpast + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_simps_modlist[['mod8']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FreD_CV30yrpast + country:FreD_CV30yrpast + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_simps_modlist[['mod9']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + country:FreD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_simps_modlist[['mod10']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + country:FreD_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_simps_modlist[['mod11']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PFreD365past + country:PFreD365past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_simps_modlist[['mod12']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast + country:FstDrE_SD10yrpast')))
  
  dia_sedi_simps_modlist[['mod13']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FreD_CV30yrpast + country:FreD_CV30yrpast')))
  
  dia_sedi_simps_modlist[['mod14']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + country:FreD_samp')))
  
  dia_sedi_simps_modlist[['mod15']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PFreD365past + country:PFreD365past')))
  
  dia_sedi_simps_modlist[['mod16']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past')))
  
  dia_sedi_simps_modlist[['mod17']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past')))
  
  dia_sedi_simps_modlist[['mod18']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast*sqrt(qsim_avg_samp) + country:FstDrE_SD10yrpast')))
  
  dia_sedi_simps_modlist[['mod19']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast + country:FstDrE_SD10yrpast + meanConD_CV30yrpast')))
  
  dia_sedi_simps_modlist[['mod20']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast + country:FstDrE_SD10yrpast + meanConD_CV30yrpast  + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_simps_modlist[['mod21']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE_SD10yrpast + country:FstDrE_SD10yrpast + meanConD_CV30yrpast  + sqrt(qsim_avg_samp) + country:sqrt(qsim_avg_samp)')))
  
  mod_perf_tab <- lapply(dia_sedi_simps_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(dia_sedi_simps_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- dia_sedi_simps_modlist[['mod20']]
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    #Check residuals, variance decomposition, summary, and predictions
    plot(chosen_mod)
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    ggplot(data=augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'meanQ3650past',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'meanQ3650past',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    #Model with RF to see
    # allvars_dt <- as.data.table(in_ssn$obs) %>%
    #   setorderv(c('country', 'site')) 
    # 
    # allvars_dt[
    #   , (in_vars) := lapply(.SD, function(x) base::scale(x, center=T, scale=T)),
    #   .SDcols = in_vars]
    #
    # check_rf <- ranger::ranger(mean_invsimpson ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_invsimpson', 'country'), with=F])
    # allvars_dt[, mean_simps_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_simps_preds, y=mean_invsimpson, color=stream_type)) + 
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 20 for continuous predictions###########
  
  # #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_dia_sedi_best_final <- quick_dia_sedi_ssn(
  #   as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type +
  #     log10(basin_area_km2) + country:log10(basin_area_km2)')),
  #   estmethod = 'reml'
  # )
  # ssn_dia_sedi_best_preds <- augment(ssn_dia_sedi_best_final,
  #                               drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_dia_sedi_best_final)
  #   SSN2::varcomp(ssn_dia_sedi_best_final)
  #   plot(ssn_dia_sedi_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_invsimpson ~ FstDrE_SD10yrpast + country:FstDrE_SD10yrpast + meanConD_CV30yrpast  + log10(meanQ3650past) + country:log10(meanQ3650past)',
    ssn.object = ssn_dia_sedi,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive)  {
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('dia_sedi_nopools'),
      formula_root = 'meanConD_CV30yrpast',
      hydro_var = 'FstDrE_SD10yrpast',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_epa_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_epa_none)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_dia_sedi_pred_final <- ssn_glm(
    formula = mean_invsimpson_log10 ~ FstDrE_SD10yrpast_z + country:FstDrE_SD10yrpast_z + 
      meanConD_CV30yrpast_z  + meanQ3650past_log10_z + country:meanQ3650past_log10_z,
    family = 'Gaussian',
    ssn.object = ssn_dia_sedi,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_dia_sedi_mod_preds <- augment(ssn_dia_sedi_pred_final, drop=F)
  
  if (interactive) {
    tidy(dia_sedi_simps_modlist[['mod20']])
    tidy(ssn_dia_sedi_pred_final)
    summary(ssn_dia_sedi_pred_final)
    plot(ssn_dia_sedi_pred_final)
    SSN2::varcomp(ssn_dia_sedi_pred_final)
    
    #Check predictions
    ggplot(data=ssn_dia_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_dia_sedi_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_dia_sedi_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    # , ssn_mod_fit_best = ssn_dia_sedi_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_dia_sedi_best_preds #Augmented data for "best" model
  ))
}


#------ model_dia_sedi_richness_yr -----------------------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

#' Diatom model for annual data
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
#' @param in_allvars_summarized A list containing merged environmental and hydrological
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
#'   
model_dia_sedi_richness_yr <- function(in_ssn_eu_summarized,
                                       in_allvars_summarized,
                                       in_cor_matrices, 
                                       ssn_covtypes,
                                       scale_predictors = T,
                                       interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_dia_sedi <- in_ssn_eu_summarized$dia_sedi_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_dia_sedi,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_dia_sedi$obs$meanQ3650past_log10 <- log10(ssn_dia_sedi$obs$meanQ3650past)
  # ssn_dia_sedi$preds$preds_hist$meanQ3650past_log10 <- log10(ssn_dia_sedi$preds$preds_hist$meanQ3650past)
  # ssn_dia_sedi$preds$preds_proj$meanQ3650past_log10 <- log10(ssn_dia_sedi$preds$preds_proj$meanQ3650past)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'meanQ3650past_log10')
  
  ssn_dia_sedi <- scale_ssn_predictors(in_ssn = ssn_dia_sedi,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_richness'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_dia_sedi$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_dia_sedi$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_dia_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_sqrt() 
    
    ggplot(ssn_dia_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_dia_sedi$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_dia_sedi$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_dia_sedi$obs, aes(x=sqrt(STcon_m10_directed_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    data.table::melt(ssn_dia_sedi$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'dia_sedi_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('dia_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'PDurD365past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_pois_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",  #Gamma(link = "log"),
      organism = c('dia_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'PDurD365past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_nbin_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('dia_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'PDurD365past',
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
    
    ssn_cov_glance[family=='poisson' & which.min(AICc),]
    plot(ssn_norm_ini[[1]][[1]])
    plot(ssn_pois_ini[[1]][[1]])
    plot(ssn_nbin_ini[[1]][[1]])
    
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose negative binomial based on the AIC scores
    ssn_pois_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",
      organism = c('dia_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'PDurD365past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_pois_cov$ssn_glance)
    
    SSN2::loocv(ssn_pois_cov$ssn_list$none_exponential_none)
    SSN2::loocv(ssn_pois_cov$ssn_list$none_linear_none)
    SSN2::loocv(ssn_pois_cov$ssn_list$none_none_spherical)
    SSN2::loocv(ssn_pois_cov$ssn_list$none_none_none)
    
    lapply(list(
      ssn_pois_cov$ssn_list$none_exponential_none,
      ssn_pois_cov$ssn_list$none_linear_none,
      ssn_pois_cov$ssn_list$none_epa_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_dia_sedi_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'poisson',
      ssn.object = ssn_dia_sedi,
      taildown_type = 'none',
      tailup_type = 'none',
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
  dia_sedi_rich_modlist <- list()
  
  #Null model
  dia_sedi_rich_modlist[['null']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  dia_sedi_rich_modlist[['mod1']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod2']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod3']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod4']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod5']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod6']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  dia_sedi_rich_modlist[['mod7']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod8']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod9']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod10']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod11']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod12']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PFreD365past + country:PFreD365past + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod13']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod14']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod15']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ Fdist_mean_10past_undirected_avg_samp + country:Fdist_mean_10past_undirected_avg_samp + log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_sedi_rich_modlist[['mod16']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod17']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PFreD365past + country:PFreD365past + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod18']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod19']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod20']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ Fdist_mean_10past_undirected_avg_samp + country:Fdist_mean_10past_undirected_avg_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_sedi_rich_modlist[['mod21']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + sqrt(STcon_m10_directed_avg_samp) + log10(meanQ3650past) + country:log10(meanQ3650past) ')))
  
  mod_perf_tab <- lapply(dia_sedi_rich_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(dia_sedi_rich_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- dia_sedi_rich_modlist[['mod11']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'PDurD365past',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'PDurD365past',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'meanQ3650past',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
  }
  
  #Add environmental variables and "stream_type"
  dia_sedi_rich_modlist[['mod22']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(meanQ3650past) + country:log10(meanQ3650past) + env_PC1')))
  
  dia_sedi_rich_modlist[['mod23']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(meanQ3650past) + country:log10(meanQ3650past) + env_PC2')))
  
  dia_sedi_rich_modlist[['mod24']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(meanQ3650past) + country:log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_sedi_rich_modlist[['mod25']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(meanQ3650past) + country:log10(meanQ3650past) + stream_type')))
  
  dia_sedi_rich_modlist[['mod26']] <- quick_dia_sedi_ssn(as.formula(paste(alpha_var,  '~ PDurD365past + country:PDurD365past + log10(meanQ3650past) + country:log10(meanQ3650past) + stream_type + log10(meanQ3650past):stream_type')))
  
  # check_rf <- ranger::ranger(mean_richness ~ ., 
  #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
  # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
  # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) + 
  #   geom_point() +
  #   facet_wrap(~country)
  
  #Summary table
  mod_perf_tab <- lapply(dia_sedi_rich_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(dia_sedi_rich_modlist)]
  
  #####CHOOSE MOD 25 for "absolute" best model###############"
  #####CHOOSE MOD 11 for continuous predictions###########
  
  #------ For "best" model -------------------------------------------------
  # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_dia_sedi_best_final <- quick_dia_sedi_ssn(
  #   as.formula(paste(alpha_var, '~ PDurD365past + country:PDurD365past + 
  #                    log10(meanQ3650past) + country:log10(meanQ3650past) +
  #                    env_PC1 + env_PC2')),
  #   estmethod = 'reml'
  # )
  # ssn_dia_sedi_best_preds <- augment(ssn_dia_sedi_best_final,
  #                               drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_dia_sedi_best_final)
  #   SSN2::varcomp(ssn_dia_sedi_best_final)
  #   plot(ssn_dia_sedi_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ PDurD365past + country:PDurD365past + log10(meanQ3650past) + country:log10(meanQ3650past)',
    ssn.object = ssn_dia_sedi,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'poisson',
      organism = c('dia_sedi_nopools'),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'PDurD365past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none)
    # plot(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none, which=1:6)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_dia_sedi_pred_final <- ssn_glm(
    formula = mean_richness ~ PDurD365past_z + country:PDurD365past_z + meanQ3650past_log10_z + country:meanQ3650past_log10_z,
    family = 'poisson',
    ssn.object = ssn_dia_sedi,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_dia_sedi_mod_preds <- augment(ssn_dia_sedi_pred_final, drop=F, type.predict = 'response')
  
  if (interactive) {
    tidy(dia_sedi_rich_modlist[['mod11']])
    tidy(ssn_dia_sedi_pred_final)
    summary(ssn_dia_sedi_pred_final)
    plot(ssn_dia_sedi_pred_final)
    SSN2::varcomp(ssn_dia_sedi_pred_final)
    
    #Check predictions
    ggplot(data=ssn_dia_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    ggplot(data=ssn_dia_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() 
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_dia_sedi_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_dia_sedi_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    #,  ssn_mod_fit_best = ssn_dia_sedi_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_dia_sedi_best_preds #Augmented data for "best" model
  ))
}


#------ model_dia_biof_richness_yr ---------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

#' Diatom model for annual data
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
#' @param in_allvars_summarized A list containing merged environmental and hydrological
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
#'   
model_dia_biof_richness_yr <- function(in_ssn_eu_summarized,
                                       in_allvars_summarized,
                                       in_cor_matrices, 
                                       ssn_covtypes,
                                       scale_predictors = T,
                                       interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_dia_biof <- in_ssn_eu_summarized$dia_biof_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_dia_biof,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_dia_biof$obs$meanQ3650past_log10 <- log10(ssn_dia_biof$obs$meanQ3650past)
  # ssn_dia_biof$preds$preds_hist$meanQ3650past_log10 <- log10(ssn_dia_biof$preds$preds_hist$meanQ3650past)
  # ssn_dia_biof$preds$preds_proj$meanQ3650past_log10 <- log10(ssn_dia_biof$preds$preds_proj$meanQ3650past)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        "meanQ3650past_log10", "particle_size")
  
  ssn_dia_biof <- scale_ssn_predictors(in_ssn = ssn_dia_biof,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_richness'
  
  
  #Exploratory plots-------
  if (interactive) {
    #Check correlations
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'dia_biof_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    ggplot(ssn_dia_biof$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_dia_biof$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_dia_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_biof$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_biof$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_dia_biof$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_biof$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_dia_biof$obs, aes(x=particle_size, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    data.table::melt(ssn_dia_biof$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with ----
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('dia_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'PDurD365past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_pois_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",  #Gamma(link = "log"),
      organism = c('dia_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'PDurD365past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_nbin_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('dia_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'PDurD365past',
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
    
    ssn_cov_glance[family=='poisson' & which.min(AICc),]
    plot(ssn_norm_ini[[1]][[1]])
    plot(ssn_pois_ini[[1]][[1]])
    plot(ssn_nbin_ini[[1]][[1]])
    
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose negative binomial based on the AIC scores
    ssn_pois_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",
      organism = c('dia_biof_nopools'),
      formula_root = ' log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_pois_cov$ssn_glance)
    
    SSN2::loocv(ssn_pois_cov$ssn_list$linear_none_none)
    SSN2::loocv(ssn_pois_cov$ssn_list$epa_none_none)
    SSN2::loocv(ssn_pois_cov$ssn_list$spherical_none_none)
    lapply(list(
      ssn_pois_cov$ssn_list$linear_none_none,
      ssn_pois_cov$ssn_list$epa_none_none,
      ssn_pois_cov$ssn_list$spherical_none_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_dia_biof_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'poisson',
      ssn.object = ssn_dia_biof,
      taildown_type = 'none',
      tailup_type = 'none',
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
  dia_biof_rich_modlist <- list()
  
  #Null model
  dia_biof_rich_modlist[['null']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  dia_biof_rich_modlist[['mod1']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  dia_biof_rich_modlist[['mod2']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  dia_biof_rich_modlist[['mod3']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  dia_biof_rich_modlist[['mod4']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_biof_rich_modlist[['mod5']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_biof_rich_modlist[['mod6']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  dia_biof_rich_modlist[['mod7']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod8']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod9']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod10']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod11']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past) + particle_size')))
  
  dia_biof_rich_modlist[['mod12']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past) + particle_size + particle_size:country')))
  
  dia_biof_rich_modlist[['mod13']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod14']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + FreD_samp:country + country*log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod15']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + DurD3650past:country + country*log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod16']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_CV30yrpast + DurD_CV30yrpast:country + country*log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod17']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ meanConD_CV30yrpast + meanConD_CV30yrpast:country + country*log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod18']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + particle_size')))
  
  dia_biof_rich_modlist[['mod19']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC1')))
  
  dia_biof_rich_modlist[['mod20']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC2')))
  
  dia_biof_rich_modlist[['mod21']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_rich_modlist[['mod22']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_rich_modlist[['mod23']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_rich_modlist[['mod24']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + stream_type:country + DurD3650past + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_rich_modlist[['mod25']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + stream_type:country + DurD3650past + country*log10(meanQ3650past)')))
  
  dia_biof_rich_modlist[['mod26']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + particle_size')))
  
  dia_biof_rich_modlist[['mod27']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + particle_size + particle_size:country')))
  
  mod_perf_tab <- lapply(dia_biof_rich_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(dia_biof_rich_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- dia_biof_rich_modlist[['mod18']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    augment(chosen_mod, drop=F, type.predict = 'response') %>%
      setDT %>%
      .[!which.max(`.resid`),] %>%
      ggplot(data=., 
             aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'meanQ3650past',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'meanQ3650past',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # check_rf <- ranger::ranger(mean_richness ~ ., 
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
    # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) + 
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 23 for "absolute" best model###############"
  #####CHOOSE MOD 18 for continuous predictions###########
  
  #------ For "best" model -------------------------------------------------
  # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_dia_biof_best_final <- quick_dia_biof_ssn(
  #   as.formula(paste(alpha_var, '~ stream_type + DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC1 + env_PC2')),
  #   estmethod = 'reml'
  # )
  # ssn_dia_biof_best_preds <- augment(ssn_dia_biof_best_final,
  #                                    drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_dia_biof_best_final)
  #   SSN2::varcomp(ssn_dia_biof_best_final)
  #   plot(ssn_dia_biof_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + particle_size',
    ssn.object = ssn_dia_biof,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'poisson',
      organism = c('dia_biof_nopools'),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    # 
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    plot(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian, which=1:6)
    print('plot')
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_dia_biof_pred_final <- ssn_glm(
    formula = mean_richness ~ DurD_samp_z + DurD_samp_z:country + country*meanQ3650past_z + particle_size,
    family = 'poisson',
    ssn.object = ssn_dia_biof,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_dia_biof_mod_preds <- augment(ssn_dia_biof_pred_final, drop=F, type.predict = 'response')
  
  if (interactive) {
    tidy(dia_biof_rich_modlist[['mod18']])
    tidy(ssn_dia_biof_pred_final)
    summary(ssn_dia_biof_pred_final)
    plot(ssn_dia_biof_pred_final, which=1:6)
    print('plot')
    SSN2::varcomp(ssn_dia_biof_pred_final)
    
    #Check predictions
    ggplot(data=ssn_dia_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    ggplot(data=ssn_dia_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() 
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_dia_biof_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_dia_biof_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    # , ssn_mod_fit_best = ssn_dia_biof_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_dia_biof_best_preds #Augmented data for "best" model
  ))
}

#------ model_dia_biof_invsimpson_yr --------------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

model_dia_biof_invsimpson_yr <- function(in_ssn_eu_summarized,
                                         in_allvars_summarized,
                                         in_cor_matrices, 
                                         ssn_covtypes,
                                         scale_predictors = T,
                                         interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_dia_biof <- in_ssn_eu_summarized$dia_biof_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_dia_biof,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_dia_biof$obs$meanQ3650past_log10 <- log10(ssn_dia_biof$obs$meanQ3650past)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'meanQ3650past_log10', 'particle_size')
  
  ssn_dia_biof <- scale_ssn_predictors(in_ssn = ssn_dia_biof,
                                       in_vars = c(hydro_candidates,
                                                   'meanQ3650past_log10'),
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_invsimpson'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_dia_biof$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_dia_biof$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_dia_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_biof$obs, aes(x=basin_area_km2, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_biof$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_dia_biof$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_dia_biof$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_dia_biof$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    data.table::melt(ssn_dia_biof$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'dia_biof_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'dia_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('dia_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_gamma_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gamma',
      organism = c('dia_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_lognorm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('dia_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD3650past',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      ssn_norm_ini$ssn_glance,
      ssn_gamma_ini$ssn_glance,
      ssn_lognorm_ini$ssn_glance
    ))
    
    norm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='mean_invsimpson',][which.min(AICc),]
    lognorm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='log10(mean_invsimpson)',][which.min(AICc),]
    
    plot(ssn_lognorm_ini[[1]][[1]])
    plot(ssn_norm_ini[[1]][[1]])
    
    # plot(ssn_nbin_ini[[1]][[1]])
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose lognormal
    ssn_lognorm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('dia_biof_nopools'),
      formula_root = ' sqrt(meanQ3650past)',
      hydro_var = 'DurD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_lognorm_cov$ssn_glance)
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_dia_biof_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_dia_biof,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'spherical',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = in_random,
      estmethod = estmethod
    )
  }
  
  #Then test multiple models -----
  # This is the main model selection loop, testing various combinations of
  # predictors.
  dia_biof_simps_modlist <- list()
  
  #Initial model
  #Null model
  dia_biof_simps_modlist[['null']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  dia_biof_simps_modlist[['mod1']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  dia_biof_simps_modlist[['mod3']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  dia_biof_simps_modlist[['mod4']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_biof_simps_modlist[['mod5']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  dia_biof_simps_modlist[['mod6']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  dia_biof_simps_modlist[['mod7']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod8']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod9']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod10']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod11']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past) + particle_size')))
  
  dia_biof_simps_modlist[['mod12']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past) + particle_size + particle_size:country')))
  
  dia_biof_simps_modlist[['mod13']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod14']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ FreD_samp + FreD_samp:country + country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod15']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + DurD3650past:country + country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod16']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_CV30yrpast + DurD_CV30yrpast:country + country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod17']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ meanConD_CV30yrpast + meanConD_CV30yrpast:country + country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod18']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ STcon_m10_undirected_avg_samp + STcon_m10_undirected_avg_samp:country + country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod19']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + particle_size')))
  
  dia_biof_simps_modlist[['mod20']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC1')))
  
  dia_biof_simps_modlist[['mod21']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC2')))
  
  dia_biof_simps_modlist[['mod22']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_simps_modlist[['mod23']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_simps_modlist[['mod24']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + DurD_samp + DurD_samp:country + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_simps_modlist[['mod25']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + stream_type:country + DurD3650past + country*log10(meanQ3650past) + env_PC1 + env_PC2')))
  
  dia_biof_simps_modlist[['mod26']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ stream_type + stream_type:country + DurD3650past + country*log10(meanQ3650past)')))
  
  dia_biof_simps_modlist[['mod27']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + particle_size')))
  
  dia_biof_simps_modlist[['mod28']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + particle_size + particle_size:country')))
  
  dia_biof_simps_modlist[['mod29']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + DurD_samp:country + country*log10(basin_area_km2) + particle_size')))
  
  dia_biof_simps_modlist[['mod30']] <- quick_dia_biof_ssn(as.formula(paste(alpha_var,  '~ DurD_samp*country + particle_size + particle_size:country')))
  
  
  mod_perf_tab <- lapply(dia_biof_simps_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(dia_biof_simps_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- dia_biof_simps_modlist[['mod28']]
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    #Check residuals, variance decomposition, summary, and predictions
    plot(chosen_mod)
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    ggplot(data=augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'particle_size',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    #Model with RF to see
    # allvars_dt <- as.data.table(in_ssn$obs) %>%
    #   setorderv(c('country', 'site')) 
    # 
    # allvars_dt[
    #   , (in_vars) := lapply(.SD, function(x) base::scale(x, center=T, scale=T)),
    #   .SDcols = in_vars]
    #
    # check_rf <- ranger::ranger(mean_invsimpson ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_invsimpson', 'country'), with=F])
    # allvars_dt[, mean_simps_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_simps_preds, y=mean_invsimpson, color=stream_type)) + 
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD  for continuous predictions###########
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_invsimpson ~ DurD_samp + DurD_samp:country + particle_size + particle_size:country',
    ssn.object = ssn_dia_biof,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive)  {
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('dia_biof_nopools'),
      formula_root = 'particle_size',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_spherical)

    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_dia_biof_pred_final <- ssn_glm(
    formula = mean_invsimpson ~ DurD_samp_z + DurD_samp_z:country + 
      particle_size_z + particle_size_z:country,
    family = 'Gaussian',
    ssn.object = ssn_dia_biof,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'spherical',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_dia_biof_mod_preds <- augment(ssn_dia_biof_pred_final, drop=F)
  
  if (interactive) {
    tidy(dia_biof_simps_modlist[['mod28']])
    tidy(ssn_dia_biof_pred_final)
    summary(ssn_dia_biof_pred_final)
    plot(ssn_dia_biof_pred_final)
    SSN2::varcomp(ssn_dia_biof_pred_final)
    
    #Check predictions
    ggplot(data=ssn_dia_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_dia_biof_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_dia_biof_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
  ))
}


#------ model_fun_sedi_invsimpson_yr -------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

model_fun_sedi_invsimpson_yr <- function(in_ssn_eu_summarized,
                                         in_allvars_summarized,
                                         in_cor_matrices, 
                                         ssn_covtypes,
                                         scale_predictors = T,
                                         interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_fun_sedi <- in_ssn_eu_summarized$fun_sedi_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_fun_sedi,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_fun_sedi$obs$STcon_m10_directed_avg_samp_log10 <- log10(ssn_fun_sedi$obs$STcon_m10_directed_avg_samp)
  # ssn_fun_sedi$preds$preds_hist$STcon_m10_directed_avg_samp_log10 <- log10(ssn_fun_sedi$preds$preds_hist$STcon_m10_directed_avg_samp)
  # ssn_fun_sedi$preds$preds_proj$STcon_m10_directed_avg_samp_log10 <- log10(ssn_fun_sedi$preds$preds_proj$STcon_m10_directed_avg_samp)
  
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'STcon_m10_directed_avg_samp_log10')
  
  ssn_fun_sedi <- scale_ssn_predictors(in_ssn = ssn_fun_sedi,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_invsimpson'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_fun_sedi$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_fun_sedi$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_fun_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_fun_sedi$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_fun_sedi$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_fun_sedi$obs, aes(x= sqrt(DurD_CV30yrpast), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    data.table::melt(ssn_fun_sedi$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'fun_sedi_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('fun_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'sqrt(DurD_CV30yrpast)',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_gamma_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gamma',
      organism = c('fun_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'sqrt(DurD_CV30yrpast)',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_lognorm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('fun_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'sqrt(DurD_CV30yrpast)',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      ssn_norm_ini$ssn_glance,
      ssn_gamma_ini$ssn_glance,
      ssn_lognorm_ini$ssn_glance
    ))
    
    View(ssn_cov_glance )
    
    norm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='mean_invsimpson',][which.min(AICc),]
    lognorm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='log10(mean_invsimpson)',][which.min(AICc),]
    
    plot(ssn_lognorm_ini[[1]][[1]])
    plot(ssn_norm_ini[[1]][[1]])
    
    # plot(ssn_nbin_ini[[1]][[1]])
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose lognormal
    ssn_lognorm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('fun_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD_CV30yrpast',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_lognorm_cov$ssn_glance)
    SSN2::loocv(ssn_lognorm_cov$ssn_list$mariah_none_none)
    plot(ssn_lognorm_cov$ssn_list$mariah_none_none)
    SSN2::loocv(ssn_lognorm_cov$ssn_list$none_none_spherical)
    plot(ssn_lognorm_cov$ssn_list$none_none_spherical)
    
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_fun_sedi_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_fun_sedi,
      taildown_type = 'mariah',
      tailup_type = 'none',
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
  ssn_fun_sedi$obs$mean_invsimpson_log10 <- log10(ssn_fun_sedi$obs$mean_invsimpson)
  alpha_var <- 'mean_invsimpson_log10'
  fun_sedi_simps_modlist <- list()
  
  #Null model
  fun_sedi_simps_modlist[['null']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  fun_sedi_simps_modlist[['mod1']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_sedi_simps_modlist[['mod2']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_sedi_simps_modlist[['mod3']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  fun_sedi_simps_modlist[['mod4']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_sedi_simps_modlist[['mod5']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_sedi_simps_modlist[['mod6']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  fun_sedi_simps_modlist[['mod7']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  fun_sedi_simps_modlist[['mod8']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  fun_sedi_simps_modlist[['mod9']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  fun_sedi_simps_modlist[['mod10']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  fun_sedi_simps_modlist[['mod11']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  fun_sedi_simps_modlist[['mod12']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  fun_sedi_simps_modlist[['mod13']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(DurD_CV30yrpast)')))
  
  fun_sedi_simps_modlist[['mod14']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(meanConD_CV30yrpast)')))
  
  fun_sedi_simps_modlist[['mod15']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp')))
  
  fun_sedi_simps_modlist[['mod16']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past')))
  
  fun_sedi_simps_modlist[['mod17']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + log10(STcon_m10_directed_avg_samp)')))
  
  fun_sedi_simps_modlist[['mod18']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + sqrt(DurD_CV10yrpast) + country:sqrt(DurD_CV10yrpast)')))
  
  mod_perf_tab <- lapply(fun_sedi_simps_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(fun_sedi_simps_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- fun_sedi_simps_modlist[['mod17']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    augment(chosen_mod, drop=F, type.predict = 'response') %>%
      setDT %>%
      .[!which.max(`.resid`),] %>%
      ggplot(data=., 
             aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # check_rf <- ranger::ranger(log10(mean_invsimpson) ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_invsimpson', 'country'), with=F])
    # allvars_dt[, mean_simps_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_simps_preds, y=mean_invsimpson, color=stream_type)) +
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 17 for continuous predictions###########
  
  # #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_fun_sedi_best_final <- quick_fun_sedi_ssn(
  #   as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type +
  #     log10(basin_area_km2) + country:log10(basin_area_km2)')),
  #   estmethod = 'reml'
  # )
  # ssn_fun_sedi_best_preds <- augment(ssn_fun_sedi_best_final,
  #                               drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_fun_sedi_best_final)
  #   SSN2::varcomp(ssn_fun_sedi_best_final)
  #   plot(ssn_fun_sedi_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_invsimpson_log10 ~ country*DurD_samp + log10(STcon_m10_directed_avg_samp)',
    ssn.object = ssn_fun_sedi,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    in_ssn_eu_summarized$fun_sedi_nopools$ssn$obs$mean_invsimpson_log10 <- 
      log10(in_ssn_eu_summarized$fun_sedi_nopools$ssn$obs$mean_invsimpson)
    
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('fun_sedi_nopools'),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_epa_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_epa_none)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_fun_sedi_pred_final <- ssn_glm(
    formula = mean_invsimpson_log10 ~ country*DurD_samp_z + STcon_m10_directed_avg_samp_log10_z,
    family = 'Gaussian',
    ssn.object = ssn_fun_sedi,
    taildown_type = 'mariah',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_fun_sedi_mod_preds <- augment(ssn_fun_sedi_pred_final, drop=F)
  
  if (interactive) {
    tidy(fun_sedi_simps_modlist[['mod17']])
    tidy(ssn_fun_sedi_pred_final)
    summary(ssn_fun_sedi_pred_final)
    plot(ssn_fun_sedi_pred_final)
    SSN2::varcomp(ssn_fun_sedi_pred_final)
    
    #Check predictions
    ggplot(data=ssn_fun_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_fun_sedi_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_fun_sedi_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    # , ssn_mod_fit_best = ssn_fun_sedi_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_fun_sedi_best_preds #Augmented data for "best" model
  ))
}
#------ model_fun_sedi_richness_yr ---------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

#' Fungi model for annual data
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
#' @param in_allvars_summarized A list containing merged environmental and hydrological
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
#'   
model_fun_sedi_richness_yr <- function(in_ssn_eu_summarized,
                                       in_allvars_summarized,
                                       in_cor_matrices, 
                                       ssn_covtypes,
                                       scale_predictors = T,
                                       interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_fun_sedi <- in_ssn_eu_summarized$fun_sedi_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_fun_sedi,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_fun_sedi$obs$Fdist_mean_10past_undirected_avg_samp_log10 <- 
    log10(ssn_fun_sedi$obs$Fdist_mean_10past_undirected_avg_samp + 1)
  ssn_fun_sedi$obs$basin_area_km2_log10 <- log10(ssn_fun_sedi$obs$basin_area_km2)
  
  # ssn_fun_sedi$preds$preds_hist$Fdist_mean_10past_undirected_avg_samp_log10 <- 
  #   log10(ssn_fun_sedi$preds$preds_hist$Fdist_mean_10past_undirected_avg_samp + 1)
  # ssn_fun_sedi$preds$preds_hist$basin_area_km2_log10 <- log10(ssn_fun_sedi$preds$preds_hist$basin_area_km2)
  # 
  # ssn_fun_sedi$preds$preds_proj$Fdist_mean_10past_undirected_avg_samp_log10 <- 
  #   log10(ssn_fun_sedi$preds$preds_proj$Fdist_mean_10past_undirected_avg_samp + 1)
  # ssn_fun_sedi$preds$preds_proj$basin_area_km2_log10 <- log10(ssn_fun_sedi$preds$preds_proj$basin_area_km2)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        "Fdist_mean_10past_undirected_avg_samp_log10",
                        "basin_area_km2_log10")
  
  ssn_fun_sedi <- scale_ssn_predictors(in_ssn = ssn_fun_sedi,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_richness'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_fun_sedi$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_fun_sedi$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_fun_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_sqrt() 
    
    ggplot(ssn_fun_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_fun_sedi$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_fun_sedi$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_fun_sedi$obs, aes(x= Fdist_mean_10past_undirected_avg_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_sedi$obs, aes(x=DurD_CV30yrpast, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    
    
    data.table::melt(ssn_fun_sedi$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'fun_sedi_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('fun_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_pois_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",  #Gamma(link = "log"),
      organism = c('fun_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_nbin_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('fun_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
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
    
    ssn_cov_glance[family=='poisson' & which.min(AICc),]
    plot(ssn_norm_ini[[1]][[1]])
    plot(ssn_pois_ini[[1]][[1]])
    plot(ssn_nbin_ini[[1]][[1]])
    
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose negative binomial based on the AIC scores
    ssn_nbin_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('fun_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_nbin_cov$ssn_glance)
    
    SSN2::loocv(ssn_nbin_cov$ssn_list$mariah_none_none)
    SSN2::loocv(ssn_nbin_cov$ssn_list$exponential_none_none)
    SSN2::loocv(ssn_nbin_cov$ssn_list$spherical_none_none)
    lapply(list(
      ssn_nbin_cov$ssn_list$mariah_none_none,
      ssn_nbin_cov$ssn_list$exponential_none_none,
      ssn_nbin_cov$ssn_list$spherical_none_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_fun_sedi_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'nbinomial',
      ssn.object = ssn_fun_sedi,
      taildown_type = 'mariah',
      tailup_type = 'none',
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
  fun_sedi_rich_modlist <- list()
  
  #Null model
  fun_sedi_rich_modlist[['null']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  fun_sedi_rich_modlist[['mod1']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod2']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod3']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod4']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod5']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod6']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  fun_sedi_rich_modlist[['mod7']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  fun_sedi_rich_modlist[['mod8']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  #fun_sedi_rich_modlist[['mod9']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  fun_sedi_rich_modlist[['mod10']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  fun_sedi_rich_modlist[['mod11']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod12']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  fun_sedi_rich_modlist[['mod13']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod14']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(Fdist_mean_10past_undirected_avg_samp+1) + country:log10(Fdist_mean_10past_undirected_avg_samp+1) + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod15']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD3650past + country:DurD3650past + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod16']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ FstDrE + country:FstDrE + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod17']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(Fdist_mean_10past_undirected_avg_samp+1) + country:log10(Fdist_mean_10past_undirected_avg_samp+1) + DurD_CV10yrpast + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod18']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(Fdist_mean_10past_undirected_avg_samp+1) + country:log10(Fdist_mean_10past_undirected_avg_samp+1) + FreD_CV10yrpast + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod19']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(Fdist_mean_10past_undirected_avg_samp+1) + country:log10(Fdist_mean_10past_undirected_avg_samp+1) + meanConD_CV10yrpast + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod20']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(Fdist_mean_10past_undirected_avg_samp+1) + meanConD_CV10yrpast + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod21']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + country*log10(meanQ3650past)')))
  
  fun_sedi_rich_modlist[['mod22']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + DurD_CV10yrpast + country*log10(basin_area_km2)')))
  
  fun_sedi_rich_modlist[['mod23']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + DurD_CV10yrpast')))
  
  fun_sedi_rich_modlist[['mod24']] <- quick_fun_sedi_ssn(as.formula(paste(alpha_var,  '~ DurD_samp + country:DurD_samp + country*log10(STcon_m10_directed_avg_samp)')))
  
  
  mod_perf_tab <- lapply(fun_sedi_rich_modlist, function(x) {
    # print(x)
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(fun_sedi_rich_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- fun_sedi_rich_modlist[['mod22']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'basin_area_km2',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'basin_area_km2',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # check_rf <- ranger::ranger(mean_richness ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
    # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) +
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 22 for continuous predictions###########
  
  #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_fun_sedi_best_final <- quick_fun_sedi_ssn(
  #   as.formula(paste(alpha_var, '~ PDurD365past + country:PDurD365past + 
  #                    log10(meanQ3650past) + country:log10(meanQ3650past) +
  #                    env_PC1 + env_PC2')),
  #   estmethod = 'reml'
  # )
  # ssn_fun_sedi_best_preds <- augment(ssn_fun_sedi_best_final,
  #                                    drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_fun_sedi_best_final)
  #   SSN2::varcomp(ssn_fun_sedi_best_final)
  #   plot(ssn_fun_sedi_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ DurD_samp + country:DurD_samp + DurD_CV10yrpast + country*log10(basin_area_km2)',
    ssn.object = ssn_fun_sedi,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 20
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    # This section focuses on creating a model suitable for prediction to new reaches,
    # using only variables available for all reaches.
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'poisson',
      organism = c('fun_sedi_nopools'),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    #Test covariance structures again, this time with REML 
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'nbinomial',
      organism = c('fun_sedi_nopools'),
      formula_root = 'log10(basin_area_km2) + DurD_CV10yrpast',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none)
    # plot(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none, which=1:6)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_fun_sedi_pred_final <- ssn_glm(
    formula = mean_richness ~ DurD_samp_z + country:DurD_samp_z + 
      DurD_CV10yrpast_z + country*basin_area_km2_log10_z,
    family = 'nbinomial',
    ssn.object = ssn_fun_sedi,
    taildown_type = 'mariah',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_fun_sedi_mod_preds <- augment(ssn_fun_sedi_pred_final, drop=F, type.predict = 'response')
  
  if (interactive) {
    tidy(fun_sedi_rich_modlist[['mod22']])
    tidy(ssn_fun_sedi_pred_final)
    summary(ssn_fun_sedi_pred_final)
    plot(ssn_fun_sedi_pred_final, which=1:6)
    print('plot')
    SSN2::varcomp(ssn_fun_sedi_pred_final)
    
    #Check predictions
    ggplot(data=ssn_fun_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point() +
      # geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    ggplot(data=ssn_fun_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      # geom_smooth(method='lm') +
      geom_abline() 
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_fun_sedi_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_fun_sedi_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot
    # , #Torgegram for prediction model 
    # ssn_mod_fit_best = ssn_fun_sedi_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_fun_sedi_best_preds #Augmented data for "best" model
  ))
}



#------ model_fun_biof_invsimpson_yr -------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

model_fun_biof_invsimpson_yr <- function(in_ssn_eu_summarized,
                                         in_allvars_summarized,
                                         in_cor_matrices, 
                                         ssn_covtypes,
                                         scale_predictors = T,
                                         interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_fun_biof <- in_ssn_eu_summarized$fun_biof_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_fun_biof,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_fun_biof$obs$DurD_CV10yrpast_sqrt <- sqrt(ssn_fun_biof$obs$DurD_CV10yrpast)
  # ssn_fun_biof$preds$preds_hist$DurD_CV10yrpast_sqrt <- sqrt(ssn_fun_biof$preds$preds_hist$DurD_CV10yrpast)
  # ssn_fun_biof$preds$preds_proj$DurD_CV10yrpast_sqrt <- sqrt(ssn_fun_biof$preds$preds_proj$DurD_CV10yrpast)
  
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        "DurD_CV10yrpast_sqrt")
  
  ssn_fun_biof <- scale_ssn_predictors(in_ssn = ssn_fun_biof,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_invsimpson'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_fun_biof$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_fun_biof$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_fun_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_biof$obs, aes(x=basin_area_km2, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_biof$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_biof$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_fun_biof$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_fun_biof$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_fun_biof$obs, aes(x= sqrt(DurD_CV30yrpast), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    data.table::melt(ssn_fun_biof$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'fun_biof_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    # ssn_norm_ini <- model_ssn_hydrowindow(
    #   in_ssn = in_ssn_eu_summarized,
    #   family = 'Gaussian',
    #   organism = c('fun_biof_nopools'),
    #   formula_root = 'log10(basin_area_km2)',
    #   hydro_var = 'DurD_samp',
    #   response_var = alpha_var,
    #   ssn_covtypes = ssn_covtypes,
    #   partition_formula = as.formula('~ country'),
    #   random_formula = NULL,
    #   estmethod='ml'
    # )
    # 
    ssn_gamma_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gamma',
      organism = c('fun_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_lognorm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('fun_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      # ssn_norm_ini$ssn_glance,
      ssn_gamma_ini$ssn_glance,
      ssn_lognorm_ini$ssn_glance
    ))
    
    View(ssn_cov_glance )
    
    gamma_best <- ssn_cov_glance[family=='Gamma' & response_var=='mean_invsimpson',][which.min(AICc),]
    lognorm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='log10(mean_invsimpson)',][which.min(AICc),]
    
    plot(ssn_lognorm_ini[[1]][[1]], which=1:6)
    plot(ssn_gamma_ini[[1]][[1]], which=1:6)
    
    # plot(ssn_nbin_ini[[1]][[1]])
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose lognormal
    ssn_lognorm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('fun_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_lognorm_cov$ssn_glance)
    SSN2::loocv(ssn_lognorm_cov$ssn_list$exponential_none_none)
    plot(ssn_lognorm_cov$ssn_list$exponential_none_none)
    print('plot')
    SSN2::loocv(ssn_lognorm_cov$ssn_list$spherical_none_none)
    plot(ssn_lognorm_cov$ssn_list$spherical_none_none)
    print('plot')
    SSN2::loocv(ssn_lognorm_cov$ssn_list$none_none_gaussian)
    plot(ssn_lognorm_cov$ssn_list$none_none_gaussian)
    print('plot')
    
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_fun_biof_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_fun_biof,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'gaussian',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = in_random,
      estmethod = estmethod
    )
  }
  
  #Then test multiple models -----
  # This is the main model selection loop, testing various combinations of
  # predictors.
  ssn_fun_biof$obs$mean_invsimpson_log10 <- log10(ssn_fun_biof$obs$mean_invsimpson)
  alpha_var <- 'mean_invsimpson_log10'
  fun_biof_simps_modlist <- list()
  
  #Null model
  fun_biof_simps_modlist[['null']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  fun_biof_simps_modlist[['mod1']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_biof_simps_modlist[['mod2']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_biof_simps_modlist[['mod3']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  fun_biof_simps_modlist[['mod4']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_biof_simps_modlist[['mod5']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_biof_simps_modlist[['mod6']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  fun_biof_simps_modlist[['mod7']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  fun_biof_simps_modlist[['mod8']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  fun_biof_simps_modlist[['mod9']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  fun_biof_simps_modlist[['mod10']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  fun_biof_simps_modlist[['mod11']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  fun_biof_simps_modlist[['mod12']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  fun_biof_simps_modlist[['mod13']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(DurD_CV30yrpast)')))
  
  fun_biof_simps_modlist[['mod14']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(FreD_CV30yrpast)')))
  
  fun_biof_simps_modlist[['mod15']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp')))
  
  fun_biof_simps_modlist[['mod16']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past')))
  
  fun_biof_simps_modlist[['mod17']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + sqrt(DurD_CV10yrpast)')))
  
  fun_biof_simps_modlist[['mod18']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + sqrt(DurD_CV10yrpast) + country:sqrt(DurD_CV10yrpast)')))
  
  fun_biof_simps_modlist[['mod19']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + log10(STcon_m10_directed_avg_samp)')))
  
  fun_biof_simps_modlist[['mod20']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + log10(STcon_m10_directed_avg_samp) + country:log10(STcon_m10_directed_avg_samp)')))
  
  fun_biof_simps_modlist[['mod21']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past*FreD_samp')))
  
  fun_biof_simps_modlist[['mod22']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past + FreD_samp + country:FreD_samp')))
  
  
  mod_perf_tab <- lapply(fun_biof_simps_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(fun_biof_simps_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- fun_biof_simps_modlist[['mod17']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    augment(chosen_mod, drop=F, type.predict = 'response') %>%
      setDT %>%
      .[!which.max(`.resid`),] %>%
      ggplot(data=., 
             aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # check_rf <- ranger::ranger(log10(mean_invsimpson) ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_invsimpson', 'country'), with=F])
    # allvars_dt[, mean_simps_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_simps_preds, y=mean_invsimpson, color=stream_type)) +
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 17 for continuous predictions###########
  
  # #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_fun_biof_best_final <- quick_fun_biof_ssn(
  #   as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type +
  #     log10(basin_area_km2) + country:log10(basin_area_km2)')),
  #   estmethod = 'reml'
  # )
  # ssn_fun_biof_best_preds <- augment(ssn_fun_biof_best_final,
  #                               drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_fun_biof_best_final)
  #   SSN2::varcomp(ssn_fun_biof_best_final)
  #   plot(ssn_fun_biof_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_invsimpson_log10 ~ country*DurD_samp + sqrt(DurD_CV10yrpast)',
    ssn.object = ssn_fun_biof,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    in_ssn_eu_summarized$fun_biof_nopools$ssn$obs$mean_invsimpson_log10 <- 
      log10(in_ssn_eu_summarized$fun_biof_nopools$ssn$obs$mean_invsimpson)
    
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('fun_biof_nopools'),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = ~country,
      estmethod='reml'
    )
    View(  ssn_mod_predictions_covtypes$ssn_glance)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$spherical_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$spherical_none_none)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_fun_biof_pred_final <- ssn_glm(
    formula = mean_invsimpson_log10 ~ country*DurD_samp_z + DurD_CV10yrpast_sqrt_z,
    family = 'Gaussian',
    ssn.object = ssn_fun_biof,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'spherical',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_fun_biof_mod_preds <- augment(ssn_fun_biof_pred_final, drop=F)
  
  if (interactive) {
    tidy(fun_biof_simps_modlist[['mod17']])
    tidy(ssn_fun_biof_pred_final)
    summary(ssn_fun_biof_pred_final)
    plot(ssn_fun_biof_pred_final, which=1:6)
    SSN2::varcomp(ssn_fun_biof_pred_final)
    
    #Check predictions
    ggplot(data=ssn_fun_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_fun_biof_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_fun_biof_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    #, ssn_mod_fit_best = ssn_fun_biof_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_fun_biof_best_preds #Augmented data for "best" model
  ))
}

#------ model_fun_biof_richness_yr ---------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

#' Fungi model for annual data
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
#' @param in_allvars_summarized A list containing merged environmental and hydrological
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
#'   
model_fun_biof_richness_yr <- function(in_ssn_eu_summarized,
                                       in_allvars_summarized,
                                       in_cor_matrices, 
                                       ssn_covtypes,
                                       scale_predictors = T,
                                       interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_fun_biof <- in_ssn_eu_summarized$fun_biof_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_fun_biof,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  # Define candidate hydrological variables
  hydro_candidates <- in_allvars_summarized$cols$hydro_con_summarized
  
  # Scale 
  ssn_fun_biof <- scale_ssn_predictors(in_ssn = ssn_fun_biof,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_richness'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_fun_biof$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_fun_biof$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_fun_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_fun_biof$obs, aes(x=meanQ3650past, y=get(alpha_var),)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_fun_biof$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_fun_biof$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    ggplot(ssn_fun_biof$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_fun_biof$obs, aes(x= STcon_m10_undirected_avg_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    data.table::melt(ssn_fun_biof$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'fun_biof_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'fun_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('fun_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_pois_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",  #Gamma(link = "log"),
      organism = c('fun_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_nbin_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('fun_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
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
    
    ssn_cov_glance[family=='poisson' & which.min(AICc),]
    plot(ssn_pois_ini[[1]][[1]])
    plot(ssn_nbin_ini[[1]][[1]])
    
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose poisson based on the AIC scores and residuals
    ssn_pois_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",
      organism = c('fun_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_pois_cov$ssn_glance)
    
    SSN2::loocv(ssn_pois_cov$ssn_list$epa_mariah_none)
    SSN2::loocv(ssn_pois_cov$ssn_list$linear_spherical_none)
    SSN2::loocv(ssn_pois_cov$ssn_list$none_none_gaussian)
    lapply(list(
      ssn_pois_cov$ssn_list$epa_mariah_none,
      ssn_pois_cov$ssn_list$linear_spherical_none,
      ssn_pois_cov$ssn_list$none_none_gaussian
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_fun_biof_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'poisson',
      ssn.object = ssn_fun_biof,
      taildown_type = 'none', #Opt for simpler structure to avoid bias in residuals
      tailup_type = 'none',
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
  fun_biof_rich_modlist <- list()
  
  #Null model
  fun_biof_rich_modlist[['null']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  fun_biof_rich_modlist[['mod1']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_biof_rich_modlist[['mod2']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  fun_biof_rich_modlist[['mod3']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  fun_biof_rich_modlist[['mod4']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_biof_rich_modlist[['mod5']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_biof_rich_modlist[['mod6']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  fun_biof_rich_modlist[['mod7']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  fun_biof_rich_modlist[['mod8']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  fun_biof_rich_modlist[['mod9']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  fun_biof_rich_modlist[['mod10']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  fun_biof_rich_modlist[['mod11']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  fun_biof_rich_modlist[['mod12']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  fun_biof_rich_modlist[['mod13']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*FreD3650past')))
  
  fun_biof_rich_modlist[['mod14']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp')))
  
  fun_biof_rich_modlist[['mod15']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past')))
  
  fun_biof_rich_modlist[['mod16']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*FstDrE')))
  
  fun_biof_rich_modlist[['mod17']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(STcon_m10_undirected_avg_samp)')))
  
  fun_biof_rich_modlist[['mod18']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*meanConD_yr')))
  
  fun_biof_rich_modlist[['mod19']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp')))
  
  fun_biof_rich_modlist[['mod20']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp + sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  fun_biof_rich_modlist[['mod21']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp + FstDrE_SD10yrpast + country:FstDrE_SD10yrpast')))
  
  fun_biof_rich_modlist[['mod22']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp + log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  fun_biof_rich_modlist[['mod23']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + temperature_c')))
  
  fun_biof_rich_modlist[['mod24']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + ph')))
  
  fun_biof_rich_modlist[['mod25']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + temperature_c + country:temperature_c')))
  
  fun_biof_rich_modlist[['mod26']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + ph + country:ph')))
  
  fun_biof_rich_modlist[['mod27']] <- quick_fun_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp + sqrt(meanQ3650past)')))
  
  mod_perf_tab <- lapply(fun_biof_rich_modlist, function(x) {
    # print(x)
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(fun_biof_rich_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- fun_biof_rich_modlist[['mod14']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'FstDrE_SD10yrpast',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    
    # check_rf <- ranger::ranger(mean_richness ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
    # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) +
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 14 for continuous predictions###########
  #Adding temperature leads to a change in coefficient of the relationship with
  #DurD_samp for Spain, which indicates excessive collinearity
  
  #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_fun_biof_best_final <- quick_fun_biof_ssn(
  #   as.formula(paste(alpha_var, '~ PDurD365past + country:PDurD365past + 
  #                    log10(meanQ3650past) + country:log10(meanQ3650past) +
  #                    env_PC1 + env_PC2')),
  #   estmethod = 'reml'
  # )
  # ssn_fun_biof_best_preds <- augment(ssn_fun_biof_best_final,
  #                                    drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_fun_biof_best_final)
  #   SSN2::varcomp(ssn_fun_biof_best_final)
  #   plot(ssn_fun_biof_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ country*DurD_samp',
    ssn.object = ssn_fun_biof,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 20
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'poisson',
      organism = c('fun_biof_nopools'),
      formula_root = 'country + log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$epa_linear_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$epa_linear_none)
    plot(ssn_mod_predictions_covtypes$ssn_list$epa_linear_none, which=1:6)
    print('plot')
    
    summary(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none)
    plot(ssn_mod_predictions_covtypes$ssn_list$mariah_none_none, which=1:6)
    print('plot')
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    plot(ssn_mod_predictions_covtypes$ssn_list$none_none_none, which=1:6)
    print('plot')
  }
  
  #Re-fit final model with REML
  ssn_fun_biof_pred_final <- ssn_glm(
    formula = mean_richness ~ country*DurD_samp_z,
    family = 'poisson',
    ssn.object = ssn_fun_biof,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_fun_biof_mod_preds <- augment(ssn_fun_biof_pred_final, drop=F, type.predict = 'response')
  
  if (interactive) {
    tidy(fun_biof_rich_modlist[['mod14']])
    tidy(ssn_fun_biof_pred_final)
    summary(ssn_fun_biof_pred_final)
    plot(ssn_fun_biof_pred_final, which=1:6)
    print('plot')
    SSN2::varcomp(ssn_fun_biof_pred_final)
    
    #Check predictions
    ggplot(data=ssn_fun_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point() +
      # geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    ggplot(data=ssn_fun_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      # geom_smooth(method='lm') +
      geom_abline() 
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_fun_biof_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_fun_biof_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot
    # , #Torgegram for prediction model 
    # ssn_mod_fit_best = ssn_fun_biof_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_fun_biof_best_preds #Augmented data for "best" model
  ))
}

#------ model_bac_sedi_invsimpson_yr -------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

model_bac_sedi_invsimpson_yr <- function(in_ssn_eu_summarized,
                                         in_allvars_summarized,
                                         in_cor_matrices, 
                                         ssn_covtypes,
                                         scale_predictors = T,
                                         interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_bac_sedi <- in_ssn_eu_summarized$bac_sedi_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_bac_sedi,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  in_ssn_eu_summarized$bac_sedi_nopools$ssn$obs %<>%
    mutate(STcon_m10_directed_avg_samp_log10=log10(STcon_m10_directed_avg_samp))
  ssn_bac_sedi$obs %<>%
    mutate(STcon_m10_directed_avg_samp_log10 = log10(STcon_m10_directed_avg_samp),
           STcon_m10_undirected_avg_samp_log10 = log10(STcon_m10_undirected_avg_samp)
    )
  
  # ssn_bac_sedi$preds$preds_hist$STcon_m10_directed_avg_samp_log10 <- log10(ssn_bac_sedi$preds$preds_hist$STcon_m10_directed_avg_samp)
  # ssn_bac_sedi$preds$preds_proj$STcon_m10_directed_avg_samp_log10 <- log10(ssn_bac_sedi$preds$preds_proj$STcon_m10_directed_avg_samp)
  
  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        'STcon_m10_directed_avg_samp_log10',
                        'STcon_m10_undirected_avg_samp_log10')
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_bac_sedi <- scale_ssn_predictors(in_ssn = ssn_bac_sedi,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_invsimpson'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_bac_sedi$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_bac_sedi$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_bac_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_bac_sedi$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm')  +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_bac_sedi$obs, aes(x= log10(STcon_m10_directed_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x= sqrt(Fdist_mean_10past_undirected_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    data.table::melt(ssn_bac_sedi$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'bac_sedi_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    in_ssn_eu_summarized[['bac_sedi_nopools']]$ssn$obs$STcon_m10_directed_avg_samp_log10 <- 
      log10(in_ssn_eu_summarized[['bac_sedi_nopools']]$ssn$obs$STcon_m10_directed_avg_samp)
    
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_gamma_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gamma',
      organism = c('bac_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_lognorm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      ssn_norm_ini$ssn_glance,
      ssn_gamma_ini$ssn_glance,
      ssn_lognorm_ini$ssn_glance
    ))
    
    View(ssn_cov_glance )
    
    norm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='mean_invsimpson',][which.min(AICc),]
    lognorm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='log10(mean_invsimpson)',][which.min(AICc),]
    
    plot(ssn_lognorm_ini[[1]][[1]])
    plot(ssn_norm_ini[[1]][[1]])
    
    # plot(ssn_nbin_ini[[1]][[1]])
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose lognormal
    ssn_norm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('bac_sedi_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_norm_cov$ssn_glance)
    
    SSN2::loocv(ssn_norm_cov$ssn_list$linear_none_none)
    plot(ssn_norm_cov$ssn_list$linear_none_none)
    print('plot')
    SSN2::loocv(ssn_norm_cov$ssn_list$exponential_none_none)
    plot(ssn_norm_cov$ssn_list$exponential_none_none)
    print('plot')
    SSN2::loocv(ssn_norm_cov$ssn_list$epa_none_none)
    plot(ssn_norm_cov$ssn_list$epa_none_none)
    print('plot')
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_bac_sedi_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_bac_sedi,
      taildown_type = 'linear',
      tailup_type = 'none',
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
  bac_sedi_simps_modlist <- list()
  
  #Null model
  bac_sedi_simps_modlist[['null']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  bac_sedi_simps_modlist[['mod1']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_sedi_simps_modlist[['mod2']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_sedi_simps_modlist[['mod3']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  bac_sedi_simps_modlist[['mod4']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_sedi_simps_modlist[['mod5']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_sedi_simps_modlist[['mod6']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  bac_sedi_simps_modlist[['mod7']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  bac_sedi_simps_modlist[['mod8']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  bac_sedi_simps_modlist[['mod9']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  bac_sedi_simps_modlist[['mod10']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  bac_sedi_simps_modlist[['mod11']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  bac_sedi_simps_modlist[['mod12']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  bac_sedi_simps_modlist[['mod13']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + log10(STcon_m10_undirected_avg_samp) + country:log10(STcon_m10_undirected_avg_samp)')))
  
  bac_sedi_simps_modlist[['mod14']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + log10(STcon_m10_undirected_avg_samp)')))
  
  bac_sedi_simps_modlist[['mod15']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past) + log10(STcon_m10_undirected_avg_samp) + country:log10(STcon_m10_undirected_avg_samp)')))
  
  bac_sedi_simps_modlist[['mod16']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + FreD3650past + country:FreD3650past')))
  
  bac_sedi_simps_modlist[['mod17']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + FreD3650past')))
  
  bac_sedi_simps_modlist[['mod18']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past) + FreD3650past + country:FreD3650past')))
  
  bac_sedi_simps_modlist[['mod19']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp + country:DurD_samp')))
  
  bac_sedi_simps_modlist[['mod20']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp')))
  
  bac_sedi_simps_modlist[['mod21']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past) + DurD_samp + country:DurD_samp')))
  
  bac_sedi_simps_modlist[['mod22']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2) + sqrt(STcon_m10_undirected_avg_samp) + country:sqrt(STcon_m10_undirected_avg_samp)')))
  
  bac_sedi_simps_modlist[['mod23']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2) + FreD3650past + country:FreD3650past')))
  
  bac_sedi_simps_modlist[['mod24']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2) + DurD_samp + country:DurD_samp')))
  
  bac_sedi_simps_modlist[['mod25']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + log10(STcon_m10_undirected_avg_samp)'))) 
  
  bac_sedi_simps_modlist[['mod26']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + log10(STcon_m10_undirected_avg_samp) +   country:log10(STcon_m10_undirected_avg_samp)'))) 
  
  bac_sedi_simps_modlist[['mod27']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + log10(STcon_m10_directed_avg_samp)')))    
  
  bac_sedi_simps_modlist[['mod28']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp + log10(STcon_m10_directed_avg_samp) + country:log10(STcon_m10_directed_avg_samp)')))
  
  bac_sedi_simps_modlist[['mod29']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(STcon_m10_directed_avg_samp) + DurD_samp')))    
  
  
  mod_perf_tab <- lapply(bac_sedi_simps_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(bac_sedi_simps_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- bac_sedi_simps_modlist[['mod25']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    augment(chosen_mod, drop=F, type.predict = 'response') %>%
      setDT %>%
      .[!which.max(`.resid`),] %>%
      ggplot(data=., 
             aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # check_rf <- ranger::ranger(log10(mean_invsimpson) ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_invsimpson', 'country'), with=F])
    # allvars_dt[, mean_simps_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_simps_preds, y=mean_invsimpson, color=stream_type)) +
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 27 for continuous predictions###########
  
  # #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_bac_sedi_best_final <- quick_bac_sedi_ssn(
  #   as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type +
  #     log10(basin_area_km2) + country:log10(basin_area_km2)')),
  #   estmethod = 'reml'
  # )
  # ssn_bac_sedi_best_preds <- augment(ssn_bac_sedi_best_final,
  #                               drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_bac_sedi_best_final)
  #   SSN2::varcomp(ssn_bac_sedi_best_final)
  #   plot(ssn_bac_sedi_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_invsimpson ~ country*DurD_samp + log10(STcon_m10_undirected_avg_samp)',
    ssn.object = ssn_bac_sedi,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    in_ssn_eu_summarized$bac_sedi_nopools$ssn$obs$mean_invsimpson_log10 <- 
      log10(in_ssn_eu_summarized$bac_sedi_nopools$ssn$obs$mean_invsimpson)
    
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_sedi_nopools'),
      formula_root = 'country + log10(STcon_m10_undirected_avg_samp)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    View(ssn_mod_predictions_covtypes$ssn_glance)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_epa_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_epa_none)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_bac_sedi_pred_final <- ssn_glm(
    formula = mean_invsimpson ~ country*DurD_samp_z + STcon_m10_undirected_avg_samp_log10_z,
    family = 'Gaussian',
    ssn.object = ssn_bac_sedi,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_bac_sedi_mod_preds <- augment(ssn_bac_sedi_pred_final, drop=F)
  
  if (interactive) {
    tidy(bac_sedi_simps_modlist[['mod27']])
    tidy(ssn_bac_sedi_pred_final)
    summary(ssn_bac_sedi_pred_final)
    plot(ssn_bac_sedi_pred_final)
    SSN2::varcomp(ssn_bac_sedi_pred_final)
    
    #Check predictions
    ggplot(data=ssn_bac_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_bac_sedi_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_bac_sedi_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    # , ssn_mod_fit_best = ssn_bac_sedi_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_bac_sedi_best_preds #Augmented data for "best" model
  ))
}

#------ model_bac_sedi_richness_yr ---------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

#' Bacteria model for annual data
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
#' @param in_allvars_summarized A list containing merged environmental and hydrological
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
#'   
model_bac_sedi_richness_yr <- function(in_ssn_eu_summarized,
                                       in_allvars_summarized,
                                       in_cor_matrices, 
                                       ssn_covtypes,
                                       scale_predictors = T,
                                       interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_bac_sedi <- in_ssn_eu_summarized$bac_sedi_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_bac_sedi,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )

  # Define candidate hydrological variables
  hydro_candidates <- in_allvars_summarized$cols$hydro_con_summarized
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_bac_sedi <- scale_ssn_predictors(in_ssn = ssn_bac_sedi,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  alpha_var <- 'mean_richness'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_bac_sedi$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_bac_sedi$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_bac_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x=basin_area_km2, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_sqrt() 
    
    ggplot(ssn_bac_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_bac_sedi$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_bac_sedi$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_bac_sedi$obs, aes(x= sqrt(Fdist_mean_10past_undirected_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_sedi$obs, aes(x= log10(STcon_m10_undirected_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    data.table::melt(ssn_bac_sedi$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'bac_sedi_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_sedi_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
  }
  
  #Test all possible covariance structures
  in_ssn_eu_summarized[['bac_sedi_nopools']]$ssn$obs$STcon_m10_undirected_avg_samp_log10 <- 
    log10(in_ssn_eu_summarized[['bac_sedi_nopools']]$ssn$obs$STcon_m10_undirected_avg_samp)
  
  if (interactive) {
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('bac_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'STcon_m10_undirected_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_pois_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",  #Gamma(link = "log"),
      organism = c('bac_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'STcon_m10_undirected_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_nbin_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('bac_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'STcon_m10_undirected_avg_samp_log10',
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
    
    ssn_cov_glance[family=='Gaussian' & which.min(AICc),]
    plot(ssn_norm_ini[[1]][[1]])
    plot(ssn_pois_ini[[1]][[1]])
    plot(ssn_nbin_ini[[1]][[1]])
    
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose negative binomial based on the AIC scores
    ssn_norm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('bac_sedi_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'STcon_m10_undirected_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_norm_cov$ssn_glance)
    
    SSN2::loocv(ssn_norm_cov$ssn_list$none_none_none)
    plot(ssn_norm_cov$ssn_list$none_none_none)
    print('plot')
    SSN2::loocv(ssn_norm_cov$ssn_list$none_spherical_none)
    plot(ssn_norm_cov$ssn_list$none_spherical_none)
    print('plot')
    
    lapply(list(
      ssn_norm_cov$ssn_list$none_none_none,
      ssn_norm_cov$ssn_list$none_spherical_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_bac_sedi_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_bac_sedi,
      taildown_type = 'none',
      tailup_type = 'none',
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
  bac_sedi_rich_modlist <- list()
  
  #Null model
  bac_sedi_rich_modlist[['null']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  bac_sedi_rich_modlist[['mod1']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_sedi_rich_modlist[['mod2']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_sedi_rich_modlist[['mod3']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  bac_sedi_rich_modlist[['mod4']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_sedi_rich_modlist[['mod5']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_sedi_rich_modlist[['mod6']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  bac_sedi_rich_modlist[['mod7']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  bac_sedi_rich_modlist[['mod8']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  bac_sedi_rich_modlist[['mod9']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  bac_sedi_rich_modlist[['mod10']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  bac_sedi_rich_modlist[['mod11']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  bac_sedi_rich_modlist[['mod12']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  bac_sedi_rich_modlist[['mod13']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*FreD3650past')))
  
  bac_sedi_rich_modlist[['mod14']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*log10(STcon_m10_undirected_avg_samp)')))
  
  # bac_sedi_rich_modlist[['mod15']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*sd6_10yrpast')))
  
  bac_sedi_rich_modlist[['mod15']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past')))
  
  bac_sedi_rich_modlist[['mod16']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country*sqrt(Fdist_mean_10past_undirected_avg_samp)')))
  
  bac_sedi_rich_modlist[['mod17']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + FreD3650past')))
  
  bac_sedi_rich_modlist[['mod18']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + FreD3650past + log10(STcon_m10_undirected_avg_samp) + country:log10(STcon_m10_undirected_avg_samp)')))
  
  bac_sedi_rich_modlist[['mod19']] <- quick_bac_sedi_ssn(as.formula(paste(alpha_var,  '~ country + FreD3650past + FstDrE_SD10yrpast + country:FstDrE_SD10yrpast')))
  
  mod_perf_tab <- lapply(bac_sedi_rich_modlist, function(x) {
    # print(x)
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(bac_sedi_rich_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- bac_sedi_rich_modlist[['mod13']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'FreD3650past',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'FreD3650past',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # check_rf <- ranger::ranger(mean_richness ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
    # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) +
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 13 for continuous predictions###########
  
  #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_bac_sedi_best_final <- quick_bac_sedi_ssn(
  #   as.formula(paste(alpha_var, '~ PDurD365past + country:PDurD365past + 
  #                    log10(meanQ3650past) + country:log10(meanQ3650past) +
  #                    env_PC1 + env_PC2')),
  #   estmethod = 'reml'
  # )
  # ssn_bac_sedi_best_preds <- augment(ssn_bac_sedi_best_final,
  #                                    drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_bac_sedi_best_final)
  #   SSN2::varcomp(ssn_bac_sedi_best_final)
  #   plot(ssn_bac_sedi_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ country*DurD3650past',
    ssn.object = ssn_bac_sedi,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 20
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    # This section focuses on creating a model suitable for prediction to new reaches,
    # using only variables available for all reaches.
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_sedi_nopools'),
      formula_root = 'log10(meanQ3650past)',
      hydro_var = 'FreD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = ~ country,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    
    
    summary(ssn_mod_predictions_covtypes$ssn_list$epa_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$epa_none_none)
    # plot(ssn_mod_predictions_covtypes$ssn_list$epa_none_none, which=1:6)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  
  #Re-fit final model with REML
  ssn_bac_sedi_pred_final <- ssn_glm(
    formula = mean_richness ~ country*FreD3650past_z,
    family = 'Gaussian',
    ssn.object = ssn_bac_sedi,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_bac_sedi_mod_preds <- augment(ssn_bac_sedi_pred_final, drop=F, type.predict = 'response')
  
  if (interactive) {
    tidy(bac_sedi_rich_modlist[['mod13']])
    tidy(ssn_bac_sedi_pred_final)
    summary(ssn_bac_sedi_pred_final)
    plot(ssn_bac_sedi_pred_final, which=1:6)
    SSN2::varcomp(ssn_bac_sedi_pred_final)
    
    #Check predictions
    ggplot(data=ssn_bac_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point() +
      # geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    ggplot(data=ssn_bac_sedi_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      # geom_smooth(method='lm') +
      geom_abline() 
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_bac_sedi_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_bac_sedi_mod_preds, #Augmented data with prediction model
    torgegram_mod_pred = tg_plot
    # , #Torgegram for prediction model 
    # ssn_mod_fit_best = ssn_bac_sedi_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_bac_sedi_best_preds #Augmented data for "best" model
  ))
}

#------ model_bac_biof_invsimpson_yr -------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = F
# scale_predictors = T

model_bac_biof_invsimpson_yr <- function(in_ssn_eu_summarized,
                                         in_allvars_summarized,
                                         in_cor_matrices, 
                                         ssn_covtypes,
                                         scale_predictors = T,
                                         interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_bac_biof <- in_ssn_eu_summarized$bac_biof_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_bac_biof,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  ssn_bac_biof$obs$STcon_m10_directed_avg_samp_log10 <- log10(ssn_bac_biof$obs$STcon_m10_directed_avg_samp)
  # ssn_bac_biof$preds$preds_hist$STcon_m10_directed_avg_samp_log10 <- log10(ssn_bac_biof$preds$preds_hist$STcon_m10_directed_avg_samp)
  # ssn_bac_biof$preds$preds_proj$STcon_m10_directed_avg_samp_log10 <- log10(ssn_bac_biof$preds$preds_proj$STcon_m10_directed_avg_samp)

  # Define candidate hydrological variables
  hydro_candidates <- c(in_allvars_summarized$cols$hydro_con_summarized,
                        "STcon_m10_directed_avg_samp_log10")
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_bac_biof <- scale_ssn_predictors(in_ssn = ssn_bac_biof,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  alpha_var <- 'mean_invsimpson'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_bac_biof$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_bac_biof$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_bac_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_biof$obs, aes(x=basin_area_km2, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_biof$obs, aes(x=meanQ3650past, y=get(alpha_var))) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_biof$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_bac_biof$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_bac_biof$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_bac_biof$obs, aes(x= log10(STcon_m10_undirected_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country, scales='free_y')
    
    ggplot(ssn_bac_biof$obs, aes(x= sqrt(Fdist_mean_10past_undirected_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      facet_wrap(~country)
    
    data.table::melt(ssn_bac_biof$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'bac_biof_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
  }
  
  #Test all possible covariance structures
  in_ssn_eu_summarized[['bac_biof_nopools']]$ssn$obs$STcon_m10_directed_avg_samp_log10 <- 
    log10(in_ssn_eu_summarized[['bac_biof_nopools']]$ssn$obs$STcon_m10_directed_avg_samp)
  
  if (interactive) {
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_gamma_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gamma',
      organism = c('bac_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_lognorm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_cov_glance <- rbindlist(list(
      ssn_norm_ini$ssn_glance,
      ssn_gamma_ini$ssn_glance,
      ssn_lognorm_ini$ssn_glance
    ))
    
    View(ssn_cov_glance )
    
    norm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='mean_invsimpson',][which.min(AICc),]
    lognorm_best <- ssn_cov_glance[family=='Gaussian' & response_var=='log10(mean_invsimpson)',][which.min(AICc),]
    
    plot(ssn_lognorm_ini[[1]][[1]])
    plot(ssn_norm_ini[[1]][[1]])
    
    # plot(ssn_nbin_ini[[1]][[1]])
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose lognormal
    ssn_lognorm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('bac_biof_nopools'),
      formula_root = 'log10(basin_area_km2)',
      hydro_var = 'STcon_m10_directed_avg_samp_log10',
      response_var = paste0('log10(', alpha_var, ')'),
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_lognorm_cov$ssn_glance)
    SSN2::loocv(ssn_lognorm_cov$ssn_list$none_none_gaussian)
    plot(ssn_lognorm_cov$ssn_list$none_none_gaussian)
    print('plot')
    print('plot')
    
    SSN2::loocv(ssn_lognorm_cov$ssn_list$none_exponential_gaussian)
    plot(ssn_lognorm_cov$ssn_list$none_epa_gaussian)
    print('plot')
    print('plot')
    
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_bac_biof_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_bac_biof,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'gaussian',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = in_random,
      estmethod = estmethod
    )
  }
  
  #Then test multiple models -----
  # This is the main model selection loop, testing various combinations of
  # predictors.
  ssn_bac_biof$obs$mean_invsimpson_log10 <- log10(ssn_bac_biof$obs$mean_invsimpson)
  alpha_var <- 'mean_invsimpson_log10'
  bac_biof_simps_modlist <- list()
  
  #Null model
  bac_biof_simps_modlist[['null']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  bac_biof_simps_modlist[['mod1']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_biof_simps_modlist[['mod2']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_biof_simps_modlist[['mod3']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  bac_biof_simps_modlist[['mod4']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_biof_simps_modlist[['mod5']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_biof_simps_modlist[['mod6']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  bac_biof_simps_modlist[['mod7']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  bac_biof_simps_modlist[['mod8']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  bac_biof_simps_modlist[['mod9']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  bac_biof_simps_modlist[['mod10']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  bac_biof_simps_modlist[['mod11']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  bac_biof_simps_modlist[['mod12']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  bac_biof_simps_modlist[['mod13']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + PDurD365past')))
  
  bac_biof_simps_modlist[['mod14']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*PDurD365past')))
  
  bac_biof_simps_modlist[['mod15']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + FstDrE')))
  
  bac_biof_simps_modlist[['mod16']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*FstDrE')))
  
  bac_biof_simps_modlist[['mod17']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country +  PFreD365past')))
  
  bac_biof_simps_modlist[['mod18']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country* PFreD365past')))
  
  bac_biof_simps_modlist[['mod19']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + FreD_samp ')))
  
  bac_biof_simps_modlist[['mod20']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*FreD_samp ')))
  
  bac_biof_simps_modlist[['mod21']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp')))
  
  bac_biof_simps_modlist[['mod22']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp')))
  
  bac_biof_simps_modlist[['mod23']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + log10(STcon_m10_directed_avg_samp)')))
  
  bac_biof_simps_modlist[['mod24']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(STcon_m10_directed_avg_samp)')))
  
  bac_biof_simps_modlist[['mod25']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + PDurD365past*FstDrE')))
  
  bac_biof_simps_modlist[['mod26']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp*PFreD365past')))
  
  bac_biof_simps_modlist[['mod27']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + DurD_samp*log10(STcon_m10_directed_avg_samp)')))
  
  bac_biof_simps_modlist[['mod28']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + PDurD365past*log10(STcon_m10_directed_avg_samp)')))
  
  mod_perf_tab <- lapply(bac_biof_simps_modlist, function(x) {
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(bac_biof_simps_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- bac_biof_simps_modlist[['mod28']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    augment(chosen_mod, drop=F, type.predict = 'response') %>%
      setDT %>%
      .[!which.max(`.resid`),] %>%
      ggplot(data=., 
             aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
    
    # check_rf <- ranger::ranger(log10(mean_invsimpson) ~ .,
    #                            data= allvars_dt[, c(hydro_candidates, 'mean_invsimpson', 'country'), with=F])
    # allvars_dt[, mean_simps_preds := ranger::predictions(check_rf)]
    # ggplot(allvars_dt, aes(x=mean_simps_preds, y=mean_invsimpson, color=stream_type)) +
    #   geom_point() +
    #   facet_wrap(~country)
  }
  
  #####CHOOSE MOD 28 for continuous predictions###########
  
  # #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_bac_biof_best_final <- quick_bac_biof_ssn(
  #   as.formula(paste(alpha_var, '~ DurD3650past + stream_type + country:stream_type +
  #     log10(basin_area_km2) + country:log10(basin_area_km2)')),
  #   estmethod = 'reml'
  # )
  # ssn_bac_biof_best_preds <- augment(ssn_bac_biof_best_final,
  #                               drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_bac_biof_best_final)
  #   SSN2::varcomp(ssn_bac_biof_best_final)
  #   plot(ssn_bac_biof_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_invsimpson_log10 ~ country + PDurD365past*log10(STcon_m10_directed_avg_samp)',
    ssn.object = ssn_bac_biof,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 15
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    in_ssn_eu_summarized$bac_biof_nopools$ssn$obs$mean_invsimpson_log10 <- 
      log10(in_ssn_eu_summarized$bac_biof_nopools$ssn$obs$mean_invsimpson)
    
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_biof_nopools'),
      formula_root = 'log10(STcon_m10_directed_avg_samp)',
      hydro_var = 'PDurD365past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_mod_predictions_covtypes$ssn_glance)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_gaussian)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
  }
  
  #Re-fit final model with REML
  ssn_bac_biof_pred_final <- ssn_glm(
    formula = mean_invsimpson_log10 ~ country + PDurD365past_z*STcon_m10_directed_avg_samp_log10_z,
    family = 'Gaussian',
    ssn.object = ssn_bac_biof,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'gaussian',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_bac_biof_mod_preds <- augment(ssn_bac_biof_pred_final, drop=F)
  
  if (interactive) {
    tidy(bac_biof_simps_modlist[['mod28']])
    tidy(ssn_bac_biof_pred_final)
    summary(ssn_bac_biof_pred_final)
    plot(ssn_bac_biof_pred_final)
    SSN2::varcomp(ssn_bac_biof_pred_final)
    
    #Check predictions
    ggplot(data=ssn_bac_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_bac_biof_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_bac_biof_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot #Torgegram for prediction model 
    # , ssn_mod_fit_best = ssn_bac_biof_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_bac_biof_best_preds #Augmented data for "best" model
  ))
}

#------ model_bac_biof_richness_yr ---------------------------------------------
# in_allvars_summarized <- tar_read(allvars_summarized)
# in_ssn_eu_summarized <- tar_read(ssn_eu_summarized)
# in_cor_matrices <- tar_read(cor_matrices_list_summarized)
# tar_load(ssn_covtypes)
# tar_load(cor_heatmaps_summarized)
# interactive = T
# scale_predictors = T

#' Bacteria model for annual data
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
#' @param in_allvars_summarized A list containing merged environmental and hydrological
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
#'   
model_bac_biof_richness_yr <- function(in_ssn_eu_summarized,
                                       in_allvars_summarized,
                                       in_cor_matrices, 
                                       ssn_covtypes,
                                       scale_predictors = T,
                                       interactive = F) {
  
  #Subset SSN and create distance matrices-----
  ssn_bac_biof <- in_ssn_eu_summarized$bac_biof_nopools$ssn
  SSN2::ssn_create_distmat(
    ssn_bac_biof,
    predpts = c("preds_hist", "preds_proj"),
    among_predpts = TRUE,
    overwrite = TRUE
  )
  
  # Define candidate hydrological variables
  hydro_candidates <- in_allvars_summarized$cols$hydro_con_summarized
  
  #Scale predictor data (mean of 0 and SD of 1) -----
  ssn_bac_biof <- scale_ssn_predictors(in_ssn = ssn_bac_biof,
                                       in_vars = hydro_candidates,
                                       scale_ssn_preds = FALSE)
  
  alpha_var <- 'mean_richness'
  
  #Exploratory plots-------
  if (interactive) {
    ggplot(ssn_bac_biof$obs, aes(x=country, y=get(alpha_var))) +
      geom_boxplot()
    
    ggplot(ssn_bac_biof$obs, aes(x=country, y=get(alpha_var), fill=stream_type)) +
      geom_boxplot()
    
    ggplot(ssn_bac_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_biof$obs, aes(x=basin_area_km2, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_sqrt() 
    
    ggplot(ssn_bac_biof$obs, aes(x=meanQ3650past, y=get(alpha_var), color=stream_type)) +
      geom_point() +
      geom_smooth(method='lm') +
      scale_x_log10() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_biof$obs, aes(x=meanQ3650past, y=get(alpha_var), color=country)) +
      geom_point() +
      geom_smooth(method='lm', se=F) +
      scale_x_log10() 
    
    ggplot(ssn_bac_biof$obs, aes(x=DurD_samp, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_bac_biof$obs, aes(x=DurD3650past, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') 
    
    ggplot(ssn_bac_biof$obs, aes(x= sqrt(Fdist_mean_10past_undirected_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    ggplot(ssn_bac_biof$obs, aes(x= log10(STcon_m10_directed_avg_samp), y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm') +
      scale_x_sqrt() +
      facet_wrap(~country)
    
    data.table::melt(ssn_bac_biof$obs, 
                     id.vars=c('site', 'country', alpha_var), 
                     measure.vars=hydro_candidates
    ) %>%
      ggplot(aes(x=value, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(aes(color=country), method='lm', se=F) +
      geom_smooth(method='lm', se=F, color='black') +
      facet_wrap(~variable, scales='free_x') +
      theme_bw()
    
    #Settle on a basic covariance structure to start  with -----
    #Check correlation
    topcors_overall <- in_cor_matrices$div[
      variable1 == alpha_var & organism == 'bac_biof_nopools' 
      # & variable2 %in% hydro_candidates 
      & !is.na(correlation),] %>%
      .[, abs_cor := abs(correlation)] %>%
      setorder(-abs_cor)
    
    topcors_bydrn <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, `:=`(cor_order = frank(-abs(correlation),ties.method="first"),
               var_label = paste(variable2, round(correlation, 2))
      ), by=country] %>%
      dcast(cor_order~country, value.var = 'var_label', fill=NA)
    
    topcors_bydrn_avg <- in_cor_matrices$div_bydrn[
      variable1 == alpha_var & organism == 'bac_biof_nopools' &
        variable2 %in% hydro_candidates & !is.na(correlation),] %>%
      .[, mean(correlation), by=variable2] %>%
      .[order(V1),]
    
    #Test all possible covariance structures
    ssn_norm_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      organism = c('bac_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'FreD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_pois_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "poisson",  #Gamma(link = "log"),
      organism = c('bac_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'FreD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='ml'
    )
    
    ssn_nbin_ini <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "nbinomial",
      organism = c('bac_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'FreD3650past',
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
    
    ssn_cov_glance[family=='poisson' & which.min(AICc),]
    plot(ssn_norm_ini[[1]][[1]])
    plot(ssn_pois_ini[[1]][[1]])
    plot(ssn_nbin_ini[[1]][[1]])
    
    # varcomp(ssn_norm_ini[[1]][[10]])
    
    #CHoose negative binomial based on the AIC scores
    ssn_norm_cov <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = "Gaussian",
      organism = c('bac_biof_nopools'),
      formula_root = ' log10(basin_area_km2)',
      hydro_var = 'FreD3650past',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    View(ssn_norm_cov$ssn_glance)
    
    SSN2::loocv(ssn_norm_cov$ssn_list$none_exponential_exponential)
    plot(ssn_norm_cov$ssn_list$none_exponential_exponential, which=1:2)
    print('plot')
    SSN2::loocv(ssn_norm_cov$ssn_list$none_none_spherical)
    plot(ssn_norm_cov$ssn_list$none_none_spherical, which=1:2)
    print('plot')
    
    lapply(list(
      ssn_norm_cov$ssn_list$none_exponential_exponential,
      ssn_norm_cov$ssn_list$none_none_spherical,
      ssn_norm_cov$ssn_list$none_none_none
    ), glance) %>% rbindlist
  }
  
  #Set a quick helper function to train an SSN GLM model with a basic model structure
  quick_bac_biof_ssn <- function(in_formula, in_random=NULL, estmethod='ml') {
    ssn_glm(
      formula = in_formula,
      family = 'Gaussian',
      ssn.object = ssn_bac_biof,
      taildown_type = 'none',
      tailup_type = 'none',
      euclid_type = 'spherical',
      additive = "afv_qsqrt",
      partition_factor = ~ country,
      random = ~country,
      estmethod = estmethod
    )
  }
  
  #Then test multiple models -----
  # This is the main model selection loop, testing various combinations of
  # predictors.
  bac_biof_rich_modlist <- list()
  
  #Null model
  bac_biof_rich_modlist[['null']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ 1')))
  
  #Initial model
  bac_biof_rich_modlist[['mod1']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_biof_rich_modlist[['mod2']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2)')))
  
  bac_biof_rich_modlist[['mod3']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(basin_area_km2)')))
  
  bac_biof_rich_modlist[['mod4']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_biof_rich_modlist[['mod5']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(basin_area_km2) + country:log10(basin_area_km2)')))
  
  bac_biof_rich_modlist[['mod6']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past)'))) #Replace by qsim_3650past
  
  bac_biof_rich_modlist[['mod7']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past)')))
  
  bac_biof_rich_modlist[['mod8']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*log10(meanQ3650past)')))
  
  bac_biof_rich_modlist[['mod9']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ sqrt(meanQ3650past) + country:sqrt(meanQ3650past)')))
  
  bac_biof_rich_modlist[['mod10']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ log10(meanQ3650past) + country:log10(meanQ3650past)')))
  
  bac_biof_rich_modlist[['mod11']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + log10(basin_area_km2)')))
  
  bac_biof_rich_modlist[['mod12']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country')))
  
  bac_biof_rich_modlist[['mod13']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ RelF_avg_samp')))
  
  bac_biof_rich_modlist[['mod14']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country + sd6_30yrpast')))
  
  bac_biof_rich_modlist[['mod15']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*FreD3650past')))
  
  bac_biof_rich_modlist[['mod16']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD_samp')))
  
  bac_biof_rich_modlist[['mod17']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*FstDrE')))
  
  bac_biof_rich_modlist[['mod18']] <- quick_bac_biof_ssn(as.formula(paste(alpha_var,  '~ country*DurD3650past')))
  
  mod_perf_tab <- lapply(bac_biof_rich_modlist, function(x) {
    print(x)
    cbind(Reduce(paste, deparse(x$formula)), 
          glance(x), 
          ifelse(length(attr(x$terms, "term.labels"))>=2, max(as.data.frame(vif(x))[['GVIF^(1/(2*Df))']]), NA),
          SSN2::loocv(x))
  }) %>% 
    rbindlist %>%
    .[, mod := names(bac_biof_rich_modlist)]
  
  if (interactive) {
    #Check model in additional depth (applied to most models before continuing)
    chosen_mod <- bac_biof_rich_modlist[['mod16']]
    
    #Check residuals, variance decomposition, summary, and predictions
    par(mfrow=c(2,3))
    plot(chosen_mod, which=1:6)
    print('plot')
    varcomp(chosen_mod)
    summary(chosen_mod)
    
    #Check residual correlations to choose variable additions
    check_resid_corr(in_ssn_mod = chosen_mod,
                     in_idcol = 'site',
                     in_response_var = alpha_var,
                     in_candidates = hydro_candidates)
    
    
    ggplot(data=  augment(chosen_mod, drop=F, type.predict = 'response'), 
           aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=country)) +
      geom_smooth(method='lm') +
      geom_abline()
    
    get_ssn_emmeans(in_mod=chosen_mod , 
                    in_pred_var = 'DurD_samp',
                    interaction_var = 'country', 
                    in_drn_dt = drn_dt,
                    plot=T, label_pred_var=T)
    
    get_ssn_emtrends(in_mod=chosen_mod , 
                     in_pred_var = 'DurD_samp',
                     interaction_var = 'country', 
                     in_drn_dt = drn_dt,
                     plot=T)
  }
  
  # check_rf <- ranger::ranger(mean_richness ~ .,
  #                            data= allvars_dt[, c(hydro_candidates, 'mean_richness', 'country'), with=F])
  # allvars_dt[, mean_rich_preds := ranger::predictions(check_rf)]
  # ggplot(allvars_dt, aes(x=mean_rich_preds, y=mean_richness, color=stream_type)) +
  #   geom_point() +
  #   facet_wrap(~country)
  
  
  #####CHOOSE MOD 16 for continuous predictions###########
  
  #------ For "best" model -------------------------------------------------
  # # Refits the selected "best" model using the REML method for better variance estimation.
  # ssn_bac_biof_best_final <- quick_bac_biof_ssn(
  #   as.formula(paste(alpha_var, '~ PDurD365past + country:PDurD365past + 
  #                    log10(meanQ3650past) + country:log10(meanQ3650past) +
  #                    env_PC1 + env_PC2')),
  #   estmethod = 'reml'
  # )
  # ssn_bac_biof_best_preds <- augment(ssn_bac_biof_best_final,
  #                                    drop=F)
  # 
  # if (interactive) {
  #   summary(ssn_bac_biof_best_final)
  #   SSN2::varcomp(ssn_bac_biof_best_final)
  #   plot(ssn_bac_biof_best_final)
  # }
  
  #------ For prediction model -------------------------------------------------
  # This section focuses on creating a model suitable for prediction to new reaches,
  # using only variables available for all reaches.
  #Check Torgegram
  tg_selected <- SSN2::Torgegram(
    formula = 'mean_richness ~ country*DurD_samp',
    ssn.object = ssn_bac_biof,
    type = c("flowcon", "flowuncon", "euclid"),
    partition_factor = as.formula('~ country'),
    bins = 20
  ) %>%
    rbindlist(idcol='dist_type')
  
  tg_plot <- ggplot(tg_selected, 
                    aes(x=dist, y=gamma, color = dist_type)) +
    geom_point(aes(size=np)) +
    geom_smooth(method='lm', linewidth=2, se=F) +
    scale_size_continuous(range=c(0.5, 7)) +
    scale_x_sqrt(breaks=c(1000, 5000, 10000, 20000)) +
    theme_classic() +
    facet_wrap(~dist_type)
  
  if (interactive) {
    ssn_mod_predictions_covtypes <- model_ssn_hydrowindow(
      in_ssn = in_ssn_eu_summarized,
      family = 'Gaussian',
      organism = c('bac_biof_nopools'),
      formula_root = 'country + log10(meanQ3650past)',
      hydro_var = 'DurD_samp',
      response_var = alpha_var,
      ssn_covtypes = ssn_covtypes,
      partition_formula = as.formula('~ country'),
      random_formula = NULL,
      estmethod='reml'
    )
    
    selected_glance <- ssn_mod_predictions_covtypes$ssn_glance
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_linear_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_linear_none)
    plot(ssn_mod_predictions_covtypes$ssn_list$none_linear_none, which=1:2)
    
    summary(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    varcomp(ssn_mod_predictions_covtypes$ssn_list$none_none_none)
    plot(ssn_mod_predictions_covtypes$ssn_list$none_none_none, which=1:2)
  }
  
  #Re-fit final model with REML
  ssn_bac_biof_pred_final <- ssn_glm(
    formula = mean_richness ~ country*DurD_samp_z,
    family = 'Gaussian',
    ssn.object = ssn_bac_biof,
    taildown_type = 'none',
    tailup_type = 'none',
    euclid_type = 'none',
    additive = "afv_qsqrt",
    partition_factor = ~ country,
    random = NULL,
    estmethod = 'reml'
  )
  
  ssn_bac_biof_mod_preds <- augment(ssn_bac_biof_pred_final, drop=F, type.predict = 'response')
  
  if (interactive) {
    tidy(bac_biof_rich_modlist[['mod16']])
    tidy(ssn_bac_biof_pred_final)
    summary(ssn_bac_biof_pred_final)
    par(mfrow=c(2,3))
    plot(ssn_bac_biof_pred_final, which=1:6)
    print('plot')
    SSN2::varcomp(ssn_bac_biof_pred_final)
    
    #Check predictions
    ggplot(data=ssn_bac_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point() +
      # geom_smooth(method='lm') +
      geom_abline() +
      facet_wrap(~country)
    
    ggplot(data=ssn_bac_biof_mod_preds, aes(x=.fitted, y=get(alpha_var))) +
      geom_point(aes(color=stream_type)) +
      # geom_smooth(method='lm') +
      geom_abline() 
  }
  
  #------ Write out results ----------------------------------------------------
  return(list(
    model_selection_table = mod_perf_tab,
    ssn_mod_fit = ssn_bac_biof_pred_final, #Prediction model fit with REML
    ssn_pred_final = ssn_bac_biof_mod_preds, #AUgmented data with prediction model
    torgegram_mod_pred = tg_plot
    # , #Torgegram for prediction model 
    # ssn_mod_fit_best = ssn_bac_biof_best_final, #"Best" model fit with REML
    # ssn_pred_best = ssn_bac_biof_best_preds #Augmented data for "best" model
  ))
}

#------ get_perf_table_multiorganism-------------------------------------------
# in_mod_list = tar_read(ssn_mod_yr_fit_multiorganism)

get_perf_table_multiorganism <- function(in_mod_list) {
  
  out_tab <- lapply(in_mod_list, function(in_mod_fit) {
    if (!inherits(in_mod_fit, c('lm','glm','ssn_lm', 'ssn_glm', 'lme', 'merMod'))) {
      return(data.table(formula=NULL))
    }
    
    varcomp_cols <- SSN2::varcomp(in_mod_fit) %>%
      as.data.table %>%
      dcast(formula=.~varcomp, value.var='proportion') %>%
      .[, `.`:=NULL] %>%
      setnames(names(.), paste0('varcomp_', names(.)))
    
    #If there is a country fixed effect, 
    # re-fit with only that and get fixed-effect R2
    if (any((grepl('country\\s*([*]|[+])', 
              as.character(formula(in_mod_fit)))))) {
      
      country_formula <- as.formula(
        paste(as.character(in_mod_fit$call$formula)[2], '~ country'))
      
      country_fit <- stats::update(in_mod_fit, 
                                   formula.=country_formula, 
                                   evaluate=TRUE,
                                   ssn.object = in_mod_fit$ssn.object)
      country_r2 <- SSN2::varcomp(country_fit) %>%
        filter(varcomp=="Covariates (PR-sq)") %>%
        .[['proportion']]
      
      
      varcomp_cols[,varcomp_hydroenv_covariates := 
                     `varcomp_Covariates (PR-sq)`-country_r2]
      
    } else {
      country_r2 <- NA
      varcomp_cols[,varcomp_hydroenv_covariates := 
                     `varcomp_Covariates (PR-sq)`]
    }
    
    varcomp_cols[,`varcomp_Covariates (PR-sq)`:=NULL]
    
    augmented_dat <- as.data.table(augment(in_mod_fit, type.predict = 'response')) 
    augmented_dat[, cv_predict := loocv(in_mod_fit,
                                        cv_predict = TRUE,
                                        type = 'response')$cv_predict
                  ]
    
    # print(in_mod_fit)
    perf_table <- cbind(
      formula = format_ssn_glm_equation(in_mod_fit, greek = TRUE), 
      glance(in_mod_fit), 
      GVIF = ifelse(length(attr(in_mod_fit$terms, "term.labels"))>=2,
                    max(as.data.frame(vif(in_mod_fit))[['GVIF^(1/(2*Df))']]), 
                    NA),
      varcomp_country_only = country_r2,
      varcomp_cols,
      mape_loocv = augmented_dat[get(names(augmented_dat)[[1]])>0,
                                 Metrics::mape(get(names(augmented_dat)[[1]]), cv_predict)],
      pseudo_r2_loocv = pseudoR2(in_mod_fit, adjust=T),
      SSN2::loocv(in_mod_fit))
    return(perf_table)
  }) %>% 
    rbindlist(fill=T, idcol='organism') 
  
  return(out_tab)
}

#------ diagnose_ssn_mod -----------------------------------------------------
# in_mod_fit <- tar_read(ssn_mods_och_richness_yr)$ssn_mod_fit
# write_plots=T
# out_dir = figdir
# response_var_label = 'Mean richness'
# in_drn_dt = drn_dt
# in_hydro_vars_dt = tar_read(hydro_vars_dt)
# plot_path_prefix <- 'ssn_mod_yr_miv_diagplot'

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
plot_ssn_mod_diagplot <- function(in_mod_fit,
                                  in_drn_dt,
                                  in_hydro_vars_dt,
                                  write_plots=T,
                                  response_var_label,
                                  plot_path_prefix=NULL,
                                  out_dir) {
  response_var <- all.vars(in_mod_fit$formula)[[1]]
  
  p_obs_pred <- plot_ssn_obs_pred(in_mod_fit, 
                                  in_drn_dt,
                                  response_var_label)
  
  p_emtrends <- plot_formula_emtrends(in_mod_fit,
                                      in_drn_dt,
                                      in_hydro_vars_dt,
                                      plot = TRUE,
                                      verbose = FALSE) 
  
  out_p <- ((p_emtrends$plot + theme(legend.position = 'none') 
             | p_obs_pred)) +
    plot_layout(guides = "collect") 
  
  
  if (write_plots) {
    ggsave(
      filename = file.path(
        out_dir, 
        paste0(plot_path_prefix, '_', response_var, '.pdf')),
      plot = out_p,
      width = 8,
      height = 4,
      units='in',
      dpi=600
    )
  }
  
  return(out_p)
}

#------ predict_ssn_mod -------------------------------------------------------
# in_ssn_mod_fit = tar_read(ssn_mods_dia_sedi_richness_yr)$ssn_mod_fit
# in_ssn_mod_fit <- tar_read(ssn_mods_bac_sedi_invsimpson_yr)$ssn_mod_fit
# in_ssn_mod_fit <- tar_read(ssn_mods_fun_sedi_invsimpson_yr)$ssn_mod_fit
# in_ssn_mod_fit <- tar_read(ssn_mods_fun_sedi_richness_yr)$ssn_mod_fit

# in_ssn_mod_fit = tar_read(ssn_mods_miv_richness_yr)$ssn_mod_fit
# in_hydrocon_sites_proj = rbindlist(tar_read(hydrocon_sites_proj_gcm))
# type_predict='link'
# predict_years=c(seq(1991, 2020), seq(2041, 2100))
# verbose = TRUE
# overwrite=TRUE

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
predict_ssn_mod <- function(in_ssn_mod_fit, in_hydrocon_sites_proj, 
                            type_predict, predict_years, 
                            overwrite=FALSE, verbose=TRUE) {
  
  
  if (!(class(in_ssn_mod_fit) %in% c('ssn_lm', 'ssn_glm'))) {
    return(data.table(formula=NULL))
  } else if (class(in_ssn_mod_fit)=='ssn_lm') {
    type_predict <- NULL
  }
  
  if (verbose) {
    print('Re-assembling SSN')
  }
  
  #Get basin area and particle size
  formatted_proj_dt <- in_ssn_mod_fit$ssn.object$obs %>%
    st_drop_geometry %>%
    .[, c('site', 'rid', 'basin_area_km2', 'particle_size', 'organism')] %>%
    merge(in_hydrocon_sites_proj, ., by=c('site')) %>%
    .[year %in% predict_years,]
  
  #Format projected variables
  names(formatted_proj_dt)
  setnames(formatted_proj_dt,
           c('DurD_yr', 'FreD_yr', 'PDurD_yr', 'PFreD_yr',
             'STcon_directed_mean_yr', 
             'STcon_undirected_mean_yr', 'Fdist_undmean_yr'), 
           c('DurD_samp', 'FreD_samp', 'PDurD365past', 'PFreD365past',
             'STcon_m10_directed_avg_samp',
             'STcon_m10_undirected_avg_samp',
             'Fdist_mean_10past_undirected_avg_samp')
  )
  
  formatted_proj_dt[, `:=`(
    meanQ3650past_sqrt = sqrt(meanQ3650past),
    meanQ3650past_log10 = log10(meanQ3650past),
    STcon_m10_directed_avg_samp_log10 = log10(STcon_m10_directed_avg_samp),
    STcon_m10_undirected_avg_samp_log10 = log10(STcon_m10_undirected_avg_samp), 
    DurD_CV10yrpast_sqrt = sqrt(DurD_CV10yrpast),
    Fdist_mean_10past_undirected_avg_samp_log10 = log10(Fdist_mean_10past_undirected_avg_samp + 1),
    basin_area_km2_log10 = log10(basin_area_km2)
  )]
  
  #Get LSN preds
  lsn_path <- gsub(paste0('_', formatted_proj_dt$organism[[1]], '.ssn'), 
                   '_lsn',
                   in_ssn_mod_fit$ssn.object$path)
  
  # Create a clean prediction layer (e.g. sf object)
  preds_proj_lsn <- st_read(file.path(lsn_path, 'preds_proj.gpkg'))
  cols_to_keep <- c('country', 'rid', 'UID', 
                    setdiff(names(formatted_proj_dt), names(preds_proj_lsn)))
  
  new_preds_proj_lsn <- merge(
    preds_proj_lsn,
    formatted_proj_dt[, cols_to_keep, with=F], 
    by = c("country", "rid", "UID"))
  
  # Re-assemble an SSN object
  proj_ssn_path <- gsub(formatted_proj_dt$organism[[1]],
                  paste0(formatted_proj_dt$organism[[1]], '_', 
                         
                         '_proj_', 
                         format(Sys.Date(), '%Y%m%d')),
                  in_ssn_mod_fit$ssn.object$path)
  
  in_ssn_mod_fit$ssn.object <- SSNbler::ssn_assemble(
    edges = in_ssn_mod_fit$ssn.object$edges,
    lsn_path = gsub(paste0('_', formatted_proj_dt$organism[[1]], '.ssn'), 
                    '_lsn',
                    in_ssn_mod_fit$ssn.object$path),
    obs =   in_ssn_mod_fit$ssn.object$obs,
    preds = list(preds_proj=new_preds_proj_lsn),
    ssn_path =  proj_ssn_path,
    import = TRUE,
    check = TRUE,
    afv_col = "afv_qsqrt",
    overwrite = TRUE
  )

  #Re-compute distance matrices
  if (verbose) {
    print('Creating distance matrix')
  }
  SSN2::ssn_create_distmat(
    ssn.object = in_ssn_mod_fit$ssn.object,
    predpts = "preds_proj",
    among_predpts = FALSE,
    overwrite = TRUE
  )
  
  #Scale variables
  if (verbose) {
    print('Scaling variables')
  }
  in_ssn_mod_fit$ssn.object <- scale_ssn_predictors(
    in_ssn = in_ssn_mod_fit$ssn.object, 
    in_vars = setdiff(names(formatted_proj_dt), 
                      c('reach_id', 'year', 'UID', 'date', 'site', 
                        'country', 'gcm', 'scenario', 'organism', 'rid')), 
    scale_ssn_preds = TRUE)
  
  #Generate predictions means and CI
  if (verbose) {
    print('Producing predictions')
  }
  preds_proj_dt <- SSN2::augment(in_ssn_mod_fit,  
                                 , newdata = 'preds_proj'
                                 , drop = FALSE
                                 , type.predict = type_predict
                                 , interval = 'prediction'
                                 , se_fit = TRUE
  ) %>% 
    .[, c('rid','organism', 'site', 'country', 'scenario', 'gcm', 'year',
          '.fitted', '.lower', '.upper', '.se.fit')] %>%
    st_drop_geometry %>%
    as.data.table 
  
  preds_proj_dt[, response_var := all.vars(in_ssn_mod_fit$formula)[[1]]]
  
  # Check that values make sense and are comparable between those used for 
  # training models and those from GCMs
  # check <- merge(in_hydrocon_sites_proj[year==2021,],
  #                st_drop_geometry(ssn_obj$obs[, c('site', 'country',
  #                                                 'STcon_m10_directed_avg_samp',
  #                                                 'Fdist_mean_10past_undirected_avg_samp',
  #                                                 'meanQ3650past')]),
  #                by=c('site', 'country'))
  # 
  # ggplot(check, aes(x=STcon_mean_yr, y=STcon_m10_directed_avg_samp, color=interaction(gcm, scenario))) +
  #   geom_point()+
  #   geom_abline() +
  #   scale_x_sqrt() +
  #   scale_y_sqrt()
  # 
  # ggplot(check[country='Hungary',], 
  #        aes(x=Fdist_undmean_yr, y=Fdist_mean_10past_undirected_avg_samp, 
  #            color=interaction(gcm, scenario))) +
  #   geom_point()+
  #   geom_abline() +
  #   scale_x_sqrt() +
  #   scale_y_sqrt()
  # 
  # ggplot(check, aes(x=meanQ3650past.x, y=meanQ3650past.y, color=interaction(gcm, scenario))) +
  #   geom_point()+
  #   geom_abline() +
  #   scale_x_log10() +
  #   scale_y_log10()
  
  #Write out results  ---------------------------------------------------------
  return(preds_proj_dt)
}


#------ compute_future_change --------------------------------------------------
# tar_load(ssn_mod_yr_fit_multiorganism)
# in_mod_fit <- names(ssn_mod_yr_fit_multiorganism)[[1]]
# in_ssn_mod_fit = ssn_mod_yr_fit_multiorganism[[in_mod_fit]]
# in_ssn_proj_dt = rbindlist(ssn_proj_dt, fill=T)[mod==in_mod_fit,]
# reference_years = seq(1991, 2020)
# future_years = list(mid_century=seq(2041, 2070), 
#                     late_century=seq(2071, 2100))
# n_sim = 100

compute_div_change <- function(in_ssn_mod_fit, in_ssn_proj_dt,
                               reference_years, future_years, n_sim) {
  
  if (!(class(in_ssn_mod_fit) %in% c('ssn_lm', 'ssn_glm'))) {
    return(data.table(formula=NULL))
  } else if (class(in_ssn_mod_fit)=='ssn_glm') {
    # Back transform from link space
    inv_link <- get_inverse_link_function(in_ssn_mod_fit$family)   # change if your model uses a different link
  } else {
    inv_link <- identity
  }
  
  #Simulate data and assign to corresponding period
  years_dt <- as.data.table(future_years) %>% 
    melt(variable.name="period", value.name="year") %>%
    rbind(data.table(year=reference_years, period='reference'))
  
  sim_dt <- in_ssn_proj_dt[,  
                    list(
                      n = seq(n_sim),
                      sim = inv_link(
                        rnorm(n = n_sim, mean = .fitted, sd = .se.fit))),
                    by=.(organism, site, country, scenario,
                         gcm, year, response_var, mod)] %>%
    merge(years_dt, by='year')
  
  #Compute change statistics
  sim_change_dt <- sim_dt[, list(mean_div = mean(sim)),
                          by=.(organism, site, country, scenario, 
                               gcm, period, n, response_var, mod)] %>%
    dcast(organism+site+country+scenario+gcm+n+response_var+mod~period, 
          value.var='mean_div') %>%
    .[, `:=`(div_change_mid = 100*(mid_century - reference)/reference,
             div_change_late = 100*(late_century - reference)/reference
    )]
  
  
  stats_change_dt <- melt(
    sim_change_dt,
    id.vars=c('organism', 'site', 'country',
              'scenario', 'gcm', 'n', 'response_var', 'mod'),
    measure.vars = c('div_change_mid', 'div_change_late')) %>%
    .[, list(mean_change = mean(value),
             lower_change = quantile(value, probs=0.025),
             upper_change = quantile(value, probs=0.975)
    ),
    by=.(organism, site, country, scenario, gcm, variable, response_var, mod)]
  
  
  return(list(
    sims_dt = sim_change_dt,
    stats_dt = stats_change_dt)
  )
}

#------ plot_ssn_proj------------------------------------------------------------
# in_future_sims_dt <- tar_read(future_change_dt) %>%
#   lapply(function(x) x[['sims_dt']]) %>%
#   rbindlist(fill=T)
# in_future_stats_dt <- tar_read(future_change_dt) %>%
#   lapply(function(x) x[['stats_dt']]) %>%
#   rbindlist(fill=T)
# in_drn_dt = drn_dt
# in_env_summarized <- tar_read(env_summarized)
# in_mod_fit_list <- tar_read(ssn_mod_yr_fit_multiorganism)
# in_organism_dt <- tar_read(organism_dt)
# mod_sub = c('fun_sedi_richness','dia_sedi_richness', 'miv_richness', 'ept_richness', 'och_richness')
# write_plots=T,
# out_dir = figdir

#' Map SSN model projections and changes
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
plot_ssn_proj <- function(
    in_future_sims_dt,
    in_future_stats_dt,
    in_drn_dt,
    in_env_summarized,
    in_mod_fit_list,
    in_organism_dt,
    mod_sub = c('miv_richness', 'ept_richness', 'och_richness', 
                'dia_sedi_richness', 'fun_sedi_richness'),
    write_plots=T,
    out_dir=NULL
    ) {
  #Only projects models for which hydro-env variables
  # explain at least 30% of variance*
  
  build_path <- function(prefix, ext='.pdf') {
    file.path(out_dir,
              paste0(prefix, '_', format(Sys.Date(), '%Y%m%d'), ext)
    )
  }
  
  #------------------ Make plots -----------------------------------------------
  #For a given organism, make a plot of scenario x country x period + fill=GCM
  sub_dt <- in_future_stats_dt[response_var=='mean_richness' &
                                 organism=='miv_nopools' &
                                 scenario %in% c('ssp126', 'ssp585'),]
  
  plot_change_miv_boxplot <- ggplot(sub_dt, 
                            aes(x=variable, y=mean_change, fill=gcm)) +
    geom_hline(yintercept=0, color='darkgrey') + 
    geom_boxplot(outliers=TRUE, alpha=0.8) +
    # geom_boxplot(aes(y=lower_change), outliers=FALSE, alpha=0.5) +
    # geom_boxplot(aes(y=upper_change), outliers=FALSE, alpha=0.5) +
    facet_grid(country~scenario) + 
    scale_x_discrete(labels=c('2041-2070', '2071-2100')) +
    scale_y_continuous(breaks=seq(-60, 20, 20)) +
    labs(y='Predicted change in species richness (%)') +
    theme_bw()
  
  plot_change_miv_ridges <- in_future_stats_dt[response_var=='mean_richness' &
                                          organism=='miv_nopools' &
                                          scenario %in% c('ssp126', 'ssp585'),] %>%
    ggplot(aes(y=gcm, fill=gcm)) +
    # ggridges::geom_density_ridges(aes(x=lower_change), color=NA,
    #                               alpha=0.3, stat="binline", bins=10) +
    # ggridges::geom_density_ridges(aes(x=upper_change), color=NA,
    #                               alpha=0.3, stat="binline", bins=10) +
    ggridges::geom_density_ridges(aes(x=mean_change), 
                                  color=NA, alpha=0.6,
                                  stat="binline", bins=10, panel_scaling = FALSE) +
    coord_flip() +
    facet_nested(country~scenario+variable) +
    ggridges::theme_ridges() +
    theme(axis.text.x = element_blank())
  
  
  #Make a cross-organism plot of distribution of mean change across GCMs  -------
  sub_mean_dt <- in_future_stats_dt[
    response_var=='mean_richness' & 
      scenario %in% c('ssp126', 'ssp585') &
      mod %in% mod_sub, 
    list(
      mean_change = mean(mean_change),
      SNR_change = mean(mean_change)/abs(max(mean_change)-min(mean_change)),
      gcm = 'Mean'), 
    by=.(site, variable, scenario, country, organism, mod)] %>%
    merge(in_env_summarized$dt_nopools[!duplicated(site), 
                                       .(stream_type, site)],
          by='site') %>%
    merge(in_drn_dt, by='country')  %>%
    .[, country := factor(
      country,
      levels = c("Finland", "France",  "Hungary", "Czechia", "Croatia", "Spain" ),
      ordered=T)
    ] %>%
    merge(in_organism_dt, by='organism') %>%
    .[, mod := factor(mod, levels=mod_sub)] %>%
    .[order(mod),] %>%
    .[, organism_label := factor(organism_label, levels=unique(organism_label))]

  color_vec <- sub_mean_dt[!duplicated(country),
                          setNames(color, country)]
  
  plot_list <- list() 
  
  plot_list[['plot_multigcm_avg_crossorganism']] <- ggplot(
    sub_mean_dt, aes(x=variable, y=mean_change, color=country, fill=country)
  ) +
    geom_hline(yintercept=0, color='darkgrey') + 
    geom_boxplot(alpha=0.5, outliers=TRUE) +
    scale_x_discrete(labels=c('2041-2070', '2071-2100')) +
    # scale_y_continuous(breaks=seq(-60, 20, 20)) +
    scale_fill_manual(values=color_vec) +
    scale_color_manual(values=color_vec) +
    labs(y='Predicted change in species richness (%)') +
    facet_grid(organism_label~scenario, scales='free_y') + 
    theme_bw()
  
  plot_list[['plot_multigcm_avg_crossorganism_type']] <- ggplot(
    sub_mean_dt,
    aes(x=variable, y=mean_change, fill=country, color=country, linetype=stream_type)
  ) +
    geom_hline(yintercept=0, color='darkgrey') + 
    geom_boxplot(alpha=0.5, outliers=TRUE) +
    scale_x_discrete(labels=c('2041-2070', '2071-2100')) +
    scale_linetype_discrete(name='Stream type',
                            labels = c('Perennial', 'Non-perennial')) +
    # scale_y_continuous(breaks=seq(-60, 20, 20)) +
    scale_fill_manual(values=color_vec) +
    scale_color_manual(values=color_vec) +
    labs(y='Predicted change in species richness (%)') +
    facet_grid(organism_label~scenario, scales='free_y') + 
    theme_bw()
  
  plot_list[['plot_SNR_crossorganism']] <- ggplot(
    sub_mean_dt, aes(x=variable, y=SNR_change, color=country, fill=country)
  ) +
    geom_hline(yintercept=0, color='darkgrey') + 
    geom_boxplot(alpha=0.5, outliers=FALSE) +
    scale_x_discrete(labels=c('2041-2070', '2071-2100')) +
    # scale_y_continuous(breaks=seq(-60, 20, 20)) +
    scale_fill_manual(values=color_vec) +
    scale_color_manual(values=color_vec) +
    labs(y='Mean/Range of % change in species richness across GCMs') +
    facet_grid(organism_label~scenario, scales='free_y') + 
    theme_bw()
  
  plot_list[['plot_SNR_crossorganism_type']] <- ggplot(
    sub_mean_dt,
    aes(x=variable, y=SNR_change, fill=country, color=country, linetype=stream_type)
  ) +
    geom_hline(yintercept=0, color='darkgrey') + 
    geom_boxplot(alpha=0.5, outliers=FALSE) +
    scale_x_discrete(labels=c('2041-2070', '2071-2100')) +
    # scale_y_continuous(breaks=seq(-60, 20, 20)) +
    scale_fill_manual(values=color_vec) +
    scale_color_manual(values=color_vec) +
    labs(y='Mean/Range of % change in species richness across GCMs') +
    facet_grid(organism_label~scenario, scales='free_y') + 
    theme_bw()
  
  #write to PDF
  if (write_plots) {
    lapply(names(plot_list), function(p_name) {
      ggsave(filename = build_path(p_name, ext='.png'),
             plot = plot_list[[p_name]],
             width = 8, height = 9, units='in', dpi=600)
    })
  }

  #------------------ Make maps -----------------------------------------------
  # Make map of predicted 1991-2020 richness -----------------------------------
  ref_multigcm_dt <- in_future_sims_dt[
    , list(ref_richness_pred_mean = mean(reference)),
    by=.(organism, site, country, scenario, gcm, mod)] %>%
    .[, list(ref_richness_pred_multigcm_mean = mean(ref_richness_pred_mean),
             ref_richness_pred_multigcm_range = (max(ref_richness_pred_mean)
                                                 -min(ref_richness_pred_mean))),
      by=.(organism, site, country, scenario, mod)]
  
  
  # mod_name <-  mod_sub[[3]]
  maps_list <- lapply(mod_sub, function(mod_name) {
    ssn_preds_hist <-in_mod_fit_list[[mod_name]]$ssn.object 
    
    ssn_preds_hist$obs <- ssn_preds_hist$obs %>%
      merge(ref_multigcm_dt[mod==mod_name,], 
            by=c('site', 'country', 'organism')) %>%
      mutate(perc_diff = 100*(ref_richness_pred_multigcm_mean-mean_richness)/mean_richness) %>%
      merge(in_organism_dt, by='organism')
    
    map_ref_pred <- map_ssn_facets(in_ssn=ssn_preds_hist, 
                                   in_pts='obs',
                                   facet_col='country',
                                   ptcolor_col = 'ref_richness_pred_multigcm_mean',
                                   shape_col = 'stream_type',
                                   linewidth_col='qsim_avg',       
                                   nbreaks=8,
                                   page_title=str_wrap(
                                     paste0('Predicted species richness - ',
                                                     unique(ssn_preds_hist$obs$organism_label)[1],
                                                     ' (1991-2020)'),
                                     60)
    )
    
    map_ref_diff <- map_ssn_facets(in_ssn=ssn_preds_hist, 
                                   in_pts='obs',
                                   facet_col='country',
                                   ptcolor_col = 'perc_diff',
                                   shape_col = 'stream_type',
                                   linewidth_col='qsim_avg', 
                                   nbreaks=10,
                                   page_title=str_wrap(
                                   paste0('Species richness (% difference 
                                                   predicted mean 1991-2020 vs observed 2021) - ',
                                                     unique(ssn_preds_hist$obs$organism_label)[1]),
                                   60)
    ) 
    
    ggplot(ssn_preds_hist$obs, aes(x=country, y=perc_diff)) +
      geom_boxplot(outliers=F)
    
    
    # Make map of future changes --------------------------------------------------
    ssn_preds_hist <-in_mod_fit_list[[mod_name]]$ssn.object 
    
    #For a specific projection scenario (2070-2099, ssp585)
    change_topmap <- sub_mean_dt[mod==mod_name &
                                   variable=='div_change_late' &
                                   scenario=='ssp585',]
    
    ssn_preds_hist$obs <- ssn_preds_hist$obs %>%
      merge(change_topmap, 
            by=c('site', 'country', 'organism', 'stream_type')) %>%
      merge(in_organism_dt, by='organism')
    
    map_change_pred <- map_ssn_facets(in_ssn=ssn_preds_hist, 
                                      in_pts='obs',
                                      facet_col='country',
                                      ptcolor_col = 'mean_change',
                                      shape_col = 'stream_type',
                                      linewidth_col='qsim_avg',       
                                      nbreaks=8,
                                      page_title=str_wrap(
                                        paste0('Predicted % change in species richness - ',
                                               unique(ssn_preds_hist$obs$organism_label)[1],
                                               ' (2071-2100 vs 1991-2020) - SSP5-8.5'),
                                        60)
                                      )
  
    out_map_list <- list(
      map_ref_pred = map_ref_pred,
      map_ref_diff = map_ref_diff,
      map_change_pred = map_change_pred
    )
    
    #write to PDF
    if (write_plots) {
      lapply(names(out_map_list), function(p_name) {
        ggsave(filename = build_path(paste0(mod_name, '_', p_name)),
               plot = out_map_list[[p_name]],
               width = 8, height = 8, units='in', dpi=600)
      })
    }
    
    return(out_map_list)
    
  }) %>% setNames(mod_sub)

  
  return(list(
    plots = plot_list,
    maps = maps_list
  ))
}



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
                                    "Czechia" = "#f78c6b", 
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


#------ plot_ssn2_marginal_effects ---------------------------------------------
# tar_load(ssn_mods_miv_richness_yr)
# write_plots=T
# out_dir = figdir
# in_mod <- ssn_mods_miv_richness_yr$ssn_mod_fit

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
  
  emm_best_stream_type <- emmeans(in_mod, 
                                  data=in_mod$ssn.object$obs, 
                                  specs=~FreD3650past_z|country,
                                  type = "response") %>%
    data.frame %>%
    ggplot(aes(x=country, y=prob, color=stream_type)) +
    geom_point(size=2) +
    geom_segment(aes(y=asymp.LCL, yend=asymp.UCL), linewidth=2, alpha=0.5) +
    scale_color_manual(
      name = 'Stream type',
      labels = c('Perennial', 'Non-perennial'),
      values=c('#2b8cbe', '#feb24c')) +
    scale_y_continuous(name='Estimated marginal mean: mean invsimpson') +
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
                            ~ log10(basin_area_km2) | country, 
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
                               ~ log10(basin_area_km2) | country, 
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
