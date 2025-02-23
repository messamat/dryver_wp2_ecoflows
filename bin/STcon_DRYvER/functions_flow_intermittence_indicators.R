#########################################################################################
## Import data
#########################################################################################
# watershed
import_watershed <- function(catchment)
{
  watershed_shp_PATH <- paste0("data/",catchment,"/watershed_small_catchment.shp")
  watershed <- st_read(watershed_shp_PATH) 
  
  watershed
}

# reaches
import_reaches <- function(catchment, watershed)
{
  reach_shp_PATH <- paste0("data/",catchment,"/river_network.shp")
  reaches <- st_read(reach_shp_PATH) 
  
  # intersection of reaches and watershed
  reaches$cnt = stringr::str_count(reaches$geometry, ",")
  reaches = filter(reaches, cnt>1)
  reaches <- st_intersection(reaches, watershed)
  
  reaches
}

get_results_file_present <- function(catchment)
{
  directory <- paste0("data/",catchment,"/Results_present_period/")
  files <- list.files(directory)
  files <- files[!(files == "readme.txt")]
  if (length(files) == 1)
  {
    paste0(directory,files)
  }
  else if (length(files) == 0)
  {
    print(paste0("Error : no file(s) in ","data/",catchment,"/Results_present_period/"))
  }else{
    day <- strsplit(files,"_")
    day <- unlist(lapply(day, function(x) substr(x[2],1,10)))
    day <- as.Date(day, format="%Y-%m-%d")
    k <- which(day == max(day))
    paste0(directory,files[k])
  }
}

get_nc_var_present <- function(nc, varname, reachID, dates)
{
  data <- ncvar_get(nc, varname)
  data <- as.data.frame(data)
  colnames(data) <- reachID
  data$dates <- dates
  
  hy <- hsaHydroYearSeasons(data$dates, month = 10, day_of_the_year = TRUE, # compute hydrological years, and other time periods (seasons, month, civil year, day of the year)
                              seasons = list("SummerFall" = 5:10, "WinterSpring" = c(11:12, 1:4)),
                              seasons_by_day = FALSE, minimal = FALSE)
  
  data$doy <- hy$j # month
  data$month <- hy$m # month
  data$season <- hy$s # season
  data$year <- hy$y # year
  data$yearHy <- hy$hy # hydrological year
  
  data <- melt(data, id.vars=c("dates","doy","month","season","year","yearHy"))
  colnames(data)[colnames(data) == "variable"] <- "reachID"

  data
}

# -------------------------------------------------------------------
# Compute a data.frame containing the hydrological year, and 
# optionnally the days of the hydrological year and the season which
# each day belongs to. (developped by Ivan Horner)
# -------------------------------------------------------------------
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

#########################################################################################
## Flow intermittence indicators
#########################################################################################
# Indicator 1 : Number of days with dry/flowing conditions per timescale
# flow_intermittence : dataframe with flow intermittence data
# condition : "dry" or "flowing"
# time_scale : "yearHy" (hydrological year) or "season" or "month"
# reaches : stream network (sf object)
indicator_1 <- function(flow_intermittence, condition, time_scale, reaches)
{
  condition_abr <- "F"
  if(condition == "dry")
  {
    flow_intermittence$value <- ifelse(flow_intermittence$value == 0,1,0)
    condition_abr <- "D"
  }
  
  k <- which(colnames(flow_intermittence) == time_scale)
  flow_intermittence$time_scale <- flow_intermittence[,k]
  df <- flow_intermittence %>%
    group_by(time_scale, yearHy, reachID) %>%
    summarise(nb_days=sum(value)) %>%
    ungroup()
  
  if(time_scale == "yearHy")
  {
    df <- df %>%
      group_by(reachID) %>%
      summarise(mean_nb_days=mean(nb_days, na.rm=FALSE))
    
    name <- paste0("Con",condition_abr)
    reaches$indicator <- sapply(reaches$cat, function(x){if(x %in% df$reachID){df$mean_nb_days[df$reachID == x]}else{NA}})
    colnames(reaches)[colnames(reaches) == "indicator"] <- name
    
  }else{
    df <- df %>%
      group_by(time_scale, reachID) %>%
      summarise(mean_nb_days=mean(nb_days, na.rm=FALSE))
    df <- df %>% dcast(reachID ~ time_scale, value.var = "mean_nb_days")
    #df$reachID <- as.numeric(as.character(df$reachID))
    
    names <- paste0("Con",condition_abr,"_",colnames(df)[-1])
    for (i in 1:length(names))
    {
      reaches <-cbind(reaches,sapply(reaches$cat, function(x){if(x %in% df$reachID){df[df$reachID == x, i+1]}else{NA}}))
      colnames(reaches)[grepl("sapply",colnames(reaches))] <- names[i]
    }
  }
  
  reaches
}




