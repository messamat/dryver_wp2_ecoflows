#' Calculate Spatiotemporal Connectivity Index
#'
#' This function calculates spatiotemporal connectivity indices based on 
#' spatiotemporal graph networks. The aim is to provide a methodological framework 
#' for calculating these indices based on high-frequency information obtained 
#' from natural ecosystems.
#'
#' The system does not necessarily need to quantify water absence/presence 
#' but each cell value must represent a feature defining connectivity "on" or "off" 
#' and that can be transmitted to build ecologically meaningful links in a
#' spatiotemporal graph. 
#'
#' @param sites_status_matrix A matrix representing the status of the site (wet/dry, 
#' active/inactive). 
#' The dataset should have columns for each site and rows for each monitored day.
#' Warning: no other columns should be included (e.g., date, other IDs)
#' @param network_structure A square matrix representing the basic connections among 
#' sites (adjacency matrix): for a given site (row), each adjacent connected 
#' site (column) is given a value of 1, all others 0. Must have the same number
#' of rows and columns as there are columns in sites_status_matrix.
#' @param direction The direction of the graph, either "directed" or "undirected".
#' @param routing_mode The direction for graph connectivity when directed, 
#' can be "in" (routing from upstream if directed), "out" (routing from 
#' downstream if directed), or "all". See ?igraph or ?igraph::closeness
# for a better understanding. 
#' @param weighting Logical; whether to weight the connectivity based on distances.
#' @param dist_matrix A distance matrix representing the distances between sites.
#' Can be any type of distance (euclidean, environmental, topographic, ...) 
#' between pairs of sites
#' @param weighting_links Logical; whether to use specific data for each time 
#' unit to quantify connections.
#' @param link_weights A list of data frames with the same structure as 
#' `sites_status_matrix`, providing daily flow data.
#' @param indirect_dispersal Logical; whether to consider indirect dispersal 
#' for lost nodes (i.e., a bump in dispersal before disconnection).
#' @param standardize_neighbors Logical; whether to standardize STCon based 
#' on the number of possible connections (i.e., number of upstream neighbors
#' for downstream nodes when routing_mode==in, as downstream nodes will
#' necessarily have higher values everything being equal).
#' @param value_s_link Numeric; value attributed to each spatial link.
#' @param value_t_link Numeric; value attributed to each temporal link.
#' @param value_no_s_link Numeric; value attributed to spatial disconnections.
#' @param value_no_t_link Numeric; value attributed to temporal disconnections.
#' @param legacy_effect Numeric; weight given to the link/no-link values for each time unit.
#' It is a vector of length legacy_length with values ranging from 0 to 1 to
#' modulate the relevance of a connection through time.
#' @param legacy_length Numeric; number of time units considered for legacy effects.
#' Set to 1 by default (i.e., the connectivity of time T1 will only be considered for T2).
#' If value is set at legacy_length=3, the connectivity of T1 will be considered 
#' for T2, T3, and T4.
#' @param convert_to_integer Logical; whether to convert the output to integers 
#' to minimize storage of outputs.
#' @param rounding_factor Numeric; factor used for rounding the output before 
#' converting to integer (can be a fraction of 1 when distances are large so that
#' STCon < 2*10^9 to be converted).
#' @param verbose Logical; whether to print messages during the function execution.
#' @param output A character vector specifying the output components: 
#' 'Main_matrix', 'STconmat', 'STcon'.
#'
#' @return A list containing the spatiotemporal connectivity matrix, 
#' the collapsing matrix (STconmat), and the standardized connectivity (STcon).
#' @export
#'
#' @references
#' The current script is a modification of the original function introduced in:
#' Cunillera-Montcusi, D., Fernandez-Calero, J. M., Polsterl, S., Valera, J., 
#' Argelich, R., Cid, N., Bonada, N., & Canedo-Arguelles, M. (2023). 
#' Navigating through space and time: a methodological approach to quantify 
#' spatiotemporal connectivity using flow intermittence data as a study case.
#'
#' @authors 
#' David Cunillera-Montcusi (\email{david.cunillera@dcm.cat}), original developer
#' Mathis Messager (\email{mathis.messager@mail.mcgill.ca}), code optimization
#'
#' @examples
#' \dontrun{
#' site_status_mat <- data.frame(
#'  MonitoredDays = c("Day1", "Day2", "Day3", "Day4"),
#'  StreamSite1 = c(0, 1, 1, 1),
#'  StreamSite2 = c(0, 1, 1, 1),
#'  StreamSite3 = c(0, 1, 1, 1),
#'  StreamSite4 = c(0, 0, 0, 1),
#'  StreamSite5 = c(0, 1, 0, 1)
#' )
#' network_structure <- matrix(1, ncol = 5, nrow = 5)
#' compute_stcon(site_status_mat[,2:ncol(site_status_mat)],
#'              network_structure,
#'              direction = "directed")

compute_stcon <- function(sites_status_matrix, 
                          network_structure,
                          direction,
                          routing_mode = "in",
                          weighting = FALSE,
                          dist_matrix,
                          weighting_links = FALSE,
                          link_weights,
                          indirect_dispersal = TRUE,
                          standardize_neighbors = FALSE,
                          value_s_link = 1L,
                          value_t_link = 1L,
                          value_no_s_link = 0L,
                          value_no_t_link = 0L,
                          legacy_effect = 1L,
                          legacy_length = 1L,
                          convert_to_integer = TRUE, 
                          rounding_factor = 100L,
                          verbose = T,
                          output = c('Main_matrix', 
                                     'STconmat',
                                     'STcon')) {
  
  # Validate `direction`
  assert_that(direction %in% c("directed", "undirected"), 
              msg = "ERROR: 'direction' must be either 'directed' or 'undirected'.")
  
  if (verbose) {
    message("Your river will be considered as a ", 
            ifelse(direction == "directed", "directed", "undirected"), " graph.")
  }
  
  # Validate `weighting_links`
  assert_that(is.logical(weighting_links), 
              msg = "ERROR: 'weighting_links' must be a logical (TRUE or FALSE).")
  
  # Validate `link_weights` if weighting is enabled
  if (isTRUE(weighting_links)) {
    assert_that(is.list(link_weights), 
                msg = "ERROR: 'link_weights' must be a list when 'weighting_links' is TRUE.")
    
    assert_that(nrow(link_weights) == nrow(sites_status_matrix) && 
                  ncol(link_weights) == ncol(sites_status_matrix),
                msg = "ERROR: 'link_weights' and 'sites_status_matrix' must have the same dimensions.")
    
    if (verbose) message("Your links will be weighted with daily data entered in 'link_weights'.")
  } else {
    if (verbose) message("Your links will be unweighted, as defined in the LINK/NO_LINK structure.")
  }
  
  # Validate `weighting`
  assert_that(is.logical(weighting), 
              msg = "ERROR: 'weighting' must be a logical (TRUE or FALSE).")
  
  if (isTRUE(weighting)) {
    if (verbose) message("Your connectivity will be Weighted; connections will be multiplied by 'dist_matrix'.")
    assert_that(is.matrix(dist_matrix), 
                msg = "ERROR: 'dist_matrix' must be a matrix object when 'weighting' is TRUE.")
  } else {
    if (verbose) message("Your connectivity will be NON-weighted; connections will not be multiplied by any distance matrix.")
  }
  
  # Validate `legacy_length` and `legacy_effect`
  assert_that(is.numeric(legacy_length) && length(legacy_length) == 1, 
              msg = "ERROR: 'legacy_length' must be a single numeric value.")
  
  assert_that(length(legacy_effect) == legacy_length, 
              msg = paste("ERROR: The length of 'legacy_effect' is", length(legacy_effect), 
                          "but 'legacy_length' is", legacy_length, ". They must be the same."))
  
  assert_that(nrow(network_structure) == ncol(sites_status_matrix),
              msg = "ERROR: number of rows in network structure must match number of columns in sites_status_matrix.")
  
  ####_______________________________________________________________________
  # River network ####
  ####_______________________________________________________________________
  # Calculate the number of nodes of our network (used along the function)  
  numn_nodes <- ncol(sites_status_matrix)
  nsteps <- nrow(sites_status_matrix)
  
  if(weighting==TRUE){dist_matr <- dist_matrix}
  
  # Build the matrix corresponding to the num. of nodes multiplied by the number
  #of time steps
  ### This matrix is the template where we will put all the values.
  ST_matrix_raw <- matrix(nrow = numn_nodes, ncol = numn_nodes*2, 
                          data = value_no_s_link, 
                          dimnames = list(colnames(sites_status_matrix),
                                          rep(colnames(sites_status_matrix), 2)))
  ST_matrix_netwGraph_raw <- matrix(nrow = numn_nodes, ncol = numn_nodes, 
                                    data = 0L,
                                    dimnames = list(colnames(sites_status_matrix),
                                                    colnames(sites_status_matrix)))
  # Define the spatial connections of the matrix
  ### i.e., the rows or columns where we have to add the values of the connections 
  spa_connections <- seq_len(numn_nodes)
  temp_connections <- spa_connections + numn_nodes
  
  spa_temp_index_daily <- function(ST_matrix, ST_matrix_netwGraph, day) {
    #print(ST_matrix[1,])
    
    # Get site status data for current and next time steps
    ## time_step_1 is the present
    ## time_step_2 is the following step (the close future)
    time_step_1 <- sites_status_matrix[day, ]
    time_step_2 <- sites_status_matrix[day+1,]
    if(weighting_links==T){day_link_weights <- link_weights[day,2:interm_ncols]}
    
    #Simple fluvial network_____________________________________________________
    #Create an adjacency matrix for time step 1 whereby:
    #for sites that are wet, get the normal structure (direct connection to sites)
    #for sites that are dry, 0s to all sites
    ST_matrix_netwGraph[time_step_1 == 1,] <- network_structure[time_step_1 == 1,]
    diag_backup <- diag(ST_matrix_netwGraph)
    ST_matrix_netwGraph[time_step_1 == 0,] <- 0L
    diag(ST_matrix_netwGraph) <- diag_backup #####REMARK Mathis: not sure why this is needed
    #It's equivalent to [-site_step] in the original code
    #ST_matrix_netwGraph[spa_connections[site_step],
    #c(spa_connections[1]:spa_connections[numn_nodes])[-site_step]] <- 0
    
    # FLuvial SPATIAL links ____________________________________________________
    ## Fill the matrix section corresponding to the time_step based on river graph (igraph::)
    a <- graph_from_adjacency_matrix(ST_matrix_netwGraph, 
                                     mode = direction, 
                                     diag = FALSE)
    
    # Compute shortest path distances for all node pairs (igraph::)
    dist_matrix_day <- distances(a, mode = routing_mode,
                                 algorithm = "unweighted")
    
    # PROCESS SPATIAL STEP -----------------------------------------------------
    # Convert distances into binary connectivity (1 if connected, 0 if not)
    All_river_paths <- fifelse(!is.infinite(dist_matrix_day), 
                               value_s_link, value_no_s_link)
    
    # Weigh the links based on daily information of flow or strength of the link.
    if (weighting_links == TRUE) {
      All_river_paths <- All_river_paths * as.numeric(day_link_weights)
    }
    # Weigh the sites based on the distances between them (a pairwise matrix)
    if (weighting == TRUE) {
      All_river_paths <- All_river_paths * dist_matr
    }
    
    # Add the "All_river_paths" filled for each node in the "big" matrix
    ST_matrix[, spa_connections] <- ST_matrix[, spa_connections] + All_river_paths
    
    # PROCESS TEMPORAL STEPS ---------------------------------------------------
    # Create the matrix to drop the information of the shortest paths.
    All_river_paths <- matrix(nrow = length(time_step_1),
                              ncol = length(time_step_1),
                              data = value_no_t_link)
    
    #Calculate temporal changes in one step
    temp_changes <- time_step_1 - time_step_2
    
    #Determine stable (can be stable 1-1 or 0-0), lost, and gained indices
    stable_indices_1 <- which(temp_changes == 0 & time_step_1 == 1)
    stable_indices_0 <- which(temp_changes == 0 & time_step_1 == 0)
    lost_indices <- which(temp_changes == 1)
    gained_indices <- which(temp_changes == -1)
    
    #Compute direct spatial links for stable connected nodes (1->1)
    All_river_paths[stable_indices_1, ] <- fifelse(
      !is.infinite(dist_matrix_day[stable_indices_1, ]),
      value_t_link, value_no_t_link)
    
    #Indirect dispersal for lost nodes (1->0) 
    if (indirect_dispersal & (weighting == FALSE)) {
      All_river_paths[lost_indices, ] <- 
        All_river_paths[lost_indices, ] + 
        fifelse(!is.infinite(dist_matrix_day[lost_indices, ]), 
                value_t_link, value_no_t_link)
    }
    
    #Apply weights 
    if (weighting_links) {All_river_paths <- All_river_paths * day_link_weights}
    if (weighting) {All_river_paths <- All_river_paths * dist_matr}
    
    #Apply legacy effects and update ST_matrix
    for (leg_eff in seq_len(legacy_length)) {
      All_river_paths_legacy <- All_river_paths * legacy_effect[leg_eff]
      ST_matrix[, temp_connections] <- 
        ST_matrix[, temp_connections] + All_river_paths_legacy
      
      #Self-connection for stable connected nodes (diagonal assignment)
      value_t_link_modif <- if (weighting_links) { 
        value_t_link * day_link_weights 
      } else value_t_link
      ST_matrix[cbind(spa_connections[stable_indices_1], 
                      temp_connections[stable_indices_1])] <- 
        ST_matrix[cbind(spa_connections[stable_indices_1], 
                        temp_connections[stable_indices_1])] + value_t_link_modif
    }
    
    return(ST_matrix)
  }
  
  # Run function for every day
  ### For time step -1 because the last day does not have a "future" from 
  ### which to extract values. 
  ST_matrix <- lapply(1:(nsteps-1), function(day) {
    if (verbose) cat("We are at time unit", day, "of", (nsteps-1), "\n")
    spa_temp_index_daily(ST_matrix = ST_matrix_raw, 
                         ST_matrix_netwGraph = ST_matrix_netwGraph_raw,
                         day = day)
  }) %>% purrr::reduce(`+`) 
  
  ####_______________________________________________________________________
  # STconmat calculation ####
  ####_______________________________________________________________________
  # Find below the lines to calculate the "collapsing" matrix that just sums 
  # all the values of all the SPATIOTEMPORAL matrix
  # These pairwise matrix is called the STconmat.   
  ST_matrix_collapsed <- ST_matrix[,spa_connections] + ST_matrix[,temp_connections]
  
  if (any(c('all', 'STconmat', 'STcon') %in% output)) {
    #Divide by the number of days so we obtain the "per day" values
    ST_matrix_collapsed_standardized <- ST_matrix_collapsed/c(nsteps-1)
  }
  
  if (any(c('all', 'STconmat') %in% output)) {
    if (convert_to_integer) {
      ST_matrix_collapsed_standardized <- round(
        rounding_factor*ST_matrix_collapsed_standardized)
      storage.mode(ST_matrix_collapsed_standardized) <- "integer"
    }
    STconmat <- ST_matrix_collapsed_standardized
  } else {
    STconmat <- NULL
  }
  
  ####_______________________________________________________________________
  # STcon calculation ####
  ####_______________________________________________________________________
  if (any(c('all', 'STcon') %in% output)) {
    spt_conn <- apply(ST_matrix_collapsed_standardized, 1, sum)
    
    # "leng_correct" corrects the baseline connectivity of nodes due to their
    # position in the network. For example, with "out" routing mode, upstream
    # nodes will have higher values when considering their number of connections.
    if (standardize_neighbors) {
      leng_correct <- network_structure %>%
        graph_from_adjacency_matrix(mode = direction) %>%
        neighborhood_size(order = numn_nodes, mode = routing_mode) - 1 #to remove the connection to itself
      spt_conn <- spt_conn/leng_correct
    }
    
    if (convert_to_integer) {
      spt_conn <- as.integer(
        round(rounding_factor*spt_conn))
    }
  } else {
    spt_conn <- NULL
  }
  
  # OUTPUTS _______________________####
  if (!(any(c('all', 'Main_matrix') %in% output))) {
    ST_matrix <- NULL
  } else {
    if (convert_to_integer) {
      ST_matrix <- round(rounding_factor*ST_matrix)
      storage.mode(ST_matrix) <- "integer"
    }
  }
  
  main_output <- list(
    Main_matrix = ST_matrix,
    STconmat = STconmat,
    STcon = spt_conn  
  )
  
  return(main_output)
  ####_______________________________________________________________________
  
}# Function

