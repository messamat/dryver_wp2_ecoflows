#_______________________________________________________________________________#
#_______________________________________________________________________________#
#______________ spat_temp_index function _______________________________________#
#_______________________________________________________________________________#
#_______________________________________________________________________________#

# We present the function "spat_temp_index" a function to calculate Spatiotemporal
#connectiviy indices based on spatiotemporal graph networks. The aim of this 
#function is to provide a methodological framework from which to calculate these 
#indices based on high-frequency (or frequency-based) information obtained from natural ecosystems. 

# The current script is a modification of the original function introduced in
# the article titled: 
#Navigating through space and time: a methodological approach to quantify 
#spatiotemporal connectivity using flow intermittence data as a study case.
#David Cunillera-Montcusi, Jose Maria Fernandez-Calero, Sebastian Polsterl, 
#Julia Valera, Roger Argelich, Nuria Cid, Nuria Bonada, Miguel Canedo-Arguelles

# The original script was written by David Cunillera-Montcusi and modified
# by Mathis Messager
# Feedback: david.cunillera@dcm.cat and mathis.messager@mail.mcgill.ca

# Intermitence_dataset ____________#
# See below, an example of the type of matrix that must be entered in the function 
#as interm_dataset
#data.frame(
#MonitoredDays= c("Day1","Day2","Day3","Day4"), # An identifier for the monitored 
#days from ORDERED from the "oldest" to the "newest"
#StreamSite1=c(0,1,1,1), # More upstream site water presence record (1= water presence, 0= water absence)
#StreamSite2=c(0,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite3=c(0,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite4=c(0,0,0,1), # Following monitored site water presence record in descending order 
#StreamSite5=c(0,1,0,1), # Following monitored site water presence record in descending order 
#)
# This information is KEY to preserve the structure and functioning of the function.
#The system does not necessarily need to quantify water absence/presence 
#but each cell value must represent a feature defining connectivity "on" or "off" 
#and that can be transmitted to build ecologically meaningful links in a
#spatiotemporal graph. 

# Sites coordinates ____________#
# The coordinates of the sites must be located at columns 3 and 4 of the 
# data.frame in order to be read properly.

# LINK and NO-LINK values ____________#
# The meaning of a "LINK" (Input matrix value= 1) or a "NO-LINK" (Input matrix value= 0) 
# can be modified inside the function according
#to the different interpretations that one once to give to the values (e.g. do 
#we want to quantify connectivity? dispersal resistance?)
#see the above mentioned paper to see an specific example. 

# Set direction ____________#
# direction can be either "directed" or "undirected". This feature will modify 
# the way the graph is being controlled. 

# sense ____________#
# Sense is referring to the direction that will be considered for the graph 
# connectivity when directed. Specially for centrality metrics. 
# It can be either "in", "out", "all". See ?igraph or ?igraph::closeness
# for a better understanding. 

# Weighting ____________#
# weighting either FALSE or TRUE --> The value of the weight is defined by the 
# matrix added (must be in a distance matrix format) and 
#can consider any type of distance between pairs of monitored sites 
# (euclidean, environmental, topographic, ...)
#dist_matrices--> corresponds to the attached matrix representing the 
# "distances" between sites. 

# network_structurecture ____________#
# Network structure is a matrix that corresponds to the "basic" connections that 
# can be possible. An adjacency matrix where all sites are connected among them. 
# Should be equivalent to the network that you would expect from a "fully"
# connected network. 

# weighting_links=FALSE & link_weights ____________#
# In case that instead of using the values used as "LINKS/NO-LINKS" that are 
# "qualitative" we have specific data for every time unit that can be used to 
# quantify the connections between sites (e.g. flow data, wind strenght etcetera). 
# This information can be incorporated by selecting weighting_links=TRUE and then 
# providing a list of data.frames with exactly the same information as the 
# "Intermitence_dataset", having the same amount of rows and columns (e.g. dauly 
# flow data in each site).

# LINKS / NO-LINKS  ____________#
# value_LINK is the value attributed for each "effective link", which is a
# connection between two nodes (a line)
# - value_s_link for spatial links
# - value_t_link for temporal links
# value_NO_link is the value for each "effective disconnection", which is a 
# "connection"void" between (a white space)
# - value_NO_s_link for spatial links
# - value_NO_t_link for temporal links

# Legacy effects & Legacy length  ____________#
#Legacy effects are a way to quantify spatiotemporal connectivity during several
#time units at the same time. This means that we  can quantify temporal links for 
#more than 1 time unit and modify the weight or relevance that each time unit has. 
#For example: 
# - Legacy length defines the "number" of time units considered. By default it is 
#set to 1 which means that the function will consider only 1 time unit 
#(the connectivity of time T1 will only be considered for T2). If value is set 
# at legacy length=3 the connectivity of 
#T1 will be considered for T2, T3, and T4 (the connections present at T1 will 
# be also considered for the 3 following T).
# - Legacy effect defines the weight given to the LINK / NO-LINK values for each 
#time unit. It is a vector with values ranging from
#0 to 1 that will modulate the relevance of a connection through time. 

spat_temp_index_edit <- function(interm_dataset, 
                                 network_structure,
                                 direction,
                                 sense = "out",
                                 weighting = FALSE,
                                 dist_matrices,
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
    
    assert_that(nrow(link_weights) == nrow(interm_dataset) && 
                  ncol(link_weights) == ncol(interm_dataset),
                msg = "ERROR: 'link_weights' and 'interm_dataset' must have the same dimensions.")
    
    if (verbose) message("Your links will be weighted with daily data entered in 'link_weights'.")
  } else {
    if (verbose) message("Your links will be unweighted, as defined in the LINK/NO_LINK structure.")
  }
  
  # Validate `weighting`
  assert_that(is.logical(weighting), 
              msg = "ERROR: 'weighting' must be a logical (TRUE or FALSE).")
  
  if (isTRUE(weighting)) {
    if (verbose) message("Your connectivity will be Weighted; connections will be multiplied by 'dist_matrices'.")
    assert_that(is.matrix(dist_matrices), 
                msg = "ERROR: 'dist_matrices' must be a matrix object when 'weighting' is TRUE.")
  } else {
    if (verbose) message("Your connectivity will be NON-weighted; connections will not be multiplied by any distance matrix.")
  }
  
  # Validate `legacy_length` and `legacy_effect`
  assert_that(is.numeric(legacy_length) && length(legacy_length) == 1, 
              msg = "ERROR: 'legacy_length' must be a single numeric value.")
  
  assert_that(length(legacy_effect) == legacy_length, 
              msg = paste("ERROR: The length of 'legacy_effect' is", length(legacy_effect), 
                          "but 'legacy_length' is", legacy_length, ". They must be the same."))

  assert_that(nrow(network_structure) == ncol(interm_dataset),
              msg = "ERROR: number of rows in network structure must match number of columns in interm_dataset.")
  
  ####_______________________________________________________________________
  # River network ####
  ####_______________________________________________________________________
  # We calculate the number of nodes of our network (used along the function)  
  numn_nodes <- ncol(interm_dataset)
  nsteps <- nrow(interm_dataset)
  
  if(weighting==TRUE){dist_matr <- dist_matrices}
  
  # We built the matrix corresponding to the num. of nodes multiplied by the
  # DAYS of HOBOS that we have
  ### This matrix is the "giant" template where we will put all the values.
  ST_matrix_raw <- matrix(nrow = numn_nodes, ncol = numn_nodes*2, 
                          data = value_no_s_link, 
                          dimnames = list(colnames(interm_dataset),
                                          rep(colnames(interm_dataset), 2)))
  ST_matrix_netwGraph_raw <- matrix(nrow = numn_nodes, ncol = numn_nodes, 
                                    data = 0L,
                                    dimnames = list(colnames(interm_dataset),
                                                    colnames(interm_dataset)))
  # First we define the spatial connections of the matrix
  ### Also known as the rows or columns at which we have to add the values of 
  #the connections 
  spa_connections <- seq_len(numn_nodes)
  temp_connections <- spa_connections + numn_nodes
  
  spa_temp_index_daily <- function(ST_matrix, ST_matrix_netwGraph, day) {
    #print(ST_matrix[1,])
    
    # We obtain the time steps:
    ## time_step_1 is the present
    ## time_step_2 is the following step (the close future)
    time_step_1 <- interm_dataset[day, ]
    time_step_2 <- interm_dataset[day+1,]
    if(weighting_links==T){day_link_weights <- link_weights[day,2:interm_ncols]}
    
    #Simple fluvial network_____________________________________________________
    #Create an adjacancy matrix for time step 1 whereby:
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
    ## Fill the matrix section corresponding to the time_step based on river graph 
    a <- igraph::graph_from_adjacency_matrix(ST_matrix_netwGraph, 
                                             mode = direction, 
                                             diag = FALSE)
    
    # Compute shortest path distances for all node pairs
    dist_matrix_day <- igraph::distances(a, mode = sense,
                                         algorithm = "unweighted")
    
    # PROCESS SPATIAL STEP -----------------------------------------------------
    # Convert distances into binary connectivity (1 if connected, 0 if not)
    All_river_paths <- fifelse(!is.infinite(dist_matrix_day), 
                               value_s_link, value_no_s_link)
    
    # Weight the links base on daily information of flow or strength of the link.
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
    # We create the matrix where we will drop the information of the shortest paths.
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
  
  if (any(c('all', 'STconmat') %in% output)) {
    ST_matrix_collapsed_standardized <- ST_matrix_collapsed/c(nsteps-1)
    if (convert_to_integer) {
      ST_matrix_collapsed_standardized <- round(
        rounding_factor*ST_matrix_collapsed_standardized)
      storage.mode(ST_matrix_collapsed_standardized) <- "integer"
    }
  } else {
    ST_matrix_collapsed_standardized <- NULL
  }

  ####_______________________________________________________________________
  # STcon calculation ####
  ####_______________________________________________________________________
  if (any(c('all', 'STcon') %in% output)) {
    spt_conn <- apply(ST_matrix_collapsed, 1, sum)
    
    # "leng_correct" is a reverse vector (from big to small) used to correct
    # the fact that uperstream nodes will have higher values when considering its 
    #number of connections. As I am "node 1" my number of connections will be higher
    #than "node 10". IF WE FOLLOW THE RIVER DOWNSTREAM!
    if (standardize_neighbors) {
      aa <- graph_from_adjacency_matrix(network_structure, mode = direction)
      leng_correct <- neighborhood_size(aa, order = numn_nodes, mode = sense) - 1 #to remove the connection to itself
      spt_conn <- spt_conn/leng_correct
    }
    # We divide by the number of days so we obtain the "per day" values
    spt_conn <- spt_conn/c(nsteps-1)
    
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
  
  Main_output <- list(
    Main_matrix = ST_matrix,
    STconmat = ST_matrix_collapsed_standardized,
    STcon = spt_conn  
  )
  
  return(Main_output)
  ####_______________________________________________________________________
  
}# Function

