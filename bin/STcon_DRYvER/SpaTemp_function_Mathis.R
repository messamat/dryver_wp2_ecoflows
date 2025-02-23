#__________________________________________________________________________________________________________#
#__________________________________________________________________________________________________________#
#______________ spat_temp_index function __________________________________________________________________#
#__________________________________________________________________________________________________________#
#__________________________________________________________________________________________________________#

# We present the function "spat_temp_index" a function to calculate Spatiotemporal connectiviy indices based on 
#spatiotemporal graph networks. The aim of this function is to provide a methodological framework from which to 
#calculate these indices based on high-frequency (or frequency-based) information obtained from natural ecosystems. 

# The current script as well as all the information contained in it and its tutorial are part of the article titled: 
#Navigating through space and time: a methodological approach to quantify spatiotemporal connectivity using flow intermittence data as a study case.
#David Cunillera-Montcusi1,2,3,4*, Jose Maria Fernandez-Calero1,2, Sebastian Polsterl5, Julia Valera1, Roger Argelich1, Nuria Cid1,6, Nuria Bonada1,2,
#Miguel Canedo-Arguelles1,7,8
#1- FEHM-Lab (Freshwater Ecology, Hydrology and Management), Departament de Biologia Evolutiva, Ecologia i Ciencies Ambientals, Facultat de Biologia,Universitat de Barcelona (UB), Diagonal 643, 08028 Barcelona, Spain.
#2- Institut de Recerca de la Biodiversitat (IRBio), Universitat de Barcelona (UB), Diagonal 643, 08028 Barcelona, Spain.
#3. Departamento de Ecologia y Gestion Ambiental, Centro Universitario Regional del Este (CURE), Universidad de la Repulica, Tacuarembo? s/n, Maldonado, Uruguay.
#4. GRECO, Institute of Aquatic Ecology, University of Girona, Girona, Spain
#5- The Lab for Artificial Intelligence in Medical Imaging (AI-Med), Department of Child and Adolescent Psychiatry, Ludwig-Maximilians-Universitat, Waltherstrae 23, 80337 Munich, Germany
#6- IRTA  Marine and Continental Waters Programme, Ctra de Poble Nou Km 5.5, E43540, La Rapita, Catalonia, Spain 
#7- Institut de Recerca de l'Aigua (IdRA), Universitat de Barcelona (UB), Diagonal 643, 08028 Barcelona, Catalonia, Spain.
#8- Institute of Environmental Assessment and Water Research (IDAEA-CSIC), Carrer de Jordi Girona, 18-26, 08034 Barcelona

#*Corresponding author: david.cunillera@dcm.cat

# The script and the follwouing functions have been written by David Cunillera-Montcusi. For any question, comment or 
#feedback contact him at: david.cunillera@dcm.cat 

# Mainly, the information that should be contained in the datasets must represent some sort of habitat availability 
# like habitat presence/absence (e.g. aquatic habitats wet/dry phases) or habitat connectivity strengths (e.g. flows 
# in rivers) or any other type of interpretation that has ecological meaning and that can be related with spatial and 
# temporal patterns. 

# Intermitence_dataset ____________#
# See below, an example of the type of matrix that must be entered in the function as interm_dataset
#data.frame(
#MonitoredDays= c("Day1","Day2","Day3","Day4"), # An identifier for the monitored days from ORDERED from the "oldest" to the "newest"
#StreamSite1=c(0,1,1,1), # More upstream site water presence record (1= water presence, 0= water absence)
#StreamSite2=c(0,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite3=c(0,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite4=c(0,0,0,1), # Following monitored site water presence record in descending order 
#StreamSite5=c(0,1,0,1), # Following monitored site water presence record in descending order 
#StreamSite6=c(0,0,0,1), # Following monitored site water presence record in descending order 
#StreamSite7=c(0,1,0,1), # Following monitored site water presence record in descending order 
#StreamSite8=c(0,0,0,1), # Following monitored site water presence record in descending order 
#StreamSite9=c(1,1,1,1), # Following monitored site water presence record in descending order 
#StreamSite10=c(1,1,1,1),# Following monitored site water presence record in descending order 
#)
# This information is KEY to preserve the structure and functioning of the function. The system does not necessarily need to quantify 
#water absence/presence but each cell value must represent a feature defining connectivity "on" or "off" and that can be transmitted 
#to built ecologically meaningful links in a spatiotemporal graph. 

# Sites coordinates ____________#
# The coordinates of the sites must be located at columns 3 and 4 of the data.frame in order to be read properly.

# LINK and NO-LINK values ____________#
# The meaning of a "LINK" (Input matrix value= 1) or a "NO-LINK" (Input matrix value= 0) can be modified inside the function according
#to the different interpretations that one once to give to the values (e.g. do we want to quantify connectivity? dispersal resistance?)
#see the above mentioned paper to see an specific example. 

# Set direction ____________#
# direction can be either "directed" or "undirected". This feature will modify the way the graph is being controlled. 

# sense ____________#
# Sense is referring to the direction that will be considered for the graph connectivity when directed. Specially for centrality metrics. 
# It can be either "in", "out", "all". See ?igraph or ?igraph::closeness for a better understanding. 

# Weighting ____________#
# weighting either FALSE or TRUE --> The value of the weight is defined by the matrix added (must be in a distance matrix format) and 
#can consider any type of distance between pairs of monitored sites (euclidean, environmental, topographic, ...)
#dist_matrices--> corresponds to the attached matrix representing the "distances" between sites. 

# Network_structure ____________#
# Network structure is a matrix that corresponds to the "basic" connections that can be possible. An adjacency matrix where all 
#sites are connected among them. Should be equivalent to the network that you would expect from a "fully" connected network. 

# weighting_links=FALSE & link_weights ____________#
# In case that instead of using the values used as "LINKS/NO-LINKS" that are "qualitative" we have specific data for every time unit 
#that can be used ot quantify the connections between sites (e.g. flow data, wind strenght etcetera). This information can be
#incorporated by selecting weighting_links=TRUE and then providing a list of data.frames with exactly the same information as
#the "Intermitence_dataset", having the same amount of rows and columns (e.g. dauly flow data in each site).

# LINKS / NO-LINKS  ____________#
# value_LINK is the value attributed for each "effective link", which is a connection between two nodes (a line)
# - value_s_link for spatial links
# - value_t_link for temporal links
# value_NO_link is the value for each "effective disconnection", which is a "connection"void" between (a white space)
# - value_NO_s_link for spatial links
# - value_NO_t_link for temporal links

# Legacy effects & Legacy lenght  ____________#
#Legacy effects are a way to quantify spatiotemporal connectivity during several time units at the same time. This means that we 
#can quantify temporal links for more than 1 time unit and modify the weight or relevance that each time unit has. 
#For example: 
# - Legacy lenght defines the "number" of time units considered. By default it is set to 1 which means that the function will consider 
#only 1 time unit (the connectivity of time T1 will only be considered for T2). If value is set at legacy lenght=3 the connectivity of 
#T1 will be considered for T2, T3, and T4 (the connections present at T1 will be also considered for the 3 following T).
# - Legacy effect defines the weight given to the LINK / NO-LINK values for each time unit. It is a vector with values ranging from
#0 to 1 that will modulate the relevance of a connection through time. 

# WARNING !# WARNING !# WARNING !# WARNING !# WARNING !# WARNING !# WARNING !
# WARNING !# WARNING !# WARNING !# WARNING !# WARNING !# WARNING !# WARNING !
# The values of interm_dataset,Sites_coordinates,Network_stru must be entered in a list format, even if only 1 matrix is used. 
#This is done in order to allow the possibility to enter several streams in one call and obtain the results according to that.
# In case of using dist_matrices and link_weights they must also be entered as list(). 

spat_temp_index <- function(interm_dataset, 
                            Sites_coordinates,
                            direction,
                            sense="out",
                            weighting=FALSE,
                            dist_matrices,
                            weighting_links=FALSE,
                            link_weights,
                            Network_stru,
                            value_s_link=1,
                            value_t_link=1,
                            value_no_s_link=0,
                            value_no_t_link=0,
                            legacy_effect=1,
                            legacy_length=1){
  
  require(gridExtra);require(igraph);require(shp2graph);require(sna,quietly = T,warn.conflicts = F);
  require(tidyverse);require(viridis);require(doParallel);require(ggnetwork)
  
  if(direction=="directed"){
    cat("Your river will be considered as a directed graph","\n")}
  if(direction=="undirected"){
    cat("Your river will be considered as an undirected graph","\n")}
  
  # We select the corresponding distance matrix
  if(weighting_links==T){
    cat("Your links will be weighted with daily data entered in the 'link_weights'","\n")
    if(nrow(link_weights[[1]])!=nrow(interm_dataset[[1]]) & ncol(link_weights[[1]])!=ncol(interm_dataset[[1]])){
      return(cat("!!!ERROR: Intermitence dataset and link_weights must have the same dimensions!","\n"))}
  }
  if(weighting_links==F){cat("Your links will be normal, as defined in the LINK/NO_LINK","\n")}
  
  # We select the corresponding distance matrix
  if(weighting==T){cat("Your connectivity will be Weighted, connections will be multiplied by 'dist_matrices'","\n")}
  if(weighting==F){cat("Your connectivity will be NON weighted, connections will not be multiplied by any distance matrix","\n")}
  
  if(legacy_length!=length(legacy_effect)){
    return(cat("!!!ERROR: The length of your legacy effects is", length(legacy_effect), "and your legacy length is", legacy_length,"! They must be the same!", "\n"))}
  
  if(length(which(c(is.list(interm_dataset),is.list(Sites_coordinates),is.list(Network_stru))==F))>0){
    return(cat("Your interm_dataset,Sites_coordinates or Network_stru must be list objects", "\n"))}
  if(weighting==T & is.list(dist_matrices)==F){return(cat("!!!ERROR: Your distance matrix must be a list object"))}
  
  # Simple river network________________####
  # Building river networks based on just directional network.
  Simple_river_network <- list()
  Simple_river_network_maps <- list()
  for (river in 1:length(interm_dataset)) {
    print("Building river networks based on just directional network")
    
    if(is.numeric(Sites_coordinates[[river]][,3])==F & is.numeric(Sites_coordinates[[river]][,4])==F){
      return(cat("!!!ERROR: X and Y coordinates must be at columns 3 and 4 of the Sites_coordinates"))}
    
    ST_matrix_out <- matrix(nrow = ncol(interm_dataset[[river]])-1,ncol = ncol(interm_dataset[[river]])-1, 
                            data=value_no_s_link)
    spa_connections <-seq(1,ncol(interm_dataset[[river]])-1,1)
    time_step_1 <- rep(1,ncol(interm_dataset[[river]])-1)
    
    for (site_step in 1:c(ncol(interm_dataset[[river]])-1)) {
      #Simple spatial links _______________________
      if(time_step_1[site_step]==1){
        ST_matrix_out[spa_connections[site_step],which(Network_stru[[river]][site_step,]==1)] <- value_s_link
      }else{
        ST_matrix_out[spa_connections[site_step],which(Network_stru[[river]][site_step,]==1)] <- value_no_s_link
      }
    }
    Simple_river_network[[river]] <- ST_matrix_out
  
  }
  ####_______________________________________________________________________
  
  ####_______________________________________________________________________
  # River network ####
  ####_______________________________________________________________________
  # River network matrix BUILDING ####
  ####_______________________________________________________________________
  pack_check <- search()
  pack_check_val <- length(which(pack_check=="package:sna"))
  if(pack_check_val>0){detach("package:sna", unload = TRUE)}
  # Below there is the function who builds the MATRIX ponderating SPATIAL lINKS=1 and TEMPORAL LINKS=1
  ST_matrix_rivers <- list()

  #Parallelization parameters
  # cores <- detectCores() #Number of cores in computer
  # cl <- makeCluster(cores[1]-1) #not to overload your computer
  # registerDoParallel(cl)
  
  out_Matrix_LIST <- list()
  pack_check <- search()
  pack_check_val <- length(which(pack_check=="package:sna"))
  if(pack_check_val>0){detach("package:sna", unload = TRUE)}
  for (river in 1:length(interm_dataset)) {
  #for (river in 1:length(interm_dataset)) { #- With this it takes 6'26''
  # We calculate the number of nodes of our network (used along the function)  
    numn_nodes <- ncol(interm_dataset[[river]])-1
    
    if(weighting==TRUE){dist_matr <- dist_matrices[[river]]}
    
    # We built the matrix corresponding to the num. of nodes multiplied by the DAYS of HOBOS that we have
    ### This matrix is the "giant" themplate where we will put all the values.
    ST_matrix <- matrix(nrow = numn_nodes,ncol = numn_nodes*2, data=value_no_s_link)
    ST_matrix_netwGraph <- matrix(nrow = numn_nodes,ncol = numn_nodes, data=0)
    
    # Once created the template we start to fill it for every day
    ### We fill it for Days (or time)-1 because the last day does not have a "future" from which to extract values. 
    for (days in 1:c(length(interm_dataset[[river]][,1])-1)) {
      cat("We are at time unit", days, "of", (length(interm_dataset[[river]][,1])-1), "and at river", river,"\n")
      # First we define the spatial connections of the matrix
      ### Also known as the rows or columns at which we have to add the values of the connections 
      spa_connections <-seq(1,length(colnames(interm_dataset[[river]]))-1,1)#+((days-1)*numn_nodes)
      
      # We obtain the time steps:
      ## time_step_1 is the present
      ## time_step_2 is the following step (the close future)
      time_step_1 <- interm_dataset[[river]][days,2:ncol(interm_dataset[[river]])]
      time_step_2 <- interm_dataset[[river]][days+1,2:ncol(interm_dataset[[river]])]
      if(weighting_links==T){day_link_weights <- link_weights[[river]][days,2:ncol(interm_dataset[[river]])]}
      
      #Simple fluvial network_______________________
      ## This step fills "the diagonal" of each time_step following the direction of the river
      ## it basically connects the river in a dendritic structure.
      for (site_step in 1:c(length(time_step_1))) {
        #print(site_step)
        if(time_step_1[site_step]==1){
          ST_matrix_netwGraph[spa_connections[site_step],
                              c(spa_connections[1]:spa_connections[numn_nodes])] <- as.numeric(Network_stru[[river]][site_step,])
        }else{
          ST_matrix_netwGraph[spa_connections[site_step],
                              c(spa_connections[1]:spa_connections[numn_nodes])[-site_step]] <- 0
        }
      }
      
      # FLuvial SPATIAL links ___________________________________________________________________________________________________________________
      # Now the party begins. 
      ## Here we fill the matrix section corresponding to the time_step based on the river graph based on a dendritic. 
      require(igraph)
      # We create the graph
      a <- graph.adjacency(ST_matrix_netwGraph[spa_connections[1]:spa_connections[numn_nodes],
                                               spa_connections[1]:spa_connections[numn_nodes]], 
                           mode=direction,diag = FALSE)
      
      # We create the matrix where we will drop the information of the shortest paths.
      ## We will fill "1" or "0" according to the shortest paths. 
      All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_no_s_link)
      #All_river_paths[upper.tri(All_river_paths)] <- value_no_s_link
      
      # For each path (e.g., from node 1 to node 7) we "check" the length of the shortest path. 
      ## check = 0 means that the graph is disconnected.
      ## check bigger than 0 means that the graph is connected.
      for (every_path in 1:c(length(time_step_1))){
        check <- length(all_shortest_paths(a, every_path, 1:numn_nodes, mode = sense)$res)
        if (check==0) {
          site <- 0 
        }else{
          # If bigger than 0. we create a sequence from the path to downstream.
          neigh <- (c(1:numn_nodes)[-every_path])
          connect_loc <- all_shortest_paths(a, every_path,neigh, mode = sense)$nrgeo
          connect_loc[every_path] <- 0
          site <- which(connect_loc>0)
        }
        # We fill the "All_river_paths" with 1 on the connections concerning to each "row" or node.
        ## Site is the vector with the connections (follwing the river downstream).
        ## When "0" site does not correspond to any row... so the "1" does not go anywhere. 
        All_river_paths[every_path,site] <- value_s_link
        # We weight the links base on daily information of flow or strength of the link.
        if(weighting_links==T){All_river_paths[every_path,] <-  as.numeric(All_river_paths[every_path,]*as.numeric(day_link_weights[every_path]))}
        # We weight the sites for the distances between them (a pairwise matrix)
        if(weighting==T){All_river_paths[every_path,] <-  as.numeric(All_river_paths[every_path,]*dist_matr[every_path,])}
      }
      
      # We add the "All_river_paths" filled for each node in the "big" matrix specific sites
      ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                spa_connections[1]:spa_connections[numn_nodes]] <- All_river_paths+ST_matrix[spa_connections[1]:spa_connections[numn_nodes],
                                                                                             spa_connections[1]:spa_connections[numn_nodes]]
      
      # In the following lines we continue the party towards temporal steps
      for (site_step in 1:length(time_step_1)) {
        # We created "All_river_paths" for temporal
        All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = value_no_t_link)
        
        # FLuvial TEMPORAL DIRECT links ___________________________________________________________________________________________________________________
        ## We generate the temporal connectins
        temp_connections <-seq(1+numn_nodes,length(colnames(interm_dataset[[river]]))-1+numn_nodes,1)#+((days)*numn_nodes) 
        
        ## We then evaluate the difference between the two time steps and therefore we quantify:
        # - Stable links: 0 (WARNING: stable links can be stable 1-1 or 0-0!)
        # - Lost links: 1
        # - Gained links: -1
        ## - This values 0,1,-1 define what we will do with the links that match such pattern
        temp_change <- time_step_1[site_step]-time_step_2[site_step]
        
        #Stable links (when temp_change=0)
        ## The mechanics is the same as previously.
        if(temp_change==0){# Temporal change is constant 
          if(time_step_1[site_step]==1){# This temporal change implies going from 1 to 1 (so a real stable connected link)
            # We then do the same as before, check, substitute and add "1" or 0 depending if the connection 
            # following the river is flow.
            #for (every_path in 1:length(time_step_1)){
            check <- length(all_shortest_paths(a, site_step, 1:numn_nodes, mode = sense)$res)
            if (check==0) {
              site <-0
            }else{
              neigh <- (c(1:numn_nodes)[-site_step])
              connect_loc <- all_shortest_paths(a, site_step,neigh, mode = sense)$nrgeo
              connect_loc[site_step] <- 0
              site <- which(connect_loc>0)
            }
            All_river_paths[site_step,site] <- value_t_link
            
            # We weight the links base on daily information of flow or strength of the link.
            if(weighting_links==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*as.numeric(day_link_weights[site_step]))}
            # We weight
            if(weighting==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*dist_matr[site_step,])}
            
            for (leg_eff in 1:legacy_length) {
              All_river_paths_legacy <- All_river_paths[site_step,]*legacy_effect[leg_eff]
              # TEMPORAL LINKS are filled in the "future" of our current matrix. This means that we are filling the matrix in 
              # in the diagonal of our "time step" for spatial links but we add the temporal links in the following time step. 
              # so, we evaluate here the present (time step 1) and the future (time step 2) but we register it as the past of the future (at time step 2)
              ST_matrix[spa_connections[site_step],
                        temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths_legacy+ST_matrix[spa_connections[site_step],
                                                                                                              temp_connections[1]:temp_connections[numn_nodes]]
              # Here we add the temporal "link" between "himself". If the link is stable and connected (from 1 to 1), we fill the 
              # diagonal value accordingly. Therefore, we will be able to evaluate the relationship between "himself". Kind of 
              # Tot_Num indicator.
              value_t_link_modif <- value_t_link
              if(weighting_links==T){value_t_link_modif <- value_t_link_modif*as.numeric(day_link_weights[site_step])}
              ST_matrix[spa_connections[site_step],
                        temp_connections[site_step]] <- (value_t_link_modif*legacy_effect[leg_eff])+ST_matrix[spa_connections[site_step],temp_connections[site_step]]
            }
          }else{# Here we check if the temporal change implies going from 0 to 0 (so a stable disconnected link). Then we put 0
            # We weight the links base on daily information of flow or strength of the link.
            if(weighting_links==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*as.numeric(day_link_weights[site_step]))}
            # We weight
            if(weighting==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*dist_matr[site_step,])}
            for (leg_eff in 1:legacy_length) {
              All_river_paths_legacy <- All_river_paths[site_step,]*legacy_effect[leg_eff]
              ST_matrix[spa_connections[site_step],
                        temp_connections[1]:temp_connections[numn_nodes]] <-All_river_paths_legacy+ST_matrix[spa_connections[site_step],
                                                                                                             temp_connections[1]:temp_connections[numn_nodes]]
            }
          }
        }
        
        #Lost links (when temp_change=1)
        ## This just needs to be filled with zeros... so no need to use "All_river_paths"
        if(temp_change==1){
          # We weight the links base on daily information of flow or strength of the link.
          if(weighting_links==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*as.numeric(day_link_weights[site_step]))}
          # We weight
          if(weighting==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*dist_matr[site_step,])}
          for (leg_eff in 1:legacy_length) {
            All_river_paths_legacy <- All_river_paths[site_step,]*legacy_effect[leg_eff]
            ST_matrix[spa_connections[site_step],
                      temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths_legacy+ST_matrix[spa_connections[site_step],
                                                                                                            temp_connections[1]:temp_connections[numn_nodes]]
          }
        }
        #Gained links (when temp_change=-1)
        ## It is a "gain" but it means that "in the present" (time step 1), the node is still disconnected. So it =0
        if(temp_change==-1){
          # We weight the links base on daily information of flow or strength of the link.
          if(weighting_links==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*as.numeric(day_link_weights[site_step]))}
          # We weight
          if(weighting==T){All_river_paths[site_step,] <- as.numeric(All_river_paths[site_step,]*dist_matr[site_step,])}
          for (leg_eff in 1:legacy_length) {
            All_river_paths_legacy <- All_river_paths[site_step,]*legacy_effect[leg_eff]
            ST_matrix[spa_connections[site_step],
                      temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths_legacy+ST_matrix[spa_connections[site_step],
                                                                                                            temp_connections[1]:temp_connections[numn_nodes]]
          }
        }
        
        # FLuvial TEMPORAL INDIRECT links ___________________________________________________________________________________________________________________
        ## Indirect links are those links defined here as the ones that are "lost for the first time" (so temp_change=1).
        ## When this occurs we asign an "extra" 1 in that particular case. Assuming that when the node dries for the first time
        ## there is an increase in "dispersal" (downstream directed).
        ### This only occurs when there is a loss of a previously wet node (from 1 in the present to 0 in the future).
        if(temp_change==1){
          ######################################################################
          #MODIF MATHIS: NEED TO RESET all_river_paths. Otherwise,
          #when using weighting == T, multiply 0.1*distances by 0.1*distances for those 
          #that are not connected.
          All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),
                                    data = value_no_t_link)
          #######################################################################
          #All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),data = 0) 
          #All_river_paths[upper.tri(All_river_paths)] <- value_no_t_link
          check <- length(all_shortest_paths(a, site_step, 1:numn_nodes, mode = sense)$res)
          if (check==0) {
            site <-0
          }else{
            neigh <- (c(1:numn_nodes)[-site_step])
            connect_loc <- all_shortest_paths(a, site_step,neigh, mode = sense)$nrgeo
            connect_loc[site_step] <- 0
            site <- which(connect_loc>0)
          }
          # We fill the sites with the value
          All_river_paths[site_step,site] <- value_t_link
          # We weight the links base on daily information of flow or strength of the link.
          if(weighting_links==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*as.numeric(day_link_weights[site_step]))}
          # We weight
          if(weighting==T){All_river_paths[site_step,] <-  as.numeric(All_river_paths[site_step,]*dist_matr[site_step,])}
          # We pass it to the main matrix
          for (leg_eff in 1:legacy_length) {
            All_river_paths_legacy <- All_river_paths[site_step,]*legacy_effect[leg_eff]
            ST_matrix[spa_connections[site_step],
                      temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths_legacy+ST_matrix[spa_connections[site_step],
                                                                                                            temp_connections[1]:temp_connections[numn_nodes]]
          }
          
        }# End of the if
        
      }# Site_step closing
    }# Days closing
    
    out_Matrix <- list(ST_matrix)
    out_Matrix_LIST[[river]] <- out_Matrix
  }# Loop for every river entered in the lists
  
  #parallel::stopCluster(cl) #Stop parallel computing
  
  # Exctracring the results into different lists
  for (river in 1:length(interm_dataset)) {ST_matrix_rivers[[river]] <- out_Matrix_LIST[[river]][[1]]}
  
  ####_______________________________________________________________________
  # STconmat calculaiton ####
  ####_______________________________________________________________________
  # Find below the lines to calculate the "collapsing" matrix that just summs all the values of all the SPATIOTEMPORAL matrix
  # These pairwise matrix is called the STconmat.   
  ST_matrix_out_out <- list()
  
  for (river in 1:length(interm_dataset)) {
    numn_nodes <- ncol(interm_dataset[[river]])-1
    temp_connections <-seq(1+numn_nodes,length(colnames(interm_dataset[[river]]))-1+numn_nodes,1)
    spa_connections <-seq(1,length(colnames(interm_dataset[[river]]))-1,1)
    
    # We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
    out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = value_no_s_link)
    
    Spatial_matrix <-  ST_matrix_rivers[[river]][,spa_connections]
    Temporal_matrix <- ST_matrix_rivers[[river]][,temp_connections]
    
    out_out <- Spatial_matrix+Temporal_matrix
    out_out <- out_out/c(length(interm_dataset[[river]][,1])-1)
    
    # We save the collapsed matrix
    ST_matrix_out_out[[river]] <- out_out
  }
  ####_______________________________________________________________________
  
  ####_______________________________________________________________________
  # STcon calculaiton ####
  ####_______________________________________________________________________
  ST_connectivity_value <- list()
  for (river in 1:length(interm_dataset)) {
    # We already know this value
    numn_nodes <- ncol(interm_dataset[[river]])-1
    temp_connections <-seq(1+numn_nodes,length(colnames(interm_dataset[[river]]))-1+numn_nodes,1)
    spa_connections <-seq(1,length(colnames(interm_dataset[[river]]))-1,1)
    
    # We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
    out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = value_no_s_link)
    
    Spatial_matrix <-  ST_matrix_rivers[[river]][,spa_connections]
    Temporal_matrix <- ST_matrix_rivers[[river]][,temp_connections]
    
    out_out <- Spatial_matrix+Temporal_matrix
    
    # "leng_correct" is a reverse vector (from big to small) used to correct the fact that uperstream nodes will have higher values when 
    # considering its number of connections. As I am "node 1" my number of connections will be higher tan "node 10". IF WE FOLLOW THE RIVER DOWNSTREAM!
    leng_correct <- c()
    aa <- graph_from_adjacency_matrix(as.matrix(Network_stru[[river]]),mode = "directed")
    for (neigg in 1:numn_nodes) {
      neighbour<- all_shortest_paths(aa, from = neigg, to = c(1:numn_nodes)[-neigg],mode = sense)$nrgeo
      neighbour[neigg] <- 0
      leng_correct[neigg] <- length(which(neighbour>0))
    }
    
    spt_conn <-apply(out_out,1,sum)/leng_correct
    # We divide by the number of days so we obtain the "per day" values
    spt_conn<- spt_conn/c(length(interm_dataset[[river]][,1])-1)
    
    ST_connectivity_value[[river]] <- spt_conn  
  }
  
  
  ####_______________________________________________________________________
  # OUTPUTS _______________________####
  # Global matrix
  NonW_ST_matrix_rivers <- ST_matrix_rivers
  
  # Spatiotemporal matrix 
  NonW_ST_matrix_out_out <- ST_matrix_out_out
  
  # Spatiotemporal connectivity 
  NonW_ST_connectivity_value <- ST_connectivity_value
  
  Main_output <- list(Main_matrix=  NonW_ST_matrix_rivers,
                      STconmat=    NonW_ST_matrix_out_out,
                      STcon=       NonW_ST_connectivity_value
                      )
  
  return(Main_output)
  ####_______________________________________________________________________
  
}# Function

