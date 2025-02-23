############################### PREP ############################################
# load libraries
library(tidyverse)
library(ncdf4)
library(lubridate)
library(reshape2)
library(sf)
#library(devtools)
library(igraph)
#library(riverconn)

working_directory <- file.path(getwd(), "bin", "STcon_DRYvER")
#setwd(working_directory)

# load functions
source("functions_flow_intermittence_indicators.R")
# The code from GitHUB includes plots and network centrality calculations
#source("https://raw.github.com/Cunillera-Montcusi/Quantifyinig-SpaTem-connectivity/main/SpaTemp_function.R")
# You are interested in the "depurated" version of the code
source("SpaTemp_function_Mathis.R")
source("SpaTemp_function_Mathis_edit.R")

DryRiverNet <- c("Albarine", "Bukkosdi", "Lepsamaanjoki", "Genal","Butiznica","Velicka")

DRN_position=5 # For looping between catchments 

#_____________________________________________________________________________________________________________________
# Import data following WP1 codes and functions to upload the flow_intermittence 
catchment <- DryRiverNet[DRN_position] # select the catchment name 
nc_file <- get_results_file_present(catchment)
nc <- nc_open(nc_file) # open netcdf file
reachID <- ncvar_get(nc, "reachID") # get list of reaches IDs
dates <- ncvar_get(nc, "date") # get dates of simulation period
dates <- as.Date(dates, origin="1950-01-01") # convert dates into R date format
flow_intermittence <- get_nc_var_present(nc, "isflowing", reachID, dates) # 0=dry, 1=flowing

#___________________________________________________________________
# WARNING 1
# Some flow intermittence values for some dates are NA, it is not necessarily for all dates but it is something 
# to take into account. I assign them 1 in case that they exist. 
flow_intermittence$value[is.na(flow_intermittence$value)] <- 1 
#__________________________________________________________________

# Shape names & uploadign
Shape_country <- c("france","hungary","finland","spain","croatia","czech")
shape_DRN <- sf::st_read(paste("D:\\DRYvER_wp2\\dryver_wp2_ecoflows\\results\\gis\\",
                               Shape_country[DRN_position],
                               "_river_network_nocomplexconf20250222_reided20250222.shp",
                               sep=""))

#___________________________________________________________________
# WARNING 2
# Croatian stream has two weird reaches.
# Reach 9001 is actually not going anywhere, its from and to are not connected to anything. I eliminate it. 
# Reach 1176 does not exist in the flow intermittence and there is no "1100" cat for the entire catchment. I suspect, 
#that it is a "typo" of some kind. The from and to from this reach connect to reaches with cat 1776, which made me think
#on this as a typo. I edited that reach to 1776 to get the corresponding assignation. 
# if(DRN_position==5){
#   shape_DRN <- shape_DRN %>% filter(cat!=9001)
#   shape_DRN[which(shape_DRN$cat==1176),4] <- 1776 }
#___________________________________________________________________

#___________________________________________________________________
# WARNING 3
# Some reaches have the same "cat" so I change their name adding a LETTER for the duplicated ones. 
#this is just ot later assign the outlet as igraph does not want duplicated reaches.
cat_to_correct <- shape_DRN %>% arrange(from) %>% group_by(cat) %>%
  mutate(count=n()) %>%  filter(count>1)%>% pull(cat)
if (length(cat_to_correct)>0) {
  shape_DRN[which(shape_DRN$cat%in%cat_to_correct),4]$cat <-
    paste(cat_to_correct,rep(LETTERS[1:length(which(cat_to_correct-cat_to_correct[1]==0))]),sep="_")
}
#___________________________________________________________________

# We built the adjacency list to built the graph
DRN_adj_table <- shape_DRN  %>% 
  arrange(from) %>% 
  select(from, to, length_m, to_cat_shp) %>%
  sf::st_drop_geometry()

# We identify where the outlet is (the NA and to_cat_shp)
Outlet_from_to <- DRN_adj_table[which(is.na(DRN_adj_table$to_cat_shp)==T),c("to")] %>%
  sf::st_drop_geometry() %>%
  as.numeric

# We create the graph from the DRN_adj_list (from - to)
igr1 <- igraph::graph_from_edgelist(as.matrix(DRN_adj_table[1:2]), directed = TRUE)

# We need to incorporate a new vertex corresponding the outlet (that we have detecte in the previous step). 
# We name the cat 11111 and it will be our last vertex and where the outlet will be directed.  
nodes_cat <- shape_DRN %>% arrange(from) %>% select(cat,from) %>% 
  sf::st_drop_geometry() %>% 
  mutate(cat=as.character(cat)) %>% 
  bind_rows(data.frame("cat"="11111",from=Outlet_from_to)) %>% arrange(from)

# We assign vertex attributes according to the cat_shape, ordered following the "from" value from the cat_shape, 
# which is not exactly "from anymore" it is an ID of the vertex
V(igr1)$cat <- as.character(nodes_cat$cat)

#Compute edge distances for distance matrices
E(igr1)$weight <- DRN_adj_table$length_m

# Using riverconn package (Damiano is the best) we set the directionality towards 11111 vertex. Our outlet.
# This is a key step to ensure that the network ends where it needs to end! Otherwise MISTAKE! 
#igr1 <- set_graph_directionality(igr1, field_name = "cat",outlet_name =  "11111")

# To ease the plotting and posterior management we create an "edge_DaFr" where we put the coordinates of the  
# starting and ending points of the segments and we edit the x and y coordinates. 
edges_DaFr <- left_join(
  as.data.frame(st_line_sample(shape_DRN, sample = 0) %>% st_coordinates())%>% 
    left_join(shape_DRN %>% select(UID,from,to) %>% st_drop_geometry(),by=c("L1"="UID")) %>% 
    rename(X1_coord=X,Y1_coord=Y),
  as.data.frame(st_line_sample(shape_DRN, sample = 1) %>% st_coordinates()) %>% 
    left_join(shape_DRN %>% select(UID,from,to) %>% st_drop_geometry(),by=c("L1"="UID")) %>% 
    rename(X2_coord=X,Y2_coord=Y),
  by=c("from","to"))

# Second step is to built a data.frame with the starting coordinates of all the segments and add 
# the initial row of the ending coordinates (which corresponds to the last point of the graph so the endpoint). 
nodes_df <- as.data.frame(st_line_sample(shape_DRN, sample = 0) %>% st_coordinates()) %>% 
  bind_rows(as.data.frame(st_line_sample(shape_DRN, sample = 1) %>% st_coordinates())[1,]) %>% 
  select(L1,"x"=X,"y"=Y) %>% mutate(L1_Copy=L1,.after=L1)

igraph_p <- ggplot()+
  geom_sf(data=shape_DRN, color='blue') +
  geom_segment(data=edges_DaFr,
               aes(x=X1_coord,y=Y1_coord, xend=X2_coord, yend=Y2_coord),
               arrow =arrow(length=unit(0.4,"cm"), ends="last"), linewidth=0.5, 
               colour="grey50", alpha=1)+
  geom_point(data=nodes_df, aes(x=x, y=y, text=L1), shape=21)+
  theme_classic()
plotly::ggplotly(igraph_p, tooltip = "text")


# We create the "End_point" site that will correspond to the 111111 in the flow_intermittence dataset
# this point will be added with the same frequency of any other reach
End_point <- flow_intermittence %>% 
  filter(reachID==flow_intermittence$reachID[1]) %>%# We select whatever reach
  rename("cat"="reachID") %>% 
  mutate(cat="11111",value=1) %>% #Modify its values as we wish 
  select(dates, cat, value) # select only this three values to merge them later

# Now last step to "merge" everything together into the Intermittence dataset
# We merge the shape_DRN dataset with the flow intermittence throw splitting the repeated "cat"
# and making sure that we are assigning the corresponding value of intermittence to the reach and 
# at the same time we maintain the "cat ids" that were used to built the graph. 
# Two options exists because not all DRN have repeated values
# We merge the Endpoint site and we "pivot_wide" the table to obtain the TRUE intermittence table, 
# where each row corresponds to a day (dates as factors) and columns to all nodes of the network.
if (length(cat_to_correct)>0) {
  interm_dataset <- 
    shape_DRN %>% select(UID,from,to,cat) %>% st_drop_geometry() %>% 
    mutate(TRUE_cat=unlist(lapply(strsplit(cat, split = "_"), `[`, 1))) %>% 
    left_join(flow_intermittence,by=c("TRUE_cat"="reachID"),relationship = "many-to-many") %>%
    select(dates, cat, value) %>% 
    bind_rows(End_point) %>% 
    pivot_wider(id_cols = dates,names_from = cat,values_from = value) %>% 
    mutate(dates=as.factor(dates))
}else{
  interm_dataset <- 
    shape_DRN %>% select(UID,from,to,cat) %>% st_drop_geometry() %>% 
    mutate(cat=as.character(cat)) %>% 
    left_join(flow_intermittence,by=c("cat"="reachID"),relationship = "many-to-many")%>%
    select(dates, cat, value) %>% 
    bind_rows(End_point) %>% 
    pivot_wider(id_cols = dates,names_from = cat,values_from = value) %>% 
    mutate(dates=as.factor(dates))
}

# CONGRATS! All is ready to STcon things! 

# You cut the length of the intermittence according to whatever you want. In this case we select the first 
# 30 days of the whole dataset. 
FL_intermitence_cut <- as.data.frame(interm_dataset[1:30,]) # THIS IS THE MOST IMPORTANT POINT! WHERE YOU DEFINE THE TIME WINDOW!!!! 

# We built the matrix of the network structure for the STcon, which is the "base" on which connectivity will be assessed. 
Network_structure <- as.data.frame(as.matrix(as_adjacency_matrix(igr1)))

# Last comprobation with some indicators. FL_intermittence has to have 1 more (dates) column
cat("There are", nrow(nodes_df), "sites.",
    ncol(FL_intermitence_cut), "columns in the intermitence dataset &",
    dim(Network_structure), "network rows and columns.")

# We put everything in lists as the funcion parallelizes things according to the number of things inside 
# the number of objects inside the lists! So, everything must be entered as a list! 
# For example: 
# In case of calculating STcon for a DRN in different dates you could create a list with 
# the DRN information into lists where each element of the list corresopnds to a date or whatever. 
# In this examples we just have 1 element per list.
interm_dataset_campaigns_To_Run <- list(FL_intermitence_cut)
Sites_coordinates_campaigns_To_Run <- list(nodes_df)
Network_stru_campaigns_To_Run <- list(Network_structure)

# We add an extra "reference river" into the lists with a REFerence intermitence where all sites are permanent (so = 1)
REF_FL_intermitence<-FL_intermitence_cut %>% replace(.==0,1)
# We add the seventh "reference river into the list
interm_dataset_campaigns_To_Run[[
  length(interm_dataset_campaigns_To_Run)+1]] <- REF_FL_intermitence
Sites_coordinates_campaigns_To_Run[[
  length(Sites_coordinates_campaigns_To_Run)+1]] <- nodes_df
Network_stru_campaigns_To_Run[[
  length(Network_stru_campaigns_To_Run)+1]] <- Network_structure

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
river_dist_mat <- igraph::distances(igr1)

interm_dataset = interm_dataset_campaigns_To_Run[[1]]
Sites_coordinates=Sites_coordinates_campaigns_To_Run[[1]]
Network_stru = Network_stru_campaigns_To_Run[[1]]
direction="directed"
sense= "out"
weighting=T
dist_matrices = river_dist_mat # Weighting pairs
weighting_links = F
link_weights = NULL # Weighting links
legacy_effect = 1
legacy_length = 1 # Legacy effects
indirect_dispersal = TRUE
value_s_link=1
value_t_link=1 # Values to links
value_no_s_link=0.1
value_no_t_link=0.1 #
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################
####_______________________________________________________________________
# River network ####
####_______________________________________________________________________
ST_matrix_rivers_edit <- list()
out_Matrix_LIST_edit <- list()

# We calculate the number of nodes of our network (used along the function)  
interm_ncols_edit <- ncol(interm_dataset)
numn_nodes_edit <- interm_ncols_edit - 1
nsteps_edit <- length(interm_dataset[,1])

if(weighting==TRUE){dist_matr <- dist_matrices}

# We built the matrix corresponding to the num. of nodes multiplied by the DAYS of HOBOS that we have
### This matrix is the "giant" template where we will put all the values.
ST_matrix_edit <- matrix(nrow = numn_nodes_edit, ncol = numn_nodes_edit*2, 
                         data = value_no_s_link)
ST_matrix_netwGraph_edit <- matrix(nrow = numn_nodes_edit, ncol = numn_nodes_edit, 
                                   data = 0)

# Once created the template we start to fill it for every day
### We fill it for Days (or time)-1 because the last day does not have a "future" from which to extract values. 
for (days in 1:(nsteps-1)) {
  cat("We are at time unit", days, "of", (nsteps-1), "\n")
  # First we define the spatial connections of the matrix
  ### Also known as the rows or columns at which we have to add the values of the connections 
  spa_connections_edit <- seq(1, numn_nodes_edit)#+((days-1)*numn_nodes) 
  
  # We obtain the time steps:
  ## time_step_1 is the present
  ## time_step_2 is the following step (the close future)
  time_step_1_edit <- interm_dataset[days, 2:interm_ncols_edit] 
  time_step_2_edit <- interm_dataset[days+1, 2:interm_ncols_edit] 
  
  #Simple fluvial network_______________________
  #Create an adjacancy matrix for time step 1 whereby:
  #for sites that are wet, get the normal structure (direct connection to sites)
  #for sites that are dry, 0s to all sites
  ST_matrix_netwGraph_edit[time_step_1_edit==1,] <- 
    as.matrix(Network_stru)[time_step_1_edit==1,]
  diag_backup_edit <- diag(ST_matrix_netwGraph_edit)
  ST_matrix_netwGraph_edit[time_step_1_edit==0,] <- 0
  diag(ST_matrix_netwGraph_edit) <- diag_backup_edit #####REMARK Mathis: not sure why this is needed
  #It's equivalent to [-site_step] in the original code
  #ST_matrix_netwGraph[spa_connections[site_step],
  #c(spa_connections[1]:spa_connections[numn_nodes])[-site_step]] <- 0
  
  # FLuvial SPATIAL links ___________________________________________________________________________________________________________________
  ## Here we fill the matrix section corresponding to the time_step based on the river graph based on a dendritic. 
  # We create the graph
  a_edit <- igraph::graph_from_adjacency_matrix(ST_matrix_netwGraph_edit, 
                                                mode=direction, 
                                                diag = FALSE)
  
  # Compute shortest path distances for all node pairs
  dist_matrix <- igraph::distances(a_edit, mode = sense)
  dist_matrix[is.infinite(dist_matrix)] <- 0
  
  # Convert distances into binary connectivity (1 if connected, 0 if not)
  All_river_paths_edit <- ifelse(dist_matrix > 0, value_s_link, value_no_s_link)
  
  # We weight the sites for the distances between them (a pairwise matrix)
  if (weighting == TRUE) {
    All_river_paths_edit <- All_river_paths_edit * dist_matr
  }
  
  # We add the "All_river_paths" filled for each node in the "big" matrix specific sites
  ST_matrix_edit[, spa_connections_edit] <- ST_matrix_edit[, spa_connections_edit] +
    All_river_paths_edit
  
  # # In the following lines we continue the party towards temporal steps
  # We create the matrix where we will drop the information of the shortest paths.
  All_river_paths_edit <- matrix(nrow = length(time_step_1_edit),
                                 ncol = length(time_step_1_edit),
                                 data = value_no_t_link)
  
  temp_connections_edit <- seq(1 + numn_nodes_edit, 2 * numn_nodes_edit)
  
  #Calculate temporal changes in one step
  temp_changes_edit <- time_step_1_edit - time_step_2_edit
  
  #Determine stable (can be stable 1-1 or 0-0), lost, and gained indices
  stable_indices_1 <- which(temp_changes_edit == 0 & time_step_1_edit == 1)
  stable_indices_0 <- which(temp_changes_edit == 0 & time_step_1_edit == 0)
  lost_indices <- which(temp_changes_edit == 1)
  gained_indices <- which(temp_changes_edit == -1)
  
  #Compute direct spatial links for stable connected nodes (1->1)
  All_river_paths_edit[stable_indices_1, ] <- ifelse(
    dist_matrix[stable_indices_1, ] > 0,
    value_t_link, value_no_t_link)
  
  #Indirect dispersal for lost nodes (1->0) 
  if (indirect_dispersal) {
    All_river_paths_edit[lost_indices, ] <- 
      All_river_paths_edit[lost_indices, ] + 
      ifelse(dist_matrix[lost_indices, ] > 0, 
             value_t_link, value_no_t_link)
  }

  #Apply weights 
  if (weighting) {All_river_paths_edit <- All_river_paths_edit * dist_matr}
  
  #Apply legacy effects and update ST_matrix
  for (leg_eff in seq_len(legacy_length)) {
    All_river_paths_legacy_edit <- All_river_paths_edit * legacy_effect[leg_eff]
    ST_matrix_edit[, temp_connections_edit] <- 
      ST_matrix_edit[, temp_connections_edit] + All_river_paths_legacy_edit
    
    #Self-connection for stable connected nodes (diagonal assignment)
    value_t_link_modif_edit <- if (weighting_links) value_t_link * day_link_weights else value_t_link
    ST_matrix_edit[cbind(spa_connections_edit[stable_indices_1], 
                         temp_connections_edit[stable_indices_1])] <- 
      ST_matrix_edit[cbind(spa_connections_edit[stable_indices_1], 
                           temp_connections_edit[stable_indices_1])] + value_t_link_modif_edit
  }
}# Days closing

out_Matrix_edit <- list(ST_matrix_edit)
out_Matrix_LIST_edit <- out_Matrix_edit

# Extracting the results into different lists
ST_matrix_rivers_edit <- out_Matrix_LIST_edit[[1]]

####_______________________________________________________________________
# STconmat calculation ####
####_______________________________________________________________________
# Find below the lines to calculate the "collapsing" matrix that just sums all the values of all the SPATIOTEMPORAL matrix
# These pairwise matrix is called the STconmat.   
ST_matrix_out_out_edit <- list()

# We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
out_out_edit <- matrix(nrow = numn_nodes_edit,ncol = numn_nodes_edit, 
                       data = value_no_s_link)

Spatial_matrix_edit <-  ST_matrix_rivers_edit[,spa_connections_edit]
Temporal_matrix_edit <- ST_matrix_rivers_edit[,temp_connections_edit]

out_out_edit <- Spatial_matrix_edit+Temporal_matrix_edit
out_out_edit <- out_out_edit/c(nsteps-1)

# We save the collapsed matrix
ST_matrix_out_out_edit <- out_out_edit

####_______________________________________________________________________

####_______________________________________________________________________
# STcon calculation ####
####_______________________________________________________________________
ST_connectivity_value <- list()

# We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, 
                  data = value_no_s_link)

Spatial_matrix <-  ST_matrix_rivers[,spa_connections]
Temporal_matrix <- ST_matrix_rivers[,temp_connections]

out_out <- Spatial_matrix+Temporal_matrix

# "leng_correct" is a reverse vector (from big to small) used to correct the fact that uperstream nodes will have higher values when 
# considering its number of connections. As I am "node 1" my number of connections will be higher tan "node 10". IF WE FOLLOW THE RIVER DOWNSTREAM!
leng_correct <- c()
aa <- graph_from_adjacency_matrix(as.matrix(Network_stru),mode = "directed")
for (neigg in 1:numn_nodes) {
  neighbour<- all_shortest_paths(aa, from = neigg, to = c(1:numn_nodes)[-neigg],mode = sense)$nrgeo
  neighbour[neigg] <- 0
  leng_correct[neigg] <- length(which(neighbour>0))
}

spt_conn <-apply(out_out,1,sum)/leng_correct
# We divide by the number of days so we obtain the "per day" values
spt_conn<- spt_conn/c(nsteps-1)

ST_connectivity_value <- spt_conn  



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


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

ST_matrix_rivers <- list()
ST_directed_Ocloseness_rivers <- list()
ST_directed_Allcloseness_rivers <- list()
ST_directed_betweennes_rivers <- list()

#Parallelization parameters
# cores <- detectCores() #Number of cores in computer
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)

out_Matrix_LIST <- list()
pack_check <- search()
pack_check_val <- length(which(pack_check=="package:sna"))
if(pack_check_val>0){detach("package:sna", unload = TRUE)}
#for (river in 1:length(interm_dataset)) { #- With this it takes 6'26''
# We calculate the number of nodes of our network (used along the function)  
numn_nodes <- ncol(interm_dataset)-1

if(weighting==TRUE){dist_matr <- dist_matrices}

# We built the matrix corresponding to the num. of nodes multiplied by the DAYS of HOBOS that we have
### This matrix is the "giant" themplate where we will put all the values.
ST_matrix <- matrix(nrow = numn_nodes,ncol = numn_nodes*2, data=value_no_s_link)
ST_matrix_netwGraph <- matrix(nrow = numn_nodes,ncol = numn_nodes, data=0)

# Once created the template we start to fill it for every day
### We fill it for Days (or time)-1 because the last day does not have a "future" from which to extract values. 
for (days in 1:c(length(interm_dataset[,1])-1)) {
  cat("We are at time unit", days, "of", (length(interm_dataset[,1])-1), "and at river", river,"\n")
  # First we define the spatial connections of the matrix
  ### Also known as the rows or columns at which we have to add the values of the connections 
  spa_connections <-seq(1,length(colnames(interm_dataset))-1,1)#+((days-1)*numn_nodes)
  
  # We obtain the time steps:
  ## time_step_1 is the present
  ## time_step_2 is the following step (the close future)
  time_step_1 <- interm_dataset[days,2:ncol(interm_dataset)]
  time_step_2 <- interm_dataset[days+1,2:ncol(interm_dataset)]
  
  #Simple fluvial network_______________________
  ## This step fills "the diagonal" of each time_step following the direction of the river
  ## it basically connects the river in a dendritic structure.
  for (site_step in 1:c(length(time_step_1))) {
    #print(site_step)
    if(time_step_1[site_step]==1){
      ST_matrix_netwGraph[spa_connections[site_step],
                          c(spa_connections[1]:spa_connections[numn_nodes])] <- as.numeric(Network_stru[site_step,])
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
  All_river_paths <- matrix(nrow =length(time_step_1),ncol = length(time_step_1),
                            data = value_no_s_link)
  #All_river_paths[upper.tri(All_river_paths)] <- value_no_s_link
  
  # For each path (e.g., from node 1 to node 7) we "check" the length of the shortest path. 
  ## check = 0 means that the graph is disconnected.
  ## check bigger than 0 means that the graph is connected.
  for (every_path in 1:c(length(time_step_1))){
    check <- length(all_shortest_paths(a, every_path, 1:numn_nodes, mode = sense)$res)
    if (check==0) {
      site <-0
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
    # We weight the sites for the distances between them (a pairwise matrix)
    if(weighting==T){All_river_paths[every_path,] <-  
      as.numeric(All_river_paths[every_path,]*dist_matr[every_path,])}
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
    temp_connections <-seq(1+numn_nodes,length(colnames(interm_dataset))-1+numn_nodes,1)#+((days)*numn_nodes) 
    
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
      All_river_paths[site_step, site] <- value_t_link
      # We weight the links base on daily information of flow or strength of the link.
      if(weighting_links==T){All_river_paths[site_step,] <- 
        as.numeric(All_river_paths[site_step,]*as.numeric(day_link_weights[site_step]))}
      # We weight
      if(weighting==T){All_river_paths[site_step,] <-  
        as.numeric(All_river_paths[site_step,]*dist_matr[site_step,])}
      # We pass it to the main matrix
      for (leg_eff in 1:legacy_length) {
        All_river_paths_legacy <- All_river_paths[site_step,]*legacy_effect[leg_eff]
        ST_matrix[spa_connections[site_step],
                  temp_connections[1]:temp_connections[numn_nodes]] <- All_river_paths_legacy+ST_matrix[spa_connections[site_step],
                                                                                                        temp_connections[1]:temp_connections[numn_nodes]]
      }
      
    }# End of the if
    
  }# Site_step closing
  
  #############################################################################################
  #CHECK CORRESPONDENCE
  # which(!sapply(seq_len(nrow(ST_matrix)), 
  #               function(i) all(round(ST_matrix_edit[i,], 5) == 
  #                                 round(ST_matrix[i,], 5))))
  ###################################################################################################
}# Days closing

out_Matrix <- list(ST_matrix,
                   ST_Oclosenness_matrix,
                   ST_Allclosenness_matrix,
                   ST_betweennes_matrix)
out_Matrix_LIST <- out_Matrix


#parallel::stopCluster(cl) #Stop parallel computing

# Exctracring the results into different lists
ST_matrix_rivers <- out_Matrix_LIST[[1]]

####_______________________________________________________________________
# STconmat calculaiton ####
####_______________________________________________________________________
# Find below the lines to calculate the "collapsing" matrix that just summs all the values of all the SPATIOTEMPORAL matrix
# These pairwise matrix is called the STconmat.   
ST_matrix_out_out <- list()


numn_nodes <- ncol(interm_dataset)-1
temp_connections <-seq(1+numn_nodes,length(colnames(interm_dataset))-1+numn_nodes,1)
spa_connections <-seq(1,length(colnames(interm_dataset))-1,1)

# We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)

Spatial_matrix <-  ST_matrix_rivers[,spa_connections]
Temporal_matrix <- ST_matrix_rivers[,temp_connections]

out_out <- Spatial_matrix+Temporal_matrix
out_out <- out_out/c(length(interm_dataset[,1])-1)

# We save the collapsed matrix
ST_matrix_out_out <- out_out

####_______________________________________________________________________

####_______________________________________________________________________
# STcon calculaiton ####
####_______________________________________________________________________
ST_connectivity_value <- list()

# We already know this value
numn_nodes <- ncol(interm_dataset)-1
temp_connections <-seq(1+numn_nodes,length(colnames(interm_dataset))-1+numn_nodes,1)
spa_connections <-seq(1,length(colnames(interm_dataset))-1,1)

# We create the out matrix which match the size of our "simple" matrix num_nodes*num_nodes
out_out <- matrix(nrow = numn_nodes,ncol = numn_nodes, data = 0)

Spatial_matrix <-  ST_matrix_rivers[,spa_connections]
Temporal_matrix <- ST_matrix_rivers[,temp_connections]

out_out <- Spatial_matrix+Temporal_matrix

# "leng_correct" is a reverse vector (from big to small) used to correct the fact that uperstream nodes will have higher values when 
# considering its number of connections. As I am "node 1" my number of connections will be higher tan "node 10". IF WE FOLLOW THE RIVER DOWNSTREAM!
leng_correct <- c()
aa <- graph_from_adjacency_matrix(as.matrix(Network_stru),mode = "directed")
for (neigg in 1:numn_nodes) {
  neighbour<- all_shortest_paths(aa, from = neigg, to = c(1:numn_nodes)[-neigg],mode = sense)$nrgeo
  neighbour[neigg] <- 0
  leng_correct[neigg] <- length(which(neighbour>0))
}

spt_conn <-apply(out_out,1,sum)/leng_correct
# We divide by the number of days so we obtain the "per day" values
spt_conn<- spt_conn/c(length(interm_dataset[,1])-1)

ST_connectivity_value <- spt_conn  



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

##################################################

