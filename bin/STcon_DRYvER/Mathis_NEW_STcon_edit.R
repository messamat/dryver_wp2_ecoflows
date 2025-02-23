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
FL_intermitence_cut <- as.data.frame(interm_dataset[1:50,]) # THIS IS THE MOST IMPORTANT POINT! WHERE YOU DEFINE THE TIME WINDOW!!!! 

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

interm_ncols <- ncol(interm_dataset_campaigns_To_Run[[1]])


############################## TEST NO WEIGHTS #################################
# Here we STCON things! 
interm_dataset = interm_dataset_campaigns_To_Run[1]
Sites_coordinates=Sites_coordinates_campaigns_To_Run[1]
Network_stru = Network_stru_campaigns_To_Run[1]
direction="directed"
sense= "out"
weighting=F
dist_matrices = NULL # Weighting pairs
weighting_links =F
link_weights = NULL # Weighting links
legacy_effect = 1
legacy_length = 1 # Legacy effects
value_s_link=1
value_t_link=1 # Values to links
value_no_s_link=0
value_no_t_link=0 #

tictoc::tic()
DirNonW_original <- spat_temp_index(
  interm_dataset = interm_dataset_campaigns_To_Run[1],
  Sites_coordinates=Sites_coordinates_campaigns_To_Run[1],
  Network_stru = Network_stru_campaigns_To_Run[1], 
  direction="directed", 
  sense= "out",
  weighting=F,
  dist_matrices = NULL, # Weighting pairs
  weighting_links =F,
  link_weights = NULL, # Weighting links
  legacy_effect = 1, 
  legacy_length = 1, # Legacy effects
  value_s_link=1,
  value_t_link=1, # Values to links
  value_no_s_link=0,
  value_no_t_link=0 # Values to links
) # Last parameters information
tictoc::toc()

tictoc::tic()
DirNonW <- spat_temp_index_edit(
  interm_dataset = as.matrix(interm_dataset_campaigns_To_Run[[1]][, 2:interm_ncols]),
  Network_stru = as.matrix(Network_stru_campaigns_To_Run[[1]]), 
  direction="directed", 
  sense= "out",
  weighting= F,
  dist_matrices = NULL, # Weighting pairs
  weighting_links =F,
  link_weights = NULL, # Weighting links
  legacy_effect = 1, 
  legacy_length = 1, # Legacy effects
  value_s_link=1,
  value_t_link=1, # Values to links
  value_no_s_link=0,
  value_no_t_link=0 # Values to links
) ast parameters information
tictoc::toc()

all(DirNonW_original$STcon[[1]] == DirNonW$STcon)

STcon <- DirNonW$STcon
#[[1]]/DirNonW$STcon[[2]]
STcon <- STcon[!is.na(STcon)]

ggplot()+
  geom_segment(data=edges_DaFr %>%mutate(STcon=STcon), 
               aes(x=X1_coord,y=Y1_coord, xend=X2_coord, yend=Y2_coord,colour=STcon), 
               arrow =arrow(length=unit(0.01,"cm"), ends="last"), linewidth=1.3, alpha=1)+
  geom_point(data=nodes_df,aes(x=x, y=y),fill="grey",shape=21,alpha=0.5)+
  scale_color_viridis()+
  theme_classic()


############################## TEST DISTANCE WEIGHTS ###########################
river_dist_mat <- igraph::distances(igr1)

tictoc::tic()
DirNonW_original <- spat_temp_index(
  interm_dataset = interm_dataset_campaigns_To_Run[1],
  Sites_coordinates=Sites_coordinates_campaigns_To_Run[1],
  Network_stru = Network_stru_campaigns_To_Run[1], 
  direction="directed", 
  sense= "out",
  weighting= T,
  dist_matrices = list(river_dist_mat), # Weighting pairs
  weighting_links =F,
  link_weights = NULL, # Weighting links
  legacy_effect = 1, 
  legacy_length = 1, # Legacy effects
  value_s_link=1,
  value_t_link=1, # Values to links
  value_no_s_link=0.1,
  value_no_t_link=0.1 # Values to links
) # Last parameters information
tictoc::toc()

# interm_dataset = as.matrix(interm_dataset_campaigns_To_Run[[1]][, 2:interm_ncols])
# Sites_coordinates=Sites_coordinates_campaigns_To_Run[[1]]
# Network_stru = Network_stru_campaigns_To_Run[[1]]
# direction="directed"
# sense= "out"
# weighting=T
# dist_matrices = river_dist_mat # Weighting pairs
# weighting_links = F
# link_weights = NULL # Weighting links
# legacy_effect = 1
# legacy_length = 1 # Legacy effects
# indirect_dispersal = TRUE
# value_s_link=1
# value_t_link=1 # Values to links
# value_no_s_link=0.1
# value_no_t_link=0.1 #

#profvis({
tictoc::tic()
DirNonW <- spat_temp_index_edit(
  interm_dataset = as.matrix(interm_dataset_campaigns_To_Run[[1]][, 2:interm_ncols]),
  Network_stru = as.matrix(Network_stru_campaigns_To_Run[[1]]), 
  direction="directed", 
  sense= "out",
  weighting= T,
  dist_matrices = river_dist_mat, # Weighting pairs
  weighting_links =F,
  link_weights = NULL, # Weighting links
  legacy_effect = 1, 
  legacy_length = 1, # Legacy effects
  value_s_link=1,
  value_t_link=1, # Values to links
  value_no_s_link=0.1,
  value_no_t_link=0.1 # Values to links
) # Last
tictoc::toc()
#   interval = 0.005
# })

all(round(DirNonW_original$STcon[[1]], 3) == round(DirNonW$STcon, 3))
all(round(DirNonW_original$STcon[[1]]/DirNonW$STcon, 3) == 1, na.rm=T)

