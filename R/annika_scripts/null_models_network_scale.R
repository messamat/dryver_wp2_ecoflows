# DRYvER 

# Network-scale assembly processes during each campaign
# - How does network-scale dying affect DRNs assembly processes?


rm(list = ls())
library(data.table)

# ENV data

env <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/ENV_all_NAs_as_blank.csv") # running_id used for merging

# FUNGI - biofilm
# ------------

fungi_biof <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/fungi_dna_Biof.csv")
fungi_biof[1:10, 1:10]
# make same running_id
fungi_biof$running_id <- paste(fungi_biof$Site, fungi_biof$Campaign, sep = "_") 
fungi_biof[1:10, 3870:3873]
# combine state_of_flow
fungi_biof <- merge(fungi_biof, env[,.(running_id, state_of_flow)], by="running_id", all.x=T, sort=F) 
fungi_biof[1:10, 3870:3874]

# make response and explanatory tables
resp <- fungi_biof[,7:3873]
exp <- fungi_biof[,c(1:6,3874)]


#remove singletons i.e. taxa that are present in one sample only
# resp.nos <- data.frame(resp[,colSums(resp>0)>1,drop=FALSE])
# resp[1:10, 1:10]
# data_cleaned = resp[colSums(resp)> 2,]

colSums(resp)
resp.2 <- resp[, which(colSums(resp) != 0)]
