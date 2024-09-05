# Temporal mean alpha diversity calculations

rm(list = ls())
library(data.table)
library(vegan)

#### Sediment diatoms ####
dia_sed <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/diatoms_dna_sediment_removed_pools_and_zero_rows_and_columns.csv")
dia_sed[1:10, 1:10]

com <- dia_sed[,c(3,8:1301)]
meta <- dia_sed[,1:7]
# first column as row names
com <- data.frame(com[,-1], row.names=com$running_id)

# diversity indices for running_id
dia_sed_S <- data.frame(specnumber(com))
dia_sed_S <- cbind(dia_sed_S, meta)
# mean diversity for each site
dia_sed_S_mean<- dia_sed_S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
dia_sed_S_mean$Organism <- "Sediment_diatoms"

rm(dia_sed, dia_sed_S)
rm(com2)

#### Biofilm diatoms ####
dia_biof <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/diatoms_dna_biofilm_removed_pools_and_zero_rows_and_columns.csv")
dia_biof[1:10, 1:10]

dia_biof$running_id <- paste(dia_biof$Site, dia_biof$Campaign, sep = "_")

dia_biof[1:10, 1:10]
com <- dia_biof[,6:1237]
meta <- dia_biof[,c(1:5, 1237)]
# first/last column as row names
com <- data.frame(com[,-1232], row.names=com$running_id)
com[1:5, 1230:1231]

# diversity indices for running_id
dia_biof_S <- data.frame(specnumber(com))
dia_biof_S <- cbind(dia_biof_S, meta)
# mean diversity for each site
dia_biof_S_mean<- dia_biof_S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
dia_biof_S_mean$Organism <- "Biofilm_diatoms"

rm(dia_biof, dia_biof_S)
#### Sediment fungi ####
data <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/fungi_dna_sediment_removed_pools_zero_rows_and_columns.csv")
data[1:10, 1:10]

data$running_id <- paste(data$Site, data$Campaign, sep = "_")

data[1:5, 1:7]
com <- data[,6:3064]
meta <- data[,c(1:5, 3064)]

# first/last column as row names
com <- data.frame(com[,-3059], row.names=com$running_id)
com[1:5, 3057:3058]

# diversity indices for running_id
fun_sed_S <- data.frame(specnumber(com))
fun_sed_S <- cbind(fun_sed_S, meta)
# mean diversity for each site
fun_sed_S_mean<- fun_sed_S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
fun_sed_S_mean$Organism <- "Sediment_fungi"

rm(fun_sed_S)
rm(S)

#### Biofilm fungi ####
data <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/fungi_dna_Biof_removed_pools_zero_rows_and_columns.csv")
data[1:5, 1:8]

data$running_id <- paste(data$Site, data$Campaign, sep = "_")

com <- data[,6:1983]
meta <- data[,c(1:5, 1983)]

# first/last column as row names
com <- data.frame(com[,-1978], row.names=com$running_id)

# diversity indices for running_id
S <- data.frame(specnumber(com))
S <- cbind(S, meta)
# mean diversity for each site
fun_biof_S_mean<- S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
fun_biof_S_mean$Organism <- "Biofilm_fungi"

rm(S)

#### Macroinvertebrates, with pools ####

data <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/DRYvER_MIV_data_genus_reduced_taxa.csv")
data[1:5, 1:8]

data$running_id <- paste(data$Site, data$Campaign, sep = "_")

com <- data[,6:300]
meta <- data[,c(1:5, 300)]

# first/last column as row names
com <- data.frame(com[,-295], row.names=com$running_id)

# diversity indices for running_id
S <- data.frame(specnumber(com))
S <- cbind(S, meta)

# mean diversity for each site
mac_pools_S_mean<- S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
mac_pools_S_mean$Organism <- "Macroinvertebrates_pools"


#### Macroinvertebrates, with no pools ####

data <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/DRYvER_MIV_data_genus_reduced_taxa_removed_pools.csv")
data[1:5, 1:8]

data$running_id <- paste(data$Site, data$Campaign, sep = "_")

com <- data[,6:300]
meta <- data[,c(1:5, 300)]

# first/last column as row names
com <- data.frame(com[,-295], row.names=com$running_id)

# diversity indices for running_id
S <- data.frame(specnumber(com))
S <- cbind(S, meta)

# mean diversity for each site
mac_S_mean<- S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
mac_S_mean$Organism <- "Macroinvertebrates"


#### Sediment bacteria ####
data <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/bacteria_DNA_Soil.csv")
dnameta <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/metadata DRYvER eDNA updated.csv")
data[1:5, 1:8] 

data$running_id <- paste(data$ENV_ID, data$Campaign, sep = "_")
dnameta$running_id <- paste(dnameta$Site_code, dnameta$Campaign, sep="_")

# remove pool sites
data2 <- merge(data, dnameta[,.(Sample.ID, Habitat)], by="Sample.ID", all.x=T, sort=F)
data3 <- droplevels(data2[!data2$Habitat == 'pool',]) 

com1 <- data3[,8:36859]
com1[1:5, 36849:36851]

# first/last column as row names
com1 <- data.frame(com1[,-36851], row.names=com1$running_id)
com1$Habitat <- NULL

# removing zero species by transposing and removing zero rows
com1t <- t(com1)
com1t[1:5, 1:5]

tot <- rowSums(com1t) 
com2t <- com1t[tot > 0, ]

com3 <- t(com2t)
com3 <- as.data.frame(com3)
com3[1:5, 1:5] # final version, transposed back, with no zero species

# write.csv(com3, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Data/Sediment_bacteria_removed_pools_and_zero_species.csv")

# näyttää OK:lta
# test <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/Sediment_bacteria_removed_pools_and_zero_species.csv")

meta <- data3[,1:6]
# diversity indices for running_id
S <- data.frame(specnumber(com3))
S <- cbind(S, meta)
# mean diversity for each site
library(dplyr)
bac_sed_S_mean<- S%>%group_by(ENV_ID)%>%summarise(mean_S=mean(specnumber.com3., na.rm = T))
bac_sed_S_mean$Organism <- "Sediment_bacteria"

bac_sed_S_mean$Site <- bac_sed_S_mean$ENV_ID
bac_sed_S_mean$ENV_ID <- NULL


#### Biofilm bacteria ####

data <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/bacteria_DNA_Biof.csv")
dnameta <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/metadata DRYvER eDNA updated.csv")
data[1:5, 1:8] 

data$running_id <- paste(data$Site, data$Campaign, sep = "_")
dnameta$running_id <- paste(dnameta$Site_code, dnameta$Campaign, sep="_")

# remove pool sites
data2 <- merge(data, dnameta[,.(Sample.ID, Habitat)], by="Sample.ID", all.x=T, sort=F)
data3 <- droplevels(data2[!data2$Habitat == 'pool',]) 

com1 <- data3[,6:36857]
com1[1:5, 36849:36852]

# first/last column as row names
com1 <- data.frame(com1[,-36851], row.names=com1$running_id)
com1$Habitat <- NULL

# removing zero species by transposing and removing zero rows
com1t <- t(com1)
com1t[1:5, 1:5]

tot <- rowSums(com1t) 
com2t <- com1t[tot > 0, ]

com3 <- t(com2t)
com3 <- as.data.frame(com3)
com3[1:5, 1:5] # final version, transposed back, with no zero species

# write.csv(com3, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Data/Biofilm_bacteria_removed_pools_and_zero_species.csv")

# näyttää OK:lta
# test <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/Sediment_bacteria_removed_pools_and_zero_species.csv")

meta <- data3[,1:4]
# diversity indices for running_id
S <- data.frame(specnumber(com3))
S <- cbind(S, meta)
# mean diversity for each site
library(dplyr)
bac_biof_S_mean<- S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com3., na.rm = T))
bac_biof_S_mean$Organism <- "Biofilm_bacteria"



#### Flying macroinvertebrates, with no pools ####

data <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_flying.csv")
data[1:5, 1:8]

data$running_id <- paste(data$Site, data$Campaign, sep = "_")

com <- data[,6:215]
meta <- data[,c(1:5, 215)]
head(meta)
# first/last column as row names
com <- data.frame(com[,-210], row.names=com$running_id)
head(com)
# diversity indices for running_id
S <- data.frame(specnumber(com))
S <- cbind(S, meta)

# mean diversity for each site
flymac_S_mean<- S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
flymac_S_mean$Organism <- "Flying_macroinvertebrates"

#### Non-flying macroinvertebrates, with no pools ####

data <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/DRYvER_MIV_data_genus_reduced_taxa_removed_pools_for_traits_nonflying.csv")
data[1:5, 1:8]

data$running_id <- paste(data$Site, data$Campaign, sep = "_")

com <- data[,6:85]
meta <- data[,c(1:5, 85)]
head(meta)
# first/last column as row names
com <- data.frame(com[,-80], row.names=com$running_id)

# diversity indices for running_id
S <- data.frame(specnumber(com))
S <- cbind(S, meta)

# mean diversity for each site
nonflymac_S_mean<- S%>%group_by(Site)%>%summarise(mean_S=mean(specnumber.com., na.rm = T))
nonflymac_S_mean$Organism <- "Nonflying_macroinvertebrates"


#### Finally combine as rbind and save as csv ####

mean_S <- rbind(dia_biof_S_mean, dia_sed_S_mean, fun_biof_S_mean, fun_sed_S_mean, mac_pools_S_mean, mac_S_mean)

#save as csv
write.csv(mean_S, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity.csv")

# UPDATE 1: to update bacteria, open old diversity data
divdata <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity.csv")


names(bac_biof_S_mean)

names(divdata)
divdata$V1 <- NULL
names(divdata)

divupdate <- rbind(divdata, bac_biof_S_mean)

# reorder columns in bacteria sediment diversity dataframes
names(bac_sed_S_mean)
bac_sed_S_mean2 <- bac_sed_S_mean[,c(3,1,2)]
names(bac_sed_S_mean2)

divupdate <- rbind(divupdate, bac_sed_S_mean2)

#save as csv
write.csv(divupdate, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity_updated.csv")



# UPDATE 2, 15.4.2024: update flying and non-flying macroinvertebrates, with no pools data

divdata <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity_updated.csv")

names(flymac_S_mean)
names(nonflymac_S_mean)
names(divdata)
divdata$V1 <- NULL
names(divdata)


divupdate2 <- rbind(divdata, flymac_S_mean, nonflymac_S_mean)


#save as csv
write.csv(divupdate2, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity_updated2.csv", row.names = F)


