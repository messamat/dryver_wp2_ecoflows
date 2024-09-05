# DRYvER 

# Site-specific  assembly processes
# - How does site-specific dying affect local assembly processes? (cf. Huttunen et al. 2017)


rm(list = ls())
library(data.table)
library(vegan)

# BACTERIA - SEDIMENT 
# ------------

bac_sed <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/Sediment_bacteria_removed_pools_and_zero_species.csv")
bac_sed[1:5, 1:5]
bac_sed$running_id <- bac_sed$V1
bac_sed[1:5, 26246:26249] 
env <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/ENV_all_NAs_as_blank.csv")
names(env)
# open env data, combine country and site info with running id

bac_sed <- merge(bac_sed, env[,.(running_id, site, DRN)], by="running_id", all.x=T, sort=F)


levels(factor(bac_sed$DRN))

#### Croatia ####
com <- bac_sed[DRN == "Croatia"] 
com[1:5, 26247:26251]
tot <- rowSums(com[,3:26249]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!site %in% rare_levels)

groups <- as.factor(com$site)
# com[1:10, 1:10]
com <- com[,3:26249] ### MUOKATTAVA !!!!

start_time_total <- Sys.time()
### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Croatia <- oecosimu(com, foo, "quasiswap", nsimul = 99, groups = groups)
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/CRO_sediment_bac_null_models_99perm.txt") ### MUOKATTAVA !!!!
paste("-----")
paste("jaccard, null model")
output_Croatia
sink()

#### Czech ####
com <- bac_sed[DRN == "Czech"] 
com[1:5, 26247:26251]
tot <- rowSums(com[,3:26249]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!site %in% rare_levels)

groups <- as.factor(com$site)
# com[1:10, 1:10]
com <- com[,3:26249] ### MUOKATTAVA !!!!

start_time_total <- Sys.time()
### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Czech <- oecosimu(com, foo, "quasiswap", nsimul = 99, groups = groups) ### MUOKATTAVA !!!!
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/CZE_sediment_bac_null_models_99perm.txt") ### MUOKATTAVA !!!!
paste("-----")
paste("jaccard, null model")
output_Czech ### MUOKATTAVA !!!!
sink()

#### Finland ####
com <- bac_sed[DRN == "Finland"] 
com[1:5, 26247:26251]
tot <- rowSums(com[,3:26249]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!site %in% rare_levels)

groups <- as.factor(com$site)
# com[1:10, 1:10]
com <- com[,3:26249] ### MUOKATTAVA !!!!

start_time_total <- Sys.time()
### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Finland <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups) ### MUOKATTAVA !!!!
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/FIN_sediment_bac_null_models_999perm.txt") ### MUOKATTAVA !!!!
paste("-----")
paste("jaccard, null model")
output_Finland ### MUOKATTAVA !!!!
sink()

#### France ####
com <- bac_sed[DRN == "France"] 
com[1:5, 26247:26251]
tot <- rowSums(com[,3:26249]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!site %in% rare_levels)

groups <- as.factor(com$site)
# com[1:10, 1:10]
com <- com[,3:26249] ### MUOKATTAVA !!!!

start_time_total <- Sys.time()
### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_France <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups) ### MUOKATTAVA !!!!
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/FRA_sediment_bac_null_models_999perm.txt") ### MUOKATTAVA !!!!
paste("-----")
paste("jaccard, null model")
output_France ### MUOKATTAVA !!!!
sink()
#### Hungary ####
com <- bac_sed[DRN == "Hungary"] 
com[1:5, 26247:26251]
tot <- rowSums(com[,3:26249]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!site %in% rare_levels)

groups <- as.factor(com$site)
# com[1:10, 1:10]
com <- com[,3:26249] ### MUOKATTAVA !!!!


### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Hungary <- oecosimu(com, foo, "quasiswap", nsimul = 99, groups = groups) ### MUOKATTAVA !!!!
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/HUN_sediment_bac_null_models_99perm.txt") ### MUOKATTAVA !!!!
paste("-----")
paste("jaccard, null model")
output_Hungary ### MUOKATTAVA !!!!
sink()
#### Spain ####
com <- bac_sed[DRN == "Spain"] 
com[1:5, 26247:26251]
tot <- rowSums(com[,3:26249]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!site %in% rare_levels)

groups <- as.factor(com$site)
# com[1:10, 1:10]
com <- com[,3:26249] ### MUOKATTAVA !!!!

start_time_total <- Sys.time()
### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Spain <- oecosimu(com, foo, "quasiswap", nsimul = 99, groups = groups) ### MUOKATTAVA !!!!
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/SPA_sediment_bac_null_models_99perm.txt") ### MUOKATTAVA !!!!
paste("-----")
paste("jaccard, null model")
output_Spain ### MUOKATTAVA !!!!
sink()

end_time_total <- Sys.time()
end_time_total - start_time_total
#### next ####
# Next : save all outputs to same dataframe !! 

res1 <- cbind(as.data.frame(output_Croatia$oecosimu$statistic), as.data.frame(output_Croatia$oecosimu$z), as.data.frame(output_Croatia$oecosimu$pval))
colnames(res1) <- c('Statistic','SES','P') 
res1$Country <- "Croatia"

res2 <- cbind(as.data.frame(output_Czech$oecosimu$statistic), as.data.frame(output_Czech$oecosimu$z), as.data.frame(output_Czech$oecosimu$pval))
colnames(res2) <- c('Statistic','SES','P') 
res2$Country <- "Czech Republic"

res3 <- cbind(as.data.frame(output_Finland$oecosimu$statistic), as.data.frame(output_Finland$oecosimu$z), as.data.frame(output_Finland$oecosimu$pval))
colnames(res3) <- c('Statistic','SES','P') 
res3$Country <- "Finland"

res4 <- cbind(as.data.frame(output_France$oecosimu$statistic), as.data.frame(output_France$oecosimu$z), as.data.frame(output_France$oecosimu$pval))
colnames(res4) <- c('Statistic','SES','P') 
res4$Country <- "France"

res5 <- cbind(as.data.frame(output_Hungary$oecosimu$statistic), as.data.frame(output_Hungary$oecosimu$z), as.data.frame(output_Hungary$oecosimu$pval))
colnames(res5) <- c('Statistic','SES','P') 
res5$Country <- "Hungary"

res6 <- cbind(as.data.frame(output_Spain$oecosimu$statistic), as.data.frame(output_Spain$oecosimu$z), as.data.frame(output_Spain$oecosimu$pval))
colnames(res6) <- c('Statistic','SES','P') 
res6$Country <- "Spain"

res <- rbind(res1, res2, res3, res4, res5, res6)
res$Organism <- "Sediment_bacteria" #### MUOKATTAVA !!! 

save(res, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Sediment_bac_jaccard_oecosimu_999p.Rdata") #### MUOKATTAVA !!! 
write.csv(res, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Sediment_bac_jaccard_oecosimu_999p.csv") #### MUOKATTAVA !!! 
# Next : save R workspace
save.image(file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Sediment_bac_jaccard_oecosimu_999p_workspace.Rdata") #### MUOKATTAVA !!! 


