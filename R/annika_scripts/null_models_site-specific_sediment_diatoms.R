# DRYvER 

# Site-specific  assembly processes
# - How does site-specific dying affect local assembly processes? (cf. Huttunen et al. 2017)


rm(list = ls())
library(data.table)
library(vegan)

# DIATOMS - sediment
# ------------

dia_sed <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/diatoms_dna_sediment_removed_pools_and_zero_rows_and_columns.csv")
dia_sed[1:10, 1:10]

levels(factor(dia_sed$Country))
dia_sed$running_id <- NULL
dia_sed$summa_sample <- NULL
dia_sed[1:10, 1:10]
#### Croatia ####
com <- dia_sed[Country == "Croatia"] 

tot <- rowSums(com[,6:1299]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$Site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!Site %in% rare_levels)

groups <- as.factor(com$Site)
# com[1:10, 1:10]
com <- com[,6:1299] ### MUOKATTAVA !!!!

start_time_total <- Sys.time()
### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Croatia <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups)
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/CRO_sediment_dia_null_models_999perm.txt") ### MUOKATTAVA !!!!
paste("-----")
paste("jaccard, null model")
output_Croatia
sink()

#### Czech ####
com <- dia_sed[Country == "Czech"] ### MUOKATTAVA !!!

tot <- rowSums(com[,6:1299]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$Site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!Site %in% rare_levels)

groups <- as.factor(com$Site)
# com[1:10, 1:10]
com <- com[,6:1299] ### MUOKATTAVA !!!

### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Czech <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups)
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/CZE_sediment_dia_null_models_999perm.txt") ### MUOKATTAVA !!!
paste("-----")
paste("jaccard, null model")
output_Czech
sink()

#### Finland ####
com <- dia_sed[Country == "Finland"]  ### MUOKATTAVA !!!

tot <- rowSums(com[,6:1299]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$Site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!Site %in% rare_levels)

groups <- as.factor(com$Site)
# com[1:10, 1:10]
com <- com[,6:1299] ### MUOKATTAVA !!!

### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Finland <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups)
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/FIN_sediment_dia_null_models_999perm.txt") ### MUOKATTAVA !!!
paste("-----")
paste("jaccard, null model")
output_Finland
sink()

#### France ####
com <- dia_sed[Country == "France"]  ### MUOKATTAVA !!!

tot <- rowSums(com[,6:1299]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$Site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!Site %in% rare_levels)

groups <- as.factor(com$Site)
# com[1:10, 1:10]
com <- com[,6:1299] ### MUOKATTAVA !!!

### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_France <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups)
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/FRA_sediment_dia_null_models_999perm.txt") ### MUOKATTAVA !!!
paste("-----")
paste("jaccard, null model")
output_France
sink()

#### Hungary ####
com <- dia_sed[Country == "Hungary"]  ### MUOKATTAVA !!!

tot <- rowSums(com[,6:1299]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$Site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!Site %in% rare_levels)

groups <- as.factor(com$Site)
# com[1:10, 1:10]
com <- com[,6:1299] ### MUOKATTAVA !!!

### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Hungary <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups)
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/HUN_sediment_dia_null_models_999perm.txt") ### MUOKATTAVA !!!
paste("-----")
paste("jaccard, null model")
output_Hungary
sink()

#### Spain ####
com <- dia_sed[Country == "Spain"] ### MUOKATTAVA !!!

tot <- rowSums(com[,6:1299]) ### MUOKATTAVA !!!!
com <- com[tot > 0, ]

tt <- table(com$Site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!Site %in% rare_levels)

groups <- as.factor(com$Site)
com[1:10, 1:10]
com <- com[,6:1299] ### MUOKATTAVA !!!

### The models start here
start_time <- Sys.time()
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_Spain <- oecosimu(com, foo, "quasiswap", nsimul = 999, groups = groups)
end_time <- Sys.time()
end_time - start_time
# Saving
sink("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/SPA_sediment_dia_null_models_999perm.txt") ### MUOKATTAVA !!!
paste("-----")
paste("jaccard, null model")
output_Spain
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
res$Organism <- "Sediment_diatoms" #### MUOKATTAVA !!! 

save(res, file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Sediment_dia_jaccard_oecosimu_999p.Rdata") #### MUOKATTAVA !!! 
# Next : save R workspace
save.image(file = "C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Sediment_dia_jaccard_oecosimu_999p_workspace.Rdata") #### MUOKATTAVA !!! 


