# DRYvER 

# Site-specific  assembly processes
# - How does site-specific dying affect local assembly processes? (cf. Huttunen et al. 2017)


rm(list = ls())
library(data.table)
library(vegan)

# ENV data

# env <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/ENV_all_NAs_as_blank.csv") # running_id used for merging

# FUNGI - biofilm
# ------------

fungi_biof <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/fungi_dna_Biof_removed_zero_rows_and_columns.csv")
fungi_biof[1:10, 1:10]
# make same running_id
# fungi_biof$running_id <- paste(fungi_biof$Site, fungi_biof$Campaign, sep = "_") 
# fungi_biof[1:10, 3870:3873]
# # combine state_of_flow
# fungi_biof <- merge(fungi_biof, env[,.(running_id, state_of_flow)], by="running_id", all.x=T, sort=F) 
# fungi_biof[1:10, 3870:3874]

# erikseen DRN:t
levels(factor(fungi_biof$Country))
com <- fungi_biof[Country == "Croatia"] 
# remove zero/empty rows
tot <- rowSums(com[,6:1982])
com <- com[tot > 0, ]

# remove sites with only one sample # Croatia and Spain produced an error -- there is only one sample per site in some occasions

tt <- table(com$Site)
rare_levels <- names(tt)[tt<2]
com <- subset(com,!Site %in% rare_levels)

# make group and response tables
# groups<- as.factor(com$Campaign)

groups <- as.factor(com$Site)
com[1:10, 1:10]
com <- com[,6:1982]

### The models start here
start_time <- Sys.time()

# null model, jaccard
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output <- oecosimu(com, foo, "quasiswap", nsimul = 99, groups = groups)

# # havaitulla jaccardilla, ilman nollamallia
# jac.obs <-vegdist(com, "jaccard", binary =TRUE)
# beta.jac.obs <-betadisper(jac.obs, groups)
# beta.jac.obs
# 
# # null model, bray
# foo <-function(x, groups, ...) diag(meandist(vegdist(x, "bray", binary =TRUE), grouping = groups))
# output_bray <- oecosimu(com, foo, "quasiswap", nsimul = 99, groups = groups)

end_time <- Sys.time()
end_time - start_time

# Saving
sink("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/CRO_biofilm_fungi_null_models.txt")
paste("-----")
paste("jaccard, null model")
output
# paste("-----")
# paste("jaccard, observed")
# beta.jac.obs
# paste("-----")
# paste("bray, null model")
# output_bray
sink()





# MI - Finland
# ------------

mi_fi <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/MIV_Fin_zero_columns_removed.csv")
mi_fi[1:10, 1:10]

# # remove zero/empty rows
# tot2 <- rowSums(mi_fi[,7:158])
# mi_fi <- mi_fi[tot2 > 0, ]

# make group and response tables
# groups<- as.factor(mi_fi$Campaign)

groups<- as.factor(mi_fi$Site)
mi_fi[1:10, 1:10]
mi_fi <- mi_fi[,7:158]

# null model
foo <-function(x, groups, ...) diag(meandist(vegdist(x, "bray", binary =TRUE), grouping = groups))
output_bray_mi_fi <- oecosimu(mi_fi, foo, "quasiswap", nsimul = 99, groups = groups)

foo <-function(x, groups, ...) diag(meandist(vegdist(x, "jac", binary =TRUE), grouping = groups))
output_mi_fi <- oecosimu(mi_fi, foo, "quasiswap", nsimul = 99, groups = groups)

# havaitulla jaccardilla, ilman nollamallia
jac.obs_mi_fi <-vegdist(mi_fi, "jaccard", binary =TRUE)
beta.jac.obs_mi_fi <-betadisper(jac.obs_mi_fi, groups)
beta.jac.obs_mi_fi

sink("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/macinv_fin_null_models.txt")
paste("jaccard, nollamalli")
output_mi_fi
paste("bray, nollamalli")
output_bray_mi_fi
paste("jaccard")
beta.jac.obs_mi_fi
sink()

# tsekkaa nmds-plotti

# hellinger transformation and euclidean distance in NMDS
bio <- decostand(mi_fi, "hellinger")
# (log1p(bio), method = 'hellinger'))
nmds = metaMDS(bio, distance = "euclidean", k=2)
nmds

plot(nmds, "sites")



env <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/ENV_all_NAs_as_blank.csv") # running_id used for merging
env$Running_ID <- env$running_id
env_2 <- merge(mi_fi2, env[,.(Running_ID, stream_type)], by="Running_ID", all.x=T, sort=F) 

# extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds)$sites)

mi_fi2[1:10, 1:10]

# data.scores$tyyppi <- exp$PaikanGeoTyyppi
data.scores$type <- env_2$stream_type
levels(factor(data.scores$type))

library(ggplot2)
# plotti ilman envfittiä
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, aes(colour = type), size = 3, alpha = 0.5) +
  scale_colour_manual(values = c("steelblue","orange")) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"),
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "type")
gg

pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/Macinv_Fin_nmds_plot_with_hellinger_transformation_euclidean_distances.pdf"), height=6.4, width=8) # max height=9.65, width=6.4
gg
dev.off()