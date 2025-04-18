




SES_TotDur90 <- all[, cor.test(z, TotDur90), by=.(Country, organism)] 
SES_discharge <- all[, cor.test(z, discharge), by=.(Country, organism)] 


<<<<<<< HEAD
# TotLeng90
for (in_organism in unique(in_env_null_models_dt$organism)) {
  print(in_organism)
  ss <- in_env_null_models_dt[organism == in_organism] 
  ss2 <- ss[significance == "sig",]
  p.vals = sapply(unique(ss$Country), function(i) {
    coef(summary(lm(z~TotLeng90, data=ss[Country==i, ])))[2,4] 
  })
  p1 <- ggplot() + 
    geom_point(aes(TotLeng90, z, colour=Country), data = ss, size = 1, alpha=0.2) +
    geom_point(aes(TotLeng90, z, colour=Country), data = ss2, size = 1, alpha=0.5) +
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotLeng90, z, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotLeng90, z, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
    geom_smooth(aes(TotLeng90, z), data=ss, colour="black", method = "lm", linewidth = 1.1, se = F) + # , linetype="dashed"
    theme_classic() +
    ggtitle("Sediment diatoms") + # MUOKKAA
    xlab("Duration of dry periods") +
    ylab("z") +
    scale_color_manual(values = c("Croatia" = "#ef476f",
                                  "Czech Republic" = "#f78c6b", 
                                  "Finland" = "#ffd166", 
                                  "France" = "#06d6a0", 
                                  "Hungary" = "#118ab2",
                                  "Spain" = "#073b4c")) +
    theme(legend.position = "none")
  
  pdf(file.path(resdir, "Null_models", 
                paste0("All_", in_organism, "_SES_vs_TotLeng90_lm2.pdf")), 
      height=3, width=4)
  print(p1)
  dev.off()
}


#Discharge
for (in_organism in unique(in_env_null_models_dt$organism)) {
  print(in_organism)
  ss <- in_env_null_models_dt[organism == in_organism] 
=======
int90 <- int[, list(TotDur90 = mean(TotDur, na.rm=T),
                    TotLeng90 = mean(TotLeng, na.rm=T)), by=.(Site, Country)]
env_mean <- env[, list(discharge = mean(discharge_l_s, na.rm=T),
                       moss = mean(moss_cover, na.rm=T),
                       particle_size = mean(particle_size, na.rm=T)
                       )
                , by=.(Site, Country, stream_type)]

all <- merge(res, int90, by=c("Site", "Country"), all.x=T, sort=F) %>%
  merge(env_mean, by=c("Site", "Country"), all.x=T, sort=F) %>%
  .[Site != 'BUK52', ] %>% #Intermittence indicators are missing here 
  .[, Country := as.factor(Country)]

SES_TotDur90 <- all[, cor.test(z, TotDur90), by=.(Country, organism)] 
SES_discharge <- all[, cor.test(z, discharge), by=.(Country, organism)] 


plt_z_by_stream_type <- function(in_dt) {
  plots <- list()
  for(i in levels(in_dt$Country)){
    d <- subset(in_dt, Country == i)
    plots[[paste0(i)]] <- ggplot(d, aes(x=z, y=Site, color=stream_type)) + 
      scale_colour_manual(values = c("steelblue","orange")) +
      geom_boxplot() + 
      coord_flip() +
      ggtitle(paste(i)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  
  pdf(file.path(resdir, "Null_models", 
                paste0(unique(in_dt$organism), "_SES_vs_streamtype.pdf")), 
      height=10, width=10)
  do.call('grid.arrange', c(plots))
  dev.off()
}

for (in_organism in unique(all$organism)) {
  print(in_organism)
  plt_z_by_stream_type(in_dt = all[organism==in_organism,])
}

#### scatterplot, all countries in the same plot, separately for each organism group ####
all[, significance := fifelse(pval <= 0.05, "sig", "nonsig")]

# TotDur90
for (in_organism in unique(all$organism)) {
  print(in_organism)
  ss <- all[organism == in_organism] 
  ss2 <- ss[significance == "sig",]
  p.vals = sapply(unique(ss$Country), function(i) {
    coef(summary(lm(z~TotDur90, data=ss[Country==i, ])))[2,4] 
  })
  p1 <- ggplot() + 
    geom_point(aes(TotDur90, z, colour=Country), data = ss, size = 1, alpha=0.2) +
    geom_point(aes(TotDur90, z, colour=Country), data = ss2, size = 1, alpha=0.5) +
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, z, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, z, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
    geom_smooth(aes(TotDur90, z), data=ss, colour="black", method = "lm", linewidth = 1.1, se = F) + # , linetype="dashed"
    theme_classic() +
    ggtitle("Sediment diatoms") + # MUOKKAA
    xlab("Number of dry days") +
    ylab("z") +
    scale_color_manual(values = c("Croatia" = "#ef476f",
                                  "Czech Republic" = "#f78c6b", 
                                  "Finland" = "#ffd166", 
                                  "France" = "#06d6a0", 
                                  "Hungary" = "#118ab2",
                                  "Spain" = "#073b4c")) +
    theme(legend.position = "none")

  pdf(file.path(resdir, "Null_models", 
                paste0("All_", in_organism, "_SES_vs_TotDur90_lm2.pdf")), 
      height=3, width=4)
  print(p1)
  dev.off()
}

# TotLeng90
for (in_organism in unique(all$organism)) {
  print(in_organism)
  ss <- all[organism == in_organism] 
  ss2 <- ss[significance == "sig",]
  p.vals = sapply(unique(ss$Country), function(i) {
    coef(summary(lm(z~TotLeng90, data=ss[Country==i, ])))[2,4] 
  })
  p1 <- ggplot() + 
    geom_point(aes(TotLeng90, z, colour=Country), data = ss, size = 1, alpha=0.2) +
    geom_point(aes(TotLeng90, z, colour=Country), data = ss2, size = 1, alpha=0.5) +
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotLeng90, z, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotLeng90, z, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
    geom_smooth(aes(TotLeng90, z), data=ss, colour="black", method = "lm", linewidth = 1.1, se = F) + # , linetype="dashed"
    theme_classic() +
    ggtitle("Sediment diatoms") + # MUOKKAA
    xlab("Duration of dry periods") +
    ylab("z") +
    scale_color_manual(values = c("Croatia" = "#ef476f",
                                  "Czech Republic" = "#f78c6b", 
                                  "Finland" = "#ffd166", 
                                  "France" = "#06d6a0", 
                                  "Hungary" = "#118ab2",
                                  "Spain" = "#073b4c")) +
    theme(legend.position = "none")
  
  pdf(file.path(resdir, "Null_models", 
                paste0("All_", in_organism, "_SES_vs_TotLeng90_lm2.pdf")), 
      height=3, width=4)
  print(p1)
  dev.off()
}


#Discharge
for (in_organism in unique(all$organism)) {
  print(in_organism)
  ss <- all[organism == in_organism] 
>>>>>>> ccdc49c786ea61bd62b5f2be33108a7f6e707e04
  ss2 <- ss[significance == "sig",]
  p.vals = sapply(unique(ss$Country), function(i) {
    coef(summary(lm(z~discharge, data=ss[Country==i, ])))[2,4] 
  })
  p1 <- ggplot() + 
    geom_point(aes(discharge, z, colour=Country), data = ss, size = 1, alpha=0.2) +
    geom_point(aes(discharge, z, colour=Country), data = ss2, size = 1, alpha=0.5) +
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(discharge, z, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
    geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(discharge, z, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
    geom_smooth(aes(discharge, z), data=ss, colour="black", method = "lm", linewidth = 1.1, se = F) + # , linetype="dashed"
    theme_classic() +
    ggtitle("Sediment diatoms") + # MUOKKAA
    xlab("Discharge") +
    ylab("z") +
    scale_color_manual(values = c("Croatia" = "#ef476f",
                                  "Czech Republic" = "#f78c6b", 
                                  "Finland" = "#ffd166", 
                                  "France" = "#06d6a0", 
                                  "Hungary" = "#118ab2",
                                  "Spain" = "#073b4c")) +
    theme(legend.position = "none")
  
  pdf(file.path(resdir, "Null_models", 
                paste0("All_", in_organism, "_SES_vs_discharge_lm2.pdf")), 
      height=3, width=4)
  print(p1)
  dev.off()
}
<<<<<<< HEAD
=======

>>>>>>> ccdc49c786ea61bd62b5f2be33108a7f6e707e04

#### mixed models ####
lmer_int <- all[,
  list(
    TotDur90_full = list(lmer(z ~ TotDur90 + (1|Country), data=.SD)),
    TotDur90_null = list(lmer(z ~ (1|Country), data=.SD)),
    TotDur90_ML = list(anova(lmer(z ~ TotDur90 + (1|Country), data=.SD),
                        lmer(z ~ (1|Country), data=.SD))),
    TotLeng90_full = list(lmer(z ~ TotLeng90 + (1|Country), data=.SD)),
    TotLeng90_null = list(lmer(z ~ (1|Country), data=.SD)),
    TotLeng90_ML = list(anova(lmer(z ~ TotLeng90 + (1|Country), data=.SD),
                         lmer(z ~ (1|Country), data=.SD)))
  )
  , by=organism] 


#### jitter plots to show deviation of SES from zero, all countries in the same plot, separately for each organism group ####
z_jitter_by_organism <- function(in_dt) {
  jitter_p <- ggplot(in_dt) + 
    scale_y_continuous() +
    geom_jitter(aes(Country, z, colour=significance), 
                shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
    theme_classic() +
    theme(axis.title.x = element_blank()) +
    labs(y = "z") +
    ggtitle(unique(in_dt$organism)) +
    scale_color_manual(values = c("sig" = "mediumblue",
                                  "nonsig" = "lightslateblue")) +
    labs(colour = "Departure from zero") +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1),
          legend.position = 'non')
  return(jitter_p)
}

plots_jitter <- list()
for (in_organism in unique(all$organism)) {
  print(in_organism)
  plots_jitter[[in_organism]] <- z_jitter_by_organism(in_dt = all[organism==in_organism,])
}

pdf(file.path(resdir, "Null_models", "Jitterplots_SES_significance.pdf"),
    height=10, width=6)
do.call("grid.arrange", c(plots_jitter, ncol=2))
dev.off()



###############################################################################################
# rm(list = ls())
# library(data.table)
# library(vegan)
# 
# # closeAllConnections()
# res <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/null_model_results_reduced.csv")
# env <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/ENV_all_NAs_as_blank.csv")
# int <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/Intermittence_all.csv") 
# 
# # TotDur= Number of dry days (how many zeros there are)
# # TotNum= Number of DRY PERIODS (how many times there is a change between 1 to 0
# # TotLeng= Average length of all DRY PERIODS
# int1 <- fread(file.path(datdir, 'Datasets',  "Intermittence_Data/Croatia_Local_Interm_90_d.csv"))
# int2 <- fread(file.path(datdir, 'Datasets',"Intermittence_Data/Czech Republic_Local_Interm_90_d.csv"))
# int3 <- fread(file.path(datdir, 'Datasets',"Intermittence_Data/Finland_Local_Interm_90_d.csv"))
# int4 <- fread(file.path(datdir, 'Datasets',"Intermittence_Data/France_Local_Interm_90_d_corrected.csv"))
# int5 <- fread(file.path(datdir, 'Datasets',"Intermittence_Data/Hungary_Local_Interm_90_d.csv"))
# int6 <- fread(file.path(datdir, 'Datasets',"Intermittence_Data/Spain_Local_Interm_90_d.csv"))
# 
# int2 <- rbind(int1, int2, int3, int4, int5, int6)
# # Lasketaan paikoille keskiarvo kuivien päivien lukumäärästä 90 vrk ajalla ennen sämpläystä (koska 90 vrk näyttää olevan max arvo useimmissa paikoissa)
# names(int)
# library(dplyr)
# TotDur90<- int%>%group_by(Sites)%>%summarise(TotDur90=mean(TotDur, na.rm = T))
# TotDur90$site <- TotDur90$Sites
# # Lasketaan paikoille keskiarvo kuivien kausien kestosta 90 vrk ajalla ennen sämpläystä
# TotLeng90<- int%>%group_by(Sites)%>%summarise(TotLeng90=mean(TotLeng, na.rm = T))
# TotLeng90$site <- TotLeng90$Sites
#
# Calculate the average of the locations of the environmental variables, 1-6 samples per location
# names(env)
# library(dplyr)
# discharge<- env%>%group_by(site)%>%summarise(discharge=mean(discharge_l_s, na.rm = T))
# moss<- env%>%group_by(site)%>%summarise(moss=mean(moss_cover, na.rm = T))
# particlesize <- env%>%group_by(site)%>%summarise(particlesize=mean(particle_size, na.rm = T))
#
# names(int)
# ConD<- int%>%group_by(Sites)%>%summarise(ConD=mean(ConD, na.rm = T))
# ConD$site <- ConD$Sites
# 
# ConDmax<- int%>%group_by(Sites)%>%summarise(ConDmax=max(ConD, na.rm = T))
# ConDmax$site <- ConD$Sites
#
#
# type<- env[,c(3,6)]
# type <- unique(type)
# res$site <- res$Site
# 
# all <- merge(res, TotDur90, by="site", all.x=T, sort=F)
# all <- merge(all, TotLeng90, by="site", all.x=T, sort=F) 
# all <- merge(all, moss, by="site", all.x=T, sort=F) 
# all <- merge(all, particlesize, by="site", all.x=T, sort=F)
# # all <- merge(all, ConD, by="site", all.x=T, sort=F)
# all <- merge(all, type, by="site", all.x=T, sort=F)
# # all <- merge(all, ConDmax, by="site", all.x=T, sort=F)
# all <- merge(all, discharge, by="site", all.x=T, sort=F)
# 
# 
# all <- droplevels(all[!all$Site == 'BUK52',]) # täältä puuttuu ainakin intermittence indicators
#
#
#
#### correlation test ####
# correlation test, countries
# ss <- all[Country == "Finland"] 
# cor.test(all[Country == "Finland"]$SES, all[Country == "Finland"]$ConDmax)
#
# levels(factor(all$Country))
# levels(factor(all$Organism))
# 
# ss <- all[Organism == "Sediment_fungi"] 
# 
# 
# # TotDur90
# sink("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_fungi_SES_TotDur90_correlations.txt")
# paste("-----")
# paste("Finland, correlation")
# cor.test(ss[Country == "Finland"]$SES, ss[Country == "Finland"]$TotDur90)
# paste("-----")
# paste("Croatia, correlation")
# cor.test(ss[Country == "Croatia"]$SES, ss[Country == "Croatia"]$TotDur90)
# paste("-----")
# paste("France, correlation")
# cor.test(ss[Country == "France"]$SES, ss[Country == "France"]$TotDur90)
# paste("-----")
# paste("Hungary, correlation")
# cor.test(ss[Country == "Hungary"]$SES, ss[Country == "Hungary"]$TotDur90)
# paste("-----")
# paste("Spain, correlation")
# cor.test(ss[Country == "Spain"]$SES, ss[Country == "Spain"]$TotDur90)
# paste("-----")
# paste("Czech, correlation")
# cor.test(ss[Country == "Czech Republic"]$SES, ss[Country == "Czech Republic"]$TotDur90)
# sink()
# 
# # discharge
# sink("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_fungi_SES_Discharge_correlations.txt")
# paste("-----")
# paste("Finland")
# cor.test(ss[Country == "Finland"]$SES, ss[Country == "Finland"]$discharge)
# paste("-----")
# paste("Croatia")
# cor.test(ss[Country == "Croatia"]$SES, ss[Country == "Croatia"]$discharge)
# paste("-----")
# paste("Czech")
# cor.test(ss[Country == "Czech Republic"]$SES, ss[Country == "Czech Republic"]$discharge)
# paste("-----")
# paste("France")
# cor.test(ss[Country == "France"]$SES, ss[Country == "France"]$discharge)
# paste("-----")
# paste("Hungary")
# cor.test(ss[Country == "Hungary"]$SES, ss[Country == "Hungary"]$discharge)
# paste("-----")
# paste("Spain")
# cor.test(ss[Country == "Spain"]$SES, ss[Country == "Spain"]$discharge)
# sink()
#
# # boxplot for a country
# ss <- all[Country == "Finland"]
# p <- ggplot(ss, aes(x=SES, y=site, color=stream_type)) +
#   geom_boxplot() + coord_flip()
# p
# 
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/FIN_biof_fungi_SES_vs_streamtype.pdf"), height=6.4, width=10) # max height=9.65, width=6.4
# p
# dev.off()
#
# # all in figure
# library(ggplot2)
# 
# levels(factor(all$Organism))
# 
# # tehty 18032024
# 
# #### boxplot: SES vs. stream type ####
# plots <- list()
# for(i in unique(all$Country)){
#   d <- subset(all, Country == i)
#   plots[[paste0(i)]] <- ggplot(d[Organism=="Macroinvertebrates"], aes(x=SES, y=site, color=stream_type)) + 
#     scale_colour_manual(values = c("steelblue","orange")) +
#     geom_boxplot() + coord_flip() +
#     ggtitle(paste(i)) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# }
# 
# library(gridExtra)
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Macroinvertebrates_SES_vs_streamtype.pdf"), height=10, width=10)
# do.call('grid.arrange', c(plots))
# dev.off()
# 
# 
# #### scatterplot, all countries in the same plot, separately for each organism group ####
# 
# # make new column based on another column's values
# library(dplyr)
# levels <- c(-Inf, 0.05,  Inf)
# labels <- c("sig", "nonsig")
# all2 <- all %>% mutate(Significance2 = cut(Significance, levels, labels = labels))
# 
# library(ggplot2)
# 
# # levels(factor(all$Country))
# levels(factor(all$Organism))
# # min(all$mean_S)
# # max(all$mean_S)
# 
# ss <- all2[Organism == "Sediment_diatoms"] 
# ss2 <- ss[Significance2 == "sig"]
# # TotDur90
# # summary(lm(ss$SES~ ss$TotDur90))
# 
# p.vals = sapply(unique(ss$Country), function(i) {
#   coef(summary(lm(SES~TotDur90, data=ss[Country==i, ])))[2,4] 
# })
# p.vals
# p1 <- ggplot() + 
#   geom_point(aes(TotDur90, SES, colour=Country), data = ss, size = 1, alpha=0.2) +
#   geom_point(aes(TotDur90, SES, colour=Country), data = ss2, size = 1, alpha=0.5) +
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, SES, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, SES, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
#   geom_smooth(aes(TotDur90, SES), data=ss, colour="black", method = "lm", linewidth = 1.1, se = F) + # , linetype="dashed"
#   theme_classic() +
#   ggtitle("Sediment diatoms") + # MUOKKAA
#   xlab("Number of dry days") +
#   ylab("SES") +
#   scale_color_manual(values = c("Croatia" = "#ef476f",
#                                 "Czech Republic" = "#f78c6b", 
#                                 "Finland" = "#ffd166", 
#                                 "France" = "#06d6a0", 
#                                 "Hungary" = "#118ab2",
#                                 "Spain" = "#073b4c")) +
#   theme(legend.position = "none")
# p1
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/All_Sediment_diatoms_SES_vs_TotDur90_lm2.pdf"), height=3, width=4) # MUOKKAA
# p1
# dev.off()
# 
# # TotLeng90
# summary(lm(ss$SES~ ss$TotLeng90))
# 
# p.vals = sapply(unique(ss$Country), function(i) {
#   coef(summary(lm(SES~TotLeng90, data=ss[Country==i, ])))[2,4] 
# })
# p.vals
# p1 <- ggplot() +
#   geom_point(aes(TotLeng90, SES, colour=Country), data = ss, size = 1, alpha=0.2) +
#   geom_point(aes(TotLeng90, SES, colour=Country), data = ss2, size = 1, alpha=0.5) +
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotLeng90, SES, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotLeng90, SES, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
#   geom_smooth(aes(TotLeng90, SES), data=ss, colour="black", method = "lm", linewidth = 1.1, se = F) + # , linetype="dashed"
#   theme_classic() +
#   ggtitle("Sediment diatoms") + # MUOKKAA
#   xlab("Duration of dry periods") +
#   ylab("SES") +
#   scale_color_manual(values = c("Croatia" = "#ef476f",
#                                 "Czech Republic" = "#f78c6b", 
#                                 "Finland" = "#ffd166", 
#                                 "France" = "#06d6a0", 
#                                 "Hungary" = "#118ab2",
#                                 "Spain" = "#073b4c")) +
#   theme(legend.position = "none")
# p1
# 
# 
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/All_Sediment_diatoms_SES_vs_TotLeng90_lm2.pdf"), height=3, width=4) # MUOKKAA
# p1
# dev.off()
#
# # testaus: polynomial
# 
# 
# ss <- all2[Organism == "Sediment_diatoms"] 
# ss2 <- ss[Significance2 == "sig"]
# 
# p.vals = sapply(unique(ss$Country), function(i) {
#   coef(summary(lm(SES~TotDur90 + I(TotDur90^2), data=ss[Country==i, ])))[2,4] 
# })
# p.vals
# p1 <- ggplot() + 
#   geom_point(aes(TotDur90, SES, colour=Country), data = ss, size = 1, alpha=0.2) +
#   geom_point(aes(TotDur90, SES, colour=Country), data = ss2, size = 1, alpha=0.5) +
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, SES, colour=Country), method = "loess", linewidth = 0.75, se = F) + 
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, SES, colour=Country), method = "loess", linewidth = 0.75, se = F, linetype="dashed") +
#   geom_smooth(aes(TotDur90, SES), data=ss, colour="black", method = "loess", linewidth = 1.1, se = F) + # , linetype="dashed"
#   theme_classic() +
#   ggtitle("Sediment diatoms") + # MUOKKAA
#   xlab("Number of dry days") +
#   ylab("SES") +
#   scale_color_manual(values = c("Croatia" = "#ef476f",
#                                 "Czech Republic" = "#f78c6b", 
#                                 "Finland" = "#ffd166", 
#                                 "France" = "#06d6a0", 
#                                 "Hungary" = "#118ab2",
#                                 "Spain" = "#073b4c")) +
#   theme(legend.position = "none")
# p1
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/All_Sediment_diatoms_SES_vs_TotDur90_lm2.pdf"), height=3, width=4) # MUOKKAA
# p1
# dev.off()
# 
# 
# 




#### ylimääräisiä plotteja ####

# # scatterplot, discharge -- linear models with significance as colour
# levels(factor(all$Organism))
# ss <- all2[Organism == "Sediment_fungi"] 
# plots3 <- list()
# for(i in unique(ss$Country)){
#   d <- subset(ss, Country == i)
#   p.vals = sapply(unique(d$Country), function(i) {
#     coef(summary(lm(SES~discharge, data=d[Country==i, ])))[2,4] 
#   })
#   plots3[[paste0(i)]] <- ggplot(d[Organism=="Sediment_fungi"], aes(x=discharge, y=SES)) + 
#     geom_point(aes(x=discharge, y=SES, colour = Significance2), size = 2) + 
#     geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals < 0.05],], aes(discharge, SES), method = "lm", linewidth = 0.5, se = F) + # vaihda kem
#     geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals > 0.05],], aes(discharge, SES), method = "lm", linewidth = 0.5, se = F, linetype="dashed") + # vaihda kem
#     scale_colour_manual(values = c("steelblue","orange")) +
#     ggtitle(paste(i)) +
#     theme_classic()
# }
# library(gridExtra)
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_fungi_SES_vs_Discharge_lm_sig.pdf"), height=10, width=10)
# do.call('grid.arrange', c(plots3))
# dev.off()


# # TÄSTÄ LÄHTEE
# ss <- all2[Organism == "Flying_macroinvertebrates"] 
# # scatterplot, TotDur90 -- linear models with significance as colour
# plots2 <- list()
# for(i in unique(ss$Country)){
#   d <- subset(ss, Country == i)
#   p.vals = sapply(unique(d$Country), function(i) {
#     coef(summary(lm(SES~TotDur90, data=d[Country==i, ])))[2,4] 
#   })
#   plots2[[paste0(i)]] <- ggplot(d[Organism=="Biofilm_diatoms"], aes(x=TotDur90, y=SES)) + 
#     geom_point(aes(x=TotDur90, y=SES, colour = Significance2), size = 2) + 
#     geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, SES), method = "lm", linewidth = 0.5, se = F) + # vaihda kem
#     geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, SES), method = "lm", linewidth = 0.5, se = F, linetype="dashed") + # vaihda kem
#     scale_colour_manual(values = c("steelblue","orange")) +
#     ggtitle(paste(i)) +
#     theme_classic()
# }
# library(gridExtra)
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Biofilm_diatoms_SES_vs_TotDur90_lm_sig.pdf"), height=10, width=10)
# do.call('grid.arrange', c(plots2))
# dev.off()

# # scatterplot discharge quadratic
# plots3 <- list()
# for(i in unique(all$Country)){
#   d <- subset(all, Country == i)
#   p.vals = sapply(unique(d$Country), function(i) {
#     coef(summary(lm(SES~discharge+I(discharge^2), data=d[Country==i, ])))[2,4]
#   })
#   plots3[[paste0(i)]] <- ggplot(d[Organism=="Macroinvertebrates"], aes(x=discharge, y=SES)) +
#     geom_point(size = 2) +
#     geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals < 0.05],], aes(discharge, SES), method = "lm",formula = y~x+I(x^2), linewidth = 0.5, se = F) + # vaihda kem
#     geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals > 0.05],], aes(discharge, SES), method = "lm",formula = y~x+I(x^2), linewidth = 0.5, se = F, linetype="dashed") + # vaihda kem
#     # scale_colour_manual(values = c("steelblue","orange"))  +
#     ggtitle(paste(i)) +
#     theme_classic()
# }

# # scatterplot, country-specific, all organisms in the same plot
# levels(factor(all$Country))
# levels(factor(all$Organism))
# min(all$SES)
# max(all$SES)
# ss <- all[Country=="Finland"] # MUOKKAA
# 
# p.vals = sapply(unique(ss$Organism), function(i) {
#     coef(summary(lm(SES~TotDur90, data=ss[Organism==i, ])))[2,4] 
#   })
# p.vals
# p1 <- ggplot(ss) + 
#   # geom_point(size = 2, alpha=0.5) + 
#   geom_smooth(data=ss[ss$Organism %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, SES, colour=Organism), method = "lm", linewidth = 0.5, se = F) + 
#   geom_smooth(data=ss[ss$Organism %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, SES, colour=Organism), method = "lm", linewidth = 0.5, se = F, linetype="dashed") +
#   theme_classic() +
#   ggtitle("Finland") + # MUOKKAA
#   scale_color_manual(values = c("Biofilm_fungi" = "cornflowerblue",
#                                 "Macroinvertebrates" = "mediumseagreen", 
#                                 "Sediment_fungi" = "darkblue", 
#                                 "Biofilm_diatoms" = "bisque3", 
#                                 "Sediment_diatoms" = "rosybrown4"))
# p1
# 
# 
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Finland_SES_vs_TotDur90_lm.pdf"), height=3, width=6) # MUOKKAA
# p1
# dev.off()

#### mixed models ####
# 
# levels(factor(all$Organism))
# library(lme4)
# # install.packages("lmerTest")
# library(MuMIn)
# library(lmerTest)
# 
# org <- "Sediment_bacteria" ### MUOKKAA
# sink(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/New_mixed_models_SES_", org, ".txt", sep="")) 
# paste(org)
# paste("---------------------------------------------------------")
# paste("TotDur90 full model")
# print(summary(lmer(SES ~ TotDur90 + (1|Country), data=all[Organism==org]))) 
# paste("---------------------------------------------------------")
# paste("TotDur90 null model")
# print(summary(lmer(SES ~ (1|Country), data=all[Organism==org]))) 
# paste("---------------------------------------------------------")
# paste("Anova, refitting model(s) with ML (instead of REML)")
# anova(lmer(SES ~ TotDur90 + (1|Country), data=all[Organism==org]), lmer(SES ~ (1|Country), data=all[Organism==org])) 
# paste("---------------------------------------------------------")
# paste("TotLeng90 full model")
# print(summary(lmer(SES ~ TotLeng90 + (1|Country), data=all[Organism==org]))) 
# paste("---------------------------------------------------------")
# paste("TotLeng90 null model")
# print(summary(lmer(SES ~ (1|Country), data=all[Organism==org]))) 
# paste("---------------------------------------------------------")
# paste("Anova, refitting model(s) with ML (instead of REML)")
# anova(lmer(SES ~ TotLeng90 + (1|Country), data=all[Organism==org]), lmer(SES ~ (1|Country), data=all[Organism==org])) 
# sink()
# 
# #### jitter plots to show deviation of SES from zero, all countries in the same plot, separately for each organism group ####
# 
# 
# # make new column based on another column's values
# library(dplyr)
# levels <- c(-Inf, 0.05,  Inf)
# labels <- c("Significant", "Non-significant")
# all2 <- all %>% mutate(Significance2 = cut(Significance, levels, labels = labels))
# 
# library(ggplot2)
# 
# names(all2)
# levels(factor(all2$Organism))
# 
# ss <- all2[Organism == "Biofilm_bacteria"] 
# jitter1 <- ggplot(ss) + 
#   scale_y_continuous() +
#   geom_jitter(aes(Country, SES, colour=Significance2), shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
#   theme_classic() +
#   theme(axis.title.x = element_blank()) +
#   labs(y = "SES") +
#   ggtitle("Biofilm bacteria") +
#   scale_color_manual(values = c("Significant" = "mediumblue",
#                                 "Non-significant" = "lightslateblue")) +
#   labs(colour = "Departure from zero") +
#   theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
#   # theme(legend.position = "bottom")
#   # theme(legend.position = "none")
# plot(jitter1)
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Jitterplot_with_legend.pdf"), height=3, width=4) # MUOKKAA
# jitter1
# dev.off()
# 
# ss <- all2[Organism == "Sediment_bacteria"] 
# jitter2 <- ggplot(ss) + 
#   scale_y_continuous() +
#   geom_jitter(aes(Country, SES, colour=Significance2), shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
#   theme_classic() +
#   theme(axis.title.x = element_blank()) +
#   labs(y = "SES") +
#   ggtitle("Sediment bacteria") +
#   scale_color_manual(values = c("Significant" = "mediumblue",
#                                 "Non-significant" = "lightslateblue")) +
#   labs(colour = "Departure from zero") +
#   theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
#   # theme(legend.position = "bottom")
#   theme(legend.position = "none")
# plot(jitter2)
# 
# ss <- all2[Organism == "Biofilm_diatoms"] 
# jitter3 <- ggplot(ss) + 
#   scale_y_continuous() +
#   geom_jitter(aes(Country, SES, colour=Significance2), shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
#   theme_classic() +
#   theme(axis.title.x = element_blank()) +
#   labs(y = "SES") +
#   ggtitle("Biofilm diatoms") +
#   scale_color_manual(values = c("Significant" = "mediumblue",
#                                 "Non-significant" = "lightslateblue")) +
#   labs(colour = "Departure from zero") +
#   theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
#   # theme(legend.position = "bottom")
#   theme(legend.position = "none")
# plot(jitter3)
# 
# levels(factor(all2$Organism))
# 
# ss <- all2[Organism == "Sediment_diatoms"] 
# jitter4 <- ggplot(ss) + 
#   scale_y_continuous() +
#   geom_jitter(aes(Country, SES, colour=Significance2), shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
#   theme_classic() +
#   theme(axis.title.x = element_blank()) +
#   labs(y = "SES") +
#   ggtitle("Sediment diatoms") +
#   scale_color_manual(values = c("Significant" = "mediumblue",
#                                 "Non-significant" = "lightslateblue")) +
#   labs(colour = "Departure from zero") +
#   theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
#   # theme(legend.position = "bottom")
#   theme(legend.position = "none")
# plot(jitter4)
# 
# levels(factor(all2$Organism))
# 
# ss <- all2[Organism == "Biofilm_fungi"] 
# jitter5 <- ggplot(ss) + 
#   scale_y_continuous() +
#   geom_jitter(aes(Country, SES, colour=Significance2), shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
#   theme_classic() +
#   theme(axis.title.x = element_blank()) +
#   labs(y = "SES") +
#   ggtitle("Biofilm fungi") +
#   scale_color_manual(values = c("Significant" = "mediumblue",
#                                 "Non-significant" = "lightslateblue")) +
#   labs(colour = "Departure from zero") +
#   theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
#   # theme(legend.position = "bottom")
#   theme(legend.position = "none")
# plot(jitter5)
# 
# ss <- all2[Organism == "Sediment_fungi"] 
# jitter6 <- ggplot(ss) + 
#   scale_y_continuous() +
#   geom_jitter(aes(Country, SES, colour=Significance2), shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
#   theme_classic() +
#   theme(axis.title.x = element_blank()) +
#   labs(y = "SES") +
#   ggtitle("Sediment fungi") +
#   scale_color_manual(values = c("Significant" = "mediumblue",
#                                 "Non-significant" = "lightslateblue")) +
#   labs(colour = "Departure from zero") +
#   theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
#   # theme(legend.position = "bottom")
#   theme(legend.position = "none")
# plot(jitter6)
# 
# levels(factor(all2$Organism))
# 
# ss <- all2[Organism == "Macroinvertebrates"] 
# jitter7 <- ggplot(ss) + 
#   scale_y_continuous() +
#   geom_jitter(aes(Country, SES, colour=Significance2), shape=16, size = 2, alpha=0.7, position=position_jitter(0.2)) +
#   theme_classic() +
#   theme(axis.title.x = element_blank()) +
#   labs(y = "SES") +
#   ggtitle("Macroinvertebrates") +
#   scale_color_manual(values = c("Significant" = "mediumblue",
#                                 "Non-significant" = "lightslateblue")) +
#   labs(colour = "Departure from zero") +
#   theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) +
#   # theme(legend.position = "bottom")
#   theme(legend.position = "none")
# plot(jitter7)
# 
# library(gridExtra)
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Jitterplots_SES_significance.pdf"), height=10, width=6)
# grid.arrange(jitter1, 
#              jitter2, 
#              jitter3,
#              jitter4,
#              jitter5,
#              jitter6, 
#              jitter7, 
#              ncol=2)
# dev.off()
