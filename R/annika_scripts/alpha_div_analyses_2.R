## Diversity analyses ##
## - same as for SES values
rm(list = ls())
library(data.table)
library(vegan)
# closeAllConnections()

div <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/Mean_diversity_updated2.csv") # includes bacteria, non-fly and fly macroinvertebrates
env <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/ENV_all_NAs_as_blank.csv")

int1 <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/Intermittence_Data/Croatia_Local_Interm_90_d.csv")
int2 <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/Intermittence_Data/Czech Republic_Local_Interm_90_d.csv")
int3 <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/Intermittence_Data/Finland_Local_Interm_90_d.csv")
int4 <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/Intermittence_Data/France_Local_Interm_90_d_corrected.csv")
int5 <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/Intermittence_Data/Hungary_Local_Interm_90_d.csv")
int6 <- fread("M:/gDRYvER/WP 2/Metacommunity models/Datasetit/Intermittence_Data/Spain_Local_Interm_90_d.csv")

int <- rbind(int1, int2, int3, int4, int5, int6)
# Lasketaan paikkojen keskiarvo kuivan kauden pituudesta 90 vrk ajalla ennen sämpläystä (koska 90 vrk näyttää olevan max arvo useimmissa paikoissa)
names(int)
library(dplyr)
TotDur90<- int%>%group_by(Sites)%>%summarise(TotDur90=mean(TotDur, na.rm = T))
TotDur90$site <- TotDur90$Sites

# Lasketaan paikoille keskiarvo kuivien kausien kestosta 90 vrk ajalla ennen sämpläystä
TotLeng90<- int%>%group_by(Sites)%>%summarise(TotLeng90=mean(TotLeng, na.rm = T))
TotLeng90$site <- TotLeng90$Sites

# Lasketaan ympäristömuuttujien paikkojen keskiarvo, näytteitä 1-6 kpl per paikka
names(env)
discharge<- env%>%group_by(site)%>%summarise(discharge=mean(discharge_l_s, na.rm = T))
moss<- env%>%group_by(site)%>%summarise(moss=mean(moss_cover, na.rm = T))
particlesize <- env%>%group_by(site)%>%summarise(particlesize=mean(particle_size, na.rm = T))


typedrn<- env[,c(2,3,6)]
typedrn <- unique(typedrn)


div$site <- div$Site

all <- merge(div, TotDur90, by="site", all.x=T, sort=F)
all <- merge(all, TotLeng90, by="site", all.x=T, sort=F) 
all <- merge(all, moss, by="site", all.x=T, sort=F) 
all <- merge(all, particlesize, by="site", all.x=T, sort=F)
# all <- merge(all, ConD, by="site", all.x=T, sort=F)
all <- merge(all, typedrn, by="site", all.x=T, sort=F)
# all <- merge(all, ConDmax, by="site", all.x=T, sort=F)
all <- merge(all, discharge, by="site", all.x=T, sort=F)

# data <- merge(yrs, sites[,.(MASTERID_bqe, nyears)], by="MASTERID_bqe", all.x=T, sort=F)

all <- droplevels(all[!all$Site == 'BUK52',]) # täältä puuttuu ainakin intermittence indicators

all <- transform(all, DRN=plyr::revalue(DRN, c("Czech"="Czech Republic")))
colnames(all)[11] <- "Country"


#### Correlation test ####
levels(factor(all$Organism))
ss <- all[Organism == "Sediment_bacteria"] 
# TotDur90
sink("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_bacteria_mean_S_TotDur90_correlations.txt")
paste("-----")
paste("Finland, correlation")
cor.test(ss[Country == "Finland"]$mean_S, ss[Country == "Finland"]$TotDur90)
paste("-----")
paste("Croatia, correlation")
cor.test(ss[Country == "Croatia"]$mean_S, ss[Country == "Croatia"]$TotDur90)
paste("-----")
paste("France, correlation")
cor.test(ss[Country == "France"]$mean_S, ss[Country == "France"]$TotDur90)
paste("-----")
paste("Hungary, correlation")
cor.test(ss[Country == "Hungary"]$mean_S, ss[Country == "Hungary"]$TotDur90)
paste("-----")
paste("Spain, correlation")
cor.test(ss[Country == "Spain"]$mean_S, ss[Country == "Spain"]$TotDur90)
paste("-----")
paste("Czech Republic, correlation")
cor.test(ss[Country == "Czech Republic"]$mean_S, ss[Country == "Czech Republic"]$TotDur90)
sink()

# discharge
sink("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_bacteria_mean_S_Discharge_correlations.txt")
paste("-----")
paste("Finland")
cor.test(ss[Country == "Finland"]$mean_S, ss[Country == "Finland"]$discharge)
paste("-----")
paste("Croatia")
cor.test(ss[Country == "Croatia"]$mean_S, ss[Country == "Croatia"]$discharge)
paste("-----")
paste("Czech Republic")
cor.test(ss[Country == "Czech Republic"]$mean_S, ss[Country == "Czech Republic"]$discharge)
paste("-----")
paste("France")
cor.test(ss[Country == "France"]$mean_S, ss[Country == "France"]$discharge)
paste("-----")
paste("Hungary")
cor.test(ss[Country == "Hungary"]$mean_S, ss[Country == "Hungary"]$discharge)
paste("-----")
paste("Spain")
cor.test(ss[Country == "Spain"]$mean_S, ss[Country == "Spain"]$discharge)
sink()



#### scatterplot, separately organisms and countries ####

#### y-akselilla on alpha diversity, x-akselilla on discharge/totdur90, kullekin eliöryhmälle ja maalle erikseen

levels(factor(all$Organism))
library(ggplot2)

# TÄSTÄ LÄHTEE
ss <- all[Organism == "Sediment_bacteria"] 
# scatterplot, TotDur90 -- linear models with significance as colour
plots2 <- list()
for(i in unique(ss$Country)){
  d <- subset(ss, Country == i)
  p.vals = sapply(unique(d$Country), function(i) {
    coef(summary(lm(mean_S~TotDur90, data=d[Country==i, ])))[2,4] 
  })
  plots2[[paste0(i)]] <- ggplot(d[Organism=="Sediment_bacteria"], aes(x=TotDur90, y=mean_S)) + 
    geom_point(aes(x=TotDur90, y=mean_S), size = 2) + 
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, mean_S), method = "lm", linewidth = 0.5, se = F) + # vaihda kem
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, mean_S), method = "lm", linewidth = 0.5, se = F, linetype="dashed") + # vaihda kem
    # scale_colour_manual(values = c("steelblue","orange")) +
    ggtitle(paste(i)) +
    theme_classic()
}
library(gridExtra)
pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_bacteria_mean_S_vs_TotDur90_lm_sig.pdf"), height=10, width=10)
do.call('grid.arrange', c(plots2))
dev.off()
# scatterplot, discharge
plots3 <- list()
for(i in unique(ss$Country)){
  d <- subset(ss, Country == i)
  p.vals = sapply(unique(d$Country), function(i) {
    coef(summary(lm(mean_S~discharge, data=d[Country==i, ])))[2,4] 
  })
  plots3[[paste0(i)]] <- ggplot(d[Organism=="Sediment_bacteria"], aes(x=discharge, y=mean_S)) + 
    geom_point(aes(x=discharge, y=mean_S), size = 2) + 
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals < 0.05],], aes(discharge, mean_S), method = "lm", linewidth = 0.5, se = F) + # vaihda kem
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals > 0.05],], aes(discharge, mean_S), method = "lm", linewidth = 0.5, se = F, linetype="dashed") + # vaihda kem
    # scale_colour_manual(values = c("steelblue","orange")) +
    ggtitle(paste(i)) +
    theme_classic()
}
library(gridExtra)
pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_bacteria_mean_S_vs_Discharge_lm_sig.pdf"), height=10, width=10)
do.call('grid.arrange', c(plots3))
dev.off()


#### mixed models ####

levels(factor(all$Organism))
library(lme4)
# install.packages("lmerTest")
library(MuMIn)
library(lmerTest)

org <- "Biofilm_diatoms" ### MUOKKAA
sink(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Mixed_models_Mean_S_2_", org, ".txt", sep="")) 
paste(org)
paste("---------------------------------------------------------")
paste("TotDur90 full model")
print(summary(lmer(mean_S ~ TotDur90 + (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("TotDur90 null model")
print(summary(lmer(mean_S ~ (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("Anova, refitting model(s) with ML (instead of REML)")
anova(lmer(mean_S ~ TotDur90 + (1|Country), data=all[Organism==org]), lmer(mean_S ~ (1|Country), data=all[Organism==org])) 
paste("---------------------------------------------------------")
paste("TotLeng90 full model")
print(summary(lmer(mean_S ~ TotLeng90 + (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("TotLeng90 null model")
print(summary(lmer(mean_S ~ (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("Anova, refitting model(s) with ML (instead of REML)")
anova(lmer(mean_S ~ TotLeng90 + (1|Country), data=all[Organism==org]), lmer(mean_S ~ (1|Country), data=all[Organism==org])) 
sink()

#### scatter plot, all countries in the same plot ####

levels(factor(all$Country))

min(all$mean_S)
max(all$mean_S)
library(ggplot2)

levels(factor(all$Organism))

ss <- all[Organism=="Sediment_fungi"] # MUOKKAA


# testataan vielä voiko piirtää mustan viivan lm:n p-arvon perusteella, vaikka mixed model ko. muuttujalla olisi parempi kuin nollamalli
summary(lm(ss$mean_S~ ss$TotDur90))
# TotDur
p.vals = sapply(unique(ss$Country), function(i) {
  coef(summary(lm(mean_S~TotDur90, data=ss[Country==i, ])))[2,4] 
})
p.vals
p1 <- ggplot(ss) + 
  geom_point(aes(TotDur90, mean_S, colour=Country), size = 1, alpha=0.5) +
  geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, mean_S, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
  geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, mean_S, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
  geom_smooth(aes(TotDur90, mean_S), colour="black", method = "lm", linewidth = 1.1, se = F) + 
  theme_classic() +
  ggtitle("Sediment fungi") + # MUOKKAA
  xlab("Number of dry days") +
  ylab("Taxonomic richness") +
  scale_color_manual(values = c("Croatia" = "#ef476f",
                                "Czech Republic" = "#f78c6b", 
                                "Finland" = "#ffd166", 
                                "France" = "#06d6a0", 
                                "Hungary" = "#118ab2",
                                "Spain" = "#073b4c")) +
  theme(legend.position = "none")
p1


pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/All_Sediment_fungi_mean_S_vs_TotDur90_lm2.pdf"), height=3, width=4) # MUOKKAA
p1
dev.off()


# TotLeng90
summary(lm(ss$mean_S~ ss$TotLeng90))

p.vals = sapply(unique(ss$Country), function(i) {
  coef(summary(lm(mean_S~TotLeng90, data=ss[Country==i, ])))[2,4] 
})
p.vals
p1 <- ggplot(ss) + 
  geom_point(aes(TotLeng90, mean_S, colour=Country), size = 1, alpha=0.5) + 
  geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotLeng90, mean_S, colour=Country), method = "lm", linewidth = 0.75, se = F) + 
  geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotLeng90, mean_S, colour=Country), method = "lm", linewidth = 0.75, se = F, linetype="dashed") +
  geom_smooth(aes(TotLeng90, mean_S), colour="black", method = "lm", linewidth = 1.1, se = F) + 
  theme_classic() +
  ggtitle("Sediment fungi") + # MUOKKAA
  xlab("Duration of dry periods") +
  ylab("Taxonomic richness") +
  scale_color_manual(values = c("Croatia" = "#ef476f",
                                "Czech Republic" = "#f78c6b", 
                                "Finland" = "#ffd166", 
                                "France" = "#06d6a0", 
                                "Hungary" = "#118ab2",
                                "Spain" = "#073b4c")) +
  theme(legend.position = "none")
p1


pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/All_Sediment_fungi_mean_S_vs_TotLeng90_lm2.pdf"), height=3, width=4) # MUOKKAA
p1
dev.off()


# 
# # discharge
# p.vals = sapply(unique(ss$Country), function(i) {
#   coef(summary(lm(mean_S~discharge, data=ss[Country==i, ])))[2,4] 
# })
# p.vals
# p1 <- ggplot(ss) + 
#   geom_point(aes(discharge, mean_S, colour=Country), size = 2, alpha=0.5) + 
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(discharge, mean_S, colour=Country), method = "lm", linewidth = 0.5, se = F) + 
#   geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(discharge, mean_S, colour=Country), method = "lm", linewidth = 0.5, se = F, linetype="dashed") +
#   theme_classic() +
#   ggtitle("Sediment_bacteria") + # MUOKKAA
#   scale_color_manual(values = c("Croatia" = "cornflowerblue",
#                                 "Czech Republic" = "mediumseagreen", 
#                                 "Finland" = "darkblue", 
#                                 "France" = "bisque3", 
#                                 "Hungary" = "orchid4",
#                                 "Spain" = "sienna3"))
# p1
# 
# 
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/All_Sediment_bacteria_mean_S_vs_discharge_lm.pdf"), height=3, width=5) # MUOKKAA
# p1
# dev.off()
