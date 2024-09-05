# DRYvER 

# Site-specific  assembly processes
# - How does drying affect local assembly processes? (cf. Huttunen et al. 2017)

rm(list = ls())
library(data.table)
library(vegan)

# closeAllConnections()


res <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Results/null_model_results_reduced.csv")
env <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/ENV_all_NAs_as_blank.csv")
# int <- fread("C:/Users/vilmi/Desktop/DRYVERR/2024/Data/Intermittence_all.csv")

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

TotLeng90<- int%>%group_by(Sites)%>%summarise(TotLeng90=mean(TotLeng, na.rm = T))
TotLeng90$site <- TotLeng90$Sites

# Lasketaan ympäristömuuttujien paikkojen keskiarvo, näytteitä 1-6 kpl per paikka
names(env)
library(dplyr)
discharge<- env%>%group_by(site)%>%summarise(discharge=mean(discharge_l_s, na.rm = T))
moss<- env%>%group_by(site)%>%summarise(moss=mean(moss_cover, na.rm = T))
particlesize <- env%>%group_by(site)%>%summarise(particlesize=mean(particle_size, na.rm = T))

# names(int)
# ConD<- int%>%group_by(Sites)%>%summarise(ConD=mean(ConD, na.rm = T))
# ConD$site <- ConD$Sites
# 
# ConDmax<- int%>%group_by(Sites)%>%summarise(ConDmax=max(ConD, na.rm = T))
# ConDmax$site <- ConD$Sites


type<- env[,c(3,6)]
type <- unique(type)
res$site <- res$Site

all <- merge(res, TotDur90, by="site", all.x=T, sort=F)
all <- merge(all, TotLeng90, by="site", all.x=T, sort=F) 
all <- merge(all, moss, by="site", all.x=T, sort=F) 
all <- merge(all, particlesize, by="site", all.x=T, sort=F)
# all <- merge(all, ConD, by="site", all.x=T, sort=F)
all <- merge(all, type, by="site", all.x=T, sort=F)
# all <- merge(all, ConDmax, by="site", all.x=T, sort=F)
all <- merge(all, discharge, by="site", all.x=T, sort=F)


all <- droplevels(all[!all$Site == 'BUK52',]) # täältä puuttuu ainakin intermittence indicators


#### correlation test ####
# correlation test, countries
# ss <- all[Country == "Finland"] 
# cor.test(all[Country == "Finland"]$SES, all[Country == "Finland"]$ConDmax)

levels(factor(all$Country))
levels(factor(all$Organism))

ss <- all[Organism == "Sediment_fungi"] 


# TotDur90
sink("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_fungi_SES_TotDur90_correlations.txt")
paste("-----")
paste("Finland, correlation")
cor.test(ss[Country == "Finland"]$SES, ss[Country == "Finland"]$TotDur90)
paste("-----")
paste("Croatia, correlation")
cor.test(ss[Country == "Croatia"]$SES, ss[Country == "Croatia"]$TotDur90)
paste("-----")
paste("France, correlation")
cor.test(ss[Country == "France"]$SES, ss[Country == "France"]$TotDur90)
paste("-----")
paste("Hungary, correlation")
cor.test(ss[Country == "Hungary"]$SES, ss[Country == "Hungary"]$TotDur90)
paste("-----")
paste("Spain, correlation")
cor.test(ss[Country == "Spain"]$SES, ss[Country == "Spain"]$TotDur90)
paste("-----")
paste("Czech, correlation")
cor.test(ss[Country == "Czech Republic"]$SES, ss[Country == "Czech Republic"]$TotDur90)
sink()

# discharge
sink("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_fungi_SES_Discharge_correlations.txt")
paste("-----")
paste("Finland")
cor.test(ss[Country == "Finland"]$SES, ss[Country == "Finland"]$discharge)
paste("-----")
paste("Croatia")
cor.test(ss[Country == "Croatia"]$SES, ss[Country == "Croatia"]$discharge)
paste("-----")
paste("Czech")
cor.test(ss[Country == "Czech Republic"]$SES, ss[Country == "Czech Republic"]$discharge)
paste("-----")
paste("France")
cor.test(ss[Country == "France"]$SES, ss[Country == "France"]$discharge)
paste("-----")
paste("Hungary")
cor.test(ss[Country == "Hungary"]$SES, ss[Country == "Hungary"]$discharge)
paste("-----")
paste("Spain")
cor.test(ss[Country == "Spain"]$SES, ss[Country == "Spain"]$discharge)
sink()

# # boxplot for a country
# ss <- all[Country == "Finland"]
# p <- ggplot(ss, aes(x=SES, y=site, color=stream_type)) +
#   geom_boxplot() + coord_flip()
# p
# 
# pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/FIN_biof_fungi_SES_vs_streamtype.pdf"), height=6.4, width=10) # max height=9.65, width=6.4
# p
# dev.off()

# all in figure
library(ggplot2)

levels(factor(all$Organism))

# tehty 18032024

#### boxplot: SES vs. stream type ####
plots <- list()
for(i in unique(all$Country)){
  d <- subset(all, Country == i)
  plots[[paste0(i)]] <- ggplot(d[Organism=="Macroinvertebrates"], aes(x=SES, y=site, color=stream_type)) + 
    scale_colour_manual(values = c("steelblue","orange")) +
    geom_boxplot() + coord_flip() +
    ggtitle(paste(i)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

library(gridExtra)
pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Macroinvertebrates_SES_vs_streamtype.pdf"), height=10, width=10)
do.call('grid.arrange', c(plots))
dev.off()





#### scatterplot, separately organisms and countries ####

#### TÄMÄ KESKEN 19033024: tee niin että p-arvo on väri, ja x-akselilla on discharge/totdur90, kullekin eliöryhmälle ja maalle erikseen

# how to make new column based on another column's values
library(dplyr)
levels <- c(-Inf, 0.05,  Inf)
labels <- c("sig", "nonsig")
all2 <- all %>% mutate(Significance2 = cut(Significance, levels, labels = labels))

levels(factor(all$Organism))
library(ggplot2)


#### scatter plot, all countries in the same plot ####

levels(factor(all$Country))
levels(factor(all$Organism))
# min(all$mean_S)
# max(all$mean_S)








# tähän jäi pe 12.4.2024 eli miten saa viimeisen rivin geom pointin laittamaan reunat vain merkitseville pisteille?







ss <- all2[Organism == "Flying_macroinvertebrates"] 
library(ggplot2)
ss2 <- ss[Significance2 == "sig"]
# TotDur
p.vals = sapply(unique(ss$Country), function(i) {
  coef(summary(lm(SES~TotDur90, data=ss[Country==i, ])))[2,4] 
})
p.vals
p1 <- ggplot() + 
  geom_point(aes(TotDur90, SES, colour=Country), data = ss, size = 2.5, alpha=0.3) +
  geom_point(aes(TotDur90, SES, colour=Country), data = ss2, size = 2.5, alpha=0.9) +
  geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, SES, colour=Country), method = "lm", linewidth = 0.8, se = F) + 
  geom_smooth(data=ss[ss$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, SES, colour=Country), method = "lm", linewidth = 0.8, se = F, linetype="dashed") +
  theme_classic() +
  ggtitle("Flying macroinvertebrates") + # MUOKKAA
  scale_color_manual(values = c("Croatia" = "#ef476f",
                                "Czech Republic" = "#f78c6b", 
                                "Finland" = "#ffd166", 
                                "France" = "#06d6a0", 
                                "Hungary" = "#118ab2",
                                "Spain" = "#073b4c"))
  # geom_point(aes(TotDur90, SES, data=ss[ss$Signicance2=="sig"]), fill=NA, colour = "black", size = 2.5, alpha=0.0)
p1


pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/All_flying_macroninvs_SES_vs_TotDur90_lm.pdf"), height=3, width=5) # MUOKKAA
p1
dev.off()

# TÄSTÄ LÄHTEE
ss <- all2[Organism == "Flying_macroinvertebrates"] 
# scatterplot, TotDur90 -- linear models with significance as colour
plots2 <- list()
for(i in unique(ss$Country)){
  d <- subset(ss, Country == i)
  p.vals = sapply(unique(d$Country), function(i) {
    coef(summary(lm(SES~TotDur90, data=d[Country==i, ])))[2,4] 
  })
  plots2[[paste0(i)]] <- ggplot(d[Organism=="Biofilm_diatoms"], aes(x=TotDur90, y=SES)) + 
    geom_point(aes(x=TotDur90, y=SES, colour = Significance2), size = 2) + 
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, SES), method = "lm", linewidth = 0.5, se = F) + # vaihda kem
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, SES), method = "lm", linewidth = 0.5, se = F, linetype="dashed") + # vaihda kem
    scale_colour_manual(values = c("steelblue","orange")) +
    ggtitle(paste(i)) +
    theme_classic()
}
library(gridExtra)
pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Biofilm_diatoms_SES_vs_TotDur90_lm_sig.pdf"), height=10, width=10)
do.call('grid.arrange', c(plots2))
dev.off()

# scatterplot, discharge -- linear models with significance as colour
levels(factor(all$Organism))
ss <- all2[Organism == "Sediment_fungi"] 
plots3 <- list()
for(i in unique(ss$Country)){
  d <- subset(ss, Country == i)
  p.vals = sapply(unique(d$Country), function(i) {
    coef(summary(lm(SES~discharge, data=d[Country==i, ])))[2,4] 
  })
  plots3[[paste0(i)]] <- ggplot(d[Organism=="Sediment_fungi"], aes(x=discharge, y=SES)) + 
    geom_point(aes(x=discharge, y=SES, colour = Significance2), size = 2) + 
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals < 0.05],], aes(discharge, SES), method = "lm", linewidth = 0.5, se = F) + # vaihda kem
    geom_smooth(data=d[d$Country %in% names(p.vals)[p.vals > 0.05],], aes(discharge, SES), method = "lm", linewidth = 0.5, se = F, linetype="dashed") + # vaihda kem
    scale_colour_manual(values = c("steelblue","orange")) +
    ggtitle(paste(i)) +
    theme_classic()
}
library(gridExtra)
pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Sediment_fungi_SES_vs_Discharge_lm_sig.pdf"), height=10, width=10)
do.call('grid.arrange', c(plots3))
dev.off()

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

# scatterplot, all organisms in the same plot
levels(factor(all$Country))
levels(factor(all$Organism))
min(all$SES)
max(all$SES)
ss <- all[Country=="Finland"] # MUOKKAA

p.vals = sapply(unique(ss$Organism), function(i) {
    coef(summary(lm(SES~TotDur90, data=ss[Organism==i, ])))[2,4] 
  })
p.vals
p1 <- ggplot(ss) + 
  # geom_point(size = 2, alpha=0.5) + 
  geom_smooth(data=ss[ss$Organism %in% names(p.vals)[p.vals < 0.05],], aes(TotDur90, SES, colour=Organism), method = "lm", linewidth = 0.5, se = F) + 
  geom_smooth(data=ss[ss$Organism %in% names(p.vals)[p.vals > 0.05],], aes(TotDur90, SES, colour=Organism), method = "lm", linewidth = 0.5, se = F, linetype="dashed") +
  theme_classic() +
  ggtitle("Finland") + # MUOKKAA
  scale_color_manual(values = c("Biofilm_fungi" = "cornflowerblue",
                                "Macroinvertebrates" = "mediumseagreen", 
                                "Sediment_fungi" = "darkblue", 
                                "Biofilm_diatoms" = "bisque3", 
                                "Sediment_diatoms" = "rosybrown4"))
p1


pdf(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Finland_SES_vs_TotDur90_lm.pdf"), height=3, width=6) # MUOKKAA
p1
dev.off()

#### mixed models ####

levels(factor(all$Organism))
library(lme4)
# install.packages("lmerTest")
library(MuMIn)
library(lmerTest)

org <- "Nonflying_macroinvertebrates" ### MUOKKAA
sink(paste("M:/gDRYvER/WP 2/Metacommunity models/Tulokset/Null models/Mixed_models_SES_", org, ".txt", sep="")) 
paste(org)
paste("---------------------------------------------------------")
paste("TotDur90 full model")
print(summary(lmer(SES ~ TotDur90 + (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("TotDur90 null model")
print(summary(lmer(SES ~ (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("Anova, refitting model(s) with ML (instead of REML)")
anova(lmer(SES ~ TotDur90 + (1|Country), data=all[Organism==org]), lmer(SES ~ (1|Country), data=all[Organism==org])) 
paste("---------------------------------------------------------")
paste("discharge full model")
print(summary(lmer(SES ~ discharge + (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("discharge null model")
print(summary(lmer(SES ~ (1|Country), data=all[Organism==org]))) 
paste("---------------------------------------------------------")
paste("Anova, refitting model(s) with ML (instead of REML)")
anova(lmer(SES ~ discharge + (1|Country), data=all[Organism==org]), lmer(SES ~ (1|Country), data=all[Organism==org])) 
sink()
