# SUPPLEMENTARY R-SCRIPT TO:
# Long-term data reveal: Recent declines in carabid beetles in a temperate forest are linked to severe droughts
# F.Weiss, H.von Wehrden & A.Linde
# F.Weiss: ORCID 0000-0003-1078-1528

# PART III: Modelling diversity

# Data available at: 

#### Package list ####

library(data.table)
library(vegan)
library(dplyr)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(lme4)
library(ggeffects)
library(gridExtra)
library(sjPlot)
library(gamm4)
library(cols4all)
library(bbmle)


colors <- c4a("berlin", 7)


#### Loading data ####

# loading processed data from previous script
beetle_samples <- read.csv("carabid_samples_raw2.csv")

# loading sampling meta data
sampling_meta <- read.csv2("sampling_meta_dwd_spei.csv")


#### Preparing data ####

# remove plots with sampling effort != 4
beetle_samples$sampling_length <- sampling_meta$Trapping_length[match(beetle_samples$sample_id, sampling_meta$sample_id)]
beetle_samples$sampling_effort <- sampling_meta$sampling_effort[match(beetle_samples$sample_id, sampling_meta$sample_id)]

# remove all samples of intervals with more or less than 4 traps
beetle_samples <- beetle_samples[beetle_samples$sampling_effort == 4,]

# update sample_id2
beetle_samples$sample_id2 <- paste(beetle_samples$year, beetle_samples$plot, sep = "_")

# check which traps have full 3 intervals!
data_table <- data.table(beetle_samples)   
consistency <- aggregate(cbind(beetle_samples$interval), by= list(beetle_samples$sample_id2), FUN=  function(x){length(unique(x))}) 
incomplete_samples <- consistency[consistency$V1 < 3,]
beetle_samples2 <- beetle_samples[!(beetle_samples$sample_id2 %in% c(incomplete_samples$Group.1)), ]

# plot data availability 
ggplot(beetle_samples2, aes(year, plot)) +
  geom_point(shape = 15, color = colors[2], size = 3)+
  ylab("Plot")+
  xlab("Year")+
  ggtitle("Data availability for diversity metrics")+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20), 
        axis.text.x  = element_text(size=10), axis.text.y  = element_text(size=20) )


# aggregate per traps over all 3 intervals for all 4 traps
beetle_samples2$sample_id2 <- paste(beetle_samples2$year, beetle_samples2$plot, sep = "_")

# aggregate by sampleID (year and plot)
beetle_samples3 <- aggregate(cbind(beetle_samples2$sampling_abundance, beetle_samples2$corr_abundance2), by= list( beetle_samples2$sample_id2, beetle_samples2$year, beetle_samples2$plot, beetle_samples2$site ,beetle_samples2$species), FUN=sum  )

beetle_samples3<-beetle_samples3 %>% 
  rename(
    sample_id2 = Group.1,
    year = Group.2,
    plot= Group.3,
    site= Group.4,
    species=Group.5,
    sampling_abundance = V1,
    corr_abundance = V2
  )

# remove zeros
beetle_samples3 <- beetle_samples3[beetle_samples3$species != "no carabids",]

# reformat to wide table
# based on corrected abundances
wide_samples <- beetle_samples3[, c(1,5,7)]

wide_samples  <- reshape(wide_samples, idvar = "sample_id2", timevar = "species", direction = "wide")
wide_samples[is.na(wide_samples)] <- 0

wide_samples2 <- wide_samples[,-1] 
row.names(wide_samples2) <- wide_samples$sample_id2




#### Compute diversity metrics ####


# observed species richness
diversity_samples <- aggregate(cbind(beetle_samples3$species), by= list(beetle_samples3$sample_id2, beetle_samples3$year, beetle_samples3$plot, beetle_samples3$site), function(x) length(unique(x)))

diversity_samples<-diversity_samples %>% 
  rename(
    sample_id2 = Group.1,
    year = Group.2,
    plot = Group.3,
    site = Group.4,
    obs_species = V1
  )


# Pielou's Evenness 
shannon <- diversity(wide_samples2, index="shannon")
shannon <- data.frame(sample_id2 = wide_samples$sample_id2, shannon = shannon)
diversity_samples$shannon <- shannon$shannon[match(diversity_samples$sample_id2, shannon$sample_id2)]
diversity_samples$evenness <- diversity_samples$shannon/log(diversity_samples$obs_species)


# turnover (Jaccard Index)

# exclude plots with <= 4 years of data
exclude <- c("2b", "09", "05", "04")

beetle_samples4 <- beetle_samples3[!(beetle_samples3$plot %in% exclude),]


# calculate jaccard as mean dissimilarity to first 2 years of each plot
# set up function for jaccard index
my_jaccard <- function(t1, t2){
  length(intersect(t1,t2)) / length(unique(c(t1,t2)))}

# set up data frame to fill 
turnover_data <- data.frame(year = as.numeric(c()), plot = as.factor(c()),jaccard = as.numeric(c())  )


# run nested loop for years in plots  
plots <- c(unique(beetle_samples4$plot))

for(i in plots){
  
  loopdata <- beetle_samples4[beetle_samples4$plot == i,]
  
  years1 <- c(sort(unique(loopdata$year))) 
  
  ref1<- c(unique(loopdata$species[loopdata$year == min(years1)]))
  
  years2 <- years1[years1 != min(years1)] 
  
  ref2<- c(unique(loopdata$species[loopdata$year == min(years2)]))
  
  years3 <- years2[years2 != min(years2)] 
  
  jacc <- data.frame(year= years3, plot = rep(i, length(years3), jaccard = rep(NA, length(years3))))
  
  
  for(j in years3){
    
    species <- c(unique(loopdata$species[loopdata$year == j]))
    
    jaccard1 <- my_jaccard (ref1, species)
    jaccard2 <- my_jaccard (ref2, species)
    
    jacc$jaccard[jacc$year == j] <- mean(jaccard1, jaccard2)
    
  }
  
  turnover_data <- rbind(turnover_data, jacc)
  
}

# merge diversity metrics in one dataframe
turnover_data$sample_id2 <- paste(turnover_data$year, turnover_data$plot, sep = "_")
diversity_samples$jaccard <- turnover_data$jaccard[match(diversity_samples$sample_id2, turnover_data$sample_id2)]

# merge with SPEI data
sampling_meta2 <- unique(sampling_meta[,c(2,13:30)])
diversity_samples2 <- left_join(diversity_samples, sampling_meta2, by= "year")

LTdiversity <- diversity_samples2

# year variable
LTdiversity$year <- as.integer(LTdiversity$year)
LTdiversity$year2 <- as.factor(LTdiversity$year)
LTdiversity$scaled_year <- c(scale(LTdiversity$year))

LTdiversity$plot <- as.factor(LTdiversity$plot)
LTdiversity$site <-as.factor( LTdiversity$site )



#### Species richness #### 

# GLMM
species_glmm <- glmmTMB(obs_species ~ scaled_year + (1|year2) + (1|site/plot) , family = poisson, data=LTdiversity)

dharma_sim1 <- simulateResiduals(fittedModel = species_glmm, re.form= NULL)
plot(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(species_glmm, dv.labels = c("species richness (observed) GLMM"))



## gamm
species_gamm <- gamm4(obs_species ~ s(scaled_year),
                      random=  ~ (1|site/plot), 
                      data= LTdiversity, 
                      family=poisson,
                      REML=TRUE)

plot(species_gamm$gam)
# GAMM estimates linear trend -> discard GAMM

# SPEI GLMM
m1 <- update(species_glmm,. ~ . + spei12)
m2<- update(species_glmm,. ~ . + spei12_1)
m3 <- update(species_glmm,. ~ . + spei12_2)
m4 <- update(species_glmm,. ~ . + spei24)
m5 <- update(species_glmm,. ~ . + spei24_1)
m6 <- update(species_glmm,. ~ . + spei24_2)
m7 <- update(species_glmm,. ~ . + spei36)
m8 <- update(species_glmm,. ~ . + spei36_1)
m9 <- update(species_glmm,. ~ . + spei36_2)
m10 <- update(species_glmm,. ~ . + spei48)
m11 <- update(species_glmm,. ~ . + spei48_1)
m12 <- update(species_glmm,. ~ . + spei48_2)
m13 <- update(species_glmm,. ~ . + spei60)
m14 <- update(species_glmm,. ~ . + spei60_1)
m15 <- update(species_glmm,. ~ . + spei60_2)
m16 <- update(species_glmm,. ~ . + spei72)
m17 <- update(species_glmm,. ~ . + spei72_1)
m18 <- update(species_glmm,. ~ . + spei72_2)

AICtab(species_glmm,m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18)

# SPEI not a meaningful addition



#### Evenness ####

# values between 0 and 1 -> beta distribution? check Geissinger et al. 2022

# GLMM
evenness_glmm <- glmmTMB(evenness ~ scaled_year + (1|year2) + (1|site/plot), family = beta_family(link="logit"), data=LTdiversity)

dharma_sim1 <- simulateResiduals(fittedModel = evenness_glmm, re.form= NULL)
plot(dharma_sim1)
dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(evenness_glmm, dv.labels = "(Pilou's) Evenness GLMM")

# No GAMM -> not beta distribution in gamm4!

# SPEI GLMM
m1 <- update(evenness_glmm,. ~ . + spei12)
m2 <- update(evenness_glmm,. ~ . + spei12_1)
m3 <- update(evenness_glmm,. ~ . + spei12_2)
m4 <- update(evenness_glmm,. ~ . + spei24)
m5 <- update(evenness_glmm,. ~ . + spei24_1)
m6 <- update(evenness_glmm,. ~ . + spei24_2)
m7 <- update(evenness_glmm,. ~ . + spei36)
m8 <- update(evenness_glmm,. ~ . + spei36_1)
m9 <- update(evenness_glmm,. ~ . + spei36_2)
m10 <- update(evenness_glmm,. ~ . + spei48)
m11 <- update(evenness_glmm,. ~ . + spei48_1)
m12 <- update(evenness_glmm,. ~ . + spei48_2)
m13 <- update(evenness_glmm,. ~ . + spei60)
m14 <- update(evenness_glmm,. ~ . + spei60_1)
m15 <- update(evenness_glmm,. ~ . + spei60_2)
m16 <- update(evenness_glmm,. ~ . + spei72)
m17 <- update(evenness_glmm,. ~ . + spei72_1)
m18 <- update(evenness_glmm,. ~ . + spei72_2)

AICtab(evenness_glmm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18)

# no meaningfull addition



#### Turnover/jaccard ####

# reorganize data
turnover_data <- LTdiversity[!(is.na(LTdiversity$jaccard) == TRUE),]

hist(turnover_data$jaccard)
unique(turnover_data$plot)

plot( turnover_data$year, turnover_data$jaccard)

# values between 0 and 1 -> beta distribution? check Geissinger et al. 2022

jaccard_glmm <- glmmTMB(jaccard ~ scaled_year + (1|year2) + (1|site/plot), family = beta_family(link="logit"), data=turnover_data)

dharma_sim1 <- simulateResiduals(fittedModel = jaccard_glmm, re.form= NULL)
plot(dharma_sim1)
dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))
# fine!

tab_model(jaccard_glmm,dv.labels = "Turnover (Jaccard index) GLMM")

summary(jaccard_glmm)

## SPEI glmm


m1 <- update(jaccard_glmm,. ~ . + spei12)
m2- update(jaccard_glmm,. ~ . + spei12_1)
m3 <- update(jaccard_glmm,. ~ . + spei12_2)
m4 <- update(jaccard_glmm,. ~ . + spei24)
m5 <- update(jaccard_glmm,. ~ . + spei24_1)
m6 <- update(jaccard_glmm,. ~ . + spei24_2)
m7 <- update(jaccard_glmm,. ~ . + spei36)
m8 <- update(jaccard_glmm,. ~ . + spei36_1)
m9 <- update(jaccard_glmm,. ~ . + spei36_2)
m10 <- update(jaccard_glmm,. ~ . + spei48)
m11 <- update(jaccard_glmm,. ~ . + spei48_1)
m12 <- update(jaccard_glmm,. ~ . + spei48_2)
m13 <- update(jaccard_glmm,. ~ . + spei60)
m14 <- update(jaccard_glmm,. ~ . + spei60_1)
m15 <- update(jaccard_glmm,. ~ . + spei60_2)
m16 <- update(jaccard_glmm,. ~ . + spei72)
m17 <- update(jaccard_glmm,. ~ . + spei72_1)
m18 <- update(jaccard_glmm,. ~ . + spei72_2)

AICtab(jaccard_glmm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18)





#### Figure ####



# species richness 
ggpreds1 <- ggpredict(species_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds1$year <- sort(unique(LTdiversity$year))  

preds_richness <- data.frame(  scaled_year= ggpreds1$x,
                               year = ggpreds1$year,
                               fit  =  ggpreds1$predicted,
                               lwr   = ggpreds1$conf.low,
                               upr   = ggpreds1$conf.high   )

tsdata <- ts(data.frame(x1=c(preds_richness$fit)), start=c(2001), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)


plot1 <- ggplot(preds_richness, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "Species (observed)", x = "Year")+
  ggtitle("A")+
  theme(plot.title = element_text(hjust =0.05,vjust=-2 ))+
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  geom_jitter(data=LTdiversity ,aes(scaled_year, obs_species),alpha=0.3, color = colors[4], size=2)+
  geom_ribbon(data = preds_richness, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[4])+
  geom_line(data=preds_richness ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=2)+
  theme(legend.position="none")


# evenness 
ggpreds1 <- ggpredict(evenness_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds1$year <- sort(unique(LTdiversity$year))  

preds_evenness <- data.frame(  scaled_year= ggpreds1$x,
                               year = ggpreds1$year,
                               fit  =  ggpreds1$predicted,
                               lwr   = ggpreds1$conf.low,
                               upr   = ggpreds1$conf.high   )

tsdata <- ts(data.frame(x1=c(preds_evenness$fit)), start=c(2001), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)

plot2 <- ggplot(preds_evenness, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "Pielou's evenness", x = "Year")+
  ggtitle("B")+
  theme(plot.title = element_text(hjust =0.05,vjust=-2))+
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_evenness$scaled_year)) )+
  geom_jitter(data=LTdiversity ,aes(scaled_year, evenness),alpha=0.3, color = colors[4], size=2)+
  geom_ribbon(data = preds_evenness, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[2])+
  geom_line(data=preds_evenness ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=2)+
  theme(legend.position="none")


# Turnover 
ggpreds1 <- ggpredict(jaccard_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds1$year <- sort(unique(turnover_data$year))  

preds_jaccard <- data.frame(  scaled_year= ggpreds1$x,
                              year = ggpreds1$year,
                              fit  =  ggpreds1$predicted,
                              lwr   = ggpreds1$conf.low,
                              upr   = ggpreds1$conf.high   )

tsdata <- ts(data.frame(x1=c(preds_jaccard$fit)), start=c(2003), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)

plot3 <- ggplot(preds_jaccard, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "Jaccard similarity", x = "Year")+
  ggtitle("C")+
  theme(plot.title = element_text(hjust =0.05,vjust=-2))+
  scale_x_continuous(labels= c(  "2003", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_jaccard$scaled_year)) )+
  geom_jitter(data=turnover_data ,aes(scaled_year, jaccard),alpha=0.3, color = colors[4], size=2)+
  geom_ribbon(data = preds_jaccard, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[6])+
  geom_line(data=preds_jaccard ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=2)+
  theme(legend.position="none")

grid.arrange(plot1, plot3, plot2, ncol=3)


# scaled plot 

# load abundance GLMM (corr_ab_glmm) from previous script
ggpreds1 <- ggpredict(corr_ab_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_length = -0.1443603, interval= "2", scaled_sampling_rain=-0.2290402,scaled_sampling_temp= 0.1761928 ))
ggpreds1$year <- unique(LTabundance$year2)  

preds_abundance <- data.frame(  scaled_year= ggpreds1$x,
                                year = ggpreds1$year,
                                fit  =  ggpreds1$predicted,
                                lwr   = ggpreds1$conf.low,
                                upr   = ggpreds1$conf.high   )


preds_abundance$scaled_fit <- preds_abundance$fit / max(preds_abundance$fit)
preds_richness$scaled_fit <- preds_richness$fit / max(preds_richness$fit)
preds_evenness$scaled_fit <- preds_evenness$fit / max(preds_evenness$fit)
preds_jaccard$scaled_fit <- preds_jaccard$fit / max(preds_jaccard$fit)


plot4 <- ggplot(preds_abundance, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "", x = "Time")+
  ggtitle("D")+
  theme(plot.title = element_text(hjust =0.05,vjust=-2))+
  scale_y_continuous(limits = c(0.5, 1.05), breaks = c(0.5, 1), labels= c("low", "high"))+
  scale_x_continuous(limits = c(-1.66,3.2), labels= c(  ), breaks = c() )+
  geom_line(data=preds_abundance ,aes(scaled_year, scaled_fit), color = colors[4],alpha=0.5, lty=2, lwd=3)+
  geom_line(data=preds_evenness ,aes(scaled_year, scaled_fit), color = colors[2], lty=1, lwd=3)+
  geom_line(data=preds_richness ,aes(scaled_year, scaled_fit), color = colors[4], lty=1, lwd=3)+
  geom_line(data=preds_jaccard ,aes(scaled_year, scaled_fit), color = colors[6], lty=1, lwd=3)+
  annotate("text", x = 2.5, y = 0.9, label = "Evenness", color = colors[2], size=6, fontface =2)+
  annotate("text", x = 2.5, y = 0.79, label = "Richness", color = colors[4], size=6, fontface =2)+
  annotate("text", x = 2.5, y = 0.74, label = "Turnover", color = colors[6], size=6, fontface =2)+
  annotate("text", x = 2.65, y = 0.69, label = "Abundance", color = colors[4],alpha=0.5, size=6, fontface =2)+
  annotate("text", x = 2.65, y = 0.66, label = "(n.s.)", color = colors[4],alpha=0.5, size=6, fontface =2)+
  theme(legend.position="none")

grid.arrange(plot1, plot2, plot3, plot4, ncol=4, widths = c(1,1,1,1.2))


#### sensitivity ####

# species richness 

coefs <- coef(summary(species_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[2,1],
                            lwr   = estimates[2,1] + 1.96*estimates[2,2],
                            upr   = estimates[2,1] - 1.96*estimates[2,2],
                            p     = estimates[2,4])

# OOS years
years <- c(sort(unique(LTdiversity$year)))

oosyear_species_GLMM <- data.frame(excluded_year= c( 2001:2022))

for(i in years){
  loop_data <- LTdiversity[LTdiversity$year != i,]
  loop_mod<- glmmTMB(obs_species ~ scaled_year + (1|year2) + (1|site/plot) , family = poisson, data=loop_data)
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  oosyear_species_GLMM$estimate[oosyear_species_GLMM$excluded_year == i] <- estimates[2,1]
  oosyear_species_GLMM$estimate_upr[oosyear_species_GLMM$excluded_year == i] <- estimates[2,1] + 1.96*estimates[2,2]
  oosyear_species_GLMM$estimate_lwr[oosyear_species_GLMM$excluded_year == i] <- estimates[2,1] - 1.96*estimates[2,2]
  oosyear_species_GLMM$p[oosyear_species_GLMM$excluded_year == i] <- estimates[2,4]
  print(i)
}

oosyear_species_GLMM_SIG <- oosyear_species_GLMM[oosyear_species_GLMM$p < 0.05,]


plot1 <- ggplot(oosyear_species_GLMM, aes(excluded_year, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "Trend coefficient", x = "Excluded year")+
  ggtitle("Species richness trend sensitivity")+
  theme(plot.title = element_text(hjust = 0.05,vjust=0))+  
  scale_x_continuous(labels= c(  "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(2001:2022) )+
  geom_hline(yintercept = full_coefs$effect, color = berlin[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = berlin[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = berlin[6], lty=2, lwd=1.5)+
  geom_errorbar(data=oosyear_species_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = berlin[4])+
  geom_point(data=oosyear_species_GLMM ,aes(excluded_year, estimate, color = berlin[4]),pch=20 ,size=5) +
  geom_point(data=oosyear_species_GLMM ,aes(excluded_year, estimate, color = berlin[2]),pch=20 ,size=4) +
  geom_point(data=oosyear_species_GLMM_SIG ,aes(excluded_year, estimate, color = berlin[4]),pch=8 ,size=4) +
  scale_color_manual(values = c(berlin[4], berlin[2]), name = "random year effect", labels=c("without", "with" ))+
  theme(legend.position="none")


# OOS plot
plots <- unique(LTdiversity$plot)

oosplot_species_GLMM <- data.frame(excluded_plot= as.factor(unique(LTdiversity$plot)))

for(i in plots){
  loop_data <- LTdiversity[LTdiversity$plot != i,]
  
  loop_mod <- glmmTMB(obs_species ~ scaled_year + (1|year2) + (1|site/plot) , family = poisson, data=loop_data)
  
  coefs<-coef(summary(loop_mod))
  
  estimates <- coefs$cond
  
  oosplot_species_GLMM$estimate[oosplot_species_GLMM$excluded_plot == i] <- estimates[2,1]
  
  oosplot_species_GLMM$estimate_upr[oosplot_species_GLMM$excluded_plot == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosplot_species_GLMM$estimate_lwr[oosplot_species_GLMM$excluded_plot == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosplot_species_GLMM$p[oosplot_species_GLMM$excluded_plot == i] <- estimates[2,4]
  
  print(i)
}

oosplot_species_GLMM_SIG <- oosplot_species_GLMM[oosplot_species_GLMM$p < 0.05,]

# plot
plot2 <- ggplot(oosplot_species_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  geom_hline(yintercept = full_coefs$effect, color = berlin[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = berlin[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = berlin[6], lty=2, lwd=1.5)+
  geom_errorbar(data=oosplot_species_GLMM, aes( x=excluded_plot , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = berlin[4])+
  geom_point(data=oosplot_species_GLMM ,aes(excluded_plot, estimate, color = berlin[4]),pch=20 ,size=5) +
  geom_point(data=oosplot_species_GLMM ,aes(excluded_plot, estimate, color = berlin[2]),pch=20 ,size=4) +
  geom_point(data=oosplot_species_GLMM_SIG ,aes(excluded_plot, estimate, color = berlin[4]),pch=8 ,size=4) +
  scale_color_manual(values = c(berlin[4], berlin[2]), name = "random year effect", labels=c("without", "with" ))+
  theme(legend.position="none")

grid.arrange(plot1, plot2, ncol=2)



# evenness 
coefs <- coef(summary(evenness_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[2,1],
                            lwr   = estimates[2,1] + 1.96*estimates[2,2],
                            upr   = estimates[2,1] - 1.96*estimates[2,2],
                            p     = estimates[2,4])

#OOS years
years <- c(sort(unique(LTdiversity$year)))

oosyear_evenness_GLMM <- data.frame(excluded_year= c( 2001:2022))

for(i in years){
  loop_data <- LTdiversity[LTdiversity$year != i,]
  
  loop_mod<- glmmTMB(evenness ~ scaled_year + (1|year2) + (1|site/plot) + ar1(year2 +0|plot ), family = beta_family(link="logit"), data=loop_data)
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosyear_evenness_GLMM$estimate[oosyear_evenness_GLMM$excluded_year == i] <- estimates[2,1]
  
  oosyear_evenness_GLMM$estimate_upr[oosyear_evenness_GLMM$excluded_year == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosyear_evenness_GLMM$estimate_lwr[oosyear_evenness_GLMM$excluded_year == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosyear_evenness_GLMM$p[oosyear_evenness_GLMM$excluded_year == i] <- estimates[2,4]
  
  print(i)
}

oosyear_evenness_GLMM_SIG <- oosyear_evenness_GLMM[oosyear_evenness_GLMM$p < 0.05,]

plot1 <- ggplot(oosyear_evenness_GLMM, aes(excluded_year, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "Trend coefficient", x = "Excluded year")+
  ggtitle("Evenness trend sensitivity")+
  theme(plot.title = element_text(hjust = 0.05,vjust=0))+  
  scale_x_continuous(labels= c(  "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(2001:2022) )+
  geom_hline(yintercept = full_coefs$effect, color = berlin[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = berlin[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = berlin[6], lty=2, lwd=1.5)+
  geom_errorbar(data=oosyear_evenness_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = berlin[4])+
  geom_point(data=oosyear_evenness_GLMM ,aes(excluded_year, estimate, color = berlin[4]),pch=20 ,size=5) +
  geom_point(data=oosyear_evenness_GLMM ,aes(excluded_year, estimate, color = berlin[2]),pch=20 ,size=4) +
  geom_point(data=oosyear_evenness_GLMM_SIG ,aes(excluded_year, estimate, color = berlin[4]),pch=8 ,size=4) +
  scale_color_manual(values = c(berlin[4], berlin[2]), name = "random year effect", labels=c("without", "with" ))+
  theme(legend.position="none")


# OOS plot
plots <- unique(LTdiversity$plot)

oosplot_evenness_GLMM <- data.frame(excluded_plot= as.factor(unique(LTdiversity$plot)))

for(i in plots){
  loop_data <- LTdiversity[LTdiversity$plot != i,]
  
  loop_mod<- glmmTMB(evenness ~ scaled_year + (1|year2) + (1|site/plot) + ar1(year2 +0|plot ), family = beta_family(link="logit"), data=loop_data)
  
  coefs<-coef(summary(loop_mod))
  
  estimates <- coefs$cond
  
  oosplot_evenness_GLMM$estimate[oosplot_evenness_GLMM$excluded_plot == i] <- estimates[2,1]
  
  oosplot_evenness_GLMM$estimate_upr[oosplot_evenness_GLMM$excluded_plot == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosplot_evenness_GLMM$estimate_lwr[oosplot_evenness_GLMM$excluded_plot == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosplot_evenness_GLMM$p[oosplot_evenness_GLMM$excluded_plot == i] <- estimates[2,4]
  
  print(i)
}

oosplot_evenness_GLMM_SIG <- oosplot_evenness_GLMM[oosplot_evenness_GLMM$p < 0.05,]

plot2 <- ggplot(oosplot_evenness_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  geom_hline(yintercept = full_coefs$effect, color = berlin[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = berlin[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = berlin[6], lty=2, lwd=1.5)+
  geom_errorbar(data=oosplot_evenness_GLMM, aes( x=excluded_plot , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = berlin[4])+
  geom_point(data=oosplot_evenness_GLMM ,aes(excluded_plot, estimate, color = berlin[4]),pch=20 ,size=5) +
  geom_point(data=oosplot_evenness_GLMM ,aes(excluded_plot, estimate, color = berlin[2]),pch=20 ,size=4) +
  geom_point(data=oosplot_evenness_GLMM_SIG ,aes(excluded_plot, estimate, color = berlin[4]),pch=8 ,size=4) +
  scale_color_manual(values = c(berlin[4], berlin[2]), name = "random year effect", labels=c("without", "with" ))+
  theme(legend.position="none")


grid.arrange(plot1, plot2, ncol=2)


# turnover 
coefs <- coef(summary(jaccard_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[2,1],
                            lwr   = estimates[2,1] + 1.96*estimates[2,2],
                            upr   = estimates[2,1] - 1.96*estimates[2,2],
                            p     = estimates[2,4])
#OOS years
years <- c(sort(unique(turnover_data$year)))

oosyear_turnover_GLMM <- data.frame(excluded_year= c( 2003:2022))

for(i in years){
  loop_data <- turnover_data[turnover_data$year != i,]
  
  loop_mod<- glmmTMB(jaccard ~ scaled_year + (1|year2) + (1|site/plot), family = beta_family(link="logit"), data=loop_data)
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosyear_turnover_GLMM$estimate[oosyear_turnover_GLMM$excluded_year == i] <- estimates[2,1]
  
  oosyear_turnover_GLMM$estimate_upr[oosyear_turnover_GLMM$excluded_year == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosyear_turnover_GLMM$estimate_lwr[oosyear_turnover_GLMM$excluded_year == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosyear_turnover_GLMM$p[oosyear_turnover_GLMM$excluded_year == i] <- estimates[2,4]
  
  print(i)
}

oosyear_turnover_GLMM_SIG <- oosyear_turnover_GLMM[oosyear_turnover_GLMM$p < 0.05,]

# plot
plot1 <- ggplot(oosyear_turnover_GLMM, aes(excluded_year, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "Trend coefficient", x = "Excluded year")+
  ggtitle("Turnover trend sensitivity")+
  theme(plot.title = element_text(hjust = 0.05,vjust=0))+  
  scale_x_continuous(labels= c(   "2003", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(2003:2022) )+
  geom_hline(yintercept = full_coefs$effect, color = berlin[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = berlin[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = berlin[6], lty=2, lwd=1.5)+
  geom_errorbar(data=oosyear_turnover_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = berlin[4])+
  geom_point(data=oosyear_turnover_GLMM ,aes(excluded_year, estimate, color = berlin[4]),pch=20 ,size=5) +
  geom_point(data=oosyear_turnover_GLMM ,aes(excluded_year, estimate, color = berlin[2]),pch=20 ,size=4) +
  geom_point(data=oosyear_turnover_GLMM_SIG ,aes(excluded_year, estimate, color = berlin[4]),pch=8 ,size=4) +
  scale_color_manual(values = c(berlin[4], berlin[2]), name = "random year effect", labels=c("without", "with" ))+
  theme(legend.position="none")


# OOS plot
plots <- unique(turnover_data$plot)

oosplot_turnover_GLMM <- data.frame(excluded_plot= as.factor(unique(turnover_data$plot)))

for(i in plots){
  loop_data <- turnover_data[turnover_data$plot != i,]
  
  loop_mod<- glmmTMB(jaccard ~ scaled_year + (1|year2) + (1|site/plot), family = beta_family(link="logit"), data=loop_data)
  
  coefs<-coef(summary(loop_mod))
  
  estimates <- coefs$cond
  
  oosplot_turnover_GLMM$estimate[oosplot_turnover_GLMM$excluded_plot == i] <- estimates[2,1]
  
  oosplot_turnover_GLMM$estimate_upr[oosplot_turnover_GLMM$excluded_plot == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosplot_turnover_GLMM$estimate_lwr[oosplot_turnover_GLMM$excluded_plot == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosplot_turnover_GLMM$p[oosplot_turnover_GLMM$excluded_plot == i] <- estimates[2,4]
  
  print(i)
}

oosplot_turnover_GLMM_SIG <- oosplot_turnover_GLMM[oosplot_turnover_GLMM$p < 0.05,]


plot2 <- ggplot(oosplot_turnover_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  geom_hline(yintercept = full_coefs$effect, color = berlin[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = berlin[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = berlin[6], lty=2, lwd=1.5)+
  geom_errorbar(data=oosplot_turnover_GLMM, aes( x=excluded_plot , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = berlin[4])+
  geom_point(data=oosplot_turnover_GLMM ,aes(excluded_plot, estimate, color = berlin[4]),pch=20 ,size=5) +
  geom_point(data=oosplot_turnover_GLMM ,aes(excluded_plot, estimate, color = berlin[2]),pch=20 ,size=4) +
  geom_point(data=oosplot_turnover_GLMM_SIG ,aes(excluded_plot, estimate, color = berlin[4]),pch=8 ,size=4) +
  scale_color_manual(values = c(berlin[4], berlin[2]), name = "random year effect", labels=c("without", "with" ))+
  theme(legend.position="none")

grid.arrange(plot1, plot2, ncol=2)
















