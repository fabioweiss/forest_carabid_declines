# SUPPLEMENTARY R-SCRIPT TO:
# Long-term data reveal: Recent declines in carabid beetles in a temperate forest are linked to severe droughts
# F.Weiss, H.von Wehrden & A.Linde
# F.Weiss: ORCID 0000-0003-1078-1528

# PART II: Modelling abundance and biomass

# Data available at: 

#### Package list #####

library(dplyr)
library(bbmle)
library(sjPlot)
library(gridExtra)
library(DHARMa)
library(gamm4)
library(glmmTMB)
library(cols4all)
library(ggplot2)
library(ggeffects)

colors <- c4a("berlin", 7)

#### Loading data ####

# loading raw data
beetle_samples <- read.csv("EWcarabids1999-2022_rawdata.csv")

# loading sampling meta data
sampling_meta <- read.csv2("sampling_meta_dwd_spei.csv")

# laoding carabid weights
# carabid weights were calculated using the approach by Weiss & Linde (2022, https://doi.org/10.1007/s10841-022-00391-6)
carabid_weights <- read.csv("EWcarabids1999-2022_species_biomass.csv")


#### Pitfall trap bias correction ####


# assign weights to data based on species
# quickly correct one name
beetle_samples$species[beetle_samples$species== "Lonicera pilicornis"] <- "Loricera pilicornis"

beetle_samples$species_weight <- carabid_weights$weight[match(beetle_samples$species, carabid_weights$species)]

beetle_samples$species_weight[beetle_samples$species =="no carabids"] <- 0


## pitfall trap bias abundance correction sensu Engel et al. (2017, https://doi.org/10.1002/ecs2.1790)

# preparing a second species weight which is levelled off at 300mg for calculating Engel`s delta
beetle_samples$species_weight2 <- beetle_samples$species_weight
beetle_samples$species_weight2[beetle_samples$species_weight2 > 300] <- 300 
beetle_samples$species_weight2[beetle_samples$species_weight2 > 0 & beetle_samples$species_weight <1] <- 1

# Adding mean temperature from trappind DWD data:


sampling_meta$sample_id <- as.factor(sampling_meta$sample_id)
beetle_samples$sample_id <- as.factor(beetle_samples$sample_id)

beetle_samples$mean_temp <- sampling_meta$sampling_temp[match(beetle_samples$sample_id, sampling_meta$sample_id)]
beetle_samples$mean_temp <- round(beetle_samples$mean_temp, digits = 0)

# modelling bias slopes according to Engel et al. (2017, Appendix S9)

# correction factors for 4trap/grid 
corr_fact_Engel <- data.frame(temp = c(15,18,21,24,27,30), slope =c( -0.51, -0.46, -0.43, -0.41, -0.38, -0.32))
mcorr <- lm(slope ~ temp + I(temp^2), data=corr_fact_Engel)
newdata <- data.frame(temp = c(0:40))
corr_fact_grid <- data.frame(temp=newdata$temp, slope= predict(mcorr, newdata))
plot(corr_fact_grid$temp,corr_fact_grid$slope)

# correction factors for 4trap/transect 
# below 24degree almost constant at -0.495
corr_fact_Engel <- data.frame(temp = c(24,27,30), slope =c( -0.49, -0.46, -0.4))
mcorr <- lm(slope ~ temp, data=corr_fact_Engel)
newdata <- data.frame(temp = c(0:40))
corr_fact_transect <- data.frame(temp=newdata$temp, slope= predict(mcorr, newdata))
corr_fact_transect$slope[corr_fact_transect$temp < 24] <- -0.495
plot(corr_fact_transect$temp,corr_fact_transect$slope)



# correction factor slope (ß) 4trap/grid
beetle_samples$beta <-  corr_fact_grid$slope[match(beetle_samples$mean_temp, corr_fact_grid$temp)]

# adapting correction factor slope (ß) 4trap/transect 
beetle_samples$beta[beetle_samples$plot =="01"] <- corr_fact_transect$slope[match(beetle_samples$mean_temp[beetle_samples$plot =="01"], corr_fact_transect$temp)]
beetle_samples$beta[beetle_samples$plot =="2s"] <- corr_fact_transect$slope[match(beetle_samples$mean_temp[beetle_samples$plot =="2s"], corr_fact_transect$temp)]
beetle_samples$beta[beetle_samples$plot =="07"] <- corr_fact_transect$slope[match(beetle_samples$mean_temp[beetle_samples$plot =="07"], corr_fact_transect$temp)]


# calculate correction factor sigma
beetle_samples$sigma <- beetle_samples$species_weight2^(beetle_samples$beta)

# deal with Zeros
beetle_samples$sigma[beetle_samples$abundance == 0] <- 1


## assign sampling abundances, corrected abundance, sampling biomass and corrected biomass

# rename abundance
beetle_samples$sampling_abundance <- beetle_samples$abundance

# corrected abundance
beetle_samples$corr_abundance <- beetle_samples$sampling_abundance * beetle_samples$sigma

# adjust corrected abundance for count data
beetle_samples$corr_abundance2 <- as.integer(ceiling(beetle_samples$corr_abundance * 10))

# check, all have assigned weights
beetle_samples[is.na(beetle_samples$species_weight),]

# sampling biomass
beetle_samples$weight_sum <- beetle_samples$species_weight * beetle_samples$sampling_abundance

# corrected biomass
beetle_samples$corr_weight_sum <- beetle_samples$species_weight * beetle_samples$corr_abundance


# Save processed data for later use in diversity analysis
write.csv(beetle_samples,"carabid_samples_raw2.csv", row.names=FALSE )


#### Aggregating abundance and biomass ####

# aggregate abundance by sampleID*trap
agg_samples <- aggregate(cbind(beetle_samples$sampling_abundance, beetle_samples$corr_abundance2, beetle_samples$weight_sum, beetle_samples$corr_weight_sum), by= list(beetle_samples$sample_id, beetle_samples$sample_id2, beetle_samples$year,  beetle_samples$interval, beetle_samples$trap, beetle_samples$plot, beetle_samples$site), FUN=sum  )

agg_samples<-agg_samples %>% 
  rename(
    sample_id = Group.1,
    sample_id2 = Group.2,
    year = Group.3,
    interval = Group.4,
    trap = Group.5,
    plot= Group.6,
    site=Group.7,
    sampling_abundance = V1,
    corr_abundance = V2,
    sampling_biomass = V3,
    corr_biomass = V4
  )


#### Add information about N.brevicollis####

# aggregate data of N.brevicollis
brevicollis <- beetle_samples[beetle_samples$species == "Nebria brevicollis",]

brevicollis_agg <- aggregate(cbind(brevicollis$sampling_abundance, brevicollis$corr_abundance2, brevicollis$weight_sum, brevicollis$corr_weight_sum), by= list(brevicollis$sample_id, brevicollis$sample_id2, brevicollis$year, brevicollis$interval,brevicollis$trap, brevicollis$plot, brevicollis$site), FUN=sum  )

brevicollis_agg<-brevicollis_agg %>% 
  rename(
    sample_id = Group.1,
    sample_id2 = Group.2,
    year = Group.3,
    interval = Group.4,
    trap = Group.5,
    plot= Group.6,
    site=Group.7,
    sampling_abundance = V1,
    corr_abundance = V2,
    sampling_biomass = V3,
    corr_biomass = V4
  )

# match n brevicollis abundance with full data

agg_samples$sampling_abundance_nb <- brevicollis_agg$sampling_abundance[match(agg_samples$sample_id2, brevicollis_agg$sample_id2)]
agg_samples$sampling_abundance_nb[is.na(agg_samples$sampling_abundance_nb)] <- 0

agg_samples$corr_abundance_nb <- brevicollis_agg$corr_abundance[match(agg_samples$sample_id2, brevicollis_agg$sample_id2)]
agg_samples$corr_abundance_nb[is.na(agg_samples$corr_abundance_nb)] <- 0

agg_samples$sampling_biomass_nb <- brevicollis_agg$sampling_biomass[match(agg_samples$sample_id2, brevicollis_agg$sample_id2)]
agg_samples$sampling_biomass_nb[is.na(agg_samples$sampling_biomass_nb)] <- 0

agg_samples$corr_biomass_nb <- brevicollis_agg$corr_biomass[match(agg_samples$sample_id2, brevicollis_agg$sample_id2)]
agg_samples$corr_biomass_nb[is.na(agg_samples$corr_biomass_nb)] <- 0


#### Match sampling metadata with aggregated samples ####

sampling_meta$sample_id <- as.factor(sampling_meta$sample_id)
agg_samples$sample_id <- as.factor(agg_samples$sample_id)

sampling_meta <- sampling_meta %>% 
  rename("sampling_length"="Trapping_length" )


# add metadata by merging
agg_samples2<- merge(agg_samples, sampling_meta[,c(7,9:30)], by.x = "sample_id", by.y ="sample_id") 

# remove all samples of intervals with more or less than 4 traps
agg_samples2 <- agg_samples2[agg_samples2$sampling_effort == 4,]


#### Restructuring data ####

#rename data
LTabundance <-agg_samples2

# year variable
LTabundance$year <- as.integer(LTabundance$year)
LTabundance$year2 <- as.factor(LTabundance$year)
LTabundance$scaled_year <- c(scale(LTabundance$year))

# sampling length and interval
LTabundance$scaled_length <- c(scale(LTabundance$sampling_length))
LTabundance$interval <- as.factor(LTabundance$interval)

# abundances and biomass excluding Nebria brevicollis
LTabundance$sampling_abundance2 <- LTabundance$sampling_abundance - LTabundance$sampling_abundance_nb
LTabundance$corr_abundance2 <- LTabundance$corr_abundance - LTabundance$corr_abundance_nb
LTabundance$sampling_biomass2 <- LTabundance$sampling_biomass - LTabundance$sampling_biomass_nb
LTabundance$corr_biomass2 <- LTabundance$corr_biomass - LTabundance$corr_biomass_nb

# new trap ID for random intercept
LTabundance$trapID <- paste(LTabundance$plot,LTabundance$year,LTabundance$trap, sep="_")
LTabundance$trapID <- as.factor(LTabundance$trapID)
LTabundance$plot <- as.factor(LTabundance$plot)
LTabundance$site <-as.factor( LTabundance$site )

# scale weather variables
LTabundance$scaled_sampling_temp <- c(scale(LTabundance$sampling_temp))
LTabundance$scaled_sampling_rain <- c(scale(LTabundance$sampling_rain))

# SPEI data is already scaled

# plot data availability
ggplot(LTabundance, aes(year, plot)) +
  geom_point(shape = 15, color = colors[2], size = 3)+
  ylab("Plot")+
  xlab("Year")+
  ggtitle("Data availability for abundance and biomass")+
  theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20), 
        axis.text.x  = element_text(size=10), axis.text.y  = element_text(size=20) )


#### Modelling abundance ####

# GLMM
corr_ab_glmm <- glmmTMB(corr_abundance2 ~ scaled_year + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp +   scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= LTabundance, ziformula = ~ 1, family=nbinom1)

# Model including N.brevicollis:
# corr_ab_glmm <- glmmTMB(corr_abundance ~ scaled_year + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp +   scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= LTabundance, ziformula = ~ 1, family=nbinom2)

# DHARMa diagnostics
dharma_sim1 <- simulateResiduals(fittedModel = corr_ab_glmm, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)
testZeroInflation(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))

tab_model(corr_ab_glmm, dv.labels = c("Corrected abundance GLMM"))


# GAMM 
# fit with unconstrained k and not random intercept for year first to fix k to appropriate dimension in a second run (Knape 2016)
# theta of the respective glmm with Ziformula=~0, family=negbinom2

corr_ab_gamm <- gamm4(corr_abundance2 ~ s(scaled_year)+ scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp + scaled_sampling_rain ,
                      random= ~ (1|site/plot/trapID), 
                      data= LTabundance, 
                      family=negbin(theta= 4.19),
                      REML=TRUE)

plot(corr_ab_gamm$gam)

corr_ab_gamm <- gamm4(corr_abundance2 ~ s(scaled_year, k=6, fx=TRUE) + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp + scaled_sampling_rain ,
                      random= ~  (1|year2) + (1|site/plot/trapID), 
                      data= LTabundance, 
                      family=negbin(theta= 4.19),
                      REML=TRUE)

# Model including N.brevicollis:
# corr_ab_gamm <- gamm4(corr_abundance ~ s(scaled_year, k=7, fx=TRUE) + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp + scaled_sampling_rain ,random= ~  (1|year2) + (1|site/plot/trapID), data= LTabundance, family=negbin(theta= 4.19), REML=TRUE)


dharma_sim2 <- simulateResiduals(fittedModel = corr_ab_gamm$mer, re.form= NULL)
plot(dharma_sim2)
testDispersion(dharma_sim2)
testZeroInflation(dharma_sim2)

dharma_sim2.1 = recalculateResiduals(dharma_sim2, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim2.1, time = unique( LTabundance$scaled_year))

plot(corr_ab_gamm$gam)
tab_model(corr_ab_gamm$mer, dv.labels = c("Corrected abundance GAMM"))


# SPEI GLMM

# SPEI modelselection
# SPEI variables are already scaled

m1 <- update(corr_ab_glmm,. ~ . + spei12)
m2 <- update(corr_ab_glmm,. ~ . + spei12_1)
m3 <- update(corr_ab_glmm,. ~ . + spei12_2)
m4 <- update(corr_ab_glmm,. ~ . + spei24)
m5 <- update(corr_ab_glmm,. ~ . + spei24_1)
m6 <- update(corr_ab_glmm,. ~ . + spei24_2)
m7 <- update(corr_ab_glmm,. ~ . + spei36)
m8 <- update(corr_ab_glmm,. ~ . + spei36_1)
m9 <- update(corr_ab_glmm,. ~ . + spei36_2)
m10 <- update(corr_ab_glmm,. ~ . + spei48)
m11 <- update(corr_ab_glmm,. ~ . + spei48_1)
m12 <- update(corr_ab_glmm,. ~ . + spei48_2)
m13 <- update(corr_ab_glmm,. ~ . + spei60)
m14 <- update(corr_ab_glmm,. ~ . + spei60_1)
m15 <- update(corr_ab_glmm,. ~ . + spei60_2)
m16 <- update(corr_ab_glmm,. ~ . + spei72)
m17 <- update(corr_ab_glmm,. ~ . + spei72_1)
m18 <- update(corr_ab_glmm,. ~ . + spei72_2)

AICtab(corr_ab_glmm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18)

ab_spei_glmm <- m18

dharma_sim1 <- simulateResiduals(fittedModel = ab_spei_glmm, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)
testZeroInflation(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))

tab_model(ab_spei_glmm, dv.labels = c("Corrected abundance GLMM (SPEI)"))

#### Abundance figure ####


# GLMM 
ggpreds1 <- ggpredict(corr_ab_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_length = -0.1443603, interval= "2", scaled_sampling_rain=-0.2290402,scaled_sampling_temp= 0.1761928 ))
ggpreds1$year <- unique(LTabundance$year2)  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  ggpreds1$predicted,
                           lwr   = ggpreds1$conf.low,
                           upr   = ggpreds1$conf.high   )


# mean percentage change
tsdata <- ts(data.frame(x1=c(preds_glmm$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)



plot1 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "Abundance (bias corrected)", x = "Year")+
  ggtitle("A")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 10, 20 , 30, 40))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=corr_abundance2, scaled_year), alpha=0.2, color= "black")+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=3)+
  
  theme(legend.position="none")


# GAMM 
ggpreds3 <- ggpredict(corr_ab_gamm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_length = -0.1443603, interval= "2", scaled_sampling_rain=-0.2290402,scaled_sampling_temp= 0.1761928 ))
ggpreds3$year <- unique(LTabundance$year2)  

preds_gamm <- data.frame(  scaled_year= ggpreds3$x,
                           year = ggpreds3$year,
                           fit  =  ggpreds3$predicted,
                           lwr   = ggpreds3$conf.low,
                           upr   = ggpreds3$conf.high   )


# mean percentage change
tsdata <- ts(data.frame(x1=c(preds_gamm$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)

plot2 <- ggplot(preds_gamm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "Year")+
  ggtitle("B")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 10, 20 , 30, 40))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=corr_abundance2, scaled_year), alpha=0.2, color= "black")+
  
  geom_ribbon(data = preds_gamm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[2])+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color =colors[2], lty=1, lwd=3)+
  
  theme(legend.position="none")

# SPEI glmm 
ggpreds1 <- ggpredict(ab_spei_glmm, terms= c("scaled_year [all]", "spei72_2 [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_length = -0.1443603, interval= "2", scaled_sampling_rain=-0.2290402,scaled_sampling_temp= 0.1761928))

ggpreds1$year <- LTabundance$year2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1$spei <- LTabundance$spei72_2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1 <- ggpreds1[ggpreds1$group == round(ggpreds1$spei, digits=3),]  

preds_glmm1 <- data.frame(  scaled_year= ggpreds1$x,
                            year = ggpreds1$year,
                            fit  =  ggpreds1$predicted,
                            lwr   = ggpreds1$conf.low,
                            upr   = ggpreds1$conf.high   )


plot3 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "Year")+
  ggtitle("C")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(0, 50), breaks = c(0, 10, 20 , 30, 40))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=corr_abundance2, scaled_year), alpha=0.2, color= "black")+
  
  geom_ribbon(data = preds_glmm1, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[6])+
  
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=3)+
  
  theme(legend.position="none")


plot4 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "Year")+
  ggtitle("D")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(2, 11), breaks = c(3,6,9))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=3)+
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=3)+
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=3)+
  
  theme(legend.position="top")

grid.arrange(plot1, plot2, plot3,  plot4,  ncol=4)


#### Abundance sensitivity ####

# GLMM 
coefs <- coef(summary(corr_ab_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[2,1],
                            lwr   = estimates[2,1] + 1.96*estimates[2,2],
                            upr   = estimates[2,1] - 1.96*estimates[2,2],
                            p     = estimates[2,4])

#OOS years

years <- sort(c(unique(LTabundance$year)))

oos_abundance_GLMM <- data.frame(excluded_year= c( 1999:2022))

for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <-  glmmTMB(corr_abundance2 ~ scaled_year + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp +   scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= loop_data, ziformula = ~1, family=nbinom1)
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oos_abundance_GLMM$estimate[oos_abundance_GLMM$excluded_year == i] <- estimates[2,1]
  
  oos_abundance_GLMM$estimate_upr[oos_abundance_GLMM$excluded_year == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oos_abundance_GLMM$estimate_lwr[oos_abundance_GLMM$excluded_year == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oos_abundance_GLMM$p[oos_abundance_GLMM$excluded_year == i] <- estimates[2,4]
  
  print(i)
}

oos_abundance_GLMM_SIG <- oos_abundance_GLMM[oos_abundance_GLMM$p < 0.05,]


plot1 <- ggplot(oos_abundance_GLMM, aes(excluded_year, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "Trend coefficient (scaled)", x = "Excluded year")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  
  geom_hline(yintercept = full_coefs$effect, color = colors[6] , lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6] , lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6] , lty=2, lwd=1.5)+
  
  geom_errorbar(data=oos_abundance_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oos_abundance_GLMM ,aes(excluded_year, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oos_abundance_GLMM ,aes(excluded_year, estimate, color = colors[2] ),pch=20 ,size=4) +
  geom_point(data=oos_abundance_GLMM_SIG ,aes(excluded_year, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2] , "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")


# OOS plot

plots <- unique(LTabundance$plot)

oosplot_abundance_GLMM <- data.frame(excluded_plot= as.factor(unique(LTabundance$plot)))

for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <-  glmmTMB(corr_abundance2 ~ scaled_year + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp +   scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= loop_data, ziformula = ~1, family=nbinom1)
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosplot_abundance_GLMM$estimate[oosplot_abundance_GLMM$excluded_plot == i] <- estimates[2,1]
  
  oosplot_abundance_GLMM$estimate_upr[oosplot_abundance_GLMM$excluded_plot == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosplot_abundance_GLMM$estimate_lwr[oosplot_abundance_GLMM$excluded_plot == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosplot_abundance_GLMM$p[oosplot_abundance_GLMM$excluded_plot == i] <- estimates[2,4]
  
  print(i)
}

oosplot_abundance_GLMM_SIG <- oosplot_abundance_GLMM[oosplot_abundance_GLMM$p < 0.05,]


plot2 <- ggplot(oosplot_abundance_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  

  geom_hline(yintercept = full_coefs$effect, color = colors[6] , lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6] , lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6] , lty=2, lwd=1.5)+
  
  geom_errorbar(data=oosplot_abundance_GLMM, aes( x=excluded_plot , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosplot_abundance_GLMM ,aes(excluded_plot, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosplot_abundance_GLMM ,aes(excluded_plot, estimate, color = colors[2] ),pch=20 ,size=4) +
  geom_point(data=oosplot_abundance_GLMM_SIG ,aes(excluded_plot, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2] , "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")

grid.arrange(plot1, plot2, ncol=2)


# GAMM 

years <- sort(c(unique(LTabundance$year)))

par(mfrow=c(6,4))

# OOS years
for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <-  gamm4(corr_abundance2 ~ s(scaled_year, k=6, fx=TRUE) + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp + scaled_sampling_rain ,
                     random= ~  (1|year2) + (1|site/plot/trapID), 
                     data= loop_data, 
                     family=negbin(theta= 4.16),
                     REML=TRUE)
  
  plot(loop_mod$gam, main = as.character(i), ylab="", xlab="", col=colors[2], shade=T, shade.col = colors[7], cex=2, ylim=c(-1.6, 1.1))
  
  print(i)
}


# oos plot
plots <- c(unique(LTabundance$plot))

par(mfrow=c(4,4))

for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <-  gamm4(corr_abundance2 ~ s(scaled_year, k=6, fx=TRUE) + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp + scaled_sampling_rain ,
                     random= ~  (1|year2) + (1|site/plot/trapID), 
                     data= loop_data, 
                     family=negbin(theta= 4.16),
                     REML=TRUE)
  
  plot(loop_mod$gam, main = as.character(i), ylab="", xlab="", col=colors[2], shade=T, shade.col = colors[7], cex=2, ylim=c(-1.6, 1.1))
  
  print(i)
}

# SPEI GLMM 

coefs <- coef(summary(ab_spei_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[9,1],
                            lwr   = estimates[9,1] + 1.96*estimates[9,2],
                            upr   = estimates[9,1] - 1.96*estimates[9,2],
                            p     = estimates[9,4])

#OOS years

years <- sort(c(unique(LTabundance$year)))

oosyears_spei_GLMM <- data.frame(excluded_year= c( 1999:2022))

for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <-  glmmTMB(corr_abundance2 ~ scaled_year + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp +   scaled_sampling_rain  + spei72_2 + (1|year2) + (1|site/plot/trapID), data= loop_data, ziformula = ~1, family=nbinom1)
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosyears_spei_GLMM$estimate[oosyears_spei_GLMM$excluded_year == i] <- estimates[9,1]
  
  oosyears_spei_GLMM$estimate_upr[oosyears_spei_GLMM$excluded_year == i] <- estimates[9,1] + 1.96*estimates[9,2]
  
  oosyears_spei_GLMM$estimate_lwr[oosyears_spei_GLMM$excluded_year == i] <- estimates[9,1] - 1.96*estimates[9,2]
  
  oosyears_spei_GLMM$p[oosyears_spei_GLMM$excluded_year == i] <- estimates[9,4]
  
  print(i)
}

oosyears_spei_GLMM_SIG <- oosyears_spei_GLMM[oosyears_spei_GLMM$p < 0.05,]


plot1 <- ggplot(oosyears_spei_GLMM, aes(excluded_year, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI coefficient (scaled)", x = "Excluded year")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  

  scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  
  geom_hline(yintercept = full_coefs$effect, color = colors[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6], lty=2, lwd=1.5)+
  
  geom_errorbar(data=oosyears_spei_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosyears_spei_GLMM ,aes(excluded_year, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosyears_spei_GLMM ,aes(excluded_year, estimate, color = colors[2]),pch=20 ,size=4) +
  geom_point(data=oosyears_spei_GLMM_SIG ,aes(excluded_year, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2], "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")


# OOS plot

plots <- unique(LTabundance$plot)

oosplot_spei_GLMM <- data.frame(excluded_plot= as.factor(unique(LTabundance$plot)))

for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <-  glmmTMB(corr_abundance2 ~ scaled_year + scaled_length+ I(scaled_length^2) + interval * scaled_sampling_temp +   scaled_sampling_rain  + spei72_2 + (1|year2) + (1|site/plot/trapID), data= loop_data, ziformula = ~1, family=nbinom1)
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosplot_spei_GLMM$estimate[oosplot_spei_GLMM$excluded_plot == i] <- estimates[9,1]
  
  oosplot_spei_GLMM$estimate_upr[oosplot_spei_GLMM$excluded_plot == i] <- estimates[9,1] + 1.96*estimates[9,2]
  
  oosplot_spei_GLMM$estimate_lwr[oosplot_spei_GLMM$excluded_plot == i] <- estimates[9,1] - 1.96*estimates[9,2]
  
  oosplot_spei_GLMM$p[oosplot_spei_GLMM$excluded_plot == i] <- estimates[9,4]
  
  print(i)
}

oosplot_spei_GLMM_SIG <- oosplot_spei_GLMM[oosplot_spei_GLMM$p < 0.05,]


plot2 <- ggplot(oosplot_spei_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  

  geom_hline(yintercept = full_coefs$effect, color = colors[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6],  lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6],  lty=2, lwd=1.5)+
  
  geom_errorbar(data=oosplot_spei_GLMM, aes( x=excluded_plot , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosplot_spei_GLMM ,aes(excluded_plot, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosplot_spei_GLMM ,aes(excluded_plot, estimate, color = colors[2]),pch=20 ,size=4) +
  geom_point(data=oosplot_spei_GLMM_SIG ,aes(excluded_plot, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2], "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")

grid.arrange(plot1, plot2, ncol=2)




#### Modelling biomass ####

# Cubic root transformation for biomass
hist(LTabundance$corr_biomass2^(1/3))
LTabundance$crt_corr_biomass2 <- (LTabundance$corr_biomass2)^(1/3)

# GLMM
biomass_glmm <- glmmTMB(crt_corr_biomass2 ~ scaled_year + scaled_length+  interval * scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= LTabundance, family=gaussian(link=identity))

dharma_sim1 <- simulateResiduals(fittedModel = biomass_glmm, re.form= NULL)
plot(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))

tab_model(biomass_glmm, dv.labels = c("Biomass GLMM"))

# GAMM 
biomass_gamm <- gamm4(crt_corr_biomass2 ~ s(scaled_year)+ scaled_length+ interval * scaled_sampling_temp + scaled_sampling_rain ,
                      random= ~ (1|site/plot/trapID), 
                      data= LTabundance, 
                      family=gaussian(link=identity),
                      REML=TRUE)

plot(biomass_gamm$gam)

biomass_gamm <- gamm4(crt_corr_biomass2 ~ s(scaled_year, k=8, fx=TRUE )+ scaled_length+ interval * scaled_sampling_temp + scaled_sampling_rain ,
                      random= ~(1|year2) + (1|site/plot/trapID), 
                      data= LTabundance, 
                      family=gaussian(link=identity),
                      REML=TRUE)


dharma_sim2 <- simulateResiduals(fittedModel = biomass_gamm$mer, re.form= NULL)
plot(dharma_sim2)

dharma_sim2.1 = recalculateResiduals(dharma_sim2, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim2.1, time = unique( LTabundance$scaled_year))

plot(biomass_gamm$gam)
tab_model(biomass_gamm$mer, dv.labels = c("Biomass GAMM"))

# SPEI GLMM 

# modelselection
m1 <- update(biomass_glmm,. ~ . + spei12)
m2 <- update(biomass_glmm,. ~ . + spei12_1)
m3 <- update(biomass_glmm,. ~ . + spei12_2)
m4 <- update(biomass_glmm,. ~ . + spei24)
m5 <- update(biomass_glmm,. ~ . + spei24_1)
m6 <- update(biomass_glmm,. ~ . + spei24_2)
m7 <- update(biomass_glmm,. ~ . + spei36)
m8 <- update(biomass_glmm,. ~ . + spei36_1)
m9 <- update(biomass_glmm,. ~ . + spei36_2)
m10 <- update(biomass_glmm,. ~ . + spei48)
m11 <- update(biomass_glmm,. ~ . + spei48_1)
m12 <- update(biomass_glmm,. ~ . + spei48_2)
m13 <- update(biomass_glmm,. ~ . + spei60)
m14 <- update(biomass_glmm,. ~ . + spei60_1)
m15 <- update(biomass_glmm,. ~ . + spei60_2)
m16 <- update(biomass_glmm,. ~ . + spei72)
m17 <- update(biomass_glmm,. ~ . + spei72_1)
m18 <- update(biomass_glmm,. ~ . + spei72_2)

AICtab(biomass_glmm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18)

biomass_spei_glmm <- m18

dharma_sim1 <- simulateResiduals(fittedModel = biomass_spei_glmm, re.form= NULL)
plot(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))

tab_model(biomass_spei_glmm, dv.labels = c("Biomass GLMM (SPEI)"))

#### Biomass figure ####


# GLMM 
ggpreds1 <- ggpredict(biomass_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_length = -0.1443603, interval= "2", scaled_sampling_rain=-0.2290402,scaled_sampling_temp= 0.1761928 ))
ggpreds1$year <- unique(LTabundance$year2)  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  (ggpreds1$predicted)^3,
                           lwr   = (ggpreds1$conf.low)^3,
                           upr   = (ggpreds1$conf.high)^3   )

# mean percentage change
tsdata <- ts(data.frame(x1=c(preds_glmm$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)

plot1 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "Biomass (mg)", x = "Year")+
  ggtitle("A")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0, 100, 200 , 300, 400))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=corr_biomass2, scaled_year), alpha=0.2, color= "black")+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=3)+
  
  theme(legend.position="none")


# GAMM 
ggpreds3 <- ggpredict(biomass_gamm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_length = -0.1443603, interval= "2", scaled_sampling_rain=-0.2290402,scaled_sampling_temp= 0.1761928 ))
ggpreds3$year <- sort(unique(LTabundance$year2))  

preds_gamm <- data.frame(  scaled_year= ggpreds3$x,
                           year = ggpreds3$year,
                           fit  =  (ggpreds3$predicted)^3,
                           lwr   = (ggpreds3$conf.low)^3,
                           upr   = (ggpreds3$conf.high)^3   )

# mean percentage change
tsdata <- ts(data.frame(x1=c(preds_gamm$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)

plot2 <- ggplot(preds_gamm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "Year")+
  ggtitle("B")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0, 100, 200 , 300, 400))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=corr_biomass2, scaled_year), alpha=0.2, color= "black")+
  
  geom_ribbon(data = preds_gamm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[2])+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color =colors[2], lty=1, lwd=3)+
  
  theme(legend.position="none")

# SPEI glmm 
ggpreds1 <- ggpredict(biomass_spei_glmm, terms= c("scaled_year [all]", "spei72_2 [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_length = -0.1443603, interval= "2", scaled_sampling_rain=-0.2290402,scaled_sampling_temp= 0.1761928))

ggpreds1$year <- LTabundance$year2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1$spei <- LTabundance$spei72_2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1 <- ggpreds1[ggpreds1$group == round(ggpreds1$spei, digits=3),]  

preds_glmm1 <- data.frame(  scaled_year= ggpreds1$x,
                            year = ggpreds1$year,
                            fit  =  (ggpreds1$predicted)^3,
                            lwr   = (ggpreds1$conf.low)^3,
                            upr   = (ggpreds1$conf.high)^3   )

plot3 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "Year")+
  ggtitle("C")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(0, 500), breaks = c(0, 100, 200 , 300, 400))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=corr_biomass2, scaled_year), alpha=0.2, color= "black")+
  
  geom_ribbon(data = preds_glmm1, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[6])+
  
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=3)+
  
  theme(legend.position="none")


plot4 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "Year")+
  ggtitle("D")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_y_continuous(limits = c(0, 120), breaks = c(30,60,90))+
  
  scale_x_continuous(labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=3)+
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=3)+
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=3)+
  
  theme(legend.position="top")


grid.arrange(plot1, plot2, plot3,  plot4,  ncol=4)

#### Biomass sensitivity ####


#GLMM

coefs <- coef(summary(biomass_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[2,1],
                            lwr   = estimates[2,1] + 1.96*estimates[2,2],
                            upr   = estimates[2,1] - 1.96*estimates[2,2],
                            p     = estimates[2,4])


# OOS years
years <- sort(c(unique(LTabundance$year)))

oosyear_biomass_GLMM <- data.frame(excluded_year= c( 1999:2022))


for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <-  glmmTMB(crt_corr_biomass2 ~ scaled_year+ scaled_length+  interval * scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosyear_biomass_GLMM$estimate[oosyear_biomass_GLMM$excluded_year == i] <- estimates[2,1]
  
  oosyear_biomass_GLMM$estimate_upr[oosyear_biomass_GLMM$excluded_year == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosyear_biomass_GLMM$estimate_lwr[oosyear_biomass_GLMM$excluded_year == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosyear_biomass_GLMM$p[oosyear_biomass_GLMM$excluded_year == i] <- estimates[2,4]
  
  print(i)
}

oosyear_biomass_GLMM_SIG <- oosyear_biomass_GLMM[oosyear_biomass_GLMM$p < 0.05,]


plot1 <- ggplot(oosyear_biomass_GLMM, aes(excluded_year, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "Trend coefficient (scaled)", x = "Excluded year")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  
  geom_hline(yintercept = full_coefs$effect, color = colors[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6], lty=2, lwd=1.5)+
  
  geom_errorbar(data=oosyear_biomass_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosyear_biomass_GLMM ,aes(excluded_year, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosyear_biomass_GLMM ,aes(excluded_year, estimate, color = colors[2]),pch=20 ,size=4) +
  geom_point(data=oosyear_biomass_GLMM_SIG ,aes(excluded_year, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2], "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")

# OOS plots

plots <- c(unique(LTabundance$plot))

oosplot_biomass_GLMM <- data.frame(excluded_plot= as.factor(unique(LTabundance$plot)))


for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <-  glmmTMB(crt_corr_biomass2 ~ scaled_year+ scaled_length+  interval * scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosplot_biomass_GLMM$estimate[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,1]
  
  oosplot_biomass_GLMM$estimate_upr[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosplot_biomass_GLMM$estimate_lwr[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosplot_biomass_GLMM$p[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,4]
  
  print(i)
}

oosplot_biomass_GLMM_SIG <- oosplot_biomass_GLMM[oosplot_biomass_GLMM$p < 0.05,]

plot2 <- ggplot(oosplot_biomass_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  scale_y_continuous(limits = c(-1.1, 0.1))+
  
  geom_hline(yintercept = full_coefs$effect, color = colors[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6], lty=2, lwd=1.5)+  
  
  geom_errorbar(data=oosplot_biomass_GLMM, aes( x=excluded_plot , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosplot_biomass_GLMM ,aes(excluded_plot, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosplot_biomass_GLMM ,aes(excluded_plot, estimate, color = colors[2]),pch=20 ,size=4) +
  geom_point(data=oosplot_biomass_GLMM_SIG ,aes(excluded_plot, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2], "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")

grid.arrange(plot1, plot2, ncol=2)


# GAMM

# oos year
years <- sort(c(unique(LTabundance$year)))

par(mfrow=c(6,4))

for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <-  gamm4(crt_corr_biomass2 ~ s(scaled_year, k=8, fx=TRUE )+ scaled_length+ interval * scaled_sampling_temp + scaled_sampling_rain ,
                     random= ~(1|year2) + (1|site/plot/trapID), 
                     data= loop_data, 
                     family=gaussian(link=identity),
                     REML=TRUE)
  
  
  plot(loop_mod$gam, main = as.character(i), ylab="", xlab="", col=colors[2], shade=T, shade.col = colors[7], cex=2, ylim=c(-3.1, 2))
  
  print(i)
}


# oos plot

plots <- c(unique(LTabundance$plot))

par(mfrow=c(4,4))

for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <-  gamm4(crt_corr_biomass2 ~ s(scaled_year, k=8, fx=TRUE )+ scaled_length+ interval * scaled_sampling_temp + scaled_sampling_rain ,
                     random= ~(1|year2) + (1|site/plot/trapID), 
                     data= loop_data, 
                     family=gaussian(link=identity),
                     REML=TRUE)
  
  
  plot(loop_mod$gam, main = as.character(i), ylab="", xlab="", col=colors[2], shade=T, shade.col = colors[7], cex=2, ylim=c(-3.1, 2))
  
  print(i)
}

par(mfrow=c(1,1))


# SPEI GLMM

coefs <- coef(summary(biomass_spei_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[8,1],
                            lwr   = estimates[8,1] + 1.96*estimates[8,2],
                            upr   = estimates[8,1] - 1.96*estimates[8,2],
                            p     = estimates[8,4])


# OOS years

years <- sort(c(unique(LTabundance$year)))

oosyears_spei_GLMM <- data.frame(excluded_year= c( 1999:2022))

for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <- glmmTMB(crt_corr_biomass2 ~ scaled_year + scaled_length+  interval * scaled_sampling_temp + scaled_sampling_rain + spei72_2  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))  
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosyears_spei_GLMM$estimate[oosyears_spei_GLMM$excluded_year == i] <- estimates[8,1]
  
  oosyears_spei_GLMM$estimate_upr[oosyears_spei_GLMM$excluded_year == i] <- estimates[8,1] + 1.96*estimates[8,2]
  
  oosyears_spei_GLMM$estimate_lwr[oosyears_spei_GLMM$excluded_year == i] <- estimates[8,1] - 1.96*estimates[8,2]
  
  oosyears_spei_GLMM$p[oosyears_spei_GLMM$excluded_year == i] <- estimates[8,4]
  
  print(i)
}

oosyears_spei_GLMM_SIG <- oosyears_spei_GLMM[oosyears_spei_GLMM$p < 0.05,]


plot1 <- ggplot(oosyears_spei_GLMM, aes(excluded_year, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI coefficient (scaled)", x = "Excluded year")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  

  scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  
  geom_hline(yintercept = full_coefs$effect, color = colors[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6], lty=2, lwd=1.5)+
  
  geom_errorbar(data=oosyears_spei_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosyears_spei_GLMM ,aes(excluded_year, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosyears_spei_GLMM ,aes(excluded_year, estimate, color = colors[2]),pch=20 ,size=4) +
  geom_point(data=oosyears_spei_GLMM_SIG ,aes(excluded_year, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2], "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")


# OOS plot

plots <- unique(LTabundance$plot)

oosplot_spei_GLMM <- data.frame(excluded_plot= as.factor(unique(LTabundance$plot)))

for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <- glmmTMB(crt_corr_biomass2 ~ scaled_year + scaled_length+  interval * scaled_sampling_temp + scaled_sampling_rain + spei72_2  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))  
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosplot_spei_GLMM$estimate[oosplot_spei_GLMM$excluded_plot == i] <- estimates[8,1]
  
  oosplot_spei_GLMM$estimate_upr[oosplot_spei_GLMM$excluded_plot == i] <- estimates[8,1] + 1.96*estimates[8,2]
  
  oosplot_spei_GLMM$estimate_lwr[oosplot_spei_GLMM$excluded_plot == i] <- estimates[8,1] - 1.96*estimates[8,2]
  
  oosplot_spei_GLMM$p[oosplot_spei_GLMM$excluded_plot == i] <- estimates[8,4]
  
  print(i)
}

oosplot_spei_GLMM_SIG <- oosplot_spei_GLMM[oosplot_spei_GLMM$p < 0.05,]


plot2 <- ggplot(oosplot_spei_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  geom_hline(yintercept = full_coefs$effect, color = colors[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6],  lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6],  lty=2, lwd=1.5)+
  
  geom_errorbar(data=oosplot_spei_GLMM, aes( x=excluded_plot , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosplot_spei_GLMM ,aes(excluded_plot, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosplot_spei_GLMM ,aes(excluded_plot, estimate, color = colors[2]),pch=20 ,size=4) +
  geom_point(data=oosplot_spei_GLMM_SIG ,aes(excluded_plot, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2], "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")

grid.arrange(plot1, plot2, ncol=2)
