# SUPPLEMENTARY R-SCRIPT TO:
# Long-term drought triggers severe declines in carabid beetles in a temperate forest
# F.Weiss, H.von Wehrden & A.Linde
# F.Weiss: ORCID 0000-0003-1078-1528

# PART II: Modelling abundance and biomass

# Data available at:  https://doi.org/10.48548/pubdata-46

#### Package list #####

library(lme4)
library(bbmle)
library(sjPlot)
library(gridExtra)
library(DHARMa)
library(gamm4)
library(glmmTMB)
library(cols4all)
library(ggplot2)
library(ggeffects)
library(ggpubr)
library(psych)


colors <- c4a("berlin", 7)

#### Loading data ####

# loading raw data
beetle_samples <- read.csv("published data/EWcarabids1999-2022_rawdata.csv")

# loading sampling meta data
sampling_meta <- read.csv2("sampling_meta_dwd_spei.csv")

# laoding carabid weights
# carabid weights were calculated using the approach by Weiss & Linde (2022, https://doi.org/10.1007/s10841-022-00391-6)
carabid_weights <- read.csv("published data/EWcarabids1999-2022_species_biomass.csv")


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

# interval
LTabundance$scaled_effort <- c(scale(LTabundance$sampling_length))

LTabundance$interval[LTabundance$interval== 1] <- "may"
LTabundance$interval[LTabundance$interval== 2] <- "june"
LTabundance$interval[LTabundance$interval== 3] <- "july"

LTabundance$interval <- as.factor(LTabundance$interval)
LTabundance$interval <- relevel(LTabundance$interval, ref = "may")

# abundances and biomass without Nebria brevicollis
LTabundance$sampling_abundance2 <- LTabundance$sampling_abundance - LTabundance$sampling_abundance_nb
LTabundance$sampling_biomass2 <- LTabundance$sampling_biomass - LTabundance$sampling_biomass_nb

# new trap ID for random intercept
LTabundance$trapID <- paste(LTabundance$plot,LTabundance$year,LTabundance$trap, sep="_")
LTabundance$trapID <- as.factor(LTabundance$trapID)
LTabundance$plot <- as.factor(LTabundance$plot)
LTabundance$site <-as.factor( LTabundance$site )

# scale precipitation variables

# scale temp by interval
intervals <- c(unique(LTabundance$interval))

LTabundance$scaled_sampling_temp <- NA
LTabundance$scaled_sampling_rain <- NA

for(i in intervals){
  LTabundance$scaled_sampling_temp[LTabundance$interval == i] <- c(scale(LTabundance$sampling_temp[LTabundance$interval == i]))
  LTabundance$scaled_sampling_rain[LTabundance$interval == i] <- c(scale(LTabundance$sampling_rain[LTabundance$interval == i]))
}

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


#### ABUNDANCE ####

### modelling ####

## GLMM ####
ab_glmm <- glmmTMB(sampling_abundance2 ~ scaled_year + interval + scaled_sampling_temp + scaled_sampling_rain + scaled_effort +I(scaled_effort^2)  + (1|year2) + (1|site/plot/trapID) , data= LTabundance, ziformula = ~ 0, family=nbinom2)


# diagnostics
dharma_sim1 <- simulateResiduals(fittedModel = ab_glmm, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)
testZeroInflation(dharma_sim1)
dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))

tab_model(ab_glmm, dv.labels = c("Abundance GLMM"))

summary(ab_glmm)

## GAMM ####

# fit with unconstrained k and not random intercept for year first to fix k to appropriate dimension in a second run (Knape 2016)
# theta of the respective glmm with Ziformula=~0, family=negbinom2

ab_gamm <- gamm4(sampling_abundance2 ~ s(scaled_year) + interval + scaled_sampling_temp + scaled_sampling_rain  + scaled_effort +I(scaled_effort^2),
                 random= ~ (1|site/plot/trapID), 
                 data= LTabundance, 
                 family=negbin(theta= 6.93),
                 REML=TRUE)

plot(ab_gamm$gam)

# k=6.33 -> 7

ab_gamm <- gamm4(sampling_abundance2 ~ s(scaled_year, k=7, fx=TRUE) + interval + scaled_sampling_temp + scaled_sampling_rain  + scaled_effort +I(scaled_effort^2),
                 random= ~ (1|year2) + (1|site/plot/trapID), 
                 data= LTabundance, 
                 family=negbin(theta= 6.93),
                 REML=TRUE)


# diagnostics
dharma_sim2 <- simulateResiduals(fittedModel = ab_gamm$mer, re.form= NULL)
plot(dharma_sim2)
testDispersion(dharma_sim2)
testZeroInflation(dharma_sim2)
#zeroinflation
dharma_sim2.1 = recalculateResiduals(dharma_sim2, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim2.1, time = unique( LTabundance$scaled_year))

tab_model(ab_gamm$mer, dv.labels = c("Abundance GAMM"))
plot(ab_gamm$gam)




## SPEI GLMM ####


# spei modelselection
# spei variables are already scaled

ab_glmm2 <- glmmTMB(sampling_abundance2 ~ scaled_year + interval + scaled_sampling_temp + scaled_sampling_rain + scaled_effort +I(scaled_effort^2)  + (1|year2) + (1|site/plot/trapID) , data= LTabundance, ziformula = ~ 0, family=nbinom2, REML = F)

m1 <- update(ab_glmm2,. ~ . + spei12)
m2 <- update(ab_glmm2,. ~ . + spei12_1)
m3 <- update(ab_glmm2,. ~ . + spei12_2)
m4 <- update(ab_glmm2,. ~ . + spei24)
m5 <- update(ab_glmm2,. ~ . + spei24_1)
m6 <- update(ab_glmm2,. ~ . + spei24_2)
m7 <- update(ab_glmm2,. ~ . + spei36)
m8 <- update(ab_glmm2,. ~ . + spei36_1)
m9 <- update(ab_glmm2,. ~ . + spei36_2)
m10 <- update(ab_glmm2,. ~ . + spei48)
m11 <- update(ab_glmm2,. ~ . + spei48_1)
m12 <- update(ab_glmm2,. ~ . + spei48_2)
m13 <- update(ab_glmm2,. ~ . + spei60)
m14 <- update(ab_glmm2,. ~ . + spei60_1)
m15 <- update(ab_glmm2,. ~ . + spei60_2)
m16 <- update(ab_glmm2,. ~ . + spei72)
m17 <- update(ab_glmm2,. ~ . + spei72_1)
m18 <- update(ab_glmm2,. ~ . + spei72_2)

AICtab(ab_glmm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18)

ab_spei_glmm <- update(ab_glmm,. ~ . + spei72_2)

dharma_sim1 <- simulateResiduals(fittedModel = ab_spei_glmm, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)
testZeroInflation(dharma_sim1)
dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))

tab_model(ab_spei_glmm, dv.labels = c("Abundance GLMM (SPEI)"))


### prediction ####


## GLMM ####

ggpreds1 <- ggpredict(ab_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
ggpreds1$year <- sort(unique(LTabundance$year2))  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  ggpreds1$predicted,
                           lwr   = ggpreds1$conf.low,
                           upr   = ggpreds1$conf.high   )


#percentage change
tsdata <- ts(data.frame(x1=c(preds_glmm$fit)), start=c(1999), frequency=1)

geometric.mean(abs(tsdata/stats::lag(tsdata,-1) - 1))                  



# bootstrap it ##########

# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model,  terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
return(ggpreds$predicted)
}

# quantile approach for CIs
abundance_boot1 <- bootMer(ab_glmm, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 128)



boot_rates <- c()


for (i in c(1:100)){
  tsdata <- ts(data.frame(x1=c(abundance_boot1$t[i,])), start=c(1999), frequency=1)
  boot_rates[i] <- mean(tsdata/stats::lag(tsdata,-1) - 1)
}

#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)










## GAMM ####

ggpreds3 <- ggpredict(ab_gamm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
ggpreds3$year <- sort(unique(LTabundance$year2)  )

preds_gamm <- data.frame(  scaled_year= ggpreds3$x,
                           year = ggpreds3$year,
                           fit  =  ggpreds3$predicted,
                           lwr   = ggpreds3$conf.low,
                           upr   = ggpreds3$conf.high   )


#percentage change
tsdata <- ts(data.frame(x1=c(preds_gamm$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)

# geometric mean
rates <- stats::lag(tsdata,+1)/tsdata
(geometric.mean(rates)-1)*100


# percent decline since 2015
(mean(preds_gamm$fit[preds_gamm$year %in% c(2016)]) - mean(preds_gamm$fit[preds_gamm$year %in% c(2022)])) /mean(preds_gamm$fit[preds_gamm$year %in% c(2016)])
# 70.5 % 



# boootstrap linear model from 2015-2022 ##########


# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model,  terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
return(ggpreds$predicted)
}

boot_data <- LTabundance[LTabundance$year >= 2015,]

ab_glmm2 <- glmmTMB(sampling_abundance2 ~ scaled_year + interval + scaled_sampling_temp + scaled_sampling_rain + scaled_effort +I(scaled_effort^2)  + (1|year2) + (1|site/plot/trapID) , data= boot_data, ziformula = ~ 0, family=nbinom2)

# quantile approach for CIs
abundance_boot3 <- bootMer(ab_glmm2, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 128)

boot_rates <- c()

#
for (i in c(1:1000)){
  boot_it <- abundance_boot3$t[i,]
  boot_rates[i] <- (boot_it[1]- boot_it[8])/boot_it[1]
}
#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)
#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)


## SPEI glmm ####

# predictions with annual spei values
ggpreds1 <- ggpredict(ab_spei_glmm, terms= c("scaled_year [all]", "spei72_2 [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))

ggpreds1$year <- LTabundance$year2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1$spei <- LTabundance$spei72_2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1 <- ggpreds1[ggpreds1$group == round(ggpreds1$spei, digits=3),]  

preds_glmm1 <- data.frame(  scaled_year= ggpreds1$x,
                            year = ggpreds1$year,
                            fit  =  ggpreds1$predicted,
                            lwr   = ggpreds1$conf.low,
                            upr   = ggpreds1$conf.high)

# predictions for background decline
ggpreds2 <- ggpredict(ab_spei_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(spei72_2 = 0, scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))

ggpreds2$year <- LTabundance$year2[match(ggpreds2$x, round(LTabundance$scaled_year, digits=2))]

preds_glmm2 <- data.frame(  scaled_year= ggpreds2$x,
                            year = ggpreds2$year,
                            fit  =  ggpreds2$predicted )



#percentage change
tsdata <- ts(data.frame(x1=c(preds_glmm2$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)


# bootstrap it ######
# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model,  terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0, spei72_2 = 0 ))
return(ggpreds$predicted)
}

# quantile approach for CIs
abundance_boot2 <- bootMer(ab_spei_glmm, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed=128)

abundance_boot2$t

boot_rates_spei_background <- c()


for (i in c(1:1000)){
  tsdata <- ts(data.frame(x1=c(abundance_boot2$t[i,])), start=c(1999), frequency=1)
  boot_rates_spei_background[i] <- mean(tsdata/stats::lag(tsdata,-1) - 1)
}

#lwr   
quantile(boot_rates_spei_background, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates_spei_background, probs=.975, na.rm=TRUE)




### plotting ####


set.seed(128)
plot1 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "Abundance", x = "")+
  ggtitle("A")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20 , 30), labels =c("0","  10", "  20", "    30"))+
  
  scale_x_continuous(limits=c(-1.7,2),labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=sampling_abundance2, scaled_year), alpha=0.1, color= "black", size=1.5)+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=4, lineend="round")+
  
  theme(legend.position="none")


set.seed(128)
plot2 <- ggplot(preds_gamm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "")+
  ggtitle("B")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20 , 30), labels =c("0","  10", "  20", "    30"))+
  
  scale_x_continuous(limits=c(-1.7,2), labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=sampling_abundance2, scaled_year), alpha=0.1, color= "black", size=1.5)+
  
  geom_ribbon(data = preds_gamm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[2])+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color =colors[2], lty=1, lwd=4, lineend="round")+
  
  theme(legend.position="none")


set.seed(128)
plot3 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "")+
  ggtitle("C")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20 , 30), labels =c("0","  10", "  20", "    30"))+
  
  scale_x_continuous(limits=c(-1.7,2),labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=sampling_abundance2, scaled_year), alpha=0.1, color= "black", size=1.5)+
  
  geom_line(data=preds_glmm2 ,aes(scaled_year, fit), color = colors[4],alpha=0.6, lty=1, lwd=3.5, lineend="round")+
  
  geom_ribbon(data = preds_glmm1, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[6])+
  
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=4, lineend="round")+
  
  theme(legend.position="none")


plot4 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "")+
  ggtitle("D")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(1,8), breaks = c(2,4,6), labels =c("   2", "   4", "      6"))+
  
  scale_x_continuous(limits=c(-1.7,2),labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=3.5, lineend="round")+
  geom_line(data=preds_glmm2 ,aes(scaled_year, fit), color = colors[4],alpha=0.5, lty=1, lwd=3.5, lineend="round")+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=3.5, lineend="round")+
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=3.5, alpha=0.8, lineend="round")+
  
  theme(legend.position="top")

#left <- text_grob("Abundance", size = 30, rot=90)
#bottom <- text_grob("Year", size = 30)
#grid.arrange(plot1, plot2, plot3,  plot4,  ncol=2, bottom=bottom, left=left)

# export 10x12 portrait



### sensitivity ####
## GLMM ####

#full model coefs

coefs <- coef(summary(ab_glmm))

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
  
  loop_mod <- glmmTMB(sampling_abundance2 ~ scaled_year + interval + scaled_sampling_temp + scaled_sampling_rain + scaled_effort +I(scaled_effort^2)  + (1|year2) + (1|site/plot/trapID) , data= loop_data, ziformula = ~ 0, family=nbinom2)
  
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
  #scale_y_continuous(limits = c(-0.43, 0.05))+
  
  scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))+
  
  
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
  
  loop_mod <- glmmTMB(sampling_abundance2 ~ scaled_year + interval + scaled_sampling_temp + scaled_sampling_rain + scaled_effort +I(scaled_effort^2)  + (1|year2) + (1|site/plot/trapID) , data= loop_data, ziformula = ~ 0, family=nbinom2)
  
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
  # scale_y_continuous(limits = c(-0.43, 0.05))+
  
  #scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))+
  
  
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


## GAMM ####

#model
ab_gamm <- gamm4(sampling_abundance2 ~ s(scaled_year, k=7, fx=TRUE) + interval + scaled_sampling_temp + scaled_sampling_rain  + scaled_effort +I(scaled_effort^2),
                 random= ~ (1|year2) + (1|site/plot/trapID), 
                 data= LTabundance, 
                 family=negbin(theta= 6.93),
                 REML=TRUE)

## oos year

years <- sort(c(unique(LTabundance$year)))

par(mfrow=c(6,4))

for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <-  gamm4(sampling_abundance2 ~ s(scaled_year, k=7, fx=TRUE) + interval + scaled_sampling_temp + scaled_sampling_rain  + scaled_effort +I(scaled_effort^2),
                     random= ~ (1|year2) + (1|site/plot/trapID), 
                     data= loop_data, 
                     family=negbin(theta= 6.93),
                     REML=TRUE)
  
  plot(loop_mod$gam, main = as.character(i), ylab="", xlab="", col=colors[2], shade=T, shade.col = colors[7], cex=2, ylim=c(-1.6, 1.1))
  
  print(i)
}

#export PDF 20x10 in

## oos plot

plots <- c(unique(LTabundance$plot))

par(mfrow=c(4,4))

for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <-  gamm4(sampling_abundance2 ~ s(scaled_year, k=7, fx=TRUE) + interval + scaled_sampling_temp + scaled_sampling_rain  + scaled_effort +I(scaled_effort^2),
                     random= ~ (1|year2) + (1|site/plot/trapID), 
                     data= loop_data, 
                     family=negbin(theta= 6.93),
                     REML=TRUE)
  
  plot(loop_mod$gam, main = as.character(i), ylab="", xlab="", col=colors[2], shade=T, shade.col = colors[7], cex=2, ylim=c(-1.6, 1.1))
  
  print(i)
}



## SPEI GLMM ####


#full model coefs

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
  
  loop_mod <- glmmTMB(sampling_abundance2 ~ scaled_year + interval + scaled_sampling_temp + scaled_sampling_rain + scaled_effort +I(scaled_effort^2)  + spei72_2  + (1|year2) + (1|site/plot/trapID) , data= loop_data, ziformula = ~ 0, family=nbinom2)
  
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
  #scale_y_continuous(limits = c(-0.43, 0.05))+
  
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
  
  loop_mod <- glmmTMB(sampling_abundance2 ~ scaled_year + interval + scaled_sampling_temp + scaled_sampling_rain + scaled_effort +I(scaled_effort^2)  + spei72_2  + (1|year2) + (1|site/plot/trapID) , data= loop_data, ziformula = ~ 0, family=nbinom2)
  
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
  #scale_y_continuous(limits = c(-0.43, 0.05))+
  
  #scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))+
  
  
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


#### BIOMASS ####

### modelling ####



# transform
LTabundance$crt_sampling_biomass2 <- (LTabundance$sampling_biomass2)^(1/3)

## GLMM ####

biomass_glmm <- glmmTMB(crt_sampling_biomass2 ~ scaled_year + scaled_effort + I(scaled_effort^2) +  interval + scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= LTabundance, family=gaussian(link=identity))

dharma_sim1 <- simulateResiduals(fittedModel = biomass_glmm, re.form= NULL)
plot(dharma_sim1)

testDispersion(dharma_sim1)
# DHARMa actually says no zero-inflation

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))
# fine

tab_model(biomass_glmm, dv.labels = c("Biomass GLMM"))


## GAMM ####

biomass_gamm <- gamm4(crt_sampling_biomass2 ~ s(scaled_year) + scaled_effort + I(scaled_effort^2) + interval + scaled_sampling_temp + scaled_sampling_rain ,
                      random= ~ (1|site/plot/trapID), 
                      data= LTabundance, 
                      family=gaussian(link=identity),
                      REML=TRUE)

plot(biomass_gamm$gam)
#7.85 -> k=8

biomass_gamm <- gamm4(crt_sampling_biomass2 ~ s(scaled_year, k=8, fx=TRUE )+ scaled_effort + I(scaled_effort^2) + interval + scaled_sampling_temp + scaled_sampling_rain ,
                      random= ~(1|year2) + (1|site/plot/trapID), 
                      data= LTabundance, 
                      family=gaussian(link=identity),
                      REML=TRUE)

plot(biomass_gamm$gam)
tab_model(biomass_gamm$mer, dv.labels = c("Biomass GAMM"))


dharma_sim2 <- simulateResiduals(fittedModel = biomass_gamm$mer, re.form= NULL)
plot(dharma_sim2)
testDispersion(dharma_sim2)

dharma_sim2.1 = recalculateResiduals(dharma_sim2, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim2.1, time = unique( LTabundance$scaled_year))
# With re.form= NULL -> significant strong negative autocorrelation
# DW = 2.7698, p-value = 0.04678
# With re.from= ~0  -> positive autocorrelation (not significant)
# DW = 2.4705, p-value = 0.2363
# but gamm4 does not allow for ar1 structures anyways 

summary(biomass_gamm$mer)
## SPEI GLMM ####

# spei modelselection
# spei variables are already scaled


biomass_glmm2 <- glmmTMB(crt_sampling_biomass2 ~ scaled_year + scaled_effort + I(scaled_effort^2) +  interval + scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= LTabundance, family=gaussian(link=identity), REML = F )

m1 <- update(biomass_glmm2,. ~ . + spei12)
m2 <- update(biomass_glmm2,. ~ . + spei12_1)
m3 <- update(biomass_glmm2,. ~ . + spei12_2)
m4 <- update(biomass_glmm2,. ~ . + spei24)
m5 <- update(biomass_glmm2,. ~ . + spei24_1)
m6 <- update(biomass_glmm2,. ~ . + spei24_2)
m7 <- update(biomass_glmm2,. ~ . + spei36)
m8 <- update(biomass_glmm2,. ~ . + spei36_1)
m9 <- update(biomass_glmm2,. ~ . + spei36_2)
m10 <- update(biomass_glmm2,. ~ . + spei48)
m11 <- update(biomass_glmm2,. ~ . + spei48_1)
m12 <- update(biomass_glmm2,. ~ . + spei48_2)
m13 <- update(biomass_glmm2,. ~ . + spei60)
m14 <- update(biomass_glmm2,. ~ . + spei60_1)
m15 <- update(biomass_glmm2,. ~ . + spei60_2)
m16 <- update(biomass_glmm2,. ~ . + spei72)
m17 <- update(biomass_glmm2,. ~ . + spei72_1)
m18 <- update(biomass_glmm2,. ~ . + spei72_2)

AICtab(biomass_glmm, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18)

biomass_spei_glmm <- update(biomass_glmm,. ~ . + spei72_2)

dharma_sim1 <- simulateResiduals(fittedModel = biomass_spei_glmm, re.form= NULL)
plot(dharma_sim1)
testDispersion(dharma_sim1)
# DHARMa actually says no zero-inflation

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTabundance$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTabundance$scaled_year))
# fine

tab_model(biomass_spei_glmm, dv.labels = c("Biomass GLMM (SPEI)"))


### prediction ####


## GLMM ####

ggpreds1 <- ggpredict(biomass_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort = 0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
ggpreds1$year <- sort(unique(LTabundance$year2))  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  (ggpreds1$predicted)^3,
                           lwr   = (ggpreds1$conf.low)^3,
                           upr   = (ggpreds1$conf.high)^3   )

#percentage change
tsdata <- ts(data.frame(x1=c(preds_glmm$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)


# bootstrap it ##########

# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model,  terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
return((ggpreds$predicted)^3)
}


# quantile approach for CIs
biomass_boot1 <- bootMer(biomass_glmm, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 128)

biomass_boot1$t

boot_rates <- c()


for (i in c(1:1000)){
  tsdata <- ts(data.frame(x1=c(biomass_boot1$t[i,])), start=c(1999), frequency=1)
  boot_rates[i] <- mean(tsdata/stats::lag(tsdata,-1) - 1)
}

#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)



## GAMM ####

ggpreds3 <- ggpredict(biomass_gamm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort = 0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
ggpreds3$year <- sort(unique(LTabundance$year2))  

preds_gamm <- data.frame(  scaled_year= ggpreds3$x,
                           year = ggpreds3$year,
                           fit  =  (ggpreds3$predicted)^3,
                           lwr   = (ggpreds3$conf.low)^3,
                           upr   = (ggpreds3$conf.high)^3   )

#percentage change
tsdata <- ts(data.frame(x1=c(preds_gamm$fit)), start=c(1999), frequency=1)

# geometric mean
rates <- stats::lag(tsdata,+1)/tsdata
(geometric.mean(rates)-1)*100



# percent decline comparable to repeated survey:
(mean(preds_gamm$fit[preds_gamm$year %in% c(1999,2001,2000)]) - mean(preds_gamm$fit[preds_gamm$year %in% c(2020,2021,2022)])) /mean(preds_gamm$fit[preds_gamm$year %in% c(1999,2001,2000)])
# 84.1 % 

# percent decline since 2015
(mean(preds_gamm$fit[preds_gamm$year %in% c(2016)]) - mean(preds_gamm$fit[preds_gamm$year %in% c(2022)])) /mean(preds_gamm$fit[preds_gamm$year %in% c(2016)])
# 88.7 % 

# boootstrap linear model from 2015-2022 ###########

# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model,  terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))
return(ggpreds$predicted)
}

boot_data <- LTabundance[LTabundance$year >= 2015,]

biomass_glmm2 <- glmmTMB(crt_sampling_biomass2 ~ scaled_year + scaled_effort + I(scaled_effort^2) +  interval + scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= boot_data, family=gaussian(link=identity))


biomass_boot3 <- bootMer(biomass_glmm2, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed= 128)

boot_rates <- c()

biomass_boot3$t

#
for (i in c(1:1000)){
  boot_it <- biomass_boot3$t[i,]^3
  boot_rates[i] <- (boot_it[2]- boot_it[8])/boot_it[2]
}
#lwr   
quantile(boot_rates, probs=.025, na.rm=TRUE)
#upr
quantile(boot_rates, probs=.975, na.rm=TRUE)

## SPEI glmm ####

ggpreds1 <- ggpredict(biomass_spei_glmm, terms= c("scaled_year [all]", "spei72_2 [all]"),ci.lvl = .95, type="fixed", condition = c( scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))

ggpreds1$year <- LTabundance$year2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1$spei <- LTabundance$spei72_2[match(ggpreds1$x, round(LTabundance$scaled_year, digits=2))]
ggpreds1 <- ggpreds1[ggpreds1$group == round(ggpreds1$spei, digits=3),]  

preds_glmm1 <- data.frame(  scaled_year= ggpreds1$x,
                            year = ggpreds1$year,
                            fit  =  (ggpreds1$predicted)^3,
                            lwr   = (ggpreds1$conf.low)^3,
                            upr   = (ggpreds1$conf.high)^3   )


# predictions for background decline
ggpreds2 <- ggpredict(biomass_spei_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(spei72_2 = 0, scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0 ))

ggpreds2$year <- LTabundance$year2[match(ggpreds2$x, round(LTabundance$scaled_year, digits=2))]

preds_glmm2 <- data.frame(  scaled_year= ggpreds2$x,
                            year = ggpreds2$year,
                            fit  =  (ggpreds2$predicted)^3 )



tsdata <- ts(data.frame(x1=c(preds_glmm2$fit)), start=c(1999), frequency=1)
mean(tsdata/stats::lag(tsdata,-1) - 1)



# bootstrap it ######
# first modify the ggpredict function to return only a numeric vector
ggpredict_boot <- function(model){  ggpreds <- ggpredict(model,  terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed", condition = c(scaled_effort=0, interval= "june", scaled_sampling_rain=0,scaled_sampling_temp= 0, spei72_2 = 0 ))
return((ggpreds$predicted)^3)
}

# quantile approach for CIs
biomass_boot2 <- bootMer(biomass_spei_glmm, ggpredict_boot  , nsim=1000, .progress = "txt", re.form = NA, seed=128)

biomass_boot2$t

boot_rates_spei_background <- c()


for (i in c(1:1000)){
  tsdata <- ts(data.frame(x1=c(biomass_boot2$t[i,])), start=c(1999), frequency=1)
  boot_rates_spei_background[i] <- mean(tsdata/stats::lag(tsdata,-1) - 1)
}

#lwr   
quantile(boot_rates_spei_background, probs=.025, na.rm=TRUE)

#upr
quantile(boot_rates_spei_background, probs=.975, na.rm=TRUE)





### plotting #####

set.seed(129)
plot5 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "Biomass (mg)", x = "")+
  ggtitle("E")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(0, 8000), breaks = c(0, 2500, 5000, 7500))+
  
  scale_x_continuous(limits=c(-1.7,2),labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=sampling_biomass2, scaled_year), alpha=0.1, color= "black", size=1.5)+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=4, lineend="round")+
  
  theme(legend.position="none")

set.seed(129)
plot6 <- ggplot(preds_gamm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "")+
  ggtitle("F")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(0, 8000), breaks = c(0, 2500, 5000, 7500))+
  
  scale_x_continuous(limits=c(-1.7,2), labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=sampling_biomass2, scaled_year), alpha=0.1, color= "black", size=1.5)+
  
  geom_ribbon(data = preds_gamm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[2])+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color =colors[2], lty=1, lwd=4, lineend="round")+
  
  theme(legend.position="none")

set.seed(129)
plot7 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "")+
  ggtitle("G")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(0, 8000), breaks = c(0, 2500, 5000, 7500))+
  
  scale_x_continuous(limits=c(-1.7,2),labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_jitter(data=LTabundance, aes(y=sampling_biomass2, scaled_year), alpha=0.1, color= "black", size=1.5)+
  
  geom_line(data=preds_glmm2 ,aes(scaled_year, fit), color = colors[4],alpha=0.6, lty=1, lwd=3.5, lineend="round")+
  
  geom_ribbon(data = preds_glmm1, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill=colors[6])+
  
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=4, lineend="round")+
  
  theme(legend.position="none")

set.seed(129)
plot8 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  
  labs(y = "", x = "")+
  ggtitle("H")+
  theme(plot.title = element_text(hjust =0.9,vjust = -4, size=30))+
  scale_y_continuous(limits = c(0, 1600), breaks = c(0,500,1000,1500))+
  
  scale_x_continuous(limits=c(-1.7,2),labels= c( "1999", "", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "","","2022"), breaks = c(unique(preds_glmm$scaled_year)) )+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=3.5, lineend="round")+
  geom_line(data=preds_glmm2 ,aes(scaled_year, fit), color = colors[4],alpha=0.5, lty=1, lwd=3.5, lineend="round")+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=3.5, lineend="round")+
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=3.5, alpha=0.8, lineend="round")+
  
  theme(legend.position="top")

#left <- text_grob("Biomass (mg)", size = 30, rot=90)
#bottom <- text_grob("Year", size = 30)
#grid.arrange(plot1, plot2, plot3,  plot4,  ncol=2, bottom=bottom, left=left)

# export 10x12 portrait

#full plot

# 20x12
bottom <- text_grob("Year", size = 30)
grid.arrange(plot1, plot2, plot3,  plot4,
             plot5, plot6, plot7, plot8, ncol=4, bottom=bottom)


### sensitivity #########
## GLMM ####

#full model coefs

coefs <- coef(summary(biomass_glmm))

estimates <- coefs$cond

full_coefs <- data.frame(   effect  = estimates[2,1],
                            lwr   = estimates[2,1] + 1.96*estimates[2,2],
                            upr   = estimates[2,1] - 1.96*estimates[2,2],
                            p     = estimates[2,4])


#OOS years
years <- sort(c(unique(LTabundance$year)))

oosyear_biomass_GLMM <- data.frame(excluded_year= c( 1999:2022))


for(i in years){
  loop_data <- LTabundance[LTabundance$year != i,]
  
  loop_mod <-  glmmTMB(crt_sampling_biomass2 ~ scaled_year + scaled_effort + I(scaled_effort^2) +  interval + scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))
  
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
  #scale_y_continuous(limits = c(-1.1, 0.1))+
  
  scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))+
  
  
  geom_hline(yintercept = full_coefs$effect, color = colors[6], lty=1, lwd=1.5)+
  geom_hline(yintercept = full_coefs$lwr, color = colors[6], lty=2, lwd=1.5)+
  geom_hline(yintercept = full_coefs$upr, color = colors[6], lty=2, lwd=1.5)+
  
  geom_errorbar(data=oosyear_biomass_GLMM, aes( x=excluded_year , ymin=estimate_lwr, ymax=estimate_upr), width=0.4, size=1, color = "black")+
  
  geom_point(data=oosyear_biomass_GLMM ,aes(excluded_year, estimate, color = "black"),pch=20 ,size=5) +
  geom_point(data=oosyear_biomass_GLMM ,aes(excluded_year, estimate, color = colors[2]),pch=20 ,size=4) +
  geom_point(data=oosyear_biomass_GLMM_SIG ,aes(excluded_year, estimate, color = "black"),pch=8 ,size=4) +
  
  scale_color_manual(values = c(colors[2], "black"), name = "random year effect", labels=c("without", "with" ))+
  
  theme(legend.position="none")

## OOS plots

plots <- c(unique(LTabundance$plot))

oosplot_biomass_GLMM <- data.frame(excluded_plot= as.factor(unique(LTabundance$plot)))


for(i in plots){
  loop_data <- LTabundance[LTabundance$plot != i,]
  
  loop_mod <-  glmmTMB(crt_sampling_biomass2 ~ scaled_year + scaled_effort + I(scaled_effort^2) +  interval + scaled_sampling_temp + scaled_sampling_rain  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))
  
  
  coefs<-coef(summary(loop_mod))
  estimates <- coefs$cond
  
  oosplot_biomass_GLMM$estimate[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,1]
  
  oosplot_biomass_GLMM$estimate_upr[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,1] + 1.96*estimates[2,2]
  
  oosplot_biomass_GLMM$estimate_lwr[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,1] - 1.96*estimates[2,2]
  
  oosplot_biomass_GLMM$p[oosplot_biomass_GLMM$excluded_plot == i] <- estimates[2,4]
  
  print(i)
}

oosplot_biomass_GLMM_SIG <- oosplot_biomass_GLMM[oosplot_biomass_GLMM$p < 0.05,]

# plot

plot2 <- ggplot(oosplot_biomass_GLMM, aes(excluded_plot, estimate)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "", x = "Excluded plot")+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.05,vjust=-4))+  
  #scale_y_continuous(limits = c(-1.1, 0.1))+
  
  #scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))+
  
  
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
# export PDF 6x12 in

## GAMM ####

## oos year

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
#export PDF A4 

## oos plot

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
#export PDF 8.3x8.3 in

par(mfrow=c(1,1))


## SPEI GLMM ####


#full model coefs

coefs <- coef(summary(biomass_spei_glmm))

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
  
  loop_mod <-  glmmTMB(crt_sampling_biomass2 ~ scaled_year + scaled_effort + I(scaled_effort^2) +  interval + scaled_sampling_temp + scaled_sampling_rain + spei72_2  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))
  
  
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
  #scale_y_continuous(limits = c(-0.43, 0.05))+
  
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
  
  loop_mod <-  glmmTMB(crt_sampling_biomass2 ~ scaled_year + scaled_effort + I(scaled_effort^2) +  interval + scaled_sampling_temp + scaled_sampling_rain + spei72_2  + (1|year2) + (1|site/plot/trapID), data= loop_data, family=gaussian(link=identity))  
  
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
  #scale_y_continuous(limits = c(-0.43, 0.05))+
  
  #scale_x_continuous(labels= c( "", "2000", "", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(1999:2022) )+
  #theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=0.4))+
  
  
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
# export PDF 6x12 in



