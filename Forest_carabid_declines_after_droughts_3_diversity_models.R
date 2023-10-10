# SUPPLEMENTARY R-SCRIPT TO:
# Long-term drought triggers severe declines in carabid beetles in a temperate forest
# F.Weiss, H.von Wehrden & A.Linde
# F.Weiss: ORCID 0000-0003-1078-1528

# PART III: Modelling diversity

# Data available at:  https://doi.org/10.48548/pubdata-46

#### Package list ####

library(data.table)
library(vegan)
library(Rarefy)
library(dplyr)
library(ggplot2)
library(iNEXT)


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

library(dplyr)

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

# remove N.brevicollis
beetle_samples3 <- beetle_samples3[beetle_samples3$species != "Nebria brevicollis",]

# reformat to wide table

# based on sampling abundances
wide_samples <- beetle_samples3[, c(1,5,6)]

wide_samples  <- reshape(wide_samples, idvar = "species", timevar = "sample_id2", direction = "wide")
wide_samples[is.na(wide_samples)] <- 0

wide_samples2 <- wide_samples[,-1] 
row.names(wide_samples2) <- wide_samples$species

#wide_samples3 <- as.list(wide_samples2)

# convert to list of vectors (type 'abundance' in iNEXT framework)
# wide_samples3 <- setNames(split(wide_samples2, seq(nrow(wide_samples2))), rownames(wide_samples2))



#### how large are sample? ####
beetle_samples4 <- aggregate(cbind(beetle_samples3$sampling_abundance), by= list( beetle_samples3$sample_id2), FUN=sum  )

# min = 6
min(beetle_samples4$V1)

# median = 108
median(beetle_samples4$V1)
nrow(beetle_samples4[beetle_samples4$V1 < 31,])

# max = 380
max(beetle_samples4$V1)








## simpson ####


stand_simpson <- iNEXT(wide_samples2, q=c(2), datatype="abundance",knots = 60, endpoint=62)

stand_simpson <- data.frame(stand_simpson$iNextEst$coverage_based)

# select rows with coverage closest to X, Standardize by coverage.
# create a loop to find coverage sweet spot for the sample, test range 0.70 to 0.99

ids <- c(unique(stand_simpson$Assemblage))
cov <- c(seq(from = 0.70, to = 0.99, by = 0.01))

for(y in cov){
  
  test <- stand_simpson[1,]
  stand_simpson$SC_2 <- abs(stand_simpson$SC - y)
  
  for(i in ids){
    loopdata <- stand_simpson[stand_simpson$Assemblage == i,]
    loopdata2 <- loopdata[loopdata$SC_2 == min(loopdata$SC_2),] 
    
    test <- rbind(test,loopdata2)
  }
  
  test <- test[-1,]
  print( paste(y, "->",max(test$SC_2)))
  
}

# repeat for coverage 0.86

simpson86 <- stand_simpson[1,]
stand_simpson$SC_2 <- abs(stand_simpson$SC - 0.86)

for(i in ids){
  loopdata <- stand_simpson[stand_simpson$Assemblage == i,]
  loopdata2 <- loopdata[loopdata$SC_2 == min(loopdata$SC_2),] 
  
  simpson86 <- rbind(simpson86,loopdata2)
}

simpson86 <- simpson86[-1,]


## richness ####


stand_richness <- iNEXT(wide_samples2, q=c(0), datatype="abundance",knots = 60, endpoint=62)

stand_richness <- data.frame(stand_richness$iNextEst$coverage_based)

# select rows with coverage closest to X, Standardize by coverage.
# create a loop to find coverage sweet spot for the sample, test range 0.70 to 0.99

ids <- c(unique(stand_richness$Assemblage))
cov <- c(seq(from = 0.70, to = 0.99, by = 0.01))

for(y in cov){
  
  test <- stand_richness[1,]
  stand_richness$SC_2 <- abs(stand_richness$SC - y)
  
  for(i in ids){
    loopdata <- stand_richness[stand_richness$Assemblage == i,]
    loopdata2 <- loopdata[loopdata$SC_2 == min(loopdata$SC_2),] 
    
    test <- rbind(test,loopdata2)
  }
  
  test <- test[-1,]
  print( paste(y, "->",max(test$SC_2)))
  
}

# repeat for coverage 0.86

richness86 <- stand_richness[1,]
stand_richness$SC_2 <- abs(stand_richness$SC - 0.86)

for(i in ids){
  loopdata <- stand_richness[stand_richness$Assemblage == i,]
  loopdata2 <- loopdata[loopdata$SC_2 == min(loopdata$SC_2),] 
  
  richness86 <- rbind(richness86,loopdata2)
}

richness86 <- richness86[-1,]




## turnover ####

# turnover

# calculated for each plot and year the Jaccard Index between first year of sampling
# maybe rather get average of Jaccard Index towards first 2 years available for each plot

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

## observed richness ####

diversity_samples <- aggregate(cbind(beetle_samples3$species), by= list(beetle_samples3$sample_id2, beetle_samples3$year, beetle_samples3$plot, beetle_samples3$site), function(x) length(unique(x)))

diversity_samples<-diversity_samples %>% 
  rename(
    sample_id2 = Group.1,
    year = Group.2,
    plot = Group.3,
    site = Group.4,
    obs_richness = V1
  )

## combine in one DF #####

turnover_data$sample_id2 <- paste(turnover_data$year,turnover_data$plot, sep = "_")

diversity_samples$jaccard <- turnover_data$jaccard[match(diversity_samples$sample_id2, turnover_data$sample_id2)]


richness86$sample_id2 <- substr(richness86$Assemblage,20,26)
diversity_samples$stand_richness <- richness86$qD[match(diversity_samples$sample_id2, richness86$sample_id2)]
diversity_samples$stand_richness_low <- richness86$qD.LCL[match(diversity_samples$sample_id2, richness86$sample_id2)]
diversity_samples$stand_richness_high <- richness86$qD.UCL[match(diversity_samples$sample_id2, richness86$sample_id2)]

simpson86$sample_id2 <- substr(simpson86$Assemblage,20,26)
diversity_samples$stand_simpson <- simpson86$qD[match(diversity_samples$sample_id2, simpson86$sample_id2)]
diversity_samples$stand_simpson_low <- simpson86$qD.LCL[match(diversity_samples$sample_id2, simpson86$sample_id2)]
diversity_samples$stand_simpson_high <- simpson86$qD.UCL[match(diversity_samples$sample_id2, simpson86$sample_id2)]

diversity_samples$stand_evenness <- diversity_samples$stand_simpson/diversity_samples$stand_richness

LTdiversity <- diversity_samples

# year variable
LTdiversity$year <- as.integer(LTdiversity$year)
LTdiversity$year2 <- as.factor(LTdiversity$year)
LTdiversity$scaled_year <- c(scale(LTdiversity$year))

LTdiversity$plot <- as.factor(LTdiversity$plot)
LTdiversity$site <-as.factor( LTdiversity$site )



#### modelling #### 

## standardized richness ####


## glmm
richness_glmm <- glmmTMB(stand_richness ~ scaled_year + (1|year2) + (1|site/plot) , family = Gamma(link="log"), data=LTdiversity)


dharma_sim1 <- simulateResiduals(fittedModel = richness_glmm, re.form= NULL)
plot(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))


tab_model(richness_glmm, dv.labels = c("Standardized Richness GLMM"))



## gamm

richness_gamm <- gamm4(stand_richness ~ s(scaled_year),
                       random=  ~ (1|site/plot), 
                       data= LTdiversity, 
                       family= gaussian(link = "identity"),
                       REML=TRUE)

plot(species_gamm$gam)

richness_gamm <- gamm4(stand_richness ~ s(scaled_year, k=3,fx=T),
                       random=  ~ (1|year2) + (1|site/plot), 
                       data= LTdiversity, 
                       family=gaussian(link = "identity"),
                       REML=TRUE)

plot(richness_gamm$gam)



dharma_sim1 <- simulateResiduals(fittedModel = richness_gamm$mer, re.form= NULL)
plot(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(richness_gamm$mer)

## SPEI glmm

m7 <- update(richness_glmm,. ~ . + spei12)
m8 <- update(richness_glmm,. ~ . + spei12_1)
m9 <- update(richness_glmm,. ~ . + spei12_2)
m10 <- update(richness_glmm,. ~ . + spei24)
m11 <- update(richness_glmm,. ~ . + spei24_1)
m12 <- update(richness_glmm,. ~ . + spei24_2)
m13 <- update(richness_glmm,. ~ . + spei36)
m14 <- update(richness_glmm,. ~ . + spei36_1)
m15 <- update(richness_glmm,. ~ . + spei36_2)
m16 <- update(richness_glmm,. ~ . + spei48)
m17 <- update(richness_glmm,. ~ . + spei48_1)
m18 <- update(richness_glmm,. ~ . + spei48_2)
m19 <- update(richness_glmm,. ~ . + spei60)
m20 <- update(richness_glmm,. ~ . + spei60_1)
m21 <- update(richness_glmm,. ~ . + spei60_2)
m22 <- update(richness_glmm,. ~ . + spei72)
m23 <- update(richness_glmm,. ~ . + spei72_1)
m24 <- update(richness_glmm,. ~ . + spei72_2)

AICtab(richness_glmm, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24)

richness_spei_glmm <- m23

dharma_sim1 <- simulateResiduals(fittedModel = richess_spei_glmm, re.form= NULL)
plot(dharma_sim1)
dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(richess_spei_glmm)




## standardized simpson #### 


## glmm
simpson_glmm <- glmmTMB(stand_simpson ~ scaled_year + (1|year2) + (1|site/plot) , family = Gamma(link="log"), data=LTdiversity)


dharma_sim1 <- simulateResiduals(fittedModel = simpson_glmm, re.form= NULL)
plot(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(simpson_glmm, dv.labels = c("Standardized Simpson GLMM"))



## gamm

simpson_gamm <- gamm4(stand_simpson ~ s(scaled_year),
                      random=  ~ (1|site/plot), 
                      data= LTdiversity, 
                      family= Gamma(link="log"),
                      REML=TRUE)

plot(simpson_gamm$gam)

# estimated as linear... no GAMM necessary

## SPEI glmm

m7 <- update(simpson_glmm,. ~ . + spei12)
m8 <- update(simpson_glmm,. ~ . + spei12_1)
m9 <- update(simpson_glmm,. ~ . + spei12_2)
m10 <- update(simpson_glmm,. ~ . + spei24)
m11 <- update(simpson_glmm,. ~ . + spei24_1)
m12 <- update(simpson_glmm,. ~ . + spei24_2)
m13 <- update(simpson_glmm,. ~ . + spei36)
m14 <- update(simpson_glmm,. ~ . + spei36_1)
m15 <- update(simpson_glmm,. ~ . + spei36_2)
m16 <- update(simpson_glmm,. ~ . + spei48)
m17 <- update(simpson_glmm,. ~ . + spei48_1)
m18 <- update(simpson_glmm,. ~ . + spei48_2)
m19 <- update(simpson_glmm,. ~ . + spei60)
m20 <- update(simpson_glmm,. ~ . + spei60_1)
m21 <- update(simpson_glmm,. ~ . + spei60_2)
m22 <- update(simpson_glmm,. ~ . + spei72)
m23 <- update(simpson_glmm,. ~ . + spei72_1)
m24 <- update(simpson_glmm,. ~ . + spei72_2)

AICtab(simpson_glmm, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24)

# SPEI not a meaningful addition



## evenness ####


## glmm
evenness_glmm <- glmmTMB(stand_evenness ~ scaled_year + (1|year2) + (1|site/plot), family = beta_family(link="logit"), data=LTdiversity)


dharma_sim1 <- simulateResiduals(fittedModel = evenness_glmm, re.form= NULL)
plot(dharma_sim1)
dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(evenness_glmm)



evenness_gamm <- gamm4(stand_evenness ~ s(scaled_year),
                       random=  ~ (1|site/plot), 
                       data= LTdiversity, 
                       family= gaussian(link="identity"),
                       REML=TRUE)

plot(evenness_gamm$gam)

evenness_gamm <- gamm4(stand_evenness ~ s(scaled_year,k=3,fx=T),
                       random=  ~ (1|year2) + (1|site/plot), 
                       data= LTdiversity, 
                       family= gaussian(link="identity"),
                       REML=TRUE)


dharma_sim1 <- simulateResiduals(fittedModel = evenness_gamm$mer, re.form= NULL)
plot(dharma_sim1)

dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(evenness_gamm$mer)


## SPEI GLMM


m7 <- update(evenness_glmm,. ~ . + spei12)
m8 <- update(evenness_glmm,. ~ . + spei12_1)
m9 <- update(evenness_glmm,. ~ . + spei12_2)
m10 <- update(evenness_glmm,. ~ . + spei24)
m11 <- update(evenness_glmm,. ~ . + spei24_1)
m12 <- update(evenness_glmm,. ~ . + spei24_2)
m13 <- update(evenness_glmm,. ~ . + spei36)
m14 <- update(evenness_glmm,. ~ . + spei36_1)
m15 <- update(evenness_glmm,. ~ . + spei36_2)
m16 <- update(evenness_glmm,. ~ . + spei48)
m17 <- update(evenness_glmm,. ~ . + spei48_1)
m18 <- update(evenness_glmm,. ~ . + spei48_2)
m19 <- update(evenness_glmm,. ~ . + spei60)
m20 <- update(evenness_glmm,. ~ . + spei60_1)
m21 <- update(evenness_glmm,. ~ . + spei60_2)
m22 <- update(evenness_glmm,. ~ . + spei72)
m23 <- update(evenness_glmm,. ~ . + spei72_1)
m24 <- update(evenness_glmm,. ~ . + spei72_2)

AICtab(evenness_glmm, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24)

evenness_spei_glmm <- m19

dharma_sim1 <- simulateResiduals(fittedModel = evenness_spei_glmm, re.form= NULL)
plot(dharma_sim1)
dharma_sim1.1 = recalculateResiduals(dharma_sim1, group = LTdiversity$scaled_year)
testTemporalAutocorrelation(dharma_sim1.1, time = unique( LTdiversity$scaled_year))

tab_model(evenness_spei_glmm)



## turnover/jaccard ####

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

#gamm
jaccard_gamm <- gamm4(jaccard ~ s(scaled_year),
                      random=  ~ (1|site/plot), 
                      data= LTdiversity, 
                      family= gaussian(link="identity"),
                      REML=TRUE)

plot(jaccard_gamm$gam)

jaccard_gamm <- gamm4(jaccard ~ s(scaled_year,k=4,fx=T),
                      random=  ~ (1|year2) + (1|site/plot), 
                      data= LTdiversity, 
                      family= gaussian(link="identity"),
                      REML=TRUE)


dharma_sim1 <- simulateResiduals(fittedModel = jaccard_gamm$mer, re.form= NULL)
plot(dharma_sim1)

tab_model(jaccard_gamm$mer)
plot(jaccard_gamm$gam)




## SPEI glmm


m7 <- update(jaccard_glmm,. ~ . + spei12)
m8 <- update(jaccard_glmm,. ~ . + spei12_1)
m9 <- update(jaccard_glmm,. ~ . + spei12_2)
m10 <- update(jaccard_glmm,. ~ . + spei24)
m11 <- update(jaccard_glmm,. ~ . + spei24_1)
m12 <- update(jaccard_glmm,. ~ . + spei24_2)
m13 <- update(jaccard_glmm,. ~ . + spei36)
m14 <- update(jaccard_glmm,. ~ . + spei36_1)
m15 <- update(jaccard_glmm,. ~ . + spei36_2)
m16 <- update(jaccard_glmm,. ~ . + spei48)
m17 <- update(jaccard_glmm,. ~ . + spei48_1)
m18 <- update(jaccard_glmm,. ~ . + spei48_2)
m19 <- update(jaccard_glmm,. ~ . + spei60)
m20 <- update(jaccard_glmm,. ~ . + spei60_1)
m21 <- update(jaccard_glmm,. ~ . + spei60_2)
m22 <- update(jaccard_glmm,. ~ . + spei72)
m23 <- update(jaccard_glmm,. ~ . + spei72_1)
m24 <- update(jaccard_glmm,. ~ . + spei72_2)

AICtab(jaccard_glmm, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23, m24)


jaccard_spei_glmm <- m24
# no support.





#### visualizing ####

berlin <- c4a("berlin", 7)


## richness ####
tab_model(richness_glmm)

ggpreds1 <- ggpredict(richness_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds1$year <- sort(unique(LTdiversity$year))  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  ggpreds1$predicted,
                           lwr   = ggpreds1$conf.low,
                           upr   = ggpreds1$conf.high   )


ggpreds2 <- ggpredict(richness_gamm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds2$year <- sort(unique(LTdiversity$year)  )

preds_gamm <- data.frame(  scaled_year= ggpreds2$x,
                           year = ggpreds2$year,
                           fit  =  ggpreds2$predicted,
                           lwr   = ggpreds2$conf.low,
                           upr   = ggpreds2$conf.high   )




ggpreds1 <- ggpredict(richness_spei_glmm, terms= c("scaled_year [all]", "spei72_1 [all]"),ci.lvl = .95, type="fixed")

ggpreds1$year <- LTdiversity$year2[match(ggpreds1$x, round(LTdiversity$scaled_year, digits=2))]
ggpreds1$spei <- LTdiversity$spei72_1[match(ggpreds1$x, round(LTdiversity$scaled_year, digits=2))]
ggpreds1 <- ggpreds1[ggpreds1$group == round(ggpreds1$spei, digits=3),]  

preds_glmm1 <- data.frame(  scaled_year= ggpreds1$x,
                            year = ggpreds1$year,
                            fit  =  ggpreds1$predicted,
                            lwr   = ggpreds1$conf.low,
                            upr   = ggpreds1$conf.high)

# predictions for background decline
ggpreds2 <- ggpredict(richness_spei_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")

ggpreds2$year <- LTdiversity$year2[match(ggpreds2$x, round(LTdiversity$scaled_year, digits=2))]

preds_glmm2 <- data.frame(  scaled_year= ggpreds2$x,
                            year = ggpreds2$year,
                            fit  =  ggpreds2$predicted,
                            lwr   = ggpreds2$conf.low,
                            upr   = ggpreds2$conf.high)




plot1 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "Standardized Richness", x = "")+
  ggtitle("A")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, stand_richness),alpha=0.3, color = colors[4], size=2)+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=2, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")


plot2 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "", x = "")+
  ggtitle("B")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, stand_richness),alpha=0.3, color = colors[4], size=2)+
  
  geom_ribbon(data = preds_gamm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[2])+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")

plot3 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "", x = "")+
  ggtitle("C")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, stand_richness),alpha=0.3, color = colors[4], size=2)+
  
  #geom_ribbon(data = preds_glmm2, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[4])+
  
  geom_line(data=preds_glmm2 ,aes(scaled_year, fit), alpha=0.8, color = colors[4], lty=1, lwd=2.5, lineend="round")+
  
  geom_ribbon(data = preds_glmm1, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[6])+
  
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")


## simpson ####
tab_model(simpson_glmm)

ggpreds1 <- ggpredict(simpson_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds1$year <- sort(unique(LTdiversity$year))  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  ggpreds1$predicted,
                           lwr   = ggpreds1$conf.low,
                           upr   = ggpreds1$conf.high   )





plot4 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "Standardized Simpson", x = "")+
  ggtitle("D")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  scale_y_continuous(limits = c(1, 12), breaks = c(3,6, 9, 12))+
  
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, stand_simpson),alpha=0.3, color = colors[4], size=2)+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")


plot5 <- ggplot() +
  annotate("text", x = 10,  y = 10,
           size = 6,
           label = "No non-linear \n trend detected") + theme_void()


plot6 <- ggplot() +
  annotate("text", x = 10,  y = 10,
           size = 6,
           label = "Drought index not \n a meaningful predictor") + theme_void()



## evenness ####
tab_model(evenness_spei_glmm)

ggpreds1 <- ggpredict(evenness_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds1$year <- sort(unique(LTdiversity$year))  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  ggpreds1$predicted,
                           lwr   = ggpreds1$conf.low,
                           upr   = ggpreds1$conf.high   )


ggpreds2 <- ggpredict(evenness_gamm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds2$year <- sort(unique(LTdiversity$year)  )

preds_gamm <- data.frame(  scaled_year= ggpreds2$x,
                           year = ggpreds2$year,
                           fit  =  ggpreds2$predicted,
                           lwr   = ggpreds2$conf.low,
                           upr   = ggpreds2$conf.high   )




ggpreds1 <- ggpredict(evenness_spei_glmm, terms= c("scaled_year [all]", "spei60 [all]"),ci.lvl = .95, type="fixed")

ggpreds1$year <- LTdiversity$year2[match(ggpreds1$x, round(LTdiversity$scaled_year, digits=2))]
ggpreds1$spei <- LTdiversity$spei60[match(ggpreds1$x, round(LTdiversity$scaled_year, digits=2))]
ggpreds1 <- ggpreds1[ggpreds1$group == round(ggpreds1$spei, digits=2),]  

preds_glmm1 <- data.frame(  scaled_year= ggpreds1$x,
                            year = ggpreds1$year,
                            fit  =  ggpreds1$predicted,
                            lwr   = ggpreds1$conf.low,
                            upr   = ggpreds1$conf.high)

# predictions for background decline
ggpreds2 <- ggpredict(evenness_spei_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")

ggpreds2$year <- LTdiversity$year2[match(ggpreds2$x, round(LTdiversity$scaled_year, digits=2))]

preds_glmm2 <- data.frame(  scaled_year= ggpreds2$x,
                            year = ggpreds2$year,
                            fit  =  ggpreds2$predicted,
                            lwr   = ggpreds2$conf.low,
                            upr   = ggpreds2$conf.high)




plot7 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "Standardized Evenness", x = "")+
  ggtitle("E")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, stand_evenness),alpha=0.3, color = colors[4], size=2)+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")


plot8 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "", x = "")+
  ggtitle("F")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, stand_evenness),alpha=0.3, color = colors[4], size=2)+
  
  geom_ribbon(data = preds_gamm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[2])+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")

plot9 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "", x = "")+
  ggtitle("G")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c( "2001", "", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(preds_richness$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, stand_evenness),alpha=0.3, color = colors[4], size=2)+
  
  #geom_ribbon(data = preds_glmm2, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[4])+
  
  geom_line(data=preds_glmm2 ,aes(scaled_year, fit), alpha=0.8, color = colors[4], lty=2, lwd=2.5, lineend="round")+
  
  geom_ribbon(data = preds_glmm1, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[6])+
  
  geom_line(data=preds_glmm1 ,aes(scaled_year, fit), color = colors[6], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")




## turnover ####
tab_model(jaccard_glmm)

ggpreds1 <- ggpredict(jaccard_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds1$year <- sort(unique(turnover_data$year))  

preds_glmm <- data.frame(  scaled_year= ggpreds1$x,
                           year = ggpreds1$year,
                           fit  =  ggpreds1$predicted,
                           lwr   = ggpreds1$conf.low,
                           upr   = ggpreds1$conf.high   )


ggpreds2 <- ggpredict(jaccard_gamm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")
ggpreds2$year <- sort(unique(turnover_data$year)  )

preds_gamm <- data.frame(  scaled_year= ggpreds2$x,
                           year = ggpreds2$year,
                           fit  =  ggpreds2$predicted,
                           lwr   = ggpreds2$conf.low,
                           upr   = ggpreds2$conf.high   )


ggpreds1 <- ggpredict(jaccard_spei_glmm, terms= c("scaled_year [all]", "spei72_2 [all]"),ci.lvl = .95, type="fixed")

ggpreds1$year <- turnover_data$year[match(ggpreds1$x, round(turnover_data$scaled_year, digits=2))]
ggpreds1$spei <- turnover_data$spei72_2[match(ggpreds1$x, round(turnover_data$scaled_year, digits=2))]
ggpreds1 <- ggpreds1[ggpreds1$group == round(ggpreds1$spei, digits=3),]  

preds_glmm1 <- data.frame(  scaled_year= ggpreds1$x,
                            year = ggpreds1$year,
                            fit  =  ggpreds1$predicted,
                            lwr   = ggpreds1$conf.low,
                            upr   = ggpreds1$conf.high)

# predictions for background decline
ggpreds2 <- ggpredict(jaccard_spei_glmm, terms= c("scaled_year [all]"),ci.lvl = .95, type="fixed")

ggpreds2$year <- turnover_data$year2[match(ggpreds2$x, round(turnover_data$scaled_year, digits=2))]

preds_glmm2 <- data.frame(  scaled_year= ggpreds2$x,
                            year = ggpreds2$year,
                            fit  =  ggpreds2$predicted,
                            lwr   = ggpreds2$conf.low)





plot10 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "Turnover (observed)", x = "")+
  ggtitle("H")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c(  "2001","", "", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(LTdiversity$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, jaccard),alpha=0.3, color = colors[4], size=2)+
  
  geom_ribbon(data = preds_glmm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[4])+
  
  geom_line(data=preds_glmm ,aes(scaled_year, fit), color = colors[4], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")


plot11 <- ggplot(preds_glmm, aes(scaled_year, fit)) + 
  theme_classic() +
  theme(text = element_text(size = 20)) +
  
  labs(y = "", x = "")+
  ggtitle("I")+
  theme(plot.title = element_text(hjust =0.9,vjust=-4 ))+
  #scale_y_continuous(limits = c(4, 25), breaks = c(5,10, 15, 20))+
  
  scale_x_continuous(labels= c(  "2001", "","", "", "", "", "", "", "", "2010", "", "", "", "", "", "", "", "", "", "2020","",""), breaks = c(unique(LTdiversity$scaled_year)) )+
  
  geom_jitter(data=LTdiversity ,aes(scaled_year, jaccard),alpha=0.3, color = colors[4], size=2)+
  
  geom_ribbon(data = preds_gamm, aes(y=fit, ymin=lwr, ymax=upr), alpha=0.5,fill= colors[2])+
  
  geom_line(data=preds_gamm ,aes(scaled_year, fit), color = colors[2], lty=1, lwd=2.5, lineend="round")+
  
  theme(legend.position="none")




plot12 <- ggplot() +
  annotate("text", x = 10,  y = 10,
           size = 6,
           label = "Drought index not \n a meaningful predictor") + theme_void()






grid.arrange(plot1,plot2,plot3,
             plot4,plot5,plot6,
             plot7, plot8, plot9,
             plot10, plot11, plot12, ncol=3)
















