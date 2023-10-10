# SUPPLEMENTARY R-SCRIPT TO:
# Long-term drought triggers severe declines in carabid beetles in a temperate forest
# F.Weiss, H.von Wehrden & A.Linde
# F.Weiss: ORCID 0000-0003-1078-1528

# PART IV: Modelling species trend and trait analysis

# Data available at: https://doi.org/10.48548/pubdata-46


#### Package list ####

library(data.table)
library(dplyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(cowplot)
library(cols4all)
library(bbmle)
library(scales)
library(foreign)
library(ggplot2)
library(MASS)
library(gridExtra)
library(ggthemes)




#### Loading data ####

# loading processed data from previous script
beetle_samples <- read.csv("carabid_samples_raw2.csv")

# loading sampling meta data
sampling_meta <- read.csv2("sampling_meta_dwd_spei.csv")


#### Re-organising data ####

# remove plots with sampling effort != 4
beetle_samples$sampling_length <- sampling_meta$Trapping_length[match(beetle_samples$sample_id, sampling_meta$sample_id)]
beetle_samples$sampling_effort <- sampling_meta$sampling_effort[match(beetle_samples$sample_id, sampling_meta$sample_id)]

# remove all samples of intervals with more or less than 4 traps
beetle_samples <- beetle_samples[beetle_samples$sampling_effort == 4,]


# check which traps have full 3 intervals!
# update sample_id2
beetle_samples$sample_id2 <- paste(beetle_samples$year, beetle_samples$plot, sep = "_")

data_table <- data.table(beetle_samples)   
consistency <- aggregate(cbind(beetle_samples$interval), by= list(beetle_samples$sample_id2), FUN=  function(x){length(unique(x))}) 
consistency2 <- consistency[consistency$V1 < 3,]
beetle_samples2 <- beetle_samples[!(beetle_samples$sample_id2 %in% c(consistency2$Group.1)), ]


## add rain and temp
beetle_samples2$temp <- sampling_meta$sampling_temp[match(beetle_samples2$sample_id, sampling_meta$sample_id)]
beetle_samples2$rain <- sampling_meta$sampling_rain[match(beetle_samples2$sample_id, sampling_meta$sample_id)]

sampling_meta2 <- unique(sampling_meta[,c(2,13:30)])


# aggregate per traps over all 3 intervals for all 4 traps
beetle_samples2$sample_id2 <- paste(beetle_samples2$year, beetle_samples2$plot, sep = "_")


# aggregate by sampleID (year and plot)
beetle_samples3 <- aggregate(cbind(beetle_samples2$sampling_abundance, beetle_samples2$corr_abundance2), by= list(beetle_samples2$sample_id2, beetle_samples2$year, beetle_samples2$plot, beetle_samples2$site ,beetle_samples2$species), FUN=sum  )

beetle_samples3B <- aggregate(cbind(beetle_samples2$temp, beetle_samples2$rain), by= list(beetle_samples2$sample_id2, beetle_samples2$year, beetle_samples2$plot, beetle_samples2$site ,beetle_samples2$species), FUN=mean  )

beetle_samples3<-beetle_samples3 %>% 
  rename(
    sample_id2 = Group.1,
    year = Group.2,
    plot= Group.3,
    site=Group.4,
    species=Group.5,
    sampling_abundance = V1,
    corr_abundance =V2
    
  )

beetle_samples3$temp <- beetle_samples3B$V1[match(beetle_samples3$sample_id2, beetle_samples3B$Group.1)]
beetle_samples3$rain <- beetle_samples3B$V2[match(beetle_samples3$sample_id2, beetle_samples3B$Group.1)]

# scale variables
beetle_samples3$scaled_year <- as.numeric(scale(beetle_samples3$year))
beetle_samples3$scaled_rain <- as.numeric(scale(beetle_samples3$rain))
beetle_samples3$scaled_temp <- as.numeric(scale(beetle_samples3$temp))


# combine species
beetle_samples3$species2 <- beetle_samples3$species
beetle_samples3$species[beetle_samples3$species == "Pterostichus nigrita"] <- "Pterostichus nigrita/rhaeticus"
beetle_samples3$species[beetle_samples3$species == "Pterostichus rhaeticus"] <- "Pterostichus nigrita/rhaeticus"
beetle_samples3$species[beetle_samples3$species == "Pterostichus diligens"] <- "Pterostichus strenuus/diligens"
beetle_samples3$species[beetle_samples3$species == "Pterostichus strenuus"] <- "Pterostichus strenuus/diligens"
beetle_samples3$species[beetle_samples3$species == "Pterostichus oblongopunctatus"] <- "Pterostichus obl./quadrif."
beetle_samples3$species[beetle_samples3$species == "Pterostichus quadrifoveolatus"] <- "Pterostichus obl./quadrif."


# remove zeros
beetle_samples4 <- beetle_samples3[beetle_samples3$species != "no carabids",]

zeros <- beetle_samples3[beetle_samples3$species == "no carabids",]


#### Species trends #####

species <- c(unique(beetle_samples4$species))

species_trends <- data.frame(
  species = species, 
  sum =  as.numeric(c(rep(NA,82))),
  corr_sum = as.numeric(c(rep(NA,82))),
  years =   as.numeric(c(rep(NA,82))),
  trend = as.numeric(c(rep(NA,82))),
  trend_AIC = as.numeric(c(rep(NA,82))),
  spei_eff = as.numeric(c(rep(NA,82))),
  spei_AIC = as.numeric(c(rep(NA,82))),
  model_family = as.character(c(rep(NA,82))),
  model_zi = as.character(c(rep(NA,82))),
  ar_test = as.character(c(rep("not_ok",82))),
  ar1 = c(rep(0,82))
)

# loop to determine number per species and number of years with data.

for(i in species){
  loopdata <- beetle_samples4[beetle_samples4$species == i,]
  species_trends$sum[species_trends$species == i] <- sum(loopdata$sampling_abundance)
  species_trends$corr_sum[species_trends$species == i] <- sum(loopdata$corr_abundance)
  species_trends$years[species_trends$species == i] <- length(unique(loopdata$year)) 
}

# only set threshhold
species_trends2 <- species_trends[species_trends$years > 3,]
species_trends2 <- species_trends2[species_trends2$species != "Bembidion spec.",]

species2 <- c(unique(species_trends2$species))

colors<- c4a("berlin", 7)
cblue <- colors[2]
corange <- colors[6]

rm(plot_lst)
plot_lst <- vector("list")

# modelselection between poisson and nbinom2, and zi ~0 and ~1
# https://stats.stackexchange.com/questions/571227/inflated-or-not-inflated-true-zero-dilemma-in-glmm

for(i in species2) {
  
  loopdata1 <- beetle_samples4[beetle_samples4$species == i,]
  
  loopdata2 <- rbind(zeros, loopdata1)
  
  loopdata3 <- aggregate(cbind(loopdata2$sampling_abundance), by= list(loopdata2$sample_id2, loopdata2$year, loopdata2$plot, loopdata2$site, loopdata2$scaled_rain, loopdata2$scaled_temp, loopdata2$scaled_year), FUN=sum  )
  
  loopdata3<-loopdata3 %>% 
    rename(
      sample_id2 = Group.1,
      year = Group.2,
      plot = Group.3,
      site = Group.4,
      scaled_rain = Group.5,
      scaled_temp = Group.6,
      scaled_year = Group.7,
      sampling_abundance=V1
    )
  
  loopdata3<- left_join(loopdata3, sampling_meta2, by= "year")
  
  loopdata3$year2 <- as.factor(as.character(loopdata3$year))
  loopdata3$plot <- as.factor(as.character(loopdata3$plot))
  loopdata3$site <- as.factor(as.character(loopdata3$site))
  
  plots <- unique(loopdata1$plot)
  loopdata3 <- loopdata3[loopdata3$plot %in% plots,]
  
  try({
    loopmodel1 <- glmmTMB(sampling_abundance ~ scaled_year + scaled_temp + scaled_rain + (1|year2) + (1|site/plot), data= loopdata3, family=poisson)
    loopmodel2 <- loopmodel1
    loopmodel3 <- loopmodel1
    loopmodel4 <- loopmodel1
  })
  try({
    loopmodel2 <- update(loopmodel1, ziformula= ~1)
  })
  try({
    loopmodel3 <- update(loopmodel1, family= nbinom2)
  })
  try({
    loopmodel4 <- update(loopmodel1, ziformula= ~1, family= nbinom2)
  })
  try({
    modList <- list( loopmodel1=loopmodel1, 
                     loopmodel2=loopmodel2,
                     loopmodel3=loopmodel3,
                     loopmodel4=loopmodel4)
    
    loopmodel <- modList[[which.min(sapply(modList, AICc))]]
    
    species_trends2$model_family[species_trends2$species == i] <- family(loopmodel)
    species_trends2$model_zi[species_trends2$species == i] <- paste( format(loopmodel$call$ziformula))
  })
  try({ 
    dharma_sim1 <- simulateResiduals(fittedModel = loopmodel, re.form= NULL)
    dharma_sim2 = recalculateResiduals(dharma_sim1, group = loopdata3$scaled_year)
    test_cor <- testTemporalAutocorrelation(dharma_sim2, time = unique(loopdata3$scaled_year), plot=F)
    species_trends2$ar_test[species_trends2$species == i] <- "ok"
  })
  
  try({
    if(test_cor$p.value < 0.05){
      loopmodel <- update(loopmodel, . ~ . + ar1(year2 +0|plot))
      species_trends2$ar1[species_trends2$species == i] <- 1
    }
  })
  
  try({
    
    coefs <- coef(summary(loopmodel))
    estimates1 <- coefs$cond
    species_trends2$trend[species_trends2$species == i] <- estimates1[2,1]
    
    loopmodel_null <- update(loopmodel, . ~ . - scaled_year )
    species_trends2$trend_AIC[species_trends2$species == i] <- AICc(loopmodel)- AICc(loopmodel_null)
  })
  
  
  try({
    
    m19 <- update(loopmodel_null,. ~ . + spei72)
    m20 <- update(loopmodel_null,. ~ . + spei72_1)
    m21 <- update(loopmodel_null,. ~ . + spei72_2)
    
    modList <- list( 
      m19=m19,
      m20=m20,
      m21=m21
    )
    
    loopmodel_spei <- modList[[which.min(sapply(modList, AICc))]]
    
    #loopmodel_spei <- update(loopmodel_null, . ~ . + spei60_2)
    
    coefs <- coef(summary(loopmodel_spei))
    estimates2 <- coefs$cond
    species_trends2$spei_eff[species_trends2$species == i] <- estimates2[4,1]
    
    species_trends2$spei_AIC[species_trends2$species == i] <- AICc(loopmodel_spei)- AICc(loopmodel_null) 
    
    pars <- rownames(estimates2)
    species_trends2$spei_type[species_trends2$species == i]<- pars[4]
  })
  
  
  
  try({
    preds <- ggpredict(loopmodel, terms= c("scaled_year [all]"),ci.lvl = 0.95, type="fixed", condition = c(scaled_temp = -0.02128506, scaled_rain = 0.02624506))
    preds$year <- sort(unique(loopdata3$year))
  })
  
  plot <- ggplot(loopdata3, aes(year, sampling_abundance)) + 
    theme_classic() +
    theme(text = element_text(size = 10)) +
    labs(y = "Sampling abundance", x = "Year")+
    ggtitle(i)+
    geom_jitter( color = "black", alpha=0.5,  size=2, height = 0.1)+
    scale_y_continuous(limits = c(0, max(loopdata3$sampling_abundance)), oob = scales::squish, labels = c(0,1,2,5,10,25,50,100,200,300,400), breaks = c(0,1,2,5,10,25,50,100,200,300,400))
  theme(legend.position="none")
  
  try({
    if( species_trends2$trend_AIC[species_trends2$species == i] < -2 &  estimates1[2,1] <0){ 
      plot <- plot + 
        geom_ribbon(data = preds, aes(x=year , y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.3,fill=corange)+
        geom_line(data=preds ,aes(year, predicted), color = corange, lty=1, lwd=2)
    }
    
    if(species_trends2$trend_AIC[species_trends2$species == i] < -2 &  estimates1[2,1] >0){       
      plot <- plot +
        geom_ribbon(data = preds, aes(x=year , y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.3,fill=cblue)+
        geom_line(data=preds ,aes(year, predicted), color = cblue, lty=1, lwd=2)
    }
    
    if(species_trends2$trend_AIC[species_trends2$species == i] > -2  &  estimates1[2,1] <0){       
      plot <- plot+
        #geom_ribbon(data = preds, aes(x=year , y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2,fill=colors[4])+
        geom_line(data=preds ,aes(year, predicted), color = colors[4], lty=1, lwd=2)
    }
    
    if(species_trends2$trend_AIC[species_trends2$species == i] > -2  &  estimates1[2,1] >0){       
      plot <- plot+
        #geom_ribbon(data = preds, aes(x=year , y=predicted, ymin=conf.low, ymax=conf.high), alpha=0.2,fill=colors[4])+
        geom_line(data=preds ,aes(year, predicted), color =  colors[4], lty=1, lwd=2)}         
  })
  
  plot_lst[[which(species2 == i)]] <- plot
  
  try({
    print(plot)
  })
  
  rm(loopmodel)
  rm(loopmodel1)
  rm(loopmodel2)
  rm(loopmodel3)
  rm(loopmodel4)
  rm(loopmodel_null)
  rm(loopmodel_spei)
  rm(modList)
  rm(coefs)
  rm(estimates2)
  rm(estimates1)
  rm(preds)
  rm( m19, m29, m21)
  #rm(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19, m29, m21)            
  print( species2 %>% { which(. == i) } )
  print(i)
  
}

# run manually for Stomis pumicatus ( AR1 borderline significant but ar1 null model does not converge, repeat without ar1 manually)

cowplot::plot_grid(plotlist = plot_lst, nrow = 8)

#### Re-organising data ####

# trend/effect classification 
species_trends2$trend2[species_trends2$trend <0 & species_trends2$trend_AIC < -2] <- "declining"
species_trends2$trend2[species_trends2$trend >0 & species_trends2$trend_AIC < -2] <- "increasing"
species_trends2$trend2[species_trends2$trend_AIC > -2] <- "uncertain"

species_trends2$spei_eff2[species_trends2$spei_eff <0 & species_trends2$spei_AIC < -2] <- "negative"
species_trends2$spei_eff2[species_trends2$spei_eff >0 & species_trends2$spei_AIC < -2] <- "positive"
species_trends2$spei_eff2[species_trends2$spei_AIC > -2] <- "uncertain"

# Traits were added based on sources (Supporting Information) outside R
species_trends <- read.csv2("published data/EWcarabids1999-2022_species_trends_traits.csv")

# order factors
species_trends$trend2 <- ordered(species_trends$trend2 , levels = c("declining", "uncertain", "increasing"))
species_trends$spei_eff2 <- ordered(species_trends$spei_eff2 , levels = c("positive", "uncertain", "negative"))
species_trends$wings <- ordered(species_trends$wings , levels = c("winged", "dimorphic", "shortwinged"))





#### Figure 5 ####

# subfigures were created in R and then compiled outside R


colors <- c4a("berlin", 7)

colors[8] <-  adjustcolor( colors[4], alpha.f = 0.5)


### traits ####

## size #####
plot1 <- ggplot(species_trends, aes(x = trend2, y = size, fill = trend2)) +
  theme_dark()+
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(trim=T, scale="width") +
  #scale_x_discrete(labels= c("decline","no trend","increase"), position = "top") +
  scale_x_discrete(labels= c("","","")) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[6],colors[8],colors[2]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank()) 

plot2 <- ggplot(species_trends, aes(x = spei_eff2, y = size, fill = spei_eff2)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(trim=T, scale="width") +
  #scale_x_discrete(labels= c("decline","no trend","increase"), position = "top") +
  scale_x_discrete(labels= c("","","")) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[2], colors[6],colors[8]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank()) 

dat <- species_trends[species_trends$spei_eff2 == "positive",]
dat[nrow(dat)+1,] <- NA
dat[11,]$spei_type <- "spei72_1"

plot3 <- ggplot(dat, aes(x = spei_type, y = size, fill = spei_type)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(trim=T, scale="width") +
  #scale_x_discrete(labels= c("decline","no trend","increase"), position = "top") +
  scale_x_discrete(labels= c("","","")) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[8],colors[2],colors[6]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())


grid.arrange(plot1, plot2, plot3, widths=c(1.2,1.2,0.8), nrow=1)
# export PDF 5x10




## wings #####
dat1 <- data.frame(table(species_trends$trend2, species_trends$wings))

plot1 <- ggplot(dat1,
                aes(x = Var1, 
                    y = Var2,
                    fill = Var1,
                    size = Freq)) +
  theme_dark()+
  theme(text = element_text(size = 30)) +
  geom_point(pch=21, color="black") +
  geom_text(aes(label = Freq), 
            colour = "black", 
            size = 7) +
  ggtitle("")+
  scale_x_discrete(labels= c("","","","")) +
  scale_size_continuous(range = c(7, 30)) + # Adjust as required.
  scale_fill_manual(values = c(colors[6],colors[8],colors[2])) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()) 

dat2 <- data.frame(table(species_trends$spei_eff2, species_trends$wings))

plot2 <- ggplot(dat2,
                aes(x = Var1, 
                    y = Var2,
                    fill = Var1,
                    size = Freq)) +
  theme_dark()+
  theme(text = element_text(size = 30)) +
  geom_point(pch=21, color="black") +
  geom_text(aes(label = Freq), 
            colour = "black", 
            size = 7) +
  ggtitle("")+
  scale_x_discrete(labels= c("","","","")) +
  #scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(7, 30)) + # Adjust as required.
  scale_fill_manual(values = c(colors[6],colors[8],colors[2])) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks = element_blank()) 

dat3 <- data.frame(table(species_trends$spei_type, species_trends$wings, species_trends$spei_eff2) )
dat3 <- dat3[dat3$Var3 == "positive",]
dat3 <- dat3[dat3$Var1 %in% c("spei72", "spei72_2"),]


plot3 <- ggplot(dat3,
                aes(x = Var1, 
                    y = Var2,
                    fill = Var1,
                    size = Freq)) +
  theme_dark()+
  theme(text = element_text(size = 30)) +
  geom_point(pch=21, color="black") +
  geom_text(aes(label = Freq), 
            colour = "black", 
            size = 7) +
  ggtitle("")+
  scale_x_discrete(labels= c("","","","")) +
  #scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(7, 30)) + # Adjust as required.
  scale_fill_manual(values = c(colors[8],colors[6])) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks = element_blank()) 


grid.arrange(plot1, plot2, plot3, widths=c(1.2,1,0.8), nrow=1)
# export PDF 5x10





## feeding guild #####

table(species_trends$feeding)

dat1 <- data.frame(table(species_trends$trend2, species_trends$feeding))

plot1 <- ggplot(dat1,
                aes(x = Var1, 
                    y = Var2,
                    fill = Var1,
                    size = Freq)) +
  theme_dark()+
  theme(text = element_text(size = 30)) +
  geom_point(pch=21, color="black") +
  geom_text(aes(label = Freq), 
            colour = "black", 
            size = 7) +
  ggtitle("")+
  scale_x_discrete(labels= c("","","","")) +
  scale_size_continuous(range = c(7, 30)) + # Adjust as required.
  scale_fill_manual(values = c(colors[6],colors[8],colors[2])) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks = element_blank()) 

dat2 <- data.frame(table(species_trends$spei_eff2, species_trends$feeding))

plot2 <- ggplot(dat2,
                aes(x = Var1, 
                    y = Var2,
                    fill = Var1,
                    size = Freq)) +
  theme_dark()+
  theme(text = element_text(size = 30)) +
  geom_point(pch=21, color="black") +
  geom_text(aes(label = Freq), 
            colour = "black", 
            size = 7) +
  ggtitle("")+
  scale_x_discrete(labels= c("","","","")) +
  #scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(7, 30)) + # Adjust as required.
  scale_fill_manual(values = c(colors[6],colors[8],colors[2])) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks = element_blank()) 


dat3 <- data.frame(table(species_trends$spei_type, species_trends$feeding, species_trends$spei_eff2) )
dat3 <- dat3[dat3$Var3 == "positive",]
dat3 <- dat3[dat3$Var1 %in% c("spei72", "spei72_2"),]


plot3 <- ggplot(dat3,
                aes(x = Var1, 
                    y = Var2,
                    fill = Var1,
                    size = Freq)) +
  theme_dark()+
  theme(text = element_text(size = 30)) +
  geom_point(pch=21, color="black") +
  geom_text(aes(label = Freq), 
            colour = "black", 
            size = 7) +
  ggtitle("")+
  scale_x_discrete(labels= c("","","","")) +
  #scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(7, 30)) + # Adjust as required.
  scale_fill_manual(values = c(colors[8],colors[6])) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y =element_blank(),
        axis.ticks = element_blank()) 



grid.arrange(plot1, plot2, plot3, widths=c(1,1,0.9), nrow=1)
# export PDF 5x10



## sustek humidity ####

species_trends$humidity_preference_SUSTEK <- as.numeric(species_trends$humidity_preference_SUSTEK)

plot1<- ggplot(species_trends, aes(x = trend2, y = humidity_preference_SUSTEK, fill = trend2)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin() +
  scale_x_discrete(labels= c("declining","uncertain","increasing")) +
  scale_y_continuous(labels= c("dry","","","","","humid"), breaks =c(3:8)) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[6],colors[8],colors[2]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank()) 


species_trends$humidity_preference_SUSTEK <- as.numeric(species_trends$humidity_preference_SUSTEK)

plot2<- ggplot(species_trends, aes(x = spei_eff2, y = humidity_preference_SUSTEK, fill = spei_eff2)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin() +
  scale_x_discrete(labels= c("positive","uncertain","negative")) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[2], colors[6],colors[8]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())


dat <- species_trends[species_trends$spei_eff2 == "positive",]
dat[nrow(dat)+1,] <- NA
dat[11,]$spei_type <- "spei72_1"

plot3 <- ggplot(dat, aes(x = spei_type, y = humidity_preference_SUSTEK, fill = spei_type)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(trim=T, scale="width") +
  scale_x_discrete(labels= c("no lag","lag 1","lag 2"), position = "top") +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[8],colors[2],colors[6]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())


grid.arrange(plot1, plot2, plot3, widths=c(1.2,1.2,0.8), nrow=1)



# export PDF 5x10


## latitude ####

plot1 <-  ggplot(species_trends, aes(x = trend2, y = latitude, fill = trend2)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(scale = "width") +
  scale_x_discrete(labels= c("declining","uncertain","increasing")) +
  scale_y_continuous(limits= c(42.5,52.5), labels= c("south","","","","north"), breaks = c(42.5,45, 47.5,50,52.5)) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[6],colors[8],colors[2]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank()) 

plot2 <- ggplot(species_trends, aes(x = spei_eff2, y = latitude, fill = spei_eff2)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(scale = "width") +
  scale_x_discrete(labels= c("declining","uncertain","increasing")) +
  scale_y_continuous(limits= c(42.5,52.5), labels= c("south","","","","north"), breaks = c(42.5,45, 47.5,50,52.5)) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[2], colors[6],colors[8]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank()) 

dat <- species_trends[species_trends$spei_eff2 == "positive",]
dat[nrow(dat)+1,] <- NA
dat[11,]$spei_type <- "spei72_1"

plot3 <- ggplot(dat, aes(x = spei_type, y = latitude, fill = spei_type)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(trim=T, scale="width") +
  scale_x_discrete(labels= c("no lag","lag 1","lag 2"), position = "top") +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[8],colors[2],colors[6]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_blank()) 

grid.arrange(plot1, plot2, plot3, widths=c(1.2,1.2,0.8), nrow=1)


## local abundance ####

plot1 <-  ggplot(species_trends, aes(x = trend2, y = sum, fill = trend2)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(scale = "width") +
  scale_x_discrete(labels= c("","","")) +
  scale_y_continuous(limits= c(0,6100), labels= c("low","high"), breaks = c(0,5800)) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[6],colors[8],colors[2]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        #axis.ticks.y = element_blank(),
        #axis.text.y = element_blank(),
        panel.background = element_blank()) 

plot2 <- ggplot(species_trends, aes(x = spei_eff2, y = sum, fill = spei_eff2)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(scale = "width") +
  scale_x_discrete(labels= c("","","")) +
  scale_y_continuous(limits= c(0,6100), labels= c("low","high"), breaks = c(0,5800)) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[2], colors[6],colors[8]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank()) 


dat <- species_trends[species_trends$spei_eff2 == "positive",]
dat[nrow(dat)+1,] <- NA
dat[11,]$spei_type <- "spei72_1"

plot3 <- ggplot(dat, aes(x = spei_type, y = sum, fill = spei_type)) +
  theme_dark() +
  theme(text = element_text(size = 30)) +
  labs(y = "", x = "")+
  ggtitle("")+
  geom_violin(trim=T, scale="width") +
  scale_x_discrete(labels= c("","","")) +
  geom_jitter(height=0.05, width=0.05, size=3) +
  scale_fill_manual(values = c(colors[8],colors[2],colors[6]))+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.background = element_blank()) 

grid.arrange(plot1, plot2, plot3, widths=c(1.2,1.2,0.8), nrow=1)