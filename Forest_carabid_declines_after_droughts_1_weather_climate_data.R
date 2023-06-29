
# SUPPLEMENTARY R-SCRIPT TO:
# Long-term data reveal: Recent declines in carabid beetles in a temperate forest are linked to severe droughts
# F.Weiss, H.von Wehrden & A.Linde
# F.Weiss: ORCID 0000-0003-1078-1528

# PART I: Compiling wheather and climate Data

# Data available at: 

#### Package list ####

library(lubridate)
library(dplyr)
library(SPEI)
library(cols4all)
library(ggplot2)

#### Loading meteorological data of two DWD stations, select relevant period and extract temperature and precipitation #####

# DWD data is available at: https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/
# Angermünde station ID 164
# Buch station ID 400

dwd_angerm <- read.csv("dwd_angerm.csv")
dwd_buch <- read.csv("dwd_buch.csv")

merge_angerm <- dwd_angerm[dwd_angerm$date >= "1961-01-01" & dwd_angerm$date < "2022-12-01",]
merge_buch <- dwd_buch[dwd_buch$date >= "1961-01-01" & dwd_buch$date < "2022-12-01",]

merge_angerm<- merge_angerm[order(as.Date(merge_angerm$date, format="%Y/%m/%d")),]
merge_buch<- merge_buch[order(as.Date(merge_buch$date, format="%Y/%m/%d")),]

merge_angerm %>% arrange(ymd(merge_angerm$date))
merge_buch %>% arrange(ymd(merge_buch$date))

dwd_merge <- data.frame(date = merge_angerm$date,
                        MESS_DATUM = merge_angerm$MESS_DATUM,
                        RSK = (merge_angerm$RSK + merge_buch$RSK)/2,
                        TMK = (merge_angerm$TMK + merge_buch$TMK)/2,
                        TXK = (merge_angerm$TXK + merge_buch$TXK)/2,
                        TNK = (merge_angerm$TNK + merge_buch$TNK)/2)


# Buch missing data, replace with Angermünde data only.

dwd_merge$RSK[dwd_merge$RSK < 0] <- merge_angerm$RSK[match(dwd_merge$MESS_DATUM[dwd_merge$RSK < 0], merge_angerm$MESS_DATUM)]

dwd_merge$TMK[dwd_merge$TMK < -200] <- merge_angerm$TMK[match(dwd_merge$MESS_DATUM[dwd_merge$TMK < -200], merge_angerm$MESS_DATUM)]

dwd_merge$TXK[dwd_merge$TXK < -200] <- merge_angerm$TXK[match(dwd_merge$MESS_DATUM[dwd_merge$TXK < -200], merge_angerm$MESS_DATUM)]

dwd_merge$TNK[dwd_merge$TNK < -200] <- merge_angerm$TNK[match(dwd_merge$MESS_DATUM[dwd_merge$TNK < -200], merge_angerm$MESS_DATUM)]

# load sampling metadata
sampling_meta <- read.csv("EWcarabids1999-2022_samplingscheme.csv")


# define date variables
sampling_meta$Trapping_Start <- as.Date(sampling_meta$Trapping_Start, format = "%Y-%m-%d", origin= "1970-01-01" )
sampling_meta$Trapping_End <- as.Date(sampling_meta$Trapping_End, format = "%Y-%m-%d", origin= "1970-01-01" )

dwd_merge$date <- as.Date(dwd_merge$date, format = "%Y-%m-%d", origin= "1970-01-01" )


# extract mean daily temperature and mean daily rainfall for sampling period

intervals <- sampling_meta$sample_id

for (i in intervals){
  
  weather_period <- as.Date(c(sampling_meta$Trapping_Start[sampling_meta$sample_id ==i] : sampling_meta$Trapping_End[sampling_meta$sample_id ==i]),format= "%Y-%m-%d", origin= "1970-01-01"  )
  
  sampling_meta$sampling_temp[sampling_meta$sample_id ==i] <- mean(dwd_merge$TMK[dwd_merge$date %in% weather_period])
  
  sampling_meta$sampling_rain[sampling_meta$sample_id ==i] <- mean(dwd_merge$RSK[dwd_merge$date %in% weather_period])
  
}

#### SPEI - drought index #####

dwd_merge$year_month <- substring(dwd_merge$date,1,7)

spei_dat1<- aggregate(cbind(dwd_merge$RSK), by= list(dwd_merge$year_month), FUN= sum)

spei_dat2<- aggregate(cbind(dwd_merge$TMK, dwd_merge$TXK, dwd_merge$TNK), by= list(dwd_merge$year_month), FUN= mean)

spei_dat<- spei_dat2 %>% 
  rename(
    year_month = Group.1,
    TMK = V1,
    TXK = V2,
    TNK = V3
  )

spei_dat$RSK <- spei_dat1$V1[match(spei_dat$year_month, spei_dat1$Group.1)]

# potential-evapotransperation sensu Thornthwaite (1948)
# lat = N 52.822110 (decimal degree Eberswalde)
spei_dat$PE <- thornthwaite(Tave = spei_dat$TMK, lat= 52.822110)

spei_dat$CWB <- spei_dat$RSK - spei_dat$PE

spei <- spei(spei_dat$CWB, scale = 12)
spei_dat$spei12 <- spei$fitted

spei <- spei(spei_dat$CWB, scale = 24)
spei_dat$spei24 <- spei$fitted

spei <- spei(spei_dat$CWB, scale = 36)
spei_dat$spei36 <- spei$fitted

spei <- spei(spei_dat$CWB, scale = 48)
spei_dat$spei48 <- spei$fitted

spei <- spei(spei_dat$CWB, scale = 60)
spei_dat$spei60 <- spei$fitted

spei <- spei(spei_dat$CWB, scale = 72)
spei_dat$spei72 <- spei$fitted

spei_dat$year <- as.numeric(substring(spei_dat$year_month, 1,4))

par(mfrow=c(1,1))


# Compile different SPEI-variables
spei_dat$year <- as.numeric(substring(spei_dat$year_month,1,4))
spei_dat$month <- as.numeric(substring(spei_dat$year_month,6,7))
str(spei_dat)


# SPEI lag0

spei_season <- spei_dat[spei_dat$month %in% c(3:7),]

spei_season <- aggregate(cbind( spei_season$spei12, spei_season$spei24, spei_season$spei36, spei_season$spei48, spei_season$spei60, spei_season$spei72), by=list(spei_season$year), FUN= mean)
spei_season<-spei_season %>% 
  rename(
    year = Group.1,
    spei12 = V1,
    spei24 = V2,
    spei36 = V3,
    spei48 = V4,
    spei60 = V5,
    spei72 = V6
)

# SPEI lag1
spei_dat$gen1[spei_dat$month %in% c(1:2)] <- spei_dat$year[spei_dat$month %in% c(1:2)]
spei_dat$gen1[spei_dat$month %in% c(3:12)] <- spei_dat$year[spei_dat$month %in% c(3:12)]+1

spei_dat$gen2 <- spei_dat$gen1+1


spei_year1 <- aggregate(cbind( spei_dat$spei12, spei_dat$spei24, spei_dat$spei36, spei_dat$spei48, spei_dat$spei60, spei_dat$spei72), by=list(spei_dat$gen1), FUN= mean)
spei_year1<-spei_year1 %>% 
  rename(
    year = Group.1,
    spei12_1 = V1,
    spei24_1 = V2,
    spei36_1 = V3,
    spei48_1 = V4,
    spei60_1 = V5,
    spei72_1 = V6
  )

# SPEI lag2
spei_year2 <- aggregate(cbind( spei_dat$spei12, spei_dat$spei24, spei_dat$spei36, spei_dat$spei48, spei_dat$spei60, spei_dat$spei72), by=list(spei_dat$gen2), FUN= mean)
spei_year2<-spei_year2 %>% 
  rename(
    year = Group.1,
    spei12_2 = V1,
    spei24_2 = V2,
    spei36_2 = V3,
    spei48_2 = V4,
    spei60_2 = V5,
    spei72_2 = V6
  )

spei_full <- merge(spei_season, spei_year1, by.x = "year", by.y ="year")
spei_full <- merge(spei_full, spei_year2, by.x = "year", by.y ="year")

# scale SPEI data
spei_scaled <- data.frame(scale(spei_full))
spei_scaled$year <- spei_full$year


sampling_meta<- merge(sampling_meta, spei_scaled, by.x = "year", by.y ="year")



# use later
write.csv2(sampling_meta,'sampling_meta_dwd_spei.csv')


#### Figures ####

# Figure 6

spei_dat$year_dec <- spei_dat$year + 1/12 * spei_dat$month

# take preds_gamm of the biomass GAMM (related Rcode file #2)

preds_gamm$scaled_fit <- c(scale(preds_gamm$fit))

preds_gamm$year_dec<- as.numeric(c(as.character(preds_gamm$year)))

spei_dat$color[spei_dat$spei72 <= 0] <- "2"
spei_dat$color[spei_dat$spei72 > 0] <- "1"

colors <- c4a("berlin", 7)

ggplot(spei_dat, aes(year_dec, spei72, color = color, group=1)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI", x = "Year")+
  ggtitle("")+
  theme(plot.title = element_text(hjust =0.05,vjust=-4))+
  scale_x_continuous(limits= c(1967,2022),labels= c(1970, 1980,1990,2000,2010,2020), breaks = c(1970, 1980,1990,2000,2010,2020)) +
  scale_color_manual(values = c(colors[2],colors[6]))+
  geom_line(data=spei_dat ,aes(year_dec, spei72), lty=1, lwd=1.5)+
  #geom_line(data=preds_gamm ,aes(year_dec, scaled_fit), color = colors[4], lty=1, lwd=3)+
  theme(legend.position="none")


# plots for Supporting Information

plot1 <- ggplot(spei_dat, aes(year_dec, spei12)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI", x = "Year")+
  ggtitle("SPEI 12")+
  scale_x_continuous(limits= c(1999,2022),labels= c(1970, 1980,1990,2000,2010,2020), breaks = c(1970, 1980,1990,2000,2010,2020)) +
  geom_line(data=spei_dat ,aes(year_dec, spei12), lty=1, lwd=3, color= colors[4], alpha =0.5)+
  geom_line(data=spei_full ,aes(year, spei12), lty=1, lwd=2, color= colors[4])+
  geom_line(data=spei_full ,aes(year, spei12_1), lty=1, lwd=2, color= colors[2])+
  geom_line(data=spei_full ,aes(year, spei12_2), lty=1, lwd=2, color= colors[6])+
  theme(legend.position="none")

plot2 <- ggplot(spei_dat, aes(year_dec, spei24)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI", x = "Year")+
  ggtitle("SPEI 24")+
  scale_x_continuous(limits= c(1999,2022),labels= c(1970, 1980,1990,2000,2010,2020), breaks = c(1970, 1980,1990,2000,2010,2020)) +
  geom_line(data=spei_dat ,aes(year_dec, spei24), lty=1, lwd=3, color= colors[4], alpha =0.5)+
  geom_line(data=spei_full ,aes(year, spei24), lty=1, lwd=2, color= colors[4])+
  geom_line(data=spei_full ,aes(year, spei24_1), lty=1, lwd=2, color= colors[2])+
  geom_line(data=spei_full ,aes(year, spei24_2), lty=1, lwd=2, color= colors[6])+
  theme(legend.position="none")

plot3 <- ggplot(spei_dat, aes(year_dec, spei36)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI", x = "Year")+
  ggtitle("SPEI 36")+
  scale_x_continuous(limits= c(1999,2022),labels= c(1970, 1980,1990,2000,2010,2020), breaks = c(1970, 1980,1990,2000,2010,2020)) +
  geom_line(data=spei_dat ,aes(year_dec, spei36), lty=1, lwd=3, color= colors[4], alpha =0.5)+
  geom_line(data=spei_full ,aes(year, spei36), lty=1, lwd=2, color= colors[4])+
  geom_line(data=spei_full ,aes(year, spei36_1), lty=1, lwd=2, color= colors[2])+
  geom_line(data=spei_full ,aes(year, spei36_2), lty=1, lwd=2, color= colors[6])+
  theme(legend.position="none")

plot4 <- ggplot(spei_dat, aes(year_dec, spei48)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI", x = "Year")+
  ggtitle("SPEI 48")+
  scale_x_continuous(limits= c(1999,2022),labels= c(1970, 1980,1990,2000,2010,2020), breaks = c(1970, 1980,1990,2000,2010,2020)) +
  geom_line(data=spei_dat ,aes(year_dec, spei48), lty=1, lwd=3, color= colors[4], alpha =0.5)+
  geom_line(data=spei_full ,aes(year, spei48), lty=1, lwd=2, color= colors[4])+
  geom_line(data=spei_full ,aes(year, spei48_1), lty=1, lwd=2, color= colors[2])+
  geom_line(data=spei_full ,aes(year, spei48_2), lty=1, lwd=2, color= colors[6])+
  theme(legend.position="none")

plot5 <- ggplot(spei_dat, aes(year_dec, spei60)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI", x = "Year")+
  ggtitle("SPEI 60")+
  scale_x_continuous(limits= c(1999,2022),labels= c(1970, 1980,1990,2000,2010,2020), breaks = c(1970, 1980,1990,2000,2010,2020)) +
  geom_line(data=spei_dat ,aes(year_dec, spei60), lty=1, lwd=3, color= colors[4], alpha =0.5)+
  geom_line(data=spei_full ,aes(year, spei60), lty=1, lwd=2, color= colors[4])+
  geom_line(data=spei_full ,aes(year, spei60_1), lty=1, lwd=2, color= colors[2])+
  geom_line(data=spei_full ,aes(year, spei60_2), lty=1, lwd=2, color= colors[6])+
  theme(legend.position="none")

plot6 <- ggplot(spei_dat, aes(year_dec, spei72)) + 
  theme_classic() +
  theme(text = element_text(size = 25)) +
  labs(y = "SPEI", x = "Year")+
  ggtitle("SPEI 72")+
  scale_x_continuous(limits= c(1999,2022),labels= c(1970, 1980,1990,2000,2010,2020), breaks = c(1970, 1980,1990,2000,2010,2020)) +
  geom_line(data=spei_dat ,aes(year_dec, spei72), lty=1, lwd=3, color= colors[4], alpha =0.5)+
  geom_line(data=spei_full ,aes(year, spei72), lty=1, lwd=2, color= colors[4])+
  geom_line(data=spei_full ,aes(year, spei72_1), lty=1, lwd=2, color= colors[2])+
  geom_line(data=spei_full ,aes(year, spei72_2), lty=1, lwd=2, color= colors[6])+
  theme(legend.position="none")


grid.arrange(plot1, plot2, plot3, nrow=3)

grid.arrange( plot4, plot5, plot6, nrow=3)














