# Libraries
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(ggplot2)
library(dplyr)
library(forcats)
library(ggpubr)
library(plyr)
library(anytime)
library(fANCOVA)
library(tweedie)
library(car)
library(locfit)
library(MuMIn)
library(tidyverse)
library(mgcv)
library(ggpubr)
library(mgcViz)
library(cplm)
library(statmod)
library(gee)
library(geepack)
library(TSA)
library(epitools)
library(lubridate)
library(survival)

#load data
site = 'GI'
saveDir = paste("G:/Baja/Seasonality/")
filename = paste(saveDir,site,"_dayData_forGLMR125.csv",sep="")
dayBinTAB = read.csv(filename) #no effort days deleted
head(dayBinTAB)
str(dayBinTAB)
dayBinTAB$Season = as.factor(dayBinTAB$Season) #change season from an integer to a factor
levels(dayBinTAB$Season)
dayBinTAB$Season = revalue(dayBinTAB$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
dayBinTAB$tbin = anytime(as.factor(dayBinTAB$tbin))

#plot data as proportion of hours per day with clicks
title1 = paste(site,"Proportion of Hours per Day with Clicks")
ggplot(dayBinTAB, aes(x=tbin,y=HoursProp))+
  ggtitle(title1)+
  labs(y="Proportion of Hours per Day with Clicks",x="Time (days)")+
  geom_line()+
  geom_point()
fig1 =paste(saveDir,site,"HoursProp_TimeSeries.png",sep="")
ggsave(fig1)

#plot data as box plot for seasons; have to plot this with no effort days deleted
title2 = paste("Seasonal Presence at",site)
ggplot(dayBinTAB, aes(x=Season, y=HoursProp, color = Season))+
  geom_boxplot()+
  ggtitle(title2)+
  labs(y="Proportion of Hours per Day with Clicks")
  scale_color_brewer(palette = "Dark2")
fig2 =paste(saveDir,site,"BoxPlot.png",sep="")
ggsave(fig2)

#groupin data according to ITS
if (site == 'PT'){
  n = 4
  GroupedDay = aggregate(dayBinTAB,list(rep(1:(nrow(dayBinTAB)%/%n+1),each=n,len=nrow(dayBinTAB))),mean)[-1];
}

if (site == 'BD'){
  n = 5
  GroupedDay = aggregate(dayBinTAB,list(rep(1:(nrow(dayBinTAB)%/%n+1),each=n,len=nrow(dayBinTAB))),mean)[-1];
}

if (site == 'QN'){
  n = 12
  GroupedDay = aggregate(dayBinTAB,list(rep(1:(nrow(dayBinTAB)%/%n+1),each=n,len=nrow(dayBinTAB))),mean)[-1];
}

if (site == 'CB'){
  n = 10
  GroupedDay = aggregate(dayBinTAB,list(rep(1:(nrow(dayBinTAB)%/%n+1),each=n,len=nrow(dayBinTAB))),mean)[-1];
}

if (site == 'AB'){
  n = 6
  GroupedDay = aggregate(dayBinTAB,list(rep(1:(nrow(dayBinTAB)%/%n+1),each=n,len=nrow(dayBinTAB))),mean)[-1];
}

if (site == 'KOA'){
  n = 3
  GroupedDay = aggregate(dayBinTAB,list(rep(1:(nrow(dayBinTAB)%/%n+1),each=n,len=nrow(dayBinTAB))),mean)[-1];
}

if (site == 'KS'){
  GroupedDay = dayBinTAB
}

if (site == 'GI'){
  GroupedDay = dayBinTAB
}

#round day, year, month, and find season for ITS data
if (site == 'KS' || site == "GI"){
}else{
GroupedDay$day = floor(GroupedDay$day)
GroupedDay$month = floor(GroupedDay$month)
GroupedDay$Year = floor(GroupedDay$Year)
GroupedDay$Season[GroupedDay$month == 1 | GroupedDay$month == 2 | GroupedDay$month == 3] = 1
GroupedDay$Season[GroupedDay$month == 4 | GroupedDay$month == 5 | GroupedDay$month == 6] = 2
GroupedDay$Season[GroupedDay$month == 7 | GroupedDay$month == 8 | GroupedDay$month == 9] = 3
GroupedDay$Season[GroupedDay$month == 10 | GroupedDay$month == 11 | GroupedDay$month == 12] = 4
GroupedDay$Season = as.factor(GroupedDay$Season) #change season from an integer to a factor
GroupedDay$Season = revalue(GroupedDay$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall")) #change the numbers in actual seasons
}

#plot data grouped with ITS as proportion of hours per day with clicks
title1 = paste(site,"Proportion of Hours per Day with Clicks")
ggplot(GroupedDay, aes(x=tbin,y=HoursProp))+
  ggtitle(title1)+
  labs(y="Proportion of Hours per Day with Clicks",x="Time (days)")+
  geom_line()+
  geom_point()
fig1 =paste(saveDir,site,"HoursProp_TimeSeriesITS.png",sep="")
ggsave(fig1)

##### grouped data by day of year - mean
filename2 = paste(saveDir,site,"_days365GroupedMean_forGLMR125.csv",sep="")
oneyear = read.csv(filename2) #bin means from days

if (nrow(oneyear) >= 365) {
names(oneyear) <- c("Day", "HoursProp", "SEM", "Std", "Variance", "Range", 'Season', "Month")
oneyear$Season = as.factor(oneyear$Season) #change season from an integer to a factor
levels(oneyear$Season)
oneyear$Season = revalue(oneyear$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
head(oneyear)
str(oneyear)
} else {
oneyear$Season = as.factor(oneyear$Season) #change season from an integer to a factor
levels(oneyear$Season)
oneyear$Season = revalue(oneyear$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
head(oneyear)
str(oneyear)
}

#groupin data according to ITS
if (site == 'PT'){
  n = 4
  GroupedYear = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'BD'){
  n = 5
  GroupedYear = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'QN'){
  n = 5
  GroupedYear = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'CB'){
  n = 21
  GroupedYear = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'AB'){
  n = 6
  GroupedYear = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'KOA'){
  n = 3
  GroupedYear = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'KS'){
  GroupedYear = oneyear
}

if (site == 'GI'){
  GroupedYear = oneyear
}

#round day, year, month, and find season for ITS data
if (site == 'KS' || site == 'GI'){
}else{
GroupedYear$Month = month(as.Date(GroupedYear$Day, origin = "2014-01-01"))
GroupedYear$Season[GroupedYear$Month == 1 | GroupedYear$Month == 2 | GroupedYear$Month == 3] = 1
GroupedYear$Season[GroupedYear$Month == 4 | GroupedYear$Month == 5 | GroupedYear$Month == 6] = 2
GroupedYear$Season[GroupedYear$Month == 7 | GroupedYear$Month == 8 | GroupedYear$Month == 9] = 3
GroupedYear$Season[GroupedYear$Month == 10 | GroupedYear$Month == 11 | GroupedYear$Month == 12] = 4
GroupedYear$Season = as.factor(GroupedYear$Season) #change season from an integer to a factor
GroupedYear$Season = revalue(GroupedYear$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall")) #change the numbers in actual seasons
}

#plot data as time series
title3 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks")
ggplot(oneyear, aes(x=Day,y=HoursProp))+
  ggtitle(title3)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_line()+
  geom_point()
fig4 =paste(saveDir,site,"AveragedHoursProp_TimeSeries.png",sep="")
ggsave(fig4)

#plot data as time series with ITS
title3 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks")
ggplot(GroupedYear, aes(x=Day,y=HoursProp))+
  ggtitle(title3)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_line()+
  geom_point()
fig4 =paste(saveDir,site,"AveragedHoursProp_TimeSeriesITS.png",sep="")
ggsave(fig4)

if (nrow(oneyear) >= 365) {
#plot data as time series with error bars
title4 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks")
ggplot(oneyear, aes(x=Day,y=HoursProp))+
  ggtitle(title4)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_errorbar(aes(ymin = HoursProp - SEM, ymax = HoursProp + SEM))+
  geom_line()+
  geom_point()
fig5 =paste(saveDir,site,"AveragedHoursProp_TimeSeries_ErrorBars.png",sep="")
ggsave(fig5)
}else{}

## GAMs ##

#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
if (site == 'AB'){
  gamTw = gam(HoursProp ~ s(day, bs = 'cc', k = 19), data = GroupedDay, family = tw, method = "REML")
  plot(gamTw, pages =1)
  summary(gamTw)
}else{
  gamTw = gam(HoursProp ~ s(day, bs = 'cc', k = 47), data = GroupedDay, family = tw, method = "REML")
  plot(gamTw, pages =1)
  summary(gamTw)
}

#GAM to check for significance between seasons
gamTwS = gam(HoursProp ~ Season, data = GroupedDay, family = tw, method = "REML")
summary(gamTwS)

#Better GAM plots
#pattern only
viz = getViz(gamTw)
print(plot(viz,allTerms=T),pages=1)

#first way to plot GAM
vizGG = plot(viz,allTerms = T) +
  labs(title = 'Sperm whales (GAM)')+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
print(vizGG,pages =1)
fig6 =paste(saveDir,site,"GAM1.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  labs(title = 'Sperm whales (GAM)')+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
fig7 =paste(saveDir,site,"GAM2.png",sep="")
ggsave(fig7)

###load all data and run GAM
filename = paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Data.csv",sep="")
AllTAB = read.csv(filename) #no effort days deleted
head(AllTAB)
str(AllTAB)
AllTAB$Season = as.factor(AllTAB$Season) #change season from an integer to a factor
levels(AllTAB$Season)
AllTAB$Season = revalue(AllTAB$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
AllTAB$tbin = anytime(as.factor(AllTAB$tbin))

#remove NaNs
AllTable = na.omit(AllTAB)
AllTable$Site = as.character(AllTable$Site)

#table with only CB, QN, PT
Central_GOA = subset(AllTable, Site!="BD" & Site!="KS" & Site!="KOA" & Site!="AB")

gamALL = gam(HoursProp ~ s(day, bs = 'cc', k = 47), data = AllTable, family = tw, method = "REML")
plot(gamALL, pages =1)
summary(gamALL)

#Better GAM plots
#pattern only
viz = getViz(gamALL)
print(plot(viz,allTerms=T),pages=1)

#first way to plot GAM
vizGG = plot(viz,allTerms = T) +
  labs(title = 'Sperm whales (GAM)')+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
print(vizGG,pages =1)
fig6 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/GAM1.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  labs(title = 'Sperm whales (GAM)')+
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
print(vizGG2,pages =1)
fig7 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/GAM2.png",sep="")
ggsave(fig7)

