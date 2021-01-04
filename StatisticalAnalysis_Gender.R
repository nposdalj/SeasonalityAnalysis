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
library(gtable)

#load data
site = 'BD'
filename = paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,site,"_binPresence.csv",sep="")
binPresence = read.csv(filename) #no effort days deleted
head(binPresence)
str(binPresence)
binPresence$Season = as.factor(binPresence$Season) #change season from an integer to a factor
levels(binPresence$Season)
binPresence$Season = revalue(binPresence$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
binPresence$tbin = anytime(as.factor(binPresence$tbin))

#plot data as proportion of hours per day with clicks
title1 = paste(site,"Proprtion of Hours/Day w/ Clicks")
plot1 = ggplot(binPresence, aes(x=tbin,y=FeHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot2 = ggplot(binPresence, aes(x=tbin,y=JuHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot3 = ggplot(binPresence, aes(x=tbin,y=MaHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
figure = ggarrange(plot1,plot2,plot3, labels = c("Social Units","  Mid-Size  ","    Males   "),align = "v",ncol = 1, nrow = 3)
annotate_figure(figure, top = text_grob(title1, face = "bold", size = 14), bottom = text_grob("Time (years)"),
                left = text_grob("Proportion of Hours/Day w/Clicks", rot = 90))
fig1 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"HoursProp_TimeSeries_StackedGroups.png",sep="")
ggsave(fig1)

#plot data as box plot for seasons; have to plot this with no effort days deleted
title2 = paste("Seasonal Presence at",site)
plot1 = ggplot(binPresence, aes(x=Season, y=FeHoursProp, color = Season))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_color_brewer(palette = "Dark2")
plot2 = ggplot(binPresence, aes(x=Season, y=JuHoursProp, color = Season))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_color_brewer(palette = "Dark2")
plot3 = ggplot(binPresence, aes(x=Season, y=MaHoursProp, color = Season))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_color_brewer(palette = "Dark2")
figure = ggarrange(plot1,plot2,plot3, labels = c("Social Units","  Mid-Size  ","    Males   "), align = "hv", ncol = 1, nrow = 3, legend = "right",common.legend = TRUE)
annotate_figure(figure, top = text_grob(title1, face = "bold", size = 14), bottom = text_grob("Seasons"),
                left = text_grob("Proportion of Hours/Day w/Clicks", rot = 90))
fig2 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"BoxPlot_StackedGroups.png",sep="")
ggsave(fig2)

#grouping data according to ITS
if (site == 'PT'){
  n = 3
  GroupedDayF = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 3
  GroupedDayJ = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 3
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}

if (site == 'BD'){
  n = 5
  GroupedDayF = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 4
  GroupedDayJ = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 4
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}

if (site == 'QN'){
  GroupedDayF = binPresence;
  n = 7
  GroupedDayJ = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 11
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}

if (site == 'CB'){
  GroupedDayF = binPresence;
  n = 9
  GroupedDayJ = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 6
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}

if (site == 'AB'){
  GroupedDayF = binPresence;
  n = 3
  GroupedDayJ = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 4
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}

if (site == 'KOA'){
  n = 2
  GroupedDayJ = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 2
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}

if (site == 'KS'){
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}

#round day, year, month, and find season for ITS data
if (nrow(binPresence) > nrow(GroupedDayF)){
if (exists('GroupedDayF')){
  GroupedDayF$day = floor(GroupedDayF$day)
  GroupedDayF$month = floor(GroupedDayF$month)
  GroupedDayF$Year = floor(GroupedDayF$Year)
  GroupedDayF$Season[GroupedDayF$month == 1 | GroupedDayF$month == 2 | GroupedDayF$month == 3] = 1
  GroupedDayF$Season[GroupedDayF$month == 4 | GroupedDayF$month == 5 | GroupedDayF$month == 6] = 2
  GroupedDayF$Season[GroupedDayF$month == 7 | GroupedDayF$month == 8 | GroupedDayF$month == 9] = 3
  GroupedDayF$Season[GroupedDayF$month == 10 | GroupedDayF$month == 11 | GroupedDayF$month == 12] = 4
  GroupedDayF$Season = as.factor(GroupedDayF$Season) #change season from an integer to a factor
  GroupedDayF$Season = revalue(GroupedDayF$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall"))
}
}

if (nrow(binPresence) > nrow(GroupedDayJ)){
  GroupedDayJ$day = floor(GroupedDayJ$day)
  GroupedDayJ$month = floor(GroupedDayJ$month)
  GroupedDayJ$Year = floor(GroupedDayJ$Year)
  GroupedDayJ$Season[GroupedDayJ$month == 1 | GroupedDayJ$month == 2 | GroupedDayJ$month == 3] = 1
  GroupedDayJ$Season[GroupedDayJ$month == 4 | GroupedDayJ$month == 5 | GroupedDayJ$month == 6] = 2
  GroupedDayJ$Season[GroupedDayJ$month == 7 | GroupedDayJ$month == 8 | GroupedDayJ$month == 9] = 3
  GroupedDayJ$Season[GroupedDayJ$month == 10 | GroupedDayJ$month == 11 | GroupedDayJ$month == 12] = 4
  GroupedDayJ$Season = as.factor(GroupedDayJ$Season) #change season from an integer to a factor
  GroupedDayJ$Season = revalue(GroupedDayJ$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall"))
}

if (nrow(binPresence) > nrow(GroupedDayM)){
  GroupedDayM$day = floor(GroupedDayM$day)
  GroupedDayM$month = floor(GroupedDayM$month)
  GroupedDayM$Year = floor(GroupedDayM$Year)
  GroupedDayM$Season[GroupedDayM$month == 1 | GroupedDayM$month == 2 | GroupedDayM$month == 3] = 1
  GroupedDayM$Season[GroupedDayM$month == 4 | GroupedDayM$month == 5 | GroupedDayM$month == 6] = 2
  GroupedDayM$Season[GroupedDayM$month == 7 | GroupedDayM$month == 8 | GroupedDayM$month == 9] = 3
  GroupedDayM$Season[GroupedDayM$month == 10 | GroupedDayM$month == 11 | GroupedDayM$month == 12] = 4
  GroupedDayM$Season = as.factor(GroupedDayM$Season) #change season from an integer to a factor
  GroupedDayM$Season = revalue(GroupedDayM$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall"))
}

#plot data grouped with ITS as proportion of hours per day with clicks
title1 = paste(site,"Proprtion of Hours/Day w/ Clicks - ITS")
plot1 = ggplot(GroupedDayF, aes(x=tbin,y=FeHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot2 = ggplot(GroupedDayJ, aes(x=tbin,y=JuHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot3 = ggplot(GroupedDayM, aes(x=tbin,y=MaHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
figure = ggarrange(plot1,plot2,plot3, labels = c("Social Units","  Mid-Size  ","    Males   "),align = "v",ncol = 1, nrow = 3)
annotate_figure(figure, top = text_grob(title1, face = "bold", size = 14), bottom = text_grob("Time (years)"),
                left = text_grob("Proportion of Hours/Day w/Clicks", rot = 90))
fig1 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"HoursProp_TimeSeriesITS_StackedGroups.png",sep="")
ggsave(fig1)

##### grouped data by day of year - mean
if (nrow(binPresence) >= 365){
filename2 = paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"_365GroupedMeanFemale.csv",sep="")
filename3 = paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"_365GroupedMeanJuvenile.csv",sep="")
filename4 = paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"_365GroupedMeanMale.csv",sep="")
oneyearF = read.csv(filename2) #bin means from days
oneyearJ = read.csv(filename3) #bin means from days
oneyearM = read.csv(filename4) #bin means from days
} else {
  filename2 = paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"_365GroupedMean.csv",sep="")
  oneyear = read.csv(filename2)
}

if (nrow(oneyearF) >= 365) {
names(oneyearF) <- c("Day", "HoursPropFE", "SEM", "Std", "Variance", "Range", 'Season', "Month")
oneyearF$Season = as.factor(oneyearF$Season) #change season from an integer to a factor
levels(oneyearF$Season)
oneyearF$Season = revalue(oneyearF$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
names(oneyearJ) <- c("Day", "HoursPropJU", "SEM", "Std", "Variance", "Range", 'Season', "Month")
oneyearJ$Season = as.factor(oneyearJ$Season) #change season from an integer to a factor
levels(oneyearJ$Season)
oneyearJ$Season = revalue(oneyearJ$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
names(oneyearM) <- c("Day", "HoursPropMA", "SEM", "Std", "Variance", "Range", 'Season', "Month")
oneyearM$Season = as.factor(oneyearM$Season) #change season from an integer to a factor
levels(oneyearM$Season)
oneyearM$Season = revalue(oneyearM$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
} else {
oneyear$Season = as.factor(oneyear$Season) #change season from an integer to a factor
levels(oneyear$Season)
oneyear$Season = revalue(oneyear$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
head(oneyear)
str(oneyear)
}

#groupin data according to ITS
if (site == 'PT'){
  n = 5
  GroupedYearF = aggregate(oneyearF,list(rep(1:(nrow(oneyearF)%/%n+1),each=n,len=nrow(oneyearF))),mean)[-1];
  n = 2
  GroupedYearJ = aggregate(oneyearJ,list(rep(1:(nrow(oneyearJ)%/%n+1),each=n,len=nrow(oneyearJ))),mean)[-1];
  n = 3
  GroupedYearM = aggregate(oneyearM,list(rep(1:(nrow(oneyearM)%/%n+1),each=n,len=nrow(oneyearM))),mean)[-1];
}

if (site == 'BD'){
  n = 7
  GroupedYearF = aggregate(oneyearF,list(rep(1:(nrow(oneyearF)%/%n+1),each=n,len=nrow(oneyearF))),mean)[-1];
  n = 5
  GroupedYearJ = aggregate(oneyearJ,list(rep(1:(nrow(oneyearJ)%/%n+1),each=n,len=nrow(oneyearJ))),mean)[-1];
  n = 4
  GroupedYearM = aggregate(oneyearM,list(rep(1:(nrow(oneyearM)%/%n+1),each=n,len=nrow(oneyearM))),mean)[-1];
}

if (site == 'QN'){
  GroupedYearF = oneyearF;
  n = 2
  GroupedYearJ = aggregate(oneyearJ,list(rep(1:(nrow(oneyearJ)%/%n+1),each=n,len=nrow(oneyearJ))),mean)[-1];
  n = 4
  GroupedYearM = aggregate(oneyearM,list(rep(1:(nrow(oneyearM)%/%n+1),each=n,len=nrow(oneyearM))),mean)[-1];
}

if (site == 'CB'){
  GroupedYearF = oneyearF;
  n = 21
  GroupedYearJ = aggregate(oneyearJ,list(rep(1:(nrow(oneyearJ)%/%n+1),each=n,len=nrow(oneyearJ))),mean)[-1];
  n = 14
  GroupedYearM = aggregate(oneyearM,list(rep(1:(nrow(oneyearM)%/%n+1),each=n,len=nrow(oneyearM))),mean)[-1];
}

if (site == 'AB'){
  GroupedYearF = oneyear;
  n = 3
  GroupedYearJ = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
  n = 4
  GroupedYearM = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'KOA'){
  n = 2
  GroupedYearJ = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
  n = 2
  GroupedYearM = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'KS'){
  GroupedYearJ = oneyear;
  GroupedYearM = oneyear;
}

if (exists('oneyearF')){
  oneyear = oneyearF
}

#round day, year, month, and find season for ITS data
if (nrow(oneyear) == nrow(GroupedYearF)){
  if (exists('GroupedYearF')){
    GroupedYearF$Day = floor(GroupedYearF$Day)
    GroupedYearF$Month = floor(GroupedYearF$Month)
    GroupedYearF$Season = as.integer(GroupedYearF$Season)
    GroupedYearF$Season[GroupedYearF$Month == 1 | GroupedYearF$Month == 2 | GroupedYearF$Month == 3] = 1
    GroupedYearF$Season[GroupedYearF$Month == 4 | GroupedYearF$Month == 5 | GroupedYearF$Month == 6] = 2
    GroupedYearF$Season[GroupedYearF$Month == 7 | GroupedYearF$Month == 8 | GroupedYearF$Month == 9] = 3
    GroupedYearF$Season[GroupedYearF$Month == 10 | GroupedYearF$Month == 11 | GroupedYearF$Month == 12] = 4
    GroupedYearF$Season = as.factor(GroupedYearF$Season) #change season from an integer to a factor
    GroupedYearF$Season = revalue(GroupedYearF$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall"))
  }
}

if (nrow(oneyear) == nrow(GroupedYearJ)){
  GroupedYearJ$Day = floor(GroupedYearJ$Day)
  GroupedYearJ$Month = floor(GroupedYearJ$Month)
  GroupedYearJ$Season = as.integer(GroupedYearJ$Season)
  GroupedYearJ$Season[GroupedYearJ$Month == 1 | GroupedYearJ$Month == 2 | GroupedYearJ$Month == 3] = 1
  GroupedYearJ$Season[GroupedYearJ$Month == 4 | GroupedYearJ$Month == 5 | GroupedYearJ$Month == 6] = 2
  GroupedYearJ$Season[GroupedYearJ$Month == 7 | GroupedYearJ$Month == 8 | GroupedYearJ$Month == 9] = 3
  GroupedYearJ$Season[GroupedYearJ$Month == 10 | GroupedYearJ$Month == 11 | GroupedYearJ$Month == 12] = 4
  GroupedYearJ$Season = as.factor(GroupedYearJ$Season) #change season from an integer to a factor
  GroupedYearJ$Season = revalue(GroupedYearJ$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall"))
}

if (nrow(oneyear) == nrow(GroupedYearM)){
  GroupedYearM$Day = floor(GroupedYearM$Day)
  GroupedYearM$Month = floor(GroupedYearM$Month)
  GroupedYearJ$Season = as.integer(GroupedYearJ$Season)
  GroupedYearM$Season[GroupedYearM$Month == 1 | GroupedYearM$Month == 2 | GroupedYearM$Month == 3] = 1
  GroupedYearM$Season[GroupedYearM$Month == 4 | GroupedYearM$Month == 5 | GroupedYearM$Month == 6] = 2
  GroupedYearM$Season[GroupedYearM$Month == 7 | GroupedYearM$Month == 8 | GroupedYearM$Month == 9] = 3
  GroupedYearM$Season[GroupedYearM$Month == 10 | GroupedYearM$Month == 11 | GroupedYearM$Month == 12] = 4
  GroupedYearM$Season = as.factor(GroupedYearM$Season) #change season from an integer to a factor
  GroupedYearM$Season = revalue(GroupedYearM$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall"))
}

#plot data as time series
#plot data as proportion of hours per day with clicks
title3 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series")
plot1 = ggplot(oneyearF, aes(x=Day,y=HoursPropFE))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot2 = ggplot(oneyearJ, aes(x=Day,y=HoursPropJU))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot3 = ggplot(oneyearM, aes(x=Day,y=HoursPropMA))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
figure = ggarrange(plot1,plot2,plot3, labels = c("Social Units","  Mid-Size  ","    Males   "),align = "v",ncol = 1, nrow = 3)
annotate_figure(figure, top = text_grob(title1, face = "bold", size = 14), bottom = text_grob("Time (years)"),
                left = text_grob("Proportion of Hours/Day w/Clicks", rot = 90))
fig1 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"AveragedHoursProp_TimeSeries_StackedGroups.png",sep="")
ggsave(fig1)

#plot data grouped with ITS as proportion of hours per day with clicks
title1 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series (ITS)")
plot1 = ggplot(GroupedYearF, aes(x=Day,y=HoursPropFE))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot2 = ggplot(GroupedYearJ, aes(x=Day,y=HoursPropJU))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot3 = ggplot(GroupedYearM, aes(x=Day,y=HoursPropMA))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
figure = ggarrange(plot1,plot2,plot3, labels = c("Social Units","  Mid-Size  ","    Males   "),align = "v",ncol = 1, nrow = 3)
annotate_figure(figure, top = text_grob(title1, face = "bold", size = 14), bottom = text_grob("Time (years)"),
                left = text_grob("Proportion of Hours/Day w/Clicks", rot = 90))
fig1 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"AveragedHoursProp_TimeSeriesITS_StackedGroups.png",sep="")
ggsave(fig1)

if (nrow(oneyear) >= 365) {
#plot data as time series with error bars
title4 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series")
plot1 = ggplot(oneyearF, aes(x=Day,y=HoursPropFE))+
  geom_errorbar(aes(ymin = HoursPropFE - SEM, ymax = HoursPropFE + SEM))+
  geom_line()+
  geom_point()+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot2 = ggplot(oneyearJ, aes(x=Day,y=HoursPropJU))+
  geom_errorbar(aes(ymin = HoursPropJU - SEM, ymax = HoursPropJU + SEM))+
  geom_line()+
  geom_point()+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
plot3 = ggplot(oneyearM, aes(x=Day,y=HoursPropMA))+
  geom_errorbar(aes(ymin = HoursPropMA - SEM, ymax = HoursPropMA + SEM))+
  geom_line()+
  geom_point()+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
figure = ggarrange(plot1,plot2,plot3, labels = c("Social Units","  Mid-Size  ","    Males   "),align = "v",ncol = 1, nrow = 3)
annotate_figure(figure, top = text_grob(title1, face = "bold", size = 14), bottom = text_grob("Time (years)"),
                left = text_grob("Proportion of Hours/Day w/Clicks", rot = 90))
fig5 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"AveragedHoursProp_TimeSeries_ErrorBars_StackedGroups.png",sep="")
ggsave(fig5)
}

## GAMs ##

#Social Units#
#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
#females
gamTw = gam(FeHoursProp ~ s(day, bs = 'cc', k = 50), data = GroupedDayF, family = tw, method = "REML")
plot(gamTw, pages =1)
summary(gamTw)

#GAM to check for significance between seasons
gamTwS = gam(FeHoursProp ~ Season, data = GroupedDayF, family = tw, method = "REML")
summary(gamTwS)

#Better GAM plots
#pattern only
viz = getViz(gamTw)
print(plot(viz,allTerms=T),pages=1)

#first way to plot GAM
vizGG = plot(viz,allTerms = T) +
  l_points() +
  l_fitLine(linetype = 3)  +
  l_fitContour()+
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciBar() +
  l_points(shape = 19, size = 1, alpha = 0.1) +
  l_rug() +
  theme_get() 
print(vizGG,pages =1)
fig6 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"GAM1_SocialUnits.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
print(vizGG2,pages =1)
fig7 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"GAM2_SocialUnits.png",sep="")
ggsave(fig7)

#Juveniles#
#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
gamTw = gam(JuHoursProp ~ s(day, bs = 'cc', k = 50), data = GroupedDayJ, family = tw, method = "REML")
plot(gamTw, pages =1)
summary(gamTw)

#GAM to check for significance between seasons
gamTwS = gam(JuHoursProp ~ Season, data = GroupedDayJ, family = tw, method = "REML")
summary(gamTwS)

#Better GAM plots
#pattern only
viz = getViz(gamTw)
print(plot(viz,allTerms=T),pages=1)

#first way to plot GAM
vizGG = plot(viz,allTerms = T) +
  l_points() +
  l_fitLine(linetype = 3)  +
  l_fitContour()+
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciBar() +
  l_points(shape = 19, size = 1, alpha = 0.1) +
  l_rug() +
  theme_get() 
print(vizGG,pages =1)
fig6 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"GAM1_Juveniles.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
print(vizGG2,pages =1)
fig7 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"GAM2_Juveniles.png",sep="")
ggsave(fig7)

#Males#
#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
gamTw = gam(MaHoursProp ~ s(day, bs = 'cc', k = 50), data = GroupedDayM, family = tw, method = "REML")
plot(gamTw, pages =1)
summary(gamTw)

#GAM to check for significance between seasons
gamTwS = gam(MaHoursProp ~ Season, data = GroupedDayM, family = tw, method = "REML")
summary(gamTwS)

#Better GAM plots
#pattern only
viz = getViz(gamTw)
print(plot(viz,allTerms=T),pages=1)

#first way to plot GAM
vizGG = plot(viz,allTerms = T) +
  l_points() +
  l_fitLine(linetype = 3)  +
  l_fitContour()+
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciBar() +
  l_points(shape = 19, size = 1, alpha = 0.1) +
  l_rug() +
  theme_get() 
print(vizGG,pages =1)
fig6 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"GAM1_Males.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
print(vizGG2,pages =1)
fig7 =paste("G:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,"GAM2_Males.png",sep="")
ggsave(fig7)






















#GAM w/ GEE

#GLMs

#GLM in days
hist(dayBinTAB$HoursProp)
glm_day = glm(HoursProp ~ day, data = dayBinTAB, family = poisson())
plot(glm_day)
summary(glm_day)

#other code
#fitting a simple GEE to view the output
gee_simp = GEE(HoursProp ~ day + Season, family = "poisson", data = oneyear)
summary(gee_simp)

#looks for statistical significane between seasons using a GLM
GLM_form = formula(dayBinTAB$HoursProp ~ dayBinTAB$Season)
gee1 = geeglm(GLM_form, data = dayBinTAB, id = Season, family=poisson("identity"))
gee1
coef(gee1)
vcov(gee1)
summary(gee1)







