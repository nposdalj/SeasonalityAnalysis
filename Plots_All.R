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

#load data from MATLAB
site = 'QC'
saveDir = paste("I:/My Drive/CentralPac_TPWS_metadataReduced/QC/Seasonality/")
filename = paste(saveDir,site,"_dayData_forGLMR125.csv",sep="")
dayBinTAB = read.csv(filename) #no effort days deleted
head(dayBinTAB)
str(dayBinTAB)
dayBinTAB$Season = as.factor(dayBinTAB$Season) #change season from an integer to a factor
levels(dayBinTAB$Season)
dayBinTAB$Season = revalue(dayBinTAB$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
dayBinTAB$tbin = anytime(as.factor(dayBinTAB$tbin))

#load data from StatisicalAnalysis_All
filenameStatAll = paste(saveDir,site,"_GroupedDay.csv",sep="")
GroupedDay = read.csv(filenameStatAll) #no effort days deleted

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

filenameGY = paste(saveDir,site,"_GroupedYear.csv",sep="")
GroupedYear = read.csv(filenameGY)

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


