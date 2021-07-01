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

#load data from MATLAB
site = 'Wake'
saveDir = paste('G:/My Drive/CentralPac_TPWS_metadataReduced/Wake/Seasonality/')
filename = paste(saveDir,site,"_binPresence.csv",sep="")
binPresence = read.csv(filename) #no effort days deleted
head(binPresence)
str(binPresence)
binPresence$Season = as.factor(binPresence$Season) #change season from an integer to a factor
levels(binPresence$Season)
binPresence$Season = revalue(binPresence$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
binPresence$tbin = anytime(as.factor(binPresence$tbin))

#load data from StatisticalAnalysis_Gender: Females
filenameStatGender_F = paste(saveDir,site,"_GroupedDayF.csv",sep="")
GroupedDayF = read.csv(filenameStatGender_F) #no effort days deleted

#load data from StatisticalAnalysis_Gender: Juveniles
filenameStatGender_J= paste(saveDir,site,"_GroupedDayJ.csv",sep="")
GroupedDayJ = read.csv(filenameStatGender_J) #no effort days deleted

#load data from StatisticalAnalysis_Gender: Males
filenameStatGender_M = paste(saveDir,site,"_GroupedDayM.csv",sep="")
GroupedDayM = read.csv(filenameStatGender_M) #no effort days deleted

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
fig1 =paste(saveDir,site,site,"HoursProp_TimeSeries_StackedGroups.png",sep="")
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
fig2 =paste(saveDir,site,"BoxPlot_StackedGroups.png",sep="")
ggsave(fig2)


#plot data grouped with ITS as proportion of hours per day with clicks
title1 = paste(site,"Proprtion of Hours/Day w/ Clicks - ITS")
theme(axis.text=element_text(size=18),
      axis.title=element_text(size=20,face="bold"))
title1 = paste("Proportion of Hours per Day with Clicks at Buldir Island, BSAI")
plot1 = ggplot(GroupedDayF, aes(x=tbin,y=FeHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank(),axis.text=element_text(size=18))+
  theme(axis.title.y = element_blank(),axis.text=element_text(size=18))
plot2 = ggplot(GroupedDayJ, aes(x=tbin,y=JuHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank(),axis.text=element_text(size=18))+
  theme(axis.title.y = element_blank(),axis.text=element_text(size=18))
plot3 = ggplot(GroupedDayM, aes(x=tbin,y=MaHoursProp))+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank(),axis.text=element_text(size=18))+
  theme(axis.title.y = element_blank(),axis.text=element_text(size=18))
figure = ggarrange(plot1,plot2,plot3, labels = c("Social Units","  Mid-Size  ","    Males   "),align = "v",ncol = 1, nrow = 3)
annotate_figure(figure, top = text_grob(title1, face = "bold", size = 20), bottom = text_grob("Time (years)", size = 20),
                left = text_grob("Proportion of Hours/Day w/Clicks", rot = 90, size = 24))
fig1 =paste(saveDir,site,"HoursProp_TimeSeriesITS_StackedGroups.png",sep="")
ggsave(fig1)

##### grouped data by day of year - mean
  filename2 = paste(saveDir,site,"_365GroupedMean.csv",sep="")
  oneyear = read.csv(filename2)

  filenameGYF = paste(saveDir,site,"_GroupedYearF.csv",sep="")
  GroupedYearF = read.csv(filenameGYF)
  
  filenameGYJ = paste(saveDir,site,"_GroupedYearJ.csv",sep="")
  GroupedYearJ = read.csv(filenameGYJ)
  
  filenameGYM = paste(saveDir,site,"_GroupedYearM.csv",sep="")
  GroupedYearM = read.csv(filenameGYM)

#plot data as time series
#plot data as proportion of hours per day with clicks
if (exists('oneyearF')){
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
  fig1 = paste(saveDir,site,"AveragedHoursProp_TimeSeries_StackedGroups.png",sep="")
  ggsave(fig1)
}

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
fig1 =paste(saveDir,site,"AveragedHoursProp_TimeSeriesITS_StackedGroups.png",sep="")
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
  fig5 =paste(saveDir,site,"AveragedHoursProp_TimeSeries_ErrorBars_StackedGroups.png",sep="")
  ggsave(fig5)
}
