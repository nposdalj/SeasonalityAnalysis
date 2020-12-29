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
site = 'KOA'
filename = paste("E:/Project_Sites/",site,"/Seasonality/",site,"_dayData_forGLMR125.csv",sep="")
dayBinTAB = read.csv(filename) #no effort days deleted
head(dayBinTAB)
str(dayBinTAB)
dayBinTAB$Season = as.factor(dayBinTAB$Season) #change season from an integer to a factor
levels(dayBinTAB$Season)
dayBinTAB$Season = revalue(dayBinTAB$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
dayBinTAB$tbin = anytime(as.factor(dayBinTAB$tbin))

#plot data as proportion of hours per day with clicks
title1 = paste(site,"Proprtion of Hours/Day w/ Clicks")
ggplot(dayBinTAB, aes(x=tbin,y=HoursProp))+
  ggtitle(title1)+
  labs(y="Proportion of Hours/Day w/Clicks",x="Time (days)")+
  geom_line()+
  geom_point()
fig1 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"HoursProp_TimeSeries.png",sep="")
ggsave(fig1)

#plot data as box plot for seasons; have to plot this with no effort days deleted
title2 = paste("Seasonal Presence at",site)
ggplot(dayBinTAB, aes(x=Season, y=HoursProp, color = Season))+
  geom_boxplot()+
  ggtitle(title2)+
  labs(y="Proportion of Hours/Day w/Clicks")
  scale_color_brewer(palette = "Dark2")
fig2 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"BoxPlot.png",sep="")
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

#round day, year, month, and find season for ITS data
GroupedDay$day = floor(GroupedDay$day)
GroupedDay$month = floor(GroupedDay$month)
GroupedDay$Year = floor(GroupedDay$Year)
GroupedDay$Season[GroupedDay$month == 1 | GroupedDay$month == 2 | GroupedDay$month == 3] = 1
GroupedDay$Season[GroupedDay$month == 4 | GroupedDay$month == 5 | GroupedDay$month == 6] = 2
GroupedDay$Season[GroupedDay$month == 7 | GroupedDay$month == 8 | GroupedDay$month == 9] = 3
GroupedDay$Season[GroupedDay$month == 10 | GroupedDay$month == 11 | GroupedDay$month == 12] = 4
GroupedDay$Season = as.factor(GroupedDay$Season) #change season from an integer to a factor
GroupedDay$Season = revalue(GroupedDay$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall")) #change the numbers in actual seasons

#plot data grouped with ITS as proportion of hours per day with clicks
title1 = paste(site,"Proprtion of Hours/Day w/ Clicks - ITS")
ggplot(GroupedDay, aes(x=tbin,y=HoursProp))+
  ggtitle(title1)+
  labs(y="Proportion of Hours/Day w/Clicks",x="Time (days)")+
  geom_line()+
  geom_point()
fig1 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"HoursProp_TimeSeriesITS.png",sep="")
ggsave(fig1)

##### grouped data by day of year - mean
filename2 = paste("E:/Project_Sites/",site,"/Seasonality/",site,"_days365GroupedMean_forGLMR125.csv",sep="")
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

#round day, year, month, and find season for ITS data
GroupedYear$Month = month(as.Date(GroupedYear$Day, origin = "2014-01-01"))
GroupedYear$Season[GroupedYear$Month == 1 | GroupedYear$Month == 2 | GroupedYear$Month == 3] = 1
GroupedYear$Season[GroupedYear$Month == 4 | GroupedYear$Month == 5 | GroupedYear$Month == 6] = 2
GroupedYear$Season[GroupedYear$Month == 7 | GroupedYear$Month == 8 | GroupedYear$Month == 9] = 3
GroupedYear$Season[GroupedYear$Month == 10 | GroupedYear$Month == 11 | GroupedYear$Month == 12] = 4
GroupedYear$Season = as.factor(GroupedYear$Season) #change season from an integer to a factor
GroupedYear$Season = revalue(GroupedYear$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall")) #change the numbers in actual seasons

#plot data as time series
title3 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series")
ggplot(oneyear, aes(x=Day,y=HoursProp))+
  ggtitle(title3)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_line()+
  geom_point()
fig4 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"AveragedHoursProp_TimeSeries.png",sep="")
ggsave(fig4)

#plot data as time series with ITS
title3 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series (ITS)")
ggplot(GroupedYear, aes(x=Day,y=HoursProp))+
  ggtitle(title3)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_line()+
  geom_point()
fig4 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"AveragedHoursProp_TimeSeriesITS.png",sep="")
ggsave(fig4)

if (nrow(oneyear) >= 365) {
#plot data as time series with error bars
title4 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series")
ggplot(oneyear, aes(x=Day,y=HoursProp))+
  ggtitle(title4)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_errorbar(aes(ymin = HoursProp - SEM, ymax = HoursProp + SEM))+
  geom_line()+
  geom_point()
fig5 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"AveragedHoursProp_TimeSeries_ErrorBars.png",sep="")
ggsave(fig5)
}

## GAMs ##

#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
gamTw = gam(HoursProp ~ s(day, bs = 'cc', k = 50), data = GroupedDay, family = tw, method = "REML")
plot(gamTw, pages =1)
summary(gamTw)

#GAM to check for significance between seasons
gamTwS = gam(HoursProp ~ Season, data = GroupedDay, family = tw, method = "REML")
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
fig6 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"GAM1.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
print(vizGG2,pages =1)
fig7 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"GAM2.png",sep="")
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







