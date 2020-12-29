# Libraries
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

#load data
site = 'KS'
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

#plot data as time series
title3 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series")
ggplot(oneyear, aes(x=Day,y=HoursProp))+
  ggtitle(title3)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_line()+
  geom_point()
fig4 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"AveragedHoursProp_TimeSeries.png",sep="")
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

#GAMs
#GAM to identify seasonal pattern
gamTw = gam(HoursProp ~ s(day, bs = 'cc', k = 50), data = dayBinTAB, family = tw, method = "REML")
plot(gamTw, pages =1)
summary(gamTw)

#GAM to check for significance between seasons
gamTwS = gam(HoursProp ~ Season, data = dayBinTAB, family = tw, method = "REML")
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



#GLMs
#fitting a simple GEE to view the output
gee_simp = GEE(HoursProp ~ day + Season, family = "poisson", data = oneyear)
summary(gee_simp)

#looks for statistical significane between seasons using a GLM
GLM_form = formula(oneyear$HoursProp ~ oneyear$Season)
gee1 = geeglm(GLM_form, data = oneyear, id = Season, family=poisson("identity"))
gee1
coef(gee1)
vcov(gee1)
summary(gee1)

#Glm in days
hist(dayBinTAB$HoursProp)
glm_day = glm(HoursProp ~ day, data = dayBinTAB, family = tweedie)
plot(glm_day)
summary(glm_day)

viz_glm = getViz(glm_day)
print(plot(viz_glm,allTerms=T),pages=1)

vizGG = plot(viz_glm,allTerms = T) +
  l_points() +
  l_fitLine(linetype = 3)  +
  l_fitContour()+
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciBar() +
  l_points(shape = 19, size = 3, alpha = 0.1) +
  l_rug() +
  theme_get() 
print(vizGG,pages =1)







