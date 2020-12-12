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
library(geepack)

#load data
site = 'CB'
filename = paste("E:/Project_Sites/",site,"/Seasonality/",site,"_dayData_forGLMR125_Boreal.csv",sep="")
dayBinTAB = read.csv(filename) #no effort days deleted
head(dayBinTAB)
head(mean365)
str(dayBinTAB)
str(mean365)
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
filename2 = paste("E:/Project_Sites/",site,"/Seasonality/",site,"_days365GroupedMean_forGLMR125_Boreal.csv",sep="")
oneyear = read.csv(filename2) #bin means from days
names(oneyear) <- c("Day", "HoursProp", "SEM", "Std", "Variance", "Range", 'Season', "Month")
oneyear$Season = as.factor(oneyear$Season) #change season from an integer to a factor
levels(oneyear$Season)
oneyear$Season = revalue(oneyear$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
head(oneyear)
str(oneyear)

#plot data as time series
title3 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series")
ggplot(oneyear, aes(x=Day,y=HoursProp))+
  ggtitle(title3)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_line()+
  geom_point()
fig4 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"AveragedHoursProp_TimeSeries.png",sep="")
ggsave(fig4)

#plot data as time series with error bars
title4 = paste(site,"Yearly Average of Proportion of Hours per Day with Clicks - Time Series")
ggplot(oneyear, aes(x=Day,y=HoursProp))+
  ggtitle(title4)+
  labs(y="Average Proportion of Hours per Day with Clicks",x="Day of the Year")+
  geom_errorbar(aes(ymin = HoursProp - SEM, ymax = HoursProp + SEM))
geom_line()+
  geom_point()
fig5 =paste("E:/Project_Sites/",site,"/Seasonality/",site,"AveragedHoursProp_TimeSeries_ErrorBars.png",sep="")
ggsave(fig5)

#Generalized Linear Models

GLM_form = formula(oneyear$HoursProp ~ oneyear$Season)
gee1 = geeglm(GLM_form, data = oneyear, id = Season, family=poisson("identity"))
gee1
coef(gee1)
vcov(gee1)
summary(gee1)

#GAMs and GLMs for Seasonality
#Gams in days
gamTw = gam(HoursProp ~ s(Day, bs = 'cc', k = 100), data = dayBinTAB, family = tw, method = "REML")
plot(gamTw, pages =1)
summary(gamTw)

gamTwS = gam(HoursProp ~ Season, data = dayBinTAB, family = tw, method = "REML")
plot(gamTwS, pages =1)
summary(gamTwS)


#Glm in days
hist(dayBinTAB$HoursProp)
glm_day = glm(HoursProp ~ day, data = dayBinTAB, family = tweedie)
plot(glm_day, pages = 1)
summary(glm_day)

viz = getViz(gamTw)
print(plot(viz,allTerms=T),pages=1)

viz_glm = getViz(glm_day)
print(plot(viz_glm,allTerms=T),pages=1)

vizGG = plot(viz,allTerms = T) +
  l_points() +
  l_fitLine(linetype = 3)  +
  l_fitContour()+
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciBar() +
  l_points(shape = 19, size = 3, alpha = 0.1) +
  l_rug() +
  theme_get() 
print(vizGG,pages =1)



main2 = paste("Grouped Loess Smoothing and Prediction for",site)
plot(oneyear$Bin, x =oneyear$Day, type="l", main=main2,xlab="Day of Year",ylab="Presence of Sperm Whales in 5 minute Bins (Median)")
lines(smoothed10, x =oneyear$Day, col="red")
lines(smoothed25, x =oneyear$Day, col="green")
lines(smoothed50, x =oneyear$Day, col="blue")
legend("topright",inset = 0.02, title = "Smoothing Span",legend=c("0.10", "0.25", "0.50"),col=c("red","green","blue"),lty =1)
fig5 =paste("F:/Seasonality/",site,"/",site,"GroupedLoess.png",sep="")
png(fig5)

#Fit A Local Polynomial Regression With Automatic Smoothing Parameter Selection
if (site == 'CORC'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 0, criterion = c("gcv"), family = "gaussian", plot = TRUE)
}
if (site == 'BD'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 1, criterion = c("gcv"), family = "symmetric", plot = TRUE)
}
if (site == 'QN'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 1, criterion = c("gcv"), family = "symmetric", plot = TRUE)
}
if (site == 'PT'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 2, criterion = c("gcv"), family = "gaussian", plot = TRUE)
}
if (site == 'CB'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 2, criterion = c("aicc"), family = "symmetric", plot = TRUE)
}
if (site == 'QC'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 2, criterion = c("gcv"), family = "gaussian", plot = TRUE)
}
if (site == 'CCE'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 2, criterion = c("gcv"), family = "gaussian", plot = TRUE)
}
if (site == 'AB'){
  AIC = loess.as(oneyear$Day,oneyear$Bin,degree = 2, criterion = c("gcv"), family = "symmetric", plot = TRUE)
}

smoothedAIC = predict(AIC)

main3 = paste("Best Grouped Loess Smoothing and Prediction for",site)
BestSpan = round(AIC$pars$span,2)
BestSpan
leg = paste("Best Span =", BestSpan)
plot(oneyear$Bin, x =oneyear$Day, type="l", main=main3,xlab="Day of Year",ylab="Presence in 5 minute Bins")
lines(smoothedAIC, x =oneyear$Day, col="red")
legend("topright",inset = 0.02,legend=c(leg),col=c("red"),lty =1)
fig7=paste("F:/Seasonality/",site,"/",site,"BestSpan_RPlot.png",sep="")
png(fig7)

title4 = paste("Best Grouped Loess Smoothing and Prediction for",site)
ggplot(oneyear, aes(x=Day,y=Bin))+
  ggtitle(title4)+
  labs(y="Daily Presence in 5 minute bins",x="Day of the Year")+
  geom_rect(data = NULL, aes(xmin = 1, xmax = 60, ymin = -Inf, ymax = Inf), fill = 'lightblue',alpha = 0.01)+
  geom_rect(data = NULL, aes(xmin = 61, xmax = 152, ymin = -Inf, ymax = Inf), fill = 'lightgreen',alpha = 0.009)+
  geom_rect(data = NULL, aes(xmin = 153, xmax = 244, ymin = -Inf, ymax = Inf), fill = 'mistyrose1',alpha = 0.025)+
  geom_rect(data = NULL, aes(xmin = 245, xmax = 335, ymin = -Inf, ymax = Inf), fill = 'navajowhite',alpha = 0.01)+
  geom_rect(data = NULL, aes(xmin = 336, xmax = 366, ymin = -Inf, ymax = Inf), fill = 'lightblue',alpha = 0.01)+
  geom_point()+
  geom_line()+
  geom_line(color = "red",aes(x=oneyear$Day,smoothedAIC),size = 2)+
  scale_x_continuous(limits = c(0,366), expand = c(0, 0)) +
  theme(panel.border = element_blank())
fig6 =paste("F:/Seasonality/",site,"/",site,"BestSpan_ggplot.png",sep="")
ggsave(fig6)

#plot data as time series with seasonal trend
dayBinTAB$julian = strftime(dayBinTAB$tbin,format = "%j")
dayBinTAB$julian = as.integer(dayBinTAB$julian)
trend = as.data.frame(smoothedAIC)
trend$julian = as.integer(1:length(smoothedAIC))
savetrend = paste("F:/Seasonality/",site,"/",site,"_onlytrend.csv",sep="")
write.csv(trend,savetrend)

GrandTab = right_join(dayBinTAB, trend, by ="julian")
GrandTab2 = GrandTab[order(GrandTab$tbin),]
saveGT2 = paste("F:/Seasonality/",site,"/",site,"_smoothedAIC.csv",sep="")
write.csv(GrandTab2,saveGT2)

title5 = paste(site,"Presence Time Series with Loess Regression")
ggplot(GrandTab2, aes(x=tbin,y=NormBin))+
  ggtitle(title5)+
  labs(y="Daily Presence (5 minute bins)",x="Time (days)")+
  geom_point()+
  geom_line()+
  geom_line(color = "red",aes(x=GrandTab2$tbin,GrandTab2$smoothedAIC),size = 2)+
  theme(panel.border = element_blank())
fig7 =paste("F:/Seasonality/",site,"/",site,"TimeSeries_Loess.png",sep="")
ggsave(fig7)

graphics.off()
