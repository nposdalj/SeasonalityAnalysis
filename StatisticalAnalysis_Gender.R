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
site = 'BP'
saveDir = paste("H:/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",site,"/",sep="")
filename = paste(saveDir,site,"_binPresence.csv",sep="")
binPresence = read.csv(filename) #no effort days deleted
head(binPresence)
str(binPresence)
binPresence$Season = as.factor(binPresence$Season) #change season from an integer to a factor
levels(binPresence$Season)
binPresence$Season = revalue(binPresence$Season, c("1"="Summer", "2"="Fall", "3"="Winter", "4"="Spring")) #change the numbers in actual seasons
binPresence$tbin = anytime(as.factor(binPresence$tbin))

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
  GroupedDayF = binPresence;
  n = 2
  GroupedDayJ = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  n = 2
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}

if (site == 'KS' || site == 'GI'|| site=='PG'){
  GroupedDayF = binPresence;
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}
if (site == 'GI'){
  n = 2
  GroupedDayF = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}
if (site == 'CORC'){
  GroupedDayF = binPresence;
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}
if (site == 'TIN'){
  GroupedDayF = binPresence;
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}
if (site == 'SAP'){
  n = 6
  GroupedDayF = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}
if (site == 'QC'){
  GroupedDayF = binPresence;
  GroupedDayJ = binPresence;
  n = 2
  GroupedDayM = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
}
if (site == 'Wake'){
  n = 9
  GroupedDayF = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}

if (site == 'BS'){
  n = 2
  GroupedDayF = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}

if (site == 'BP'){
  GroupedDayF = binPresence
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
}

#export GroupedDayF as .csv
fileName_GDF = paste(saveDir,site,"_GroupedDayF.csv",sep="")
write.csv(GroupedDayF,fileName_GDF, row.names = FALSE)

#export GroupedDayJ as .csv
fileName_GDJ = paste(saveDir,site,"_GroupedDayJ.csv",sep="")
write.csv(GroupedDayJ,fileName_GDJ, row.names = FALSE)

#export GroupedDayM as .csv
fileName_GDM = paste(saveDir,site,"_GroupedDayM.csv",sep="")
write.csv(GroupedDayM,fileName_GDM, row.names = FALSE)

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


##### grouped data by day of year - mean
if (nrow(binPresence) >= 365){
filename2 = paste(saveDir,site,"_365GroupedMeanFemale.csv",sep="")
filename3 = paste(saveDir,site,"_365GroupedMeanJuvenile.csv",sep="")
filename4 = paste(saveDir,site,"_365GroupedMeanMale.csv",sep="")
oneyearF = read.csv(filename2) #bin means from days
oneyearJ = read.csv(filename3) #bin means from days
oneyearM = read.csv(filename4) #bin means from days
} else {
  filename2 = paste(saveDir,site,"_365GroupedMean.csv",sep="")
  oneyear = read.csv(filename2)
}

if (exists('oneyearF')){
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
}} else {
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
  GroupedYearF = oneyear;
  n = 2
  GroupedYearJ = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
  n = 2
  GroupedYearM = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
}

if (site == 'KS' || site == 'GI'){
  GroupedYearF = oneyearF;
  GroupedYearJ = oneyearJ;
  GroupedYearM = oneyearM;
}
if (site == 'GI'){
  n = 2
  GroupedYearF = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
  GroupedYearJ = oneyearJ;
  GroupedYearM = oneyearM;
}
if (site == 'CORC'){
  GroupedYearF = oneyearF;
  GroupedYearJ = oneyearJ;
  GroupedYearM = oneyearM;
}
if (site == 'TIN'){
  GroupedYearF = oneyearF;
  GroupedYearJ = oneyearJ;
  GroupedYearM = oneyearM;
}
if (site == 'SAP'){
  GroupedYearF = oneyearF;
  n = 8
  GroupedYearJ = aggregate(oneyearJ,list(rep(1:(nrow(oneyearJ)%/%n+1),each=n,len=nrow(oneyearJ))),mean)[-1];
  GroupedYearM = oneyearM;
}
if (site == 'QC'){
  n = 2
  GroupedYearF = aggregate(oneyear,list(rep(1:(nrow(oneyear)%/%n+1),each=n,len=nrow(oneyear))),mean)[-1];
  GroupedYearJ = oneyearJ;
  GroupedYearM = oneyearM;
}
if (site == 'Wake'){
  GroupedYearF = oneyearF;
  n = 4
  GroupedYearJ = aggregate(oneyearJ,list(rep(1:(nrow(oneyearJ)%/%n+1),each=n,len=nrow(oneyearJ))),mean)[-1];
  GroupedYearM = oneyearM;
}
if (site == 'BS'){
  n = 2
  GroupedDayF = aggregate(binPresence,list(rep(1:(nrow(binPresence)%/%n+1),each=n,len=nrow(binPresence))),mean)[-1];
  GroupedDayJ = binPresence;
  GroupedDayM = binPresence;
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
  GroupedYearM$Season = as.integer(GroupedYearJ$Season)
  GroupedYearM$Season[GroupedYearM$Month == 1 | GroupedYearM$Month == 2 | GroupedYearM$Month == 3] = 1
  GroupedYearM$Season[GroupedYearM$Month == 4 | GroupedYearM$Month == 5 | GroupedYearM$Month == 6] = 2
  GroupedYearM$Season[GroupedYearM$Month == 7 | GroupedYearM$Month == 8 | GroupedYearM$Month == 9] = 3
  GroupedYearM$Season[GroupedYearM$Month == 10 | GroupedYearM$Month == 11 | GroupedYearM$Month == 12] = 4
  GroupedYearM$Season = as.factor(GroupedYearM$Season) #change season from an integer to a factor
  GroupedYearM$Season = revalue(GroupedYearM$Season, c("1"="Winter", "2"="Spring", "3"="Summer", "4"="Fall"))
}

#export GroupedYearF as .csv
fileName_GYF = paste(saveDir,site,"_GroupedYearF.csv",sep="")
write.csv(GroupedYearF,fileName_GYF, row.names = FALSE)

#export GroupedYearJ as .csv
fileName_GYJ = paste(saveDir,site,"_GroupedYearJ.csv",sep="")
write.csv(GroupedYearJ,fileName_GYJ, row.names = FALSE)

#export GroupedYearM as .csv
fileName_GYM = paste(saveDir,site,"_GroupedYearM.csv",sep="")
write.csv(GroupedYearM,fileName_GYM, row.names = FALSE)

## GAMs ##

#Social Units#
#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
if (sum(GroupedDayF$Female > 0)){
#females
gamTw = gam(FeHoursProp ~ s(day, bs = 'cc', k = -1), data = GroupedDayF, family = tw, method = "REML")
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
  labs(title = 'Social Units (GAM)', x = 'Day of Year')+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
  theme_get() 
print(vizGG,pages =1)
fig6 =paste(saveDir,site,"GAM1_SocialUnits.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  labs(title = 'Social Units (GAM)')+
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  #l_points(shape = 19, size = 1, alpha = 0.1) + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))+
  theme_classic()
print(vizGG2,pages =1)
fig7 =paste(saveDir,site,"GAM2_SocialUnits.png",sep="")
ggsave(fig7)
}

#Juveniles#
#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
gamTw = gam(JuHoursProp ~ s(day, bs = 'cc', k = -1), data = GroupedDayJ, family = tw, method = "REML")
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
  labs(title = 'Mid-Size (GAM)')+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
print(vizGG,pages =1)
fig6 =paste(saveDir,site,"GAM1_Juveniles.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  labs(title = 'Mid-Size (GAM)')+
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  l_points(shape = 19, size = 1, alpha = 0.1) + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))+
  theme_classic()
print(vizGG2,pages =1)
fig7 =paste(saveDir,site,"GAM2_Juveniles.png",sep="")
ggsave(fig7)

#Males#
#GAMs with appropiate ITS binning

#GAM to identify seasonal pattern
gamTw = gam(MaHoursProp ~ s(day, bs = 'cc', k = -1), data = GroupedDayM, family = tw, method = "REML")
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
  labs(title = 'Males (GAM)', x = 'Day of Year')+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
theme_get() 
fig6 =paste(saveDir,site,"GAM1_Males.png",sep="")
ggsave(fig6)

#second way to plot GAM
vizGG2 = plot(viz, allTerms = T) +
  l_fitLine(colour = "red") + l_rug(mapping = aes(x=x,y=y), alpha=0.8) +
  labs(title = 'Males (GAM)')+
  l_ciLine(mul = 5, colour = "blue", linetype = 2)+
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
print(vizGG2,pages =1)
fig7 =paste(saveDir,site,"GAM2_Males.png",sep="")
ggsave(fig7)

#### Kruskall Wallis Test + Dunn's Test
#Females
kruskal.test(FeHoursProp ~ Season, data = binPresence)
dunnTest(FeHoursProp ~ Season,
         data=binPresence,
         method="bonferroni")
pairwise.wilcox.test(binPresence$FeHoursProp, binPresence$Season,
                     p.adjust.method = "BH")
#Juveniles
kruskal.test(JuHoursProp ~ Season, data = binPresence)
dunnTest(JuHoursProp ~ Season,
         data=binPresence,
         method="bonferroni")
pairwise.wilcox.test(binPresence$JuHoursProp, binPresence$Season,
                     p.adjust.method = "BH")
#Males
kruskal.test(MaHoursProp ~ Season, data = binPresence)
dunnTest(MaHoursProp ~ Season,
         data=binPresence,
         method="bonferroni")
pairwise.wilcox.test(binPresence$MaHoursProp, binPresence$Season,
                     p.adjust.method = "BH")