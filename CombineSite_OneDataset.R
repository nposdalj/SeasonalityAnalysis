## This is script was written to take the data from each site in a region and create one table fro GAM_GEE modeling ##

#load libraries
library(tidyverse)
library(anytime)
library(lubridate)

Sites = c('CB','PT','QN','AB','KOA','BD','KS') #The GOA and BSAI Sites

fileDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",Sites, sep="")#setting the directory
filename = paste(fileDir,"/",Sites,"_dayData_forGLMR125.csv",sep="")

DayTab = data.frame(matrix(ncol = 0, nrow=0))

for (i in 1:length(Sites)){
  dayBinTAB = read.csv(filename[i])
  modDayBinTAB = dayBinTAB %>%
    dplyr::select(tbin, HoursNorm)
  name = Sites[i]
  names(modDayBinTAB)[names(modDayBinTAB) == "HoursNorm"] = name
  DayTab = merge(DayTab, modDayBinTAB, all = TRUE)
}

DayTab$tbin = anytime(as.factor(DayTab$tbin))
DayTab$tbin = format(as.POSIXct(DayTab$tbin,format='%m/%d/%Y %H:%M:%S'),format='%m/%d/%Y')
DayTab$tbin = as.Date(DayTab$tbin)
DayTab = DayTab[order(DayTab$tbin),]
