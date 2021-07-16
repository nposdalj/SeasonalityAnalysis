## This is script was written to take the data from each site in a region and create one table fro GAM_GEE modeling ##

#load libraries
library(tidyverse)
library(anytime)
library(lubridate)

Sites = c('CB','PT','QN','AB','KOA','BD','KS') #The GOA and BSAI Sites

fileDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",Sites, sep="")#setting the directory
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/")
filename = paste(fileDir,"/",Sites,"_dayData_forGLMR125.csv",sep="")

#Day table with each column a different site
DayTab = data.frame(matrix(ncol = 0, nrow=0))

for (i in 1:length(Sites)){
  dayBinTAB = read.csv(filename[i])
  modDayBinTAB = dayBinTAB %>%
    dplyr::select(tbin, HoursNorm)
  name = Sites[i]
  names(modDayBinTAB)[names(modDayBinTAB) == "HoursNorm"] = name
  DayTab = merge(DayTab, modDayBinTAB, all = TRUE)
}

#Sort table by ascending order by time
DayTab$tbin = anytime(as.factor(DayTab$tbin))
DayTab$tbin = format(as.POSIXct(DayTab$tbin,format='%m/%d/%Y %H:%M:%S'),format='%d/%m/%Y')
DayTab$tbin <- lubridate::dmy(DayTab$tbin)
DayTab = DayTab %>% arrange(ymd(DayTab$tbin))
DayTab$Julian = format(DayTab$tbin,"%j")
DayTab$Year = format(DayTab$tbin,"%Y")

#Export grouped table as .csv
fileName = paste(saveDir,"AllSitesGrouped_GAMGEE_COL.csv",sep="")
write.csv(DayTab,fileName, row.names = FALSE)

#Day table with each row that has a site name
DayTab = data.frame(matrix(ncol = 0, nrow=0))

for (i in 1:length(Sites)){
  dayBinTAB = read.csv(filename[i])
  modDayBinTAB = dayBinTAB %>%
    dplyr::select(tbin, HoursNorm)
  name = Sites[i]
  modDayBinTAB$Site = name
  DayTab = merge(DayTab, modDayBinTAB, all = TRUE)
}

#Sort table by ascending order by time
DayTab$tbin = anytime(as.factor(DayTab$tbin))
DayTab$tbin = format(as.POSIXct(DayTab$tbin,format='%m/%d/%Y %H:%M:%S'),format='%d/%m/%Y')
DayTab$tbin <- lubridate::dmy(DayTab$tbin)
DayTab = DayTab %>% arrange(ymd(DayTab$tbin))
DayTab$Julian = format(DayTab$tbin,"%j")
DayTab$Year = format(DayTab$tbin,"%Y")

#Export grouped table as .csv
fileName = paste(saveDir,"AllSitesGrouped_GAMGEE_ROW.csv",sep="")
write.csv(DayTab,fileName, row.names = FALSE)
