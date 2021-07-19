## This is script was written to take the data from each site in a region and create one table fro GAM_GEE modeling ##

#load libraries
library(tidyverse)
library(anytime)
library(lubridate)

Sites = c('CB','PT','QN','AB','KOA','BD','KS') #The GOA and BSAI Sites

fileDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",Sites, sep="")#setting the directory
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/")
filename = paste(fileDir,"/",Sites,"_dayData_forGLMR125.csv",sep="")
filename2 = paste(fileDir,"/",Sites,"_binData_forGAMGEE.csv",sep="")

#### Day Data
#Day table with each column a different site
DayTab = data.frame(matrix(ncol = 0, nrow=0))

for (i in 1:length(Sites)){
  dayBinTAB = read.csv(filename[i])
  modDayBinTAB = dayBinTAB %>%
    dplyr::select(tbin, HoursNorm)
  name = Sites[i]
  names(modDayBinTAB)[names(modDayBinTAB) == "HoursNorm"] = name
  if (i == 1){
    DayTab = modDayBinTAB
  }
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
  modDayBinTAB$ID = i
  if (i == 1){
    DayTab = modDayBinTAB
  }
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

#### Hourly Data
#Hourly table with each column a different site
HourTab = data.frame(matrix(ncol = 0, nrow=0))

for (i in 1:length(Sites)){
  HourBinTab = read.csv(filename2[i])
  modHourBinTAB = HourBinTab %>%
    dplyr::select(tbin, PreAbs)
  name = Sites[i]
  names(modHourBinTAB)[names(modHourBinTAB) == "PreAbs"] = name
  if (i == 1){
    HourTab = modHourBinTAB
  }
  HourTab = merge(HourTab, modHourBinTAB, all = TRUE)
}

#Sort table by ascending order by time
HourTab$tbin = anytime(as.factor(HourTab$tbin))
HourTab$tbin = format(as.POSIXct(HourTab$tbin,format='%m/%d/%Y %H:%M:%S'),format='%d/%m/%Y')
HourTab$tbin <- lubridate::dmy(HourTab$tbin)
HourTab = HourTab %>% arrange(ymd(HourTab$tbin))
HourTab$Julian = format(HourTab$tbin,"%j")
HourTab$Year = format(HourTab$tbin,"%Y")

#Export grouped table as .csv
fileName = paste(saveDir,"AllSitesGrouped_Binary_GAMGEE_COL.csv",sep="")
write.csv(HourTab,fileName, row.names = FALSE)

#Hourly table with each row that has a site name
HourTab = data.frame(matrix(ncol = 0, nrow=0))

for (i in 1:length(Sites)){
  HourBinTab = read.csv(filename2[i])
  modHourBinTAB = HourBinTab %>%
    dplyr::select(tbin, PreAbs)
  name = Sites[i]
  modHourBinTAB$Site = name
  modHourBinTAB$ID = i
  if (i == 1){
    HourTab = modHourBinTAB
  }
  HourTab = merge(HourTab, modHourBinTAB, all = TRUE)
}

#Sort table by ascending order by time
HourTab$tbin = anytime(as.factor(HourTab$tbin))
HourTab$tbin = format(as.POSIXct(HourTab$tbin,format='%m/%d/%Y %H:%M:%S'),format='%d/%m/%Y')
HourTab$tbin <- lubridate::dmy(HourTab$tbin)
HourTab = HourTab %>% arrange(ymd(HourTab$tbin))
HourTab$Julian = format(HourTab$tbin,"%j")
HourTab$Year = format(HourTab$tbin,"%Y")

#Export grouped table as .csv
fileName = paste(saveDir,"AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="")
write.csv(HourTab,fileName, row.names = FALSE)

