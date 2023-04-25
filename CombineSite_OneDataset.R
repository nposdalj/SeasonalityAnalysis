## This is script was written to take the data from each site in a region and create one table fro GAM_GEE modeling ##

#load libraries
library(tidyverse)
library(anytime)
library(lubridate)

GDrive = 'G'

#Sites
#Sites = c('CB','PT','QN','AB','KOA','BD','KS') #The GOA and BSAI Sites
Sites = c('HZ','OC','NC','BC','WC','NFC','GS','BP','BS','JAX')
#Sites = c('CA','CCE','CORC','DCPP01C','GI','HOKE','PS1','PS2','QC')
#Sites = c('CSM','Equator','Kauai','King','Kona','LSM','Pagan','Palmyra','PHR','Saipan','Tinian','Wake')

#Regions
#Regions = c('CentralPac')
#region = 'CentralPac'
#Regions = c('CCE')
#region = 'CCE'
#Regions = c("GOA","BSAI") #GOA AND BSAI Region
#region = 'GofAK'
region = 'WAT'
Regions = c("North","South")


fileDir = paste(GDrive,":/My Drive/",region,"_TPWS_metadataReduced/SeasonalityAnalysis/",Sites, sep="")#setting the directory
saveDir = paste(GDrive,":/My Drive/",region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/",sep="")
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
    dplyr::select(tbin, HoursNorm, Effort_Bin, Effort_Sec)
  name = Sites[i]
  modDayBinTAB$Site = name
  modDayBinTAB$ID = i
  if (between(i,1,5)){
    modDayBinTAB$Region = Regions[1]}
  if (between(i,6,9)){
    modDayBinTAB$Region = Regions[2]}
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

#PreAbs
DayTab$HoursNorm[is.na(DayTab$HoursNorm)] = 0
#DayTab$HoursNorm %>% mutate_if(is.numeric, ~1 * (. > 0)) #replaces NA values with 0
DayTab$PreAbs = ifelse(DayTab$HoursNorm>0,1,0)

#Export grouped table as .csv
fileName = paste(saveDir,"AllSitesGrouped_GAMGEE_ROW.csv",sep="")
write.csv(DayTab,fileName, row.names = FALSE)

#### Hourly Data
#Hourly table with each column a different site
HourTab = data.frame(matrix(ncol = 0, nrow=0))

for (i in 1:length(Sites)){
  HourBinTab = read.csv(filename2[i])
  modHourBinTAB = HourBinTab %>%
    dplyr::select(tbin, PreAbs, Effort_Bin, Effort_Sec)
  name = Sites[i]
  names(modHourBinTAB)[names(modHourBinTAB) == "PreAbs"] = name
  if (i == 1){
    HourTab = modHourBinTAB
  }
  HourTab = merge(HourTab, modHourBinTAB, all = TRUE)
}

#Add Julian Day and year
HourTab$tbin = anytime(as.factor(HourTab$tbin))
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
    dplyr::select(tbin, PreAbs, Effort_Bin, Effort_Sec)
  name = Sites[i]
  modHourBinTAB$Site = name
  modHourBinTAB$ID = i
  if (between(i,1,6)){
    modHourBinTAB$Region = Regions[1]}
  if (between(i,7,10)){
    modHourBinTAB$Region = Regions[2]}
  if (i == 1){
    HourTab = modHourBinTAB
  }
  HourTab = merge(HourTab, modHourBinTAB, all = TRUE)
}

#Add Julian Day and year
HourTab$tbin = anytime(as.factor(HourTab$tbin))
HourTab$Julian = format(HourTab$tbin,"%j")
HourTab$Year = format(HourTab$tbin,"%Y")

#Export grouped table as .csv
fileName = paste(saveDir,"AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="")
write.csv(HourTab,fileName, row.names = FALSE)

