#load libraries
library(ncdf4)
library(sp)
library(rgdal)
library(httr)
library(sf)
library(dplyr)
library(lubridate)
library(raster)
library(rgeos)
library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")
library(tidyverse)
library(mgcv)
library(tweedie)
library(anytime)
library(mgcViz)
library(zoo)
library(car)

#load functions
source('C:/Users/DAM1/Documents/Github/SeasonalityAnalysis/GetChlA.R')
source('C:/Users/DAM1/Documents/Github/SeasonalityAnalysis/GetSST.R')
source('C:/Users/DAM1/Documents/Github/SeasonalityAnalysis/GetAVISO.R')

#increasing memory limit
memory.limit(size=300000)

# User Defined sections
#define the lat and long of interest
#df1 = data.frame("lat" = c(19.29, 19.2, 19.2467, 19.2467), "long" = c(166.69, 166.69, 166.74, 166.64)) #Wake
#df1 = data.frame("lat" = c(15.36, 15.27, 15.3186, 15.3186), "long" = c(145.46, 145.46, 145.51, 145.41)) #Saipan
df1 = data.frame("lat" = c(15.08, 14.99, 15.0387, 15.0387), "long" = c(145.75, 145.75, 145.8, 145.7))#Tinian
#df1 = data.frame("lat" = c(29.185, 29.1, 29.1410, 29.1410), "long" = c(-118.26, -118.26, -118.31, -118.21))#GI
#df1 = data.frame("lat" = c(31.79, 31.7, 31.747, 31.747), "long" = c(-121.38, -121.38, -121.43, -121.33))#CORC
#df1 = data.frame("lat" = c(52.4, 52.31, 52.3547, 52.3547), "long" = c(-175.635, -175.635, -175.71, -175.56))#BD
#df1 = data.frame("lat" = c(47.54, 47.45, 47.4936, 47.4936), "long" = c(-125.378, -125.378, -125.44, -125.31))#QC
#df1 = data.frame("lat" = c(56.385, 56.295, 56.34, 56.34), "long" = c(-145.182, -145.182, -145.26, -145.1))#QN
#df1 = data.frame("lat" = c(56.29, 56.2, 56.2434, 56.2434), "long" = c(-142.75, -142.75, -142.83, -142.67))#PT
#df1 = data.frame("lat" = c(58.71, 58.62, 58.6668, 58.6668), "long" = c(-148.0034, -148.0034, -148.12, -147.94))#CB

#define the start and end of the data 
startTime = "2011-04-13" #this should be formatted like this: 2010-03-05
endTime = "2019-05-12"

#ITS
#ITS = 22 #Wake
#ITSF = 9 #Wake Females
#ITSJ = 7 #Wake Juveniles
#ITS = 4 #Saipan
#ITSF = 6 #Saipan Females
#ITSM = 1 #Saipan Males
ITS = 2 #Tinian
ITSJ = 7 #Tinian Juveniles
ITSM = 1 #Tinian Males

#loading the environmental data
envDir = paste("O:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory

#loading sperm whale data
site = 'TIN'
saveDir = paste("O:/My Drive/CentralPac_TPWS_metadataReduced/Tinian/Seasonality/")#setting the directory

#load data from StatisicalAnalysis_All
filenameStatAll = paste(saveDir,site,"_Day.csv",sep="")#_Day.csv created from StatisticalAnalysis_All.R
DayData = read.csv(filenameStatAll) #load files as data frame
DayTable = DayData %>%
  dplyr::select(tbin, Count_Click, Count_Bin, HoursProp, HoursNorm)
names(DayTable)[1] = "time"
DayTable$time = as.Date(DayTable$time)#converting time from character to date

#load sex specific data
filename_sex = paste(saveDir,site,"_binPresence.csv",sep="")
SexDayData = read.csv(filename_sex) #load files as data frame
SexDayData = SexDayData %>%
  dplyr::select(tbin, FemaleHoursNorm, MaleHoursNorm, JuvenileHoursNorm)
names(SexDayData)[1] = "time"

SexDayData$time = anytime(as.factor(SexDayData$time))
SexDayData$time = as.Date(SexDayData$time)

#Load chlorophyll
GetChla(envDir)

#Load SST
GetSST(envDir)
#############################
LoadAVISO(envDir)
GetSSH(AVISO)
GetDEN(AVISO)
GetSAL(AVISO)
GetTEMP(AVISO)
GetEKE(AVISO)

#Visualize Eddy Kinetic Energy
hist(EKE$EKE_cm)
plot(EKE)

#merge the data sets in chunks because of memory issues
tab <- left_join(DayTable, SST4, by = "time") %>%
  left_join(., ChlAdf, by = "time") %>%
  left_join(., SSHdf, by = "time") %>%
  left_join(., DENdf, by = "time") %>%
  left_join(., SALdf, by = "time") %>%
  left_join(., TEMPdf, by = "time") %>%
  left_join(., EKE, by = "time")

#replacing SST Nan values from satellite with resTEMP (model data)
tab$mean_SST <- ifelse(is.na(tab$mean_SST), tab$resTEMP, tab$mean_SST)
tab = tab[complete.cases(tab[ , 2:4]),]#remove any rows with lat or long as na

#Join sex specific presence in daily hours information
tab <- left_join(tab, SexDayData, by = 'time')

#Change column names to look good on plots
#TabBinned_Grouped = TabBinned_Grouped %>% 
 # dplyr::rename(
 #   'Julian Day' = 'Julian',
 # )
#TabBinned_Grouped = TabBinned_Grouped %>% 
  #dplyr::rename(
   # 'Mean SST (C)' = 'mean_SST',
  #)
#TabBinned_Grouped = TabBinned_Grouped %>% 
  #dplyr::rename(
  #  'Standard Deviation of SST (C)' = 'SD_SST',
  #)
#TabBinned_Grouped = TabBinned_Grouped %>% 
  #dplyr::rename(
   # 'Chlorophyll (mg m^-3)' = 'resChlA',
 # )
#TabBinned_Grouped = TabBinned_Grouped %>% 
  #dplyr::rename(
  #  'SSH (m)' = 'resSSH',
 # )
#TabBinned_Grouped = TabBinned_Grouped %>% 
  #dplyr::rename(
    #'Mixed Layer Thickness (m)' = 'resDEN',
 # )
#TabBinned_Grouped = TabBinned_Grouped %>% 
  #dplyr::rename(
  #  'Salinity (ppt)' = 'resSAL',
  #)
#TabBinned_Grouped = TabBinned_Grouped %>% 
  #dplyr::rename(
  #  'EKE (cm^2 s^-2)' = 'EKE_cm',
 # )
  
#Group by ITS for the general sperm whale model
startDate = tab$time[1]
endDate = tab$time[nrow(tab)]
timeseries = data.frame(date=seq(startDate, endDate, by="days"))
div = floor(nrow(timeseries)/ITS)
ITSgroups = rep(1:div, times=1, each=ITS)
divdiff = nrow(timeseries) - length(ITSgroups)
last = tail(ITSgroups, n = 1)
lastVec = rep(last,each = divdiff)
timeseries$groups = c(ITSgroups,lastVec)
TabBinned = left_join(tab,timeseries,by = c("time" = "date"))
TabBinned_Grouped = aggregate(TabBinned[, c(1:length(TabBinned))], list(TabBinned$groups), mean, na.rm = TRUE)
TabBinned_Grouped$Julian = as.numeric(format(TabBinned_Grouped$time,"%j"))
TabBinned_Grouped$Year = as.numeric(format(TabBinned_Grouped$time,"%Y"))

#Group by ITS for Female
if (exists("ITSF")){
  if (ITSF > 1){
    ITSgroupsF = rep(1:(floor(nrow(timeseries)/ITSF)), times=1, each=ITSF)
    timeseriesF = data.frame(date=seq(startDate, endDate, by="days"))
    div = floor(nrow(timeseriesF)/ITSF)
    ITSgroups = rep(1:div, times=1, each=ITSF)
    divdiff = nrow(timeseriesF) - length(ITSgroups)
    last = tail(ITSgroups, n = 1)
    lastVec = rep(last,each = divdiff)
    timeseriesF$groups = c(ITSgroups,lastVec)
    TabBinnedF = left_join(tab,timeseriesF,by = c("time" = "date"))
    TabBinned_GroupedF = aggregate(TabBinnedF[, c(1:length(TabBinnedF))], list(TabBinnedF$groups), mean, na.rm = TRUE)
    TabBinned_GroupedF$Julian = as.numeric(format(TabBinned_GroupedF$time,"%j"))
    TabBinned_GroupedF$Year = as.numeric(format(TabBinned_GroupedF$time,"%Y"))
  } else {
    TabBinned_GroupedF = tab
    TabBinned_GroupedF$Julian = as.numeric(format(TabBinned_GroupedF$time,"%j"))
    TabBinned_GroupedF$Year = as.numeric(format(TabBinned_GroupedF$time,"%Y"))
  }
}

#Group by ITS for Juvenile
if (exists("ITSJ")){
  if (ITSJ > 1){
    ITSgroupsJ = rep(1:(floor(nrow(timeseries)/ITSJ)), times=1, each=ITSJ)
    timeseriesJ = data.frame(date=seq(startDate, endDate, by="days"))
    div = floor(nrow(timeseriesJ)/ITSJ)
    ITSgroups = rep(1:div, times=1, each=ITSJ)
    divdiff = nrow(timeseriesJ) - length(ITSgroups)
    last = tail(ITSgroups, n = 1)
    lastVec = rep(last,each = divdiff)
    timeseriesJ$groups = c(ITSgroups,lastVec)
    TabBinnedJ = left_join(tab,timeseriesJ,by = c("time" = "date"))
    TabBinned_GroupedJ = aggregate(TabBinnedJ[, c(1:length(TabBinnedJ))], list(TabBinnedJ$groups), mean, na.rm = TRUE)
    TabBinned_GroupedJ$Julian = as.numeric(format(TabBinned_GroupedJ$time,"%j"))
    TabBinned_GroupedJ$Year = as.numeric(format(TabBinned_GroupedJ$time,"%Y"))
  } else {
    TabBinned_GroupedJ = tab
    TabBinned_GroupedJ$Julian = as.numeric(format(TabBinned_GroupedJ$time,"%j"))
    TabBinned_GroupedJ$Year = as.numeric(format(TabBinned_GroupedJ$time,"%Y"))
  }
}

#Group by ITS for Male
if (exists("ITSM")){
  if (ITSM > 1){
    ITSgroupsM = rep(1:(floor(nrow(timeseries)/ITSM)), times=1, each=ITSM)
    timeseriesM = data.frame(date=seq(startDate, endDate, by="days"))
    div = floor(nrow(timeseriesM)/ITSM)
    ITSgroups = rep(1:div, times=1, each=ITSM)
    divdiff = nrow(timeseriesM) - length(ITSgroups)
    last = tail(ITSgroups, n = 1)
    lastVec = rep(last,each = divdiff)
    timeseriesM$groups = c(ITSgroups,lastVec)
    TabBinnedM = left_join(tab,timeseriesM,by = c("time" = "date"))
    TabBinned_GroupedM = aggregate(TabBinnedM[, c(1:length(TabBinnedM))], list(TabBinnedM$groups), mean, na.rm = TRUE)
    TabBinned_GroupedM$Julian = as.numeric(format(TabBinned_GroupedM$time,"%j"))
    TabBinned_GroupedM$Year = as.numeric(format(TabBinned_GroupedM$time,"%Y"))
  } else {
    TabBinned_GroupedM = tab
    TabBinned_GroupedM$Julian = as.numeric(format(TabBinned_GroupedM$time,"%j"))
    TabBinned_GroupedM$Year = as.numeric(format(TabBinned_GroupedM$time,"%Y"))
  }
}


#run GAM
#If you get either of these errors, it might be because there are too many NAs in the variable, try interpolating:
#1) Error in if (theta > 0) (b + a * exp(-theta))/(1 + exp(-theta)) else (b *  : 
#missing value where TRUE/FALSE needed
#2) Error in gam.reparam(UrS, sp, grderiv) : 
#NA/NaN/Inf in foreign function call (arg 3)
#Interpolating using na.approx:
#TabBinned_GroupedJ$resChlA = na.approx(TabBinned_GroupedJ$resChlA, na.rm = FALSE)

#Test how each covariate should be used (linear, smooth, as.factor())
#empty model for comparison
GAM_empty = gam(HoursNorm ~ 1, data = TabBinned_Grouped, family = tw, method = "REML")

#Julian Day
GAM_01a = gam(HoursNorm ~ Julian, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_01b = gam(HoursNorm ~ s(Julian, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model01 = c('empty','01a','01b')
AIC01 = c(AIC(GAM_empty),AIC(GAM_01a),AIC(GAM_01b))
data.frame(rbind(model01,AIC01))
#Saipan
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   2143.75884213975 2130.74080097241 2103.62871708431
#Julian as a "cc" smooth

#Tinian
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   1634.79054479956 1623.26560883975 1616.30134957535
#Julian as a 'cc' smooth

#Wake
#Julian is always a "cc" smooth

#Chlorophyll
GAM_02a = gam(HoursNorm ~ resChlA, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_02b = gam(HoursNorm ~ s(resChlA, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_02c = gam(HoursNorm ~ s(resChlA, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model02 = c('empty','02a','02b','02c')
AIC02 = c(AIC(GAM_empty),AIC(GAM_02a),AIC(GAM_02b),AIC(GAM_02c))
data.frame(rbind(model02,AIC02))
#Saipan
#X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   2141.48369945917 1214.96924100947 1213.11010200871 1213.17495326208
#Chlorophyll as a 'cr' smooth

#Tinian
#X1               X2               X3
#model02            empty              02a              02b
#AIC02   1634.79054479956 406.874091957057 405.618387661288
#Chlorophyll as a smooth

#Wake
#X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   342.418503970956 318.213389695321 318.574762805002 318.613374773736
#Chlorophyll as linear

#EKE
GAM_03a = gam(HoursNorm ~ EKE_cm, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_03b = gam(HoursNorm ~ s(EKE_cm, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_03c = gam(HoursNorm ~ s(EKE_cm, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model03 = c('empty','03a','03b', '03c')
AIC03 = c(AIC(GAM_empty),AIC(GAM_03a),AIC(GAM_03b),AIC(GAM_03c))
data.frame(rbind(model03,AIC03))
#Saipan
#                      X1               X2               X3               X4
#model03            empty              03a              03b              03c
#AIC03   2141.48369945917 2143.42543933207 2143.42912412885 2143.42741443162
#EKE as linear

#Tinian
#X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   1634.79054479956 406.874091957057 406.503032970221 406.645195238531
#EKE as linear

#Wake
#X1               X2               X3               X4
#model03            empty              03a              03b              03c
#AIC03   342.418503970956 332.704565927828 332.704981498064 332.704994948509
#EKE as linear

#Salinity
GAM_04a = gam(HoursNorm ~ resSAL, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_04b = gam(HoursNorm ~ s(resSAL, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_04c = gam(HoursNorm ~ s(resSAL, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model04 = c('empty','04a','04b','04c')
AIC04 = c(AIC(GAM_empty),AIC(GAM_04a),AIC(GAM_04b),AIC(GAM_04c))
data.frame(rbind(model04,AIC04))
#Saipan
#                      X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   2141.48369945917 2141.98403859539 2141.80059043796 2141.83559361882
#salinity as a 'cr' smooth

#Tinian
#                     X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   1634.79054479956 1617.06202675983 1617.06560017288 1617.06606751368
#salinity as linear

#Wake
#                      X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   342.418503970956 340.984401015384 338.536236907081 338.594316591802
#Salinity as "cr" smooth

#Mean SST
GAM_05a = gam(HoursNorm ~ mean_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_05b = gam(HoursNorm ~ s(mean_SST, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_05c = gam(HoursNorm ~ s(mean_SST, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model05 = c('empty','05a','05b','05c')
AIC05 = c(AIC(GAM_empty),AIC(GAM_05a),AIC(GAM_05b),AIC(GAM_05c))
data.frame(rbind(model05,AIC05))
#Saipan
#                      X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   2141.48369945917 2136.45760434501 2131.27818522867 2131.99464816362
#mean SST as a 'cr' smooth

#Tinian
#                      X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   1634.79054479956 1632.44148060028 1632.44530087046 1632.44227753014
#mean SST as 'cr' smooth 

#Wake
#X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   342.418503970956 338.219635726337 338.221139921545 338.220324276489
#mean SST as linear

#SSH
GAM_06a = gam(HoursNorm ~ resSSH, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_06b = gam(HoursNorm ~ s(resSSH, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_06c = gam(HoursNorm ~ s(resSSH, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model06 = c('empty','06a','06b','06c')
AIC06 = c(AIC(GAM_empty),AIC(GAM_06a),AIC(GAM_06b),AIC(GAM_06c))
data.frame(rbind(model06,AIC06))
#Saipan
#                      X1               X2               X3               X4
#model06            empty              06a              06b              06c
#AIC06   2141.48369945917 2136.90561489243 2130.15296421187 2131.91936153117
#SSH as a 'cr' smooth

#Tinian
#                      X1               X2               X3              X4
#model06            empty              06a              06b             06c
#AIC06   1634.79054479956 1635.73759027985 1634.93938324326 1635.3949923932
#SSH as a 'cr' smooth

#Wake
#X1               X2              X3               X4
#model06            empty              06a             06b              06c
#AIC06   342.418503970956 339.984798947495 337.93286226822 337.624468243133
#SSH as a normal smooth

#Density
GAM_07a = gam(HoursNorm ~ resDEN, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_07b = gam(HoursNorm ~ s(resDEN, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_07c = gam(HoursNorm ~ s(resDEN, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model07 = c('empty','07a','07b','07c')
AIC07 = c(AIC(GAM_empty),AIC(GAM_07a),AIC(GAM_07b),AIC(GAM_07c))
data.frame(rbind(model07,AIC07))
#Saipan
#                      X1               X2               X3               X4
#model07            empty              07a              07b              07c
#AIC07   2141.48369945917 2139.77816579722 2141.49509148073 2140.07866063997
#density as linear

#Tinian
#                      X1               X2               X3               X4
#model07            empty              07a              07b              07c
#AIC07   1634.79054479956 1629.24706552712 1628.45393301092 1629.22563851813
#density as a 'cr' smooth

#Wake
#X1               X2               X3               X4
#model07            empty              07a              07b              07c
#AIC07   342.418503970956 343.931005901443 342.977786180689 344.064351292256
#density as "cr" smooth

#Standard deviation of SST
GAM_08a = gam(HoursNorm ~ SD_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_08b = gam(HoursNorm ~ s(SD_SST, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_08c = gam(HoursNorm ~ s(SD_SST, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model08 = c('empty','08a','08b','08c')
AIC08 = c(AIC(GAM_empty),AIC(GAM_08a),AIC(GAM_08b),AIC(GAM_08c))
data.frame(rbind(model08,AIC08))
#Saipan
#                      X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   2141.48369945917 1457.77690764797 1456.05186214008 1457.78489342236
#std dev of SST as a 'cr' smooth

#Tinian
#                      X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   1634.79054479956 941.982727021977 942.329997658581 941.986021503707
#std dev of SST as linear

#Wake
#X1              X2               X3               X4
#model08            empty             08a              08b              08c
#AIC08   342.418503970956 300.37593926682 298.573419287569 300.376783521815
#std dev of SST as "cr" smooth 

#Test which covariates we should keep
#Round 1
#Initial model
Full1 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, k=-1)+EKE_cm+
              s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
              resDEN +SD_SST
              , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ s(resChlA, k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
C = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)++EKE_cm+
           s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
           resDEN +SD_SST
         , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, k=-1)+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, k=-1)+EKE_cm+
          mean_SST +s(resSSH, k=-1)+
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, k=-1)+EKE_cm+
           s(resSAL,k=-1) + s(resSSH, k=-1)+
           resDEN +SD_SST
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST +
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, k=-1)+EKE_cm+
            s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
            resDEN 
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR1 = c('Full1','J','C','E','S','MT','H','D','SDT')
AICR1 = c(AIC(Full1),AIC(J),AIC(C),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR1,AICR1))
#Saipan
#                      X1               X2               X3               X4               X5
#modelR1            Full1                J                C                E                S
#AICR1   1074.92346880652 1099.44748857908 1425.95668084402 1075.36050706271 1072.61740230551
#X6              X7               X8               X9
#modelR1               MT               H                D              SDT
#AICR1   1074.22394361901 1084.4955849055 1075.38254469558 1179.54591701796
#std dev of SST removed

#Tinian
#                      X1               X2               X3               X4
#modelR1            Full1                J                C                E
#AICR1   384.505165917479 387.538727116988 922.043580227863 388.795226841835
#X5               X6               X7               X8
#modelR1                S               MT                H                D
#AICR1   386.862320761001 382.726979524462 385.912192947752 386.598563403263
#X9
#modelR1              SDT
#AICR1   395.566731253434
#Chlorophyll removed

#Wake
#X1              X2               X3               X4
#modelR1            Full1               J                C                E
#AICR1   251.671545884863 265.07466465574 276.401097726524 259.633713158517
#X5               X6               X7               X8
#modelR1                S               MT                H                D
#AICR1   263.322382631849 249.897600993043 255.776793612547 248.931834405921
#X9
#modelR1              SDT
#AICR1   282.921676502635
#std dev of SST removed

#Round 2
Full2 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
              s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
              resDEN +SD_SST
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ EKE_cm+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          mean_SST +s(resSSH, k=-1)+
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
           s(resSAL,k=-1) + s(resSSH, k=-1)+
           resDEN +SD_SST
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST +
          resDEN +SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          SD_SST
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
            s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
            resDEN 
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR2 = c('Full2','J','E','S','MT','H','D','SDT')
AICR2 = c(AIC(Full2),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR2,AICR2))
#Saipan
#                      X1               X2               X3               X4               X5
#modelR2            Full2                J                C                E                S
#AICR2   1179.54591701796 1203.68115687539 2094.28743583581 1178.61957724007 1176.94244251436
#X6               X7               X8
#modelR2               MT                H                D
#AICR2   1179.99647563686 1190.44263901117 1179.80940046926
#chlorophyll removed

#Tinian
#                      X1               X2               X3               X4
#modelR2            Full2                J                E                S
#AICR2   922.043580227863 935.898525770589 924.915361132531 920.541663629191
#X5               X6               X7               X8
#modelR2               MT                H                D              SDT
#AICR2   925.149902506105 924.021947974003 923.889096554455 1601.40562663682
#std dev of SST removed

#Wake
#X1               X2               X3              X4
#modelR2            Full2                J                C               E
#AICR2   282.921676502635 294.567213743167 305.858760371536 291.10470217936
#X5               X6               X7               X8
#modelR2                S               MT                H                D
#AICR2   294.824575710083 282.291931048225 287.009710075041 248.931834405921
#Chlorophyll removed

#Round 3
Full3 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
              s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
              resDEN 
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ EKE_cm+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          resDEN 
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)+
          resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          mean_SST +s(resSSH, k=-1)+
          resDEN 
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
           s(resSAL,k=-1) + s(resSSH, k=-1)+
           resDEN 
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST +
          resDEN 
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR3 = c('Full3','J','E','S','MT','H','D')
AICR3 = c(AIC(Full3),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR3,AICR3))
#Saipan
#                      X1               X2               X3               X4               X5
#modelR3            Full3                J                E                S               MT
#AICR3   2094.28743583581 2125.64438030823 2093.67229772668 2092.33947845392 2099.21516276714
#X6               X7
#modelR3                H                D
#AICR3   2098.54437778115 2092.43166977793
#mean SST removed, keep Julian despite high AIC

#Tinian
#                      X1               X2               X3               X4
#modelR3            Full3                J                E                S
#AICR3   1601.40562663682 1616.04878322308 1604.50643968103 1602.71499674374
#X5               X6               X7
#modelR3              MT                H                D
#AICR3   1602.1128395232 1602.79575326853 1607.61059972203
#density removed, keeping Julian despite high AIC

#Wake
#X1               X2               X3               X4
#modelR3            Full3                J                E                S
#AICR3   305.858760371536 312.425850162277 311.018375427713 316.020082102712
#X5               X6               X7
#modelR3               MT                H                D
#AICR3   306.683640205421 309.717116072572 274.810202653348
#Salinity removed

#Round 4
Full4 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
              s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ EKE_cm+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1) 
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL,k=-1) + mean_SST +s(resSSH, k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          mean_SST +s(resSSH, k=-1) 
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
           s(resSAL,k=-1) + s(resSSH, k=-1) 
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + mean_SST 
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR4 = c('Full4','J','E','S','MT','H')
AICR4 = c(AIC(Full4),AIC(J),AIC(E),AIC(S),AIC(MT), AIC(H))
data.frame(rbind(modelR4,AICR4))
#Saipan
#                      X1              X2               X3               X4               X5
#modelR4            Full4               J                E                S                H
#AICR4   2099.21516276714 2133.6442884902 2097.71945023963 2097.32620785073 2105.98181795106
#X6
#modelR4                D
#AICR4   2097.27847400687
#SSH removed 

#Tinian
#                      X1               X2               X3              X4
#modelR4            Full4                J                E               S
#AICR4   1607.61059972203 1615.83298282944 1608.91184657886 1612.0549582165
#X5              X6
#modelR4               MT               H
#AICR4   1607.05642493893 1608.4761722472
#salinity removed

#Wake
#                      X1               X2               X3               X4
#modelR4            Full4                J                E               MT
#AICR4   316.020082102712 322.458074308942 330.839704856908 318.712970168162
#X5               X6
#modelR4                H                D
#AICR4   315.063034763153 285.032746767053
#EKE removed

#Round 5
Full5 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
              mean_SST +s(resSSH, k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ EKE_cm+
          mean_SST +s(resSSH, k=-1) 
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+
          mean_SST +s(resSSH, k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
           s(resSSH, k=-1) 
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          mean_SST 
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR5 = c('Full5','J','MT','E','H')
AICR5 = c(AIC(Full5),AIC(J),AIC(MT), AIC(E),AIC(H))
data.frame(rbind(modelR5,AICR5))
#Saipan
#X1               X2               X3               X4               X5
#modelR5            Full5                J               MT                S                D
#AICR5   2105.98181795106 2143.28572286566 2099.21516276714 2105.04680159995 2104.87818890922
#keep remaining variables

#Tinian
#                     X1               X2               X3               X4
#modelR5           Full5                J               MT                E
#AICR5   1612.0549582165 1632.77595004857 1614.41109779104 1614.85453585155
#X5
#modelR5                H
#AICR5   1612.04281512021
#keep remaining environmental vars for full model

#Wake
#X1               X2               X3               X4
#modelR5            Full5                J               MT                H
#AICR5   330.839704856908 332.630113616472 334.171518038998 331.847779269344
#X5
#modelR5                D
#AICR5   295.826168046641
#keep remaining variables

#Final full model
FinalGAM = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
                 mean_SST +s(resSSH, k=-1)
               , data = TabBinned_Grouped, family = tw, method = "REML")
summary(FinalGAM)
viz = getViz(FinalGAM)
vizGG = plot(viz,allTerms = T) +
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
print(vizGG,pages =1)
fig1 =paste(saveDir,site,"_ENV_GAM.png",sep="")
ggsave(fig1)

#running Sex Specific GAM
  ##Females
#Test how each covariate should be used (linear, smooth, as.factor())
#empty model for comparison
GAM_empty = gam(FemaleHoursNorm ~ 1, data = TabBinned_GroupedF, family = tw, method = "REML")

#Julian
GAM_01a = gam(FemaleHoursNorm ~ Julian, data = TabBinned_GroupedF, family = tw, method = "REML")
GAM_01b = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k =-1), data = TabBinned_GroupedF, family = tw, method = "REML")

model01 = c('empty','01a','01b')
AIC01 = c(AIC(GAM_empty),AIC(GAM_01a),AIC(GAM_01b))
data.frame(rbind(model01,AIC01))
#Saipan
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   2665.37688613595 2658.91448130588 2633.63404062137
#Julian as a 'cc' smooth

#Wake
#Julian is always a smooth

#Chlorophyll
GAM_02a = gam(FemaleHoursNorm ~ resChlA, data = TabBinned_GroupedF, family = tw, method = "REML")
GAM_02b = gam(FemaleHoursNorm ~ s(resChlA, bs = "cr", k =-1), data = TabBinned_GroupedF, family = tw, method = "REML")
GAM_02c = gam(FemaleHoursNorm ~ s(resChlA, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model02 = c('empty','02a','02b','02c')
AIC02 = c(AIC(GAM_empty),AIC(GAM_02a),AIC(GAM_02b),AIC(GAM_02c))
data.frame(rbind(model02,AIC02))
#Saipan
#                      X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   2141.48369945917 1375.88544364518 1376.06454424454 1526.95091237974
#Chlorophyll as linear

#Wake
#X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   342.418503970956 552.853095505039 553.370769822935 306.198743533215
#Chlorophyll as a normal smooth

#EKE
GAM_03a = gam(FemaleHoursNorm ~ EKE_cm, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_03b = gam(FemaleHoursNorm ~ s(EKE_cm, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_03c = gam(FemaleHoursNorm ~ s(EKE_cm, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model03 = c('empty','03a','03b','03c')
AIC03 = c(AIC(GAM_empty),AIC(GAM_03a),AIC(GAM_03b),AIC(GAM_03c))
data.frame(rbind(model03,AIC03))
#Saipan
#                      X1              X2               X3               X4
#model03            empty             03a              03b              03c
#AIC03   2141.48369945917 2665.9152563983 2666.11409520255 2666.17020430492
#EKE as linear

#Wake
#X1               X2               X3               X4
#model03            empty              03a              03b              03c
#AIC03   342.418503970956 322.013220493546 322.013766797164 322.014377920677
#EKE as linear

#Salinity
GAM_04a = gam(FemaleHoursNorm ~ resSAL, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_04b = gam(FemaleHoursNorm ~ s(resSAL, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_04c = gam(FemaleHoursNorm ~ s(resSAL, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model04 = c('empty','04a','04b','04c')
AIC04 = c(AIC(GAM_empty),AIC(GAM_04a),AIC(GAM_04b),AIC(GAM_04c))
data.frame(rbind(model04,AIC04))
#Saipan
#                      X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   2141.48369945917 2665.09446554882 2665.09912189751 2665.09741643719
#salinity as a normal smooth

#Wake
#X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   342.418503970956 316.925946704085 316.468756083078 316.467172651007
#Salinity as a normal smooth

#Mean SST
GAM_05a = gam(FemaleHoursNorm ~ mean_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_05b = gam(FemaleHoursNorm ~ s(mean_SST, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_05c = gam(FemaleHoursNorm ~ s(mean_SST, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model05 = c('empty','05a','05b','05c')
AIC05 = c(AIC(GAM_empty),AIC(GAM_05a),AIC(GAM_05b),AIC(GAM_05c))
data.frame(rbind(model05,AIC05))
#Saipan
#                      X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   2141.48369945917 2661.14560511157 2661.15210797022 2661.15350515159
#Mean SST as a normal smooth

#Wake
#X1               X2              X3               X4
#model05            empty              05a             05b              05c
#AIC05   342.418503970956 325.407383810289 317.50494400982 317.379127234032
#Mean SST as a normal smooth

#SSH
GAM_06a = gam(FemaleHoursNorm ~ resSSH, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_06b = gam(FemaleHoursNorm ~ s(resSSH, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_06c = gam(FemaleHoursNorm ~ s(resSSH, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model06 = c('empty','06a','06b','06c')
AIC06 = c(AIC(GAM_empty),AIC(GAM_06a),AIC(GAM_06b),AIC(GAM_06c))
data.frame(rbind(model06,AIC06))
#Saipan
#                      X1               X2               X3               X4
#model06            empty              06a              06b              06c
#AIC06   2141.48369945917 2660.44911654165 2651.43967180873 2651.14127853984
#SSH as a normal smooth

#Wake
#                      X1               X2               X3               X4
#model06            empty              06a              06b              06c
#AIC06   342.418503970956 332.636187968703 326.278595735254 326.267615888791
#SSH as a normal smooth

#Density
GAM_07a = gam(FemaleHoursNorm ~ resDEN, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_07b = gam(FemaleHoursNorm ~ s(resDEN, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_07c = gam(FemaleHoursNorm ~ s(resDEN, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model07 = c('empty','07a','07b','07c')
AIC07 = c(AIC(GAM_empty),AIC(GAM_07a),AIC(GAM_07b),AIC(GAM_07c))
data.frame(rbind(model07,AIC07))
#Saipan
#                      X1              X2               X3               X4
#model07            empty             07a              07b              07c
#AIC07   2141.48369945917 2658.5385365849 2658.68048855064 2658.63697297238
#Density as linear

#Wake
#X1               X2               X3               X4
#model07            empty              07a              07b              07c
#AIC07   342.418503970956 330.046034889351 330.311408687101 330.359112573299
#Density as linear

#Standard deviation of SST
GAM_08a = gam(FemaleHoursNorm ~ SD_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_08b = gam(FemaleHoursNorm ~ s(SD_SST, bs = "cr", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")
GAM_08c = gam(FemaleHoursNorm ~ s(SD_SST, k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model08 = c('empty','08a','08b','08c')
AIC08 = c(AIC(GAM_empty),AIC(GAM_08a),AIC(GAM_08b),AIC(GAM_08c))
data.frame(rbind(model08,AIC08))
#Saipan
#                      X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   2141.48369945917 1728.41431180364 1728.41894537372 1728.42036375173
#std dev of SST as a 'cr' smooth

#Wake
#X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   342.418503970956 302.614778713009 302.446890123018 302.495229772987
#std dev of SST as a cr smooth

#Test which covariates we should keep for Female Specific GAM
#Round 1
#Initial model
Full1 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+EKE_cm+
              s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ resChlA+EKE_cm+
           s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
         , data = TabBinned_Grouped, family = tw, method = "REML")
C = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+EKE_cm+
          s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+EKE_cm+
           s(resSAL, k=-1) + s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+EKE_cm+
            s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR1 = c('Full1','J','C','E','S','MT','H','D','SDT')
AICR1 = c(AIC(Full1),AIC(J),AIC(C),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR1,AICR1))
#Saipan
#                     X1               X2               X3              X4               X5
#modelR1           Full1                J                C               E                S
#AICR1   1335.7612010297 1354.34133850303 1707.67131139023 1333.7611456684 1333.89913476741
#X6               X7               X8               X9
#modelR1               MT                H                D              SDT
#AICR1   1333.72118706082 1339.47402313534 1335.24786560955 1507.38984825146
#remove chlorophyll

#Wake
#X1               X2              X3               X4
#modelR1            Full1                J               C                E
#AICR1   242.033040895924 253.267392343263 277.64557938092 254.876012601097
#X5              X6               X7               X8
#modelR1                S              MT                H                D
#AICR1   253.998739616413 241.61445530537 240.060981624594 240.706846162488
#X9
#modelR1            SDT
#AICR1   271.5183301578
#remove chlorophyll

#Round 2
Full2 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
              s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
           s(resSAL, k=-1) + s(resSSH, k=-1)+resDEN+s(SD_SST, bs="cr", k=-1)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +resDEN+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+s(SD_SST, bs="cr", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
            s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR2 = c('Full2','J','E','S','MT','H','D','SDT')
AICR2 = c(AIC(Full2),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR2,AICR2))
#Saipan 
#                      X1              X2               X3               X4               X5
#modelR2            Full2               J                E                S               MT
#AICR2   1707.67131139023 1724.9201707516 1704.48561078599 1706.10618471288 1703.45131635063
#X6               X7               X8
#modelR2                H                D              SDT
#AICR2   1709.86865476309 1706.23837538322 2628.04810164109
#remove std dev SST

#Wake
#X1               X2               X3               X4
#modelR2           Full2                J                E                S
#AICR2   277.64557938092 277.645206045198 282.797109140055 282.080287469881
#X5               X6               X7               X8
#modelR2               MT                H                D              SDT
#AICR2   280.668741424683 279.042117417182 276.138284284723 301.775364697152
#remove std dev SST

#Round 3
Full3 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
              s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
           s(resSAL, k=-1) + s(resSSH, k=-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +s(resSSH, k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR3 = c('Full3','J','E','S','MT','H','D')
AICR3 = c(AIC(Full3),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR3,AICR3))
#Saipan
#                      X1               X2               X3               X4              X5
#modelR3            Full3                J                E                S              MT
#AICR3   2628.04810164109 2650.21832476474 2626.16281587374 2626.01032201449 2625.9105095972
#X6               X7
#modelR3                H                D
#AICR3   2632.91121516797 2626.26883246402
#remove SSH, but keeping Julian

#Wake
#X1               X2               X3              X4
#modelR3            Full3                J                E               S
#AICR3   301.775364697152 301.774822285016 307.812576354195 309.26853813405
#X5               X6               X7
#modelR3               MT                H                D
#AICR3   307.530052538192 305.443773245884 300.581016975958
#remove Salinity

#Round 4
Full4 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
                s(resSAL, k=-1) + s(mean_SST, k=-1) +resDEN
              , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, k=-1) + s(mean_SST, k=-1) +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(mean_SST, k=-1) +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
           s(resSAL, k=-1) +resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
          s(resSAL, k=-1) + s(mean_SST, k=-1) 
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR4 = c('Full4','J','E','MT','H','D')
AICR4 = c(AIC(Full4),AIC(J),AIC(E),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR4,AICR4))
#Saipan
#                 X1               X2               X3               X4               X5
#modelR4            Full4                J                E               MT                H
#AICR4   2632.91121516797 2663.59874059349 2631.15089084986 2635.96592645301 2632.91121516797
#X6
#modelR4                D
#AICR4   2632.76745478379
#keep remaining variables

#Wake
#X1               X2              X3               X4
#modelR4           Full4                J               E               MT
#AICR4   309.26853813405 310.339492358624 318.94767697495 317.770894843397
#X5               X6
#modelR4                H                D
#AICR4   306.616527605559 308.397900360968
#remove EKE

#Round 5
Full5 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
              s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ s(mean_SST, k=-1) +s(resSSH, k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
           s(resSSH, k=-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(mean_SST, k=-1) +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(mean_SST, k=-1) +s(resSSH, k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR5 = c('Full5','J','MT','H','D')
AICR5 = c(AIC(Full5),AIC(J),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR5,AICR5))
#Saipan
#X1               X2              X3               X4
#modelR5           Full5                J               E                S
#AICR5   2635.9519212653 2659.77120739666 2635.3997284668 2635.73518860949
#X5
#modelR5                D
#AICR5   2634.76393388367
#keep round 5 remaining variables

#Wake
#X1               X2               X3               X4
#modelR5           Full5                J               MT                H
#AICR5   318.94767697495 318.945305508075 327.518768362659 317.370411071748
#X5
#modelR5                D
#AICR5   317.964484697383
#keep or remove MT??

#Final Female Sex Specific full model
FinalFemaleGAM = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+EKE_cm+
                       s(resSAL, k=-1) + s(mean_SST, k=-1) +resDEN
                     , data = TabBinned_Grouped, family = tw, method = "REML")
summary(FinalFemaleGAM)
viz = getViz(FinalFemaleGAM)
vizGG = plot(viz,allTerms = T) +
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
print(vizGG,pages =1)
fig1 =paste(saveDir,site,"_ENV_FemaleGAM.png",sep="")
ggsave(fig1)

#running Sex Specific GAM
  ##Juveniles
#Test how each covariate should be used (linear, smooth, as.factor())
#empty model for comparison
GAM_empty = gam(JuvenileHoursNorm ~ 1, data = TabBinned_GroupedJ, family = tw, method = "REML")

#Julian
GAM_01a = gam(JuvenileHoursNorm ~ Julian, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_01b = gam(JuvenileHoursNorm ~ s(Julian, bs = "cc", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model01 = c('empty','01a','01b')
AIC01 = c(AIC(GAM_empty),AIC(GAM_01a),AIC(GAM_01b))
data.frame(rbind(model01,AIC01))
#Tinian
#X1               X2               X3
#model01          empty              01a              01b
#AIC01   1860.496588857 1863.16815196655 1863.96680521664
#julian as a 'cc' smooth

#Wake
#Julian is always a cc smooth

#Chlorophyll
GAM_02a = gam(JuvenileHoursNorm ~ resChlA, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_02b = gam(JuvenileHoursNorm ~ s(resChlA, bs = "cr", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_02c = gam(JuvenileHoursNorm ~ s(resChlA, k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model02 = c('empty','02a','02b','02c')
AIC02 = c(AIC(GAM_empty),AIC(GAM_02a),AIC(GAM_02b),AIC(GAM_02c))
data.frame(rbind(model02,AIC02))
#Tinian
#                      X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   1634.79054479956 176.769468172997 177.956593246178 178.058378339903
#Chlorophyll as linear

#Wake
#X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   342.418503970956 362.517349918233 366.213503411699 365.488549830909
#Chlorophyll as linear

#EKE
GAM_03a = gam(JuvenileHoursNorm ~ EKE_cm, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_03b = gam(JuvenileHoursNorm ~ s(EKE_cm, bs = "cr", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_03c = gam(JuvenileHoursNorm ~ s(EKE_cm, k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model03 = c('empty','03a','03b','03c')
AIC03 = c(AIC(GAM_empty),AIC(GAM_03a),AIC(GAM_03b),AIC(GAM_03c))
data.frame(rbind(model03,AIC03))
#Tinian
#                      X1               X2               X3               X4
#model03            empty              03a              03b              03c
#AIC03   1634.79054479956 301.429044782972 301.429904903718 301.429858701043
#EKE as linear

#Wake
#X1               X2               X3               X4
#model03            empty              03a              03b              03c
#AIC03   342.418503970956 1817.27695909936 1817.27832551855 1817.27699248743
#EKE as linear

#Salinity
GAM_04a = gam(JuvenileHoursNorm ~ resSAL, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_04b = gam(JuvenileHoursNorm ~ s(resSAL, bs = "cr", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_04c = gam(JuvenileHoursNorm ~ s(resSAL, k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model04 = c('empty','04a','04b','04c')
AIC04 = c(AIC(GAM_empty),AIC(GAM_04a),AIC(GAM_04b),AIC(GAM_04c))
data.frame(rbind(model04,AIC04))
#Tinian
#                      X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   1634.79054479956 301.436877254286 301.437109572235 301.437046935855
#salinity as a normal smooth

#Wake
#X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   342.418503970956 1819.82363864342 1819.82493569296 1819.82384791739
#salinity as a normal smooth

#Mean SST
GAM_05a = gam(JuvenileHoursNorm ~ mean_SST, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_05b = gam(JuvenileHoursNorm ~ s(mean_SST, bs = "cr", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_05c = gam(JuvenileHoursNorm ~ s(mean_SST,  k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model05 = c('empty','05a','05b','05c')
AIC05 = c(AIC(GAM_empty),AIC(GAM_05a),AIC(GAM_05b),AIC(GAM_05c))
data.frame(rbind(model05,AIC05))
#Tinian
#                      X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   1634.79054479956 299.905811865276 299.906104922116 299.906121684052
#Mean SST as 'cr' smooth

#Wake
#X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   342.418503970956 1819.14886738336 1819.14967892526 1819.15022461166
#Mean SST as a normal smooth

#SSH
GAM_06a = gam(JuvenileHoursNorm ~ resSSH, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_06b = gam(JuvenileHoursNorm ~ s(resSSH, bs = "cr", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_06c = gam(JuvenileHoursNorm ~ s(resSSH, k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model06 = c('empty','06a','06b','06c')
AIC06 = c(AIC(GAM_empty),AIC(GAM_06a),AIC(GAM_06b),AIC(GAM_06c))
data.frame(rbind(model06,AIC06))
#Tinian
#                      X1               X2               X3               X4
#model06            empty              06a              06b              06c
#AIC06   1634.79054479956 300.405347919073 300.405765286041 300.405809888613
#SSH as a 'cr' smooth

#Wake
#X1               X2               X3               X4
#model06            empty              06a              06b              06c
#AIC06   342.418503970956 1820.63298939036 1823.09958711937 1823.15456634509
#SSH as a smooth

#Density
GAM_07a = gam(JuvenileHoursNorm ~ resDEN, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_07b = gam(JuvenileHoursNorm ~ s(resDEN, bs = "cr", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_07c = gam(JuvenileHoursNorm ~ s(resDEN, k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model07 = c('empty','07a','07b','07c')
AIC07 = c(AIC(GAM_empty),AIC(GAM_07a),AIC(GAM_07b),AIC(GAM_07c))
data.frame(rbind(model07,AIC07))
#Tinian
#                      X1               X2               X3               X4
#model07            empty              07a              07b              07c
#AIC07   1634.79054479956 300.598651566171 300.598957781108 300.598833322288
#density as linear

#Wake
#X1               X2              X3              X4
#model07            empty              07a             07b             07c
#AIC07   342.418503970956 1820.51669375678 1823.5649468905 1823.7214779085
#density as linear

#Standard deviation of SST
GAM_08a = gam(JuvenileHoursNorm ~ SD_SST, data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_08b = gam(JuvenileHoursNorm ~ s(SD_SST, bs = "cr", k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")
GAM_08c = gam(JuvenileHoursNorm ~ s(SD_SST, k =-1), data = TabBinned_GroupedJ, family = tw, method = "REML")

model08 = c('empty','08a','08b','08c')
AIC08 = c(AIC(GAM_empty),AIC(GAM_08a),AIC(GAM_08b),AIC(GAM_08c))
data.frame(rbind(model08,AIC08))
#Tinian
#                      X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   1634.79054479956 252.232290608707 252.232594826277 252.232499923249
#std dev SST as normal smooth

#Wake
#X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   342.418503970956 833.437434329325 833.437719248008 833.438109830244
#std dev SST as normal smooth

#Test which covariates we should keep for Juvenile Specific GAM
#Round 1
#Initial model
Full1 = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
              s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
            +resDEN+s(SD_SST,k=-1)
            , data = TabBinned_GroupedJ, family = tw, method = "REML")
J = gam(JuvenileHoursNorm ~ resChlA+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
C = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
E = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
S = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
           s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
MT = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
           s(resSAL,k=-1) + s(resSSH,bs="cr",k=-1)
         +resDEN+s(SD_SST,k=-1)
         , data = TabBinned_GroupedJ, family = tw, method = "REML")
H = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) 
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
D = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
SDT = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
            s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
          +resDEN
          , data = TabBinned_GroupedJ, family = tw, method = "REML")
modelR1 = c('Full1','J','C','E','S','MT','H','D','SDT')
AICR1 = c(AIC(Full1),AIC(J),AIC(C),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR1,AICR1))
#Tinian
#                      X1               X2               X3               X4               X5
#modelR1            Full1                J                C                E                S
#AICR1   238.182297822619 229.109452423101 265.406542894899 201.224306549872 213.169466037939
#X6              X7               X8               X9
#modelR1               MT               H                D              SDT
#AICR1   214.950428094598 222.00639349114 219.945216977607 212.385066933658
#remove chlorophyll

#Wake
#X1               X2               X3               X4
#modelR1           Full1                J                C                E
#AICR1   227.85489372494 227.855281539179 278.551753705583 226.395177598489
#X5               X6               X7               X8
#modelR1                S               MT                H                D
#AICR1   226.759249276606 225.767507088368 223.235356678859 225.561289355681
#X9
#modelR1              SDT
#AICR1   266.865821200058
#remove chlorophyll

#Round 2
Full2= gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
             s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
           +resDEN+s(SD_SST,k=-1)
           , data = TabBinned_GroupedJ, family = tw, method = "REML")
J = gam(JuvenileHoursNorm ~ EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
E = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
S = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
MT = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
           s(resSAL,k=-1) + s(resSSH,bs="cr",k=-1)
         +resDEN+s(SD_SST,k=-1)
         , data = TabBinned_GroupedJ, family = tw, method = "REML")
H = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) 
        +resDEN+s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
D = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +s(SD_SST,k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
SDT = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
            s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
          +resDEN
          , data = TabBinned_GroupedJ, family = tw, method = "REML")
modelR2 = c('Full2','J','E','S','MT','H','D','SDT')
AICR2 = c(AIC(Full2),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR2,AICR2))
#Tinian
#                      X1               X2               X3               X4               X5
#modelR2            Full2                J                E                S               MT
#AICR2   265.406542894899 265.406288575829 261.471600529843 261.444080799133 262.374395342803
#X6               X7               X8
#modelR2               H                D              SDT
#AICR2   261.35632287089 261.017244276263 312.040727525366
#remove std dev SST

#Wake
#X1               X2               X3               X4
#modelR2            Full2                J                E                S
#AICR2   278.551753705583 278.552387495673 279.869689617303 276.793229635488
#X5               X6               X7               X8
#modelR2               MT                H                D              SDT
#AICR2   276.183619254203 274.243219909598 276.166098128741 346.272736933212
#remove std dev of SST

#Round 3
Full3 = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
              s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
            +resDEN
            , data = TabBinned_GroupedJ, family = tw, method = "REML")
J = gam(JuvenileHoursNorm ~ EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
E = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
S = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        +resDEN
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
MT = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
           s(resSAL,k=-1) + s(resSSH,bs="cr",k=-1)
         +resDEN
         , data = TabBinned_GroupedJ, family = tw, method = "REML")
H = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) 
        +resDEN
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
D = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
        , data = TabBinned_GroupedJ, family = tw, method = "REML")
modelR3 = c('Full3','J','E','S','MT','H','D')
AICR3 = c(AIC(Full3),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR3,AICR3))
#Tinian
#                      X1               X2               X3               X4               X5
#modelR3            Full3                J                E                S               MT
#AICR3   312.040727525366 312.040578287567 308.602220949753 308.599782815027 308.593933562717
#X6               X7
#modelR3                H                D
#AICR3   308.360393991731 308.144069393953
#keep remaining variables

#Wake
#X1               X2               X3               X4
#modelR3            Full3                J                E                S
#AICR3   346.272736933212 345.910184358419 350.152489489062 343.464873754557
#X5              X6               X7
#modelR3               MT               H                D
#AICR3   347.011165230058 341.19088573123 343.817152989437
#keep remaining variables

#Final Juvenile Sex Specific full model
FinalJuvenileGAM = gam(JuvenileHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
                         s(resSAL,k=-1) + s(mean_SST,bs="cr",k=-1) +s(resSSH,bs="cr",k=-1)
                       +resDEN
                       , data = TabBinned_GroupedJ, family = tw, method = "REML")
summary(FinalJuvenileGAM)
viz = getViz(FinalJuvenileGAM)
vizGG = plot(viz,allTerms = T) +
  labs(ylab = "Sperm Whale Presence (Hours/Day)")+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
print(vizGG,pages =1)
fig1 =paste(saveDir,site,"_ENV_JuvenileGAM.png",sep="")
ggsave(fig1)

#running Sex Specific GAM
  ##Males
#Test how each covariate should be used (linear, smooth, as.factor())
#empty model for comparison
GAM_empty = gam(MaleHoursNorm ~ 1, data = TabBinned_GroupedM, family = tw, method = "REML")

#Julian
GAM_01a = gam(MaleHoursNorm ~ Julian, data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_01b = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model01 = c('empty','01a','01b')
AIC01 = c(AIC(GAM_empty),AIC(GAM_01a),AIC(GAM_01b))
data.frame(rbind(model01,AIC01))
#Saipan
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   887.529203151121 878.265022630491 873.011160556566
#Julian as a 'cc' smooth

#Tinian
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   2110.00713779915 2109.03174260952 2109.74342449606
#Julian as 'cc' smooth

#Chlorophyll
GAM_02a = gam(MaleHoursNorm ~ resChlA, data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_02b = gam(MaleHoursNorm ~ s(resChlA, bs = "cr", k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_02c = gam(MaleHoursNorm ~ s(resChlA, k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model02 = c('empty','02a','02b','02c')
AIC02 = c(AIC(GAM_empty),AIC(GAM_02a),AIC(GAM_02b),AIC(GAM_02c))
data.frame(rbind(model02,AIC02))
#Saipan
#                      X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   2141.48369945917 607.706638891198 609.502715316711 609.601044701568
#chlorophyll as linear

#Tinian
#X1               X2               X3               X4
#model02            empty              02a              02b              02c
#AIC02   2110.00713779915 272.564454835334 274.423469134624 274.494240318998
#Chlorophyll as linear

#EKE
GAM_03a = gam(MaleHoursNorm ~ EKE_cm, data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_03b = gam(MaleHoursNorm ~ s(EKE_cm, bs = "cr", k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_03c = gam(MaleHoursNorm ~ s(EKE_cm, k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model03 = c('empty','03a','03b','03c')
AIC03 = c(AIC(GAM_empty),AIC(GAM_03a),AIC(GAM_03b),AIC(GAM_03c))
data.frame(rbind(model03,AIC03))
#Saipan
#                      X1               X2               X3               X4
#model03            empty              03a              03b              03c
#AIC03   2141.48369945917 2786.26146152255 2787.03966288725 2787.24552034257
#EKE as linear

#Tinian
#                      X1               X2              X3               X4
#model03            empty              03a             03b              03c
#AIC03   2110.00713779915 2109.17637455186 2109.1786916091 2109.17688323346
#EKE as linear

#Salinity
GAM_04a = gam(MaleHoursNorm ~ resSAL, data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_04b = gam(MaleHoursNorm ~ s(resSAL, bs = "cr", k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_04c = gam(MaleHoursNorm ~ s(resSAL, k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model04 = c('empty','04a','04b','04c')
AIC04 = c(AIC(GAM_empty),AIC(GAM_04a),AIC(GAM_04b),AIC(GAM_04c))
data.frame(rbind(model04,AIC04))
#Saipan
#                      X1               X2               X3               X4
#model04            empty              04a              04b              04c
#AIC04   2141.48369945917 2784.34360193609 2784.63871105464 2784.65889373016
#Salinity as a 'cr' smooth

#Tinian
#                      X1              X2               X3               X4
#model04            empty             04a              04b              04c
#AIC04   2110.00713779915 2108.1701495849 2109.89348283305 2109.90740688014
#Salinity as a 'cr' smooth

#Mean SST
GAM_05a = gam(MaleHoursNorm ~ mean_SST, data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_05b = gam(MaleHoursNorm ~ s(mean_SST, bs = "cr", k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_05c = gam(MaleHoursNorm ~ s(mean_SST, k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model05 = c('empty','05a','05b','05c')
AIC05 = c(AIC(GAM_empty),AIC(GAM_05a),AIC(GAM_05b),AIC(GAM_05c))
data.frame(rbind(model05,AIC05))
#Saipan
#                     X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   2141.48369945917 2776.67511827778 2776.67701022956 2776.67732123249
#Mean SST as a 'cr' smooth

#Tinian
#                      X1               X2               X3               X4
#model05            empty              05a              05b              05c
#AIC05   2110.00713779915 2110.55989643511 2112.92759367526 2114.06987391218
#Mean SST as linear

#SSH
GAM_06a = gam(MaleHoursNorm ~ resSSH, data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_06b = gam(MaleHoursNorm ~ s(resSSH, bs = "cr", k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_06c = gam(MaleHoursNorm ~ s(resSSH, k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model06 = c('empty','06a','06b','06c')
AIC06 = c(AIC(GAM_empty),AIC(GAM_06a),AIC(GAM_06b),AIC(GAM_06c))
data.frame(rbind(model06,AIC06))
#Saipan
#                      X1               X2               X3               X4
#model06            empty              06a              06b              06c
#AIC06   2141.48369945917 2785.62034294709 2785.62568962918 2785.62534403544
#SSH as a normal smooth

#Tinian
#                      X1               X2               X3               X4
#model06            empty              06a              06b              06c
#AIC06   2110.00713779915 2112.19464402909 2114.44309409673 2114.94189041267
#SSH as a linear

#Density
GAM_07a = gam(MaleHoursNorm ~ s(resDEN, bs="cc", k=-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_07b = gam(MaleHoursNorm ~ s(resDEN, bs="cr", k=-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_07c = gam(MaleHoursNorm ~ s(resDEN, k=-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model07 = c('empty','07a','07b','07c')
AIC07 = c(AIC(GAM_empty),AIC(GAM_07a),AIC(GAM_07b),AIC(GAM_07c))
data.frame(rbind(model07,AIC07))
#Saipan
#                      X1               X2               X3               X4
#model07            empty              07a              07b              07c
#AIC07   2141.48369945917 2783.49768165688 2779.51574405995 2779.51525231626
#density as a normal smooth

#Tinian
#                      X1               X2               X3               X4
#model07            empty              07a              07b              07c
#AIC07   2110.00713779915 2111.20528552966 2109.29214328276 2109.29214180054
#density as a normal smooth

#Standard deviation of SST
GAM_08a = gam(MaleHoursNorm ~ SD_SST, data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_08b = gam(MaleHoursNorm ~ s(SD_SST, bs = "cr", k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")
GAM_08c = gam(MaleHoursNorm ~ s(SD_SST, k =-1), data = TabBinned_GroupedM, family = tw, method = "REML")

model08 = c('empty','08a','08b','08c')
AIC08 = c(AIC(GAM_empty),AIC(GAM_08a),AIC(GAM_08b),AIC(GAM_08c))
data.frame(rbind(model08,AIC08))
#Saipan
#                      X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   2141.48369945917 998.414840258139 998.426127575516 998.525085238502
#std dev SST as a cr smooth

#Tinian
#                      X1               X2               X3               X4
#model08            empty              08a              08b              08c
#AIC08   2110.00713779915 787.106333769798 787.106972961245 787.107663034754
#std dev SST as 'cr' smooth

#Test which covariates we should keep for Male Specific GAM
#Round 1
#Initial model
Full1 = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
              s(resSAL, bs = "cr", k=-1) + mean_SST +s(resSSH,k=-1)+
              s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
            , data = TabBinned_GroupedM, family = tw, method = "REML")
J = gam(MaleHoursNorm ~ resChlA+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
C = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
E = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+
         s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
         s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
       , data = TabBinned_GroupedM, family = tw, method = "REML")
S = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
          mean_SST  +s(resSSH,k=-1)+
           s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
         , data = TabBinned_GroupedM, family = tw, method = "REML")
MT = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
           s(resSAL, bs = "cr", k=-1) +s(resSSH,k=-1)+
           s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
         , data = TabBinned_GroupedM, family = tw, method = "REML")
H = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +
          s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
D = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
SDT = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+resChlA+EKE_cm+
            s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
            s(resDEN,k=-1)
          , data = TabBinned_GroupedM, family = tw, method = "REML")
modelR1 = c('Full1','J','C','E','S','MT','H','D','SDT')
AICR1 = c(AIC(Full1),AIC(J),AIC(C),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR1,AICR1))
#Saipan
#X1               X2               X3               X4              X5
#modelR1            Full1                J                C                E               S
#AICR1   538.211424869672 538.210457885585 1009.49796194627 535.472545978177 535.45423807791
#X6               X7               X8               X9
#modelR1               MT                H                D              SDT
#AICR1   532.726874701811 538.072420030499 533.804832445118 624.481470935852
#remove chlorophyll

#Tinian
#                      X1               X2               X3               X4               X5
#modelR1            Full1                J                C                E                S
#AICR1   305.515939839273 321.169959519976 794.538582896943 271.428432911797 269.883837343044
#X6               X7               X8               X9
#modelR1               MT                H                D              SDT
#AICR1   266.781623446821 267.721694311764 268.244460480538 281.745841887533
#remove chlorophyll

#Round 2
Full2 = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
              s(resSAL, bs = "cr", k=-1) + mean_SST +s(resSSH,k=-1)+
              s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
            , data = TabBinned_GroupedM, family = tw, method = "REML")
J = gam(MaleHoursNorm ~ EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
E = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
S = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
MT = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
           s(resSAL, bs = "cr", k=-1) +s(resSSH,k=-1)+
           s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
         , data = TabBinned_GroupedM, family = tw, method = "REML")
H = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +
          s(resDEN,k=-1)+s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
D = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(SD_SST, bs="cr",k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
SDT = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
            s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
            s(resDEN,k=-1)
          , data = TabBinned_GroupedM, family = tw, method = "REML")
modelR2 = c('Full2','J','E','S','MT','H','D','SDT')
AICR2 = c(AIC(Full2),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR2,AICR2))
#Saipan
#                      X1               X2               X3               X4               X5
#modelR2            Full2                J                E                S               MT
#AICR2   1009.49796194627 1013.01689582351 1007.01440508224 1007.15754772359 1007.97221926775
#X6               X7               X8
#modelR2                H                D              SDT
#AICR2   1008.21151304873 1006.87024650936 2781.55704035258
#remove std dev of SST

#Tinian
#                      X1               X2               X3               X4               X5
#modelR2            Full2                J                E                S               MT
#AICR2   794.538582896943 794.537692662672 792.664926001935 795.977605921437 795.774852633116
#X6               X7               X8
#modelR2                H                D              SDT
#AICR2   792.110464371626 793.713095132094 2120.02391853154
#remove std dev of SST

#Round 3
Full3 = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
              s(resSAL, bs = "cr", k=-1) + mean_SST +s(resSSH,k=-1)+
              s(resDEN,k=-1)
            , data = TabBinned_GroupedM, family = tw, method = "REML")
J = gam(MaleHoursNorm ~ EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
E = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
S = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          mean_SST  +s(resSSH,k=-1)+
          s(resDEN,k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
MT = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
           s(resSAL, bs = "cr", k=-1) +s(resSSH,k=-1)+
           s(resDEN,k=-1)
         , data = TabBinned_GroupedM, family = tw, method = "REML")
H = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +
          s(resDEN,k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
D = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
          s(resSAL, bs = "cr", k=-1) + mean_SST  +s(resSSH,k=-1)
        , data = TabBinned_GroupedM, family = tw, method = "REML")
modelR3 = c('Full3','J','E','S','MT','H','D')
AICR3 = c(AIC(Full3),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR3,AICR3))
#Saipan
#                      X1               X2               X3             X4              X5
#modelR3            Full3                J                E              S              MT
#AICR3   2781.55704035258 2783.16566069915 2780.51352187552 2780.285906716 2779.3983711904
#X6               X7
#modelR3                H                D
#AICR3   2781.79858098144 2779.32542286913
#keep round 3 remaining variables

#Tinian
#                      X1              X2               X3               X4               X5
#modelR3            Full3               J                E                S               MT
#AICR3   2120.02391853154 2117.9815376294 2121.89786924475 2116.26769502815 2118.20134842514
#X6               X7
#modelR3                H                D
#AICR3   2114.70242695529 2118.97781436002
#keep round 3 remaining variables

#Final Male Sex Specific full model
FinalMaleGAM = gam(MaleHoursNorm ~ s(Julian, bs="cc", k=-1)+EKE_cm+
                     s(resSAL, bs = "cr", k=-1) + mean_SST +s(resSSH,k=-1)+
                     s(resDEN,k=-1)
                   , data = TabBinned_GroupedM, family = tw, method = "REML")
summary(FinalMaleGAM)
viz = getViz(FinalMaleGAM)
vizGG = plot(viz,allTerms = T) +
  labs(ylab = "Sperm Whale Presence (Hours/Day)")+
  l_fitLine(linetype = 1, size = 2)  +
  l_fitContour()+
  #l_ciLine(mul = 5, colour = "blue", linetype = 2) +
  l_ciPoly(level = 0.95, alpha = 1/2)+
  l_ciBar() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20,face="bold"))
print(vizGG,pages =1)
fig1 =paste(saveDir,site,"_ENV_MaleGAM.png",sep="")
ggsave(fig1)

#GAM_trial = gam(HoursNorm ~ s(Julian, bs = "cc", k = -1) + +s(resChlA, bs = "cc", k = -1) +
                 # s(EKE_cm, bs = "cc", k = -1) + s(resSAL, bs = "cc", k = -1) + 
                 # s(mean_SST, bs = "cc", k = -1) + s(resSSH, bs = "cc", k = -1) + 
                 # s(resDEN, bs = "cc", k = -1) + s(SD_SST, bs = "cc", k=-1),
          #  data = TabBinned_Grouped, family = tw, method = "REML")

#summary(GAM_trial)
#plot(GAM_trial, pages = 1)
#viz = getViz(GAM_trial)
#print(plot(viz,allTerms=T),pages=1)

#gam.check(GAM_trial)
#concurvity(GAM_trial, full= TRUE)
#concurvity(GAM_trial, full= FALSE)

#do we want 2D Gams? 
GAM_2D = gam(HoursNorm ~ s(Julian, resChlA, bs = "fs", k = -1) +
                  s(EKE_cm, bs = "cc", k = -1) + s(resSAL, bs = "cc", k = -1) + 
                  s(mean_SST, bs = "cc", k = -1) + s(resSSH, bs = "cc", k = -1) + 
                  s(resDEN, bs = "cc", k = -1) + s(SD_SST, bs = "cc", k=-1),
                data = TabBinned_Grouped, family = tw, method = "REML")
plot(GAM_2D, pages = 1)
plot(GAM_2D, scheme= 2)#change scheme for dif plots

vis.gam(x = GAM_trial, view = c("EKE_cm", "resSSH"), plot.type = "contour") 
      #plot types include: persp, contour
      #not displaying 3D plot from tutorial
      #add too.far function to specify which data should not be included
      #se function to create high and low predictions surfaces
      #alter orientation: theta = horizontal, phi = vertical, r = zoom

#GAM_tensor = gam using te()
plot(GAM_tensor)

#plot GAMs

