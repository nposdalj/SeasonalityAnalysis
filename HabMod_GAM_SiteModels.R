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

#increasing memory limit
memory.limit(size=300000)

# User Defined sections
#define the lat and long of interest
#df1 = data.frame("lat" = c(19.29, 19.2, 19.2467, 19.2467), "long" = c(-166.69, -166.69, -166.74, -166.64)) #Wake
#df1 = data.frame("lat" = c(15.36, 15.27, 15.3186, 15.3186), "long" = c(145.46, 145.46, 145.51, 145.41)) #Saipan
df1 = data.frame("lat" = c(15.08, 14.99, 15.0387, 15.0387), "long" = c(145.75, 145.75, 145.8, 145.7))#Tinian
#df1 = data.frame("lat" = c(29.185, 29.1, 29.1410, 29.1410), "long" = c(-118.26, -118.26, -118.31, -118.21))#GI
#df1 = data.frame("lat" = c(31.79, 31.7, 31.747, 31.747), "long" = c(-121.38, -121.38, -121.43, -121.33))#CORC
#df1 = data.frame("lat" = c(52.4, 52.31, 52.3547, 52.3547), "long" = c(-175.635, -175.635, -175.71, -175.56))#BD
#df1 = data.frame("lat" = c(47.54, 47.45, 47.4936, 47.4936), "long" = c(-125.378, -125.378, -125.44, -125.31))#QC
#df1 = data.frame("lat" = c(56.385, 56.295, 56.34, 56.34), "long" = c(-145.182, -145.182, -145.26, -145.1))#QN
#df1 = data.frame("lat" = c(56.29, 56.2, 56.2434, 56.2434), "long" = c(-142.75, -142.75, -142.83, -142.67))#PT
#df1 = data.frame("lat" = c(58.71, 58.62, 58.6668, 58.6668), "long" = c(-148.0034, -148.0034, -148.12, -147.94))#CB

#define the start and end of the datame 
startTime = "2011-04-13" #this should be formatted like this: 2010-03-05
endTime = "2019-05-13" 

#ITS
#ITS = 4 #Saipan
#ITSF = 6 #Saipan Females
#ITSM = 1 #Saipan Males
ITS = 2 #Tinian
ITSJ = 1 #Tinian Juveniles
ITSM = 1 #Tinian Males

#loading the environmental data
envDir = paste("O:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory)

#spatial polygon for area of interest
ch <- chull(df1$long, df1$lat)
coords <- df1[c(ch, ch[1]), ]#creating convex hull
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = 1)))#converting convex hull to spatial polygon

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

#removing unnecessary variables
rm("DayData")

#clear memory 
gc()


#chlorophyll data
#filenameStatAll = paste(envDir,"Chl2.csv",sep="")#load files as data frame
#chl = read.csv(filenameStatAll)
#chl = chl[-1,] #delete first row
#chl$latitude = as.numeric(chl$latitude)
#chl$longitude = as.numeric(chl$longitude)
#chl = chl[complete.cases(chl[ , 2:3]),]#remove any rows with lat or long as na

#loading as .ncfile
filenameStatAll = paste(envDir,"Chl2.nc",sep="")#load files as data frame
ChlA = nc_open(filenameStatAll)
v1=ChlA$var[[1]]
ChlAvar=ncvar_get(ChlA,v1)
ChlA_lon=v1$dim[[1]]$vals
ChlA_lat=v1$dim[[2]]$vals
ChlA_dates= as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT')
ChlA_dates = as.Date(ChlA_dates)

#CHLA
#plotting in ggplot
r = raster(t(ChlAvar[,,1]),xmn = min(ChlA_lon),xmx = max(ChlA_lon),ymn=min(ChlA_lat),ymx=max(ChlA_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="Chl"
mid = mean(df$Chl)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = Chl)) + 
  ggtitle(paste("Daily Chl on", ChlA_dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting timeseries
I=which(ChlA_lon>=min(df1$long) & ChlA_lon<= max(df1$long)) #only extract the region we care about
J=which(ChlA_lat>=min(df1$lat) & ChlA_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(ChlA_dates>= startTime & ChlA_dates<= endTime) #extract only the dates we care about
ChlA2=ChlAvar[I,J,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(ChlA2)[3] #find the length of time

#take the mean
resChlA=rep(NA,n) 
for (i in 1:n) 
  resChlA[i]=mean(ChlA2[,,i],na.rm=TRUE) 

#plot the time series
plot(1:n,resChlA,axes=FALSE,type='o',pch=20,xlab='',ylab='SSH',las = 3) 
axis(2) 
axis(1,1:n,format(resChlA[K]),las = 3) 
box()

#subset the dataframe based on the area of interest
#coordinates(chl) <- c("latitude", "longitude")
#coords <- over(chl, sp_poly)
#chl2 <- chl[coords == 1 & !is.na(coords),]
#chl2$chlorophyll = as.numeric(chl2$chlorophyll)#converting SST from character to numeric
#chl2$time = as.Date(chl2$time)#converting time from character to date
#chl3 = as.data.frame(chl2)#converting SPDF back to DF
#chl3$chlorophyll[chl3$chlorophyll < 0] <- NA #making anything<0 NA

#average the environmental variable based on the ITS over the area of interest
#chl4 = chl3 %>%
 #mutate(time = floor_date(time)) %>%
  #group_by(time) %>%
  #summarize(mean_chl = mean(chlorophyll), SD_chl = sd(chlorophyll)) #finding daily mean

#data exploration
#mean(ChlA$chlorophyll, na.rm = TRUE)#finding overall mean
#save standard deviation
#sd(chl2$chlorophyll, na.rm = TRUE)
#SST histogram
#hist(chl2$chlorophyll)

#plot time series
#plot(chl4$time, chl4$mean_chl)#exploratory plot
#title1 = paste(site,"Chlorophyll Plot")
#ggplot(chl4, aes(x=time,y=mean_chl))+
  #ggtitle(title1)+
  #labs(y="Mean Chlorophyll (mg m-3)",x="Time (days)")+
  #geom_line()+
  #geom_point()

#rm(chl)
#rm(chl2)
#rm(chl3)

#clear memory 
gc()

#SST data
filenameStatAll = paste(envDir,"AQUASST_SAPTIN.csv",sep="")#load files as data frame
SST = read.csv(filenameStatAll)
SST = SST[-1,] #delete first row
SST$latitude = as.numeric(SST$latitude)
SST$longitude = as.numeric(SST$longitude)
SST = SST[complete.cases(SST[ , 2:3]),]#remove any rows with lat or long as na

#subset the dataframe based on the area of interest
coordinates(SST) <- c("latitude", "longitude")
coords <- over(SST, sp_poly)
SST2 <- SST[coords == 1 & !is.na(coords),]
SST2$sstMasked = as.numeric(SST2$sstMasked)#converting SST from character to numeric
SST2$time = as.Date(SST2$time)#converting time from character to date
SST3 = as.data.frame(SST2)#converting SPDF back to DF

#average the environmental variable based on the ITS over the area of interest
SST4 = SST3 %>%
  mutate(time = floor_date(time)) %>%
  group_by(time) %>%
  summarize(mean_SST = mean(sstMasked), SD_SST = sd(sstMasked)) #finding daily mean

#data exploration
mean(SST2$sstMasked, na.rm = TRUE)#finding overall mean
#save standard deviation
sd(SST2$sstMasked, na.rm = TRUE)
#SST histogram
hist(SST2$sstMasked)

#plot time series
plot(SST4$time, SST4$mean_SST)#exploratory plot
title1 = paste(site,"Sea Surface Temperature Plot")
ggplot(SST4, aes(x=time,y=mean_SST))+
  ggtitle(title1)+
  labs(y="Mean SST (C)",x="Time (days)")+
  geom_line()+
  geom_point()

rm(SST)
rm(SST2)
rm(SST3)
###############################
#SSH anomaly data
filenameStatAll = paste(envDir,"AVISOglobalvars.nc",sep="")#load files as data frame
AVISO = nc_open(filenameStatAll)
names(AVISO$var)

#zos - SSH
v6=AVISO$var[[6]]
SSHvar=ncvar_get(AVISO,v6)
SSH_lon=v6$dim[[1]]$vals
SSH_lat=v6$dim[[2]]$vals
SSH_dates=as.POSIXlt(v6$dim[[3]]$vals*60*60,origin='1950-01-01') #extract the date/time
SSH_dates = as.Date(SSH_dates, format = "%m/%d/%y") #get rid of the time

#mlotst - density ocean mixed layer thickness
v1=AVISO$var[[1]]
DENvar=ncvar_get(AVISO,v1)
DEN_lon=v1$dim[[1]]$vals
DEN_lat=v1$dim[[2]]$vals
DEN_dates=as.POSIXlt(v1$dim[[3]]$vals*60*60,origin='1950-01-01') #extract the date/time
DEN_dates = as.Date(DEN_dates, format = "%m/%d/%y") #get rid of the time

#so - salinity
v2=AVISO$var[[5]]
SALvar=ncvar_get(AVISO,v2)
SAL_lon=v2$dim[[1]]$vals
SAL_lat=v2$dim[[2]]$vals
SAL_dates=as.POSIXlt(v2$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
SAL_dates = as.Date(SAL_dates, format = "%m/%d/%y") #get rid of the time

#thetao - temperature
v3=AVISO$var[[3]]
TEMPvar=ncvar_get(AVISO,v3)
TEMP_lon=v3$dim[[1]]$vals
TEMP_lat=v3$dim[[2]]$vals
TEMP_dates=as.POSIXlt(v3$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
TEMP_dates = as.Date(TEMP_dates, format = "%m/%d/%y") #get rid of the time

#uo - eastward velocity
v4=AVISO$var[[4]]
EASTVvar=ncvar_get(AVISO,v4)
EASTV_lon=v4$dim[[1]]$vals
EASTV_lat=v4$dim[[2]]$vals
EAST_dates=as.POSIXlt(v4$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
EAST_dates = as.Date(EAST_dates, format = "%y/%m/%d") #get rid of the time
EASTdf <- as.data.frame(EASTVvar)

#vo - northward velocity
v5=AVISO$var[[2]]
NORVvar=ncvar_get(AVISO,v5)
NORV_lon=v5$dim[[1]]$vals
NORV_lat=v5$dim[[2]]$vals
NOR_dates=as.POSIXlt(v5$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
NOR_dates = as.Date(NOR_dates, format = "%m/%d/%y") #get rid of the time
NORdf <- as.data.frame(NORVvar)

#removing unnecessary variables
rm("v1","v2","v3","v4","v5","v6")

#Plotting maps and time series
#loading the world
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#SSH
#plotting in ggplot
r = raster(t(SSHvar[,,1]),xmn = min(SSH_lon),xmx = max(SSH_lon),ymn=min(SSH_lat),ymx=max(SSH_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="SSH"
mid = mean(df$SSH)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = SSH)) + 
  ggtitle(paste("Daily SSH on", SSH_dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting timeseries
I=which(SSH_lon>=min(df1$long) & SSH_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(SSH_lat>=min(df1$lat) & SSH_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K=which(SSH_dates>= startTime & SSH_dates<= endTime) #extract only the dates we care about

SSH2=SSHvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(SSH2)[3] #find the length of time

#take the mean
resSSH=rep(NA,n) 
for (i in 1:n) 
  resSSH[i]=mean(SSH2[,,i],na.rm=TRUE) 

#plot the time series
plot(1:n,resSSH,axes=FALSE,type='o',pch=20,xlab='',ylab='SSH',las = 3) 
axis(2) 
axis(1,1:n,format(SSH_dates[K]),las = 3) 
box()

#remove unnecessary variables
rm("SSHvar","SSH2", "SSH_lon","SSH_lat")

#Density ocean mixed layer thickness
#Plotting in ggplot
r = raster(t(DENvar[,,1]),xmn = min(DEN_lon),xmx = max(DEN_lon),ymn=min(DEN_lat),ymx=max(DEN_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="DEN"
mid = mean(df$DEN)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = DEN)) + 
  ggtitle(paste("Daily Density Ocean Mized Layer Thickness on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting time series SAPTIN 
I=which(DEN_lon>=min(df1$long) & DEN_lon<= max(df1$long)) #only extract the region we care about
J=which(DEN_lat>=min(df1$lat) & DEN_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(DEN_dates>= startTime & DEN_dates<= endTime) #extract only the dates we care about
DEN2=DENvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(DEN2)[3] #find the length of time

#take the mean
resDEN=rep(NA,n) 
for (i in 1:n) 
  resDEN[i]=mean(DEN2[,,i],na.rm=TRUE) 

#plot the time series
plot(1:n,resDEN,axes=FALSE,type='o',pch=20,xlab='',ylab='Density',las = 3) 
axis(2) 
axis(1,1:n,format(DEN_dates[K]),las = 3) 
box()

#remove unnecessary variables
rm("DENvar","DEN2", "DEN_lon","DEN_lat")

#Salinity
#Plotting in ggplot
r = raster(t(SALvar[,,1]),xmn = min(SAL_lon),xmx = max(SAL_lon),ymn=min(SAL_lat),ymx=max(SAL_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="SAL"
mid = mean(df$SAL)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = SAL)) + 
  ggtitle(paste("Daily Salinity on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting time series SAPTIN 
I=which(SAL_lon>=min(df1$long) & SAL_lon<= max(df1$long)) #only extract the region we care about
J=which(SAL_lat>=min(df1$lat) & SAL_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(SAL_dates>= startTime & SAL_dates<= endTime) #extract only the dates we care about
SAL2=SALvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(SAL2)[3] #find the length of time

#take the mean
resSAL=rep(NA,n) 
for (i in 1:n) 
  resSAL[i]=mean(SAL2[,,i],na.rm=TRUE) 

#plot the time series
plot(1:n,resSAL,axes=FALSE,type='o',pch=20,xlab='',ylab='Salinity',las = 3) 
axis(2) 
axis(1,1:n,format(SAL_dates[K]),las = 3) 
box()

#remove unnecessary variables
rm("SALvar","SAL2", "SAL_lon","SAL_lat")

#Temperature
#Plotting in ggplot
r = raster(t(TEMPvar[,,1]),xmn = min(TEMP_lon),xmx = max(TEMP_lon),ymn=min(TEMP_lat),ymx=max(TEMP_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="TEMP"
mid = mean(df$TEMP)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = TEMP)) + 
  ggtitle(paste("Daily Temperature on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting time series SAPTIN 
I=which(TEMP_lon>=min(df1$long) & TEMP_lon<= max(df1$long)) #only extract the region we care about
J=which(TEMP_lat>=min(df1$lat) & TEMP_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(TEMP_dates>= startTime & TEMP_dates<= endTime) #extract only the dates we care about
TEMP2=TEMPvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(TEMP2)[3] #find the length of time

#take the mean
resTEMP=rep(NA,n) 
for (i in 1:n) 
  resTEMP[i]=mean(TEMP2[,,i],na.rm=TRUE)

#plot the time series
plot(1:n,resTEMP,axes=FALSE,type='o',pch=20,xlab='',ylab='Temperature',las = 3) 
axis(2) 
axis(1,1:n,format(TEMP_dates[K]),las = 3) 
box()

#remove unnecessary variables
rm("TEMPvar","TEMP2", "TEMP_lon","TEMP_lat")

#Eastward Velocity
#Plotting in ggplot
r = raster(t(EASTVvar[,,1]),xmn = min(EASTV_lon),xmx = max(EASTV_lon),ymn=min(EASTV_lat),ymx=max(EASTV_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="EASTV"
mid = mean(df$EASTV)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = EASTV)) + 
  ggtitle(paste("Daily Eastward Velocity on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting time series SAPTIN 
I=which(EASTV_lon>=min(df1$long) & EASTV_lon<= max(df1$long)) #only extract the region we care about
J=which(EASTV_lat>=min(df1$lat) & EASTV_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(EAST_dates>= startTime & EAST_dates<= endTime) #extract only the dates we care about
EASTV2=EASTVvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(EASTV2)[3] #find the length of time

#take the mean
resEV=rep(NA,n) 
for (i in 1:n) 
  resEV[i]=mean(EASTV2[,,i],na.rm=TRUE) 

#plot the time series
plot(1:n,resEV,axes=FALSE,type='o',pch=20,xlab='',ylab='Eastward Velocity',las = 3) 
axis(2) 
axis(1,1:n,format(EAST_dates[K]),las = 3) 
box()

#remove unnecessary variables
rm("EASTVvar","EASTV2", "EASTV_lon","EASTV_lat")

#Northward Velocity
#Plotting in ggplot
r = raster(t(NORVvar[,,1]),xmn = min(NORV_lon),xmx = max(NORV_lon),ymn=min(NORV_lat),ymx=max(NORV_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="NORV"
mid = mean(df$NORV)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = NORV)) + 
  ggtitle(paste("Daily Northward Velocity on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting time series SAPTIN 
I=which(NORV_lon>=min(df1$long) & NORV_lon<= max(df1$long)) #only extract the region we care about
J=which(NORV_lat>=min(df1$lat) & NORV_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(NOR_dates>= startTime & NOR_dates<= endTime) #extract only the dates we care about
NORV2=NORVvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(NORV2)[3] #find the length of time

#take the mean
resNV=rep(NA,n) 
for (i in 1:n) 
  resNV[i]=mean(NORV2[,,i],na.rm=TRUE)

#plot the time series
plot(1:n,resNV,axes=FALSE,type='o',pch=20,xlab='',ylab='Northward Velocity',las = 3) 
axis(2) 
axis(1,1:n,format(NOR_dates[K]),las = 3) 
box()

#remove unnecessary variables
rm("NORVvar","NORV2", "NORV_lon","NORV_lat")

#remove universal variables
rm("I","JJ","J")

#subset the dataframe based on the area of interest
#average the environmental variable based on the ITS over the area of interest
#save standard deviation 
 

####

#Calculate EKE
u <- (resEV)^2
v <- (resNV)^2
velocity = u + v
EKE_meters = 0.5 * velocity
EKE_cm = EKE_meters * 10000 

#converting _dates to data frames and renaming column to 'time'
SSH_ddf <- as.data.frame(SSH_dates[K])
SSH_ddf = SSH_ddf %>% 
  dplyr::rename(
    time = 'SSH_dates[K]',
  )
DEN_ddf <- as.data.frame(DEN_dates[K])
DEN_ddf = DEN_ddf %>% 
  dplyr::rename(
    time = 'DEN_dates[K]',
  )
SAL_ddf <- as.data.frame(SAL_dates[K])
SAL_ddf = SAL_ddf %>% 
  dplyr::rename(
    time = 'SAL_dates[K]',
  )
TEMP_ddf <- as.data.frame(TEMP_dates[K])
TEMP_ddf = TEMP_ddf %>% 
  dplyr::rename(
    time = 'TEMP_dates[K]',
  )
EV_ddf <- as.data.frame(EAST_dates[K])
EV_ddf = EV_ddf %>% 
  dplyr::rename(
    time = 'EAST_dates[K]',
  )
NV_ddf <- as.data.frame(NOR_dates[K])
NV_ddf = NV_ddf %>% 
  dplyr::rename(
    time = 'NOR_dates[K]',
  )
ChlA_ddf <- as.data.frame(ChlA_dates[K])
ChlA_ddf = ChlA_ddf %>% 
  dplyr::rename(
    time = 'ChlA_dates[K]',
  )
ChlA_ddf$time=as.Date(ChlA_ddf$time)

#merge res dataframes with dates
SSHdf<- bind_cols(SSH_ddf,as.data.frame(resSSH))
DENdf<- bind_cols(DEN_ddf,as.data.frame(resDEN))
SALdf<- bind_cols(SAL_ddf,as.data.frame(resSAL))
TEMPdf<- bind_cols(TEMP_ddf,as.data.frame(resTEMP))
EVdf<- bind_cols(EV_ddf,as.data.frame(resEV))
NVdf<- bind_cols(NV_ddf,as.data.frame(resNV))
EKE <- bind_cols(SSH_ddf,as.data.frame(EKE_cm))
ChlAdf <- bind_cols(ChlA_ddf, as.data.frame(resChlA))

#clear memory and increase memory limit size
gc()
rm("SSH_ddf","DEN_ddf","SAL_ddf","TEMP_ddf","EV_ddf","NV_ddf")
rm("points","sp_poly")

#merge the data sets in chunks because of memory issues
tab <- left_join(DayTable, SST4, by = "time") %>%
  left_join(., ChlAdf, by = "time") %>%
  left_join(., SSHdf, by = "time") %>%
  left_join(., DENdf, by = "time") %>%
  left_join(., SALdf, by = "time") %>%
  left_join(., TEMPdf, by = "time") %>%
  left_join(., EVdf, by = "time") %>%
  left_join(., NVdf, by = "time") %>%
  left_join(., EKE, by = "time")

#replacing SST Nan values from satellite with resTEMP (model data)
tab$mean_SST <- ifelse(is.na(tab$mean_SST), tab$resTEMP, tab$mean_SST)
tab = tab[complete.cases(tab[ , 2:4]),]#remove any rows with lat or long as na

#Join sex specific presence in daily hours information
tab <- left_join(tab, SexDayData, by = 'time')
  
#Group by ITS for the general sperm whale model
startDate = tab$time[1]
endDate = tab$time[nrow(tab)]
timeseries = data.frame(date=seq(startDate, endDate, by="days"))
#timeseries$groups = rep(1:(nrow(timeseries)/ITS), times=1, each=ITS)
ITSgroups = rep(1:(floor(nrow(timeseries)/ITS)), times=1, each=ITS)
timeseries$groups = c(ITSgroups,ITSgroups[3180]+1)
TabBinned = left_join(tab,timeseries,by = c("time" = "date"))
TabBinned_Grouped = aggregate(TabBinned[, c(1:length(TabBinned))], list(TabBinned$groups), mean, na.rm = TRUE)
TabBinned_Grouped$Julian = as.numeric(format(TabBinned_Grouped$time,"%j"))
TabBinned_Grouped$Year = as.numeric(format(TabBinned_Grouped$time,"%Y"))

#Group by ITS for Female
if (exists("ITSF")){
  if (ITSF > 1){
    ITSgroupsF = rep(1:(floor(nrow(timeseries)/ITSF)), times=1, each=ITSF)
    timeseriesF = data.frame(date=seq(startDate, endDate, by="days"))
    timeseriesF$groups = c(ITSgroupsF,ITSgroupsF[3180]+1) #FIGURE THIS OUT
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

#Group by ITS for Male


#run GAM
#Test how each covariate should be used (linear, smooth, as.factor())
#empty model for comparison
GAM_empty = gam(HoursNorm ~ 1, data = TabBinned_Grouped, family = tw, method = "REML")

#Julian day
GAM_01a = gam(HoursNorm ~ Julian, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_01b = gam(HoursNorm ~ s(Julian, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model01 = c('empty','01a','01b')
AIC01 = c(AIC(GAM_empty),AIC(GAM_01a),AIC(GAM_01b))
data.frame(rbind(model01,AIC01))
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   2143.75884213975 2130.74080097241 2103.62871708431
#Julian day as a smooth

#Chlorophyll
GAM_02a = gam(HoursNorm ~ resChlA, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_02b = gam(HoursNorm ~ s(resChlA, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model02 = c('empty','02a','02b')
AIC02 = c(AIC(GAM_empty),AIC(GAM_02a),AIC(GAM_02b))
data.frame(rbind(model02,AIC02))
#X1               X2               X3
#model02            empty              02a              02b
#AIC02   2143.75884213975 1032.98688073045 1033.55960791435
#Chlorophyll as linear

#EKE
GAM_03a = gam(HoursNorm ~ EKE_cm, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_03b = gam(HoursNorm ~ s(EKE_cm, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model03 = c('empty','03a','03b')
AIC03 = c(AIC(GAM_empty),AIC(GAM_03a),AIC(GAM_03b))
data.frame(rbind(model03,AIC03))
#X1               X2              X3
#model03            empty              03a             03b
#AIC03   2143.75884213975 2145.70388331462 2143.7678528048
#EKE as a smooth

#Salinity
GAM_04a = gam(HoursNorm ~ resSAL, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_04b = gam(HoursNorm ~ s(resSAL, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model04 = c('empty','04a','04b')
AIC04 = c(AIC(GAM_empty),AIC(GAM_04a),AIC(GAM_04b))
data.frame(rbind(model04,AIC04))
#X1               X2               X3
#model04            empty              04a              04b
#AIC04   2143.75884213975 2144.18447053367 2143.55635615522
#salinity as a smooth

#Mean SST
GAM_05a = gam(HoursNorm ~ mean_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_05b = gam(HoursNorm ~ s(mean_SST, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model05 = c('empty','05a','05b')
AIC05 = c(AIC(GAM_empty),AIC(GAM_05a),AIC(GAM_05b))
data.frame(rbind(model05,AIC05))
#X1               X2               X3
#model05            empty              05a              05b
#AIC05   2143.75884213975 2137.89747496607 2138.79817539631
#mean SST as linear

#SSH
GAM_06a = gam(HoursNorm ~ resSSH, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_06b = gam(HoursNorm ~ s(resSSH, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model06 = c('empty','06a','06b')
AIC06 = c(AIC(GAM_empty),AIC(GAM_06a),AIC(GAM_06b))
data.frame(rbind(model06,AIC06))
#X1               X2              X3
#model06            empty              06a             06b
#AIC06   2143.75884213975 2139.08068098461 2132.5023262054
#SSH as a smooth

#Density
GAM_07a = gam(HoursNorm ~ resDEN, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_07b = gam(HoursNorm ~ s(resDEN, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model07 = c('empty','07a','07b')
AIC07 = c(AIC(GAM_empty),AIC(GAM_07a),AIC(GAM_07b))
data.frame(rbind(model07,AIC07))
#X1               X2               X3
#model07            empty              07a              07b
#AIC07   2143.75884213975 2142.00922461761 2143.75905922045
#density as linear

#Standard deviation of SST
GAM_08a = gam(HoursNorm ~ SD_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_08b = gam(HoursNorm ~ s(SD_SST, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model08 = c('empty','08a','08b')
AIC08 = c(AIC(GAM_empty),AIC(GAM_08a),AIC(GAM_08b))
data.frame(rbind(model08,AIC08))
#X1               X2               X3
#model08            empty              08a              08b
#AIC08   2143.75884213975 1457.77690764797 1456.05186214008
#std dev of SST as a smooth

#Test which covariates we should keep
#Round 1
#Initial model
Full1 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
              , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ resChlA+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
C = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+s(EKE_cm, bs = "cc", k = -1)
           + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1) +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+resChlA+s(EKE_cm, bs = "cc", k = -1)+
            s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR1 = c('Full1','J','C','E','S','MT','H','D','SDT')
AICR1 = c(AIC(Full1),AIC(J),AIC(C),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR1,AICR1))
#X1              X2               X3               X4
#modelR1            Full1               J                C                E
#AICR1   830.227882831821 832.56314687881 1425.18811746228 830.229480278946
#X5               X6               X7               X8
#modelR1                S               MT                H                D
#AICR1   829.379105033852 830.511182156163 835.736161297421 829.268449913528
#X9
#modelR1              SDT
#AICR1   1016.80934718085
#Chlorophyll removed

#Round 2
Full2 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)  +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
         s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN+s(SD_SST, bs="cc", k=-1),
        data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
modelR2 = c('Full2','J','E','S','MT','H','D','SDT')
AICR2 = c(AIC(Full2),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR2,AICR2))
#X1               X2               X3               X4
#modelR2            Full2                J                E                S
#AICR2   1425.18811746228 1446.14871965395 1425.18858469934 1424.29137767899
#X5               X6               X7               X8
#modelR2              MT                H                D              SDT
#AICR2   1426.6336193432 1433.34011165727 1425.80106718315 2096.24862054496
#std dev SST removed

#Round 3
Full3 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)  +s(resSSH, bs="cc", k=-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN,
        data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR3 = c('Full3','J','E','S','MT','H','D')
AICR3 = c(AIC(Full3),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR3,AICR3))
#X1               X2               X3               X4
#modelR3            Full3                J                E                S
#AICR3   1425.18811746228 2128.97009012411 2096.24772996916 2096.24713217825
#X5               X6               X7
#modelR3             MT                H                D
#AICR3   2097.304971053 2107.85285812374 2094.33356083903
#SSH removed because want to keep Julian Day despite high AIC

#Round 4
Full4 = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(HoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          mean_SST +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST 
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR4 = c('Full4','J','E','S','MT','D')
AICR4 = c(AIC(Full4),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(D))
data.frame(rbind(modelR4,AICR4))
#X1               X2               X3              X4
#modelR4            Full4                J                E               S
#AICR4   2107.85285812374 2139.59020978204 2107.86893670364 2106.8874359167
#X5               X6
#modelR4               MT                D
#AICR4   2105.16883787426 2105.44559132224
#keeping round 4 remaining variables


#Final full model
FinalGAM = gam(HoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
                 s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN
               , data = TabBinned_Grouped, family = tw, method = "REML")
summary(FinalGAM)
viz = getViz(FinalGAM)
vizGG = plot(viz,allTerms = T) +
  labs(title = 'Sperm whales (GAM)')+
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
#Test how each covariate should be used (linear, smooth, as.factor())
#empty model for comparison
GAM_empty = gam(FemaleHoursNorm ~ 1, data = TabBinned_Grouped, family = tw, method = "REML")

#Julian day
GAM_01a = gam(FemaleHoursNorm ~ Julian, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_01b = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model01 = c('empty','01a','01b')
AIC01 = c(AIC(GAM_empty),AIC(GAM_01a),AIC(GAM_01b))
data.frame(rbind(model01,AIC01))
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   2665.37688613595 2658.91448130588 2633.63404062137
#Julian day as a smooth

#Chlorophyll
GAM_02a = gam(FemaleHoursNorm ~ resChlA, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_02b = gam(FemaleHoursNorm ~ s(resChlA, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model02 = c('empty','02a','02b')
AIC02 = c(AIC(GAM_empty),AIC(GAM_02a),AIC(GAM_02b))
data.frame(rbind(model02,AIC02))
#X1              X2               X3
#model02            empty             02a              02b
#AIC02   2665.37688613595 1519.2450330591 1518.06630610966
#Chlorophyll as a smooth

#EKE
GAM_03a = gam(FemaleHoursNorm ~ EKE_cm, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_03b = gam(FemaleHoursNorm ~ s(EKE_cm, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model03 = c('empty','03a','03b')
AIC03 = c(AIC(GAM_empty),AIC(GAM_03a),AIC(GAM_03b))
data.frame(rbind(model03,AIC03))
#X1               X2               X3
#model03            empty              03a              03b
#AIC03   2665.37688613595 2667.05372330502 2666.03714500583
#EKE as a smooth

#Salinity
GAM_04a = gam(FemaleHoursNorm ~ resSAL, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_04b = gam(FemaleHoursNorm ~ s(resSAL, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model04 = c('empty','04a','04b')
AIC04 = c(AIC(GAM_empty),AIC(GAM_04a),AIC(GAM_04b))
data.frame(rbind(model04,AIC04))
#X1               X2               X3
#model04            empty              04a              04b
#AIC04   2665.37688613595 2666.16486167496 2666.00015422402
#salinity as a smooth

#Mean SST
GAM_05a = gam(FemaleHoursNorm ~ mean_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_05b = gam(FemaleHoursNorm ~ s(mean_SST, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model05 = c('empty','05a','05b')
AIC05 = c(AIC(GAM_empty),AIC(GAM_05a),AIC(GAM_05b))
data.frame(rbind(model05,AIC05))
#X1               X2               X3
#model05            empty              05a              05b
#AIC05   2665.37688613595 2661.73743691337 2662.81896789659
#Mean SST as linear

#SSH
GAM_06a = gam(FemaleHoursNorm ~ resSSH, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_06b = gam(FemaleHoursNorm ~ s(resSSH, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model06 = c('empty','06a','06b')
AIC06 = c(AIC(GAM_empty),AIC(GAM_06a),AIC(GAM_06b))
data.frame(rbind(model06,AIC06))
#X1              X2               X3
#model06            empty             06a              06b
#AIC06   2665.37688613595 2661.5099005417 2651.75497601258
#SSH as a smooth

#Density
GAM_07a = gam(FemaleHoursNorm ~ resDEN, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_07b = gam(FemaleHoursNorm ~ s(resDEN, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model07 = c('empty','07a','07b')
AIC07 = c(AIC(GAM_empty),AIC(GAM_07a),AIC(GAM_07b))
data.frame(rbind(model07,AIC07))
#X1               X2               X3
#model07            empty              07a              07b
#AIC07   2665.37688613595 2659.61205283345 2665.50705315705
#Density as linear

#Standard deviation of SST
GAM_08a = gam(FemaleHoursNorm ~ SD_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_08b = gam(FemaleHoursNorm ~ s(SD_SST, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model08 = c('empty','08a','08b')
AIC08 = c(AIC(GAM_empty),AIC(GAM_08a),AIC(GAM_08b))
data.frame(rbind(model08,AIC08))
#X1               X2               X3
#model08            empty              08a              08b
#AIC08   2665.37688613595 1728.41431180364 1726.45921466348
#std dev of SST as a smooth

#Test which covariates we should keep for Female Specific GAM
#Round 1
#Initial model
Full1 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k=-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
C = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)
        + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1) +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
            s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR1 = c('Full1','J','C','E','S','MT','H','D','SDT')
AICR1 = c(AIC(Full1),AIC(J),AIC(C),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR1,AICR1))
#X1               X2               X3               X4
#modelR1            Full1                J                C                E
#AICR1   1333.91533897873 1353.23573935008 1704.28395302219 1333.29836015401
#X5              X6               X7               X8
#modelR1                S              MT                H                D
#AICR1   1333.91856692961 1334.3639124352 1333.83177576131 1332.61648595199
#X9
#modelR1              SDT
#AICR1   1498.74673993798
#remove chlorophyll

#Round 2
Full2 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)  +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN+s(SD_SST, bs="cc", k=-1),
        data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
            s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR2 = c('Full2','J','E','S','MT','H','D','SDT')
AICR2 = c(AIC(Full2),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR2,AICR2))
#X1               X2               X3              X4
#modelR2            Full2                J                E               S
#AICR2   1704.28395302219 1718.93507024444 1704.28606643639 1703.8859863176
#X5               X6              X7               X8
#modelR2               MT                H               D              SDT
#AICR2   1702.08751504329 1707.67467869235 1703.0252692576 2625.47011965109
#remove std dev SST

#Round 3
Full3 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)  +s(resSSH, bs="cc", k=-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN,
        data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR3 = c('Full3','J','E','S','MT','H','D')
AICR3 = c(AIC(Full3),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR3,AICR3))
#X1               X2               X3               X4
#modelR3            Full3                J                E                S
#AICR3   2625.47011965109 2647.45354184231 2625.47880711001 2625.47677247455
#X5              X6               X7
#modelR3               MT               H                D
#AICR3   2623.49914510541 2634.7100193122 2623.89766622477
#remove SSH, but keeping Julian Day

#Round 4
Full4 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          mean_SST +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST 
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR4 = c('Full4','J','E','S','MT','D')
AICR4 = c(AIC(Full4),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(D))
data.frame(rbind(modelR4,AICR4))
#                     X1               X2               X3               X4
#modelR4           Full4                J                E                S
#AICR4   2634.7100193122 2660.85670133633 2633.91829869743 2634.78131232744
#X5               X6
#modelR4              MT                D
#AICR4   2635.9519212653 2634.07268645227
#remove mean SST

#Round 5
Full5 = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) +resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(FemaleHoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) +resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
        +resDEN , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR5 = c('Full5','J','E','S','D')
AICR5 = c(AIC(Full5),AIC(J),AIC(E),AIC(S),AIC(D))
data.frame(rbind(modelR5,AICR5))
#X1               X2              X3               X4
#modelR5           Full5                J               E                S
#AICR5   2635.9519212653 2659.77120739666 2635.3997284668 2635.73518860949
#X5
#modelR5                D
#AICR5   2634.76393388367
#keep round 5 remaining variables

#Final Female Sex Specific full model
FinalFemaleGAM = gam(FemaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
                       s(resSAL, bs = "cc", k =-1) +resDEN
                     , data = TabBinned_Grouped, family = tw, method = "REML")
summary(FinalFemaleGAM)
viz = getViz(FinalFemaleGAM)
vizGG = plot(viz,allTerms = T) +
  labs(title = 'Sperm whales (GAM)')+
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

#running Male Sex Specific GAM
#Test how each covariate should be used (linear, smooth, as.factor())
#empty model for comparison
GAM_empty = gam(MaleHoursNorm ~ 1, data = TabBinned_Grouped, family = tw, method = "REML")

#Julian day
GAM_01a = gam(MaleHoursNorm ~ Julian, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_01b = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model01 = c('empty','01a','01b')
AIC01 = c(AIC(GAM_empty),AIC(GAM_01a),AIC(GAM_01b))
data.frame(rbind(model01,AIC01))
#X1               X2               X3
#model01            empty              01a              01b
#AIC01   887.529203151121 878.265022630491 873.011160556566
#Julian Day as a smooth

#Chlorophyll
GAM_02a = gam(MaleHoursNorm ~ resChlA, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_02b = gam(MaleHoursNorm ~ s(resChlA, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model02 = c('empty','02a','02b')
AIC02 = c(AIC(GAM_empty),AIC(GAM_02a),AIC(GAM_02b))
data.frame(rbind(model02,AIC02))
#X1               X2               X3
#model02            empty              02a              02b
#AIC02   887.529203151121 522.333605324521 520.782661963009
#chlorophyll as a smooth

#EKE
GAM_03a = gam(MaleHoursNorm ~ EKE_cm, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_03b = gam(MaleHoursNorm ~ s(EKE_cm, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model03 = c('empty','03a','03b')
AIC03 = c(AIC(GAM_empty),AIC(GAM_03a),AIC(GAM_03b))
data.frame(rbind(model03,AIC03))
#X1              X2               X3
#model03            empty             03a              03b
#AIC03   887.529203151121 889.43488311573 888.033470215873
#EKE as a smooth

#Salinity
GAM_04a = gam(MaleHoursNorm ~ resSAL, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_04b = gam(MaleHoursNorm ~ s(resSAL, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model04 = c('empty','04a','04b')
AIC04 = c(AIC(GAM_empty),AIC(GAM_04a),AIC(GAM_04b))
data.frame(rbind(model04,AIC04))
#X1               X2               X3
#model04            empty              04a              04b
#AIC04   887.529203151121 888.063737473537 885.106938358472
#Salinity as a smooth

#Mean SST
GAM_05a = gam(MaleHoursNorm ~ mean_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_05b = gam(MaleHoursNorm ~ s(mean_SST, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model05 = c('empty','05a','05b')
AIC05 = c(AIC(GAM_empty),AIC(GAM_05a),AIC(GAM_05b))
data.frame(rbind(model05,AIC05))
#X1               X2              X3
#model05            empty              05a             05b
#AIC05   887.529203151121 878.203163215944 878.71817204281
#Mean SST as linear

#SSH
GAM_06a = gam(MaleHoursNorm ~ resSSH, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_06b = gam(MaleHoursNorm ~ s(resSSH, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model06 = c('empty','06a','06b')
AIC06 = c(AIC(GAM_empty),AIC(GAM_06a),AIC(GAM_06b))
data.frame(rbind(model06,AIC06))
#X1               X2               X3
#model06            empty              06a              06b
#AIC06   887.529203151121 888.818483448203 887.531172601541
#SSH as a smooth

#Density
GAM_07a = gam(MaleHoursNorm ~ resDEN, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_07b = gam(MaleHoursNorm ~ s(resDEN, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model07 = c('empty','07a','07b')
AIC07 = c(AIC(GAM_empty),AIC(GAM_07a),AIC(GAM_07b))
data.frame(rbind(model07,AIC07))
#X1               X2              X3
#model07            empty              07a             07b
#AIC07   887.529203151121 882.063834607968 884.37513611068
#density as linear

#Standard deviation of SST
GAM_08a = gam(MaleHoursNorm ~ SD_SST, data = TabBinned_Grouped, family = tw, method = "REML")
GAM_08b = gam(MaleHoursNorm ~ s(SD_SST, bs = "cc", k =-1), data = TabBinned_Grouped, family = tw, method = "REML")

model08 = c('empty','08a','08b')
AIC08 = c(AIC(GAM_empty),AIC(GAM_08a),AIC(GAM_08b))
data.frame(rbind(model08,AIC08))
#X1               X2               X3
#model08            empty              08a              08b
#AIC08   887.529203151121 640.444848766138 638.156146161628
#std dev SST as a smooth

#Test which covariates we should keep for Male Specific GAM
#Round 1
#Initial model
Full1 = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k=-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(MaleHoursNorm ~ s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
C = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)
        + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1) +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(resChlA, bs ="cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
            s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR1 = c('Full1','J','C','E','S','MT','H','D','SDT')
AICR1 = c(AIC(Full1),AIC(J),AIC(C),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR1,AICR1))
#X1               X2               X3               X4
#modelR1            Full1                J                C                E
#AICR1   476.794286020769 477.021708264104 633.124767065043 476.041484985759
#X5               X6               X7               X8
#modelR1                S               MT                H                D
#AICR1   476.794660205568 477.561228754239 476.794329114588 474.924271986226
#X9
#modelR1              SDT
#AICR1   517.850796842214
#remove chlorophyll

#Round 2
Full2 = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(MaleHoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          mean_SST +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)  +s(resSSH, bs="cc", k=-1)+resDEN+s(SD_SST, bs="cc", k=-1)
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN+s(SD_SST, bs="cc", k=-1),
        data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+s(SD_SST, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
SDT = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
            s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
          , data = TabBinned_Grouped, family = tw, method = "REML")
modelR2 = c('Full2','J','E','S','MT','H','D','SDT')
AICR2 = c(AIC(Full2),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D),AIC(SDT))
data.frame(rbind(modelR2,AICR2))
#X1               X2               X3               X4
#modelR2            Full2                J                E                S
#AICR2   633.124767065043 631.022583067217 632.259454185958 631.670886890875
#X5               X6               X7               X8
#modelR2               MT                H                D              SDT
#AICR2   631.141841046402 634.106134166604 630.724512863152 877.860989589944
#remove std dev of SST

#Round 3
Full3 = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
              s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
            , data = TabBinned_Grouped, family = tw, method = "REML")
J = gam(MaleHoursNorm ~ s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
E = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
S = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
        , data = TabBinned_Grouped, family = tw, method = "REML")
MT = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
           s(resSAL, bs = "cc", k =-1)  +s(resSSH, bs="cc", k=-1)+resDEN
         , data = TabBinned_Grouped, family = tw, method = "REML")
H = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +resDEN,
        data = TabBinned_Grouped, family = tw, method = "REML")
D = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
          s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)
        , data = TabBinned_Grouped, family = tw, method = "REML")
modelR3 = c('Full3','J','E','S','MT','H','D')
AICR3 = c(AIC(Full3),AIC(J),AIC(E),AIC(S),AIC(MT),AIC(H),AIC(D))
data.frame(rbind(modelR3,AICR3))
#X1               X2               X3               X4
#modelR3            Full3                J                E                S
#AICR3   877.860989589944 877.351357655869 877.806456886425 877.168759749715
#X5               X6               X7
#modelR3               MT                H                D
#AICR3   876.701137422412 876.777925410025 875.713479764058
#keep round 3 remaining variables

#Final Male Sex Specific full model
FinalMaleGAM = gam(MaleHoursNorm ~ s(Julian, bs = "cc", k=-1)+s(EKE_cm, bs = "cc", k = -1)+
                       s(resSAL, bs = "cc", k =-1) + mean_SST +s(resSSH, bs="cc", k=-1)+resDEN
                     , data = TabBinned_Grouped, family = tw, method = "REML")
summary(FinalMaleGAM)
viz = getViz(FinalFemaleGAM)
vizGG = plot(viz,allTerms = T) +
  labs(title = 'Sperm whales (GAM)')+
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

GAM_trial = gam(HoursNorm ~ s(Julian, bs = "cc", k = -1) + +s(resChlA, bs = "cc", k = -1) +
                  s(EKE_cm, bs = "cc", k = -1) + s(resSAL, bs = "cc", k = -1) + 
                  s(mean_SST, bs = "cc", k = -1) + s(resSSH, bs = "cc", k = -1) + 
                  s(resDEN, bs = "cc", k = -1) + s(SD_SST, bs = "cc", k=-1),
            data = TabBinned_Grouped, family = tw, method = "REML")

summary(GAM_trial)
plot(GAM_trial, pages = 1)
viz = getViz(GAM_trial)
print(plot(viz,allTerms=T),pages=1)

gam.check(GAM_trial)
concurvity(GAM_trial, full= TRUE)
concurvity(GAM_trial, full= FALSE)

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

