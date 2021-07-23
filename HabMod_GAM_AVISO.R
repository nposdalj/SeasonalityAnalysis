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

#define the lat and long of interest
#df1 = data.frame("lat" = c(19.29, 19.2, 19.2467, 19.2467), "long" = c(-166.69, -166.69, -166.74, -166.64)) #Wake
df1 = data.frame("lat" = c(15.36, 15.27, 15.3186, 15.3186), "long" = c(145.46, 145.46, 145.51, 145.41)) #Saipan
#df1 = data.frame("lat" = c(15.08, 14.99, 15.0387, 15.0387), "long" = C(145.75, 145.75, 145.8, 145.7))#Tinian
#df1 = data.frame("lat" = c(29.185, 29.1, 29.1410, 29.1410), "long" = C(-118.26, -118.26, -118.31, -118.21))#GI
#df1 = data.frame("lat" = c(31.79, 31.7, 31.747, 31.747), "long" = C(-121.38, -121.38, -121.43, -121.33))#CORC
#df1 = data.frame("lat" = c(52.4, 52.31, 52.3547, 52.3547), "long" = C(-175.635, -175.635, -175.71, -175.56))#BD
#df1 = data.frame("lat" = c(47.54, 47.45, 47.4936, 47.4936), "long" = C(-125.378, -125.378, -125.44, -125.31))#QC
#df1 = data.frame("lat" = c(56.385, 56.295, 56.34, 56.34), "long" = C(-145.182, -145.182, -145.26, -145.1))#QN
#df1 = data.frame("lat" = c(56.29, 56.2, 56.2434, 56.2434), "long" = C(-142.75, -142.75, -142.83, -142.67))#PT
#df1 = data.frame("lat" = c(58.71, 58.62, 58.6668, 58.6668), "long" = C(-148.0034, -148.0034, -148.12, -147.94))#CB

#define the start and end of the datame 
startTime = "2010-03-05" #this should be formatted like this: 2010-03-05 00:05:00 PDT
endTime = "2019-02-02" 

#ITS
ITS = 4

#loading the environmental data
envDir = paste("O:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory)

#spatial polygon for area of interest
ch <- chull(df1$long, df1$lat)
coords <- df1[c(ch, ch[1]), ]#creating convex hull
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = 1)))#converting convex hull to spatial polygon

#loading sperm whale data
site = 'SAP'
saveDir = paste("O:/My Drive/CentralPac_TPWS_metadataReduced/Saipan/Seasonality/")#setting the directory

#load data from StatisicalAnalysis_All
filenameStatAll = paste(saveDir,site,"_Day.csv",sep="")
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
  ggtitle(paste("Daily SSH on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting timeseries
I=which(SSH_lon>=min(df1$long) & SSH_lon<= max(df1$long)) #only extract the region we care about
J=which(SSH_lat>=min(df1$lat) & SSH_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(SSH_dates>= startTime & SSH_dates<= endTime) #extract only the dates we care about
SSH2=SSHvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

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
  rename(
    time = 'SSH_dates[K]',
  )
DEN_ddf <- as.data.frame(DEN_dates[K])
DEN_ddf = DEN_ddf %>% 
  rename(
    time = 'DEN_dates[K]',
  )
SAL_ddf <- as.data.frame(SAL_dates[K])
SAL_ddf = SAL_ddf %>% 
  rename(
    time = 'SAL_dates[K]',
  )
TEMP_ddf <- as.data.frame(TEMP_dates[K])
TEMP_ddf = TEMP_ddf %>% 
  rename(
    time = 'TEMP_dates[K]',
  )
EV_ddf <- as.data.frame(EAST_dates[K])
EV_ddf = EV_ddf %>% 
  rename(
    time = 'EAST_dates[K]',
  )
NV_ddf <- as.data.frame(NOR_dates[K])
NV_ddf = NV_ddf %>% 
  rename(
    time = 'NOR_dates[K]',
  )
ChlA_ddf <- as.data.frame(ChlA_dates[K])
ChlA_ddf = ChlA_ddf %>% 
  rename(
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
