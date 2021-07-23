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
startTime = "2010-03-05" #this should be formatted like this: 2010-03-05
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
  ggtitle(paste("Daily Chl on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
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