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

#spatial polygon for area of interest
ch <- chull(df1$long, df1$lat)
coords <- df1[c(ch, ch[1]), ]#creating convex hull
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = 1)))#converting convex hull to spatial polygon

#loading sperm whale data
site = 'SAP'
saveDir = paste("I:/My Drive/CentralPac_TPWS_metadataReduced/Saipan/Seasonality/")#setting the directory

#load data from StatisicalAnalysis_All
filenameStatAll = paste(saveDir,site,"_GroupedDay.csv",sep="")
GroupedDay = read.csv(filenameStatAll) #load files as data frame

#load sex specific data
#Females
filename_GDF = paste(saveDir,site,"_GroupedDayF.csv",sep="")
GroupedDayF = read.csv(filename_GDF) #load files as data frame
#Juveniles
filename_GDJ = paste(saveDir,site,"_GroupedDayJ.csv",sep="")
GroupedDayJ = read.csv(filename_GDJ) #load files as data frame
#Males
filename_GDM = paste(saveDir,site,"_GroupedDayM.csv",sep="")
GroupedDayM = read.csv(filename_GDM) #load files as data frame

#clear memory 
gc()

#loading the environmental data
envDir = paste("I:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory

#chlorophyll data
filenameStatAll = paste(envDir,"Chl2.csv",sep="")#load files as data frame
chl = read.csv(filenameStatAll)
chl = chl[-1,] #delete first row
chl$latitude = as.numeric(chl$latitude)
chl$longitude = as.numeric(chl$longitude)
chl = chl[complete.cases(chl[ , 2:3]),]#remove any rows with lat or long as na


#subset the dataframe based on the area of interest
coordinates(chl) <- c("latitude", "longitude")
coords <- over(chl, sp_poly)
chl2 <- chl[coords == 1 & !is.na(coords),]
chl2$chlorophyll = as.numeric(chl2$chlorophyll)#converting SST from character to numeric
chl2$time = as.Date(chl2$time)#converting time from character to date
chl3 = as.data.frame(chl2)#converting SPDF back to DF
chl3$chlorophyll[chl3$chlorophyll < 0] <- NA #making anything<0 NA

#average the environmental variable based on the ITS over the area of interest
chl4 = chl3 %>%
  mutate(time = floor_date(time)) %>%
  group_by(time) %>%
  summarize(mean_chl = mean(chlorophyll), SD_chl = sd(chlorophyll)) #finding daily mean

#data exploration
mean(chl2$chlorophyll, na.rm = TRUE)#finding overall mean
#save standard deviation
sd(chl2$chlorophyll, na.rm = TRUE)
#SST histogram
hist(chl2$chlorophyll)

#plot time series
plot(chl4$time, chl4$mean_chl)#exploratory plot
title1 = paste(site,"Chlorophyll Plot")
ggplot(chl4, aes(x=time,y=mean_chl))+
  ggtitle(title1)+
  labs(y="Mean Chlorophyll (mg m-3)",x="Time (days)")+
  geom_line()+
  geom_point()

rm(chl)
rm(chl2)
rm(chl3)