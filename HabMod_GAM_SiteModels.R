#load libraries
library(ncdf4)
library(sp)
library(rgdal)
library(httr)
library(sf)
library(dplyr)
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

ch <- chull(df1$long, df1$lat)
coords <- df1[c(ch, ch[1]), ]
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = 1)))


#loading sperm whale data
site = 'Wake'
saveDir = paste("O:/My Drive/CentralPac_TPWS_metadataReduced/Wake/Seasonality/")#setting the directory
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
envDir = paste("O:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory
#chlorophyll data
filenameStatAll = paste(envDir,"Chlorophyll.csv",sep="")#load files as data frame
chl = read.csv(filenameStatAll)
#subset the dataframe based on the area of interest
coordinates(chl) <- ~lat + long
coords <- over(chl, sp_poly)
chl[coords == 1 & !is.na(coords),]
#average the environmental variable based on the ITS over the area of interest
  #save standard deviation

#clear memory 
gc()

#SST data
filenameStatAll = paste(envDir,"SST_SAPTIN2.csv",sep="")#load files as data frame
SST = read.csv(filenameStatAll)
SST = SST[-1,] #delete first row
SST$latitude = as.numeric(SST$latitude)
SST$longitude = as.numeric(SST$longitude)
SST = SST[complete.cases(SST[ , 2:3]),]#remove any rows with lat or long as na

#data exploration
#plot time series
  #SST 
#load files
SST = nc_open("SST_SAPTIN.nc")
names(SST$var)
v1=SST$var[[1]]
SSTvar=ncvar_get(SST,v1)
SST_lon=v1$dim[[1]]$vals
SST_lat=v1$dim[[2]]$vals
dates=as.POSIXlt(v1$dim[[3]]$vals,origin='2010-01-01',tz='GMT')

#loading the world
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#setting color breaks
h=hist(SSTvar[,,1],100,plot=FALSE)
breaks=h$breaks
n=length(breaks)-1
jet.colors <-colorRampPalette(c("blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
c=jet.colors(n)

#plotting in base R
layout(matrix(c(1,2,3,0,4,0), nrow=1, ncol=2), widths=c(5,1), heights=4) 
layout.show(2) 
par(mar=c(3,3,3,1))
r = raster(t(SSTvar[,,1]),xmn = min(SST_lon),xmx = max(SST_lon),ymn=min(SST_lat),ymx=max(SST_lat))
image(r,col=c,breaks=breaks,xlab='',ylab='',axes=TRUE,xlim = c(min(SST_lon),max(SST_lon)),xaxs='i',asp=0, main=paste("Monthly SST", dates[1]))
#points(-66.35, rep(41.06165),pch=20,cex=2) #do we need these points
#points(-76, rep(33.6699),pch=20,cex=2)

#adding color scale
par(mar=c(3,1,3,3))
source('scale.R') 
image.scale(sst[,,1], col=c, breaks=breaks, horiz=FALSE, yaxt="n",xlab='',ylab='',main='SST') 
axis(4, las=1) 
box()

#plotting time series SAPTIN 
I=which(SST_lon>=-144.45 & SST_lon<=-146,76) #change lon to SST_lon values to match ours, use max and min function
J=which(SST_lat>=14.03 & SST_lat<=16.3) #change ""
sst2=SSTvar[I,J,]

n=dim(sst2)[3] 

res=rep(NA,n) 
for (i in 1:n) 
  res[i]=mean(sst2[,,i],na.rm=TRUE) 

plot(1:n,res,axes=FALSE,type='o',pch=20,xlab='',ylab='SST (ºC)') 
axis(2) 
axis(1,1:n,format(dates,'%m')) 
box()

#plotting time series WAKE #add code to load SST_WAKE 
I=which(SST_lon>=-76.25 & SST_lon<=-75.75) #change lon to SST_lon values to match ours, use max and min function
J=which(SST_lat>=33.41991667 & SST_lat<=33.91991667) #change ""
sst2=SSTvar[I,J,] 

n=dim(sst2)[3] 

res=rep(NA,n) 
for (i in 1:n) 
  res[i]=mean(sst2[,,i],na.rm=TRUE) 

plot(1:n,res,axes=FALSE,type='o',pch=20,xlab='',ylab='SST (ºC)') 
axis(2) 
axis(1,1:n,format(dates,'%m')) 
box()
#find min and max

#subset the dataframe based on the area of interest
coordinates(SST) <- c("latitude","longitude")
coords <- over(SST, sp_poly)
SST = SST[coords == 1 & !is.na(coords),]

#SSH anomaly data
filenameStatAll = paste(envDir,"AVISOglobalvars.nc",sep="")#load files as data frame
SSH = nc_open(filenameStatAll)
names(SSH$var)

#zos - SSH
v6=SSH$var[[6]]
SSHvar=ncvar_get(SSH,v6)
SSH_lon=v6$dim[[1]]$vals
SSH_lat=v6$dim[[2]]$vals

#mlotst - density ocean mixed layer thickness
v1=DEN$var[[1]]
DENvar=ncvar_get(DEN,v1)
DEN_lon=v1$dim[[1]]$vals
DEN_lat=v1$dim[[2]]$vals

#so - salinity
v2=SAL$var[[2]]
SALvar=ncvar_get(SAL,v2)
SAL_lon=v2$dim[[1]]$vals
SAL_lat=v2$dim[[2]]$vals

#thetao - temperature
v3=TEMP$var[[3]]
TEMPvar=ncvar_get(TEMP,v3)
TEMP_lon=v3$dim[[1]]$vals
TEMP_lat=v3$dim[[2]]$vals

#uo - eastward velocity
v4=EASTV$var[[4]]
EASTVvar=ncvar_get(EASTV,v4)
EASTV_lon=v4$dim[[1]]$vals
EASTV_lat=v4$dim[[2]]$vals

#vo - northward velocity
v5=NORV$var[[5]]
NORVvar=ncvar_get(NORV,v5)
NORV_lon=v5$dim[[1]]$vals
NORV_lat=v5$dim[[2]]$vals



#subset the dataframe based on the area of interest
#average the environmental variable based on the ITS over the area of interest
#save standard deviation 



#merge the data sets




#run GAM




#plot GAMs

