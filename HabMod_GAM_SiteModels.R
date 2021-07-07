#load libraries
library(ncdf4)
library(sp)
library(rgdal)

#define the lat and long of interest
#df1 = data.frame("lat" = c(19.29, 19.2, 19.2467, 19.2467), "long" = c(-166.69, -166.69, -166.74, -166.64)) #Wake
df1 = data.frame("lat" = c(15.36, 15.27, 15.3186, 15.3186), "long" = c(145.46, 145.46, 145.51, 145.41)) #Saipan
#Tinian
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

#SST data
filenameStatAll = paste(envDir,"SST_WAKE.csv",sep="")#load files as data frame
SST = read.csv(filenameStatAll)
SST = SST[-1,] #delete first row
SST$latitude = as.numeric(SST$latitude)
SST$longitude = as.numeric(SST$longitude)
SST = SST[complete.cases(SST[ , 2:3]),]#remove any rows with lat or long as na

#data exploration
#plot time series
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