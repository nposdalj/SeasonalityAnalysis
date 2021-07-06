#load libraries
library(ncdf4)
library(sp)
library(rgdal)

#define the lat and long of interest
df1 = data.frame("lat" = c(19.29, 19.2, 19.2467, 19.2467), "long" = c(-166.69, -166.69, -166.74, -166.64))
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = 1)))


#loading sperm whale data
site = 'Wake'
saveDir = paste("G:/My Drive/CentralPac_TPWS_metadataReduced/Wake/Seasonality/")#setting the directory
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
envDir = paste("G:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory
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
SST = SST[-1,]

#subset the dataframe based on the area of interest
coordinates(SST) <- c("latitude","longitude")
coords <- over(SST, sp_poly)
SST[coords == 1 & !is.na(coords),]

#SSH anomaly data
filenameStatAll = paste(envDir,"AVISOglobalvars.nc",sep="")#load files as data frame
SSH = nc_open(filenameStatAll)
names(SSH$var)

#zos - SSH
v6=SSH$var[[6]]
SSHvar=ncvar_get(SSH,v6)
SSH_lon=v6$dim[[1]]$vals
SSH_lat=v6$dim[[2]]$vals

#mlotst - density ocean mixed layer thickeness
#so - salinity
#uo - eastward velocity
#vo - northward velocity
#thetao - temperature


#subset the dataframe based on the area of interest
#average the environmental variable based on the ITS over the area of interest
#save standard deviation 



#merge the data sets




#run GAM




#plot GAMs