#load libraries
library(ncdf4)

#define the lat and long of interest

#loading sperm whale data
site = 'Wake'
saveDir = paste("D:/My Drive/CentralPac_TPWS_metadataReduced/Wake/Seasonality/")#setting the directory
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

#loading the environmental data
envDir = paste("D:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory
#chlorophyll data
filenameStatAll = paste(envDir,"Chlorophyll.csv",sep="")#load files as data frame
chl = read.csv(filenameStatAll)
#subset the dataframe based on the area of interest
#average the environmental variable based on the ITS over the area of interest
  #save standard deviation 

#SSH anomaly data
filenameStatAll = paste(envDir,"SSHAnomaly.nc",sep="")#load files as data frame
SSH = nc_open(filenameStatAll)
names(SSH$var)
v1=SSH$var[[1]]
SSHvar=ncvar_get(SSH,v1)
SSH_lon=v1$dim[[1]]$vals
SSH_lat=v1$dim[[2]]$vals

#subset the dataframe based on the area of interest
#average the environmental variable based on the ITS over the area of interest
#save standard deviation 



#merge the data sets




#run GAM




#plot GAMs