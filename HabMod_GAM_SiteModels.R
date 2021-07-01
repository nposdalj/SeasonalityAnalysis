#load libraries

#define the lat and long of interest

#loading sperm whale data
site = 'Wake'
saveDir = paste("G:/My Drive/CentralPac_TPWS_metadataReduced/Wake/Seasonality/")#setting the directory
#load data from StatisicalAnalysis_All
filenameStatAll = paste(saveDir,site,"_GroupedDay.csv",sep="")
GroupedDay = read.csv(filenameStatAll) #load files as data frame



#loading the environmental data
envDir = paste("G:/My Drive/Gaia_EnvironmentalData/CentralPac/")#setting the directory
filenameStatAll = paste(envDir,"Chlorophyll.csv",sep="")#load files as data frame
chl = read.csv(filenameStatAll)
#subset the dataframe based on the area of interest
#average the environmental variable based on the ITS over the area of interest
  #save standard deviation 





#merge the data sets




#run GAM




#plot GAMs