## This code plots histograms of julian day for each month (binary presence in each month)
## NP 03242022

#Load libraries
library(ggplot2)
library(lubridate)
library(plyr)

#User Directory
GDrive = 'H'
FileDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/",sep="")
fileName = paste(FileDir,'BigModel_gamgeeOutput.RData',sep="")
load(fileName)
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",sep="")

Sites = c('CB','PT','QN','AB','KOA','BD','KS') #The GOA and BSAI Sites
Regions = c("GOA","BSAI") #GOA AND BSAI Region

#Add month
SiteHourTableB$Month = as.factor(month(SiteHourTableB$date))

#SITE
#Subset in loop and plot histograms
for (i in 1:length(Sites)){
  site = Sites[[1]]
  #Subset Table
  TempTable = subset(SiteHourTableB, Site == site)#subset the table for the site only
  
  #Average daily values
  TempTableAgg = ddply(TempTable,~date,summarise,mean=mean(PreAbs),sd=sd(PreAbs))
  TempTableAgg$Month = as.factor(month(TempTableAgg$date))
  
  #Average monthly values
  TempTableAggMonth = ddply(TempTable,~Month,summarise,mean=mean(PreAbs),sd=sd(PreAbs),sum=sum(PreAbs))
  
  #Calculate the # of 1s and the proportion of 1s
  BinN = sum(TempTable$PreAbs)
  PropBin = sum(TempTable$PreAbs)/length(TempTable)
  
  filename = paste(saveDir,site,'/',site,'_Ns.txt',sep="")
  sink(filename)
  paste(BinN,' bins with 1s')
  paste(PropBin,' is the proportion of bins with 1s')
  sink(file = NULL)
  
  #Plot violin plots
  ggplot(TempTableAgg)+geom_bar(aes(x=Month,y=mean),stat="identity",fill="skyblue",alpha=0.7)
  
  ggplot(TempTableAggMonth)+geom_bar(aes(x=Month,y=mean),stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(x=Month, ymin=mean-sd, ymax=mean+sd), width=0.4, colour="orange", alpha=0.9, size=1.3)
  


}
