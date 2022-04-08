## This code plots histograms of julian day for each month (binary presence in each month)
## NP 03242022

#Load libraries
library(ggplot2)
library(lubridate)
library(plyr)
library(plotrix)

#User Directory
GDrive = 'I'
FileDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/",sep="")
fileName = paste(FileDir,'BigModel_gamgeeOutput.RData',sep="")
fileNameSex = paste(FileDir,'BigModel_gamgeeOutput_sexClasses.RData',sep="")
load(fileName)
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",sep="")
PlotDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",sep="")

Sites = c('CB','PT','QN','AB','KOA','BD','KS') #The GOA and BSAI Sites
Regions = c("GOA","BSAI") #GOA AND BSAI Region

#Add month
SiteHourTableB$Month = as.factor(month(SiteHourTableB$date))

#### GENERAL MODELS######
#SITE
#Subset in loop and plot histograms
for (i in 1:length(Sites)){
  site = Sites[[i]]
  #Subset Table
  TempTable = subset(SiteHourTableB, Site == site)#subset the table for the site only
  
  #Average daily values
  TempTableAgg = ddply(TempTable,~date,summarise,mean=mean(PreAbs),sd=sd(PreAbs),N=sum(PreAbs))
  TempTableAgg$Month = as.factor(month(TempTableAgg$date))
  
  #Average monthly values
  TempTableAggMonth = ddply(TempTable,~Month,summarise,Mean=mean(PreAbs),sd=sd(PreAbs),Presence=sum(PreAbs),SEM=std.error(PreAbs),Absence=sum(PreAbs==0))
  TempTableAggMonth$Count = TempTableAggMonth$Presence + TempTableAggMonth$Absence
  TempTableAggMonth$CI = TempTableAggMonth$SEM * qt((1-0.05)/2 + 0.5,TempTableAggMonth$Count-1)
  
  #Calculate the # of 1s and the proportion of 1s
  BinN = sum(TempTable$PreAbs)
  PropBin = sum(TempTable$PreAbs)/length(TempTable)
  
  filename = paste(saveDir,site,'/',site,'_Ns.txt',sep="")
  sink(filename)
  print(paste(nrow(TempTable),'total number of bins'),sep="")
  print(paste(BinN,' bins with 1s'),sep="")
  print(paste(PropBin,' is the proportion of bins with 1s'),sep="")
  sink(file = NULL)
  
  PlotsaveDir = paste(PlotDir,site,sep='')
  
  #Plot histograms
  ggplot(TempTableAggMonth,aes(x=Month,y=Mean))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_',site,'.jpeg',sep=''))

  ggplot(TempTableAggMonth,aes(x=Month,y=Mean))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=Mean-CI, ymax=Mean+CI), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=Count))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=Presence))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_',site,'.jpeg',sep=''))
}

#Region
#Subset in loop and plot histograms
for (i in 1:length(Regions)){
  region = Regions[[i]]
  #Subset Table
  TempTable = subset(SiteHourTableB, Region == region)#subset the table for the site only
  
  #Average daily values
  TempTableAgg = ddply(TempTable,~date,summarise,mean=mean(PreAbs),sd=sd(PreAbs),N=sum(PreAbs))
  TempTableAgg$Month = as.factor(month(TempTableAgg$date))
  
  #Average monthly values
  TempTableAggMonth = ddply(TempTable,~Month,summarise,Mean=mean(PreAbs),sd=sd(PreAbs),Presence=sum(PreAbs),SEM=std.error(PreAbs),Absence=sum(PreAbs==0))
  TempTableAggMonth$Count = TempTableAggMonth$Presence + TempTableAggMonth$Absence
  TempTableAggMonth$CI = TempTableAggMonth$SEM * qt((1-0.05)/2 + 0.5,TempTableAggMonth$Count-1)
  
  #Calculate the # of 1s and the proportion of 1s
  BinN = sum(TempTable$PreAbs)
  PropBin = sum(TempTable$PreAbs)/length(TempTable)
  
  filename = paste(saveDir,region,'/',region,'_Ns.txt',sep="")
  sink(filename)
  print(paste(nrow(TempTable),'total number of bins'),sep="")
  print(paste(BinN,' bins with 1s'),sep="")
  print(paste(PropBin,' is the proportion of bins with 1s'),sep="")
  sink(file = NULL)
  
  PlotsaveDir = paste(PlotDir,region,sep='')
  
  #Plot histograms
  ggplot(TempTableAggMonth,aes(x=Month,y=Mean))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=Mean))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=Mean-CI, ymax=Mean+CI), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=Count))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=Presence))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_',region,'.jpeg',sep=''))
}

#Big Model
#Average daily values
TempTableAgg = ddply(SiteHourTableB,~date,summarise,mean=mean(PreAbs),sd=sd(PreAbs),N=sum(PreAbs))
TempTableAgg$Month = as.factor(month(TempTableAgg$date))
  
#Average monthly values
TempTableAggMonth = ddply(TempTable,~Month,summarise,Mean=mean(PreAbs),sd=sd(PreAbs),Presence=sum(PreAbs),
                          SEM=std.error(PreAbs),Absence=sum(PreAbs==0))
TempTableAggMonth$Count = TempTableAggMonth$Presence + TempTableAggMonth$Absence
TempTableAggMonth$CI = TempTableAggMonth$SEM * qt((1-0.05)/2 + 0.5,TempTableAggMonth$Count-1)
  
#Calculate the # of 1s and the proportion of 1s
BinN = sum(TempTable$PreAbs)
PropBin = sum(TempTable$PreAbs)/length(TempTable)
  
filename = paste(saveDir,'BigModel/BigModel_Ns.txt',sep="")
sink(filename)
print(paste(nrow(TempTable),'total number of bins'),sep="")
print(paste(BinN,' bins with 1s'),sep="")
print(paste(PropBin,' is the proportion of bins with 1s'),sep="")
sink(file = NULL)

PlotsaveDir = paste(PlotDir,'BigModel',sep='')
  
#Plot histograms
ggplot(TempTableAggMonth,aes(x=Month,y=Mean))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=Mean-SEM, ymax=Mean+SEM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Big Moddel",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=Mean))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=Mean-CI, ymax=Mean+CI), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=Count))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Distribution of Raw Data for Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=Presence))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Monthly Binary Presence for Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_BigModel.jpeg',sep=''))

#### SEX SPECIFIC MODELS######
load(fileNameSex)
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",sep="")
#SITE
#Subset in loop and plot histograms
for (i in 1:length(Sites)){
  site = Sites[[i]]
  #Subset Table
  TempTable = subset(SiteHourTableB, Site == site)#subset the table for the site only
  TempTable$Month = month(TempTable$date)
  
  #Average monthly values
  TempTableAggMonth = ddply(TempTable,~Month,summarise,MeanF=mean(PreAbsF),sdF=sd(PreAbsF),PresenceF=sum(PreAbsF),SEMF=std.error(PreAbsF),AbsenceF=sum(PreAbsF==0),
                            MeanJ=mean(PreAbsJ),sdJ=sd(PreAbsJ),PresenceJ=sum(PreAbsJ),SEMJ=std.error(PreAbsJ),AbsenceJ=sum(PreAbsJ==0),
                            MeanM=mean(PreAbsM),sdM=sd(PreAbsM),PresenceM=sum(PreAbsM),SEMM=std.error(PreAbsM),AbsenceM=sum(PreAbsM==0))
  TempTableAggMonth$CountF = TempTableAggMonth$PresenceF + TempTableAggMonth$AbsenceF
  TempTableAggMonth$CountJ = TempTableAggMonth$PresenceJ + TempTableAggMonth$AbsenceJ
  TempTableAggMonth$CountM = TempTableAggMonth$PresenceM + TempTableAggMonth$AbsenceM
  TempTableAggMonth$CIF = TempTableAggMonth$SEMF * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountF-1)
  TempTableAggMonth$CIJ = TempTableAggMonth$SEMJ * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountJ-1)
  TempTableAggMonth$CIM = TempTableAggMonth$SEMM * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountM-1)
  
  #Calculate the # of 1s and the proportion of 1s
  BinNF = sum(TempTable$PreAbsF)
  PropBinF = sum(TempTable$PreAbsF)/length(TempTable$PreAbsF)
  BinNJ = sum(TempTable$PreAbsJ)
  PropBinJ = sum(TempTable$PreAbsJ)/length(TempTable$PreAbsJ)
  BinNM = sum(TempTable$PreAbsM)
  PropBinM = sum(TempTable$PreAbsM)/length(TempTable$PreAbsM)
  
  filename = paste(saveDir,site,'/',site,'_sexSpecific_Ns.txt',sep="")
  sink(filename)
  print("Social Groups")
  print(paste(nrow(TempTable),'total number of bins'),sep="")
  print(paste(BinNF,' bins with 1s'),sep="")
  print(paste(PropBinF,' is the proportion of bins with 1s'),sep="")
  print("Mid-Size")
  print(paste(BinNJ,' bins with 1s'),sep="")
  print(paste(PropBinJ,' is the proportion of bins with 1s'),sep="")
  print("Males")
  print(paste(BinNM,' bins with 1s'),sep="")
  print(paste(PropBinM,' is the proportion of bins with 1s'),sep="")
  sink(file = NULL)
  
  PlotsaveDir = paste(PlotDir,site,sep='')
  
  #Plot histograms for Social Groups
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanF-SEMF, ymax=MeanF+SEMF), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Social Groups at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_SocialGroups_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanF-CIF, ymax=MeanF+CIF), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Social Groups at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_SocialGroups_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=CountF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data for Social Groups at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_SocialGroups_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=PresenceF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence for Social Groups at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_SocialGroups_',site,'.jpeg',sep=''))
  
  #Plot histograms for Mid-Size
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanJ-SEMJ, ymax=MeanJ+SEMJ), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Mid-Size at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_MidSize_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanJ-CIJ, ymax=MeanJ+CIJ), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Mid-Size at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_MidSize_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=CountJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data for Mid-Size at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_MidSize_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=PresenceJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence for Mid-Size at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_MidSize_',site,'.jpeg',sep=''))
  
  #Plot histograms for Males
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanM-SEMM, ymax=MeanM+SEMM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Males at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_Males_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanM-CIM, ymax=MeanM+CIM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Males at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_Males_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=CountM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data for Males at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_Males_',site,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=PresenceM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence for Males at ",site,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_Males_',site,'.jpeg',sep=''))
}

#Region
#Subset in loop and plot histograms
for (i in 1:length(Regions)){
  region = Regions[[i]]
  #Subset Table
  TempTable = subset(SiteHourTableB, Region == region)#subset the table for the site only
  TempTable$Month = month(TempTable$date)
  
  #Average monthly values
  TempTableAggMonth = ddply(TempTable,~Month,summarise,MeanF=mean(PreAbsF),sdF=sd(PreAbsF),PresenceF=sum(PreAbsF),SEMF=std.error(PreAbsF),AbsenceF=sum(PreAbsF==0),
                            MeanJ=mean(PreAbsJ),sdJ=sd(PreAbsJ),PresenceJ=sum(PreAbsJ),SEMJ=std.error(PreAbsJ),AbsenceJ=sum(PreAbsJ==0),
                            MeanM=mean(PreAbsM),sdM=sd(PreAbsM),PresenceM=sum(PreAbsM),SEMM=std.error(PreAbsM),AbsenceM=sum(PreAbsM==0))
  TempTableAggMonth$CountF = TempTableAggMonth$PresenceF + TempTableAggMonth$AbsenceF
  TempTableAggMonth$CountJ = TempTableAggMonth$PresenceJ + TempTableAggMonth$AbsenceJ
  TempTableAggMonth$CountM = TempTableAggMonth$PresenceM + TempTableAggMonth$AbsenceM
  TempTableAggMonth$CIF = TempTableAggMonth$SEMF * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountF-1)
  TempTableAggMonth$CIJ = TempTableAggMonth$SEMJ * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountJ-1)
  TempTableAggMonth$CIM = TempTableAggMonth$SEMM * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountM-1)
  
  #Calculate the # of 1s and the proportion of 1s
  BinNF = sum(TempTable$PreAbsF)
  PropBinF = sum(TempTable$PreAbsF)/length(TempTable$PreAbsF)
  BinNJ = sum(TempTable$PreAbsJ)
  PropBinJ = sum(TempTable$PreAbsJ)/length(TempTable$PreAbsJ)
  BinNM = sum(TempTable$PreAbsM)
  PropBinM = sum(TempTable$PreAbsM)/length(TempTable$PreAbsM)
  
  filename = paste(saveDir,region,'/',region,'_sexSpecific_Ns.txt',sep="")
  sink(filename)
  print("Social Groups")
  print(paste(nrow(TempTable),'total number of bins'),sep="")
  print(paste(BinNF,' bins with 1s'),sep="")
  print(paste(PropBinF,' is the proportion of bins with 1s'),sep="")
  print("Mid-Size")
  print(paste(BinNJ,' bins with 1s'),sep="")
  print(paste(PropBinJ,' is the proportion of bins with 1s'),sep="")
  print("Males")
  print(paste(BinNM,' bins with 1s'),sep="")
  print(paste(PropBinM,' is the proportion of bins with 1s'),sep="")
  sink(file = NULL)
  
  PlotsaveDir = paste(PlotDir,region,sep='')
  
  #Plot histograms for Social Groups
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanF-SEMF, ymax=MeanF+SEMF), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Social Groups at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_SocialGroups_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanF-CIF, ymax=MeanF+CIF), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Social Groups at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_SocialGroups_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=CountF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data for Social Groups at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_SocialGroups_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=PresenceF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence for Social Groups at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_SocialGroups_',region,'.jpeg',sep=''))
  
  #Plot histograms for Mid-Size
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanJ-SEMJ, ymax=MeanJ+SEMJ), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Mid-Size at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_MidSize_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanJ-CIJ, ymax=MeanJ+CIJ), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Mid-Size at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_MidSize_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=CountJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data for Mid-Size at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_MidSize_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=PresenceJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence for Mid-Size at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_MidSize_',region,'.jpeg',sep=''))
  
  #Plot histograms for Males
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanM-SEMM, ymax=MeanM+SEMM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Males at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_Males_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=MeanM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    geom_errorbar(aes(ymin=MeanM-CIM, ymax=MeanM+CIM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
    ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Males at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_Males_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=CountM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Distribution of Raw Data for Males at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_Males_',region,'.jpeg',sep=''))
  
  ggplot(TempTableAggMonth,aes(x=Month,y=PresenceM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
    ggtitle(paste("Monthly Binary Presence for Males at ",region,sep=""))
  ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_Males_',region,'.jpeg',sep=''))
}

#Big Model
#Average monthly values
TempTableAggMonth = ddply(TempTable,~Month,summarise,MeanF=mean(PreAbsF),sdF=sd(PreAbsF),PresenceF=sum(PreAbsF),SEMF=std.error(PreAbsF),AbsenceF=sum(PreAbsF==0),
                          MeanJ=mean(PreAbsJ),sdJ=sd(PreAbsJ),PresenceJ=sum(PreAbsJ),SEMJ=std.error(PreAbsJ),AbsenceJ=sum(PreAbsJ==0),
                          MeanM=mean(PreAbsM),sdM=sd(PreAbsM),PresenceM=sum(PreAbsM),SEMM=std.error(PreAbsM),AbsenceM=sum(PreAbsM==0))
TempTableAggMonth$CountF = TempTableAggMonth$PresenceF + TempTableAggMonth$AbsenceF
TempTableAggMonth$CountJ = TempTableAggMonth$PresenceJ + TempTableAggMonth$AbsenceJ
TempTableAggMonth$CountM = TempTableAggMonth$PresenceM + TempTableAggMonth$AbsenceM
TempTableAggMonth$CIF = TempTableAggMonth$SEMF * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountF-1)
TempTableAggMonth$CIJ = TempTableAggMonth$SEMJ * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountJ-1)
TempTableAggMonth$CIM = TempTableAggMonth$SEMM * qt((1-0.05)/2 + 0.5,TempTableAggMonth$CountM-1)

#Calculate the # of 1s and the proportion of 1s
BinNF = sum(TempTable$PreAbsF)
PropBinF = sum(TempTable$PreAbsF)/length(TempTable$PreAbsF)
BinNJ = sum(TempTable$PreAbsJ)
PropBinJ = sum(TempTable$PreAbsJ)/length(TempTable$PreAbsJ)
BinNM = sum(TempTable$PreAbsM)
PropBinM = sum(TempTable$PreAbsM)/length(TempTable$PreAbsM)

filename = paste(saveDir,'BigModel/BigModel_sexSpecific_Ns.txt',sep="")
sink(filename)
print("Social Groups")
print(paste(nrow(TempTable),'total number of bins'),sep="")
print(paste(BinNF,' bins with 1s'),sep="")
print(paste(PropBinF,' is the proportion of bins with 1s'),sep="")
print("Mid-Size")
print(paste(BinNJ,' bins with 1s'),sep="")
print(paste(PropBinJ,' is the proportion of bins with 1s'),sep="")
print("Males")
print(paste(BinNM,' bins with 1s'),sep="")
print(paste(PropBinM,' is the proportion of bins with 1s'),sep="")
sink(file = NULL)

PlotsaveDir = paste(PlotDir,'BigModel',sep='')

#Plot histograms for Social Groups
ggplot(TempTableAggMonth,aes(x=Month,y=MeanF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=MeanF-SEMF, ymax=MeanF+SEMF), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Social Groups in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_SocialGroups_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=MeanF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=MeanF-CIF, ymax=MeanF+CIF), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Social Groups in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_SocialGroups_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=CountF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Distribution of Raw Data for Social Groups in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_SocialGroups_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=PresenceF))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Monthly Binary Presence for Social Groups in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_SocialGroups_BigModel.jpeg',sep=''))

#Plot histograms for Mid-Size
ggplot(TempTableAggMonth,aes(x=Month,y=MeanJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=MeanJ-SEMJ, ymax=MeanJ+SEMJ), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Mid-Size in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_MidSize_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=MeanJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=MeanJ-CIJ, ymax=MeanJ+CIJ), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Mid-Size in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_MidSize_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=CountJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Distribution of Raw Data for Mid-Size in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_MidSize_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=PresenceJ))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Monthly Binary Presence for Mid-Size in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_MidSize_BigModel.jpeg',sep=''))

#Plot histograms for Males
ggplot(TempTableAggMonth,aes(x=Month,y=MeanM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=MeanM-SEMM, ymax=MeanM+SEMM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Standard Error of the Mean for Males in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_SEM_Males_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=MeanM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  geom_errorbar(aes(ymin=MeanM-CIM, ymax=MeanM+CIM), width=0.4, colour="orange", alpha=0.9, size=1.3)+
  ggtitle(paste("Mean Monthly Binary Presence with Confidence Intervals for Males in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MeanMonthlyBinaryPresence_CI_Males_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=CountM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Distribution of Raw Data for Males in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/RawDataDistribution_Males_BigModel.jpeg',sep=''))

ggplot(TempTableAggMonth,aes(x=Month,y=PresenceM))+geom_bar(stat="identity",fill="skyblue",alpha=0.7)+
  ggtitle(paste("Monthly Binary Presence for Males in Big Model",sep=""))
ggsave(file = paste(PlotsaveDir,'/MonthlyBinaryPresence_Males_Big Model.jpeg',sep=''))