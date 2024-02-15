#Load libraries
library("eulerr")
library("tidyverse")
library("dplyr")

#load data
GDrive =  'L'
Region = c('GofAK')
dir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites",sep="")

#Site Names
SiteNames = c('CB','PT','QN','BD','AB','KOA','KS')
#SiteNames = c('CA','CCE','CORC','DCPP01C','GI','HOKE','PS1','PS2','QC')
#SiteNames = c('CSM','Equator','Kauai','King','Kona','LSM','Pagan','Palmyra','PHR','Saipan','Tinian','Wake')
#SiteNames = c('HZ','OC','NC','BC','WC','NFC','HAT_A','HAT_B','GS','BP','BS','JAX')

saveDir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/Plots/",sep="")

#General Data
fileName1 = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/AllSitesGrouped_GAMGEE_ROW.csv",sep="")#setting the directory
DayTable = read.csv(fileName1) #no effort days deleted
#DayTable$Region = 'CCE'
DayTable = na.omit(DayTable)
DayTable$tbin = as.Date(DayTable$tbin)

#Sex Specific Data
fileName2 = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/AllSitesGrouped_GAMGEE_ROW_sexClasses.csv",sep="")#setting the directory
SexDayTable = read.csv(fileName2) #no effort days deleted
#SexDayTable$Region = 'CCE'
SexDayTable = na.omit(SexDayTable)
SexDayTable$tbin = as.Date(SexDayTable$tbin)

#Make new columns to find ratio
SexDayTable$F = SexDayTable$PreAbsF
SexDayTable$J = SexDayTable$PreAbsJ
SexDayTable$M = SexDayTable$PreAbsM
SexDayTable$FJ = SexDayTable$F + SexDayTable$J
SexDayTable$JM = SexDayTable$M + SexDayTable$J
SexDayTable$FM = SexDayTable$F + SexDayTable$M
SexDayTable$FJM = SexDayTable$F + SexDayTable$J + SexDayTable$M

SexDayTable$FJM = replace(SexDayTable$FJM, which(SexDayTable$FJM <3 ), 0) #delete rows that don't have all sexes
SexDayTable$FJM = replace(SexDayTable$FJM, which(SexDayTable$FJM == 3), 1) #only keep rows that have both sexes and make it equal to 1

SexDayTable$FJ = replace(SexDayTable$FJ, which(SexDayTable$FJ <2), 0) #delete rows that don't have both sexes
SexDayTable$FJ = replace(SexDayTable$FJ, which(SexDayTable$FJM == 1), 0) #delete rows that have all sexes
SexDayTable$FJ = replace(SexDayTable$FJ, which(SexDayTable$FJ == 2), 1) #only keep rows that have both sexes and make it equal to 1

SexDayTable$JM = replace(SexDayTable$JM, which(SexDayTable$JM <2 ), 0) #delete rows that don't have both sexes
SexDayTable$JM = replace(SexDayTable$JM, which(SexDayTable$FJM == 1), 0) #delete rows that have all sexes
SexDayTable$JM = replace(SexDayTable$JM, which(SexDayTable$JM == 2), 1) #only keep rows that have both sexes and make it equal to 1

SexDayTable$FM = replace(SexDayTable$FM, which(SexDayTable$FM <2 ), 0) #delete rows that don't have both sexes
SexDayTable$FM = replace(SexDayTable$FM, which(SexDayTable$FJM == 1), 0) #delete rows that have all sexes
SexDayTable$FM = replace(SexDayTable$FM, which(SexDayTable$FM == 2), 1) #only keep rows that have both sexes and make it equal to 1

SexDayTable$F = replace(SexDayTable$F, which(SexDayTable$FJ ==1 | SexDayTable$FJM ==1 | SexDayTable$FM ==1),0) #delete F only rows when it's being accounted for in another group
SexDayTable$J = replace(SexDayTable$J, which(SexDayTable$FJ ==1 | SexDayTable$JM ==1 | SexDayTable$FJM ==1),0) #delete J only rows when it's being accounted for in another group
SexDayTable$M = replace(SexDayTable$M, which(SexDayTable$JM ==1 | SexDayTable$FM ==1 | SexDayTable$FJM ==1),0) #delete M only rows when it's being accounted for in another group

#loop through each site and pull out relevant data
for (i in 1:length(SiteNames)){
  site = SiteNames[i]
  SiteDayTable = dplyr::filter(SexDayTable, grepl(site,Site))
  SiteDayTableGen = dplyr::filter(DayTable, grepl(site,Site))
  
  #Find Total Days and Days with Presence in General
  TotalDays = nrow(SiteDayTableGen)
  PresentDays = apply(SiteDayTable,2,function(x) sum(x > 0))
    if (i == 1){
    SumTable = data.frame("site" = site, "TotalDays" = TotalDays, "PresentDays" = sum(SiteDayTableGen$PreAbs))
  }else{
    SumTableTemp = data.frame("site" = site, "TotalDays" = TotalDays, "PresentDays" = sum(SiteDayTableGen$PreAbs))
    SumTable = rbind(SumTable, SumTableTemp)
  }
  
  #Find Ratio of F/M/J Days
  if (i == 1){
  SiteRatio = c(F = sum(SiteDayTable$F), J = sum(SiteDayTable$J), M = sum(SiteDayTable$M), "F&J" = sum(SiteDayTable$FJ),
                          "F&M" = sum(SiteDayTable$FM), "J&M" = sum(SiteDayTable$JM), "F&J&M" = sum(SiteDayTable$FJM))
  SuMSiteRatio = sum(SiteRatio)
  SiteRatioPercent = ceiling((SiteRatio/SuMSiteRatio) * 100)
  SiteRatioPercent = list(SiteRatioPercent)
  }else{
    SiteRatioTemp = c(F = sum(SiteDayTable$F), J = sum(SiteDayTable$J), M = sum(SiteDayTable$M), "F&J" = sum(SiteDayTable$FJ),
                               "F&M" = sum(SiteDayTable$FM), "J&M" = sum(SiteDayTable$JM), "F&J&M" = sum(SiteDayTable$FJM))
    SuMSiteRatioTEMP = sum(SiteRatioTemp)
    SiteRatioPercentTEMP = ceiling((SiteRatioTemp/SuMSiteRatioTEMP) * 100)
    SiteRatioPercent = append(SiteRatioPercent,list(SiteRatioPercentTEMP))
  }
}

#For some reason, I need to run this manually
# Pie Charts --------------------------------------------------------------
for (i in 1:length(SiteNames)){
  site = SiteNames[i]
  title = paste(saveDir,"/",site,"/PropVennDaily", site,".pdf",sep="")
  pdf(title)
  plot(euler(SiteRatioPercent[[i]]),fills = c('#66c2a5','#fc8d62','#8da0cb'), labels = FALSE, quantities = list(cex=3),
       main = paste("Proportional Venn Diagram for ",site,sep = ""))
  print(paste('Plot for ',site,' Saved',sep=""))
  dev.off()
}

# Relative Bars & Effort -----------------------------------------------------------

#Relative Effort
MaxEffort = max(SumTable$TotalDays)
SumTable$RelativeEffort = SumTable$TotalDays/MaxEffort

#Relative Presence
MaxPresence = max(SumTable$PresentDays)
SumTable$Relative_toEachOther_Presence = SumTable$PresentDays/MaxPresence
SumTable$Relative_toEffort_Presence = SumTable$PresentDays/SumTable$TotalDays

#Stacked bar plot for relative effort
title = paste(saveDir,"/RelativeEffortDays",".pdf",sep="")
pdf(title)
p = ggplot(SumTable, aes(x=site, y=RelativeEffort)) +
  geom_bar(stat="identity", fill="lightgrey")+
  theme_minimal()
p+ theme(
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "grey")
)
dev.off()

#Stacked bar plot for relative to each other presence
title = paste(saveDir,"/Relative_toEachOther_PresenceDays",".pdf",sep="")
pdf(title)
p = ggplot(SumTable, aes(x=site, y=Relative_toEachOther_Presence)) +
  geom_bar(stat="identity", fill="darkgray")+
  theme_minimal()
p+ theme(
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "grey")
)
dev.off()

#Stacked bar plot for relative to effort presence
title = paste(saveDir,"/Relative_toEffort_PresenceDays",".pdf",sep="")
pdf(title)
p = ggplot(SumTable, aes(x=site, y=Relative_toEffort_Presence)) +
  geom_bar(stat="identity", fill="darkgray")+
  theme_minimal()
p+ theme(
  # Remove panel border
  panel.border = element_blank(),  
  # Remove panel grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Remove panel background
  panel.background = element_blank(),
  # Add axis line
  axis.line = element_line(colour = "grey")
)
dev.off()
