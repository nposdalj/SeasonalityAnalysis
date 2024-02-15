#Load libraries
library("eulerr")
library("tidyverse")
library("dplyr")

#load data
GDrive =  'L'
Region = c('WAT')
dir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites",sep="")

#Site Names
#SiteNames = c('CB','PT','QN','BD','AB','KOA','KS')
#SiteNames = c('CA','CCE','CORC','DCPP01C','GI','HOKE','PS1','PS2','QC')
#SiteNames = c('CSM','Equator','Kauai','King','Kona','LSM','Pagan','Palmyra','PHR','Saipan','Tinian','Wake')
SiteNames = c('HZ','OC','NC','BC','WC','NFC','HAT_A','HAT_B','GS','BP','BS','JAX')

saveDir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/Plots/",sep="")

#General Data
fileName1 = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="")#setting the directory
HourTable = read.csv(fileName1) #no effort days deleted
#HourTable$Region = 'CCE'
HourTable = na.omit(HourTable)
HourTable$tbin = as.Date(HourTable$tbin)

#Sex Specific Data
fileName2 = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="")#setting the directory
SexHourTable = read.csv(fileName2) #no effort days deleted
#SexHourTable$Region = 'CCE'
SexHourTable = na.omit(SexHourTable)
SexHourTable$tbin = as.Date(SexHourTable$tbin)

#Make new columns to find ratio
SexHourTable$F = SexHourTable$PreAbsF
SexHourTable$J = SexHourTable$PreAbsJ
SexHourTable$M = SexHourTable$PreAbsM
SexHourTable$FJ = SexHourTable$F + SexHourTable$J
SexHourTable$JM = SexHourTable$M + SexHourTable$J
SexHourTable$FM = SexHourTable$F + SexHourTable$M
SexHourTable$FJM = SexHourTable$F + SexHourTable$J + SexHourTable$M

SexHourTable$FJM = replace(SexHourTable$FJM, which(SexHourTable$FJM <3 ), 0) #delete rows that don't have all sexes
SexHourTable$FJM = replace(SexHourTable$FJM, which(SexHourTable$FJM == 3), 1) #only keep rows that have both sexes and make it equal to 1

SexHourTable$FJ = replace(SexHourTable$FJ, which(SexHourTable$FJ <2), 0) #delete rows that don't have both sexes
SexHourTable$FJ = replace(SexHourTable$FJ, which(SexHourTable$FJM == 1), 0) #delete rows that have all sexes
SexHourTable$FJ = replace(SexHourTable$FJ, which(SexHourTable$FJ == 2), 1) #only keep rows that have both sexes and make it equal to 1

SexHourTable$JM = replace(SexHourTable$JM, which(SexHourTable$JM <2 ), 0) #delete rows that don't have both sexes
SexHourTable$JM = replace(SexHourTable$JM, which(SexHourTable$FJM == 1), 0) #delete rows that have all sexes
SexHourTable$JM = replace(SexHourTable$JM, which(SexHourTable$JM == 2), 1) #only keep rows that have both sexes and make it equal to 1

SexHourTable$FM = replace(SexHourTable$FM, which(SexHourTable$FM <2 ), 0) #delete rows that don't have both sexes
SexHourTable$FM = replace(SexHourTable$FM, which(SexHourTable$FJM == 1), 0) #delete rows that have all sexes
SexHourTable$FM = replace(SexHourTable$FM, which(SexHourTable$FM == 2), 1) #only keep rows that have both sexes and make it equal to 1

SexHourTable$F = replace(SexHourTable$F, which(SexHourTable$FJ ==1 | SexHourTable$FJM ==1 | SexHourTable$FM ==1),0) #delete F only rows when it's being accounted for in another group
SexHourTable$J = replace(SexHourTable$J, which(SexHourTable$FJ ==1 | SexHourTable$JM ==1 | SexHourTable$FJM ==1),0) #delete J only rows when it's being accounted for in another group
SexHourTable$M = replace(SexHourTable$M, which(SexHourTable$JM ==1 | SexHourTable$FM ==1 | SexHourTable$FJM ==1),0) #delete M only rows when it's being accounted for in another group

#loop through each site and pull out relevant data
for (i in 1:length(SiteNames)){
  site = SiteNames[i]
  SiteHourTable = dplyr::filter(SexHourTable, grepl(site,Site))
  siteHourTableGen = dplyr::filter(HourTable, grepl(site,Site))
  
  #Find Total Hours and Hours with Presence in General
  TotalHours = nrow(SiteHourTable)
  PresentHours = apply(SiteHourTable,2,function(x) sum(x > 0))
    if (i == 1){
    SumTable = data.frame("site" = site, "TotalHours" = TotalHours, "PresentHours" = sum(siteHourTableGen$PreAbs))
  }else{
    SumTableTemp = data.frame("site" = site, "TotalHours" = TotalHours, "PresentHours" = sum(siteHourTableGen$PreAbs))
    SumTable = rbind(SumTable, SumTableTemp)
  }
  
  #Find Ratio of F/M/J Days
  if (i == 1){
  SiteRatio = c(F = sum(SiteHourTable$F), J = sum(SiteHourTable$J), M = sum(SiteHourTable$M), "F&J" = sum(SiteHourTable$FJ),
                          "F&M" = sum(SiteHourTable$FM), "J&M" = sum(SiteHourTable$JM), "F&J&M" = sum(SiteHourTable$FJM))
  SuMSiteRatio = sum(SiteRatio)
  SiteRatioPercent = ceiling((SiteRatio/SuMSiteRatio) * 100)
  SiteRatioPercent = list(SiteRatioPercent)
  }else{
    SiteRatioTemp = c(F = sum(SiteHourTable$F), J = sum(SiteHourTable$J), M = sum(SiteHourTable$M), "F&J" = sum(SiteHourTable$FJ),
                               "F&M" = sum(SiteHourTable$FM), "J&M" = sum(SiteHourTable$JM), "F&J&M" = sum(SiteHourTable$FJM))
    SuMSiteRatioTEMP = sum(SiteRatioTemp)
    SiteRatioPercentTEMP = ceiling((SiteRatioTemp/SuMSiteRatioTEMP) * 100)
    SiteRatioPercent = append(SiteRatioPercent,list(SiteRatioPercentTEMP))
  }
}

#For some reason, I need to run this manually
# Pie Charts --------------------------------------------------------------
for (i in 1:length(SiteNames)){
  site = SiteNames[i]
  title = paste(saveDir,"/",site,"/PropVennHourly", site,".pdf",sep="")
  pdf(title)
  plot(euler(SiteRatioPercent[[i]]),fills = c('#66c2a5','#fc8d62','#8da0cb'), labels = FALSE, quantities = list(cex=3),
       main = paste("Proportional Hourly Venn Diagram for ",site,sep = ""))
  print(paste('Plot for ',site,' Saved',sep=""))
  dev.off()
}

# Relative Bars & Effort -----------------------------------------------------------

#Relative Effort
MaxEffort = max(SumTable$TotalHours)
SumTable$RelativeEffort = SumTable$TotalHours/MaxEffort

#Relative Presence
MaxPresence = max(SumTable$PresentHours)
SumTable$Relative_toEachOther_Presence = SumTable$PresentHours/MaxPresence
SumTable$Relative_toEffort_Presence = SumTable$PresentHours/SumTable$TotalHours

#Stacked bar plot for relative effort
title = paste(saveDir,"/RelativeEffortHours",".pdf",sep="")
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
title = paste(saveDir,"/Relative_toEachOther_PresenceHours",".pdf",sep="")
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
title = paste(saveDir,"/Relative_toEffort_PresenceHours",".pdf",sep="")
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

