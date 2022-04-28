#Load libraries
library("eulerr")
library("tidyverse")

#load data
GDrive =  'I'
dir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites",sep="")

#General Data
fileName1 = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/AllSitesGrouped_GAMGEE_ROW.csv",sep="")#setting the directory
DayTable = read.csv(fileName1) #no effort days deleted
DayTable = na.omit(DayTable)
DayTable$tbin = as.Date(DayTable$tbin)

#Sex Specific Data
fileName2 = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/AllSitesGrouped_GAMGEE_ROW_sexClasses.csv",sep="")#setting the directory
SexDayTable = read.csv(fileName2) #no effort days deleted
SexDayTable = na.omit(SexDayTable)
SexDayTable$tbin = as.Date(SexDayTable$tbin)

#Make new columns to find ratio
SexDayTable$F = SexDayTable$PreAbsF
SexDayTable$J = SexDayTable$PreAbsJ
SexDayTable$M = SexDayTable$PreAbsM
SexDayTable$FJ = SexDayTable$F + SexDayTable$J
SexDayTable$FJ = replace(SexDayTable$FJ, which(SexDayTable$FJ <2 ), 0)
SexDayTable$FJ = replace(SexDayTable$FJ, which(SexDayTable$FJ == 2), 1)
SexDayTable$JM = SexDayTable$M + SexDayTable$J
SexDayTable$JM = replace(SexDayTable$JM, which(SexDayTable$JM <2 ), 0)
SexDayTable$JM = replace(SexDayTable$JM, which(SexDayTable$JM == 2), 1)
SexDayTable$FM = SexDayTable$F + SexDayTable$M
SexDayTable$FM = replace(SexDayTable$FM, which(SexDayTable$FM <2 ), 0)
SexDayTable$FM = replace(SexDayTable$FM, which(SexDayTable$FM == 2), 1)
SexDayTable$FJM = SexDayTable$F + SexDayTable$J + SexDayTable$M
SexDayTable$FJM = replace(SexDayTable$FJM, which(SexDayTable$FJM <3 ), 0)
SexDayTable$FJM = replace(SexDayTable$FJM, which(SexDayTable$FJM == 3), 1)

saveDir = paste(GDrive,":/My Drive/Manuscripts/GOA/Figures",sep="")

#Site Names
SiteNames = c('CB','PT','QN','BD','AB','KOA','KS')

#loop through each site and pull out relevant data
for (i in 1:length(SiteNames)){
  site = SiteNames[i]
  SiteDayTable = dplyr::filter(SexDayTable, grepl(site,Site))
  
  #Find Total Days and Days with Presence in General
  TotalDays = nrow(SiteDayTable)
  PresentDays = apply(SiteDayTable,2,function(x) sum(x > 0))
    if (i == 1){
    SumTable = data.frame("site" = site, "TotalDays" = TotalDays, "PresentDays" = PresentDays[10])
  }else{
    SumTableTemp = data.frame("site" = site, "TotalDays" = TotalDays, "PresentDays" = PresentDays[10])
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
  title = paste(saveDir,"/PropVenn", site,".pdf",sep="")
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
SumTable$RelativePresence = SumTable$PresentDays/SumTable$TotalDays

#Stacked bar plot for relative effort
title = paste(saveDir,"/RelativeEffort",".pdf",sep="")
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

#Stacked bar plot for relative presence
title = paste(saveDir,"/RelativePresence",".pdf",sep="")
pdf(title)
p = ggplot(SumTable, aes(x=site, y=RelativePresence)) +
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

