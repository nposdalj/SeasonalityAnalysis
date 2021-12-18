#Load libraries
library("eulerr")
library("tidyverse")

#Load data from PropVenn Matlab Output (Manually put in the proportions into this .csv)
Prop = read.csv("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/Proportion.csv")
saveDir = 'I:/My Drive/Manuscripts/GOA/Figures'

#Plot Proportions
PropPlot = Prop #replicate table
PropPlot$Site <- NULL

#Crate data frames
CB <-  c(F = 1, J = 23, M = 14,"F&J" = 27, "J&M" = 313, "F&M" = 26, "F&J&M" = 26)
PT <- c(F = 19, J = 54, M = 21,"F&J" = 22, "J&M" = 46, "F&M" = 11, "F&J&M" = 9)
QN <- c(F = 16, J = 20, M = 38,"F&J" = 22, "J&M" = 65, "F&M" = 12, "F&J&M" = 9)
BD <- c(F = 5, J = 46, M = 130,"F&J" = 13, "J&M" = 134, "F&M" = 16, "F&J&M" = 7)
AB <- c(F = 2, J = 24, M = 24,"F&J" = 1, "J&M" = 14, "F&M" = 1, "F&J&M" = 1)
KOA <- c(F = 0, J = 32, M = 50,"F&J" = 0, "J&M" = 22, "F&M" = 0, "F&J&M" = 0)
KS <- c(F = 0, J = 22, M = 19,"F&J" = 0, "J&M" = 16, "F&M" = 0, "F&J&M" = 0)
AllSites = list(CB,PT,QN,BD,AB,KOA,KS)

#Site Names
SiteNames = c('CB','PT','QN','BD','AB','KOA','KS')

for (i in 1:length(SiteNames)){
  site = SiteNames[i]
  title = paste(saveDir,"/PropVenn", site,".pdf",sep="")
  pdf(title)
  plot(euler(AllSites[[i]]),quantities = list(cex=1.5,type = "percent"),labels = list(labels= c("Social Groups", "Mid-Size", "Adult Males"),cex=1.5),
       main = paste("Proportional Venn Diagram for ",site,sep = ""))
  print(paste('Plot for ',site,' Saved',sep=""))
  dev.off()
}


