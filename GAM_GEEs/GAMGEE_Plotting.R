#load libraries
library(boot)
library(pracma)
library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
library(splines)         # to construct the B-splines within a GEE-GLM
library(tidyverse)       # because it literally does everything
library(rjags)           # replacement for geeglm which is out of date
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2)         # to build the partial residual plots
library(mvtnorm)         # to build the partial residual plots
library(gridExtra)       # to build the partial residual plots
library(SimDesign)
library(lubridate)
library(regclass)
library(mgcv)
library(ChemoSpecUtils)
library(car)            # to run an ANOVA
library(splines2)       # to use mSpline for the GEEs
#load functions
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_Plotting_Functions.R')

# Load Workspace --------------------------------------------------
#site = 'Big'
region = 'GOA'
if (exists("site")){
  if (site == 'Big'){
    saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/")
    fileName = paste(saveWorkspace,'BigModel_gamgeeOutput.RData',sep="")
    load(fileName)
  }else{
    saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
    fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
    load(fileName)
  }
}

if (exists("region")){
  saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
  fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput.RData',sep="")
  load(fileName)
  site = region
}

# Plot Julian Day ---------------------------------------------------------
if (exists("site")){
  if (site == 'CB'){
    BasePlot_JD_Year(PODFinal,SiteHourTableB)
    ggPlot_JD_Year(PODFinal, SiteHourTableB)
  }else if (site == 'Big'){
    BasePlot_JD_Year(PODFinal,SiteHourTableB)
    ggPlot_JD_Year(PODFinal, SiteHourTableB)
  }else{
    BasePlot_JD(PODFinal,SiteHourTableB)
    ggPlot_JD(PODFinal,SiteHourTableB)
  }
}

if (exists("region")){
  if (region == 'GOA'){
    BasePlot_JD_Year(PODFinal,SiteHourTableB)
    ggPlot_JD_Year(PODFinal, SiteHourTableB)
  }else{
    BasePlot_JD(PODFinal,SiteHourTableB)
    ggPlot_JD(PODFinal,SiteHourTableB)
}
} 


# Plot Year ---------------------------------------------------------------
if (site == 'CB'){
  ggPlot_Year(PODFinal,SiteHourTableB,site)
}else if (site == "Big"){
  ggPlot_Year_Big(PODFinal,SiteHourTableB,site)
  # ggPlot_Year_Big(PODFinal_Region,SiteHourTableB,site) #model looked bad with region
}else if (site == 'GOA'){
  ggPlot_Year(PODFinal,SiteHourTableB,site)
  # ggPlot_Year(PODFinal_Site,SiteHourTableB,site) #model looked bad with region
}

# Plot Site ---------------------------------------------------------------
if (exists("region")){
   if (region == 'GOA'){
     #do nothing since it wasn't included in the model
# ggPlot_Site_Year(PODFinal_Site,SiteHourTableB,site)
# ggPlot_Site_asFactor_Year(PODFinal_Site,SiteHourTableB,site)
   }else{
ggPlot_Site(PODFinal,SiteHourTableB,site)
ggPlot_Site_asFactor(PODFinal,SiteHourTableB,site)
  }
}

# Plot Region ---------------------------------------------------------------
if (site == 'Big'){
  ggPlot_Region_Big(PODFinal_Region,SiteHourTableB,site)
  ggPlot_Region_asFactor_Big(PODFinal_Region,SiteHourTableB,site)
}

