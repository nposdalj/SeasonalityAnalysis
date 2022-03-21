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
site = 'Big'
GDrive = 'H'
#region = 'GOA'
if (exists("site")){
  if (site == 'Big'){
    saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/",sep="")
    fileName = paste(saveWorkspace,'BigModel_gamgeeOutput.RData',sep="")
    load(fileName)
  }else{
    saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
    fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
    load(fileName)
  }
}

if (exists("region")){
  saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
  fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput.RData',sep="")
  load(fileName)
}

#If it's a leap year, delete julian day 366 for plotting
SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),]

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
    ggPlot_JD(PODFinal,SiteHourTableB,site)
  }
}

if (exists("region")){
  site = region
  if (region == 'GOA'){
    BasePlot_JD_AfterSite(PODFinal,SiteHourTableB)
    ggPlot_JD_AfterSite(PODFinal, SiteHourTableB)
  }else{
    BasePlot_JD(PODFinal,SiteHourTableB)
    ggPlot_JD(PODFinal,SiteHourTableB,site)
}
} 


# Plot Year ---------------------------------------------------------------
if (site == 'CB'){
  ggPlot_Year(PODFinal,SiteHourTableB,site)
}else if (site == "Big"){
  ggPlot_Year_Big(PODFinal,SiteHourTableB,site)
}

# Plot Site ---------------------------------------------------------------
if (exists("region")){
   if (region == 'GOA'){
ggPlot_Site_GOA(PODFinal,SiteHourTableB,site)
ggPlot_Site_asfactor_GOA(PODFinal,SiteHourTableB,site)
   }else{
ggPlot_Site(PODFinal,SiteHourTableB,site)
ggPlot_Site_asFactor(PODFinal,SiteHourTableB,site)
  }
}

# Plot Region ---------------------------------------------------------------
if (site == 'Big'){
  ggPlot_Region_Big(PODFinal,SiteHourTableB,site)
  ggPlot_Region_asFactor_Big(PODFinal,SiteHourTableB,site)
}

