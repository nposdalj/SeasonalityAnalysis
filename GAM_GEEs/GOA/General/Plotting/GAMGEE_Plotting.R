rm(list = ls()) #clear environment

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
library(scales)
library(magick)
library(cowplot)
#load functions
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_Plotting_Functions_RealProbs_HistAbove.R')

# Load Workspace --------------------------------------------------
site = 'Big'
GDrive = 'I'
#region = 'GOA'
if (exists("site")){
  if (site == 'Big'){
    saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/",sep="")
    fileName = paste(saveWorkspace,'BigModel_gamgeeOutput_modified.RData',sep="")
    load(fileName)
    saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/AllSites/",sep="")
  }else{
    saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
    fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_modified.RData',sep="")
    load(fileName)
    saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",site,'/',sep="")
  }
}

if (exists("region")){
  saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
  fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_modified.RData',sep="")
  load(fileName)
  saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",region,'/',sep="")
}

SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting

SiteHourTableB$Year = as.factor(SiteHourTableB$Year) #Change year to categorical variable for plotting histograms

# Plot Julian Day ---------------------------------------------------------
if (exists("site")){
  if (site == 'Big'){
    ggPlot_JD_Year(PODFinal, SiteHourTableB)
  }else if (site == 'CB'){
  }else{
    ggPlot_JD(PODFinal,SiteHourTableB,site)
  }
}

if (exists("region")){
  site = region
  if (region == 'GOA'){
    ggPlot_JD_AfterYear(PODFinal, SiteHourTableB,site)
  }else{
    ggPlot_JD(PODFinal,SiteHourTableB,site)
}
} 


# Plot Year ---------------------------------------------------------------
if (site == 'CB' | site =='GOA'){
  ggPlot_Year(PODFinal,SiteHourTableB,site)
}else if (site == "Big"){
  ggPlot_Year_Big(PODFinal,SiteHourTableB,site)
}

# Plot Site ---------------------------------------------------------------
if (exists("region")){
   if (region == 'BSAI'){
ggPlot_Site(PODFinal,SiteHourTableB,site)
  }
}

