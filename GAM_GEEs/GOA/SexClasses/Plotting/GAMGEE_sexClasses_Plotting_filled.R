rm(list = ls()) #clear environment

#load libraries
library(boot)
library(pracma)
library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
library(splines)         # to construct the B-splines within a GEE-GLM
library(tidyverse)       # because it literally does everything
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2)         # to build the partial residual plots
library(mvtnorm)         # to build the partial residual plots
library(gridExtra)       # to build the partial residual plots
library(lubridate)
library(regclass)
library(mgcv)
library(ChemoSpecUtils)
library(car)            # to run an ANOVA
library(splines2)       # to use mSpline for the GEEs
library(scales)
library(magick)
library(cowplot)
library(ggExtra)
library(plyr)

#load functions
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_sexClasses_Plotting_Functions_RealProbs_HistAbove_filled.R')

# Load Workspace --------------------------------------------------
site = 'QN'
#region = 'BSAI'
GDrive = 'G'
sexGroups = c('Social Groups','Mid-Size','Males')

if (exists("site")){
  if (site == 'Big'){
    saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/",sep="")
    fileName = paste(saveWorkspace,'BigModel_gamgeeOutput_sexClasses_modified.RData',sep="")
    load(fileName)
    GDrive = 'G'
    saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/AllSites/",sep="")
  }else{
saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses_modified.RData',sep="")
load(fileName)
GDrive = 'G'
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",site,'/',sep="")
  }
}

if (exists("region")){
  saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
  fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_sexClasses_modified.RData',sep="")
  load(fileName)
  site = region
  GDrive = 'G'
  saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",region,'/',sep="")
}

SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting

SiteHourTableB$Year = as.factor(SiteHourTableB$Year) #Change year to categorical variable for plotting histograms

for (i in 1:length(sexGroups)){
  GDrive = 'G'
  
  # if (exists("region")){
  # saveDir = saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
  # }else{
  # saveDir = saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
  # }
  sex = sexGroups[i]

# Set up Sex Classes --------------------------------------------------
if (sex == 'Social Groups'){
    PODFinal = PODFinalF
    pr = prf
    COL = '#66c2a5' 
}else if(sex == 'Mid-Size'){
    PODFinal = PODFinalJ
    pr = prj
    COL = '#fc8d62'
}else if(sex == 'Males'){
    PODFinal = PODFinalM
    pr = prm
    COL = '#8da0cb'
}

# Plot Julian Day ---------------------------------------------------------
if (site == 'CB' & sex == 'Social Groups' | site == 'GOA' & sex == 'Mid-Size'){
  ggPlot_JD_AfterYear(PODFinal, SiteHourTableB,site,sex,COL)
}else if (site == 'GOA' & sex == 'Males'){
  ggPlot_JD_AfterSite(PODFinal, SiteHourTableB,site,sex,COL)
}else if (site == 'Big'){
  ggPlot_JD_AfterYearB(PODFinal, SiteHourTableB,site,sex,COL)
}else if (site == 'GOA' & sex == 'Social Groups'){
  ggPlot_JD_AfterYearSite(PODFinal,SiteHourTableB,site,sex,COL)
}else{
  ggPlot_JD_sex(PODFinal,SiteHourTableB,sex,COL)
}

# Plot Year ---------------------------------------------------------------
if (site == 'CB'){
  if (sex=='Males'){
  ggPlot_Year_AfterJD(PODFinal,SiteHourTableB,sex,COL)
  }else if (sex == 'Social Groups'){
    ggPlot_Year(PODFinal,SiteHourTableB,sex,COL)
  }
}
  
if (exists("region")){
if (region == 'GOA' & sex == 'Mid-Size' | region == 'GOA' & sex == 'Social Groups'){
      ggPlot_Year(PODFinal,SiteHourTableB,sex,COL)
}}
  
if (site == 'Big'){
    ggPlot_Year_Big(PODFinal,SiteHourTableB,sex,COL)
} 

# Plot Site ---------------------------------------------------------------
if (exists("region")){
  if (region == 'GOA' & sex == 'Males'){
    ggPlot_Site_Year(PODFinal,SiteHourTableB,sex,COL)
  }else if (region == 'BSAI' & sex == 'Social Groups' | region == 'BSAI' & sex == 'Males'){
    ggPlot_Site(PODFinal,SiteHourTableB,sex,COL)
  }else if (region == 'GOA' & sex == 'Social Groups'){
    ggPlot_Site_AfterYear(PODFinal,SiteHourTableB,sex,COL)
  }
}
}

