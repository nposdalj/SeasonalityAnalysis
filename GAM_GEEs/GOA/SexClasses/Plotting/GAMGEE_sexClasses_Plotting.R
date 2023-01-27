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
library(ggExtra)
#load functions
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_sexClasses_Plotting_Functions_RealProbs_HistAbove.R')

# Load Workspace --------------------------------------------------
#site = 'Big'
region = 'GOA'
GDrive = 'I'
sexGroups = c('Social Groups','Mid-Size','Males')

if (exists("site")){
  if (site == 'Big'){
    saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/",sep="")
    fileName = paste(saveWorkspace,'BigModel_gamgeeOutput_sexClasses_modified.RData',sep="")
    load(fileName)
  }else{
saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses_modified.RData',sep="")
load(fileName)
  }
}

if (exists("region")){
  saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
  fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_sexClasses_modified.RData',sep="")
  load(fileName)
  site = region
}

#If it's a leap year, delete julian day 366 for plotting
SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),]

for (i in 1:length(sexGroups)){
  sex = sexGroups[i]

# Set up Sex Classes --------------------------------------------------
if (sex == 'Social Groups'){
    PODFinal = PODFinalF
    pr = prf
}else if(sex == 'Mid-Size'){
    PODFinal = PODFinalJ
    pr = prj
}else if(sex == 'Males'){
    PODFinal = PODFinalM
    pr = prm
}

# Plot Julian Day ---------------------------------------------------------
if (site == 'CB' & sex == 'Social Groups' | site == 'GOA' & sex == 'Mid-Size'){
  ggPlot_JD_AfterYear(PODFinal, SiteHourTableB,sex)
}else if (site == 'GOA' & sex == 'Males'){
  ggPlot_JD_AfterSite(PODFinal, SiteHourTableB,sex)
}else if (site == 'Big'){
  ggPlot_JD_AfterYearB(PODFinal, SiteHourTableB,sex)
}else if (site == 'GOA' & sex == 'Social Groups'){
  ggPlot_JD_AfterYearSite(PODFinal,SiteHourTableB,sex)
}else{
  ggPlot_JD_sex(PODFinal,SiteHourTableB,sex)
}

# Plot Year ---------------------------------------------------------------
if (site == 'CB'){
  if (sex=='Males'){
  ggPlot_Year_AfterJD(PODFinal,SiteHourTableB,sex)
  }else if (sex == 'Social Groups'){
    ggPlot_Year(PODFinal,SiteHourTableB,sex)
  }
}
  
if (exists("region")){
if (region == 'GOA' & sex == 'Mid-Size' | region == 'GOA' & sex == 'Social Groups'){
      ggPlot_Year(PODFinal,SiteHourTableB,sex)
}}
  
if (site == 'Big'){
    ggPlot_Year_Big(PODFinal,SiteHourTableB,sex)
} 

# Plot Site ---------------------------------------------------------------
if (exists("region")){
  if (region == 'GOA' & sex == 'Males'){
    ggPlot_Site_Year(PODFinal,SiteHourTableB,sex)
  }else if (region == 'BSAI' & sex == 'Social Groups' | region == 'BSAI' & sex == 'Males'){
    ggPlot_Site(PODFinal,SiteHourTableB,sex)
  }else if (region == 'GOA' & sex == 'Social Groups'){
    ggPlot_Site_AfterYear(PODFinal,SiteHourTableB,sex)
  }
}
}

