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
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_sexClasses_Plotting_Functions.R')

# Load Workspace --------------------------------------------------
site = 'Big'
#region = 'GOA'
sexGroups = c('Social Groups','Mid-Size','Males')

if (exists("site")){
  if (site == 'Big'){
    saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/")
    fileName = paste(saveWorkspace,'BigModel_gamgeeOutput_sexClasses.RData',sep="")
    load(fileName)
  }else{
saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses.RData',sep="")
load(fileName)
  }
}

if (exists("region")){
  saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
  fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_sexClasses.RData',sep="")
  load(fileName)
  site = region
}

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
if (site == 'CB' & sex == 'Social Groups'){
  BasePlot_JD_Year(PODFinal,SiteHourTableB,sex)
  ggPlot_JD_Year(PODFinal, SiteHourTableB,sex)
}else if (site == 'GOA' & !sex == 'Males'){
  BasePlot_JD_Year(PODFinal,SiteHourTableB,sex)
  ggPlot_JD_Year(PODFinal, SiteHourTableB,sex)
}else if (site == 'GOA' & sex == 'Males' | site == 'Big' & sex == 'Mid-Size'){
  BasePlot_JD_Year_M(PODFinalM,SiteHourTableB,sex)
  ggPlot_JD_Year_M(PODFinalM, SiteHourTableB,sex)
}else if (site == 'Big' & sex == 'Social Groups' | site == 'Big' & sex == 'Males'){
  BasePlot_JD_Year(PODFinal,SiteHourTableB,sex)
  ggPlot_JD_Year(PODFinal, SiteHourTableB,sex)
}else{
BasePlot_JD_sex(PODFinal,SiteHourTableB,sex)
ggPlot_JD_sex(PODFinal,SiteHourTableB,sex)
}

# Plot Year ---------------------------------------------------------------
if (site == 'CB'){
  if (sex == 'Social Groups'){
  ggPlot_Year(PODFinal,SiteHourTableB,sex)
  }else{
    ggPlot_Year_MM(PODFinal,SiteHourTableB,sex)
  }
}
  
if (region == 'GOA'){
    if (sex == 'Social Groups' | sex == 'Mid-Size'){
      ggPlot_Year(PODFinal,SiteHourTableB,sex)
}}
  
if (site == 'Big'){
  if (sex == 'Social Groups' | sex == 'Males'){
    ggPlot_Year_Big(PODFinal,SiteHourTableB,sex)
  }else{
    ggPlot_Year_Big_Mid(PODFinal,SiteHourTableB,sex)
  }
  } 

# Plot Site ---------------------------------------------------------------
if (exists("region")){
  if (region == 'GOA' & sex == 'Males'){
    ggPlot_Site_Year(PODFinal,SiteHourTableB,sex)
    ggPlot_Site_asFactor_Year(PODFinal,SiteHourTableB,sex)
  }else if (region == 'BSAI'){
    ggPlot_Site(PODFinal,SiteHourTableB,sex)
    ggPlot_Site_asFactor(PODFinal,SiteHourTableB,sex)
  }
}
}
