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
site = 'CB'
sex = 'Social Groups'
saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses.RData',sep="")
load(fileName)

# Set up Sex Classes --------------------------------------------------
if (sex == 'Social Groups'){
  PODFinal = PODFinalF
}

if (sex == 'Mid Size'){
  PODFinal = PODFinalJ
}

if (sex == 'Male'){
  PODFinal = PODFinalM
}

# Plot Julian Day ---------------------------------------------------------
if (site == 'CB'){
  BasePlot_JD_Year(PODFinal,SiteHourTableB)
  ggPlot_JD_Year(PODFinal, SiteHourTableB)
} else {
BasePlot_JD_sex(PODFinal,SiteHourTableB,sex)
ggPlot_JD_sex(PODFinal,SiteHourTableB,sex)
}

# Plot Year ---------------------------------------------------------------
if (site == 'CB'){
  ggPlot_Year(PODFinal,SiteHourTableB)
}
