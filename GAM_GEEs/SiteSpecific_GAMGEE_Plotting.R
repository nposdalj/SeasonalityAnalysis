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
source('C:/Users/Alba/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_Plotting_Functions.R')

# Load Workspace --------------------------------------------------
site = 'PT'
saveWorkspace = paste("D:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
load(fileName)
saveDir = paste("D:/My Drive/GofAK_TPWS_metadataReduced/Plots/",site,'/',sep="")

# Plot Julian Day ---------------------------------------------------------
if (site == 'CB'){
  BasePlot_JD_Year(PODFinal,SiteHourTableB)
  ggPlot_JD_Year(PODFinal, SiteHourTableB)
} else {
BasePlot_JD(PODFinal,SiteHourTableB)
ggPlot_JD(PODFinal,SiteHourTableB,site)
}

# Plot Year ---------------------------------------------------------------
if (site == 'CB'){
  ggPlot_Year(PODFinal,SiteHourTableB)
}
