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
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_Plotting_Functions_RealProbs_HistAbove.R')  #on Nat's computer

# Load Workspace --------------------------------------------------
region = 'North'
GDrive = 'I'
varOrder = cbind('Site','Julian Day','Year') #Variables in the final model and their order (ex: (1)'Julian Day','Year' (2)'Year,'Julian Day', (3)'Julian Day', (4)'Year')
saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_modified.RData',sep="")
load(fileName)
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/",region,'/',sep="")

SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting

SiteHourTableB$Year = as.factor(SiteHourTableB$Year) #Change year to categorical variable for plotting histograms

if (region == 'North'){
#Site
ggPlot_Site_WATN(PODFinal,SiteHourTableB,region)
#JulianDay
ggPlot_JD_AfterYear_WAT_2015(PODFinal,SiteHourTableB,region)
#Year
ggPlot_Year_Regional(PODFinal,SiteHourTableB,region)}

if (region == 'South'){
#Site
  ggPlot_Site_WATS(PODFinal,SiteHourTableB,region)
#JulianDay
  ggPlot_JD_AfterYear_WAT(PODFinal,SiteHourTableB,region)
#Year
ggPlot_Year_WATS(PODFinal,SiteHourTableB,region)}
