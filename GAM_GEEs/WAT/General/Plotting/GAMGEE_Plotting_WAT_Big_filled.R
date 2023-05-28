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
library(plyr)           # for histograms above plots

#load functions
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_Plotting_Functions_RealProbs_HistAbove_filled.R')  #on Nat's computer

# Load Workspace --------------------------------------------------
region = 'AllSites'
GDrive = 'G'
COL = '#D3D3D3'
varOrder = cbind('Region','Julian Day','Year') #Variables in the final model and their order (ex: (1)'Julian Day','Year' (2)'Year,'Julian Day', (3)'Julian Day', (4)'Year')
saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/",sep="")
fileName = paste(saveWorkspace,'AllSites_gamgeeOutput.RData',sep="")
load(fileName)
GDrive = 'G'
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/AllSites",sep="")

SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting
SiteHourTableB$Year = as.factor(SiteHourTableB$Year) #Change year to categorical variable for plotting histograms

#All Models include Effort above
#Region
ggPlot_Region(PODFinal,SiteHourTableB,region,COL)
#Julian Day
ggPlot_JD_BigModel_WAT(PODFinal,SiteHourTableB,region,COL)
#Year
ggPlot_Year_WAT_BigM(PODFinal,SiteHourTableB,region,COL)

