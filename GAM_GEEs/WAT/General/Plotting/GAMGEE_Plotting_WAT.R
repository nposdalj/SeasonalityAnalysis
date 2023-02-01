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
site = 'OC'
GDrive = 'I'
varOrder = cbind('Julian Day','Year') #Variables in the final model and their order (ex: (1)'Julian Day','Year' (2)'Year,'Julian Day', (3)'Julian Day', (4)'Year')
saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
load(fileName)
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/",site,'/',sep="")

SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting

SiteHourTableB$Year = as.factor(SiteHourTableB$Year) #Change year to categorical variable for plotting histograms

#Only includes 2016-2019
if (site == 'BC' || site == 'GS' || site == 'BP' || site == 'BS' || site == 'WC' || site == 'JAX'){
if (length(varOrder) == 2){ #If the model has both variables
  if (varOrder[1] == "Julian Day"){ #If JD is the first variable in the model
  # Plot Julian Day First---------------------------------------------------------
  ggPlot_JD(PODFinal,SiteHourTableB,site)
  # Plot Year Second---------------------------------------------------------------
  ggPlot_Year_WAT(PODFinal,SiteHourTableB,site)}
  if (varOrder[1] == "Year"){ #If Year is the first variable in the model
  # Plot Year First---------------------------------------------------------
    ggPlot_Year_First_WAT(PODFinal,SiteHourTableB,site)
  # Plot Julian Day Second---------------------------------------------------------
    ggPlot_JD_AfterYear_WAT(PODFinal,SiteHourTableB,site)}}
if (length(varOrder) == 1){ #If the model only has one variable
  if (varOrder[1] == "Julian Day"){ #If the model only has JD
    ggPlot_JD(PODFinal,SiteHourTableB,site)}
  if (varOrder[1] == "Year"){ #If the model only has Year
    ggPlot_Year_First_WAT(PPODFinal,SiteHourTableB,site)}}}

#Includes 2015 - 2019
if (site == 'OC' || site == 'NC' || site == 'HZ'){
  if (length(varOrder) == 2){ #If the model has both variables
    if (varOrder[1] == "Julian Day"){ #If JD is the first variable in the model
      # Plot Julian Day First---------------------------------------------------------
      ggPlot_JD(PODFinal,SiteHourTableB,site)
      # Plot Year Second---------------------------------------------------------------
      ggPlot_Year_WAT_2015(PODFinal,SiteHourTableB,site)}
    if (varOrder[1] == "Year"){ #If Year is the first variable in the model
      # Plot Year First---------------------------------------------------------
      ggPlot_Year_First_WAT_2015(PODFinal,SiteHourTableB,site)
      # Plot Julian Day Second---------------------------------------------------------
      ggPlot_JD_AfterYear_WAT_2015(PODFinal,SiteHourTableB,site)}}
  if (length(varOrder) == 1){ #If the model only has one variable
    if (varOrder[1] == "Julian Day"){ #If the model only has JD
      ggPlot_JD(PODFinal,SiteHourTableB,site)}
    if (varOrder[1] == "Year"){ #If the model only has Year
      ggPlot_Year_First_WAT_2015(PODFinal,SiteHourTableB,site)}}}
