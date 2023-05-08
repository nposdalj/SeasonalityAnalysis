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
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_sexClasses_Plotting_Functions_RealProbs_HistAbove.R')  #on Nat's computer

# Load Workspace --------------------------------------------------
GDrive = 'G'
model = 'Big'
saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/",sep="")
fileName = paste(saveWorkspace,'_Big_gamgeeOutput_sexClasses.RData',sep="")
load(fileName)
GDrive = 'G'
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/AllSites",sep="")
sexGroups = c('Social Groups','Mid-Size','Males')

SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting
SiteHourTableB$Year = as.factor(SiteHourTableB$Year) #Change year to categorical variable for plotting histograms

for (i in 1:length(sexGroups)){
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
  if (sex == 'Males'){
    ggPlot_JD_AfterSite(PODFinal,SiteHourTableB,model,sex,COL)
  }else if (sex == 'Mid-Size'){
    ggPlot_JD_WATBIG(PODFinal,SiteHourTableB,model,sex,COL)
  }else{
    ggPlot_JD_SG_WATBIG(PODFinal,SiteHourTableB,model,sex,COL)
  }
  
  # Plot Year ---------------------------------------------------------------
  if (sex == 'Males'){
    ggPlot_Year_WAT_Big(PODFinal,SiteHourTableB,sex,COL)
  }else if (sex == 'Mid-Size'){
    ggPlot_Year_WATTT_Big(PODFinal,SiteHourTableB,sex,COL)
  }else{
    ggPlot_Year_WATT_Big(PODFinal,SiteHourTableB,sex,COL)
  }
  
  # Plot Region ---------------------------------------------------------------
  if (sex == 'Males'){
    ggPlot_Region_WAT_last(PODFinal,SiteHourTableB,model,sex,COL)
  }else{
  ggPlot_Region_WAT(PODFinal,SiteHourTableB,model,sex,COL)
  }
}