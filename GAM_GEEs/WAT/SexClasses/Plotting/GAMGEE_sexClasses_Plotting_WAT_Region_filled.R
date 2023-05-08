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
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_sexClasses_Plotting_Functions_RealProbs_HistAbove_filled.R')  #on Nat's computer

# Load Workspace --------------------------------------------------
region = 'South'
site = "South"
GDrive = 'G'

saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_sexClasses.RData',sep="")
load(fileName)
GDrive = 'G'
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/",region,'/',sep="")
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
  if (region == 'South' & sex == 'Males'){
    ggPlot_JD_AfterYear_WAT(PODFinal,SiteHourTableB,sex,COL)
  }else if (region == 'North' & sex == 'Social Groups'){
    ggPlot_JD_AfterYearB(PODFinal,SiteHourTableB,sex,COL)
  }else if (region == 'North' & sex == 'Mid-Size' | region == 'North' & sex == 'Males'){
    ggPlot_JD_sex(PODFinal,SiteHourTableB,sex,COL)
  }else{
    ggPlot_JD_Last(PODFinal,SiteHourTableB,region,sex,COL)
  }
  
  # Plot Year ---------------------------------------------------------------
  if (region == 'South' & sex == 'Males'){
    print('Skip this') #skip this
  }else{    
  if (length(unique(SiteHourTableB$Year)) > 4){
    if (sex == 'Social Groups'){
    ggPlot_Year_WATT_north(PODFinal,SiteHourTableB,sex,COL)
    }else{
      ggPlot_Year_WATTT_north(PODFinal,SiteHourTableB,sex,COL)}
    }else{
      ggPlot_Year_WAT(PODFinal,SiteHourTableB,sex,COL)
    }
  }
  
  # Plot Site ---------------------------------------------------------------
  if (region == 'South'){
  ggPlot_Site_WAT(PODFinal,SiteHourTableB,region,sex,COL)
  }else if (region == 'North' & sex == 'Mid-Size' | region == 'North' & sex == 'Males'){
    ggPlot_Site_WATTT(PODFinal,SiteHourTableB,region,sex,COL)
  }else{
    ggPlot_Site_WATT(PODFinal,SiteHourTableB,region,sex,COL)
  }
}