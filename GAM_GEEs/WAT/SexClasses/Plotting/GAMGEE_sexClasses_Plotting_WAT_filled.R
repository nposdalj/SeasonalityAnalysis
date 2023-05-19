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
site = 'JAX'
GDrive = 'G'
saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses.RData',sep="")
load(fileName)
GDrive = 'G'
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/",site,'/',sep="")
sexGroups = c('Social Groups','Mid-Size','Males')
SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting
SiteHourTableB$Year = as.factor(SiteHourTableB$Year)

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
  if (site == 'BC' & sex == 'Social Groups' | site == 'WC' & sex == 'Social Groups'){
    ggPlot_JD_AfterYear_WAT(PODFinal,SiteHourTableB,sex,COL)
  }else if (site == 'BP' & sex == 'Social Groups' | site == 'BP' & sex == 'Mid-Size' | site == 'BP' & sex == 'Males' |
            site == 'JAX' & sex == 'Mid-Size' | site == 'JAX' & sex == 'Males'){
    print('Skip this')
  }else{
    ggPlot_JD_sex(PODFinal,SiteHourTableB,sex,COL)
  }

# Plot Year ---------------------------------------------------------------
  if (site == 'HZ' & sex == 'Social Groups' | site == 'HZ' & sex == 'Males' | site == 'OC' & sex == 'Social Groups' | 
      site == 'BC' & sex == 'Mid-Size' | site == 'WC' & sex == 'Males' | site == 'BP' & sex == 'Males' | site == 'BP' & sex == 'Mid-Size'){
    #skip this
    print('Skip this')
  }else if (site == 'BC' & sex == 'Social Groups' | site == 'WC' & sex == 'Social Groups' | site == 'BP' & sex == 'Social Groups' | 
            site == 'BP' & sex == 'Mid-Size' | site == 'JAX' & sex == 'Mid-Size' | site == 'JAX' & sex == 'Males'){
    ggPlot_Year_WAT_first(PODFinal,SiteHourTableB,sex,COL)
  }else{    
    if (length(unique(SiteHourTableB$Year)) > 4){
    ggPlot_Year_WATT(PODFinal,SiteHourTableB,sex,COL)
    } else if{ 
    }else{
    ggPlot_Year_WAT(PODFinal,SiteHourTableB,sex,COL)
    }
  }
}

