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
library(scales)
library(magick)
library(cowplot)
#load functions
source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/GAMGEE_sexClasses_Plotting_Functions_RealProbs_HistAbove.R')

# Load Workspace --------------------------------------------------
site = 'BP'
GDrive = 'I'
saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses.RData',sep="")
load(fileName)
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/",site,'/',sep="")

sexGroups = c('Social Groups','Mid-Size','Males')

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
  if (site == 'BP'){
    print('Skip this')
  }else{
    ggPlot_JD_sex(PODFinal,SiteHourTableB,sex)
  }

# Plot Year ---------------------------------------------------------------
  if (site == 'BS'& sex == 'Mid-Size' | site == 'BP' & sex == 'Social Groups' | site == 'BP' & sex == 'Males'){
    #skip this
    print('Skip this')
  }else if (site == 'BP' & sex == 'Mid-Size'){
    ggPlot_Year_WAT_first(PODFinal,SiteHourTableB,sex)
  }else{
    ggPlot_Year_WAT(PODFinal,SiteHourTableB,sex)
  }
}


