# This script will load all of the workspaces from the models, save the final model,
# and eventually plot all of the residuals for each model

#load libraries
library(sure)
library(PResiduals)
library(arm)
library(performance)
library(parameters)
library(see)

#load functions
#source('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/GAM_GEEs/Binned_Residuals_Functions.R')

#Specify Directories and find subdirectories
GDrive = 'H'
MainDir = paste(GDrive,':/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis',sep='')
PlotDir = paste(GDrive,':/My Drive/GofAK_TPWS_metadataReduced/Plots',sep='')
AllFolders = list.dirs(MainDir, full.names = TRUE, recursive = TRUE) #find all subfolders

#SITES - General
#Loop through each folder and load R workspace, extract final model for residual plotting
for (i in 2:length(AllFolders)){
  #Seperate site name
  SiteFolder = AllFolders[i]
  Site = gsub(".*SeasonalityAnalysis/","",SiteFolder)
  
  FileName = paste(SiteFolder,'/',Site,'_SiteSpecific_gamGeeOutput.RData',sep='')#FileName
  SaveDir = paste(PlotDir,'/',Site,sep='')
  
  #Load workspace (if it exists)
  if (file.exists(FileName)){
  load(FileName)
    jpeg(file = paste(SaveDir,'/Residuals_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin, SiteHourTableB$PreAbs)
    plot(PODFinal$residuals)
    plot(PODFinal$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinal),resid(PODFinal))
    qqnorm(resid(PODFinal))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(predict(PODFinal),resid(PODFinal,type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinal))
    dev.off()
    
    filename = paste(MainDir,'/',Site,'/',Site,'_ModelEval.txt',sep="")
    sink(filename)
    binned_residuals(PODFinal)
    print(r2(PODFinal))
    sink(file = NULL)
  }
} 

#SITES - Sex Specific
#Loop through each folder and load R workspace, extract final model for residual plotting
for (i in 2:length(AllFolders)){
  #Seperate site name
  SiteFolder = AllFolders[i]
  Site = gsub(".*SeasonalityAnalysis/","",SiteFolder)
  
  FileName = paste(SiteFolder,'/',Site,'_SiteSpecific_gamGeeOutput_sexClasses.RData',sep='')#FileName
  SaveDir = paste(PlotDir,'/',Site,sep='')
  
  #Load workspace (if it exists)
  if (file.exists(FileName)){
    load(FileName)
    
    #Social Groups
    jpeg(file = paste(SaveDir,'/Residuals_SocialGroups',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsF)
    plot(PODFinalF$residuals)
    plot(PODFinalF$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_SocialGroups_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalF),resid(PODFinalF))
    qqnorm(resid(PODFinalF))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_SocialGroups_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalF),resid(PODFinalF, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_SocialGroups_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalF))
    dev.off()
    
    #Mid-Size
    jpeg(file = paste(SaveDir,'/Residuals_MidSize',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsJ)
    plot(PODFinalJ$residuals)
    plot(PODFinalJ$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_MidSize_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalJ),resid(PODFinalJ))
    qqnorm(resid(PODFinalJ))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_MidSize_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalJ),resid(PODFinalJ, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_MidSize_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalJ))
    dev.off()
    
    #Males
    jpeg(file = paste(SaveDir,'/Residuals_Males',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsM)
    plot(PODFinalM$residuals)
    plot(PODFinalM$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_Males_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalM),resid(PODFinalM))
    qqnorm(resid(PODFinalM))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_Males_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalM),resid(PODFinalM, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_Males_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalM))
    dev.off()
    
    filename = paste(MainDir,'/',Site,'/',Site,'_ModelEvalSexGroups.txt',sep="")
    sink(filename)
    binned_residuals(PODFinalF)
    print(r2(PODFinalF))
    binned_residuals(PODFinalJ)
    print(r2(PODFinalJ))
    binned_residuals(PODFinalM)
    print(r2(PODFinalM))
    sink(file = NULL)
  }
}

#Regions - General
#Loop through each folder and load R workspace, extract final model for residual plotting
for (i in 2:length(AllFolders)){
  #Seperate site name
  SiteFolder = AllFolders[i]
  Site = gsub(".*SeasonalityAnalysis/","",SiteFolder)
  
  FileName = paste(SiteFolder,'/',Site,'_RegionSpecific_gamGeeOutput.RData',sep='')#FileName
  SaveDir = paste(PlotDir,'/',Site,sep='')
  
  #Load workspace (if it exists)
  if (file.exists(FileName)){
    load(FileName)
    jpeg(file = paste(SaveDir,'/Residuals_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbs)
    plot(PODFinal$residuals)
    plot(PODFinal$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinal),resid(PODFinal))
    qqnorm(resid(PODFinal))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinal),resid(PODFinal, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinal))
    dev.off()
    
    filename = paste(MainDir,'/',Site,'/',Site,'_ModelEval.txt',sep="")
    sink(filename)
    binned_residuals(PODFinal)
    print(r2(PODFinal))
    sink(file = NULL)
  }
} 

#Region - sex specific
#Loop through each folder and load R workspace, extract final model for residual plotting
for (i in 2:length(AllFolders)){
  #Seperate site name
  SiteFolder = AllFolders[i]
  Site = gsub(".*SeasonalityAnalysis/","",SiteFolder)
  
  FileName = paste(SiteFolder,'/',Site,'_RegionSpecific_gamGeeOutput_sexClasses.RData',sep='')#FileName
  SaveDir = paste(PlotDir,'/',Site,sep='')
  
  #Load workspace (if it exists)
  if (file.exists(FileName)){
    load(FileName)
    
    #Social Groups
    jpeg(file = paste(SaveDir,'/Residuals_SocialGroups',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsF)
    plot(PODFinalF$residuals)
    plot(PODFinalF$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_SocialGroups_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalF),resid(PODFinalF))
    qqnorm(resid(PODFinalF))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_SocialGroups',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalF),resid(PODFinalF, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_SocialGroups_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalF))
    dev.off()
    
    #Mid-Size
    jpeg(file = paste(SaveDir,'/Residuals_MidSize',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsJ)
    plot(PODFinalJ$residuals)
    plot(PODFinalJ$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_MidSize_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalJ),resid(PODFinalJ))
    qqnorm(resid(PODFinalJ))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_MidSize_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalJ),resid(PODFinalJ, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_MidSize_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalJ))
    dev.off()
    
    #Males
    jpeg(file = paste(SaveDir,'/Residuals_Males',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsM)
    plot(PODFinalM$residuals)
    plot(PODFinalM$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_Males_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalM),resid(PODFinalM))
    qqnorm(resid(PODFinalM))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_Males',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalM),resid(PODFinalM,type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_Males_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalM))
    dev.off()
    
    filename = paste(MainDir,'/',Site,'/',Site,'_ModelEvalSexGroups.txt',sep="")
    sink(filename)
    binned_residuals(PODFinalF)
    print(r2(PODFinalF))
    binned_residuals(PODFinalJ)
    print(r2(PODFinalJ))
    binned_residuals(PODFinalM)
    print(r2(PODFinalM))
    sink(file = NULL)
  }
}

#Big Model
#Loop through each folder and load R workspace, extract final model for residual plotting
for (i in 2:length(AllFolders)){
  #Seperate site name
  SiteFolder = AllFolders[i]
  Site = gsub(".*SeasonalityAnalysis/","",SiteFolder)
  
  FileName = paste(SiteFolder,'/',Site,'_gamGeeOutput.RData',sep='')#FileName
  SaveDir = paste(PlotDir,'/',Site,sep='')
  
  #Load workspace (if it exists)
  if (file.exists(FileName)){
    load(FileName)
    jpeg(file = paste(SaveDir,'/Residuals_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbs)
    plot(PODFinal$residuals)
    plot(PODFinal$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinal),resid(PODFinal))
    qqnorm(resid(PODFinal))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinal),resid(PODFinal, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinal))
    dev.off()
    
    filename = paste(MainDir,'/',Site,'/',Site,'_ModelEval.txt',sep="")
    sink(filename)
    binned_residuals(PODFinal)
    print(r2(PODFinal))
    sink(file = NULL)
  }
} 

#Big Model - Sex Specific
#Loop through each folder and load R workspace, extract final model for residual plotting
for (i in 2:length(AllFolders)){
  #Seperate site name
  SiteFolder = AllFolders[i]
  Site = gsub(".*SeasonalityAnalysis/","",SiteFolder)
  
  FileName = paste(SiteFolder,'/',Site,'_gamGeeOutput_sexClasses.RData',sep='')#FileName
  SaveDir = paste(PlotDir,'/',Site,sep='')
  
  #Load workspace (if it exists)
  if (file.exists(FileName)){
    load(FileName)
    
    #Social Groups
    jpeg(file = paste(SaveDir,'/Residuals_SocialGroups',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsF)
    plot(PODFinalF$residuals)
    plot(PODFinalF$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_SocialGroups_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalF),resid(PODFinalF))
    qqnorm(resid(PODFinalF))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_SocialGroups',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalF),resid(PODFinalF, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_SocialGroups_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalF))
    dev.off()
    
    #Mid-Size
    jpeg(file = paste(SaveDir,'/Residuals_MidSize',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsJ)
    plot(PODFinalJ$residuals)
    plot(PODFinalJ$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_MidSize_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalJ),resid(PODFinalJ))
    qqnorm(resid(PODFinalJ))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_MidSize_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalJ),resid(PODFinalJ, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_MidSize_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalJ))
    dev.off()
    
    #Males
    jpeg(file = paste(SaveDir,'/Residuals_Males',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,3))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsM)
    plot(PODFinalM$residuals)
    plot(PODFinalM$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_Males_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalM),resid(PODFinalM))
    qqnorm(resid(PODFinalM))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedPlots_Males_',Site,'.jpeg',sep=''),width = 500, height = 400)
    binnedplot(fitted(PODFinalM),resid(PODFinalM, type="response"))
    dev.off()
    
    jpeg(file = paste(SaveDir,'/BinnedResiduals_Males_',Site,'.jpeg',sep=''),width = 500, height = 400)
    print(binned_residuals(PODFinalM))
    dev.off()
    
    filename = paste(MainDir,'/',Site,'/',Site,'_ModelEvalSexGroups.txt',sep="")
    sink(filename)
    binned_residuals(PODFinalF)
    print(r2(PODFinalF))
    binned_residuals(PODFinalJ)
    print(r2(PODFinalJ))
    binned_residuals(PODFinalM)
    print(r2(PODFinalM))
    sink(file = NULL)
  }
}

