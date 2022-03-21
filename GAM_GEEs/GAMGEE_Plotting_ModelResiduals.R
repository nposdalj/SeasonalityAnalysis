# This script will load all of the workspaces from the models, save the final model,
# and eventually plot all of the residuals for each model

#load libraries

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
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin, SiteHourTableB$PreAbs)
    plot(PODFinal$residuals)
    plot(PODFinal$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinal),resid(PODFinal))
    qqnorm(resid(PODFinal))
    dev.off()
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
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsF)
    plot(PODFinalF$residuals)
    plot(PODFinalF$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_SocialGroups_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalF),resid(PODFinalF))
    qqnorm(resid(PODFinalF))
    dev.off()
    
    #Mid-Size
    jpeg(file = paste(SaveDir,'/Residuals_MidSize',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsJ)
    plot(PODFinalJ$residuals)
    plot(PODFinalJ$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_MidSize_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalJ),resid(PODFinalJ))
    qqnorm(resid(PODFinalJ))
    dev.off()
    
    #Males
    jpeg(file = paste(SaveDir,'/Residuals_Males',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsM)
    plot(PODFinalM$residuals)
    plot(PODFinalM$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_Males_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalM),resid(PODFinalM))
    qqnorm(resid(PODFinalM))
    dev.off()
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
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbs)
    plot(PODFinal$residuals)
    plot(PODFinal$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinal),resid(PODFinal))
    qqnorm(resid(PODFinal))
    dev.off()
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
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsF)
    plot(PODFinalF$residuals)
    plot(PODFinalF$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_SocialGroups_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalF),resid(PODFinalF))
    qqnorm(resid(PODFinalF))
    dev.off()
    
    #Mid-Size
    jpeg(file = paste(SaveDir,'/Residuals_MidSize',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsJ)
    plot(PODFinalJ$residuals)
    plot(PODFinalJ$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_MidSize_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalJ),resid(PODFinalJ))
    qqnorm(resid(PODFinalJ))
    dev.off()
    
    #Males
    jpeg(file = paste(SaveDir,'/Residuals_Males',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsM)
    plot(PODFinalM$residuals)
    plot(PODFinalM$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_Males_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalM),resid(PODFinalM))
    qqnorm(resid(PODFinalM))
    dev.off()
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
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbs)
    plot(PODFinal$residuals)
    plot(PODFinal$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinal),resid(PODFinal))
    qqnorm(resid(PODFinal))
    dev.off()
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
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsF)
    plot(PODFinalF$residuals)
    plot(PODFinalF$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_SocialGroups_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalF),resid(PODFinalF))
    qqnorm(resid(PODFinalF))
    dev.off()
    
    #Mid-Size
    jpeg(file = paste(SaveDir,'/Residuals_MidSize',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsJ)
    plot(PODFinalJ$residuals)
    plot(PODFinalJ$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_MidSize_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalJ),resid(PODFinalJ))
    qqnorm(resid(PODFinalJ))
    dev.off()
    
    #Males
    jpeg(file = paste(SaveDir,'/Residuals_Males',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(3,1))
    plot(SiteHourTableB$tbin,SiteHourTableB$PreAbsM)
    plot(PODFinalM$residuals)
    plot(PODFinalM$fitted.values)
    dev.off()
    
    jpeg(file = paste(SaveDir,'/ResQQ_Males_',Site,'.jpeg',sep=''),width = 1000, height = 400)
    par(mfrow=c(1,2))
    plot(fitted(PODFinalM),resid(PODFinalM))
    qqnorm(resid(PODFinalM))
    dev.off()
  }
}