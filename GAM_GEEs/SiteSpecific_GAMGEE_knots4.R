### Modelling habitat preference of sperm whales using Generalized Additive Models with Generalized Estimating Equations ###
### Script adapted from Pirotta et al. (2011) and Benjamins ###  
### Example from the GofAK + BSAI ###
### 7 Models total:
        #Site specific models: CB, PT, QN, BD (more than 270 d of recording)
        #Region specific models: BSAI + GOA
        #Big model: all 7 sites

# All the libraries have to be installed prior to their utilization (see R help on library installation)

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
library(ggfortify)      # extract confidence interval for ACF plots

site = 'PT' #specify the site of interest

# Step 1: Load the Data -----------------------------------------------------------
dir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/Plots/",site, sep="")
saveWorkspace = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(site,Site))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable = subset(SiteHourTable, Site == site)#subset the table for the site only

#If it's a leap year, delete julian day 366
SiteHourTable = SiteHourTable[!(SiteHourTable$Julian==366),]

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region (GOA or BSAI).
# The column Julian represents the day of the year and the column year represents the year of recording.


# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
if (site == 'CB'){
BlockMod<-glm(PreAbs~
               bs(Julian,k=4)+
               TimeLost+
               as.factor(Year)
             ,data=SiteHourTable,family=binomial)

}else{
BlockMod<-glm(PreAbs~
                bs(Julian,k=4)+
                TimeLost, data=SiteHourTable,family=binomial)
}

ACF = acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1))
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval = ACFidx[1]

#create the blocks based on the full timesereies
startDate = SiteHourTable$tbin[1]
endDate = SiteHourTable$tbin[nrow(SiteHourTable)]
timeseries = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseries)/ACFval)), times=1, each=ACFval)
divdiff = nrow(timeseries) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseries$block = c(preBlock,lastVec)
timeseries = rename(timeseries, tbin = date)
SiteHourTableB = left_join(SiteHourTable,timeseries,by="tbin")

#Make blocks continuous 
gaps = check4Gaps(SiteHourTableB$block,tol=1)
UnBlock = as.data.frame(unique(SiteHourTableB$block))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$Blocks = UnBlock$sequence[match(SiteHourTableB$block,UnBlock$`unique(SiteHourTableB$block)`)]

difference = diff(SiteHourTableB$Blocks) #find difference between rows
gapsCont = check4Gaps(SiteHourTableB$Blocks)
SiteHourTableB$Waves = rep(0,nrow(SiteHourTableB)) #make space for waves

for (i in 1:nrow(gapsCont)){
  StartB = gapsCont$beg.indx[i]
  EndB = gapsCont$end.indx[i]
  Num = length(StartB:EndB)
  SiteHourTableB$Waves[StartB:EndB] = seq.int(1,Num)
}

# Step 3: ANOVA to Check for Significance of Variables --------------------
#ANOVA with car package
Anova(BlockMod)
summary(BlockMod)

#PT 
# bs(Julian, k = 4)  195.177  4     <2e-16 ***
#   TimeLost             1.282  1     0.2576   

#QN
# bs(Julian, k = 4)      654  4     <2e-16 ***
#   TimeLost                 8  1     0.0052 ** 

#BD
# bs(Julian, k = 4)      695  4     <2e-16 ***
#   TimeLost                 2  1       0.22  

#CB
# bs(Julian, k = 4)     1682  4     <2e-16 ***
#   TimeLost                 2  1       0.19    
# as.factor(Year)       1121  7     <2e-16 ***


# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
if (site == 'CB'){
  GLM1 = glm(PreAbs ~ bs(Julian) + TimeLost + as.factor(Year), family = binomial, data = SiteHourTableB)
  VIF(GLM1)
#CB
  # bs(Julian)      1.62  3            1.08
  # TimeLost        1.01  1            1.01
  # as.factor(Year) 1.62  7            1.04
} else {
  #Other sites
  GLM1 = glm(PreAbs~bs(Julian)+TimeLost,family=binomial,data=SiteHourTableB)
  #VIF scores in GLM to work out collinearity:
  VIF(GLM1)
}

#PT
# bs(Julian) 1.000182  3        1.000030
# TimeLost   1.000182  1        1.000091
#QN
# bs(Julian)    1  3               1
# TimeLost      1  1               1
#BD
# bs(Julian)    1  3               1
# TimeLost      1  1               1


# Step 5: Model Selection - Covariate Preparation -------------------------

# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,365,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.
POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

#Julian Day
POD0a = geeglm(PreAbs ~ mSpline(Julian,
                                knots=quantile(Julian, probs=c(0.333,0.666)),
                                Boundary.knots=c(1,365),
                                periodic=T), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD0b = geeglm(PreAbs ~ AvgDayMat, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model0A<-c("POD0", "POD0a", "POD0b")
QIC0A<-c(QIC(POD0)[1],QIC(POD0a)[1],QIC(POD0b)[1])
QICmod0A<-data.frame(rbind(model0A,QIC0A))
QICmod0A

#CB
#QIC          QIC.1            QIC.2
# model0A             POD0           POD0a            POD0b
# QIC0A   62743.5048448833 61894.0018326248 61896.3522746218
#Julian day as mspline, but I'll use a variance covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
# QIC0A   11571.2264645272 11410.2528505801 11410.1658116943
#Julian day as a variance covariance matrix

#QN
# QIC            QIC.1            QIC.2
# model0A             POD0            POD0a            POD0b
# QIC0A   16378.2460093917 15734.6417060365 15734.1555437428
#Julian day as a variance covariance matrix

#BD
# QIC            QIC.1            QIC.2
# model0A             POD0            POD0a            POD0b
# QIC0A   23650.3441943795 23015.7537566308 23011.7382084554
#Julian day as mspline has lower QIC but I'll still go with a variance covariance matrix

if (site == 'CB'){
  #Year
  POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  POD1c = geeglm(PreAbs ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),
                                  family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model1A<-c("POD0", "POD1a", "POD1b","POD1c")
  QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
  QICmod1A<-data.frame(rbind(model1A,QIC1A))
  QICmod1A
}

#CB
# QIC            QIC.1            QIC.2           QIC.3
# model1A             POD0            POD1a            POD1b           POD1c
# QIC1A   62743.5048448833 61291.184631157 62727.1764111448 61263.0577232449
#Year as mspline but I'll use factor.

#TimeLost
POD2a = geeglm(PreAbs ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2b = geeglm(PreAbs ~ TimeLost, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2c = geeglm(PreAbs ~ mSpline(TimeLost,
                                knots=4),
                                family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model2A<-c("POD0", "POD2a", "POD2b","POD2c")
QIC2A<-c(QIC(POD0)[1],QIC(POD2a)[1],QIC(POD2b)[1],QIC(POD2c)[1])
QICmod2A<-data.frame(rbind(model2A,QIC2A))
QICmod2A

#CB
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   62743.5048448833 63025.9560730867 62827.9224182023 62870.1796619116
#TimeLost as linear

#PT
# QIC          QIC.1            QIC.2            QIC.3
# model2A             POD0          POD2a            POD2b            POD2c
# QIC2A   11571.2264645272 11578.9417632632 11585.353504737 135053.474264298
#TimeLost as factor

#QN
# QIC            QIC.1            QIC.2            QIC.3
# model2A             POD0            POD2a            POD2b            POD2c
# QIC2A   16378.2460093917 16427.7867543528 16382.0099195888 16396.6153574563
#TimeLost as linear

#BD
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
# QIC2A   23650.3441943795 23875.7697393579 23662.5320430141 23723.9228871438
#TimeLost as linear

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
if (site == "CB"){
  #CB (with Year as factor)
  #The initial full model is:
  POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without AvgDayMat
  POD3b = geeglm(PreAbs ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Year
  POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Timelost
  POD3d = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model3A = c("POD0","POD3a","POD3b","POD3c","POD3d")
  QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1])
  QICmod3A<-data.frame(rbind(model3A,QIC3A))
  QICmod3A
  #CB
  #QIC            QIC.1            QIC.2            QIC.3            QIC.4
  # model3A             POD0            POD3a            POD3b           POD3c            POD3d
  #QIC3A   62743.5048448833 60915.8541569011 61277.9740202374 61977.6718987442 60920.4684434373
  #Remove TimeLost
  
  #The  full model without TimeLost is:
  POD3e = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without AvgDayMat
  POD3f = geeglm(PreAbs ~ as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Year
  POD3g = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model3B = c("POD0","POD3e","POD3f","POD3g")
  QIC3B = c(QIC(POD0)[1],QIC(POD3e)[1],QIC(POD3f)[1],QIC(POD3g)[1])
  QICmod3B<-data.frame(rbind(model3B,QIC3B))
  QICmod3B
  #CB
  #QIC            QIC.1           QIC.2            QIC.3
  #model3B             POD0            POD3e           POD3f            POD3g
  #QIC3B   62743.5048448833 60920.4684434373 61291.184631157 61896.3522746218
  #Full model is the best, year than avgdaymat
}

if (site == 'BD' | site == 'QN' | site == 'PT'){
#Other sites (without year, timeLost as linear)
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3c = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#QN
# QIC            QIC.1            QIC.2            QIC.3
# model3A             POD0            POD3a            POD3b            POD3c
# QIC3A   16378.2460093917 15736.326454506 16382.0099195888 15734.1555437428
#Remove Time Lost as a variable. Final model is POD3C.

#BD
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
# QIC3A   23650.3441943795 23021.3485765059 23662.5320430141 23011.7382084554
#Remove Time Lost as a variable. Final model is POD3C.

#PT
# QIC           QIC.1           QIC.2            QIC.3
# model3A             POD0           POD3a           POD3b            POD3c
# QIC3A   11571.2264645272 11422.908619339 11585.353504737 11410.1658116943
#Remove Time Lost as variable. Final model is POD3C.
}

# Step 7: Finalize Model --------------------------------------------------
#In descending order:
#CB
#Year
#AvgDaMat
#For PT,QN,BD only AvgDayMat was a significant variable, so the order doesn't matter...

#Year as factor
if (site == 'CB'){
  #CB
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ as.factor(Year)+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
} else {
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}


# STEP 8: Interpreting the summary of the model --------------------------
# How to interpret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reported by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero
anova(PODFinal)
#CB
# as.factor(Year)  7 41.9   5.4e-07 ***
#   AvgDayMat        2  2.6      0.28    

#BD
# Df     X2 P(>|Chi|)  
# AvgDayMat  2 16.6   0.00024 ***

#PT
# Df   X2 P(>|Chi|)   
# AvgDayMat  2 7.3137   0.02581 *
  
#QN
# Df   X2 P(>|Chi|)  
# AvgDayMat  2 37.3   7.8e-09 ***

#Remove variables from model that were not significant after running anova
#CB
if (site == 'CB'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
anova(PODFinal)
#CB
#as.factor(Year)  7 41.9   5.4e-07 ***

#Save model output
filename = paste(saveWorkspace,site,'_SiteSpecificModelSummary.txt',sep="")
sink(filename)
summary(PODFinal)
anova(PODFinal)
sink(file = NULL)

# Step 9: Construction of the ROC curve    --------------------------------
pr <- predict(PODFinal, type="response")  
pred <- prediction(pr,SiteHourTableB$PreAbs) 
perf <- performance(pred, measure="tpr", x.measure="fpr")   
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5))
#This creates a ROC plot
#Interpretting ROC curves -- 

# Choice of the best cut-off probability

y<-as.data.frame(perf@y.values)
x<-as.data.frame(perf@x.values)
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve

#Find max and position of max
max = max(d,na.rm=TRUE)
names(d)[1] <- 'Col1'
position = which.max(d$Col1)

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
cutoff = alpha[position,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve

# This value can now be used to build the confusion matrix:

DATA<-matrix(0,dim(SiteHourTableB)[1] ,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim()) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:dim(SiteHourTableB)[1]                                    # the first column is filled with an ID value that is unique for each row
DATA$Observed<-SiteHourTableB$PreAbs                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinal,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:

#CB
#observed
#predicted     1     0
#1 13494 10046
#0  7008 14806

#PT
#observed
#predicted     1     0
# 1 1250 6961
# 0  612 7405

#QN
# observed
# predicted    1    0
# 1 2150 5696
# 0 1065 7815

#BD
#observed
#predicted    1    0
# 1 7916 6167
# 0  758 2218

# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")


# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
save.image(file = fileName)