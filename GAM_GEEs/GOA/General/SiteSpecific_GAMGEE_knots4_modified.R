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
library(ggfortify)      # extract confidence interval for ACF plots

site = 'CB' #specify the site of interest
GDrive = 'I'

# Step 1: Load the Data -----------------------------------------------------------
dir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel")
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",site, sep="")
saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
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
names(timeseries)[names(timeseries) == 'date'] = 'tbin'
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
#   TimeLost             1.282  1     0.25782  

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
#Skipped this step for the modified code because I know I'm going to use Julian day as a covariance-variance matrix, year as a factor...
#when applicable and TimeLost is not going to be included because it's not significant.

# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,365,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
if (site == "CB"){
  #CB (with Year as factor)
  #The initial full model is:
  POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without AvgDayMat
  POD3b = geeglm(PreAbs ~ as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Year
  POD3c = geeglm(PreAbs ~ AvgDayMat ,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model3A = c("POD0","POD3a","POD3b","POD3c")
  QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1])
  QICmod3A<-data.frame(rbind(model3A,QIC3A))
  QICmod3A
  #CB
  # QIC            QIC.1            QIC.2            QIC.3
  # model3A           POD0            POD3a            POD3b            POD3c
  # QIC3A   62711.82863537 60895.8941267794 61246.3370753183 61885.8728145065
  #Full model is best
  #Year first, then AvgDayMat
}

if (site == 'BD' | site == 'QN' | site == 'PT'){
#Other sites (without year)
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A

#QN
# QIC            QIC.1
# model3A             POD0            POD3a
# QIC3A   16378.2460093917 15734.1555437428
#Full model is best

#BD
# QIC            QIC.1
# model3A             POD0            POD3a
# QIC3A   23650.3441943795 23011.7382084554
#Full model is best

#PT
# QIC            QIC.1
# model3A             POD0            POD3a
# QIC3A   11564.5271762722 11406.8167171976
#Full model is best

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
# AvgDayMat  2 7.3137   0.02605 *
  
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
#as.factor(Year)  7 32.7     3e-05 ***

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
# 1 1245 7158
# 0  617 7184

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
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_modified.RData',sep="")
save.image(file = fileName)