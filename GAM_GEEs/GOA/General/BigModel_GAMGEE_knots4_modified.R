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

GDrive = 'I'

# Step 1: Load the Data -----------------------------------------------------------
dir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel",sep="")
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots",sep="")
saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel",sep="")
fileName = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = HourTable
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable$Region = as.factor(SiteHourTable$Region)
SiteHourTable$Site = as.factor(SiteHourTable$Site)

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
BlockMod<-glm(PreAbs~
               bs(Julian,k=4)+
               as.factor(Year)+
               as.factor(Region)+
                TimeLost
             ,data=SiteHourTable,family=binomial)

ACF = acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1))
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval = ACFidx[1]
ACFval = 463 #average of all of the sites

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
#Quick ANOVA to check for significance of variables - using the car packageAnova(BlockMod)
Anova(BlockMod)

# LR Chisq Df Pr(>Chisq)    
# bs(Julian, k = 4)     2069  4    < 2e-16 ***
#   as.factor(Year)       4120  8    < 2e-16 ***
#   as.factor(Region)      203  1    < 2e-16 ***
#   TimeLost                32  1    1.5e-08 ***

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
GLM1 = glm(PreAbs ~ bs(Julian) + TimeLost + as.factor(Year) + as.factor(Region), family = binomial, data = SiteHourTableB)
#VIF scores in GLM to work out collinearity:
VIF(GLM1)
# bs(Julian)        1.56  3            1.08
# TimeLost          1.01  1            1.00
# as.factor(Year)   3.72  8            1.09
# as.factor(Region) 2.76  1            1.66

# Step 5: Model Selection - Covariate Preparation -------------------------
# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.
POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The Initial full model:
POD3aa = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3bb = geeglm(PreAbs ~ as.factor(Year)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3cc = geeglm(PreAbs ~ AvgDayMat+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Region
POD3dd = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3AA = c("POD0","POD3aa","POD3bb","POD3cc","POD3dd")
QIC3AA = c(QIC(POD0)[1],QIC(POD3aa)[1],QIC(POD3bb)[1],QIC(POD3cc)[1],QIC(POD3dd)[1])
QICmod3AA<-data.frame(rbind(model3AA,QIC3AA))
QICmod3AA

# QIC            QIC.1            QIC.2            QIC.3            QIC.4
# model3AA             POD0           POD3aa           POD3bb           POD3cc           POD3dd
# QIC3AA   136332.865499456 129076.604470597 130573.869559692 135121.900263247 128284.703049996
#Remove region

#The Initial full model without region is:
POD3aa = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3bb = geeglm(PreAbs ~ as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3cc = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3AA = c("POD0","POD3aa","POD3bb","POD3cc")
QIC3AA = c(QIC(POD0)[1],QIC(POD3aa)[1],QIC(POD3bb)[1],QIC(POD3cc)[1])
QICmod3AA<-data.frame(rbind(model3AA,QIC3AA))
QICmod3AA

# QIC            QIC.1           QIC.2            QIC.3
# model3AA             POD0           POD3aa          POD3bb           POD3cc
# QIC3AA   136332.865499456 128284.703049996 129623.47419353 134495.514519045

#In descending order:
#as.factor(Year)
#AvgDayMat

# Step 7: Finalize Model --------------------------------------------------
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ as.factor(Year)+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero
anova(PODFinal)
# as.factor(Year)  8 119.913 < 2.2e-16 ***
#   AvgDayMat        2  23.728  7.04e-06 ***

#Save model output
filename = paste(saveWorkspace,'_ModelSummary.txt',sep="")
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
cmx(DATA, threshold = cutoff)      
# observed
# predicted     1     0
# 1 22189 22623
# 0 13936 40558
  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,'/BigModel_gamgeeOutput_modified.RData',sep="")
save.image(file = fileName)