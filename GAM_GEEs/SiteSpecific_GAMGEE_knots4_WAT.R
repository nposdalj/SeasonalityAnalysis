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
library(ChemoSpec2D)
library(car)            # to run an ANOVA
library(splines2)       # to use mSpline for the GEEs
library(ggfortify)      # extract confidence interval for ACF plots

site = 'BS' #specify the site of interest

# Step 1: Load the Data -----------------------------------------------------------
dir = paste("H:/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("H:/My Drive/WAT_TPWS_metadataReduced/Plots/",site, sep="")
saveWorkspace = paste("H:/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste("H:/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
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
# The column TimeLost represents the amount of time that we didn't record in each one hour bin


# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
BlockMod<-glm(PreAbs~
               bs(Julian,k=4)+
               TimeLost+
               as.factor(Year)
             ,data=SiteHourTable,family=binomial)

ACF = acf(residuals(BlockMod))
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
names(timeseries)[1] = 'tbin'
SiteHourTableB = left_join(SiteHourTable,timeseries,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$block))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$Blocks = UnBlock$sequence[match(SiteHourTableB$block,UnBlock$`unique(SiteHourTableB$block)`)]
difference = diff(SiteHourTableB$Blocks) #find difference between rows

# Step 3: ANOVA to Check for Significance of Variables --------------------
#ANOVA with car package
Anova(BlockMod) #run this to get p-values
summary(BlockMod)

#BP
# bs(Julian, k = 4)   91.532  4     <2e-16 ***
# TimeLost             1.879  1     0.1705    
# as.factor(Year)     83.482  3     <2e-16 ***
  
#BS
#bs(Julian, k = 4)   321.76  4     <2e-16 ***
#TimeLost              0.06  1      0.802    
#as.factor(Year)     164.22  3     <2e-16 ***

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
  GLM1 = glm(PreAbs ~ bs(Julian) + TimeLost + as.factor(Year), family = binomial, data = SiteHourTableB)
  VIF(GLM1)
  
 # BP
  # bs(Julian)      1.248195  3        1.037641
  # TimeLost        1.003835  1        1.001916
  # as.factor(Year) 1.252221  3        1.038198
  
#BS
#bs(Julian)      1.231767  3        1.035352
#TimeLost        1.007290  1        1.003638
#as.factor(Year) 1.238047  3        1.036230


# Step 5: Model Selection - Covariate/variable Preparation -------------------------

# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,365,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.
POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
  #Year as factor)
  #The initial full model is:
  POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without AvgDayMat
  POD3b = geeglm(PreAbs ~ as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Year
  POD3c = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

  model3A = c("POD0","POD3a","POD3b","POD3c")
  QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1])
  QICmod3A<-data.frame(rbind(model3A,QIC3A))
  QICmod3A
  
  #BS
  #QIC            QIC.1            QIC.2            QIC.3            
  # model3A             POD0            POD3a            POD3b           POD3c      
  #QIC3A   62743.5048448833 60915.8541569011 61277.9740202374 61977.6718987442 
  # model order: year, then avg julian day
  
  #BP
  # QIC            QIC.1            QIC.2            QIC.3
  # model3A             POD0            POD3a            POD3b            POD3c
  # QIC3A   12735.9210928996 12539.6262418463 12615.8897572215 12615.5422645209
  # model order: day, year

# Step 7: Finalize Model --------------------------------------------------

#Year as factor
if (site == 'BP'){
  #BS
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ as.factor(Year)+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
} else {
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}


# STEP 8: Interpreting the summary of the model --------------------------
# How to interpret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reported by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero
anova(PODFinal)
  
  # BS
  # Df   X2 P(>|Chi|)    
  # AvgDayMat        2 33.9   4.4e-08 ***
  #   as.factor(Year)  3 33.0   3.2e-07 ***
  # BP
  # Df     X2 P(>|Chi|)    
  # AvgDayMat        2 38.548 4.261e-09 ***
  # as.factor(Year)  3 25.575 1.170e-05 ***

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

# BS
# observed
# predicted     1     0
# 1   857  6364
# 0  1030 19147

# BP
# observed
# predicted     1     0
# 1  1053 11484
# 0   653 13877

# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")


# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
save.image(file = fileName)