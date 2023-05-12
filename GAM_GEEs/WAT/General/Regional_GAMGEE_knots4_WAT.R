### Modelling habitat preference of sperm whales using Generalized Additive Models with Generalized Estimating Equations ###
### Script adapted from Pirotta et al. (2011) and Benjamins ###  
### 9 Models total:
  #Site specific models: 'HZ','OC','NC','BC','WC','GS','BP','BS','JAX'
  #Region specific models: North + South
  #Big model: all 9 sites

# All the libraries have to be installed prior to their utilization (see R help on library installation)

library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
library(splines)         # to construct the B-splines within a GEE-GLM
library(tidyverse)       # because it literally does everything
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

region = 'South' #specify the region of interest
Region = 'WAT'
GDrive= 'G'

# Step 1: Load the Data -----------------------------------------------------------
dir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites", sep="")
saveDir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/Plots/",region, sep="")
saveWorkspace = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
fileName = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(region,Region))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12

#If it's a leap year, delete julian day 366
SiteHourTable = SiteHourTable[!(SiteHourTable$Julian==366),]

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region
# The column Julian represents the day of the year and the column year represents the year of recording.

# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
BlockMod<-glm(PreAbs~ bs(Julian, k = 4)+
                as.factor(Year)+
                as.factor(Site)+
                TimeLost, data=SiteHourTable,family=binomial)

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
#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockMod)

#North
#bs(Julian, k = 4)   4205.8  4  < 2.2e-16 ***
#as.factor(Year)     2438.6  4  < 2.2e-16 ***
#as.factor(Site)     6293.7  4  < 2.2e-16 ***
#TimeLost             118.7  1  < 2.2e-16 ***

#South
#bs(Julian, k = 4)      221  4     <2e-16 ***
#as.factor(Year)         76  3     <2e-16 ***
#as.factor(Site)       1115  3     <2e-16 ***
#TimeLost                 7  1     0.0079 ** 

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
GLM = glm(PreAbs~bs(Julian)+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
VIF(GLM)

#North
#bs(Julian)      1.356390  3        1.052117
#TimeLost        1.012317  1        1.006140
#as.factor(Site) 1.117421  4        1.013975
#as.factor(Year) 1.433952  4        1.046085

#South
#bs(Julian)      1.41  3            1.06
#TimeLost        1.01  1            1.00
#as.factor(Site) 1.02  3            1.00
#as.factor(Year) 1.41  3            1.06

# Step 5: Model Selection - Covariate Preparation -------------------------
#Skipped this step for the modified code because I know I'm going to use Julian day as a covariance-variance matrix, year as a factor...
#when applicable and TimeLost is not going to be included because it's not significant.
# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ as.factor(Site)+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3c = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#Without Year
POD3d = geeglm(PreAbs ~ AvgDayMat+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c","POD3d")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A

#North
#QIC            QIC.1            QIC.2           QIC.3            QIC.4
#model3A             POD0            POD3a            POD3b            POD3c            POD3d
#QIC3A   259944.799358274 237873.379451355 245365.13585481 250197.394630712 241294.183321952
#Full model is the best
#Model Order - Site, Julian Day, Year

#South
#QIC           QIC.1           QIC.2            QIC.3            QIC.4
#model3A             POD0           POD3a           POD3b            POD3c            POD3d
#QIC3A   41361.5688974939 40021.468102432 40201.598831496 41120.7347941663 40091.3963344378
#Full model is the best
#Model Order - Site, Julian Day, Year

# Step 7: Finalize Model --------------------------------------------------
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ as.factor(Site)+AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinal)

#North
#as.factor(Site)  4 3287 < 2.2e-16 ***
#AvgDayMat        2 1485 < 2.2e-16 ***
#as.factor(Year)  4  937 < 2.2e-16 ***

#South
#as.factor(Site)  3 577   < 2e-16 ***
#AvgDayMat        2  78   < 2e-16 ***
#as.factor(Year)  3  35   1.1e-07 ***

filename = paste(saveWorkspace,region,'_RegionalModelSummary.txt',sep="")
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
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45? line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45? line and the ROC curve

#Find max and position of max
max = max(d,na.rm=TRUE)
names(d)[1] <- 'Col1'
position = which.max(d$Col1)


alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
cutoff = alpha[position,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45? line and the curve

# This value can now be used to build the confusion matrix:

DATA<-matrix(0,dim(SiteHourTableB)[1] ,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim()) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:dim(SiteHourTableB)[1]                                    # the first column is filled with an ID value that is unique for each row
DATA$Observed<-SiteHourTableB$PreAbs                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinal,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_modified.RData',sep="")
save.image(file = fileName)