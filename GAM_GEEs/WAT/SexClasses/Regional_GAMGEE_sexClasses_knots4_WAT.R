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
library(car) #for ANOVA
library(splines2)       # to use mSpline for the GEEs
library(ggfortify)      # extract confidence interval for ACF plots
library(dplyr)

region = 'North' #specify the region of interest
Region = 'WAT'
GDrive = 'G'

# Step 1: Load the data ---------------------------------------------------
#Hourly data
dir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites", sep="")
saveDir = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/Plots/",region, sep="")
saveWorkspace = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
fileName = paste(GDrive,":/My Drive/",Region,"_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(region,Region))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable = subset(SiteHourTable, Region == region)#subset the table for the site only
SiteHourTable = SiteHourTable[!(SiteHourTable$Julian==366),] #If it's a leap year, delete julian day 366

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region
# The column Julian represents the day of the year and the column year represents the year of recording.


# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
#Females
BlockModF<-glm(PreAbsF~
               bs(Julian, k = 4)+
               TimeLost+
               as.factor(Year)+
                 as.factor(Site)
             ,data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian, k = 4)+
                 TimeLost+
                 as.factor(Year)+
                 as.factor(Site)
               ,data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                 as.factor(Year)+
                 as.factor(Site)+
                 bs(Julian, k = 4),
               data=SiteHourTable,family=binomial)

#Social Groups
ACFF = acf(residuals(BlockModF), lag.max = 2000, ylim=c(-0.1,0.1))
CIF = ggfortify:::confint.acf(ACFF)
CIF_neg = CIF*-1
ACFidxF = which(ACFF[["acf"]] < CIF & ACFF[["acf"]] > CIF_neg , arr.ind=TRUE)
ACFvalF = ACFidxF[1]

#Mid Size
ACFJ = acf(residuals(BlockModJ), lag.max = 2000, ylim=c(-0.1,0.1))
CIJ = ggfortify:::confint.acf(ACFJ)
CIJ_neg = CIJ*-1
ACFidxJ = which(ACFJ[["acf"]] < CIJ & ACFJ[["acf"]] > CIJ_neg, arr.ind=TRUE)
ACFvalJ = ACFidxJ[1]

#Males
ACFM = acf(residuals(BlockModM), lag.max = 200, ylim=c(-0.1,0.1))
CIM = ggfortify:::confint.acf(ACFM)
CIM_neg = CIM*-1
ACFidxM = which(ACFM[["acf"]] < CIM & ACFM[["acf"]] > CIM_neg, arr.ind=TRUE)
if (region =='North'){
  ACFvalM = 172 #Closest to the CI intervals without having a ridiculous number over 1600+
}else{
ACFvalM = ACFidxM[1]}

#create the blocks based on the full timesereies
startDate = SiteHourTable$tbin[1]
endDate = SiteHourTable$tbin[nrow(SiteHourTable)]

#Females
timeseriesF = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseriesF)/ACFvalF)), times=1, each=ACFvalF)
divdiff = nrow(timeseriesF) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseriesF$blockF = c(preBlock,lastVec)
timeseriesF = dplyr::rename(timeseriesF, tbin = date)
SiteHourTableB = left_join(SiteHourTable,timeseriesF,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$blockF))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksF = UnBlock$sequence[match(SiteHourTableB$blockF,UnBlock$`unique(SiteHourTableB$blockF)`)]
difference = diff(SiteHourTableB$BlocksF) #find difference between rows

#Juveniles
timeseriesJ = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseriesJ)/ACFvalJ)), times=1, each=ACFvalJ)
divdiff = nrow(timeseriesJ) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseriesJ$blockJ = c(preBlock,lastVec)
timeseriesJ = dplyr::rename(timeseriesJ, tbin = date)
SiteHourTableB = left_join(SiteHourTableB,timeseriesJ,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$blockJ))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksJ = UnBlock$sequence[match(SiteHourTableB$blockJ,UnBlock$`unique(SiteHourTableB$blockJ)`)]
difference = diff(SiteHourTableB$BlocksJ) #find difference between rows

#Males
timeseriesM = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseriesM)/ACFvalM)), times=1, each=ACFvalM)
divdiff = nrow(timeseriesM) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseriesM$blockM = c(preBlock,lastVec)
timeseriesM = dplyr::rename(timeseriesM, tbin = date)
SiteHourTableB = left_join(SiteHourTableB,timeseriesM,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$blockM))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksM = UnBlock$sequence[match(SiteHourTableB$blockM,UnBlock$`unique(SiteHourTableB$blockM)`)]
difference = diff(SiteHourTableB$BlocksM) #find difference between rows

# Step 3: ANOVA to Check for Significance of Variables --------------------
#ANOVA with car package
Anova(BlockModF)
#North
#bs(Julian, k = 4)      489  4     <2e-16 ***
#TimeLost               145  1     <2e-16 ***
#as.factor(Year)       2653  4     <2e-16 ***
#as.factor(Site)       3237  4     <2e-16 ***

#South
#bs(Julian, k = 4)       29  4    9.5e-06 ***
#TimeLost                 5  1      0.028 *  
#as.factor(Year)         74  3    6.7e-16 ***
#as.factor(Site)        344  3    < 2e-16 ***

Anova(BlockModJ)
#North
#bs(Julian, k = 4)      489  4     <2e-16 ***
#TimeLost               145  1     <2e-16 ***
#as.factor(Year)       2653  4     <2e-16 ***
#as.factor(Site)       3237  4     <2e-16 ***

#South
#bs(Julian, k = 4)       60  4    3.2e-12 ***
#TimeLost                 8  1     0.0061 ** 
#as.factor(Year)         51  3    6.0e-11 ***
#as.factor(Site)        577  3    < 2e-16 ***
  
Anova(BlockModM)
#North
#bs(Julian, k = 4)     2783  4     <2e-16 ***
#TimeLost                 4  1      0.057 .  
#as.factor(Year)        136  4     <2e-16 ***
#as.factor(Site)       1089  4     <2e-16 ***

#South
#bs(Julian, k = 4)     73.1  4    5.1e-15 ***
#TimeLost               2.2  1       0.14    
#as.factor(Year)        6.2  3       0.10    
#as.factor(Site)      154.0  3    < 2e-16 ***

# Step 4: Data Exploration and Initial Analysis ------------------------------------------------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:

  #Females
  GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
  #Males
  GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
  #Juveniles
  GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)
#North
#Julian          1.32  1            1.15
#TimeLost        1.00  1            1.00
#as.factor(Year) 1.32  4            1.03

#Julian          1.22  1            1.10
#TimeLost        1.00  1            1.00
#as.factor(Year) 1.22  4            1.03

#Julian          1.18  1            1.08
#TimeLost        1.00  1            1.00
#as.factor(Year) 1.18  4            1.02

#South
#Julian           1.3  1            1.14
#TimeLost         1.0  1            1.00
#as.factor(Year)  1.3  3            1.04


#Julian          1.35  1            1.16
#TimeLost        1.00  1            1.00
#as.factor(Year) 1.35  3            1.05

#Julian           1.3  1            1.14
#TimeLost         1.0  1            1.00
#as.factor(Year)  1.3  3            1.04

# Step 5: Model Selection - Covariate Preparation -------------------------
# Construct variance-covariance matrices for cyclic covariates:
#Females
AvgDayBasisF <- gam(PreAbsF~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMatF = as.matrix(AvgDayBasisF)

#Juveniles
AvgDayBasisJ <- gam(PreAbsJ~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMatJ = as.matrix(AvgDayBasisJ)

#Males
AvgDayBasisM <- gam(PreAbsM~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMatM = as.matrix(AvgDayBasisM)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.

#Female
POD0f<-geeglm(PreAbsF ~ 1, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)

#Juvenile
POD0j<-geeglm(PreAbsJ ~ 1, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Male
POD0m<-geeglm(PreAbsM ~ 1, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)


# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#Females
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fc = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Site
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)

model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA

#South
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fA            POD0f           POD3fa           POD3fb           POD3fc           POD3fd
#QIC3fA   12806.7635613707 12386.4482332141 12407.9468021846 12440.8229626073 12711.7594856545
#Full model is best - Site, Year, Julian Day

#North
# QIC            QIC.1            QIC.2          QIC.3           QIC.4
# model3fA            POD0f           POD3fa           POD3fb      POD3fc          POD3fd
# QIC3fA   152072.337986146 146451.889676579 146649.198626759  148841.712674249 149655.60754625
#Full model is best - Site, Year, Julian Day

#Juveniles
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Site
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA

#South
#QIC            QIC.1           QIC.2            QIC.3           QIC.4
#model3jA            POD0j           POD3ja          POD3jb           POD3jc          POD3jd
#QIC3jA   15879.2960879808 15270.7800575956 15275.132341962 15311.0636492576 15830.158334417
#Full model - Site, Year, Julian Day

#North
# QIC            QIC.1            QIC.2            QIC.3            QIC.4
# model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd
# QIC3jA   81971.4861594802 74491.1781115864 79771.8295944595 75227.2290293803 76301.4072099103
#Full model - Julian Day, Site, Year

#Males
#The initial full model is:
if (region == 'South'){
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Site
POD3mc = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

model3mA = c("POD0m","POD3ma","POD3mb","POD3mc")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#South
#QIC            QIC.1            QIC.2            QIC.3
#model3mA            POD0m           POD3ma           POD3mb           POD3mc
#QIC3mA   5771.70385452481 5621.78286706419 5632.15957850428 5762.88942847682
#Full model is best - Site, Julian Day
}else{
  POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without AvgDayMat
  POD3mb = geeglm(PreAbsM ~ as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without Year
  POD3mc = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without Site
  POD3md = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  
  model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md")
  QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1])
  QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
  QICmod3mA
  
  #North
  # QIC            QIC.1            QIC.2            QIC.3           QIC.4
  # model3mA            POD0m           POD3ma           POD3mb           POD3mc          POD3md
  # QIC3mA   22781.3047776683 19038.2586405958 21822.8623096594 19131.8809396324 20143.457783147
  #Full model - Julian Day, Site, Year
}

# Step 7: Finalize Model --------------------------------------------------
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalF = geeglm(PreAbsF ~ as.factor(Site)+as.factor(Year)+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)

#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
if (region == 'South'){
  PODFinalJ = geeglm(PreAbsJ ~  as.factor(Site)+as.factor(Year)+AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}else{
  PODFinalJ = geeglm(PreAbsJ ~  AvgDayMatJ+as.factor(Site)+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
if (region == 'South'){
  PODFinalM = geeglm(PreAbsM ~ as.factor(Site)+AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
}else{
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site)+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
}
# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinalF)
#South
#as.factor(Site)  3 100.0   < 2e-16 ***
#as.factor(Year)  3  34.2   1.8e-07 ***
#AvgDayMatF       2  12.1    0.0024 ** 

#North
# as.factor(Site)  4 318.32 < 2.2e-16 ***
# as.factor(Year)  4 191.06 < 2.2e-16 ***
# AvgDayMatF       2  21.96 7.7e-14 ***
  
anova(PODFinalJ)
#South
#as.factor(Site)  3 342   < 2e-16 ***
#as.factor(Year)  3  28   2.9e-06 ***
#AvgDayMatJ       2   6     0.042 *  

#North
# AvgDayMatJ       2 1609.35 < 2.2e-16 ***
# as.factor(Site)  4  473.14 < 2.2e-16 ***
# as.factor(Year)  4  191.72 < 2.2e-16 ***

anova(PODFinalM)
#South
#as.factor(Site)  3 74.6   4.4e-16 ***
#AvgDayMatM       2 11.2    0.0038 ** 

#North
#AvgDayMatM       2 150.2   < 2e-16 ***
#as.factor(Site)  7 199.4   < 2e-16 ***
#as.factor(Year)  4  32.2   1.7e-06 ***

filename = paste(saveWorkspace,region,'_Regional_sexClasses_ModelSummary.txt',sep="")
sink(filename)
summary(PODFinalF)
anova(PODFinalF)
summary(PODFinalJ)
anova(PODFinalJ)
summary(PODFinalM)
anova(PODFinalM)
sink(file = NULL)
  
# Step 9: Construction of the ROC curve    --------------------------------
#Females
prf <- predict(PODFinalF, type="response")  
pred <- prediction(prf,SiteHourTableB$PreAbsF) 
perf <- performance(pred, measure="tpr", x.measure="fpr")   
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5))
#This creates a ROC plot
#Interpretting ROC curves -- 

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
DATA$Observed<-SiteHourTableB$PreAbsF                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalF,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here
  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

#Juveniles
prj <- predict(PODFinalJ, type="response")  
pred <- prediction(prj,SiteHourTableB$PreAbsJ) 
perf <- performance(pred, measure="tpr", x.measure="fpr")   
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5))
#This creates a ROC plot
#Interpretting ROC curves -- 

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
DATA$Observed<-SiteHourTableB$PreAbsJ                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalJ,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

#Males
prm <- predict(PODFinalM, type="response")  
pred <- prediction(prm,SiteHourTableB$PreAbsM) 
perf <- performance(pred, measure="tpr", x.measure="fpr")   
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5))
#This creates a ROC plot
#Interpretting ROC curves -- 

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
DATA$Observed<-SiteHourTableB$PreAbsM                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalM,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here


# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,region,'_Regional_gamgeeOutput_sexClasses.RData',sep="")
save.image(file = fileName)