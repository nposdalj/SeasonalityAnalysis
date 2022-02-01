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

region = 'GOA' #specify the region of interest

# Step 1: Load the Data -----------------------------------------------------------
dir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/Plots/",region, sep="")
fileName = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="")
saveWorkspace = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(region,Region))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region (GOA or BSAI).
# The column Julian represents the day of the year and the column year represents the year of recording.

# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals

#Social Groups
if (region == 'BSAI'){
  BlockModF<-glm(PreAbsF~
                   bs(Julian,k=4)+
                   TimeLost+
                   as.factor(Site)
                 ,data=SiteHourTable,family=binomial)
}else{
  BlockModF<-glm(PreAbsF~
                   bs(Julian,k=4)+bs(Year,k=4)+
                   as.factor(Site)+TimeLost
                 ,data=SiteHourTable,family=binomial)
}

#Juveniles
if (region == 'BSAI'){
  BlockModJ<-glm(PreAbsJ~
                   bs(Julian,k=4)+
                   TimeLost+
                   as.factor(Site)
                 ,data=SiteHourTable,family=binomial)
}else{
BlockModJ<-glm(PreAbsJ~
                 bs(Julian)+
                 TimeLost+
                 as.factor(Site)+
                 bs(Year,k=4)
               ,data=SiteHourTable,family=binomial)
}

#Males
if (region == 'BSAI'){
  BlockModM<-glm(PreAbsM~
                   bs(Julian,k=4)+
                   TimeLost+
                   as.factor(Site)
                 ,data=SiteHourTable,family=binomial)
}else{
  BlockModM<-glm(PreAbsM~
                  bs(Julian,k=4)+
                  TimeLost+
                  as.factor(Site)+
                  bs(Year,k=4)
                ,data=SiteHourTable,family=binomial)
}

#Social Groups
ACFF = acf(residuals(BlockModF), lag.max = 2000, ylim=c(0,0.1))
CIF = ggfortify:::confint.acf(ACFF)
ACFidxF = which(ACFF[["acf"]] < CIF, arr.ind=TRUE)
ACFvalF = ACFidxF[1]

#Mid Size
ACFJ = acf(residuals(BlockModJ), lag.max = 2000, ylim=c(0,0.1))
CIJ = ggfortify:::confint.acf(ACFJ)
ACFidxJ = which(ACFJ[["acf"]] < CIJ, arr.ind=TRUE)
ACFvalJ = ACFidxJ[1]

#Males
ACFM = acf(residuals(BlockModM), lag.max = 2000, ylim=c(0,0.1))
CIM = ggfortify:::confint.acf(ACFM)
ACFidxM = which(ACFM[["acf"]] < CIM, arr.ind=TRUE)
ACFvalM = ACFidxM[1]

# #Females
# acf(residuals(BlockModF), lag.max = 2000, ylim=c(0,0.1))
# acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1), xlim =c(0,50)) 
# if (region == 'BSAI'){
#   ACFvalF = 220
# }else{
#   ACFvalF = 35
# }
# 
# #Juveniles
# acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
# acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1), xlim =c(50,100))
# if (region == 'BSAI'){
#   ACFvalJ = 162
# }else{
#   ACFvalJ = 53
# }
# 
# #Males
# acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
# acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1), xlim =c(20,40))
# if (region == 'BSAI'){
#   ACFvalM = 139
# }else{
#   ACFvalM = 31
# }

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
timeseriesF = rename(timeseriesF, tbin = date)
SiteHourTableB = left_join(SiteHourTable,timeseriesF,by="tbin")

#Make blocks continuous 
gaps = check4Gaps(SiteHourTableB$blockF,tol=1)
UnBlock = as.data.frame(unique(SiteHourTableB$blockF))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksF = UnBlock$sequence[match(SiteHourTableB$blockF,UnBlock$`unique(SiteHourTableB$blockF)`)]

difference = diff(SiteHourTableB$BlocksF) #find difference between rows
gapsCont = check4Gaps(SiteHourTableB$BlocksF)
SiteHourTableB$WavesF = rep(0,nrow(SiteHourTableB)) #make space for waves

for (i in 1:nrow(gapsCont)){
  StartB = gapsCont$beg.indx[i]
  EndB = gapsCont$end.indx[i]
  Num = length(StartB:EndB)
  SiteHourTableB$WavesF[StartB:EndB] = seq.int(1,Num)
}

#Juveniles
timeseriesJ = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseriesJ)/ACFvalJ)), times=1, each=ACFvalJ)
divdiff = nrow(timeseriesJ) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseriesJ$blockJ = c(preBlock,lastVec)
timeseriesJ = rename(timeseriesJ, tbin = date)
SiteHourTableB = left_join(SiteHourTableB,timeseriesJ,by="tbin")

#Make blocks continuous 
gaps = check4Gaps(SiteHourTableB$blockJ,tol=1)
UnBlock = as.data.frame(unique(SiteHourTableB$blockJ))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksJ = UnBlock$sequence[match(SiteHourTableB$blockJ,UnBlock$`unique(SiteHourTableB$blockJ)`)]

difference = diff(SiteHourTableB$BlocksJ) #find difference between rows
gapsCont = check4Gaps(SiteHourTableB$BlocksJ)
SiteHourTableB$WavesJ = rep(0,nrow(SiteHourTableB)) #make space for waves

for (i in 1:nrow(gapsCont)){
  StartB = gapsCont$beg.indx[i]
  EndB = gapsCont$end.indx[i]
  Num = length(StartB:EndB)
  SiteHourTableB$WavesJ[StartB:EndB] = seq.int(1,Num)
}

#Males
timeseriesM = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseriesM)/ACFvalM)), times=1, each=ACFvalM)
divdiff = nrow(timeseriesM) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseriesM$blockM = c(preBlock,lastVec)
timeseriesM = rename(timeseriesM, tbin = date)
SiteHourTableB = left_join(SiteHourTableB,timeseriesM,by="tbin")

#Make blocks continuous 
gaps = check4Gaps(SiteHourTableB$blockM,tol=1)
UnBlock = as.data.frame(unique(SiteHourTableB$blockM))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksM = UnBlock$sequence[match(SiteHourTableB$blockM,UnBlock$`unique(SiteHourTableB$blockM)`)]

difference = diff(SiteHourTableB$BlocksM) #find difference between rows
gapsCont = check4Gaps(SiteHourTableB$BlocksM)
SiteHourTableB$WavesM = rep(0,nrow(SiteHourTableB)) #make space for waves

for (i in 1:nrow(gapsCont)){
  StartB = gapsCont$beg.indx[i]
  EndB = gapsCont$end.indx[i]
  Num = length(StartB:EndB)
  SiteHourTableB$WavesM[StartB:EndB] = seq.int(1,Num)
}

# Step 3: ANOVA to Check for Significance of Variables --------------------
#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockModF)
#GOA
# bs(Julian, k = 4)   218.14  4  < 2.2e-16 ***
#   bs(Year, k = 4)      55.96  3  4.285e-12 ***
#   as.factor(Site)     524.01  4  < 2.2e-16 ***
#   TimeLost              3.66  1    0.05557 .  

#BSAI
# bs(Julian, k = 4)     1097  4     <2e-16 ***
#   TimeLost                 1  1       0.34    
# as.factor(Site)          0  1       0.72 

Anova(BlockModJ)
#GOA
# bs(Julian)        1513.5  3    < 2e-16 ***
#   TimeLost             5.6  1    0.01776 *  
#   as.factor(Site)   4022.6  4    < 2e-16 ***
#   bs(Year, k = 4)    305.1  3    < 2e-16 ***

#BSAI
# bs(Julian, k = 4)    316.1  4    < 2e-16 ***
#   TimeLost               0.4  1       0.52    
# as.factor(Site)       34.4  1    4.4e-09 ***

Anova(BlockModM)
#GOA
# bs(Julian, k = 4)   1392.5  4  < 2.2e-16 ***
#   TimeLost               8.6  1   0.003397 ** 
#   as.factor(Site)     3480.0  4  < 2.2e-16 ***
#   bs(Year, k = 4)      553.9  3  < 2.2e-16 ***

#BSAI
# bs(Julian, k = 4)      338  4    < 2e-16 ***
#   TimeLost                 1  1       0.37    
# as.factor(Site)         31  1    2.5e-08 ***

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
#Social Groups
if (region == 'BSAI'){
  GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Site),family=binomial,data=SiteHourTableB)
}else{
  GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)
}
VIF(GLMF)
#BSAI
# Julian        TimeLost as.factor(Site) 
# 1               1               1 

#GOA
# GVIF Df GVIF^(1/(2*Df))
# Julian          1.111702  1        1.054373
# TimeLost        1.009662  1        1.004819
# as.factor(Site) 2.438748  4        1.117882
# Year            2.569866  1        1.603080

#Mid-Size
if (region == 'BSAI'){
  GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
}else{
  GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)
}
VIF(GLMJ)
#BSAI
# GVIF Df GVIF^(1/(2*Df))
# Julian          1.47  1            1.21
# TimeLost        1.00  1            1.00
# as.factor(Site) 2.25  1            1.50
# as.factor(Year) 2.95  2            1.31

#GOA
# GVIF Df GVIF^(1/(2*Df))
# Julian          1.091565  1        1.044780
# TimeLost        1.006564  1        1.003277
# as.factor(Site) 1.275524  4        1.030887
# Year            1.379046  1        1.174328

#Male
if (region == 'BSAI'){
  GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Site),family=binomial,data=SiteHourTableB)
}else{
  GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)
}
VIF(GLMM)
#BSAI
# Julian        TimeLost as.factor(Site) 
# 1               1               1 

#GOA
# GVIF Df GVIF^(1/(2*Df))
# Julian          1.068123  1        1.033500
# TimeLost        1.006136  1        1.003064
# as.factor(Site) 1.207081  4        1.023804
# Year            1.281787  1        1.132160

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

#Julian Day
POD0fa = geeglm(PreAbsF ~ mSpline(Julian,
                                  knots=quantile(Julian, probs=c(0.333,0.666)),
                                  Boundary.knots=c(1,366),
                                  periodic=T), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD0fb = geeglm(PreAbsF ~ AvgDayMatF, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model0fA<-c("POD0f", "POD0fa", "POD0fb")
QIC0fA<-c(QIC(POD0f)[1],QIC(POD0fa)[1],QIC(POD0fb)[1])
QICmod0fA<-data.frame(rbind(model0fA,QIC0fA))
QICmod0fA
#GOA
#QIC           QIC.1            QIC.2
#model0fA            POD0f          POD0fa           POD0fb
# QIC0fA   9105.67439809655 9036.80148694459 9039.43051192032
#Julian day as a covariance matrix

#BSAI
#QIC            QIC.1           QIC.2
#model0fA            POD0f           POD0fa          POD0fb
#QIC0fA   3403.79550617836 2371.14317724232 2370.83782002882
#Julian day as a covariance matrix

if (region == 'GOA'){
  #Year
  POD1fa = geeglm(PreAbsF ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  POD1fb = geeglm(PreAbsF ~ Year, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  POD1fc = geeglm(PreAbsF ~ mSpline(Year,
                                    knots=quantile(Year, probs=c(0.333,0.666)),
                                    Boundary.knots=c(2011,2019)), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  model1fA<-c("POD0f", "POD1fa", "POD1fb","POD1fc")
  QIC1fA<-c(QIC(POD0f)[1],QIC(POD1fa)[1],QIC(POD1fb)[1],QIC(POD1fc)[1])
  QICmod1fA<-data.frame(rbind(model1fA,QIC1fA))
  QICmod1fA
  #GOA
  #QIC           QIC.1            QIC.2            QIC.3
  #model1fA            POD0f          POD1fa           POD1fb           POD1fc
  # QIC1fA   9105.67439809655 9096.510787458 9104.6130344458 9111.44877706165
  #Year as mSpline.
}

#TimeLost
POD2fa = geeglm(PreAbsF ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fb = geeglm(PreAbsF ~ TimeLost, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fc = geeglm(PreAbsF ~ mSpline(TimeLost,
                                  knots=4), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model2fA<-c("POD0f", "POD2fa", "POD2fb","POD2fc")
QIC2fA<-c(QIC(POD0f)[1],QIC(POD2fa)[1],QIC(POD2fb)[1],QIC(POD2fc)[1])
QICmod2fA<-data.frame(rbind(model2fA,QIC2fA))
QICmod2fA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
# QIC2fA   9105.67439809655 9192.95914527949 9134.26189400808 57861.0556253454
#TimeLost as linear

#BSAI
#QIC            QIC.1           QIC.2            QIC.3
#model2fA            POD0f           POD2fa          POD2fb           POD2fc
#
#TimeLost as linear

#Juvenile
POD0j<-geeglm(PreAbsJ ~ 1, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Julian Day
POD0ja = geeglm(PreAbsJ ~ mSpline(Julian,
                                  knots=quantile(Julian, probs=c(0.333,0.666)),
                                  Boundary.knots=c(1,366),
                                  periodic=T), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD0jb = geeglm(PreAbsJ ~ AvgDayMatJ, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model0jA<-c("POD0j", "POD0ja", "POD0jb")
QIC0jA<-c(QIC(POD0j)[1],QIC(POD0ja)[1],QIC(POD0jb)[1])
QICmod0jA<-data.frame(rbind(model0jA,QIC0jA))
QICmod0jA
#GOA
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   68226.5136442113 66959.5429913374 66981.0738821536
#Julian day as a covariance matrix.

#BSAI
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   14901.4969508935 14559.7509115965 14572.6290773741
#Julian day as a covariance matrix.

if (region == 'GOA'){
#Year
POD1ja = geeglm(PreAbsJ ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jb = geeglm(PreAbsJ ~ Year, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jc = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model1jA<-c("POD0j", "POD1ja", "POD1jb","POD1jc")
QIC1jA<-c(QIC(POD0j)[1],QIC(POD1ja)[1],QIC(POD1jb)[1],QIC(POD1jc)[1])
QICmod1jA<-data.frame(rbind(model1jA,QIC1jA))
QICmod1jA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model1jA            POD0j           POD1ja           POD1jb           POD1jc
# QIC1jA   68226.5136442113 66342.4628844564 68224.1712143521 66453.7626840511
#Year as mSpline
}

#TimeLost
POD2ja = geeglm(PreAbsJ ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jb = geeglm(PreAbsJ ~ TimeLost, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jc = geeglm(PreAbsJ ~ mSpline(TimeLost,
                                  knots=4), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model2jA<-c("POD0j", "POD2ja", "POD2jb","POD2jc")
QIC2jA<-c(QIC(POD0j)[1],QIC(POD2ja)[1],QIC(POD2jb)[1],QIC(POD2jc)[1])
QICmod2jA<-data.frame(rbind(model2jA,QIC2jA))
QICmod2jA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
# QIC2jA   68226.5136442113 68314.5554993268 68224.2885803012 68229.5597270708
#TimeLost as linear

#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
#QIC2jA   14901.4969508935 180933.116344926 14914.3869437389 181007.820011848
#TimeLost as linear.

#Male
POD0m<-geeglm(PreAbsM ~ 1, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

#Julian Day
POD0ma = geeglm(PreAbsM ~ mSpline(Julian,
                                  knots=quantile(Julian, probs=c(0.333,0.666)),
                                  Boundary.knots=c(1,366),
                                  periodic=T), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD0mb = geeglm(PreAbsM ~ AvgDayMatJ, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model0mA<-c("POD0m", "POD0ma", "POD0mb")
QIC0mA<-c(QIC(POD0m)[1],QIC(POD0ma)[1],QIC(POD0mb)[1])
QICmod0mA<-data.frame(rbind(model0mA,QIC0mA))
QICmod0mA
#GOA
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
# QIC0mA   58066.9036081687 56460.4232524141 56461.4027284332
#Julian day as a covariance matrix.

#BSAI
#QIC            QIC.1           QIC.2
#model0mA            POD0m           POD0ma          POD0mb
# QIC0mA   16089.97000828 15873.8160740996 15870.3286691748
##Julian day as a covariance matrix.

if (region == 'GOA'){
#Year
POD1ma = geeglm(PreAbsM ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mb = geeglm(PreAbsM ~ Year, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mc = geeglm(PreAbsM ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model1mA<-c("POD0m", "POD1ma", "POD1mb","POD1mc")
QIC1mA<-c(QIC(POD0m)[1],QIC(POD1ma)[1],QIC(POD1mb)[1],QIC(POD1mc)[1])
QICmod1mA<-data.frame(rbind(model1mA,QIC1mA))
QICmod1mA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model1mA            POD0m           POD1ma           POD1mb           POD1mc
# QIC1mA   58066.9036081687 55717.5796140073 57701.7790312048 55862.6272306594
#Year as a mSpline
}

#TimeLost
POD2ma = geeglm(PreAbsM ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mb = geeglm(PreAbsM ~ TimeLost, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mc = geeglm(PreAbsM ~ mSpline(TimeLost,
                                  knots=4), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model2mA<-c("POD0m", "POD2ma", "POD2mb","POD2mc")
QIC2mA<-c(QIC(POD0m)[1],QIC(POD2ma)[1],QIC(POD2mb)[1],QIC(POD2mc)[1])
QICmod2mA<-data.frame(rbind(model2mA,QIC2mA))
QICmod2mA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
# QIC2mA   58066.9036081687 58155.3041102724 58077.5294636147 58090.1142780062
#TimeLost as linear.

#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
#QIC2mA   16089.97000828 16125.5706744776 16098.2949857906 16116.4850683141
#TimeLost as linear.

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.

if (region == 'GOA'){
#Social Groups
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                            knots=quantile(Year, probs=c(0.333,0.666)),
                                            Boundary.knots=c(2011,2019))+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fc = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                                              Boundary.knots=c(2011,2019))+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Site
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                                      knots=quantile(Year, probs=c(0.333,0.666)),
                                                      Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#Without Year
POD3fdd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd","POD3fdd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1],QIC(POD3fdd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#GOA
#QIC            QIC.1           QIC.2            QIC.3            QIC.4            QIC.5
#model3fA            POD0f           POD3fa          POD3fb           POD3fc           POD3fd          POD3fdd
# QIC3fA   9105.67439809655 56007.6079723919 9306.62063602652 60304.5931527526 9071.42459773896 8885.10600320385
#Remove site.

#The  full model without Site is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fg = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without TimeLost
POD3fgg = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fB = c("POD0f","POD3fe","POD3ff","POD3fg","POD3fgg")
QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1],QIC(POD3fgg)[1])
QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
QICmod3fB
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fB            POD0f           POD3fe           POD3ff           POD3fg          POD3fgg
# QIC3fB   9105.67439809655 9071.42459773896 9140.39772395276 9060.57604515016 9048.22689310987
#Remove TimeLost

#The full model without TimeLost is:
POD3fh = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fi = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fj = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3fC = c(QIC(POD0f)[1],QIC(POD3fh)[1],QIC(POD3fi)[1],QIC(POD3fj)[1])
QICmod3fC<-data.frame(rbind(model3fC,QIC3fC))
QICmod3fC
#GOA
#QIC            QIC.1           QIC.2            QIC.3
#model3fC            POD0f           POD3fh          POD3fi           POD3fj
# QIC3fC   9105.67439809655 9048.22689310987 9111.44877706165 9039.43051192032
#Full model is best

#Model Order
#Year
#AvgDayMat

}else{
  #The initial full model is:
  POD3fa = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  #without AvgDayMat
  POD3fb = geeglm(PreAbsF ~ TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  #without Timelost
  POD3fc = geeglm(PreAbsF ~ AvgDayMatF +as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  #without Site
  POD3fd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd")
  QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1])
  QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
  QICmod3fA
  #BSAI
  #QIC            QIC.1            QIC.2            QIC.3            QIC.4
  #model3fA            POD0f           POD3fa           POD3fb           POD3fc           POD3fd
  # QIC3fA   3403.79550617836 2370.80426520256 3338.49682826181 2370.83708484994 2370.80426520256
  #Remove Time Lost
  
  #The  full model without TimeLost is:
  POD3fe = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  #without AvgDayMat
  POD3ff = geeglm(PreAbsF ~ as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  #without Site
  POD3fg = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  model3fB = c("POD0f","POD3fe","POD3ff","POD3fg")
  QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1])
  QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
  QICmod3fB
  #QIC            QIC.1            QIC.2           QIC.3
  #model3fB            POD0f           POD3fe           POD3ff          POD3fg
  # QIC3fB   3403.79550617836 2370.83708484994 3338.52819557516 2370.83782002882
  #The full model is the best.
  #AvgDayMat first, then site.
}

#Mid-Size
if (region == 'GOA'){
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Site
POD3jdd = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd","POD3jdd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1],QIC(POD3jdd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#GOA
# QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
# model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd          POD3jdd
# QIC3jA   66581.7159280241 909612.557218144 1259861.02410151 845128.884936515 1214164.40247372 63760.1730377111
#Remove Site.

#The  full model without Site is:
POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jf = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jg = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#Without TimeLost
POD3jee = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jB = c("POD0j","POD3je","POD3jf","POD3jg","POD3jee")
QIC3jB = c(QIC(POD0j)[1],QIC(POD3je)[1],QIC(POD3jf)[1],QIC(POD3jg)[1],QIC(POD3jee)[1])
QICmod3jB<-data.frame(rbind(model3jB,QIC3jB))
QICmod3jB
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3jB            POD0j           POD3je           POD3jf           POD3jg          POD3jee
#QIC3jB   66581.7159280241 63760.1730377111 65085.3002445263 65218.7715841914 63741.470492651
#Remove TimeLost

#The full model without TimeLost is:
POD3jh = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3ji = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jj = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3jC = c(QIC(POD0j)[1],QIC(POD3jh)[1],QIC(POD3ji)[1],QIC(POD3jj)[1])
QICmod3jC<-data.frame(rbind(model3jC,QIC3jC))
QICmod3jC
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model3jC            POD0f           POD3fh           POD3fi           POD3fj
#QIC3jC   66581.7159280241 63741.470492651 65050.1224949519 65205.196089869
#Full model is the best

#Model Order
#Year
#AvgDayMat

}else{
  #The initial full model is:
  POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  #without AvgDayMat
  POD3jb = geeglm(PreAbsJ ~ TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  #without Timelost
  POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  #without Site
  POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd")
  QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1])
  QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
  QICmod3jA
  #BSAI
  #QIC            QIC.1            QIC.2            QIC.3            QIC.4
  #model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd
  #QIC3jA   14901.4969508935 14561.3121248878 14828.3226892144 14551.3819047153 14582.967436424
  #Remove TimeLost.
  
  #The  full model without TimeLost is:
  POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  #without AvgDayMat
  POD3jf = geeglm(PreAbsJ ~ as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  #without Site
  POD3jg = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  model3jB = c("POD0j","POD3je","POD3jf","POD3jg")
  QIC3jB = c(QIC(POD0j)[1],QIC(POD3je)[1],QIC(POD3jf)[1],QIC(POD3jg)[1])
  QICmod3jB<-data.frame(rbind(model3jB,QIC3jB))
  QICmod3jB
  #BSAI
  #QIC            QIC.1            QIC.2            QIC.3
  #model3jB            POD0j           POD3je           POD3jf           POD3jg
  #QIC3jB   14901.4969508935 14551.3819047153 14816.6138985244 14572.6290773741
  #Full model is the best.
  #AvgDayMat first, then site.
}

if (region == 'GOA'){
#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Site
POD3mdd = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md","POD3mdd")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1],QIC(POD3mdd)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#GOA
#QIC            QIC.1            QIC.2            QIC.3           QIC.4            QIC.5
#model3mA            POD0m           POD3ma           POD3mb           POD3mc          POD3md          POD3mdd
#QIC3mA   57273.6041843083 1137626.03605746 667077.747579427 53708.9487038409 833290.259843814 54450.6165691385
#Without Year

#The  full model without Year is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Site
POD3mg = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#Without TimeLost
POD3mgg = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3me","POD3mf","POD3mg","POD3mgg")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1],QIC(POD3mgg)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
QICmod3mB
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3mB            POD0m           POD3me           POD3mf           POD3mg          POD3mgg
#QIC3mB   57273.6041843083 53708.9487038409 55085.4284209363 55959.3710282636 53672.2327165976
#Remove TimeLost

#The full model without TimeLost is:
POD3mh = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mi = geeglm(PreAbsM ~ as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Site
POD3mj = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3mh","POD3mi","POD3mj")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3mh)[1],QIC(POD3mi)[1],QIC(POD3mj)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
QICmod3mB
#GOA
# QIC            QIC.1            QIC.2            QIC.3
# model3mB            POD0m           POD3mh           POD3mi           POD3mj
# QIC3mB   57273.6041843083 53672.2327165976 55043.9665372996 55939.0467556249
#Full model is best.

#Model order
#Site
#AvgDayMat

}else{
  #The initial full model is:
  POD3ma = geeglm(PreAbsM ~ AvgDayMatM+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without AvgDayMat
  POD3mb = geeglm(PreAbsM ~ TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without Site
  POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without Timelost
  POD3md = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md")
  QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1])
  QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
  QICmod3mA
  #BSAI
  #QIC            QIC.1            QIC.2            QIC.3            QIC.4
  #model3mA            POD0m           POD3ma           POD3mb           POD3mc           POD3md
  #QIC3mA   16089.97000828 15842.4680660985 16060.1507409164 15878.8854855284 15831.2782164232
  #Remove TimeLost.
  
  #The  full model without Timelost is:
  POD3me = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without AvgDayMat
  POD3mf = geeglm(PreAbsM ~ as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without Site
  POD3mg = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  model3mB = c("POD0m","POD3me","POD3mf","POD3mg")
  QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1])
  QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
  QICmod3mB
  #BSAI
  #QIC            QIC.1            QIC.2           QIC.3
  #model3mB            POD0m           POD3me           POD3mf          POD3mg
  #QIC3mB   16089.97000828 15831.2782164232 16048.3274797541 15870.3286691748
  #Full model is best 
  #AvgDayMat first, then site.
}

# Step 7: Finalize Model --------------------------------------------------
#Social Groups
if (region == 'BSAI'){
  dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalF = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  } else {
  dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalF = geeglm(PreAbsF ~ mSpline(Year,
                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                     Boundary.knots=c(2011,2019))+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  PODFinalF_Site = geeglm(PreAbsF ~ mSpline(Year,
                                       knots=quantile(Year, probs=c(0.333,0.666)),
                                       Boundary.knots=c(2011,2019))+AvgDayMatF+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
}

#Mid-Size
if (region == 'BSAI'){
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~ mSpline(Year,
                                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                                  Boundary.knots=c(2011,2019))+AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  PODFinalJ_Site = geeglm(PreAbsJ ~ mSpline(Year,
                                                knots=quantile(Year, probs=c(0.333,0.666)),
                                                Boundary.knots=c(2011,2019))+AvgDayMatJ+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ,data=SiteHourTableB)
}

#Males
if (region == 'BSAI'){
  dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalM = geeglm(PreAbsM ~ as.factor(Site) +AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  PODFinalM_Year = geeglm(PreAbsM ~ as.factor(Site)+AvgDayMatM+mSpline(Year,
                                                                       knots=quantile(Year, probs=c(0.333,0.666)),
                                                                       Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
}

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinalF)
#GOA
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5 22.843 0.0003617 ***
#   AvgDayMatF                                                                                      2 18.963 7.624e-05 ***

#BSAI
# AvgDayMatF       2  4.74     0.093 .  
# as.factor(Site)  1 15.61   7.8e-05 ***

if (region == 'GOA'){
  anova(PODFinalF_Site)
  # mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5 2.3000e+01 0.0003617 ***
  #   AvgDayMatF                                                                                      2 1.9000e+01 7.624e-05 ***
  #   as.factor(Site)                                                                                 4 8.2783e+10 < 2.2e-16 ***
}

anova(PODFinalJ)
#GOA
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5 101.056 < 2.2e-16 ***
#   AvgDayMatJ                                                                                      2  57.234  3.73e-13 ***

#BSAI
# AvgDayMatJ       2 9.21      0.01 *
#   as.factor(Site)  1 1.58      0.21  

if (region == 'GOA'){
anova(PODFinalJ_Site)
  # mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5  101.06 < 2.2e-16 ***
  #   AvgDayMatJ                                                                                      2   57.23  3.73e-13 ***
  #   as.factor(Site)                                                                                 4 1920.28 < 2.2e-16 ***
}

anova(PODFinalM)
#GOA
# as.factor(Site)  4 180.10   < 2e-16 ***
#   AvgDayMatM       2  13.23   0.00134 ** 

#BSAI
# AvgDayMatM       2 18.5276  9.48e-05 ***
#   as.factor(Site)  1  6.3368   0.01183 *  

if (region == 'GOA'){
  anova(PODFinalM_Year)
  # as.factor(Site)                                                                                 4 180.101 < 2.2e-16 ***
  #   AvgDayMatM                                                                                      2  13.230   0.00134 ** 
  #   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5  73.175 2.232e-14 ***
}

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
DATA$Observed<-SiteHourTableB$PreAbsF                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalF,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:
#BSAI - Social Group
# observed
# predicted     1     0
# 1   298  2053
# 0    40 16793

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
DATA$Observed<-SiteHourTableB$PreAbsJ                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalJ,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:
#BSAI - Mid-Size
# observed
# predicted     1     0
# 1  1518  6559
# 0   995 10112

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
DATA$Observed<-SiteHourTableB$PreAbsM                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalM,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:
#BSAI - Males
# observed
# predicted    1    0
# 1 2098 8722
# 0  741 7623

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput_sexClasses.RData',sep="")
save.image(file = fileName)
