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
GDrive= 'I'
dir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel",sep="")
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",region, sep="")
fileName = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="")
saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
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
# The column "Region" corresponds to the recording region (GOA or BSAI).
# The column Julian represents the day of the year and the column year represents the year of recording.

# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals

#Social Groups
if (region == 'BSAI'){
  BlockModF<-glm(PreAbsF~
                   bs(Julian)+
                   TimeLost+
                   as.factor(Site)
                 ,data=SiteHourTable,family=binomial)
}else{
  BlockModF<-glm(PreAbsF~
                   bs(Julian,k=4)+
                   as.factor(Year)+
                   as.factor(Site)+
                   TimeLost
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
                 bs(Julian,k=4)+
                 as.factor(Year)+
                 TimeLost+
                 as.factor(Site)
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
                   as.factor(Year)+
                  TimeLost+
                  as.factor(Site)
                ,data=SiteHourTable,family=binomial)
}

#Social Groups
ACFF = acf(residuals(BlockModF), lag.max = 2000, ylim=c(0,0.1))
CIF = ggfortify:::confint.acf(ACFF)
ACFidxF = which(ACFF[["acf"]] < CIF, arr.ind=TRUE)
ACFvalF = ACFidxF[1]
#ACFval calculated from averaging all of the GOA sites (PT, QN, CB)
ACFvalF = 443

#Mid Size
ACFJ = acf(residuals(BlockModJ), lag.max = 2000, ylim=c(0,0.1))
CIJ = ggfortify:::confint.acf(ACFJ)
ACFidxJ = which(ACFJ[["acf"]] < CIJ, arr.ind=TRUE)
ACFvalJ = ACFidxJ[1]
#ACFval calculated from averaging all of the GOA sites (PT, QN, CB)
ACFvalJ = 290

#Males
ACFM = acf(residuals(BlockModM), lag.max = 2000, ylim=c(0,0.1))
CIM = ggfortify:::confint.acf(ACFM)
ACFidxM = which(ACFM[["acf"]] < CIM, arr.ind=TRUE)
ACFvalM = ACFidxM[1]
#ACFval calculated from averaging all of the GOA sites (PT, QN, CB)
ACFvalM = 189

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
names(timeseriesF)[names(timeseriesF) == 'date'] = 'tbin'
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
names(timeseriesJ)[names(timeseriesJ) == 'date'] = 'tbin'
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
names(timeseriesM)[names(timeseriesM) == 'date'] = 'tbin'
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
# bs(Julian, k = 4)      462  4     <2e-16 ***
#   as.factor(Year)        297  7     <2e-16 ***
#   as.factor(Site)        702  4     <2e-16 ***
#   TimeLost                 3  1       0.11    

#BSAI
# bs(Julian)           776  3    < 2e-16 ***
#   TimeLost               0  1       0.98    
# as.factor(Site)       58  1    2.2e-14 ***

Anova(BlockModJ)
#GOA
# bs(Julian, k = 4)     1301  4     <2e-16 ***
#   as.factor(Year)        411  7     <2e-16 ***
#   TimeLost                 7  1     0.0089 ** 
#   as.factor(Site)       2906  4     <2e-16 ***

#BSAI
# bs(Julian, k = 4)    259.8  4    < 2e-16 ***
#   TimeLost               0.1  1       0.78    
# as.factor(Site)       24.2  1    8.7e-07 ***

Anova(BlockModM)
#GOA
# bs(Julian, k = 4)     1417  4     <2e-16 ***
#   as.factor(Year)        706  7     <2e-16 ***
#   TimeLost                10  1     0.0015 ** 
#   as.factor(Site)       2704  4     <2e-16 ***

#BSAI
# bs(Julian, k = 4)      340  4     <2e-16 ***
#   TimeLost                 2  1       0.17    
# as.factor(Site)         79  1     <2e-16 ***

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
#Social Groups
if (region == 'BSAI'){
  GLMF = glm(PreAbsF~bs(Julian)+TimeLost+as.factor(Site),family=binomial,data=SiteHourTableB)
}else{
  GLMF = glm(PreAbsF~bs(Julian)+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
}
VIF(GLMF)
#BSAI
# bs(Julian)      1.780039  3        1.100876
# TimeLost        1.001390  1        1.000695
# as.factor(Site) 1.779350  1        1.333923

#GOA
# GVIF Df GVIF^(1/(2*Df))
# bs(Julian)      1.92e+00  3            1.11
# TimeLost        1.02e+00  1            1.01
# as.factor(Site) 6.68e+06  4            7.13
# as.factor(Year) 9.95e+06  7            3.16

#Mid-Size
if (region == 'BSAI'){
  GLMJ = glm(PreAbsJ~bs(Julian)+TimeLost+as.factor(Site),family=binomial,data=SiteHourTableB)
}else{
  GLMJ = glm(PreAbsJ~bs(Julian)+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
}
VIF(GLMJ)
#BSAI
# bs(Julian)      1.141818  3        1.022350
# TimeLost        1.001430  1        1.000715
# as.factor(Site) 1.142469  1        1.068863

#GOA
# GVIF Df GVIF^(1/(2*Df))
# bs(Julian)      1.67  3            1.09
# TimeLost        1.01  1            1.01
# as.factor(Site) 2.45  4            1.12
# as.factor(Year) 3.78  7            1.10

#Male
if (region == 'BSAI'){
  GLMM = glm(PreAbsM~bs(Julian)+TimeLost+as.factor(Site),family=binomial,data=SiteHourTableB)
}else{
  GLMM = glm(PreAbsM~bs(Julian)+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
}
VIF(GLMM)
#BSAI
# bs(Julian)      1.077366  3        1.012497
# TimeLost        1.000812  1        1.000406
# as.factor(Site) 1.077038  1        1.037804

#GOA
# GVIF Df GVIF^(1/(2*Df))
# bs(Julian)      1.69  3            1.09
# TimeLost        1.01  1            1.01
# as.factor(Site) 2.07  4            1.10
# as.factor(Year) 3.34  7            1.09

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
# QIC0fA   11133.6721960309 10992.4329362201 10995.2210020837
#Julian day as mspline but I'll use a covariance matrix

#BSAI
#QIC            QIC.1           QIC.2
#model0fA            POD0f           POD0fa          POD0fb
# QIC0fA   5202.930974178 4670.3899174931 4662.37409589423
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
  # QIC1fA   11133.6721960309 10969.0658696822 11114.8248943987 11035.5646618752
  #Year as mspline but I'll use factor
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
# QIC2fA   11133.6721960309 11248.1934428678 11166.0957422099 74035.5480750523
#TimeLost as linear

#BSAI
#QIC            QIC.1           QIC.2            QIC.3
#model2fA            POD0f           POD2fa          POD2fb           POD2fc
#QIC2fA   5202.930974178 35832.4334580555 5204.15239629091 5207.78863235595
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
# QIC0jA   69688.050030344 68356.8796353347 68385.0932462919
#Julian day as mspline but I'll use a covariance matrix 

#BSAI
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   15899.9327468631 15510.2345937398 15512.0963512478
#Julian day as mSpline, but I'll use a covariance matrix.

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
# QIC1jA   69688.050030344 67907.6071348806 69692.9425047134 68021.415828091
#Year as factor
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
# QIC2jA   69688.050030344 69775.2476789701 69687.2661175405 69691.0256109908
#TimeLost as linear

#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
# QIC2jA   15899.9327468631 16409.5733666624 15908.1559606344 1168358.62564298
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
# QIC0mA   58482.1073959251 56925.6754971063 56926.3902946418
#Julian day as mspline, but I'll use a covariance matrix.

#BSAI
#QIC            QIC.1           QIC.2
#model0mA            POD0m           POD0ma          POD0mb
# QIC0mA   15553.7948708722 15378.0956045986 15374.5776213622
#Julian day as a covariance matrix.

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
# QIC1mA   58482.1073959251 56197.1446834397 58154.3735691095 56335.9862956417
#Year as a factor
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
# QIC2mA   58482.1073959251 58567.989103007 58492.2768116643 58504.2169092167
#TimeLost as linear.

#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
#QIC2mA   15553.7948708722 15854.0710966028 15558.0403266 15575.2139703477
#TimeLost as linear.

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.

if (region == 'GOA'){
#Social Groups
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ as.factor(Year)+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fc = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Site
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#Without Year
POD3fdd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd","POD3fdd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1],QIC(POD3fdd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#GOA
#QIC            QIC.1           QIC.2            QIC.3            QIC.4            QIC.5
#model3fA            POD0f           POD3fa          POD3fb           POD3fc           POD3fd          POD3fdd
# QIC3fA   11133.6721960309 345318.014426275 74529.0718370052 74765.8132701579 10752.8129276296 10752.5195226224
#Remove Site

#The  full model without Site is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ TimeLost+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fg = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without TimeLost
POD3fgg = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fB = c("POD0f","POD3fe","POD3ff","POD3fg","POD3fgg")
QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1],QIC(POD3fgg)[1])
QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
QICmod3fB
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fB            POD0f           POD3fe           POD3ff           POD3fg          POD3fgg
# QIC3fB   11133.6721960309 10752.8129276297 11006.823877059 11016.8798947422 10726.6559134207
#Remove TimeLost

#The full model without TimeLost is:
POD3fh = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fi = geeglm(PreAbsF ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fj = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3fC = c(QIC(POD0f)[1],QIC(POD3fh)[1],QIC(POD3fi)[1],QIC(POD3fj)[1])
QICmod3fC<-data.frame(rbind(model3fC,QIC3fC))
QICmod3fC
#GOA
#QIC            QIC.1           QIC.2            QIC.3
#model3fC            POD0f           POD3fh          POD3fi           POD3fj
# QIC3fC   11133.6721960309 10726.6559134207 10969.0658696822 10995.2210020837
#Final model

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
  POD3fd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd")
  QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1])
  QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
  QICmod3fA
  #BSAI
  #QIC            QIC.1            QIC.2            QIC.3            QIC.4
  #model3fA            POD0f           POD3fa           POD3fb           POD3fc           POD3fd
  # QIC3fA   5202.930974178 4576.81649417978 5203.55125985791 4576.13482489699 4663.3931245821
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
  # QIC3fB   5202.930974178 4576.13482489699 5202.24732511049 4662.37409589423
  #The full model is the best.
  #AvgDayMat first, then site.
}

#Mid-Size
if (region == 'GOA'){
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ as.factor(Year)+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Site
POD3jdd = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd","POD3jdd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1],QIC(POD3jdd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#GOA
# QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
# model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd          POD3jdd
# QIC3jA   69688.050030344 1396458.43345325 872936.896678686 1651747.77720791 1379203.56370439 66859.0410480842
#Remove Site.

#The  full model without Site is:
POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jf = geeglm(PreAbsJ ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jg = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#Without TimeLost
POD3jee = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jB = c("POD0j","POD3je","POD3jf","POD3jg","POD3jee")
QIC3jB = c(QIC(POD0j)[1],QIC(POD3je)[1],QIC(POD3jf)[1],QIC(POD3jg)[1],QIC(POD3jee)[1])
QICmod3jB<-data.frame(rbind(model3jB,QIC3jB))
QICmod3jB
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3jB            POD0j           POD3je           POD3jf           POD3jg          POD3jee
#QIC3jB   69688.050030344 66859.0410480842 67914.6600993375 68376.8992287065 66857.129060385
#Remove TimeLost

#The full model without TimeLost is:
POD3jh = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3ji = geeglm(PreAbsJ ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jj = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3jC = c(QIC(POD0j)[1],QIC(POD3jh)[1],QIC(POD3ji)[1],QIC(POD3jj)[1])
QICmod3jC<-data.frame(rbind(model3jC,QIC3jC))
QICmod3jC
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model3jC            POD0f           POD3fh           POD3fi           POD3fj
#QIC3jC   69688.050030344 66857.129060385 67907.6071348806 68385.0932462919
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
  #QIC3jA   15899.9327468631 15506.0771134112 15788.7915373773 15500.9728562968 15517.7002658635
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
  #QIC3jB   15899.9327468631 15500.9728562968 15782.6817047258 15512.0963512478
  #Full model is the best.
  #AvgDayMat first, then site.
}

if (region == 'GOA'){
#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ as.factor(Year)+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Site
POD3mdd = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md","POD3mdd")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1],QIC(POD3mdd)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#GOA
#QIC            QIC.1            QIC.2            QIC.3           QIC.4            QIC.5
#model3mA            POD0m           POD3ma           POD3mb           POD3mc          POD3md          POD3mdd
#QIC3mA   58482.1073959251 3601698.03868641 1399599.04881838 53723.8068781844 1003420.38337421 54976.7469709464
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
#QIC3mB  558482.1073959251 53723.8068781844 55357.1068106415 56935.9503789595 53709.992609199
#Remove TimeLost

#The full model without TimeLost is:
POD3mh = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mi = geeglm(PreAbsM ~ as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Site
POD3mj = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mBB = c("POD0m","POD3mh","POD3mi","POD3mj")
QIC3mBB = c(QIC(POD0m)[1],QIC(POD3mh)[1],QIC(POD3mi)[1],QIC(POD3mj)[1])
QICmod3mBB<-data.frame(rbind(model3mBB,QIC3mBB))
QICmod3mBB
#GOA
# QIC            QIC.1            QIC.2            QIC.3
# model3mB            POD0m           POD3mh           POD3mi           POD3mj
# QIC3mB   58482.1073959251 53709.992609199 55332.7701418434 56926.3902946418
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
  #QIC3mA   15553.7948708722 15288.4136686694 15473.2400610951 15378.1832935421 15283.4303511309
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
  #QIC3mB   15553.7948708722 15283.4303511309 15467.319522993 15374.5776213622
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
  PODFinalF = geeglm(PreAbsF ~ as.factor(Year)+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
}

#Mid-Size
if (region == 'BSAI'){
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~ as.factor(Year)+AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}

#Males
if (region == 'BSAI'){
  dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalM = geeglm(PreAbsM ~ as.factor(Site) +AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  }

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinalF)
#GOA
# as.factor(Year)  7 24.2     0.001 ** 
#   AvgDayMatF       2 21.2   2.5e-05 ***

#BSAI
# AvgDayMatF       2 27.33   1.2e-06 ***
#   as.factor(Site)  1  7.38    0.0066 ** 

anova(PODFinalJ)
#GOA
# as.factor(Year)  7 229.1    <2e-16 ***
#   AvgDayMatJ       2  98.9    <2e-16 ***

#BSAI
# AvgDayMatJ       2 31.22   1.7e-07 ***
#   as.factor(Site)  1  1.46      0.23    
# Remove site

anova(PODFinalM)
#GOA
# as.factor(Site)  4 303.2   < 2e-16 ***
#   AvgDayMatM       2  31.8   1.3e-07 ***

#BSAI
# AvgDayMatM       2 21.9   1.8e-05 ***
#   as.factor(Site)  1 14.4   0.00015 ***

#Remove variables from model that were not significant after running anova
#Mid-Size
if (region == 'BSAI'){
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}
anova(PODFinalJ)
#BSAI
#AvgDayMatJ  2 31.215 1.666e-07 ***

filename = paste(saveWorkspace,region,'_Regional_SexSpecific_ModelSummary.txt',sep="")
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
