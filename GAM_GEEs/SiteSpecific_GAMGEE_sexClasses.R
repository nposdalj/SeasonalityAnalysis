### Modelling habitat preference of sperm whales using Generalized Additive Models with Generalized Estimating Equations ###
### Script adapted from Pirotta et al. (2011) and Benjamins ###  
### Example from the GofAK + BSAI ###
### 7 Models total:
        #Site specific models: CB, PT, QN, BD (more than 270 d of recording)
        #Region specific models: BSAI + GOA
        #Big model: all 7 sites

## STEP 1: require the libraries needed ##

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
library(car) #for ANOVA

## STEP 1: the data ##
#Hourly data
site = 'CB'
dir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
fileName = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(site,Site))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable = subset(SiteHourTable, Site == site)#subset the table for the site only

#Daily data - for block calculations
fileName2 = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_GAMGEE_ROW_sexClasses.csv")#setting the directory
DayTable = read.csv(fileName2) #no effort days deleted
DayTable = na.omit(DayTable)
DayTable$tbin = as.Date(DayTable$tbin)
SiteDayTable = dplyr::filter(DayTable,grepl(site,Site))
SiteDayTable$Effort_Bin[SiteDayTable$Effort_Bin > 12] = 12
SiteDayTable = subset(SiteDayTable, Site == site)#subset the table for the site only

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin
SiteDayTable$TimeLost = max(SiteDayTable$Effort_Bin) - SiteDayTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region (GOA or BSAI).
# The column Julian represents the day of the year and the column year represents the year of recording.

## Step 3: identify the best blocking structure
#ON MODEL RESIDUALS
#CB ONLY
#Females
BlockModF<-glm(PreAbsF~
               bs(Julian)+
               TimeLost+
               bs(Year)
             ,data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                bs(Julian)+
                TimeLost+
                bs(Year)
              ,data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                bs(Julian)+
                TimeLost+
                bs(Year)
              ,data=SiteHourTable,family=binomial)

#Other sites
#Females
BlockModF<-glm(PreAbsF~
                bs(Julian)+
                TimeLost, data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                bs(Julian)+
                TimeLost, data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                bs(Julian)+
                TimeLost, data=SiteHourTable,family=binomial)

#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockModF)
#BD
#LR Chisq Df Pr(>Chisq)    
#bs(Julian)   926.19  3     <2e-16 ***
  #TimeLost       0.80  1     0.3708   

#PT
#bs(Julian)      410  3     <2e-16 ***
  #TimeLost          2  1       0.18  

#QN
#bs(Julian)     88.5  3     <2e-16 ***
  #TimeLost        3.8  1      0.053 .  

#CB
#bs(Julian)          79.7  3     <2e-16 ***
  #TimeLost             2.4  1       0.12    
#as.factor(Year)    212.1  7     <2e-16 ***

Anova(BlockModJ)
#BD
#bs(Julian)  227.746  3     <2e-16 ***
  #TimeLost      0.004  1     0.9496 

#PT
#bs(Julian)    307.8  3     <2e-16 ***
  #TimeLost        0.7  1       0.42  

#QN
#bs(Julian)      140  3     <2e-16 ***
  #TimeLost          0  1       0.83 

#CB
#bs(Julian)          1647  3     <2e-16 ***
  #TimeLost               6  1      0.019 *  
  #as.factor(Year)      207  7     <2e-16 ***

Anova(BlockModM)
#BD
#bs(Julian)   455.22  3     <2e-16 ***
  #TimeLost       0.04  1     0.8337 

#PT
#bs(Julian)     76.6  3     <2e-16 ***
  #TimeLost        0.1  1        0.8 

#QN
#bs(Julian)      232  3     <2e-16 ***
  #TimeLost          2  1       0.16   

#CB
#bs(Julian)          1336  3    < 2e-16 ***
  #TimeLost              11  1    0.00084 ***
  #as.factor(Year)      852  7    < 2e-16 ***

#Females
acf(residuals(BlockModF), lag.max = 2000, ylim=c(0,0.1))
acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1), xlim =c(0,50)) 
ACFvalF = 34
#BD Females - 207
#PT Females - 34
#QN Females - 23
#CB Females - 34

#Juveniles
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1), xlim =c(500,600)) 
ACFvalJ = 584
#BD Juveniles - ??
#PT Juveniles - 110
#QN Juveniles - 30
#CB Juveniles - 584

#Males
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1), xlim =c(400,500))
ACFvalM = 276
#BD Males - 143
#PT Males - 221
#QN Males - 379
#CB Males - 276 (Got close enough); 422 (got close enough when Year was a smooth instead of a factor)

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

## Step 4: Data exploration and initial analysis ##
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
#CB ONLY
#Females
GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)

#OTHER SITES
#Females
GLMF = glm(PreAbsF~Julian+TimeLost,family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~Julian+TimeLost,family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~Julian+TimeLost,family=binomial,data=SiteHourTableB)

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)

#BD
#VIF(GLMF)
#Julian TimeLost 
#1.000109 1.000109 
#VIF(GLMJ)
#Julian TimeLost 
#1.000081 1.000081 
#VIF(GLMM)
#Julian TimeLost 
#1.000133 1.000133 

#PT/QN
#VIF(GLMF)
#Julian TimeLost 
#1        1 
#VIF(GLMJ)
#Julian TimeLost 
#1        1 
#VIF(GLMM)
#Julian TimeLost 
#1        1 

#CB
#VIF(GLMF)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.35  1            1.16
#TimeLost        1.02  1            1.01
#as.factor(Year) 1.38  7            1.02
#VIF(GLMJ)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.42  1            1.19
#TimeLost        1.01  1            1.00
#as.factor(Year) 1.43  7            1.03
#VIF(GLMM)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.51  1            1.23
#TimeLost        1.01  1            1.00
#as.factor(Year) 1.52  7            1.03

## STEP 4: Model selection - covariate preparation ##
# Construct variance-covariance matrices for cyclic covariates:
#Females
AvgDayBasisF <- gam(PreAbsF~s(Julian, bs ="cc", k=-1), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=6)))$X[,2:5]
AvgDayMatF = as.matrix(AvgDayBasisF)

#Juveniles
AvgDayBasisJ <- gam(PreAbsJ~s(Julian, bs ="cc", k=-1), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=6)))$X[,2:5]
AvgDayMatJ = as.matrix(AvgDayBasisJ)

#Males
AvgDayBasisM <- gam(PreAbsM~s(Julian, bs ="cc", k=-1), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=6)))$X[,2:5]
AvgDayMatM = as.matrix(AvgDayBasisM)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.

#Female
POD0f<-geeglm(PreAbsF ~ 1, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)

#Julian Day
POD0fa = geeglm(PreAbsF ~ bs(Julian, knots=6), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD0fb = geeglm(PreAbsF ~ AvgDayMatF, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model0fA<-c("POD0f", "POD0fa", "POD0fb")
QIC0fA<-c(QIC(POD0f)[1],QIC(POD0fa)[1],QIC(POD0fb)[1])
QICmod0fA<-data.frame(rbind(model0fA,QIC0fA))
QICmod0fA
#BD
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   3124.96557501026 2178.50992322419 2137.61197218677
#Julian day as a covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   3005.53153922808 26592.8902621707 2402.53818766155
#Julian day as a covariance matrix

#QN
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   2372.04767788591 17018.8715640433 2317.84035612086
#Julian day as a covariance matrix

#CB
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   2329.98680269474 12838.7326423582 2303.18308193766
#Julian day as a covariance matrix

#CB ONLY
#Year
POD1fa = geeglm(PreAbsF ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD1fb = geeglm(PreAbsF ~ Year, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD1fc = geeglm(PreAbsF ~ bs(Year), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model1fA<-c("POD0f", "POD1fa", "POD1fb","POD1fc")
QIC1fA<-c(QIC(POD0f)[1],QIC(POD1fa)[1],QIC(POD1fb)[1],QIC(POD1fc)[1])
QICmod1fA<-data.frame(rbind(model1fA,QIC1fA))
QICmod1fA
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model1fA            POD0f           POD1fa           POD1fb           POD1fc
#QIC1fA   2329.98680269474 2155.65331278052 2302.74218858981 2285.18592339013
#Year as factor.

#TimeLost
POD2fa = geeglm(PreAbsF ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fb = geeglm(PreAbsF ~ TimeLost, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fc = geeglm(PreAbsF ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model2fA<-c("POD0f", "POD2fa", "POD2fb","POD2fc")
QIC2fA<-c(QIC(POD0f)[1],QIC(POD2fa)[1],QIC(POD2fb)[1],QIC(POD2fc)[1])
QICmod2fA<-data.frame(rbind(model2fA,QIC2fA))
QICmod2fA
#BD
#QIC           QIC.1            QIC.2            QIC.3
#model2fA            POD0f          POD2fa           POD2fb           POD2fc
#QIC2fA   3124.96557501026 3124.3042051176 3124.95351379776 22193.2887665125
#TimeLost as linear. (Even though as.factor had lowest QIC)

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
#QIC2fA   3005.53153922808 3002.31969255705 22421.2152203148 22421.2153627057
#TimeLost as factor.

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
#QIC2fA   2372.04767788591 2372.86655692344 2370.72248860827 2374.59938113841
#TimeLost as linear.

#CB
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
#QIC2fA   2329.98680269474 12674.5142129543 2338.85032706174 2334.55305261024
#TimeLost as linear.

#Juvenile
POD0j<-geeglm(PreAbsJ ~ 1, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Julian Day
POD0ja = geeglm(PreAbsJ ~ bs(Julian, knots = 6), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD0jb = geeglm(PreAbsJ ~ AvgDayMatJ, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model0jA<-c("POD0j", "POD0ja", "POD0jb")
QIC0jA<-c(QIC(POD0j)[1],QIC(POD0ja)[1],QIC(POD0jb)[1])
QICmod0jA<-data.frame(rbind(model0jA,QIC0jA))
QICmod0jA
#BD
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   13212.5443149178 12996.1539988939 12711.4222500691
#Julian day as covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   6351.80655303413 59755.6613908351 6304.11269776008
#Julian day as covariance matrix

#QN
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   5202.24909252701 50601.1581030481 5047.82935323086
#Julian day as a covariance matrix

#CB                      QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   46282.9243919683 43935.1745235464 45438.1446277905
#Julian day as a covariance matrix

#CB ONLY
#Year
POD1ja = geeglm(PreAbsJ ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jb = geeglm(PreAbsJ ~ Year, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jc = geeglm(PreAbsJ ~ bs(Year), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model1jA<-c("POD0j", "POD1ja", "POD1jb","POD1jc")
QIC1jA<-c(QIC(POD0j)[1],QIC(POD1ja)[1],QIC(POD1jb)[1],QIC(POD1jc)[1])
QICmod1jA<-data.frame(rbind(model1jA,QIC1jA))
QICmod1jA
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model1jA            POD0j           POD1ja           POD1jb           POD1jc
#QIC1jA   46282.9243919683 45352.2087787943 46188.4957092451 45804.8895179568
#Year as factor.

#TimeLost
POD2ja = geeglm(PreAbsJ ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jb = geeglm(PreAbsJ ~ TimeLost, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jc = geeglm(PreAbsJ ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model2jA<-c("POD0j", "POD2ja", "POD2jb","POD2jc")
QIC2jA<-c(QIC(POD0j)[1],QIC(POD2ja)[1],QIC(POD2jb)[1],QIC(POD2jc)[1])
QICmod2jA<-data.frame(rbind(model2jA,QIC2jA))
QICmod2jA
#BD
#QIC           QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
#QIC2jA   13212.5443149178 158647.006516794 13223.1141241812 13237.9776447194
#TimeLost as linear.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
#QIC2jA   6351.80655303413 6384.51733176367 6358.76854561236 6358.05905313353
#TimeLost as linear.

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
#QIC2jA   5202.24909252701 47040.5691370706 5205.25607128381 5200.86795837766
#TimeLost as linear. (Even though smooth had the lowest QIC)

#CB
#QIC           QIC.1            QIC.2            QIC.3
#model2jA            POD0j          POD2ja           POD2jb           POD2jc
#QIC2jA   46282.9243919683 46427.700926742 46297.0392742771 46314.1783902309
#TimeLost as linear.

#Male
POD0m<-geeglm(PreAbsM ~ 1, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

#Julian Day
POD0ma = geeglm(PreAbsM ~ bs(Julian, knots=6), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD0mb = geeglm(PreAbsM ~ AvgDayMatJ, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model0mA<-c("POD0m", "POD0ma", "POD0mb")
QIC0mA<-c(QIC(POD0m)[1],QIC(POD0ma)[1],QIC(POD0mb)[1])
QICmod0mA<-data.frame(rbind(model0mA,QIC0mA))
QICmod0mA
#BD
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
#QIC0mA   14995.1532685815 14565.7063029356 14512.7636741952
#Julian Day as covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
#QIC0mA   3070.61673630417 3009.99059215269 3070.06982622696
#Julian Day as covariance matrix

#QN
#QIC            QIC.1           QIC.2
#model0mA            POD0m           POD0ma          POD0mb
#QIC0mA   5784.45236600162 57590.2357444109 5688.3973482198
#Julian Day as covariance matrix

#CB
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
#QIC0mA   40421.7085304264 38926.3807141555 38586.8059227603
#Julian Day as covariance matrix

#CB ONLY
#Year
POD1ma = geeglm(PreAbsM ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mb = geeglm(PreAbsM ~ Year, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mc = geeglm(PreAbsM ~ bs(Year), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model1mA<-c("POD0m", "POD1ma", "POD1mb","POD1mc")
QIC1mA<-c(QIC(POD0m)[1],QIC(POD1ma)[1],QIC(POD1mb)[1],QIC(POD1mc)[1])
QICmod1mA<-data.frame(rbind(model1mA,QIC1mA))
QICmod1mA
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model1mA            POD0m           POD1ma           POD1mb           POD1mc
#QIC1mA   40421.7085304264 39564.2836988386 40398.2543820786 39636.0229647741
#Year as factor

#TimeLost
POD2ma = geeglm(PreAbsM ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mb = geeglm(PreAbsM ~ TimeLost, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mc = geeglm(PreAbsM ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model2mA<-c("POD0m", "POD2ma", "POD2mb","POD2mc")
QIC2mA<-c(QIC(POD0m)[1],QIC(POD2ma)[1],QIC(POD2mb)[1],QIC(POD2mc)[1])
QICmod2mA<-data.frame(rbind(model2mA,QIC2jA))
QICmod2mA
#BD
#QIC           QIC.1            QIC.2            QIC.3
#model2mA            POD0m          POD2ma           POD2mb           POD2mc
#QIC2jA   13212.5443149178 3124.3042051176 3124.95351379776 22193.2887665125
#TimeLost as factor.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
#QIC2jA   6351.80655303413 6384.51733176367 6358.76854561236 6358.05905313353
#TimeLost as linear.

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
#QIC2jA   5202.24909252701 47040.5691370706 5205.25607128381 5200.86795837766
#TimeLost as linear (Even though smooth had a lower QIC)

#CB
#QIC           QIC.1            QIC.2            QIC.3
#model2mA            POD0m          POD2ma           POD2mb           POD2mc
#QIC2jA   46282.9243919683 46427.700926742 46297.0392742771 46314.1783902309
#TimeLost as linear.

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#CB (with Year as factor)
#Females
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fc = geeglm(PreAbsF ~ AvgDayMatF +TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#CB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fA            POD0f           POD3fa           POD3fb           POD3fc           POD3fd
#QIC3fA   2329.98680269474 2092.14542796175 2162.16652690524 2309.00848851999 2083.11460422949
#Remove TimeLost

#The  full model without TimeLost is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fg = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fB = c("POD0f","POD3fe","POD3ff","POD3fg")
QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1])
QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
QICmod3fB
PODFinalf = POD3fe
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model3fB            POD0f           POD3fe           POD3ff           POD3fg
#QIC3fB   2329.98680269474 2083.11460422949 2155.65331278052 2303.18308193766
#Full model is the best.

#Juveniles
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#CB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd
#QIC3jA   46282.9243919683 44804.4291555885 45372.5285549423 45466.8386199378 44777.747938863
#Remove TimeLost

#The  full model without TimeLost is:
POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jf = geeglm(PreAbsJ ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jg = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jB = c("POD0j","POD3je","POD3jf","POD3jg")
QIC3jB = c(QIC(POD0j)[1],QIC(POD3je)[1],QIC(POD3jf)[1],QIC(POD3jg)[1])
QICmod3jB<-data.frame(rbind(model3jB,QIC3jB))
QICmod3jB
PODFinalj = POD3je
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model3jB            POD0j           POD3je           POD3jf           POD3jg
#QIC3jB   46282.9243919683 44777.747938863 45352.2087787943 45438.1446277905
#Full model is the best.

#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#CB
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3mA            POD0m           POD3ma           POD3mb           POD3mc          POD3md
#QIC3mA   40421.7085304264 37766.669659805 39578.0703832276 38585.9961341015 37778.0580942321
#Remove TimeLost.

#The  full model without TimeLost is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mg = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3me","POD3mf","POD3mg")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
QICmod3mB
PODFinalm = POD3me
#CB
#QIC           QIC.1            QIC.2            QIC.3
#model3mB            POD0m          POD3me           POD3mf           POD3mg
#QIC3mB   40421.7085304264 37778.0580942321 39564.2836988386 38586.8059227603
#Full model is the best.

#Other sites (without year)
#Females - TimeLost as linear (BD, PT)
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without TimeLost
POD3fc = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("PODf0","POD3fa","POD3fb","POD3fc")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
PODFinalf = POD3fc
#BD
#QIC            QIC.1           QIC.2            QIC.3
#model3fA            PODf0           POD3fa           POD3fb           POD3fc
#QIC3fA   3124.96557501026 2137.60559920481 3124.95351379776 2137.61197218677
#Remove TimeLost. POD3fc is final model.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model3fA            PODf0           POD3fa           POD3fb           POD3fc
#QIC3fA   3005.53153922808 47078.6219076777 22421.2152203148 2402.53818766155
#Remove TimeLost. POD3fc is final model.

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model3fA            PODf0           POD3fa           POD3fb           POD3fc
#QIC3fA   2372.04767788591 2316.58036715227 2370.72248860827 2317.84035612086
#Remove TimeLost. POD3fc is final model.


#Juveniles - TimeLost as linear (BD, PT)
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without TimeLost
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3fb","POD3fc")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
PODFinalj = POD3jc
#BD
#QIC            QIC.1           QIC.2            QIC.3
#model3jA            POD0j           POD3ja           POD3fb           POD3fc
#QIC3jA   13212.5443149178 12723.1480150494 13223.1141241812 12711.4222500691
#Remove TimeLost. POD3jc

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model3jA            POD0j           POD3ja           POD3fb           POD3fc
#QIC3jA   6351.80655303413 6311.90947090316 6358.76854561236 6304.11269776008
#Remove TimeLost. POD3jc

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model3jA            POD0j           POD3ja           POD3fb           POD3fc
#QIC3jA   5202.24909252701 5050.14106517662 5205.25607128381 5047.82935323086
#Remove TimeLost. POD3jc

#Males - TimeLost as factor (BD, )
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(TimeLost),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ as.factor(TimeLost),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without TimeLost
POD3mc = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
PODFinalm = POD3mc
#BD
#QIC            QIC.1            QIC.2            QIC.3
#model3mA            POD0m           POD3ma           POD3mb           POD3mc
#QIC3mA   14995.1532685815 14528.2397167441 15027.6497125557 14512.7636741952
#Remove TimeLost. POD3mc

#Males - TimeLost as linear (PT, QN)
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without TimeLost
POD3mc = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
PODFinalm = POD3mc
#PT
#QIC            QIC.1            QIC.2            QIC.3
#model3mA            POD0m           POD3ma           POD3mb           POD3mc
#QIC3mA   3070.61673630417 3075.51183535275 3076.21481422971 3070.06982622696
#Remove TimeLost. POD3mc

#QN
#QIC            QIC.1            QIC.2           QIC.3
#model3mA            POD0m           POD3ma           POD3mb          POD3mc
#QIC3mA   5784.45236600162 5709.72503471671 5816.46401715114 5688.3973482198
#Remove TimeLost. POD3mc

# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

anova(PODFinalf)
#BD
#AvgDayMatF  4 33.200 1.087e-06 ***
#PT
#AvgDayMatF  4 71.2   1.2e-14 ***
#QN
#AvgDayMatF  4 15.3     0.004 **
#CB
#AvgDayMatF  4  9.7109   0.04559 *  
#as.factor(Year)  7 39799   < 2e-16 ***

anova(PODFinalj)
#BD
#AvgDayMatJ  4 9.177   0.05682 .
#PT
#AvgDayMatJ  4 11.4     0.022 *
#QN
#AvgDayMatJ  4 25.6   3.8e-05 ***
#CB
#AvgDayMatJ       4 15.782  0.003326 **
  #as.factor(Year)  7 19.843  0.005918 **

anova(PODFinalm)
#BD
#AvgDayMatM  4 26.788 2.194e-05 ***
#PT
#AvgDayMatM  4 1.43      0.84
#QN
#AvgDayMatM  4 10      0.04 *
#CB
#AvgDayMatM       4 61.331 1.523e-12 ***
  #as.factor(Year)  7 39.957 1.283e-06 ***

# STEP 6: Interpretting the summary of the model
#CB ONLY
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)
#(Intercept)          -5.0055   0.4208  141.480   <2e-16 ***
  #AvgDayMatFADBM1      -0.9295   0.9992    0.865   0.3523    
#AvgDayMatFADBM2       2.5181   1.1546    4.756   0.0292 *  
  #AvgDayMatFADBM3       1.4632   0.7352    3.961   0.0466 *  
  #AvgDayMatFADBM4       0.8436   0.4581    3.392   0.0655 .  
#as.factor(Year)2012   0.6048   0.5267    1.319   0.2508    
#as.factor(Year)2013  -0.5426   0.7941    0.467   0.4945    
#as.factor(Year)2014 -16.7198   1.4728  128.873   <2e-16 ***
  #as.factor(Year)2015  -2.8886   1.5412    3.513   0.0609 .  
#as.factor(Year)2017  -0.2074   0.6262    0.110   0.7405    
#as.factor(Year)2018  -1.8173   0.8546    4.522   0.0335 *  
  #as.factor(Year)2019 -42.8734   0.4591 8719.492   <2e-16 ***
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)   0.6775   3.866
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha   0.7566    1.25
#Number of clusters:   1344  Maximum cluster size: 35 

#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj)
#(Intercept)           -0.734   0.260  7.98  0.00472 ** 
  #AvgDayMatJADBM1       -0.808   0.297  7.40  0.00653 ** 
  #AvgDayMatJADBM2        0.306   0.395  0.60  0.43939    
#AvgDayMatJADBM3        0.282   0.217  1.70  0.19218    
#AvgDayMatJADBM4       -0.403   0.218  3.41  0.06492 .  
#as.factor(Year)2012   -0.678   0.352  3.71  0.05411 .  
#as.factor(Year)2013   -0.544   0.407  1.79  0.18052    
#as.factor(Year)2014   -0.796   0.430  3.43  0.06413 .  
#as.factor(Year)2015   -1.206   0.322 14.03  0.00018 ***
  #as.factor(Year)2017   -0.287   0.343  0.70  0.40336    
#as.factor(Year)2018   -1.063   0.367  8.41  0.00373 ** 
  #as.factor(Year)2019   -0.712   0.313  5.17  0.02292 *  
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)     0.99   0.108
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.931  0.0108
#Number of clusters:   83  Maximum cluster size: 585 

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm)
#(Intercept)           -1.414   0.141 100.09  < 2e-16 ***
  #AvgDayMatMADBM1       -2.276   0.414  30.26  3.8e-08 ***
  #AvgDayMatMADBM2       -0.717   0.269   7.10  0.00772 ** 
  #AvgDayMatMADBM3        0.660   0.178  13.72  0.00021 ***
  #AvgDayMatMADBM4        0.640   0.171  14.03  0.00018 ***
  #as.factor(Year)2012   -0.242   0.199   1.48  0.22327    
#as.factor(Year)2013   -0.408   0.210   3.77  0.05232 .  
#as.factor(Year)2014   -0.814   0.187  18.94  1.4e-05 ***
  #as.factor(Year)2015   -0.470   0.245   3.68  0.05496 .  
#as.factor(Year)2017   -1.007   0.285  12.52  0.00040 ***
  #as.factor(Year)2018    0.128   0.272   0.22  0.63875    
#as.factor(Year)2019    0.120   0.328   0.13  0.71445    
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)     1.07   0.595
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.806  0.0955
#Number of clusters:   170  Maximum cluster size: 277 

#Other sites
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)
#BD
#Estimate Std.err   Wald Pr(>|W|)    
#(Intercept)      -6.8317  0.9493 51.792 6.17e-13 ***
  #AvgDayMatFADBM1   3.8403  0.7203 28.421 9.76e-08 ***
  #AvgDayMatFADBM2  -7.6374  6.5473  1.361    0.243    
#AvgDayMatFADBM3   0.9036  0.7042  1.647    0.199    
#AvgDayMatFADBM4  -5.3541  7.8051  0.471    0.493    
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)    482.3 2360897
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha   0.8862   531.9
#Number of clusters:   85  Maximum cluster size: 208 

#PT
#(Intercept)       -5.401   0.463 136.30  < 2e-16 ***
  #AvgDayMatFADBM1   -6.280   1.960  10.27   0.0014 ** 
  #AvgDayMatFADBM2   -0.948   0.615   2.37   0.1236    
#AvgDayMatFADBM3    3.396   0.434  61.30  4.9e-15 ***
  #AvgDayMatFADBM4    0.483   0.653   0.55   0.4594    
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)      2.9     649
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha     0.25    49.6
#Number of clusters:   427  Maximum cluster size: 35 

#QN
#(Intercept)       -4.201   0.180 541.74   <2e-16 ***
  #AvgDayMatFADBM1   -0.667   0.363   3.37    0.066 .  
#AvgDayMatFADBM2   -1.163   0.496   5.51    0.019 *  
  #AvgDayMatFADBM3    0.611   0.349   3.06    0.080 .  
#AvgDayMatFADBM4    0.916   0.388   5.56    0.018 *  
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)     1.03    1.35
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.777   0.288
#Number of clusters:   579  Maximum cluster size: 24 

dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj)
#BD
#Estimate Std.err   Wald Pr(>|W|)    
#(Intercept)       -2.005   0.147 185.58   <2e-16 ***
  #AvgDayMatJADBM1    0.500   0.461   1.18    0.278    
#AvgDayMatJADBM2   -0.519   0.437   1.41    0.235    
#AvgDayMatJADBM3    0.568   0.388   2.15    0.143    
#AvgDayMatJADBM4    1.119   0.443   6.38    0.012 *  
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)        1   0.347
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.922  0.0317
#Number of clusters:   53  Maximum cluster size: 332 

#PT
#(Intercept)       -2.841   0.157 326.46   <2e-16 ***
  #AvgDayMatJADBM1   -0.289   0.304   0.91    0.341    
#AvgDayMatJADBM2    0.273   0.329   0.69    0.407    
#AvgDayMatJADBM3    0.550   0.252   4.76    0.029 *  
  #AvgDayMatJADBM4   -0.378   0.295   1.64    0.200    
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)     1.01   0.467
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.863  0.0735
#Number of clusters:   134  Maximum cluster size: 111 

#QN
#(Intercept)       -3.194   0.120 710.99  < 2e-16 ***
  #AvgDayMatJADBM1   -3.179   0.788  16.26  5.5e-05 ***
  #AvgDayMatJADBM2   -0.205   0.395   0.27    0.604    
#AvgDayMatJADBM3    0.714   0.299   5.69    0.017 *  
  #AvgDayMatJADBM4    0.487   0.241   4.08    0.043 *  
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)     1.05    2.54
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.738   0.558
#Number of clusters:   445  Maximum cluster size: 31 

dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm)
#BD
#Estimate Std.err   Wald Pr(>|W|)    
#(Intercept)      -1.7937  0.0859 436.22  < 2e-16 ***
  #AvgDayMatMADBM1  -0.8917  0.3423   6.79  0.00919 ** 
  #AvgDayMatMADBM2  -1.2530  0.3493  12.87  0.00033 ***
  #AvgDayMatMADBM3  -0.3113  0.2173   2.05  0.15193    
#AvgDayMatMADBM4  -0.2153  0.2279   0.89  0.34483    
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)    0.987   0.132
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.794  0.0293
#Number of clusters:   123  Maximum cluster size: 144 

#PT
#(Intercept)      -3.8319  0.2821 184.49   <2e-16 ***
  #AvgDayMatMADBM1  -0.2072  0.3721   0.31     0.58    
#AvgDayMatMADBM2  -0.2463  0.4878   0.25     0.61    
#AvgDayMatMADBM3  -0.0227  0.5326   0.00     0.97    
#AvgDayMatMADBM4  -0.3580  0.3916   0.84     0.36    
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)     1.03    1.65
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.763   0.361
#Number of clusters:   68  Maximum cluster size: 222 

#QN
#(Intercept)      -2.9390  0.2272 167.32   <2e-16 ***
  #AvgDayMatMADBM1  -2.3395  0.9040   6.70   0.0097 ** 
  #AvgDayMatMADBM2  -0.6073  0.4962   1.50   0.2210    
#AvgDayMatMADBM3   0.2475  0.5331   0.22   0.6424    
#AvgDayMatMADBM4  -0.0146  0.5659   0.00   0.9794    
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)    0.994   0.847
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.826   0.142
#Number of clusters:   38  Maximum cluster size: 380 

# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

# STEP 7: Construction of the ROC curve  
#Females
prf <- predict(PODFinalf, type="response")  
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
DATA$Predicted<-predict(PODFinalf,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:
#BD - Female
#observed
#predicted     1     0
#1   266  1701
#0    45 15416

#PT - Female
#observed
#predicted     1     0
#1   251  2077
#0    60 12063

#QN - Female
#observed
#predicted    1    0
#1  212 6999
#0   24 6025
  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

#Juveniles
prj <- predict(PODFinalj, type="response")  
pred <- prediction(prf,SiteHourTableB$PreAbsJ) 
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
DATA$Predicted<-predict(PODFinalj,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:
#BD - Juvenile
#observed
#predicted     1     0
#1  2198 15230
#0     0     0

#PT - Juvenile
#observed
#predicted     1     0
#1   829 13622
#0     0     0

#QN - Juvenile
#observed
#predicted     1     0
#1   649 11291
#0     3  1317

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

#Males
prm <- predict(PODFinalj, type="response")  
pred <- prediction(prf,SiteHourTableB$PreAbsM) 
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
DATA$Predicted<-predict(PODFinalm,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:
#BD - Males
#observed
#predicted     1     0
#1  2689 14739
#0     0     0

#PT - Males
#observed
#predicted     1     0
#1   319 14132
#0     0     0

#QN - Males
#observed
#predicted     1     0
#1   741 11607
#0    11   901

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# STEP 7: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
library(boot)
library(pracma)

#Probability of covariate #1: AvgDayBasisMat:
#Female
BootstrapParameters3<-rmvnorm(10000, coef(PODFinalf),summary(PODFinalf)$cov.unscaled)
start=2; finish=5; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinalf)[,start:finish]*coef(PODFinalf)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc",k=6), fit=F, family=binomial, knots=list(PlottingVar3=seq(1,366,length=6)))$X[,2:5]
RealFit3<-Basis3%*%coef(PODFinalf)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main ="Julian Day" , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)

#Juvenile
BootstrapParameters3<-rmvnorm(10000, coef(PODFinalj),summary(PODFinalj)$cov.unscaled)
start=2; finish=5; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinalj)[,start:finish]*coef(PODFinalj)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc",k=6), fit=F, family=binomial, knots=list(PlottingVar3=seq(1,366,length=6)))$X[,2:5]
RealFit3<-Basis3%*%coef(PODFinalj)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main ="Julian Day" , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)

#Male
BootstrapParameters3<-rmvnorm(10000, coef(PODFinalm),summary(PODFinalm)$cov.unscaled)
start=2; finish=5; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinalm)[,start:finish]*coef(PODFinalm)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc",k=6), fit=F, family=binomial, knots=list(PlottingVar3=seq(1,366,length=6)))$X[,2:5]
RealFit3<-Basis3%*%coef(PODFinalm)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main ="Julian Day" , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)

#as.factor(Year) (CB ONLY)
SiteHourTableB$prf = prf
ggplot(SiteHourTableB, aes(x = Year, y = prf)) +
  geom_boxplot(aes(fill = factor(Year)), alpha = .2)

SiteHourTableB$prj = prj
ggplot(SiteHourTableB, aes(x = Year, y = prj)) +
  geom_boxplot(aes(fill = factor(Year)), alpha = .2)

SiteHourTableB$prm = prm
ggplot(SiteHourTableB, aes(x = Year, y = prm)) +
  geom_boxplot(aes(fill = factor(Year)), alpha = .2)

### IGNORE BELOW THIS LINE

#Probability of covariate #2: TimeLost:
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalf),summary(PODFinalf)$cov.unscaled)
val = 6; Variable=SiteHourTableB$TimeLost; xlabel="Time Lost"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalf)[,val]*coef(PODFinalf)[val]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
#Basis1<-gam(rbinom(5000,1,0.5)~PlottingVar1, fit=F, family=binomial)
#Basis1<-glm(rbinom(5000,1,0.5)~PlottingVar1, family=binomial)
Basis1<-gam(rbinom(5000,1,0.5)~PlottingVar1, fit=F, family=binomial)$X[,2]
RealFit1<-Basis1*coef(PODFinalf)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)
plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(0,10), main ="Time Lost" , cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(PlottingVar1)
