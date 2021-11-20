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
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/Plots/",site, sep="")
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
#bs(Julian)  1053.54  3     <2e-16 ***
  #TimeLost       0.86  1     0.3537  

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
#bs(Julian)  207.040  3     <2e-16 ***
#TimeLost      0.503  1      0.478 

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
#bs(Julian)   319.38  3     <2e-16 ***
#TimeLost       0.07  1     0.7882  

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
acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1), xlim =c(200,300)) 
ACFvalF = 226
#BD Females - 226
#PT Females - 34
#QN Females - 23
#CB Females - 34

#Juveniles
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1), xlim =c(260,280)) 
ACFvalJ = 272
#BD Juveniles - 272
#PT Juveniles - 110
#QN Juveniles - 30
#CB Juveniles - 584

#Males
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1), xlim =c(100,150))
ACFvalM = 139
#BD Males - 139
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
#1.000118 1.000118 
#VIF(GLMJ)
#Julian TimeLost 
#1.000128 1.000128 
#VIF(GLMM)
#Julian TimeLost 
#1.000119 1.000119 

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
POD0fa = geeglm(PreAbsF ~ bs(Julian, knots=6), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD0fb = geeglm(PreAbsF ~ AvgDayMatF, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model0fA<-c("POD0f", "POD0fa", "POD0fb")
QIC0fA<-c(QIC(POD0f)[1],QIC(POD0fa)[1],QIC(POD0fb)[1])
QICmod0fA<-data.frame(rbind(model0fA,QIC0fA))
QICmod0fA
#BD
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   3344.50244232974 2252.53294008571 2234.0221551908
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
#QIC2fA   3344.50244232974 24362.7152164628 3344.49647341016 24531.7725534563
#TimeLost as linear.

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
#QIC0jA   13025.0450626805 12797.3182476176 12578.8138401894
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
#QIC2jA   13025.0450626805 156315.980205119 13043.8050485186 13074.6655418734
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
#QIC0mA   14925.7545361902 14638.8964876502 14464.0713236508
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
#QIC2jA   13025.0450626805 156315.980205119 13043.8050485186 13074.6655418734
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

#CB (with Year as smooth)
#Females
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fc = geeglm(PreAbsF ~ AvgDayMatF +TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#CB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fA            POD0f           POD3fa           POD3fb           POD3fc           POD3fd
#QIC3fA   2329.98680269474 2236.91997288269 2298.00492186197 2309.00848851999 2229.10419514264
#Remove TimeLost

#The  full model without TimeLost is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
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
#QIC3fB   2329.98680269474 2229.10419514264 2285.18592339013 2303.18308193766
#Full model is the best.

#Juveniles
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#CB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd
#QIC3jA   46282.9243919683 45161.7305048944 45831.3847182987 45466.8386199378 45128.6191087935
#Remove TimeLost

#The  full model without TimeLost is:
POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jf = geeglm(PreAbsJ ~ bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
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
#QIC3jB   46282.9243919683 45128.6191087935 45804.8895179568 45438.1446277905
#Full model is the best.

#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+bs(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#CB
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3mA            POD0m           POD3ma           POD3mb           POD3mc          POD3md
#QIC3mA   
#Remove TimeLost.

#The  full model without TimeLost is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+bs(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ bs(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
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
#QIC3mB   40421.7085304264 38014.3681731602 39648.2399971516 38585.9961341015 38021.988662637
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
#QIC3fA   3344.50244232974 2234.00018665068 3344.49647341016 2234.0221551908
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
#QIC3jA   13025.0450626805 12596.3873961207 13043.8050485186 12578.8138401894
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
#QIC3mA   14925.7545361902 14485.3626400531 14961.5206171963 14464.0713236508
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
#AvgDayMatF  4 20.181 0.0004599 ***
#PT
#AvgDayMatF  4 71.2   1.2e-14 ***
#QN
#AvgDayMatF  4 15.3     0.004 **
#CB
#AvgDayMatF  4  9.7109   0.04559 *  
#as.factor(Year)  7 39799   < 2e-16 ***

#CB
#AvgDayMatF  4  9.71     0.046 *  
  #bs(Year)    3 28.65   2.6e-06 ***

anova(PODFinalj)
#BD
#AvgDayMatJ  4 12.431   0.01442 *
#PT
#AvgDayMatJ  4 11.4     0.022 *
#QN
#AvgDayMatJ  4 25.6   3.8e-05 ***
#CB
#AvgDayMatJ       4 15.782  0.003326 **
  #as.factor(Year)  7 19.843  0.005918 **

#CB
#AvgDayMatJ  4 15.78    0.0033 **
  #bs(Year)    3  5.33    0.1489   

anova(PODFinalm)
#BD
#AvgDayMatM  4 20.11 0.0004751 ***
#PT
#AvgDayMatM  4 1.43      0.84
#QN
#AvgDayMatM  4 10      0.04 *
#CB
#AvgDayMatM       4 61.331 1.523e-12 ***
  #as.factor(Year)  7 39.957 1.283e-06 ***

#CB
#AvgDayMatM       4 61.3   1.5e-12 ***
  #as.factor(Year)  7 40.0   1.3e-06 ***

# STEP 6: Interpretting the summary of the model
#CB ONLY
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)

dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF +bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)
# Call:
#   geeglm(formula = PreAbsF ~ AvgDayMatF + as.factor(Year), family = binomial, 
#          data = SiteHourTableB, id = BlocksF, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err    Wald Pr(>|W|)    
# (Intercept)           -4.918   0.531   85.66   <2e-16 ***
#   AvgDayMatFADBM1        0.893   0.405    4.87    0.027 *  
#   AvgDayMatFADBM2       -0.372   0.395    0.89    0.347    
# as.factor(Year)2012    0.456   0.611    0.56    0.455    
# as.factor(Year)2013   -0.804   0.876    0.84    0.359    
# as.factor(Year)2014  -13.776   1.162  140.54   <2e-16 ***
#   as.factor(Year)2015   -2.323   1.156    4.04    0.044 *  
#   as.factor(Year)2017   -0.183   0.695    0.07    0.792    
# as.factor(Year)2018   -1.157   0.952    1.48    0.224    
# as.factor(Year)2019  -41.847   0.501 6967.36   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     0.76    4.87
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.781    1.29
# Number of clusters:   1344  Maximum cluster size: 35

#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj)

dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj)
# Call:
#   geeglm(formula = PreAbsJ ~ AvgDayMatJ + as.factor(Year), family = binomial, 
#          data = SiteHourTableB, id = BlocksJ, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err  Wald Pr(>|W|)    
# (Intercept)          -0.9974  0.2245 19.73  8.9e-06 ***
#   AvgDayMatJADBM1      -0.8737  0.2288 14.58  0.00013 ***
#   AvgDayMatJADBM2       0.0582  0.1754  0.11  0.73981    
# as.factor(Year)2012  -0.5933  0.3192  3.45  0.06310 .  
# as.factor(Year)2013  -0.5341  0.3537  2.28  0.13102    
# as.factor(Year)2014  -0.5371  0.3631  2.19  0.13904    
# as.factor(Year)2015  -0.6976  0.3505  3.96  0.04654 *  
#   as.factor(Year)2017  -0.0936  0.2696  0.12  0.72831    
# as.factor(Year)2018  -0.5175  0.4944  1.10  0.29520    
# as.factor(Year)2019  -0.2402  0.3069  0.61  0.43385    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.03   0.172
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.918  0.0145
# Number of clusters:   83  Maximum cluster size: 585 

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm)

dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM+bs(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm)
# Call:
#   geeglm(formula = PreAbsM ~ AvgDayMatM + as.factor(Year), family = binomial, 
#          data = SiteHourTableB, id = BlocksM, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err  Wald Pr(>|W|)    
# (Intercept)          -1.2269  0.2530 23.51  1.2e-06 ***
#   AvgDayMatMADBM1       0.5282  0.1390 14.45  0.00014 ***
#   AvgDayMatMADBM2       0.6184  0.1581 15.30  9.2e-05 ***
#   as.factor(Year)2012  -0.2734  0.2834  0.93  0.33478    
# as.factor(Year)2013  -0.5071  0.2550  3.96  0.04672 *  
#   as.factor(Year)2014  -0.9062  0.2805 10.44  0.00124 ** 
#   as.factor(Year)2015  -0.7251  0.3607  4.04  0.04440 *  
#   as.factor(Year)2017  -0.9963  0.3536  7.94  0.00484 ** 
#   as.factor(Year)2018  -0.0583  0.4050  0.02  0.88548    
# as.factor(Year)2019  -0.0460  0.3425  0.02  0.89315    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.02   0.168
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha     0.83  0.0285
# Number of clusters:   170  Maximum cluster size: 277 

#Other sites
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)
#BD
#Call:
# geeglm(formula = PreAbsF ~ AvgDayMatF, family = binomial, data = SiteHourTableB, 
#id = BlocksF, corstr = "ar1")

#Coefficients:
#  Estimate Std.err Wald Pr(>|W|)   
#(Intercept)       -10.15    3.85 6.97   0.0083 **
#  AvgDayMatFADBM1    -1.73    1.31 1.73   0.1885   
#AvgDayMatFADBM2   -10.76    5.47 3.87   0.0491 * 
#Estimate Std.err
#(Intercept)     25.8   26596
#Link = identity 

#Estimated Correlation Parameters:
# Estimate Std.err
#alpha    0.794     186
#Number of clusters:   78  Maximum cluster size: 227

#PT
# Call:
#   geeglm(formula = PreAbsF ~ AvgDayMatF, family = binomial, data = SiteHourTableB, 
#          id = BlocksF, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)       -4.189   0.240 305.54  < 2e-16 ***
#   AvgDayMatFADBM1    1.508   0.434  12.08  0.00051 ***
#   AvgDayMatFADBM2   -0.312   0.211   2.18  0.13961    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.22    4.67
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.794   0.739
# Number of clusters:   427  Maximum cluster size: 35 

#QN
# Call:
#   geeglm(formula = PreAbsF ~ AvgDayMatF, family = binomial, data = SiteHourTableB, 
#          id = BlocksF, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)       -4.153   0.212 383.09   <2e-16 ***
#   AvgDayMatFADBM1   -0.825   0.503   2.69      0.1    
# AvgDayMatFADBM2    0.154   0.184   0.70      0.4    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.14    2.99
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.801   0.511
# Number of clusters:   579  Maximum cluster size: 24 

dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj)
#BD
#Call:
 # geeglm(formula = PreAbsJ ~ AvgDayMatJ, family = binomial, data = SiteHourTableB, 
  #       id = BlocksJ, corstr = "ar1")

#Coefficients:
 # Estimate Std.err   Wald Pr(>|W|)    
#(Intercept)      -1.9958  0.1385 207.74   <2e-16 ***
#  AvgDayMatJADBM1   0.6356  0.2738   5.39     0.02 *  
#  AvgDayMatJADBM2   0.0886  0.2412   0.13     0.71    
#Estimate Std.err
#(Intercept)    0.988   0.254
#Link = identity 

#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha    0.913  0.0251
#Number of clusters:   65  Maximum cluster size: 273 

#PT
# Call:
#   geeglm(formula = PreAbsJ ~ AvgDayMatJ, family = binomial, data = SiteHourTableB, 
#          id = BlocksJ, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)       -3.005   0.145 427.44  < 2e-16 ***
#   AvgDayMatJADBM1    0.352   0.217   2.64      0.1    
# AvgDayMatJADBM2    1.222   0.309  15.62  7.8e-05 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.05   0.711
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.822   0.115
# Number of clusters:   134  Maximum cluster size: 111 

#QN
# Call:
#   geeglm(formula = PreAbsJ ~ AvgDayMatJ, family = binomial, data = SiteHourTableB, 
#          id = BlocksJ, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)       -3.101   0.112 772.94   <2e-16 ***
#   AvgDayMatJADBM1    0.919   0.247  13.80   0.0002 ***
#   AvgDayMatJADBM2    0.584   0.215   7.35   0.0067 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.03   0.537
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha     0.75   0.123
# Number of clusters:   445  Maximum cluster size: 31 

dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm)
#BD
#Call:
#  geeglm(formula = PreAbsM ~ AvgDayMatM, family = binomial, data = SiteHourTableB, 
#         id = BlocksM, corstr = "ar1")

#Coefficients:
#  Estimate Std.err   Wald Pr(>|W|)    
#(Intercept)       -1.733   0.115 228.32   <2e-16 ***
#  AvgDayMatMADBM1   -0.388   0.181   4.61    0.032 *  
#  AvgDayMatMADBM2    0.395   0.232   2.90    0.088 .  
#Estimate Std.err
#(Intercept)    0.987   0.155
#Link = identity 

#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha    0.806  0.0364
#Number of clusters:   65  Maximum cluster size: 273

#PT
# Call:
#   geeglm(formula = PreAbsM ~ AvgDayMatM, family = binomial, data = SiteHourTableB, 
#          id = BlocksM, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)      -3.9240  0.1901 425.90   <2e-16 ***
#   AvgDayMatMADBM1  -0.0253  0.2662   0.01     0.92    
# AvgDayMatMADBM2   0.9024  0.6332   2.03     0.15    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.06     1.4
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.727   0.296
# Number of clusters:   68  Maximum cluster size: 222

#QN
# Call:
#   geeglm(formula = PreAbsM ~ AvgDayMatM, family = binomial, data = SiteHourTableB, 
#          id = BlocksM, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)       -3.016   0.187 261.14   <2e-16 ***
#   AvgDayMatMADBM1    1.088   0.339  10.32   0.0013 ** 
#   AvgDayMatMADBM2    0.983   0.482   4.17   0.0412 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.04    1.02
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.808    0.17
# Number of clusters:   38  Maximum cluster size: 380 

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
#1   298  2053
#0    40 15037

#PT - Female
# observed
# predicted     1     0
# 1   216  1704
# 0    95 12436

#QN - Female
# observed
# predicted     1     0
# 1    95  1249
# 0   141 11775

#CB - Female
# observed
# predicted     1     0
# 1   149 19102
# 0    29 26074
  
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
#1  2150 15278
#0     0     0

#PT - Juvenile
# observed
# predicted     1     0
# 1   829 13622
# 0     0     0

#QN - Juvenile
# observed
# predicted    1    0
# 1  578 9106
# 0   74 3502

#CB
# observed
# predicted     1     0
# 1  9396 35958
# 0     0     0

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
#1  2668 14760
#0     0     0

#PT - Males
# observed
# predicted     1     0
# 1   319 14132
# 0     0     0

#QN - Males
# observed
# predicted     1     0
# 1   752 12508
# 0     0     0

#CB - Males
# observed
# predicted     1     0
# 1  7419 37935
# 0     0     0

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# STEP 7: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
library(boot)
library(pracma)

#Probability of covariate #1: AvgDayBasisMat:
#Female
BootstrapParameters3<-rmvnorm(10000, coef(PODFinalf),summary(PODFinalf)$cov.unscaled)
start=2; finish=3; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinalf)[,start:finish]*coef(PODFinalf)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,366,length=4)))$X[,2:3]
RealFit3<-Basis3%*%coef(PODFinalf)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)

#Base R Plotting
title = paste(saveDir,"/BaseR_Julian Day_SocialGroups - ", site,".png",sep="")
png(title)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main = title , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)
dev.off()

#ggplot
# Calculate kernel density of Jday observations
dJday = stats::density(Variable,na.rm = TRUE,n=5000,from=1,to=366)
dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
colnames(dens) = c("Day", "Density")
dens$Density = dens$Density / 0.15 #max(dens$Density) # normalize kernel density
if (min(cis3a[1,])<0){ # set kernel density at bottom of y axis
  dens$Density = dens$Density - abs(min(cis3a[1,])) 
} else {
  dens$Density = dens$Density + min(cis3a[1,])
}

plotDF = data.frame(PlottingVar3, RealFitCenter3a)
colnames(plotDF) = c("Jday", "Fit")

ggplot(plotDF, aes(Jday, Fit),
) + geom_polygon(data=dens,
                 aes(Day,Density),
                 fill=4,
                 alpha=0.2
) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                stat ="identity"
) + labs(x = "Julian Day",
         y = "Probability",
         title = paste('Julian Day at',site),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_SocialGroups - ", site,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalf),summary(PODFinalf)$cov.unscaled)
start=4; finish=10; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalf)[,start:finish]*coef(PODFinalf)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~as.factor(PlottingVar1), fit=F, family=binomial, knots=list(PlottingVar1=seq(2011,2019,length=3)))$X[,4:10]
RealFit1<-Basis1%*%coef(PODFinalf)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Year as factor
ggtitle = paste(saveDir,"/Probability of Year (SocialGroups) - ", site,".png",sep="")
SiteHourTableB$prf = prf
ggplot(SiteHourTableB, aes(x = Year, y = prf)) +
  geom_boxplot(aes(fill = factor(Year)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedYearCoefs = data.frame(
  c(
    BootstrapParameters1[, 1] - mean(BootstrapParameters1[, 1]),
    BootstrapParameters1[, 2],
    BootstrapParameters1[, 3],
    BootstrapParameters1[, 4],
    BootstrapParameters1[, 5],
    BootstrapParameters1[, 6],
    BootstrapParameters1[, 7],
    BootstrapParameters1[, 8],
    BootstrapParameters1[, 9]
  ),
  as.factor(rep(2011:2019, each = 10000))
)
colnames(AdjustedYearCoefs) = c("Coefficient", "Year")

ggtitle = paste(saveDir,"/Year_SocialGroups - ", site,".png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year at',site))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

#Juvenile
BootstrapParameters3<-rmvnorm(10000, coef(PODFinalj),summary(PODFinalj)$cov.unscaled)
start=2; finish=3; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinalj)[,start:finish]*coef(PODFinalj)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,366,length=4)))$X[,2:3]
RealFit3<-Basis3%*%coef(PODFinalj)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)

#Base R Plotting
title = paste(saveDir,"/BaseR_Julian Day_MidSize - ", site,".png",sep="")
png(title)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main = title , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)
dev.off()

#ggplot
# Calculate kernel density of Jday observations
dJday = stats::density(Variable,na.rm = TRUE,n=5000,from=1,to=366)
dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
colnames(dens) = c("Day", "Density")
dens$Density = dens$Density / 0.15 #max(dens$Density) # normalize kernel density
if (min(cis3a[1,])<0){ # set kernel density at bottom of y axis
  dens$Density = dens$Density - abs(min(cis3a[1,])) 
} else {
  dens$Density = dens$Density + min(cis3a[1,])
}

plotDF = data.frame(PlottingVar3, RealFitCenter3a)
colnames(plotDF) = c("Jday", "Fit")

ggplot(plotDF, aes(Jday, Fit),
) + geom_polygon(data=dens,
                 aes(Day,Density),
                 fill=4,
                 alpha=0.2
) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                stat ="identity"
) + labs(x = "Julian Day",
         y = "Probability",
         title = paste('Julian Day at',site),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_MidSize - ", site,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalj),summary(PODFinalj)$cov.unscaled)
start=4; finish=10; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalj)[,start:finish]*coef(PODFinalj)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~as.factor(PlottingVar1), fit=F, family=binomial, knots=list(PlottingVar1=seq(2011,2019,length=3)))$X[,4:10]
RealFit1<-Basis1%*%coef(PODFinalj)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Year as factor
ggtitle = paste(saveDir,"/Probability of Year (MidSize) - ", site,".png",sep="")
SiteHourTableB$prj = prj
ggplot(SiteHourTableB, aes(x = Year, y = prj)) +
  geom_boxplot(aes(fill = factor(Year)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedYearCoefs = data.frame(
  c(
    BootstrapParameters1[, 1] - mean(BootstrapParameters1[, 1]),
    BootstrapParameters1[, 2],
    BootstrapParameters1[, 3],
    BootstrapParameters1[, 4],
    BootstrapParameters1[, 5],
    BootstrapParameters1[, 6],
    BootstrapParameters1[, 7],
    BootstrapParameters1[, 8],
    BootstrapParameters1[, 9]
  ),
  as.factor(rep(2011:2019, each = 10000))
)
colnames(AdjustedYearCoefs) = c("Coefficient", "Year")

ggtitle = paste(saveDir,"/Year_MidSize - ", site,".png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year at',site))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device


#Male
BootstrapParameters3<-rmvnorm(10000, coef(PODFinalm),summary(PODFinalm)$cov.unscaled)
start=2; finish=3; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinalm)[,start:finish]*coef(PODFinalm)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,366,length=4)))$X[,2:3]
RealFit3<-Basis3%*%coef(PODFinalm)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)

#Base R Plotting
title = paste(saveDir,"/BaseR_Julian Day_Male - ", site,".png",sep="")
png(title)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main = title , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)
dev.off()

#ggplot
# Calculate kernel density of Jday observations
dJday = stats::density(Variable,na.rm = TRUE,n=5000,from=1,to=366)
dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
colnames(dens) = c("Day", "Density")
dens$Density = dens$Density / 0.15 #max(dens$Density) # normalize kernel density
if (min(cis3a[1,])<0){ # set kernel density at bottom of y axis
  dens$Density = dens$Density - abs(min(cis3a[1,])) 
} else {
  dens$Density = dens$Density + min(cis3a[1,])
}

plotDF = data.frame(PlottingVar3, RealFitCenter3a)
colnames(plotDF) = c("Jday", "Fit")

ggplot(plotDF, aes(Jday, Fit),
) + geom_polygon(data=dens,
                 aes(Day,Density),
                 fill=4,
                 alpha=0.2
) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                stat ="identity"
) + labs(x = "Julian Day",
         y = "Probability",
         title = paste('Julian Day at',site),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_Males - ", site,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalj),summary(PODFinalj)$cov.unscaled)
start=4; finish=10; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalj)[,start:finish]*coef(PODFinalj)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~as.factor(PlottingVar1), fit=F, family=binomial, knots=list(PlottingVar1=seq(2011,2019,length=3)))$X[,4:10]
RealFit1<-Basis1%*%coef(PODFinalj)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Year as factor
ggtitle = paste(saveDir,"/Probability of Year (Male) - ", site,".png",sep="")
SiteHourTableB$prj = prm
ggplot(SiteHourTableB, aes(x = Year, y = prm)) +
  geom_boxplot(aes(fill = factor(Year)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedYearCoefs = data.frame(
  c(
    BootstrapParameters1[, 1] - mean(BootstrapParameters1[, 1]),
    BootstrapParameters1[, 2],
    BootstrapParameters1[, 3],
    BootstrapParameters1[, 4],
    BootstrapParameters1[, 5],
    BootstrapParameters1[, 6],
    BootstrapParameters1[, 7],
    BootstrapParameters1[, 8],
    BootstrapParameters1[, 9]
  ),
  as.factor(rep(2011:2019, each = 10000))
)
colnames(AdjustedYearCoefs) = c("Coefficient", "Year")

ggtitle = paste(saveDir,"/Year_Male - ", site,".png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year at',site))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device
