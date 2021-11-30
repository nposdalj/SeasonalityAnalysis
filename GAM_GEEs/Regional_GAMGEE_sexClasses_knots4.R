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
region = 'GOA'
dir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/Plots/",region, sep="")
fileName = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="")
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(region,Region))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12

#Daily data - for block calculations
fileName2 = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_GAMGEE_ROW_sexClasses.csv")#setting the directory
DayTable = read.csv(fileName2) #no effort days deleted
DayTable = na.omit(DayTable)
DayTable$tbin = as.Date(DayTable$tbin)
SiteDayTable = dplyr::filter(DayTable,grepl(region,Region))
SiteDayTable$Effort_Bin[SiteDayTable$Effort_Bin > 12] = 12

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
#Region Specific
#BSAI
#Females
BlockModF<-glm(PreAbsF~
                 bs(Julian)+
                 TimeLost+
                 as.factor(Site)+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian)+
                 TimeLost+
                 as.factor(Site)+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian)+
                 TimeLost+
                 as.factor(Site)+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

#GOA
#Females
BlockModF<-glm(PreAbsF~
                 bs(Julian)+
                 TimeLost+
                 as.factor(Site)+
                 bs(Year)
               ,data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian)+
                 TimeLost+
                 as.factor(Site)+
                 bs(Year)
               ,data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian)+
                 TimeLost+
                 as.factor(Site)+
                 bs(Year)
               ,data=SiteHourTable,family=binomial)


#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockModF)
#GOA
#bs(Julian)        378.93  3    < 2e-16 ***
#TimeLost            4.43  1    0.03541 *  
#as.factor(Site)   433.66  4    < 2e-16 ***
#bs(Year)          112.41  3    < 2e-16 ***

#BSAI
#bs(Julian)           793  3     <2e-16 ***
#TimeLost               1  1       0.38    
#as.factor(Site)        0  1       1.00    
#as.factor(Year)        4  2       0.11  

Anova(BlockModJ)
#GOA
#bs(Julian)       1636.21  3    < 2e-16 ***
#TimeLost            3.99  1    0.04569 *  
#as.factor(Site)  3023.25  4    < 2e-16 ***
#bs(Year)          321.33  3    < 2e-16 ***

#BSAI
#bs(Julian)         165.7  3    < 2e-16 ***
#TimeLost             0.1  1       0.78    
#as.factor(Site)     55.0  1    1.2e-13 ***
#as.factor(Year)    179.3  2    < 2e-16 ***

Anova(BlockModM)
#GOA
#bs(Julian)       1159.37  3  < 2.2e-16 ***
#TimeLost            7.38  1   0.006586 ** 
#as.factor(Site)  2574.30  4  < 2.2e-16 ***
#bs(Year)          585.36  3  < 2.2e-16 ***

#BSAI
#bs(Julian)           424  3     <2e-16 ***
#TimeLost               1  1     0.3597    
#as.factor(Site)       10  1     0.0017 ** 
#as.factor(Year)      151  2     <2e-16 ***

#Females
acf(residuals(BlockModF), lag.max = 2000, ylim=c(0,0.1))
acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1), xlim =c(0,50)) 
ACFvalF = 31
#GOA - 31
#BSAI - 217

#Juveniles
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1), xlim =c(0,70))
ACFvalJ = 53
#GOA - 53
#BSAI - 274

#Males
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1), xlim =c(0,50))
ACFvalM = 31
#GOA - 31
#BSAI - 131

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
#Region Specific
#GOA
#Females
GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)

#BSAI
#Females
GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)

#GOA
#Julian          1.115622  1        1.056230
#TimeLost        1.011201  1        1.005585
#as.factor(Site) 1.291714  4        1.032514
#Year            1.412620  1        1.188537
#VIF(GLMJ)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.093228  1        1.045576
#TimeLost        1.006406  1        1.003198
#as.factor(Site) 1.277982  4        1.031135
#Year            1.384001  1        1.176436
#VIF(GLMM)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.069528  1        1.034180
#TimeLost        1.005952  1        1.002971
#as.factor(Site) 1.242515  4        1.027514
#Year            1.320004  1        1.148914

#BSAI
#Julian          1.00  1            1.00
#TimeLost        1.00  1            1.00
#as.factor(Site) 1.54  1            1.24
#as.factor(Year) 1.54  2            1.11
#VIF(GLMJ)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.55  1            1.24
#TimeLost        1.00  1            1.00
#as.factor(Site) 2.02  1            1.42
#as.factor(Year) 2.74  2            1.29
#VIF(GLMM)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.42  1            1.19
#TimeLost        1.00  1            1.00
#as.factor(Site) 1.40  1            1.19
#as.factor(Year) 1.84  2            1.17

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
#GOA
#QIC           QIC.1            QIC.2
#model0fA            POD0f          POD0fa           POD0fb
#QIC0fA   8456.24692284754 8151.4973958272 8032.27430646315
#Julian day as a covariance matrix

#BSAI
#QIC            QIC.1           QIC.2
#model0fA            POD0f           POD0fa          POD0fb
#QIC0fA   3184.88930066242 2160.28051775351 2121.0881909267
#Julian day as a covariance matrix

#GOA ONLY
#Year
POD1fa = geeglm(PreAbsF ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD1fb = geeglm(PreAbsF ~ Year, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD1fc = geeglm(PreAbsF ~ bs(Year), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model1fA<-c("POD0f", "POD1fa", "POD1fb","POD1fc")
QIC1fA<-c(QIC(POD0f)[1],QIC(POD1fa)[1],QIC(POD1fb)[1],QIC(POD1fc)[1])
QICmod1fA<-data.frame(rbind(model1fA,QIC1fA))
QICmod1fA
#GOA
#QIC           QIC.1            QIC.2            QIC.3
#model1fA            POD0f          POD1fa           POD1fb           POD1fc
#QIC1fA   8456.24692284754 8154.7264917688 8337.53585000172 8325.14077719338
#Year as factor

#TimeLost
POD2fa = geeglm(PreAbsF ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fb = geeglm(PreAbsF ~ TimeLost, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fc = geeglm(PreAbsF ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model2fA<-c("POD0f", "POD2fa", "POD2fb","POD2fc")
QIC2fA<-c(QIC(POD0f)[1],QIC(POD2fa)[1],QIC(POD2fb)[1],QIC(POD2fc)[1])
QICmod2fA<-data.frame(rbind(model2fA,QIC2fA))
QICmod2fA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
#QIC2fA   8456.24692284754 8521.81291766263 8464.69313539556 8502.37071568201
#TimeLost as linear

#BSAI
#QIC            QIC.1           QIC.2            QIC.3
#model2fA            POD0f           POD2fa          POD2fb           POD2fc
#QIC2fA   3184.88930066242 22421.9335840961 3184.8825198751 22397.4094063534
#TimeLost as linear

#Juvenile
POD0j<-geeglm(PreAbsJ ~ 1, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Julian Day
POD0ja = geeglm(PreAbsJ ~ bs(Julian, knots = 6), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD0jb = geeglm(PreAbsJ ~ AvgDayMatJ, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model0jA<-c("POD0j", "POD0ja", "POD0jb")
QIC0jA<-c(QIC(POD0j)[1],QIC(POD0ja)[1],QIC(POD0jb)[1])
QICmod0jA<-data.frame(rbind(model0jA,QIC0jA))
QICmod0jA
#GOA
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   66581.7159280241 64606.6957774796 65569.0171679479
#Julian day as a covariance matrix.

#BSAI
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   15081.9977983076 14827.5166886239 14573.7082927332
#Julian day as a covariance matrix.

#GOA ONLY
#Year
POD1ja = geeglm(PreAbsJ ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jb = geeglm(PreAbsJ ~ Year, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jc = geeglm(PreAbsJ ~ bs(Year), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model1jA<-c("POD0j", "POD1ja", "POD1jb","POD1jc")
QIC1jA<-c(QIC(POD0j)[1],QIC(POD1ja)[1],QIC(POD1jb)[1],QIC(POD1jc)[1])
QICmod1jA<-data.frame(rbind(model1jA,QIC1jA))
QICmod1jA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model1jA            POD0j           POD1ja           POD1jb           POD1jc
#QIC1jA   66581.7159280241 64914.0002469341 66586.0907020832 65172.0081373343
#Year as factor.

#TimeLost
POD2ja = geeglm(PreAbsJ ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jb = geeglm(PreAbsJ ~ TimeLost, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jc = geeglm(PreAbsJ ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model2jA<-c("POD0j", "POD2ja", "POD2jb","POD2jc")
QIC2jA<-c(QIC(POD0j)[1],QIC(POD2ja)[1],QIC(POD2jb)[1],QIC(POD2jc)[1])
QICmod2jA<-data.frame(rbind(model2jA,QIC2jA))
QICmod2jA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
#QIC2jA   66581.7159280241 66748.0688240194 66602.0567455407 66621.7408878972
#TimeLost as linear

#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
#QIC2jA   15081.9977983076 185051.245308155 15088.5464462827 15104.0831880007
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
#GOA
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
#QIC0mA   57273.6041843083 55723.6653887203 55666.3718955172
#Julian day as a covariance matrix.

#BSAI
#QIC            QIC.1           QIC.2
#model0mA            POD0m           POD0ma          POD0mb
#QIC0mA   16161.6666238992 15711.0154842336 15669.191016971
##Julian day as a covariance matrix.

#GOA ONLY
#Year
POD1ma = geeglm(PreAbsM ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mb = geeglm(PreAbsM ~ Year, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mc = geeglm(PreAbsM ~ bs(Year), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model1mA<-c("POD0m", "POD1ma", "POD1mb","POD1mc")
QIC1mA<-c(QIC(POD0m)[1],QIC(POD1ma)[1],QIC(POD1mb)[1],QIC(POD1mc)[1])
QICmod1mA<-data.frame(rbind(model1mA,QIC1mA))
QICmod1mA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model1mA            POD0m           POD1ma           POD1mb           POD1mc
#QIC1mA   57273.6041843083 55091.8989342703 56876.9193908313 55403.0108099303
#Year as a factor.

#TimeLost
POD2ma = geeglm(PreAbsM ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mb = geeglm(PreAbsM ~ TimeLost, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mc = geeglm(PreAbsM ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model2mA<-c("POD0m", "POD2ma", "POD2mb","POD2mc")
QIC2mA<-c(QIC(POD0m)[1],QIC(POD2ma)[1],QIC(POD2mb)[1],QIC(POD2mc)[1])
QICmod2mA<-data.frame(rbind(model2mA,QIC2mA))
QICmod2mA
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
#QIC2jA   66581.7159280241 66748.0688240194 66602.0567455407 66621.7408878972
#TimeLost as linear.

#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
#QIC2mA   16161.6666238992 16176.1972135311 16168.9108400892 16170.8901930977
#TimeLost as linear.

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work).
#The initial full model is:
#GOA (Year as a factor)
#Females
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
#GOA
#QIC            QIC.1           QIC.2            QIC.3            QIC.4            QIC.5
#model3fA            POD0f           POD3fa          POD3fb           POD3fc           POD3fd          POD3fdd
#QIC3fA   8456.24692284754 54581.1064863982 52676.242533809 7783.01585606247 84558.0624891759 7556.73377769405
#Remove site.

#The  full model without Site is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
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
#QIC3fB   8456.24692284754 7556.73377769405 8167.08900751181 8038.74409420976 7546.99867982881
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
PODFinalf = POD3fh
#GOA
#Full model is final model
#QIC            QIC.1           QIC.2            QIC.3
#model3fC            POD0f           POD3fh          POD3fi           POD3fj
#QIC3fC   8456.24692284754 7546.99867982881 8154.7264917688 8032.27430646315

#Juveniles
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
#QIC           QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3jA            POD0j          POD3ja           POD3jb           POD3jc           POD3jd          POD3jdd
#QIC3jA   66581.7159280241 2135636.3721179 1902888.03677668 2152683.45570096 1857714.68363474 64158.4767610496
#Remove Site

#The  full model without TimeLost is:
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
#QIC3jB   66581.7159280241 64158.4767610496 64952.8551634191 65590.4614504138 64124.0402900555
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
#QIC3jC   66581.7159280241 64124.0402900555 64914.0002469341 65569.0171679479
#Full model is the best
PODFinalj = POD3jh

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
#QIC3mA   57273.6041843083 1707857.12191111 1446439.48446054 54368.0074282296 1047081.2450656 53755.8909065657
#Without Site

#The  full model without Site is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mg = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#Without TimeLost
POD3mgg = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3me","POD3mf","POD3mg","POD3mgg")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1],QIC(POD3mgg)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
QICmod3mB
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3mB            POD0m           POD3me           POD3mf           POD3mg          POD3mgg
#QIC3mB   57273.6041843083 53755.8909065657 55105.0864548706 55683.8142192286 55683.8142192286
#Model without TimeLost is the best
PODFinalm = POD3mgg


#GOA (Year as a smooth)
#Females
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Site)+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ TimeLost+as.factor(Site)+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fc = geeglm(PreAbsF ~ AvgDayMatF +as.factor(Site)+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Site
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#Without Year
POD3fdd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd","POD3fdd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1],QIC(POD3fdd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#GOA
#QIC            QIC.1           QIC.2            QIC.3            QIC.4            QIC.5
#model3fA            POD0f           POD3fa          POD3fb           POD3fc           POD3fd          POD3fdd
#QIC3fA   8456.24692284754 54581.1064863982 52676.242533809 7783.01585606247 84558.0624891759 7556.73377769405
#Remove site.

#The  full model without Site is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
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
#QIC3fB   8456.24692284754 7556.73377769405 8167.08900751181 8038.74409420976 7546.99867982881
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
PODFinalf = POD3fh
#GOA
#Full model is final model
#QIC            QIC.1           QIC.2            QIC.3
#model3fC            POD0f           POD3fh          POD3fi           POD3fj
#QIC3fC   8456.24692284754 7546.99867982881 8154.7264917688 8032.27430646315

#Juveniles
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
#QIC           QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3jA            POD0j          POD3ja           POD3jb           POD3jc           POD3jd          POD3jdd
#QIC3jA   66581.7159280241 2135636.3721179 1902888.03677668 2152683.45570096 1857714.68363474 64158.4767610496
#Remove Site

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
#QIC3jB   66581.7159280241 64158.4767610496 64952.8551634191 65590.4614504138 64124.0402900555
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
#QIC3jC   66581.7159280241 64124.0402900555 64914.0002469341 65569.0171679479
#Full model is the best
PODFinalj = POD3jh

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
#QIC3mA   57273.6041843083 1707857.12191111 1446439.48446054 54368.0074282296 1047081.2450656 53755.8909065657
#Without Site

#The  full model without Site is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mg = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#Without TimeLost
POD3mgg = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3me","POD3mf","POD3mg","POD3mgg")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1],QIC(POD3mgg)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
QICmod3mB
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3mB            POD0m           POD3me           POD3mf           POD3mg          POD3mgg
#QIC3mB   57273.6041843083 53755.8909065657 55105.0864548706 55683.8142192286 55683.8142192286
#Model without TimeLost is the best
PODFinalm = POD3mgg

#BSAI
#Females
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
#QIC3fA   3184.88930066242 2109.96325853725 3125.00592702048 2109.96989916953 2109.96325853725
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
PODFinalf=POD3fe
#QIC            QIC.1            QIC.2           QIC.3
#model3fB            POD0f           POD3fe           POD3ff          POD3fg
#QIC3fB   3184.88930066242 2109.96989916953 3125.01012024612 2121.0881909267
#The full model is the best.

#Juveniles
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
#QIC3jA   15081.9977983076 14565.3267700078 15009.0771874134 14561.0008757966 14578.7139682383
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
PODFinalj = POD3je
#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model3jB            POD0j           POD3je           POD3jf           POD3jg
#QIC3jB  15081.9977983076 14561.0008757966 15003.5385343253 14573.7082927332
#Full model is the best.

#Males
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
#QIC3mA   16161.6666238992 15625.2542677672 16131.9455258738 15674.8963198763 15615.0088244095
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
PODFinalm = POD3me
#BSAI
#QIC            QIC.1            QIC.2           QIC.3
#model3mB            POD0m           POD3me           POD3mf          POD3mg
#QIC3mB   16161.6666238992 15615.0088244095 16118.4547702909 15669.191016971
#Full model is best

# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

anova(PODFinalf)
#GOA
#AvgDayMatF       4    51 2.122e-10 ***
#as.factor(Year)  7 87105 < 2.2e-16 ***

#BSAI
#AvgDayMatF       4 26.000 3.165e-05 ***
#as.factor(Site)  1 71.837 < 2.2e-16 ***

anova(PODFinalj)
#GOA
#AvgDayMatJ       4 87.987 < 2.2e-16 ***
#as.factor(Year)  7 94.477 < 2.2e-16 ***

#BSAI
#AvgDayMatJ       4 11.8647   0.01839 *
#as.factor(Site)  1  1.1442   0.28476  

anova(PODFinalm)
#GOA
#AvgDayMatM       4 170.68 < 2.2e-16 ***
#as.factor(Year)  7 233.85 < 2.2e-16 ***

#BSAI
#AvgDayMatM       4 30.9612 3.118e-06 ***
#as.factor(Site)  1  8.6132  0.003337 ** 

# STEP 6: Interpretting the summary of the model
#GOA ONLY (Year as factor)
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)
# Call:
#   geeglm(formula = PreAbsF ~ AvgDayMatF + as.factor(Year), family = binomial, 
#          data = SiteHourTableB, id = BlocksF, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err    Wald Pr(>|W|)    
# (Intercept)           -4.453   0.481   85.63  < 2e-16 ***
#   AvgDayMatFADBM1        1.084   0.253   18.33  1.9e-05 ***
#   AvgDayMatFADBM2       -0.387   0.163    5.62    0.018 *  
#   as.factor(Year)2012    0.455   0.532    0.73    0.393    
# as.factor(Year)2013    0.415   0.530    0.61    0.434    
# as.factor(Year)2014   -1.061   0.577    3.39    0.066 .  
# as.factor(Year)2015   -1.006   0.615    2.68    0.102    
# as.factor(Year)2017   -0.408   0.552    0.55    0.460    
# as.factor(Year)2018   -1.629   0.836    3.80    0.051 .  
# as.factor(Year)2019  -41.400   0.480 7452.95  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)    0.881    1.26
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.686   0.362
# Number of clusters:   1558  Maximum cluster size: 96 

#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj)
# Call:
#   geeglm(formula = PreAbsJ ~ AvgDayMatJ + as.factor(Year), family = binomial, 
#          data = SiteHourTableB, id = BlocksJ, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err  Wald Pr(>|W|)    
# (Intercept)          -1.0004  0.1375 52.96  3.4e-13 ***
#   AvgDayMatJADBM1      -0.7593  0.1003 57.33  3.7e-14 ***
#   AvgDayMatJADBM2       0.1233  0.0829  2.21  0.13690    
# as.factor(Year)2012  -0.9220  0.1664 30.72  3.0e-08 ***
#   as.factor(Year)2013  -1.4906  0.1574 89.64  < 2e-16 ***
#   as.factor(Year)2014  -0.8760  0.1740 25.35  4.8e-07 ***
#   as.factor(Year)2015  -0.7551  0.1946 15.05  0.00010 ***
#   as.factor(Year)2017  -0.5534  0.1671 10.97  0.00093 ***
#   as.factor(Year)2018  -0.5633  0.2524  4.98  0.02564 *  
#   as.factor(Year)2019  -0.5048  0.1865  7.33  0.00679 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.01  0.0906
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.821  0.0168
# Number of clusters:   913  Maximum cluster size: 162 

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm)
# Call:
#   geeglm(formula = PreAbsM ~ AvgDayMatM + as.factor(Year), family = binomial, 
#          data = SiteHourTableB, id = BlocksM, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err  Wald Pr(>|W|)    
# (Intercept)          -1.3046  0.1308 99.42  < 2e-16 ***
#   AvgDayMatMADBM1       0.3719  0.0736 25.54  4.3e-07 ***
#   AvgDayMatMADBM2       0.7112  0.0780 83.21  < 2e-16 ***
#   as.factor(Year)2012  -0.4829  0.1532  9.93   0.0016 ** 
#   as.factor(Year)2013  -1.4349  0.1496 92.00  < 2e-16 ***
#   as.factor(Year)2014  -1.1294  0.1544 53.48  2.6e-13 ***
#   as.factor(Year)2015  -0.9613  0.1652 33.86  5.9e-09 ***
#   as.factor(Year)2017  -1.0010  0.1565 40.89  1.6e-10 ***
#   as.factor(Year)2018   0.1366  0.1917  0.51   0.4761    
# as.factor(Year)2019  -0.3311  0.1574  4.43   0.0354 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)    0.995   0.102
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.703   0.027
# Number of clusters:   1558  Maximum cluster size: 96 

#GOA ONLY (Year as smooth)
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalf_smooth = geeglm(PreAbsF ~ AvgDayMatF+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf_smooth)
#Call:
#  geeglm(formula = PreAbsF ~ AvgDayMatF + bs(Year), family = binomial, 
#         data = SiteHourTableB, id = BlocksF, corstr = "ar1")

#Coefficients:
#  Estimate Std.err    Wald Pr(>|W|)    
#(Intercept)      -4.3392  0.2619 274.448  < 2e-16 ***
#  AvgDayMatFADBM1  -2.8977  1.1727   6.106  0.01347 *  
#  AvgDayMatFADBM2  -0.3902  0.4628   0.711  0.39909    
#AvgDayMatFADBM3   2.0965  0.2752  58.029 2.59e-14 ***
#  AvgDayMatFADBM4   0.9193  0.2941   9.770  0.00177 ** 
#  bs(Year)1        -0.4379  0.8837   0.246  0.62026    
#bs(Year)2        -1.0367  0.8840   1.375  0.24092    
#bs(Year)3        -2.7779  0.4128  45.293 1.70e-11 ***
#Estimate Std.err
#(Intercept)    1.547    46.3
#Link = identity 

#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha   0.7659   6.161
#Number of clusters:   1558  Maximum cluster size: 96 

#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalj_smooth = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj_smooth)
#Call:
#  geeglm(formula = PreAbsJ ~ AvgDayMatJ + bs(Year), family = binomial, 
#         data = SiteHourTableB, id = BlocksJ, corstr = "ar1")

#Coefficients:
#  Estimate Std.err  Wald Pr(>|W|)    
#(Intercept)       -0.875   0.129 45.84  1.3e-11 ***
#  AvgDayMatJADBM1   -0.966   0.140 47.52  5.5e-12 ***
#  AvgDayMatJADBM2    0.190   0.166  1.31    0.252    
#AvgDayMatJADBM3    0.345   0.112  9.54    0.002 ** 
#  AvgDayMatJADBM4   -0.489   0.109 20.24  6.8e-06 ***
#  bs(Year)1         -2.897   0.355 66.63  3.3e-16 ***
#  bs(Year)2          0.220   0.278  0.63    0.428    
#bs(Year)3         -0.905   0.183 24.36  8.0e-07 ***
#Estimate Std.err
#(Intercept)    0.994  0.0687
#Link = identity 

#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha    0.838  0.0131
#Number of clusters:   913  Maximum cluster size: 162 

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalm_smooth = geeglm(PreAbsM ~ AvgDayMatM+bs(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm_smooth)
#Call:
#  geeglm(formula = PreAbsM ~ AvgDayMatM + bs(Year), family = binomial, 
#         data = SiteHourTableB, id = BlocksM, corstr = "ar1")

#Coefficients:
#  Estimate Std.err   Wald Pr(>|W|)    
#(Intercept)      -1.4084  0.1014 192.93  < 2e-16 ***
#  AvgDayMatMADBM1  -1.8515  0.2054  81.21  < 2e-16 ***
#  AvgDayMatMADBM2  -0.4271  0.1419   9.06   0.0026 ** 
#  AvgDayMatMADBM3   0.4239  0.1003  17.85  2.4e-05 ***
#  AvgDayMatMADBM4   0.2610  0.1012   6.65   0.0099 ** 
#  bs(Year)1        -2.3565  0.2741  73.91  < 2e-16 ***
#  bs(Year)2        -0.5633  0.2214   6.47   0.0110 *  
#  bs(Year)3        -0.0979  0.1440   0.46   0.4963    
# Estimate Std.err
#(Intercept)     1.02   0.213
#Link = identity 

#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha    0.691  0.0525
#Number of clusters:   1558  Maximum cluster size: 96 

#BSAI
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Site),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)
# Call:
#   geeglm(formula = PreAbsF ~ AvgDayMatF + as.factor(Site), family = binomial, 
#          data = SiteHourTableB, id = BlocksF, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err  Wald Pr(>|W|)    
# (Intercept)          -9.96    3.82  6.78   0.0092 ** 
#   AvgDayMatFADBM1      -1.82    1.31  1.94   0.1638    
# AvgDayMatFADBM2      -9.58    5.10  3.53   0.0602 .  
# as.factor(Site)KS   -30.03    6.78 19.61  9.5e-06 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate  Std.err
# (Intercept)     1144 39605298
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha   0.0105     365
# Number of clusters:   91  Maximum cluster size: 218 

#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Site),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinalj)
# Call:
#   geeglm(formula = PreAbsJ ~ AvgDayMatJ + as.factor(Site), family = binomial, 
#          data = SiteHourTableB, id = BlocksJ, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)         -1.976   0.154 165.54   <2e-16 ***
#   AvgDayMatJADBM1      0.556   0.314   3.14    0.076 .  
# AvgDayMatJADBM2     -0.114   0.240   0.23    0.635    
# as.factor(Site)KS    0.487   0.330   2.18    0.140    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.01   0.251
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.924  0.0253
# Number of clusters:   72  Maximum cluster size: 275

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Site),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalm)
# Call:
#   geeglm(formula = PreAbsM ~ AvgDayMatM + as.factor(Site), family = binomial, 
#          data = SiteHourTableB, id = BlocksM, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)        -1.7392  0.0929 350.45   <2e-16 ***
#   AvgDayMatMADBM1    -0.4087  0.1515   7.28   0.0070 ** 
#   AvgDayMatMADBM2     0.5711  0.1991   8.23   0.0041 ** 
#   as.factor(Site)KS  -0.5821  0.2008   8.40   0.0037 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)    0.994   0.135
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.803  0.0329
# Number of clusters:   149  Maximum cluster size: 132 

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

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

#Males
prm <- predict(PODFinalm, type="response")  
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

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# STEP 8: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
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
title = paste(saveDir,"/BaseR_Julian Day_SocialGroups - ", region,".png",sep="")
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
         title = paste('Julian Day at',region),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_SocialGroups - ", region,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Site
BootstrapParameters2<-rmvnorm(10000, coef(PODFinalf),summary(PODFinalf)$cov.unscaled)
val=4; Variable=SiteHourTableB$Site; xlabel="Site"; ylabel="Probability"  
BootstrapCoefs2<-BootstrapParameters2[, c(1, val)]

ggtitle = paste(saveDir,"/Probability of Site_SocialGroups) - ", region,".png",sep="")
SiteHourTableB$prf = prf
ggplot(SiteHourTableB, aes(x = Site, y = prf)) +
  geom_boxplot(aes(fill = factor(Site)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedSiteCoefs = data.frame(c(BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),BootstrapCoefs2[,2]),
                               as.factor(strrep(c('KS','BD'),times=1)))
colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")

ggtitle = paste(saveDir,"/Site_SocialGroups - ", region,".png",sep="")
ggplot(AdjustedSiteCoefs, aes(Site, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Site at',region))

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
ggtitle = paste(saveDir,"/Probability of Year (SocialGroups) - ", region,".png",sep="")
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

ggtitle = paste(saveDir,"/Year_SocialGroups - ", region,".png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year at',region))

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
title = paste(saveDir,"/BaseR_Julian Day_MidSize - ", region,".png",sep="")
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
         title = paste('Julian Day at',region),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_MidSize - ", region,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device


#Site as factor
ggtitle = paste(saveDir,"/Probability of Site_MidSize) - ", region,".png",sep="")
SiteHourTableB$prj = prj
ggplot(SiteHourTableB, aes(x = Site, y = prj)) +
  geom_boxplot(aes(fill = factor(Site)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedSiteCoefs = data.frame(c(BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),BootstrapCoefs2[,2]),
                               as.factor(strrep(c('KS','BD'),times=1)))
colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")

ggtitle = paste(saveDir,"/Site_MidSize - ", region,".png",sep="")
ggplot(AdjustedSiteCoefs, aes(Site, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Site at',region))

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
ggtitle = paste(saveDir,"/Probability of Year (MidSize) - ", region,".png",sep="")
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

ggtitle = paste(saveDir,"/Year_MidSize - ", region,".png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year at',region))

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
title = paste(saveDir,"/BaseR_Julian Day_Male - ", region,".png",sep="")
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
         title = paste('Julian Day at',region),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_Males - ", region,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

#Site as factor
ggtitle = paste(saveDir,"/Probability of Site_Male) - ", region,".png",sep="")
SiteHourTableB$prm = prm
ggplot(SiteHourTableB, aes(x = Site, y = prm)) +
  geom_boxplot(aes(fill = factor(Site)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedSiteCoefs = data.frame(c(BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),BootstrapCoefs2[,2]),
                               as.factor(strrep(c('KS','BD'),times=1)))
colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")

ggtitle = paste(saveDir,"/Site_Male - ", region,".png",sep="")
ggplot(AdjustedSiteCoefs, aes(Site, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Site at',region))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalm),summary(PODFinalm)$cov.unscaled)
start=4; finish=10; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalm)[,start:finish]*coef(PODFinalm)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~as.factor(PlottingVar1), fit=F, family=binomial, knots=list(PlottingVar1=seq(2011,2019,length=3)))$X[,4:10]
RealFit1<-Basis1%*%coef(PODFinalm)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Year as factor
ggtitle = paste(saveDir,"/Probability of Year (Male) - ", region,".png",sep="")
SiteHourTableB$prm = prm
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

ggtitle = paste(saveDir,"/Year_Male - ", region,".png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year at',region))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device