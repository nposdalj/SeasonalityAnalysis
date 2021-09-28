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
site = 'BD'
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
BlockMod<-glm(PreAbsF~
               bs(Julian)+
               TimeLost+
               as.factor(Year)
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

Anova(BlockModJ)
#BD
#bs(Julian)  227.746  3     <2e-16 ***
  #TimeLost      0.004  1     0.9496 

Anova(BlockModM)
#BD
#bs(Julian)   455.22  3     <2e-16 ***
  #TimeLost       0.04  1     0.8337 

#Females
acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1), xlim =c(200,300)) 
ACFvalF = 207
#BD Females - 207

#Juveniles
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1), xlim =c(300,400)) 
ACFvalJ = 331
#BD Juveniles - 331

#Males
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1), xlim =c(100,200)) 
ACFvalM = 143
#BD Males - 143

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
GLM1_CB = glm(PreAbs ~ Julian + TimeLost + as.factor(Year), family = binomial, data = SiteHourTableB)

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
POD0fa = geeglm(PreAbsF ~ bs(Julian, knots=10), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
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

#CB ONLY
#Year
POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1c = geeglm(PreAbs ~ bs(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model1fA<-c("POD0", "POD1a", "POD1b","POD1c")
QIC1fA<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
QICmod1A<-data.frame(rbind(model1A,QIC1A))
QICmod1A
#CB
#QIC            QIC.1            QIC.2           QIC.3
#model1A             POD0            POD1a            POD1b           POD1c
#QIC1A   62461.7399063735 61014.3232753975 62441.7567814493 60931.500813237
#Year as a smooth has the lower QIC but it should really be a factor, right?
#Kept it as a smooth

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

#Juvenile
POD0j<-geeglm(PreAbsJ ~ 1, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Julian Day
POD0ja = geeglm(PreAbsJ ~ bs(Julian), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
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
#TimeLost as linear

#Male
POD0m<-geeglm(PreAbsM ~ 1, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

#Julian Day
POD0ma = geeglm(PreAbsM ~ bs(Julian, knots=10), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
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

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#CB (with Year)
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ bs(Year)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3d = geeglm(PreAbs ~ AvgDayMat+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c","POD3d")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#CB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3A             POD0            POD3a            POD3b            POD3c            POD3d
#QIC3A   62461.7399063735 58529.9166816421 60981.3235686168 59792.9947933176 58489.4165947094
#Remove TimeLost

#The  full model without TimeLost is:
POD3e = geeglm(PreAbs ~ AvgDayMat+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3f = geeglm(PreAbs ~ bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3g = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3B = c("POD0","POD3e","POD3f","POD3g")
QIC3B = c(QIC(POD0)[1],QIC(POD3e)[1],QIC(POD3f)[1],QIC(POD3g)[1])
QICmod3B<-data.frame(rbind(model3B,QIC3B))
QICmod3B
#CB
#QIC            QIC.1           QIC.2            QIC.3
#model3B             POD0            POD3e           POD3f            POD3g
#QIC3B   62461.7399063735 58489.4165947094 60931.500813237 59735.0202180173
#Full model is the best

#Other sites (without year)
#Females - TimeLost as linear (BD, )
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
PODFinalf = POD3fa
#BD
#QIC            QIC.1           QIC.2            QIC.3
#model3fA            PODf0           POD3fa           POD3fb           POD3fc
#QIC3fA   3124.96557501026 2137.60559920481 3124.95351379776 2137.61197218677
#Full model has the lowest QIC. POD3fa


#Juveniles - TimeLost as linear (BD, )
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

# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

#BD FEMALES HAD 2 VARIABLES THAT WERE SIGNIFICANT
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD3faa = geeglm(PreAbsF ~ TimeLost+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fAA = c("PODf0","POD3fa","POD3faa")
QIC3fAA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3faa)[1])
QICmod3fAA<-data.frame(rbind(model3fAA,QIC3fAA))
QICmod3fAA
#QIC            QIC.1            QIC.2
#model3fAA            PODf0           POD3fa          POD3faa
#QIC3fAA   3124.96557501026 2137.60559920481 2137.60559920481
#QIC is exactly the same so we'll keep the original model order.
PODFinalf = POD3fa

anova(PODFinalf)
#BD
#AvgDayMatF  4 33.200 1.087e-06 ***
  #TimeLost    1  5.379   0.02039 * 

anova(PODFinalj)
#BD
#AvgDayMatJ  4 9.177   0.05682 .

anova(PODFinalm)
#AvgDayMatM  4 26.788 2.194e-05 ***

# STEP 6: Interpretting the summary of the model
#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalf)
#BD
#Coefficients:
 # Estimate    Std.err   Wald Pr(>|W|)    
#(Intercept)     -6.8316791  0.9492835 51.792 6.17e-13 ***
 # AvgDayMatFADBM1  3.8402428  0.7203610 28.420 9.77e-08 ***
#  AvgDayMatFADBM2 -7.6374837  6.5474283  1.361   0.2434    
#AvgDayMatFADBM3  0.9035787  0.7041580  1.647   0.1994    
#AvgDayMatFADBM4 -5.3540664  7.8054434  0.471   0.4928    
#TimeLost        -0.0020792  0.0008965  5.379   0.0204 *  
  ---
#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)    482.3 2360989
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha   0.8862     532
#Number of clusters:   85  Maximum cluster size: 208 

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
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc",k=10), fit=F, family=binomial, knots=list(PlottingVar3=seq(1,366,length=10)))$X[,2:5]
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

#Probability of covariate #2: TimeLost:
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalf),summary(PODFinalf)$cov.unscaled)
start=6; finish=6; Variable=SiteHourTableB$TimeLost; xlabel="Time Lost"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalf)[,start:finish]*coef(PODFinalf)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~PlottingVar1, fit=F, family=binomial)
RealFit1<-Basis1%*%coef(PODFinalf)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)
plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(2011,2019), main ="Year" , cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(PlottingVar1)
