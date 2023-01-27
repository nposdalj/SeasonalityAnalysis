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
library(car) #for ANOVA
library(splines2)       # to use mSpline for the GEEs
library(ggfortify)      # extract confidence interval for ACF plots

site = 'CB' #specify the site of interest

# Step 1: Load the data ---------------------------------------------------
#Hourly data
dir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/Plots/",site, sep="")
saveWorkspace = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="") #setting the directory
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


# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
if (site == 'CB'){
#CB
#Females
BlockModF<-glm(PreAbsF~
               bs(Julian, k = 4)+
               TimeLost+
               as.factor(Year)
             ,data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian, k = 4)+
                 TimeLost+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian, k = 4)+
                 TimeLost+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)
}else{
#Other sites
#Females
BlockModF<-glm(PreAbsF~
                bs(Julian, k = 4)+
                TimeLost, data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian, k = 4)+
                 TimeLost, data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian, k = 4)+
                 TimeLost, data=SiteHourTable,family=binomial)
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
#ANOVA with car package
Anova(BlockModF)
#BD
#LR Chisq Df Pr(>Chisq)    
# bs(Julian, k = 4)   791.30  4     <2e-16 ***
#   TimeLost              0.23  1     0.6294  

#PT
# bs(Julian, k = 4)      452  4     <2e-16 ***
#   TimeLost                 2  1       0.17    

#QN
# bs(Julian, k = 4)     98.3  4     <2e-16 ***
#   TimeLost               2.4  1       0.12 

#CB
# bs(Julian, k = 4)     79.8  4     <2e-16 ***
#   TimeLost               2.4  1       0.12    
# as.factor(Year)      212.1  7     <2e-16 ***

Anova(BlockModJ)
#BD
# bs(Julian, k = 4)  272.724  4     <2e-16 ***
#   TimeLost             0.104  1     0.7465     

#PT
# bs(Julian, k = 4)    235.4  4     <2e-16 ***
#   TimeLost               0.7  1       0.41    

#QN
# bs(Julian, k = 4)  293.079  4     <2e-16 ***
#   TimeLost             0.028  1      0.13

#CB
# bs(Julian, k = 4)     1648  4     <2e-16 ***
#   TimeLost                 5  1      0.019 *  
#   as.factor(Year)        207  7     <2e-16 ***

Anova(BlockModM)
#BD
# bs(Julian, k = 4)   348.15  4     <2e-16 ***
#   TimeLost              0.80  1     0.3717     

#PT
# bs(Julian, k = 4)     43.9  4    6.8e-09 ***
#   TimeLost               0.0  1       0.96    

#QN
# bs(Julian, k = 4)   482.34  4     <2e-16 ***
#   TimeLost              2.34  1     0.27  

#CB
# bs(Julian, k = 4)     1430  4    < 2e-16 ***
#   TimeLost                12  1    0.00063 ***
#   as.factor(Year)        854  7    < 2e-16 ***


# Step 4: Data Exploration and Initial Analysis ------------------------------------------------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
if (site == 'CB'){
  #Females
  GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
  #Males
  GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
  #Juveniles
  GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
}else{
#Females
GLMF = glm(PreAbsF~bs(Julian)+TimeLost,family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~bs(Julian)+TimeLost,family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~bs(Julian)+TimeLost,family=binomial,data=SiteHourTableB)
}

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)

#BD
# bs(Julian) 1.000635  3        1.000106
# TimeLost   1.000635  1        1.000317
# 
# bs(Julian) 1.000431  3        1.000072
# TimeLost   1.000431  1        1.000215
# 
# bs(Julian) 1.000654  3        1.000109
# TimeLost   1.000654  1        1.000327

#QN
# bs(Julian)    1  3               1
# TimeLost      1  1               1

# bs(Julian)    1  3               1
# TimeLost      1  1               1

# bs(Julian)    1  3               1
# TimeLost      1  1               1

#CB
# Julian          1.35  1            1.16
# TimeLost        1.02  1            1.01
# as.factor(Year) 1.38  7            1.02

# Julian          1.42  1            1.19
# TimeLost        1.01  1            1.00
# as.factor(Year) 1.43  7            1.03

# Julian          1.51  1            1.23
# TimeLost        1.01  1            1.00
# as.factor(Year) 1.52  7            1.03

#PT
# bs(Julian)    1  3               1
# TimeLost      1  1               1
# 
# bs(Julian)    1  3               1
# TimeLost      1  1               1
# 
# bs(Julian) 1.000172  3        1.000029
# TimeLost   1.000172  1        1.000086


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
#BD
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
# QIC0fA   4751.49393940578 4105.54862731724 4112.1147291363
#Julian day as mSpline, but I'll use a covariance matrix.

#PT
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
# QIC0fA   3463.64294970002 3202.95499728672 3197.46426507714
#Julian day as a variance covariance matrix

#QN
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
# QIC0fA   3607.55084614072 3570.51616432596 3571.65922008805
#Julian day as a covariance matrix

#CB
#QIC            QIC.1            QIC.2
# model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   2332.35543716108 2328.50417490589 2329.59626247959
#Julian day as mSpline, but I'll use covariance matrix.

if (site == 'CB'){
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
}
#CB
# QIC            QIC.1            QIC.2            QIC.3
# model1fA            POD0f           POD1fa           POD1fb           POD1fc
#QIC1fA   2332.35543716108 2159.73362663339 2306.12102082096 2168.63342439418
#Year as factor

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
#QIC2fA   4751.49393940578 5194.535449834 4753.28094461306 4757.54511720217
#TimeLost as linear.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
# QIC2fA   3463.64294970002 3460.48758474816 26025.5039983983 26025.5041290322
#TimeLost as factor, but I'll still use linear.

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
# QIC2fA   3607.55084614072 3613.09066119389 3607.49408899963 3610.21564902294
#TimeLost as linear.

#CB
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
#QIC2fA   2332.35543716108 3248753.19517323 2339.63553588597 2338.88322454545
#TimeLost as linear.

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
#BD
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
#QIC0jA   13611.1246028419 13311.7593271007 13310.4043227396
#Julian day covariance matrix.

#PT
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   6885.14193980651 6696.44634897979 6698.84212655132
#Julian day as mspline but I'll use covariance matrix

#QN
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   7319.30794378126 7098.72481561668 7098.54952080784
#Julian day as a covariance matrix

#CB                      QIC            QIC.1            QIC.2
# QIC            QIC.1            QIC.2
# model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   46285.3794152219 44418.3259438001 44465.9619632154
#Julian day as mSpline, but I'll use covariance matrix.

if (site == 'CB'){
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
}
#CB
#QIC            QIC.1            QIC.2            QIC.3
# QIC            QIC.1            QIC.2            QIC.3
# model1jA            POD0j           POD1ja           POD1jb           POD1jc
# QIC1jA   46285.3794152219 45357.0773456882 46194.2458972643 45731.3178289939
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
#QIC2jA   13611.1246028419 14086.9382444343 13620.7413000152 13630.2790024408
#TimeLost as linear.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
# QIC2jA   6885.14193980651 6891.2804607924 6888.91636373419 6892.03041952167
#TimeLost as linear

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
# QIC2jA   7319.30794378126 7318.10375567827 7318.75386723759 7318.9981127154
#TimeLost as linear.

#CB
#QIC           QIC.1            QIC.2            QIC.3
#model2jA            POD0j          POD2ja           POD2jb           POD2jc
#QIC2jA   46285.3794152219 46425.4439192836 46297.1440587426 46316.4902442723
#TimeLost as linear.

#Male
POD0m<-geeglm(PreAbsM ~ 1, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

#Julian Day
POD0ma = geeglm(PreAbsM ~  mSpline(Julian,
                                   knots=quantile(Julian, probs=c(0.333,0.666)),
                                   Boundary.knots=c(1,366),periodic=T), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD0mb = geeglm(PreAbsM ~ AvgDayMatJ, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model0mA<-c("POD0m", "POD0ma", "POD0mb")
QIC0mA<-c(QIC(POD0m)[1],QIC(POD0ma)[1],QIC(POD0mb)[1])
QICmod0mA<-data.frame(rbind(model0mA,QIC0mA))
QICmod0mA
#BD
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
# QIC0mA   14385.8910442873 14195.1027938734 14194.3609863831
#Julian Day as covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
# QIC0mA   2809.16464589342 2778.86703030812 2778.58654398979
#Julian Day as covariance matrix

#QN
#QIC            QIC.1           QIC.2
#model0mA            POD0m           POD0ma          POD0mb
# QIC0mA   6866.81696443204 6470.76820440076 6469.65140317566
#Julian Day as covariance matrix

#CB
# QIC            QIC.1            QIC.2
# model0mA            POD0m           POD0ma           POD0mb
# QIC0mA   40423.9662520198 39369.8966996129 39378.2074085191
#Julian Day as mSpline, but I'll use covariance matrix

if (site == 'CB'){
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
}
#CB
#QIC            QIC.1            QIC.2            QIC.3
# model1mA            POD0m           POD1ma           POD1mb           POD1mc
# QIC1mA   40423.9662520198 39581.9007019536 40403.7675021959 39652.3872439494
#Year as factor

#TimeLost
POD2ma = geeglm(PreAbsM ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mb = geeglm(PreAbsM ~ TimeLost, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mc = geeglm(PreAbsM ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model2mA<-c("POD0m", "POD2ma", "POD2mb","POD2mc")
QIC2mA<-c(QIC(POD0m)[1],QIC(POD2ma)[1],QIC(POD2mb)[1],QIC(POD2mc)[1])
QICmod2mA<-data.frame(rbind(model2mA,QIC2mA))
QICmod2mA
#BD
#QIC           QIC.1            QIC.2            QIC.3
#model2mA            POD0m          POD2ma           POD2mb           POD2mc
# QIC2mA   14385.8910442873 14685.1523295243 14392.6841956061 14410.0607944982
#TimeLost as linear.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
# QIC2mA   2809.16464589342 2839.16272510763 2817.40671935769 19933.3272719903
#TimeLost as linear

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
# QIC2mA   QIC2mA   6866.81696443204 62829.6417888161 6877.85977627095 63252.8177234494
#TimeLost as linear.

#CB
# QIC           QIC.1            QIC.2            QIC.3
# model2mA            POD0m          POD2ma           POD2mb           POD2mc
#QIC2mA   40423.9662520198 40549.654177705 40448.449207252 40472.8872071309
#TimeLost as linear.


# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
if (site == 'CB'){
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
# QIC            QIC.1            QIC.2            QIC.3           QIC.4
# model3fA            POD0f           POD3fa           POD3fb           POD3fc          POD3fd
# QIC3fA   2332.35543716108 2119.09464555095 2170.23029684459 2335.68477317228 2109.8464014643
#Remove TimeLost.

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
#CB
# QIC           QIC.1            QIC.2            QIC.3
# model3fB            POD0f          POD3fe           POD3ff           POD3fg
# QIC3fB   2332.35543716108 2109.8464014643 2159.73362663339 2329.59626247959
#Full model is the best.
#Year is first.
#Then AvgdayMat
}else{
  #Females
  #The initial full model is:
  POD3fa = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  #without AvgDayMat
  POD3fb = geeglm(PreAbsF ~ TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  #without Timelost
  POD3fc = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  model3fA = c("POD0f","POD3fa","POD3fb","POD3fc")
  QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1])
  QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
  QICmod3fA
  #BD
  # QIC            QIC.1            QIC.2            QIC.3
  # model3fA            POD0f           POD3fa           POD3fb           POD3fc
  #QIC3fA   4751.49393940578 4112.98865807994 4753.28094461306 4112.1147291363
  #Remove Time Lost.
  
  #PT
  #BD
  # QIC            QIC.1            QIC.2            QIC.3
  # model3fA            POD0f           POD3fa           POD3fb           POD3fc
  # QIC3fA   3463.64294970002 26030.1183253369 26025.5039983983 3197.46426507714
  #Remove Time Lost.
  
  #QN
  #QIC            QIC.1            QIC.2            QIC.3
  #model3fA            POD0f           POD3fa           POD3fb           POD3fc
  #QIC3fA   3607.55084614072 3572.83896989477 3607.49408899963 3571.65922008805
  #Remove Time Lost.
  
}

#Juveniles
if (site == 'CB'){
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
# QIC           QIC.1            QIC.2            QIC.3            QIC.4
# model3jA            POD0j          POD3ja           POD3jb           POD3jc           POD3jd
# QIC3jA   46285.3794152219 44186.9227494936 45368.5906566386 44478.9750496734 44174.2397229725
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
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model3jB            POD0j           POD3je           POD3jf           POD3jg
#QIC3jB   46285.3794152219 44174.2397229725 45357.0773456882 44465.9619632154
#Full model is the best.
#AvgDayMat first in model then year.
}else{
  POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  #without AvgDayMat
  POD3jb = geeglm(PreAbsJ ~ TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  #without Timelost
  POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
  model3jA = c("POD0j","POD3ja","POD3jb","POD3jc")
  QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1])
  QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
  QICmod3jA
  #BD
  # QIC            QIC.1            QIC.2            QIC.3
  # model3jA            POD0j           POD3ja           POD3jb           POD3jc
  # QIC3jA   12935.4258268309 12684.3860185494 12954.2204673188 12667.9161993789
  #Remove TimeLost
  
  #PT
  # QIC            QIC.1            QIC.2            QIC.3
  # model3jA            POD0j           POD3ja           POD3jb           POD3jc
  # QIC3jA   6885.14193980651 6701.90896058155 6888.91636373419 6698.84212655132
  #Remove TimeLost.
  
  #QN
  #QIC            QIC.1            QIC.2            QIC.3
  #model3jA            POD0j           POD3ja           POD3jb           POD3jc
  #QIC3jA   7319.30794378126 7099.03692182995 7318.75386723759 7098.54952080784
  #Remvoe Time Lost.
}

#Males
if (site == 'CB'){
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#CB
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3mA            POD0m           POD3ma           POD3mb           POD3mc          POD3md
#QIC3mA   40423.9662520198 38652.8736932768 39596.8818068662 39389.3660536625 38649.9534972491
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
#CB
#QIC           QIC.1            QIC.2            QIC.3
#model3mB            POD0m          POD3me           POD3mf           POD3mg
#QIC3mB   40423.9662520198 38649.9534972491 39581.9007019536 39378.2074085191
#Full model is the best.
#AvgDayMat first, then year.
}else{
  #The initial full model is:
  POD3ma = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without AvgDayMat
  POD3mb = geeglm(PreAbsM ~ TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  #without Timelost
  POD3mc = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
  model3mA = c("POD0m","POD3ma","POD3mb","POD3mc")
  QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1])
  QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
  QICmod3mA
  #BD
  # QIC            QIC.1            QIC.2            QIC.3
  # model3mA            POD0m           POD3ma           POD3mb           POD3mc
  # QIC3mA   14385.8910442873 14200.0927735405 14392.6841956061 14194.3609863831
  #Remove Timelost.
  
  #PT
  # QIC            QIC.1            QIC.2           QIC.3
  # model3mA            POD0m           POD3ma           POD3mb          POD3mc
  # QIC3mA   2809.16464589342 2785.7245829409 2817.40671935769 2778.58654398979
  #Remove Time Lost.
  
  #QN
  #QIC            QIC.1            QIC.2            QIC.3
  #model3mA            POD0m           POD3ma           POD3mb           POD3mc
  #QIC3mA   6866.81696443204 6475.44170865174 6877.85977627095 6469.65140317566
  #Remove Time Lost.
}

# Step 7: Finalize Model --------------------------------------------------
#Females
if (site == 'CB'){
  #CB
  dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalF = geeglm(PreAbsF ~ as.factor(Year)+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalF = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
}

#Juveniles
if (site == 'CB'){
  #CB
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~  AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}

#Males
if (site == 'CB'){
  dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM + as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
}
# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero
anova(PODFinalF)
#BD
#Df   X2 P(>|Chi|)  
#AvgDayMatF  2 39.4   2.8e-09 ***
#PT
# AvgDayMatF  2 14.1   0.00088 ***
#QN
#AvgDayMatF  2 11.5    0.0032 **
#CB
#as.factor(Year)  7 5497    <2e-16 ***
#AvgDayMatF       2    7     0.039 *  

anova(PODFinalJ)
#BD
#Df   X2 P(>|Chi|)   
#AvgDayMatJ  2 19.1   7.2e-05 ***
#PT
# AvgDayMatJ  2 8.3     0.016 *
#QN
#AvgDayMatJ  2 41.3   1.1e-09 ***
#CB
#AvgDayMatJ       2 22.73   1.2e-05 ***
#as.factor(Year)  7  7.32       0.4    
#Remove year
  
anova(PODFinalM)
#BD
#Df   X2 P(>|Chi|)    
#AvgDayMatM  2 24.1     6e-06 ***
#PT
#AvgDayMatM  2 2.62      0.27
#QN
#AvgDayMatM  2 33.3     6e-08 ***
#CB
#AvgDayMatM       2 18.1   0.00012 ***
#as.factor(Year)  7 25.9   0.00053 ***

#Remove variables from model that were not significant after running anova
if (site == 'CB'){
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~  AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}
anova(PODFinalJ)
#CB
#AvgDayMatJ  2 22.7   1.2e-05 ***

filename = paste(saveWorkspace,site,'_SiteSpecific_sexClasses_ModelSummary.txt',sep="")
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
#BD - Female
# observed
# predicted     1     0
# 1   298  2053
# 0    40 15037

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
# 1   149 17766
# 0    29 27410

#QN - Female
# observed
# predicted     1     0
# 1   178  4719
# 0   180 11649
  
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
#BD - Juvenile
# observed
# predicted     1     0
# 1  1131  4950
# 0  1019 10328

#PT - Juvenile
# observed
# predicted     1     0
# 1   395  2870
# 0   434 10752

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

#QN - MidSize
# observed
# predicted    1    0
# 1  656 6896
# 0  236 8938

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
#BD - Males
# observed
# predicted    1    0
# 1 2091 8489
# 0  577 6271

#PT - Males
# observed
# predicted     1     0
# 1   130  1743
# 0   189 12389

#QN - Males
# observed
# predicted    1    0
# 1  824 7262
# 0  164 8476

#CB - Males
# observed
# predicted     1     0
# 1  5767 20413
# 0  1652 17522

#QN - Males
#observed
#predicted    1    0
#1  592 6582
#0  160 5926

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses.RData',sep="")
save.image(file = fileName)