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

site = 'PT' #specify the site of interest

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
               bs(Year, k = 4)
             ,data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian, k = 4)+
                 TimeLost+
                 bs(Year, k = 4)
               ,data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian, k = 4)+
                 TimeLost+
                 bs(Year, k = 4)
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
if (site == 'PT'){ #The first value that goes below the CI intervals is WAY too high; this is the first lowest number that is reasonable
  ACFvalF = 346
}

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

# ## Manualy doing it
# #Females
# acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1), xlim =c(0,40)) 
# if (site == 'CB'){
#   ACFvalF = 34
# } else if (site == 'PT'){
#   ACFvalF = 34
# } else if (site == 'QN'){
#   ACFvalF = 23
# } else if (site == 'BD'){
#   ACFvalF = 219
# }
# 
# #Juveniles
# acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
# acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1), xlim =c(50,100)) 
# if (site == 'CB'){
#   ACFvalJ = 584
# } else if (site == 'PT'){
#   ACFvalJ = 110
# } else if (site == 'QN'){
#   ACFvalJ = 98
# } else if (site == 'BD'){
#   ACFvalJ = 163
# }
# 
# #Males
# acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
# acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1), xlim =c(800,850))
# if (site == 'CB'){
#   ACFvalM = 422
# } else if (site == 'PT'){
#   ACFvalM = 126
# } else if (site == 'QN'){
#   ACFvalM = 808
# } else if (site == 'BD'){
#   ACFvalM = 140
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
#ANOVA with car package
Anova(BlockModF)
#BD
#LR Chisq Df Pr(>Chisq)    
# bs(Julian, k = 4)  1082.46  4     <2e-16 ***
#   TimeLost              0.89  1     0.3444   

#PT
# bs(Julian, k = 4)  180.854  4     <2e-16 ***
#   TimeLost             0.629  1     0.4277    

#QN
# bs(Julian, k = 4)  123.745  4    < 2e-16 ***
#   TimeLost             5.319  1    0.02109 *  

#CB
#bs(Julian, k = 4)   77.618  4  5.565e-16 ***
#  TimeLost             2.767  1    0.09622 .  
#bs(Year, k = 4)     69.007  3  6.965e-15 ***

Anova(BlockModJ)
#BD
#bs(Julian, k = 4)  310.426  4     <2e-16 ***
#  TimeLost             0.417  1     0.5182    

#PT
# bs(Julian, k = 4)  240.160  4    < 2e-16 ***
#   TimeLost             3.468  1    0.06255 .  

#QN
# bs(Julian, k = 4)  293.079  4     <2e-16 ***
#   TimeLost             0.028  1      0.867 

#CB
#bs(Julian, k = 4)   2096.3  4    < 2e-16 ***
#  TimeLost               5.8  1    0.01598 *  
#  bs(Year, k = 4)      182.6  3    < 2e-16 ***

Anova(BlockModM)
#BD
# bs(Julian, k = 4)   369.96  4     <2e-16 ***
#   TimeLost              0.16  1     0.6887    

#PT
# bs(Julian, k = 4)  116.204  4     <2e-16 ***
#   TimeLost             0.017  1     0.8952    

#QN
# bs(Julian, k = 4)   482.34  4     <2e-16 ***
#   TimeLost              2.34  1     0.1257  

#CB
#bs(Julian, k = 4)  1261.76  4  < 2.2e-16 ***
#TimeLost              9.04  1   0.002638 ** 
#bs(Year, k = 4)     586.41  3  < 2.2e-16 ***


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
GLMF = glm(PreAbsF~Julian+TimeLost,family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~Julian+TimeLost,family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~Julian+TimeLost,family=binomial,data=SiteHourTableB)
}

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)

#BD
# > VIF(GLMF)
# Julian TimeLost 
# 1.000118 1.000118 
# > VIF(GLMJ)
# Julian TimeLost 
# 1.000101 1.000101 
# > VIF(GLMM)
# Julian TimeLost 
# 1.000095 1.000095 

#QN
#VIF(GLMF)
# Julian TimeLost 
# 1.000118 1.000118 
#VIF(GLMJ)
#Julian TimeLost 
# 1.000173 1.000173 
#VIF(GLMM)
#Julian TimeLost 
# 1.000409 1.000409 

#CB
#VIF(GLMF)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.349637  1        1.161739
#TimeLost        1.022190  1        1.011034
#as.factor(Year) 1.378006  7        1.023167
#VIF(GLMJ)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.420667  1        1.191918
#TimeLost        1.009089  1        1.004534
#as.factor(Year) 1.432957  7        1.026029
#VIF(GLMM)
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.505839  1        1.227126
#TimeLost        1.008383  1        1.004183
#as.factor(Year) 1.517947  7        1.030260

#PT
# > VIF(GLMF)
# Julian TimeLost 
# 1        1 
# > VIF(GLMJ)
# Julian TimeLost 
# 1        1 
# > VIF(GLMM)
# Julian TimeLost 
# 1.000003 1.000003 


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
# QIC0fA   3323.88506665023 2364.31222292268 2370.76211786386
#Julian day as a covariance matrix.

#PT
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
# QIC0fA   1586.20038808768 1432.98136422689 1430.67903144549
#Julian day as a covariance matrix

#QN
#QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
# QIC0fA   3467.65822569374 3403.42761798235 3404.98754944958
#Julian day as a covariance matrix

#CB
#QIC            QIC.1            QIC.2
# model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   2331.10672465949 2324.03165188438 2325.15845280872
#Julian day as a covariance matrix

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
#QIC1fA   2331.10672465949 2152.76566182616 2304.5708319038 2162.61274606088
#Year as mSpline smooth

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
#QIC2fA   3323.88506665023 3341.10236013668 3323.87324469837 24561.9008706672
#TimeLost as linear.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
# QIC2fA   1586.20038808768 1581.3568138426 9776.42488661546 10114.8336214469
#TimeLost as factor.

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
# QIC2fA   3467.65822569374 3482.68376078001 3464.8441112808 3466.07554836999
#TimeLost as linear.

#CB
#QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
#QIC2fA   2331.10672465949 13893.6758356065 2341.74526788362 2355.64130793398
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
#QIC0jA   12935.4258268309 12671.6597864823 12667.9161993789
#Julian day as a spline, but I'm using a covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   4371.25803284889 4258.86373742367 4254.23191727409
#Julian day as covariance matrix

#QN
#QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   6967.08452447368 6713.52727309581 6713.54813826824
#Julian day as a covariance matrix

#CB                      QIC            QIC.1            QIC.2
# QIC            QIC.1            QIC.2
# model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   46285.3794152219 44418.3259438001 44465.9619632154
#Julian day as a covariance matrix

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
#Year as mSpline

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
#QIC2jA   12935.4258268309 157581.557900694 12954.2204673188 12979.328579183
#TimeLost as linear.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
# QIC2jA   4371.25803284889 4378.17937272286 35028.3091037713 35027.1767401884
#TimeLost as Factor

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
# QIC2jA   46285.3794152219 46425.4439192836 46297.1440587426 46316.4902442723
#TimeLost as linear.

#CB
#QIC           QIC.1            QIC.2            QIC.3
#model2jA            POD0j          POD2ja           POD2jb           POD2jc
#QIC2jA   46282.9243919683 46427.700926742 46297.0392742771 46314.1783902309
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
# QIC0mA   14799.8337920185 14543.559604596 14542.0721427205
#Julian Day as covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0mA            POD0m           POD0ma           POD0mb
# QIC0mA   1973.75752317758 1861.83637418158 1859.71972457202
#Julian Day as covariance matrix

#QN
#QIC            QIC.1           QIC.2
#model0mA            POD0m           POD0ma          POD0mb
# QIC0mA   7511.96093198324 7077.02603965591 7075.85912254909
#Julian Day as covariance matrix

#CB
# QIC            QIC.1            QIC.2
# model0mA            POD0m           POD0ma           POD0mb
# QIC0mA   40424.2730583234 39370.4552933131 39378.7643820587
#Julian Day as covariance matrix

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
# QIC1mA   40424.2730583234 39583.1224200198 40404.0029235511 39652.7056596262
#Year as an mSpline

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
# QIC2mA   14799.8337920185 14836.8194965385 14809.3228897158 14818.6670360153
#TimeLost as linear.

#PT
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
# QIC2mA   1973.75752317758 1968.53209502516 1974.50869640061 12893.4512396524
#TimeLost as factor

#QN
#QIC            QIC.1            QIC.2            QIC.3
#model2mA            POD0m           POD2ma           POD2mb           POD2mc
# QIC2mA   7511.96093198324 7712.36487446723 7526.62242348286 71069.6809481596
#TimeLost as linear.

#CB
# QIC           QIC.1            QIC.2            QIC.3
# model2mA            POD0m          POD2ma           POD2mb           POD2mc
#QIC2mA   40424.2730583234 40551.5854890866 40448.6077629622 40474.2397506081
#TimeLost as linear.


# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
if (site == 'CB'){
#Females
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fc = geeglm(PreAbsF ~ AvgDayMatF +TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#CB
# QIC            QIC.1            QIC.2            QIC.3           QIC.4
# model3fA            POD0f           POD3fa           POD3fb           POD3fc          POD3fd
# QIC3fA   2331.10672465949 2117.97929396374 2160.64847763955 2334.15262246906 2117.74937039568
#Remove TimeLost.

#The  full model without TimeLost is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fg = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fB = c("POD0f","POD3fe","POD3ff","POD3fg")
QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1])
QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
QICmod3fB
#CB
# QIC           QIC.1            QIC.2            QIC.3
# model3fB            POD0f          POD3fe           POD3ff           POD3fg
# QIC3fB   2331.10672465949 2117.74937039568 2162.61274606088 2325.15845280872
#Full model is the best.
#Year is first.
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
  #QIC3fA   3323.88506665023 2370.7287736748 3323.87324469837 2370.76211786386
  #Remove Time Lost.
  
  #PT
  #BD
  # QIC            QIC.1            QIC.2            QIC.3
  # model3fA            POD0f           POD3fa           POD3fb           POD3fc
  # QIC3fA   1586.20038808768 9171.36715076586 9776.42488661546 1430.67903144549
  #Remove Time Lost.
  
  #QN
  #QIC            QIC.1            QIC.2            QIC.3
  #model3fA            POD0f           POD3fa           POD3fb           POD3fc
  #QIC3fA   3467.65822569374 3402.27379674434 3464.8441112808 3404.98754944958
  #Full model is best.But Since we remove TimeLost for majority of the other models..
  
}

#Juveniles
if (site == 'CB'){
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#CB
# QIC           QIC.1            QIC.2            QIC.3            QIC.4
# model3jA            POD0j          POD3ja           POD3jb           POD3jc           POD3jd
# QIC3jA   46285.3794152219 44246.0909921872 45743.1089644356 44478.9750496734 44233.7851739266
#Remove TimeLost

#The  full model without TimeLost is:
POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jf = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jg = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jB = c("POD0j","POD3je","POD3jf","POD3jg")
QIC3jB = c(QIC(POD0j)[1],QIC(POD3je)[1],QIC(POD3jf)[1],QIC(POD3jg)[1])
QICmod3jB<-data.frame(rbind(model3jB,QIC3jB))
QICmod3jB
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model3jB            POD0j           POD3je           POD3jf           POD3jg
#QIC3jB   46285.3794152219 44233.7851739266 45731.3178289939 44465.9619632154
#Full model is the best.
#AvgDayMat first in model.
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
  # QIC3jA   4371.25803284889 35019.0661240705 35028.3091037713 4254.23191727409
  #Remove TimeLost.
  
  #QN
  #QIC            QIC.1            QIC.2            QIC.3
  #model3jA            POD0j           POD3ja           POD3jb           POD3jc
  #QIC3jA   6967.08452447368 6714.7105686191 6968.0639960482 6713.54813826824
  #Remvoe Time Lost.
}

#Males
if (site == 'CB'){
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ mSpline(Year,
                                   knots=quantile(Year, probs=c(0.333,0.666)),
                                   Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#CB
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3mA            POD0m           POD3ma           POD3mb           POD3mc          POD3md
#QIC3mA   40424.2730583234 38787.8224486274 39672.6091313743 39390.1077164302 38780.5554482199
#Remove TimeLost.

#The  full model without TimeLost is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mg = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3me","POD3mf","POD3mg")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
QICmod3mB
#CB
#QIC           QIC.1            QIC.2            QIC.3
#model3mB            POD0m          POD3me           POD3mf           POD3mg
#QIC3mB   40424.2730583234 38780.5554482199 39652.7056596262 39378.7643820587
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
  # QIC3mA   14799.8337920185 14551.3905845884 14809.3228897158 14542.0721427205
  #Remove Timelost.
  
  #PT
  # QIC            QIC.1            QIC.2           QIC.3
  # model3mA            POD0m           POD3ma           POD3mb          POD3mc
  # QIC3mA   1973.75752317758 1860.33843239953 1974.50869640061 1859.71972457202
  #Remove Time Lost.
  
  #QN
  #QIC            QIC.1            QIC.2            QIC.3
  #model3mA            POD0m           POD3ma           POD3mb           POD3mc
  #QIC3mA   7511.96093198324 7085.60099657606 7526.62242348286 7075.85912254909
  #Remove Time Lost.
}

# Step 7: Finalize Model --------------------------------------------------
#Females
if (site == 'CB'){
  #CB
  dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalF = geeglm(PreAbsF ~ mSpline(Year,
                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                     Boundary.knots=c(2011,2019))+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalF = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
}

#Juveniles
if (site == 'CB'){
  #CB
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~  AvgDayMatJ+mSpline(Year,
                                       knots=quantile(Year, probs=c(0.333,0.666)),
                                       Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
} else {
  dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}

#Males
if (site == 'CB'){
  dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM + mSpline(Year,
                                       knots=quantile(Year, probs=c(0.333,0.666)),
                                       Boundary.knots=c(2011,2019)) ,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
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
# Df   X2 P(>|Chi|)  
#AvgDayMatF  2 5.0355   0.08064 .
#PT
# AvgDayMatF  2 12.78  0.001678 **
  #QN
#AvgDayMatF  4 15.3     0.004 **
#CB
#mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5 17.8677  0.003117 **
 # AvgDayMatF                                                                                      2  5.0595  0.079677 . 
#QN
#Df     X2 P(>|Chi|)
# AvgDayMatF  2 8.6922   0.01296 *

anova(PODFinalJ)
#BD
# Df   X2 P(>|Chi|)   
#AvgDayMatJ  2 8.8342   0.01207 *
#PT
# Df     X2 P(>|Chi|)    
# AvgDayMatJ  2 3.6419    0.1619  
#QN
# AvgDayMatJ  2 42.553 5.751e-10 ***
  #CB
# AvgDayMatJ                                                                                        2 22.7263 1.162e-05***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5  5.0074     0.415  
#QN
#AvgDayMatJ  2 13.277  0.001309 **
  
anova(PODFinalM)
#BD
# Df   X2 P(>|Chi|)    
#AvgDayMatM  2 22.601 1.237e-05 ***
  #PT
# AvgDayMatM  2 4.3566    0.06076 
#QN
#AvgDayMatM  4 10      0.04 *
#CB
#AvgDayMatM                                                                                      2 18.301 0.0001061 ***
#  mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5 15.059 0.0101116 *  
#QN
# AvgDayMatM  2 33 6.825e-08 ***
  
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