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

# Step 1: Load the Data -----------------------------------------------------------
dir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("H:/My Drive/GofAK_TPWS_metadataReduced/Plots")
saveWorkspace = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
fileName = paste("H:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="")
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = HourTable
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable$Region = as.factor(SiteHourTable$Region)
SiteHourTable$Site = as.factor(SiteHourTable$Site)

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
#Females
BlockModF<-glm(PreAbsF~
                 bs(Julian,k=4)+
                 TimeLost+
                 as.factor(Region)+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

ACF = acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1))
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFvalF = ACFidx[1]

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian,k=4)+
                 TimeLost+
                 as.factor(Region)+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

ACF = acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFvalJ = ACFidx[1]

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian,k=4)+
                 TimeLost+
                 as.factor(Region)+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

ACF = acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFvalM = ACFidx[1]


#Females
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
# bs(Julian, k = 4)    156.9  4     <2e-16 ***
#   TimeLost               1.7  1       0.19    
# as.factor(Region)    201.4  1     <2e-16 ***
#   as.factor(Year)      141.8  8     <2e-16 ***

Anova(BlockModJ)
# bs(Julian, k = 4)     1090  4    < 2e-16 ***
#   TimeLost                18  1    1.9e-05 ***
#   as.factor(Region)       17  1    3.3e-05 ***
#   as.factor(Year)       1049  8    < 2e-16 ***

Anova(BlockModM)
# bs(Julian, k = 4)      920  4     <2e-16 ***
#   TimeLost                 1  1      0.346    
# as.factor(Region)       11  1      0.001 ** 
#   as.factor(Year)       2274  8     <2e-16 ***

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
#Social Groups
GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)

#Females
# Julian          1.32e+00  1            1.15
# TimeLost        1.01e+00  1            1.00
# as.factor(Site) 1.51e+08  6            4.80
# as.factor(Year) 1.64e+08  8            3.26
# 
# Julian           1.37  1            1.17
# TimeLost         1.01  1            1.00
# as.factor(Site) 10.87  6            1.22
# as.factor(Year) 13.17  8            1.17
# 
# Julian          1.46  1            1.21
# TimeLost        1.01  1            1.00
# as.factor(Site) 6.34  6            1.17
# as.factor(Year) 8.14  8            1.14

# Step 5: Model Selection - Covariate Preparation -------------------------
# Construct variance-covariance matrices for cyclic covariates:
#Social Groups
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
#                      QIC            QIC.1            QIC.2
#model0fA            POD0f           POD0fa           POD0fb
#QIC0fA   14728.6846557085 14511.2345484173 14518.0936207307
#Julian day as mspline, but I'll use a covariance matrix.

#Year
POD1fa = geeglm(PreAbsF ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD1fb = geeglm(PreAbsF ~ Year, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD1fc = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019)), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model1fA<-c("POD0f", "POD1fa", "POD1fb","POD1fc")
QIC1fA<-c(QIC(POD0f)[1],QIC(POD1fa)[1],QIC(POD1fb)[1],QIC(POD1fc)[1])
QICmod1fA<-data.frame(rbind(model1fA,QIC1fA))
QICmod1fA
#QIC           QIC.1            QIC.2            QIC.3
#model1fA            POD0f          POD1fa           POD1fb           POD1fc
# QIC1fA   14728.6846557085 14501.2436001281 14644.8359111158 14517.9817015324
#Year as factor

#TimeLost
POD2fa = geeglm(PreAbsF ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fb = geeglm(PreAbsF ~ TimeLost, family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
POD2fc = geeglm(PreAbsF ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model2fA<-c("POD0f", "POD2fa", "POD2fb","POD2fc")
QIC2fA<-c(QIC(POD0f)[1],QIC(POD2fa)[1],QIC(POD2fb)[1],QIC(POD2fc)[1])
QICmod2fA<-data.frame(rbind(model2fA,QIC2fA))
QICmod2fA
#                      QIC            QIC.1            QIC.2            QIC.3
#model2fA            POD0f           POD2fa           POD2fb           POD2fc
#QIC2fA   14728.6846557085 14980.4420950885 14762.0695509189 14785.4541298785
#TimeLost as linear

#Region is a factor, no need to check them.

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
#                      QIC            QIC.1            QIC.2
#model0jA            POD0j           POD0ja           POD0jb
# QIC0jA   84123.4855413207 83245.3084134772 83242.5258626712
#Julian day as a covariance-variance matrix.

#Year
POD1ja = geeglm(PreAbsJ ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jb = geeglm(PreAbsJ ~ Year, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD1jc = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019)), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model1jA<-c("POD0j", "POD1ja", "POD1jb","POD1jc")
QIC1jA<-c(QIC(POD0j)[1],QIC(POD1ja)[1],QIC(POD1jb)[1],QIC(POD1jc)[1])
QICmod1jA<-data.frame(rbind(model1jA,QIC1jA))
QICmod1jA
#QIC            QIC.1            QIC.2           QIC.3
#model1jA            POD0j           POD1ja           POD1jb          POD1jc
#QIC1jA   83462.272117909 82487.3993723317 83454.2312664553 82656.7466525561
#Year as mSpline.

#TimeLost
POD2ja = geeglm(PreAbsJ ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jb = geeglm(PreAbsJ ~ TimeLost, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
POD2jc = geeglm(PreAbsJ ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model2jA<-c("POD0j", "POD2ja", "POD2jb","POD2jc")
QIC2jA<-c(QIC(POD0j)[1],QIC(POD2ja)[1],QIC(POD2jb)[1],QIC(POD2jc)[1])
QICmod2jA<-data.frame(rbind(model2jA,QIC2jA))
QICmod2jA
#QIC            QIC.1            QIC.2            QIC.3
#model2jA            POD0j           POD2ja           POD2jb           POD2jc
#QIC2jA   84123.4855413207 83049.6567212038 84128.4795626085 83200.1690705311
#Time Lost as linear. 

#Region is a factor, no need to check them.

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
#QIC            QIC.1          QIC.2
#model0mA           POD0m           POD0ma         POD0mb
#QIC0mA   73793.3201663778 72613.9822574439 72616.622652085
#Julian day as mSpline, but I'll use a Variance-covariance matrix

#Year
POD1ma = geeglm(PreAbsM ~ as.factor(Year), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mb = geeglm(PreAbsM ~ Year, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD1mc = geeglm(PreAbsM ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019)), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model1mA<-c("POD0m", "POD1ma", "POD1mb","POD1mc")
QIC1mA<-c(QIC(POD0m)[1],QIC(POD1ma)[1],QIC(POD1mb)[1],QIC(POD1mc)[1])
QICmod1mA<-data.frame(rbind(model1mA,QIC1mA))
QICmod1mA
#QIC            QIC.1            QIC.2            QIC.3
#model1mA           POD0m           POD1ma           POD1mb           POD1mc
# QIC1mA   73793.3201663778 71188.1973557482 73742.081249784 71323.6660997178
#Year as factor

#TimeLost
POD2ma = geeglm(PreAbsM ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mb = geeglm(PreAbsM ~ TimeLost, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
POD2mc = geeglm(PreAbsM ~ bs(TimeLost), family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model2mA<-c("POD0m", "POD2ma", "POD2mb","POD2mc")
QIC2mA<-c(QIC(POD0m)[1],QIC(POD2ma)[1],QIC(POD2mb)[1],QIC(POD2mc)[1])
QICmod2mA<-data.frame(rbind(model2mA,QIC2mA))
QICmod2mA
#QIC            QIC.1            QIC.2            QIC.3
#model2mA           POD0m           POD2ma           POD2mb           POD2mc
# QIC2mA   73793.3201663778 74517.536057113 73809.9706155258 73818.8156969022
#TimeLost as linear.

#Region is a factor, no need to check them.

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The Initial full model:
#Social Groups
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Region)+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ TimeLost+as.factor(Region)+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fc = geeglm(PreAbsF ~ AvgDayMatF +as.factor(Region)+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Region
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#Without Year
POD3fdd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd","POD3fdd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1],QIC(POD3fdd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#                      QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3fA            POD0f           POD3fa           POD3fb           POD3fc           POD3fd          POD3fdd
#QIC3fA   14728.6846557085 113968.890098732 14816.2975264557 99844.9920706722 14245.34629925 15522.2901131574
#Remove region

#The initial full model without Region is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fg = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fgg = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fB = c("POD0f","POD3fe","POD3ff","POD3fg","POD3fgg")
QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1],QIC(POD3fgg)[1])
QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
QICmod3fB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fB            POD0f           POD3fe           POD3ff           POD3fg          POD3fgg
#QIC3fB   14728.6846557085 14245.34629925 14547.353663868 14550.3534311585 14207.1224707961
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
#QIC            QIC.1            QIC.2            QIC.3
#model3fC            POD0f           POD3fh           POD3fi           POD3fj
#QIC3fC   14728.6846557085 14207.1224707961 14501.2436001281 14518.0936207307
#Full model is best.

#Model order
#AvgDayMat
#Year

#Juveniles
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ as.factor(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Region
POD3jdd = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd","POD3jdd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1],QIC(POD3jdd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd          POD3jdd
#QIC3jA   84123.4855413207 3403765.58485822 1647251.56567596 89755.4870580481 2180961.15804936 82197.9596154759
#Remove Region.

#The  full model without Region is:
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
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3jB            POD0j           POD3je           POD3jf           POD3jg          POD3jee
# QIC3jB   84123.4855413207 82197.9596154759 83049.1216526588 83234.137549429 82201.6643985028
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
#                      QIC            QIC.1           QIC.2            QIC.3
#model3jC            POD0f           POD3fh          POD3fi           POD3fj
#QIC3jC   84123.4855413207 82201.6643985028 83049.6567212038 83242.5258626712
#Final model is best.

#Model Order
#Year
#AvgDayMat

#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ as.factor(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Region
POD3mdd = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md","POD3mdd")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1],QIC(POD3mdd)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#                     QIC            QIC.1            QIC.2            QIC.3            QIC.4           QIC.5
#model3mA           POD0m           POD3ma           POD3mb           POD3mc           POD3md         POD3mdd
#QIC3mA   73793.3201663778 70362.4979029885 71277.0492352081 73465.0849160167 70351.9393414044 70248.6676445438
#Remove region.

#The  full model without Region is:
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
#                     QIC           QIC.1            QIC.2            QIC.3           QIC.4
#model3mB           POD0m          POD3me           POD3mf           POD3mg         POD3mgg
#QIC3mB   73793.3201663778 70248.6676445438 71199.3681886654 72632.5051403468 70237.7735718734
#Model without TimeLost is best

#The  full model without TimeLost is:
POD3mh = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mi = geeglm(PreAbsM ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mj = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mC = c("POD0m","POD3mh","POD3mi","POD3mj")
QIC3mC = c(QIC(POD0m)[1],QIC(POD3mh)[1],QIC(POD3mi)[1],QIC(POD3mj)[1])
QICmod3mC<-data.frame(rbind(model3mC,QIC3mC))
QICmod3mC
#QIC            QIC.1            QIC.2            QIC.3
#model3mC            POD0m           POD3mh           POD3mi           POD3mj
#QIC3mC   73793.3201663778 70237.7735718734 71188.1973557482 72616.622652085
#Full model is best

#Model Order
#Year
#AvgDayMat

# Step 7: Finalize Model --------------------------------------------------
#Social Groups
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalF = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)

#Mid-Size
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalM = geeglm(PreAbsM ~ as.factor(Year)+AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinalF)
# AvgDayMatF       2 24.1   5.7e-06 ***
#   as.factor(Year)  8 32.4   7.8e-05 ***

anova(PODFinalJ)
#AvgDayMatJ       2  78.6    <2e-16 ***
#as.factor(Year)  8 123.5    <2e-16 ***
  
anova(PODFinalM)
# as.factor(Year)  8 207.7   < 2e-16 ***
#   AvgDayMatM       2  57.2   3.9e-13 ***

#Save model output
filename = paste(saveWorkspace,'_ModelSummary.txt',sep="")
sink(filename)
summary(PODFinalF)
anova(PODFinalF)
summary(PODFinalJ)
anova(PODFinalJ)
summary(PODFinalM)
anova(PODFinalM)
sink(file = NULL)

# Step 9: Construction of the ROC curve    --------------------------------
#Social Groups
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

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

#Mid-Size
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

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,'/BigModel_gamgeeOutput_sexClasses.RData',sep="")
save.image(file = fileName)