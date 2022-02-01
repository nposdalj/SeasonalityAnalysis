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
dir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/Plots")
saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
fileName = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="")
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = HourTable
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable$Region = as.factor(SiteHourTable$Region)
SiteHourTable$Site = as.factor(SiteHourTable$Site)

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
                 bs(Year,k=4)
               ,data=SiteHourTable,family=binomial)

#Females
acf(residuals(BlockModF), lag.max = 2000, ylim=c(0,0.1))
acf(residuals(BlockModF), lag.max = 1000, ylim=c(0,0.1), xlim =c(150,200))
ACFvalF = 187

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian,k=4)+
                 TimeLost+
                 as.factor(Region)+
                 bs(Year,k=4)
               ,data=SiteHourTable,family=binomial)

#Juveniles
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModJ), lag.max = 1000, ylim=c(0,0.1), xlim =c(0,70))
ACFvalJ = 49

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian,k=4)+
                 TimeLost+
                 as.factor(Region)+
                 bs(Year,k=4)
               ,data=SiteHourTable,family=binomial)

#Males
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockModM), lag.max = 1000, ylim=c(0,0.1), xlim =c(260,280))
ACFvalM = 275

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
# bs(Julian, k = 4)    218.7  4    < 2e-16 ***
#   TimeLost               2.9  1      0.088 .  
# as.factor(Region)     31.7  1    1.8e-08 ***
#   bs(Year, k = 4)      182.7  3    < 2e-16 ***

Anova(BlockModJ)
# bs(Julian, k = 4)     1097  4    < 2e-16 ***
#   TimeLost                15  1    0.00011 ***
#   as.factor(Region)       23  1    1.5e-06 ***
#   bs(Year, k = 4)        176  3    < 2e-16 ***

Anova(BlockModM)
# bs(Julian, k = 4)      597  4     <2e-16 ***
#   TimeLost                 0  1       0.86    
# as.factor(Region)        0  1       0.75    
# bs(Year, k = 4)       1055  3     <2e-16 ***

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
#Social Groups
GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)
#Males
GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)
#Juveniles
GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Site)+Year,family=binomial,data=SiteHourTableB)

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)

#Females
# Julian          1.08  1            1.04
# TimeLost        1.01  1            1.00
# as.factor(Site) 2.39  6            1.08
# Year            2.42  1            1.56

#Juvenile
# Julian          1.09  1            1.04
# TimeLost        1.01  1            1.00
# as.factor(Site) 1.85  6            1.05
# Year            1.94  1            1.39

#Male
# Julian          1.07  1            1.04
# TimeLost        1.01  1            1.00
# as.factor(Site) 2.00  6            1.06
# Year            2.07  1            1.44

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
#QIC0fA   11946.7214698826 11614.9135891646 11623.4765208528
#Julian day as a covariance matrix.

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
# QIC1fA   11946.7214698826 11419.427277472 11808.1407068822 11436.517206088
#Year as mSpline.

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
#QIC2fA   11946.7214698826 12057.8011789883 11954.8727135589 11988.1926486955
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
# QIC0jA   81505.4986388886 80785.7363174957 80770.8731497392
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
#QIC1jA   81505.4986388886 80699.9487504421 81501.3180472977 80877.7526931465
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
#QIC2jA   81505.4986388886 81678.3544482711 81532.8315256467 81554.0566221052
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
#QIC0mA   73511.3121477284 72491.7074305815 72492.2857059766
#Julian day as a Variance-covariance matrix

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
# QIC1mA   73511.3121477284 71130.0258574058 73450.0443632068 71237.0538709762
#Year as mSpline.

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
# QIC2mA   73511.3121477284 73648.9095969306 73538.9064443902 73554.552182489
#TimeLost as linear.

#Region is a factor, no need to check them.

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The Initial full model:
#Social Groups
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Region)+mSpline(Year,
                                                                        knots=quantile(Year, probs=c(0.333,0.666)),
                                                                        Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ TimeLost+as.factor(Region)+mSpline(Year,
                                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                                             Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Timelost
POD3fc = geeglm(PreAbsF ~ AvgDayMatF +as.factor(Region)+mSpline(Year,
                                                                knots=quantile(Year, probs=c(0.333,0.666)),
                                                                Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Region
POD3fd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+mSpline(Year,
                                                      knots=quantile(Year, probs=c(0.333,0.666)),
                                                      Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#Without Year
POD3fdd = geeglm(PreAbsF ~ AvgDayMatF+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fA = c("POD0f","POD3fa","POD3fb","POD3fc","POD3fd","POD3fdd")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1],QIC(POD3fd)[1],QIC(POD3fdd)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA
#                      QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3fA            POD0f           POD3fa           POD3fb           POD3fc           POD3fd          POD3fdd
#QIC3fA   11946.7214698826 12399.8472416558 76046.4247991459 12385.5303963656 11114.6390856615 12893.0578674885
#Remove region

#The initial full model without Region is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fg = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without TimeLost
POD3fgg = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fB = c("POD0f","POD3fe","POD3ff","POD3fg","POD3fgg")
QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1],QIC(POD3fgg)[1])
QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
QICmod3fB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fB            POD0f           POD3fe           POD3ff           POD3fg          POD3fgg
#QIC3fB   11946.7214698826 11114.6390856615 11449.901792683 11630.9255920233 11104.3543241628
#Remove TimeLost

#The full model without TimeLost is:
POD3fh = geeglm(PreAbsF ~ AvgDayMatF+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fi = geeglm(PreAbsF ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fj = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3fC = c(QIC(POD0f)[1],QIC(POD3fh)[1],QIC(POD3fi)[1],QIC(POD3fj)[1])
QICmod3fC<-data.frame(rbind(model3fC,QIC3fC))
QICmod3fC
#QIC            QIC.1            QIC.2            QIC.3
#model3fC            POD0f           POD3fh           POD3fi           POD3fj
#QIC3fC   11946.7214698826 11104.3543241628 11436.517206088 11623.4765208528
#Full model is best.

#Model order
#Year
#AvgDayMat

#Juveniles
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019))+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019))+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019))+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Region
POD3jdd = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd","POD3jdd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1],QIC(POD3jdd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd          POD3jdd
#QIC3jA   81505.4986388886 1582936.15999493 1882839.37807817 1029034.64270664 1889294.29972503 80055.2025546555
#Remove Region.

#The  full model without Region is:
POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jf = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jg = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#Without TimeLost
POD3jee = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jB = c("POD0j","POD3je","POD3jf","POD3jg","POD3jee")
QIC3jB = c(QIC(POD0j)[1],QIC(POD3je)[1],QIC(POD3jf)[1],QIC(POD3jg)[1],QIC(POD3jee)[1])
QICmod3jB<-data.frame(rbind(model3jB,QIC3jB))
QICmod3jB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3jB            POD0j           POD3je           POD3jf           POD3jg          POD3jee
# QIC3jB   81505.4986388886 80055.2025546555 80918.4973930406 80793.8668132026 80020.5259301117
#Remove TimeLost

#The full model without TimeLost is:
POD3jh = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3ji = geeglm(PreAbsJ ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jj = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3jC = c(QIC(POD0j)[1],QIC(POD3jh)[1],QIC(POD3ji)[1],QIC(POD3jj)[1])
QICmod3jC<-data.frame(rbind(model3jC,QIC3jC))
QICmod3jC
#                      QIC            QIC.1           QIC.2            QIC.3
#model3jC            POD0f           POD3fh          POD3fi           POD3fj
#QIC3jC   81505.4986388886 80020.5259301117 80877.7526931465 80770.8731497392
#Final model is best.

#Model Order
#AvgDayMat
#Year

#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019))+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019))+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019))+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Region
POD3mdd = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md","POD3mdd")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1],QIC(POD3mdd)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#                     QIC            QIC.1            QIC.2            QIC.3            QIC.4           QIC.5
#model3mA           POD0m           POD3ma           POD3mb           POD3mc           POD3md         POD3mdd
#QIC3mA   73511.3121477284 70874.62270983 71459.8275170118 73616.537936738 70856.0214748948 70608.0260214968
#Remove region.

#The  full model without Region is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mg = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#Without TimeLost
POD3mgg = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                              knots=quantile(Year, probs=c(0.333,0.666)),
                                              Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3me","POD3mf","POD3mg","POD3mgg")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1],QIC(POD3mgg)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
#                     QIC           QIC.1            QIC.2            QIC.3           QIC.4
#model3mB           POD0m          POD3me           POD3mf           POD3mg         POD3mgg
#QIC3mB   73511.3121477284 70608.0260214968 71254.8393665503 72518.403885725 70589.9884181587
#Model without TimeLost is best

#The  full model without TimeLost is:
POD3mh = geeglm(PreAbsM ~ AvgDayMatM+mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mi = geeglm(PreAbsM ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mj = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mC = c("POD0m","POD3mh","POD3mi","POD3mj")
QIC3mC = c(QIC(POD0m)[1],QIC(POD3mh)[1],QIC(POD3mi)[1],QIC(POD3mj)[1])
QICmod3mC<-data.frame(rbind(model3mC,QIC3mC))
QICmod3mC
#QIC            QIC.1            QIC.2            QIC.3
#model3mC            POD0m           POD3mh           POD3mi           POD3mj
#QIC3mC   73511.3121477284 70589.9884181587 71237.0538709762 72492.2857059766
#Full model is best

#Model Order
#Year
#AvgDayMat

# Step 7: Finalize Model --------------------------------------------------
#Social Groups
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalF = geeglm(PreAbsF ~ mSpline(Year,
                                    knots=quantile(Year, probs=c(0.333,0.666)),
                                    Boundary.knots=c(2010,2019))+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
PODFinalF_Region = geeglm(PreAbsF ~ mSpline(Year,
                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                     Boundary.knots=c(2010,2019))+AvgDayMatF+as.factor(Region),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#Mid-Size
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                     Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
PODFinalJ_Region = geeglm(PreAbsJ ~ AvgDayMatJ+mSpline(Year,
                                                knots=quantile(Year, probs=c(0.333,0.666)),
                                                Boundary.knots=c(2010,2019))+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalM = geeglm(PreAbsM ~ mSpline(Year,
                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                     Boundary.knots=c(2010,2019))+AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
PODFinalM_Region = geeglm(PreAbsM ~ mSpline(Year,
                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                     Boundary.knots=c(2010,2019))+AvgDayMatM+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinalF)
# Df     X2 P(>|Chi|)    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 30.331 1.269e-05 ***
#   AvgDayMatF                                                                                      2 19.347 6.292e-05 ***
anova(PODFinalF_Region)
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 30.3310 1.269e-05 ***
#   AvgDayMatF                                                                                      2 19.3472 6.292e-05 ***
#   as.factor(Region)                                                                               1  0.5512    0.4578    

anova(PODFinalJ)
# AvgDayMatJ                                                                                      2 39.245 3.007e-09 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 58.023 3.112e-11 ***
anova(PODFinalJ_Region)
# AvgDayMatJ                                                                                      2 39.245 3.007e-09 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 58.023 3.112e-11 ***
#   as.factor(Region)                                                                               1  6.256   0.01238 *  

anova(PODFinalM)
# Df     X2 P(>|Chi|)    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 92.540 < 2.2e-16 ***
#   AvgDayMatM                                                                                      2 18.772 8.388e-05 ***
anova(PODFinalM_Region)
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 92.540 < 2.2e-16 ***
#   AvgDayMatM                                                                                      2 18.772 8.388e-05 ***
#   as.factor(Region)                                                                               1  3.254   0.07126 .  

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
DATA$Predicted<-predict(PODFinalJ,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

#Males
prm <- predict(PODFinalM, type="response")  
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
DATA$Predicted<-predict(PODFinalM,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,'/BigModel_gamgeeOutput_sexClasses.RData',sep="")
save.image(file = fileName)