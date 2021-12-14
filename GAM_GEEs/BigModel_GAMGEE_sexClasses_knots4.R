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
#QIC3fA   11946.7214698826 691666.993088642 76058.4857555332 806367.420530073 11205.3627115433 1542627.48900464
#Remove region

#The initial full model without Region is:
POD3fe = geeglm(PreAbsF ~ AvgDayMatF+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3ff = geeglm(PreAbsF ~ bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fg = geeglm(PreAbsF ~ AvgDayMatF+TimeLost,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without TimeLost
POD3fgg = geeglm(PreAbsF ~ AvgDayMatF+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fB = c("POD0f","POD3fe","POD3ff","POD3fg","POD3fgg")
QIC3fB = c(QIC(POD0f)[1],QIC(POD3fe)[1],QIC(POD3ff)[1],QIC(POD3fg)[1],QIC(POD3fgg)[1])
QICmod3fB<-data.frame(rbind(model3fB,QIC3fB))
QICmod3fB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3fB            POD0f           POD3fe           POD3ff           POD3fg          POD3fgg
#QIC3fB   11946.7214698826 11205.3627115433 11700.7932397054 11530.6469765503 11199.6058208388
#Remove TimeLost

#The full model without TimeLost is:
POD3fh = geeglm(PreAbsF ~ AvgDayMatF+bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fi = geeglm(PreAbsF ~ bs(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fj = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
model3fC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3fC = c(QIC(POD0f)[1],QIC(POD3fh)[1],QIC(POD3fi)[1],QIC(POD3fj)[1])
QICmod3fC<-data.frame(rbind(model3fC,QIC3fC))
QICmod3fC
PODFinalf = POD3fh
#QIC            QIC.1            QIC.2            QIC.3
#model3fC            POD0f           POD3fh           POD3fi           POD3fj
#QIC3fC   11946.7214698826 11199.6058208388 11691.5078651858 11523.7052388388
#Full model is best

#Juveniles
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ bs(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ +TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Timelost
POD3jd = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year)+as.factor(Region),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Region
POD3jdd = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jA = c("POD0j","POD3ja","POD3jb","POD3jc","POD3jd","POD3jdd")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1],QIC(POD3jd)[1],QIC(POD3jdd)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3jA            POD0j           POD3ja           POD3jb           POD3jc           POD3jd          POD3jdd
#QIC3jA   81512.6839888135 1160427.37539302 1024564.45275153 2175241.01334051 1025336.74174502 80802.4748768065
#Remove Region.


#The  full model without Region is:
POD3je = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jf = geeglm(PreAbsJ ~ bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jg = geeglm(PreAbsJ ~ AvgDayMatJ+TimeLost,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#Without TimeLost
POD3jee = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jB = c("POD0j","POD3je","POD3jf","POD3jg","POD3jee")
QIC3jB = c(QIC(POD0j)[1],QIC(POD3je)[1],QIC(POD3jf)[1],QIC(POD3jg)[1],QIC(POD3jee)[1])
QICmod3jB<-data.frame(rbind(model3jB,QIC3jB))
QICmod3jB
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3jB            POD0j           POD3je           POD3jf           POD3jg          POD3jee
#QIC3jB   81512.6839888135 80802.4748768065 81357.0537542215 80902.9614559011 80779.6699658956
#Remove TimeLost

#The full model without TimeLost is:
POD3jh = geeglm(PreAbsJ ~ AvgDayMatJ+bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3ji = geeglm(PreAbsJ ~ bs(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jj = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
model3jC = c("POD0f","POD3fh","POD3fi","POD3fj")
QIC3jC = c(QIC(POD0j)[1],QIC(POD3jh)[1],QIC(POD3ji)[1],QIC(POD3jj)[1])
QICmod3jC<-data.frame(rbind(model3jC,QIC3jC))
QICmod3jC
PODFinalj = POD3jh
#                      QIC            QIC.1           QIC.2            QIC.3
#model3jC            POD0f           POD3fh          POD3fi           POD3fj
#QIC3jC   81512.6839888135 80779.6699658956 81336.066557015 80884.6310603188
#Final model is best.

#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+bs(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ bs(Year)+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM +TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Timelost
POD3md = geeglm(PreAbsM ~ AvgDayMatM+bs(Year)+as.factor(Region),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Region
POD3mdd = geeglm(PreAbsM ~ AvgDayMatM+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mA = c("POD0m","POD3ma","POD3mb","POD3mc","POD3md","POD3mdd")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1],QIC(POD3md)[1],QIC(POD3mdd)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#                     QIC            QIC.1            QIC.2            QIC.3            QIC.4           QIC.5
#model3mA           POD0m           POD3ma           POD3mb           POD3mc           POD3md         POD3mdd
#QIC3mA   73504.382473807 70975.9591813263 72115.8070266027 72632.4869098921 70961.9706161557 70668.241453907
#Remove region.

#The  full model without Region is:
POD3me = geeglm(PreAbsM ~ AvgDayMatM+bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mf = geeglm(PreAbsM ~ bs(Year)+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mg = geeglm(PreAbsM ~ AvgDayMatM+TimeLost,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#Without TimeLost
POD3mgg = geeglm(PreAbsM ~ AvgDayMatM+bs(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
model3mB = c("POD0m","POD3me","POD3mf","POD3mg","POD3mgg")
QIC3mB = c(QIC(POD0m)[1],QIC(POD3me)[1],QIC(POD3mf)[1],QIC(POD3mg)[1],QIC(POD3mgg)[1])
QICmod3mB<-data.frame(rbind(model3mB,QIC3mB))
QICmod3mB
PODFinalm = POD3mgg
#                     QIC           QIC.1            QIC.2            QIC.3           QIC.4
#model3mB           POD0m          POD3me           POD3mf           POD3mg         POD3mgg
#QIC3mB   73504.382473807 70668.241453907 71910.3164044811 71879.2503940143 70654.526638768
#Model without TimeLost is best


# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

#Decided to still include Region since I'm curious if the two regions are different or not

#In descending order:
#AvgDayMat
#as.factor(Year)

#Females
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
dimnames(AvgYrBasisF)<-list(NULL,c("AYBM1","AYBM2"))
PODFinalf = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
PODFinalfY = geeglm(PreAbsF ~ AvgDayMatF+AvgYrBasisF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
summary(PODFinalfY)
# Call:
#   geeglm(formula = PreAbsF ~ AvgDayMatF + AvgYrBasisF, family = binomial, 
#          data = SiteHourTableB, id = BlocksF, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)        -4.804   0.171 790.28   <2e-16 ***
#   AvgDayMatFADBM1     0.430   0.271   2.51   0.1131    
# AvgDayMatFADBM2    -1.020   0.313  10.61   0.0011 ** 
#   AvgYrBasisFAYBM1    0.668   0.338   3.90   0.0484 *  
#   AvgYrBasisFAYBM2   -0.948   0.389   5.95   0.0147 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.12       5
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.739   0.959
# Number of clusters:   322  Maximum cluster size: 564


#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
dimnames(AvgYrBasisJ)<-list(NULL,c("AYBM1","AYBM2"))
PODFinalj = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
PODFinaljY = geeglm(PreAbsJ ~ AvgDayMatJ+AvgYrBasisJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
summary(PODFinaljY)
# Call:
#   geeglm(formula = PreAbsJ ~ AvgDayMatJ + AvgYrBasisJ, family = binomial, 
#          data = SiteHourTableB, id = BlocksJ, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err    Wald Pr(>|W|)    
# (Intercept)       -1.8008  0.0295 3727.30  < 2e-16 ***
#   AvgDayMatJADBM1   -0.4933  0.0657   56.43  5.8e-14 ***
#   AvgDayMatJADBM2    0.2018  0.0513   15.46  8.4e-05 ***
#   AvgYrBasisJAYBM1  -0.2911  0.0506   33.14  8.6e-09 ***
#   AvgYrBasisJAYBM2   0.0395  0.0697    0.32     0.57    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     0.99  0.0491
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.809  0.0109
# Number of clusters:   2583  Maximum cluster size: 72 

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
dimnames(AvgYrBasisM)<-list(NULL,c("AYBM1","AYBM2"))
PODFinalm = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
PODFinalmY = geeglm(PreAbsM ~ AvgDayMatM+AvgYrBasisM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
summary(PODFinalmY)
# Call:
#   geeglm(formula = PreAbsM ~ AvgDayMatM + AvgYrBasisM, family = binomial, 
#          data = SiteHourTableB, id = BlocksM, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)       -2.0772  0.0368 3184.4  < 2e-16 ***
#   AvgDayMatMADBM1    0.2198  0.0699    9.9   0.0017 ** 
#   AvgDayMatMADBM2    0.4873  0.0752   41.9  9.4e-11 ***
#   AvgYrBasisMAYBM1  -0.4800  0.0592   65.8  4.4e-16 ***
#   AvgYrBasisMAYBM2  -0.5665  0.0812   48.7  3.0e-12 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)    0.999  0.0781
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.727    0.02
# Number of clusters:   1082  Maximum cluster size: 168 

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
BootstrapParameters3<-rmvnorm(10000, coef(PODFinalfY),summary(PODFinalfY)$cov.unscaled)
start=2; finish=3; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinalfY)[,start:finish]*coef(PODFinalfY)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,366,length=4)))$X[,2:3]
RealFit3<-Basis3%*%coef(PODFinalfY)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)

#Base R Plotting
title = paste(saveDir,"/BaseR_Julian Day_SocialGroups.png",sep="")
png(title)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(2011,2019), main = title , cex.lab = 1.5, cex.axis=1.5)    
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
         title = paste('Julian Day'),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_SocialGroups.png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year as factor
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

ggtitle = paste(saveDir,"/Probability of Year (SocialGroups).png",sep="")
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

ggtitle = paste(saveDir,"/Year_SocialGroups.png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year'))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year as smooth
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalfY),summary(PODFinalfY)$cov.unscaled)
start=4; finish=5; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalfY)[,start:finish]*coef(PODFinalfY)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~s(PlottingVar1, bs="cc",k=4), fit=F, family=binomial, knots=list(PlottingVar1=seq(2010,2019,length=4)))$X[,2:3]
RealFit1<-Basis1%*%coef(PODFinalfY)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Base R Plotting
title = paste(saveDir,"/BaseR_YearSmooth_SocialGroups.png",sep="")
png(title)
plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(2010,2019), main = title , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(PlottingVar1)
dev.off()

#ggplot
# Calculate kernel density of Jday observations
dYear = stats::density(Variable,na.rm = TRUE,n=5000,from=2010,to=2019)
dens = data.frame(c(dYear$x, rev(dYear$x)), c(dYear$y, rep(0, length(dYear$y))))
colnames(dens) = c("Year", "Density")
dens$Density = dens$Density / 2 #max(dens$Density) # normalize kernel density
if (min(cis1a[1,])<0){ # set kernel density at bottom of y axis
  dens$Density = dens$Density - abs(min(cis1a[1,])) 
} else {
  dens$Density = dens$Density + min(cis1a[1,])
}

plotDF = data.frame(PlottingVar1, RealFitCenter1a)
colnames(plotDF) = c("Year", "Fit")

ggplot(plotDF, aes(Year, Fit),
) + geom_polygon(data=dens,
                 aes(Year,Density),
                 fill=4,
                 alpha=0.2
) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=cis1a[1,], ymax=cis1a[2,]),
                stat ="identity"
) + labs(x = "Year",
         y = "Probability",
         title = paste('Year'),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Year_SocialGroups.png",sep="")

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
title = paste(saveDir,"/BaseR_Julian Day_MidSize.png",sep="")
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
         title = paste('Julian Day'),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_MidSize.png",sep="")

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
ggtitle = paste(saveDir,"/Probability of Year (MidSize).png",sep="")
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

ggtitle = paste(saveDir,"/Year_MidSize.png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year'))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year as smooth
BootstrapParameters1<-rmvnorm(10000, coef(PODFinaljY),summary(PODFinaljY)$cov.unscaled)
start=4; finish=5; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinaljY)[,start:finish]*coef(PODFinaljY)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~s(PlottingVar1, bs="cc",k=4), fit=F, family=binomial, knots=list(PlottingVar1=seq(2010,2019,length=4)))$X[,2:3]
RealFit1<-Basis1%*%coef(PODFinaljY)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Base R Plotting
title = paste(saveDir,"/BaseR_YearSmooth_MidSize.png",sep="")
png(title)
plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(2010,2019), main = title , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(PlottingVar1)
dev.off()

#ggplot
# Calculate kernel density of Jday observations
dYear = stats::density(Variable,na.rm = TRUE,n=5000,from=2010,to=2019)
dens = data.frame(c(dYear$x, rev(dYear$x)), c(dYear$y, rep(0, length(dYear$y))))
colnames(dens) = c("Year", "Density")
dens$Density = dens$Density / 6 #max(dens$Density) # normalize kernel density
if (min(cis1a[1,])<0){ # set kernel density at bottom of y axis
  dens$Density = dens$Density - abs(min(cis1a[1,])) 
} else {
  dens$Density = dens$Density + min(cis1a[1,])
}

plotDF = data.frame(PlottingVar1, RealFitCenter1a)
colnames(plotDF) = c("Year", "Fit")

ggplot(plotDF, aes(Year, Fit),
) + geom_polygon(data=dens,
                 aes(Year,Density),
                 fill=4,
                 alpha=0.2
) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=cis1a[1,], ymax=cis1a[2,]),
                stat ="identity"
) + labs(x = "Year",
         y = "Probability",
         title = paste('Year'),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/YearSmooth_Midsize.png",sep="")

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
title = paste(saveDir,"/BaseR_Julian Day_Male.png",sep="")
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
         title = paste('Julian Day'),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/Julian Day_Males.png",sep="")

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
ggtitle = paste(saveDir,"/Probability of Year (Male).png",sep="")
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

ggtitle = paste(saveDir,"/Year_Male.png",sep="")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Year'))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year as smooth
BootstrapParameters1<-rmvnorm(10000, coef(PODFinalmY),summary(PODFinalmY)$cov.unscaled)
start=4; finish=5; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinalmY)[,start:finish]*coef(PODFinalmY)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~s(PlottingVar1, bs="cc",k=4), fit=F, family=binomial, knots=list(PlottingVar1=seq(2010,2019,length=4)))$X[,2:3]
RealFit1<-Basis1%*%coef(PODFinalmY)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Base R Plotting
title = paste(saveDir,"/BaseR_YearSmooth_Male.png",sep="")
png(title)
plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(2010,2019), main = title , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(PlottingVar1)
dev.off()

#ggplot
# Calculate kernel density of Jday observations
dYear = stats::density(Variable,na.rm = TRUE,n=5000,from=2010,to=2019)
dens = data.frame(c(dYear$x, rev(dYear$x)), c(dYear$y, rep(0, length(dYear$y))))
colnames(dens) = c("Year", "Density")
dens$Density = dens$Density / 6 #max(dens$Density) # normalize kernel density
if (min(cis1a[1,])<0){ # set kernel density at bottom of y axis
  dens$Density = dens$Density - abs(min(cis1a[1,])) 
} else {
  dens$Density = dens$Density + min(cis1a[1,])
}

plotDF = data.frame(PlottingVar1, RealFitCenter1a)
colnames(plotDF) = c("Year", "Fit")

ggplot(plotDF, aes(Year, Fit),
) + geom_polygon(data=dens,
                 aes(Year,Density),
                 fill=4,
                 alpha=0.2
) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=cis1a[1,], ymax=cis1a[2,]),
                stat ="identity"
) + labs(x = "Year",
         y = "Probability",
         title = paste('Year'),
) + theme(axis.line = element_line(),
          panel.background = element_blank()
)

ggtitle = paste(saveDir,"/YearSmooth_Male.png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device