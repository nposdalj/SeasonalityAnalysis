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
fileName = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(region,Region))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12

#Daily data - for block calculations
fileName2 = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_GAMGEE_ROW.csv")#setting the directory
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
#NOT ON MODEL RESIDUALS
#Day Table
#acf on the 1 day binned data
acf(SiteDayTable$PreAbs, lag.max = 50)
#BSAI - at day 41
#GOA - at day 6

#Hour Table
#acf on the 1 hour binned data
acf(SiteHourTable$PreAbs, lag.max = 2000)
acf(SiteHourTable$PreAbs, lag.max = 2000, ylim=c(0,0.1), xlim =c(300,400)) 
#BSAI - at hour 746
#GOA - at hour 349

#ON MODEL RESIDUALS
#Region Specific
#BSAI
BlockMod<-glm(PreAbs~
                bs(Julian)+
                TimeLost + 
                as.factor(Site), data=SiteHourTable,family=binomial)

#GOA
BlockMod<-glm(PreAbs~
                bs(Julian)+
                TimeLost +
                bs(Year) + 
                as.factor(Site), data=SiteHourTable,family=binomial)

#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockMod)

#BSAI
#bs(Julian)        453.19  3  < 2.2e-16 ***
#TimeLost            1.49  1     0.2218    
#as.factor(Site)    37.71  1  8.213e-10 ***

#GOA
#bs(Julian)        1754.5  3    < 2e-16 ***
#TimeLost             4.1  1    0.04179 *  
#bs(Year)          1003.7  3    < 2e-16 ***
#as.factor(Site)   6314.6  4    < 2e-16 ***

summary(BlockMod)
acf(residuals(BlockMod), lag.max = 600, ylim=c(0,0.1))
acf(residuals(BlockMod), lag.max = 1000, ylim=c(-0.1,0.1), xlim =c(500,600)) 
ACFval = 519
#BSAI - at hour 653
#GOA - at hour 519

#create the blocks based on the full timesereies
startDate = SiteHourTable$tbin[1]
endDate = SiteHourTable$tbin[nrow(SiteHourTable)]
timeseries = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseries)/ACFval)), times=1, each=ACFval)
divdiff = nrow(timeseries) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseries$block = c(preBlock,lastVec)
timeseries = rename(timeseries, tbin = date)
SiteHourTableB = left_join(SiteHourTable,timeseries,by="tbin")

#Make blocks continuous 
gaps = check4Gaps(SiteHourTableB$block,tol=1)
UnBlock = as.data.frame(unique(SiteHourTableB$block))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$Blocks = UnBlock$sequence[match(SiteHourTableB$block,UnBlock$`unique(SiteHourTableB$block)`)]

difference = diff(SiteHourTableB$Blocks) #find difference between rows
gapsCont = check4Gaps(SiteHourTableB$Blocks)
SiteHourTableB$Waves = rep(0,nrow(SiteHourTableB)) #make space for waves

for (i in 1:nrow(gapsCont)){
  StartB = gapsCont$beg.indx[i]
  EndB = gapsCont$end.indx[i]
  Num = length(StartB:EndB)
  SiteHourTableB$Waves[StartB:EndB] = seq.int(1,Num)
}

## Step 4: Data exploration and initial analysis ##
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
#Region Specific
GLM_BD = glm(PreAbs~Julian+TimeLost+Site,family=binomial,data=SiteHourTableB)
#VIF scores in GLM to work out collinearity:
VIF(GLM_BD)
#BSAI
#Julian TimeLost     Site 
#1.000441 1.001597 1.001921 

GLM_GOA = glm(PreAbs~Julian+TimeLost+Site+Year,family=binomial,data=SiteHourTableB)
#VIF scores in GLM to work out collinearity:
VIF(GLM_GOA)
#GOA
#Julian   1.081902  1        1.040145
#TimeLost 1.006228  1        1.003109
#Site     1.345288  4        1.037772
#Year     1.430608  1        1.196080

## STEP 4: Model selection - covariate preparation ##
# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.
POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

#Julian Day
POD0a = geeglm(PreAbs ~ bs(Julian, knots=10), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD0b = geeglm(PreAbs ~ AvgDayMat, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model0A<-c("POD0", "POD0a", "POD0b")
QIC0A<-c(QIC(POD0)[1],QIC(POD0a)[1],QIC(POD0b)[1])
QICmod0A<-data.frame(rbind(model0A,QIC0A))
QICmod0A
#BSAI
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   26598.7716063937 26221.8264368644 25615.2847735921
#Julian day as a covariance matrix

#GOA
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   101941.659918885 99560.1332277305 99133.0064315274
#Julian day as a covariance matrix

#Year
POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1c = geeglm(PreAbs ~ bs(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model1A<-c("POD0", "POD1a", "POD1b","POD1c")
QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
QICmod1A<-data.frame(rbind(model1A,QIC1A))
QICmod1A
#Not enough years for BSAI
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model1A             POD0            POD1a            POD1b            POD1c
#QIC1A   101941.659918885 97461.0521473272 101712.657305758 97839.0336235966
#Year as a smooth and factor?

#TimeLost
POD2a = geeglm(PreAbs ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2b = geeglm(PreAbs ~ TimeLost, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2c = geeglm(PreAbs ~ bs(TimeLost, k=3), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model2A<-c("POD0", "POD2a", "POD2b","POD2c")
QIC2A<-c(QIC(POD0)[1],QIC(POD2a)[1],QIC(POD2b)[1],QIC(POD2c)[1])
QICmod2A<-data.frame(rbind(model2A,QIC2A))
QICmod2A
#BSAI
#QIC            QIC.1            QIC.2          QIC.3
#model2A             POD0            POD2a            POD2b          POD2c
#QIC2A   26598.7716063937 26636.0005421372 26610.4495685621 26651.64796247
#TimeLost as linear

#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   101941.659918885 102336.206365645 102051.894489098 102090.903903032
#TimeLost as linear

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work).
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3d = geeglm(PreAbs ~ AvgDayMat+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c","POD3d")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#BSAI
#QIC          QIC.1            QIC.2            QIC.3            QIC.4
#model3A             POD0          POD3a            POD3b            POD3c            POD3d
#QIC3A   26598.7716063937 25552.50349009 26599.4583801071 25624.2308193703 25543.2040816703
#Remove timelost as variable.

#Second round of model testing without timelost
#The initial full model is:
POD3e = geeglm(PreAbs ~ AvgDayMat+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3f = geeglm(PreAbs ~ as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3g = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3B = c("POD0","POD3e","POD3f","POD3g")
QIC3B = c(QIC(POD0)[1],QIC(POD3e)[1],QIC(POD3f)[1],QIC(POD3g)[1])
QICmod3B<-data.frame(rbind(model3B,QIC3B))
QICmod3B
#BSAI
#QIC            QIC.1            QIC.2            QIC.3
#model3B             POD0            POD3e            POD3f            POD3g
#QIC3B   26598.7716063937 25543.2040816703 26587.3271097547 25615.2847735921
#The full model has the lowest QIC.

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work).
#For GOA with Year
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+TimeLost+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ as.factor(Site)+TimeLost+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3d = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#Without Year
POD3e = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c","POD3d","POD3e")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1],QIC(POD3e)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4           QIC.5
#model3A             POD0            POD3a            POD3b            POD3c            POD3d           POD3e
#QIC3A   101941.659918885 2581678.37093259 1762959.93139538 95959.3609289755 2579481.02680015 2670888.8339794
#Remove Site.

#Second round of model testing without Site
#The initial full model is:
POD3f = geeglm(PreAbs ~ AvgDayMat +TimeLost+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3g = geeglm(PreAbs ~ TimeLost+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3h = geeglm(PreAbs ~ AvgDayMat +bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3i = POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3B = c("POD0","POD3f","POD3g","POD3h","POD3i")
QIC3B = c(QIC(POD0)[1],QIC(POD3f)[1],QIC(POD3g)[1],QIC(POD3h)[1],QIC(POD3i)[1])
QICmod3B<-data.frame(rbind(model3B,QIC3B))
QICmod3B
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3B             POD0            POD3f            POD3g            POD3h            POD3i
#QIC3B   101941.659918885 95959.3609289755 97904.5271985591 95908.1063053551 99205.6791394299
#Remove TimeLost
PODFinal = POD3h

# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

#In descending order:
#AvgDayMat
#as.factor(Year)

anova(PODFinal)
#BSAI
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: PreAbs
#Terms added sequentially (first to last)

#Df     X2 P(>|Chi|)    
#AvgDayMat        4 38.370   9.4e-08 ***
  #as.factor(Site)  1  6.576   0.01033 *  
  ---
  #Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  
#GOA
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: PreAbs
#Terms added sequentially (first to last)

#Df     X2 P(>|Chi|)    
#AvgDayMat  4 54.256 4.653e-11 ***
#bs(Year)   3 64.879 5.318e-14 ***

# STEP 6: Interpretting the summary of the model
#BSAI
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ AvgDayMat+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
summary(PODFinal)
# Call:
#   geeglm(formula = PreAbs ~ AvgDayMat + as.factor(Site), family = binomial, 
#          data = SiteHourTableB, id = Blocks, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err Wald Pr(>|W|)  
# (Intercept)         0.0269  0.1386 0.04    0.846  
# AvgDayMatADBM1     -0.1527  0.2539 0.36    0.547  
# AvgDayMatADBM2      0.6205  0.2825 4.82    0.028 *
#   as.factor(Site)KS  -0.5217  0.2153 5.87    0.015 *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)        1  0.0269
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha     0.92  0.0124
# Number of clusters:   30  Maximum cluster size: 654 

#GOA
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ AvgDayMat + as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
summary(PODFinal)
# Call:
#   geeglm(formula = PreAbs ~ AvgDayMat + as.factor(Year), family = binomial, 
#          data = SiteHourTableB, id = Blocks, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err  Wald Pr(>|W|)    
# (Intercept)           0.4339  0.2779  2.44  0.11852    
# AvgDayMatADBM1       -0.0899  0.1409  0.41  0.52345    
# AvgDayMatADBM2        0.3969  0.1279  9.63  0.00191 ** 
#   as.factor(Year)2012  -0.8352  0.3050  7.50  0.00618 ** 
#   as.factor(Year)2013  -1.7570  0.2948 35.52  2.5e-09 ***
#   as.factor(Year)2014  -1.4173  0.3227 19.29  1.1e-05 ***
#   as.factor(Year)2015  -1.3085  0.3608 13.15  0.00029 ***
#   as.factor(Year)2017  -1.1412  0.3234 12.45  0.00042 ***
#   as.factor(Year)2018  -0.5007  0.5577  0.81  0.36929    
# as.factor(Year)2019  -0.5631  0.3189  3.12  0.07749 .  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)    0.996  0.0322
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.846  0.0103
# Number of clusters:   98  Maximum cluster size: 1560

# STEP 6: Construction of the ROC curve     
pr <- predict(PODFinal, type="response")  
pred <- prediction(pr,SiteHourTableB$PreAbs) 
perf <- performance(pred, measure="tpr", x.measure="fpr")   
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5))
#This creates a ROC plot

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
DATA$Observed<-SiteHourTableB$PreAbs                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinal,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:
#BSAI
# observed
# predicted    1    0
# 1 8651 7365
# 0  815 2353

#GOA
# observed
# predicted     1     0
# 1 14451 16254
# 0 12208 37209

  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# STEP 7: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
library(boot)
library(pracma)

#Probability of covariate #1: AvgDayBasisMat:
BootstrapParameters3<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
start=2; finish=3; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinal)[,start:finish]*coef(PODFinal)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,366,length=4)))$X[,2:3]
RealFit3<-Basis3%*%coef(PODFinal)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)

#Base R Plotting
title = paste(saveDir,"/BaseR_Julian Day - ", region,".png",sep="")
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

ggtitle = paste(saveDir,"/Julian Day - ", region,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Site
BootstrapParameters2<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
val=4; Variable=SiteHourTableB$Site; xlabel="Site"; ylabel="Probability"  
BootstrapCoefs2<-BootstrapParameters2[, c(1, val)]

#Site as factor
ggtitle = paste(saveDir,"/Probability of Site - ", region,".png",sep="")
SiteHourTableB$pr = pr
ggplot(SiteHourTableB, aes(x = Site, y = pr)) +
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

ggtitle = paste(saveDir,"/Site - ", region,".png",sep="")
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

#FOR GOA ONLY
#Probability of covariate #2: as.smooth(Year):
BootstrapParameters1<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
start=4; finish=10; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(PODFinal)[,start:finish]*coef(PODFinal)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~as.factor(PlottingVar1), fit=F, family=binomial, knots=list(PlottingVar1=seq(2011,2019,length=3)))$X[,4:10]
RealFit1<-Basis1%*%coef(PODFinal)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)

#Year as factor
ggtitle = paste(saveDir,"/Probability of Year - ", region,".png",sep="")
SiteHourTableB$pr = pr
ggplot(SiteHourTableB, aes(x = Year, y = pr)) +
  geom_boxplot(aes(fill = factor(Year)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year
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

ggtitle = paste(saveDir,"/Year - ", region,".png",sep="")
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
