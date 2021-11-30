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
dir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/Plots")
fileName = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = HourTable
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable$Region = as.factor(SiteHourTable$Region)
SiteHourTable$Site = as.factor(SiteHourTable$Site)

#Daily data - for block calculations
fileName2 = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_GAMGEE_ROW.csv")#setting the directory
DayTable = read.csv(fileName2) #no effort days deleted
DayTable = na.omit(DayTable)
DayTable$tbin = as.Date(DayTable$tbin)
SiteDayTable = DayTable
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
BlockMod<-glm(PreAbs~
               bs(Julian)+
               bs(Year)+
               as.factor(Site)+
                TimeLost
             ,data=SiteHourTable,family=binomial)

#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockMod)
summary(BlockMod)
acf(residuals(BlockMod),lag.max = 1000, ylim=c(-0.1,0.1))
acf(residuals(BlockMod),lag.max = 1000, ylim=c(-0.1,0.1), xlim=c(0,20))
ACFval = 17

#bs(Julian)            1901  3     <2e-16 ***
#TimeLost                 8  1     0.0041 ** 
#bs(Year)               656  3     <2e-16 ***
#as.factor(Site)       7336  5     <2e-16 ***
#as.factor(Region)           0               

#Without region
#bs(Julian)        1901.3  3  < 2.2e-16 ***
#TimeLost             8.2  1   0.004077 ** 
#bs(Year)           655.8  3  < 2.2e-16 ***
#as.factor(Site)   7712.7  6  < 2.2e-16 ***

#without site
#bs(Julian)         1487.43  3  < 2.2e-16 ***
#TimeLost             36.70  1   1.38e-09 ***
#bs(Year)           1596.37  3  < 2.2e-16 ***
#as.factor(Region)   376.34  1  < 2.2e-16 ***

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
GLM1 = glm(PreAbs ~ Julian + TimeLost + Year + as.factor(Site), family = binomial, data = SiteHourTableB)
#VIF scores in GLM to work out collinearity:
VIF(GLM1)
#Julian          1.09  1            1.04
#TimeLost        1.01  1            1.00
#Year            2.16  1            1.47
#as.factor(Site) 2.10  6            1.06

## STEP 4: Model selection - covariate preparation ##
# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.
POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

#Julian Day
POD0a = geeglm(PreAbs ~ bs(Julian, knots=4), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD0b = geeglm(PreAbs ~ AvgDayMat, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model0A<-c("POD0", "POD0a", "POD0b")
QIC0A<-c(QIC(POD0)[1],QIC(POD0a)[1],QIC(POD0b)[1])
QICmod0A<-data.frame(rbind(model0A,QIC0A))
QICmod0A
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   130222.387183502 128122.315690403 127402.382736781
#Julian day as a covariance matrix

#Year
POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1c = geeglm(PreAbs ~ bs(Year,k=4), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model1A<-c("POD0", "POD1a", "POD1b","POD1c")
QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
QICmod1A<-data.frame(rbind(model1A,QIC1A))
QICmod1A
#QIC            QIC.1            QIC.2            QIC.3
#model1A             POD0            POD1a            POD1b            POD1c
#QIC1A   130222.387183502 124067.379854311 130129.233073581 126225.320432525
#Year as factor.

#TimeLost
POD2a = geeglm(PreAbs ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2b = geeglm(PreAbs ~ TimeLost, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2c = geeglm(PreAbs ~ bs(TimeLost, k=3), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model2A<-c("POD0", "POD2a", "POD2b","POD2c")
QIC2A<-c(QIC(POD0)[1],QIC(POD2a)[1],QIC(POD2b)[1],QIC(POD2c)[1])
QICmod2A<-data.frame(rbind(model2A,QIC2A))
QICmod2A
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   130222.387183502 130322.506135796 130229.745853173 130205.595844546
#TimeLost as a smooth.

#Region is a factor, no need to check them.

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The initial full model (with year as a factor) is:
POD3aa = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+bs(TimeLost, k = 3)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3bb = geeglm(PreAbs ~ as.factor(Year)+bs(TimeLost, k = 3)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3cc = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost, k = 3)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3dd = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Region
POD3ee = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+bs(TimeLost, k = 3),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3AA = c("POD0","POD3aa","POD3bb","POD3cc","POD3dd","POD3ee")
QIC3AA = c(QIC(POD0)[1],QIC(POD3aa)[1],QIC(POD3bb)[1],QIC(POD3cc)[1],QIC(POD3dd)[1],QIC(POD3ee)[1][1])
QICmod3AA<-data.frame(rbind(model3AA,QIC3AA))
QICmod3AA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3AA             POD0           POD3aa           POD3bb           POD3cc           POD3dd           POD3ee
#QIC3AA   130222.387183502 121650.160197515 124167.868211663 126794.710541971 121685.484447061 121579.628156675
#Remove Region

#The  full model without region is:
POD3ff = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost, k = 3)+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3gg = geeglm(PreAbs ~ bs(TimeLost, k = 3)+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3hh = geeglm(PreAbs ~ AvgDayMat + as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3ii = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost, k = 3),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3BB = c("POD0","POD3ff","POD3gg",'POD3hh',"POD3ii")
QIC3BB = c(QIC(POD0)[1],QIC(POD3ff)[1],QIC(POD3gg)[1],QIC(POD3hh)[1],QIC(POD3ii)[1])
QICmod3BB<-data.frame(rbind(model3BB,QIC3BB))
QICmod3BB
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3BB             POD0           POD3ff           POD3gg           POD3hh          POD3ii
#QIC3BB   130222.387183502 121579.628156675 124021.107274206 121615.81017731 127397.886268327
#Remove timelost.

#The full model without time lost is:
POD3jj = geeglm(PreAbs ~ AvgDayMat + as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3kk = geeglm(PreAbs ~ as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#Without Year
POD3ll = geeglm(PreAbs ~ AvgDayMat ,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3CC = c("POD0","POD3jj","POD3kk",'POD3ll')
QIC3CC = c(QIC(POD0)[1],QIC(POD3jj)[1],QIC(POD3kk)[1],QIC(POD3ll)[1])
QICmod3CC<-data.frame(rbind(model3CC,QIC3CC))
QICmod3CC
#QIC            QIC.1            QIC.2            QIC.3
#model3CC             POD0           POD3jj           POD3kk           POD3ll
#QIC3CC   130222.387183502 121615.81017731 124067.379854311 127402.382736781
#full model is the best.

#The initial full model (with year as a smooth) is:
POD3aa = geeglm(PreAbs ~ AvgDayMat+bs(Year)+bs(TimeLost)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3bb = geeglm(PreAbs ~ bs(Year)+bs(TimeLost)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3cc = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3dd = geeglm(PreAbs ~ AvgDayMat+bs(Year)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Region
POD3ee = geeglm(PreAbs ~ AvgDayMat+bs(Year)+bs(TimeLost),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3AA = c("POD0","POD3aa","POD3bb","POD3cc","POD3dd","POD3ee")
QIC3AA = c(QIC(POD0)[1],QIC(POD3aa)[1],QIC(POD3bb)[1],QIC(POD3cc)[1],QIC(POD3dd)[1],QIC(POD3ee)[1][1])
QICmod3AA<-data.frame(rbind(model3AA,QIC3AA))
QICmod3AA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3AA             POD0           POD3aa           POD3bb           POD3cc           POD3dd           POD3ee
#QIC3AA   130222.387183502 124057.706684998 126242.619639652 126808.950916833 124079.816109347 124045.237123099
#Remove Timlost

#The  full model without TimeLost is:
POD3ff = geeglm(PreAbs ~ AvgDayMat+bs(Year)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3gg = geeglm(PreAbs ~ bs(Year)+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3hh = geeglm(PreAbs ~ AvgDayMat+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Region
POD3ii = geeglm(PreAbs ~ AvgDayMat+bs(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3BB = c("POD0","POD3ff","POD3gg",'POD3hh',"POD3ii")
QIC3BB = c(QIC(POD0)[1],QIC(POD3ff)[1],QIC(POD3gg)[1],QIC(POD3hh)[1],QIC(POD3ii)[1])
QICmod3BB<-data.frame(rbind(model3BB,QIC3BB))
QICmod3BB
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3BB             POD0           POD3ff           POD3gg           POD3hh          POD3ii
#QIC3BB   130222.387183502 124079.816109347 126279.662155278 126814.070624723 124071.325165137
#Sticking with the full model.

# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

#Decided to still include Region since I'm curious if the two regions are different or not

#In descending order:
#bs(Year)
#AvgDayMat
#as.factor(Region)

dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ as.factor(Year)+AvgDayMat+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
summary(PODFinal)
# Call:
#   geeglm(formula = PreAbs ~ as.factor(Year) + AvgDayMat + as.factor(Region), 
#          family = binomial, data = SiteHourTableB, id = Blocks, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)           -0.3878  0.0753  26.52  2.6e-07 ***
#   as.factor(Year)2011    0.6502  0.0973  44.64  2.4e-11 ***
#   as.factor(Year)2012    0.2981  0.1083   7.57  0.00593 ** 
#   as.factor(Year)2013   -0.7764  0.1265  37.70  8.2e-10 ***
#   as.factor(Year)2014   -0.3868  0.1337   8.37  0.00382 ** 
#   as.factor(Year)2015   -0.2854  0.1451   3.87  0.04925 *  
#   as.factor(Year)2017   -0.1276  0.1340   0.91  0.34090    
# as.factor(Year)2018    0.5581  0.1695  10.84  0.00100 ***
#   as.factor(Year)2019    0.4181  0.1391   9.04  0.00264 ** 
#   AvgDayMatADBM1        -0.1682  0.0434  15.02  0.00011 ***
#   AvgDayMatADBM2         0.4633  0.0401 133.32  < 2e-16 ***
#   as.factor(Region)GOA  -0.1655  0.0907   3.33  0.06821 .  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)    0.989 0.00937
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.762 0.00465
# Number of clusters:   3492  Maximum cluster size: 54 

# STEP 7: Construction of the ROC curve     
pr <- predict(PODFinal, type="response")  
pred <- prediction(pr,SiteHourTableB$PreAbs) 
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
DATA$Observed<-SiteHourTableB$PreAbs                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinal,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

# observed
# predicted     1     0
# 1 22270 22675
# 0 13855 40506
  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# STEP 7: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
library(boot)
library(pracma)

#Probability of covariate #1: AvgDayBasisMat:
BootstrapParameters3<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
start=10; finish=11; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
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
title = paste(saveDir,"/BaseR_Julian Day.png",sep="")
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

ggtitle = paste(saveDir,"/Julian Day.png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

#Probability of covariate #2: as.smooth(Year):
BootstrapParameters1<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
start=2; finish=8; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
BootstrapCoefs1<-BootstrapParameters1[,c(1,start:finish)]

#Year as factor
ggtitle = paste(saveDir,"/Probability of Year.png",sep="")
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

ggtitle = paste(saveDir,"/Year.png",sep="")
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

#Region as factor
BootstrapParameters2<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
val=12; Variable=SiteHourTableB$Region; xlabel="Region"; ylabel="Probability"  
BootstrapCoefs2<-BootstrapParameters2[, c(1, val)]

ggtitle = paste(saveDir,"/Probability of Region.png",sep="")
SiteHourTableB$pr = pr
ggplot(SiteHourTableB, aes(x = Region, y = pr)) +
  geom_boxplot(aes(fill = factor(Region)), alpha = .2)

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedSiteCoefs = data.frame(c(BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),BootstrapCoefs2[,2]),
                               as.factor(strrep(c('BSAI','GOA'),times=1)))
colnames(AdjustedSiteCoefs) = c("Coefficient", "Region")

ggtitle = paste(saveDir,"/Region.png",sep="")
ggplot(AdjustedSiteCoefs, aes(Region, Coefficient)
) + geom_boxplot(
) + theme(axis.line = element_line(),
          panel.background = element_blank()
) + labs(title = paste('Region'))

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device