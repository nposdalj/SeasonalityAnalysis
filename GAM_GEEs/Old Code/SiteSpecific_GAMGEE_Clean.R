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
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/Plots/",site, sep="")
fileName = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
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
fileName2 = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_GAMGEE_ROW.csv")#setting the directory
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
#NOT ON MODEL RESIDUALS
#Day Table
#acf on the 1 day binned data
acf(SiteDayTable$PreAbs, lag.max = 50)
#CB - at day 30 it gets close enough to the CI intervals
#PT - at day 12
#QN - at day 6
#BD - at day 39

#Hour Table
#acf on the 1 hour binned data
acf(SiteHourTable$PreAbs, lag.max = 2000)
acf(SiteHourTable$PreAbs, lag.max = 2000, ylim=c(0,0.1), xlim =c(700,800)) 
#CB - at hour 1135 it goes below the CI intervals (approximately 47 days)
#PT - at hour 345
#QN - at hour 894
#BD - at hour 896

#ON MODEL RESIDUALS
#CB ONLY
BlockMod<-glm(PreAbs~
               bs(Julian)+
               TimeLost+
               as.factor(Year)
             ,data=SiteHourTable,family=binomial)

#Other sites
BlockMod<-glm(PreAbs~
                bs(Julian)+
                TimeLost, data=SiteHourTable,family=binomial)

#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockMod)

#CB
#bs(Julian)           1671.58  3  < 2.2e-16 ***
#TimeLost            1.49  1     0.2219  
#as.factor(Year)      1104.34  7  < 2.2e-16 ***

#PT 
#bs(Julian)  288.960  3     <2e-16 ***
#TimeLost      1.383  1     0.2397   

#QN
#bs(Julian)  152.421  3  < 2.2e-16 ***
#TimeLost     10.449  1   0.001227 ** 

#BD
#bs(Julian)            438.21  3     <2e-16 ***
#TimeLost       1.57  1     0.2102 

summary(BlockMod)
acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1), xlim =c(600,700)) 
ACFval = 682
#CB - at hour 602 it drops below the CI intervals (approximately 25 days)
#PT - at hour 286
#QN - at hour 867
#BD - at hour 682

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
#CB ONLY
GLM1_CB = glm(PreAbs ~ Julian + TimeLost + as.factor(Year), family = binomial, data = SiteHourTableB)
#VIF scores in GLM to work out collinearity:
VIF(GLM1_CB)
#CB
#GVIF Df GVIF^(1/(2*Df))
#Julian          1.388517  1        1.178353
#TimeLost        1.008421  1        1.004202
#as.factor(Year) 1.399990  7        1.024324

#Other sites
GLM1 = glm(PreAbs~Julian+TimeLost,family=binomial,data=SiteHourTableB)
#VIF scores in GLM to work out collinearity:
VIF(GLM1)
#PT
#Julian TimeLost 
#1.000027 1.000027
#QN
#Julian TimeLost 
#1.000364 1.000364 
#BD
#Julian TimeLost 
#1.000101 1.000101

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

#CB
#QIC          QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   62461.7399063735 60235.9376468338 59735.0202180173
#Julian day as a variance covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   10602.0472061315 10325.8214073948 10340.9814184713
#Julian day as a variance covariance matrix

#QN
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   12308.4810495205 12159.0682034341 12131.6772312595
#Julian day as a variance covariance matrix

#BD
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   24170.2243879127 23789.7701233741 23115.5325908096
#Julian day as a variance-covariance matrix.

#Year
POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1c = geeglm(PreAbs ~ bs(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model1A<-c("POD0", "POD1a", "POD1b","POD1c")
QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
QICmod1A<-data.frame(rbind(model1A,QIC1A))
QICmod1A
#CB
#QIC            QIC.1            QIC.2           QIC.3
#model1A             POD0            POD1a            POD1b           POD1c
#QIC1A   62461.7399063735 61014.3232753975 62441.7567814493 60931.500813237
#Year as a smooth.

#TimeLost
POD2a = geeglm(PreAbs ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2b = geeglm(PreAbs ~ TimeLost, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2c = geeglm(PreAbs ~ bs(TimeLost, k=3), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model2A<-c("POD0", "POD2a", "POD2b","POD2c")
QIC2A<-c(QIC(POD0)[1],QIC(POD2a)[1],QIC(POD2b)[1],QIC(POD2c)[1])
QICmod2A<-data.frame(rbind(model2A,QIC2A))
QICmod2A
#CB
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   62461.7399063735 62751.4240991442 62547.0152954542 62582.0658630441
#TimeLost as linear

#PT
#QIC          QIC.1            QIC.2            QIC.3
#model2A             POD0          POD2a            POD2b            POD2c
#QIC2A   10602.0472061315 10613.44085275 10616.4265146059 125143.790899462
#TimeLost as factor

#QN
#QIC            QIC.1           QIC.2            QIC.3
#model2A             POD0            POD2a           POD2b            POD2c
#QIC2A   12308.4810495205 12376.5022314965 12313.076207006 12326.2495230743
#TimeLost as linear

#BD
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   24170.2243879127 24446.4762123491 24184.1694914502 24255.3050259412
#TimeLost as linear

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
PODFinal = POD3e
#CB
#QIC            QIC.1           QIC.2            QIC.3
#model3B             POD0            POD3e           POD3f            POD3g
#QIC3B   62461.7399063735 58489.4165947094 60931.500813237 59735.0202180173
#Full model is the best

#Other sites (without year, timeLost as linear)
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3c = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#QN
#QIC            QIC.1           QIC.2            QIC.3
#model3A             POD0            POD3a           POD3b            POD3c
#QIC3A   12308.4810495205 12135.2016697477 12313.076207006 12131.6772312595
#Remove Time Lost as a variable. Final model is POD3C.

#BD
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   24170.2243879127 23128.3155446675 24184.1694914502 23115.5325908096
#Remove Time Lost as a variable. Final model is POD3C.

#Other sites (without year, timeLost as factor) for PT
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(TimeLost),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ as.factor(TimeLost),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3c = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#PT
#QIC          QIC.1          QIC.2            QIC.3
#model2A             POD0          POD2a          POD2b            POD2c
#QIC2A   10602.0472061315 234782.6515914 10613.44085275 10340.9814184713
#Remove TimeLost as a variable. Final model is POD2C.

#Final model for other sites
PODFinal = POD3c

# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

#In descending order:
#CB
#AvgDayMat
#Year

anova(PODFinal)

#CB
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: PreAbs
#Terms added sequentially (first to last)

#Df     X2 P(>|Chi|)    
#AvgDayMat  4 40.046 4.235e-08 ***
#bs(Year)   3 27.305 5.081e-06 ***
# Retain all covariates. This is the final model.

#For PT,QN,BD only AvgDayMat was a significant variable, so the order doesn't matter...

# STEP 6: Interpretting the summary of the model
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
summary(PODFinal)

#CB
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ AvgDayMat + as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
summary(PODFinal)

# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

#QN
#Call:
 # geeglm(formula = PreAbs ~ AvgDayMat, family = binomial, data = SiteHourTableB, 
         #id = Blocks, corstr = "ar1")

#Coefficients:
#  Estimate Std.err  Wald Pr(>|W|)    
#(Intercept)      -1.592   0.182 76.28   <2e-16 ***
#  AvgDayMatADBM1    0.411   0.245  2.82    0.093 .  
#AvgDayMatADBM2    0.457   0.362  1.60    0.206    
#Estimate Std.err
#(Intercept)     1.01   0.173
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha     0.86  0.0223
#Number of clusters:   18  Maximum cluster size: 868

#PT
#Call:
 # geeglm(formula = PreAbs ~ AvgDayMat, family = binomial, data = SiteHourTableB, 
         #id = Blocks, corstr = "ar1")

#Coefficients:
 # Estimate Std.err   Wald Pr(>|W|)    
#(Intercept)      -2.068   0.125 273.33   <2e-16 ***
 # AvgDayMatADBM1    0.451   0.223   4.10   0.0428 *  
#  AvgDayMatADBM2    0.770   0.296   6.75   0.0094 ** 
#  Estimate Std.err
#(Intercept)    0.997   0.216
#Link = identity 

#Estimated Correlation Parameters:
#  Estimate Std.err
#alpha    0.847  0.0337
#Number of clusters:   52  Maximum cluster size: 287 

#BD
#Call:
  #geeglm(formula = PreAbs ~ AvgDayMat, family = binomial, data = SiteHourTableB, 
   #      id = Blocks, corstr = "ar1")

#Coefficients:
 # Estimate Std.err Wald Pr(>|W|)  
#(Intercept)     -0.0374  0.1436 0.07    0.795  
#AvgDayMatADBM1  -0.1108  0.2761 0.16    0.688  
#AvgDayMatADBM2   0.6156  0.2531 5.92    0.015 *
#Estimate Std.err
#(Intercept)        1  0.0293
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha     0.94  0.0106
#Number of clusters:   26  Maximum cluster size: 683 

#CB
#Call:
 # geeglm(formula = PreAbs ~ AvgDayMat + as.factor(Year), family = binomial, 
         #data = SiteHourTableB, id = Blocks, corstr = "ar1")

#Coefficients:
 # Estimate Std.err Wald Pr(>|W|)   
#(Intercept)           0.7046  0.4057 3.02   0.0824 . 
#AvgDayMatADBM1        0.0163  0.1902 0.01   0.9316   
#AvgDayMatADBM2        0.3098  0.1785 3.01   0.0827 . 
#as.factor(Year)2012  -0.8651  0.4799 3.25   0.0714 . 
#as.factor(Year)2013  -0.8953  0.4552 3.87   0.0492 * 
  #as.factor(Year)2014  -1.2585  0.4488 7.86   0.0050 **
  #as.factor(Year)2015  -1.4030  0.4836 8.42   0.0037 **
  #as.factor(Year)2017  -0.8823  0.4725 3.49   0.0619 . 
#as.factor(Year)2018  -1.0036  0.7794 1.66   0.1979   
#as.factor(Year)2019  -0.3235  0.4257 0.58   0.4474   

#Correlation structure = ar1 
#Estimated Scale Parameters:
  
  #Estimate Std.err
#(Intercept)     1.01  0.0448
#Link = identity 

#Estimated Correlation Parameters:
  #Estimate Std.err
#alpha    0.939 0.00471
#Number of clusters:   80  Maximum cluster size: 603 

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
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45� line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45� line and the ROC curve

#Find max and position of max
max = max(d,na.rm=TRUE)

if (site == 'BD') {
position = which.max(d$c.0..0.00137080191912269..0.00468357322366918..0.00742517706191455..)
} else if (site == 'CB'){
  position = which.max(d$c.0..0.000523096732657331..0.000845002414292612..0.00116690809592789..)
} else if (site == 'PT'){
  position = which.max(d$c.0..0.00188694079723249..0.00377388159446497..0.00550357732526142..)
} else if (site == 'QN'){
  position = which.max(d$c.0..0.00438756855575868..0.00941499085923217..0.0144424131627057..)
}

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
cutoff = alpha[position,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45� line and the curve

# This value can now be used to build the confusion matrix:

DATA<-matrix(0,dim(SiteHourTableB)[1] ,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim()) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:dim(SiteHourTableB)[1]                                    # the first column is filled with an ID value that is unique for each row
DATA$Observed<-SiteHourTableB$PreAbs                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinal,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

#Confusion matrix:

#CB
#observed
#predicted     1     0
#1 13494 10046
#0  7008 14806

#PT
#observed
#predicted     1     0
#1 1124 5334
#0  608 7385

#QN
#observed
#predicted    1    0
#1 2055 8253
#0  265 2687

#BD
#observed
#predicted    1    0
#1 6755 5538
#0 1919 3216

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
title = paste(saveDir,"/BaseR_Julian Day - ", site,".png",sep="")
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

ggtitle = paste(saveDir,"/Julian Day - ", site,".png",sep="")

ggsave(
  ggtitle,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

#FOR CB ONLY
#Probability of covariate #2: as.smooth(Year):
BootstrapParameters1<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
start=4; finish=10; Variable=SiteHourTableB$Year; xlabel="Year"; ylabel="Probability"  
BootstrapCoefs1<-BootstrapParameters1[,c(1,start:finish)]

#Year as factor
ggtitle = paste(saveDir,"/Probability of Year - ", site,".png",sep="")
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

ggtitle = paste(saveDir,"/Year - ", site,".png",sep="")
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
