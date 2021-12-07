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

site = 'CB' #specify the site of interest

# Step 1: Load the Data -----------------------------------------------------------
#Hourly data
dir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/Plots/",site, sep="")
saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
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
BlockMod<-glm(PreAbs~
               bs(Julian,k=4)+
               TimeLost+
               bs(Year,k=4)
             ,data=SiteHourTable,family=binomial)

}else{

#Other sites
BlockMod<-glm(PreAbs~
                bs(Julian,k=4)+
                TimeLost, data=SiteHourTable,family=binomial)
}

acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1), xlim =c(600,640)) 
#CB - at hour 628 it drops below the CI intervals (approximately 25 days)
#PT - at hour 286
#QN - at hour 867
#BD - at hour 682

if (site == 'CB'){
  ACFval = 628
} else if (site == 'PT'){
  ACFval = 286
} else if (site == 'QN'){
  ACFval = 286
} else if (site == 'BD'){
  ACFval = 682
}

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



# Step 3: ANOVA to Check for Significance of Variables --------------------
#ANOVA with car package
Anova(BlockMod)
summary(BlockMod)

#CB
# Call:
#   glm(formula = PreAbs ~ bs(Julian, k = 4) + TimeLost + bs(Year, 
#                                                            k = 4), family = binomial, data = SiteHourTable)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.560  -1.087  -0.709   1.178   1.996  
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -0.1060     0.1945   -0.55    0.586    
# bs(Julian, k = 4)1  -0.7921     0.1994   -3.97  7.1e-05 ***
#   bs(Julian, k = 4)2   1.1081     0.2100    5.28  1.3e-07 ***
#   bs(Julian, k = 4)3   0.5918     0.2108    2.81    0.005 ** 
#   bs(Julian, k = 4)4   0.9497     0.1957    4.85  1.2e-06 ***
#   TimeLost             0.0223     0.0152    1.46    0.143    
# bs(Year, k = 4)1    -0.1182     0.0474   -2.49    0.013 *  
#   bs(Year, k = 4)2    -1.0626     0.0613  -17.34  < 2e-16 ***
#   bs(Year, k = 4)3    -1.5080     0.0962  -15.67  < 2e-16 ***
#   bs(Year, k = 4)4         NA         NA      NA       NA    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 62456  on 45353  degrees of freedom
# Residual deviance: 59172  on 45345  degrees of freedom
# AIC: 59190
# 
# Number of Fisher Scoring iterations: 4

#PT 
# bs(Julian)    289.0  3     <2e-16 ***
#   TimeLost        1.4  1       0.24    

#QN
# bs(Julian)    152.4  3     <2e-16 ***
#   TimeLost       10.4  1     0.0012 ** 

#BD
#bs(Julian)     438.21  3     <2e-16 ***
#TimeLost       1.57  1     0.2102 


# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
if (site == 'CB'){
  GLM1 = glm(PreAbs ~ Julian + TimeLost + bs(Year), family = binomial, data = SiteHourTableB)
  VIF(GLM1)
#CB
  # GVIF Df GVIF^(1/(2*Df))
  # Julian   1.10  1            1.05
  # TimeLost 1.01  1            1.00
  # bs(Year) 1.11  3            1.02
} else {
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
}


# Step 5: Model Selection - Covariate Preparation -------------------------

# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.
POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

#Julian Day
POD0a = geeglm(PreAbs ~ mSpline(Julian,
                                knots=quantile(Julian, probs=c(0.333,0.666)),
                                Boundary.knots=c(1,366),
                                periodic=T), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD0b = geeglm(PreAbs ~ AvgDayMat, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model0A<-c("POD0", "POD0a", "POD0b")
QIC0A<-c(QIC(POD0)[1],QIC(POD0a)[1],QIC(POD0b)[1])
QICmod0A<-data.frame(rbind(model0A,QIC0A))
QICmod0A

#CB
#QIC          QIC.1            QIC.2
# model0A             POD0           POD0a            POD0b
# QIC0A   62461.7399063735 61609.192902932 61611.5439649097
#Julian day as a variance covariance matrix

#PT
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#IC0A   62461.7399063735 60235.9376468338 59735.0202180173
#Julian day as a variance covariance matrix

#QN
# QIC            QIC.1            QIC.2
# model0A             POD0            POD0a            POD0b
# QIC0A   12302.1776102301 12154.7001195908 12152.7261988382
#Julian day as a variance covariance matrix

#BD
# QIC            QIC.1            QIC.2
# model0A             POD0            POD0a            POD0b
# QIC0A   24170.2243879127 23757.9446246292 23757.1611343918
#Julian day as a variance-covariance matrix.

if (site == 'CB'){
  #Year
  POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  POD1c = geeglm(PreAbs ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),
                                  family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model1A<-c("POD0", "POD1a", "POD1b","POD1c")
  QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
  QICmod1A<-data.frame(rbind(model1A,QIC1A))
  QICmod1A
}

#CB
# QIC            QIC.1            QIC.2           QIC.3
# model1A             POD0            POD1a            POD1b           POD1c
# QIC1A   62461.7399063735 61014.3232753975 62441.7567814493 61024.349983411
#Year as a smooth.

#TimeLost
POD2a = geeglm(PreAbs ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2b = geeglm(PreAbs ~ TimeLost, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2c = geeglm(PreAbs ~ mSpline(TimeLost,
                                knots=4),
                                family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model2A<-c("POD0", "POD2a", "POD2b","POD2c")
QIC2A<-c(QIC(POD0)[1],QIC(POD2a)[1],QIC(POD2b)[1],QIC(POD2c)[1])
QICmod2A<-data.frame(rbind(model2A,QIC2A))
QICmod2A

#CB
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   62461.7399063735 62751.4240991442 62547.0152954542 62585.5738287031
#TimeLost as linear

#PT
# QIC          QIC.1            QIC.2            QIC.3
# model2A             POD0          POD2a            POD2b            POD2c
# QIC2A   10602.0472061315 10613.44085275 10616.4265146059 913773.471987913
#TimeLost as factor

#QN
# QIC            QIC.1            QIC.2            QIC.3
# model2A             POD0            POD2a            POD2b            POD2c
# QIC2A   12302.1776102301 12376.4452089374 12307.0001363939 12327.6245819357
#TimeLost as linear

#BD
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   24170.2243879127 24446.4762123491 24184.1694914502 24261.4981967489
#TimeLost as linear

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
if (site == "cB"){
  #CB (with Year as mSpline)
  #The initial full model is:
  POD3a = geeglm(PreAbs ~ AvgDayMat+mSpline(Year,
                                            knots=quantile(Year, probs=c(0.333,0.666)),
                                            Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without AvgDayMat
  POD3b = geeglm(PreAbs ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019))+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Year
  POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Timelost
  POD3d = geeglm(PreAbs ~ AvgDayMat+mSpline(Year,
                                            knots=quantile(Year, probs=c(0.333,0.666)),
                                            Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model3A = c("POD0","POD3a","POD3b","POD3c","POD3d")
  QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1])
  QICmod3A<-data.frame(rbind(model3A,QIC3A))
  QICmod3A
  #CB
  #QIC            QIC.1            QIC.2            QIC.3            QIC.4
  #model3A             POD0            POD3a            POD3b            POD3c            POD3d
  #QIC3A   62461.7399063735 60594.5985528871 61098.2989861884 61696.2499631834 60528.878840482
  #Remove TimeLost
  
  
  #The  full model without TimeLost is:
  POD3e = geeglm(PreAbs ~ AvgDayMat+mSpline(Year,
                                            knots=quantile(Year, probs=c(0.333,0.666)),
                                            Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without AvgDayMat
  POD3f = geeglm(PreAbs ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Year
  POD3g = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model3B = c("POD0","POD3e","POD3f","POD3g")
  QIC3B = c(QIC(POD0)[1],QIC(POD3e)[1],QIC(POD3f)[1],QIC(POD3g)[1])
  QICmod3B<-data.frame(rbind(model3B,QIC3B))
  QICmod3B
  #CB
  #QIC            QIC.1           QIC.2            QIC.3
  #model3B             POD0            POD3e           POD3f            POD3g
  #QIC3B   62461.7399063735 60528.878840482 61024.349983411 61611.5439649097
  #Full model is the best  
}

if (site == 'BD' | site == 'QN'){
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
# QIC            QIC.1            QIC.2            QIC.3
# model3A             POD0            POD3a            POD3b            POD3c
# QIC3A   12302.1776102301 12156.2162940731 12307.0001363939 12152.7261988382
#Remove Time Lost as a variable. Final model is POD3C.

#BD
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   24170.2243879127 23766.9898540658 24184.1694914502 23757.1611343918
#Remove Time Lost as a variable. Final model is POD3C.
}

if (site == 'PT'){
#Other sites (without year, timeLost as factor)
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
# QIC            QIC.1          QIC.2            QIC.3
# model3A             POD0            POD3a          POD3b            POD3c
# QIC3A   10602.0472061315 10323.1582706194 10613.44085275 10324.6611762466
#Remove TimeLost as a variable. Final model is POD2C.
}


# Step 7: Finalize Model --------------------------------------------------
if (site == 'CB'){
  #CB
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ mSpline(Year,
                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                     Boundary.knots=c(2011,2019)) + AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
} else {
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}

#Testing covariate significance
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

#In descending order:
#CB
#Year
#AvgDaMat
#For PT,QN,BD only AvgDayMat was a significant variable, so the order doesn't matter...

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero
anova(PODFinal)

#CB
# Df      X2
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5 25.9462
# AvgDayMat                                                                                       2  4.7955
# P(>|Chi|)    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019)) 9.141e-05 ***
#   AvgDayMat                                                                                        0.09092 .  

#BD
# Df     X2 P(>|Chi|)  
# AvgDayMat  2 7.0489   0.02947 *

#PT
# Df   X2 P(>|Chi|)   
# AvgDayMat  2 12.2    0.0023 **

#QN
# Df   X2 P(>|Chi|)  
# AvgDayMat  2 5.93     0.052 .

summary(PODFinal)

#CB
# Call:
#   geeglm(formula = PreAbs ~ mSpline(Year, knots = quantile(Year, 
#                                                            probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019)) + 
#            AvgDayMat, family = binomial, data = SiteHourTableB, id = Blocks, 
#          corstr = "ar1")
# 
# Coefficients:
#   Estimate
# (Intercept)                                                                                       0.6553
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))1  -0.9283
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))2  -1.6183
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))3  -3.4238
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))4  -1.0889
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))5  -0.3674
# AvgDayMatADBM1                                                                                   -0.0412
# AvgDayMatADBM2                                                                                    0.3936
# Std.err
# (Intercept)                                                                                      0.4083
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))1  0.7045
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))2  1.1886
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))3  1.7208
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))4  1.1547
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))5  0.4491
# AvgDayMatADBM1                                                                                   0.1753
# AvgDayMatADBM2                                                                                   0.1805
# Wald
# (Intercept)                                                                                     2.58
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))1 1.74
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))2 1.85
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))3 3.96
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))4 0.89
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))5 0.67
# AvgDayMatADBM1                                                                                  0.06
# AvgDayMatADBM2                                                                                  4.75
# Pr(>|W|)
# (Intercept)                                                                                        0.108
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))1    0.188
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))2    0.173
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))3    0.047
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))4    0.346
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))5    0.413
# AvgDayMatADBM1                                                                                     0.814
# AvgDayMatADBM2                                                                                     0.029
# 
# (Intercept)                                                                                      
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))1  
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))2  
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))3 *
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))4  
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))5  
# AvgDayMatADBM1                                                                                   
# AvgDayMatADBM2                                                                                  *
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     1.01  0.0392
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.938 0.00479
# Number of clusters:   80  Maximum cluster size: 603 

#QN
# Call:
#   geeglm(formula = PreAbs ~ AvgDayMat, family = binomial, data = SiteHourTableB, 
#          id = Blocks, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)      -1.586   0.130 149.45   <2e-16 ***
#   AvgDayMatADBM1    0.374   0.229   2.67    0.102    
# AvgDayMatADBM2    0.506   0.252   4.02    0.045 *  
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)        1   0.135
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.859  0.0203
# Number of clusters:   48  Maximum cluster size: 287 

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

# Step 9: Construction of the ROC curve    --------------------------------
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

#Confusion matrix:

#CB
#observed
#predicted     1     0
#1 13494 10046
#0  7008 14806

#PT
#observed
#predicted     1     0
# 1 1124 5334
# 0  608 7385

#QN
# observed
# predicted    1    0
# 1 2028 8016
# 0  292 2924

#BD
#observed
#predicted    1    0
#1 6755 5538
#0 1919 3216

# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

#save workspace
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
save.image(file = fileName)



