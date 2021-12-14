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

region = 'GOA' #specify the region of interest

# Step 1: Load the Data -----------------------------------------------------------
dir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites")
saveDir = paste("I:/My Drive/GofAK_TPWS_metadataReduced/Plots/",region, sep="")
saveWorkspace = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
fileName = paste("I:/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(region,Region))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region (GOA or BSAI).
# The column Julian represents the day of the year and the column year represents the year of recording.

# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
if (region == 'BSAI'){
BlockMod<-glm(PreAbs~
                bs(Julian, k=4)+
                TimeLost + 
                as.factor(Site), data=SiteHourTable,family=binomial)
}else{

BlockMod<-glm(PreAbs~
                bs(Julian)+
                bs(Year) + 
                TimeLost, data=SiteHourTable,family=binomial)
}

acf(residuals(BlockMod), lag.max = 700, ylim=c(0,0.1))
acf(residuals(BlockMod), lag.max = 1000, ylim=c(-0.1,0.1), xlim =c(0,30)) 

if (region == 'BSAI'){
  ACFval = 653
}else{
  ACFval = 23
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
#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockMod)

#BSAI
# bs(Julian, k = 4)      529  4    < 2e-16 ***
#   TimeLost                 2  1       0.21    
# as.factor(Site)         40  1    2.4e-10 ***

#GOA
# bs(Julian, k = 4)     1763  4     <2e-16 ***
#   TimeLost                 4  1      0.043 *  
#   bs(Year, k = 4)       1003  3     <2e-16 ***
#   as.factor(Site)       6317  4     <2e-16 ***

summary(BlockMod)

#BSAI
# Call:
#   glm(formula = PreAbs ~ bs(Julian, k = 4) + TimeLost + as.factor(Site), 
#       family = binomial, data = SiteHourTable)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -3.35   -1.11   -0.99    1.18    1.45  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)          5.6003     1.1780    4.75  2.0e-06 ***
#   bs(Julian, k = 4)1  -5.4922     1.1870   -4.63  3.7e-06 ***
#   bs(Julian, k = 4)2  -7.4281     1.1739   -6.33  2.5e-10 ***
#   bs(Julian, k = 4)3  -3.3485     1.1919   -2.81    0.005 ** 
#   bs(Julian, k = 4)4  -6.2190     1.1778   -5.28  1.3e-07 ***
#   TimeLost             0.0667     0.0535    1.25    0.213    
# as.factor(Site)KS   -0.3360     0.0532   -6.32  2.7e-10 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 26591  on 19183  degrees of freedom
# Residual deviance: 26047  on 19177  degrees of freedom
# AIC: 26061
# 
# Number of Fisher Scoring iterations: 6

#GOA
# Call:
#   glm(formula = PreAbs ~ bs(Julian, k = 4) + TimeLost + bs(Year, 
#                                                            k = 4) + as.factor(Site), family = binomial, data = SiteHourTable)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.540  -0.932  -0.597   1.131   2.558  
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -1.7815     0.1900   -9.38  < 2e-16 ***
#   bs(Julian, k = 4)1  -0.5607     0.1872   -3.00   0.0027 ** 
#   bs(Julian, k = 4)2   1.1451     0.1943    5.89  3.8e-09 ***
#   bs(Julian, k = 4)3   1.0088     0.1952    5.17  2.4e-07 ***
#   bs(Julian, k = 4)4   0.7805     0.1838    4.25  2.2e-05 ***
#   TimeLost             0.0280     0.0139    2.02   0.0432 *  
#   bs(Year, k = 4)1     0.0633     0.0462    1.37   0.1707    
# bs(Year, k = 4)2    -1.5632     0.0565  -27.65  < 2e-16 ***
#   bs(Year, k = 4)3    -0.4976     0.0895   -5.56  2.7e-08 ***
#   bs(Year, k = 4)4         NA         NA      NA       NA    
# as.factor(Site)CB    1.4117     0.0476   29.64  < 2e-16 ***
#   as.factor(Site)KOA   0.5053     0.0663    7.62  2.4e-14 ***
#   as.factor(Site)PT   -0.1722     0.0560   -3.07   0.0021 ** 
#   as.factor(Site)QN    0.1612     0.0514    3.14   0.0017 ** 
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 101930  on 80121  degrees of freedom
# Residual deviance:  90242  on 80109  degrees of freedom
# AIC: 90268
# 
# Number of Fisher Scoring iterations: 4

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
if (region == 'BSAI'){
  GLM = glm(PreAbs~Julian+TimeLost+Site,family=binomial,data=SiteHourTableB)
}else{
  GLM = glm(PreAbs~Julian+TimeLost+Site+Year,family=binomial,data=SiteHourTableB)
}
VIF(GLM)

#BSAI
# Julian TimeLost     Site 
# 1        1        1 

#GOA
# GVIF Df GVIF^(1/(2*Df))
# Julian   1.08  1            1.04
# TimeLost 1.01  1            1.00
# Site     1.35  4            1.04
# Year     1.43  1            1.20

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
#BSAI
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
# QIC0A   26598.7716063937 26231.2071664586 26223.7668717875
#Julian day as a covariance matrix

#GOA
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
# QIC0A   101973.193288013 100785.805498286 100764.708449876
#Julian day as a covariance matrix

if (region == 'GOA'){
  #Year
  POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  POD1c = geeglm(PreAbs ~ mSpline(Year,
                                  knots=quantile(Year, probs=c(0.333,0.666)),
                                  Boundary.knots=c(2011,2019)), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  model1A<-c("POD0", "POD1a", "POD1b","POD1c")
  QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
  QICmod1A<-data.frame(rbind(model1A,QIC1A))
  QICmod1A
}
#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model1A             POD0            POD1a            POD1b            POD1c
#QIC1A   101973.193288013 97401.7545591589 101737.310584738 97528.4566324262
#Year as an mSpline.

#TimeLost
POD2a = geeglm(PreAbs ~ as.factor(TimeLost), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2b = geeglm(PreAbs ~ TimeLost, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD2c = geeglm(PreAbs ~ mSpline(TimeLost,
                                knots=4), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model2A<-c("POD0", "POD2a", "POD2b","POD2c")
QIC2A<-c(QIC(POD0)[1],QIC(POD2a)[1],QIC(POD2b)[1],QIC(POD2c)[1])
QICmod2A<-data.frame(rbind(model2A,QIC2A))
QICmod2A
#BSAI
#QIC            QIC.1            QIC.2          QIC.3
#model2A             POD0            POD2a            POD2b          POD2c
# QIC2A   26598.7716063937 26636.0005421372 26610.4495685621 26654.5141226769
#TimeLost as linear

#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
#QIC2A   101973.193288013 102049.559092936 101981.723519973 101955.841619011
#TimeLost as linear

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
if (region == 'BSAI'){
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
# QIC3A   26598.7716063937 26173.6683253399 26599.4583801071 26233.0799088394 26162.9487085645
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
# QIC3B   26598.7716063937 26162.9487085645 26587.3271097547 26223.7668717875
#model3B             POD0            POD3e            POD3f            POD3g
#The full model has the lowest QIC. AvgDayMat then Site.
}else{
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+TimeLost+mSpline(Year,
                                                                   knots=quantile(Year, probs=c(0.333,0.666)),
                                                                   Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ as.factor(Site)+TimeLost+mSpline(Year,
                                                         knots=quantile(Year, probs=c(0.333,0.666)),
                                                         Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost+mSpline(Year,
                                                    knots=quantile(Year, probs=c(0.333,0.666)),
                                                    Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3d = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+mSpline(Year,
                                                          knots=quantile(Year, probs=c(0.333,0.666)),
                                                          Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#Without Year
POD3e = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c","POD3d","POD3e")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1],QIC(POD3e)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4           QIC.5
#model3A             POD0            POD3a            POD3b            POD3c            POD3d           POD3e
#QIC3A   101973.193288013 2014594.95955435 2530294.12534769 97106.7063896764 2032709.30182762 2609370.35408955
#Remove Site.

#Second round of model testing without Site
#The initial full model is:
POD3f = geeglm(PreAbs ~ AvgDayMat +TimeLost+mSpline(Year,
                                                    knots=quantile(Year, probs=c(0.333,0.666)),
                                                    Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3g = geeglm(PreAbs ~ TimeLost+mSpline(Year,
                                         knots=quantile(Year, probs=c(0.333,0.666)),
                                         Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3h = geeglm(PreAbs ~ AvgDayMat + mSpline(Year,
                                          knots=quantile(Year, probs=c(0.333,0.666)),
                                          Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3i = POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3B = c("POD0","POD3f","POD3g","POD3h","POD3i")
QIC3B = c(QIC(POD0)[1],QIC(POD3f)[1],QIC(POD3g)[1],QIC(POD3h)[1],QIC(POD3i)[1])
QICmod3B<-data.frame(rbind(model3B,QIC3B))
QICmod3B
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3B             POD0            POD3f            POD3g            POD3h            POD3i
#QIC3B   101973.193288013 97106.7063896764 97538.8701419481 97098.4175104023 100771.303533705
#Remove TimeLost

#Third round of model testing without TimeLost
POD3i = geeglm(PreAbs ~ AvgDayMat +mSpline(Year,
                                                    knots=quantile(Year, probs=c(0.333,0.666)),
                                                    Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3j = geeglm(PreAbs ~ mSpline(Year,
                                         knots=quantile(Year, probs=c(0.333,0.666)),
                                         Boundary.knots=c(2011,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3k = POD3c = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3C = c("POD0","POD3i","POD3j","POD3k")
QIC3C = c(QIC(POD0)[1],QIC(POD3i)[1],QIC(POD3j)[1],QIC(POD3k)[1])
QICmod3C<-data.frame(rbind(model3C,QIC3C))
QICmod3C
#GOA
# QIC            QIC.1            QIC.2            QIC.3
# model3C             POD0            POD3i            POD3j            POD3k
# QIC3C   101973.193288013 97098.4175104023 97528.4566324262 100764.708449876
#Full model is the best. Year goes first, then AvgDayMat.
}

# Step 7: Finalize Model --------------------------------------------------
if (region == 'BSAI'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat + as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
} else {
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ mSpline(Year,
                                               knots=quantile(Year, probs=c(0.333,0.666)),
                                               Boundary.knots=c(2011,2019))+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinal)

#BSAI
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: PreAbs
#Terms added sequentially (first to last)

#Df     X2 P(>|Chi|)    
# AvgDayMat        2 6.52     0.038 *
#   as.factor(Site)  1 5.87     0.015 *
  
#GOA
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: PreAbs
#Terms added sequentially (first to last)

#Df     X2 P(>|Chi|)    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))  5 462   < 2e-16 ***
#   AvgDayMat                                                                                       2  50   1.2e-11 ***

summary(PODFinal)

#BSAI
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
# Call:
#   geeglm(formula = PreAbs ~ mSpline(Year, knots = quantile(Year, 
#                                                            probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019)) + 
#            AvgDayMat, family = binomial, data = SiteHourTableB, id = Blocks, 
#          corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)                                                                                       0.4957  0.1114  19.78  8.7e-06 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))1  -0.5622  0.1895   8.80  0.00301 ** 
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))2  -4.5773  0.3132 213.57  < 2e-16 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))3  -1.3814  0.4688   8.68  0.00321 ** 
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))4  -1.9020  0.2864  44.11  3.1e-11 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2011, 2019))5  -0.5177  0.1349  14.73  0.00012 ***
#   AvgDayMatADBM1                                                                                   -0.0679  0.0550   1.53  0.21657    
# AvgDayMatADBM2                                                                                    0.3190  0.0453  49.48  2.0e-12 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)     0.99  0.0126
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.792 0.00527
# Number of clusters:   2099  Maximum cluster size: 72 

# Step 9: Construction of the ROC curve    --------------------------------
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
# 1 14556 16541
# 0 12103 36922

  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput.RData',sep="")
save.image(file = fileName)