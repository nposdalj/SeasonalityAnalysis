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

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region (GOA or BSAI).
# The column Julian represents the day of the year and the column year represents the year of recording.

# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
BlockMod<-glm(PreAbs~
               bs(Julian,k=4)+
               bs(Year,k=4)+
               as.factor(Region)+
                TimeLost
             ,data=SiteHourTable,family=binomial)

acf(residuals(BlockMod),lag.max = 1000, ylim=c(-0.1,0.1))
acf(residuals(BlockMod),lag.max = 1000, ylim=c(-0.1,0.1), xlim=c(20,30))
ACFval = 29

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
#Quick ANOVA to check for significance of variables - using the car packageAnova(BlockMod)
Anova(BlockMod)

# LR Chisq Df Pr(>Chisq)    
# bs(Julian, k = 4)  1520.96  4  < 2.2e-16 ***
#   bs(Year, k = 4)    1591.85  3  < 2.2e-16 ***
#   as.factor(Region)   374.94  1  < 2.2e-16 ***
#   TimeLost             36.56  1  1.484e-09 ***

summary(BlockMod)
# Call:
#   glm(formula = PreAbs ~ bs(Julian, k = 4) + bs(Year, k = 4) + 
#         as.factor(Region) + TimeLost, family = binomial, data = SiteHourTable)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.6861  -0.9395  -0.7666   1.2766   1.9180  
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)           0.71954    0.13535   5.316 1.06e-07 ***
#   bs(Julian, k = 4)1   -0.81665    0.13630  -5.991 2.08e-09 ***
#   bs(Julian, k = 4)2   -0.45488    0.14346  -3.171  0.00152 ** 
#   bs(Julian, k = 4)3    0.65089    0.14460   4.501 6.76e-06 ***
#   bs(Julian, k = 4)4   -0.13094    0.13398  -0.977  0.32842    
# bs(Year, k = 4)1     -0.40425    0.04686  -8.628  < 2e-16 ***
#   bs(Year, k = 4)2     -0.84662    0.04853 -17.447  < 2e-16 ***
#   bs(Year, k = 4)3     -1.80590    0.06748 -26.761  < 2e-16 ***
#   bs(Year, k = 4)4           NA         NA      NA       NA    
# as.factor(Region)GOA -0.51772    0.02674 -19.361  < 2e-16 ***
#   TimeLost              0.07999    0.01326   6.031 1.63e-09 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 130202  on 99305  degrees of freedom
# Residual deviance: 124539  on 99296  degrees of freedom
# AIC: 124559
# 
# Number of Fisher Scoring iterations: 4

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
GLM1 = glm(PreAbs ~ Julian + TimeLost + Year + as.factor(Region), family = binomial, data = SiteHourTableB)
#VIF scores in GLM to work out collinearity:
VIF(GLM1)
# Julian          TimeLost              Year as.factor(Region) 
# 1.068787          1.004482          1.594076          1.552767 

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
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
#QIC0A   130220.476099582 128868.381114667 128832.835812775
#Julian day as a covariance matrix

#Year
POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1c = geeglm(PreAbs ~ mSpline(Year,
                                knots=quantile(Year, probs=c(0.333,0.666)),
                                Boundary.knots=c(2010,2019)), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model1A<-c("POD0", "POD1a", "POD1b","POD1c")
QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
QICmod1A<-data.frame(rbind(model1A,QIC1A))
QICmod1A
#QIC            QIC.1            QIC.2            QIC.3
#model1A             POD0            POD1a            POD1b            POD1c
#QIC1A   130220.476099582 124065.478597119 130127.05877415 124216.392022745
#Year as mSpline

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
#QIC2A   130220.476099582 130382.415066292 130250.286587733 130235.333647841
#TimeLost as linear.

#Region is a factor, no need to check them.

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The Initial full model:
POD3aa = geeglm(PreAbs ~ AvgDayMat+mSpline(Year,
                                           knots=quantile(Year, probs=c(0.333,0.666)),
                                           Boundary.knots=c(2010,2019))+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3bb = geeglm(PreAbs ~ mSpline(Year,
                                 knots=quantile(Year, probs=c(0.333,0.666)),
                                 Boundary.knots=c(2010,2019))+TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3cc = geeglm(PreAbs ~ AvgDayMat +TimeLost+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3dd = geeglm(PreAbs ~ AvgDayMat+mSpline(Year,
                                           knots=quantile(Year, probs=c(0.333,0.666)),
                                           Boundary.knots=c(2010,2019))+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Region
POD3ee = geeglm(PreAbs ~ AvgDayMat+mSpline(Year,
                                           knots=quantile(Year, probs=c(0.333,0.666)),
                                           Boundary.knots=c(2010,2019))+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3AA = c("POD0","POD3aa","POD3bb","POD3cc","POD3dd","POD3ee")
QIC3AA = c(QIC(POD0)[1],QIC(POD3aa)[1],QIC(POD3bb)[1],QIC(POD3cc)[1],QIC(POD3dd)[1],QIC(POD3ee)[1][1])
QICmod3AA<-data.frame(rbind(model3AA,QIC3AA))
QICmod3AA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3AA             POD0           POD3aa           POD3bb           POD3cc           POD3dd           POD3ee
#                      QIC            QIC.1            QIC.2            QIC.3           QIC.4
# model3AA             POD0           POD3aa           POD3bb           POD3cc          POD3dd
# QIC3AA   130220.476099582 123505.317964113 124458.710476665 128400.600209527 123486.83952005
# QIC.5
# model3AA          POD3ee
# QIC3AA   123342.43807106
#Remove TimeLost.

#The  full model without TimeLost is:
POD3ff = geeglm(PreAbs ~ AvgDayMat +as.factor(Region)+mSpline(Year,
                                                                knots=quantile(Year, probs=c(0.333,0.666)),
                                                                Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3gg = geeglm(PreAbs ~ as.factor(Region)+mSpline(Year,
                                                     knots=quantile(Year, probs=c(0.333,0.666)),
                                                     Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Region
POD3hh = geeglm(PreAbs ~ AvgDayMat + mSpline(Year,
                                             knots=quantile(Year, probs=c(0.333,0.666)),
                                             Boundary.knots=c(2010,2019)),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3ii = geeglm(PreAbs ~ AvgDayMat + as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3BB = c("POD0","POD3ff","POD3gg",'POD3hh',"POD3ii")
QIC3BB = c(QIC(POD0)[1],QIC(POD3ff)[1],QIC(POD3gg)[1],QIC(POD3hh)[1],QIC(POD3ii)[1])
QICmod3BB<-data.frame(rbind(model3BB,QIC3BB))
QICmod3BB
#                      QIC           QIC.1           QIC.2            QIC.3            QIC.4
# model3BB             POD0          POD3ff          POD3gg           POD3hh           POD3ii
# QIC3BB   130220.476099582 123486.83952005 124437.84094239 123326.698228703 128387.081896753
#Decided to still include Region since I'm curious if the two regions are different or not

#In descending order:
#mSpline(Year)
#AvgDayMat
#as.factor(Region)

# Step 7: Finalize Model --------------------------------------------------
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinal = geeglm(PreAbs ~ mSpline(Year,
                                   knots=quantile(Year, probs=c(0.333,0.666)),
                                   Boundary.knots=c(2010,2019))+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
PODFinal_Region = geeglm(PreAbs ~ mSpline(Year,
                                   knots=quantile(Year, probs=c(0.333,0.666)),
                                   Boundary.knots=c(2010,2019))+AvgDayMat+as.factor(Region),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero
anova(PODFinal)
# Df     X2 P(>|Chi|)    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 475.70 < 2.2e-16 ***
#   AvgDayMat                                                                                       2  91.76 < 2.2e-16 ***

anova(PODFinal_Region)

# Df  X2 P(>|Chi|)    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))  5 476    <2e-16 ***
#   AvgDayMat                                                                                       2  92    <2e-16 ***
#   as.factor(Region)                                                                               1   2      0.13    

summary(PODFinal)
# Call:
#   geeglm(formula = PreAbs ~ mSpline(Year, knots = quantile(Year, 
#                                                            probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019)) + 
#            AvgDayMat + as.factor(Region), family = binomial, data = SiteHourTableB, 
#          id = Blocks, corstr = "ar1")
# 
# Coefficients:
#   Estimate Std.err   Wald Pr(>|W|)    
# (Intercept)                                                                                      -0.4338  0.0857  25.64  4.1e-07 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))1   2.0273  0.1962 106.82  < 2e-16 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))2  -2.1051  0.3885  29.37  6.0e-08 ***
#   mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))3  -0.1236  0.5120   0.06  0.80926    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))4  -0.1449  0.3251   0.20  0.65575    
# mSpline(Year, knots = quantile(Year, probs = c(0.333, 0.666)), Boundary.knots = c(2010, 2019))5   0.7064  0.2037  12.03  0.00052 ***
#   AvgDayMatADBM1                                                                                   -0.1414  0.0529   7.15  0.00752 ** 
#   AvgDayMatADBM2                                                                                    0.3739  0.0461  65.89  4.4e-16 ***
#   as.factor(Region)GOA                                                                             -0.1669  0.1106   2.28  0.13132    
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation structure = ar1 
# Estimated Scale Parameters:
#   
#   Estimate Std.err
# (Intercept)    0.987  0.0104
# Link = identity 
# 
# Estimated Correlation Parameters:
#   Estimate Std.err
# alpha    0.803 0.00477
# Number of clusters:   2050  Maximum cluster size: 90 

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
cmx(DATA, threshold = cutoff)      
# observed
# predicted     1     0
# 1 22189 22623
# 0 13936 40558
  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,'/BigModel_gamgeeOutput.RData',sep="")
save.image(file = fileName)