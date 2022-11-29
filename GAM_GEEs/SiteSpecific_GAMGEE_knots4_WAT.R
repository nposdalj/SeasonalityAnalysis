### Modelling habitat preference of sperm whales using Generalized Additive Models with Generalized Estimating Equations ###
### Script adapted from Pirotta et al. (2011) and Benjamins ###  
### Example from the GofAK + BSAI ###
### 7 Models total:
        #Site specific models: CB, PT, QN, BD (more than 270 d of recording)
        #Region specific models: BSAI + GOA
        #Big model: all 7 sites
        #Change

# All the libraries have to be installed prior to their utilization (see R help on library installation)

library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
library(splines)         # to construct the B-splines within a GEE-GLM
library(tidyverse)       # because it literally does everything
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2)         # to build the partial residual plots
library(mvtnorm)         # to build the partial residual plots
library(gridExtra)       # to build the partial residual plots
library(SimDesign)
library(lubridate)
library(regclass)
library(mgcv)
library(ChemoSpec2D)
library(car)            # to run an ANOVA
library(splines2)       # to use mSpline for the GEEs
library(ggfortify)      # extract confidence interval for ACF plots
library(lubridate)      # to adjust dates

site = 'JAX' #specify the site of interest 

# Step 1: Load the Data -----------------------------------------------------------
GDir = 'G'
dir = paste(GDir,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites",sep="")
saveDir = paste(GDir,":/My Drive/WAT_TPWS_metadataReduced/Plots/",site, sep="")
saveWorkspace = paste(GDir,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(GDir,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/All_Sites/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin, "%m/%d/%Y")
HourTable$tbin = as.POSIXct(strptime(HourTable$tbin, "%m/%d/%Y %H:%M"))
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(site,Site))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable = subset(SiteHourTable, Site == site)#subset the table for the site only

#If it's a leap year, delete julian day 366
SiteHourTable = SiteHourTable[!(SiteHourTable$Julian==366),]

#Time lost variable
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region (WAT)
# The column Julian represents the day of the year and the column year represents the year of recording.
# The column TimeLost represents the amount of time that we didn't record in each one hour bin


# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
BlockMod<-glm(PreAbs~
               bs(Julian,k=4)+
               TimeLost+
               as.factor(Year)
             ,data=SiteHourTable,family=binomial)

ACF = acf(residuals(BlockMod),lag.max = 10000)
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval = ACFidx[1] #<---- check this value to make sure it's reasonable (<600)

#create the blocks based on the full timeseries
startDate = SiteHourTable$tbin[1]
endDate = SiteHourTable$tbin[nrow(SiteHourTable)]
timeseries = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseries)/ACFval)), times=1, each=ACFval)
divdiff = nrow(timeseries) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff) 
timeseries$block = c(preBlock,lastVec)
names(timeseries)[1] = 'tbin'
SiteHourTableB = left_join(SiteHourTable,timeseries,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$block))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$Blocks = UnBlock$sequence[match(SiteHourTableB$block,UnBlock$`unique(SiteHourTableB$block)`)]
difference = diff(SiteHourTableB$Blocks) #find difference between rows

# Step 3: ANOVA to Check for Significance of Variables --------------------
#ANOVA with car package
Anova(BlockMod) #run this to get p-values
summary(BlockMod)

#BP
# bs(Julian, k = 4)   91.532  4     <2e-16 ***
# TimeLost             1.879  1     0.1705    
# as.factor(Year)     83.482  3     <2e-16 ***
  
#BS
#bs(Julian, k = 4)   321.76  4     <2e-16 ***
#TimeLost              0.06  1      0.802    
#as.factor(Year)     164.22  3     <2e-16 ***

#BC
# bs(Julian, k = 4)  1382.12  4  < 2.2e-16 ***
# TimeLost             27.76  1  1.376e-07 ***
# as.factor(Year)    1003.91  3  < 2.2e-16 ***

# NC
# bs(Julian, k = 4)   743.72  4  < 2.2e-16 ***
# TimeLost             35.22  1  2.938e-09 ***
# as.factor(Year)     557.67  4  < 2.2e-16 ***

# HZ
# bs(Julian, k = 4)  1472.36  4  < 2.2e-16 ***
# TimeLost             10.22  1   0.001388 ** 
# as.factor(Year)     390.93  4  < 2.2e-16 ***

# OC
# bs(Julian, k = 4)   955.14  4  < 2.2e-16 ***
#   TimeLost             13.16  1  0.0002861 ***
#   as.factor(Year)     355.66  4  < 2.2e-16 ***

# WC
# bs(Julian, k = 4)      792  4    < 2e-16 ***
#   TimeLost                14  1    0.00019 ***
#   as.factor(Year)       1296  3    < 2e-16 ***

# GS
# bs(Julian, k = 4)      545  4     <2e-16 ***
#   TimeLost                 5  1      0.031 *  
#   as.factor(Year)        290  3     <2e-16 ***

# JAX
# bs(Julian, k = 4)    118.4  4     <2e-16 ***
# TimeLost               0.1  1       0.76    
# as.factor(Year)      103.2  3     <2e-16 ***

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
  GLM1 = glm(PreAbs ~ bs(Julian) + TimeLost + as.factor(Year), family = binomial, data = SiteHourTableB)
  VIF(GLM1)
  
#BP
#bs(Julian)      1.248195  3        1.037641
#TimeLost        1.003835  1        1.001916
#as.factor(Year) 1.252221  3        1.038198
  
#BS
#bs(Julian)      1.231767  3        1.035352
#TimeLost        1.007290  1        1.003638
#as.factor(Year) 1.238047  3        1.036230
  
#BC
# bs(Julian)      1.428575  3        1.061249
# TimeLost        1.005951  1        1.002971
# as.factor(Year) 1.426685  3        1.061015
  
# NC
# bs(Julian)      1.468314  3        1.066113
# TimeLost        1.003124  1        1.001561
# as.factor(Year) 1.466258  4        1.049002
  
# HZ
# bs(Julian)      1.560773  3        1.077019
# TimeLost        1.002440  1        1.001219
# as.factor(Year) 1.560166  4        1.057174
  
# OC
# bs(Julian)      1.303200  3        1.045126
# TimeLost        1.005726  1        1.002859
# as.factor(Year) 1.298717  4        1.033212
  
# WC
# bs(Julian)      1.434  3           1.062
# TimeLost        1.040  1           1.020
# as.factor(Year) 1.411  3           1.059
  
# GS
# bs(Julian)      1.65  3            1.09
# TimeLost        1.00  1            1.00
# as.factor(Year) 1.65  3            1.09
  
# JAX
# bs(Julian)      2.55  3            1.17
# TimeLost        1.00  1            1.00
# as.factor(Year) 2.55  3            1.17

# Step 5: Model Selection - Covariate/variable Preparation -------------------------

# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=4), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,365,length=4)))$X[,2:3]
AvgDayMat = as.matrix(AvgDayBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.
POD0<-geeglm(PreAbs ~ 1, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 

#The model with the lowest QIC is the final model.
#Model order - the variable, when removed, that increases the QIC the most goes first

  #The initial full model is:
  POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without AvgDayMat
  POD3b = geeglm(PreAbs ~ as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
  #without Year
  POD3c = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

  model3A = c("POD0","POD3a","POD3b","POD3c")
  QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1])
  QICmod3A<-data.frame(rbind(model3A,QIC3A))
  QICmod3A
  
  #BS
  #QIC            QIC.1            QIC.2            QIC.3            
  # model3A             POD0            POD3a            POD3b           POD3c      
  #QIC3A   62743.5048448833 60915.8541569011 61277.9740202374 61977.6718987442
  #Final model includes year and j day
  # model order: year, then avg julian day
  
  #BP
  # QIC            QIC.1            QIC.2            QIC.3
  # model3A             POD0            POD3a            POD3b            POD3c
  # QIC3A   12735.9210928996 12539.6262418463 12615.8897572215 12615.5422645209
  #Final model includes year and j day
  # model order adjusted to match BS, JD  and Year have same values
  
  # BC
  # QIC            QIC.1           QIC.2            QIC.3
  # model3A             POD0            POD3a           POD3b            POD3c
  # QIC3A   33127.7674172908 30806.6812723813 32296.702835807 31803.1418968153
  # Final Model: jday then year
  
  # NC
  # QIC            QIC.1           QIC.2            QIC.3
  # model3A            POD0            POD3a           POD3b            POD3c
  # QIC3A   40824.170090849 40138.8454691975 40552.406139129 40399.5999343466
  #Final model includes year and j day, order: j day then year
  
  # HZ
  # QIC            QIC.1            QIC.2            QIC.3
  # model3A            POD0            POD3a            POD3b            POD3c
  # QIC3A   40374.245194265 39475.8025007485 40765.7570671783 39352.1721116676
  #Final model includes only j day
  
  # OC
  # QIC            QIC.1            QIC.2            QIC.3
  #  model3A             POD0            POD3a            POD3b            POD3c
  #  QIC3A   42055.0086155732 41070.7674227592 41843.4065567286 41181.3824881721
  #Final model includes year and j day, , order: j day then year
  
  # WC
  # QIC            QIC.1            QIC.2            QIC.3
  # model3A             POD0            POD3a            POD3b            POD3c
  # QIC3A   30230.6737362667 27906.1944266229 28659.4816960935 29190.2991358086
  #Final model includes year and j day, order: year then j day
  
  # GS
  # QIC            QIC.1            QIC.2            QIC.3
  # model3A             POD0            POD3a            POD3b            POD3c
  # QIC3A   16364.2045633377 15909.9719423346 16271.4510039233 16103.6826382639
  #Final model includes year and j day, order: j day then year
  
  # JAX
  # QIC            QIC.1            QIC.2            QIC.3
  # model3A             POD0            POD3a            POD3b            POD3c
  # QIC3A   6000.25899926762 5860.23861537633 5953.00008800914 5925.57746064399
  # Final model: jday then year
  

# Step 7: Finalize Model ----------------------------------------
  

#Year as factor
if (site == 'BP'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ as.factor(Year)+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
} 
if (site == 'BS') {
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ as.factor(Year)+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
if (site == 'BC'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
if (site == 'NC'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
if (site == ('HZ')){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
if (site == 'OC'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
if (site == 'WC'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ as.factor(Year)+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
if (site == 'GS'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
if (site == 'JAX'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}
# STEP 8: Interpreting the summary of the model --------------------------
# How to interpret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reported by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero
anova(PODFinal)
#re-run *
#redo any models that have insignificant values  
  # BS 
  # Df   X2 P(>|Chi|)    
  # as.factor(Year)  3 17.1   0.00066 ***
  # AvgDayMat        2 46.0     1e-10 ***
  
  # BP 
  # Df   X2 P(>|Chi|)    
  # as.factor(Year)  3 34.8   1.4e-07 ***
  # AvgDayMat        2 25.1   3.6e-06 ***
  
  # BC
  # AvgDayMat        2  67.5   2.1e-15 ***
  # as.factor(Year)  3 109.2   < 2e-16 ***
  
  # NC 
  # Df    X2 P(>|Chi|)   
  # AvgDayMat        2  8.98    0.0112 * 
  # as.factor(Year)  4 15.52    0.0037 **
  
  # HZ 
  # Df   X2 P(>|Chi|)    
  # AvgDayMat  2 21.9   1.7e-05 ***

  # OC
  # Df     X2 P(>|Chi|)    
  # AvgDayMat        2 28.050 8.108e-07 ***
  #   as.factor(Year)  4 10.987   0.02671 * 
  
  # WC 
  # Df   X2 P(>|Chi|)    
  # as.factor(Year)  3 46.0   5.7e-10 ***
  # AvgDayMat        2 34.7   2.9e-08 ***
  
  # GS
  # Df   X2 P(>|Chi|)    
  # AvgDayMat        2 21.8   1.9e-05 ***
  #   as.factor(Year)  3 25.1   1.4e-05 ***
  
  # JAX
  # AvgDayMat        2 20.5   3.6e-05 ***
  # as.factor(Year)  3 17.3   0.00062 ***
  
#Save model output
filename = paste(saveWorkspace,site,'_SiteSpecificModelSummary.txt',sep="")
sink(filename)
summary(PODFinal)
anova(PODFinal)
 sink(file = NULL)

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
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45? line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45? line and the ROC curve

#Find max and position of max
max = max(d,na.rm=TRUE)
names(d)[1] <- 'Col1'
position = which.max(d$Col1)

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
cutoff = alpha[position,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45? line and the curve

# This value can now be used to build the confusion matrix:

DATA<-matrix(0,dim(SiteHourTableB)[1] ,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim()) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:dim(SiteHourTableB)[1]                                    # the first column is filled with an ID value that is unique for each row
DATA$Observed<-SiteHourTableB$PreAbs                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinal,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")


# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput.RData',sep="")
save.image(file = fileName)
