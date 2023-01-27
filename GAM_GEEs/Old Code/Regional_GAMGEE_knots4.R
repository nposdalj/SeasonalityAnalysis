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
GDrive= 'I'
dir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel", sep="")
saveDir = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/Plots/",region, sep="")
saveWorkspace = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/",region,'/',sep="")
fileName = paste(GDrive,":/My Drive/GofAK_TPWS_metadataReduced/SeasonalityAnalysis/BigModel/AllSitesGrouped_Binary_GAMGEE_ROW.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(region,Region))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12

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
if (region == 'BSAI'){
BlockMod<-glm(PreAbs~
                bs(Julian, k=4)+
                TimeLost + 
                as.factor(Site), data=SiteHourTable,family=binomial)
}else{

BlockMod<-glm(PreAbs~ bs(Julian, k = 4)+
                as.factor(Year)+
                as.factor(Site)+
                TimeLost, data=SiteHourTable,family=binomial)
}

ACF = acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1))
CI = ggfortify:::confint.acf(ACF)
ACFidx = which(ACF[["acf"]] < CI, arr.ind=TRUE)
ACFval = ACFidx[1]

#ACFval calculated from averaging all of the GOA sites (PT, QN, CB)
ACFval = 390

#create the blocks based on the full timesereies
startDate = SiteHourTable$tbin[1]
endDate = SiteHourTable$tbin[nrow(SiteHourTable)]
timeseries = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseries)/ACFval)), times=1, each=ACFval)
divdiff = nrow(timeseries) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseries$block = c(preBlock,lastVec)
names(timeseries)[names(timeseries) == 'date'] = 'tbin'
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
# bs(Julian, k = 4)   719.87  4     <2e-16 ***
#   TimeLost              0.62  1     0.4316    
# as.factor(Site)      71.62  1     <2e-16 ***

#GOA
# bs(Julian, k = 4)     1579  4     <2e-16 ***
#   as.factor(Year)       1122  7     <2e-16 ***
#   as.factor(Site)       5806  4     <2e-16 ***
#   TimeLost                 5  1      0.031 * 

# Step 4: Data Exploration and Initial Analysis ---------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:
if (region == 'BSAI'){
  GLM = glm(PreAbs~bs(Julian)+TimeLost+as.factor(Site),family=binomial,data=SiteHourTableB)
}else{
  GLM = glm(PreAbs~bs(Julian)+TimeLost+as.factor(Site)+as.factor(Year),family=binomial,data=SiteHourTableB)
}
VIF(GLM)

#BSAI
# bs(Julian)      1.133165  3        1.021054
# TimeLost        1.001865  1        1.000932
# as.factor(Site) 1.134216  1        1.064996

#GOA
# bs(Julian)      1.597  3           1.081
# TimeLost        1.010  1           1.005
# as.factor(Site) 3.111  4           1.152
# as.factor(Year) 4.594  7           1.115

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
# QIC0A   26417.3656260167 25860.248483494 25856.1463739235
#Julian day as a covariance matrix

#GOA
#QIC            QIC.1            QIC.2
#model0A             POD0            POD0a            POD0b
# QIC0A   107873.224599673 106330.17458963 106321.525271609
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
# QIC1A   107873.224599673 103022.472613358 107662.928732495 103058.77687551
#Year as a factor.

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
# QIC2A   26417.3656260167 26458.5347150241 26426.0172751744 26462.5544299619
#TimeLost as linear

#GOA
#QIC            QIC.1            QIC.2            QIC.3
#model2A             POD0            POD2a            POD2b            POD2c
# QIC2A   107873.224599673 107869.834236126 107854.824216274 107803.257314075
#TimeLost as mSpline, but I'm going to use linear

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
# QIC3A   26417.3656260167 25753.8401404939 26417.8686943176 25863.0210955829 25746.0005325001
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
# QIC3B   26417.3656260167 25746.0005325001 26408.8121046668 25856.1463739235
#model3B             POD0            POD3e            POD3f            POD3g
#The full model has the lowest QIC. AvgDayMat then Site.
}else{
#The initial full model is:
POD3a = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+TimeLost+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ as.factor(Site)+TimeLost+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3c = geeglm(PreAbs ~ AvgDayMat +TimeLost+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3d = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+as.factor(Year),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#Without Year
POD3e = geeglm(PreAbs ~ AvgDayMat+as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c","POD3d","POD3e")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1],QIC(POD3e)[1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4           QIC.5
#model3A             POD0            POD3a            POD3b            POD3c            POD3d           POD3e
#QIC3A   107873.224599673 2268321.88261158 2294195.98522917 102233.183647118 2014807.0576355 101815.271827164
#Remove Year

#Second round of model testing without Time Lost
#The initial full model is:
POD3f = geeglm(PreAbs ~ AvgDayMat +as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3g = geeglm(PreAbs ~ as.factor(Site)+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3h = geeglm(PreAbs ~ AvgDayMat+TimeLost,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3i = POD3c = geeglm(PreAbs ~ AvgDayMat+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3B = c("POD0","POD3f","POD3g","POD3h","POD3i")
QIC3B = c(QIC(POD0)[1],QIC(POD3f)[1],QIC(POD3g)[1],QIC(POD3h)[1],QIC(POD3i)[1])
QICmod3B<-data.frame(rbind(model3B,QIC3B))
QICmod3B
#GOA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4
#model3B             POD0            POD3f            POD3g            POD3h            POD3i
#QIC3B   107873.224599673 101815.271827164 100118.804637723 106304.976135675 101807.980067339
#Remove TimeLost

#Third round of model testing without TimeLost
POD3i = geeglm(PreAbs ~ AvgDayMat+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3j = geeglm(PreAbs ~ as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3k = POD3c = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3C = c("POD0","POD3i","POD3j","POD3k")
QIC3C = c(QIC(POD0)[1],QIC(POD3i)[1],QIC(POD3j)[1],QIC(POD3k)[1])
QICmod3C<-data.frame(rbind(model3C,QIC3C))
QICmod3C
#GOA
# QIC            QIC.1            QIC.2            QIC.3
# model3C             POD0            POD3i            POD3j            POD3k
#QIC3C   107873.224599673 101807.980067339 100078.781796158 106321.525271609
#Full model is the best. Site goes first, then AvgDayMat.
}

# Step 7: Finalize Model --------------------------------------------------
if (region == 'BSAI'){
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ AvgDayMat + as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
} else {
  dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinal = geeglm(PreAbs ~ as.factor(Site)+AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
}

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinal)

#BSAI
#AvgDayMat        2 13.2    0.0014 ** 
#as.factor(Site)  1 16.9     4e-05 ***

#GOA
#as.factor(Site)  4 587   < 2e-16 ***
#AvgDayMat        2  40   1.7e-09 ***

filename = paste(saveWorkspace,region,'_RegionalModelSummary.txt',sep="")
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
# 1 15604 31995
# 0  5515 32462

  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,region,'_RegionSpecific_gamgeeOutput.RData',sep="")
save.image(file = fileName)