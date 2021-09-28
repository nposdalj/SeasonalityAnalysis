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
#NOT ON MODEL RESIDUALS
#Day Table
#acf on the 1 day binned data
acf(SiteDayTable$PreAbs, lag.max = 50)
#AllSites = 11

#Hour Table
#acf on the 1 hour binned data
acf(SiteHourTable$PreAbs, lag.max = 2000)
acf(SiteHourTable$PreAbs, lag.max = 2000, ylim=c(0,0.1), xlim =c(400,500)) 
#AllSites = 457

#ON MODEL RESIDUALS
BlockMod<-glm(PreAbs~
               bs(Julian)+
               TimeLost+
               bs(Year)+
                as.factor(Site),
                #as.factor(Region)
             ,data=SiteHourTable,family=binomial)

#Quick ANOVA to check for significance of variables - using the car package
Anova(BlockMod)
#bs(Julian)            1901  3     <2e-16 ***
  #TimeLost                 8  1     0.0041 ** 
  #bs(Year)               656  3     <2e-16 ***
  #as.factor(Site)       7336  5     <2e-16 ***
  #as.factor(Region)           0               

#Without region
#bs(Julian)          1901  3     <2e-16 ***
  #TimeLost               8  1     0.0041 ** 
  #bs(Year)             656  3     <2e-16 ***
  #as.factor(Site)     7713  6     <2e-16 ***
  
summary(BlockMod)
#Call:
  #glm(formula = PreAbs ~ bs(Julian) + TimeLost + bs(Year) + as.factor(Site) + 
        #as.factor(Region), family = binomial, data = SiteHourTable)

#Deviance Residuals: 
  #Min      1Q  Median      3Q     Max  
#-1.698  -1.014  -0.587   1.155   2.410  

#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)           -2.2011     0.0717  -30.70  < 2e-16 ***
  #bs(Julian)1            0.5414     0.0904    5.99  2.2e-09 ***
  #bs(Julian)2            1.7347     0.0543   31.98  < 2e-16 ***
  #bs(Julian)3            0.7701     0.0495   15.55  < 2e-16 ***
  #TimeLost               0.0386     0.0135    2.87  0.00414 ** 
  #bs(Year)1              0.3129     0.0938    3.33  0.00086 ***
  #bs(Year)2             -1.2188     0.0666  -18.31  < 2e-16 ***
  #bs(Year)3              0.6676     0.0571   11.69  < 2e-16 ***
  #as.factor(Site)BD      1.4087     0.0546   25.80  < 2e-16 ***
  #as.factor(Site)CB      1.2277     0.0472   26.02  < 2e-16 ***
  #as.factor(Site)KOA     0.1208     0.0651    1.85  0.06367 .  
#as.factor(Site)KS      1.0348     0.0787   13.14  < 2e-16 ***
  #as.factor(Site)PT     -0.4656     0.0549   -8.48  < 2e-16 ***
  #as.factor(Site)QN     -0.0549     0.0510   -1.08  0.28164    
#as.factor(Region)GOA       NA         NA      NA       NA    

#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 130202  on 99305  degrees of freedom
#Residual deviance: 117236  on 99292  degrees of freedom
#AIC: 117264

#Number of Fisher Scoring iterations: 4

acf(residuals(BlockMod), lag.max = 1000, ylim=c(0,0.1))
acf(residuals(BlockMod), lag.max = 1000, ylim=c(-0.1,0.1), xlim =c(0,50)) 
ACFval = 18

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
#Julian          1.086  1           1.042
#TimeLost        1.007  1           1.003
#Year            2.155  1           1.468
#as.factor(Site) 2.103  6           1.064

## STEP 4: Model selection - covariate preparation ##
# Construct variance-covariance matrices for cyclic covariates:
AvgDayBasis <- gam(PreAbs~s(Julian, bs ="cc", k=-1), fit=F, data = SiteHourTableB, family =binomial, knots = list(HOUR=seq(1,366,length=6)))$X[,2:5]
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
#QIC           QIC.1            QIC.2
#model0A             POD0           POD0a            POD0b
#QIC0A   130217.029932687 128114.47179077 127396.454394609
#Julian day as a covariance matrix

#Year
POD1a = geeglm(PreAbs ~ as.factor(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1b = geeglm(PreAbs ~ Year, family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
POD1c = geeglm(PreAbs ~ bs(Year), family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model1A<-c("POD0", "POD1a", "POD1b","POD1c")
QIC1A<-c(QIC(POD0)[1],QIC(POD1a)[1],QIC(POD1b)[1],QIC(POD1c)[1])
QICmod1A<-data.frame(rbind(model1A,QIC1A))
QICmod1A
#QIC            QIC.1            QIC.2            QIC.3
#model1A             POD0            POD1a            POD1b            POD1c
#QIC1A   130217.029932687 124060.408660561 130120.757199256 126225.320432525
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
#QIC2A   130217.029932687 130325.965379635 130227.908717307 130207.287074268
#TimeLost as a smooth

#Site and region are factors, no need to check them.

#Make Region a binary variable
#SiteHourTableB$RegionBinary = grepl("BSAI", SiteHourTableB$Region)
#SiteHourTableB$RegionBinary = as.integer(SiteHourTableB$RegionBinary)

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The initial full model (with year as a factor) is:
POD3aa = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3bb = geeglm(PreAbs ~ as.factor(Year)+bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3cc = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3dd = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3ee = geeglm(PreAbs ~ AvgDayMat+as.factor(Year)+bs(TimeLost),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3AA = c("POD0","POD3aa","POD3bb","POD3cc","POD3dd","POD3ee")
QIC3AA = c(QIC(POD0)[1],QIC(POD3aa)[1],QIC(POD3bb)[1],QIC(POD3cc)[1],QIC(POD3dd)[1],QIC(POD3ee)[1][1])
QICmod3AA<-data.frame(rbind(model3AA,QIC3AA))
QICmod3AA
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3AA             POD0           POD3aa           POD3bb           POD3cc           POD3dd           POD3ee
#QIC3AA   130217.029932687 2629570.87540136 2654690.83506322 119153.902424341 2708976.49373487 121589.401307196
#Remove year

#The  full model without Year is:
POD3ff = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3gg = geeglm(PreAbs ~ bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3hh = geeglm(PreAbs ~ AvgDayMat + as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3ii = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3BB = c("POD0","POD3ff","POD3gg",'POD3hh',"POD3ii")
QIC3BB = c(QIC(POD0)[1],QIC(POD3ff)[1],QIC(POD3gg)[1],QIC(POD3hh)[1],QIC(POD3ii)[1])
QICmod3BB<-data.frame(rbind(model3BB,QIC3BB))
QICmod3BB
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3BB             POD0           POD3ff           POD3gg           POD3hh          POD3ii
#QIC3BB   130217.029932687 119153.902424341 121591.489964346 119083.872722853 127406.07596043
#Remove timelost.

#The full model without time lost is:
POD3jj = geeglm(PreAbs ~ AvgDayMat + as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3kk = geeglm(PreAbs ~ as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#Without Site
POD3ll = geeglm(PreAbs ~ AvgDayMat ,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3CC = c("POD0","POD3jj","POD3kk",'POD3ll')
QIC3CC = c(QIC(POD0)[1],QIC(POD3jj)[1],QIC(POD3kk)[1],QIC(POD3ll)[1])
QICmod3CC<-data.frame(rbind(model3CC,QIC3CC))
QICmod3CC
#QIC            QIC.1            QIC.2            QIC.3
#model3CC             POD0           POD3jj           POD3kk           POD3ll
#QIC3CC   130217.029932687 119083.872722853 121548.427290765 127396.454394609
#full model is the best.

## STEP 5: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#The initial full model (with year as a smooth) is:
POD3a = geeglm(PreAbs ~ AvgDayMat+bs(Year)+bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3b = geeglm(PreAbs ~ bs(Year)+bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Year
POD3c = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Timelost
POD3d = geeglm(PreAbs ~ AvgDayMat+bs(Year)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3e = geeglm(PreAbs ~ AvgDayMat+bs(Year)+bs(TimeLost),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3A = c("POD0","POD3a","POD3b","POD3c","POD3d","POD3e")
QIC3A = c(QIC(POD0)[1],QIC(POD3a)[1],QIC(POD3b)[1],QIC(POD3c)[1],QIC(POD3d)[1],QIC(POD3e)[1][1])
QICmod3A<-data.frame(rbind(model3A,QIC3A))
QICmod3A
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model3A             POD0            POD3a            POD3b            POD3c            POD3d            POD3e
#QIC3A   130217.029932687 2912982.78637363 2556654.14552418 119153.902424341 2604791.84918335 124046.040537498
#Model without year has the lowest QIC.

#The  full model without Year is:
POD3f = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3g = geeglm(PreAbs ~ bs(TimeLost)+as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without TimeLost
POD3h = geeglm(PreAbs ~ AvgDayMat + as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without Site
POD3i = geeglm(PreAbs ~ AvgDayMat +bs(TimeLost),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3B = c("POD0","POD3f","POD3g",'POD3h',"POD3i")
QIC3B = c(QIC(POD0)[1],QIC(POD3f)[1],QIC(POD3g)[1],QIC(POD3h)[1],QIC(POD3i)[1])
QICmod3B<-data.frame(rbind(model3B,QIC3B))
QICmod3B
#QIC            QIC.1            QIC.2            QIC.3           QIC.4
#model3B             POD0            POD3f            POD3g            POD3h           POD3i
#QIC3B   130217.029932687 119153.902424341 121591.489964346 119083.872722853 127406.07596043
#Remove timelost.

#The full model without time lost is:
POD3j = geeglm(PreAbs ~ AvgDayMat + as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#without AvgDayMat
POD3k = geeglm(PreAbs ~ as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
#Without Site
POD3l = geeglm(PreAbs ~ AvgDayMat ,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
model3C = c("POD0","POD3j","POD3k",'POD3l')
QIC3C = c(QIC(POD0)[1],QIC(POD3j)[1],QIC(POD3k)[1],QIC(POD3l)[1])
QICmod3C<-data.frame(rbind(model3C,QIC3C))
QICmod3C
#QIC            QIC.1            QIC.2            QIC.3
#model3C             POD0            POD3j            POD3k            POD3l
#QIC3C   130217.029932687 119083.872722853 121548.427290765 127396.454394609
#Full model is the best.


# STEP 6: Testing covariate significance.3
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

#In descending order:
#AvgDayMat
#Site

anova(POD3j)
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: PreAbs
#Terms added sequentially (first to last)

#Df  X2 P(>|Chi|)    
#AvgDayMat        4 314    <2e-16 ***
  #as.factor(Site)  6 407    <2e-16 ***

PODFinal = POD3j

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

#observed
#predicted     1     0
#1 26914 25513
#0  9211 37668
  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

# STEP 7: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
dimnames(AvgDayMat)<-list(NULL,c("ADBM1", "ADBM2", "ADBM3", "ADBM4"))

PODFinal = geeglm(PreAbs ~ AvgDayMat + as.factor(Site),family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)

library(boot)
library(pracma)

#Probability of covariate #1: AvgDayBasisMat:
BootstrapParameters3<-rmvnorm(10000, coef(PODFinal),summary(PODFinal)$cov.unscaled)
start=2; finish=5; Variable=SiteHourTableB$Julian; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(PODFinal)[,start:finish]*coef(PODFinal)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=6), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:5]
RealFit3<-Basis3%*%coef(PODFinal)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main ="Julian Day" , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)

#Probability of covariate #2: as.factor(Site)
#https://stackoverflow.com/questions/41518638/graph-glm-in-ggplot2-where-x-variable-is-categorical
SiteHourTableB$pr = pr
ggplot(SiteHourTableB, aes(x = Site, y = pr)) +
     geom_boxplot(aes(fill = factor(Site)), alpha = .2)



