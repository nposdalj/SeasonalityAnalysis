### Modelling habitat preference of sperm whales using Generalized Additive Models with Generalized Estimating Equations ###
### Script adapted from Pirotta et al. (2011) and Benjamins ###  
### 9 Models total:
  #Site specific models: 'HZ','OC','NC','BC','WC','GS','BP','BS','JAX'
  #Region specific models: North + South
  #Big model: all 9 sites

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
library(ChemoSpecUtils)
library(car) #for ANOVA
library(splines2)       # to use mSpline for the GEEs
library(ggfortify)      # extract confidence interval for ACF plots
library(dplyr)

site = 'WC' #specify the site of interest
GDrive = 'G'

# Step 1: Load the data ---------------------------------------------------
#Hourly data
dir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/AllSites",sep="")
saveDir = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/Plots/",site, sep="")
saveWorkspace = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/",site,'/',sep="")
fileName = paste(GDrive,":/My Drive/WAT_TPWS_metadataReduced/SeasonalityAnalysis/AllSites/AllSitesGrouped_Binary_GAMGEE_ROW_sexClasses.csv",sep="") #setting the directory
HourTable = read.csv(fileName)
HourTable = na.omit(HourTable)
HourTable$date = as.Date(HourTable$tbin)
HourTable$tbin = as.POSIXct(HourTable$tbin)
HourTable = HourTable[ order(HourTable$tbin , decreasing = FALSE ),]
SiteHourTable = dplyr::filter(HourTable, grepl(site,Site))
SiteHourTable$Hour = hour(SiteHourTable$tbin)
SiteHourTable$Effort_Bin[SiteHourTable$Effort_Bin > 12] = 12
SiteHourTable = subset(SiteHourTable, Site == site)#subset the table for the site only
SiteHourTable = SiteHourTable[!(SiteHourTable$Julian==366),] #If it's a leap year, delete julian day 366
SiteHourTable$TimeLost = max(SiteHourTable$Effort_Bin) - SiteHourTable$Effort_Bin #Time lost variable

# Each record in the dataset corresponds to an hour of recording, which constitutes the unit of analysis. 
# The column "PreAbs" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).
# The column "Site" corresponds to the recording site.
# The column "Region" corresponds to the recording region (GOA or BSAI).
# The column Julian represents the day of the year and the column year represents the year of recording.


# Step 2: Identify the Best Blocking Structure ----------------------------
#Calculate acf on model residuals
#Females
BlockModF<-glm(PreAbsF~
               bs(Julian, k = 4)+
               TimeLost+
               as.factor(Year)
             ,data=SiteHourTable,family=binomial)

#Juveniles
BlockModJ<-glm(PreAbsJ~
                 bs(Julian, k = 4)+
                 TimeLost+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

#Males
BlockModM<-glm(PreAbsM~
                 bs(Julian, k = 4)+
                 TimeLost+
                 as.factor(Year)
               ,data=SiteHourTable,family=binomial)

#Social Groups
ACFF = acf(residuals(BlockModF), lag.max = 2000, ylim=c(0,0.1))
CIF = ggfortify:::confint.acf(ACFF)
ACFidxF = which(ACFF[["acf"]] < CIF, arr.ind=TRUE)
ACFvalF = ACFidxF[1]

#Mid Size
ACFJ = acf(residuals(BlockModJ), lag.max = 2000, ylim=c(0,0.1))
CIJ = ggfortify:::confint.acf(ACFJ)
ACFidxJ = which(ACFJ[["acf"]] < CIJ, arr.ind=TRUE)
ACFvalJ = ACFidxJ[1]

#Males
ACFM = acf(residuals(BlockModM), lag.max = 2000, ylim=c(0,0.1))
CIM = ggfortify:::confint.acf(ACFM)
ACFidxM = which(ACFM[["acf"]] < CIM, arr.ind=TRUE)
ACFvalM = ACFidxM[1]

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
timeseriesF = dplyr::rename(timeseriesF, tbin = date)
SiteHourTableB = left_join(SiteHourTable,timeseriesF,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$blockF))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksF = UnBlock$sequence[match(SiteHourTableB$blockF,UnBlock$`unique(SiteHourTableB$blockF)`)]
difference = diff(SiteHourTableB$BlocksF) #find difference between rows

#Juveniles
timeseriesJ = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseriesJ)/ACFvalJ)), times=1, each=ACFvalJ)
divdiff = nrow(timeseriesJ) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseriesJ$blockJ = c(preBlock,lastVec)
timeseriesJ = dplyr::rename(timeseriesJ, tbin = date)
SiteHourTableB = left_join(SiteHourTableB,timeseriesJ,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$blockJ))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksJ = UnBlock$sequence[match(SiteHourTableB$blockJ,UnBlock$`unique(SiteHourTableB$blockJ)`)]
difference = diff(SiteHourTableB$BlocksJ) #find difference between rows

#Males
timeseriesM = data.frame(date=seq(startDate, endDate, by="hours"))
preBlock = rep(1:(floor(nrow(timeseriesM)/ACFvalM)), times=1, each=ACFvalM)
divdiff = nrow(timeseriesM) - length(preBlock)
last = tail(preBlock, n = 1)+1
lastVec = rep(last,each = divdiff)
timeseriesM$blockM = c(preBlock,lastVec)
timeseriesM = dplyr::rename(timeseriesM, tbin = date)
SiteHourTableB = left_join(SiteHourTableB,timeseriesM,by="tbin")

#Make blocks continuous 
UnBlock = as.data.frame(unique(SiteHourTableB$blockM))
UnBlock$sequence = seq.int(nrow(UnBlock))
SiteHourTableB$BlocksM = UnBlock$sequence[match(SiteHourTableB$blockM,UnBlock$`unique(SiteHourTableB$blockM)`)]
difference = diff(SiteHourTableB$BlocksM) #find difference between rows

# Step 3: ANOVA to Check for Significance of Variables --------------------
#ANOVA with car package
Anova(BlockModF)
#HZ
#bs(Julian, k = 4)      886  4    < 2e-16 ***
#TimeLost                19  1    1.1e-05 ***
#as.factor(Year)        482  4    < 2e-16 ***

#OC
# bs(Julian, k = 4)  155.671  4  < 2.2e-16 ***
# TimeLost            27.626  1  1.472e-07 ***
# as.factor(Year)    161.145  4  < 2.2e-16 ***

#NC
# bs(Julian, k = 4)      286  4    < 2e-16 ***
# TimeLost                31  1    2.9e-08 ***
# as.factor(Year)        636  4    < 2e-16 ***

#BC
# bs(Julian, k = 4)   375.98  4  < 2.2e-16 ***
# TimeLost             18.81  1  1.443e-05 ***
# as.factor(Year)    1186.99  3  < 2.2e-16 ***

#WC
# bs(Julian, k = 4)      348  4     <2e-16 ***
# TimeLost                 7  1     0.0066 ** 
# as.factor(Year)       1224  3     <2e-16 ***

Anova(BlockModJ)
#HZ
#bs(Julian, k = 4)      919  4     <2e-16 ***
#TimeLost                17  1      3e-05 ***
#as.factor(Year)        238  4     <2e-16 ***

#OC
# bs(Julian, k = 4)  2042.34  4  < 2.2e-16 ***
# TimeLost             15.56  1  7.995e-05 ***
# as.factor(Year)     542.68  4  < 2.2e-16 ***

#NC
# bs(Julian, k = 4)     1846  4    < 2e-16 ***
# TimeLost                12  1    0.00062 ***
# as.factor(Year)        209  4    < 2e-16 ***

#BC
# bs(Julian, k = 4)   871.40  4  < 2.2e-16 ***
#   TimeLost             19.53  1  9.898e-06 ***
#   as.factor(Year)      10.24  3    0.01664 *  

#WC
# bs(Julian, k = 4)      677  4    < 2e-16 ***
# TimeLost                17  1    4.5e-05 ***
# as.factor(Year)         24  3    2.0e-05 ***
  
Anova(BlockModM)
#HZ
#bs(Julian, k = 4)      840  4    < 2e-16 ***
#TimeLost                 8  1    0.00361 ** 
#as.factor(Year)         23  4    0.00015 ***

#OC
# bs(Julian, k = 4)   717.50  4     <2e-16 ***
# TimeLost              1.31  1     0.2525    
# as.factor(Year)     128.16  4     <2e-16 ***

#NC
# bs(Julian, k = 4)     1015  4    < 2e-16 ***
# TimeLost                 1  1        0.3    
# as.factor(Year)         57  4    1.5e-11 ***

#BC
# bs(Julian, k = 4)  137.732  4  < 2.2e-16 ***
# TimeLost             2.076  1     0.1497    
# as.factor(Year)     28.444  3   2.93e-06 ***

#WC
# bs(Julian, k = 4)     73.4  4    4.4e-15 ***
# TimeLost               0.1  1       0.74    
# as.factor(Year)       28.5  3    2.9e-06 ***

# Step 4: Data Exploration and Initial Analysis ------------------------------------------------------------------
# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).
# Basic model for VIF analysis:

  #Females
  GLMF = glm(PreAbsF~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
  #Males
  GLMJ = glm(PreAbsJ~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)
  #Juveniles
  GLMM = glm(PreAbsM~Julian+TimeLost+as.factor(Year),family=binomial,data=SiteHourTableB)

#VIF scores in GLM to work out collinearity:
VIF(GLMF)
VIF(GLMJ)
VIF(GLMM)

#HZ
#Julian          1.44  1            1.20
#TimeLost        1.00  1            1.00
#as.factor(Year) 1.44  4            1.05

#Julian          1.54  1            1.24
#TimeLost        1.00  1            1.00
#as.factor(Year) 1.54  4            1.06

#Julian          1.26  1            1.12
#TimeLost        1.00  1            1.00
#as.factor(Year) 1.26  4            1.03

#OC
# Julian          1.232183  1        1.110037
# TimeLost        1.002017  1        1.001008
# as.factor(Year) 1.231875  4        1.026410
# 
# Julian          1.173532  1        1.083297
# TimeLost        1.003268  1        1.001633
# as.factor(Year) 1.172871  4        1.020132
# 
# Julian          1.143291  1        1.069248
# TimeLost        1.004415  1        1.002205
# as.factor(Year) 1.141632  4        1.016695

#NC
# Julian          1.31  1            1.14
# TimeLost        1.00  1            1.00
# as.factor(Year) 1.31  4            1.03

# Julian          1.24  1            1.11
# TimeLost        1.00  1            1.00
# as.factor(Year) 1.23  4            1.03

# Julian          1.23  1            1.11
# TimeLost        1.00  1            1.00
# as.factor(Year) 1.23  4            1.03

#BC
# Julian          1.403842  1        1.184838
# TimeLost        1.001959  1        1.000979
# as.factor(Year) 1.405638  3        1.058390

# Julian          1.243672  1        1.115201
# TimeLost        1.001864  1        1.000932
# as.factor(Year) 1.244674  3        1.037152

# Julian          1.159798  1        1.076939
# TimeLost        1.001586  1        1.000793
# as.factor(Year) 1.159718  3        1.025004

#WC
# Julian          1.44  1            1.20
# TimeLost        1.02  1            1.01
# as.factor(Year) 1.46  3            1.07

# Julian          1.28  1            1.13
# TimeLost        1.02  1            1.01
# as.factor(Year) 1.29  3            1.04

# Julian          1.51  1            1.23
# TimeLost        1.04  1            1.02
# as.factor(Year) 1.53  3            1.07

# Step 5: Model Selection - Covariate Preparation -------------------------
# Construct variance-covariance matrices for cyclic covariates:
#Females
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


#Juvenile
POD0j<-geeglm(PreAbsJ ~ 1, family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)


#Male
POD0m<-geeglm(PreAbsM ~ 1, family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)


# Step 6: Determine which covariates are most relevant --------------------
#and which can be removed (on the basis of previous collinearity work). 
#The reduced model with the lowest QIC is the one to use in the following step.
#Females
#The initial full model is:
POD3fa = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without AvgDayMat
POD3fb = geeglm(PreAbsF ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
#without Year
POD3fc = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)

model3fA = c("POD0f","POD3fa","POD3fb","POD3fc")
QIC3fA = c(QIC(POD0f)[1],QIC(POD3fa)[1],QIC(POD3fb)[1],QIC(POD3fc)[1])
QICmod3fA<-data.frame(rbind(model3fA,QIC3fA))
QICmod3fA

#HZ
#QIC            QIC.1            QIC.2            QIC.3
#model3fA           POD0f           POD3fa           POD3fb           POD3fc
#QIC3fA   35534.540037402 34652.3272649365 35308.9594610671 34995.1320247913
#Full model is best
#Model Order - Julian Day, Year

#OC
# QIC            QIC.1            QIC.2            QIC.3
# model3fA            POD0f           POD3fa           POD3fb           POD3fc
# QIC3fA   33540.8939891145 33357.6299158409 33457.3594848748 33288.0598266738
#Full model is best
#Model Order - Julian Day, Year

#NC
# QIC            QIC.1            QIC.2            QIC.3
# model3fA            POD0f           POD3fa           POD3fb           POD3fc
# QIC3fA   34424.5812332557 33816.3590822841 34055.8707388389 34327.7762331342
#Full model is best
#Model Order - Julian Day, Year

#BC
# QIC            QIC.1            QIC.2            QIC.3
# model3fA            POD0f           POD3fa           POD3fb           POD3fc
# QIC3fA   24868.1322744819 23317.9471990927 23654.0263091128 24457.1150825552
#Full model is best
#Model Order -Year, Julian Day

#WC
# QIC            QIC.1            QIC.2            QIC.3
# model3fA            POD0f           POD3fa           POD3fb           POD3fc
# QIC3fA   20975.4099442765 18927.7256130296 19166.2312160498 20094.3636360685
#Full model is best
#Model Order -Year, Julian Day

#Juveniles
#The initial full model is:
POD3ja = geeglm(PreAbsJ ~ AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without AvgDayMat
POD3jb = geeglm(PreAbsJ ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
#without Year
POD3jc = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

model3jA = c("POD0j","POD3ja","POD3jb","POD3jc")
QIC3jA = c(QIC(POD0j)[1],QIC(POD3ja)[1],QIC(POD3jb)[1],QIC(POD3jc)[1])
QICmod3jA<-data.frame(rbind(model3jA,QIC3jA))
QICmod3jA
#HZ
#QIC            QIC.1            QIC.2            QIC.3
#model3jA            POD0j           POD3ja           POD3jb           POD3jc
#QIC3jA   16161.0497192109 15108.4668660763 15926.9380495206 15356.2833301027
#Full model is best
#Model Order - Julian Day, Year

#OC
# QIC            QIC.1            QIC.2            QIC.3
# model3jA            POD0j           POD3ja           POD3jb           POD3jc
# QIC3jA   18227.6382129555 15985.1146704017 17916.2396399579 16535.1316953002
#Full model best
#order-jday then year

#NC
# QIC            QIC.1            QIC.2            QIC.3
# model3jA            POD0j           POD3ja           POD3jb           POD3jc
# QIC3jA   23447.3543208099 21956.2110368831 23483.3995161646 22028.5210385183
#Full model is best
#Model Order - Julian Day, Year

#BC
# QIC            QIC.1            QIC.2           QIC.3
# model3jA            POD0j           POD3ja           POD3jb          POD3jc
# QIC3jA   13093.5024108312 12286.4489312997 13110.6283467241 12273.093187511
#Full model is best
#Model Order - Julian Day, Year

#WC
# QIC            QIC.1            QIC.2            QIC.3
# model3jA            POD0j           POD3ja           POD3jb           POD3jc
# QIC3jA   9528.54058792815 8895.97317038746 9501.02493689036 8905.42534429207
#Full model is best
#Model Order - Julian Day, Year

#Males
#The initial full model is:
POD3ma = geeglm(PreAbsM ~ AvgDayMatM+as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without AvgDayMat
POD3mb = geeglm(PreAbsM ~ as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
#without Year
POD3mc = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

model3mA = c("POD0m","POD3ma","POD3mb","POD3mc")
QIC3mA = c(QIC(POD0m)[1],QIC(POD3ma)[1],QIC(POD3mb)[1],QIC(POD3mc)[1])
QICmod3mA<-data.frame(rbind(model3mA,QIC3mA))
QICmod3mA
#HZ
#QIC            QIC.1            QIC.2            QIC.3
#model3mA            POD0m           POD3ma           POD3mb           POD3mc
#QIC3mA   7185.54933472434 5989.75261850721 6772.72891517581 5999.97795529014
#Full model is best
#Model Order - Julian Day, Year

#OC
# QIC            QIC.1            QIC.2            QIC.3
# model3mA            POD0m           POD3ma           POD3mb           POD3mc
# QIC3mA   4555.64015796564 3801.42159750823 4578.63152492947 3894.85147323616
#Full model is best
#Model Order - Julian Day, Year

#NC
# QIC            QIC.1            QIC.2            QIC.3
# model3mA            POD0m           POD3ma           POD3mb           POD3mc
# QIC3mA   6662.91974809815 5548.08846446073 6531.07402564486 5589.48961523317
#Full model is best
#Model Order - Julian Day, Year

#BC
# QIC            QIC.1           QIC.2            QIC.3
# model3mA            POD0m           POD3ma          POD3mb           POD3mc
# QIC3mA   2047.35391249467 1933.00630211378 2044.4887775271 1945.35480754113
#Full model is best
#Model Order - Julian Day, Year

#WC
# QIC            QIC.1            QIC.2            QIC.3
# model3mA            POD0m           POD3ma           POD3mb           POD3mc
# QIC3mA   1555.65040286281 1528.20056341222 1553.37911093667 1529.69529680764
#Full model is best
#Model Order - Julian Day, Year

# Step 7: Finalize Model --------------------------------------------------
#Females
if(site == 'BC' || site =='WC'){
  dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
  PODFinalF = geeglm(PreAbsF ~ as.factor(Year)+AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
}else{
dimnames(AvgDayMatF)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalF = geeglm(PreAbsF ~ AvgDayMatF+as.factor(Year),family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
}
#Juveniles
dimnames(AvgDayMatJ)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalJ = geeglm(PreAbsJ ~  AvgDayMatJ+as.factor(Year),family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)

#Males
dimnames(AvgDayMatM)<-list(NULL,c("ADBM1", "ADBM2"))
PODFinalM = geeglm(PreAbsM ~ AvgDayMatM + as.factor(Year),family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)

# STEP 8: Interpreting the summary of the model --------------------------
# How to intepret model results
# Standard error - robust estimate - provides reasonable variance estimates even when the specified correlation model is in correct
# Wald - square of the z-statistic reproted by the gee function
# P - values are the upper tailed probabilities from the chi-squared random variable with 1 degree of freedom...
# distribution and test whether the true parameter value is different from zero

anova(PODFinalF)
#HZ
#AvgDayMatF       2 14.90   0.00058 ***
#as.factor(Year)  4  6.48   0.16598    

#OC
# AvgDayMatF       2 14.7717   0.00062 ***
# as.factor(Year)  4  5.9087   0.20607    

#NC
# AvgDayMatF       2 11.1    0.0039 ** 
# as.factor(Year)  4 60.7   2.1e-12 ***

#BC
# as.factor(Year)  3 109.44 < 2.2e-16 ***
#   AvgDayMatF       2  23.45  8.09e-06 ***

#WC
# AvgDayMatF       2 22.8   1.1e-05 ***
# as.factor(Year)  3 59.5   7.6e-13 ***
anova(PODFinalJ)
#HZ
#AvgDayMatJ       2 193.5   < 2e-16 ***
#as.factor(Year)  4  34.5   5.8e-07 ***

#OC
# AvgDayMatJ       2 117.031 < 2.2e-16 ***
# as.factor(Year)  4  41.201 2.443e-08 ***

#NC
# AvgDayMatJ       2 94.2    <2e-16 ***
# as.factor(Year)  4 12.3     0.015 *  

#BC
# AvgDayMatJ       2 255.554    <2e-16 ***
# as.factor(Year)  3   1.791     0.617    

#WC
# AvgDayMatJ       2 200.7    <2e-16 ***
#   as.factor(Year)  3   9.2     0.027 *  
anova(PODFinalM)
#HZ
#AvgDayMatM       2 160.2    <2e-16 ***
#as.factor(Year)  4   6.4      0.17    

#OC
# AvgDayMatM       2 49.657 1.649e-11 ***
# as.factor(Year)  4 32.667 1.397e-06 ***

#NC
# AvgDayMatM       2 52.6   3.8e-12 ***
# as.factor(Year)  4 13.8    0.0079 ** 

#BC
# AvgDayMatM       2 68.594 1.221e-15 ***
# as.factor(Year)  3  8.425     0.038 *  

#WC
# AvgDayMatM       2 14.76   0.00062 ***
# as.factor(Year)  3  6.07   0.10838    

#Remove any insignificant variables
if (site == 'HZ'){
  #Remove Year from Female and Male models
  PODFinalF = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
}
if (site == 'OC'){
  #Remove Year from Female
  PODFinalF = geeglm(PreAbsF ~ AvgDayMatF,family = binomial, corstr="ar1", id=BlocksF, data=SiteHourTableB)
}
if (site == 'BC'){
  #Remove Year from Juvenile
  PODFinalJ = geeglm(PreAbsJ ~ AvgDayMatJ,family = binomial, corstr="ar1", id=BlocksJ, data=SiteHourTableB)
}

if (site == 'WC'){
  #Remove Year from Male
  PODFinalM = geeglm(PreAbsM ~ AvgDayMatM,family = binomial, corstr="ar1", id=BlocksM, data=SiteHourTableB)
}
filename = paste(saveWorkspace,site,'_SiteSpecific_sexClasses_ModelSummary.txt',sep="")
sink(filename)
summary(PODFinalF)
anova(PODFinalF)
summary(PODFinalJ)
anova(PODFinalJ)
summary(PODFinalM)
anova(PODFinalM)
sink(file = NULL)
  
# Step 9: Construction of the ROC curve    --------------------------------
#Females
prf <- predict(PODFinalF, type="response")  
pred <- prediction(prf,SiteHourTableB$PreAbsF) 
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
DATA$Observed<-SiteHourTableB$PreAbsF                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalF,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here
  
# The area under the curve (auc) can also be used as an rough indication of model performance:
  
auc <- performance(pred, measure="auc")

#Juveniles
prj <- predict(PODFinalJ, type="response")  
pred <- prediction(prj,SiteHourTableB$PreAbsJ) 
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
DATA$Observed<-SiteHourTableB$PreAbsJ                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalJ,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here

# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

#Males
prm <- predict(PODFinalM, type="response")  
pred <- prediction(prm,SiteHourTableB$PreAbsM) 
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
DATA$Observed<-SiteHourTableB$PreAbsM                                          # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(PODFinalM,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = cutoff)                                   # the identified cut-off must be used here


# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# Step 10: Save Workspace -------------------------------------------------
fileName = paste(saveWorkspace,site,'_SiteSpecific_gamgeeOutput_sexClasses.RData',sep="")
save.image(file = fileName)