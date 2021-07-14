### Modelling habitat preference of sperm whales using Generalized Additive Models with Generalized Estimating Equations ###

### Example with Balearic data (2003-2008) ###

## STEP 1: the data ##

# Open the R workspace provided
load("~/GitHub/SeasonalityAnalysis/Pirotta/m436p257_supp2-dataset.RData") #load dataset
head(dat); dim(dat)  # to visualise the structure of the data
attach(dat)          # to attach the dataset

# Each record in the dataset corresponds to a GPS fix, which constitutes the unit of analysis. Each fix has been associated to the value of each environmental 
# covariate in that spatial position, including multiple spatial or temporal scales for some of the variables. The column "Line_Id" specifies the block to 
# which each point belongs. The column "Pres" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact). In the 
# column "Group", the encounters with single animals are classified as 0s, while the encounters with groups as 1s; 2 identifies the absences (this is 
# just to be able to easily exclude the encounters with single animals or the encounters with groups when carrying out the analysis by grouping behaviour).  

## STEP 2: require the libraries needed ##

# All the libraries have to be installed prior to their utilization (see R help on library installation)

library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
devtools::install_github("vjcitn/yags")
library(yags)            # for the GEEs (QIC provided)
library(rjags)           # replacement for yags which is out of date
library(splines)         # to construct the B-splines within a GEE-GLM
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2)         # to build the partial residual plots
library(mvtnorm)         # to build the partial residual plots
library(gridExtra)       # to build the partial residual plots

## STEP 3: identify the best temporal or spatial scale for the covariates available at multiple scales ##

# The library yags is used to carry out model selection because it automatically provides the QIC score in the model output

# An empty model is fitted: the binary response "Pres" is modelled as a function of Latitude and Longitude only. These are expressed as B-splines with 
# one knot positioned at the average value. The independence working correlation model is used and the block is defined on the basis of the "Line_Id" values.
empty_form = formula(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long)), data = dat)
empty<-yags(empty_form,family="binomial", corstruct="independence",id=Line_Id)#, data = dat) 

# A series of models is fitted, each containing Latitude, Longitude and chlorophyll-a at one of the scales under examination. Because the package splines 
# does not allow the 'shrinkage', the inclusion of each covariate as a linear term is also tested.  

# Spatial scale: 0.05 degrees

feb1<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_feb1,knots=mean(chla_feb1)),data = dat,family=binomial, corstruct="independence",id=Line_Id)
feb1l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_feb1,family=binomial, corstruct="independence",id=Line_Id)
feb2<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_feb2,knots=mean(chla_feb2)),family=binomial, corstruct="independence",id=Line_Id)
feb2l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_feb2,family=binomial, corstruct="independence",id=Line_Id)
feb3<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_feb3,knots=mean(chla_feb3)),family=binomial, corstruct="independence",id=Line_Id)
feb3l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_feb3,family=binomial, corstruct="independence",id=Line_Id)

apr1<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_apr1,knots=mean(chla_apr1)),family=binomial, corstruct="independence",id=Line_Id)
apr1l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_apr1,family=binomial, corstruct="independence",id=Line_Id)
apr2<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_apr2,knots=mean(chla_apr2)),family=binomial, corstruct="independence",id=Line_Id)
apr2l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_apr2,family=binomial, corstruct="independence",id=Line_Id)
apr3<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_apr3,knots=mean(chla_apr3)),family=binomial, corstruct="independence",id=Line_Id)
apr3l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_apr3,family=binomial, corstruct="independence",id=Line_Id)

jun1<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_jun1,knots=mean(chla_jun1)),family=binomial, corstruct="independence",id=Line_Id)
jun1l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_jun1,family=binomial, corstruct="independence",id=Line_Id)
jun2<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_jun2,knots=mean(chla_jun2)),family=binomial, corstruct="independence",id=Line_Id)
jun2l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_jun2,family=binomial, corstruct="independence",id=Line_Id)
jun3<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_jun3,knots=mean(chla_jun3)),family=binomial, corstruct="independence",id=Line_Id)
jun3l<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_jun3,family=binomial, corstruct="independence",id=Line_Id)

# Spatial scale: 0.5 degrees

feb1x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_feb1x10,knots=mean(chla_feb1x10)),family=binomial, corstruct="independence",id=Line_Id)
feb1lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_feb1x10,family=binomial, corstruct="independence",id=Line_Id)
feb2x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_feb2x10,knots=mean(chla_feb2x10)),family=binomial, corstruct="independence",id=Line_Id)
feb2lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_feb2x10,family=binomial, corstruct="independence",id=Line_Id)
feb3x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_feb3x10,knots=mean(chla_feb3x10)),family=binomial, corstruct="independence",id=Line_Id)
feb3lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_feb3x10,family=binomial, corstruct="independence",id=Line_Id)

apr1x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_apr1x10,knots=mean(chla_apr1x10)),family=binomial, corstruct="independence",id=Line_Id)
apr1lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_apr1x10,family=binomial, corstruct="independence",id=Line_Id)
apr2x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_apr2x10,knots=mean(chla_apr2x10)),family=binomial, corstruct="independence",id=Line_Id)
apr2lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_apr2x10,family=binomial, corstruct="independence",id=Line_Id)
apr3x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_apr3x10,knots=mean(chla_apr3x10)),family=binomial, corstruct="independence",id=Line_Id)
apr3lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)

jun1x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_jun1x10,knots=mean(chla_jun1x10)),family=binomial, corstruct="independence",id=Line_Id)
jun1lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_jun1x10,family=binomial, corstruct="independence",id=Line_Id)
jun2x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_jun2x10,knots=mean(chla_jun2x10)),family=binomial, corstruct="independence",id=Line_Id)
jun2lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_jun2x10,family=binomial, corstruct="independence",id=Line_Id)
jun3x10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(chla_jun3x10,knots=mean(chla_jun3x10)),family=binomial, corstruct="independence",id=Line_Id)
jun3lx10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+chla_jun3x10,family=binomial, corstruct="independence",id=Line_Id)

# Extract the QIC scores from the yags object to compare the empty model with the others and select the best temporal and spatial scale to use.

model<-c("empty","feb1","feb1l","feb1x10","feb1lx10","feb2","feb2l","feb2x10","feb2lx10","feb3","feb3l","feb3x10","feb3lx10","apr1","apr1l","apr1x10","apr1lx10","apr2","apr2l","apr2x10","apr2lx10","apr3","apr3l","apr3x10","apr3lx10","jun1","jun1l","jun1x10","jun1lx10","jun2","jun2l","jun2x10","jun2lx10","jun3","jun3l","jun3x10","jun3lx10")
QIC<-c(empty@pan.aic,feb1@pan.aic,feb1l@pan.aic,feb1x10@pan.aic,feb1lx10@pan.aic,feb2@pan.aic,feb2l@pan.aic,feb2x10@pan.aic,feb2lx10@pan.aic,feb3@pan.aic,feb3l@pan.aic,feb3x10@pan.aic,feb3lx10@pan.aic,apr1@pan.aic,apr1l@pan.aic,apr1x10@pan.aic,apr1lx10@pan.aic,apr2@pan.aic,apr2l@pan.aic,apr2x10@pan.aic,apr2lx10@pan.aic,apr3@pan.aic,apr3l@pan.aic,apr3x10@pan.aic,apr3lx10@pan.aic,jun1@pan.aic,jun1l@pan.aic,jun1x10@pan.aic,jun1lx10@pan.aic,jun2@pan.aic,jun2l@pan.aic,jun2x10@pan.aic,jun2lx10@pan.aic,jun3@pan.aic,jun3l@pan.aic,jun3x10@pan.aic,jun3lx10@pan.aic)
QICchla<-data.frame(rbind(model,QIC))
QICchla

# model            empty             feb1            feb1l          feb1x10         feb1lx10             feb2            feb2l          feb2x10
# QIC   5026.34793615389 5187.15778286387 5114.56776261041 5174.82012046918 5094.69737877542 5128.32066946709 5093.71074321849 5187.05781045695

# model         feb2lx10             feb3            feb3l          feb3x10        feb3lx10             apr1            apr1l          apr1x10
# QIC   5091.24275816031 5123.19122565112 5079.63471700079 5102.99654102128 5068.0938861138 5227.89966708911 5112.30470325738 5051.15079409273

# model         apr1lx10             apr2            apr2l          apr2x10        apr2lx10            apr3            apr3l          apr3x10
# QIC   5061.60129276699 5178.13160903504 5057.75716607626 5072.43365973927 5005.8862255677 9600.0971098554 5024.86832535829 5079.86270230891

# model         apr3lx10            jun1            jun1l          jun1x10         jun1lx10             jun2            jun2l          jun2x10
# QIC   5000.12095826579 5151.3618764083 5091.47581326649 5165.56464281556 5072.20966807794 5161.86614789761 5083.67817220354 5165.79048250028

# model         jun2lx10             jun3            jun3l          jun3x10         jun3lx10
# QIC   5072.16031111493 5171.37258715258 5094.32682733815 5167.00226478057 5071.12590748435

# The scale selected for chlorophyll is the April peak, on the 20x20 nm scale, with standard deviation of the weights equal to 1/10 of the lag and in 
# linear form (i.e. apr3lx10), as this model corresponds to the lowest QIC score.

# The same procedure is repeated to select one single scale for slope, in the correct form (linear or smooth).

sl1x<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Slope1x,knots=mean(Slope1x)),family=binomial, corstruct="independence",id=Line_Id)
sl1xl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope1x,family=binomial, corstruct="independence",id=Line_Id)
sl5x<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Slope5x,knots=mean(Slope5x)),family=binomial, corstruct="independence",id=Line_Id)
sl5xl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope5x,family=binomial, corstruct="independence",id=Line_Id)
sl10x<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Slope10x,knots=mean(Slope10x)),family=binomial, corstruct="independence",id=Line_Id)
sl10xl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x,family=binomial, corstruct="independence",id=Line_Id)

model<-c("empty","slope1x","slope1x_linear","slope5x","slope5x_linear","slope10x","slope10x_linear")
QIC<-c(empty@pan.aic,sl1x@pan.aic,sl1xl@pan.aic,sl5x@pan.aic,sl5xl@pan.aic,sl10x@pan.aic,sl10xl@pan.aic)
QIC_slope<-data.frame(rbind(model,QIC))
QIC_slope

# model            empty         slope1x   slope1x_linear          slope5x   slope5x_linear         slope10x  slope10x_linear
# QIC   5026.34793615389 5014.0674801405 5038.57972112403 5068.51205021298 5025.94329894353 5045.63901019836 5009.06984242041

# Slope on the 20x20 nm scale, in a linear form, must be retained (i.e. slope10x_linear).

# The same procedure is repeated to select either monthly or weekly values of SST (linear or smooth).

monthly<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(SST_monthly,knots=mean(SST_monthly)),family=binomial, corstruct="independence",id=Line_Id)
monthlyl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+SST_monthly,family=binomial, corstruct="independence",id=Line_Id)
weekly<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(SST_weekly,knots=mean(SST_weekly)),family=binomial, corstruct="independence",id=Line_Id)
weeklyl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+SST_weekly,family=binomial, corstruct="independence",id=Line_Id)

model<-c("empty","monthly","monthly_linear","weekly","weekly_linear")
QIC<-c(empty@pan.aic,monthly@pan.aic,monthlyl@pan.aic,weekly@pan.aic,weeklyl@pan.aic)
QIC_sst<-data.frame(rbind(model,QIC))
QIC_sst

# model            empty          monthly  monthly_linear           weekly    weekly_linear
# QIC   5026.34793615389 5180.13860197143 5044.1039424879 5074.17086176602 5013.68200876173

# The QIC score suggests that weekly values of the SST must be used. The form is again linear (i.e. weekly_linear). 

## STEP 4: fit the full model ##

fit1<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  

# Exploratory data analysis: summarise the values of the covariates for the response "Pres". Boxplots could also be used for a visual inspection of the data.

tapply(Depth,Pres,summary)
tapply(Slope10x,Pres,summary)
tapply(Aspect,Pres,summary)
tapply(wind,Pres,summary)
tapply(SSH,Pres,summary)
tapply(SST_weekly,Pres,summary)
tapply(SST_weekly_slope,Pres,summary)
tapply(SST_mdeviation,Pres,summary)
tapply(chla_apr3x10,Pres,summary)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term. 

fitdl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Depth+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
fital<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+Aspect+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
fitwl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+wind+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
fithl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
fitwsl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
fitmdl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+SST_mdeviation+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)

model<-c("full","dl","al","wl","hl","wsl","mdl")
QIC<-c(fit1@pan.aic,fitdl@pan.aic,fital@pan.aic,fitwl@pan.aic,fithl@pan.aic,fitwsl@pan.aic,fitmdl@pan.aic)
data.frame(rbind(model,QIC))

# model             full               dl               al              wl               hl              wsl              mdl
# QIC   5271.74422487162 9529.79249364125 5308.29218593693 9328.4962017818 5236.25469756788 5250.70542492964 5266.88171399961

# SSH should be retained as linear term.

fit2<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  # i.e. fithl
fitdl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Depth+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
fital<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+Aspect+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
fitwl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+wind+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
fitwsl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
fitmdl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  

model<-c("full","dl","al","wl","wsl","mdl")
QIC<-c(fit2@pan.aic,fitdl@pan.aic,fital@pan.aic,fitwl@pan.aic,fitwsl@pan.aic,fitmdl@pan.aic)
data.frame(rbind(model,QIC))

# model             full               dl               al               wl              wsl              mdl
# QIC   5236.25469756788 5359.68898267533 5266.71325970759 5257.48584805803 5214.92397199749 9334.61689900357

# SST_weekly_slope should be retained as a linear term

fit3<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  # i.e. fitwsl
fitdl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Depth+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
fital<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+Aspect+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
fitwl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+wind+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
fitmdl<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  

model<-c("full","dl","al","wl","mdl")
QIC<-c(fit3@pan.aic,fitdl@pan.aic,fital@pan.aic,fitwl@pan.aic,fitmdl@pan.aic)
data.frame(rbind(model,QIC))

# model             full               dl               al               wl              mdl
# QIC   5214.92397199749 5339.35752588423 5244.56148984484 5236.70700508862 5221.49444337552

# All the QIC scores are higher than the score for the full model, therefore all the remaining covariates are considered in the best form (linear or smooth).

## STEP 5: perform a manual stepwise variable selection to identify the best subset of covariates ##
# The initial full model is:

fit3<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id) 

# A series of reduced models is fitted. Each contains all the covariates but one. The reduced model with the lowest QIC is the one to use in the following step.

y<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
d<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
sl10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
a<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
w<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
ssh<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
sst<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
mdev<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
ws<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
chla<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope,family=binomial, corstruct="independence",id=Line_Id)  

model<-c("full","y","d","sl10","a","w","ssh","sst","mdev","ws","chla")
QIC<-c(fit3@pan.aic,y@pan.aic,d@pan.aic,sl10@pan.aic,a@pan.aic,w@pan.aic,ssh@pan.aic,sst@pan.aic,mdev@pan.aic,ws@pan.aic,chla@pan.aic)
data.frame(rbind(model,QIC))

# model             full                y                d             sl10                a                w              ssh              sst
# QIC   5214.92397199749 5001.85722884214 5314.56599253774 5268.89315991936 5259.84066647466 5218.35811737537 5174.54748704375 5206.05871992886

# model             mdev               ws             chla
# QIC   5192.10072346801 5208.74717055123 5221.74009915839

# Year has to be removed (the model without Year is the one with the lowest QIC).
# The procedure must be continued until none of the remaining variables can be removed without increasing QIC.

fit4<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  # i.e. the model y
d<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
sl10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
a<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
w<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
ssh<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id) 
sst<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
mdev<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
ws<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)  
chla<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope,family=binomial, corstruct="independence",id=Line_Id)  

model<-c("full","d","sl10","a","w","ssh","sst","mdev","ws","chla")
QIC<-c(fit4@pan.aic,d@pan.aic,sl10@pan.aic,a@pan.aic,w@pan.aic,ssh@pan.aic,sst@pan.aic,mdev@pan.aic,ws@pan.aic,chla@pan.aic)
data.frame(rbind(model,QIC))

# model             full               d             sl10                a                w              ssh              sst             mdev
# QIC   5001.85722884214 5077.2429394218 5060.53097754001 5039.22898604842 5001.21999935888 4976.20576420094 5027.24395026128 4959.08214495535

# model               ws             chla
# QIC   4989.29621591137 5003.10646963412

# SST_mdeviation can be removed.

fit5<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id) # i.e. the model mdev
d<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sl10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
a<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
w<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
ssh<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sst<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
ws<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
chla<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope,family=binomial, corstruct="independence",id=Line_Id)

model<-c("full","d","sl10","a","w","ssh","sst","ws","chla")
QIC<-c(fit5@pan.aic,d@pan.aic,sl10@pan.aic,a@pan.aic,w@pan.aic,ssh@pan.aic,sst@pan.aic,ws@pan.aic,chla@pan.aic)
data.frame(rbind(model,QIC))

# model             full                d             sl10                a                w            ssh              sst              ws
# QIC   4959.08214495535 5034.41802080704 5005.09475080114 4993.23054794395 4956.33987701265 4931.638982768 4975.94044305331 4945.8335460295

# model             chla
# QIC   4960.72196547842

# SSH can be removed.

fit6<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id) # i.e. the model ssh
d<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sl10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
a<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
w<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sst<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
ws<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
chla<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope,family=binomial, corstruct="independence",id=Line_Id)

model<-c("full","d","sl10","a","w","sst","ws","chla")
QIC<-c(fit6@pan.aic,d@pan.aic,sl10@pan.aic,a@pan.aic,w@pan.aic,sst@pan.aic,ws@pan.aic,chla@pan.aic)
data.frame(rbind(model,QIC))

# model           full               d             sl10                a                w              sst             ws             chla
# QIC   4931.638982768 4997.0163528512 4976.04665219566 4960.98551817401 4907.16918526904 4952.86373685549 4919.284606216 4939.77543563068

# wind has to be removed.

fit7<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id) # i.e. the model w
d<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sl10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
a<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sst<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly_slope+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
ws<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
chla<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_weekly_slope,family=binomial, corstruct="independence",id=Line_Id)

model<-c("full","d","sl10","a","sst","ws","chla")
QIC<-c(fit7@pan.aic,d@pan.aic,sl10@pan.aic,a@pan.aic,sst@pan.aic,ws@pan.aic,chla@pan.aic)
data.frame(rbind(model,QIC))

# model             full               d             sl10                a              sst               ws             chla
# QIC   4907.16918526904 4997.0163528512 4935.05818149937 4930.88179159234 4923.56155384573 4896.03099630217 4917.26036903389

# SST_weekly_slope has to be removed.

fit8<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id) # i.e. the model ws
d<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sl10<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+SST_weekly+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
a<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+SST_weekly+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
sst<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+chla_apr3x10,family=binomial, corstruct="independence",id=Line_Id)
chla<-yags(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly,family=binomial, corstruct="independence",id=Line_Id)

model<-c("full","d","sl10","a","sst","chla")
QIC<-c(fit8@pan.aic,d@pan.aic,sl10@pan.aic,a@pan.aic,sst@pan.aic,chla@pan.aic)
data.frame(rbind(model,QIC))

# model             full                d            sl10                a              sst             chla
# QIC   4896.03099630217 4988.67266756013 4921.1484004708 4917.37344944322 4907.89429351439 4900.65861429135

# All the covariates, when removed, cause the QIC score to increase.
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

fit8.g<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id)

# The significance of the covariates is tested using Wald's tests (i.e. the anova.geeglm function of the library geeglm). Non-significant covariates are removed one at a time.

anova.geeglm(fit8.g)

# - chla_apr3x10 (p=0.3340)

fit9<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly,family=binomial, corstr="independence",id=Line_Id)
anova.geeglm(fit9)

# - SST_weekly (p=0.1674)

fit10<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect)),family=binomial, corstr="independence",id=Line_Id)
anova.geeglm(fit10)

# - Slope10x (p=0.0917)

fit11<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect)),family=binomial, corstr="independence",id=Line_Id)
anova.geeglm(fit11)

# fit11 represents the final model.

## STEP 6: construct the confusion matrix ##

# Construction of the ROC curve

pr <- predict(fit11,dat, type="response")                          # the final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
pred <- prediction(pr,dat$Pres)                                    # to specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
perf <- performance(pred, measure="tpr", x.measure="fpr")          # to assess model performance in the form of the true positive rate and the false positive rate
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5)) # to plot the ROC curve

# Choice of the best cut-off probability

y<-as.data.frame(perf@y.values)
x<-as.data.frame(perf@x.values)
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
write.table(d,"C:\\distances.txt")                                 # to write a table with the computed distances

# The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
# MAX d= 0.277124605038953 --> position 2098

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
alpha[2098,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve

# Best cutoff:   0.2516195
# This value can now be used to build the confusion matrix:

DATA<-matrix(0,4973,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim(dat)) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:4973                                                # the first column is filled with an ID value that is unique for each row
DATA$Observed<-dat$Pres                                            # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(fit11,dat,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = 0.2516195)                                   # the identified cut-off must be used here

#             observed
# predicted     1    0
#          1  815 1283
#          0  310 2565

# The confusion matrix can then be transformed into percentages:

#             observed
# predicted     1    0
#          1  72%  33%
#          0  24%  67%

# The area under the curve can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# auc = 0.7723559

## STEP 7: construct a prediction map ##

pred<-read.table("C:\\Grid_prediction.txt", header=TRUE)            # import a grid of points created in a GIS. Each point should be associated with the values of the covariates that are retained in the final model
pred<-cbind(pred,predict(fit11,pred, type="response"))              # predict animal occurrence using the final model
write.table(pred, "C:\\prediction.txt")                             # save the results in a new text file, that can then be imported in a GIS to draw the prediction map

## STEP 8: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the estimated relationship between the response (on the link scale) and each predictor ## 

LatForPlotting<- seq(min(Lat), max(Lat), length=50)
LongForPlotting<- seq(min(Long), max(Long), length=50)
DepthForPlotting<- seq(min(Depth), max(Depth), length=50)
AspectForPlotting<- seq(min(Aspect), max(Aspect), length=50)
PredictionData<- expand.grid(LatForPlotting,LongForPlotting,DepthForPlotting,AspectForPlotting)
names(PredictionData)<- c("LatForPlotting","LongForPlotting","DepthForPlotting","AspectForPlotting")
BootstrapParameters<-rmvnorm(10000, coef(fit11), summary(fit11)$cov.unscaled)
test<- glm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect)),family=binomial)
x1<-model.matrix(test)[,2:5]%*%coef(fit11)[c(2:5)]
x2<-model.matrix(test)[,6:9]%*%coef(fit11)[c(6:9)]
x3<-model.matrix(test)[,10:13]%*%coef(fit11)[c(10:13)]
x4<-model.matrix(test)[,14:17]%*%coef(fit11)[c(14:17)]

# Partial plot for Latitude

BootstrapCoefs<- BootstrapParameters[,1:5]
Basis<- cbind( rep(1, 10),bs(LatForPlotting, knots=mean(Lat), Boundary.knots=range(Lat)))
RealFit<- Basis%*%coef(fit11)[c(1:5)]
RealFitCenter1<- RealFit-mean(x1)-coef(fit11)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim1<- min(cis-mean(x1)-coef(fit11)[1])
MaximumYlim1<- max(cis-mean(x1)-coef(fit11)[1])
cil1<-cis[1,]-mean(x1)-coef(fit11)[1]
ciu1<-cis[2,]-mean(x1)-coef(fit11)[1]
plot1<-qplot(LatForPlotting,RealFitCenter1,xlab="Latitude", ylab="s(Latitude)", main="a",ylim=c(MinimumYlim1, MaximumYlim1),geom="line")+opts(axis.line = theme_segment(),panel.background=theme_blank(),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil1,ymax=ciu1),stat="identity")+geom_rug(aes(x = Lat, y=-10000))

# Partial plot for Longitude

BootstrapCoefs<- BootstrapParameters[,c(1,6:9)]
Basis<- cbind(rep(1,10), bs(LongForPlotting, knots=mean(Long), Boundary.knots=range(Long)))
RealFit<- Basis%*%coef(fit11)[c(1,6:9)]
RealFitCenter2<- RealFit-mean(x2)-coef(fit11)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim2<- min(cis-mean(x2)-coef(fit11)[1])
MaximumYlim2<- max(cis-mean(x2)-coef(fit11)[1])
cil2<-cis[1,]-mean(x2)-coef(fit11)[1]
ciu2<-cis[2,]-mean(x2)-coef(fit11)[1]
plot2<-qplot(LongForPlotting,RealFitCenter2,xlab="Longitude", ylab="s(Longitude)", main="b",ylim=c(MinimumYlim2, MaximumYlim2),geom="line")+opts(axis.line = theme_segment(),panel.background=theme_blank(),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil2,ymax=ciu2),stat="identity")+geom_rug(aes(x = Long, y=-10000))

# Partial plot for Depth

BootstrapCoefs<- BootstrapParameters[,c(1,10:13)]
Basis<- cbind(rep(1,10), bs(DepthForPlotting, knots=mean(Depth), Boundary.knots=range(Depth)))
RealFit<- Basis%*%coef(fit11)[c(1,10:13)]
RealFitCenter3<- RealFit-mean(x3)-coef(fit11)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim3<- min(cis-mean(x3)-coef(fit11)[1])
MaximumYlim3<- max(cis-mean(x3)-coef(fit11)[1])
cil3<-cis[1,]-mean(x3)-coef(fit11)[1]
ciu3<-cis[2,]-mean(x3)-coef(fit11)[1]
plot3<-qplot(DepthForPlotting,RealFitCenter3,xlab="Depth", ylab="s(Depth)", main="c",ylim=c(MinimumYlim3, MaximumYlim3),geom="line")+opts(axis.line = theme_segment(),panel.background=theme_blank(),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil3,ymax=ciu3),stat="identity")+geom_rug(aes(x = Depth, y=-10000))

# Partial plot for Aspect

BootstrapCoefs<- BootstrapParameters[,c(1,14:17)]
Basis<- cbind(rep(1,10), bs(AspectForPlotting, knots=mean(Aspect), Boundary.knots=range(Aspect)))
RealFit<- Basis%*%coef(fit11)[c(1,14:17)]
RealFitCenter4<- RealFit-mean(x4)-coef(fit11)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim4<- min(cis-mean(x4)-coef(fit11)[1])
MaximumYlim4<- max(cis-mean(x4)-coef(fit11)[1])
cil4<-cis[1,]-mean(x4)-coef(fit11)[1]
ciu4<-cis[2,]-mean(x4)-coef(fit11)[1]
plot4<-qplot(AspectForPlotting,RealFitCenter4,xlab="Aspect", ylab="s(Aspect)", main="d",ylim=c(MinimumYlim4, MaximumYlim4),geom="line")+opts(axis.line = theme_segment(),panel.background=theme_blank(),panel.grid.major=theme_blank(),panel.grid.minor=theme_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil4,ymax=ciu4),stat="identity")+geom_rug(aes(x = Aspect, y=-10000))

# Plot the four plots altogether

sidebysideplot <- grid.arrange(plot1,plot2,plot3,plot4,ncol=2,nrow=2)

## STEP 9: analysis by grouping behaviour ##

# The two subsets of the data set can be created as follows:

groups  <- subset(dat,Group!=0) # to select only the encounters with groups
singles <- subset(dat,Group!=1) # to select only the encounters with single males

# All the steps from 1 to 8 should now be repeated using each of the two subsets, in order to assess habitat use by groups and singletons.

{"mode":"full","isActive":false}