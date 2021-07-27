### Assessing influence of environmental covariates on harbour porpoise detections in a tidally active site ###
### Script adapted from Pirotta et al. (2011) ###  
### Example from Scarba array data (mooring M7) ###

## STEP 1: the data ##
# Open the R workspace provided
M7 = read.csv('C:/Users/nposd/Documents/GitHub/SeasonalityAnalysis/Benjamins/M7.csv')
M7 = zap_labels(M7)
head(M7); dim(M7) # to visualise the structure of the data
attach(M7) # to attach the dataset

# Each record in the dataset corresponds to a single minute, which constitutes the unit of analysis. Each minute has been associated to the value of each environmental
# covariate at that time. The column "Panel" specifies the block to which each point belongs. The column "DPM" represents the binary response variable (0: no animal in acoustic contact; 1: acoustic contact).

## STEP 2: require the libraries needed ##
# All the libraries have to be installed prior to their utilization (see R help on library installation)
library(geepack) # for the GEEs (Wald's hypothesis tests allowed)
library(geeglm) # for the GEEs (QIC provided)
library(splines) # to construct the B-splines within a GEE-GLM
library(ROCR) # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2) # to build the partial residual plots
library(mvtnorm) # to build the partial residual plots
library(regclass)

## STEP 3: Data exploration and initial analysis ##
# Investigate autocorrelation to work out block sizes: 
acf(M7$DPM, lag.max=300)
acf(M7$DPM, lag.max=310, ylim=c(0,0.1), xlim =c(300,310)) 
#acf_M7 = 302 minutes

# Follow data exploration protocols suggested by Zuur et al. (2010), Zuur (2012), to generate pair plots, box plots, and assess collinearity between covariates in the dataset using Varinace Inflation Factors (vif).

# Basic model for VIF analysis:
GLM1<-glm(DPM ~ DATENO + HOUR + PODAngle + X.TimeLost + Nall.m + MINUTE, family = binomial, data=M7)
#VIF scores in GLM to work out collinearity:
VIF(GLM1)
#Their Results
#DATENO       HOUR   PODAngle X.TimeLost     Nall.m MINUTE 
#1.009453   1.019137   1.208167   1.256448   1.260450   1.191046 

#My Results
#DATENO       HOUR   PODAngle X.TimeLost     Nall.m     MINUTE 
#1.007072   1.009321   1.010812   1.264334   1.267977   1.000300 

## STEP 4: Model selection - covariate preparation ##
# Construct variance-covariance matrices for cyclic covariates:
TideBasis<-gam(DPM~s(MINUTE, bs="cc", k=6), fit=F, data=M7, family=binomial, knots=list(MINUTE=seq(0,1,length=6)))$X[,2:5]
AvgHrBasis<-gam(DPM~s(HOUR, bs="cc", k=6), fit=F, data=M7, family=binomial, knots=list(HOUR=seq(0,23,length=6)))$X[,2:5]

TideBasisMat<-as.matrix(TideBasis)
AvgHrBasisMat<-as.matrix(AvgHrBasis)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term vs. a smoother and compared against an empty model.
# Extract the QIC scores from the geeglm object to compare the empty model with the others and select how to treat the covariate.

POD0<-geeglm(DPM ~ 1, family = binomial, corstr="independence", id=Panel, data=M7)

#DATENO
POD0a<-geeglm(DPM ~ bs(DATENO , knots=mean(DATENO)), family = binomial, corstr="independence", id=Panel, data=M7)
POD0b<-geeglm(DPM ~ as.factor(DATENO), family = binomial, corstr="independence", id=Panel, data=M7)
POD0c<-geeglm(DPM ~ DATENO, family = binomial, corstr="independence", id=Panel, data=M7)

model0A<-c("POD0", "POD0a", "POD0b", "POD0c")
QIC0A<-c(QIC(POD0)[1],QIC(POD0a)[1],QIC(POD0b)[1],QIC(POD0c)[1])
QICmod0A<-data.frame(rbind(model0A,QIC0A))
QICmod0A
#Their results
#X1               X2               X3               X4
#model0A             POD0            POD0a            POD0b            POD0c
#QIC0A   15401.3578427153 14173.6778372064 202025.047204914 14494.8958234142
#So DATENO should be treated as a smoother
#My results showed that DATENO should be treated as a factor

#PODAngle:
POD0f<-geeglm(DPM ~ bs(PODAngle, knots=mean(PODAngle)) , family = binomial, corstr="independence", id=Panel, data=M7)
POD0g<-geeglm(DPM ~ PODAngle, family = binomial, corstr="independence", id=Panel, data=M7)

model0C<-c("POD0", "POD0f", "POD0g")
QIC0C<-c(QIC(POD0)[1],QIC(POD0f)[1],QIC(POD0g)[1])
QICmod0C<-data.frame(rbind(model0C,QIC0C))
QICmod0C
#Their results
#X1               X2               X3
#model0C             POD0            POD0f            POD0g
#QIC0C   15401.3578427153 15080.8965136641 15109.7867620991
#Angle should be treated as a smoother.
#My results showed that angled should also be treated as a smoother.

#X.TimeLost:
POD0h<-geeglm(DPM ~ bs(X.TimeLost , knots=mean(X.TimeLost)), family = binomial, corstr="independence", id=Panel, data=M7)
POD0i<-geeglm(DPM ~ X.TimeLost, family = binomial, corstr="independence", id=Panel, data=M7)

model0D<-c("POD0", "POD0h", "POD0i")
QIC0D<-c(QIC(POD0)[1],QIC(POD0h)[1],QIC(POD0i)[1])
QICmod0D<-data.frame(rbind(model0D,QIC0D))
QICmod0D
#Their results
#X1               X2               X3
#model0D             POD0            POD0h            POD0i
#QIC0D   15401.3578427153 30736.3848027354 15402.0479608358
#X.TimeLost should be treated as linear.
#My results is that X.Timelost should be treated as a smoother.

## STEP 4: Determine which covariates are most relevant, and which can be removed (on the basis of previous collinearity work).         
# The empty model is:
POD0<-geeglm(DPM ~ 1, family = binomial, corstr="independence", id=Panel, data=M7)
# The initial full model is:
POD1<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO)) + AvgHrBasisMat + bs(PODAngle, knots=mean(PODAngle)) + X.TimeLost + TideBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

#A series of reduced models is fitted. Each contains all the covariates but one. The reduced model with the lowest QIC is the one to use in the following step.
POD1a<-geeglm(DPM ~ AvgHrBasisMat + bs(PODAngle, knots=mean(PODAngle)) + X.TimeLost + TideBasisMat, family = binomial, corstr="independence", id=Panel, data=M7)

POD1b<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO))  + bs(PODAngle, knots=mean(PODAngle)) + X.TimeLost + TideBasisMat, family = binomial, corstr="independence", id=Panel, data=M7)

POD1c<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO)) + AvgHrBasisMat + X.TimeLost + TideBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

POD1d<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO)) + AvgHrBasisMat + bs(PODAngle, knots=mean(PODAngle))  + TideBasisMat, family = binomial, corstr="independence", id=Panel, data=M7)

POD1e<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO)) + AvgHrBasisMat + bs(PODAngle, knots=mean(PODAngle)) + X.TimeLost, family = binomial, corstr="independence", id=Panel, data=M7)

model1<-c("POD0", "POD1", "POD1a", "POD1b", "POD1c", "POD1d", "POD1e")
QIC1<-c(POD0@pan.aic, POD1@pan.aic, POD1a@pan.aic, POD1b@pan.aic, POD1c@pan.aic, POD1d@pan.aic, POD1e@pan.aic)
QICmod1<-data.frame(rbind(model1,QIC1))
QICmod1
#X1               X2               X3               X4
#model1             POD0             POD1            POD1a            POD1b
#QIC1   15401.3578427153 13357.7779995713 14548.9472235935 13363.3777457815
#X5               X6               X7
#model1            POD1c            POD1d            POD1e
#QIC1   13349.9646158432 13358.5697983154 13687.7385781282

#Remove PODANGLE:
POD2<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO)) + AvgHrBasisMat + X.TimeLost + TideBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

POD2a<-geeglm(DPM ~AvgHrBasisMat + X.TimeLost + TideBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

POD2b<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO))  + X.TimeLost + TideBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

POD2c<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO)) + AvgHrBasisMat  + TideBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

POD2d<-geeglm(DPM ~bs(DATENO , knots=mean(DATENO)) + AvgHrBasisMat + X.TimeLost  , family = binomial, corstr="independence", id=Panel, data=M7)

model2<-c("POD0", "POD2", "POD2a", "POD2b", "POD2c", "POD2d")
QIC2<-c(POD0@pan.aic, POD2@pan.aic, POD2a@pan.aic, POD2b@pan.aic, POD2c@pan.aic, POD2d@pan.aic)
QICmod2<-data.frame(rbind(model2,QIC2))
QICmod2
X1               X2               X3               X4
model2             POD0             POD2            POD2a            POD2b
QIC2   15401.3578427153 13349.9646158432 14607.0850138747 13353.3619129114
X5               X6
model2            POD2c            POD2d
QIC2   13350.6788475452 14220.9059321465

# Retain all other covariates.

# STEP 5: Testing covariate significance.
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

# In descending order: 
AvgHrBasisMat
X.TimeLost
bs(DATENO , knots=mean(DATENO))
TideBasisMat 

# So: 
POD5<-geeglm(DPM ~  bs(DATENO , knots=mean(DATENO)) + TideBasisMat  + AvgHrBasisMat + X.TimeLost  , family = binomial, corstr="independence", id=Panel, data=M7) 
# The significance of the covariates is tested using Wald's tests (i.e. the anova.geeglm function of the library geeglm). Non-significant covariates are removed one at a time.
anova(POD5)
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: DPM
#Terms added sequentially (first to last)

#Df     X2 P(>|Chi|)    
#bs(DATENO, knots = mean(DATENO))  4 92.535   < 2e-16 ***
  #TideBasisMat                      4 99.797   < 2e-16 ***
  #AvgHrBasisMat                     4 10.084   0.03904 *  
  #X.TimeLost                        1  1.643   0.19990    
#---
  #Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Remove X.TimeLost:
POD6<-geeglm(DPM ~  bs(DATENO , knots=mean(DATENO)) + TideBasisMat  + AvgHrBasisMat , family = binomial, corstr="independence", id=Panel, data=M7) 
# The significance of the covariates is tested using Wald's tests (i.e. the anova.geeglm function of the library geeglm). Non-significant covariates are removed one at a time.
anova(POD6)
#Analysis of 'Wald statistic' Table
#Model: binomial, link: logit
#Response: DPM
#Terms added sequentially (first to last)

#Df     X2 P(>|Chi|)    
#bs(DATENO, knots = mean(DATENO))  4 92.535   < 2e-16 ***
  #TideBasisMat                      4 99.797   < 2e-16 ***
  #AvgHrBasisMat                     4 10.084   0.03904 *  
  #---
  #Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Retain all covariates. This is the final model.

# STEP 6: Construction of the ROC curve     
pr <- predict(POD6, type="response")  
pred <- prediction(pr,M7$DPM) 
perf <- performance(pred, measure="tpr", x.measure="fpr")   
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5))
#This creates a ROC plot

# Choice of the best cut-off probability

y<-as.data.frame(perf@y.values)
x<-as.data.frame(perf@x.values)
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
write.table(d,"D:\\distancesM7d.txt")                                 # to write a table with the computed distances

# The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
# MAX d= 0.322334462026791 --> position 24662

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
alpha[24662,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve
[1] 0.01823497
# Best cutoff:     0.01823497
# This value can now be used to build the confusion matrix:

DATA<-matrix(0,77863 ,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim()) 

dim(M7)
[1] 77863    52
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:77863                                    # the first column is filled with an ID value that is unique for each row
DATA$Observed<-M7$DPM                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(POD6,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = 0.01823497)                                   # the identified cut-off must be used here

#Confusion matrix:

observed
predicted     1     0
1  1196 23464
0   372 52831

# The confusion matrix can then be transformed into percentages:

observed
predicted     1    0
1  76.3%  30.8%
  0  23.7%  69.2%
  
  # The area under the curve (auc) can also be used as an rough indication of model performance:
  
  auc <- performance(pred, measure="auc")
auc = 0.7957429

# STEP 7: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
# Re-run the final model:
POD8<-geeglm(DPM ~  bs(DATENO , knots=mean(DATENO)) + TideBasisMat  + AvgHrBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

TideBasis<-gam(DPM~s(MINUTE, bs="cc", k=6), fit=F, data=M7, family=binomial, knots=list(MINUTE=seq(0,1,length=6)))$X[,2:5]
AvgHrBasis<-gam(DPM~s(HOUR, bs="cc", k=6), fit=F, data=M7, family=binomial, knots=list(HOUR=seq(0,23,length=6)))$X[,2:5]

TideBasisMat<-as.matrix(TideBasis)
AvgHrBasisMat<-as.matrix(AvgHrBasis)

dimnames(TideBasisMat)<-list(NULL,c("TBM1", "TBM2", "TBM3", "TBM4"))
dimnames(AvgHrBasisMat)<-list(NULL,c("AHBM1", "AHBM2", "AHBM3", "AHBM4"))

POD8<-geeglm(DPM ~  bs(DATENO , knots=mean(DATENO)) + TideBasisMat  + AvgHrBasisMat , family = binomial, corstr="independence", id=Panel, data=M7)

#Probability of covariate #1: bs(DATENO , knots=mean(DATENO)):
BootstrapParameters1<-rmvnorm(10000, coef(POD8),summary(POD8)$cov.unscaled, method="chol")
start=2; finish=5; Variable=DATENO; xlabel="# Days since deployment"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(POD8)[,start:finish]*coef(POD8)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~s(PlottingVar1, bs="tp", k=6), fit=F, family=binomial, knots=list(PlottingVar1=seq(0,61,length=6)))$X[,2:5]
RealFit1<-Basis1%*%coef(POD8)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)
plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(0,61), main ="M7 covariate 1: # Days since deployment" , cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey", main = "Influence of DATENO")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(PlottingVar1)




#Probability of covariate #2: TideBasisMat:  
BootstrapParameters2<-rmvnorm(10000, coef(POD8),summary(POD8)$cov.unscaled, method="chol")
start=6; finish=9; Variable=MINUTE; xlabel="Tidal cycle (0 = 1 = Low Tide at Oban)"; ylabel="Probability"   
PlottingVar2<-seq(min(Variable), max(Variable), length=5000)
CenterVar2<-model.matrix(POD8)[,start:finish]*coef(POD8)[c(start:finish)]
BootstrapCoefs2<-BootstrapParameters2[,c(start:finish)]
Basis2<-gam(rbinom(5000,1,0.5)~s(PlottingVar2, bs="cc", k=6), fit=F, family=binomial, knots=list(PlottingVar2=seq(0,1,length=6)))$X[,2:5]
RealFit2<-Basis2%*%coef(POD8)[c(start:finish)]
RealFitCenter2<-RealFit2-mean(CenterVar2)
RealFitCenter2a<-inv.logit(RealFitCenter2)
BootstrapFits2<-Basis2%*%t(BootstrapCoefs2)
quant.func2<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis2<-apply(BootstrapFits2, 1, quant.func2)-mean(CenterVar2)
cis2a<-inv.logit(cis2)
plot(PlottingVar2,(RealFitCenter2a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(0,1),main ="M7 covariate 2: Tidal cycle" , cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar2,(cis2a[1,]),PlottingVar2,(cis2a[2,]), col="grey")
lines(PlottingVar2,(RealFitCenter2a),lwd=2, col=1)
rug(PlottingVar2)



#Probability of covariate #3: AvgHrBasisMat:
BootstrapParameters3<-rmvnorm(10000, coef(POD8),summary(POD8)$cov.unscaled, method="chol")
start=10; finish=13; Variable=HOUR; xlabel="Diel Hour"; ylabel="Probability"  
PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
CenterVar3<-model.matrix(POD8)[,start:finish]*coef(POD8)[c(start:finish)]
BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=6), fit=F, family=binomial, knots=list(PlottingVar2=seq(0,23,length=6)))$X[,2:5]
RealFit3<-Basis3%*%coef(POD8)[c(start:finish)]
RealFitCenter3<-RealFit3-mean(CenterVar3)
RealFitCenter3a<-inv.logit(RealFitCenter3)
BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
cis3a<-inv.logit(cis3)
plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(0,23), main ="M7 covariate 3: Diel Hour" , cex.lab = 1.5, cex.axis=1.5)    
segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
rug(PlottingVar3)



## Finished ## 
