#sent by Ally Rice on 8/19/2021, was written to deal with Pm data. Modified by MAZ as a starting point for temporal models for ch. 2

rm(list=ls())
graphics.off()

library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
library(nlme)
library(mgcv)
library(mgcViz)
library(parsedate)
library(geepack)
library(splines)
library(pracma)  #used for linspace and size 
library(mctest)
library(statmod)
library(tweedie)
library(lubridate)
library(gridBase)
library(grid)
library(tidyverse)       # because it literally does everything
# library(rjags)           # replacement for geeglm which is out of date
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(mvtnorm)         # to build the partial residual plots
library(gridExtra)       # to build the partial residual plots
library(SimDesign)
library(regclass)
library(ChemoSpecUtils)
library(car) #for ANOVA


####################################################lets try running model stuff for each count data file we have
site = 'Kauai'
infolder = file.path('E:/ch2/countdata',site)
outfolder = file.path('E:/ch2/gamResults',site)
allfiles = list.files(infolder,pattern = "*.csv",full.names = TRUE)


####make directory if needed
if (!dir.exists(outfolder)){
dir.create(outfolder)
}

##########run for each file 
for (i in seq_along(allfiles)){
#setup our file 
Daily <- read.table(file = allfiles[i],sep = ',',header = TRUE)

#extract name
namesh = basename(allfiles[i])
outnamesh = str_replace(namesh,'.csv','_gam.RData')
outname = file.path(outfolder,outnamesh)

#######if site is Kona, remove rows from bad deployments
if (site == 'Kona'){
  rmvdeps = c('HAWAII01','HAWAII02','Hawaii03','Hawaii13','Hawaii14','Hawaii15','Hawaii_K_25')
  #find the rows for bad deployments and get rid of them
for (ir in 1:size(rmvdeps,2)){
  badrow = c(grep(rmvdeps[ir],Daily$depsave))
  #only do this if not empty because otherwise it makes Daily completely empty 
  if (!isempty(badrow)){
  Daily = Daily[-badrow,]
  
}}}

#################get a diel column 
#can just use the actual hours values for diel stuff, just need to split them into 0-24. Which is easy enough.
hrofdaytemp = as.POSIXct((Daily$hours - 719529)*86400, origin = "1970-01-01", tz = "UTC")
Daily$hrofday = hour(round(hrofdaytemp,units = "hours"))


#get moon frac to % illuminated  
Daily$moonfrac = Daily$moonfrac*100 #set moonfrac to % illuminated

#get pres/absence column
Daily$pres = Daily$counts
Daily$pres[Daily$counts>0] = 1


########look at autocorrelation
acf1 = acf(Daily$pres,lag.max = 100)
acf(Daily$pres, ylim=c(0,0.5), xlim =c(0,100))

#run basic glm model to get residuals of model and acf of residuals for final binning
#from Alba: using residuals for acf allows for accounting for multiple variables in the autocorrelation.
#for my purposes, maybe this will make the stenellid version better 
basicMod<-glm(pres~
                bs(juldays)+
                hrofday+
                as.factor(year)
              ,data=Daily,family=binomial)

#anova to check for variable significance, summary to look at other stuff
anova(basicMod)
summary(basicMod)
acfRes = acf(residuals(basicMod), lag.max = 1000, ylim=c(0,0.5))

#find acf value below 0.1 to use
acuse = Position(function(x) x < 0.01,acfRes$acf)
#ex: autocorrelation drops below 0.1 after 2 hours, so clusters for id should be 2 observations long
#this is used to set id in geeglm(); each row in id is of observations that might be correlated,
#observations in different rows are assumed to be not correlated
Daily$id = 0
deps = data.frame(unique(Daily$depsave))
#set up blocks at deployment level so that no block spans multiple deployments 
for (idep in 1:nrow(deps)){
  nrowdep = c(grep(deps$unique.Daily.depsave[idep],Daily$depsave))
Daily$id[nrowdep] = round(linspace(1,size(nrowdep,2)/acuse,size(nrowdep,2)))
}


#save plots
plotname1 = str_replace(outname,'.RData','acf.png')
dev.copy(png,plotname1)
dev.off()

#########simple gam- how to account for autocorrelation??-> using corstruct() in gam
#after talking to Annebelle, who recommended paper, switching to multinomial logistic regression family of distributions. 
#more appropriate for proportions arising from count data, which is what bins/hr is 
#paper: Douma & Weedon 2018 (stored in Mendeley ch 2)
#tweedie was bad because assumes upper limit can be infinity, which it actually can't with my data 
#switching again to try poisson(). multinomial requires the number of linear predictors to equal k, which causes 
# a problem since number of categories isn't related to number of variables in this case... 
#for now, switching over to using pres/absence and binomial() family.. not as detailed but 
#additionally, after talking to Natalie, switching to geeglm() because it works better with autocorrelated data. 
#can see what she did in her code, but seems like "independence" corstructure is what people are using (used in Harbor porpoise paper)

#####make variance/covariance matrix to use for cyclic covariates- julday and hr of day
hrgam = gam(pres~s(hrofday, bs="cc", k=6), fit=F, data=Daily,
                        family=binomial, knots=list(HOUR=seq(0,23,length=6)))$X[,2:5]
julgam = gam(pres~s(juldays, bs="cc", k=6), fit=F, data=Daily,
              family=binomial, knots=list(days=seq(1,366,length=6)))$X[,2:5]
hrmat = as.matrix(hrgam)
julmat = as.matrix(julgam)

##data exploration and colinearity analsis- skip for now, but will need to do for full environmental models


#base model
mod0 = geeglm(pres~1,family=binomial,corstr="ar1",id=id,data = Daily)
##select linear or smooth for each variable by fitting various models with only one term included as either a linear or smooth
#example for now with moonfrac
#can do this with our cyclical data by doing julmat vs bs(juldays)
#moonfrac
#not running as factor because that wouldn't make sense 
mfmod1 = geeglm(pres ~ moonfrac, family = binomial, corstr="ar1", id=id,data = Daily)
mfmod2 = geeglm(pres ~ bs(moonfrac,knots=6), family = binomial, corstr="ar1",id=id,data = Daily)
modelmf<-c("mod0", "mfmod1", "mfmod2")
QICmf<-c(QIC(mod0)[1],QIC(mfmod1)[1],QIC(mfmod2)[1])
QICmodmf<-data.frame(rbind(modelmf,QICmf))
QICmodmf

#create an if statement to pull out final way to use moonfrac, so can still run it in loop
minloc = which.min(QICmf[2:tail(QICmf,n=1)])
rm("moonuse")
if (minloc==1){
#choose first model, include moonfrac as a non-smooth
  moonuse = Daily$moonfrac
}else if (minloc==2){
  #otherwise if second model is better, use smooth for variable
  moonuse = bs(Daily$moonfrac,knots = 6)
}
#when we run this with more stuff, could make above into a little function


##test a variety of models to see what variables are most significant in order to determine final variable order 
#full model
m1full = geeglm(pres~hrmat + moonuse + julmat + as.factor(year),
                family = binomial(),data = Daily,corstr="ar1",id=id)
#remove hrmat
m1a = geeglm(pres~moonuse + julmat + as.factor(year),
             family = binomial(),data = Daily,corstr="ar1",id=id)
#remove moonuse
m1b = geeglm(pres~hrmat + julmat + as.factor(year),
                family = binomial(),data = Daily,corstr="ar1",id=id)
#remove julmat
m1c = geeglm(pres~hrmat + moonuse + as.factor(year),
                family = binomial(),data = Daily,corstr="ar1",id=id)
#remove year
m1d = geeglm(pres~hrmat + moonuse + julmat,
                family = binomial(),data = Daily,corstr="ar1",id=id)

#get QICs into one and look at full 
model1<-c("mod0", "m1full", "m1a", "m1b", "m1c", "m1d")
QICs<-c(QIC(mod0)[1],QIC(m1full)[1],QIC(m1a)[1],QIC(m1b)[1],QIC(m1c)[1],QIC(m1d)[1])
QICmod1<-data.frame(rbind(model1,QICs))
QICmod1

#choose final order of variables based on above results
#The order in which the covariates enter the model is determined by the QIC score
#(the ones that, if removed, determine the biggest increase in QIC enter the model first).
#remove null model and full model
QICdiffs = QICs[-c(1,2)] - QICs[2] #get differences of all models by comparing them to the full model
#add a row for indexing purposes
QICdiffall = rbind(QICdiffs,letters[1:size(QICdiffs,2)])
#rank them correctly 
sortedQIC = QICdiffall[,order(QICdiffs)]
#put all variables into one structure together
allvars = list(a=hrmat,b=moonuse,c=julmat,d=Daily$year,pres=Daily$pres,id = Daily$id)
termlist = c("hourofday","moonfrac","julianday","year")

#do final model with variables in order from sortedQIC
#WILL NEED TO MODIFY THIS DIRECTLY WHEN HAVE MORE VARIABLES, or find a cleaner way 
m1a = geeglm(allvars[["pres"]]~allvars[[sortedQIC[2,4]]]+allvars[[sortedQIC[2,3]]]+
               allvars[[sortedQIC[2,2]]]+allvars[[sortedQIC[2,1]]],family = binomial(),
             corstr = "ar1",id = allvars[["id"]], data = allvars)

QICfinal = QIC(m1a)
m1asum = anova(m1a)
m1asum

# # summary + plots
# plot(m1a,pages=1,all.terms=T,scale=0,shade=TRUE,rug = TRUE,shade.col = 'gray',residuals=F)
# title(namesh)
# #save smooth plots  
plotname2 = str_replace(outname,'.RData','smooths.png')
# dev.copy(png,plotname2)
# dev.off()

#redo model with changed dimensions for matrices, otherwise it bugs up in the plotting
dimnames(hrmat) = list(NULL,c("hr1","hr2","hr3","hr4"))
dimnames(julmat) = list(NULL,c("ju1","ju2","ju3","ju4"))
allvars = list(a=hrmat,b=moonuse,c=julmat,d=Daily$year,pres=Daily$pres,id = Daily$id)

#redo final model with these names
m1a = geeglm(allvars[["pres"]]~allvars[[sortedQIC[2,4]]]+allvars[[sortedQIC[2,3]]]+
               allvars[[sortedQIC[2,2]]]+allvars[[sortedQIC[2,1]]],family = binomial(),
             corstr = "ar1",id = allvars[["id"]], data = allvars)

#######function for plotting covariate plot adapted from Pirotta type plots 
plot_covariate = function(mod,variable,varname,range,type,st,fn){
  ####use code from Pirotta to plot partials 
  library(boot)
  library(SimDesign)
  #Probability of covariate
  BootstrapParameters3<-rmvnorm(10000, coef(mod),summary(mod)$cov.unscaled)
  start=st; finish=fn; Variable=variable; xlabel=varname; ylabel="Probability"  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(mod)[,start:finish]*coef(mod)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  if (grepl("factor",type)){
    Basis3<-gam(rbinom(5000,1,0.5)~as.factor(PlottingVar3), 
                fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=6)))$X[,2:5]
  }else if (grepl("linear",type)){
    Basis3<-gam(rbinom(5000,1,0.5)~PlottingVar3, 
                fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=6)))$X[,2:5]
  }else if (grepl("smooth",type)){
    Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3,k=6), 
                fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=6)))$X[,2:5]
  }else if (grepl("cyclical",type)){
    Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=6), 
                fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=6)))$X[,2:5]
  }
  if (start!=finish){ #i.e. if theres multiple terms
    RealFit3<-Basis3%*%coef(mod)[c(start:finish)]
  }else{
    #repmat the coefficient to needed length
    coefuse = coef(mod)[c(start:finish)]
    coefmat = repmat(coefuse,1,size(Basis3,2))
    RealFit3 <- Basis3%*%t(coefmat)
  }
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter3a<-inv.logit(RealFitCenter3)
  if (start==finish){ #i.e. if theres only one term, repeat bootstrap
    #repeat bootstrap columns to be equal
    bootrep = repmat(BootstrapCoefs3,size(Basis3,2),1)
    BootstrapFits3 = Basis3%*%bootrep
    BootstrapFits3 = t(BootstrapFits3)
  }else{
    BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  }
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  cis3a<-inv.logit(cis3)
  plot(PlottingVar3,RealFitCenter3a, col = 1,type="l",ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=range, main =NULL , cex.lab = 1.5, cex.axis=1.5)    
  segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
  lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
  rug(PlottingVar3)
  
}


###plot hour of day
par(mfrow=c(2,2))

# plotbase = str_replace(plotname2,'smooths.png','')
# plotname3 = paste(plotbase,paste(termlist[1],'_partial.png',sep=""),sep="")
#need to use coeff() to get starts and finishes for all variables
coeforder = names(coef(m1a))

#since variable was a cyclic mat, pass in the ORIGINAL variable, pre-cycling
hrind = str_which(coeforder,"hr")
plot_covariate(m1a,Daily$hrofday,termlist[1],range = c(0,23),type="cyclical",st=min(hrind),fn=max(hrind))
title(main = namesh)
###plot moonuse
#figure out which one is moon
QICloc = str_which(sortedQIC[2,],"b")
moonind = str_which(coeforder,fixed(paste("[2, ",QICloc,"]",sep = "")))
plot_covariate(m1a,Daily$moonfrac,termlist[2],range=c(0,100),type="smooth",st=min(moonind),fn=max(moonind))
###plot jul
julind = str_which(coeforder,"ju")
plot_covariate(m1a,Daily$juldays,termlist[3],range=c(0,365),type="cyclical",st=min(julind),fn=max(julind))
###plot year
#figure out which one is year
QICloc2 = str_which(sortedQIC[2,],"d")
yearind = str_which(coeforder,fixed(paste("[2, ",QICloc2,"]",sep = "")))
plot_covariate(m1a,Daily$year,termlist[4],range=c(min(Daily$year),max(Daily$year)),type="factor",st=min(yearind),fn=max(yearind))

#saveplot
dev.copy(png,plotname2)
dev.off()

#save model
save(list = c("m1a","m1asum","acuse","Daily","basicMod","QICmodmf","QICmod1","QICfinal","termlist"),file=outname)

}


##############notes on understanding residual plots
# gam.check() will also generate four plots. Each of these gives a different way of looking at your model residuals. 
# These plots show the results from the original, poorly fit model. On the top-left is a Q-Q plot, which compares the model 
# residuals to a normal distribution. A well-fit model's residuals will be close to a straight line. On bottom left is a 
# histogram of residuals. We would expect this to have a symmetrical bell shape. On top-right is a plot of residual values. 
# These should be evenly distributed around zero. Finally, on the bottom-right is plot of response against fitted values. 
# A perfect model would form a straight line. We don't expect a perfect model, but we do expect the pattern to cluster around 
# the 1-to-1 line.

# test whether the basis dimension (basis size) for a smooth is adquate (Wood, 2017):
#- k-index: k-index is a measure of the remaining pattern in the residuals. The further below 1, the more likely it is that there is missed pattern 
#           left in the residuals.
#- p-value: calculated based on the distribution of the k-index after randomizing the order of the residuals. Low values may indicate k to low, especially
#           if the edf is close to k. But p-values can be low for reasons other than a too low k.


###########skipping this for now, because using corstruct() argument instead
# ####use acf blocks to rebin data to account for autocorrelation
# Daily.acf = list()
# rowidx = 1
# irowuse = 1
# for (irow in 1:nrow(Daily)){
#   #if irowuse is greater than irow, then skip this row
#   if (irowuse < irow){
#   #go to first row
#   firstid = Daily$id[irow]
#   thisid = which(Daily$id == firstid)
#   #get only consecutive values 
#   findgap = diff(thisid)
#   lastval = which(findgap>1)
#   if (!isempty(lastval)){
#   userows = thisid[1:lastval[1]]}
#   else{
#     userows = thisid
# }
#   #get our shortened data
#   Daily.acf$depsave[rowidx] = Daily$depsave[userows[1]]
#   #average of counts
#   Daily.acf$counts[rowidx] = ceil(mean(Daily$counts[userows])) #some counts have half bins due to duty cycle. 
#   #These should be rounded up as any half bin would be given a 1 for presence
#   #use first hour value
#   Daily.acf$hours[rowidx] = Daily$hours[userows[1]]
#   #get mean for moonfrac
#   Daily.acf$moonfrac[rowidx] = mean(Daily$moonfrac[userows])
#   #get first julday value
#   Daily.acf$juldays[rowidx] = Daily$juldays[userows[1]]
#   #get first hrperday value
#   Daily.acf$hrofday[rowidx] = Daily$hrofday[userows[1]]
#   #get year
#   Daily.acf$year[rowidx] = Daily$year[userows[1]]
#   #get diel value
#   Daily.acf$dielval[rowidx] = Daily$dielval[userows[1]]
#   #get id
#   Daily.acf$id[rowidx] = Daily$id[userows[1]]
#   
#   #skip other rows of this type and move to next id value
#   irowuse = tail(userows,n=1)+1
#   rowidx = rowidx + 1
# }}
#  


######### calcualte tweedie parameter- old vers. Alba has said she's had better luck with tw() than Tweedie()
# y<-Daily$counts
# mean(y,na.rm=TRUE)

# out<-tweedie.profile(y~1,p.vec=seq(1,2,by=0.1))
# tweedie.plot(seq(0, 120, length=120), mu=mean(y,na.rm=TRUE),xi=out$p.max, phi=out$phi.max)
# par(new=TRUE)
# hist(y)
# str(out)



# m1a<-gamm(form1,family=Tweedie(p=out$p.max),correlation = corAR1(form =~1|id),data=Daily)
# ###try using gam() instead with tw() family. K = number of categories-1
# #find number of categories
# Kmax = 12 # most bins any could have is 12, so total categories is 13 

