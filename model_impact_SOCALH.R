
##########################################################################################################################
### CORRELATION PLOT

library(corrplot) #correlation map
#devtools::install_git("https://github.com/lindesaysh/MRSea.git")  ###if you haven't got MRSea installed'
library(MRSea) 
library (splines) #required for bs()
library (car)
library (geepack)
library(splines2)

getPvalues<-function(model, varlist=NULL, factorlist=NULL){
  #-----------------------------------------------------------------------------
  #' Calculate marginal p-values from a \code{model}. 
  #' 
  #' An ANOVA is fitted repeatedly with each covariate being the last so that the output is marginal.  \code{varlist} and \code{factorlist} are optional and shorten the variable names in the output.
  #' 
  #' @param model Fitted model object of class \code{gee}.
  #' @param varlist (default =\code{NULL}). Vector of covariate names (continous covariates only) used to make the output table names shorter.  Useful if spline parameters are specified in the model.
  #' @param factorlist (default =\code{NULL}). Vector of covariate names (factor covariates only) used to make the output table names shorter. Useful if spline parameters are specified in the model.
  #' 
  #' @return
  #' Print out table of each variable and its associated marginal p-value.
  #' 
  #' @examples 
  #' 
  #' # load data
  #' data(ns.data.re)
  #' 
  #' # make blocking structure
  #' ns.data.re$blockid<-paste(ns.data.re$GridCode, ns.data.re$Year, ns.data.re$MonthOfYear, 
  #'                     ns.data.re$DayOfMonth, sep='')
  #' ns.data.re$blockid<-as.factor(ns.data.re$blockid)
  #' 
  #' initialModel<- geeglm(birds ~ as.factor(floodebb) + as.factor(impact) + observationhour + x.pos + 
  #'               y.pos + offset(log(area)), family='poisson',data=ns.data.re, id=blockid)
  #' 
  #' getPvalues(initialModel, varlist=c('observationhour', 'x.pos', 'y.pos'), 
  #'             factorlist=c('floodebb', 'impact'))
  #' 
  #' getPvalues(initialModel)
  #' 
  #' @export
  #' 
  if(class(model)[1]!='geeglm'){stop('Class of model not geeglm. Fit model as GEE or use an appropriate ANOVA function (e.g. Anova or anova)')}
  
  print("Getting marginal p-values")
  
  # make list of terms if varlist and factorlist are spcified
  if(is.null(varlist)==FALSE){
    termlist<-rep(NA, length=length(labels(terms(model))))
    cressterm<-grep('LocalRadial', labels(terms(model)))
    
    if(length(varlist)+length(factorlist) + length(cressterm)!=length(termlist)) stop('not all one dimensional terms specified in varlist or factorlist')
    
    if(length(varlist)>0){
      for(v in 1:length(varlist)){
        termlist[grep(varlist[v], labels(terms(model)))]<- varlist[v]  
      }
    }
    
    if(length(factorlist)>0){
      for(f in 1:length(factorlist)){
        termlist[grep(factorlist[f], labels(terms(model)))]<- factorlist[f]  
      }
    }
    
    if(length(cressterm)>0){
      for(a in 1:length(cressterm)){
        if(is.na(termlist[cressterm[a]])==FALSE){
          termlist[cressterm[a]]<- paste('s(x.pos, y.pos):', termlist[cressterm[a]], sep='')
        }else{termlist[cressterm[a]]<- 's(x.pos, y.pos)'  }
      }  
    }
    
  }
  
  
  #get marginal p-values for each term:
  store<- matrix(0, ncol=2, nrow=length(labels(terms(model))))
  for(i in 1:length(labels(terms(model)))){
    covariate<- labels(terms(model))[i]  
    
    text1<-paste("anova(update(model, . ~ . - ",  covariate, "+",  covariate, "))", sep="")
    test<-eval(parse(text=text1))
    
    if(is.null(varlist)){
      store[i,1] <- covariate
    }else{      
      store[i,1] <- termlist[i]
    }
    
    store[i,2]<- round(test[which(labels(test)[[1]]==covariate),3],6)
    if(as.numeric(store[i,2])<0.0001){store[i,2]<-'<0.0001'}
  }
  
  # check for zero p-values
  which(store[,2]==0)
  #print results with p-values
  store<- as.data.frame(store)
  names(store)<-c("Variable", "p-value")
  print(store)
  
}

plotCols = c(6, 10, 12, 14:21)
covarList = names(binData[c(plotCols)])
varUnits = c("Sonar presence","Sonar lag (min)","Min. since sunrise","Min. since sunset","Day-Night","Julian date",
             "Year","Time of day (min)","Sonar prop.","Cum. SEL","Max. RLpp")

names(varUnits) = covarList

# omit nan rows
data = na.omit(binData)

## create correlation matrix
corData = cor(data[covarList])

## plot correlation map in AOE order
windows()
corrplot(corData,type = "lower", order = "AOE",tl.col = "black",tl.srt=45,diag = FALSE)

###########################################################################################################################
#Eliminate NA
zcdata<- binData[is.na(binData$sLag)==F,] 

###########################################################################################################################
### SELECT VARIABLES: GLM & VIF
# Apply a basic GLM assuming all points are independent. Calculate variance inflation factor (VIF) to assess which covariates are collinear. 
# To find which ones do not contain collinearity, variables are removed one at a time recalculating VIF values. A VIF value above 5 indicates 
# high correlation and is cause for concern.

zc.glm<-glm(zcPres~
              DN+
              bs(sLag)+
              bs(sProp)+
              bs(jd)+
              bs(timeofd)+
              bs(sunriseLag)+
              bs(sunsetLag)+
              as.factor(year)+
              sPres:bs(maxRLpp)+
              sPres:bs(cumSEL)
            ,data=binData,family=binomial)

summary(zc.glm)

###Sequentially remove terms with a high vif. >5

vif1 <- 5000
bitsmain <- attributes(terms (zc.glm))$term.labels
threshVIF = 5
while (vif1>sqrt(threshVIF)){
  
  
  for (x in 1: length (bitsmain)){
    if (x==1){newformula <- bitsmain[1]}else{newformula <- paste (newformula, bitsmain[x], sep="+")}
  }
  
  newformula <- paste ("zcPres~", newformula,sep="")
  modeltemp <- glm (as.formula (newformula), data=zcdata,family=binomial) 
  
  viftemp <- vif(modeltemp)
  print (viftemp)
  
  vif2 <- match (max(viftemp[,3]), viftemp[,3])
  vif1 <- viftemp[vif2, 3]
  if (vif1>sqrt(threshVIF)){   bitsmain <- bitsmain[-vif2] }
}

bitsmain ##so one can see what kept

### RESULTS
# GVIF Df GVIF^(1/(2*Df))
# DN                1.241908e+03  1       35.240720
# bs(sLag)          1.780155e+00  3        1.100888
# bs(sProp)         1.341097e+01  3        1.541381
# bs(jd)            2.612927e+01  3        1.722614
# bs(timeofd)       2.100564e+03  3        3.578678
# bs(sunriseLag)    1.915467e+04  3        5.172642
# bs(sunsetLag)     6.909152e+04  3        6.405762
# as.factor(year)   1.957418e+00  9        1.038017
# sPres:bs(maxRLpp) 1.734001e+08  3       23.614291
# sPres:bs(cumSEL)  2.120660e+08  3       24.419972
# GVIF Df GVIF^(1/(2*Df))
# bs(sLag)          1.545448e+00  3        1.075249
# bs(sProp)         1.339930e+01  3        1.541157
# bs(jd)            1.675834e+00  3        1.089863
# bs(timeofd)       2.044564e+03  3        3.562598
# bs(sunriseLag)    1.342753e+02  3        2.262903
# bs(sunsetLag)     6.405376e+02  3        2.936009
# as.factor(year)   1.805425e+00  9        1.033367
# sPres:bs(maxRLpp) 1.760661e+08  3       23.674418
# sPres:bs(cumSEL)  2.153219e+08  3       24.482065
# GVIF Df GVIF^(1/(2*Df))
# bs(sLag)             1.545499  3        1.075255
# bs(sProp)           10.005826  3        1.467942
# bs(jd)               1.675473  3        1.089824
# bs(timeofd)       2044.941024  3        3.562707
# bs(sunriseLag)     134.276645  3        2.262907
# bs(sunsetLag)      640.585210  3        2.936046
# as.factor(year)      1.803803  9        1.033315
# sPres:bs(maxRLpp)   10.121230  3        1.470750
# GVIF Df GVIF^(1/(2*Df))
# bs(sLag)           1.545483  3        1.075253
# bs(sProp)         10.005828  3        1.467942
# bs(jd)             1.672147  3        1.089463
# bs(sunriseLag)    19.980557  3        1.647282
# bs(sunsetLag)     20.046628  3        1.648189
# as.factor(year)    1.803758  9        1.033314
# sPres:bs(maxRLpp) 10.120652  3        1.470736

###variables kept: 
#[1] "bs(sLag)"          "bs(sProp)"         "bs(jd)"            "bs(sunriseLag)"    "bs(sunsetLag)"     "as.factor(year)"   "sPres:bs(maxRLpp)"
#removed variables: DN, cumSEL and timeofd

zc.glm2<-glm(zcPres~
              bs(sLag)+
              bs(sProp)+
              bs(jd)+
              bs(sunriseLag)+
              bs(sunsetLag)+
              as.factor(year)+
              sPres:bs(maxRLpp)
            ,data=binData,family=binomial)

summary(zc.glm2)

# Call:
#   glm(formula = zcPres ~ bs(sLag) + bs(sProp) + bs(jd) + bs(sunriseLag) + 
#         bs(sunsetLag) + as.factor(year) + sPres:bs(maxRLpp), family = binomial, 
#       data = binData)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -0.4796  -0.1791  -0.1486  -0.1248   3.7994  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)         -5.59301    0.06416 -87.173  < 2e-16 ***
#   bs(sLag)1           -0.50945    0.05488  -9.284  < 2e-16 ***
#   bs(sLag)2           -0.14519    0.10150  -1.431  0.15257    
# bs(sLag)3            1.03885    0.09408  11.042  < 2e-16 ***
#   bs(sProp)1           1.06595    0.36935   2.886  0.00390 ** 
#   bs(sProp)2          -3.36315    0.83858  -4.011 6.06e-05 ***
#   bs(sProp)3           2.32619    1.70682   1.363  0.17292    
# bs(jd)1              2.86064    0.06464  44.258  < 2e-16 ***
#   bs(jd)2             -2.05582    0.03773 -54.493  < 2e-16 ***
#   bs(jd)3              1.63093    0.03331  48.963  < 2e-16 ***
#   bs(sunriseLag)1      2.35376    0.10397  22.639  < 2e-16 ***
#   bs(sunriseLag)2      0.22253    0.07990   2.785  0.00535 ** 
#   bs(sunriseLag)3      0.52235    0.03871  13.494  < 2e-16 ***
#   bs(sunsetLag)1       0.10384    0.09820   1.057  0.29032    
# bs(sunsetLag)2       0.56544    0.07765   7.282 3.29e-13 ***
#   bs(sunsetLag)3      -0.13196    0.03311  -3.986 6.73e-05 ***
#   as.factor(year)2008 -0.66071    0.04388 -15.056  < 2e-16 ***
#   as.factor(year)2009 -0.45931    0.04046 -11.351  < 2e-16 ***
#   as.factor(year)2010 -0.09523    0.03879  -2.455  0.01409 *  
#   as.factor(year)2011 -0.22442    0.03872  -5.797 6.77e-09 ***
#   as.factor(year)2012  0.08015    0.04698   1.706  0.08798 .  
# as.factor(year)2013 -0.72872    0.04635 -15.721  < 2e-16 ***
#   as.factor(year)2014 -0.61838    0.03964 -15.599  < 2e-16 ***
#   as.factor(year)2015 -0.37612    0.03937  -9.553  < 2e-16 ***
#   as.factor(year)2016  0.49162    0.03832  12.830  < 2e-16 ***
#   sPres:bs(maxRLpp)1  -8.54791    1.49600  -5.714 1.10e-08 ***
#   sPres:bs(maxRLpp)2   9.15894    1.33733   6.849 7.45e-12 ***
#   sPres:bs(maxRLpp)3  -7.80483    0.85248  -9.155  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 426212  on 3082667  degrees of freedom
# Residual deviance: 415084  on 3082640  degrees of freedom
# (131005 observations deleted due to missingness)
# AIC: 415140
# 
# Number of Fisher Scoring iterations: 8


##########################################################################################################

##########################################################################################################
### IDENTIFY BLOCKS OF DATA: check residuals & apply Autocorrelation function

Anova(zc.glm2)  #calculates p values using type 3 sums of squares so each p represents the term as placed last in the model.

# Analysis of Deviance Table (Type II tests)
# 
# Response: zcPres
# LR Chisq Df Pr(>Chisq)    
# bs(sLag)             336.7  3  < 2.2e-16 ***
#   bs(sProp)             30.6  3  1.045e-06 ***
#   bs(jd)              4028.5  3  < 2.2e-16 ***
#   bs(sunriseLag)       841.5  3  < 2.2e-16 ***
#   bs(sunsetLag)        124.2  3  < 2.2e-16 ***
#   as.factor(year)     5053.6  9  < 2.2e-16 ***
#   sPres:bs(maxRLpp)    196.9  3  < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Evaluate residuals with the auto correlation function
windows()
acf (residuals (zc.glm2))

#plot zoomed in
acf (residuals (zc.glm2),ylim=c(0,0.05),xlim=c(0,50), lag.max=100)

## Create blocks for GEEs
diff <- zcdata$Startnum[2:dim (zcdata)[1]]-zcdata$Startnum[1:c(dim (zcdata)[1]-1)]
diff <- round (diff, 6)
zcdata$Diff <- c(NA, diff)
InitialID <- 1      ####Panel ID for each continguous bit of data
for (i in 2:dim (zcdata)[1]){
  if (zcdata$Diff[i]>0.000694){InitialID[i] <- InitialID[i-1]+1}else{InitialID[i] <- InitialID[i-1]}
} 
#######
NoInitPanels <- length (unique (InitialID))
for (i in 1:NoInitPanels){
  temp <- zcdata[InitialID==i,]
  timesincestart <- (temp$Startnum-temp$Startnum[1])
  ID40 <-  i*10000 +floor (timesincestart/0.000694*40 )
  DayID <- i*10000 +ceiling (timesincestart)
  WeekID  <- i*10000 +ceiling (timesincestart/7)
  
  if (i==1){finalDayID <- DayID
  finalID40 <- ID40
  
  }else{finalDayID <- c(finalDayID, DayID)
  finalID40<- c(finalID40, ID40)
  }
  
}
zcdata$ID40 <- finalID40

#########################################################################################################
### GEE MODEL SELECTION

# transform variables
#zcdata$newjd <- ifelse (zcdata$jd>186,zcdata$jd-365,zcdata$jd)    ####not quite sure why I did this. 
zcdata$year <- as.factor (zcdata$year)
zcdata$sLag = zcdata$sLag/1440
zcdata$sunriseLag = zcdata$sunriseLag/60
zcdata$sunsetLag = zcdata$sunsetLag/60

quantile(zcdata$jd,c(0.25,0.50,0.75))
quantile(zcdata$sunriseLag,c(1/3,2/3))
quantile(zcdata$sunsetLag,c(1/3,2/3))

zc.gee1 = geeglm(zcPres~
                   bs(sLag)+
                   bs(sProp)+
                   bs(jd)+
                   bs(sunriseLag)+
                   bs(sunsetLag)+
                   as.factor(year)+
                   sPres:bs(maxRLpp)
                 ,data=zcdata,family=binomial,id=ID40) 

getPvalues(zc.gee1)
# [1] "Getting marginal p-values"
# Variable p-value
# 1          bs(sLag) <0.0001
# 2         bs(sProp) <0.0001
# 3            bs(jd) <0.0001
# 4    bs(sunriseLag) <0.0001
# 5     bs(sunsetLag) <0.0001
# 6   as.factor(year) <0.0001
# 7 sPres:bs(maxRLpp) <0.0001
# QICb: 415503

zc.gee3 = geeglm(zcPres~
                   bs(sLag)+
                   bs(sProp)+
                   mSpline(jd,knots = c(80, 160, 240), Boundary.knots = c(1,365),periodic = T)+
                   bs(sunriseLag)+
                   bs(sunsetLag)+
                   as.factor(year)+
                   sPres:bs(maxRLpp)
                 ,data=zcdata,family=binomial,id=ID40) 

getPvalues(zc.gee3)
# [1] "Getting marginal p-values"
# Variable p-value
# 1                                                                       bs(sLag) <0.0001
# 2                                                                      bs(sProp) <0.0001
# 3 mSpline(jd, knots = c(80, 160, 240), Boundary.knots = c(1, 365), periodic = T) <0.0001
# 4                                                                 bs(sunriseLag) <0.0001
# 5                                                                  bs(sunsetLag) <0.0001
# 6                                                                as.factor(year) <0.0001
# 7                                                              sPres:bs(maxRLpp) <0.0001
# QICb: 416077

zc.gee2 = geeglm(zcPres~
                   bs(sLag)+
                   bs(sProp)+
                   mSpline(jd,knots = c(160, 240), Boundary.knots = c(1,365),periodic = T)+
                   mSpline(sunriseLag,knots = c(8,16), Boundary.knots = c(0,24),periodic = T)+
                   mSpline(sunsetLag,knots = c(12,19), Boundary.knots = c(0,24),periodic = T)+
                   as.factor(year)+
                   sPres:bs(maxRLpp)
                 ,data=zcdata,family=binomial,id=ID40) 