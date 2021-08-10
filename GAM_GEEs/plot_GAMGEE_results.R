 ##########################Graphing for Alba

library (rlang)
#devtools::install_git("https://github.com/lindesaysh/MRSea.git")  ###if you haven't got MRSea installed'
library(MRSea) 
library (splines)
library (carData)
library (car)
library (geepack)
library(MASS)

#load("G:/Shared drives/SOCAL_Sonar_Impact/Impact_analysis/4.Model_results/SOCAL_H/SOCAL_H_GEE_results.RData")
####Some functions

myviolin<-function(xvals,miny,maxy,rel=10){
###produces a violin plot of the the observed frequency of a particular variable
##Inputs: relevant variable and the min and max of the graph. 
  densx<-density(xvals)
  x <- c(densx$x[1], densx$x)
  y <- c(0, densx$y)
  y1<- y/max(y)*((maxy-miny)/rel)
 
  newy<-miny+y1
  data.frame(x=x,y=newy)
}

predgridmaker <- function (model1, thingstoshow,namethingstoshow ,otherstuff=otherstuff, otherstuffvalue=otherstuffvalue,axisLabel=axisLabel, maxy=maxyVal, maxviolin=maxViolinVal){
#########Produces an effect plot from a variable, typically plotted on the response scale with all variables at given levels. 
##Inputs
#model1 Model object
#thingstoshow: the relevant input data (to get a min and max)
#namethingstoshow: the name of the variable (for the xaxis label)
#otherstuff: the other variables in the model (I could probably automate this)
##otherstuffvalue  (the values of the other variables that will be set for the plot (typically mean values))
##maxy. Can be set to NA in which case the maxium of the upper 97.5% bound is used. Sets the limit of the y axis. 

  bootN = 1000
  interaction <- F    ####only works if there is not an interaction
  for (i in 1:length (otherstuff)){
      eval (parse (text= paste (otherstuff[i], "=", otherstuffvalue[i]))   )
      }
  if (is.factor (thingstoshow)==F){
     if (interaction==T){
###
   }else{
   thingstoshow2 <- seq(min (thingstoshow), max(thingstoshow), length=bootN)}}else{
   thingstoshow2 <-sort (unique (thingstoshow))
   }
   phrase1 <- c("expand.grid(thingstoshow2")
   for (j in 1: length (otherstuff)){
       phrase1 <- paste (phrase1,otherstuff[j], sep=",")
       if (j==length (otherstuff)){
          phrase1 <- paste (phrase1,")")
       }
   }
   names2 <- c("thingstoshow2", otherstuff)
   gridpred = eval(parse (text=phrase1)) ######creates a grid to predict over
   names (gridpred)[1] <- namethingstoshow
   names (gridpred)[2:dim (gridpred)[2]]<- otherstuff
   gridpred$Pred <- predict (POD3, gridpred, type="response")
   bspreds <- matrix (NA, bootN, dim (gridpred)[1])
   bootstrapmodel <- model1
   set.seed (101)
   #############Bootstrapping to obtain confidence interval 
   rcoefs <- try(MASS::mvrnorm(bootN, coef(model1), summary(model1)$cov.scaled), silent = T)
   if (is.null(rcoefs) || length(which(is.na(rcoefs) ==  T)) > 0) {
       rcoefs <- MASS::mvrnorm(bootN, coef(model1), as.matrix(nearPD(summary(model1)$cov.scaled)$mat))
            }   #####bootstrap replicate of model coefficient. 
   for (i in 1:bootN){ 
       bootstrapmodel$coefficients <- rcoefs[i,]    
       bspreds[i,] <- predict (bootstrapmodel, gridpred, type="response")
   }
  lower <- apply (bspreds, 2, quantile, 0.025) 
  upper <- apply (bspreds, 2, quantile, 0.975)  
  if (is.na(maxy)==T){maxy=max(upper, na.rm=T)}
  if (is.factor (thingstoshow)==F){
     if (thingstoshow=="newjd"){tempjd <- ifelse (gridpred[,1]<0, gridpred[,1]+365,gridpred[,1] )
        plot (tempjd[1:489], gridpred$Pred[1:489], type="l" ,ylab="Probability" , ylim=c(min (lower), maxy), xlim=c(0,366), xlab="Dayofyear")
       lines (tempjd[490:1000], gridpred$Pred[490:1000])
        lines(tempjd[1:489] , lower[1:489], col = "red", lty = 2)
        lines(tempjd[490:1000] , lower[490:1000], col = "red", lty = 2)
       lines(tempjd[490:1000] , upper[490:1000], col = "red", lty = 2)  
        lines(tempjd[1:489] , upper[1:489], col = "red", lty = 2) 
     }else {
        plot (gridpred[,1], gridpred$Pred, type="l" ,ylab="Probability" ,xlab = axisLabel, ylim=c(min (lower), maxy) )
        lines(gridpred[,1] , lower, col = "orange", lty = 2)
        lines(gridpred[,1] , upper, col = "orange", lty = 2) } }
      
     
     else{
     plot (gridpred[,1], gridpred$Pred, pch=20 ,ylab="Probability" ,xlab = axisLabel, ylim=c(min (lower),maxy) )
    arrows(seq(1, length (gridpred[,1])), gridpred$Pred  , seq(1, length (gridpred[,1])),   lower, angle=90, col = "orange", lty = 2)
    arrows(seq(1, length (gridpred[,1])) ,gridpred$Pred,  seq(1, length (gridpred[,1])),  upper, angle=90, col = "orange", lty = 2) 
    }
  
    ##adds violin plot
    if (is.factor(thingstoshow)==F){
    if (thingstoshow=="newjd"){thingstoshow<- "jd"}###plots on original scale. 
    violinrug<-myviolin(xvals = thingstoshow,miny = min(lower,na.rm=T),maxy = maxviolin,rel=5)
    polygon(violinrug,col=rgb(0,1,1,alpha=.3)) }
     return (gridpred )             
}

#When there is only one thing in the model
predgridmakerONE <- function (model1, thingstoshow,namethingstoshow, axisLabel=axisLabel, maxy=maxyVal, maxviolin=maxViolinVal){
  #########Produces an effect plot from a variable, typically plotted on the response scale with all variables at given levels. 
  ##Inputs
  #model1 Model object
  #thingstoshow: the relevant input data (to get a min and max)
  #namethingstoshow: the name of the variable (for the xaxis label)
  ##maxy. Can be set to NA in which case the maxium of the upper 97.5% bound is used. Sets the limit of the y axis. 
  
  bootN = 1000
  interaction <- F    ####only works if there is not an interaction
  thingstoshow = SiteHourTableB$Julian
  if (is.factor (thingstoshow)==F){
    if (interaction==T){
      ###
    }else{
      thingstoshow2 <- seq(min (thingstoshow), max(thingstoshow), length=bootN)}}else{
        thingstoshow2 <-sort (unique (thingstoshow))
      }
  
  phrase1 <- c("expand.grid(thingstoshow2)")

  names2 <- c("thingstoshow2")
  
  gridpred <- eval(parse (text=phrase1)) ######creates a grid to predict over
  namethingstoshow = namethingstoshow= "jd"
  names (gridpred)[1] <- namethingstoshow
  
  gridpred$Pred <- predict(POD3, gridpred, type="response")
  bspreds <- matrix (NA, bootN, dim (gridpred)[1])
  bootstrapmodel <- POD3
  set.seed (101)
  #############Bootstrapping to obtain confidence interval 
  rcoefs <- try(MASS::mvrnorm(bootN, coef(model1), summary(model1)$cov.scaled), silent = T)
  if (is.null(rcoefs) || length(which(is.na(rcoefs) ==  T)) > 0) {
    rcoefs <- MASS::mvrnorm(bootN, coef(model1), as.matrix(nearPD(summary(model1)$cov.scaled)$mat))
  }   #####bootstrap replicate of model coefficient. 
  for (i in 1:bootN){ 
    bootstrapmodel$coefficients <- rcoefs[i,]    
    bspreds[i,] <- predict (bootstrapmodel, gridpred, type="response")
  }
  lower <- apply (bspreds, 2, quantile, 0.025) 
  upper <- apply (bspreds, 2, quantile, 0.975)  
  if (is.na(maxy)==T){maxy=max(upper, na.rm=T)}
  if (is.factor (thingstoshow)==F){
    if (thingstoshow=="newjd"){tempjd <- ifelse (gridpred[,1]<0, gridpred[,1]+365,gridpred[,1] )
    plot (tempjd[1:489], gridpred$Pred[1:489], type="l" ,ylab="Probability" , ylim=c(min (lower), maxy), xlim=c(0,366), xlab="Dayofyear")
    lines (tempjd[490:1000], gridpred$Pred[490:1000])
    lines(tempjd[1:489] , lower[1:489], col = "red", lty = 2)
    lines(tempjd[490:1000] , lower[490:1000], col = "red", lty = 2)
    lines(tempjd[490:1000] , upper[490:1000], col = "red", lty = 2)  
    lines(tempjd[1:489] , upper[1:489], col = "red", lty = 2) 
    }else {
      plot (gridpred[,1], gridpred$Pred, type="l" ,ylab="Probability" ,xlab = axisLabel, ylim=c(min (lower), maxy) )
      lines(gridpred[,1] , lower, col = "orange", lty = 2)
      lines(gridpred[,1] , upper, col = "orange", lty = 2) } }
  
  
  else{
    plot (gridpred[,1], gridpred$Pred, pch=20 ,ylab="Probability" ,xlab = axisLabel, ylim=c(min (lower),maxy) )
    arrows(seq(1, length (gridpred[,1])), gridpred$Pred  , seq(1, length (gridpred[,1])),   lower, angle=90, col = "orange", lty = 2)
    arrows(seq(1, length (gridpred[,1])) ,gridpred$Pred,  seq(1, length (gridpred[,1])),  upper, angle=90, col = "orange", lty = 2) 
  }
  
  ##adds violin plot
  if (is.factor(thingstoshow)==F){
    if (thingstoshow=="newjd"){thingstoshow<- "jd"}###plots on original scale. 
    violinrug<-myviolin(xvals = thingstoshow,miny = min(lower,na.rm=T),maxy = maxviolin,rel=5)
    polygon(violinrug,col=rgb(0,1,1,alpha=.3)) }
  return (gridpred )             
}
 


windows() ###### Plots timeofd
otherstuff = c("sPres","sLag","sProp","maxRLpp","jd", "year")
otherstuffvalue <- c(1,0,mean (zcdata$sProp), mean (zcdata$maxRLpp),183,2011)       ###perhaps better to set jd=183
namethingstoshow= "timeofd"
axisLabel = "Normalized time of day"
predgridmaker (zc.gee1, zcdata$timeofd, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.003, maxviolin=0.003)

windows() ##### Plots maxRLpp
otherstuff = c("sPres","sProp", "sLag", "jd", "timeofd", "year")
otherstuffvalue <- c(1,mean (zcdata$sProp), 0, 183, 0, 2011)      ###check these values are suitable 
namethingstoshow="maxRLpp"
axisLabel = "Max. RLpp (dB)"
predgridmaker (zc.gee1, zcdata$maxRLpp[zcdata$sPres==1], namethingstoshow, otherstuff, otherstuffvalue,axisLabel, 0.02, 0.02) 

windows() ##### Plots sProp
otherstuff = c("sPres","maxRLpp", "sLag", "jd", "timeofd", "year")
otherstuffvalue <- c(1,mean (zcdata$maxRLpp), 0, 183, 0,2011)       ###perhaps better to set jd=183
namethingstoshow="sProp"
axisLabel = "Sonar prop./min"
predgridmaker (zc.gee1, zcdata$sProp[zcdata$sPres==1], namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.015, maxviolin=0.015)

windows() ###### Plots sLag
otherstuff = c("sPres","sProp","maxRLpp", "jd", "timeofd","sunsetLag", "year")
otherstuffvalue <- c(1,mean (zcdata$sProp), mean (zcdata$maxRLpp), 183, 0,2011)       ###perhaps better to set jd=183
namethingstoshow= "sLag"
axisLabel = "Sonar lag (days)"
predgridmaker (zc.gee3, zcdata$sLag, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.010, maxviolin=0.010)          


################Plotting 
     #####Think hard about the appropriate level for the other variables

  #####new model maxppl. Not quite the same as the poster. 
  ###The violin plot could be fiddled with

windows() ##### Plots maxRLpp
otherstuff = c("sPres","sProp", "sLag", "jd", "sunriseLag","sunsetLag", "year")
otherstuffvalue <- c(1,mean (zcdata$sProp), 0, 183, 12, 12, 2011)      ###check these values are suitable 
namethingstoshow="maxRLpp"
axisLabel = "Max. RLpp (dB)"
predgridmaker (zc.gee3, zcdata$maxRLpp[zcdata$sPres==1], namethingstoshow, otherstuff, otherstuffvalue,axisLabel, 0.02, 0.02)  
      
windows() ##### Plots sProp
otherstuff = c("sPres","maxRLpp", "sLag", "jd", "sunriseLag","sunsetLag", "year")
otherstuffvalue <- c(1,mean (zcdata$maxRLpp), 0, 183, 12,12,2011)       ###perhaps better to set jd=183
namethingstoshow="sProp"
axisLabel = "Sonar prop./min"
predgridmaker (zc.gee3, zcdata$sProp[zcdata$sPres==1], namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.015, maxviolin=0.015)

windows() ###### Plots sLag
otherstuff = c("sPres","sProp","maxRLpp", "jd", "sunriseLag","sunsetLag", "year")
otherstuffvalue <- c(1,mean (zcdata$sProp), mean (zcdata$maxRLpp), 183, 12,12,2011)       ###perhaps better to set jd=183
namethingstoshow= "sLag"
axisLabel = "Sonar lag (days)"
predgridmaker (zc.gee3, zcdata$sLag, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.010, maxviolin=0.010)          
 
windows() ###### Plots jd
otherstuff = c("Year")
otherstuffvalue <- c(2015)       ###perhaps better to set jd=183
namethingstoshow= "jd"
axisLabel = "Julian Date"
predgridmaker (POD3, SiteHourTable$Julian, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.015, maxviolin = 0.015)
 
windows() ###### Plots sunriseLag
otherstuff = c("sPres","sLag","sProp","maxRLpp","jd","sunsetLag", "year")
otherstuffvalue <- c(1,0,mean (zcdata$sProp), mean (zcdata$maxRLpp),183, 12,2011)       ###perhaps better to set jd=183
namethingstoshow= "sunriseLag"
axisLabel = "Hours since sunrise"
predgridmaker (zc.gee3, zcdata$sunriseLag, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.010, maxviolin=0.010) 

windows() ###### Plots sunriseLag
otherstuff = c("sPres","sLag","sProp","maxRLpp","jd","sunriseLag", "year")
otherstuffvalue <- c(1,0,mean (zcdata$sProp), mean (zcdata$maxRLpp),183, 12,2011)       ###perhaps better to set jd=183
namethingstoshow= "sunsetLag"
axisLabel = "Hours since sunset"
predgridmaker (zc.gee3, zcdata$sunsetLag, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.010, maxviolin=0.010) 
 
windows() ###### Plots year
otherstuff = c("sPres","sLag","sProp","maxRLpp","jd", "sunriseLag","sunsetLag")
otherstuffvalue <- c(1,0,mean (zcdata$sProp), mean (zcdata$maxRLpp),183,12,12)       ###perhaps better to set jd=183
namethingstoshow= "year"
axisLabel = "Year"
predgridmaker (zc.gee3, zcdata$year, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=NA, maxviolin=0.02) 
 
 
windows() ###### Plots timeofd
otherstuff = c("sPres","sLag","sProp","maxRLpp","jd", "year")
otherstuffvalue <- c(1,0,mean (zcdata$sProp), mean (zcdata$maxRLpp),183,2011)       ###perhaps better to set jd=183
namethingstoshow= "timeofd"
axisLabel = "Time of day (hours)"
predgridmaker (zc.gee3, zcdata$sunsetLag, namethingstoshow, otherstuff, otherstuffvalue,axisLabel, maxy=0.02, maxviolin=0.02)      
      