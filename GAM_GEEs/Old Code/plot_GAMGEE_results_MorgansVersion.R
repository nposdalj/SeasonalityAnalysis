###violin function - from Alba
#notes from Morgan
#note that this version uses k=4 throughout instead of k=6
#this did cause a weird error in the plotting stuff for some reason where dimensions weren't matching up properly for cyclic covariates. 
#so the input coefficients from the model for that variable would have 4 terms, but the 'Basis3' variable would only have 2... 
#I worked around this by just using a subset of the coefficients from the model in the plotting BUT I'm not sure how sketchy that is to do
#long story short, if you're sticking with k=6 elsewhere, then change it back in here too and it just won't hit those if statements so the other stuff shouldn't come into play

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
#######function for plotting covariate plot adapted from Pirotta type plots
plot_covariate = function(mod,variable,varname,range,type,st,fn,data){
  ####use code from Pirotta to plot partials
  library(boot)
  library(SimDesign)
  start=st; finish=fn; Variable=variable; xlabel=varname; ylabel="Probability"
  #if you have a factor, do a different type of plotting
  if(grepl("factor",type)){
    pr = predict(mod,type = "response")
    prdf = data.frame(pr,Variable)
    boxplot(pr~Variable,data = prdf,xlab = xlabel,ylab=ylabel,cex.lab = 1.5)
    #plot outliers?
  }else{
    #Probability of covariate
    BootstrapParameters3<-rmvnorm(10000, coef(mod),summary(mod)$cov.unscaled)
    PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
    CenterVar3<-model.matrix(mod)[,start:finish]*coef(mod)[c(start:finish)]
    BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
    # if (grepl("factor",type)){
    # Basis3<-gam(rbinom(5000,1,0.5)~as.factor(PlottingVar3),
    #             fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=4)))$X[,-1]
    if (grepl("linear",type)){
      Basis3<-gam(rbinom(5000,1,0.5)~PlottingVar3,
                  fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=4)))$X[,-1]
    }else if (grepl("smooth",type)){
      Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3,k=4),
                  fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=4)))$X[,-1]
    }else if (grepl("cyclical",type)){
      Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4),
                  fit=F, family=binomial, knots=list(PlottingVar2=seq(range[1],range[2],length=4)))$X[,-1]
    }
    if (start!=finish){ #i.e. if theres multiple terms
      # #see if there's a mismatch of terms, in which case only use 3 of the 4 coef values. idk why this is coming up
      # #and this maybe isn't a very good solution. But...
      if (size(Basis3,2)!= size(coef(mod)[c(start:finish)],2)){
        shortcoef = coef(mod)[c(start:finish)]
        rmvcoef = size(coef(mod)[c(start:finish)],2) - size(Basis3,2)
        shortcoef = shortcoef[-c(1:rmvcoef)]
        RealFit3 = Basis3%*%shortcoef
      }else{
        RealFit3<-Basis3%*%coef(mod)[c(start:finish)]}
    }else{
      #repmat the coefficient to needed length
      coefuse = coef(mod)[c(start:finish)]
      #repeat it correctly, so using smaller dimension of Basis3
      coefmat = repmat(coefuse,1,min(size(Basis3)))
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
      if (size(Basis3,2)!=size(BootstrapCoefs3,2)){
        bootstrapshort = BootstrapCoefs3[,-c(1:rmvcoef)]
        BootstrapFits3<-Basis3%*%t(bootstrapshort)
      }else{
        BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)}
    }
    quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
    cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
    cis3a<-inv.logit(cis3)
    plot(PlottingVar3,RealFitCenter3a, col = 1,type="l",ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=range, main =NULL , cex.lab = 1.5, cex.axis=1.5)
    segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
    lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
    #add violin plot to plots
    ###EDIT TO CORRECTLY MAKE VIOLIN USING VALUES, AND PLOT WITH MY STUFF
    violinrug<-myviolin(xvals = Variable,miny = 0,maxy = max(RealFitCenter3a),rel=5)
    polygon(violinrug,col=rgb(0,1,1,alpha=.3))
  }}
###plot hour of day
par(mfrow=c(2,2))
# plotbase = str_replace(plotname2,'smooths.png','')
# plotname3 = paste(plotbase,paste(termlist[1],'_partial.png',sep=""),sep="")
#need to use coeff() to get starts and finishes for all variables
coeforder = names(coef(m1a))
termlist = c('hrofday','moonfrac','julday','year')
#since variable was a cyclic mat, pass in the ORIGINAL variable, pre-cycling
hrind = str_which(coeforder,"hr")
plot_covariate(m1a,Daily$hrofday,termlist[1],range = c(0,23),type="cyclical",st=min(hrind),fn=max(hrind),data=Daily)
title(main = namesh)
###plot moonuse
#figure out which one is moon
moonind = str_which(coeforder,"moon")
plot_covariate(m1a,Daily$moonfrac,termlist[2],range=c(0,100),type=moontype,st=min(moonind),fn=max(moonind),data=Daily)
###plot jul
julind = str_which(coeforder,"ju")
plot_covariate(m1a,Daily$juldays,termlist[3],range=c(0,365),type="cyclical",st=min(julind),fn=max(julind),data=Daily)
###plot year
#figure out which one is year
yearind = str_which(coeforder,"year")
plot_covariate(m1a,Daily$year,termlist[4],range=c(min(Daily$year),max(Daily$year)),type="factor",st=min(yearind),fn=max(yearind),data=Daily)
#saveplot
plotname2 = str_replace(plotnameres,'_res.png','_smooths.png')
dev.copy(png,plotname2)
dev.off()