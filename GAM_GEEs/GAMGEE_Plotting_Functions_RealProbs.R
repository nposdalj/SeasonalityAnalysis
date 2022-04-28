################General Patterns########################
# Plot Julian Day with a Base R Plot ---------------------------------------------------------
BasePlot_JD <- function(model, table){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; finish=3; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)]
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter4<- RealFitCenter3 + coef(model)[1]
  RealFitCenter3a<- inv.logit(RealFitCenter4)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  BootstrapFit3a<- BootstrapFits3 +coef(model)[1]
  cis3<-apply(BootstrapFit3a, 1, quant.func3)
  cis3a<-inv.logit(cis3)
  
  #Base R Plotting
  title = paste(saveDir,"/BaseR_Julian Day - ", site,".png",sep="")
  png(title)
  xlabel="Julian Day"; ylabel="Probability" 
  plot(PlottingVar3,(RealFitCenter3a), type='l', col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,365), main = title , cex.lab = 1.5, cex.axis=1.5)    
  segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
  lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
  rug(PlottingVar3)
  dev.off()
  print('Plot Saved')
}

# Plot Julian Day with a Base R Plot (CB/Big) ---------------------------------------------------------
BasePlot_JD_Year <- function(model, table){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled) #probability of covariate
  start=9; finish=10; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000) #use version of variable with range of original, linearly spaced
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)] #center the output so nothing is below zero
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)] #coefficients of 10000 bootstrapped models to give us gray area when multiplied with our splines; what would 10000 of these look like
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)]
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter4<- RealFitCenter3 + coef(model)[1]
  RealFitCenter3a<- inv.logit(RealFitCenter4)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  BootstrapFit3a<- BootstrapFits3 +coef(model)[1]
  cis3<-apply(BootstrapFit3a, 1, quant.func3)
  cis3a<-inv.logit(cis3)
  
  #Base R Plotting
  title = paste(saveDir,"/BaseR_Julian Day - ", site,".png",sep="")
  png(title)
  xlabel="Julian Day"; ylabel="Probability" 
  plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,365), main = title , cex.lab = 1.5, cex.axis=1.5)    
  segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
  lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
  rug(PlottingVar3)
  dev.off()
  print('Plot Saved')
}


# Plot Julian Day with ggplot ---------------------------------------------
ggPlot_JD <- function(model, table, site){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; finish=3; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter3a<-inv.logit(RealFitCenter3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  cis3a<-inv.logit(cis3)
  
  # Calculate kernel density of Jday observations
  dJday = stats::density(Variable,na.rm = TRUE,n=5000,from=1,to=365)
  dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
  colnames(dens) = c("Day", "Density")
  dens$Density = dens$Density / 0.15 #max(dens$Density) # normalize kernel density
  if (min(cis3[1,])<0){ # set kernel density at bottom of y axis
    dens$Density = dens$Density - abs(min(cis3[1,])) 
  } else {
    dens$Density = dens$Density + min(cis3[1,])
  }
  
  plotDF = data.frame(PlottingVar3, RealFitCenter3a)
  colnames(plotDF) = c("Jday", "Fit")
  
    #Histogram for Jday observations
    x = seq(from = 0,to = 365,by = 1)
    fullhist = hist(SiteHourTableB$Julian,x)
    yhist = fullhist$counts
    #normalize it to be @ most half as high as y scale
    ymod = yhist/(max(yhist)-abs(min(cis3a[1,])))
    ymodd = yhist/(max(yhist)-min(cis3[1,]))
    newy = ymod/((max(ymod)/min(plotDF$Fit))/.5) #used to use 3.5 but it was too generic
    plothist = data.frame(x=x[-1],y=newy)
    
    #Figure out y-scale labels based on plotting julian day as a coefficient
    plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
    colnames(plotDF_labels) = c("Jday", "Fit")
    breaksy = seq(min(cis3a[1,]),max(cis3a[2,]),length.out = 3)
    diffy = breaksy[3]-breaksy[1]
    breaksP = seq(min(cis3[1,]),max(cis3[2,]),length.out = 3)
    diffP = breaksP[3]+abs(breaksP[1])
    breaksH = seq(min(newy),max(cis3a[2,]),length.out=3)
    diffH = breaksH[3]-breaksH[1]
    diffNP = (diffH*diffP)/diffy
    BottomL = breaksP[3]-diffNP
    breaksPH = rev(round(seq(breaksP[3],BottomL,length.out=3),digits = 1))
    breaksPH_labels = as.character(breaksPH)
  
  ggplot(plotDF, aes(Jday, Fit),
  ) + geom_smooth(fill = "grey",
                  colour = "black",
                  aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                  stat ="identity"
  ) + geom_bar(data = plothist,
            aes(x,y),
            fill = 4,
            alpha = 0.2,
            position = "dodge",
            stat = "identity",
            width = 1
  ) + coord_cartesian(ylim = c(min(newy), max(cis3a[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)",
           title = paste('Julian Day at',site),
  ) + scale_y_continuous(breaks = seq(min(newy),max(cis3a[2,]),length.out = 3), labels = breaksPH_labels
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) +theme(axis.line = element_line(),
            panel.background = element_blank(),
           text = element_text(size = 25)
  )
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,".pdf",sep="")
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}
# Plot Julian Day with ggplot (Big) ---------------------------------------------
ggPlot_JD_Year <- function(model, table){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=9; finish=10; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter3a<-inv.logit(RealFitCenter3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  cis3a<-inv.logit(cis3)
  
  # Calculate kernel density of Jday observations
  dJday = stats::density(Variable,na.rm = TRUE,n=5000,from=1,to=365)
  dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
  colnames(dens) = c("Day", "Density")
  dens$Density = dens$Density / 0.15 #max(dens$Density) # normalize kernel density
  if (min(cis3[1,])<0){ # set kernel density at bottom of y axis
    dens$Density = dens$Density - abs(min(cis3[1,])) 
  } else {
    dens$Density = dens$Density + min(cis3[1,])
  }
  
  plotDF = data.frame(PlottingVar3, RealFitCenter3a)
  colnames(plotDF) = c("Jday", "Fit")
  
  #Histogram for Jday observations
  x = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,x)
  yhist = fullhist$counts
  #normalize it to be @ most half as high as y scale
  ymod = yhist/(max(yhist)-abs(min(cis3a[1,])))
  ymodd = yhist/(max(yhist)-min(cis3[1,]))
  newy = ymod/((max(ymod)/min(plotDF$Fit))/.5) #used to use 3.5 but it was too generic
  plothist = data.frame(x=x[-1],y=newy)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  breaksy = seq(min(cis3a[1,]),max(cis3a[2,]),length.out = 3)
  diffy = breaksy[3]-breaksy[1]
  breaksP = seq(min(cis3[1,]),max(cis3[2,]),length.out = 3)
  diffP = breaksP[3]+abs(breaksP[1])
  breaksH = seq(min(newy),max(cis3a[2,]),length.out=3)
  diffH = breaksH[3]-breaksH[1]
  diffNP = (diffH*diffP)/diffy
  BottomL = breaksP[3]-diffNP
  breaksPH = rev(round(seq(breaksP[3],BottomL,length.out=3),digits = 1))
  breaksPH_labels = as.character(breaksPH)
  
  ggplot(plotDF, aes(Jday, Fit),
  ) + geom_smooth(fill = "grey",
                  colour = "black",
                  aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                  stat ="identity"
  ) + geom_bar(data = plothist,
               aes(x,y),
               fill = 4,
               alpha = 0.2,
               position = "dodge",
               stat = "identity",
               width = 1
  ) + coord_cartesian(ylim = c(min(newy), max(cis3a[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)",
           title = paste('Julian Day at',site),
  ) + scale_y_continuous(breaks = seq(min(newy),max(cis3a[2,]),length.out = 3), labels = breaksPH_labels
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) +theme(axis.line = element_line(),
           panel.background = element_blank(),
           text = element_text(size = 25)
  )
  
  ggplot(plotDF, aes(Jday, Fit),
  ) + geom_smooth(fill = "grey",
                  colour = "black",
                  aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                  stat ="identity"
  ) + gg
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,".pdf",sep="")
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}
# Plot Julian Day with a Base R Plot (GOA) ---------------------------------------------------------
BasePlot_JD_AfterSite <- function(model, table){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=9; finish=10; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter3a<-inv.logit(RealFitCenter3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  cis3a<-inv.logit(cis3)
  
  #Base R Plotting
  title = paste(saveDir,"/BaseR_Julian Day - ", site,".png",sep="")
  png(title)
  xlabel="Julian Day"; ylabel="Probability" 
  plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,365), main = title , cex.lab = 1.5, cex.axis=1.5)    
  segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
  lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
  rug(PlottingVar3)
  dev.off()
  print('Plot Saved')
}

# Plot Julian Day with ggplot (GOA) ---------------------------------------------
ggPlot_JD_AfterSite <- function(model, table){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=6; finish=7; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter3a<-inv.logit(RealFitCenter3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  cis3a<-inv.logit(cis3)
  
  plotDF = data.frame(PlottingVar3, RealFitCenter3a)
  colnames(plotDF) = c("Jday", "Fit")
  
  # Calculate kernel density of Jday observations
  dJday = stats::density(Variable,na.rm = TRUE,n=5000,from=1,to=365)
  dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
  colnames(dens) = c("Day", "Density")
  dens$Density = dens$Density / 0.15 #max(dens$Density) # normalize kernel density
  if (min(cis3a[1,])<0){ # set kernel density at bottom of y axis
    dens$Density = dens$Density - abs(min(cis3a[1,])) 
  } else {
    dens$Density = dens$Density + min(cis3a[1,])
  }
  
  #Histogram for Jday observations
  x = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,x)
  yhist = fullhist$counts
  #normalize it to be @ most half as high as y scale
  ymod = yhist/(max(yhist)-abs(min(cis3a[1,])))
  ymodd = yhist/(max(yhist)-min(cis3[1,]))
  newy = ymod/((max(ymod)/min(plotDF$Fit))/.5) #used to use 3.5 but it was too generic
  plothist = data.frame(x=x[-1],y=newy)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  breaksy = seq(min(cis3a[1,]),max(cis3a[2,]),length.out = 3)
  diffy = breaksy[3]-breaksy[1]
  breaksP = seq(min(cis3[1,]),max(cis3[2,]),length.out = 3)
  diffP = breaksP[3]+abs(breaksP[1])
  breaksH = seq(min(newy),max(cis3a[2,]),length.out=3)
  diffH = breaksH[3]-breaksH[1]
  diffNP = (diffH*diffP)/diffy
  BottomL = breaksP[3]-diffNP
  breaksPH = rev(round(seq(breaksP[3],BottomL,length.out=3),digits = 1))
  breaksPH_labels = as.character(breaksPH)
  
  ggplot(plotDF, aes(Jday, Fit),
  ) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                stat ="identity"
    ) + geom_bar(data = plothist,
                 aes(x,y),
                 fill = 4,
                 alpha = 0.2,
                 stat = "identity",
                 width=1
    ) + coord_cartesian(ylim = c(min(newy), abs(max(cis3a[2,])))
    ) + labs(x = "Julian Day",
             y = "Probability",
             title = paste('Julian Day at',site),
    ) + scale_y_continuous(breaks = seq(min(newy),max(cis3a[2,]),length.out = 3), labels = breaksPH_labels
    ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
    ) +theme(axis.line = element_line(),
             panel.background = element_blank(),
             text = element_text(size = 25)
    )
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,".pdf",sep="")
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Plot Year as factor with ggplot --------------------------------------------------
ggPlot_Year <- function(model, table,site){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=8; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1],
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5],
    BootstrapCoefs2[, 6],
    BootstrapCoefs2[, 7],
    BootstrapCoefs2[, 8]),
  as.factor(rep(1:8, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2011","2012","2013","2014","2015","2017","2018","2019")
  names(trans) = c(1,2,3,4,5,6,7,8)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,".pdf",sep="")
  ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(
  ) + theme(axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + labs(title = paste('Year at',site), 
           x = "Year",
           y = "s(Year)")
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}

#Plot Year as factor with ggplot for Big Model --------------------------------------------------
ggPlot_Year_Big <- function(model, table,site){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=9; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1],
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5],
    BootstrapCoefs2[, 6],
    BootstrapCoefs2[, 7],
    BootstrapCoefs2[, 8],
    BootstrapCoefs2[, 9]),
    as.factor(rep(1:9, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2010","2011","2012","2013","2014","2015","2017","2018","2019")
  names(trans) = c(1,2,3,4,5,6,7,8,9)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,".pdf",sep="")
  ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(
  ) + theme(axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + labs(title = paste('Year at',site), 
           x = "Year",
           y = "s(Year)")
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}

#Plot Site with ggplot --------------------------------------------------
ggPlot_Site <- function(model, table,site){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  val=4; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,val)]
  
  # Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(c(BootstrapCoefs2[, 1],
                                   BootstrapCoefs2[,2]),
                                 as.factor(rep(1:2, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("KS","BD")
  names(trans) = c(1,2)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]

  ggtitle = paste(saveDir,"/Site - ", region,".pdf",sep="")
  ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(
  ) + theme(axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + labs(title = paste('Site at',region), 
           x = "Site",
           y = "s(Site)")
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}

#Plot Site with ggplot as a factor --------------------------------------------------
ggPlot_Site_asFactor <- function(model,table,region){
  
  ggtitle = paste(saveDir,"/Probability of Site - ", region,".pdf",sep="")
  table$pr = pr
  ggplot(table, aes(x = Site, y = pr)) +
    geom_boxplot(aes(fill = factor(Site)), alpha = .2)
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}
#Plot Site with ggplot (GOA) --------------------------------------------------
ggPlot_Site_GOA <- function(model, table,region){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=5; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]

  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1],
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]
  ),
                                 as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("AB","CB","KOA","PT","QN")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]

  ggtitle = paste(saveDir,"/Site - ", region,".pdf",sep="")
  ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(
  ) + theme(axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + labs(title = paste('Site at',region), 
           x = "Site",
           y = "s(Site)")

  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}

#Plot Site with ggplot as a factor (GOA) --------------------------------------------------
ggPlot_Site_asfactor_GOA <- function(model,table,region){

  ggtitle = paste(saveDir,"/Probability of Site - ", region,".pdf",sep="")
  table$pr = pr
  ggplot(table, aes(x = Site, y = pr)) +
    geom_boxplot(aes(fill = factor(Site)), alpha = .2)

  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}
#Plot Region with ggplot --------------------------------------------------
ggPlot_Region_Big <- function(model, table,region){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  val = 12; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,val)]
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1],
    BootstrapCoefs2[, 2]
  ),
  as.factor(rep(1:2, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("BSAI","GOA")
  names(trans) = c(1,2)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]
  
  ggtitle = paste(saveDir,"/Site - ", region,".pdf",sep="")
  ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(
  ) + theme(axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + labs(title = paste('Site at',region), 
           x = "Site",
           y = "s(Site)")
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}

#Plot Region with ggplot as a factor --------------------------------------------------
ggPlot_Region_asFactor_Big <- function(model,table,region){
  
  ggtitle = paste(saveDir,"/Probability of Region - ", region,".pdf",sep="")
  table$pr = pr
  ggplot(table, aes(x = Region, y = pr)) +
    geom_boxplot(aes(fill = factor(Region)), alpha = .2)
  
  ggsave(
    ggtitle,
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
}
