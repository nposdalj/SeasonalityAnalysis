# Plot Julian Day with a Base R Plot ---------------------------------------------------------
BasePlot_JD <- function(model, table){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; finish=3; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,366,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)]
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
  plot(PlottingVar3,(RealFitCenter3a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(1,366), main = title , cex.lab = 1.5, cex.axis=1.5)    
  segments(PlottingVar3,(cis3a[1,]),PlottingVar3,(cis3a[2,]), col="grey")
  lines(PlottingVar3,(RealFitCenter3a),lwd=2, col=1)
  rug(PlottingVar3)
  dev.off()
  print('Plot Saved')
}


# Plot Julian Day with ggplot ---------------------------------------------
ggPlot_JD <- function(model, table){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; finish=3; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,366,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)]
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  RealFitCenter3a<-inv.logit(RealFitCenter3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  cis3a<-inv.logit(cis3)
  
  # Calculate kernel density of Jday observations
  dJday = stats::density(Variable,na.rm = TRUE,n=5000,from=1,to=366)
  dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
  colnames(dens) = c("Day", "Density")
  dens$Density = dens$Density / 0.15 #max(dens$Density) # normalize kernel density
  if (min(cis3a[1,])<0){ # set kernel density at bottom of y axis
    dens$Density = dens$Density - abs(min(cis3a[1,])) 
  } else {
    dens$Density = dens$Density + min(cis3a[1,])
  }
  
  plotDF = data.frame(PlottingVar3, RealFitCenter3a)
  colnames(plotDF) = c("Jday", "Fit")
  
  ggplot(plotDF, aes(Jday, Fit),
  ) + geom_polygon(data=dens,
                   aes(Day,Density),
                   fill=4,
                   alpha=0.2
  ) + geom_smooth(fill = "grey",
                  colour = "black",
                  aes(ymin=cis3a[1,], ymax=cis3a[2,]),
                  stat ="identity"
  ) + labs(x = "Julian Day",
           y = "Probability",
           title = paste('Julian Day at',site),
  ) + theme(axis.line = element_line(),
            panel.background = element_blank()
  )
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,".png",sep="")
  
  ggsave(
    ggtitle,
    device = "png") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

# Plot Year ---------------------------------------------------------------


