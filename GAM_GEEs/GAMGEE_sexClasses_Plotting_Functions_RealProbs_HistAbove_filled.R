################Sex Specific###########################
# Plot Julian Day with ggplot ---------------------------------------------
ggPlot_JD_sex <- function(model, table,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; finish=3; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

ggPlot_JD_AfterYear_WAT <- function(model, table,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=5; finish=6; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

ggPlot_JD_SG_WATBIG <- function(model, table,site,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=7; finish=8; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

ggPlot_JD_Last <- function(model,table,site,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=8; finish=9; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

ggPlot_JD_WATBIG <- function(model,table,site,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=3; finish=4; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}



# Plot Julian Day with ggplot (When Julian day is after year for GOA (Mid-Size) or CB)  ---------------------------------------------
ggPlot_JD_AfterYear <- function(model,table,site,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=9; finish=10; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}
# Plot Julian Day with ggplot (when julian day is after site for GOA (Males)) ---------------------------------------------
ggPlot_JD_AfterSite <- function(model,table,site,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=6; finish=7; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}
# Plot Julian Day with ggplot (when julian day is after year and site for GOA (Social Groups)) ---------------------------------------------
ggPlot_JD_AfterYearSite <- function(model,table,site,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=13; finish=14; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

# Plot Julian Day with ggplot (when julian day is after Year for Big Model) ---------------------------------------------
ggPlot_JD_AfterYearB <- function(model,table,site,sex,COL){
  BootstrapParameters3<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=10; finish=11; Variable=table$Julian;  
  PlottingVar3<-seq(min(Variable), max(Variable), length=5000)
  CenterVar3<-model.matrix(model)[,start:finish]*coef(model)[c(start:finish)]
  BootstrapCoefs3<-BootstrapParameters3[,c(start:finish)]
  Basis3<-gam(rbinom(5000,1,0.5)~s(PlottingVar3, bs="cc", k=4), fit=F, family=binomial, knots=list(PlottingVar2=seq(1,365,length=4)))$X[,2:3]
  RealFit3<-Basis3%*%coef(model)[c(start:finish)] #COEFFICIENTS TO PLOT
  RealFitCenter3<-RealFit3-mean(CenterVar3)
  BootstrapFits3<-Basis3%*%t(BootstrapCoefs3)
  quant.func3<-function(x){quantile(x,probs=c(0.025, 0.975))}
  cis3<-apply(BootstrapFits3, 1, quant.func3)-mean(CenterVar3)
  
  #Histogram for Jday observations
  Jday = seq(from = 0,to = 365,by = 1)
  fullhist = hist(SiteHourTableB$Julian,Jday)
  yhist = fullhist$counts
  plothist = data.frame(x=Jday[-1],y=yhist)
  
  #Figure out y-scale labels based on plotting julian day as a coefficient
  plotDF_labels = data.frame(PlottingVar3, RealFitCenter3)
  colnames(plotDF_labels) = c("Jday", "Fit")
  
  pmain = ggplot(plotDF_labels, aes(Jday, Fit),
  ) + geom_smooth(fill = COL, colour = "black", aes(ymin=cis3[1,], ymax=cis3[2,]), stat ="identity",size = 1
  ) + coord_cartesian(ylim = c(min(cis3[1,]), max(cis3[2,]))
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)"
  ) + scale_x_continuous(breaks = seq(20,350,length.out = 12),labels = c('J','F','M','A','M','J','J','A','S','O','N','D')
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  )
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            position = "dodge",
  #            stat = "identity",
  #            width = 1
  #   )
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggtitle = paste(saveDir,"/Julian Day - ", site,'_',sex,"filled.pdf",sep="")
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

#Plot Year as factor with ggplot --------------------------------------------------
ggPlot_Year <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=8; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' '
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
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
  
  ggtitle = paste(saveDir,"/Year - ", site,'_',sex,"filled.pdf",sep="")
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2011","2012","2013","2014","2015","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2011","2012","2013","2014","2015","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

#Plot Year as factor with ggplot (males at CB) --------------------------------------------------
ggPlot_Year_AfterJD <- function(model,table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=4; end=10; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' '
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
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
  
  ggtitle = paste(saveDir,"/Year - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2011","2012","2013","2014","2015","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2011","2012","2013","2014","2015","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Plot Year as factor with ggplot for Big Model --------------------------------------------------
ggPlot_Year_Big <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=9; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' '
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
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
  
  ggtitle = paste(saveDir,"/Year - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2010","2011","2012","2013","2014","2015","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2010","2011","2012","2013","2014","2015","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}


#Plot Site with ggplot --------------------------------------------------
ggPlot_Site <- function(model,table,sex,COL){
BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
val=4; Variable=table$Site
BootstrapCoefs2<-BootstrapParameters2[, c(1,val)]

#Histogram for Site observations
plothist = data.frame(table(SiteHourTableB$Site))
colnames(plothist) = c("x","y")

# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedSiteCoefs = data.frame(c(BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
                                 BootstrapCoefs2[,2]),
                               as.factor(rep(1:2, each = 10000)))
colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
trans = c("KS","BD")
names(trans) = c(1,2)
AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]

ggtitle = paste(saveDir,"/Site - ", site,'_',sex,"filled.pdf",sep="")

pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
) + geom_boxplot(fill=COL
) + theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.line = element_line(),
          panel.background = element_blank(),
          text = element_text(size = 25)
) + scale_x_discrete(labels = c("BD","KS")
) + labs(x = "Site",
         y = "s(Site)")

# xens = axis_canvas(pmain, axis = "x")+
#   geom_bar(data = plothist,
#            aes(x,y),
#            fill = 4,
#            alpha = 0.2,
#            #position = "dodge",
#            stat = "identity",
#            width = 1) + scale_x_discrete(labels = c("BD","KS"))
# 
# p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
ggdraw(pmain)

ggsave(
  ggtitle, width = 7, height = 6, units = "in",
  device = "pdf") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device
print("Plot Saved")

}

#Plot Site with ggplot (GOA - Males) --------------------------------------------------
ggPlot_Site_Year <- function(model,table,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=5; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Site observations
  plothist = data.frame(table(SiteHourTableB$Site))
  colnames(plothist) = c("x","y")
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
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
  
  ggtitle = paste(saveDir,"/Site - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("AB","CB","KOA","PT","QN")
  ) + labs(x = "Site",
           y = "s(Site)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("AB","CB","KOA","PT","QN"))
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

#Plot Site for WAT Regional model (first term in south model) --------------------------------------------------
ggPlot_Site_WAT <- function(model,table,site,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=4; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Site observations
  plothist = data.frame(table(SiteHourTableB$Site))
  colnames(plothist) = c("x","y")
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4]
  ),
  as.factor(rep(1:4, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("GS","BP","BS","JAX")
  names(trans) = c(1,2,3,4)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]
  
  ggtitle = paste(saveDir,"/Site - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("GS","BP","BS","JAX")
  ) + labs(x = "Site",
           y = "s(Site)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("GS","BS","BP","JAX"))
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Plot region for Wat --------------------------------------------------
ggPlot_Region_WAT <- function(model,table,site,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=2; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Site observations
  plothist = data.frame(table(SiteHourTableB$Site))
  colnames(plothist) = c("x","y")
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2]
  ),
  as.factor(rep(1:2, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("North","South")
  names(trans) = c(1,2)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]
  
  ggtitle = paste(saveDir,"/Region - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("North","South")
  ) + labs(x = "Region",
           y = "s(Region)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("North","South"))
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Plot region for Wat (last term) --------------------------------------------------
ggPlot_Region_WAT_last <- function(model,table,site,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=8; end=8; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Site observations
  plothist = data.frame(table(SiteHourTableB$Site))
  colnames(plothist) = c("x","y")
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2]
  ),
  as.factor(rep(1:2, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("North","South")
  names(trans) = c(1,2)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]
  
  ggtitle = paste(saveDir,"/Region - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("North","South")
  ) + labs(x = "Region",
           y = "s(Region)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("North","South"))
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Plot Site for WAT Regional model (first term in northern model) --------------------------------------------------
ggPlot_Site_WATT <- function(model,table,site,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=5; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Site observations
  plothist = data.frame(table(SiteHourTableB$Site))
  colnames(plothist) = c("x","y")
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]
  ),
  as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("HZ","OC","NC","BC","WC")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]
  
  ggtitle = paste(saveDir,"/Site - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("HZ","OC","NC","BC","WC")
  ) + labs(x = "Site",
           y = "s(Site)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("HZ","OC","NC","BC","WC"))
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Plot Site for WAT Regional model (second term in northern model) --------------------------------------------------
ggPlot_Site_WATTT <- function(model,table,site,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=4; end=7; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Site observations
  plothist = data.frame(table(SiteHourTableB$Site))
  colnames(plothist) = c("x","y")
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]
  ),
  as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Coefficient", "Site")
  trans = c("HZ","OC","NC","BC","WC")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$SiteName = trans[as.character(AdjustedSiteCoefs$Site)]
  
  ggtitle = paste(saveDir,"/Site - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("HZ","OC","NC","BC","WC")
  ) + labs(x = "Site",
           y = "s(Site)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("HZ","OC","NC","BC","WC"))
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Plot Site with ggplot (GOA - Social Groups) --------------------------------------------------
ggPlot_Site_AfterYear <- function(model,table,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=9; end=12; Variable=table$Site
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Site observations
  plothist = data.frame(table(SiteHourTableB$Site))
  colnames(plothist) = c("x","y")
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
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
  
  ggtitle = paste(saveDir,"/Site - ", site,'_',sex,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(SiteName, Coefficient)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("AB","CB","KOA","PT","QN")
  ) + labs(x = "Site",
           y = "s(Site)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = plothist,
  #            aes(x,y),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("AB","CB","KOA","PT","QN"))
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
  
}

#For WAT
#Year is second term
ggPlot_Year_WAT <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=4; end=6; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' '
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4]),
    as.factor(rep(1:4, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2016","2017","2018","2019")
  names(trans) = c(1,2,3,4)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Year is second term for regional model
ggPlot_Year_WAT_regional <- function(model, table,site,sex,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=5; end=7; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4]),
    as.factor(rep(1:4, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2016","2017","2018","2019")
  names(trans) = c(1,2,3,4)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ",sex,'_', site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

ggPlot_Year_WATT <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=4; end=7; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(table$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]),
    as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2015", "2016","2017","2018","2019")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2015", "2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

ggPlot_Year_WATTT_north <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=8; end=11; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(table$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]),
    as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2015", "2016","2017","2018","2019")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2015", "2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

ggPlot_Year_WATT_north <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=6; end=9; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(table$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]),
    as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2015", "2016","2017","2018","2019")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2015", "2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Year is first term
ggPlot_Year_WAT_first <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=4; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4]),
    as.factor(rep(1:4, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2016","2017","2018","2019")
  names(trans) = c(1,2,3,4)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Year is first term for Big Model
ggPlot_Year_WAT_Big <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=2; end=5; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]),
    as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2015","2016","2017","2018","2019")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Year is first term for Big Model
ggPlot_Year_WATT_Big <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=3; end=6; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 
  
  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]),
    as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2015","2016","2017","2018","2019")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}

#Year is first term for Big Model
ggPlot_Year_WATTT_Big <- function(model, table,site,COL){
  BootstrapParameters2<-rmvnorm(10000, coef(model),summary(model)$cov.unscaled)
  start=5; end=8; Variable=table$Year
  BootstrapCoefs2<-BootstrapParameters2[, c(1,start:end)]
  
  #Histogram for Year observations
  counts = count(SiteHourTableB$Year)
  counts$x = as.character(counts$x)
  counts$freq = as.numeric(counts$freq)
  counts$days = round(counts$freq/12)
  counts$label = '*'
  counts$label[counts$days > 365] <- ' ' 

  #Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedSiteCoefs = data.frame(  c(
    BootstrapCoefs2[, 1] - mean(BootstrapCoefs2[, 1]),
    BootstrapCoefs2[, 2],
    BootstrapCoefs2[, 3],
    BootstrapCoefs2[, 4],
    BootstrapCoefs2[, 5]),
    as.factor(rep(1:5, each = 10000)))
  colnames(AdjustedSiteCoefs) = c("Probability", "Year")
  trans = c("2015","2016","2017","2018","2019")
  names(trans) = c(1,2,3,4,5)
  AdjustedSiteCoefs$YearVal = trans[as.character(AdjustedSiteCoefs$Year)]
  
  ggtitle = paste(saveDir,"/Year - ", site,"filled.pdf",sep="")
  
  pmain = ggplot(AdjustedSiteCoefs, aes(YearVal, Probability)
  ) + geom_boxplot(fill=COL
  ) + theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),

            axis.line = element_line(),
            panel.background = element_blank(),
            text = element_text(size = 25)
  ) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  # xens = axis_canvas(pmain, axis = "x")+
  #   geom_bar(data = counts,
  #            aes(x,freq),
  #            fill = 4,
  #            alpha = 0.2,
  #            #position = "dodge",
  #            stat = "identity",
  #            width = 1) + scale_x_discrete(labels = c("2015","2016","2017","2018","2019")
  #            )
  # 
  # p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(pmain)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")
}