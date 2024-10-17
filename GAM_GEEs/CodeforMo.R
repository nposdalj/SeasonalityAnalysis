rm(list = ls()) #clear environment

#load libraries
library(boot)
library(pracma)
library(geepack)         # for the GEEs (Wald's hypothesis tests allowed)
library(splines)         # to construct the B-splines within a GEE-GLM
library(tidyverse)       # because it literally does everything
library(ROCR)            # to build the ROC curve
library(PresenceAbsence) # to build the confusion matrix
library(ggplot2)         # to build the partial residual plots
library(mvtnorm)         # to build the partial residual plots
library(gridExtra)       # to build the partial residual plots
library(lubridate)
library(regclass)
library(mgcv)
library(ChemoSpecUtils)
library(car)            # to run an ANOVA
library(splines2)       # to use mSpline for the GEEs
library(scales)
library(magick)
library(cowplot)
library(ggExtra)
library(plyr)

SiteHourTableB = SiteHourTable[!(SiteHourTableB$Julian==366),] #If it's a leap year, delete julian day 366 for plotting
SiteHourTableB$Year = as.factor(SiteHourTableB$Year) #Change year to categorical variable for plotting histograms
COL = '#D3D3D3'
model = PODFinal
table = SiteHourTableB

  ggPlot_Year_WAT(PODFinal,SiteHourTableB,site,COL)


## I had the following part as a function that I then ran with this:
  #ggPlot_Year_WAT(PODFinal,SiteHourTableB,site,COL)
  
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
            text = element_text(size = 35)
  ) + scale_x_discrete(labels = c("2016","2017","2018","2019")
  ) + labs(x = "Year",
           y = "s(Year)")
  
  xens = axis_canvas(pmain, axis = "x")+
    geom_bar(data = counts,
             aes(x,freq),
             fill = 4,
             alpha = 0.2,
             #position = "dodge",
             stat = "identity",
             width = 1) + scale_x_discrete(labels = c("2016","2017","2018","2019")
             )
  
  p1 = insert_xaxis_grob(pmain,xens,grid::unit(.2, "null"), position = "top")
  ggdraw(p1)
  
  ggsave(
    ggtitle, width = 7, height = 6, units = "in",
    device = "pdf") # save figure
  while (dev.cur() > 1) {
    dev.off()
  } # close graphics device
  print("Plot Saved")



