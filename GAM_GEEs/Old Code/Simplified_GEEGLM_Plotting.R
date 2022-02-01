library(tidyverse)
library(geepack)
library(mgcv)
library(pracma)
library(splines2)
library(SimDesign)


# CONSTRUCT GEEGLM
# ThisSite is daily counts of 5-min bins with presence, Jday is Julian day, and yearGroup is year
modJday = geeglm(
  thisSite ~ mSpline(
    Jday,
    knots = quantile(Jday, probs = c(0.333, 0.666)),
    Boundary.knots = c(1, 365),
    periodic = T
  )
  + as.factor(yearGroup),
  family = poisson(link = "log"),
  id = reducedClustID,
  # grouping variable to account for autocorrelation
  corstr = "ar1"
)


# BOOTSTRAP GEEGLM PARAMETER ESTIMATES based on estimated means and variances,
# for later construction of confidence intervals
BootstrapParameters <-
  rmvnorm(10000, coef(modJday), summary(modJday)$cov.unscaled)
JDayBootstrapCoefs <- BootstrapParameters[, 2:3]
YearBootstrapCoefs <- BootstrapParameters[, c(1, 4:6)]

# Predict presence at each X value based on one covariate (I think; not sure why we need this)
Jx1 <- model.matrix(modJday)[, 2:3] %*% coef(modJday)[c(2:3)]
# Jx2 <- model.matrix(testJday)[, c(1, 4:6)] %*% coef(modJday)[c(1, 4:6)]


# PLOT GEEGLM PARTIAL RESIDUALS
# Julian Day
JDayForPlotting <- seq(min(Jday), max(Jday), length = 256)
# get basis functions for smooth of Julian day
Basis <-
  mSpline(
    JDayForPlotting,
    knots = c(120, 250),
    Boundary.knots = c(1, 365),
    periodic = T
  )
# multiply basis functions by model coefficients to get values of spline at each X
RealFit <- Basis %*% coef(modJday)[c(2:3)]
# adjust offset
RealFitCenterJ <- RealFit - mean(Jx1) - coef(modJday)[1] # again, not sure why the Jx1 adjustment is needed
# create data frame for plotting
plotDF = data.frame(JDayForPlotting, RealFitCenterJ)
colnames(plotDF) = c("Jday", "Fit")

# get spread of spline values based on distributions of each coefficient
JDayBootstrapFits <- Basis %*% t(JDayBootstrapCoefs)
# get quantiles for confidence interval of smooth function estimate
quant.func <- function(x) {
  quantile(x, probs = c(0.025, 0.975))
}
cisJ <- apply(JDayBootstrapFits, 1, quant.func)
Jcil <- cisJ[1, ] - mean(Jx1) - coef(modJday)[1] # upper CI bound
Jciu <- cisJ[2, ] - mean(Jx1) - coef(modJday)[1] # lower CI bound

# Calculate kernel density of Jday observations
dJday = stats::density(Jday,na.rm = TRUE,n=256,from=1,to=365)
dens = data.frame(c(dJday$x, rev(dJday$x)), c(dJday$y, rep(0, length(dJday$y))))
colnames(dens) = c("Day", "Density")
dens$Density = dens$Density / max(dens$Density) # normalize kernel density
if (min(Jcil)<0){ # set kernel density at bottom of y axis
  dens$Density = dens$Density - abs(min(Jcil)) 
} else {
  dens$Density = dens$Density + min(Jcil)
}

saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDayPlot.png",sep="")
ggplot(plotDF, aes(Jday, Fit),
  ) + geom_polygon(data=dens,
                 aes(Day,Density),
                 fill=4,
                 alpha=0.2
  ) + geom_smooth(fill = "grey",
                colour = "black",
                aes(ymin=Jcil, ymax=Jciu),
                stat ="identity"
  ) + labs(x = "Julian Day",
         y = "s(Julian Day)",
         title = paste(CTname, 'at',sites[j]),
  ) + theme(axis.line = element_line(),
          panel.background = element_blank()
  )

ggsave(
  saveName,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device

# Year
# Center intercept (1st level of year factor) at 0 and show other levels relative to it
AdjustedYearCoefs = data.frame(
  c(
    YearBootstrapCoefs[, 1] - mean(YearBootstrapCoefs[, 1]),
    YearBootstrapCoefs[, 2],
    YearBootstrapCoefs[, 3],
    YearBootstrapCoefs[, 4]
  ),
  as.factor(rep(2016:2019, each = 10000))
)
colnames(AdjustedYearCoefs) = c("Coefficient", "Year")

saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot.png",sep = "")
ggplot(AdjustedYearCoefs, aes(Year, Coefficient)
       ) + geom_boxplot(
       ) + theme(axis.line = element_line(),
                panel.background = element_blank()
       ) + labs(title = paste(CTname, 'at', sites[j]))
ggsave(
  saveName,
  device = "png") # save figure
while (dev.cur() > 1) {
  dev.off()
} # close graphics device