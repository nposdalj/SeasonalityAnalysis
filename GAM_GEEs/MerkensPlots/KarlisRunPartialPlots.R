## The code from Karli in the main script should be...
#library(ggplot2)
#runPartialPlots(model = gee_f, data = subdata_palt,varlist= c("date"), showKnots = T)
#source('~/R/MRSea/KarlisRunPartialPlots.R')
#KarlisRunPartialPlots(model = gee_j, data = subdata_palt, varlist = c('date'),
                      #save = T, showKnots = F)

KarlisRunPartialPlots <- function (model, data, factorlist = NULL, varlist = NULL, showKnots = FALSE, 
          save = FALSE, savedata = F) 
{
  print("Making partial plots")
  require(mvtnorm)
  require(splines)
  require(Matrix)
  if (save == T) {
    png("PartialFitsLink%i.png", height = 550, width = 600)
  }
  else {
    devAskNewPage(ask = TRUE)
  }
  par(mfrow = c(1, 1))
  if (is.null(factorlist) == F) {
    for (i in 1:length(factorlist)) {
      coeffac <- c(grep(factorlist[i], colnames(model.matrix(model))))
      coefradial <- c(grep("LocalRadialFunction", colnames(model.matrix(model))))
      coefpos <- coeffac[which(is.na(match(coeffac, coefradial)))]
      xvals <- data[, which(names(data) == factorlist[i])]
      newX <- sort(unique(xvals))
      #newX <- newX[2:length(newX)] #try removing this truncation
      #partialfit <- coef(model)[coefpos]
      partialfit <- c(1E-20,coef(model)[coefpos]) #Added 0 for the reference level
      rcoefs <- NULL
      try(rcoefs <- rmvnorm(1000, coef(model), summary(model)$cov.scaled), 
          silent = T)
      if (is.null(rcoefs) || length(which(is.na(rcoefs) == 
                                            T)) > 0) {
        rcoefs <- rmvnorm(1000, coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))
      }
      rpreds <- rcoefs[, coefpos]
      quant.func <- function(x) {
        quantile(x, probs = c(0.025, 0.975)) #Calculates the probability of 0.025 tails
      }
      ###This first part is when there is only one x-value to calc CIs for (e.g. two factors, one being reference)
#       if (is.null(dim(rpreds))) {
#         cis <- quant.func(rpreds)
#         
#         referencecis <- c(0,0)
#         names(referencecis) = c("2.5%", "97.5%")
#         cis <- rbind(referencecis,cis,deparse.level = 3)#Add a 0 0 to the start of cis as the confidence interval for the reference level.
#         
#         plot(newX, partialfit, pch = 20, xlab = factorlist[i], 
#              ylab = "Partial Fit", lwd = 2, xaxt = "n", 
#              ylim = c(range(0, cis)), xlim=c(0,length(newX)),
#              cex.lab = 1.3, cex.axis = 1.3)
#         #axis(1, at = newX, labels = colnames(model.matrix(model))[coefpos])
#         axis(1, at = (1:length(newX)), labels = newX)
#         segments(newX, cis[1], newX, cis[2], lwd = 2)
#         abline(h = 0, lwd = 2)
#       }
      if (is.null(dim(rpreds))) {
        cis <- quant.func(rpreds)
      }
      else {
        cis <- t(apply(rpreds, 2, quant.func))
      }
        referencecis <- c(0,0)
        names(referencecis) = c("2.5%", "97.5%")
        cis <- rbind(referencecis,cis,deparse.level = 3)#Add a 0 0 to the start of cis as the confidence interval for the reference level.
        
        #Format x labels
        if (factorlist[i] == 'preamp'){
          #xlab <- c('Instrument ID')
          newX <- gsub("[^[:digit:]]", "", newX)
        }
        else if (factorlist[i] == 'diel'){
          #xlab <- c('Diel Phase')
          newX <- c('Night','Dawn','Day','Dusk')
        }
        else if (factorlist[i] == 'subreg'){
          #xlab <- c('Sub-region')
        }
        else if (factorlist[i] == 'regsite'){
          #xlab <- c('Site')
          #take the first 3 characters of each regsite name, except LSMS, D, PALW, N
          levels(newX)[levels(newX)=="CSMA"] <- "CSM"
          levels(newX)[levels(newX)=="EQUA"] <- "EQU"
          levels(newX)[levels(newX)=="HAWK"] <- "HAW"
          levels(newX)[levels(newX)=="KAUA"] <- "KAU"
          levels(newX)[levels(newX)=="KINA"] <- "KIN"
          levels(newX)[levels(newX)=="LSMD"] <- "LSM-D"
          levels(newX)[levels(newX)=="LSMS"] <- "LSM-S"
          levels(newX)[levels(newX)=="PALS"] <- "PAL-N"
          levels(newX)[levels(newX)=="PALT"] <- "PAL-W"
          levels(newX)[levels(newX)=="PHRA"] <- "PHR"
          levels(newX)[levels(newX)=="SAIA"] <- "SAI"
          levels(newX)[levels(newX)=="TINA"] <- "TIN"
          levels(newX)[levels(newX)=="WAKS"] <- "WAK"
          newX <- c(levels(newX))
        }
        else {
          xlab <- factorlist[i]
        }

        #Change ylim range to 0.5 if needed to zoom out, 0.05 if zoom in
        plot(1:length(newX), partialfit, pch = 15, cex = 2,
             xlab = "", ylab = "",
             lwd = 2, xlim = c(0.5,length(newX)+0.5), xaxt = "n",
             ylim = c(range(0.05, cis)), cex.lab = 1.3, cex.axis = 1.3)
        axis(1, at = c(0,length(newX)+1), labels = c("",""), lwd.ticks=0)
        axis(1, at = c(1:length(newX)), labels = c(newX))
        segments(1:length(newX), cis[, 1], 1:length(newX), 
                 cis[, 2], lwd = 2)
        abline(h = 0, lwd = 1)
        ##cut ,xlab = xlab, ylab = "Partial Fit (95% CI)"
      }
      if (savedata == T) {
        partialdata <- data.frame(newX, partialfit, cis[1], 
                                  cis[2])
        save(partialdata, file = paste("PartialData_", 
                                       factorlist[i], ".RData", sep = ""), compress = "bzip2")
      }
  }
  if (is.null(varlist) == F) {
    n <- length(varlist)
    for (i in 1:n) {
      coefpos <- c(1, grep(varlist[i], colnames(model.matrix(model))))
      xvals <- data[, which(names(data) == varlist[i])]
      newX <- seq(min(xvals), max(xvals), length = 500)
      eval(parse(text = paste(varlist[i], "<- newX", sep = "")))
      response <- rep(1, 500)
      newBasis <- eval(parse(text = labels(terms(model))[grep(varlist[i], 
                                                              labels(terms(model)))]))
      partialfit <- cbind(rep(1, 500), newBasis) %*% coef(model)[coefpos]
      rcoefs <- NULL
      try(rcoefs <- rmvnorm(1000, coef(model), summary(model)$cov.scaled), 
          silent = T)
      if (is.null(rcoefs) || length(which(is.na(rcoefs) == 
                                            T)) > 0) {
        rcoefs <- rmvnorm(1000, coef(model), as.matrix(nearPD(summary(model)$cov.scaled)$mat))
      }
      rpreds <- cbind(rep(1, 500), newBasis) %*% t(rcoefs[, 
                                                          coefpos])
      quant.func <- function(x) {
        quantile(x, probs = c(0.025, 0.975))
      }
      cis <- t(apply(rpreds, 1, quant.func))

      #Format x labels
      if (varlist[i] == 'doy'){
        #xlab <- c('Day of the Year')
        ticks <- c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC')
        xlength <- c(0,366)
        tickspace <- seq(15,366,by=30)
      }
      else if (varlist[i] == 'date'){
        #xlab <- c('Date')
        if (data$regsite[1]=='HAWK') {  
          ##8/1/07 733255, 8/15/07 733269, 1/1/09 733774, 1/1/11 734504, 1/1/13 735235, 10/15/13 735522, 12/15/16 735583
          ticks <- c("Aug 2007","Jan 2009","Jan 2011", "Jan 2013","Oct 2013" )
          xlength <- c(733255,735583)
          tickspace <- c(733269,733774,734504,735235,735522)
        }
        else if (data$regsite[1]=='PALT') {
          ##10/01/06 732951, 10/15/06 732965, 1/1/07 733043, 1/1/08 733408, 1/1/09 733774, 4/1/09 733864, 4/15/09 733878
          ticks <- c("Oct 2006","Jan 2007","Jan 2008", "Jan 2009","Apr 2009" )
          xlength <- c(732951,733878)
          tickspace <- c(732965,733043,733408,733774,733864)
        } 
        else if (data$regsite[1]=='PHRA') {
          ##10/01/09 , 10/15/06 732965, 1/1/12 , 1/1/08 733408, 1/1/09 733774, 4/1/09 733864, 4/15/09 733878
          ##Start 734066, End 734875
          ticks <- c("Oct 2009","Jan 2011","Jan 2012")
          xlength <- c(734066,734875)
          tickspace <- c(734080,734504,734860)
        } 
      }
      else if(varlist[i] == 'lunar_dy_F'){
        ticks <- c('FULL','NEW','FULL') #seq(1,29, by = 1)
        xlength <- c(0,30)
        tickspace <- c(1,15,30)
      }
      else {
        xlab <- varlist[i]
#         ylength <- 
      }
      plot(newX, partialfit, type = "l", lwd = 2, ylim = range(cis),
           xlab = "", ylab = "",
           xlim = xlength, xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
      axis(1, at = xlength, labels = c("",""), lwd.ticks=0)
      axis(1, at = tickspace, labels = ticks, lwd.ticks=0)
      ##Cut , xlab = xlab,ylab = "Partial Fit (95% CI)"

      if(varlist[i] == 'lunar_dy_F'){
      rug(jitter(xvals))
      }
      else {
        rug(xvals)
      }

      # add fill
      polygon(c(rev(newX), newX), c(rev(cis[, 1]), cis[, 2]), col = 'grey80', border = NA)
      # intervals
      # lines(newX, cis[, 1], lty = 'dashed', col = 'red')
      # lines(newX, cis[, 2], lty = 'dashed', col = 'red')
      lines(newX, partialfit, type = "l", lwd = 2, ylim = range(cis),
           xlab = "", ylab = "",
           xlim = xlength, xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
      
      #old dotted error lines
      #lines(newX, cis[, 1], col = "black", lty = 4)
      #lines(newX, cis[, 2], col = "black", lty = 4)
      if (showKnots == "TRUE") {
        if (length(splineParams[[(i + 1)]]$knots) == 
              1 & is.character(splineParams[[(i + 1)]]$knots)) {
          XX <- FALSE
        }
        else {
          XX <- TRUE
        }
        if (XX == TRUE) {
          abline(v = splineParams[[(i + 1)]]$knots, lty = 4, 
                 col = "grey")
        }
      }
      eval(parse(text = paste("rm(", varlist[i], ")", sep = "")))
      rm(response)
      if (savedata == T) {
        partialdata <- data.frame(newX, partialfit, cis[, 
                                                        1], cis[, 2])
        save(partialdata, file = paste("PartialData_", 
                                       varlist[i], ".RData", sep = ""), compress = "bzip2")
      }
    }
    
  }
  if (save == T) {
    dev.off()
    #If there are 3 plots, uncomment this
    ##dev.off()
  }
  else {
    devAskNewPage(ask = FALSE)
  }
if (save == T) {
  print('If you see an error, it is fine, the script is done, the plots should have been made successfully')
  dev.off()
  
}
}