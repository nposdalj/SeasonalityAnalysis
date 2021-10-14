#Pm_GEE_smooth_cycledata_160205
#Running the final set of models through SALSA. 

library(geepack)
require(MRSea)
library(calibrate)
library(splines)

#The raw data
cycledata <- read.csv('C:/Users/KMerkens/Documents/SpermWhales/ModelingData/CycleData_170613.csv', header = TRUE)

#####
#Generate the doy vector
days = cycledata$date - 719529   # 719529 = days from 1-1-0000 to 1-1-1970
secs = days * 86400 # 86400 seconds in a day
Rdate = as.POSIXct(secs, origin = '1970-1-1', tz = 'UTC')

# ##Extract year
# year <- strftime(Rdate, format = "%Y")
# cycledata$year <- as.numeric(year)

##Then covernt to doy
doy <- strftime(Rdate, format = "%j")
cycledata$doy <- as.numeric(doy)
# doystd <- cycledata$doy/366 #used 366, just to be able to deal with leap years
# cycledata$doy <- as.numeric(doystd)

#####
#Generate at time of day (tod) vector
dateonly <- floor(cycledata$date)
timeonly <- cycledata$date - dateonly
cycledata$tod <- timeonly

#####
#add colum with response called "response" 
cycledata$response <- cycledata$presabs

#####
#Set foldid for cross validation calculation 
cycledata$blockid<-as.factor(cycledata$x12_day)
#cycledata$foldid<-getCVids(cycledata, folds=5, block='blockid')

#########
##Change the order of the levels in the diel parameter so that 4, night, 
##is used as the reference. Does not alter the actual data.
cycledata$diel <- factor(cycledata$diel, levels = c(4,1,2,3))

########
##Subsample data
subdata_palt <- subset(cycledata, cycledata$regsite == "PALT")

#Remove the extra regsite levels that don't go away
subdata_palt$regsite <- droplevels(subdata_palt$regsite) 
subdata_palt$regsitedep <- droplevels(subdata_palt$regsitedep) 


#Make a foldid for this subset of data
subdata_palt$foldid<-getCVids(subdata_palt, folds=5, block='blockid')


#########
#Adjust lunar_dy by 15 days, to make start/end at day 15/16
#First, for values less than 15, add 16 days (that's 30-14, to bring it around to the end)
#start by making a new vector, Calling it Q1 b/c start is now at first quarter
subdata_palt$lunar_dy_F <- subdata_palt$lunar_dy
beforeF <-which(subdata_palt$lunar_dy < 15)
subdata_palt$lunar_dy_F[beforeF] <- subdata_palt$lunar_dy_F[beforeF]+16
#For values 7 and larger, just subtract 6, to bring everything down
afterF <-which(subdata_palt$lunar_dy >= 15)
subdata_palt$lunar_dy_F[afterF] <- subdata_palt$lunar_dy_F[afterF]-14




########################################
########################################

#Testing all possible models, including specific salsa for each
#Initial models:
#a NA
initialModel_b <- glm(response ~ as.factor(diel), data=subdata_palt, family = 'binomial')                      #diel, doy
#c NA
initialModel_d <- glm(response ~ 1, data=subdata_palt, family = 'binomial')                                    #doy only
#E, F NA

#Salsa lists:
# salsa1dlist <- list(fitnessMeasure = "QICb", minKnots_1d = c(1),
#                     maxKnots_1d = c(5), startKnots_1d = c(1), degree = c(1), maxIterations = 10,
#                     gaps = c(5))


salsa1dlist_doy_lu <- list(fitnessMeasure = "QICb", minKnots_1d = c(1,1),
                           maxKnots_1d = c(5,5), startKnots_1d = c(1,1), degree = c(1,1), maxIterations = 10,
                           gaps = c(5,5))



###################################
# ###################################
# #A - diel, regsite, doy, date
# #NA for single sites
# 
# 
###################################
#B - diel, doy, lunar (removal = F for PLOTTING ONLY)
salsa1dOutput_b <- runSALSA1D_withremoval(initialModel_b, salsa1dlist_doy_lu, varlist = c('doy','lunar_dy_F'),
                                          varlist_cyclicSplines = c('doy','lunar_dy_F'),factorlist=c('diel'),
                                          predictionData=NULL,datain=subdata_palt, removal=F)

splineParams<-salsa1dOutput_b$splineParams
baseModel_b <- salsa1dOutput_b$bestModel
baseModel_b$formula
# response ~ as.factor(diel)

gee_b <- geeglm(formula(baseModel_b),
                data = subdata_palt, family = binomial,id = blockid, corstr = 'ar1')
QICb(gee_b)
# 1612.295

# ###
# #Plotting and p-values for explanatory model
#
library(ggplot2)
source('~/R/MRSea/KarlisRunPartialPlots.R')
KarlisRunPartialPlots(model = gee_b, data = subdata_palt, varlist = c('doy','lunar_dy_F'),
                      factorlist = c('diel'), save = T, showKnots = F)


# 
# 
# 
# #Get p values, and save the output to a txt file, in case r crashes or the computer shuts down.
# sink('PALT-Pvalues-dieldoy_180402.txt')
# getPvalues(gee_b, factorlist = c('diel'), varlist = c('doy'))
# sink()
# # [1] "Getting marginal p-values"
# # Variable  p-value
# # 1     diel 0.407382
# # 2      doy 0.786679



# ###################################
# #C - regsite, doy, date
# #NA for single sites

# 
# ###################################
# #D - doy, lunar
salsa1dOutput_d <- runSALSA1D_withremoval(initialModel_d, salsa1dlist_doy_lu, varlist = c('doy','lunar_dy_F'),
                                          varlist_cyclicSplines = c('doy','lunar_dy_F'),
                                          predictionData=NULL,datain=subdata_palt, removal=TRUE)

splineParams<-salsa1dOutput_d$splineParams
baseModel_d <- salsa1dOutput_d$bestModel
baseModel_d$formula
# response ~ 1

gee_d <- geeglm(formula(baseModel_d),
                data = subdata_palt, family = binomial,id = blockid, corstr = 'ar1')
QICb(gee_d)
# 1582.555

# 



###################################
###################################
##Winning model: F date 

#check
#plot
#other output

plotCumRes(gee_j, varlist= c("date"), splineParams=splineParams)

# #Not plotting the runs because I know there's correlation
# library(lawstat)
# runs.test(residuals(gee_e, type = "pearson"), alternative = c("two.sided"))
# plotRunsProfile(gee_e, varlist= c("doy","date"))


library(ggplot2)
#runPartialPlots(model = gee_f, data = subdata_palt,varlist= c("date"), showKnots = T)
source('~/R/MRSea/KarlisRunPartialPlots.R')
KarlisRunPartialPlots(model = gee_j, data = subdata_palt, varlist = c('date'),
                      save = T, showKnots = F)


runDiagnostics(gee_j)



#timeInfluenceCheck(gee_f, subdata_palt$blockid, splineParams, d2k = NULL)
influence<-runInfluence(gee_j, subdata_palt$blockid, splineParams, d2k = NULL)


#Get p values, and save the output to a txt file, in case r crashes or the computer shuts down. 
sink('PALT-Pvalues_160804.txt')
getPvalues(gee_j,varlist= c("date"))
sink()
# 
# [1] "Getting marginal p-values"
# Variable  p-value
# 1     date 0.039533



