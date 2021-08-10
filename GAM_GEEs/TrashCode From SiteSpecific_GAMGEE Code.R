##Trash code from SiteSpecific_GAMGEE.R code

#Calculating Block Size the Merkens Way
startDate = SiteDayTable$tbin[1]
endDate = SiteDayTable$tbin[nrow(SiteDayTable)]
timeseries = data.frame(date=seq(startDate, endDate, by="days"))
timeseries$one = 1:nrow(timeseries)
oneday = left_join(SiteHourTable,timeseries,by="date")
onedaygrouped = aggregate(oneday[, c(2,8)], list(oneday$one), mean)
onedaygrouped$Group.1 = as.factor(onedaygrouped$Group.1)
onedaygrouped = onedaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))

timeseries$two = rep(1:(nrow(timeseries)/2), times=1, each=2)
twoday = left_join(SiteHourTable,timeseries,by="date")
twodaygrouped = aggregate(twoday[, c(2,8)], list(twoday$two), mean)
twodaygrouped$Group.1 = as.factor(twodaygrouped$Group.1)
twodaygrouped = twodaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))
acf(twodaygrouped$PreAbs, plot = FALSE)
acf(twodaygrouped$PreAbs)
timeseries$two = NULL
#autocorrelated

three = rep(1:(floor(nrow(timeseries)/3)), times=1, each=3)
timeseries$three = c(three,three[3402]+1,three[3402]+1)
threeday = left_join(HourTable,timeseries,by="date")
threedaygrouped = aggregate(threeday[, c(2,8)], list(threeday$three), mean)
threedaygrouped$Group.1 = as.factor(threedaygrouped$Group.1)
threedaygrouped = threedaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))
acf(threedaygrouped$PreAbs, plot = FALSE)
acf(threedaygrouped$PreAbs)
timeseries$three = NULL
#autocorrelated

timeseries$four = rep(1:(floor(nrow(timeseries)/4)), times=1, each=4)
fourday = left_join(HourTable,timeseries,by="date")
fourdaygrouped = aggregate(fourday[, c(2,8)], list(fourday$four), mean)
fourdaygrouped$Group.1 = as.factor(fourdaygrouped$Group.1)
fourdaygrouped = fourdaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))
acf(fourdaygrouped$PreAbs, plot = FALSE)
acf(fourdaygrouped$PreAbs)
timeseries$four = NULL
#autocorrelated

five = rep(1:(floor(nrow(timeseries)/5)), times=1, each=5)
timeseries$five = c(five,five[3400]+1,five[3400]+1,five[3400]+1,five[3400]+1)
fiveday = left_join(HourTable,timeseries,by="date")
fivedaygrouped = aggregate(fiveday[, c(2,8)], list(fiveday$five), mean)
fivedaygrouped$Group.1 = as.factor(fivedaygrouped$Group.1)
fivedaygrouped = fivedaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))
acf(fivedaygrouped$PreAbs, plot = FALSE)
acf(fivedaygrouped$PreAbs)
timeseries$five = NULL
#autocorrelated

six = rep(1:(floor(nrow(timeseries)/6)), times=1, each=6)
timeseries$six = c(six,six[3402]+1,six[3402]+1)
sixday = left_join(HourTable,timeseries,by="date")
sixdaygrouped = aggregate(sixday[, c(2,8)], list(sixday$six), mean)
sixdaygrouped$Group.1 = as.factor(sixdaygrouped$Group.1)
sixdaygrouped = sixdaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))
acf(sixdaygrouped$PreAbs, plot = FALSE)
acf(sixdaygrouped$PreAbs)
timeseries$six = NULL
#autocorrelated

seven = rep(1:(floor(nrow(timeseries)/7)), times=1, each=7)
timeseries$seven = c(seven,seven[3402]+1,seven[3402]+1)
sevenday = left_join(HourTable,timeseries,by="date")
sevendaygrouped = aggregate(sevenday[, c(2,8)], list(sevenday$seven), mean)
sevendaygrouped$Group.1 = as.factor(sevendaygrouped$Group.1)
sevendaygrouped = sevendaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))
acf(sevendaygrouped$PreAbs, plot = FALSE)
acf(sevendaygrouped$PreAbs)
timeseries$seven = NULL
#autocorrelated

eight = rep(1:(floor(nrow(timeseries)/8)), times=1, each=8)
timeseries$eight = c(eight,eight[3400]+1,eight[3400]+1,eight[3400]+1,eight[3400]+1)
eightday = left_join(HourTable,timeseries,by="date")
eightdaygrouped = aggregate(eightday[, c(2,8)], list(eightday$eight), mean)
eightdaygrouped$Group.1 = as.factor(eightdaygrouped$Group.1)
eightdaygrouped = eightdaygrouped %>% mutate_if(is.numeric, ~1 * (. != 0))
acf(eightdaygrouped$PreAbs, plot = FALSE)
acf(eightdaygrouped$PreAbs)
timeseries$eight = NULL
#no autocorrrelation at this point

nine = rep(1:(floor(nrow(timeseries)/9)), times=1, each=9)
timeseries$nine = c(nine,nine[3402]+1,nine[3402]+1)

ten = rep(1:(floor(nrow(timeseries)/10)), times=1, each=10)
timeseries$ten = c(ten,ten[3400]+1,ten[3400]+1,ten[3400]+1,ten[3400]+1)

eleven = rep(1:(floor(nrow(timeseries)/11)), times=1, each=11)
timeseries$eleven = c(eleven,eleven[3399]+1,eleven[3399]+1,eleven[3399]+1,eleven[3399]+1,eleven[3399]+1)

twelve = rep(1:(floor(nrow(timeseries)/12)), times=1, each=12)
timeseries$twelve = c(twelve,twelve[3396]+1,twelve[3396]+1,twelve[3396]+1,twelve[3396]+1,twelve[3396]+1,twelve[3396]+1,twelve[3396]+1,twelve[3396]+1)

thirteen = rep(1:(floor(nrow(timeseries)/13)), times=1, each=13)
timeseries$thirteen = c(thirteen,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1,thirteen[3393]+1)

fourteen = rep(1:(floor(nrow(timeseries)/14)), times=1, each=14)
timeseries$fourteen = c(fourteen, fourteen[3402]+1, fourteen[3402]+1)

fifteen = rep(1:(floor(nrow(timeseries)/15)), times=1, each=15)
timeseries$fifteen = c(fifteen, fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1,fifteen[3390]+1)

# Group everything in blocks from 1d to 15d (15 different data sets)
HourTableBinned = left_join(HourTable,timeseries,by = "date")
HourTableBinned = HourTableBinned[ order(HourTableBinned$tbin , decreasing = FALSE ),]

#Step 7: Construction of the ROC curve
# STEP 6: Construction of the ROC curve     
pr <- predict(POD2c, type="response")  
pred <- prediction(pr,SiteHourTableB$PreAbs) 
perf <- performance(pred, measure="tpr", x.measure="fpr")   
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5))
#This creates a ROC plot

# Choice of the best cut-off probability

y<-as.data.frame(perf@y.values)
x<-as.data.frame(perf@x.values)
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
fileName2 = paste(dir,'/',site,'distances.txt',sep="")
write.table(d,fileName2)                                 # to write a table with the computed distances

# The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
# MAX d= 0.147907576861441 --> position 999

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
alpha[999,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve
#[1] 0.4422233
# Best cutoff:     0.4422233
# This value can now be used to build the confusion matrix:

DATA<-matrix(0,45354 ,9)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim())
dim(SiteHourTableB)
#CB
#[1] 45354     9
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:45354                                    # the first column is filled with an ID value that is unique for each row
DATA$Observed<-SiteHourTableB$PreAbs                                           # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(POD3,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = 0.4422233)                                   # the identified cut-off must be used here

#Confusion matrix:
#observed
#predicted     1     0
#1          13115 10693
#0           7387 14159

# The confusion matrix can then be transformed into percentages:
# The area under the curve (auc) can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")
#auc = 0.6623965

#Step 8: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the relationship between the response (on the response scale) and each predictor ##
POD3 = geeglm(PreAbs ~ AvgDayMat,family = binomial, corstr="ar1", id=Blocks, data=SiteHourTableB)
dimnames(AvgDayMat)<-list(NULL,c("AHBM1", "AHBM2", "AHBM3", "AHBM4"))

#Probability of covariate #1: AvgDayMat:
BootstrapParameters1<-rmvnorm(10000, coef(POD3),summary(POD3)$cov.unscaled, method="chol")
fstart=2; finish=5; Variable=AvgDayMat; xlabel="Julian Day"; ylabel="Probability"  
PlottingVar1<-seq(min(Variable), max(Variable), length=5000)
CenterVar1<-model.matrix(POD3)[,start:finish]*coef(POD3)[c(start:finish)]
BootstrapCoefs1<-BootstrapParameters1[,c(start:finish)]
Basis1<-gam(rbinom(5000,1,0.5)~s(PlottingVar1, bs="tp", k=6), fit=F, family=binomial, knots=list(PlottingVar1=seq(0,61,length=6)))$X[,2:5]
RealFit1<-Basis1%*%coef(POD8)[c(start:finish)]
RealFitCenter1<-RealFit1-mean(CenterVar1)
RealFitCenter1a<-inv.logit(RealFitCenter1)
BootstrapFits1<-Basis1%*%t(BootstrapCoefs1)
quant.func1<-function(x){quantile(x,probs=c(0.025, 0.975))}
cis1<-apply(BootstrapFits1, 1, quant.func1)-mean(CenterVar1)
cis1a<-inv.logit(cis1)
plot(PlottingVar1,(RealFitCenter1a), type="l", col=1,ylim=c(0, 1),xlab=xlabel, ylab=ylabel, xlim=c(0,61), main ="M7 covariate 1: # Days since deployment" , cex.lab = 1.5, cex.axis=1.5)
segments(PlottingVar1,(cis1a[1,]),PlottingVar1,(cis1a[2,]), col="grey", main = "Influence of DATENO")
lines(PlottingVar1,(RealFitCenter1a),lwd=2, col=1)
rug(PlottingVar1)


## STEP 4: fit the full model ##

fit1<-geeglm(PreAbs ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  

# Exploratory data analysis: summarise the values of the covariates for the response "Pres". Boxplots could also be used for a visual inspection of the data.

tapply(dat$Depth,dat$Pres,summary)
boxplot(dat$Depth~dat$Pres)
tapply(dat$Slope10x,dat$Pres,summary)
boxplot(dat$Slope10x~dat$Pres)
tapply(dat$Aspect,dat$Pres,summary)
boxplot(dat$Aspect~dat$Pres)
tapply(dat$wind,dat$Pres,summary)
boxplot(dat$wind~dat$Pres)
tapply(dat$SSH,dat$Pres,summary)
boxplot(dat$SSH~dat$Pres)
tapply(dat$SST_weekly,dat$Pres,summary)
boxplot(dat$SST_weekly~dat$Pres)
tapply(dat$SST_weekly_slope,dat$Pres,summary)
boxplot(dat$SST_weekly_slope~dat$Pres)
tapply(dat$SST_mdeviation,dat$Pres,summary)
boxplot(dat$SST_mdeviation~dat$Pres)
tapply(dat$chla_apr3x10,dat$Pres,summary)
boxplot(dat$chla_apr3x10~dat$Pres)

# Selection of the correct form (linear or smooth) for the covariates available at a single scale.
# A series of models is fitted where each of these covariate is tested as a linear term. 

fitdl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Depth+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fital<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+Aspect+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fitwl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+wind+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fithl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fitwsl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fitmdl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+bs(SSH,knots=mean(SSH))+SST_weekly+SST_mdeviation+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)

model<-c("full","dl","al","wl","hl","wsl","mdl")
QIC<-c(QIC(fit1)[1],QIC(fitdl)[1],QIC(fital)[1],QIC(fitwl)[1],QIC(fithl)[1],QIC(fitwsl)[1],QIC(fitmdl)[1])
data.frame(rbind(model,QIC))

#Their results
# model             full               dl               al              wl               hl              wsl              mdl
# QIC   5271.74422487162 9529.79249364125 5308.29218593693 9328.4962017818 5236.25469756788 5250.70542492964 5266.88171399961

# SSH should be retained as linear term.

#My results
#QIC            QIC.1            QIC.2           QIC.3            QIC.4            QIC.5            QIC.6
#model             full               dl               al              wl               hl              wsl              mdl
#QIC   5490.52655845423 5671.06108188895 5513.68282434987 5504.4825619765 5405.55929711371 5462.45394279189 5464.39159151921
# SSH should be retained as linear term. -- same as them

fit2<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fithl
fitdl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Depth+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fital<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+Aspect+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fitwl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+wind+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fitwsl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)
fitmdl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+bs(SST_weekly_slope,knots=mean(SST_weekly_slope))+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)

model<-c("full","dl","al","wl","wsl","mdl")
QIC<-c(QIC(fit2)[1],QIC(fitdl)[1],QIC(fital)[1],QIC(fitwl)[1],QIC(fitwsl)[1],QIC(fitmdl)[1])
data.frame(rbind(model,QIC))

#Their results
# model             full               dl               al               wl              wsl              mdl
# QIC   5236.25469756788 5359.68898267533 5266.71325970759 5257.48584805803 5214.92397199749 9334.61689900357

# SST_weekly_slope should be retained as a linear term

#My results
#QIC            QIC.1            QIC.2            QIC.3           QIC.4            QIC.5
#model             full               dl               al               wl             wsl              mdl
#QIC   5405.55929711371 5550.02753936587 5424.85324129154 5426.86240317372 5380.6578960802 5390.37850303909

# SST_weekly_slope should be retained as a linear term -- same as them

fit3<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
fitdl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Depth+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat) 
fital<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+Aspect+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat) 
fitwl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+wind+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat) 
fitmdl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat) 

model<-c("full","dl","al","wl","mdl")
QIC<-c(QIC(fit3)[1],QIC(fitdl)[1],QIC(fital)[1],QIC(fitwl)[1],QIC(fitmdl)[1])
data.frame(rbind(model,QIC))

#Their results
# model             full               dl               al               wl              mdl
# QIC   5214.92397199749 5339.35752588423 5244.56148984484 5236.70700508862 5221.49444337552

# All the QIC scores are higher than the score for the full model, therefore all the remaining covariates are considered in the best form (linear or smooth).

#My results
#QIC           QIC.1            QIC.2            QIC.3            QIC.4
#model            full              dl               al               wl              mdl
#QIC   5380.6578960802 5530.1650306503 5398.60488348528 5404.34671482125 5364.88586866761

# SST_mdeviation should be retained as a linear term -- different results form them

fit4<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
fitdl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Depth+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat) 
fital<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+Aspect+bs(wind,knots=mean(wind))+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat) 
fitwl<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+wind+SSH+SST_weekly+bs(SST_mdeviation,knots=mean(SST_mdeviation))+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat) 

model<-c("full","dl","al","wl")
QIC<-c(QIC(fit4)[1],QIC(fitdl)[1],QIC(fital)[1],QIC(fitwl)[1])
data.frame(rbind(model,QIC))

#My results
#QIC           QIC.1            QIC.2            QIC.3
#model             full              dl               al               wl
#QIC   5364.88586866761 5530.1650306503 5398.60488348528 5404.34671482125
# All the QIC scores are higher than the score for the full model, therefore all the remaining covariates are considered in the best form (linear or smooth).


## STEP 5: perform a manual stepwise variable selection to identify the best subset of covariates ##
# The initial full model is:

fit4<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl

# A series of reduced models is fitted. Each contains all the covariates but one. The reduced model with the lowest QIC is the one to use in the following step.

y<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
d<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sl10<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
a<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
w<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
ssh<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sst<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
mdev<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
ws<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
chla<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ as.factor(Year)+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope,family=binomial, corstr="independence",id=Line_Id, data = dat)  

model<-c("full","y","d","sl10","a","w","ssh","sst","mdev","ws","chla")
QIC<-c(QIC(fit4)[1],QIC(y)[1],QIC(d)[1],QIC(sl10)[1],QIC(a)[1],QIC(w)[1],QIC(ssh)[1],QIC(sst)[1],QIC(mdev)[1],QIC(ws)[1],QIC(chla)[1])
data.frame(rbind(model,QIC))

#Their results

# model             full                y                d             sl10                a                w              ssh              sst
# QIC   5214.92397199749 5001.85722884214 5314.56599253774 5268.89315991936 5259.84066647466 5218.35811737537 5174.54748704375 5206.05871992886

# model             mdev               ws             chla
# QIC   5192.10072346801 5208.74717055123 5221.74009915839

# Year has to be removed (the model without Year is the one with the lowest QIC).

#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5            QIC.6            QIC.7           QIC.8            QIC.9           QIC.10
#model             full                y                d             sl10                a                w              ssh              sst            mdev               ws             chla
#QIC   5364.88586866761 5112.24106987531 5457.06404160795 5396.92416579089 5365.81175128658 5372.97544962279 5315.54385893167 5295.49937820926 5333.6105941239 5346.08586150156 5395.69076709248
# Year has to be removed (the model without Year is the one with the lowest QIC). -- same as them

# The procedure must be continued until none of the remaining variables can be removed without increasing QIC.

fit5<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
d<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sl10<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
a<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ +bs(Depth,knots=mean(Depth))+Slope10x+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
w<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
ssh<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sst<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
mdev<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
ws<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
chla<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SSH+SST_weekly+SST_mdeviation+SST_weekly_slope,family=binomial, corstr="independence",id=Line_Id, data = dat)  

model<-c("full","d","sl10","a","w","ssh","sst","mdev","ws","chla")
QIC<-c(QIC(fit5)[1],QIC(d)[1],QIC(sl10)[1],QIC(a)[1],QIC(w)[1],QIC(ssh)[1],QIC(sst)[1],QIC(mdev)[1],QIC(ws)[1],QIC(chla)[1])
data.frame(rbind(model,QIC))

#Their results
# model             full               d             sl10                a                w              ssh              sst             mdev
# QIC   5001.85722884214 5077.2429394218 5060.53097754001 5039.22898604842 5001.21999935888 4976.20576420094 5027.24395026128 4959.08214495535
# model               ws             chla
# QIC   4989.29621591137 5003.10646963412
# SST_mdeviation can be removed.

#My results
#QIC           QIC.1            QIC.2            QIC.3            QIC.4            QIC.5            QIC.6            QIC.7            QIC.8           QIC.9
#model             full               d             sl10                a                w              ssh              sst             mdev               ws            chla
#QIC   5112.24106987531 5183.2864605278 5145.84004343482 5116.60420433129 5098.58793709458 5075.33979097764 5121.58047104614 5075.88354271824 5090.15955879618 5142.8832718297
# SSH can be removed.

fit6<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
d<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sl10<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
a<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ +bs(Depth,knots=mean(Depth))+Slope10x+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
w<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sst<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
mdev<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
ws<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
chla<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+bs(wind,knots=mean(wind))+SST_weekly+SST_mdeviation+SST_weekly_slope,family=binomial, corstr="independence",id=Line_Id, data = dat)  

model<-c("full","d","sl10","a","w","sst","ws","chla")
QIC<-c(QIC(fit6)[1],QIC(d)[1],QIC(sl10)[1],QIC(a)[1],QIC(w)[1],QIC(sst)[1],QIC(ws)[1],QIC(chla)[1])
data.frame(rbind(model,QIC))

#Their results
# model             full                d             sl10                a                w            ssh              sst              ws
# QIC   4959.08214495535 5034.41802080704 5005.09475080114 4993.23054794395 4956.33987701265 4931.638982768 4975.94044305331 4945.8335460295
# model             chla
# QIC   4960.72196547842
# SSH can be removed.

#My results
#QIC           QIC.1            QIC.2            QIC.3            QIC.4            QIC.5            QIC.6            QIC.7
#model             full               d             sl10                a                w              sst               ws             chla
#QIC   5075.33979097764 5138.7771780659 5107.35866056379 5078.33952821448 5039.77913132179 5085.13565572294 5053.60505186332 5122.59696306518
#Wind can be removed

fit7<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
d<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sl10<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
a<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ +bs(Depth,knots=mean(Depth))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sst<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
mdev<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
ws<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
chla<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+SST_weekly_slope,family=binomial, corstr="independence",id=Line_Id, data = dat)  

model<-c("full","d","sl10","a","sst","ws","chla")
QIC<-c(QIC(fit7)[1],QIC(d)[1],QIC(sl10)[1],QIC(a)[1],QIC(sst)[1],QIC(ws)[1],QIC(chla)[1])
data.frame(rbind(model,QIC))

#Their results
# model           full               d             sl10                a                w              sst             ws             chla
# QIC   4931.638982768 4997.0163528512 4976.04665219566 4960.98551817401 4907.16918526904 4952.86373685549 4919.284606216 4939.77543563068
# wind has to be removed.

#My results
#QIC            QIC.1           QIC.2            QIC.3            QIC.4            QIC.5            QIC.6
#model             full                d            sl10                a              sst               ws             chla
#QIC   5039.77913132179 5085.36628582922 5056.3274668178 5064.06083966443 5042.96710917149 5020.56360670406 5083.15870314743
# SST_weekly_slope has to be removed.

fit8<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
d<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+SST_weekly_slope+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sl10<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
a<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ +bs(Depth,knots=mean(Depth))+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
sst<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
mdev<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  
chla<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+ bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation,family=binomial, corstr="independence",id=Line_Id, data = dat)  

model<-c("full","d","sl10","a","sst","chla")
QIC<-c(QIC(fit8)[1],QIC(d)[1],QIC(sl10)[1],QIC(a)[1],QIC(sst)[1],QIC(chla)[1])
data.frame(rbind(model,QIC))

#Their results
# model             full               d             sl10                a              sst               ws             chla
# QIC   4907.16918526904 4997.0163528512 4935.05818149937 4930.88179159234 4923.56155384573 4896.03099630217 4917.26036903389
# SST_weekly_slope has to be removed.

#My results
#QIC            QIC.1            QIC.2            QIC.3            QIC.4            QIC.5
#model             full                d             sl10                a              sst             chla
#QIC   5020.56360670406 5085.36628582922 5036.01789213509 5044.21984440794 5023.11432858729 5062.28986132861

# All the covariates, when removed, cause the QIC score to increase.
# At this point, the resulting model is fitted using the library geeglm. The order in which the covariates enter the model is determined by the QIC score
# (the ones that, if removed, determine the biggest increase in QIC enter the model first).

fit8<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation+chla_apr3x10,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl

# The significance of the covariates is tested using Wald's tests (i.e. the anova.geeglm function of the library geeglm). Non-significant covariates are removed one at a time.

anova(fit8)

# - chla_apr3x10 (p=0.3340)
#My results - chla_apr3x10 (p=0.3343)

fit9<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly+SST_mdeviation,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
anova(fit9)

# - SST_weekly (p=0.1674)
#My results - SST_mdeviation (p = 0.6447016)

fit10<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect))+SST_weekly,family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
anova(fit10)

# - Slope10x (p=0.0917)
#My results - SST_weekly (p = 0.167)

fit11<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+Slope10x+bs(Aspect,knots=mean(Aspect)),family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
anova(fit11)

# Their results - fit11 represents the final model.
#My results - Slope10x (p = 0.0917)

fit12<-geeglm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect)),family=binomial, corstr="independence",id=Line_Id, data = dat)  # i.e. fitwsl
anova(fit12)

# My results - fit12 represents the final model.

## STEP 6: construct the confusion matrix ##

# Construction of the ROC curve

pr <- predict(fit12,dat, type="response")                          # the final model is used to predict the data on the response scale (i.e. a value between 0 and 1)
pred <- prediction(pr,dat$Pres)                                    # to specify the vector of predictions (pr) and the vector of labels (i.e. the observed values "Pres")
perf <- performance(pred, measure="tpr", x.measure="fpr")          # to assess model performance in the form of the true positive rate and the false positive rate
plot(perf, colorize=TRUE, print.cutoffs.at=c(0.1,0.2,0.3,0.4,0.5)) # to plot the ROC curve

# Choice of the best cut-off probability

y<-as.data.frame(perf@y.values)
x<-as.data.frame(perf@x.values)
fi <- atan(y/x) - pi/4                                             # to calculate the angle between the 45° line and the line joining the origin with the point (x;y) on the ROC curve
L <- sqrt(x^2+y^2)                                                 # to calculate the length of the line joining the origin to the point (x;y) on the ROC curve
d <- L*sin(fi)                                                     # to calculate the distance between the 45° line and the ROC curve
#write.table(d,"C:\\distances.txt")                                 # to write a table with the computed distances

#What they did
# The table should then be opened in Microsoft Excel to find the maximum distance with the command "Sort", and the relative position (i.e. the number of the corresponding record)
# MAX d= 0.277124605038953 --> position 2098

#What I did
max(d, na.rm = TRUE)
which(d$c.0..0..0..0..0..0..0..0..0..0..0..0..0..0..0..0..0..0..0..0.. == max(d, na.rm = TRUE))
#Max d = 0.2771246 --> position 2098

alpha<-as.data.frame(perf@alpha.values)                            # the alpha values represent the corresponding cut-offs
alpha[2098,]                                                       # to identify the alpha value (i.e. the cut-off) that corresponds to the maximum distance between the 45° line and the curve

# Best cutoff:   0.2516195 (I got the same cutoff as them)
# This value can now be used to build the confusion matrix:

DATA<-matrix(0,4973,3)                                             # to build a matrix with 3 columns and n rows, where n is the dimension of the data set (here 4973 - the number of rows can be checked with dim(dat)) 
DATA<-as.data.frame(DATA)
names(DATA)<-c("plotID","Observed","Predicted")
DATA$plotID<-1:4973                                                # the first column is filled with an ID value that is unique for each row
DATA$Observed<-dat$Pres                                            # the second column reports the observed response (0s and 1s)
DATA$Predicted<-predict(fit11,dat,type="response")                 # the third column reports the predictions
cmx(DATA, threshold = 0.2516195)                                   # the identified cut-off must be used here

#Their results
#             observed
# predicted     1    0
#          1  815 1283
#          0  310 2565

#My results
#             observed
#predicted    1    0
#             1  762 1119
#             0  363 2729

# The confusion matrix can then be transformed into percentages:

#             observed
# predicted     1    0
#          1  72%  33%
#          0  24%  67%

# The area under the curve can also be used as an rough indication of model performance:

auc <- performance(pred, measure="auc")

# auc = 0.7723559 (same value as)

## STEP 7: construct a prediction map ##

pred<-read.table("C:\\Grid_prediction.txt", header=TRUE)            # import a grid of points created in a GIS. Each point should be associated with the values of the covariates that are retained in the final model
pred<-cbind(pred,predict(fit12,pred, type="response"))              # predict animal occurrence using the final model
write.table(pred, "C:\\prediction.txt")                             # save the results in a new text file, that can then be imported in a GIS to draw the prediction map

## STEP 8: visualise the contribution of the explanatory variables by means of the partial residual plots, which plot the estimated relationship between the response (on the link scale) and each predictor ## 

LatForPlotting<- seq(min(dat$Lat), max(dat$Lat), length=50)
LongForPlotting<- seq(min(dat$Long), max(dat$Long), length=50)
DepthForPlotting<- seq(min(dat$Depth), max(dat$Depth), length=50)
AspectForPlotting<- seq(min(dat$Aspect), max(dat$Aspect), length=50)
PredictionData<- expand.grid(LatForPlotting,LongForPlotting,DepthForPlotting,AspectForPlotting)
names(PredictionData)<- c("LatForPlotting","LongForPlotting","DepthForPlotting","AspectForPlotting")
BootstrapParameters<-rmvnorm(10000, coef(fit12), summary(fit12)$cov.unscaled)
test<- glm(Pres ~ bs(Lat,knots=mean(Lat))+bs(Long,knots=mean(Long))+bs(Depth,knots=mean(Depth))+bs(Aspect,knots=mean(Aspect)),family=binomial, data = dat)
x1<-model.matrix(test)[,2:5]%*%coef(fit12)[c(2:5)]
x2<-model.matrix(test)[,6:9]%*%coef(fit12)[c(6:9)]
x3<-model.matrix(test)[,10:13]%*%coef(fit12)[c(10:13)]
x4<-model.matrix(test)[,14:17]%*%coef(fit12)[c(14:17)]

# Partial plot for Latitude

BootstrapCoefs<- BootstrapParameters[,1:5]
Basis<- cbind( rep(1, 10),bs(LatForPlotting, knots=mean(dat$Lat), Boundary.knots=range(dat$Lat)))
RealFit<- Basis%*%coef(fit12)[c(1:5)]
RealFitCenter1<- RealFit-mean(x1)-coef(fit12)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim1<- min(cis-mean(x1)-coef(fit12)[1])
MaximumYlim1<- max(cis-mean(x1)-coef(fit12)[1])
cil1<-cis[1,]-mean(x1)-coef(fit12)[1]
ciu1<-cis[2,]-mean(x1)-coef(fit12)[1]
plot1<-qplot(LatForPlotting,RealFitCenter1,xlab="Latitude", ylab="s(Latitude)", main="a",ylim=c(MinimumYlim1, MaximumYlim1),geom="line")+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil1,ymax=ciu1),stat="identity")+geom_rug(aes(x = dat$Lat, y=-10000))

# Partial plot for Longitude

BootstrapCoefs<- BootstrapParameters[,c(1,6:9)]
Basis<- cbind(rep(1,10), bs(LongForPlotting, knots=mean(dat$Long), Boundary.knots=range(dat$Long)))
RealFit<- Basis%*%coef(fit12)[c(1,6:9)]
RealFitCenter2<- RealFit-mean(x2)-coef(fit12)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim2<- min(cis-mean(x2)-coef(fit12)[1])
MaximumYlim2<- max(cis-mean(x2)-coef(fit12)[1])
cil2<-cis[1,]-mean(x2)-coef(fit12)[1]
ciu2<-cis[2,]-mean(x2)-coef(fit12)[1]
plot2<-qplot(LongForPlotting,RealFitCenter2,xlab="Longitude", ylab="s(Longitude)", main="b",ylim=c(MinimumYlim2, MaximumYlim2),geom="line")+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil2,ymax=ciu2),stat="identity")+geom_rug(aes(x = dat$Long, y=-10000))

# Partial plot for Depth

BootstrapCoefs<- BootstrapParameters[,c(1,10:13)]
Basis<- cbind(rep(1,10), bs(DepthForPlotting, knots=mean(dat$Depth), Boundary.knots=range(dat$Depth)))
RealFit<- Basis%*%coef(fit12)[c(1,10:13)]
RealFitCenter3<- RealFit-mean(x3)-coef(fit12)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim3<- min(cis-mean(x3)-coef(fit12)[1])
MaximumYlim3<- max(cis-mean(x3)-coef(fit12)[1])
cil3<-cis[1,]-mean(x3)-coef(fit12)[1]
ciu3<-cis[2,]-mean(x3)-coef(fit12)[1]
plot3<-qplot(DepthForPlotting,RealFitCenter3,xlab="Depth", ylab="s(Depth)", main="c",ylim=c(MinimumYlim3, MaximumYlim3),geom="line")+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil3,ymax=ciu3),stat="identity")+geom_rug(aes(x = dat$Depth, y=-10000))

# Partial plot for Aspect

BootstrapCoefs<- BootstrapParameters[,c(1,14:17)]
Basis<- cbind(rep(1,10), bs(AspectForPlotting, knots=mean(dat$Aspect), Boundary.knots=range(dat$Aspect)))
RealFit<- Basis%*%coef(fit12)[c(1,14:17)]
RealFitCenter4<- RealFit-mean(x4)-coef(fit12)[1]
BootstrapFits<- Basis%*%t(BootstrapCoefs)
quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
cis<-apply(BootstrapFits, 1, quant.func)
MinimumYlim4<- min(cis-mean(x4)-coef(fit12)[1])
MaximumYlim4<- max(cis-mean(x4)-coef(fit12)[1])
cil4<-cis[1,]-mean(x4)-coef(fit12)[1]
ciu4<-cis[2,]-mean(x4)-coef(fit12)[1]
plot4<-qplot(AspectForPlotting,RealFitCenter4,xlab="Aspect", ylab="s(Aspect)", main="d",ylim=c(MinimumYlim4, MaximumYlim4),geom="line")+theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+geom_smooth(fill="grey",colour="black",aes(ymin=cil4,ymax=ciu4),stat="identity")+geom_rug(aes(x = dat$Aspect, y=-10000))

# Plot the four plots altogether

sidebysideplot <- grid.arrange(plot1,plot2,plot3,plot4,ncol=2,nrow=2)

## STEP 9: analysis by grouping behaviour ##

# The two subsets of the data set can be created as follows:

groups  <- subset(dat,Group!=0) # to select only the encounters with groups
singles <- subset(dat,Group!=1) # to select only the encounters with single males

# All the steps from 1 to 8 should now be repeated using each of the two subsets, in order to assess habitat use by groups and singletons.