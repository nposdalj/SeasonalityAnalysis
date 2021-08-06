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