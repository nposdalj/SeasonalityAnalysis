GetSST <- function(envDir){
  
startTime = as.Date(startTime) #this should be formatted like this: 2010-03-05
endTime = as.Date(endTime) 

#SST data
#spatial polygon for area of interest
ch <- chull(df1$long, df1$lat)
coords <- df1[c(ch, ch[1]), ]#creating convex hull
sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID = 1)))#converting convex hull to spatial polygon

filenameStatAll = paste(envDir,"AQUASST_SAPTIN.csv",sep="")#load files as data frame
SST = read.csv(filenameStatAll)
SST = SST[-1,] #delete first row
SST$latitude = as.numeric(SST$latitude)
SST$longitude = as.numeric(SST$longitude)
SST = SST[complete.cases(SST[ , 2:3]),]#remove any rows with lat or long as na

#subset the dataframe based on the area of interest
coordinates(SST) <- c("latitude", "longitude")
coords <- over(SST, sp_poly)
SST2 <- SST[coords == 1 & !is.na(coords),]
SST2$sstMasked = as.numeric(SST2$sstMasked)#converting SST from character to numeric
SST2$time = as.Date(SST2$time)#converting time from character to date
SST3 = as.data.frame(SST2)#converting SPDF back to DF

#average the environmental variable based on the ITS over the area of interest
library(dplyr)
SST4 = SST3 %>%
  mutate(time = floor_date(time, "day")) %>%
  group_by(time) %>%
  summarize(mean_SST = mean(sstMasked, na.rm = TRUE), SD_SST = sd(sstMasked, na.rm = TRUE)) #finding daily mean

SST4 <<- SST4

#data exploration
mean(SST2$sstMasked, na.rm = TRUE)#finding overall mean
#save standard deviation
sd(SST2$sstMasked, na.rm = TRUE)
#SST histogram
hist(SST2$sstMasked)

#plot time series
plot(SST4$time, SST4$mean_SST)#exploratory plot
title1 = paste(site,"Sea Surface Temperature Plot")
ggplot(SST4, aes(x=time,y=mean_SST))+
  ggtitle(title1)+
  labs(y="Mean SST (C)",x="Time (days)")+
  geom_line()+
  geom_point()
}
