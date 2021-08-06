GetChla <- function(envDir){

#loading as .ncfile
filenameStatAll = paste(envDir,"Chl2.nc",sep="")#load files as data frame
ChlA = nc_open(filenameStatAll)
v1=ChlA$var[[1]]
ChlAvar=ncvar_get(ChlA,v1)
ChlA_lon=v1$dim[[1]]$vals
ChlA_lat=v1$dim[[2]]$vals
ChlA_dates= as.POSIXlt(v1$dim[[3]]$vals,origin='1970-01-01',tz='GMT')
ChlA_dates <<- as.Date(ChlA_dates)

#plotting timeseries
I=which(ChlA_lon>=min(df1$long) & ChlA_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(ChlA_lat>=min(df1$lat) & ChlA_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K<<-which(ChlA_dates>=startTime & ChlA_dates<=endTime) #extract only the dates we care about

ChlA2=ChlAvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(ChlA2)[3] #find the length of time

#take the mean
resChlA=rep(NA,n) 
for (i in 1:n) 
  resChlA[i] = mean(ChlA2[,,i],na.rm=TRUE) 

resChlA <<- resChlA

ChlA_ddf <- as.data.frame(ChlA_dates[K])
ChlA_ddf = ChlA_ddf %>% 
  dplyr::rename(
    time = 'ChlA_dates[K]',
  )
ChlA_ddf$time=as.Date(ChlA_ddf$time)

ChlAdf <<- bind_cols(ChlA_ddf, as.data.frame(resChlA))

#plot the time series
plot(1:n,resChlA,axes=FALSE,type='o',pch=20,xlab='',ylab='ChlA',las = 3)
axis(2)
axis(1,1:n,format(resChlA[K]),las = 3)
box()
}
