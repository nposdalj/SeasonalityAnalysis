LoadAVISO <- function(envDir){

filenameStatAll = paste(envDir,"AVISOglobalvars.nc",sep="")#load files as data frame
AVISO <<- nc_open(filenameStatAll)
names(AVISO$var)
}

GetSSH <- function(AVISO){
#zos - SSH
v6=AVISO$var[[6]]
SSHvar=ncvar_get(AVISO,v6)
SSH_lon=v6$dim[[1]]$vals
SSH_lat=v6$dim[[2]]$vals
SSH_dates=as.POSIXlt(v6$dim[[3]]$vals*60*60,origin='1950-01-01') #extract the date/time
SSH_dates = as.Date(SSH_dates, format = "%m/%d/%y") #get rid of the time

SSH_ddf <- as.data.frame(SSH_dates[K])
SSH_ddf = SSH_ddf %>% 
  rename(
    time = 'SSH_dates[K]',
  )

#plotting timeseries
I=which(SSH_lon>=min(df1$long) & SSH_lon<= max(df1$long)) #only extract the region we care about
J=which(SSH_lat>=min(df1$lat) & SSH_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(SSH_dates>= startTime & SSH_dates<= endTime) #extract only the dates we care about
SSH2=SSHvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(SSH2)[3] #find the length of time

#take the mean
resSSH=rep(NA,n) 
for (i in 1:n) 
  resSSH[i]=mean(SSH2[,,i],na.rm=TRUE)

SSHdf<<- bind_cols(SSH_ddf,as.data.frame(resSSH))

#plot the time series
plot(1:n,resSSH,axes=FALSE,type='o',pch=20,xlab='',ylab='SSH',las = 3) 
axis(2) 
axis(1,1:n,format(SSH_dates[K]),las = 3) 
box()
}

GetDEN <- function(AVISO){
#mlotst - density ocean mixed layer thickness
v1=AVISO$var[[1]]
DENvar=ncvar_get(AVISO,v1)
DEN_lon=v1$dim[[1]]$vals
DEN_lat=v1$dim[[2]]$vals
DEN_dates=as.POSIXlt(v1$dim[[3]]$vals*60*60,origin='1950-01-01') #extract the date/time
DEN_dates = as.Date(DEN_dates, format = "%m/%d/%y") #get rid of the time

#plotting time series SAPTIN 
I=which(DEN_lon>=min(df1$long) & DEN_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(DEN_lat>=min(df1$lat) & DEN_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K=which(DEN_dates>= startTime & DEN_dates<= endTime) #extract only the dates we care about
DEN2=DENvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(DEN2)[3] #find the length of time

#take the mean
resDEN=rep(NA,n) 
for (i in 1:n) 
  resDEN[i]=mean(DEN2[,,i],na.rm=TRUE) 

DEN_ddf <- as.data.frame(DEN_dates[K])
DEN_ddf = DEN_ddf %>% 
  rename(
    time = 'DEN_dates[K]',
  )

DENdf<<- bind_cols(DEN_ddf,as.data.frame(resDEN))

#plot the time series
plot(1:n,resDEN,axes=FALSE,type='o',pch=20,xlab='',ylab='Density',las = 3) 
axis(2) 
axis(1,1:n,format(DEN_dates[K]),las = 3) 
box()
}

#so - salinity
v2=AVISO$var[[5]]
SALvar=ncvar_get(AVISO,v2)
SAL_lon=v2$dim[[1]]$vals
SAL_lat=v2$dim[[2]]$vals
SAL_dates=as.POSIXlt(v2$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
SAL_dates = as.Date(SAL_dates, format = "%m/%d/%y") #get rid of the time

#thetao - temperature
v3=AVISO$var[[3]]
TEMPvar=ncvar_get(AVISO,v3)
TEMP_lon=v3$dim[[1]]$vals
TEMP_lat=v3$dim[[2]]$vals
TEMP_dates=as.POSIXlt(v3$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
TEMP_dates = as.Date(TEMP_dates, format = "%m/%d/%y") #get rid of the time

GetEKE <- function(AVISO){
#uo - eastward velocity
v4=AVISO$var[[4]]
EASTVvar=ncvar_get(AVISO,v4)
EASTV_lon=v4$dim[[1]]$vals
EASTV_lat=v4$dim[[2]]$vals
EAST_dates=as.POSIXlt(v4$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
EAST_dates = as.Date(EAST_dates, format = "%y/%m/%d") #get rid of the time
EASTdf <- as.data.frame(EASTVvar)

#plotting time series SAPTIN 
I=which(EASTV_lon>=min(df1$long) & EASTV_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(EASTV_lat>=min(df1$lat) & EASTV_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K=which(EAST_dates>= startTime & EAST_dates<= endTime) #extract only the dates we care about
EASTV2=EASTVvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(EASTV2)[3] #find the length of time

#take the mean
resEV=rep(NA,n) 
for (i in 1:n) 
  resEV[i]=mean(EASTV2[,,i],na.rm=TRUE) 

EV_ddf <- as.data.frame(EAST_dates[K])
EV_ddf = EV_ddf %>% 
  rename(
    time = 'EAST_dates[K]',
  )

EVdf<- bind_cols(EV_ddf,as.data.frame(resEV))

#plot the time series
plot(1:n,resEV,axes=FALSE,type='o',pch=20,xlab='',ylab='Eastward Velocity',las = 3) 
axis(2) 
axis(1,1:n,format(EAST_dates[K]),las = 3) 
box()

#vo - northward velocity
v5=AVISO$var[[2]]
NORVvar=ncvar_get(AVISO,v5)
NORV_lon=v5$dim[[1]]$vals
NORV_lat=v5$dim[[2]]$vals
NOR_dates=as.POSIXlt(v5$dim[[4]]$vals*60*60,origin='1950-01-01') #extract the date/time
NOR_dates = as.Date(NOR_dates, format = "%m/%d/%y") #get rid of the time
NORdf <- as.data.frame(NORVvar)

#plotting time series SAPTIN 
I=which(NORV_lon>=min(df1$long) & NORV_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(NORV_lat>=min(df1$lat) & NORV_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K=which(NOR_dates>= startTime & NOR_dates<= endTime) #extract only the dates we care about
NORV2=NORVvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(NORV2)[3] #find the length of time

#take the mean
resNV=rep(NA,n) 
for (i in 1:n) 
  resNV[i]=mean(NORV2[,,i],na.rm=TRUE)

NV_ddf <- as.data.frame(NOR_dates[K])
NV_ddf = NV_ddf %>% 
  rename(
    time = 'NOR_dates[K]',
  )

NVdf<- bind_cols(NV_ddf,as.data.frame(resNV))

#plot the time series
plot(1:n,resNV,axes=FALSE,type='o',pch=20,xlab='',ylab='Northward Velocity',las = 3) 
axis(2) 
axis(1,1:n,format(NOR_dates[K]),las = 3) 
box()



####

#Calculate EKE
u <- (resEV)^2
v <- (resNV)^2
velocity = u + v
EKE_meters = 0.5 * velocity
EKE_cm = EKE_meters * 10000 

EKE <<- bind_cols(SSH_ddf,as.data.frame(EKE_cm))
}


#Salinity
#Plotting in ggplot
r = raster(t(SALvar[,,1]),xmn = min(SAL_lon),xmx = max(SAL_lon),ymn=min(SAL_lat),ymx=max(SAL_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="SAL"
mid = mean(df$SAL)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = SAL)) + 
  ggtitle(paste("Daily Salinity on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting time series SAPTIN 
I=which(SAL_lon>=min(df1$long) & SAL_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(SAL_lat>=min(df1$lat) & SAL_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K=which(SAL_dates>= startTime & SAL_dates<= endTime) #extract only the dates we care about
SAL2=SALvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(SAL2)[3] #find the length of time

#take the mean
resSAL=rep(NA,n) 
for (i in 1:n) 
  resSAL[i]=mean(SAL2[,,i],na.rm=TRUE) 

#plot the time series
plot(1:n,resSAL,axes=FALSE,type='o',pch=20,xlab='',ylab='Salinity',las = 3) 
axis(2) 
axis(1,1:n,format(SAL_dates[K]),las = 3) 
box()

#remove unnecessary variables
rm("SALvar","SAL2", "SAL_lon","SAL_lat")

#Temperature
#Plotting in ggplot
r = raster(t(TEMPvar[,,1]),xmn = min(TEMP_lon),xmx = max(TEMP_lon),ymn=min(TEMP_lat),ymx=max(TEMP_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="TEMP"
mid = mean(df$TEMP)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = TEMP)) + 
  ggtitle(paste("Daily Temperature on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")

#plotting time series SAPTIN 
I=which(TEMP_lon>=min(df1$long) & TEMP_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(TEMP_lat>=min(df1$lat) & TEMP_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K=which(TEMP_dates>= startTime & TEMP_dates<= endTime) #extract only the dates we care about
TEMP2=TEMPvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(TEMP2)[3] #find the length of time

#take the mean
resTEMP=rep(NA,n) 
for (i in 1:n) 
  resTEMP[i]=mean(TEMP2[,,i],na.rm=TRUE)

#plot the time series
plot(1:n,resTEMP,axes=FALSE,type='o',pch=20,xlab='',ylab='Temperature',las = 3) 
axis(2) 
axis(1,1:n,format(TEMP_dates[K]),las = 3) 
box()


#converting _dates to data frames and renaming column to 'time'


SAL_ddf <- as.data.frame(SAL_dates[K])
SAL_ddf = SAL_ddf %>% 
  rename(
    time = 'SAL_dates[K]',
  )
TEMP_ddf <- as.data.frame(TEMP_dates[K])
TEMP_ddf = TEMP_ddf %>% 
  rename(
    time = 'TEMP_dates[K]',
  )

#merge res dataframes with dates

SALdf<- bind_cols(SAL_ddf,as.data.frame(resSAL))
TEMPdf<- bind_cols(TEMP_ddf,as.data.frame(resTEMP))
EVdf<- bind_cols(EV_ddf,as.data.frame(resEV))

