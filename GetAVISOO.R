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
SSH_dates <<- as.Date(SSH_dates, format = "%m/%d/%y") #get rid of the time

#plotting timeseries
I=which(SSH_lon>=min(df1$long) & SSH_lon<= max(df1$long)) #only extract the region we care about
if (length(I) == 1){ #if the longitude only has 1 value, add a second
  II = I:(I+1)
}else{
  II = I
}
J=which(SSH_lat>=min(df1$lat) & SSH_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(I+1)
}else{
  JJ = J
}
K=which(SSH_dates>= startTime & SSH_dates<= endTime) #extract only the dates we care about
SSH2=SSHvar[II,JJ,K] #index the original data frame to extract the lat, long, dates we care about

n=dim(SSH2)[3] #find the length of time

#take the mean
resSSH=rep(NA,n) 
for (i in 1:n) 
  resSSH[i]=mean(SSH2[,,i],na.rm=TRUE)

SSH_ddf <- as.data.frame(SSH_dates[K])
SSH_ddf <- SSH_ddf %>% 
  rename(
    time = 'SSH_dates[K]',
  )

SSHdf<<- bind_cols(SSH_ddf,as.data.frame(resSSH))

#plot the time series
plot(1:n,resSSH,axes=FALSE,type='o',pch=20,xlab='',ylab='SSH',las = 3) 
axis(2) 
axis(1,1:n,format(SSH_dates[K]),las = 3) 
box()
}

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
J=which(EASTV_lat>=min(df1$lat) & EASTV_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(EAST_dates>= startTime & EAST_dates<= endTime) #extract only the dates we care about
EASTV2=EASTVvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

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
J=which(NORV_lat>=min(df1$lat) & NORV_lat<=max(df1$lat)) #only extract the region we care about
if (length(J) == 1){ #if the latitude only has 1 value, add a second
  JJ = J:(J+1)
}
K=which(NOR_dates>= startTime & NOR_dates<= endTime) #extract only the dates we care about
NORV2=NORVvar[I,JJ,K] #index the original data frame to extract the lat, long, dates we care about

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

EKE <<- bind_cols(NV_ddf,as.data.frame(EKE_cm))
}