#plotting in ggplot
r = raster(t(ChlAvar[,,1]),xmn = min(ChlA_lon),xmx = max(ChlA_lon),ymn=min(ChlA_lat),ymx=max(ChlA_lat))
points = rasterToPoints(r, spatial = TRUE)
df = data.frame(points)
names(df)[names(df)=="layer"]="Chl"
mid = mean(df$Chl)
ggplot(data=world) +  geom_sf()+coord_sf(xlim= c(min(df1$long),max(df1$long)),ylim= c(min(df1$lat),max(df1$lat)),expand=FALSE)+
  geom_raster(data = df , aes(x = x, y = y, fill = Chl)) + 
  ggtitle(paste("Daily Chl on", dates[1]))+geom_point(x = 145.46, y = 15.3186, color = "black",size=3)+
  xlab("Latitude")+ylab("Longitude")+
  scale_fill_gradient2(midpoint = mid, low="yellow", mid = "orange",high="red")