library(marmap)
library(ggplot2)
library(oce)
library(ocedata)
library(mapproj)
library(maps)
library(mapdata)
library("ggplot2")
library("sf")
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(GISTools)
library(prettymapr)
data(coastlineWorldFine)
data(topoWorld)

Inlet = getNOAA.bathy(lon1 = -90, lon2 = -70, lat1 = 65, lat2 = 80, resolution = 1)
summary(Inlet)

#make HARP points
HARP = data.frame(lon = c(-76.2), lat = c(72.7))

#Option 1

#convert bathymetry to data frame
bf = fortify.bathy(Inlet)

# get regional polygons
reg = map_data("world2")
reg = subset(reg, region %in% c('Canada', 'USA'))

# convert lat longs
reg$long = (360 - reg$long)*-1

#set map limits
lons = c(-81, -74)
lats = c(71.5, 73.5)

#make a plot
ggplot()+
  
  #add 100m countour
  geom_contour(data = bf,
               aes(x = x, y = y, z = z),
               breaks = c(100),
               size = c(0.3),
               colours = "grey")+
  
  # add 250m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-250),
               size=c(0.6),
               colour="grey")+
  
  #add coastline
  geom_polygon(data = reg, aes(x = long, y = lat, group = group), 
               fill = "darkgrey", color = NA) + 
  
  
  # add points
  geom_point(data = HARP, aes(x = lon, y = lat),
             colour = "black", fill = "grey", 
             stroke = .5, size = 2, 
             alpha = 1, shape = 21)+
  
  # configure projection and plot domain
  coord_map(xlim = lons, ylim = lats)+
  
  # formatting
  ylab("")+xlab("")+
  theme_bw()


#Option 2

plot(Inlet, image = TRUE, land = TRUE, axes = FALSE, lwd = 0.1,
     bpal = list(c(0, max(Inlet), grey(.7), grey(.9), grey(.95)), 
                 c(min(Inlet), 0, "darkblue","lightblue")))
plot(Inlet, n = 1, lwd = 0.5, add = TRUE)

scaleBathy(Inlet, deg = 2, x = "bottomleft",inset = 5)

blues = colorRampPalette(c("red","purple","blue",
                            "cadetblue1","white"))
plot(Inlet, image = TRUE, bpal = blues(100))


# Creating a custom palette of blues
blues <- c("lightsteelblue4", "lightsteelblue3",
           "lightsteelblue2", "lightsteelblue1")
# Plotting the bathymetry with different colors for land and sea
plot(Inlet, image = TRUE, land = TRUE, lwd = 0.1,
     bpal = list(c(0, max(Inlet), "grey"),
                 c(min(Inlet),0,blues)))
# Making the coastline more visible
plot(Inlet, deep = 0, shallow = 0, step = 0,
     lwd = 0.4, add = TRUE)

#Option 3
map("worldHires","Canada", xlim=c(-141,-53), ylim=c(40,85), col="gray90", fill=TRUE)

#Option 4
mp <- function() {
  mapPlot(coastlineWorldFine, projection="+proj=merc",
          longitudelim = c(-80, -74),
          latitudelim = c(72, 73), col='grey', drawBox =TRUE)
}
mp()

b = as.topo(Inlet)

mp()
mapImage(b, col=oceColorsGebco, breaks=seq(-4000, 0, 500))
mapPolygon(coastlineWorldFine, col='grey')
mapPoints(-76.2, 72.7,  col = "blue", bg = "blue", lwd = 2) #HARP
mapText(-76.2, 72.7, "Pond Inlet", pos = 4, col = "blue")
mapPoints(-76.55, 72.68,  col = "red", bg = "red", lwd = 2) #Other recorder
mapText(-76.55, 72.68, "Guys Bight", pos = 4, col = "red", offset = -4.5)
mapScalebar(x = "topleft")
addnortharrow(pos = "bottomright")
mapGrid()

#Larger Plot using Option #4
par(mar=c(1.5, 1.5, 0.5, 0.5))
mapPlot(coastlineWorldFine,
        longitudelim=c(-130,-70), latitudelim=c(40, 80),
        projection="+proj=lcc +lat_0=30 +lat_1=60 +lon_0=-100", col='gray')

theme_set(theme_bw())
world = ne_countries(scale = "medium", returnclass = "sf")
crs_longlat <- "+proj=longlat +ellps=GRS80 +no_defs"
ggplot(data = world) +
  geom_sf() +
  coord_sf(expand = FALSE, crs = crs_longlat)
  coord_sf(crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ")



