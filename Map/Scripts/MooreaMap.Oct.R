###############################################################################################
##Making map of moorea for October 2019
##created by Danielle Becker 03/02/20
##edited by Danielle Becker 04/02/20

#clear dataframe#####
rm(list=ls())

##Install packages
# load packages
library(tidyverse)
library(lubridate)
library(data.table)
library(ggrepel)
library(ggpubr)
library(ggsn)
library(magick)
library(maps)
library(maptools)
library(ggplot2)
library(grid)
library(ggimage)
library(ggmap)
library(here)


#set wd
here()

#zoomed out map version
#make map and set sources
myMap <- get_googlemap(center = c(-149.84, -17.54), zoom = 12, maptype = 'satellite' ) 

#Get the bounding box for the map and reformat it into a data.frame for scalebar
bb <- attr(myMap, "bb")
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))


#Add the scalebar to a ggmap, need ggsn package
map <- ggmap(myMap) + 
  labs(x = 'Longitude', y = 'Latitude') +
  scalebar(data = bb2, dist = 2, st.color = "white",transform = TRUE, model  = "GRS67",  dist_unit = "km", #data are bounding box of map, model gives dataum from google maps, transform recognizes as GPS points as decimal units, location sets location on map, anchor sets bounds of locatio non map
           location = "topright", st.dist = 0.036, box.fill = "white", #sets scalebar bottom left, st.dist is distance bewteen box and numbers,box.fill fills box
           anchor = c( x = bb$ll.lon + 0.05, y = bb$ll.lat +0.03)) + #sets scalebar in specific position long x lat
  geom_point(aes(x = -149.817167, y = -17.473000),
             alpha = .5, color="green3", size = 3, shape = 20)  #geom_point adds specific plot points to map


#add north symbol to map
north2 (map, x=.27, y=.29, symbol=4) 



#zoomed in map version
#make map and set sources
myMap <- get_googlemap(center = c(-149.81, -17.48), maptype = "satellite", zoom = 15) 

#Get the bounding box for the map and reformat it into a data.frame for scalebar
bb <- attr(myMap, "bb")
bb2 <- data.frame(long = unlist(bb[c(2, 4)]), lat = unlist(bb[c(1,3)]))


#Add the scalebar to a ggmap, need ggsn package
map <- ggmap(myMap) + 
  labs(x = 'Longitude', y = 'Latitude') +
  scalebar(data = bb2, dist = 300, st.color = "white",transform = TRUE, model  = "GRS67",  dist_unit = "m", #data are bounding box of map, model gives dataum from google maps, transform recognizes as GPS points as decimal units, location sets location on map, anchor sets bounds of locatio non map
           location = "topright", st.dist = 0.036, box.fill = "white", #sets scalebar bottom left, st.dist is distance bewteen box and numbers,box.fill fills box
           anchor = c( x = bb$ll.lon + 0.025, y = bb$ll.lat +0.022)) + #sets scalebar in specific position long x lat
  geom_point(aes(x = -149.817167, y = -17.473000),
             alpha = .5, color="yellow", size = 3, shape = 19) + #geom_point adds specific plot points to map
  geom_point(aes(x = -149.8178, y = -17.473100),
             alpha = .5, color="green4", size = 3, shape = 19) 
 
#add north symbol to map
north2 (map, x=.784, y=.91, symbol=4) 


png(filename = "Output/MooreaMap.Oct.zoom.png")
plot(map)
dev.off














