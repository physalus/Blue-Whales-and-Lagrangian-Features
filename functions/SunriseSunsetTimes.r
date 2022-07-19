
# create Spatial Points object of Sakhalvasho coordinates
crds <- cbind(data$long, data$lat)
data$crds <- SpatialPoints(crds,proj4string=CRS("+proj=longlat +datum=WGS84"))

# create Date object to pass to sunriset function
dates <- gsub('-', '/', as.character(data$date))
data$dates <- as.POSIXct(dates,tz=tzOffset)  

#first time need to install package
#install.packages("maptools")
library(maptools)
# calculate sunrise times
data$srise <- sunriset(data$crds, data$dates, direction=c("sunrise"),POSIXct.out=TRUE)[,2]

# Calculate Dawn,Sunrise, Sunset and Dusk times
data$sset <- sunriset(data$crds, data$dates, direction=c("sunset"),POSIXct.out=TRUE)[,2]
data$dawn <- crepuscule(data$crds, data$dates, solarDep=12, direction=c("dawn"),POSIXct.out=TRUE)[,2]
data$dusk <- crepuscule(data$crds, data$dates, solarDep=12, direction=c("dusk"),POSIXct.out=TRUE)[,2]

# Calculate solar position (altitude)
data$solarpos <- solarpos(data$crds, data$dttz)[,2]
data$solarnoon <- solarnoon(data$crds, data$dttz,POSIXct.out=TRUE)[,2]

# determine if point is during night or day
data$night <- ifelse(data$dttz < data$srise | data$dttz > data$sset,'night','day')

# determine if day night dusk or dawn
data$astronomical <- ifelse(hms::as.hms(data$dttz,tz = tzOffset) > hms::as.hms(data$srise,tz = tzOffset) & hms::as.hms(data$dttz,tz = tzOffset) < hms::as.hms(data$sset,tz = tzOffset),'day', 
                            ifelse(hms::as.hms(data$dttz,tz = tzOffset) < hms::as.hms(data$srise,tz = tzOffset) & hms::as.hms(data$dttz,tz = tzOffset) > hms::as.hms(data$dawn,tz = tzOffset),'dawn',
                                   ifelse(hms::as.hms(data$dttz,tz = tzOffset) < hms::as.hms(data$dusk,tz = tzOffset) & hms::as.hms(data$dttz,tz = tzOffset) > hms::as.hms(data$sset,tz = tzOffset) ,'dusk','night')))

drops <- c("dates","crds","season")
data <- data[ , !(names(data) %in% drops)]
rm(crds)
rm(drops)
rm(dates)
