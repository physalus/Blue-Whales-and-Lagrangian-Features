#first time need to install package
#install.packages("maptools")
library(maptools)
library(lubridate)


## Function to find solar data based on lat/long and time. 
# Note: this uses the Astronomical definition of Sunrise/Sunset (i.e. -18 degrees below horizon)
# Input: 
#   longitude (numeric, decimal degrees) - Can be a single value or a vector of length(datetime_Local)
#   latitude (numeric, decimal degrees) - Can be a single value or a vector of length(datetime_Local)
#   dateTime_Local (posixCT, Local Timezone) - vector of datetimes in local timezone
# Output:
#   Dataframe containing:
# 
#   solarpos - sun angle
#   solarnoon - posixct time of solar noon
#   srise, sset,dawn, dusk - posixct times of these events
#   night - (day or night) - Crepuscular is included in day
#   astronomical - (day, night, dusk, dawn)
find_Astronomical <- function(longitude, latitude, datetime_Local) {
  # create a check if length datetime_local is the same as long and lat
  # if long lat are 1, repeat that value, else use full length
  if(length(latitude) == length(longitude) & length(latitude) == length(datetime_Local)){
    data <- data.frame(dttz = datetime_Local, long = longitude, lat = latitude)  
  } else{
    data <- data.frame(dttz = datetime_Local, 
                       long = rep(longitude, length(datetime_Local)), 
                       lat = rep(latitude, length(datetime_Local)))
  }
  # create Spatial Points object of Sakhalvasho coordinates
  crds <- cbind(data$long, data$lat)
  data$crds <- SpatialPoints(crds,proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # create Date object to pass to sunriset function
  dates <- gsub('-', '/', as.character(lubridate::date(data$dttz)))
  data$dates <- as.POSIXct(dates,tz=attr(data$dttz, "tzone") )  
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
  # Remove coordinates
  drops <- c("dates","crds")
  data <- data[ , !(names(data) %in% drops)]
  # Determine if point is during night or day
  # Note: this includes crepuscular with daytime
  data$night <- ifelse(hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) < hms::as.hms(data$dawn, tz=attr(data$dttz, "tzone")) | hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) > hms::as.hms(data$dusk, tz=attr(data$dttz, "tzone")),'night','day')
  # Determine if day night dusk or dawn
  data$astronomical <- ifelse(hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) > hms::as.hms(data$srise, tz=attr(data$dttz, "tzone")) & hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) < hms::as.hms(data$sset, tz=attr(data$dttz, "tzone")),'day', 
                              ifelse(hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) < hms::as.hms(data$srise, tz=attr(data$dttz, "tzone")) & hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) > hms::as.hms(data$dawn, tz=attr(data$dttz, "tzone")),'dawn',
                                     ifelse(hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) < hms::as.hms(data$dusk, tz=attr(data$dttz, "tzone")) & hms::as.hms(data$dttz, tz=attr(data$dttz, "tzone")) > hms::as.hms(data$sset, tz=attr(data$dttz, "tzone")) ,'dusk','night')))
  return(data)
}
