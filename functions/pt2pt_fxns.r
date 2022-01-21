pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
pkgTest("fossil")

## A function to calculate the westing/easting between consecutive points ##
## Output in meters ##
pt2pt.easting <- function(longitude.m){
	require(fossil)
	easting <- c(long.m[2:length(long.m)] - long.m[1:(length(long.m)-1)],"NA")
	return(easting)
}

## A function to calculate the distances between consecutive points ##
## Output in meters ##
pt2pt.distance <- function(latitude, longitude){
	require(fossil)
	distance <- c(deg.dist(lat1=latitude[1:(length(latitude)-1)], long1=longitude[1:(length(longitude)-1)], lat2=latitude[2:length(latitude)], long2=longitude[2:length(longitude)]),NA)*1000
	return(distance)
}

## A function to calculate the backward distances between consecutive points ##
## Output in meters #
pt2pt.back.distance <- function(latitude, longitude){
	require(fossil)
	distance <- c(NA,deg.dist(lat1=latitude[1:(length(latitude)-1)], long1=longitude[1:(length(longitude)-1)], lat2=latitude[2:length(latitude)], long2=longitude[2:length(longitude)]))*1000
	return(distance)
}

## A function to calculate the time increment between consecutive points ##
## Default output in seconds, other options are "auto", "mins", "hours","days", or "weeks" ##
pt2pt.duration <- function(datetime, output.units='secs'){
	duration <- c(difftime(datetime[2:length(datetime)], datetime[1:(length(datetime)-1)], units=output.units), 'NA')
	return(duration)
}

## A function to calculate the backward time increment between consecutive points ##
## Default output in seconds, other options are "auto", "mins", "hours","days", or "weeks" ##
pt2pt.back.duration <- function(datetime, output.units='secs'){
	duration <- c(NA,abs(difftime(datetime[1:(length(datetime)-1)], datetime[2:length(datetime)], units=output.units)))
	return(duration)
}

## A function to calculate the movement direction beteen consecutive points ##
## Output in degrees ##
pt2pt.direction <- function(latitude,longitude){
	require(fossil)
	## First calculate direction ##
	direction <- c(earth.bear(long1=longitude[1:(length(longitude)-1)], lat1=latitude[1:(length(latitude)-1)], long2=longitude[2:length(longitude)], lat2=latitude[2:length(latitude)]),NA)
	## Adjust direction (on scale from 0 to 360) to GPS scale (on scale from -180 to 180) ##
	direction <- ifelse(direction > 180, direction-360, direction)
	return(direction)
}

## A function to calculate the backward movement direction beteen consecutive points ##
## Output in degrees ##
pt2pt.back.direction <- function(latitude,longitude){
  require(fossil)
  ## First calculate direction ##
  direction <- c(NA,earth.bear(long1=longitude[1:(length(longitude)-1)], lat1=latitude[1:(length(latitude)-1)], long2=longitude[2:length(longitude)], lat2=latitude[2:length(latitude)]))
  ## Adjust direction (on scale from 0 to 360) to GPS scale (on scale from -180 to 180) ##
  direction <- ifelse(direction > 180, direction-360, direction)
  return(direction)
}

## A function to calculate the speed of movement between consecutive points ##
pt2pt.speed <- function(distance, duration){
	return(distance/duration)
}

## A function to calculate 1) the distance from a point to a static point and 2) a flag indicating whether the distance is within a threshold ##
pt2pt.range <- function(latitude, longitude, ptofinterest.lat, ptofinterest.lon, threshold=1){
	require(fossil)
	range.distance <- deg.dist(lat1=latitude, long1=longitude, lat2=ptofinterest.lat, long2=ptofinterest.lon)*1000
	is.in.range <- ifelse(range.distance < threshold, 1, 0)

return(data.frame(range.distance, is.in.range))
}

## A function to calculate a new lat and long coordinate from an original lat and long a bearing and a distance ##
mk.new.lat.lon <- function (lat, long, bearing, distance) 
{
    rad <- pi/180
    a1 <- lat * rad
    a2 <- long * rad
    tc <- bearing * rad
    d <- distance/6378.1 # distance / radius of earth in KM
    nlat <- asin(sin(a1) * cos(d) + cos(a1) * sin(d) * cos(tc))
    dlon <- atan2(sin(tc) * sin(d) * cos(a1), cos(d) - sin(a1) * 
        sin(nlat))
    nlon <- ((a2 + dlon + pi)%%(2 * pi)) - pi
    npts <- data.frame(nlat/rad, nlon/rad)
    return(npts)
}
