
# Parses the DateTime from a netCDF File 
# Default timezone is UTC unless provided
# Returns:
#   A vector of posixCt values
ncdate <- function(nc, tz = 'GMT') {
  ncdims <- names(nc$dim) #Extract dimension names
  timevar <- ncdims[which(ncdims %in% c("time", "Time", "datetime", "Datetime",
                                        "date", "Date"))[1]] # Pick the time dimension
  ntstep <-nc$dim[[timevar]]$len
  t <- ncvar_get(nc, timevar) # Extract the timestep count
  tunits <- ncatt_get(nc, timevar, "units") # Extract the long name of units
  tspace <- t[2] - t[1] # Calculate time period between two timesteps, for the "by" argument 
  tstr <- strsplit(tunits$value, " ") # Extract string components of the time unit
  a<-unlist(tstr[1]) # Isolate the unit .i.e. seconds, hours, days etc.
  uname <- a[which(a %in% c("seconds","hours","days"))[1]] # Check unit
  startd <- as.POSIXct(gsub(paste(uname,'since '),'',tunits$value),format="%Y-%m-%d %H:%M:%S", tz=tz) ## Extract the start / origin date
  tmulti <- 3600 # Declare hourly multiplier for date
  if (uname == "days") tmulti =86400 # Declare daily multiplier for date
  ## Rename "seconds" to "secs" for "by" argument and change the multiplier.
  if (uname == "seconds") {
    uname <- "secs"
    tmulti <- 1 }
  byt <- paste(tspace,uname) # Define the "by" argument
  if (byt == "0.0416666679084301 days") { ## If the unit is "days" but the "by" interval is in hours
    byt= "1 hour"                       ## R won't understand "by < 1" so change by and unit to hour.
    uname = "hours"}
  datev <- seq(from=as.POSIXct(startd+t[1]*tmulti),by= byt, units=uname,length=ntstep)
}

