#### Download Netcdf from HFRnet ####
## This function downloads a netCDF file from the hfrnet thredds server:
##     https://hfrnet-tds.ucsd.edu/thredds/catalog.html
## The downloaded fill will contain a subsetted portion of the dataset using the 
## paramaters identified below. The function creates a URL that queries the THREDDS 
## server and downloads the resulting nc file into the working directory.

## Inputs:
  # fname will be the beginning of the name of the saved nc file 
  #    i.e.: fname + "_HFRADAR_US_West_Coast_" + grid + "km_Resolution_Hourly.nc"
  # dStart is the start date in format 2019-10-01
  # tStart is the UTC start time in format 00:00:00
  # dEnd is the end date in format 2019-10-01
  # tEnd is the UTC end time in format 00:00:00
  # bbox is a vector containing the Lat/Longs of the extent to download c(N,S,W,E)
  # grid is the Grid size in KM (Can be either 2km or 6km)

## Returns the filename of the downloaded NC file and the URL that was used to create the file

## To test, run: (this will download a netCDF for Monterey Bay containing data from the last 24 hours)
# hfRadar_Download("test",dStart = Sys.Date()-1, bbox = c(37.25,36.25,-122.5,-121.75),grid = 6)

hfRadar_Download <- function(fname,dStart = '2011-10-01',tStart = '00:00:00',dEnd = Sys.Date(), tEnd = '00:00:00', 
                             bbox = c(49.9920,30.2500,-130.3600,-115.8055), grid = 6) {
  # url <- paste0("http://hfrnet-tds.ucsd.edu/thredds/ncss/HFR/USWC/", grid, 
  #               "km/hourly/RTV/HFRADAR_US_West_Coast_", grid,
  #               "km_Resolution_Hourly_RTV_best.ncd?var=hdop&var=number_of_radials&var=number_of_sites&var=u&var=v&north=",
  #               bbox[1], "&west=", bbox[3], "&east=", bbox[4], "&south=", bbox[2],
  #               "&disableProjSubset=on&horizStride=1&time_start=", dStart, 
  #               "T", tStart, "Z&time_end=", dEnd, "T", tEnd,"Z&timeStride=1&addLatLon=true&accept=netcdf")
  url <- paste0("http://hfrnet-tds.ucsd.edu/thredds/ncss/HFR/USWC/", grid,
                "km/hourly/RTV/HFRADAR_US_West_Coast_", grid,
                "km_Resolution_Hourly_RTV_best.ncd?var=hdop&var=number_of_radials&var=number_of_sites&var=u&var=v&north=",
                #"km_Resolution_Hourly_RTV_best.ncd?var=DOPx&var=DOPy&var=u&var=v&north=",
                bbox[1], "&west=", bbox[3], "&east=", bbox[4], "&south=", bbox[2],
                "&disableProjSubset=on&horizStride=1&time_start=", dStart,
                "T", tStart, "Z&time_end=", dEnd, "T", tEnd,"Z&timeStride=1&addLatLon=true&accept=netcdf" )
  ncname <-  paste0("./", fname,"_HFRADAR_US_West_Coast_", grid, "km_Resolution_Hourly.nc")
  download.file(url,ncname, method = "auto",
                quiet = TRUE, mode="wb", cacheOK = TRUE)
  
  cat(paste0("\nFile downloaded to: ",  ncname))
  cat(paste0("\nURL: ",url))
  return(ncname)
}


