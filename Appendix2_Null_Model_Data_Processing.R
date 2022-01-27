#### Null Model Overview ####
## This script produces several types of Null Model animal distributions for comparison to 
#   the real animal tracks

#### Load Libraries and Functions ####
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
pkgTest("tidyverse")
pkgTest("ncdf4")
pkgTest("chron") 
pkgTest("lattice")
pkgTest("RColorBrewer")
pkgTest("oce")
pkgTest("raster")
pkgTest("ggpubr")
pkgTest("lubridate")
pkgTest("rnaturalearth")
pkgTest("maptools")
pkgTest("geosphere")
pkgTest("fossil")
pkgTest("ggmap")
pkgTest("maps")
pkgTest("rnaturalearthdata")
pkgTest("sf")
pkgTest("ggspatial")
pkgTest("marmap")
pkgTest("mapdata")
pkgTest("metR")
pkgTest("rgdal")
pkgTest("stats")
pkgTest("adehabitatHR")

# Download HF Radar data
source('./functions/hfRadar_Download.R')
# Process netCDF file dates and times
source("./functions/ncdate.R")

set.seed(1516)

## Load Pre-processed Data
load("./dataProcessed/datasetFTLE.RData") # Load the processed dive dataframe
# Create a dataframe with no missing FTLE vales and columns for grouping and plotting
data_Gamm_ALL <- datasetFTLE %>% 
  # Only include Dives with FTLE data
  dplyr::filter(is.nan(ftle48)==0, is.na(ftle48)==0) %>% 
  group_by(depid) %>% 
  # Add a time Variable for the AR1 correlation structure 
  mutate(time = startI-min(startI),
         timeHrs = time/3600,
         # Add a binary Foraging field for Logistic Regression. 
         Foraging = if_else(LungeCount > 0,1,0),
         # Add a hourly feeding rate for use in HMM
         total_Duration = dive_Duration + surf_Duration,
         feedingRateHR = LungeCount/total_Duration*60) %>% 
  ungroup()
## Add 1 meter jitter to the 5 duplicate Lat/Long pairs for Spatial Autocorrelation structure
data_Gamm_ALL$unique <- row.names(data_Gamm_ALL)
# look for duplicate Lat/Longs
duplicateLocs <-data_Gamm_ALL[duplicated(data_Gamm_ALL[,c("Lat","Long")]),] 
# add a 1 meter jitter (.00001 degrees of latitude)
data_Gamm_ALL$Lat[data_Gamm_ALL$unique %in% duplicateLocs$unique] <- data_Gamm_ALL$Lat[data_Gamm_ALL$unique %in% duplicateLocs$unique] + 0.00001
data_Gamm_ALL <- data_Gamm_ALL %>% dplyr::select(-unique)
rm(duplicateLocs)

# convert bathymetry to data frame
b = getNOAA.bathy(lon1 = min(data_Gamm_ALL$Long-3,na.rm=TRUE), lon2 = max(data_Gamm_ALL$Long+3,na.rm=TRUE), lat1 =max(data_Gamm_ALL$Lat+3,na.rm=TRUE), lat2 = min(data_Gamm_ALL$Lat-3,na.rm=TRUE),resolution = 1)
bf = fortify.bathy(b)
world <- ne_countries(scale = "large", returnclass = "sf")
FTLEdatapath <- "./dataRaw/FTLE/"
FTLE_Integrations <- c("48")

#### Create Background Points ####
## In order to arrive at a 10:1 ratio of background points, we choose 10 random points for every animal location. However, since 
#  every animal location corresponds to a specific time period for the FTLE data (in hourly increments), we only extract the FTLE 
#  data of the background points for the time-period of their corresponding animal location. 

# Create empty dataframes to hold all the data
datasetBackgroundFTLE <- data.frame()
for(i in 1:length(unique(data_Gamm_ALL$depid))){
  # Get the depid and time zone for this deployment
  depid <- unique(data_Gamm_ALL$depid)[i]
  # find tz_Offset (from external file)
  tzOffset <- dep_metadata %>% filter(ID == depid) %>% .$tzOffset %>% as.character()
  cat(paste0("\nProcessing: ",  unique(data_Gamm_ALL$depid)[i]))
  cat(paste0("\nDeployment ",  i, " of ", length(unique(data_Gamm_ALL$depid))))
  ptm <- proc.time()
  # filter by ID and remove dives without location data
  data_Sub <- data_Gamm_ALL %>% dplyr::filter(depid == unique(data_Gamm_ALL$depid)[i],
                                            !is.na(Lat),!is.na(Long))
  for(ii in 1:length(FTLE_Integrations)){
    ptm1<-proc.time()
    gc()
    cat(paste0("\nFTLE: ",  FTLE_Integrations[ii],"hrs"))
    ncname <- paste0(FTLEdatapath,depid,"_FTLE_",FTLE_Integrations[ii],"Hrs.nc")
    ncin_FTLE <- nc_open(ncname, verbose=FALSE)
    # FTLE Time
    time_FTLE <- ncdate(ncin_FTLE)
    # convert to local time 
    time_FTLE_Local <- time_FTLE
    attr(time_FTLE_Local, "tzone") <- tzOffset # change the timezone to tzOffset
    # Create a FTLE Raster Stack
    FTLE_stack <- stack(ncname, varname = "FTLE")
    # Create a very large dataframe with all divergence values
    ftlestack_df <- as.data.frame(FTLE_stack, xy = TRUE) %>% #View()
      tidyr::gather(.,layer,FTLE,seq(3,length(time_FTLE)+2)) %>% 
      dplyr::filter(!is.na(FTLE)) %>% 
      mutate(depid = unique(data_Sub$depid),
             layer = str_remove(layer,pattern="X"),
             layer = parse_datetime(layer, format= "%Y.%m.%d.%H.%M.%S", locale = locale(tz = "UTC")) )
    for(j in 1:length(time_FTLE)){ 
      # subset by the time of FTLE observation (+- timeDifs)
      tempSP <- data_Sub[(data_Sub$dtStart <= time_FTLE[j]+timeDifS & data_Sub$dtStart > time_FTLE[j]-timeDifS),]
      if(dim(tempSP)[1] != 0 ){ # make sure there is data for this period
        tempFTLE <- ftlestack_df[ftlestack_df$layer == time_FTLE[j],]
        tempFTLE$layerNum <- seq(1,length(tempFTLE$layer),1)
        # Each dive gets 10 random background points
        for(jj in 1:length(tempSP$DiveNum)){
          tempDF <- tempFTLE %>% 
            dplyr::filter(layerNum %in% sample(tempFTLE$layerNum, 10, replace=FALSE)) %>% 
            rename("Long"="x", "Lat"="y","time_FTLE"="layer") %>% 
            mutate(DiveNum = tempSP$DiveNum[jj],
                   Distribution = "Background",
                   Integration = paste0('ftle',FTLE_Integrations[ii])) %>% 
            dplyr::select(depid, DiveNum, Distribution,Integration,everything(), -layerNum)
          # Add to Results DF
          datasetBackgroundFTLE <- rbind(datasetBackgroundFTLE,tempDF)
          rm(tempDF)
        }
        rm(tempFTLE)
      }
      rm(tempSP)
    }
    nc_close(ncin_FTLE)
    rm(ncin_FTLE, ncname,FTLE_stack,time_FTLE,time_FTLE_Local,ftlestack_df)
    tPTM1 <-  proc.time() - ptm1
    cat(paste0("\nElapsed time: ", round(tPTM1[3],1), " seconds"))
    rm(ptm1,tPTM1)   
  }
  rm(depid, tzOffset)  
  tPTM <-  proc.time() - ptm
  cat(paste0("\nOverall elapsed time: ", round(tPTM[3]/60,1), " minutes"))
  rm(ptm, tPTM)
  gc()
}
## Save Background Points
save(datasetBackgroundFTLE, file="./dataProcessed/datasetBackgroundFTLE.RData")


#### Create Simulated Tracks #### 
### Calculate H-value and R-value for all individuals ### 
# R-value:
# Wrapped normal distribution is a wrapped probability distribution that results from the "wrapping" of the normal distribution 
# around the unit circle. It finds application in the theory of Brownian motion and is a solution to the heat equation for periodic 
# boundary conditions. It is closely approximated by the von Mises distribution, which, due to its mathematical simplicity and 
# tractability, is the most commonly used distribution in directional statistics
# Summary Table for Simulations
simmSum <- data_Gamm_ALL %>% 
  dplyr::filter(!is.na(Long)) %>% 
  arrange(depid,dttzStart) %>% 
  group_by(depid) %>% 
  summarize(startLat = first(Lat),
            startLong = first(Long),
            startTime = first(dttzStart),
            endLat = last(Lat),
            endLong = last(Long),
            endTime = last(dttzStart),
            hValue=NA,
            rValue=NA)
for(j in 1:length(simmSum$depid)){
  dataSub <- data_Gamm_ALL[data_Gamm_ALL$depid==simmSum$depid[j],] %>% dplyr::filter(!is.na(Long))
  xy <- coordinates(dataSub[,c("Long","Lat")])
  templTraj <- as.ltraj(xy = xy, date = dataSub$dttzStart,
                        id = dataSub$depid,  typeII = TRUE, proj4string = CRS("+proj=longlat +datum=WGS84"))
  simmSum$hValue[j] <- hbrown(templTraj)   # Calculate H-Value
  tempDF <- ld(templTraj)$abs.angle
  simmSum$rValue[j] <- est.kappa(tempDF[!is.na(tempDF)])/pi
  rm(xy,dataSub,templTraj,tempDF)
}

## Create 10 starting locations sampled randomly from within the MCP of each deployment
simStarts <- data.frame()
for(j in 1:length(simmSum$depid)){
  dataSub <- data_Gamm_ALL[data_Gamm_ALL$depid==simmSum$depid[j],] %>% dplyr::filter(!is.na(Long))
  xy <- coordinates(dataSub[,c("Long","Lat")])
  xySP <- SpatialPointsDataFrame(xy, dataSub, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"), match.ID = FALSE)
  depMCP <- mcp(xySP, percent=100)
  samples <- st_sample( st_as_sf(depMCP), 1000) #%>% st_drop_geometry()
  samplesDF<-do.call(rbind, st_geometry(samples)) %>% 
    as_tibble() %>% setNames(c("Long","Lat"))
  xy <- coordinates(samplesDF[,c("Long","Lat")])
  depth_Values <- as.data.frame(raster::extract(as.raster(b),xy, cellnumbers = TRUE))
  # Remove any sampled start points on land
  samplesDF <- samplesDF[depth_Values$layer<0,] 
  # Take 10 Randomly from the remaining start locations with the MCP
  samplesDF <- samplesDF %>% 
    dplyr::filter(row_number() %in% sample.int(length(samplesDF$Lat),10,replace = FALSE)) %>% 
    mutate(depid = paste0(simmSum$depid[j],"_SimmRS-",row_number()),
           FTLE_ID = simmSum$depid[j]) %>% 
    dplyr::select(depid,Lat,Long, everything())                     
  xy <- coordinates(samplesDF[,c("Long","Lat")])
  xySP <- SpatialPoints(xy, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"))
  
  plot(depMCP)
  plot(samples, add=T)
  plot(xySP,col='red',add=T)
  
  if(j==1){simStarts <- samplesDF}
  else{simStarts <- rbind(simStarts,samplesDF)}
  #cleanup
  rm(xy,xySP,dataSub,samples,samplesDF,depth_Values,depMCP)
}

### Create Simulations ###
simTracks_All <- data.frame()
for(j in 1:length(simmSum$depid)){
  dataSub <- data_Gamm_ALL[data_Gamm_ALL$depid==simmSum$depid[j],] %>% dplyr::filter(!is.na(Long)) %>% arrange(dttzStart)
  xy <- coordinates(dataSub[,c("Long","Lat")])
  templTraj <- as.ltraj(xy = xy, date = dataSub$dttzStart,
                        id = dataSub$depid,  typeII = TRUE, proj4string = CRS("+proj=longlat +datum=WGS84"))
  simStartSub <- simStarts[simStarts$FTLE_ID==simmSum$depid[j],]
  
  simTracks <- data.frame()
  for(k in 1:length(simStartSub$depid)){
    ## Correlated Random Walk - Individual H, Individual R
    depthCheck <- FALSE
    while(!depthCheck){
      # retain the same times as original data, use randomized start locations for CRW
      tempTraj <- simm.crw(date = dataSub$dttzStart,
                          x0 = c(simStartSub$Long[k],simStartSub$Lat[k]),h=hbrown(templTraj), r=simmSum$rValue[j],
                          id= paste0(simmSum$depid[j],"_SimmCRW-",k), proj4string = CRS("+proj=longlat +datum=WGS84"))  
      # Convert to dataframe
      temp_DF_CRW <- ld(tempTraj) %>% 
        rename("Long"= "x", "Lat" = "y", "dttzStart"="date","depid"="id") %>%
        mutate(method = "CRW", 
               FTLE_ID = simmSum$depid[j]) %>% 
        dplyr::select(depid,FTLE_ID,dttzStart,Lat,Long,method) %>% 
        left_join(dplyr::select(dataSub,depid,dttzStart,LungeCount),by=c("FTLE_ID"= "depid","dttzStart"="dttzStart"))
      # Determine if any points go on land
      xy1 <- coordinates(cbind(temp_DF_CRW$Long, temp_DF_CRW$Lat))
      tempdCRWSP <- SpatialPointsDataFrame(xy1, temp_DF_CRW, proj4string = CRS("+proj=longlat"), match.ID = FALSE)
      ### Extract Values ### 
      depth_Values <- as.data.frame(raster::extract(as.raster(b),tempdCRWSP, cellnumbers = TRUE))
      if(max(depth_Values$layer)< 0) depthCheck <- TRUE
      rm(depth_Values,tempdCRWSP,xy1,tempTraj)
    }
    if(k==1){
      simTracks <- temp_DF_CRW
    }else{
      simTracks <- rbind(simTracks, temp_DF_CRW)
    }
    rm(temp_DF_CRW)
  }
  
  pT1<- ggplot()+
    geom_sf(data=world) +
    coord_sf(xlim = c(min(simTracks$Long,na.rm = TRUE)-.25, max(simTracks$Long,na.rm = TRUE)+.25), 
             ylim = c(min(simTracks$Lat,na.rm = TRUE)-.25, max(simTracks$Lat,na.rm = TRUE)+.25), expand = FALSE) +
    geom_path(data=simTracks,aes(x=Long,y=Lat,color=depid,group=depid))+    
    geom_path(data=dataSub,aes(x=Long,y=Lat),color='black')+
    ggtitle(simmSum$depid[j], subtitle = "CRW") + 
    theme(legend.position = "none")
  pT1
  ggsave(paste0(simmSum$depid[j],"-Tracks.png"), width=16, height=10, dpi=300)
  # Save values
  if(k==1){
    simTracks_All <- rbind(simTracks) 
  }else{
    simTracks_All <- rbind(simTracks_All, simTracks)
  }
  rm(simTracks)
}
## Save simTracks_All
save(simTracks_All, file="./dataProcessed/simTracks_All.RData")

#### Create Time-Shifted Tracks ####
## Shift all timestamps by a certain amount to test the temporal persistence of FTLE features

timeShifts <- c(24,48,96,192) # time shift in hours
timeshiftDataset <- data.frame()
for(i in 1:length(timeShifts)){
  tempDF <- data_Gamm_ALL %>% 
    dplyr::filter(!is.na(Long)) %>% 
    mutate(dttzStart = dttzStart+(timeShifts[i]*60*60),
           FTLE_ID = paste0(depid, "-TS")) 
  shft <- if_else(timeShifts[i]>0, paste0("_plus", timeShifts[i],"hrs"),paste0("_minus",abs(timeShifts[i]), "hrs")  )
  tempDF$depid <- paste0(tempDF$depid,shft)
  tempDF$method <- shft
  tempDF <- tempDF %>% dplyr::select(depid,FTLE_ID,dttzStart,Lat,Long,method,LungeCount) 
  if(i==1){
    timeshiftDataset <- tempDF
  }else{
    timeshiftDataset <- rbind(timeshiftDataset,tempDF)
  }         
  rm(tempDF, shft)
}
rm(i)
## Save timeshiftDataset
save(timeshiftDataset, file="./dataProcessed/timeshiftDataset.RData")

# Combine Timeshifted and CRW data for FTLE Processing
datasetSims <- rbind(simTracks_All,timeshiftDataset)

#### Download Timeshift HF Radar Data ####
useOldNC <- FALSE # if true, script will search for a local NC file before downloading a new one
numDays <- 7 # number of days of HF Radar data before and after deployment
degreeBuffer <- 1 # number of degrees around the min/max lat/longs of deployment
ncPattern <- paste0("_HFRADAR_US_West_Coast_", gridSize, "km_Resolution_Hourly.nc")
for(i in 1:length(unique(timeshiftDataset$FTLE_ID))){
  # Get the depid and time zone for this deployment
  depid <- unique(timeshiftDataset$FTLE_ID)[i]
  # find tz_Offset, region and location (from external file)
  tzOffset <- dep_metadata %>% filter(ID == depid) %>% .$tzOffset %>% as.character()
  cat(paste0("\nProcessing: ",  unique(timeshiftDataset$FTLE_ID)[i]))
  cat(paste0("\nDeployment ",  i, " of ", length(unique(timeshiftDataset$FTLE_ID))))
  # filter by ID and remove dives without location data
  data_Sub <- timeshiftDataset %>% filter(FTLE_ID == unique(timeshiftDataset$FTLE_ID)[i],
                                      !is.na(Lat),!is.na(Long)) %>%
    arrange(dttzStart)
  data_Sub$dtStart <- data_Sub$dttzStart
  attr(data_Sub$dtStart, "tzone") <- "GMT" # Set TZ to GMT
  # Determine the start and end times of deployment, +- 1 hour
  start <- first(data_Sub$dtStart)-(3600)
  end <- last(data_Sub$dtStart)+(3600)
  # Only download new NC file if it doesn't exist locally
  if(useOldNC & file.exists(paste0("./dataRaw/RawHFRadarData/",depid,ncPattern))){
    ncname <- paste0("./dataRaw/RawHFRadarData/",depid,ncPattern)
    cat(paste0("Using previously downloaded HF Radar data:\n",ncname))
  }else{ # Download the appropriate NC file using parameters from tag
    ncname <- hfRadar_Download(paste0(unique(timeshiftDataset$FTLE_ID)[i],"-TS"),
                               dStart = lubridate::date(start)-numDays,tStart = paste0(lubridate::hour(start),':00:00'),
                               dEnd = lubridate::date(end)+numDays, tEnd = paste0(lubridate::hour(end),':00:00'),
                               bbox = c(round(max(data_Sub$Lat, na.rm = TRUE)+degreeBuffer,2),round(min(data_Sub$Lat, na.rm = TRUE)-degreeBuffer,2),
                                        round(min(data_Sub$Long, na.rm = TRUE)-degreeBuffer,2),round(max(data_Sub$Long, na.rm = TRUE)+degreeBuffer,2) ),
                               grid = gridSize)
  }
}


#### Process FTLE  ####
## Import the externally processed FTLE data (netCDF) and Extract Values at Whale Locations
## Restoration Settings:
#   Data was restored using automatic land detection and restored up to shoreline. 
#   Concave Hull around available points with 10km alpha
## FTLE Settings: 
#   Include Land and Apply a Free-Slip boundary condition
#   Initial Times Increment: 1 Hour
#   Integration Duration: 2 days
#   Reverse Direction of Time 
## FTLE is calculated with a grid of tracers having 10x the resolution of the input data.
#   The final resolution varies slightly by deployment, but is roughly 590m or 6km HF radar resolution/10.

FTLEdatapath <- "./dataRaw/FTLE/"
FTLE_Integrations <- c("48")
# Create empty dataframes to hold all the data
datasetSimsFTLE <- data.frame()
for(i in 1:length(unique(datasetSims$FTLE_ID))){
  # Get the depid and time zone for this deployment
  depid <- unique(datasetSims$FTLE_ID)[i]
  # find tz_Offset, region and location (from external file)
  tzOffset <- dep_metadata %>% filter(ID == depid) %>% .$tzOffset %>% as.character()
  cat(paste0("\nProcessing: ",  unique(datasetSims$FTLE_ID)[i]))
  cat(paste0("\nDeployment ",  i, " of ", length(unique(datasetSims$FTLE_ID))))
  ptm <- proc.time()
  # filter by ID and remove dives without location data
  data_Sub <- datasetSims %>% dplyr::filter(FTLE_ID == unique(datasetSims$FTLE_ID)[i],
                                          !is.na(Lat),!is.na(Long))
  # Create an empty dataframe for extracted FTLE values
  ftle_ExtractAll <- data.frame()
  for(ii in 1:length(FTLE_Integrations)){
    ptm1<-proc.time()
    gc()
    cat(paste0("\nFTLE: ",  FTLE_Integrations[ii],"hrs"))
    ncname <- paste0(FTLEdatapath,depid,"_FTLE_",FTLE_Integrations[ii],"Hrs.nc")
    ncin_FTLE <- nc_open(ncname, verbose=FALSE)
    # FTLE Time
    time_FTLE <- ncdate(ncin_FTLE)
    # convert to local time 
    time_FTLE_Local <- time_FTLE
    attr(time_FTLE_Local, "tzone") <- tzOffset # change the timezone to tzOffset
    # Create a FTLE Raster Stack and save for later use
    FTLE_stack <- stack(ncname, varname = "FTLE")
    ### Raster Extraction ###
    # Create a spatial points dataframe of whale movement
    xy <- coordinates(cbind(data_Sub$Long, data_Sub$Lat))
    dataSP <- SpatialPointsDataFrame(xy, data_Sub, proj4string = CRS("+proj=longlat +datum=WGS84"), match.ID = FALSE)
    cat("\nExtracting FTLE Data")
    # Create an empty dataframe for extracted FTLE values
    ftle_Extract <- data.frame()
    for(j in 1:length(time_FTLE_Local)){ 
      # subset by the time of FTLE observation (+- timeDifs)
      tempSP <- dataSP[(dataSP$dttzStart <= time_FTLE_Local[j]+timeDifS & dataSP$dttzStart > time_FTLE_Local[j]-timeDifS),]
      if(dim(tempSP)[1] != 0 ){ # make sure there is data for this period
        # extract the FTLE values from the focal layer
        ftle_Values <- raster::extract(FTLE_stack,tempSP,layer = j, nl=1, cellnumbers = TRUE)
        colnames(ftle_Values) <- c('rasterCellFTLE',paste0('ftle',FTLE_Integrations[ii]))          
        # extract the dcell number of each position
        ftle_Cells <- as.data.frame(raster::cellFromXY(FTLE_stack[[j]],tempSP))
        colnames(ftle_Cells) <- "cell"
        ftle_Cells <- cbind(ftle_Cells, raster::rowColFromCell(FTLE_stack[[j]],ftle_Cells$cell))
        ftle_meanLayer <- cellStats(FTLE_stack[[j]],'mean',na.rm=TRUE)
        # combine into a single dataframe
        tempSP <- as.data.frame(tempSP)
        tempSP$extractNumFTLE <- j          
        tempSP <- cbind(tempSP,ftle_Values)
        tempSP$ftle_meanLayer <- ftle_meanLayer
        names(tempSP)[names(tempSP)=='ftle_meanLayer'] <- paste0('ftle_meanLayer',FTLE_Integrations[ii])
        # append it to an output dataframe
        ftle_Extract <- rbind(ftle_Extract, tempSP)
        rm(ftle_Values, ftle_Cells,ftle_meanLayer) # cleanup
      }
      rm(tempSP)
    }
    rm(j)
    
    if(ii==1){
      ftle_ExtractAll <- ftle_Extract 
    } else{
      ftle_ExtractAll <- cbind(ftle_ExtractAll,dplyr::select(ftle_Extract, starts_with("ftle"),-FTLE_ID))
    } 
    nc_close(ncin_FTLE)
    rm(ncin_FTLE, ncname,FTLE_stack,time_FTLE,time_FTLE_Local,ftle_Extract,xy,dataSP)
    tPTM1 <-  proc.time() - ptm1
    cat(paste0("\nElapsed time: ", round(tPTM1[3],1), " seconds"))
    rm(ptm1,tPTM1)
  }
  # Save to an output dataframe 
  if(i==1){
    datasetSimsFTLE <- ftle_ExtractAll 
  } else{
    datasetSimsFTLE <- rbind(datasetSimsFTLE, ftle_ExtractAll)
  } 
  tPTM <-  proc.time() - ptm
  cat(paste0("\nOverall elapsed time: ", round(tPTM[3]/60,1), " minutes"))
  rm(ptm, tPTM)
  gc()
}
# clean up unused variables
rm(i,ii,tzOffset, depid, ftle_ExtractAll, FTLE_Integrations,FTLEdatapath,data_Sub)

# Remove the remnants of the Spatial Points DF conversion
datasetSimsFTLE <- datasetSimsFTLE %>% dplyr::select(-c(coords.x1, coords.x2))   
# Check for INF
sapply(datasetSimsFTLE, function(x) sum(is.infinite(x))>0)

##Save Data for later import
save(datasetSimsFTLE,file="./dataProcessed/datasetSimsFTLE.RData")


