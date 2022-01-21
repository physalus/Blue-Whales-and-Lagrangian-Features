#### Surface Current Tracking Overview ####
# This script 
# 1) Loads 1Hz dive data exported from Matlab, incorporates 
#    foraging and movement data (by dive)
# 2) Downloads and processes HF Radar data concurrent to the whale data in space 
#    and time
# 3) Imports FTLE Data (processed externally from the data downloaded in step 2)
#    and Extracts FTLE data for all whale positions

# James Fahlbusch
# Goldbogen Lab, Stanford University
# Version 20220111

## Hypothesis:
#   Surface current features indicative of aggregative processes will improve predictions of both where 
#    and when foraging will occur

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
  pkgTest("raster")  
  pkgTest("lubridate")  
  pkgTest("adehabitatHR")
  pkgTest("move")  
  pkgTest("R.matlab")
  pkgTest("rstudioapi")
  
  # Helper functions specific to this analysis
  `%notin%` <- Negate(`%in%`)
  # Process netCDF file dates and times
  source("./functions/ncdate.R")
  # Finds Point to Point metrics
  source("./functions/pt2pt_fxns.R")
  # Finds Vertical Velocity
  source("./functions/tagtools_fxns.R")
  # pkgTest("metR") # note: getNOAA.bathy was not working, but an alternate version on github is used here
  source('./functions/getNOAAbathy.R')
  
#### Global Variables ####
  
  # Find locations withinTime minutes (15) of the start of the dive
  # Selects the nearest to the start of a dive
  withinTime <- 15 # number of minutes +- to search for locations
  minDiveDepth <- 10 # depth in meters to be considered a dive
 
  # HF Radar Grid Size (2km or 6km)
  gridSize <- 6
  # HF Radar extraction - number of seconds to center around timeP (+-) 
  timeDifS <- 30*60 
  # Load Metadata for each deployment
  # Note: this file makes assumptions of whaleNames
  dep_metadata <- read.csv(paste0(getwd(),"./dataRaw/DeploymentMetadata.csv") , sep=",",header=TRUE)

#### Process Whale Data ####

  ### Import Dive Data (prhSC.mat file)###
  #filename <- file.choose() 
  filename <- "./dataRaw/Bm/Bm160523-A20 1HzprhSC.mat"
  species <- substr(basename(filename),1,2) # extract the species
  pathChoice <- dirname(filename)
  # This imports all GPS.csv files in the directory chosen
  filenames <- list.files(path = pathChoice, pattern = "* 1HzprhSC.mat", all.files = FALSE, full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  # Create an empty dataframe to hold all the 1 hz data for all animals
  dataset1hz = data.frame() 
  # Create an empty dataframe to hold all the dive data for all animals
  diveDataset = data.frame() 
  # Create an empty dataframe to hold all the location data for all animals
  locDataset = data.frame() 
  # Create an empty dataframe to hold all the location data for all animals
  lungeDataset = data.frame() 
  # Iteratively import files from each deployment and create columns for sorting and processing.
  for(i in 1:length(filenames)) {
    # determine deployment ID from filename (used throughout)
    depid <- str_split(basename(filenames[i])," 1HzprhSC.mat")[[1]][1] 
    cat(paste0("\nProcessing: ",  depid))
    # find tz_Offset, region and animalID (from external file)
    tzOffset <- dep_metadata %>% filter(ID == depid) %>% .$tzOffset %>% as.character()
    # Deployment Times for filtering (Note: These are in Local Time)
    tagon <- dep_metadata %>% filter(ID == depid) %>% .$tagon %>% as.character()
    tagon <- as.POSIXct(strptime(tagon,format="%m/%d/%Y %H:%M"), tz = tzOffset)
    tagoff <- dep_metadata %>% filter(ID == depid) %>% .$tagoff %>% as.character()  
    tagoff <- as.POSIXct(strptime(tagoff,format="%m/%d/%Y %H:%M"), tz = tzOffset)
    cat(paste0("\nTagon: ",  tagon))
    cat(paste0("\nTagoff: ",  tagoff))    
    
    ## Import 1Hz Dive Data ## 
    # 1Hz prh exported from matlab
    prh1hzData <- readMat(filenames[i]) 
    prh1hzData_DF <- data.frame(p=prh1hzData$p,DN=prh1hzData$DN)
                                # ,LL=prh1hzData$Light,Temp=prh1hzData$T)
    prh1hzData_DF$depid <- depid # set the deployment ID
    prh1hzData_DF$tzOffset <- tzOffset
    prh1hzData_DF$region <- dep_metadata %>% filter(ID == depid) %>% .$Region %>% as.character()
    prh1hzData_DF$animalID <- dep_metadata %>% filter(ID == depid) %>% .$Animal_ID %>% as.character()
    # Create proper datetime fields (dttz is the datetime in local time)
    ## print fractional seconds
    options(digits.secs=2)
    prh1hzData_DF$dttz <- round_date(as.POSIXct((prh1hzData_DF$DN - 719529)*86400, origin = "1970-01-01", tz = 'GMT'),unit="second")
    prh1hzData_DF$dttz <- force_tz(prh1hzData_DF$dttz, tzOffset)
    # make a UTC Datetime column
    prh1hzData_DF$dt <- prh1hzData_DF$dttz
    attr(prh1hzData_DF$dt, "tzone") <- "GMT" # Set TZ to GMT
    # Check to see if a corrected depth was included in the prhfile
    if(exists("Depth",prh1hzData)){
      prh1hzData_DF$p_Orig <- prh1hzData_DF$p
      pnew <- data.frame(p=prh1hzData$Depth)
      prh1hzData_DF$p <- round(pnew$p,1)
      rm(pnew)
    }
    # Reorder the columns
    prh1hzData_DF <- prh1hzData_DF %>% dplyr::select(depid, region, animalID,dttz,dt, everything(),-DN)
    
    ## Dives and Surfacings (dives >10m) ##
    # find_dives in tagtools returns 4 values per dive
    # startI, endI, maxDepth, maxDepthI
    dives <- do.call(cbind.data.frame, prh1hzData$dives)
    colnames(dives) <- c('startI', 'endI', 'maxDepth', 'maxDepthI')
    #determine datetime of dives
    dives$dtStart <- prh1hzData_DF$dt[dives$startI]
    dives$dtEnd <- prh1hzData_DF$dt[dives$endI]
    dives$dttzStart <- prh1hzData_DF$dttz[dives$startI]
    dives$dttzEnd <- prh1hzData_DF$dttz[dives$endI]
    dives$DiveNum <- seq.int(nrow(dives))
    dives <- dives %>% 
      # dive duration in min
      mutate(dive_Duration = (endI-startI)/60,
             depid = depid) 
    dives$region <- dep_metadata %>% dplyr::filter(ID == depid) %>% .$Region %>% as.character()
    dives$animalID <- dep_metadata %>% dplyr::filter(ID == depid) %>% .$Animal_ID %>% as.character()
    # findsurfs2 in tagtools returns 2 variables
    # startI, endI (NOTE: these are post-dive surfacings)
    surfs <- as.data.frame(prh1hzData$surfs)
    colnames(surfs) <- c('surf_startI', 'surf_endI')
    surfs <- surfs %>% 
      mutate(surf_Duration = (surf_endI-surf_startI)/60) # surface duration in min
    # determine datetime of surfacings
    surfs$surf_dtStart <- prh1hzData_DF$dt[surfs$surf_startI]
    surfs$surf_dtEnd <- prh1hzData_DF$dt[surfs$surf_endI]
    # depending on tag type, pre-deployment time may be identified as a surfacing
    if(length(surfs$surf_startI)>length(dives$DiveNum)){
      surfs <- surfs[-c(1),]
    }
    surfs$surf_DiveNum <- seq.int(nrow(surfs))
    # join Surface interval to dives
    dives <- left_join(dives,surfs, by = c("DiveNum"= "surf_DiveNum"))
    dives <- dives %>% dplyr::filter(dttzStart>=tagon,dttzStart<=tagoff,!is.na(dttzEnd))
    # Reorder the columns
    dives <- dives %>% dplyr::select(depid, DiveNum, region, animalID, everything())
    prh1hzData_DF <- prh1hzData_DF %>%  dplyr::filter(dttz >= tagon ,dttz <= tagoff)
    
    ## Import Movement Data ##
    # Note: The location nearest to the dive start time is selected,
    # but the search time for positions is 10 min, meaning some locations
    # could be associated with multiple dives. 10 min seems to be a reasonable span
    fnLocs <- paste0(str_split(filenames[i]," 1HzprhSC.mat")[[1]][1],"-GPS.csv" )#file.choose()
    whaleLocs <- read_csv(fnLocs,col_types=cols(col_double(),DateTimeUTC = col_datetime(format = ""),
                                      Lat = col_double(),Long = col_double(),lc = col_double()))
    # check that dt imported properly to GMT time
    if(is.null(attr(whaleLocs$DateTimeUTC, "tzone"))){
      whaleLocs$DateTimeUTC <- as.POSIXct(strptime(whaleLocs$DateTimeUTC,format="%m/%d/%Y %H:%M:%S"), tz = 'GMT')
    }
    whaleLocs$locNum <- seq(1:length(whaleLocs$DateTimeUTC)) # create a column for joining later
    whaleLocs <- subset(whaleLocs, select = -c(1, lc)) # remove unused columns
    whaleLocs$depid <- depid
    # find tz_Offset
    whaleLocs$tzOffset <- tzOffset
    # make a Local Datetime column
    whaleLocs$dttz <- whaleLocs$DateTimeUTC
    attr(whaleLocs$dttz, "tzone") <- tzOffset # Set TZ to Local
    # filter to include positions that overlap with dive data (tagon, tagoff)
    whaleLocs <- whaleLocs %>% 
      dplyr::filter(DateTimeUTC >= tagon - minutes(withinTime),
                    DateTimeUTC <= tagoff + minutes(withinTime))
    ## Calculate time differences between locations (lead and lag)
    whaleLocs <- whaleLocs %>% 
      arrange(DateTimeUTC) %>% 
      # time difference between consecutive points
      mutate(time_Diff_lag = pt2pt.back.duration(dttz,output.units='mins'),
             time_Diff_lead = pt2pt.duration(dttz,output.units='mins'))
    
    # Reorder the columns
    whaleLocs <- whaleLocs %>% dplyr::select(depid, locNum,dttz, everything())
    # Attach GPS Data to Dives
    dives$Lat <- NA
    dives$Long <- NA
    dives$locNum <- NA # To be used later as a unique identifier
    dives$locClosest <- NA # Nearest location number
    dives$locTimeDif <- NA # Time diff in minutes to the nearest location
    dives$knownPosGap <- NA # Gap between lead and lag positions
    
    # Find locations withinTime minutes (15) of the start of the dive
    # Note: Can be before or after start of dive (since positions can occur on Initial surfacing
    # as well and terminal dives)
    for(k in 1:dim(dives)[1]){
      temp<- whaleLocs %>% 
        # filter locations withinTime min, before the start of the dive
        dplyr::filter(dives$dtStart[k] >= min(whaleLocs$DateTimeUTC),
                      abs(difftime(dives$dtStart[k],
                                   whaleLocs$DateTimeUTC, units = "secs")) <= (60*withinTime)) %>% 

        # calculate the time difference and arrange smallest first
        mutate(tDif = abs(as.numeric(difftime(dives$dtStart[k], .$DateTimeUTC, units = "secs")/60))) %>% 
        arrange(tDif) 
      if(dim(temp)[1] != 0){
        # select the location nearest to the dive start time
        dives$Lat[k] <- temp$Lat[1]
        dives$Long[k] <- temp$Long[1]
        dives$locNum[k] <- as.character(temp$locNum[1])
        dives$locClosest[k] <- as.character(temp$locNum[1])
        dives$locTimeDif[k] <- temp$tDif[1]
        if(dives$dtStart[k] > temp$DateTimeUTC){
          dives$knownPosGap[k] <- temp$time_Diff_lead[1]
        }else{
          dives$knownPosGap[k] <- temp$time_Diff_lag[1]
        }
        
      }
    }
    rm(k)
    dives <- dives %>% 
      group_by(locNum) %>% 
      arrange(locTimeDif) %>% 
      mutate(locRep = if_else(n()==1 ,0,1),
             locPriority = seq(1,n(),by=1)) %>% 
        ungroup()   %>% 
      arrange(dttzStart) 
    dives$locRep[is.na(dives$Lat)] <- as.numeric(NA)
    dives$locPriority[is.na(dives$Lat)] <- as.numeric(NA)
        
    # Interpolate between locations that are < withintime minutes apart (i.e. 15min gap)
    if(sum(dives$locRep,na.rm=TRUE)>0 & length(dives$dttzStart[!is.na(dives$locRep) & 
                                                               dives$locRep==1 & 
                                                               dives$locPriority > 1 &
                                                               dives$knownPosGap < withinTime &
                                                               dives$dttzStart > min(whaleLocs$dttz)& 
                                                               dives$dttzStart < max(whaleLocs$dttz)])>0){
      ## Create interpolated locations for missing dives
      #  Note: local time used as the timestamps
      locsMove <- move(x=whaleLocs$Long, y=whaleLocs$Lat, time=whaleLocs$dttz, proj=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"),
                       data=whaleLocs, animal=whaleLocs$depid, sensor='unknown')
      ## Create an interpolated location dataframe using only the time for repeated locations
      #  from the previous step. This means that all of these locations were within (withinTime) minutes
      #   Note: Times must be in between the first and last for the interpolation to work
      #         local time used as the timestamps
      locsInterpolated <- interpolateTime(locsMove,
                                          time = dives$dttzStart[!is.na(dives$locRep) & 
                                                                   dives$locRep==1 & 
                                                                   dives$locPriority > 1 &
                                                                   dives$knownPosGap < withinTime &
                                                                   dives$dttzStart > min(whaleLocs$dttz)& 
                                                                   dives$dttzStart < max(whaleLocs$dttz)],
                                          spaceMethod='greatcircle')
      # plot to make sure it looks ok
      plot(locsMove, col="red",pch=20, main="Interpolated Positions")
      points(locsInterpolated)
      #Convert to Dataframe
      locsInterpolated_DF <- as.data.frame(locsInterpolated)
      locsInterpolated_DF$locNum <- as.character(as.numeric(row.names(locsInterpolated_DF))+.1)
      # Reorder the columns
      locsInterpolated_DF <- locsInterpolated_DF %>% 
        dplyr::select(depid,locNum,timestamps,coords.x1, coords.x2)
      colnames(locsInterpolated_DF) <- c("depid","locNum2","dttzStart","longInt","latInt")
      # join and select interpolated position for dives missing a location
      dives <- dives %>% 
        left_join(dplyr::select(locsInterpolated_DF,depid,locNum2,dttzStart,latInt,longInt),by=c("depid","dttzStart")) %>% 
        mutate(Lat = if_else(!is.na(latInt),latInt,Lat),
               Lat = if_else(!is.na(locRep) & locRep==1 & locPriority > 1 & is.na(latInt),as.numeric(NA),Lat),
               Long = if_else(!is.na(longInt),longInt,Long),
               Long = if_else(!is.na(locRep) & locRep==1 & locPriority >1 & is.na(longInt),as.numeric(NA),Long),
               locNum = if_else(is.na(locNum2),locNum,paste0(locNum,".",as.character(locPriority)))) %>% 
        dplyr::select(-longInt, -latInt,-locNum2) 
      rm(locsInterpolated, locsInterpolated_DF) #clean up
    }
    # This cleans up Lat/Longs for cases that do not have any interpolation and are not the priority 
    #  dive (i.e. there is a dive start that is closer to the location time)
    #  Note: No locations are used more than once.
    dives <- dives %>% 
      mutate(Lat = if_else(!is.na(locRep) & locRep==1 & locPriority > 1 & !str_detect(locNum,'\\.'),as.numeric(NA),Lat),
             Long = if_else(!is.na(locRep) & locRep==1 & locPriority >1 & !str_detect(locNum,'\\.'),as.numeric(NA),Long)) %>% 
      dplyr::select(-locRep)
    dives$locNum[is.na(dives$Lat)] <- as.character(NA)
    # red plot should look like the deployment map (grey)
    # some segments will be missing because dives are a subset (>10m)
    plot(whaleLocs$Long,whaleLocs$Lat, col="grey",type="l",
         ylim=c(min(whaleLocs$Lat),max(whaleLocs$Lat)), 
         xlim=c(min(whaleLocs$Long),max(whaleLocs$Long)),
         xlab="Longitude",ylab="Latitude")
    points(dives$Long,dives$Lat,col="red",pch = 20)
    cat(paste0("\n# Dives: ",  length(dives$depid))) 
    cat(paste0("\n# With Locations: ",  sum(!is.na(dives$Lat)))) 
    cat(paste0("\n# Interpolated Positions (< ", withinTime, " mins of 2 known pos): ",  sum(!is.na(dives$Lat)&str_detect(dives$locNum,'\\.')))) 
    
    ## Calculate Distance and Speed Metrics (only for dives with locations)
    diveDistMetrics <- dives %>% 
      filter(!is.na(locNum)) %>% 
      # differences between consecutive points (backward)
      mutate(diveTimeDif_locs = pt2pt.back.duration(dttzStart,output.units='mins'),
             diveDist_locs = pt2pt.back.distance(Lat,Long)/1000, # Dist in KM
             diveSPD_locs = pt2pt.speed(diveDist_locs,diveTimeDif_locs/60), # KM/hr
             diveDIR_locs = pt2pt.back.direction(Lat,Long)) # degrees
    # join metrics to dives
    dives <- left_join(dives,dplyr::select(diveDistMetrics,DiveNum,
                                           diveTimeDif_locs,diveDist_locs, diveSPD_locs, diveDIR_locs), by = "DiveNum")
    rm(diveDistMetrics)
    
    ## Load Lunges ##
    ## Select the lunges.mat file
    lungeFile <- paste0(str_split(filenames[i]," 1HzprhSC.mat")[[1]][1]," lunges.mat" )#file.choose()
    lunges = readMat(lungeFile)
    lungeDT <- round_date(as.POSIXct((lunges$LungeDN - 719529)*86400, origin = "1970-01-01", tz = 'GMT'),unit = "seconds")
    lunge_dttz  <- force_tz(lungeDT,tzone=tzOffset) # Force TZ
    # attr(lunge_dttz, "tzone") #Check tz
    lungeDF <- as.data.frame(lunge_dttz)
    lungeDF$lunge_dt <- lungeDF$lunge_dttz
    attr(lungeDF$lunge_dt, "tzone") <- 'GMT'
    lungeDF$lungeDepth <- lunges$LungeDepth
    lungeDF$lungeC <- lunges$LungeC
    rm(lunge_dttz, lungeDT)
    lungeDF$depid <- depid
    # find tz_Offset
    lungeDF$tzOffset <- tzOffset
    # Select only High-Confidence Lunges Detected
    lungeDF <- lungeDF %>% 
      dplyr::filter(lungeC == 3,lunge_dttz >= tagon, lunge_dttz <= tagoff) %>% 
      dplyr::select(-lungeC) # remove unneeded columns
    # Reorder the columns
    lungeDF <- lungeDF %>% dplyr::select(depid, lunge_dttz,lunge_dt, everything())
    
    # Attach Foraging Data to Dives
    dives$LungeCount <- 0
    dives$meanLungeDepth <- NA
    dives$sdLungeDepth <- NA
    dives$meanILI <- NA # Inter Lunge Interval
    dives$sdILI <- NA
    
    # Find lunges between start and end of dive and 
    for(k in 1:dim(dives)[1]){
      # Attach Lunges to Dives
      temp <- lungeDF %>% 
        filter(lunge_dt >=  dives$dtStart[k] & lunge_dt <=  dives$dtEnd[k]) 
      if(dim(temp)[1] != 0){
        dives$LungeCount[k] <- dim(temp)[1]
        dives$meanLungeDepth[k] <- mean(temp$lungeDepth)
        dives$sdLungeDepth[k] <- sd(temp$lungeDepth)
        # Inter Lunge Interval (only on dives with more than 1 lunge)
        if(dim(temp)[1] > 0){
          dives$meanILI[k] <- mean(pt2pt.back.duration(temp$lunge_dt,output.units='secs'),na.rm=TRUE) 
          dives$sdILI[k] <- sd(pt2pt.back.duration(temp$lunge_dt,output.units='secs'),na.rm=TRUE)
        }
      }
      rm(temp)
    }
    rm(k)
    # make sure 1 lunge dives have a SD of 0
    dives$sdLungeDepth[dives$LungeCount==1] <- 0

    dives %>% 
      filter(LungeCount>0) %>% 
      summarize(Foraging = n(),
                Not_Foraging = length(dives$DiveNum)-n())
    # Save individual data to overall dataset
    if(i==1){
      dataset1hz <- prh1hzData_DF # save the 1Hz data
      diveDataset <- dives # save the dives
      locDataset <- whaleLocs
      lungeDataset <- lungeDF
      } else{
        dataset1hz <- rbind(dataset1hz,prh1hzData_DF)
        diveDataset <- rbind(diveDataset,dives)
        locDataset <- rbind(locDataset,whaleLocs)
        lungeDataset <- rbind(lungeDataset,lungeDF)
      } 
    rm(lungeDF, prh1hzData_DF, prh1hzData)
    rm(tagon, tagoff, lunges, depid, tzOffset)
  } 
  rm(lungeFile,pathChoice,i, surfs, dives, whaleLocs, filename,filenames, fnLocs, locsMove)

  ## Save Intermediate Dataframe ##
  # Check for NAs
  sapply(diveDataset, function(x) sum(is.na(x)))
  ## Save diveDataset to reload later
  ## Note: This dataframe has all dives, not just those with Locations
  save(diveDataset,file="./dataProcessed/diveDataset.RData")
  ## Save dataset1hz for Cast processing 
  save(dataset1hz, file="./dataProcessed/dataset1hz.RData")
  ## Save locDataset
  save(locDataset, file="./dataProcessed/locDataset.RData")
  ## Save lungeDataset
  save(lungeDataset, file="./dataProcessed/lungeDataset.RData")
 
#### Process Surface Current Data ####
 ## Download HF Radar data, Import Raw Radar data, Extract Radar Precision Metrics (HDOP, Radials, Sites),
  # NOTE: This script downloads SC data with a numDays buffer on either end of the deployment, and
  # adds a degreeBuffer around the deployment points to avoid edge-effects
  useOldNC <- TRUE # if true, script will search for a local NC file before downloading a new one
  numDays <- 7 # number of days of HF Radar data before and after deployment
  degreeBuffer <- 1 # number of degrees around the min/max lat/longs of deployment
  ncPattern <- paste0("_HFRADAR_US_West_Coast_", gridSize, "km_Resolution_Hourly.nc")
  for(i in 1:length(unique(diveDataset$depid))){
    # Get the depid and time zone for this deployment
    depid <- unique(diveDataset$depid)[i]
    # find tz_Offset, region and animalID (from external file)
    tzOffset <- dep_metadata %>% filter(ID == depid) %>% .$tzOffset %>% as.character()
    cat(paste0("\nProcessing: ",  unique(diveDataset$depid)[i]))
    cat(paste0("\nDeployment ",  i, " of ", length(unique(diveDataset$depid))))
    # filter by ID and remove dives without location data
    data_Sub <- diveDataset %>% filter(depid == unique(diveDataset$depid)[i],
                                       !is.na(Lat),!is.na(Long))
    # Determine the start and end times of deployment, +- 1 hour
    start <- first(data_Sub$dtStart)-(3600) 
    end <- last(data_Sub$dtEnd)+(3600)
    # Only download new NC file if it doesn't exist locally
    if(useOldNC & file.exists(paste0("./dataRaw/RawHFRadarData/",depid,ncPattern))){ 
      ncname <- paste0("./dataRaw/RawHFRadarData/",depid,ncPattern)
      cat(paste0("Using previously downloaded HF Radar data:\n",ncname))
    }else{ # Download the appropriate NC file using parameters from tag 
      ncname <- hfRadar_Download(paste0("dataRaw/RawHFRadarData/", unique(diveDataset$depid)[i]),
                                 dStart = lubridate::date(start)-numDays,tStart = paste0(lubridate::hour(start),':00:00'),
                                 dEnd = lubridate::date(end)+numDays, tEnd = paste0(lubridate::hour(end),':00:00'), 
                                 bbox = c(round(max(data_Sub$Lat, na.rm = TRUE)+degreeBuffer,2),round(min(data_Sub$Lat, na.rm = TRUE)-degreeBuffer,2),
                                          round(min(data_Sub$Long, na.rm = TRUE)-degreeBuffer,2),round(max(data_Sub$Long, na.rm = TRUE)+degreeBuffer,2) ), 
                                 grid = gridSize)
    }
    rm(ncname)
  }
   
  
#### Process FTLE  ####
  ## Import and Extract (NetCDF) ##
  ##  NOTE: FTLE Calculated Externally
  
  ## Import the externally processed FTLE data (netCDF) and Extract Values at Whale Locations
  ## Restoration Settings:
  #   Data was restored using automatic land detection and restored up to shoreline. 
  #   Concave Hull around available points with 10km alpha
  ## FTLE Settings: 
  #   Include Land and Apply a Free-Slip boundary condition
  #   Initial Time: 1 day before end of data (i.e. 1 day before today)
  #   Number of Initial Times: 24 
  #   Initial Times Increment: 1 Hour
  #   Integration Duration: 4 days
  #   Reverse Direction of Time 
  ## FTLE is calculated with a grid of tracers having 10x the resolution of the input data.
  #   The final resolution varies slightly by deployment, but is roughly 590m or 6km HF radar resolution/10.
  #   All FTLE Files for a deployment (i.e. 12, 24, 48 and 96Hr) have the same temporal and spatial resolution.
  FTLEdatapath <- "./dataRaw/FTLE/"
  FTLE_Integrations <- c("48")
  # Create empty dataframes to hold all the data
  datasetFTLE <- data.frame()
  for(i in 1:length(unique(diveDataset$depid))){
    # Get the depid and time zone for this deployment
    depid <- unique(diveDataset$depid)[i]
    # find tz_Offset, region and animalID (from external file)
    tzOffset <- dep_metadata %>% filter(ID == depid) %>% .$tzOffset %>% as.character()
    cat(paste0("\nProcessing: ",  unique(diveDataset$depid)[i]))
    cat(paste0("\nDeployment ",  i, " of ", length(unique(diveDataset$depid))))
    ptm <- proc.time()
    # filter by ID and remove dives without location data
    data_Sub <- diveDataset %>% filter(depid == unique(diveDataset$depid)[i],
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
      save(FTLE_stack,file= paste0("./dataProcessed/",unique(data_Sub$depid), "_FTLEstack-",FTLE_Integrations[ii], "hrs.RData"))
      
      ### Raster Extraction ###
      ## Create a spatial points dataframe of whale movement
      xy <- coordinates(cbind(data_Sub$Long, data_Sub$Lat))
      dataSP <- SpatialPointsDataFrame(xy, data_Sub, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"), match.ID = FALSE)
      cat("\nExtracting FTLE Data")
      # Create an empty dataframe for extracted FTLE values
      ftle_Extract <- data.frame()
      for(j in 1:length(time_FTLE)){ 
        #cat(paste0("\nLayer: ",  j))
        # subset by the time of FTLE observation (+- timeDifs)
        tempSP <- dataSP[(dataSP$dtStart <= time_FTLE[j]+timeDifS & dataSP$dtStart > time_FTLE[j]-timeDifS),]
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
        ftle_ExtractAll <- cbind(ftle_ExtractAll,dplyr::select(ftle_Extract, starts_with("ftle")))
      } 
      nc_close(ncin_FTLE)
      rm(ncin_FTLE, ncname,FTLE_stack,time_FTLE,time_FTLE_Local,ftle_Extract,xy,dataSP)
      tPTM1 <-  proc.time() - ptm1
      cat(paste0("\nElapsed time: ", round(tPTM1[3],1), " seconds"))
      rm(ptm1,tPTM1)
    }
    # Save to an output dataframe 
    if(i==1){
      datasetFTLE <- ftle_ExtractAll 
    } else{
      datasetFTLE <- rbind(datasetFTLE, ftle_ExtractAll)
    } 
    tPTM <-  proc.time() - ptm
    cat(paste0("\nOverall elapsed time: ", round(tPTM[3]/60,1), " minutes"))
    rm(ptm, tPTM)
    gc()
  }
  # clean up unused variables
  rm(i,ii,tzOffset, depid, ftle_ExtractAll, FTLE_Integrations,FTLEdatapath,data_Sub)

  # Remove the remnants of the Spatial Points DF conversion
  datasetFTLE <- datasetFTLE %>% dplyr::select(-c(coords.x1, coords.x2))   
  # Check for INF
  sapply(datasetFTLE, function(x) sum(is.infinite(x))>0)
  sapply(datasetFTLE, function(x) sum(is.na(x))>0)
  ##Save Data for later import
  save(datasetFTLE,file="./dataProcessed/datasetFTLE.RData")
  ##Save Data for Analysis
  save(datasetFTLE,file="./rmdFiles/datasetFTLE.RData")
  