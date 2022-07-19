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
pkgTest("viridis")
pkgTest("scales")
pkgTest("rnaturalearth")
pkgTest("metR")
pkgTest("ggspatial")
pkgTest("sf")


# Process netCDF file dates and times
source("./functions/ncdate.R")
# pkgTest("metR") # note: getNOAA.bathy was not working, but an alternate version on github is used here
source('./functions/getNOAAbathy.R')

world <- ne_countries(scale = "large", returnclass = "sf")

#### Setup ####

load("./dataProcessed/datasetFTLE.RData") # Load the processed dive+FTLE dataframe
# Create a dataframe with no missing FTLE vales and columns for grouping and plotting
data_Gamm_ALL <- datasetFTLE %>% 
  # Only include Dives with FTLE data
  dplyr::filter(is.nan(ftle48)==0, is.na(ftle48)==0) %>% 
  group_by(depid) %>% 
  # Add a time Variable for the AR1 correlation structure 
  mutate(time = startI-min(startI),
         timeHrs = time/3600,
         # Add a binary Foraging field for Logistic Regression 
         Foraging = if_else(LungeCount > 0,1,0),
         # Add a hourly feeding rate for use in HMM
         total_Duration = dive_Duration + surf_Duration,
         feedingRateHR = LungeCount/total_Duration*60) %>% 
  ungroup()
load("./Figure_files/Figure2-MapData/bf.RData")
load("./Figure_files/Figure2-MapData/hfRadarStations.RData")
# HF Radar Coverage Area and Station Locations
hfRadarStns$name <- "Radar Station"
# filter by ID and remove dives without location data
data_Sub <- data_Gamm_ALL %>% filter(depid == "Bm170926-TDR14",
                                   !is.na(Lat),!is.na(Long))

FTLEdatapath <- "./dataRaw/FTLE/"
depid <- "Bm170926-TDR14"
tzOffset <- "Etc/GMT+7"

ncin_FTLE48 <- nc_open("./Figure_files/Figure2-FTLE/Bm170926-TDR14_FTLE_48Hrs.nc", verbose=FALSE)
# FTLE Time (All have the same time dimension)
time_FTLE <- ncdate(ncin_FTLE48)
# convert to local time 
time_FTLE_Local <- time_FTLE
attr(time_FTLE_Local, "tzone") <- tzOffset # change the timezone to tzOffset
# Create a FTLE Raster Stack for quantiles
FTLE_stack48 <- stack("./Figure_files/Figure2-FTLE/Bm170926-TDR14_FTLE_48Hrs.nc", varname = "FTLE")
# Feeding Colors
# Dark 3 green and purple
# colFeed <- c("Feeding" = "#00AA5A", "Not Feeding" = "#B675E0")
colFeed <- c("Feeding" = "#00FF00", "Not Feeding" = "#66CCFF")
# # For Sunrise Sunset Times in Plot
meanLat <- mean(data_Sub$Lat)
meanLong<- mean(data_Sub$Long)
data<- data.frame(time_FTLE_Local)
data$dttz <- time_FTLE_Local
data$date <- date(data$time_FTLE_Local)
data$lat <- meanLat
data$long <- meanLong
# Determine Solar Position for each point
source('./functions/SunriseSunsetTimes.R')
colfunc <- colorRampPalette(c("black", "white"))
# colfunc(10)
data$shade <- ifelse(data$solarpos >= 15,"#FFFFFF",
                     ifelse(data$solarpos < 15 & data$solarpos >= 10, "#DADADA",
                            ifelse(data$solarpos < 10 & data$solarpos >= 5, "#B6B6B6",
                                   ifelse(data$solarpos < 5 & data$solarpos >= 0, "#919191",
                                          ifelse(data$solarpos < 0 & data$solarpos >= -5, "#6D6D6D",
                                                 ifelse(data$solarpos < -5 & data$solarpos >= -10, "#484848",
                                                        ifelse(data$solarpos < -5 & data$solarpos >= -10, "#242424","#000000")))))))
data$textShade <- ifelse(data$solarpos >= -5,"#000000","#FFFFFF")

# Create a placeholder Dataframe for Feeding in legend
feedDF <- data.frame(Lat = rep(meanLat,2),Long = rep(meanLong,2),feed = c("Feeding","Not Feeding"))

#### Plot All ####
# FTLE Color Scaling
qFTLE<-quantile(FTLE_stack48,probs = seq(.5,.95,by=.05),na.rm=TRUE)
qSumFTLE<-as.data.frame(summary(qFTLE)) %>% dplyr::select(-Var1) %>% 
  dplyr::filter(grepl("Mean",Freq)) %>%
  group_by(Var2) %>% 
  mutate(Means = unlist(str_split(Freq,":"))[2]) %>% 
  ungroup()

degreeBuffer <- 0.25
cat("Plotting data for", depid, "from", as.character(time_FTLE_Local[70]),"to", as.character(time_FTLE_Local[130]))
tOI <- 70:130
ftle_alpha <- function(n, flatten, scrunch) {
  c(rescale(exp((1:flatten)/scrunch), to = c(0.25, 1)), rep(1, n - flatten))
}
for (iii in tOI) {
  cat(paste0("\nPlot ",  iii-(min(tOI)-1), " of ", length(tOI)))
  # create a dataframe of each integration
  FTLE48_spdf <- as(FTLE_stack48[[iii]], "SpatialPixelsDataFrame")
  FTLE48_df <- as.data.frame(FTLE48_spdf)
  colnames(FTLE48_df) <- c("FTLE", "x", "y")
  bwOI <- data_Sub %>% 
    filter(dtStart<=time_FTLE[iii]+1800 & dtStart>=time_FTLE[iii]-1800 )
  bwOI <- bwOI %>% 
    mutate(feed = ifelse(LungeCount > 0, "Feeding", "Not Feeding"))
  # FTLE PLOT
  p1<-ggplot() +
    geom_sf(data = world,color=NA,fill=NA) +
    geom_raster(data = FTLE48_df[complete.cases(FTLE48_df),], 
                mapping = aes(x, y, fill = FTLE)) +
    labs(x="", y="") +
    scale_fill_gradientn(name = 'FTLE',
                         colors = viridis(option = "plasma", n = 40, alpha = ftle_alpha(40, 8, 8)),
                         values = rescale(c(0,.1,as.numeric(qSumFTLE$Means[1:9]))),
                         limits=c(0, .8), 
                         breaks = c(0,.2,.4,.6,.8),
                         # labels = c()
                         oob=scales::squish,
                         na.value = NA) +
    # # add 100m and 200m contours
    # ggplot2::geom_contour(data = bf,
    #                       aes(x=x, y=y, z=z),
    #                       breaks=c(-100),
    #                       size=c(0.4),
    #                       colour="darkgray", show.legend = FALSE) +
    # metR::geom_text_contour(data = bf, aes(x=x,y=y,z=z),breaks = c(-100),label.placement = label_placement_n(5),
    #                         show.legend = FALSE, size = 2.2, alpha = .8, nudge_y = -.002, color="lightgrey") +
    # ggplot2::geom_contour(data = bf,
    #                       aes(x=x, y=y, z=z),
    #                       breaks=c(-200),
    #                       size=c(0.4),
    #                       colour="darkgray", show.legend = FALSE) +
    # metR::geom_text_contour(data = bf, aes(x=x,y=y,z=z),breaks = c(-200),label.placement = label_placement_n(5),
    #                         show.legend = FALSE, size = 2.2, alpha = .8, nudge_y = -.002, color="lightgrey") +
  # add 100m and 200m contours
  ggplot2::geom_contour(data = bf,
                        aes(x=x, y=y, z=z),
                        breaks=c(-100),
                        size=c(0.5),
                        colour="darkgray", show.legend = FALSE) +
    metR::geom_text_contour(data = bf, aes(x=x,y=y,z=z),breaks = c(-100),label.placement = label_placement_n(5),
                            show.legend = FALSE, size = 2.2, alpha = .9, nudge_y = -.002, color="black") +
    ggplot2::geom_contour(data = bf,
                          aes(x=x, y=y, z=z),
                          breaks=c(-200),
                          size=c(0.5),
                          colour="darkgray", show.legend = FALSE) +
    metR::geom_text_contour(data = bf, aes(x=x,y=y,z=z),breaks = c(-200),label.placement = label_placement_n(5),
                            show.legend = FALSE, size = 2.2, alpha = .9, nudge_y = -.002, color="black") +
    geom_point(data = bwOI, mapping = aes(Long,Lat,color=feed,shape=feed,stroke=feed),
               size=1.75, show.legend = FALSE, alpha = 1,#shape=21,stroke=2
               fill="black") +
    geom_point(data = feedDF, aes(Long-1000,Lat+1000, color = feed,shape=feed,stroke=feed),
               size=1.75, show.legend = TRUE, alpha = 1,#shape=21,stroke=2
               fill="black")+
    geom_sf(data = world, size=1.25,fill="lightgrey",color="grey") +
    geom_point(data=hfRadarStns, aes(st_coordinates(hfRadarStns$geometry)[,1],
                                     st_coordinates(hfRadarStns$geometry)[,2],
                                     size = name), 
               color="black",size=2.5,stroke=1.5,
               fill = 'white', shape=23, show.legend = FALSE) +
    scale_shape_manual(values = c("Feeding" = 21, "Not Feeding" = 22), name = "Blue Whale Locations") +
    scale_discrete_manual("stroke",  values = c("Feeding" = 2, "Not Feeding" = 2), name = "Blue Whale Locations") +
    scale_size_manual(values = c("Radar Station" = 3), name="High Frequency Radar") +
    coord_sf(xlim = c(-120.7, -119.5), 
             ylim = c(33.88, 34.48), expand = FALSE) +
    scale_color_manual(values = colFeed, name = "Blue Whale Locations") +
    theme_classic() +
    annotation_scale(location = "bl", width_hint = 0.3, text_face = "bold", 
                     text_col = "black",text_cex = 1.2) + 
    annotation_north_arrow(location = "bl", 
                           which_north = "true", 
                           pad_x = unit(0.25, "in"), 
                           pad_y = unit(0.25, "in"), 
                           style = north_arrow_fancy_orienteering) + 
    xlab("Longitude") + 
    ylab("Latitude") + 
    labs(colour = data$textShade[iii]) +
    guides(color = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
           shape = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
           stroke = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
           size = guide_legend(order=3,direction = 'vertical', spacing.y=unit(.1, "cm"), 
                               override.aes = list(linetype = 1,color="black",size=2.5,stroke=1.5,
                                                   fill = 'white', shape=23)),
           fill = guide_colourbar(order=2,title.position = "top",
                                  title.hjust = .5,
                                  label.position = "bottom") )+
    ggtitle("Lagrangian Features and Blue Whale Foraging",  subtitle =  sprintf("%s", time_FTLE_Local[iii])) + 
    theme(plot.background = element_rect(fill = data$shade[iii]),
          panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
          panel.background = element_rect(fill = "aliceblue"),
          axis.title = element_text(family="Times",face="bold", size=18, color = data$textShade[iii]),
          axis.text = element_text(family="Times",face="bold", size=16, color = data$textShade[iii]),
          legend.text = element_text(face="bold", size=16),
          legend.position = "bottom", legend.box = 'horizontal',
          legend.key.size = unit(.75, "cm"),
          legend.title = element_text(face="bold", size=20),
          legend.box.background = element_rect(color="white", size=1,fill="white"),
          legend.spacing.y = unit(.1, "cm"),legend.spacing.x = unit(1, "cm"),
          plot.margin=unit(c(0.5,2.5,0.5,0.25), "cm"), # adds a buffer on the sides so the plot is centered
          plot.title = element_text(face="bold", size=23,color=data$textShade[iii], hjust = 0.5),
          plot.subtitle = element_text(face="bold", size=20,color=data$textShade[iii], hjust = 0.5))
  p1
  ggsave(sprintf("Output/GIF_files/%s_FTLE_%03i.png", depid, iii ), width=12, height=9, units = "in", dpi=600)
}



