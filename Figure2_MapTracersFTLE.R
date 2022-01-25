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
pkgTest("viridis")
pkgTest("scales")
# for the BW data
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

world <- ne_countries(scale = "large", returnclass = "sf")

load("./dataProcessed/datasetFTLE.RData") # Load the processed dive+FTLE dataframe

#### Panel A - Map of Deployments ####
hourlySum <- datasetFTLE %>% 
  mutate(yday = yday(dttzStart),
         hr = hour(dttzStart)) %>% 
  group_by(depid,yday, hr) %>% 
  summarize(LungeCount = sum(LungeCount),
            Lat = mean(Lat,na.rm=TRUE),
            Long = mean(Long,na.rm=TRUE))
p1<-ggplot() +
  geom_sf(data = world,color=NA,fill=NA,aes(geometry = geometry)) +
  geom_polygon(data=fortify(HFRadarCoverageMCP),
               aes(x=long,y=lat, alpha = name),
               color = NA, fill = "blue") +
  # add 200m contour
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-200),
               size=c(0.5),
               colour="darkgrey", show.legend = FALSE) +
  geom_text_contour(data = bf, aes(x=x, y=y,z = z),breaks=c(-200),color='darkgrey', 
                    show.legend = FALSE, size = 4.2, alpha = .6, nudge_y = -.002) +
  geom_sf(data = world, fill='grey',aes(geometry = geometry)) +
  geom_polygon(aes(x = c(-120.7,-119.5,-119.5,-120.7), y = c(33.88,33.88,34.48,34.48)),
               size = 1.5,linetype=2,fill=NA, color='black') +
  # Add whale Locations
  geom_path(data=data_Gamm_ALL, aes(x=Long,y=Lat,group=depid,color=depid),
            alpha = 0.95,size = .9, show.legend = TRUE) +
  geom_point(data=hourlySum[hourlySum$LungeCount>0,],
             aes(x=Long,y=Lat,group = depid,
                 size = LungeCount, fill = depid),
             alpha = 0.75, shape=21,show.legend = TRUE) +

  geom_point(data=hfRadarStns,
             aes(st_coordinates(hfRadarStns$geometry)[,1],st_coordinates(hfRadarStns$geometry)[,2],
                 shape = name),size=2.5, fill = 'magenta') +
  scale_alpha_manual(values=c("Coverage Area"=.15))+
  scale_shape_manual(values = c("Radar Station" = 23))+
  coord_sf(xlim = c(min(data_Gamm_ALL$Long,na.rm = TRUE)-.05, max(data_Gamm_ALL$Long,na.rm = TRUE)+.05),
           ylim = c(min(data_Gamm_ALL$Lat,na.rm = TRUE)-.05, max(data_Gamm_ALL$Lat,na.rm = TRUE)+.05),
           expand = FALSE) +

  labs(size = "Hourly Lunges",alpha=NULL, shape="HF Radar", color="Deployment",fill="Deployment") +
  guides(size = guide_legend(order=4, override.aes = list(linetype = 0)),
         alpha = guide_legend(order=2, override.aes = list(linetype = 0,size=0,color=NA)),
         shape = guide_legend(order=1, override.aes = list(linetype = 0)),
         color = guide_legend(order=3, override.aes = list(size=1.5,shape=NA,fill=NA)),
         fill = guide_legend(order=3, override.aes = list(size=3,shape=21)))+
  xlab("Longitude") + 
  ylab("Latitude") + 
  theme(axis.title = element_text(family="Times",face="bold", size=20),
        axis.text = element_text(family="Times", face="bold", size=18),
        panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(face="bold", size=16),
        legend.title = element_text(face="bold", size=20,hjust = 0.5),
        plot.title = element_text(face="bold", size=24,hjust = 0.5),
        legend.position= c(.8,.7), 
        legend.key.size = unit(.5, "cm"),
        legend.box.background = element_rect(color="white", size=1,fill="white"),
        # legend.box.margin = margin(6, 6, 6, 6),
        legend.box = "vertical", legend.key = element_blank(),
        legend.margin = margin(0.25,0.15,0.25,0.15, unit="cm"),
        legend.spacing.y = unit(.1, "cm"))
p1
ggsave("./Output/FTLE_Map_Fig2a.png", width=8, height=14, units = "in", dpi=400)


#### Panel B - FTLE Tracers ####
# We are wanting the tracers for 9.28.2017 14:00 PST (Layer 77) to align with the FTLE in Figure 1

depid <- "Bm170926-TDR14"
tzOffset <- "Etc/GMT+7"

# File to be converted
ncname_Tracers <- "./Figure_files/Figure2-Tracers/Bm170926-TDR14_Tracers.nc"# file.choose()
convertName <- str_remove(basename(ncname_Tracers),"\\.nc")
ncin_Tracers <- nc_open(ncname_Tracers)
tracerLatitudes <- ncvar_get(ncin_Tracers,"Tracer_Latitudes")
tracerLongitudes <- ncvar_get(ncin_Tracers,"Tracer_Longitudes")

time_Tracers <- ncdate(ncin_Tracers)
# convert to local time 
time_Tracers_Local <- time_Tracers
attr(time_Tracers_Local, "tzone") <- tzOffset # change the timezone to tzOffset
attr(time_Tracers_Local, "tzone") # check that dttz is in local time

Lats <- as.vector(tracerLatitudes[,,length(time_Tracers_Local)])
Longs <- as.vector(tracerLongitudes[,,length(time_Tracers_Local)])
tracerLocs <- data.frame(Lats,Longs) %>% filter(!is.na(Lats))
#Find the index of the time we want to plot
i <- which(time_Tracers_Local=="2017-09-28 14:00:00 -07") #49

hfRadarStns$name <- "Radar Station"
bwOI <- data_Sub %>% 
  filter(dttzStart<=time_Tracers_Local[i]+1800 & dttzStart>=time_Tracers_Local[i]-1800 )
bwOI <- bwOI %>% 
  mutate(feed = ifelse(LungeCount > 0, "Feeding", "Not Feeding"))
bwOISum <- bwOI %>% 
  summarize(Lat = mean(Lat),
            Long = mean(Long),
            Lunges = sum(LungeCount))
pT<-ggplot() +
  geom_sf(data = world,color=NA,fill=NA) +
  geom_point(data=tracerLocs, aes(x=Longs,y=Lats  ),shape=21,size=1.75,color='black',fill='grey')+
  geom_sf(data = world, size=1.25,fill="darkgrey",color="black") +
  geom_point(data=hfRadarStns,
             aes(st_coordinates(hfRadarStns$geometry)[,1],st_coordinates(hfRadarStns$geometry)[,2],
                 shape = name),
             fill = 'magenta', shape=23, show.legend = FALSE) +
  scale_size_manual(values = c("Radar Station" = 3), name="High Frequency Radar") +
  coord_sf(xlim = c(-120.7, -119.5), 
           ylim = c(33.88, 34.48), expand = FALSE) +
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
  guides(color = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
         shape = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
         stroke = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
         size = guide_legend(order=3,direction = 'vertical', spacing.y=unit(.1, "cm"), 
                             override.aes = list(linetype = 1,shape=23)),
         #Right side colorbar
         # fill = guide_colourbar("FTLE", order=1, title.position = "top",
         #                        title.hjust = .5, direction = 'vertical',
         #                        barwidth = unit(1.5, "cm"),barheight = unit(8, "cm"))) + 
         #Bottom Colorbar
         fill = guide_colourbar("FTLE", order=1, title.position = "top",
                                title.hjust = .5,direction = 'horizontal',
                                label.position = "bottom",
                                barwidth = unit(10, "cm"),barheight = unit(1, "cm"))) + #
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_text(face="bold", size=18),
        axis.text = element_text(face="bold", size=16),
        legend.text = element_text(face="bold", size=18),
        # legend.position = c(1.1,0.5), # Right side colorbar
        legend.position = "bottom",
        legend.box = 'horizontal',
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(face="bold", size=24),
        legend.box.background = element_rect(color="white", size=1,fill="white"),
        legend.spacing.y = unit(.25, "cm"),legend.spacing.x = unit(.25, "cm"),
        plot.margin=unit(c(0,2.5,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        plot.title = element_text(face="bold", size=24, hjust = 0.5),
        plot.subtitle = element_text(face="bold", size=20, hjust = 0.5))
pT
ggsave(sprintf("Output/%s_Tracers_%03i.png", "Bm170926-TDR14", i ), width=12, height=8, units = "in", dpi=400)

#### Panel C - FTLE ####
FTLEdatapath <- "./dataRaw/FTLE/"
depid <- "Bm170926-TDR14"
# filter by ID and remove dives without location data
data_Sub <- datasetFTLE %>% filter(depid == "Bm170926-TDR14",
                                   !is.na(Lat),!is.na(Long))

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
colFeed <- c("Feeding" = "#00AA5A", "Not Feeding" = "#B675E0")
# # For Sunrise Sunset Times in Plot
meanLat <- mean(data_Sub$Lat)
meanLong<- mean(data_Sub$Long)
b = getNOAA.bathy(lon1 = min(data_Sub$Long-.5,na.rm=TRUE), lon2 = max(data_Sub$Long+.5,na.rm=TRUE), 
                  lat1 =max(data_Sub$Lat+.5,na.rm=TRUE), lat2 = min(data_Sub$Lat-.5,na.rm=TRUE),resolution = 1)
# convert bathymetry to data frame
bf = fortify.bathy(b)
# Create a placeholder Dataframe for Feeding in legend
feedDF <- data.frame(Lat = rep(meanLat,2),Long = rep(meanLong,2),feed = c("Feeding","Not Feeding"))
# FTLE Color Scaling
qFTLE<-quantile(FTLE_stack48,probs = seq(.5,.95,by=.05),na.rm=TRUE)
qSumFTLE<-as.data.frame(summary(qFTLE)) %>% dplyr::select(-Var1) %>% 
  dplyr::filter(grepl("Mean",Freq)) %>%
  group_by(Var2) %>% 
  mutate(Means = unlist(str_split(Freq,":"))[2]) %>% 
  ungroup()

degreeBuffer <- 0.25

#Find the index of the time we want to plot
ii <- which(time_FTLE_Local=="2017-09-28 14:00:00 -07") #77
FTLE48_spdf <- as(FTLE_stack48[[ii]], "SpatialPixelsDataFrame")
FTLE48_df <- as.data.frame(FTLE48_spdf)
colnames(FTLE48_df) <- c("FTLE", "x", "y")

bwOI <- data_Sub %>% 
  filter(dtStart<=time_FTLE[ii]+1800 & dtStart>=time_FTLE[ii]-1800 )
bwOI <- bwOI %>% 
  mutate(feed = ifelse(LungeCount > 0, "Feeding", "Not Feeding"))
bwOISum <- bwOI %>% 
  summarize(Lat = mean(Lat),
            Long = mean(Long),
            Lunges = sum(LungeCount))
# FTLE PLOT
ftle_alpha <- function(n, flatten, scrunch) {
  c(rescale(exp((1:flatten)/scrunch), to = c(0.25, 1)), rep(1, n - flatten))
}
pF<-ggplot() +
  geom_sf(data = world,color=NA,fill=NA) +
  geom_raster(data = FTLE48_df[complete.cases(FTLE48_df),], 
              mapping = aes(x, y, fill = FTLE)) + 
  labs(x="", y="") +
  scale_fill_gradientn(name = 'FTLE',
                       colors = viridis(option = "plasma", n = 40, alpha = ftle_alpha(40, 8, 8)),
                       values = rescale(c(0,.1,as.numeric(qSumFTLE$Means[1:9]))),
                       limits=c(0, .8), 
                       breaks = c(0,.2,.4,.6,.8),
                       oob=scales::squish,
                       na.value = NA) +
  # add 100m and 200m contours
  ggplot2::geom_contour(data = bf,
                        aes(x=x, y=y, z=z),
                        breaks=c(-100),
                        size=c(0.7),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf, aes(x=x,y=y,z=z),breaks = c(-100),label.placement = label_placement_n(5),
                          show.legend = FALSE, size = 4.2, alpha = .9, nudge_y = -.002, color="black") +
  ggplot2::geom_contour(data = bf,
                        aes(x=x, y=y, z=z),
                        breaks=c(-200),
                        size=c(0.7),
                        colour="darkgray", show.legend = FALSE) +
  metR::geom_text_contour(data = bf, aes(x=x,y=y,z=z),breaks = c(-200),label.placement = label_placement_n(5),
                          show.legend = FALSE, size = 4.2, alpha = .9, nudge_y = -.002, color="black") +
  geom_path(data=data_Sub, aes(Long,Lat),size=1.25,alpha=1, color="#00BE70")+
  geom_point(data = bwOISum, mapping = aes(Long,Lat),stroke=4, color="black", shape = 23, size = 12,
             show.legend = FALSE, alpha = 1,  
             fill="#00BE70") +
  geom_label(data = bwOISum, mapping = aes(x=Long,y=Lat, label = Lunges),label.size = NA, fill=NA, fontface="bold")+
  geom_sf(data = world, size=1.25,fill="darkgrey",color="black") +
  geom_point(data=hfRadarStns, aes(st_coordinates(hfRadarStns$geometry)[,1],
                                   st_coordinates(hfRadarStns$geometry)[,2],
                                   size = name), 
             fill = 'magenta', shape=23, show.legend = FALSE) +
  scale_shape_manual(values = c("Feeding" = 21, "Not Feeding" = 22), name = "Blue Whale Locations") +
  scale_discrete_manual("stroke",  values = c("Feeding" = 2, "Not Feeding" = 2), name = "Blue Whale Locations") +
  scale_size_manual(values = c("Radar Station" = 3), name="High Frequency Radar") +
  coord_sf(xlim = c(-120.7, -119.5), 
           ylim = c(33.88, 34.48), expand = FALSE) +
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
  guides(color = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
         shape = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
         stroke = guide_legend(order=1,direction = 'vertical', spacing.y=unit(.1, "cm")),
         size = guide_legend(order=3,direction = 'vertical', spacing.y=unit(.1, "cm"), 
                             override.aes = list(linetype = 1,shape=23)),
         #Bottom Colorbar
         fill = guide_colourbar("FTLE", order=1, title.position = "top",
                                title.hjust = .5,direction = 'horizontal',
                                label.position = "bottom",
                                barwidth = unit(10, "cm"),barheight = unit(1, "cm"))) + #
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_text(face="bold", size=18),
        axis.text = element_text(face="bold", size=16),
        legend.text = element_text(face="bold", size=18),
        # legend.position = c(1.1,0.5), # Right side colorbar
        legend.position = "bottom",
        legend.box = 'horizontal',
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(face="bold", size=24),
        legend.box.background = element_rect(color="white", size=1,fill="white"),
        legend.spacing.y = unit(.25, "cm"),legend.spacing.x = unit(.25, "cm"),
        plot.margin=unit(c(0,2.5,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        plot.title = element_text(face="bold", size=24, hjust = 0.5),
        plot.subtitle = element_text(face="bold", size=20, hjust = 0.5))
pF
ggsave(sprintf("./Output/%s_FTLE_%03i.png", depid, ii ), width=12, height=8, units = "in", dpi=400)

# Plot just the color bar 
pCB<-ggplot() +
  geom_sf(data = world,color=NA,fill=NA) +
  geom_raster(data = FTLE48_df[complete.cases(FTLE48_df),], 
              mapping = aes(x, y, fill = FTLE)) + 
  labs(x="", y="") +
  scale_fill_gradientn(name = 'FTLE',
                       colors = viridis(option = "plasma", n = 40, alpha = ftle_alpha(40, 8, 8)),
                       values = rescale(c(0,.1,as.numeric(qSumFTLE$Means[1:9]))),
                       limits=c(0, .8), 
                       breaks = c(0,.2,.4,.6,.8),
                       oob=scales::squish,
                       na.value = NA) +
  geom_sf(data = world, size=1.25,fill="darkgrey",color="black") +
  coord_sf(xlim = c(-120.7, -119.5), 
           ylim = c(33.88, 34.48), expand = FALSE) +
  theme_classic() +
  guides(fill = guide_colourbar("FTLE", order=1, title.position = "top",
                                title.hjust = .5, direction = 'vertical',
                           barwidth = unit(1, "cm"), barheight = unit(8, "cm"))) + 
  # fill = guide_colourbar("FTLE", order=1, title.position = "left",
  #                               direction = 'horizontal',#title.vjust = 'top',
  #                               label.position = "bottom",
  #                               barwidth = unit(10, "cm"),
  #                               barheight = unit(1, "cm"))) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_text(face="bold", size=18),
        axis.text = element_text(face="bold", size=16),
        legend.text = element_text(face="bold", size=18),
        # legend.position = c(1.1,0.5), # Right side colorbar
        legend.position = "right",
        legend.box = 'horizontal',
        legend.key.size = unit(1, "cm"),
        legend.title = element_text(face="bold", size=24),
        legend.box.background = element_rect(color="white", size=1,fill="white"),
        legend.spacing.y = unit(.25, "cm"),legend.spacing.x = unit(.5, "cm"),
        plot.margin=unit(c(0,2.5,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        plot.title = element_text(face="bold", size=24, hjust = 0.5),
        plot.subtitle = element_text(face="bold", size=20, hjust = 0.5))
pCB
ggsave(sprintf("./Output/%s_FTLE_Legend.png", depid), width=12, height=8, units = "in", dpi=400)