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

load("./dataProcessed/dataset1hz.RData")
lungeFile <- "./Figure_files/Figure1-BmSubset/Bm170622-TDR12 lunges.mat" 
lunges = readMat(lungeFile)
lungeDT <- round_date(as.POSIXct((lunges$LungeDN - 719529)*86400, origin = "1970-01-01", tz = 'GMT'),unit = "seconds")
lunge_dttz  <- force_tz(lungeDT,tzone="Etc/GMT+7") # Force TZ
lungeDF <- as.data.frame(lunge_dttz)
lungeDF$lunge_dt <- lungeDF$lunge_dttz
attr(lungeDF$lunge_dt, "tzone") <- 'GMT'
lungeDF$lungeDepth <- lunges$LungeDepth
lungeDF$lungeC <- lunges$LungeC
rm(lunge_dttz, lungeDT)
# Select only High-Confidence Lunges Detected
lungeDF <- lungeDF %>% 
  dplyr::filter(lungeC == 3) %>% 
  dplyr::select(-lungeC) 

diveSample <- as.data.frame(divesHMM4state) %>% 
  dplyr::filter(depid == "Bm170622-TDR12", DiveNum >= 2556, DiveNum <= 2573) 
sample1Hz <- dataset1hz %>% dplyr::filter(depid == "Bm170622-TDR12", 
                                          dttz >= min(diveSample$dttzStart), 
                                          dttz <= (max(diveSample$dttzStart)))

lungeSample <- lungeDF %>% dplyr::filter(lunge_dttz >= min(diveSample$dttzStart), 
                                        lunge_dttz <= (max(diveSample$dttzStart)+last(diveSample$dive_Duration)+last(diveSample$surf_Duration)))

ggplot() +
 geom_rect(data=diveSample %>% slice(-n()),aes(xmin=dttzStart,xmax=(dttzStart+(dive_Duration*60)),ymin=max(sample1Hz$p + 5),ymax=-5, 
                                               fill=as.factor(stateFeedingRate)),alpha=.25, show.legend = TRUE)+
 scale_fill_manual(name=NULL,
                    values = c("1"="firebrick3",
                               "2"="forestgreen",
                               "3"="darkorchid4",
                               "4"="royalblue3"),
                    labels = c("Not Feeding",
                               "Light Feeding",
                               "Moderate Feeding",
                               "Heavy Feeding")) + 
  ylab("Depth (m)") +xlab("Time") +
  geom_line(data=sample1Hz, aes(x=dttz,y=p)) +
  geom_point(data= lungeSample, aes(x=lunge_dttz,y=lungeDepth), color='red', size = 1.25)+
  scale_y_reverse(expand=c(0,0))+
  scale_x_datetime(expand=c(0,0))+
  theme_classic()+
  theme(axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18),
        axis.text.x = element_text(face="bold", size = 18),
        plot.title = element_text(face="bold", size=24,hjust = 0.5),
        plot.subtitle = element_text( face="bold", size=20,hjust = 0.5),
        plot.caption = element_text( face="bold", size=16,hjust = 0),
        legend.text = element_text(face="bold", size=18),
        legend.title = element_text(face="bold", size=22,hjust = 0.5),
        legend.margin = margin(0.15,0.15,0.15,0.15, unit="cm"),
        legend.background = element_rect(color=NA),
        legend.key = element_rect(fill=NA,color=NA),
        legend.key.size = unit(1,"cm"),
        legend.position = "top")
ggsave("Output/Figure1B.png", width=12, height=4, units = "in", dpi=400)

#### Blown up Dives ####
# 8Hz prh exported from matlab
prh8hzData <- readMat("./Figure_files/Figure1-BmSubset/Bm170622-TDR12 8HzprhSubset.mat") 
prh_DF <- data.frame(p=prh8hzData$p,DN=prh8hzData$DN,pitch=prh8hzData$pitch,roll=prh8hzData$roll,speed = prh8hzData$speedJJ)

prh_DF$dttz <- as.POSIXct((prh_DF$DN - 719529)*86400, origin = "1970-01-01", tz = 'GMT')
prh_DF$dttz  <- force_tz(prh_DF$dttz,tzone="Etc/GMT+7") # Force TZ

## Non-Feeding Dive
prhSubNF <- prh_DF %>% 
  dplyr::filter(dttz>=diveSample$dttzStart[1],dttz<=diveSample$dttzStart[2])
plot(prhSubNF$DN,prhSubNF$p)

pNF <- ggplot()+
  geom_line(data=prhSubNF,aes(x=dttz,y=p),size=1.25) +
  scale_y_reverse(expand=c(0,0), limits=c(220,0))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Depth\n(m)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(1,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
pNF
pitchNF <- ggplot()+
  geom_line(data=prhSubNF,aes(x=dttz,y=pitch*(180/pi)),size=1.25) +
  scale_y_continuous(expand=c(0,0), limits = c(-90,90))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Pitch\n(deg)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(0,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
rollNF <- ggplot()+
  geom_line(data=prhSubNF,aes(x=dttz,y=roll*(180/pi)),size=1.25) +
  scale_y_continuous(expand=c(0,0), limits = c(-90,90))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Roll\n(deg)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(0,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
speedNF <- ggplot()+
  geom_line(data=prhSubNF,aes(x=dttz,y=speedFiltered),size=1.25) +
  scale_y_continuous(expand=c(0,0),limits=c(0,4.5))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Speed\n(ms^-2)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(0,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18))

ggarrange(pNF,pitchNF,rollNF,speedNF, nrow=4,heights = c(1.5,1,1,1),align="v")
ggsave("./Output/Figure1NotFeeding.png", width=16, height=10, units = "in", dpi=400)

## Feeding Dive
prhSubF <- prh_DF %>% 
  dplyr::filter(dttz>=diveSample$dttzStart[13],dttz<=diveSample$dttzStart[14])
plot(prhSubF$DN,prhSubF$p)

lungeSampleSub <- lungeSample %>% 
  dplyr::filter(lunge_dttz>=min(prhSubF$dttz),
                lunge_dttz<=max(prhSubF$dttz)) 
lungeSampleSub$prhI <- c(882, 1807,2620,3463,4252,5085)
lungeSampleSub$pitch <- prhSubF$pitch[ c(882, 1807,2620,3463,4252,5085)]
lungeSampleSub$roll <- prhSubF$roll[ c(882, 1807,2620,3463,4252,5085)]
lungeSampleSub$speed <- prhSubF$speedFiltered[ c(882, 1807,2620,3463,4252,5085)]

pF <- ggplot()+
  geom_line(data=prhSubF,aes(x=dttz,y=p),size=1.25) +
  geom_point(data= lungeSample %>% dplyr::filter(lunge_dttz>=min(prhSubF$dttz),
                                                 lunge_dttz<=max(prhSubF$dttz)), 
             aes(x=lunge_dttz,y=lungeDepth), color='red', size = 4)+
  scale_y_reverse(expand=c(0,0), limits=c(220,0))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Depth\n(m)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(1,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
pF
pitchF <- ggplot()+
  geom_line(data=prhSubF,aes(x=dttz,y=pitch*(180/pi)),size=1.25) +
  geom_point(data= lungeSampleSub,
             aes(x=lunge_dttz,y=pitch*(180/pi)), color='red', size = 4)+
  
  scale_y_continuous(expand=c(0,0), limits = c(-90,90))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Pitch\n(deg)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(0,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
pitchF
rollF <- ggplot()+
  geom_line(data=prhSubF,aes(x=dttz,y=roll*(180/pi)),size=1.25) +
  geom_point(data= lungeSampleSub,
             aes(x=lunge_dttz,y=roll*(180/pi)), color='red', size = 4)+
  scale_y_continuous(expand=c(0,0), limits = c(-90,90))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Roll\n(deg)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(0,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
rollF
speedF <- ggplot()+
  geom_line(data=prhSubF,aes(x=dttz,y=speedFiltered),size=1.25) +
  geom_point(data=lungeSampleSub, 
             aes(x=lunge_dttz,y=speed), color='red', size = 4)+
  scale_y_continuous(expand=c(0,0),limits=c(0,4.5))+
  scale_x_datetime(expand=c(0,0))+
  ylab("Speed\n(ms^-2)") +xlab("Time") +
  theme_bw()+
  theme(plot.margin=unit(c(0,1,0,0.25), "cm"), # adds a buffer on the sides so the plot is centered
        axis.title = element_text(face="bold", size=20),
        axis.text = element_text( face="bold", size=18))

ggarrange(pF,pitchF,rollF,speedF, nrow=4,heights = c(1.5,1,1,1),align="v")
ggsave("./Output/Figure1Feeding.png", width=16, height=10, units = "in", dpi=400)


