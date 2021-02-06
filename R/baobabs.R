#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

## Set environmental variable
## For MAXENT.Phillips with JAVA to work on RStudio server
Sys.unsetenv("DISPLAY")

## Libraries
pkg <- c("curl", "here", "readr", "dplyr", "rgdal", "sp", "raster",
         "biomod2", "foreach", "doParallel", "grDevices")
load.pkg <- function(x) {
  if(!require(x, character.only = T)) {
    install.packages(x)
    require(x, character.only = T)
  }
}
loaded <- lapply(pkg,load.pkg)
rm(pkg,load.pkg,loaded)

## Directory names (with trailing slash "/")
dir_var_sdm <- "data/gisdata/sdm_variables/"
dir_baobabs <- "data/baobabs/"

## Figure size
textwidth <- 16.6 # in cm

##=======================
## Environmental data
##=======================

## Rasters of environmental variables
if (!file.exists(paste0(dir_var_sdm,"environ.tif"))) {
  dir.create(dir_var_sdm, recursive=TRUE, showWarnings=FALSE)
  list.url <- c(
    ## Environment
    "https://madaclim.cirad.fr/environ/environ.tif",
    ## Current climate
    "https://madaclim.cirad.fr/climate/current.tif",
    ## GISS-E2-R
    "https://madaclim.cirad.fr/climate/gs_45_2050.tif",
    "https://madaclim.cirad.fr/climate/gs_45_2080.tif",
    "https://madaclim.cirad.fr/climate/gs_85_2050.tif",
    "https://madaclim.cirad.fr/climate/gs_85_2080.tif",
    ## HadGEM2-ES
    "https://madaclim.cirad.fr/climate/he_45_2050.tif",
    "https://madaclim.cirad.fr/climate/he_45_2080.tif",
    "https://madaclim.cirad.fr/climate/he_85_2050.tif",
    "https://madaclim.cirad.fr/climate/he_85_2080.tif",
    ## NorESM1-M
    "https://madaclim.cirad.fr/climate/no_45_2050.tif",
    "https://madaclim.cirad.fr/climate/no_45_2080.tif",
    "https://madaclim.cirad.fr/climate/no_85_2050.tif",
    "https://madaclim.cirad.fr/climate/no_85_2080.tif"
    
  )
  for (i in 1:length(list.url)) {
    cat(paste0("Downloading ",list.url[i]," \n"))
    dest.f <- file.path(dir_var_sdm,basename(list.url[i]))
    curl::curl_download(url=list.url[i],destfile=dest.f)
  }
}

current <- raster::stack(paste0(dir_var_sdm,"current.tif"))
names(current) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),
                    paste("prec",1:12,sep=""),paste("bio",1:19,sep=""),
                    paste("pet",1:12,sep=""),"pet","cwd","ndm")
environ <- raster::stack(paste0(dir_var_sdm,"environ.tif"))
names(environ) <- c("alt","slope","asp","solar","geol","soil","veg","wshed","percfor2010")
## Stack of explicative variables
wc <- which(names(current) %in% c("bio1","bio4","bio12","cwd"))
we <- which(names(environ) %in% c("solar","alt","soil","percfor2010"))
s <- stack(current[[wc]])
names(s) <- c("tmean","tseas","prec","cwd")

## Remove data for Comoro Islands
bbCom <- extent(xmin(s),600000,8500000,ymax(s)) # bounding-box
cellsCom <- cellsFromExtent(s,bbCom)
values(s)[cellsCom,] <- NA
values(environ$alt)[cellsCom] <- NA
s <- stack(s) # Transform back from RasterBrick to RasterStack

## Create output dir
dir.create("outputs",showWarnings=FALSE)

# ## Plot environmental variables
# if(!file.exists("outputs/environ.pdf")) {
#   pdf(file="outputs/environ.pdf")
#   plot(s)
#   dev.off()
# }

##===================================
## Occurrence data for baobab species
##===================================

## Building the occurrence dataset from raw data
## You can disregard the warnings
source(here("R/data_baobabs.R"))
## Load dataset
df.orig <- read.csv(file=paste0(dir_baobabs,"data_Adansonia.csv"),header=TRUE,sep=",")
## Make a SpatialPointsDataFrame object
coords <- cbind(df.orig$Long,df.orig$Lat)
## see PROJ.6 style here - https://www.gaia-gis.it/fossil/libspatialite/wiki?name=PROJ.6
df.sp <- SpatialPointsDataFrame(coords, data=df.orig, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +type=crs"))
## Reproject into UTM 38S - SET exactly as "s" raster created above
df.sp <- spTransform(df.sp,CRS("+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +type=crs"))
## Only for Baoabab data: change species names
df.sp$Species <- gsub("A_","Adansonia ",df.sp$Species)
# Species
sp.names <- levels(as.factor(df.sp$Species)) # Sorted in alphabetical order
sp.dir <- gsub(" ",".",sp.names)
n.species <- length(sp.names)

##==================================
## Computation per species
##==================================

## Load run_species() function
source(here("R/run_species.R"))
## Path to MaxEnt (indicate your path)
path_to_maxent.jar <- here("maxent")
## Make a cluster with all possible cores
n.core <- max(1,detectCores()-2)
clust <- makeCluster(n.core)
## Register the number of parallel workers (here all CPUs)
registerDoParallel(clust)
## Return number of parallel workers
getDoParWorkers() 
## Package names for parallel computations
pkg.names.clust <- c("rgdal","sp","raster","biomod2")
## Parallel computations
t.start <- Sys.time() ## Start the clock
foreach(i=1:n.species,.packages=pkg.names.clust) %dopar% run.species(i, path_to_maxent.jar, run.models=TRUE)
## Stop the cluster
stopCluster(clust)
## Time computation
t.stop <- Sys.time() ## Stop the clock
t.diff <- difftime(t.stop,t.start,units="min")
cat(paste0("Computation time (min): ",round(t.diff,2)),file="outputs/computation_time.txt")

##==================================
## Plot anomalies
##==================================
source(here("R/plot_anomalies.R"))

##==================================
## Plot climate change in SDA
##==================================
source(here("R/plot_climate_change_SDA.R"))

# ====================================================
# Plot species distribution shift under climate change 
# ====================================================
source(here("R/plot_SDA.R"))

# ==========================================================
# Draw points in the SDA and extract environmental variables
# ==========================================================

map_result <- stack(s,environ$alt)
as.v.na <- na.omit(map_result)

as.v.na2 <- rasterToPoints(as.v.na)
dados_tcc <- as.data.frame(as.v.na2[complete.cases(as.v.na2), ] )

set.seed(20)
data_frame2 <- sample_n(dados_tcc, size = 1000, replace = F)

# Plot MAP-MAT
range(data_frame2$alt)
range(data_frame2$tmean)

### Alt x Tmean
map.mat <- ggplot(data_frame2, aes(x=alt, y=tmean)) + xlim(0,2041) + ylim(137,280) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.5) + 
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Elevation (m)",y="Mean annual temp. (°C x 10)",size=4) +
  annotate("text", x  = 0, y = 280 , size=7, label = "(a)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
ggsave(file=paste0("./outputs/meantem_elevation.png"),
       plot=map.mat,width=10,height=10,dpi='print')

#### Altitude x T. seas
range(data_frame2$alt)
range(data_frame2$tseas)

map.mat2 <- ggplot(data_frame2, aes(x=alt, y=tseas)) + xlim(0,2041) + ylim(991,3500) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.4) +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Elevation (m)",y="Temp. seasonality (°C sd x 100)",size=4) +
  annotate("text", x  = 0, y = 3500 , size=7, label = "(b)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
ggsave(file=paste0("./outputs/meantempseas_elevation.png"),
       plot=map.mat2,width=10,height=10,dpi='print')

## Ploting comparing the latitude
## Tmean x Latitude

range(data_frame2$tmean)
range(data_frame2$y) #y latitude, don't forget ;) 7172500 8673500

lat.mat1 <- ggplot(data_frame2, aes(x=data_frame2$y, y=tmean)) + ylim(137,280) + xlim(7172500,8673500) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.4) +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Latitude (UTM)",y="Mean annual temp. (°C x 10)",size=4) +
  annotate("text", x  = 7200000, y = 280 , size=7, label = "(c)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
ggsave(file=paste0("./outputs/meantemp_latitude.png"),
       plot=lat.mat1,width=10,height=10,dpi='print')

## Tseas x Latitude

range(data_frame2$tseas)
range(data_frame2$y) #y latitude, don't forget ;) 7172500 8673500

lat.mat2 <- ggplot(data_frame2, aes(x=data_frame2$y, y=tseas)) + ylim(991,3500) + xlim(7172500,8673500) +
  #geom_density2d(data=Abs.df,col=grey(0.5)) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.4) +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Latitude (UTM)",y="Temp. seasonality (°C sd x 100)",size=4) +
  annotate("text", x  = 7200000, y = 3500 , size=7, label = "(d)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
ggsave(file=paste0("./outputs/meantemp.seas_latitude.png"),
       plot=lat.mat2,width=10,height=10,dpi='print')

## Ploting comparing the longitude
## Tmean x Latitude

#range(data_frame2$tmean)
#range(data_frame2$x) #y latitude, don't forget ;) 7172500 8673500

#lon.mat1 <- ggplot(data_frame2, aes(x=data_frame2$x, y=tmean)) + ylim(121,275) + xlim(318500,1087500) +
  #geom_density2d(data=Abs.df,col=grey(0.5)) +
 # geom_point(data=data_frame,col="darkgreen",alpha=0.4) +
  #geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  #labs(x="Longitude (UTM)",y="Mean annual temp. (°C x 10)",size=4) +
  #theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
   #     axis.title.x=element_text(size=rel(2)),
    #    axis.title.y=element_text(size=rel(2)),
     #   axis.text.x=element_text(size=rel(1.75)),
      #  axis.text.y=element_text(size=rel(1.75)))
#ggsave(file=paste0("./outputs/meantemp_longitude.png"),
  #     plot=lon.mat1,width=10,height=10,dpi='print')

## Tseas x Latitude

#range(data_frame2$tseas)
#range(data_frame2$x) #y latitude, don't forget ;) 7172500 8673500

#data_frame2$tseas
#lon.mat2 <- ggplot(data_frame2, aes(x=data_frame2$x, y=tseas)) + ylim(897,3303) + xlim(318500,1087500) +
  #geom_density2d(data=Abs.df,col=grey(0.5)) +
#  geom_point(data=data_frame,col="darkgreen",alpha=0.4) +
 # geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  #labs(x="Longitude (UTM)",y="Temp. seasonality (°C sd x 100)",size=4) +
  #theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
   #     axis.title.x=element_text(size=rel(2)),
    #    axis.title.y=element_text(size=rel(2)),
     #   axis.text.x=element_text(size=rel(1.75)),
      #  axis.text.y=element_text(size=rel(1.75)))
#ggsave(file=paste0("./outputs/meantemp.seas_longitude.png"),
 #      plot=lon.mat2,width=10,height=10,dpi='print')
library(grid)

a5 <- grid.arrange( map.mat, map.mat2, lat.mat1, lat.mat2, nrow=2, ncol=2)


ggsave(file=paste0("./outputs/four_chart_lat_elevations.png"),
       plot=a5,width=12,height=12,dpi="print")

####################### Box plots according to each SDAp and SDAf compared

##### Non threatened baobab species ####
### A digitata ###

pred_dig <- stack(paste0("Adansonia.digitata/proj_current/proj_current_Adansonia.digitata_ensemble.grd"))
ca_dig <- pred_dig[[1]]
ca_test <- stack(ca_dig,environ$alt)
ca_test_na <- na.omit(ca_test)

ca_test_na2 <- rasterToPoints(ca_test_na)
ca_test_na3 <- as.data.frame(ca_test_na2[complete.cases(ca_test_na2), ] )

alegria <- ca_test_na3 %>% filter(Adansonia.digitata_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

####################### Box plots according to each SDAp and SDAf compared

####################### Box plots according to each SDAp and SDAf compared

##### Non threatened baobab species ####
### A digitata ###

pred_dig <- stack(paste0("Adansonia.digitata/proj_current/proj_current_Adansonia.digitata_ensemble.grd"))
ca_dig <- pred_dig[[1]]
ca_test <- stack(ca_dig,environ$alt)
ca_test_na <- na.omit(ca_test)

ca_test_na2 <- rasterToPoints(ca_test_na)
ca_test_na3 <- as.data.frame(ca_test_na2[complete.cases(ca_test_na2), ] )

alegria <- ca_test_na3 %>% filter(Adansonia.digitata_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

#set.seed(20)
alegria2 <- sample_n(alegria, size = 1000, replace = F)
alegria2$Proj <- rep(c("Present"),1000)
names(alegria2) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_dig <- raster(paste0("Adansonia.digitata/caFut_85_2080.tif"))
ca_test_fut <- stack(caZD_dig,environ$alt)
ca_test_na_fut <- na.omit(ca_test_fut)

ca_test_na2_fut <- rasterToPoints(ca_test_na_fut)
ca_test_na3_fut <- as.data.frame(ca_test_na2_fut[complete.cases(ca_test_na2_fut), ] )
alegria_fut <- ca_test_na3_fut %>% filter(caFut_85_2080 >= 1500)

#set.seed(20)
alegria2_fut <- sample_n(alegria_fut, size = 1000, replace = F)
alegria2_fut$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test <- rbind(alegria2,alegria2_fut)
tail(alegria_test)

### Future testing 2050 
caZD_digi_2055 <- raster(paste0("Adansonia.digitata/caFut_85_2050.tif"))
ca_test_fut_digi_2055 <- stack(caZD_digi_2055,environ$alt)
ca_test_na_fut_digi_55 <- na.omit(ca_test_fut_digi_2055)

ca_test_na2_fut_digi_55 <- rasterToPoints(ca_test_na_fut_digi_55)
ca_test_na3_fut_digi_55 <- as.data.frame(ca_test_na2_fut_digi_55[complete.cases(ca_test_na2_fut_digi_55), ] )
alegria_fut_digi_55 <- ca_test_na3_fut_digi_55 %>% filter(caFut_85_2050 >= 1500)

#set.seed(20)
alegria2_fut_digi_55 <- sample_n(alegria_fut_digi_55, size = 1000, replace = F)
alegria2_fut_digi_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_digi_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_digi_finale <- rbind(alegria_test, alegria2_fut_digi_55)


alegria_test_digi_finale$Scenario = factor(alegria_test_digi_finale$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot

latitude_digi <- ggplot(alegria_test_digi_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8650000 , size=7, label = "(b)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

## Altitude plot
altitude_digi <- ggplot(alegria_test_digi_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 800 , size=7, label = "(a)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

########
######## A grandidieri ##############

pred_grand <- stack(paste0("Adansonia.grandidieri/proj_current/proj_current_Adansonia.grandidieri_ensemble.grd"))
ca_grand <- pred_grand[[1]]
ca_test_grand <- stack(ca_grand,environ$alt)
ca_test_na_grand <- na.omit(ca_test_grand)

ca_test_na2_grand <- rasterToPoints(ca_test_na_grand)
ca_test_na3_grand <- as.data.frame(ca_test_na2_grand[complete.cases(ca_test_na2_grand), ] )

alegria_grand <- ca_test_na3_grand %>% filter(Adansonia.grandidieri_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

#set.seed(20)
alegria2_grand <- sample_n(alegria_grand, size = 1000, replace = F)
alegria2_grand$Proj <- rep(c("Present"),1000)
names(alegria2_grand) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_grand <- raster(paste0("Adansonia.grandidieri/caFut_85_2080.tif"))
ca_test_fut_grand <- stack(caZD_grand,environ$alt)
ca_test_na_fut_grand <- na.omit(ca_test_fut_grand)

ca_test_na2_fut_grand <- rasterToPoints(ca_test_na_fut_grand)
ca_test_na3_fut_grand <- as.data.frame(ca_test_na2_fut_grand[complete.cases(ca_test_na2_fut_grand), ] )
alegria_fut_grand <- ca_test_na3_fut_grand %>% filter(caFut_85_2080 >= 1500)

#set.seed(20)
alegria2_fut_grand <- sample_n(alegria_fut_grand, size = 1000, replace = F)
alegria2_fut_grand$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_grand) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_grand <- rbind(alegria2_grand,alegria2_fut_grand)


### Future testing 2050 
caZD_grand_2055 <- raster(paste0("Adansonia.grandidieri/caFut_85_2050.tif"))
ca_test_fut_grand_2055 <- stack(caZD_grand_2055,environ$alt)
ca_test_na_fut_grand_55 <- na.omit(ca_test_fut_grand_2055)

ca_test_na2_fut_grand_55 <- rasterToPoints(ca_test_na_fut_grand_55)
ca_test_na3_fut_grand_55 <- as.data.frame(ca_test_na2_fut_grand_55[complete.cases(ca_test_na2_fut_grand_55), ] )
alegria_fut_grand_55 <- ca_test_na3_fut_grand_55 %>% filter(caFut_85_2050 >= 1500)


#set.seed(20)
alegria2_fut_grand_55 <- sample_n(alegria_fut_grand_55, size = 1000, replace = F)
alegria2_fut_grand_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_grand_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_grand_finale <- rbind(alegria2_fut_grand_55,alegria_test_grand)

alegria_test_grand_finale$Scenario = factor(alegria_test_grand_finale$Scenario,
                                     levels = c("Present", "Future_2055","Future_2085"),
                                     labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_grand <- ggplot(alegria_test_grand_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8250000 , size=7, label = "(d)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

## Altitude plot
altitude_grand <- ggplot(alegria_test_grand_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1050 , size=7, label = "(c)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

######## A za ###############

pred_za <- stack(paste0("Adansonia.za/proj_current/proj_current_Adansonia.za_ensemble.grd"))
ca_za <- pred_za[[1]]
ca_test_za <- stack(ca_za,environ$alt)
ca_test_na_za <- na.omit(ca_test_za)

ca_test_na2_za <- rasterToPoints(ca_test_na_za)
ca_test_na3_za <- as.data.frame(ca_test_na2_za[complete.cases(ca_test_na2_za), ] )

alegria_za <- ca_test_na3_za %>% filter(Adansonia.za_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

#set.seed(20)
alegria2_za <- sample_n(alegria_za, size = 1000, replace = F)
alegria2_za$Proj <- rep(c("Present"),1000)
names(alegria2_za) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_za <- raster(paste0("Adansonia.za/caFut_85_2080.tif"))
ca_test_fut_za <- stack(caZD_za,environ$alt)
ca_test_na_fut_za <- na.omit(ca_test_fut_za)

ca_test_na2_fut_za <- rasterToPoints(ca_test_na_fut_za)
ca_test_na3_fut_za <- as.data.frame(ca_test_na2_fut_za[complete.cases(ca_test_na2_fut_za), ] )
alegria_fut_za <- ca_test_na3_fut_za %>% filter(caFut_85_2080 >= 1500)

#set.seed(20)
alegria2_fut_za <- sample_n(alegria_fut_za, size = 1000, replace = F)
alegria2_fut_za$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_za) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_za <- rbind(alegria2_za,alegria2_fut_za)

### Future testing 2050 
caZD_za_2055 <- raster(paste0("Adansonia.za/caFut_85_2050.tif"))
ca_test_fut_za_2055 <- stack(caZD_za_2055,environ$alt)
ca_test_na_fut_za_55 <- na.omit(ca_test_fut_za_2055)

ca_test_na2_fut_za_55 <- rasterToPoints(ca_test_na_fut_za_55)
ca_test_na3_fut_za_55 <- as.data.frame(ca_test_na2_fut_za_55[complete.cases(ca_test_na2_fut_za_55), ] )
alegria_fut_za_55 <- ca_test_na3_fut_za_55 %>% filter(caFut_85_2050 >= 1500)

#set.seed(20)
alegria2_fut_za_55 <- sample_n(alegria_fut_za_55, size = 1000, replace = F)
alegria2_fut_za_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_za_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_za_finale <- rbind(alegria2_fut_za_55,alegria_test_za)

alegria_test_za_finale$Scenario = factor(alegria_test_za_finale$Scenario,
                                         levels = c("Present", "Future_2055","Future_2085"),
                                         labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_za <- ggplot(alegria_test_za_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8650000 , size=7, label = "(f)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

## Altitude plot
altitude_za <- ggplot(alegria_test_za_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1350 , size=7, label = "(e)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

###########Threatened baobab species ###############
######## A madagascariensis ############

pred_mada <- stack(paste0("Adansonia.madagascariensis/proj_current/proj_current_Adansonia.madagascariensis_ensemble.grd"))
ca_mada <- pred_mada[[1]]
ca_test_mada <- stack(ca_mada,environ$alt)
ca_test_na_mada <- na.omit(ca_test_mada)

ca_test_na2_mada <- rasterToPoints(ca_test_na_mada)
ca_test_na3_mada <- as.data.frame(ca_test_na2_mada[complete.cases(ca_test_na2_mada), ] )

alegria_mada <- ca_test_na3_mada %>% filter(Adansonia.madagascariensis_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

#set.seed(20)
alegria2_mada <- sample_n(alegria_mada, size = 1000, replace = F)
alegria2_mada$Proj <- rep(c("Present"),1000)
names(alegria2_mada) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caFut_mada <- raster(paste0("Adansonia.madagascariensis/caFut_85_2080.tif"))
ca_test_fut_mada <- stack(caFut_mada,environ$alt)
ca_test_na_fut_mada <- na.omit(ca_test_fut_mada)

ca_test_na2_fut_mada <- rasterToPoints(ca_test_na_fut_mada)
ca_test_na3_fut_mada <- as.data.frame(ca_test_na2_fut_mada[complete.cases(ca_test_na2_fut_mada), ] )
alegria_fut_mada <- ca_test_na3_fut_mada %>% filter(caFut_85_2080 >= 1500)

#set.seed(20)
alegria2_fut_mada <- sample_n(alegria_fut_mada, size = 1000, replace = F)
alegria2_fut_mada$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_mada) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_mada <- rbind(alegria2_mada,alegria2_fut_mada)
tail(alegria_test_mada)

### Future testing 2050 
caFut_mada_2055 <- raster(paste0("Adansonia.madagascariensis/caFut_85_2050.tif"))
ca_test_fut_mada_2055 <- stack(caFut_mada_2055,environ$alt)
ca_test_na_fut_mada_55 <- na.omit(ca_test_fut_mada_2055)

ca_test_na2_fut_mada_55 <- rasterToPoints(ca_test_na_fut_mada_55)
ca_test_na3_fut_mada_55 <- as.data.frame(ca_test_na2_fut_mada_55[complete.cases(ca_test_na2_fut_mada_55), ] )
alegria_fut_mada_55 <- ca_test_na3_fut_mada_55 %>% filter(caFut_85_2050 >= 1500)


#set.seed(20)
alegria2_fut_mada_55 <- sample_n(alegria_fut_mada_55, size = 1000, replace = F)
alegria2_fut_mada_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_mada_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_mada_finale <- rbind(alegria2_fut_mada_55,alegria_test_mada)

alegria_test_mada_finale$Scenario = factor(alegria_test_mada_finale$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_mada <- ggplot(alegria_test_mada_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8650000 , size=7, label = "(b)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +    
  theme(legend.position = "none")


## Altitude plot

altitude_mada <- ggplot(alegria_test_mada_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 800 , size=7, label = "(a)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

######## A perrieri ###############

pred_perrieri <- stack(paste0("Adansonia.perrieri/proj_current/proj_current_Adansonia.perrieri_ensemble.grd"))
ca_perrieri <- pred_perrieri[[1]]
ca_test_perrieri <- stack(ca_perrieri,environ$alt)
ca_test_na_perrieri <- na.omit(ca_test_perrieri)

ca_test_na2_perrieri <- rasterToPoints(ca_test_na_perrieri)
ca_test_na3_perrieri <- as.data.frame(ca_test_na2_perrieri[complete.cases(ca_test_na2_perrieri), ] )

alegria_perrieri <- ca_test_na3_perrieri %>% filter(Adansonia.perrieri_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

#set.seed(20)
alegria2_perrieri <- sample_n(alegria_perrieri, size = 1000, replace = F)
alegria2_perrieri$Proj <- rep(c("Present"),1000)
names(alegria2_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 (RCP 4.5!!!) Atention!
caFut_perrieri <- raster(paste0("Adansonia.perrieri/caFut_45_2080.tif"))
ca_test_fut_perrieri <- stack(caFut_perrieri,environ$alt)
ca_test_na_fut_perrieri <- na.omit(ca_test_fut_perrieri)

ca_test_na2_fut_perrieri <- rasterToPoints(ca_test_na_fut_perrieri)
ca_test_na3_fut_perrieri <- as.data.frame(ca_test_na2_fut_perrieri[complete.cases(ca_test_na2_fut_perrieri), ] )
alegria_fut_perrieri <- ca_test_na3_fut_perrieri %>% filter(caFut_45_2080 >= 1500)

#set.seed(20)
alegria2_fut_perrieri <- sample_n(alegria_fut_perrieri, size = 416, replace = F)
alegria2_fut_perrieri$Proj <- rep(c("Future_2085"),416)
names(alegria2_fut_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_perrieri <- rbind(alegria2_perrieri,alegria2_fut_perrieri)
tail(alegria_test_perrieri)

### Future testing 2050 
caFut_perrieri_2055 <- raster(paste0("Adansonia.perrieri/caFut_85_2050.tif"))
ca_test_fut_perrieri_2055 <- stack(caFut_perrieri_2055,environ$alt)
ca_test_na_fut_perrieri_55 <- na.omit(ca_test_fut_perrieri_2055)

ca_test_na2_fut_perrieri_55 <- rasterToPoints(ca_test_na_fut_perrieri_55)
ca_test_na3_fut_perrieri_55 <- as.data.frame(ca_test_na2_fut_perrieri_55[complete.cases(ca_test_na2_fut_perrieri_55), ] )
alegria_fut_perrieri_55 <- ca_test_na3_fut_perrieri_55 %>% filter(caFut_85_2050 >= 1500)


#set.seed(20)
alegria2_fut_perrieri_55 <- sample_n(alegria_fut_perrieri_55, size = 1000, replace = F)
alegria2_fut_perrieri_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_perrieri_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_perrieri_finale <- rbind(alegria2_fut_perrieri_55,alegria_test_perrieri)

alegria_test_perrieri_finale$Scenario = factor(alegria_test_perrieri_finale$Scenario,
                                           levels = c("Present", "Future_2055","Future_2085"),
                                           labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_perrieri <- ggplot(alegria_test_perrieri_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8600000 , size=7, label = "(d)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +    
  theme(legend.position = "none")


## Altitude plot
altitude_perrieri <- ggplot(alegria_test_perrieri_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1450 , size=7, label = "(c)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

######## A rubrostipa ###############
pred_rubro <- stack(paste0("Adansonia.rubrostipa/proj_current/proj_current_Adansonia.rubrostipa_ensemble.grd"))
ca_rubro <- pred_rubro[[1]]
ca_test_rubro <- stack(ca_rubro,environ$alt)
ca_test_na_rubro <- na.omit(ca_test_rubro)

ca_test_na2_rubro <- rasterToPoints(ca_test_na_rubro)
ca_test_na3_rubro <- as.data.frame(ca_test_na2_rubro[complete.cases(ca_test_na2_rubro), ] )

alegria_rubro <- ca_test_na3_rubro %>% filter(Adansonia.rubrostipa_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

#set.seed(20)
alegria2_rubro <- sample_n(alegria_rubro, size = 1000, replace = F)
alegria2_rubro$Proj <- rep(c("Present"),1000)
names(alegria2_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caFut_rubro <- raster(paste0("Adansonia.rubrostipa/caFut_85_2080.tif"))
ca_test_fut_rubro <- stack(caFut_rubro,environ$alt)
ca_test_na_fut_rubro <- na.omit(ca_test_fut_rubro)

ca_test_na2_fut_rubro <- rasterToPoints(ca_test_na_fut_rubro)
ca_test_na3_fut_rubro <- as.data.frame(ca_test_na2_fut_rubro[complete.cases(ca_test_na2_fut_rubro), ] )
alegria_fut_rubro <- ca_test_na3_fut_rubro %>% filter(caFut_85_2080 >= 1500)

#set.seed(20)
alegria2_fut_rubro <- sample_n(alegria_fut_rubro, size = 1000, replace = F)
alegria2_fut_rubro$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_rubro <- rbind(alegria2_rubro,alegria2_fut_rubro)

### Future testing 2050 
caFut_rubro_2055 <- raster(paste0("Adansonia.rubrostipa/caFut_85_2050.tif"))
ca_test_fut_rubro_2055 <- stack(caFut_rubro_2055,environ$alt)
ca_test_na_fut_rubro_55 <- na.omit(ca_test_fut_rubro_2055)

ca_test_na2_fut_rubro_55 <- rasterToPoints(ca_test_na_fut_rubro_55)
ca_test_na3_fut_rubro_55 <- as.data.frame(ca_test_na2_fut_rubro_55[complete.cases(ca_test_na2_fut_rubro_55), ] )
alegria_fut_rubro_55 <- ca_test_na3_fut_rubro_55 %>% filter(caFut_85_2050 >= 1500)


#set.seed(20)
alegria2_fut_rubro_55 <- sample_n(alegria_fut_rubro_55, size = 1000, replace = F)
alegria2_fut_rubro_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_rubro_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_rubro_finale <- rbind(alegria2_fut_rubro_55,alegria_test_rubro)

alegria_test_rubro_finale$Scenario = factor(alegria_test_rubro_finale$Scenario,
                                           levels = c("Present", "Future_2055","Future_2085"),
                                           labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_rubro <- ggplot(alegria_test_rubro_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8600000 , size=7, label = "(f)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")


## Altitude plot
altitude_rubro <- ggplot(alegria_test_rubro_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1000 , size=7, label = "(e)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

######## A suarezensis ###############

pred_suare <- stack(paste0("Adansonia.suarezensis/proj_current/proj_current_Adansonia.suarezensis_ensemble.grd"))
ca_suare <- pred_suare[[1]]
ca_test_suare <- stack(ca_suare,environ$alt)
ca_test_na_suare <- na.omit(ca_test_suare)

ca_test_na2_suare <- rasterToPoints(ca_test_na_suare)
ca_test_na3_suare <- as.data.frame(ca_test_na2_suare[complete.cases(ca_test_na2_suare), ] )

alegria_suare <- ca_test_na3_suare %>% filter(Adansonia.suarezensis_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

#set.seed(20)
alegria2_suare <- sample_n(alegria_suare, size = 1000, replace = F)
alegria2_suare$Proj <- rep(c("Present"),1000)
names(alegria2_suare) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 Atention RCP 4.5 (SAME AS A. perrieri)
caFut_suare <- raster(paste0("Adansonia.suarezensis/caFut_45_2080.tif"))
ca_test_fut_suare <- stack(caFut_suare,environ$alt)
ca_test_na_fut_suare <- na.omit(ca_test_fut_suare)

ca_test_na2_fut_suare <- rasterToPoints(ca_test_na_fut_suare)
ca_test_na3_fut_suare <- as.data.frame(ca_test_na2_fut_suare[complete.cases(ca_test_na2_fut_suare), ] )
alegria_fut_suare <- ca_test_na3_fut_suare %>% filter(caFut_45_2080 >= 1500)

#set.seed(20)
alegria2_fut_suare <- sample_n(alegria_fut_suare, size = 105, replace = F)
alegria2_fut_suare$Proj <- rep(c("Future_2085"),105)
names(alegria2_fut_suare) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_suare <- rbind(alegria2_suare,alegria2_fut_suare)

### Future testing 2050 
caFut_suare_2055 <- raster(paste0("Adansonia.suarezensis/caFut_45_2050.tif"))
ca_test_fut_suare_2055 <- stack(caFut_suare_2055,environ$alt)
ca_test_na_fut_suare_55 <- na.omit(ca_test_fut_suare_2055)

ca_test_na2_fut_suare_55 <- rasterToPoints(ca_test_na_fut_suare_55)
ca_test_na3_fut_suare_55 <- as.data.frame(ca_test_na2_fut_suare_55[complete.cases(ca_test_na2_fut_suare_55), ] )
alegria_fut_suare_55 <- ca_test_na3_fut_suare_55 %>% filter(caFut_45_2050 >= 1500)

#set.seed(20)
alegria2_fut_suare_55 <- sample_n(alegria_fut_suare_55, size = 15, replace = F)
alegria2_fut_suare_55$Proj <- rep(c("Future_2055"),15)
names(alegria2_fut_suare_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_suare_finale <- rbind(alegria2_fut_suare_55,alegria_test_suare)

alegria_test_suare_finale$Scenario = factor(alegria_test_suare_finale$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_suare <- ggplot(alegria_test_suare_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8670000 , size=7, label = "(h)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +    
  theme(legend.position = "none")


## Altitude plot
altitude_suare <- ggplot(alegria_test_suare_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 950 , size=7, label = "(g)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")



## Legend names for plot
grob_digitata <- textGrob("A. digitata",rot=90, gp=gpar(cex=2,fontface="italic"),
                          hjust=0.3, vjust=2)
grob_grand <- textGrob("A. grandidieri",rot=90, gp=gpar(cex=2,fontface="italic"),
                          hjust=0.3, vjust=2)
grob_mada <- textGrob("A. madagascariensis",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.4, vjust=2)
grob_perri <- textGrob("A. perrieri",rot=90, gp=gpar(cex=2,fontface="italic"),
                     hjust=0.3, vjust=2)
grob_rubro <- textGrob("A. rubrostipa",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=2)
grob_suare <- textGrob("A. suarezensis",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=2)
grob_za <- textGrob("A. za",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=2)
## Combine plots
### Threatened species
day_1 <- rbind(c(1,rep(seq(5,6,by=1),each=5)),
              c(2,rep(seq(7,8,by=1),each=5)),
              c(3,rep(seq(9,10,by=1),each=5)),
              c(4,rep(seq(11,12,by=1),each=5)))

### Non-threatened species
day_2 <- rbind(c(1,rep(seq(5,6,by=1),each=5)),
               c(2,rep(seq(7,8,by=1),each=5)),
               c(3,rep(seq(9,10,by=1),each=5)))             

#### Threatened species
a5 <- grid.arrange(grob_mada, grob_perri,
                   grob_rubro,grob_suare,
                   altitude_mada,latitude_mada,altitude_perrieri,latitude_perrieri,
                   altitude_rubro,latitude_rubro,altitude_suare,latitude_suare,
                   layout_matrix = day_1)
#### Non-threatened species

a6 <- grid.arrange(grob_digitata, grob_grand,
                   grob_za,
                   altitude_digi,latitude_digi,altitude_grand,latitude_grand,
                   altitude_za,latitude_za,layout_matrix = day_2)

ggsave(file=paste0("./outputs/lat_ele_all_species_within_SDAcf_threat.png"),
       plot=a5,width=18,height=15,dpi="print")

ggsave(file=paste0("./outputs/lat_ele_all_species_within_SDAcf_non_threat.png"),
       plot=a6,width=18,height=15,dpi="print")

################### Alternative plot 

a_dig_lat_ele <- ggplot(alegria2, aes(x=Lat, y=Alt)) + ylim(-100,1000) + xlim(7170500,8664500) +
  #geom_density2d(data=Abs.df,col=grey(0.5)) +
  geom_density_2d(data=alegria2,col="blue",alpha=0.4) +
  geom_point() +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Latitude (UTM)",y="Elevation present (m) ",size=4) +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))

p1 <- ggMarginal(a_dig_lat_ele,type="boxplot")

ggsave(file=paste0("./outputs/meantemp.seas_latitude.png"),
       plot=lat.mat2,width=10,height=10,dpi='print')

######################################################################
### World Seasonality Map
#######################################################################

library(maptools)
r <- raster::getData("worldclim",var="bio",res=10)
r <- r[[c(4)]]
names(r) <- c("Seas")
bio4.he.2080_world <- getData('CMIP5', var='bio', res=10, rcp=85, model='HE', year=70)
bio4.he.2080_world <- bio4.he.2080_world[[c(4)]]

bio4.gs.2080_world <- getData('CMIP5', var='bio', res=10, rcp=85, model='GS', year=70)
bio4.gs.2080_world <- bio4.gs.2080_world[[c(4)]]

bio4.no.2080_world <- getData('CMIP5', var='bio', res=10, rcp=85, model='NO', year=70)
bio4.no.2080_world <- bio4.no.2080_world[[c(4)]]

Stack.bio4.2080_world <- stack(c(bio4.gs.2080_world,bio4.he.2080_world,
                                 bio4.no.2080_world))

bio4.2080_world <- mean(Stack.bio4.2080_world) ### 2080 temp seas. raster

# set the tropical range current and future 
## future
wp <-extent(c(-180, 180, -23.5, 23.5))
future_ws <- crop(bio4.2080_world,wp)
crs(future_ws) <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## current
Stack.bio4.current_worldcrop <- crop(r,wp)
crs(Stack.bio4.current_worldcrop) <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

current_ws <- Stack.bio4.current_worldcrop
Stack.bio4.current_worldcrop[] <- future_ws[]-current_ws[]

## Setting basic theme options for plot with ggplot2
theme_base_final <- theme(
  ## Axis
  axis.line=element_blank(),
  axis.text=element_blank(),
  axis.ticks=element_blank(),
  axis.title=element_blank(),
  ## Legend
  legend.position="bottom",
  legend.title=element_blank(),
  legend.text=element_text(size=12),
  legend.key.height=unit(0.5,"line"),
  legend.key.width=unit(2.5,"line"),
  legend.box.background=element_blank(),
  ## Plot
  plot.title=element_text(hjust=0.5,size=14),
  plot.background=element_rect(fill="transparent"),
  ## Panel
  panel.background=element_rect(fill="transparent"),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  panel.border=element_blank(),
  ## Margins
  plot.margin=unit(c(0,0,0,0),units="line"),
  legend.margin=margin(c(0,0,0,0),unit="line"),
  legend.box.margin=margin(c(0,0,0,0)),
  legend.box.spacing=unit(0,units="line")
)

## Function to plot anomalies
plot_anomaly_ws <- function(r, title, label="x") {
  rdf <- data.frame(rasterToPoints(r)) # To plot raster with geom_raster()
  names(rdf) <- c("x", "y", "z")
  p <- ggplot(NULL, aes(x, y)) +
    geom_tile(data=rdf, aes(fill=z)) +
    theme_bw() + theme_base_final +
    geom_hline(yintercept = 0,linetype = "dashed") +
    coord_fixed(xlim=c(-130,180), ylim=c(-23.5,23.5)) +
    annotate("text",x=-100,y=5,label=label,hjust=1,vjust=0,size=5) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    labs(title=title)
  return(p)
}

## Rescale function for legend color
rescale <- function(x,from.min,from.max,to.min=0,to.max=1) {
  a <- from.min; b <- from.max; c <- to.min; d <- to.max
  int <- (b*c-a*d)/(b-a) ; slope <- (d-c)/(b-a)
  return(int+slope*x)
}

## Color scales
## current world seasonality
col_scale_var_cur_ws <- scale_fill_gradientn(
  colours=viridis(255, option="B", direction=1),
  na.value="transparent",
  values=rescale(seq(0,7600,l=255),0,7600),
  limits=c(0,7600),
  breaks=seq(0,7600,l=6),
  labels=seq(0,7600,l=6)
)


## Color scales
## future world seasonality
col_scale_var_fut_ws <- scale_fill_gradientn(
  colours=viridis(255, option="B", direction=1),
  na.value="transparent",
  values=rescale(seq(100,8200,l=255),100,8200),
  limits=c(0,8200),
  breaks=seq(0,8200,l=6),
  labels=seq(0,8200,l=6)
)

## Color scales
## anomaly world seasonality
col_scale_var_anomws <- scale_fill_gradientn(
  colours=c(viridis(255, option="C", direction=1),rep(grey(0.5),3)),
  na.value="transparent",
  values=rescale(c(seq(0.01,1000,l=255),-700,0),-700,1000),
  limits=c(-700,1000),
  breaks=c(seq(0,1000,l=3),-700, -350),
  labels=c(seq(0,1000,l=3),-700, -350)
)

# Current ws
cur_ws <- plot_anomaly_ws(r=current_ws, label="(a)",
                      title="Current temperature seasonality (°C sd x 100)") + col_scale_var_cur_ws
                      
# Future anomaly
ano_ws <- plot_anomaly_ws(r=Stack.bio4.current_worldcrop, label="(b)",
                          title="Temperature seasonality anomaly (current vs. 2085 RCP 8.5, °C sd x 100)") + col_scale_var_anomws

# future seasonality
fut_ws <- plot_anomaly_ws(r=future_ws, label="(c)",
                          title="Future temperature seasonality (2085 RCP 8.5, °C sd x 100)") + col_scale_var_fut_ws

## Combine plots
lay_6 <- rbind(c(rep(seq(1,1,by=1),each=3)),
               c(rep(seq(1,1,by=1),each=3)),
               c(rep(seq(1,1,by=1),each=3)),
               c(rep(seq(2,2,by=1),each=3)),
               c(rep(seq(2,2,by=1),each=3)),
               c(rep(seq(2,2,by=1),each=3)),
               c(rep(seq(3,3,by=1),each=3)),
               c(rep(seq(3,3,by=1),each=3)),
               c(rep(seq(3,3,by=1),each=3)))

plot_world_anomaly <- grid.arrange(cur_ws, ano_ws, fut_ws,layout_matrix=lay_6)

ggsave(file=paste0("./outputs/anomaly_world_chart.pdf"),
       plot=plot_world_anomaly,width=11,height=8,dpi="print")

ggsave(file=paste0("./outputs/anomaly_world_chart.png"),
       plot=plot_world_anomaly,width=11,height=8,dpi="print")

##### Generating study tables

## area change
digitata <- read.table(paste0("Adansonia.digitata/sda_fut.txt"), header=T,sep="\t")
grandidieri <- read.table(paste0("Adansonia.grandidieri/sda_fut.txt"), header=T,sep="\t")
mada <- read.table(paste0("Adansonia.madagascariensis/sda_fut.txt"), header=T,sep="\t")
perrieri <- read.table(paste0("Adansonia.perrieri/sda_fut.txt"), header=T,sep="\t")
rubrostipa <- read.table(paste0("Adansonia.rubrostipa/sda_fut.txt"), header=T,sep="\t")
suare <- read.table(paste0("Adansonia.suarezensis/sda_fut.txt"), header=T,sep="\t")
za <- read.table(paste0("Adansonia.za/sda_fut.txt"), header=T,sep="\t")

species <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                 "A. perrieri", "A. rubrostipa","A. suarezensis",
                 "A. za"),each=8)
all_species <- rbind(digitata,grandidieri,mada,perrieri,rubrostipa,suare,za)
all_species <- cbind(species,all_species)

## altitudinal change
digitata_alt <- read.table(paste0("Adansonia.digitata/alt_fut.txt"), header=T,sep="\t")
grandidieri_alt <- read.table(paste0("Adansonia.grandidieri/alt_fut.txt"), header=T,sep="\t")
mada_alt <- read.table(paste0("Adansonia.madagascariensis/alt_fut.txt"), header=T,sep="\t")
perrieri_alt <- read.table(paste0("Adansonia.perrieri/alt_fut.txt"), header=T,sep="\t")
rubrostipa_alt <- read.table(paste0("Adansonia.rubrostipa/alt_fut.txt"), header=T,sep="\t")
suare_alt <- read.table(paste0("Adansonia.suarezensis/alt_fut.txt"), header=T,sep="\t")
za_alt <- read.table(paste0("Adansonia.za/alt_fut.txt"), header=T,sep="\t")

all_species_alt <- rbind(digitata_alt,grandidieri_alt,mada_alt,
                         perrieri_alt,rubrostipa_alt,suare_alt,za_alt)

all_species_cool <- cbind(all_species,all_species_alt)
tail(all_species_cool)

final_table1 <- subset(all_species_cool, select = -c(5,6,11,12,13,14,15,17,18))

colnames(final_table1)<- c("Baobab species","SDAp (km²)","RCP","Year","Dispersal",
                           "SDAf (km²)","Change SDApf (%)","Current altitude mean (m)",
                           "Future altitude mean (m)","Altitude range shift (%)")
head(final_table1)
write.table(final_table1,paste0("./outputs/table1.txt"),sep="\t")

############################
### importance variable table
digitata_vi <- read.table(paste0("Adansonia.digitata/varimp.txt"), header=T,sep="\t")
grandidieri_vi <- read.table(paste0("Adansonia.grandidieri/varimp.txt"), header=T,sep="\t")
mada_vi <- read.table(paste0("Adansonia.madagascariensis/varimp.txt"), header=T,sep="\t")
perrieri_vi <- read.table(paste0("Adansonia.perrieri/varimp.txt"), header=T,sep="\t")
rubrostipa_vi <- read.table(paste0("Adansonia.rubrostipa/varimp.txt"), header=T,sep="\t")
suare_vi <- read.table(paste0("Adansonia.suarezensis/varimp.txt"), header=T,sep="\t")
za_vi <- read.table(paste0("Adansonia.za/varimp.txt"), header=T,sep="\t")

library(data.table)
all_species_vi <- rbind(digitata_vi,grandidieri_vi,mada_vi,
                         perrieri_vi,rubrostipa_vi,suare_vi,za_vi)
all_species_vi$mean <- c("NA")
all_species_final <- subset(all_species_vi, select = -c(5,6))

setDT(all_species_final)
trying <- all_species_final[, .(Mean = rowMeans(.SD)), by = mean]
trying$mean <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                 "A. perrieri", "A. rubrostipa","A. suarezensis",
                 "A. za"),each=4)
trying <- as.data.frame(trying)
all_species_vi_done <- cbind(trying,all_species_final)
all_species_vi_done <- subset(all_species_vi_done, select = -c(7))

all_species_vi_done$clim_var <- rep(c("tmean","tseas","prec","cwd"),7)
all_species_vi_done <- all_species_vi_done[, c(1, 7, 3, 4, 5, 6, 2)]
colnames(all_species_vi_done)<- c("Baobab_species","Clim_var",
                                  "GLM","GAM","RF","Maxent","Mean")
try <- all_species_vi_done %>% group_by(Baobab_species) %>%
  mutate(Rank = row_number(max(Mean) - Mean))
View(try)

write.table(try,paste0("./outputs/table2vi.txt"),sep="\t")

#### Table 3 with performance indexes

suare_dig <- read.table(paste0("Adansonia.digitata/performance_ca.txt"), header=T,sep="\t")
suare_grandi <- read.table(paste0("Adansonia.grandidieri/performance_ca.txt"), header=T,sep="\t")
suare_mada <- read.table(paste0("Adansonia.madagascariensis/performance_ca.txt"), header=T,sep="\t")
suare_perri <- read.table(paste0("Adansonia.perrieri/performance_ca.txt"), header=T,sep="\t")
suare_rubro <- read.table(paste0("Adansonia.rubrostipa/performance_ca.txt"), header=T,sep="\t")
suare_suare <- read.table(paste0("Adansonia.suarezensis/performance_ca.txt"), header=T,sep="\t")
suare_za <- read.table(paste0("Adansonia.za/performance_ca.txt"), header=T,sep="\t")
ca_species <- rbind(suare_dig,suare_grandi,suare_mada,suare_perri,suare_rubro,
                    suare_suare,suare_za)

rownames(ca_species) <- c("A. digitata","A. grandidieri","A. madagascariensis", 
                            "A. perrieri", "A. rubrostipa","A. suarezensis","A. za")

only_sen_spe_tss <- subset(ca_species, select=c("Sen", "Spe","TSS"))

write.table(only_sen_spe_tss,paste0("./outputs/table3_performance_ready.txt"),sep="\t")


####### #### Table A2 with performance indexes

perf_dig <- read.table(paste0("Adansonia.digitata/current_model_evaluation.txt"), header=T,sep="\t")
perf_grand <- read.table(paste0("Adansonia.grandidieri/current_model_evaluation.txt"), header=T,sep="\t")
perf_mada <- read.table(paste0("Adansonia.madagascariensis/current_model_evaluation.txt"), header=T,sep="\t")
perf_perri <- read.table(paste0("Adansonia.perrieri/current_model_evaluation.txt"), header=T,sep="\t")
perf_rubro <- read.table(paste0("Adansonia.rubrostipa/current_model_evaluation.txt"), header=T,sep="\t")
perf_suare <- read.table(paste0("Adansonia.suarezensis/current_model_evaluation.txt"), header=T,sep="\t")
perf_za <- read.table(paste0("Adansonia.za/current_model_evaluation.txt"), header=T,sep="\t")

perf_species <- rbind(perf_dig,perf_grand,perf_mada,perf_perri,perf_rubro,
                    perf_suare,perf_za)

# using subset function
newdata <- subset(perf_species, Index == "Testing.data",
                  select=c(wIndex,Model, Run, Value))


newdatap <- newdata[newdata$Run != 'Full', ]
newdatap$species <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                            "A. perrieri", "A. rubrostipa","A. suarezensis","A. za")
                        ,each=60) 

partial_final <-  aggregate(Value~species + wIndex + Model, FUN=mean, data=newdatap)
## to calculate mean of the 5 runs

partial_final1 <- partial_final %>%
  filter(wIndex != 'KAPPA')#omit tempcol in output


write.table(partial_final1,paste0("./outputs/tableA2_performance_partial.txt"),sep="\t")

### To calculate mean over the full dataset
newdataf <- newdata[newdata$Run == 'Full', ]
newdataf$species <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                          "A. perrieri", "A. rubrostipa","A. suarezensis","A. za")
                        ,each=12)
full_final <-  aggregate(Value~species + wIndex + Model, FUN=mean, data=newdataf)
full_final1 <- full_final %>%
  filter(wIndex != 'KAPPA')#omit tempcol in output

partial_final1
full_final1
write.table(full_final1,paste0("./outputs/tableA2_performance_full.txt"),sep="\t")

############# Table 4

nichesf_dig <- read.table(paste0("Adansonia.digitata/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_dig <- read.table(paste0("Adansonia.digitata/niche.txt"), header=T,sep="\t")
niches_dig <- cbind(nichesf_dig,nichesc_dig)
niches_dig$species <- rep(c("A. digitata"),each=3)

nichesf_grand <- read.table(paste0("Adansonia.grandidieri/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_grand <- read.table(paste0("Adansonia.grandidieri/niche.txt"), header=T,sep="\t")
niches_grand <- cbind(nichesf_grand,nichesc_grand)
niches_grand$species <- rep(c("A. grandidieri"),each=3)

nichesf_mada <- read.table(paste0("Adansonia.madagascariensis/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_mada <- read.table(paste0("Adansonia.madagascariensis/niche.txt"), header=T,sep="\t")
niches_mada <- cbind(nichesf_mada,nichesc_mada)
niches_mada$species <- rep(c("A. madagascariensis"),each=3)

nichesf_perri <- read.table(paste0("Adansonia.perrieri/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_perri <- read.table(paste0("Adansonia.perrieri/niche.txt"), header=T,sep="\t")
niches_perri <- cbind(nichesf_perri,nichesc_perri)
niches_perri$species <- rep(c("A. perrieri"),each=3)

nichesf_rubro <- read.table(paste0("Adansonia.rubrostipa/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_rubro <- read.table(paste0("Adansonia.rubrostipa/niche.txt"), header=T,sep="\t")
niches_rubro <- cbind(nichesf_rubro,nichesc_rubro)
niches_rubro$species <- rep(c("A. rubrostipa"),each=3)

nichesf_suare <- read.table(paste0("Adansonia.suarezensis/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_suare <- read.table(paste0("Adansonia.suarezensis/niche.txt"), header=T,sep="\t")
niches_suare <- cbind(nichesf_suare,nichesc_suare)
niches_suare$species <- rep(c("A. suarezensis"),each=3)

nichesf_za <- read.table(paste0("Adansonia.za/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_za <- read.table(paste0("Adansonia.za/niche.txt"), header=T,sep="\t")
niches_za <- cbind(nichesf_za,nichesc_za)
niches_za$species <- rep(c("A. za"),each=3)

all_niches <- cbind(niches_dig,niches_grand,niches_mada,niches_perri,
                    niches_rubro,niches_suare,niches_za)

write.table(all_niches,paste0("./outputs/tableA3_niches_fut_cur_comparison.txt"),sep="\t")

# ===========
# End of file
# ===========