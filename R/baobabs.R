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
# Plot relationship between climatic, elevational and lat gradients in Mada
# ==========================================================
source(here("R/plot_gradients_Mada.R"))

# ==========================================================
# Plot box plots comparing current and future projections and the changes in Lat and Elevation
# ==========================================================
source(here("R/plot_boxplot_alt_lat_species.R"))

# ==========================================================
# Plot World Seasonality Map
# ==========================================================
source(here("R/plot_world_seasonality.R"))

# ==========================================================
# Create all study's tables
# ==========================================================
source(here("R/all_tables.R"))

# ===========
# End of file
# ===========
