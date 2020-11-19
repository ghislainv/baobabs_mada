#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com / mario.tagliari@posgrad.ufsc.br
# license         :GPLv3
# ==============================================================================

## Set environmental variable
## For MAXENT.Phillips with JAVA to work on RStudio server
Sys.unsetenv("DISPLAY")

## Libraries
pkg <- c("curl", "readr", "dplyr", "rgdal", "sp", "raster",
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
## Plot environmental variables
if(!file.exists("outputs/environ.pdf")) {
  dir.create("outputs",showWarnings=FALSE)
  pdf(file="outputs/environ.pdf")
  plot(s)
  dev.off()
}

##==================================
## Occurence data for baobab species
##==================================

## Building the occurrence dataset from raw data
source("R/data_baobabs.R")
## Load dataset
df.orig <- read.csv(file=paste0(dir_baobabs,"data_Adansonia.csv"),header=TRUE,sep=",")
## Make a SpatialPointsDataFrame object
coords <- cbind(df.orig$Long,df.orig$Lat)
## see PROJ.6 style here – thanks universe! https://www.gaia-gis.it/fossil/libspatialite/wiki?name=PROJ.6
df.sp <- SpatialPointsDataFrame(coords,data=df.orig,proj4string=CRS("+proj=longlat +south +datum=WGS84 +no_defs +type=crs")) # lat long here proj!
#df.sp <- SpatialPointsDataFrame(coords,data=df.orig,proj4string=CRS("+init=epsg:4326")) # our old version
## Reproject into UTM 38S - SET exactly as "s" raster created above
df.sp <- spTransform(df.sp,CRS("+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +type=crs")) # utm proj!
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
source("R/run_species.R")
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
foreach(i=1:n.species,.packages=pkg.names.clust) %dopar% run.species(i, run.models=TRUE)
## Stop the cluster
stopCluster(clust)
## Time computation
t.stop <- Sys.time() ## Stop the clock
t.diff <- difftime(t.stop,t.start,units="min")
cat(paste0("Computation time (min): ",round(t.diff,2)),file="outputs/computation_time.txt")

##==================================
## Compute anomalies
##==================================

library(raster)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)

# Present climate
current <- stack("data/gisdata/sdm_variables/current.tif")
names(current) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),
                    paste("prec",1:12,sep=""),paste("bio",1:19,sep=""),
                    paste("pet",1:12,sep=""),"pet","cwd","ndm")

# Future climate
he_85_2080 <- stack("data/gisdata/sdm_variables/he_85_2080.tif")
gs_85_2080 <- stack("data/gisdata/sdm_variables/gs_85_2080.tif")
no_85_2080 <- stack("data/gisdata/sdm_variables/no_85_2080.tif")
names(no_85_2080) <- names(gs_85_2080) <- names(he_85_2080) <- names(current)

# Compute future bioclimatic variable mean and anomalies
var <- c("bio1", "bio4", "bio12", "cwd")
var_85_2080 <- stack()
ano_85_2080 <- stack()
for (i in 1:length(var)) {
  v <- var[i]
  cat(paste0("Variable: ", v, "\n"))
  v_85_2080 <- stack(he_85_2080[[v]],
                     gs_85_2080[[v]],
                     no_85_2080[[v]])
  var_85_2080 <- addLayer(var_85_2080, mean(v_85_2080))
  # Note: E(X+Y) = E(X) + E(Y) so, mean of differences = difference of the means
  ano_85_2080 <- addLayer(ano_85_2080, mean(v_85_2080)-current[[v]])
}
var_name <- c("tmean", "tseas", "prec", "cwd")
names(var_85_2080) <- var_name
names(ano_85_2080) <- var_name

## Remove data for Comoro Islands
bbCom <- extent(xmin(current),600000,8500000,ymax(current)) # bounding-box
cellsCom <- cellsFromExtent(current,bbCom)
values(var_85_2080)[cellsCom,] <- NA
var_85_2080 <- stack(var_85_2080) # Transform back from RasterBrick to RasterStack
values(ano_85_2080)[cellsCom,] <- NA
ano_85_2080 <- stack(ano_85_2080) # Transform back from RasterBrick to RasterStack

## Setting basic theme options for plot with ggplot2
theme_base <- theme(
    ## Axis
    axis.line=element_blank(),
    axis.text=element_blank(),
    axis.ticks=element_blank(),
    axis.title=element_blank(),
    ## Legend
    legend.position="bottom",
    legend.title=element_blank(),
    legend.text=element_text(size=10),
    legend.key.height=unit(0.5,"line"),
    legend.key.width=unit(1.5,"line"),
    legend.box.background=element_blank(),
    ## Plot
    plot.title=element_text(hjust=0.5,size=12),
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
plot_anomaly <- function(r, title, label="x") {
  rdf <- data.frame(rasterToPoints(r)) # To plot raster with geom_raster()
  names(rdf) <- c("x", "y", "z")
  p <- ggplot(NULL, aes(x, y)) +
      geom_raster(data=rdf, aes(fill=z)) +
      theme_bw() + theme_base +
      coord_fixed(xlim=c(313000,1090000), ylim=c(7167000,8676000)) +
      annotate("text",x=600000,y=8500000,label=label,hjust=1,vjust=0,size=5) +
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
## tmean
col_scale_var_tmean <- scale_fill_gradientn(
    colours=viridis(255, option="B", direction=1),
    na.value="transparent",
    values=rescale(seq(140,320,l=255),140,320),
    limits=c(140,320),
    breaks=seq(140,320,l=4),
    labels=seq(140,320,l=4)
)
col_scale_ano_tmean <- scale_fill_gradientn(
    colours=viridis(255, option="B", direction=1),
    na.value="transparent",
    values=rescale(seq(25,40,l=255),25,40),
    limits=c(25,40),
    breaks=seq(25,40,l=4),
    labels=seq(25,40,l=4)
)
## tseas
col_scale_var_tseas <- scale_fill_gradientn(
    colours=viridis(255, option="D", direction=1),
    na.value="transparent",
    values=rescale(seq(1000,3400,l=255),1000,3400),
    limits=c(1000,3400),
    breaks=seq(1000,3400,l=4),
    labels=seq(1000,3400,l=4)
)
col_scale_ano_tseas <- scale_fill_gradientn(
    colours=c(rep(grey(0.5),2),viridis(255, option="D", direction=1)),
    na.value="transparent",
    values=rescale(c(-50,-0.01,seq(0,300,l=255)),-50,300),
    limits=c(-50,300),
    breaks=c(-50,seq(0,300,l=3)),
    labels=c(-50,seq(0,300,l=3))
)
## prec
col_scale_var_prec <- scale_fill_gradientn(
    colours=viridis(255, option="E", direction=-1),
    na.value="transparent",
    values=rescale(seq(0,3500),0,3500),
    limits=c(0,3500),
    breaks=c(0,1000,2250,3500),
    labels=c(0,1000,2250,3500)
)
col_scale_ano_prec <- scale_fill_gradientn(
    colours=c(viridis(255, option="E", direction=-1),rep(grey(0.5),2)),
    na.value="transparent",
    values=rescale(c(seq(-300,0,l=255),0.01,50),-300,50),
    limits=c(-300,50),
    breaks=c(seq(-300,0,l=3),50),
    labels=c(seq(-300,0,l=3),50)
)
## cwd
col_scale_var_cwd <- scale_fill_gradientn(
    colours=viridis(255, option="C", direction=1),
    na.value="transparent",
    values=rescale(seq(0,2400,l=255),0,2400),
    limits=c(0,2400),
    breaks=seq(0,2400,l=4),
    labels=seq(0,2400,l=4)
)
col_scale_ano_cwd <- scale_fill_gradientn(
    colours=viridis(255, option="C", direction=1),
    na.value="transparent",
    values=rescale(seq(0,1500,l=255),0,1500),
    limits=c(0,1500),
    breaks=seq(0,1500,l=4),
    labels=seq(0,1500,l=4)
)

## Plot anomalies for each climatic variable
## tmean
p1 <- plot_anomaly(r=var_85_2080[["tmean"]], label="(a)",
                   title="Temperature\n(°C x 10)") +
    col_scale_var_tmean
p2 <- plot_anomaly(r=ano_85_2080[["tmean"]], label="(a')",
                   title="") +
    col_scale_ano_tmean
p3 <- plot_anomaly(r=var_85_2080[["tseas"]], label="(b)",
                   title="T. seasonality\n(°C sd x 100)") +
    col_scale_var_tseas
p4 <- plot_anomaly(r=ano_85_2080[["tseas"]], label="(b')",
                   title="") +
    col_scale_ano_tseas
## prec
p5 <- plot_anomaly(r=var_85_2080[["prec"]], label="(c)",
                   title="Precipitation\n(mm/y)") +
    col_scale_var_prec
p6 <- plot_anomaly(r=ano_85_2080[["prec"]], label="(c')",
                   title="") +
    col_scale_ano_prec
## cwd
p7 <- plot_anomaly(r=var_85_2080[["cwd"]], label="(d)",
                   title="Climatic water deficit\n(mm/y)") +
    col_scale_var_cwd
p8 <- plot_anomaly(r=ano_85_2080[["cwd"]], label="(d')",
                   title="") +
    col_scale_ano_cwd

## Combine plots
lay <- rbind(c(1,rep(seq(3,9,by=2),each=3)),
             c(2,rep(seq(4,10,by=2),each=3)))
tgrob_pres <- textGrob("Present climate",
                       rot=90, gp=gpar(cex=1.25), hjust=0.5, vjust=0.5)
tgrob_fut <- textGrob("Future anomaly\n(RCP 8.5, 2085)",
                      rot=90, gp=gpar(cex=1.25), hjust=0.5, vjust=0.5)
plot_anomalies <- grid.arrange(tgrob_pres, tgrob_fut,
                               p1, p2, p3, p4, p5, p6, p7, p8,
                               layout_matrix=lay)

ggsave(filename="outputs/Fig1_climatic_anomalies.png", plot=plot_anomalies,
       width=8, height=7, dpi="print")
ggsave(filename="outputs/Fig1_climatic_anomalies.pdf", plot=plot_anomalies,
       width=8, height=7, dpi="print")

##==================================
## Climate change in SDA
##==================================

## Draw points in the SDA and extract future environmental variables
# In SDA

for (i in 1: length(sp.dir)) {
  pred <-  stack(paste0(sp.dir[i],"/proj_current/proj_current_",sp.dir[i],"_ensemble.grd"))
  ca <- pred[[1]]  
  wC.anomalies <- which(values(ca)>=500) 
  nC.anomalies <- length(wC.anomalies)
  Samp.anomalies <- if (nC.anomalies>1000) {sample(wC.anomalies,1000,replace=FALSE)} else {wC.anomalies}
  mapmat.df.anomalies <- as.data.frame(var_85_2080)[Samp.anomalies,] 
  
  ## table to compare current and future bioclimatic changes for each species 
  # Used this table to create density niche curves to compare future baobabs niche
  ## over current distribution (ca)
  names(mapmat.df.anomalies) <- c("tmeanf","tseasf","precf","cwdf")
  mapmat.df <- read.csv(file=paste0("./",sp.dir[i],"/","niche_graph_species.csv"),sep=";")
  head(mapmat.df)
  mapmat.final <- cbind(mapmat.df.anomalies,mapmat.df)
  mapmat.f <- mapmat.final[,c(1,2,3,4,6,7,8,9,10,11)]
  write.csv2(mapmat.final,paste0("./",sp.dir[i],"/","niche_graph_species_compared_anomaly.csv"))
  write.table(mapmat.final,paste0("./",sp.dir[i],"/niche_graph_species_compared_anomaly.txt"),sep="\t")
  
  #### Future niche over current SDA!!! NO ANOMALY
  wC.future <- which(values(ca)>=500)
  mapmat.df.future <- as.data.frame(var_85_2080)[wC.future,]
  names(mapmat.df.future) <- c("tmeanf","tseasf","precf","cwdf")
  Mean_ok_future_fut <- round(apply(mapmat.df.anomalies,2,mean,na.rm=TRUE))
  q_ok_future_fut <- round(apply(mapmat.df.anomalies,2,quantile,c(0.025,0.975),na.rm=TRUE))
  niche_ok_future <- as.data.frame(rbind(Mean_ok_future_fut,q_ok_future_fut))
  write.csv2(niche_ok_future,paste0("./",sp.dir[i],"/","mean_niche_with_future.csv"))
  write.table(niche_ok_future,paste0("./",sp.dir[i],"/mean_niche_with_future.txt"),sep="\t")
  
}

###################################################################################
######## Generating variables histograms graphs for all seven baobab species ######
###################################################################################
## read species environmental data
database1_suare <- read.csv(file=paste0("Adansonia.suarezensis/","niche_graph_species.csv"), header=T,sep=";")
database2_digi <- read.csv(file=paste0("Adansonia.digitata/","niche_graph_species.csv"), header=T,sep=";")
database3_grand <- read.csv(file=paste0("Adansonia.grandidieri/","niche_graph_species.csv"), header=T,sep=";")
database4_mada <- read.csv(file=paste0("Adansonia.madagascariensis/","niche_graph_species.csv"), header=T,sep=";")
database5_perri <- read.csv(file=paste0("Adansonia.perrieri/","niche_graph_species.csv"), header=T,sep=";")
database6_rubro <- read.csv(file=paste0("Adansonia.rubrostipa/","niche_graph_species.csv"), header=T,sep=";")
database7_za <- read.csv(file=paste0("Adansonia.za/","niche_graph_species.csv"), header=T,sep=";")
data_teste <- rbind(database2_digi,database3_grand,
              database4_mada,database5_perri,database6_rubro,database1_suare,database7_za)
levels(data_teste$species)
quantiles <- data_teste[,c(2,3,4,5,6)]
quantiles <- round(apply(quantiles,2,quantile,c(0.025,0.975),na.rm=TRUE)) #use it to set new plot

data_teste$species = factor(data_teste$species,
                            levels = c("Adansonia.digitata","Adansonia.grandidieri",
                                       "Adansonia.madagascariensis","Adansonia.perrieri",
                                       "Adansonia.rubrostipa","Adansonia.suarezensis",
                                       "Adansonia.za"),
                            labels = c("A. digitata","A. grandidieri",
                                       "A. madagascariensis","A. perrieri",
                                       "A. rubrostipa","A. suarezensis",
                                       "A. za"))
# Density plots ## Temp. Seasonality
range(data_teste$tseas)
# generate break positions
breaks <- c(round(seq(min(data_teste$tseas),max(data_teste$tseas),length=6)))
# and labels
labels = as.character(breaks) # labels must be the same as breaks, otherwise error
# plot the map
my_plot = ggplot(data_teste, aes(x=tseas, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + #specific palette
  geom_vline(xintercept = range(data_teste$tseas), show.legend = T, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$tseas), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
my_plot_tseas =  my_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 #legend.text = element_text(face= "italic",size=15),
                                 legend.title=element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100ºC)", y = "Density",size=5) +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.6),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.6), colour="black")) +
  theme(legend.position="none")
  

# Density plots ### Annual Mean Temperature
# generate break positions
library(lemon)
breaks <- c(round(seq(min(data_teste$tmean),max(data_teste$tmean),length=6)))
labels = as.character(breaks)
my_plot2 = ggplot(data_teste, aes(x=tmean, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = range(data_teste$tmean), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$tmean), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
  my_plot_tmean =  my_plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  legend.text = element_text(face= "italic",size=12),
                                  legend.title=element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10) ", y = "Density",size=5) +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.6),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.6), colour="black")) +
  theme(legend.key.size = unit(0.5, "cm"))
  
my_plot_tmean <- reposition_legend(my_plot_tmean, 'top left', offset = 0.05)

  # Density plots ### Annual Mean Precipitation
# generate break positions
breaks <- c(round(seq(min(data_teste$prec),max(data_teste$prec),length=6)))
labels = as.character(breaks)
# plot and be happy ;)
my_plot3 = ggplot(data_teste, aes(x=prec, color=species)) + geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = range(data_teste$prec), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$prec), breaks = breaks, labels = labels)
my_plot_prec =  my_plot3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             #legend.text = element_text(face= "italic",size=15),
                             legend.title=element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "Density",size=5) +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.6),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.6), colour="black")) +
  theme(legend.position="none")
  
# Density plots ### Climatic Water Deficit
# generate break positions
breaks <- c(round(seq(min(data_teste$cwd),max(data_teste$cwd),length=6)))
labels = as.character(breaks)
my_plot4 = ggplot(data_teste, aes(x=cwd, color=species)) + geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = range(data_teste$cwd), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$cwd), breaks = breaks, labels = labels)
my_plot_cwd =  my_plot4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             #legend.text = element_text(face= "italic",size=15),
                             legend.title=element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=5) +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.6),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.6), colour="black")) +
  theme(legend.position="none")
  

  
## Combine plots 
lay_2 <- cbind(c(rep(seq(1,4,by=1),each=3)),
               c(rep(seq(2,4,by=2),each=3)))
colnames(lay_2) <-c("a","b")
lay_2 <- subset( lay_2, select = -b )

plot_densities <- grid.arrange(my_plot_tmean, my_plot_tseas, my_plot_prec, my_plot_cwd,
                               layout_matrix=lay_2)

ggsave(filename="outputs/Fig_ap_sps_niche.png", plot=plot_densities,
       width=8, height=11, dpi="print")

ggsave(filename="outputs/Fig_ap1_sps_niche.pdf", plot=plot_densities,
       width=8, height=11, dpi="print")


### Comparing bioclimatic niche of each species inside current SDA according to the 
### most important variables ######
# A_suarezensis
# Import data_set
a.suare <- read.csv(file=paste0("Adansonia.suarezensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
# 1st importance variable - Tseas
range(a.suare$tseas) # 1202 1421
range(a.suare$tseasf) # 1307.000 1529.333
mean(a.suare$tseas) # 1276.642
mean(a.suare$tseasf) # 1387.751
breaks <- c(round(seq(min(a.suare$tseas),max(a.suare$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_suar = ggplot(a.suare, aes(x=tseas, y=..density..)) + 
  geom_density(aes(tseas,fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(1175,1550), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_suar =  plot.seas_suar + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (ºC sd x 100)", y = "A. suarezensis",size=5) +
  labs(col = "") +
  annotate("text", x  = 1190, y = 0.0105 , size=7, label = "(k)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


### second importance variable - Prec
range(a.suare$prec) # 1113 - 1522
range(a.suare$precf) # 922.3 1325
breaks <- c(round(seq(min(a.suare$precf),max(a.suare$prec),length=8)))
labels = as.character(breaks)

plot.prec_suar = ggplot(a.suare, aes(x=prec, y=..density..)) + 
  geom_density(aes(prec,fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(800,1600), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_suar =  plot.prec_suar + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "",size=5) +
  labs(col = "") +
  annotate("text", x  = 825, y = 0.0032 , size=7, label = "(l)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")
  
   
################################
# Import data_set A perrieri
#1st importance variable - Tseas 
a.perrieri <- read.csv(file=paste0("Adansonia.perrieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.perrieri)
range(a.perrieri$tseas) # 860 1958
range(a.perrieri$tseasf) # 1040 2177
round(mean(a.perrieri$tseas)) # 1392
round(mean(a.perrieri$tseasf)) # 1535
breaks <- c(round(seq(min(a.perrieri$tseas),max(a.perrieri$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_per = ggplot(a.perrieri, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(750, 2250), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_perrieri =  plot.seas_per + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        legend.text = element_blank(),
                                        panel.background = element_blank(), 
                                        axis.line = element_line(colour = "black"))+
  
  labs(x="Temp. Seasonality (ºC sd x 100)",  y = "A. perrieri",size=5) +
  labs(col = "") +
  annotate("text", x  = 805, y = 0.0025 , size=7, label = "(g)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

## 2nd Most IV A perrieri - climatic water deficit
breaks <- c(round(seq(min(a.perrieri$cwd),max(a.perrieri$cwdf),length=8)))
labels = as.character(breaks)

plot.cwd_perrieri = ggplot(a.perrieri, aes(x=cwd, y=..density..)) + 
  geom_density(aes(cwd, fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(30, 1700), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.cwd_perrieri =  plot.cwd_perrieri + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Climatic Water Deficit (mm)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 80, y = 0.0023 , size=7, label = "(h)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

############## A.rubrostipa
# Import dataset
a.rubrostipa <- read.csv(file=paste0("Adansonia.rubrostipa/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

## 1st importance variable - Climatic Water Deficit
range(a.rubrostipa$cwd) # 704 967
range(a.rubrostipa$cwdf) # 1289.000 2029.333
mean(a.rubrostipa$cwd) # 811.045
mean(a.rubrostipa$cwdf) # 1672.713
breaks <- c(round(seq(min(a.rubrostipa$cwd),max(a.rubrostipa$cwdf),length=6)))
labels = as.character(breaks)

plot.cwd_rubro = ggplot(a.rubrostipa, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(650, 2100), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.cwd_rubro =  plot.cwd_rubro + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_blank(),
                             panel.background = element_blank(), 
                             axis.line = element_line(colour = "black"))+
  
  labs(x="Climatic Water Deficit (mm)", y = "A. rubrostipa",size=5) +
  labs(col = "") +
  annotate("text", x  = 698, y = 0.0065 , size=7, label = "(i)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


### 2nd Importance Variable - Prec
range(a.rubrostipa$prec) # 340 1726
round(range(a.rubrostipa$precf)) # 326 1559
mean(a.rubrostipa$prec) # 1093.882
mean(a.rubrostipa$precf) # 954.4527
breaks <- c(round(seq(min(a.rubrostipa$prec),max(a.rubrostipa$precf),length=8)))
labels = as.character(breaks)

plot.prec_rubro = ggplot(a.rubrostipa, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(150, 1900), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_rubro =  plot.prec_rubro + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             legend.text = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 200, y = 0.00125 , size=7, label = "(j)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


############## A. madagascariensis 
## Import dataset
a.madagascariensis <- read.csv(file=paste0("Adansonia.madagascariensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

# 1st importance variable - Tseas
range(a.madagascariensis$tseas) # 903 2076
range(a.madagascariensis$tseasf) # 1058.667 2186.667

breaks <- c(round(seq(min(a.madagascariensis$tseas),max(a.madagascariensis$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_mada = ggplot(a.madagascariensis, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(800,2300), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_mada =  plot.seas_mada + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                         legend.text = element_blank(),
                                         panel.background = element_blank(), 
                                         axis.line = element_line(colour = "black"))+
  
  labs(x="Temp. Seasonality (ºC sd x 100)", y = "A. madagascariensis",size=5) +
  labs(col = "") +
  annotate("text", x  = 875, y = 0.0025 , size=7, label = "(e)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.5),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


### 2nd importance variable - Tmean
range(a.madagascariensis$tmean) # 243 275
range(a.madagascariensis$tmeanf) #  271.0000 311.6667
mean(a.madagascariensis$tmean) # 262
mean(a.madagascariensis$tmeanf) # 296

breaks <- c(round(seq(min(a.madagascariensis$tmean),max(a.madagascariensis$tmeanf),length=8)))
labels = as.character(breaks)

plot.mean_mada = ggplot(a.madagascariensis, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(240,320), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.mean_mada =  plot.mean_mada + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Annual Mean Temperature (ºC x 10)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 241.5, y = 0.075 , size=7, label = "(f)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

############## A. grandidieri
## Import dataset
a.grandidieri <- read.csv(file=paste0("Adansonia.grandidieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

# 1st Importance Variable - Prec
range(a.grandidieri$prec) # 381 1033
range(a.grandidieri$precf) # 347.0000 962.6667
mean(a.grandidieri$prec) # 747.566
mean(a.grandidieri$precf) # 707.908
breaks <- c(round(seq(min(a.grandidieri$precf),max(a.grandidieri$prec),length=8)))
labels = as.character(breaks)

plot.prec_grand = ggplot(a.grandidieri, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) + 
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(300,1100), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_grand =  plot.prec_grand + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           legend.text = element_blank(),
                                           panel.background = element_blank(), 
                                           axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "A. grandidieri",size=5) +
  labs(col = "") +
  annotate("text", x  = 336, y = 0.005 , size=7, label = "(c)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


## 2nd Importance Variable - Tmean
range(a.grandidieri$tmean) # 242 - 261
range(a.grandidieri$tmeanf) # 274 7
mean(a.grandidieri$tmean) # 250
mean(a.grandidieri$tmeanf) # 285
breaks <- c(round(seq(min(a.grandidieri$tmean),max(a.grandidieri$tmeanf),length=8)))
labels = as.character(breaks)

plot.mean_grand = ggplot(a.grandidieri, aes(x=tmean, y=..density..)) + 
  geom_density(aes(tmean,fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(237, 305), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.mean_grand =  plot.mean_grand + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           legend.text = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Annual Mean Temperature (ºC x 10)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 238, y = 0.08 , size=7, label = "(d)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


############## A za
## Import dataset
a.za <- read.csv(file=paste0("Adansonia.za/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.za)
# 1st Importance Variable - Prec
range(a.za$prec) # 349 1842
range(a.za$precf) # 325 1542

breaks <- c(round(seq(min(a.za$precf),max(a.za$prec),length=8)))
labels = as.character(breaks)

plot.prec_za = ggplot(a.za, aes(x=prec, y=..density..)) + 
  geom_density(aes(prec,fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(200,1950), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_za =  plot.prec_za + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                     legend.text = element_blank(),
                                     panel.background = element_blank(), 
                                     axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "A. za",size=5) +
  labs(col = "") +
  annotate("text", x  = 250, y = 0.002 , size=7, label = "(m)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

## 2nd Importance Variable - Tmean
range(a.za$tmean) # 214 275
range(a.za$tmeanf) # 250 311
mean(a.za$tmean) # 243
mean(a.za$tmeanf) # 278
breaks <- c(round(seq(min(a.za$tmean),max(a.za$tmeanf),length=8)))
labels = as.character(breaks)

plot.mean_za = ggplot(a.za, aes(x=tmean, y=..density..)) + 
  geom_density(aes(tmean,fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(200, 320), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.mean_za =  plot.mean_za + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                     legend.text = element_blank(),
                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Annual Mean Temperature (ºC x 10)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 204, y = 0.03 , size=7, label = "(n)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

##############A digitata 
## Import dataset
a.digitata <- read.csv(file=paste0("Adansonia.digitata/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.digitata)

##1st Importance Variable - Tseas

# range tseas 
range(a.digitata$tseas) # 882 2473
range(a.digitata$tseasf) # 1066 2561
mean(a.digitata$tseas) # 1526
mean(a.digitata$tseasf) # 1640
breaks <- c(round(seq(min(a.digitata$tseas),max(a.digitata$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_dig = ggplot(a.digitata, na.rm=T, aes(x=tseas, y=..density..)) + 
  
  geom_density(aes(tseas, fill=species,color= 'tseas'), alpha=.5, na.rm=T) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(800,2700), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_dig =  plot.seas_dig + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       legend.text = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Temp. Seasonality (ºC sd x 100)",y="A. digitata",size=5) +
  labs(col = "") +
  annotate("text", x  = 882, y = 0.004 , size=7, label = "(a)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


## 2nd Importance Variable - Climatic Water Deficit
range(a.digitata$cwd) # 625 - 920
range(a.digitata$cwdf) # 1218.667 2001.000
mean(a.digitata$cwd) # 796.666
mean(a.digitata$cwdf) # 1692
breaks <- c(round(seq(min(a.digitata$cwd),max(a.digitata$cwdf),length=8)))
labels = as.character(breaks)

plot.cwd_dig = ggplot(a.digitata, aes(x=cwd, y=..density..)) + 
  
  geom_density(aes(cwd,fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(600, 2100), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

  plot.cwd_dig =  plot.cwd_dig + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        legend.text = element_blank(),
                                        panel.background = element_blank(), 
                                        axis.line = element_line(colour = "black"))+
    
  labs(x="Climatic Water Deficit (mm)", y = "",size=5) +
         labs(col = "") +
         annotate("text", x  = 628, y = 0.007 , size=7, label = "(b)") +
         theme(axis.title.x = element_text(size = rel(2),colour="black")) +
         theme(axis.title.y = element_text(size = rel(2),colour="black")) +
         theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
         theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
         theme(legend.position="none")



## Combine plots 
lay_3 <- rbind(c(rep(seq(1,2,by=1),each=3)),
               c(rep(seq(3,4,by=1),each=3)),
               c(rep(seq(5,6,by=1),each=3)),
               c(rep(seq(7,8,by=1),each=3)),
               c(rep(seq(9,10,by=1),each=3)),
               c(rep(seq(11,12,by=1),each=3)),
               c(rep(seq(13,14,by=1),each=3)))

plot_densities_curves <- grid.arrange(plot.seas_dig,plot.cwd_dig,
                                      plot.prec_grand,plot.mean_grand,
                                      plot.seas_mada,plot.mean_mada,
                                      plot.seas_perrieri,plot.cwd_perrieri,
                                      plot.cwd_rubro,plot.prec_rubro,
                                      plot.seas_suar, plot.prec_suar,
                                      plot.prec_za,plot.mean_za,
                                      layout_matrix=lay_3)

ggsave(file=paste0("./outputs/all_species_current_future_niche_comparison.pdf"),
       plot=plot_densities_curves,width=18,height=15,dpi="print")

ggsave(file=paste0("./outputs/all_species_current_future_niche_comparison_2.png"),
       plot=plot_densities_curves,width=18,height=15,dpi="print")

######################################################################
###### baobaobs vulnerability to climate change 
#####################################################################
sps <- c("Adansonia.digitata","Adansonia.grandidieri","Adansonia.madagascariensis",
         "Adansonia.perrieri","Adansonia.rubrostipa","Adansonia.suarezensis","Adansonia.za")

## Setting basic theme options for plot with ggplot2
theme_base_3 <- theme(
  ## Axis
  axis.line=element_blank(),
  axis.text=element_blank(),
  axis.ticks=element_blank(),
  axis.title=element_blank(),
  ## Legend
  legend.position="bottom",
  legend.title=element_blank(), #add Vote for bottom species
  legend.text=element_text(size=11),
  legend.key.height=unit(0.5,"line"),
  legend.key.width=unit(1.2,"line"),
  legend.box.background=element_blank(),
  ## Plot
  plot.title=element_text(hjust=0.5,size=13),
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

## Function to plot climate change scenarios (current RCP 4.5 and RCP 8.5)
plot_anomaly_2 <- function(r, title, label="x") {
  rdf <- data.frame(rasterToPoints(r)) # To plot raster with geom_raster()
  names(rdf) <- c("x", "y", "z")
  p <- ggplot(NULL, aes(x, y)) +
    geom_tile(data=rdf, aes(fill=z)) +
    theme_bw() + theme_base_3 +
    coord_fixed(xlim=c(313000,1090000), ylim=c(7167000,8676000)) +
    annotate("text",x=600000,y=8500000,label=label,hjust=1,vjust=0,size=5) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    labs(title=title) 
    return(p)
}

## Function to plot altitude range with current points
plot_anomaly_alt <- function(m,n,title, label="x") {
    rdf <- data.frame(rasterToPoints(m)) # To plot raster with geom_raster()
    names(rdf) <- c("x", "y", "z")
    rdf2 <- as.data.frame(n)
     ang <-  ggplot(NULL, aes(x,y)) + 
      geom_tile(data=rdf, aes(fill=z)) +
      geom_point(data=rdf2, shape=1, aes(x,y), color="black", size=2)+
      theme_bw() + theme_base_3 +
      coord_fixed(xlim=c(313000,1090000), ylim=c(7167000,8676000)) +
      annotate("text",x=600000,y=8500000,label=label,hjust=1,vjust=0,size=5) +
      scale_y_continuous(expand=c(0,0)) +
      scale_x_continuous(expand=c(0,0)) +
      labs(title=title)
    return(ang)
}


## Color scales future
col_scale_var_test_fut <- scale_fill_gradientn(
  colours = c(grey(c(0.90,seq(0.7,0.50,-0.05))),gcolors(7)),
  na.value="transparent",
  values=rescale(seq(-125,3125),-125,3125),
  limits=c(-125,3125),
  breaks= seq(0,3000,by=500),
  labels= c(0,2,4,6,8,10,12)
)


## Color scales pres
col_scale_var_test_pres <- scale_fill_gradientn(
  colours = c(grey(c(0.90,seq(0.9,0.7,-0.2))),gcolors(3)),
  #colours = c(grey(c(0.90,seq(0.9,0.7,-0.2))),gcolors(3)),
  na.value="transparent",
  values=rescale(seq(-125,1125),-125,1125),
  limits=c(-125,1125),
  breaks= seq(0,1000,by=250),
  labels= c(0,1,2,3,4)
)

# altitude map presence points all species
## Presence points and altitude
# Legend specifications

col_scale_var_test_alt <- scale_fill_gradientn(
  colours = terrain.colors(255)[255:1],
  na.value="transparent",
  values=rescale(seq(0,3000),0,3000),
  limits=c(0,3000),
  breaks= seq(0,3000,by=1500),
  labels= seq(0,3000,by=1500)
)

# Ploting
# non threatened
a_za <- read.table(paste0("Adansonia.za/function_ready.txt"), header=T,sep="\t")
a_digi <- read.table(paste0("Adansonia.digitata/function_ready.txt"), header=T,sep="\t")
a_grand <- read.table(paste0("Adansonia.grandidieri/function_ready.txt"), header=T,sep="\t")
#threatened
a_mada <- read.table(paste0("Adansonia.madagascariensis/function_ready.txt"), header=T,sep="\t")
a_perri <- read.table(paste0("Adansonia.perrieri/function_ready.txt"), header=T,sep="\t")
a_rubro <- read.table(paste0("Adansonia.rubrostipa/function_ready.txt"), header=T,sep="\t")
a_suare <- read.table(paste0("Adansonia.suarezensis/function_ready.txt"), header=T,sep="\t")

# non threatened
a1 <- plot_anomaly_alt(m=environ$alt,n=a_digi,label="(a)",
                     title="Presence Points\nElevation(m)") + 
  col_scale_var_test_alt

e2 <- plot_anomaly_alt(m=environ$alt,n=a_grand,label="(e)",
                       title="") + col_scale_var_test_alt

i3 <- plot_anomaly_alt(m=environ$alt,n=a_za,label="(i)",
                       title="") + col_scale_var_test_alt

# threatened species
a1t <- plot_anomaly_alt(m=environ$alt,n=a_mada,label="(a)",
                       title="Presence Points\nElevation (m)") + 
  col_scale_var_test_alt

e2t <- plot_anomaly_alt(m=environ$alt,n=a_perri,label="(e)",
                       title="") + col_scale_var_test_alt

i3t <- plot_anomaly_alt(m=environ$alt,n=a_rubro,label="(i)",
                       title="") + col_scale_var_test_alt

m4t <- plot_anomaly_alt(m=environ$alt,n=a_suare,label="(m)",
                       title="") + col_scale_var_test_alt

# Load predictions and update extent present
# non threatened 
pred_dig <- stack(paste0("Adansonia.digitata/proj_current/proj_current_Adansonia.digitata_ensemble.grd"))
ca_dig <- pred_dig[[1]]
pred_gran <- stack(paste0("Adansonia.grandidieri/proj_current/proj_current_Adansonia.grandidieri_ensemble.grd"))
ca_gran <- pred_gran[[1]]
pred_za <- stack(paste0("Adansonia.za/proj_current/proj_current_Adansonia.za_ensemble.grd"))
ca_za <- pred_za[[1]]

# threatened 
pred_mada <- stack(paste0("Adansonia.madagascariensis/proj_current/proj_current_Adansonia.madagascariensis_ensemble.grd"))
ca_mada <- pred_mada[[1]]
pred_perri <- stack(paste0("Adansonia.perrieri/proj_current/proj_current_Adansonia.perrieri_ensemble.grd"))
ca_perri <- pred_perri[[1]]
pred_rubro <- stack(paste0("Adansonia.rubrostipa/proj_current/proj_current_Adansonia.rubrostipa_ensemble.grd"))
ca_rubro <- pred_rubro[[1]]
pred_suare <- stack(paste0("Adansonia.suarezensis/proj_current/proj_current_Adansonia.suarezensis_ensemble.grd"))
ca_suare <- pred_suare[[1]]

# non threatened
b1 <- plot_anomaly_2(r=ca_dig,label="(b)",
                     title="Current\nDistribution") + 
  col_scale_var_test_pres

f2 <- plot_anomaly_2(r=ca_gran,label="(f)",
                     title="") + col_scale_var_test_pres

j3 <- plot_anomaly_2(r=ca_za,label="(j)",
                     title="") + col_scale_var_test_pres
# threatened 
b1t <- plot_anomaly_2(r=ca_mada,label="(b)",
                     title="Current\nDistribution") + 
  col_scale_var_test_pres

f2t <- plot_anomaly_2(r=ca_perri,label="(f)",
                     title="") + col_scale_var_test_pres

j3t <- plot_anomaly_2(r=ca_rubro,label="(j)",
                     title="") + col_scale_var_test_pres

n4t <- plot_anomaly_2(r=ca_suare,label="(n)",
                      title="") + col_scale_var_test_pres

################## Future maps
# Non threatened
cafut_dig <- raster(paste0("Adansonia.digitata/caFut.tif"))
cafut_gran <- raster(paste0("Adansonia.grandidieri/caFut.tif"))
cafut_za <- raster(paste0("Adansonia.za/caFut.tif"))

caZD_dig <- raster(paste0("Adansonia.digitata/caZD.tif"))
caZD_gran <- raster(paste0("Adansonia.grandidieri/caZD.tif"))
caZD_za <- raster(paste0("Adansonia.za/caZD.tif"))

# threatened 
cafut_mada <- raster(paste0("Adansonia.madagascariensis/caFut.tif"))
cafut_perri <- raster(paste0("Adansonia.perrieri/caFut.tif"))
cafut_rubro <- raster(paste0("Adansonia.rubrostipa/caFut.tif"))
cafut_suare <- raster(paste0("Adansonia.suarezensis/caFut.tif"))

caZD_mada <- raster(paste0("Adansonia.madagascariensis/caZD.tif"))
caZD_perri <- raster(paste0("Adansonia.perrieri/caZD.tif"))
caZD_rubro <- raster(paste0("Adansonia.rubrostipa/caZD.tif"))
caZD_suare <- raster(paste0("Adansonia.suarezensis/caZD.tif"))

# non threatened 
c1 <- plot_anomaly_2(r=cafut_dig, label="(c)",
                     title="Full Dispersal\nRCP 8.5 2080") +
  col_scale_var_test_fut

g2 <- plot_anomaly_2(r=cafut_gran, label="(g)",
                     title="") + col_scale_var_test_fut

k3 <- plot_anomaly_2(r=cafut_za, label="(k)",
                     title="") + col_scale_var_test_fut

d1 <- plot_anomaly_2(r=caZD_dig, label="(d)",
                     title="Zero Dispersal\nRCP 8.5 2080") + 
  col_scale_var_test_fut

h2 <- plot_anomaly_2(r=caZD_gran, label="(h)",
                     title="") + col_scale_var_test_fut

l3 <- plot_anomaly_2(r=caZD_za, label="(l)",
                     title="") + col_scale_var_test_fut
# threatened

c1t <- plot_anomaly_2(r=cafut_mada, label="(c)",
                     title="Full Dispersal\nRCP 8.5 2080") + 
  col_scale_var_test_fut

g2t <- plot_anomaly_2(r=cafut_perri, label="(g)",
                     title="") + col_scale_var_test_fut

k3t <- plot_anomaly_2(r=cafut_rubro, label="(k)",
                     title="") + col_scale_var_test_fut

o4t <- plot_anomaly_2(r=cafut_suare, label="(o)",
                      title="") + col_scale_var_test_fut

### Zero Dispersal

d1t <- plot_anomaly_2(r=caZD_mada, label="(d)",
                     title="Zero Dispersal\nRCP 8.5 2080") + 
  col_scale_var_test_fut

h2t <- plot_anomaly_2(r=caZD_perri, label="(h)",
                     title="") + col_scale_var_test_fut

l3t <- plot_anomaly_2(r=caZD_rubro, label="(l)",
                     title="") + col_scale_var_test_fut

p4t <- plot_anomaly_2(r=caZD_suare, label="(p)",
                      title="") + col_scale_var_test_fut

## Non threatened legend name

tgrob_dig <- textGrob("A. digitata",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                    hjust=0.5, vjust=4)
tgrob_gran <- textGrob("A. grandidieri",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                       hjust=0.5, vjust=4)
tgrob_za <- textGrob("A. za",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                     hjust=0, vjust=4)


## Combine plots
lay_4 <- rbind(c(1,rep(seq(4,7,by=1),each=3)),
               c(1,rep(seq(4,7,by=1),each=3)),
               c(2,rep(seq(8,11,by=1),each=3)),
               c(2,rep(seq(8,11,by=1),each=3)),
               c(3,rep(seq(12,15,by=1),each=3)),
               c(3,rep(seq(12,15,by=1),each=3)))
               

plot_baobabs <- grid.arrange(tgrob_dig, tgrob_gran, tgrob_za,
                             a1,b1,c1,d1,e2,f2,g2,h2,
                             i3,j3,k3,l3,layout_matrix=lay_4)

#ggsave(file=paste0("./outputs/non_threat.pdf"),
       #plot=plot_baobabs,width=11,height=10,dpi="print")

ggsave(file=paste0("./outputs/non_threat.png"),
       plot=plot_baobabs,width=9,height=10,dpi="print")


## Threatened legend name
tgrob_mada <- textGrob("A. madagascariensis",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                     hjust=0.5, vjust=4)
tgrob_perri <- textGrob("A. perrieri",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                       hjust=0.5, vjust=4)
tgrob_rubro <- textGrob("A. rubrostipa",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                      hjust=0.5, vjust=4)
tgrob_suare <- textGrob("A. suarezensis",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                        hjust=0.5, vjust=4)
## Combine plots
lay_5 <- rbind(c(1,rep(seq(5,8,by=1),each=3)),
               c(1,rep(seq(5,8,by=1),each=3)),
               c(2,rep(seq(9,12,by=1),each=3)),
               c(2,rep(seq(9,12,by=1),each=3)),
               c(3,rep(seq(13,16,by=1),each=3)),
               c(3,rep(seq(13,16,by=1),each=3)),
               c(4,rep(seq(17,20,by=1),each=3)),
               c(4,rep(seq(17,20,by=1),each=3)))

plot_baobabs_threatened <- grid.arrange(tgrob_mada, tgrob_perri, tgrob_rubro,tgrob_suare,
                             a1t,b1t,c1t,d1t,e2t,f2t,g2t,h2t,
                             i3t,j3t,k3t,l3t,m4t,n4t,o4t,p4t,layout_matrix=lay_5)

#ggsave(file=paste0("./outputs/threat.pdf"),
       #plot=plot_baobabs_threatened,width=11,height=10,dpi="print")

ggsave(file=paste0("./outputs/threat.png"),
       plot=plot_baobabs_threatened,width=9,height=12,dpi="print")

############## New maps requested to new Results section
## Draw points in the SDA and extract environmental variables

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

### A digitata ###

pred_dig <- stack(paste0("Adansonia.digitata/proj_current/proj_current_Adansonia.digitata_ensemble.grd"))
ca_dig <- pred_dig[[1]]
ca_test <- stack(ca_dig,environ$alt)
ca_test_na <- na.omit(ca_test)

ca_test_na2 <- rasterToPoints(ca_test_na)
ca_test_na3 <- as.data.frame(ca_test_na2[complete.cases(ca_test_na2), ] )

alegria <- ca_test_na3 %>% filter(Adansonia.digitata_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

set.seed(20)
alegria2 <- sample_n(alegria, size = 1000, replace = F)
alegria2$Proj <- rep(c("Present"),1000)
names(alegria2) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_dig <- raster(paste0("Adansonia.digitata/caFD_85_2080.tif"))
ca_test_fut <- stack(caZD_dig,environ$alt)
ca_test_na_fut <- na.omit(ca_test_fut)

ca_test_na2_fut <- rasterToPoints(ca_test_na_fut)
ca_test_na3_fut <- as.data.frame(ca_test_na2_fut[complete.cases(ca_test_na2_fut), ] )
alegria_fut <- ca_test_na3_fut %>% filter(caFD_85_2080 >= 1500)

set.seed(20)
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

set.seed(20)
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
  annotate("text", x  = 0.5, y = 600 , size=7, label = "(a)") +
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

set.seed(20)
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

set.seed(20)
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


set.seed(20)
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

######## A madagascariensis ############

pred_mada <- stack(paste0("Adansonia.madagascariensis/proj_current/proj_current_Adansonia.madagascariensis_ensemble.grd"))
ca_mada <- pred_mada[[1]]
ca_test_mada <- stack(ca_mada,environ$alt)
ca_test_na_mada <- na.omit(ca_test_mada)

ca_test_na2_mada <- rasterToPoints(ca_test_na_mada)
ca_test_na3_mada <- as.data.frame(ca_test_na2_mada[complete.cases(ca_test_na2_mada), ] )

alegria_mada <- ca_test_na3_mada %>% filter(Adansonia.madagascariensis_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

set.seed(20)
alegria2_mada <- sample_n(alegria_mada, size = 1000, replace = F)
alegria2_mada$Proj <- rep(c("Present"),1000)
names(alegria2_mada) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_mada <- raster(paste0("Adansonia.madagascariensis/caZD_85_2080.tif"))
ca_test_fut_mada <- stack(caZD_mada,environ$alt)
ca_test_na_fut_mada <- na.omit(ca_test_fut_mada)

ca_test_na2_fut_mada <- rasterToPoints(ca_test_na_fut_mada)
ca_test_na3_fut_mada <- as.data.frame(ca_test_na2_fut_mada[complete.cases(ca_test_na2_fut_mada), ] )
alegria_fut_mada <- ca_test_na3_fut_mada %>% filter(caZD_85_2080 >= 1500)

set.seed(20)
alegria2_fut_mada <- sample_n(alegria_fut_mada, size = 1000, replace = F)
alegria2_fut_mada$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_mada) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_mada <- rbind(alegria2_mada,alegria2_fut_mada)
tail(alegria_test_mada)

### Future testing 2050 
caZD_mada_2055 <- raster(paste0("Adansonia.madagascariensis/caZD_85_2050.tif"))
ca_test_fut_mada_2055 <- stack(caZD_mada_2055,environ$alt)
ca_test_na_fut_mada_55 <- na.omit(ca_test_fut_mada_2055)

ca_test_na2_fut_mada_55 <- rasterToPoints(ca_test_na_fut_mada_55)
ca_test_na3_fut_mada_55 <- as.data.frame(ca_test_na2_fut_mada_55[complete.cases(ca_test_na2_fut_mada_55), ] )
alegria_fut_mada_55 <- ca_test_na3_fut_mada_55 %>% filter(caZD_85_2050 >= 1500)


set.seed(20)
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
  annotate("text", x  = 0.5, y = 450 , size=7, label = "(a)") +
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

set.seed(20)
alegria2_perrieri <- sample_n(alegria_perrieri, size = 1000, replace = F)
alegria2_perrieri$Proj <- rep(c("Present"),1000)
names(alegria2_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 (RCP 4.5!!!) Atention!
caZD_perrieri <- raster(paste0("Adansonia.perrieri/caZD_45_2080.tif"))
ca_test_fut_perrieri <- stack(caZD_perrieri,environ$alt)
ca_test_na_fut_perrieri <- na.omit(ca_test_fut_perrieri)

ca_test_na2_fut_perrieri <- rasterToPoints(ca_test_na_fut_perrieri)
ca_test_na3_fut_perrieri <- as.data.frame(ca_test_na2_fut_perrieri[complete.cases(ca_test_na2_fut_perrieri), ] )
alegria_fut_perrieri <- ca_test_na3_fut_perrieri %>% filter(caZD_45_2080 >= 1500)

set.seed(20)
alegria2_fut_perrieri <- sample_n(alegria_fut_perrieri, size = 427, replace = F)
alegria2_fut_perrieri$Proj <- rep(c("Future_2085"),427)
names(alegria2_fut_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_perrieri <- rbind(alegria2_perrieri,alegria2_fut_perrieri)
tail(alegria_test_perrieri)

### Future testing 2050 
caZD_perrieri_2055 <- raster(paste0("Adansonia.perrieri/caZD_85_2050.tif"))
ca_test_fut_perrieri_2055 <- stack(caZD_perrieri_2055,environ$alt)
ca_test_na_fut_perrieri_55 <- na.omit(ca_test_fut_perrieri_2055)

ca_test_na2_fut_perrieri_55 <- rasterToPoints(ca_test_na_fut_perrieri_55)
ca_test_na3_fut_perrieri_55 <- as.data.frame(ca_test_na2_fut_perrieri_55[complete.cases(ca_test_na2_fut_perrieri_55), ] )
alegria_fut_perrieri_55 <- ca_test_na3_fut_perrieri_55 %>% filter(caZD_85_2050 >= 1500)


set.seed(20)
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

set.seed(20)
alegria2_rubro <- sample_n(alegria_rubro, size = 1000, replace = F)
alegria2_rubro$Proj <- rep(c("Present"),1000)
names(alegria2_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_rubro <- raster(paste0("Adansonia.rubrostipa/caZD_85_2080.tif"))
ca_test_fut_rubro <- stack(caZD_rubro,environ$alt)
ca_test_na_fut_rubro <- na.omit(ca_test_fut_rubro)

ca_test_na2_fut_rubro <- rasterToPoints(ca_test_na_fut_rubro)
ca_test_na3_fut_rubro <- as.data.frame(ca_test_na2_fut_rubro[complete.cases(ca_test_na2_fut_rubro), ] )
alegria_fut_rubro <- ca_test_na3_fut_rubro %>% filter(caZD_85_2080 >= 1500)

set.seed(20)
alegria2_fut_rubro <- sample_n(alegria_fut_rubro, size = 1000, replace = F)
alegria2_fut_rubro$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_rubro <- rbind(alegria2_rubro,alegria2_fut_rubro)

### Future testing 2050 
caZD_rubro_2055 <- raster(paste0("Adansonia.rubrostipa/caZD_85_2050.tif"))
ca_test_fut_rubro_2055 <- stack(caZD_rubro_2055,environ$alt)
ca_test_na_fut_rubro_55 <- na.omit(ca_test_fut_rubro_2055)

ca_test_na2_fut_rubro_55 <- rasterToPoints(ca_test_na_fut_rubro_55)
ca_test_na3_fut_rubro_55 <- as.data.frame(ca_test_na2_fut_rubro_55[complete.cases(ca_test_na2_fut_rubro_55), ] )
alegria_fut_rubro_55 <- ca_test_na3_fut_rubro_55 %>% filter(caZD_85_2050 >= 1500)


set.seed(20)
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
  annotate("text", x  = 0.5, y = 325 , size=7, label = "(e)") +
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

set.seed(20)
alegria2_suare <- sample_n(alegria_suare, size = 1000, replace = F)
alegria2_suare$Proj <- rep(c("Present"),1000)
names(alegria2_suare) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 Atention RCP 4.5 (SAME AS A. perrieri)
caZD_suare <- raster(paste0("Adansonia.suarezensis/caZD_45_2080.tif"))
ca_test_fut_suare <- stack(caZD_suare,environ$alt)
ca_test_na_fut_suare <- na.omit(ca_test_fut_suare)

ca_test_na2_fut_suare <- rasterToPoints(ca_test_na_fut_suare)
ca_test_na3_fut_suare <- as.data.frame(ca_test_na2_fut_suare[complete.cases(ca_test_na2_fut_suare), ] )
alegria_fut_suare <- ca_test_na3_fut_suare %>% filter(caZD_45_2080 >= 1500)

set.seed(20)
alegria2_fut_suare <- sample_n(alegria_fut_suare, size = 100, replace = F)
alegria2_fut_suare$Proj <- rep(c("Future_2085"),100)
names(alegria2_fut_suare) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_suare <- rbind(alegria2_suare,alegria2_fut_suare)

### Future testing 2050 
caZD_suare_2055 <- raster(paste0("Adansonia.suarezensis/caZD_45_2050.tif"))
ca_test_fut_suare_2055 <- stack(caZD_suare_2055,environ$alt)
ca_test_na_fut_suare_55 <- na.omit(ca_test_fut_suare_2055)

ca_test_na2_fut_suare_55 <- rasterToPoints(ca_test_na_fut_suare_55)
ca_test_na3_fut_suare_55 <- as.data.frame(ca_test_na2_fut_suare_55[complete.cases(ca_test_na2_fut_suare_55), ] )
alegria_fut_suare_55 <- ca_test_na3_fut_suare_55 %>% filter(caZD_45_2050 >= 1500)

set.seed(20)
alegria2_fut_suare_55 <- sample_n(alegria_fut_suare_55, size = 13, replace = F)
alegria2_fut_suare_55$Proj <- rep(c("Future_2055"),13)
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

######## A za ###############

pred_za <- stack(paste0("Adansonia.za/proj_current/proj_current_Adansonia.za_ensemble.grd"))
ca_za <- pred_za[[1]]
ca_test_za <- stack(ca_za,environ$alt)
ca_test_na_za <- na.omit(ca_test_za)

ca_test_na2_za <- rasterToPoints(ca_test_na_za)
ca_test_na3_za <- as.data.frame(ca_test_na2_za[complete.cases(ca_test_na2_za), ] )

alegria_za <- ca_test_na3_za %>% filter(Adansonia.za_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

set.seed(20)
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

set.seed(20)
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

set.seed(20)
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

## Legend names for plot
grob_digitata <- textGrob("A. digitata",rot=90, gp=gpar(cex=2,fontface="italic"),
                          hjust=0.3, vjust=5)
grob_grand <- textGrob("A. grandidieri",rot=90, gp=gpar(cex=2,fontface="italic"),
                          hjust=0.3, vjust=5)
grob_mada <- textGrob("A. madagascariensis",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.4, vjust=5)
grob_perri <- textGrob("A. perrieri",rot=90, gp=gpar(cex=2,fontface="italic"),
                     hjust=0.3, vjust=5)
grob_rubro <- textGrob("A. rubrostipa",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=5)
grob_suare <- textGrob("A. suarezensis",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=5)
grob_za <- textGrob("A. za",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=5)
## Combine plots
### Threatened species
day_1 <- rbind(c(1,rep(seq(5,6,by=1),each=3)),
              c(2,rep(seq(7,8,by=1),each=3)),
              c(3,rep(seq(9,10,by=1),each=3)),
              c(4,rep(seq(11,12,by=1),each=3)))

### Non-threatened species
day_2 <- rbind(c(1,rep(seq(5,6,by=1),each=3)),
               c(2,rep(seq(7,8,by=1),each=3)),
               c(3,rep(seq(9,10,by=1),each=3)))             

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

ggsave(file=paste0("./outputs/lat_ele_all_species_within_SDAcf.png"),
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

##===========================================================================
## End of script - Have a Nice Day and Enjoy The View and go for a hike!
##===========================================================================
