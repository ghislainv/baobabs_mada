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
df.sp <- SpatialPointsDataFrame(coords,data=df.orig,proj4string=CRS("+init=epsg:4326"))
## Reproject into UTM 38S
df.sp <- spTransform(df.sp,CRS("+init=epsg:32738"))
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
p1 <- plot_anomaly(r=var_85_2080[["tmean"]], label="a",
                   title="Temperature\n(°C x 10)") +
    col_scale_var_tmean
p2 <- plot_anomaly(r=ano_85_2080[["tmean"]], label="a'",
                   title="") +
    col_scale_ano_tmean
p3 <- plot_anomaly(r=var_85_2080[["tseas"]], label="b",
                   title="T. seasonality\n(°C sd x 100)") +
    col_scale_var_tseas
p4 <- plot_anomaly(r=ano_85_2080[["tseas"]], label="b'",
                   title="") +
    col_scale_ano_tseas
## prec
p5 <- plot_anomaly(r=var_85_2080[["prec"]], label="c",
                   title="Precipitation\n(mm/y)") +
    col_scale_var_prec
p6 <- plot_anomaly(r=ano_85_2080[["prec"]], label="c'",
                   title="") +
    col_scale_ano_prec
## cwd
p7 <- plot_anomaly(r=var_85_2080[["cwd"]], label="d",
                   title="Climatic water deficit\n(mm/y)") +
    col_scale_var_cwd
p8 <- plot_anomaly(r=ano_85_2080[["cwd"]], label="d'",
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
  head(mapmat.f)
  mapmat.f <- mapmat.final[,c(1,2,3,4,6,7,8,9,10,11)]
  write.csv2(mapmat.final,paste0("./",sp.dir[i],"/","niche_graph_species_compared_anomaly.csv"))
  #### Future niche over current SDA!!! NO ANOMALY
  wC.future <- which(values(ca)>=500)
  mapmat.df.future <- as.data.frame(var_85_2080)[wC.future,]
  names(mapmat.df.future) <- c("tmeanf","tseasf","precf","cwdf")
  Mean_ok_future_fut <- round(apply(mapmat.df.anomalies,2,mean,na.rm=TRUE))
  q_ok_future_fut <- round(apply(mapmat.df.anomalies,2,quantile,c(0.025,0.975),na.rm=TRUE))
  niche_ok_future <- as.data.frame(rbind(Mean_ok_future_fut,q_ok_future_fut))
  write.csv2(niche_ok_future,paste0("./",sp.dir[i],"/","mean_niche_with_future.csv"))
  
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
quantiles <- data_teste[,c(2,3,4,5,6)]
quantiles <- round(apply(quantiles,2,quantile,c(0.025,0.975),na.rm=TRUE)) #use it to set new plot
head(quantiles)

# Density plots ## Temp. Seasonality
range(data_teste$tseas)
# generate break positions
#breaks = c(848,1200,1600,2000,2400,2800,3309)
breaks <- c(round(seq(min(data_teste$tseas),max(data_teste$tseas),length=6)))

# and labels
labels = as.character(breaks) # labels must be the same as breaks, otherwise error
# plot the map
my_plot = ggplot(data_teste, aes(x=tseas, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + #specific palette
  geom_vline(xintercept = range(data_teste$tseas), show.legend = T, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$tseas), breaks = breaks, labels = labels)
my_plot_tseas =  my_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 #legend.text = element_text(face= "italic",size=15),
                                 legend.title=element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100ºC)", y = "Density",size=5) +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.5),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.5), colour="black")) +
  theme(legend.position="none")

# Density plots ### Annual Mean Temperature
range(data_teste$tmean)
# generate break positions
breaks <- c(round(seq(min(data_teste$tmean),max(data_teste$tmean),length=6)))
labels = as.character(breaks)
my_plot2 = ggplot(data_teste, aes(x=tmean, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = range(data_teste$tmean), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$tmean), breaks = breaks, labels = labels)
my_plot_tmean =  my_plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  #legend.text = element_text(face= "italic",size=15),
                                  legend.title=element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10) ", y = "Density",size=5) +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.5),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.5), colour="black")) +
  theme(legend.position="none")

# Density plots ### Annual Mean Precipitation
range(data_teste$prec)
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
  theme(axis.text.x = element_text(size = rel(1.5),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.5), colour="black")) +
  theme(legend.position="none")

# Density plots ### Climatic Water Deficit
range(data_teste$cwd)
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
  theme(axis.text.x = element_text(size = rel(1.5),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.5), colour="black")) +
  theme(legend.position="none")
  
## Combine plots 
lay_2 <- cbind(c(rep(seq(1,4,by=1),each=3)),
               c(rep(seq(2,4,by=2),each=3)))
colnames(lay_2) <-c("a","b")
lay_2 <- subset( lay_2, select = -b )

plot_densities <- grid.arrange(my_plot_tmean, my_plot_tseas, my_plot_prec, my_plot_cwd,
                               layout_matrix=lay_2)
ggsave(file=paste0("./outputs/climaticwd_species.pdf"),plot=my_plot4,width=20,height=5)

### Comparing bioclimatic niche of each species inside current SDA according to the 
### most important variables ######
# A_suarezensis
# Import data_set
a.suare <- read.csv(file=paste0("Adansonia.suarezensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
# 1st importance variable - Tseas
range(a.suare$tseas) # 1194 - 1417 
range(a.suare$tseasf) # 1314- 1452
mean(a.suare$tseas) # 1275,3
mean(a.suare$tseasf) # 1374
breaks = c(1150,1200,1250,1300,1350,1400,1450,1500)
labels = as.character(breaks)
plot.seas = ggplot(a.suare, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(1175, 1450), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(1150, 1500), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas =  plot.seas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic",size=15),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/adan.suare_current_future_niche_comparison_in_sda_v1.pdf"),
       plot=plot.seas,width=20,height=5)

### second importance variable - Prec
range(a.suare$prec) # 1113 - 1516
range(a.suare$precf) # 922.3 1312.3
breaks = c(850,900,1000,1100,1200,1300,1400,1500,1600)
labels = as.character(breaks)
plot.precv2= ggplot(a.suare, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(830, 1600), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(850, 1600), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.precv2 =  plot.precv2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_text(face= "italic", size=15),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean annual precipitation (mm.y-¹)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/adan.suare_current_future_niche_comparison_in_sda_42.pdf"),
       plot=plot.precv2,width=20,height=5)

################################
# Import data_set A perrieri
#1st importance variable - Tseas 
a.perrieri <- read.csv(file=paste0("Adansonia.perrieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.perrieri)
range(a.perrieri$tseas) # 890 - 1966
range(a.perrieri$tseasf) # 1035 - 2165
mean(a.perrieri$tseas) # 1391,5
mean(a.perrieri$tseasf) # 11462,9
breaks = c(800,1000,1200,1400,1600,1800,2000,2250)
labels = as.character(breaks)
plot.seas = ggplot(a.perrieri, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(800, 1950), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(800, 2250), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas =  plot.seas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic",size=15),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.perrieri_current_future_niche_comparison_in_sdap.pdf"),
       plot=plot.seas,width=20,height=5)

## 2nd Most IV A perrieri - climatic water deficit
range(a.perrieri$cwd) # 179 - 781
range(a.perrieri$cwdf) # 372 - 1529
mean(a.perrieri$cwd) # 467
mean(a.perrieri$cwdf) # 1014
breaks = c(100,400,700,1000,1300,1600)
labels = as.character(breaks)
plot.cwd = ggplot(a.perrieri, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(100, 1600), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.cwd2 =  plot.cwd+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            legend.text = element_text(face= "italic",size=15),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=15) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.perrieri_current_future_niche_comparison_in_sdap_v2.pdf"),
       plot=plot.cwd2,width=20,height=5)

############## A.rubrostipa
# Import dataset
a.rubrostipa <- read.csv(file=paste0("Adansonia.rubrostipa/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.rubrostipa)
## 1st importance variable - Climatic Water Deficit
range(a.rubrostipa$cwd) # 684 - 967
range(a.rubrostipa$cwdf) # 1373 - 1979
mean(a.rubrostipa$cwd) # 808,8
mean(a.rubrostipa$cwdf) # 1688,6
breaks = c(650,810,970,1370,1690,2000)
labels = as.character(breaks)

plot.cwd = ggplot(a.rubrostipa, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(650, 2000), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.cwd =  plot.cwd + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic",size=15),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.rubrostipa_current_future_niche_comparison_in_sdapv1.pdf"),
       plot=plot.cwd,width=20,height=5)

### 2nd Importance Variable - Tmean
range(a.rubrostipa$tmean) # 233 272
range(a.rubrostipa$tmeanf) # 264.3 305,6
mean(a.rubrostipa$tmean) # 257.8
mean(a.rubrostipa$tmeanf) # 291.3
breaks = c(240,260,280,300)
labels = as.character(breaks)
plot.tmean2= ggplot(a.rubrostipa, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(220, 320), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.tmean2 =  plot.tmean2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_text(face= "italic",size=15),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.rubrostipa_current_future_niche_comparison_in_sda_v2.pdf"),
       plot=plot.tmean2,width=20,height=5)

############## A. madagascariensis 
## Import dataset
a.madagascariensis <- read.csv(file=paste0("Adansonia.madagascariensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.madagascariensis)
tail(a.madagascariensis)
# 1st importance variable - Tseas
range(a.madagascariensis$tseas) # 913 - 2066
range(a.madagascariensis$tseasf) # 1058 - 2066,3
mean(a.madagascariensis$tseas) # 1497
mean(a.madagascariensis$tseasf) # 1585
breaks = c(900,1200,1500,1800,2100)
labels = as.character(breaks)

plot.seas2 = ggplot(a.madagascariensis, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(800,2200), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas2 =  plot.seas2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic",size=15),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (ºC sd x 100)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.mada_current_future_niche_comparison_in_sdapv1.pdf"),
       plot=plot.seas2,width=20,height=5)

### 2nd importance variable - Tmean
range(a.madagascariensis$tmean) # 248 - 275
range(a.madagascariensis$tmeanf) # 278 - 311.6
mean(a.madagascariensis$tmean) # 263.1
mean(a.madagascariensis$tmeanf) # 297.2
breaks = c(250,265,280,295,310)
labels = as.character(breaks)

plot.mean2 = ggplot(a.madagascariensis, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(235, 330), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(235,320), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.mean2 =  plot.mean2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic",size=15),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10) ", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.mada_current_future_niche_comparison_in_sdapv2.pdf"),
       plot=plot.mean2,width=20,height=5)

############## A. grandidieri
## Import dataset
a.grandidieri <- read.csv(file=paste0("Adansonia.grandidieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.grandidieri)
# 1st Importance Variable - Prec
range(a.grandidieri$prec) # 383 - 1029
range(a.grandidieri$precf) # 366 - 920
mean(a.grandidieri$prec) # 750
mean(a.grandidieri$precf) # 711
breaks = c(350,600,750,900,1050)
labels = as.character(breaks)

plot.prec = ggplot(a.grandidieri, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(300,1100), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
  plot.prec =  plot.prec + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic",size=15),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Precipitation (mm.y-¹)", y = "Density",size=10) +
  labs(col = "") +
    theme(axis.title.x = element_text(size = rel(3))) +
    theme(axis.title.y = element_text(size = rel(3))) +
    theme(axis.text.x = element_text(size = rel(3))) +
    theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.grandi_current_future_niche_comparison_in_sdapv1.pdf"),
       plot=plot.prec,width=20,height=5)

## 2nd Importance Variable - Tmean
range(a.grandidieri$tmean) # 241 - 263
range(a.grandidieri$tmeanf) # 276 297
mean(a.grandidieri$tmean) # 250,7
mean(a.grandidieri$tmeanf) # 286,1
breaks = c(240,250,260,280,290,300)
labels = as.character(breaks)

plot.mean2 = ggplot(a.grandidieri, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(230, 310), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.mean2 =  plot.mean2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic",size=15),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10) ", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.grand_current_future_niche_comparison_in_sdapv2.pdf"),
       plot=plot.mean2,width=20,height=5)

############## A za
## Import dataset
a.za <- read.csv(file=paste0("Adansonia.za/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.za)
# 1st Importance Variable - Prec
range(a.za$prec) # 355 1635
range(a.za$precf) # 334,1 - 1320,3
mean(a.za$prec) # 823
mean(a.za$precf) # 642,1
breaks = c(300,650,850,1200,1650)
labels = as.character(breaks)

plot.prec = ggplot(a.za, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(200,1800), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.prec =  plot.prec + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic",size=15),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Precipitation (mm.y-¹)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.za_current_future_niche_comparison_in_sdapv1.pdf"),
       plot=plot.prec,width=20,height=5)

## 2nd Importance Variable - Tmean
range(a.za$tmean) # 214 275
range(a.za$tmeanf) # 252 307
mean(a.za$tmean) # 245,2
mean(a.za$tmeanf) # 276,9
breaks = c(210,240,270,310)
labels = as.character(breaks)

plot.mean2 = ggplot(a.za, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(205, 320), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(200, 320), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.mean2 =  plot.mean2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic",size=15),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10) ", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.za_current_future_niche_comparison_in_sdapv2.pdf"),
       plot=plot.mean2,width=20,height=5)

##############A digitata 
## Import dataset
a.digitata <- read.csv(file=paste0("Adansonia.digitata/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.digitata)
## 1st Importance Variable - Climatic Water Deficit
range(a.digitata$cwd) # 622 - 920
range(a.digitata$cwdf) # 1229 - 1936,3
mean(a.digitata$cwd) # 822,4
mean(a.digitata$cwdf) # 1737,32
breaks = c(600,820,1000,1230,1740,1940)
labels = as.character(breaks)

plot.cwd = ggplot(a.digitata, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(650, 1950), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(600, 2000), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.cwd =  plot.cwd + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic",size=15),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.digitata_current_future_niche_comparison_in_sdapv1.pdf"),
       plot=plot.cwd,width=20,height=5)

##2nd Importance Variable - Tseas

# range tseas 
range(a.digitata$tseas) # 848 2482
range(a.digitata$tseasf) # 1023,33 2035
mean(a.digitata$tseas) # 1531
mean(a.digitata$tseasf) # 1572,4
breaks = c(840,1100,1400,1700,2050)
labels = as.character(breaks)

plot.seas2 = ggplot(a.digitata, na.rm=T, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5, na.rm=T) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(750,2600), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas2 =  plot.seas2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic",size=15),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (ºC  sd x 100)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.digi_current_future_niche_comparison_in_sdapv2.pdf"),
       plot=plot.seas2,width=20,height=5)

#######################################################################
### World Seasonality Map
#######################################################################


library(maptools)
data("wrld_simpl")

r <- raster::getData("worldclim",var="bio",res=10)
r <- r[[c(4)]]
names(r) <- c("Seas")

bio4.gs.2080_world <- raster("data/global_future_climate/gs85bi704.tif")
bio4.he.2080_world <- raster("data/global_future_climate/he85bi704.tif")
bio4.no.2080_world <- raster("data/global_future_climate/no85bi704.tif")
Stack.bio4.2080_world <- stack(c(bio4.gs.2080_world,bio4.he.2080_world,
                                 bio4.no.2080_world))
bio4.2080_world <- mean(Stack.bio4.2080_world) ### 2080 temp seas. raster
crs(bio4.2080_world) <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(wrld_simpl) <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Plot present

bio4.anomalies.pres_world <- r
bio4.anomalies.pres_world[] <- bio4.2080_world[]-r[]

pdf("./outputs/tseas_world_3_charts.pdf",width=15,height=10) # all mapas
par(mfrow=c(3,1))
plot(r,col=viridis_pal(option ="E")(255),xlim=c(-90,160),
     ylim=c(-23.5,23.5),
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,legend=T, horizontal=T, 
     #breaks=breakpoints2,axis.args=a.arg2, zlim=c(800,3400),
     main="Temperature seasonality (sd x 100)")
plot(wrld_simpl, add=T)
abline(h=0,lty=1,col="black") #equator solid line
abline(h=-23.5, lty=2, col="black")
abline(h=23.5, lty=2, col="black")

## Anomaly future
breakpoints <- c(-640,-600,-560,-520,-480,-440,-400,-360,-320,-280,
                 -240,-200,-160,-120,-80,-40,-30,-20,-10,
                 0,
                 55,110,165,210,265,320,375,430,485,540,
                 595,650,705,760,815,870,925,980,1035,1090,
                 1145,1200)
colors <- c("#00204DFF", "#002251FF","#002860FF","#002F6FFF","#00336FFF","#07366EFF",
            "#1D3B6DFF", "#2A406CFF","#34456BFF","#3C4A6BFF","#444F6BFF","#4C546CFF",
            "#52596CFF", "#5F636EFF","#656870FF","#888579FF","#807E79FF","#848279FF",
            "#868379FF",
            "gray70",
            "#BFB06EFF","#C2B36DFF","#C6B66BFF","#CCBA69FF","#D0BE67FF","#D5C364FF",
            "#DDC95FFF","#E0CB5EFF","#E2CD5CFF","#E4CF5BFF","#E7D259FF","#EAD357FF",
            "#ECD555FF","#F0D852FF","#F2DA50FF","#F4DC4EFF","#F7DF4BFF","#FAE149FF",
            "#FDE346FF","#FFE544FF","#FFE742FF","#FFEA46FF")

a.arg <- list(at=seq(-800,1200,length.out=11), 
              labels=c("-800","-600","-400","-200","0","200",
                       "400","600","800","1000","1200"))
#pdf("./outputs/seas_world_anomaly.pdf",width=5,height=20) # all mapas

plot(bio4.anomalies.pres_world,col=colors, 
     breaks=breakpoints, axis.arg=a.arg, zlim=c(-640,1200),xlim=c(-90,160),
     ylim=c(-23.5,23.5),
     axes=FALSE,box=FALSE,legend=T,
     main="Future Climatic Anomaly",horizontal=TRUE)
plot(wrld_simpl, add=T)
abline(h=0,lty=1,col="black") #equator solid line
abline(h=-23.5, lty=2, col="black")
abline(h=23.5, lty=2, col="black")


### 2080 climate

plot(bio4.2080_world,col=viridis_pal(option ="E")(255),xlim=c(-90,160),
     ylim=c(-23.5,23.5),
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,legend=T, horizontal=T, 
     #breaks=breakpoints,axis.args=a.arg, #zlim=c(800,3400),
     main="Future Temperature Seasonality (sd x 100)")
plot(wrld_simpl, add=T)
abline(h=0,lty=1,col="black") #equator solid line
abline(h=-23.5, lty=2, col="black")
abline(h=23.5, lty=2, col="black")

dev.off()

##===========================================================================
## End of script    
##===========================================================================
