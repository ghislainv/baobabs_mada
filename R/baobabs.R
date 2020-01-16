#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com / mario.tagliari@posgrad.ufsc.br
# license         :GPLv3
# ==============================================================================

## Set environmental variable
# For MAXENT.Phillips with JAVA to work on RStudio server
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
# Note: E(X+Y) = E(X) + E(Y) so, mean of differences = difference of the means
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
  ano_85_2080 <- addLayer(ano_85_2080, mean(v_85_2080)-current[[v]])
}
var_name <- c("tmean", "tseas", "prec", "cwd")
names(var_85_2080) <- var_name
names(ano_85_2080) <- var_name

# Plot anomalies
pdf(file="outputs/Fig1_climatic_anomalies.pdf")

# Setting basic theme options for plot with ggplot2
theme_base <- theme(axis.line=element_blank(),
										axis.text.x=element_blank(),
										axis.text.y=element_blank(),
										axis.ticks=element_blank(),
										axis.title.x=element_blank(),
										axis.title.y=element_blank(),
										legend.position="bottom",
										legend.title=element_blank(),
										legend.text=element_text(size=12),
										legend.key.height=unit(0.5,"cm"),
										legend.key.width=unit(2,"cm"),
										legend.spacing=grid::unit(c(0,0,0,0),"lines"),
										plot.title=element_text(hjust=0.5),
										plot.margin=grid::unit(c(0,0,0,0),"lines"),
										panel.spacing=grid::unit(c(0,0,0,0),"null"),
										plot.background=element_rect(fill="transparent"),
										panel.background=element_rect(fill="transparent"),
										panel.grid.major=element_blank(),
										panel.grid.minor=element_blank(),
										panel.border=element_blank())

# Function to plot anomalies
plot_anomaly <- function(r, title, viridis=list(option="A", direction=1)) {
  rdf <- data.frame(rasterToPoints(r)) # To plot raster with geom_raster()
  names(rdf) <- c("x", "y", "z")
  p <- ggplot(NULL, aes(x, y)) +
    geom_raster(data=rdf, aes(fill=z)) +
    scale_fill_viridis_c(option=viridis$option, direction=viridis$direction) +
    theme_bw() + theme_base + coord_fixed() +
    labs(title=title)
  return(p)
}

p1 <- plot_anomaly(r=var_85_2080[["tmean"]],
                  viridis=list(option="B", direction=1),
                  title="Annual mean temp.\n(°C x 10)") +
      annotate("text",x=-Inf,y=Inf,label="a",hjust=0,vjust=1,size=10,fontface="bold")
p2 <- plot_anomaly(r=var_85_2080[["tmean"]],
                  viridis=list(option="B", direction=1),
                  title="Annual mean temp.\n(°C x 10)")
grid.arrange


# tmean
plot(var_85_2080[["tmean"]], col=viridis_pal(option="B")(255),
     axes=FALSE, box=FALSE, legend=TRUE, horizontal=TRUE, 
     main="Annual mean temp.\n(°C x 10)")
plot(ano_85_2080[["tmean"]], col=viridis_pal(option="B")(255),
     axes=FALSE, box=FALSE, legend=TRUE,
     main="Future anomaly\n", horizontal=TRUE)
# tseas
plot(var_85_2080[["tseas"]], col=viridis_pal(option="D")(255),
     axes=FALSE, box=FALSE, legend=TRUE, horizontal=TRUE, 
     main="Temp. seas.\n(°C sd x 100)")
plot(ano_85_2080[["tseas"]], col=viridis_pal(option="D")(255),
     axes=FALSE, box=FALSE, legend=TRUE,
     main="Future anomaly\n", horizontal=TRUE)
# precip
plot(var_85_2080[["prec"]], col=viridis_pal(option="E", direction=-1)(255),
     axes=FALSE, box=FALSE, legend=TRUE, horizontal=TRUE, 
     main="Annual precipitation\n(mm/y)")
plot(ano_85_2080[["prec"]], col=viridis_pal(option="E", direction=-1)(255),
     axes=FALSE, box=FALSE, legend=TRUE,
     main="Future anomaly\n", horizontal=TRUE)
# cwd
plot(var_85_2080[["cwd"]], col=viridis_pal(option="C")(255),
     axes=FALSE, box=FALSE, legend=TRUE, horizontal=TRUE, 
     main="Climatic water deficit\n(mm/y)")
plot(ano_85_2080[["cwd"]], col=viridis_pal(option="C")(255),
     axes=FALSE, box=FALSE, legend=TRUE,
     main="Future anomaly\n", horizontal=TRUE)
dev.off()

##==================================
## Climate change in SDA
##==================================

## Draw points in the SDA and extract future environmental variables
# In SDA
getwd()
setwd("./baobabs_mada")
sp.names <- c("Adansonia.digitata","Adansonia.grandidieri","Adansonia.madagascariensis","Adansonia.perrieri", 
              "Adansonia.rubrostipa","Adansonia.suarezensis","Adansonia.za")

for (l in 1: length(sp.names)) {
  
  
  pred <- stack(paste0("./",sp.names[l],"/proj_current/proj_current_",sp.names[l],"_ensemble.grd"))
  ca <- pred[[1]]
  
  wC.anomalies <- which(values(ca)>500) 
  nC.anomalies <- length(wC.anomalies)
  Samp.anomalies <- if (nC.anomalies>1000) {sample(wC.anomalies,1000,replace=FALSE)} else {wC.anomalies}
  mapmat.df.anomalies <- as.data.frame(anomalies_sf)[Samp.anomalies,] 
  
  ## table to compare current and future bioclimatic changes for each species 
  # Used this table to create density niche curves to compare future baobabs niche
  ## over current distribution (ca)
  
  names(mapmat.df.anomalies) <- c("tmeanf","tseasf","precf","cwdf")
  mapmat.df <- read.csv(file=paste0("./",sp.names[l],"/","niche_graph_species.csv"),sep=";")
  head(mapmat.df)
  mapmat.final <- cbind(mapmat.df.anomalies,mapmat.df)
  head(mapmat.final)
  str(mapmat.final)
  mapmat.f <- mapmat.final[,c(1,5,2,6,3,7,4,8,9,10)]
  write.csv2(mapmat.final,paste0("./",sp.names[l],"/","niche_graph_species_compared_anomaly.csv"))
  
  #### Future niche over current SDA!!! NO ANOMALY
  
  wC.future <- which(values(ca)>500)
  mapmat.df.future <- as.data.frame(anomalies_sf)[wC.future,]
  names(mapmat.df.future) <- c("tmeanf","tseasf","precf","cwdf")
  Mean_ok_future_fut <- round(apply(mapmat.df.anomalies,2,mean,na.rm=TRUE))
  q_ok_future_fut <- round(apply(mapmat.df.anomalies,2,quantile,c(0.025,0.975),na.rm=TRUE))
  niche_ok_future <- as.data.frame(rbind(Mean_ok_future_fut,q_ok_future_fut))
  write.csv2(niche_ok_future,paste0("./",sp.names[l],"/","mean_niche_with_future.csv"))
  
}


###################################################################################
######## Generating variables histograms graphs##################
###############################################################################


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


# Density plots ## Temp. Seasonality
range(data_teste$tseas)
# generate break positions
breaks = c(847,1200,1600,2000,2400,2800,3298)
# and labels
labels = as.character(breaks)

# plot

my_plot = ggplot(data_teste, aes(x=tseas, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(861, 3327), show.legend = T, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(861, 3327), breaks = breaks, labels = labels)

my_plot =  my_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           legend.text = element_text(face= "italic", size=15),
                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100ºC)", y = "Density", size=10) +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/temp_seas_species.pdf"),plot=my_plot,width=20,height=5)


# Density plots ### Annual Mean Temperature
range(data_teste$tmean)
# generate break positions
breaks = c(180,200,220,240,260,275)
# and labels
labels = as.character(breaks)
my_plot2 = ggplot(data_teste, aes(x=tmean, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(180, 275), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(180, 275), breaks = breaks, labels = labels)
my_plot2 =  my_plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic",size=15),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10) ", y = "Density",size=10) +
  labs(colour= "Species") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/annual_mean_temp_species.pdf"),plot=my_plot2,width=20,height=5)

# Density plots ### Precipitation plot
range(data_teste$prec)
# generate break positions
breaks = c(339,650,1000,1350,1700,2123)
# and labels
labels = as.character(breaks)
# plot and be happy
my_plot3 = ggplot(data_teste, aes(x=prec, color=species)) + geom_density(size=1.3)+
  #scale_color_viridis(discrete=T,option="D")+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(339, 2123), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(339, 2123), breaks = breaks, labels = labels)
my_plot3 =  my_plot3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic",size=15),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "Density",size=10) +
  labs(colour= "Species")+
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/precipitation_species.pdf"),plot=my_plot3,width=20,height=5)

# Density plots ### Climatic Water Deficit
range(data_teste$cwd)
# generate break positions
breaks = c(166,300,450,600,750,966)
# and labels
labels = as.character(breaks)
my_plot4 = ggplot(data_teste, aes(x=cwd, color=species)) + geom_density(size=1.3)+
  # scale_color_viridis(discrete=T,option="D")+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(166, 966), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(166, 966), breaks = breaks, labels = labels)
my_plot4 =  my_plot4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic",size=15),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=10) +
  labs(colours= "Species")+
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/climaticwd_species.pdf"),plot=my_plot4,width=20,height=5)

### Comparing bioclimatic niche of each species inside current SDA and the 
### two most important variables ######

# Import data_set
a.suare <- read.csv(file=paste0("Adansonia.suarezensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

range(a.suare$tseas) # 1198 - 1375 
range(a.suare$tseasf) # 1307- 1504.33

mean(a.suare$tseas) # 1270.8
mean(a.suare$tseasf) # 1384.73

breaks = c(1150,1200,1250,1300,1350,1400,1450,1500)
labels = as.character(breaks)
par(mfrow = c(1,2))
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
ggsave(file=paste0("./outputs/adan.suare_current_future_niche_comparison_in_sda.pdf"),
       plot=plot.seas,width=20,height=5)

### second IV
range(a.suare$prec) # 1112 - 1523 
range(a.suare$precf) # 922.367 1330.372

mean(a.suare$prec) 
mean(a.suare$precf) 

breaks = c(800,900,1000,1100,1200,1300,1400,1500,1600)
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
  scale_x_continuous(limits = c(830, 1600), breaks = breaks, labels = labels) +
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
# Import data_set
a.perrieri <- read.csv(file=paste0("Adansonia.perrieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.perrieri)
tail(a.perrieri)

# and label - s
range(a.perrieri$tseas) # 916 - 1964
range(a.perrieri$tseasf) # 1064,3 - 2193,3

mean(a.perrieri$tseas) # 1413,4
mean(a.perrieri$tseasf) # 1571,22
# range(mapmat.final$tseas)

breaks = c(800,1000,1200,1400,1600,1800,2000,2200)
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
  scale_x_continuous(limits = c(800, 2200), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas =  plot.seas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic",size=15),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100)", y = "Density",size=10) +
  #labs(x="Annual precipitation (mm.y-1)",y="Mean annual temp. (°C x 10)",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.perrieri_current_future_niche_comparison_in_sdap.pdf"),
       plot=plot.seas,width=20,height=5)

## 2nd Most IV a perrieri - annual mean temp

# and label - s
range(a.perrieri$tmean) # 180 - 266
range(a.perrieri$tmeanf) # 212,66 - 297,66

mean(a.perrieri$tmean) # 241,16
mean(a.perrieri$tmeanf) # 269,35
# range(mapmat.final$tseas)

breaks = c(160,180,200,220,240,260,280,300,320)
labels = as.character(breaks)
plot.cwd = ggplot(a.perrieri, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(100, 1850), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(150, 320), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.cwd =  plot.cwd+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            legend.text = element_text(face= "italic",size=15),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10)", y = "Density",size=15) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.perrieri_current_future_niche_comparison_in_sdap_v2.pdf"),plot=plot.cwd,width=10,height=5)


############## A.rubrostipa cwd comparison graphic

a.rubrostipa <- read.csv(file=paste0("Adansonia.rubrostipa/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.rubrostipa)
tail(a.rubrostipa)


# range cwd
range(a.rubrostipa$cwd) # 706 - 963
range(a.rubrostipa$cwdf) # 1387 - 1962.3

mean(a.rubrostipa$cwd) # 815,64
mean(a.rubrostipa$cwdf) # 1671.56

breaks = c(600,700,820,960,1370,1650,1975)
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
  #geom_vline(xintercept = c(675, 2050), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(675, 2050), breaks = breaks, labels = labels) +
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

### second IV tmean

range(a.rubrostipa$tmean) # 234 269
range(a.rubrostipa$tmeanf) # 264.3 306.3

mean(a.rubrostipa$tmean) # 257.48
mean(a.rubrostipa$tmeanf) # 291.21
breaks = c(220,230,240,250,260,270,280,290,300,320)
labels = as.character(breaks)
plot.precv2= ggplot(a.rubrostipa, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  #geom_vline(xintercept = c(100, 1900), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(220, 320), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.precv2 =  plot.precv2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_text(face= "italic",size=15),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (ºC x 10)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.rubrostipa_current_future_niche_comparison_in_sda_v2.pdf"),
       plot=plot.precv2,width=20,height=5)

############## A. madagascariensis seas comparison graphic

a.madagascariensis <- read.csv(file=paste0("Adansonia.madagascariensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.madagascariensis)
tail(a.madagascariensis)

# range cwd
range(a.madagascariensis$tseas) # 891 - 2014
range(a.madagascariensis$tseasf) # 1061.6 - 2066

mean(a.madagascariensis$tseas) # 1464.298
mean(a.madagascariensis$tseasf) # 1510.303

breaks = c(850,1100,1350,1600,1850,2100)
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
 # geom_vline(xintercept = c(775, 2100), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(850,2200), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas2 =  plot.seas2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic",size=15),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (Cº sd x 100)", y = "Density",size=10) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(3))) +
  theme(axis.title.y = element_text(size = rel(3))) +
  theme(axis.text.x = element_text(size = rel(3))) +
  theme(axis.text.y = element_text(size = rel(3)))
ggsave(file=paste0("./outputs/a.mada_current_future_niche_comparison_in_sdapv1.pdf"),plot=plot.seas2,width=10,height=5)

### second MIV mean annual temp a mada

range(a.madagascariensis$tmean) # 250 - 274
range(a.madagascariensis$tmeanf) # 27.3 - 311.6

mean(a.madagascariensis$tmean) # 263.45
mean(a.madagascariensis$tmeanf) # 297.06

breaks = c(235,250,265,280,290,305,320)
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
  scale_x_continuous(limits = c(235,330), breaks = breaks, labels = labels) +
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
####################

############## A. grandidieri seas comparison graphic

a.grandidieri <- read.csv(file=paste0("Adansonia.grandidieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.grandidieri)
tail(a.grandidieri)

# range cwd
range(a.grandidieri$prec) # 437 - 998
range(a.grandidieri$precf) # 367 - 924

mean(a.grandidieri$prec) # 749
mean(a.grandidieri$precf) # 712

breaks = c(350,450,550,650,750,850,950,1050)
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
  #geom_vline(xintercept = c(325,1075), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(325,1100), breaks = breaks, labels = labels) +
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

## second importance variable mean annual 

range(a.grandidieri$tmean) # 242 - 261
range(a.grandidieri$tmeanf) # 275 297

mean(a.grandidieri$tmean) # 250,7
mean(a.grandidieri$tmeanf) # 286,1

breaks = c(240,250,260,270,280,290,300)
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
  #geom_vline(xintercept = c(230, 305), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(230, 305), breaks = breaks, labels = labels) +
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

############## A za - precipitation

a.za <- read.csv(file=paste0("Adansonia.za/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.za)
tail(a.za)


# range cwd
range(a.za$prec) # 341 - 1827
range(a.za$precf) # 329.5 - 1536

mean(a.za$prec) # 751.586
mean(a.za$precf) # 703

breaks = c(300,600,900,1200,1500,1800,2000)
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
  #geom_vline(xintercept = c(225,2075), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(225,2075), breaks = breaks, labels = labels) +
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


## Second iv A za tmean


range(a.za$tmean) # 215 275
range(a.za$tmeanf) # 251.6 311

mean(a.za$tmean) # 243.3
mean(a.za$tmeanf) # 279.3

breaks = c(200,220,240,260,280,300,320)
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
  geom_vline(xintercept = c(205, 320), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(205, 320), breaks = breaks, labels = labels) +
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

############## A digitata - precipitation

a.digitata <- read.csv(file=paste0("Adansonia.digitata/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.digitata)
tail(a.digitata)
str(a.suare)

# range tseas 2nd most important variable
range(a.digitata$tseas) # 864 1861
range(a.digitata$tseasf) # 1021- 2045
a.digitata  <- subset(a.digitata, !is.na(tseasf) & !is.na(tseas))

mean(a.digitata$tseas) # 1523
mean(a.digitata$tseasf) # 1566,8

breaks = c(800,1100,1400,1700,2000,2200)
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
  #geom_vline(xintercept = c(800, 2100), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(800,2100), breaks = breaks, labels = labels) +
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

#### 1st miv A dig cwd

range(a.digitata$cwd) # 662 - 920
range(a.digitata$cwdf) # 1123 - 1940,6

mean(a.digitata$cwd) # 822,4
mean(a.digitata$cwdf) # 1731,32

breaks = c(600,700,820,960,1370,1650,1975)
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
ggsave(file=paste0("./outputs/a.digitata_current_future_niche_comparison_in_sdapv1.pdf"),
       plot=plot.cwd,width=20,height=5)

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
