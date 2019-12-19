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

##=========================================
## OK until here: corrections GV 19/12/2019
## Mario, please arrange code below
##=========================================
  
##==================================
## Calculate anomalies
##==================================

## Temperature seasonality
## Actual

setwd("~/MEGA/Artigos/baobas_article/anomalies")

Xvar <- stack("./Xvar.tif") # rasters available at Vieilledent et al. 2016 --> https://doi.org/10.1111/1365-2745.12548

bio4.region <- Xvar[[5]]

## Mean in 2080
setwd("~/MEGA/Artigos/baobas_article/")

### Tseas
bio4.gs.2080 <- raster("./anomalies/climproj/bio4_gs_85_2080_region.tif")
bio4.he.2080 <- raster("./anomalies/climproj/bio4_he_85_2080_region.tif")
bio4.no.2080 <- raster("./anomalies/climproj/bio4_no_85_2080_region.tif")
Stack.bio4.2080 <- stack(c(bio4.gs.2080,bio4.he.2080,
                           bio4.no.2080))
bio4.2080 <- mean(Stack.bio4.2080)

##= Annual mean temp
## Actual

bio1.region <- Xvar[[4]]  # Annual mean present
bio1.gs.2080 <- raster("./anomalies/climproj/bio1_gs_85_2080_region.tif")
bio1.he.2080 <- raster("./anomalies/climproj/bio1_he_85_2080_region.tif")
bio1.no.2080 <- raster("./anomalies/climproj/bio1_no_85_2080_region.tif")
Stack.bio1.2080 <- stack(c(bio1.gs.2080,bio1.he.2080,
                           bio1.no.2080))

bio1.2080 <- mean(Stack.bio1.2080) ### 2080 temp seas. raster

## Annual mean precipitation

bio12.region <- Xvar[[6]]  # Prec present

bio6.gs.2080 <- raster("./anomalies/climproj/bio12_gs_85_2080_region.tif")
bio6.he.2080 <- raster("./anomalies/climproj/bio12_he_85_2080_region.tif")
bio6.no.2080 <- raster("./anomalies/climproj/bio12_no_85_2080_region.tif")
Stack.bio12.2080 <- stack(c(bio6.gs.2080,bio6.he.2080,
                            bio6.no.2080))
bio12.2080 <- mean(Stack.bio12.2080) ### 2080 prec. raster

## Climatic water deficit

s_cwd_fut

# configuring
step1 <- resample(bio1.2080,s_cwd_fut,"bilinear")
ex_2 <- extent(s_cwd_fut)
step2 <- crop(step1,ex_2)
step3 <- mask(step1,step2)

##
stepa <- resample(bio4.2080,s_cwd_fut,"bilinear")
ex_2 <- extent(s_cwd_fut)
stepb <- crop(stepa,ex_2)
stepc <- mask(stepa,stepb)

##

stepd <- resample(bio12.2080,s_cwd_fut,"bilinear")
ex_2 <- extent(s_cwd_fut)
stepe <- crop(stepd,ex_2)
stepf <- mask(stepd,stepe)
anomalies_sf <- stack (step3,stepc,stepf,s_cwd_fut)

names(anomalies_sf) <- c("tempf","tseasf","precf","cwdf")

## Draw points in the SDA and extract future environmental variables
# In SDA

wC.anomalies <- which(values(ca)>500) 
nC.anomalies <- length(wC.anomalies)
Samp.anomalies <- if (nC.anomalies>1000) {sample(wC.anomalies,1000,replace=FALSE)} else {wC.anomalies}
mapmat.df.anomalies <- as.data.frame(anomalies_sf)[Samp.anomalies,] 

## table to compare current and future bioclimatic changes for each species 
# Used this table to create density niche curves to compare future baobabs niche
## over current distribution (ca)

names(mapmat.df.anomalies) <- c("tmeanf","tseasf","precf","cwdf")
mapmat.final <- cbind(mapmat.df.anomalies,mapmat.df)
head(mapmat.final)
mapmat.f <- mapmat.final[,c(1,5,2,6,3,7,4,8,9,10)]
head(mapmat.f)
setwd("~/MEGA/Artigos/baobas_article/BIOMOD")
write.csv2(mapmat.f,paste0(spdir,"/niche_graph_species_compared_anomaly.csv"))

#### Future niche over current SDA!!! NO ANOMALY

wC.future <- which(values(ca)>500)
mapmat.df.future <- as.data.frame(anomalies_sf)[wC.future,]
names(mapmat.df.future) <- c("tmeanf","tseasf","precf","cwdf")
Mean_ok_future <- round(apply(mapmat.df.future,2,mean,na.rm=TRUE))
q_ok_future <- round(apply(mapmat.df.future,2,quantile,c(0.025,0.975),na.rm=TRUE))
niche_ok_future <- as.data.frame(rbind(Mean_ok_future,q_ok_future))
setwd("~/MEGA/Artigos/baobas_article/BIOMOD")
write.table(niche_ok_future,paste0(spdir,"/mean_niche_with_future.txt"),sep="\t")
write.csv2(niche_ok_future,paste0(spdir,"/mean_niche_with_future.csv"))

# Maps for article

### 1 - Figure with current, delta and 2080 climatic for 4 environmental variables

setwd("~/MEGA/Artigos/baobas_article/anomalies")

##= Temperature seasonality - 
## Current 

bio4.region <- Xvar[[5]]  # Temp Seas present

## Mean in 2080  # Download 2080 temp seas rasters

bio4.gs.2080 <- raster("./climproj/bio4_gs_85_2080_region.tif")
bio4.he.2080 <- raster("./climproj/bio4_he_85_2080_region.tif")
bio4.no.2080 <- raster("./climproj/bio4_no_85_2080_region.tif")
Stack.bio4.2080 <- stack(c(bio4.gs.2080,bio4.he.2080,
                           bio4.no.2080))
bio4.2080 <- mean(Stack.bio4.2080) ### 2080 temp seas. raster

## Plot present

bio4.anomalies.pres <- bio4.region
bio4.anomalies.pres[] <- bio4.2080[]-bio4.region[]

range(bio4.anomalies.pres[],na.rm=TRUE)


pdf("./results/seas_three_chart_map.pdf",width=10,height=7)
pdf("./results/anomalies_map.pdf",width=40,height=28) # all maps

# pdf(paste0(spdir,"/results/seas_three_hart_map.pdf"),width=10,height=7)
par(mfrow=c(2,4))
par(cex=1.2,mar=c(0,0,4,0))## 2010
plot(bio4.region,col=viridis_pal(option ="D")(255),
     axes=FALSE,box=FALSE,legend=T, horizontal=T, 
     #breaks=breakpoints2,axis.args=a.arg2, zlim=c(800,3400),
     main="Temp.seas(sd x 100)")
## Anomaly future
par(cex=1.2,mar=c(0,0,4,0))
plot(bio4.anomalies.pres,col=viridis_pal(option ="D")(255),
     scale_fill_viridis_d(direction = -1),
     axes=FALSE,box=FALSE,legend=T,
     main="Future anomaly",horizontal=TRUE)
dev.off()

######################################

###########################################
### 2 - Figure with current, delta and 2080 climatic for 4 environmental variables

##= Annual mean temp - 
## Actual

bio1.region <- Xvar[[4]]  # Annual mean present
bio1.gs.2080 <- raster("./climproj/bio1_gs_85_2080_region.tif")
bio1.he.2080 <- raster("./climproj/bio1_he_85_2080_region.tif")
bio1.no.2080 <- raster("./climproj/bio1_no_85_2080_region.tif")
Stack.bio1.2080 <- stack(c(bio1.gs.2080,bio1.he.2080,
                           bio1.no.2080))
bio1.2080 <- mean(Stack.bio1.2080) ### 2080 temp seas. raster

## Plot present

bio1.anomalies.pres <- bio1.region
bio1.anomalies.pres[] <- bio1.2080[]-bio1.region[]
range(bio1.anomalies.pres[],na.rm=TRUE)


pdf("./results/tmean_three_chart_map_ok.pdf",width=10,height=7)
par(mfrow=c(1,2))
## 2010
par(cex=1.2,mar=c(0,0,4,0))
plot(bio1.region,col=viridis_pal(option ="B")(255),
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,legend=T, horizontal=T,
     main="Annual Mean Temp(?Cx10)")

## Anomaly future

par(cex=1.2,mar=c(0,0,4,0))
plot(bio1.anomalies.pres,col=viridis_pal(option ="B")(255),
     axes=FALSE,box=FALSE,legend=T,
     main="Future anomaly",horizontal=TRUE)
dev.off()

### 3 - Figure with current, delta and 2080 climatic for 3 environmental variables

## Actual

bio12.region <- Xvar[[6]]  # Prec present

## Mean in 2080  # Download 2080  prec rasters

bio6.gs.2080 <- raster("./climproj/bio12_gs_85_2080_region.tif")
bio6.he.2080 <- raster("./climproj/bio12_he_85_2080_region.tif")
bio6.no.2080 <- raster("./climproj/bio12_no_85_2080_region.tif")
Stack.bio12.2080 <- stack(c(bio6.gs.2080,bio6.he.2080,
                            bio6.no.2080))
bio12.2080 <- mean(Stack.bio12.2080) ### 2080 prec. raster

## Plot present

bio12.anomalies.pres <- bio12.region
bio12.anomalies.pres[] <- bio12.2080[]-bio12.region[]
range(bio12.anomalies.pres[],na.rm=TRUE)

##= Plots of prec charts
pdf("./results/prec_three_chart_map3.pdf",width=10,height=7)
par(mfrow=c(1,2))
## 2010
par(cex=1.2,mar=c(0,0,4,0))
plot(bio12.region,col=viridis_pal(option ="E",direction = -1)(255),
     axes=FALSE,box=FALSE,legend=T, horizontal=T,
     main="Annual Precipitation (mm.y-1)")

## new test

library(scales)

nCol <- 112
myCol <- viridis(n = nCol)
myCol
myCol <- viridis(n = nCol, option = "E")
show_col(myCol)
show_col(myCol)

breakpoints <- c(-300,-285,-270,-255,-240,-225,-210,-195,
                 -180,-165,-150,-135,-120,-105,
                 -90,-75,-55,-40,-25,-10,0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
colors <- c("#FFEA46FF","#FFE843FF","#FFE543FF","#FDE346FF","#FAE149FF","#F7DF4BFF",
            "#F4DC4EFF","#F2DA50FF","#F0D852FF","#ECD555FF","#EAD357FF",
            "#E7D259FF","#E4CF5BFF","#E2CD5CFF","#E0CB5EFF","#DDC95FFF",
            "#D5C364FF","#D0BE67FF","#CCBA69FF","#C6B66BFF",#"#C1B26DFF",
            #"#BCAF6FFF" ,"#B7AA71FF" ,"#B2A672FF" ,"#B0A473FF","#ACA174FF","#A69D75FF" 
            "white",
            "#888579FF","#848279FF","#807E79FF",
            "#797977FF","#757575FF","#6D6E72FF","#656870FF","#5F636EFF","#52596CFF","#4C546CFF",
            "#444F6BFF","#3C4A6BFF", "#34456BFF", "#2A406CFF", "#1D3B6DFF","#07366EFF",
            "#00336FFF", "#002F6FFF","#002860FF","#002251FF")

a.arg <- list(at=seq(-300,100,length.out=5), labels=c("-300","-200","-100","0","100"))

par(cex=1.2,mar=c(0,0,4,0))
plot(bio12.anomalies.pres,col=colors, #col=viridis_pal(option ="E",direction = -1)(255),
     breaks=breakpoints, axis.arg=a.arg, zlim=c(-300,100),
     axes=FALSE,box=FALSE,legend=T,#breaks=breakpoints2,axis.args=a.arg2, zlim=c(-300,100),
     main="Future Climatic Anomaly",horizontal=TRUE)

dev.off()


#### 4 - Climatic Water Deficit

## Actual
setwd("~/MEGA/Artigos/baobas_article/anomalies")

names(Stack.cwd.2080) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),
                           paste("prec",1:12,sep=""),paste("bio",1:19,sep=""),
                           paste("pet",1:12,sep=""),"pet","cwd","ndm")

## Rasters of environmental variables
if (!file.exists("data/rasters_anomalia/environ.tif")) {
  list.url <- c(
    ## Environment
    #"https://madaclim.cirad.fr/environ/environ.tif",
    ## Current climate
    #"https://madaclim.cirad.fr/climate/current.tif",
    #ACESS1-0
    #"https://madaclim.cirad.fr/climate/ac_85_2080.tif",
    ## CCSM4
    #"https://madaclim.cirad.fr/climate/cc_85_2080.tif",
    ## GISS-E2-R
    "https://madaclim.cirad.fr/climate/gs_85_2080.tif",
    ## HadGEM2-ES
    "https://madaclim.cirad.fr/climate/he_85_2080.tif",
    ## IPSL-CM54-LR
    #"https://madaclim.cirad.fr/climate/ip_85_2080.tif",
    ## MIROC5
    #"https://madaclim.cirad.fr/climate/mc_85_2080.tif",
    ## NorESM1-M (no) Picked up this one because of the cc variability
    "https://madaclim.cirad.fr/climate/no_85_2080.tif"
    
  )
  dir.create("data/rasters_anomalia", recursive=TRUE, showWarnings=FALSE)
  for (i in 1:length(list.url)) {
    dest.f <- file.path("data/rasters_anomalia",basename(list.url[i]))
    curl::curl_download(url=list.url[i],destfile=dest.f)
  }
}

getwd()
current_ano <- stack("data/rasters_anomalia/current.tif")
names(current_ano) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),
                        paste("prec",1:12,sep=""),paste("bio",1:19,sep=""),
                        paste("pet",1:12,sep=""),"pet","cwd","ndm")
wc_ano <- which(names(current_ano) %in% c("cwd"))
s_cwd <- stack(current_ano[[wc_ano]])
names(s_cwd) <- c("cwd")

## Future distribution with Future Data - MadaClim   
model <- c("gs","he","no") # For global climate models (GCMs)
rcps <- c("85") # For representative concentration pathways (RCP): RCP 8.5
yrs <- c("2080") # For ywae 2080
n.mods <- length(model)*length(rcps)*length(yrs)

if (run.models) {
  for (cm in 1:length(model)) {
    for (k in 1:length(rcps)) {
      for (m in 1:length(yrs)) {
        
        ## Message
        i.mods <- (cm-1)*length(rcps)*length(yrs) + (k-1)*length(yrs) + m
        cat(paste0("\n","Model ",i.mods,"/",n.mods,": ",model[cm],"_",rcps[k],"_",yrs[m],"\n"))
        
        ## Load climatic data
        getwd()
        setwd("C:/Users/Usu?rio/Documents/MEGA/Artigos/baobas_article/data/rasters_anomalia")
        futures <- stack(paste0(model[cm],"_",rcps[k],"_",yrs[m],".tif"))
        
        names(futures) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),
                            paste("prec",1:12,sep=""),paste("bio",1:19,sep=""),
                            paste("pet",1:12,sep=""),"pet","cwd","ndm")
      }
    }
  }
}

future_cwd <- stack("data/rasters_anomalia/current.tif")

cwd_fut_2080 <- which(names(futures) %in% c("cwd"))
s_cwd_fut <- stack(futures[[cwd_fut_2080]])

## Plot present

cwd.anomalies.pres <- s_cwd
cwd.anomalies.pres[] <- s_cwd_fut[]-s_cwd[]

range(cwd.anomalies.pres[],na.rm=TRUE)
range(s_cwd[],na.rm=TRUE)
##= Plots of prec charts
setwd("~/MEGA/Artigos/baobas_article/anomalies")

pdf("./results/cwd_three_chart_map_new.pdf",width=10,height=7)

par(mfrow=c(1,2))
## 2010
plot(s_cwd)
par(cex=1.2,mar=c(0,0,4,0))
plot(s_cwd,col=viridis_pal(option ="C")(255),
     axes=FALSE,box=FALSE,legend=T, horizontal=T,
     main="Annual Climatic Water Deficit (mm)")

## Anomaly future

par(cex=1.2,mar=c(0,0,4,0))
plot(cwd.anomalies.pres,col=viridis_pal(option ="C")(255),
     axes=FALSE,box=FALSE,legend=T,
     main="Future anomaly",horizontal=TRUE)
dev.off()

###################################################################################
######## Generating variables histograms graphs##################
###############################################################################

setwd("~/MEGA/Artigos/baobas_article/BIOMOD")

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

setwd("~/MEGA/Artigos/baobas_article/anomalies")

# Density plots ## Temp. Seasonality
range(data_teste$tseas)
# generate break positions
breaks = c(854,1200,1600,2000,2400,2800,3257)
# and labels
labels = as.character(breaks)

# plot

my_plot = ggplot(data_teste, aes(x=tseas, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(854, 3257), show.legend = T, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(854, 3257), breaks = breaks, labels = labels)

my_plot =  my_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           legend.text = element_text(face= "italic"),
                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100?C)", y = "Density", size=10) +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("./results/temp_seas_interval_3.pdf"),plot=my_plot,width=10,height=5)

# Density plots ### Annual Mean Temperature
range(data_teste$tmean)
# generate break positions
breaks = c(181,200,220,240,260,275)
# and labels
labels = as.character(breaks)
my_plot2 = ggplot(data_teste, aes(x=tmean, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(181, 275), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(181, 275), breaks = breaks, labels = labels)
my_plot2 =  my_plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic"),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (?C x 10) ", y = "Density",size=10) +
  labs(colour= "Species") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("./results/annual_mean_temp_interval.pdf"),plot=my_plot2,width=10,height=5)

# Density plots ### Precipitation plot
range(data_teste$prec)
# generate break positions
breaks = c(341,650,1000,1350,1700,2120)
# and labels
labels = as.character(breaks)
# plot and be happy
my_plot3 = ggplot(data_teste, aes(x=prec, color=species)) + geom_density(size=1.3)+
  #scale_color_viridis(discrete=T,option="D")+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(341, 2120), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(341, 2120), breaks = breaks, labels = labels)
my_plot3 =  my_plot3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic"),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Precipitation (mm-y ?)", y = "Density",size=10) +
  labs(colour= "Species")+
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("./results/precipitation_interval.pdf"),plot=my_plot3,width=10,height=5)

# Density plots ### Climatic Water Deficit
range(data_teste$cwd)
# generate break positions
breaks = c(181,300,450,600,750,961)
# and labels
labels = as.character(breaks)
my_plot4 = ggplot(data_teste, aes(x=cwd, color=species)) + geom_density(size=1.3)+
  # scale_color_viridis(discrete=T,option="D")+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = c(181, 961), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = c(181, 961), breaks = breaks, labels = labels)
my_plot4 =  my_plot4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic"),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=10) +
  labs(colours= "Species")+
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("./results/climaticwd_interval.pdf"),plot=my_plot4,width=10,height=5)

### Creating table for seasonality anomalie accordind to most important variable ######

# Import data_set
setwd("~/MEGA/Artigos/baobas_article/BIOMOD")
a.suare <- read.csv(file=paste0("Adansonia.suarezensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

head(a.suare)
tail(a.suare)
list1 <- 1:1000
list2 <- rep("A_suarezensis", length(list1))

a.suare <- cbind(list2,a.suare)
#a.suare$anom_tseas <- round(a.suare$anom_tseas,2)
# and labels
range(a.suare$tseas) # 1194 - 1375 
range(a.suare$tseasf) # 1262 - 1432

mean(a.suare$tseas) # 1271.045
mean(a.suare$tseasf) # 1327.174 

breaks = c(1150,1200,1250,1300,1350,1400,1450,1500)
labels = as.character(breaks)
par(mfrow = c(1,2))
plot.seas = ggplot(a.suare, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=Species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(1175, 1450), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(1175, 1450), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas =  plot.seas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic"),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("adan.suare_current_future_niche_comparison_in_sda.pdf"),plot=plot.seas,width=10,height=5)

### second IV
range(a.suare$prec) # 1112 - 1523 
range(a.suare$precf) # 922.367 1330.372

mean(a.suare$prec) # 1291.187
mean(a.suare$precf) # 1105.046 
# range(mapmat.final$tseas)

breaks = c(800,900,1000,1100,1200,1300,1400,1500,1600)
labels = as.character(breaks)
plot.precv2= ggplot(a.suare, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(precf, fill=Species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(830, 1600), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(830, 1600), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.precv2 =  plot.precv2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_text(face= "italic"),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean annual precipitation (mm.y-?)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))

ggsave(file=paste0("adan.suare_current_future_niche_comparison_in_sda_v2.pdf"),
       plot=plot.precv2,width=10,height=5)

################################
# Import data_set
setwd("~/MEGA/Artigos/baobas_article/BIOMOD")
a.perrieri <- read.csv(file=paste0("Adansonia.perrieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.perrieri)
list1 <- 1:1000
list2 <- rep("A_perrieri", length(list1))
a.perrieri <- cbind(list2,a.perrieri)

colnames(a.perrieri) <- c("Species","X","tmeanf","tmean","tseasf","tseas","precf", "prec", "cwdf", "cwd", "alt")

# and label - s
range(a.perrieri$tseas) # 854 - 1767
range(a.perrieri$tseasf) # 989 - 1920

mean(a.perrieri$tseas) # 1289.79
mean(a.perrieri$tseasf) # 1373.785
# range(mapmat.final$tseas)

breaks = c(800,1000,1200,1400,1600,1800,2000)
labels = as.character(breaks)
plot.seas = ggplot(a.perrieri, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=Species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(800, 1950), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(800, 1950), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas =  plot.seas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic"),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (sd x 100)", y = "Density",size=2) +
  #labs(x="Annual precipitation (mm.y-1)",y="Mean annual temp. (?C x 10)",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.perrieri_current_future_niche_comparison_in_sdap.pdf"),plot=plot.seas,width=10,height=5)

## 2nd Most IV a perrieri - climatic water deficit

# and label - s
range(a.perrieri$cwd) # 198 - 865
range(a.perrieri$cwdf) # 340 - 1746

mean(a.perrieri$cwd) # 541.92
mean(a.perrieri$cwdf) # 994.749
# range(mapmat.final$tseas)

breaks = c(0,250,500,750,1000,1250,1500,1750)
labels = as.character(breaks)
plot.cwd = ggplot(a.perrieri, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(cwdf, fill=Species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(100, 1850), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(100, 1850), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.cwd =  plot.cwd+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            legend.text = element_text(face= "italic"),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic water deficit (mm)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.perrieri_current_future_niche_comparison_in_sdap_cwd_vi2.pdf"),plot=plot.cwd,width=10,height=5)


############## A.rubrostipa cwd comparison graphic

a.rubrostipa <- read.csv(file=paste0("Adansonia.rubrostipa/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.rubrostipa)
tail(a.rubrostipa)
list1 <- 1:1000
list2 <- rep("A_rubrostipa", length(list1))
a.rubrostipa <- cbind(list2,a.rubrostipa)
colnames(a.rubrostipa) <- c("Species","X","tmeanf","tmean","tseasf","tseas","precf", "prec", "cwdf", "cwd", "alt")
# range cwd
range(a.rubrostipa$cwd) # 713 - 960
range(a.rubrostipa$cwdf) # 1372 - 1975

mean(a.rubrostipa$cwd) # 818.676
mean(a.rubrostipa$cwdf) # 1650.095

breaks = c(600,700,820,960,1370,1650,1975)
labels = as.character(breaks)

plot.cwd = ggplot(a.rubrostipa, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(cwdf, fill=Species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(675, 2050), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(675, 2050), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.cwd =  plot.cwd + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic"),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.rubrostipa_current_future_niche_comparison_in_sdapv1.pdf"),plot=plot.cwd,width=10,height=5)

### second IV

range(a.rubrostipa$prec) # 349 1687
range(a.rubrostipa$precf) # 324 1519

mean(a.rubrostipa$prec) # 1030
mean(a.rubrostipa$precf) # 924
breaks = c(300,500,700,900,1100,1300,1500,1700)
labels = as.character(breaks)
plot.precv2= ggplot(a.rubrostipa, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(precf, fill=Species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(100, 1900), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(100, 1900), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.precv2 =  plot.precv2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_text(face= "italic"),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Mean annual precipitation (mm.y-?)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))

ggsave(file=paste0("adan.rubro_current_future_niche_comparison_in_sda_v2.pdf"),
       plot=plot.precv2,width=10,height=5)

############## A. madagascariensis seas comparison graphic

a.madagascariensis <- read.csv(file=paste0("Adansonia.madagascariensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.madagascariensis)
tail(a.madagascariensis)
list1 <- 1:1000
list2 <- rep("A_madagascariensis", length(list1))
a.madagascariensis <- cbind(list2,a.madagascariensis)
colnames(a.madagascariensis) <- c("Species","X","tmeanf","tmean","tseasf","tseas","precf", "prec", "cwdf", "cwd", "alt")

# range cwd
range(a.madagascariensis$tseas) # 861 - 2013
range(a.madagascariensis$tseasf) # 972 - 1957

mean(a.madagascariensis$tseas) # 1464.298
mean(a.madagascariensis$tseasf) # 1510.303

breaks = c(850,1050,1250,1450,1650,1850,2050)
labels = as.character(breaks)

plot.seas2 = ggplot(a.madagascariensis, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=Species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(775, 2100), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(775,2100), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas2 =  plot.seas2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic"),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (?C sd x 100)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.mada_current_future_niche_comparison_in_sdapv1.pdf"),plot=plot.seas2,width=10,height=5)

### second MIV mean annual temp a mada

# range cwd
range(a.madagascariensis$tmean) # 249 - 275
range(a.madagascariensis$tmeanf) # 278 - 312

mean(a.madagascariensis$tmean) # 263.941
mean(a.madagascariensis$tmeanf) # 298.5791

breaks = c(235,250,265,280,290,305,320)
labels = as.character(breaks)

plot.mean2 = ggplot(a.madagascariensis, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=Species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(235, 320), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(235,320), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.mean2 =  plot.mean2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic"),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (?C x 10) ", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.mada_current_future_niche_comparison_in_sdapv2.pdf"),plot=plot.mean2,width=10,height=5)
####################

############## A. grandidieri seas comparison graphic

a.grandidieri <- read.csv(file=paste0("Adansonia.grandidieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.grandidieri)
tail(a.grandidieri)
list1 <- 1:1000
list2 <- rep("A_grandidieri", length(list1))
a.grandidieri <- cbind(list2,a.grandidieri)
colnames(a.grandidieri) <- c("Species","X","tmeanf","tmean",
                             "tseasf","tseas","precf", "prec", "cwdf", "cwd", "alt")
# range cwd
range(a.grandidieri$prec) # 410 - 1006
range(a.grandidieri$precf) # 375 - 930

mean(a.grandidieri$prec) # 745,5
mean(a.grandidieri$precf) # 705,4

breaks = c(350,450,550,650,750,850,950,1050)
labels = as.character(breaks)

plot.prec = ggplot(a.grandidieri, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=Species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(325,1075), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(325,1075), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.prec =  plot.prec + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic"),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Precipitation (mm.y-?)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.grandi_current_future_niche_comparison_in_sdapv1.pdf"),plot=plot.prec,width=10,height=5)

## second importance variable mean annual 

# range cwd
range(a.grandidieri$tmean) # 242 - 261
range(a.grandidieri$tmeanf) # 276 297

mean(a.grandidieri$tmean) # 250,7
mean(a.grandidieri$tmeanf) # 286,6

breaks = c(240,250,260,270,280,290,300)
labels = as.character(breaks)

plot.mean2 = ggplot(a.grandidieri, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=Species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(230, 305), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(230, 305), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.mean2 =  plot.mean2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic"),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (?C x 10) ", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.grand_current_future_niche_comparison_in_sdapv2.pdf"),plot=plot.mean2,width=10,height=5)
############## A za - precipitation

a.za<- read.csv(file=paste0("Adansonia.za/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.za)
tail(a.za)
list1 <- 1:1000
list2 <- rep("A_za", length(list1))
a.za <- cbind(list2,a.za)
colnames(a.za) <- c("Species","X","tmeanf","tmean",
                    "tseasf","tseas","precf", "prec", "cwdf", "cwd", "alt","species")
library(viridis)
# range cwd
range(a.za$prec) # 340 - 1555
range(a.za$precf) # 321.5 - 1313.3

mean(a.za$prec) # 807.351
mean(a.za$precf) # 722

breaks = c(300,500,700,900,1100,1300,1500,1700)
labels = as.character(breaks)


plot.prec = ggplot(a.za, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=Species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(225,2075), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(225,2075), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.prec =  plot.prec + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_text(face= "italic"),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Precipitation (mm.y-?)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.za_current_future_niche_comparison_in_sdapv1.pdf"),plot=plot.prec,width=10,height=5)


## Second iv A za tmean

# range cwd
range(a.za$tmean) # 215 271
range(a.za$tmeanf) # 249.6 308.2

mean(a.za$tmean) # 244
mean(a.za$tmeanf) # 279

breaks = c(200,220,240,260,280,300,320)
labels = as.character(breaks)

plot.mean2 = ggplot(a.za, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tmeanf, fill=Species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(205, 320), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(205, 320), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.mean2 =  plot.mean2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic"),
                                 panel.background = element_blank(), 
                                 axis.line = element_line(colour = "black"))+
  labs(x="Mean Annual Temperature (?C x 10) ", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.za_current_future_niche_comparison_in_sdapv2.pdf"),plot=plot.mean2,width=10,height=5)

############## A digitata - precipitation

a.digitata<- read.csv(file=paste0("Adansonia.digitata/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.digitata)
tail(a.digitata)
list1 <- 1:1000
list2 <- rep("A_digitata", length(list1))
a.digitata <- cbind(list2,a.digitata)
colnames(a.digitata) <- c("Species","X","tmeanf","tmean",
                          "tseasf","tseas","precf", "prec", "cwdf", "cwd", "alt")
# range tseas
range(a.digitata$tseas) # 864 1861
range(a.digitata$tseasf) # 999- 1952

mean(a.digitata$tseas) # 1523
mean(a.digitata$tseasf) # 1509

breaks = c(800,1000,1200,1400,1600,1800,2000)
labels = as.character(breaks)

plot.seas2 = ggplot(a.digitata, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=Species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(800, 2100), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(800,2100), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.seas2 =  plot.seas2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_text(face= "italic"),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (?C sd x 100)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.digi_current_future_niche_comparison_in_sdapv1.pdf"),plot=plot.seas2,width=10,height=5)

#### 2nd miv A dig cwd

# range cwd
range(a.digitata$cwd) # 668 - 921
range(a.digitata$cwdf) # NA NA there are one NA in the df
a.digitata$cwdf[is.na(a.digitata$cwdf)] <- sample(a.digitata$cwdf)
range(a.digitata$cwdf) # 1110 1850

mean(a.digitata$cwd) # 826
mean(a.digitata$cwdf) # 1657

breaks = c(600,700,820,960,1370,1650,1975)
labels = as.character(breaks)

plot.cwd = ggplot(a.digitata, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=Species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(cwdf, fill=Species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  geom_vline(xintercept = c(650, 1950), colour="gray70", linetype="dashed") +
  scale_x_continuous(limits = c(650, 1950), breaks = breaks, labels = labels) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(discrete=T,option="E")
plot.cwd =  plot.cwd + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_text(face= "italic"),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=2) +
  labs(col = "") +
  theme(axis.title.x = element_text(size = rel(1.8))) +
  theme(axis.title.y = element_text(size = rel(1.8))) +
  theme(axis.text.x = element_text(size = rel(1.8))) +
  theme(axis.text.y = element_text(size = rel(1.8)))
ggsave(file=paste0("a.digitata_current_future_niche_comparison_in_sdapv2.pdf"),plot=plot.cwd,width=10,height=5)


#plot.seas= ggplot(mapmat.final, aes(x=tseas, y=..density..)) + 
# geom_density(aes(fill=Species,color= 'tseas'), alpha=.6) + 
#geom_vline(aes(xintercept=mean(tseas)),
#          color="red", linetype="solid", size=0.5) +
#geom_density(aes(tseasf, fill=Species,color='tseasf'),alpha=.6) +
#geom_vline(aes(xintercept=mean(tseasf)),
#          color="green", linetype="solid", size=0.5) +
#geom_vline(xintercept = c(1143, 1567), show.legend = T, colour="gray70", linetype="dashed") +
#scale_x_continuous(limits = c(1143, 1567), breaks = breaks, labels = labels) +
#scale_color_manual(values = c('tseas' = 'red', 'tseasf' = 'green'))+
#scale_fill_viridis(discrete=T,option="A")
#plot.seas =  plot.seas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                            legend.text = element_text(face= "italic"),
#                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#labs(title="Baobab climatic niche: Future - Tseasf - (2080/RCP 8.5)/ Present Niche (Tseas)",
#    x="Temp. Seasonality (sd x 100)", y = "Density",size=6) 
#ggsave(file=paste0(spdir,"./current_future_niche_comparison_in_sda.pdf"),plot=plot.seas,width=10,height=10)

## Developing

# generate the first plot for the species 
# Density plots ## Temp. Seasonality present and future
# generate break positions
# breaks = c(1143,1200,1300,1400,1500,1567)
# and labels
#labels = as.character(breaks)
#range(mapmat.final$tseas)
#range(mapmat.final$tseasf)

# Density plots with semi-transparent fill.
# alpha is the transparency of the overlaid color
#ggplot2.density(data=weight, xName='weight', groupName='sex',
#               legendPosition="top",
#              alpha=0.5, fillGroupDensity=TRUE )
# plot
#  plot.seas= ggplot(mapmat.final, aes(x=tseas, y=..density..)) + 
#   geom_density(aes(fill=Species,color=tseas), alpha=.6) + 
#  geom_vline(aes(xintercept=mean(tseas)),
#            color="red", linetype="dashed", size=1) +
#  geom_density(aes(tseasf, fill=Species,color=tseasf),alpha=.6) +
#scale_fill_viridis(discrete=T,option="E") +
# geom_vline(aes(xintercept=mean(tseasf)),
#           color="green", linetype="dashed", size=1) +
# scale_fill_viridis(discrete=T,option="C")#+
#geom_vline(xintercept = c(1143, 1567), show.legend = T, colour="red", linetype="dashed") +
#  scale_x_continuous(limits = c(1143, 1567), breaks = breaks, labels = labels)
#plot.seas =  plot.seas + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                              legend.text = element_text(face= "italic"),
#                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#labs(title="Baobab climatic niche: Future (2080/RCP 8.5) and present Niche",x="Temp. Seasonality (sd x 100)", y = "Density",size=6) +
#labs(fill= "Species")
#ggsave(file=paste0(spdir/"current_future_niche_comparison_in_sda.pdf"),plot=my_plot,width=10,height=10)

setwd("./BIOMOD/")
##### Draw all species (prec x temp) in the niche ggplot2

# Import niche graph species

database1_suare <- read.csv(file=paste0("Adansonia.suarezensis/","niche_graph_species.csv"), header=T,sep=";")
#database2 <- read.csv(file=paste0("Adansonia.za.perrieri/","niche_graph_species.csv"), header=T,sep=";")
database2_digi <- read.csv(file=paste0("Adansonia.digitata/","niche_graph_species.csv"), header=T,sep=";")
database3_grand <- read.csv(file=paste0("Adansonia.grandidieri/","niche_graph_species.csv"), header=T,sep=";")
database4_mada <- read.csv(file=paste0("Adansonia.madagascariensis/","niche_graph_species.csv"), header=T,sep=";")
database5_perri <- read.csv(file=paste0("Adansonia.perrieri/","niche_graph_species.csv"), header=T,sep=";")
database6_rubro <- read.csv(file=paste0("Adansonia.rubrostipa/","niche_graph_species.csv"), header=T,sep=";")
database7_za <- read.csv(file=paste0("Adansonia.za/","niche_graph_species.csv"), header=T,sep=";")

data_teste <- rbind(database2_digi,database3_grand,
                    database4_mada,database5_perri,database6_rubro,database1_suare,database7_za)
# map.mat <- ggplot(data_teste, aes(x=prec, y=tmean, fill=species)) + 
#  xlim(200,2300) + ylim(200,280) +
# labs(x="Annual precipitation (mm.y-1)",y="Mean annual temp. (?C x 10)",size=2) +
#theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
#      theme(axis.title.x = element_text(size = rel(1.8))) +
#       theme(axis.title.y = element_text(size = rel(1.8))) +
#      theme(axis.text.x = element_text(size = rel(1.8))) +
#     theme(axis.text.y = element_text(size = rel(1.8)))+
#theme(title = element_text(size = rel(1.8))) +
#  theme(legend.text = element_text(size = rel(1.2))) +

# theme(legend.text = element_text(face= "italic")) +
# geom_density2d(aes(col = species)) +
#geom_density2d(aes(colour = species)) + 
#labs(col = "Baobabs species") +
#labs(title = "Baobabs climatic niche comparison") +
#theme(legend.background = element_blank()) +
#theme(legend.position = c(0.87, .25))
#theme_bw()  
#map.mat
#ggsave(file=paste0("tmean_X_prec_all_species.pdf"),plot=map.mat,width=10,height=10)

##################### plot all current sps distributions for article map


#pred_digitata <- stack(paste0("Adansonia.digitata/proj_current/proj_current_Adansonia.digitata_ensemble.grd"))
#pred_grandidieri <- stack(paste0("Adansonia.grandidieri/proj_current/proj_current_Adansonia.grandidieri_ensemble.grd"))
#pred_madagascariensis <- stack(paste0("Adansonia.madagascariensis/proj_current/proj_current_Adansonia.madagascariensis_ensemble.grd"))
#pred_perrieri <- stack(paste0("Adansonia.perrieri/proj_current/proj_current_Adansonia.perrieri_ensemble.grd"))
#pred_rubrostipa <- stack(paste0("Adansonia.rubrostipa/proj_current/proj_current_Adansonia.rubrostipa_ensemble.grd"))
#pred_suar <- stack(paste0("Adansonia.suarezensis/proj_current/proj_current_Adansonia.suarezensis_ensemble.grd"))
#pred_za <- stack(paste0("Adansonia.za/proj_current/proj_current_Adansonia.za_ensemble.grd"))
#stack_all <- stack(pred_digitata,pred_grandidieri,pred_madagascariensis,pred_perrieri,
#                  pred_rubrostipa,pred_suar,pred_za)
#names(stack_all) <- c('Adansonia digitata', 'Adansonia grandidieri', 'Adansonia madagascariensis',
##                    "Adansonia perrieri","Adansonia rubrostipa", "Adansonia suarezensis",
#                  "Adansonia za")
#plot(stack_all)

#breakpoints <- c(-100,100,300,500,700,900,1100)
#colors <- c(grey(c(0.90,seq(0.7,0.3,-0.2))),"#568203","#013220")
#a.arg <- list(at=seq(0,1000,length.out=5), labels=c("0","1","2","3","4"),cex.axis=1.5)
#l.arg <- list(text="Vote",side=2, line=0.5, cex=1.5)
# Plot plot(bio4.anomalies.pres,col=viridis_pal(option ="D")(255)
#getwd()
#pdf(paste0("./grandi.pdf"),width=6.5,height=10)
#par(mar=c(0,0,0,r.mar),cex=1.4)
#plot(stack_all,col=colors,ext=sf, breaks=breakpoints,
#    legend.width=1.5,legend.shrink=0.6,legend.mar=7,
#   axis.args=a.arg,legend.arg=l.arg,axes=FALSE, box=F, zlim=c(0,1000))
#dev.off()


#######################################################################
### World Seasonality Map
#######################################################################
setwd("/MEGA/Artigos/baobabs_article")
library(raster)
library(sp)
library(viridis)
data("wrld_simpl")

r <- raster::getData("worldclim",var="bio",res=10)
r <- r[[c(4)]]
names(r) <- c("Seas")
plot(r)
plot(r, col= viridis_pal(option = "D")(255), xlim=c(-90,160), ylim=c(-23.5,23.5))

bio4.gs.2080_world <- raster("./anomalies/global_future/gs85bi704.tif")
bio4.he.2080_world <- raster("./anomalies/global_future/he85bi704.tif")
bio4.no.2080_world <- raster("./anomalies/global_future/no85bi704.tif")
Stack.bio4.2080_world <- stack(c(bio4.gs.2080_world,bio4.he.2080_world,
                                 bio4.no.2080_world))
bio4.2080_world <- mean(Stack.bio4.2080_world) ### 2080 temp seas. raster
crs(bio4.2080_world) <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(wrld_simpl) <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

## Plot present

bio4.anomalies.pres_world <- r
bio4.anomalies.pres_world[] <- bio4.2080_world[]-r[]


range(bio4.anomalies.pres_world[],na.rm=TRUE)
plot(bio4.anomalies.pres_world, col= viridis_pal(option = "D")(255), xlim=c(-90,160), ylim=c(-23.5,23.5))
pdf("./anomalies/global_future/anomalies_map_world.pdf",width=50,height=28) # all mapas

0.6*500
300/2.1
500+142.85
# range anomalies tropical world: -642.85 - 1102.333
1102.33+642.85
1745.18/40
1102/20

library(scales) #find scales to get hex codes
nCol <- 112
myCol <- viridis(n = nCol)
myCol
myCol <- viridis(n = nCol, option = "E")
x11()
show_col(myCol)

x11()
par(mfrow=c(2,1))
plot(r,col=viridis_pal(option ="E")(255),xlim=c(-90,160),
     ylim=c(-23.5,23.5),
     ##maxpixels=50000,
     axes=FALSE,box=FALSE,legend=T, horizontal=T, 
     #breaks=breakpoints2,axis.args=a.arg2, zlim=c(800,3400),
     main="Temperature seasonality (sd x 100)")
plot(wrld_simpl, add=T)


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
#colors <- c(col=viridis_pal(option ="E",direction = -1)(18))
a.arg <- list(at=seq(-800,1200,length.out=11), 
              labels=c("-800","-600","-400","-200","0","200",
                       "400","600","800","1000","1200"))
plot(bio4.anomalies.pres_world,col=colors, #col=viridis_pal(option ="E",direction = -1)(255),
     breaks=breakpoints, axis.arg=a.arg, zlim=c(-640,1200),xlim=c(-90,160),
     ylim=c(-23.5,23.5),
     axes=FALSE,box=FALSE,legend=T,#breaks=breakpoints2,axis.args=a.arg2, zlim=c(-300,100),
     main="Future Climatic Anomaly",horizontal=TRUE)
plot(wrld_simpl, add=T)

dev.off()

###################################

par(cex=1.2,mar=c(0,0,4,0))
plot(bio4.anomalies.pres_world,col=viridis_pal(option ="D")(255),xlim=c(-90,160),
     ylim=c(-23.5,23.5),
     axes=FALSE,box=FALSE,legend=T,
     main="Future anomaly",horizontal=TRUE)
plot(wrld_simpl, add=T)
dev.off()

# Command to restore R loaded environment 
library(session)
memory.limit(size = 10000000000000)
setwd("~/MEGA/Artigos/baobas_article")
# save.session("baobabs_article_2019.Rda")
restore.session("baobabs_article_2019.Rda")


##===========================================================================
## End of script
##===========================================================================