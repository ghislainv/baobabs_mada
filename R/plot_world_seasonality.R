#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

library(maptools)
library(raster)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)
library(here)

# Using WorlClim v1.4
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

## current
Stack.bio4.current_worldcrop <- crop(r,wp)

current_ws <- Stack.bio4.current_worldcrop
Stack.bio4.current_worldcrop[] <- future_ws[]-current_ws[]

## Save rasters
tseas <- stack(current_ws, future_ws, Stack.bio4.current_worldcrop)
names(tseas) <- c("tseas_pres", "tseas_fut", "tseas_ano")
# writeRaster(tseas, file="outputs/tseas_tropics.tif", overwrite=TRUE, options=c("COMPRESS=DEFLATE", "PREDICTOR=2"))
## Load data
tseas <- stack("outputs/tseas_tropics.tif")

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
cur_ws <- plot_anomaly_ws(r=tseas[[1]], label="(a)",
                          title="Current temperature seasonality (°C sd x 1000)") + col_scale_var_cur_ws

# Future anomaly
ano_ws <- plot_anomaly_ws(r=tseas[[3]], label="(b)",
                          title="Temperature seasonality anomaly (current vs. 2085 RCP 8.5,°C sd x 1000)") + col_scale_var_anomws

# future seasonality
fut_ws <- plot_anomaly_ws(r=tseas[[2]], label="(c)",
                          title="Future temperature seasonality (2085 RCP 8.5,°C sd x 1000)") + col_scale_var_fut_ws

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

ggsave(file=paste0("./outputs/anomaly_world_chart.png"),
       plot=plot_world_anomaly,width=11,height=8,dpi=300,scale=0.8)

# =============================================
# Figure for press release
# =============================================

## Color scales
## current world seasonality
col_scale_cur <- scale_fill_gradientn(
  colours=viridis(255, option="B", direction=1),
  na.value="transparent",
  values=rescale(seq(0,7.600,l=255),0,7.600),
  limits=c(0,7.600),
  breaks=seq(0,7.600,l=6),
  labels=seq(0,7.600,l=6)
)

## Color scales
## future world seasonality
col_scale_fut <- scale_fill_gradientn(
  colours=viridis(255, option="B", direction=1),
  na.value="transparent",
  values=rescale(seq(0,8.200,l=255),0,8.200),
  limits=c(0,8.200),
  breaks=seq(0,8.200,l=6),
  labels=seq(0,8.200,l=6)
)

## Color scales
## anomaly world seasonality
col_scale_ano <- scale_fill_gradientn(
  colours=c(grey(seq(0.3, 0.7, l=3)), viridis(255, option="C", direction=1)),
  na.value="transparent",
  values=rescale(c(-0.7, -0.35, seq(0, 1, l=255)), -0.7, 1),
  limits=c(-0.7,1),
  breaks=c(-0.7, -0.35, seq(0, 1, l=3)),
  labels=c(-0.7, -0.35, seq(0, 1, l=3))
)

# Plots
cur_ws <- plot_anomaly_ws(r=tseas[[1]]/1000, label="(a)",
                          title="Current temperature seasonality (in °C)") + col_scale_cur
# Future anomaly
ano_ws <- plot_anomaly_ws(r=tseas[[3]]/1000, label="(b)",
                          title="Temperature seasonality anomaly (2085 vs. current, in °C)") + col_scale_ano
# future seasonality
fut_ws <- plot_anomaly_ws(r=tseas[[2]]/1000, label="(c)",
                          title="Future temperature seasonality (2085, in °C)") + col_scale_fut

lay_6 <- rbind(c(rep(seq(1,1,by=1),each=3)),
               c(rep(seq(1,1,by=1),each=3)),
               c(rep(seq(1,1,by=1),each=3)),
               c(rep(seq(2,2,by=1),each=3)),
               c(rep(seq(2,2,by=1),each=3)),
               c(rep(seq(2,2,by=1),each=3)),
               c(rep(seq(3,3,by=1),each=3)),
               c(rep(seq(3,3,by=1),each=3)),
               c(rep(seq(3,3,by=1),each=3)))

plot_tropics <- grid.arrange(cur_ws, ano_ws, fut_ws,layout_matrix=lay_6)

ggsave(file=paste0("outputs/tseas_tropics_degree.png"),
       plot=plot_tropics,width=11,height=8,dpi=300,scale=0.8)

# =============================================
# Tseas unit in WorlClim v1.4 or 2.1
# =============================================

# Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>

# Directory for outputs
dir.create("wc_comp")

# Version 1.4
zipfile="wc_comp/bio_10m_bil.zip"
download.file("https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip", destfile=zipfile, method = "wget")
unzip(zipfile, files=c("bio4.bil", "bio4.hdr"), exdir="wc_comp")
tseas_wc_1_4 <- raster("wc_comp/bio4.bil")
pdf("wc_comp/tseas_wc_1_4.pdf")
plot(tseas_wc_1_4, main="tseas WorldClim 1.4") # goes from 0 to 20000 (**twenty** thousand)
dev.off()

# Version 2.1
zipfile="wc_comp/wc2.1_10m_bio.zip"
download.file("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip", destfile=zipfile, method = "wget")
unzip(zipfile, files=c("wc2.1_10m_bio_4.tif"), exdir="wc_comp")
tseas_wc_2_1 <- raster("wc_comp/wc2.1_10m_bio_4.tif")
pdf("wc_comp/tseas_wc_2_1.pdf")
plot(tseas_wc_2_1, main="tseas WorldClim 2.1") # goes from 0 to 2000 (**two** thousand)
dev.off()
# crop to common extent
tseas_wc_2_1_crop <- raster::crop(tseas_wc_2_1, tseas_wc_1_4)

# Sample points and compare
tseas <- stack(tseas_wc_1_4, tseas_wc_2_1_crop)
set.seed(1234)
samp <- as.data.frame(sampleRandom(tseas, 1000))
names(samp) <- c("wc_1_4", "wc_2_1")
# Simple linear model to find correlation
mod <- lm(wc_1_4~ -1 + wc_2_1, data=samp)
print(mod)
# Call:
# lm(formula = wc_1_4 ~ -1 + wc_2_0, data = samp)
# 
# Coefficients:
# wc_2_0  
#  9.706
pdf("wc_comp/comp.pdf")
plot(samp$wc_2_1, samp$wc_1_4,
     xlab="bio4 (tseas) WorldClim 2.1",
     ylab="bio4 (tseas) WorldClim 1.4")
abline(a=0, b=1, col="red")
dev.off()

# ===========
# End of file
# ===========