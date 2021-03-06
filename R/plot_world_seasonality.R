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
                          title="Temperature seasonality anomaly (current vs. 2085 RCP 8.5,°C sd x 100)") + col_scale_var_anomws

# future seasonality
fut_ws <- plot_anomaly_ws(r=future_ws, label="(c)",
                          title="Future temperature seasonality (2085 RCP 8.5,°C sd x 100)") + col_scale_var_fut_ws

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

# ===========
# End of file
# ===========
