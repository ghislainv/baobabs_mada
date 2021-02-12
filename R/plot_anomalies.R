#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

library(raster)
library(viridis)
library(ggplot2)
library(grid)
library(gridExtra)
library(here)

# Present climate
current <- stack(here("data/gisdata/sdm_variables/current.tif"))
names(current) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),
                    paste("prec",1:12,sep=""),paste("bio",1:19,sep=""),
                    paste("pet",1:12,sep=""),"pet","cwd","ndm")

# Future climate
he_85_2080 <- stack(here("data/gisdata/sdm_variables/he_85_2080.tif"))
gs_85_2080 <- stack(here("data/gisdata/sdm_variables/gs_85_2080.tif"))
no_85_2080 <- stack(here("data/gisdata/sdm_variables/no_85_2080.tif"))
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
    legend.text=element_text(size=8),
    legend.key.height=unit(0.5,"line"),
    legend.key.width=unit(1.10,"line"),
    legend.box.background=element_blank(),
    ## Plot
    plot.title=element_text(hjust=0.5,size=10),
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
      annotate("text",x=600000,y=8500000,label=label,hjust=1,vjust=0,size=4) +
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
                       rot=90, gp=gpar(cex=1), hjust=0.5, vjust=0.5)
tgrob_fut <- textGrob("Future anomaly\n(RCP 8.5, 2085)",
                      rot=90, gp=gpar(cex=1), hjust=0.5, vjust=0.5)
plot_anomalies <- grid.arrange(tgrob_pres, tgrob_fut,
                               p1, p2, p3, p4, p5, p6, p7, p8,
                               layout_matrix=lay)

ggsave(filename=here("outputs/climatic_anomalies.png"), plot=plot_anomalies,
       width=textwidth, height=0.8*textwidth, dpi=300, units="cm")

# ggsave(filename=here("outputs/Fig1_climatic_anomalies.pdf"), plot=plot_anomalies,
#       width=8, height=7, dpi="print")

# ===========
# End of file
# ===========