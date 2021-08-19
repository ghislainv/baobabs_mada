#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

library(grid)
require(ggplot2)

### Stack variables through Madagascar

map_result <- stack(s,environ$alt)
as.v.na <- na.omit(map_result)

as.v.na2 <- rasterToPoints(as.v.na)
dados_tcc <- as.data.frame(as.v.na2[complete.cases(as.v.na2), ] )

set.seed(20)
data_frame2 <- sample_n(dados_tcc, size = 1000, replace = F)

# Plot 
range(data_frame2$alt)
range(data_frame2$tmean)

### Alt x Tmean
map.mat <- ggplot(data_frame2, aes(x=alt, y=tmean)) + xlim(0,2041) + ylim(137,280) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.5) + 
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Elevation (m)",y="Mean annual temp. (Â°C x 10)",size=4) +
  annotate("text", x  = 0, y = 280 , size=7, label = "(a)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
#ggsave(file=paste0("./outputs/meantem_elevation.png"),
      # plot=map.mat,width=10,height=10,dpi='print')

#### Altitude x T. seas
range(data_frame2$alt)
range(data_frame2$tseas)

map.mat2 <- ggplot(data_frame2, aes(x=alt, y=tseas)) + xlim(0,2041) + ylim(991,3500) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.4) +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Elevation (m)",y="Temp. seasonality (Â°C sd x 1000)",size=4) +
  annotate("text", x  = 0, y = 3500 , size=7, label = "(b)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
#ggsave(file=paste0("./outputs/meantempseas_elevation.png"),
      # plot=map.mat2,width=10,height=10,dpi='print')

## Ploting comparing the latitude
## Tmean x Latitude

range(data_frame2$tmean)
range(data_frame2$y) # y latitude, don't forget ;) 7172500 8673500

lat.mat1 <- ggplot(data_frame2, aes(x=data_frame2$y, y=tmean)) + ylim(137,280) + xlim(7172500,8673500) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.4) +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Latitude (UTM)",y="Mean annual temp. (Â°C x 10)",size=4) +
  annotate("text", x  = 7200000, y = 280 , size=7, label = "(c)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
#ggsave(file=paste0("./outputs/meantemp_latitude.png"),
      # plot=lat.mat1,width=10,height=10,dpi='print')

## Tseas x Latitude

range(data_frame2$tseas)
range(data_frame2$y) #y latitude, don't forget ;) 7172500 8673500

lat.mat2 <- ggplot(data_frame2, aes(x=data_frame2$y, y=tseas)) + ylim(991,3500) + xlim(7172500,8673500) +
  #geom_density2d(data=Abs.df,col=grey(0.5)) +
  geom_point(data=data_frame2,col="darkgreen",alpha=0.4) +
  geom_smooth(method=loess , color="red", fill="#69b3a2", se=TRUE) +
  labs(x="Latitude (UTM)",y="Temp. seasonality (Â°C sd x 1000)",size=4) +
  annotate("text", x  = 7200000, y = 3500 , size=7, label = "(d)") +
  theme(plot.margin=unit(c(0.5,1,1,0.5), "lines"), 
        axis.title.x=element_text(size=rel(2)),
        axis.title.y=element_text(size=rel(2)),
        axis.text.x=element_text(size=rel(1.75)),
        axis.text.y=element_text(size=rel(1.75)))
#ggsave(file=paste0("./outputs/meantemp.seas_latitude.png"),
       #plot=lat.mat2,width=10,height=10,dpi='print')

a5 <- grid.arrange( map.mat, map.mat2, lat.mat1, lat.mat2, nrow=2, ncol=2)


ggsave(file=paste0("./outputs/four_chart_lat_elevations.png"),
       plot=a5,width=12,height=12,dpi=300)
ggsave(filename=here("outputs/four_chart_lat_elevation.pdf"), plot=a5,
       width=12, height=12, dpi=300)

# ===========
# End of file
# ===========
