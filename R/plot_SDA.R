#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

## Libraries
require(ggplot2)

## Colors for legend
gcolors <- colorRampPalette(c("#568203","#013220"))

## Setting basic theme options for plot with ggplot2
theme_base <- theme(
  ## Axis
  axis.line=element_blank(),
  axis.text=element_blank(),
  axis.ticks=element_blank(),
  axis.title=element_blank(),
  ## Legend
  legend.position="bottom",
  legend.title=element_blank(), #add Vote for bottom species
  legend.text=element_text(size=10),
  legend.key.height=unit(0.5,"line"),
  legend.spacing.x=unit(0.0, "line"),
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

## Color altitude
col_alt <- scale_fill_gradientn(
  colours = terrain.colors(255)[255:1],
  na.value="transparent",
  values=rescale(seq(0,3000),0,3000),
  limits=c(0,3000),
  breaks= seq(0,3000,by=1500),
  labels= seq(0,3000,by=1500)
)

## Color scales pres
col_present_SDA <- scale_fill_manual(
  na.value="transparent",
  values=c(grey(seq(0.9,0.7,-0.2)), gcolors(3)),
  labels=as.character(0:4),
  drop=FALSE
)

## Color scales future
col_future_SDA <- scale_fill_manual(
  na.value="transparent",
  values=c(grey(c(0.90,seq(0.7,0.50,-0.05))), gcolors(7)),
  labels=c("0","","2","","4","","6","","8","","10","","12"),
  drop=FALSE
)

## Function to plot altitude and occurrence points
plot_alt <- function(m, n, title, label="x", xlim=c(313000,1090000), ylim=c(7167000,8676000)) {
  # Data
  rdf <- data.frame(rasterToPoints(m)) # To plot raster with geom_raster()
  names(rdf) <- c("x", "y", "z")
  rdf2 <- as.data.frame(n)
  # Text position
  x_text <- xlim[1]+0.10*(xlim[2]-xlim[1])
  y_text <- ylim[1]+0.90*(ylim[2]-ylim[1])
  # Plot
  p <-  ggplot(NULL, aes(x,y)) + 
    geom_tile(data=rdf, aes(fill=z)) +
    col_alt +
    geom_point(data=rdf2, shape=1, aes(x,y), color="black", size=2)+
    theme_bw() + theme_base +
    coord_fixed(xlim=xlim, ylim=ylim) +
    annotate("text",x=x_text,y=y_text,label=label,hjust=0,vjust=1,size=5) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    labs(title=title)
  return(p)
}

## Function to plot current species distribution
plot_current_SDA <- function(r, title, label="x", xlim=c(313000,1090000), ylim=c(7167000,8676000)) {
  # Data
  rdf <- data.frame(rasterToPoints(r)) # To plot raster with geom_raster()
  names(rdf) <- c("x", "y", "z")
  rdf$z_factor <- factor(rdf$z, levels=seq(0,1000,250), labels=as.character(0:4))
  # Text position
  x_text <- xlim[1]+0.10*(xlim[2]-xlim[1])
  y_text <- ylim[1]+0.90*(ylim[2]-ylim[1])
  # Plot
  p <- ggplot(NULL, aes(x, y)) +
    geom_tile(data=rdf, aes(fill=z_factor)) +
    col_present_SDA +
    theme_bw() + theme_base +
    coord_fixed(xlim=xlim, ylim=ylim) +
    annotate("text",x=x_text,y=y_text,label=label,hjust=0,vjust=1,size=5) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    labs(title=title) +
    guides(fill=guide_legend(label.position="bottom", nrow=1, keywidth=unit(1.5,"line")))
  return(p)
}

## Function to plot future species distribution (RCP 4.5 and RCP 8.5)
plot_future_SDA <- function(r, title, label="x", xlim=c(313000,1090000), ylim=c(7167000,8676000)) {
  # Data
  rdf <- data.frame(rasterToPoints(r)) # To plot raster with geom_raster()
  names(rdf) <- c("x", "y", "z")
  rdf$z_factor <- factor(rdf$z, levels=seq(0,3000,250), labels=as.character(0:12))
  # Text position
  x_text <- xlim[1]+0.10*(xlim[2]-xlim[1])
  y_text <- ylim[1]+0.90*(ylim[2]-ylim[1])
  # Plot
  p <- ggplot(NULL, aes(x, y)) +
    geom_tile(data=rdf, aes(fill=z_factor)) +
    col_future_SDA +
    theme_bw() + theme_base + 
    coord_fixed(xlim=xlim, ylim=ylim) +
    annotate("text",x=x_text,y=y_text,label=label,hjust=0,vjust=1,size=5) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    labs(title=title) +
    theme(legend.text=element_text(size=7)) +
    guides(fill=guide_legend(label.position="bottom", nrow=1, keywidth=unit(0.5,"line")))
  return(p)
}

# =============================================
# Zoom regions
# =============================================

# Zoom region for A. perrieri
xmin_per <- 776286; xmax_per <- 1066742
ymin_per <- 8344529; ymax_per <- 8676100
rect_perrieri <- geom_rect(aes(xmin=xmin_per, xmax=xmax_per, ymin=ymin_per, ymax=ymax_per),
                           inherit.aes=FALSE, fill="transparent", color="black")

# Zoom region for A. suarezensis
xmin_sua <- 900000; xmax_sua <- 1040000
ymin_sua <- 8510000; ymax_sua <- 8676100
rect_suarezensis <- geom_rect(aes(xmin=xmin_sua, xmax=xmax_sua, ymin=ymin_sua, ymax=ymax_sua),
                              inherit.aes=FALSE, fill="transparent", color="black")

# =============================================
# Elevation
# =============================================

## Non-threatened species
a_za <- read.table(paste0("Adansonia.za/function_ready.txt"), header=T,sep="\t")
a_digi <- read.table(paste0("Adansonia.digitata/function_ready.txt"), header=T,sep="\t")
a_grand <- read.table(paste0("Adansonia.grandidieri/function_ready.txt"), header=T,sep="\t")

## Threatened species
a_mada <- read.table(paste0("Adansonia.madagascariensis/function_ready.txt"), header=T,sep="\t")
a_perri <- read.table(paste0("Adansonia.perrieri/function_ready.txt"), header=T,sep="\t")
a_rubro <- read.table(paste0("Adansonia.rubrostipa/function_ready.txt"), header=T,sep="\t")
a_suare <- read.table(paste0("Adansonia.suarezensis/function_ready.txt"), header=T,sep="\t")

## Non-threatened species
a1 <- plot_alt(m=environ$alt,n=a_digi,label="(a)", title="Presence Points\nElevation(m)")
e2 <- plot_alt(m=environ$alt,n=a_grand,label="(e)", title="")
i3 <- plot_alt(m=environ$alt,n=a_za,label="(i)", title="")

## Threatened species
a1t <- plot_alt(m=environ$alt,n=a_mada,label="(a)", title="Presence Points\nElevation (m)")
e2t <- plot_alt(m=environ$alt,n=a_perri,label="(e)", title="") + rect_perrieri
i3t <- plot_alt(m=environ$alt,n=a_rubro,label="(i)", title="") 
m4t <- plot_alt(m=environ$alt,n=a_suare,label="(m)", title="") + rect_suarezensis

# =============================================
# Current SDA
# =============================================

## Non-threatened species
ca_dig <- raster(here("Adansonia.digitata/committee_averaging.tif"))
ca_gran <- raster(here("Adansonia.grandidieri/committee_averaging.tif"))
ca_za <- raster(here("Adansonia.za/committee_averaging.tif"))

## Threatened species
ca_mada <- raster(here("Adansonia.madagascariensis/committee_averaging.tif"))
ca_perri <- raster(here("Adansonia.perrieri/committee_averaging.tif"))
ca_rubro <- raster(here("Adansonia.rubrostipa/committee_averaging.tif"))
ca_suare <- raster(here("Adansonia.suarezensis/committee_averaging.tif"))

## Non-threatened species
b1 <- plot_current_SDA(r=ca_dig,label="(b)", title="Current\nDistribution")
f2 <- plot_current_SDA(r=ca_gran,label="(f)", title="")
j3 <- plot_current_SDA(r=ca_za,label="(j)", title="")

## Threatened species
b1t <- plot_current_SDA(r=ca_mada,label="(b)", title="Current\nDistribution")
f2t <- plot_current_SDA(r=ca_perri,label="(f)", title="", xlim=c(xmin_per, xmax_per), ylim=c(ymin_per, ymax_per))
j3t <- plot_current_SDA(r=ca_rubro,label="(j)", title="")
n4t <- plot_current_SDA(r=ca_suare,label="(n)", title="", xlim=c(xmin_sua, xmax_sua), ylim=c(ymin_sua, ymax_sua))

# =============================================
# Future SDA
# =============================================

## Non-threatened species
cafut_dig <- raster(paste0("Adansonia.digitata/caFut_85_2080.tif"))
cafut_gran <- raster(paste0("Adansonia.grandidieri/caFut_85_2080.tif"))
cafut_za <- raster(paste0("Adansonia.za/caFut_85_2080.tif"))
# Zero-dispersal hypothesis
caZD_dig <- raster(paste0("Adansonia.digitata/caZD_85_2080.tif"))
caZD_gran <- raster(paste0("Adansonia.grandidieri/caZD_85_2080.tif"))
caZD_za <- raster(paste0("Adansonia.za/caZD_85_2080.tif"))

## Threatened species
cafut_mada <- raster(paste0("Adansonia.madagascariensis/caFut_85_2080.tif"))
cafut_perri <- raster(paste0("Adansonia.perrieri/caFut_85_2080.tif"))
cafut_rubro <- raster(paste0("Adansonia.rubrostipa/caFut_85_2080.tif"))
cafut_suare <- raster(paste0("Adansonia.suarezensis/caFut_85_2080.tif"))
# Zero-dispersal hypothesis
caZD_mada <- raster(paste0("Adansonia.madagascariensis/caZD_85_2080.tif"))
caZD_perri <- raster(paste0("Adansonia.perrieri/caZD_85_2080.tif"))
caZD_rubro <- raster(paste0("Adansonia.rubrostipa/caZD_85_2080.tif"))
caZD_suare <- raster(paste0("Adansonia.suarezensis/caZD_85_2080.tif"))

## Non-threatened species
c1 <- plot_future_SDA(r=cafut_dig, label="(c)", title="Full Dispersal\nRCP 8.5 2085")
g2 <- plot_future_SDA(r=cafut_gran, label="(g)", title="")
k3 <- plot_future_SDA(r=cafut_za, label="(k)", title="")
# Zero-dispersal hypothesis
d1 <- plot_future_SDA(r=caZD_dig, label="(d)", title="Zero Dispersal\nRCP 8.5 2085")
h2 <- plot_future_SDA(r=caZD_gran, label="(h)", title="")
l3 <- plot_future_SDA(r=caZD_za, label="(l)", title="")

## Threatened species
c1t <- plot_future_SDA(r=cafut_mada, label="(c)", title="Full Dispersal\nRCP 8.5 2085")
g2t <- plot_future_SDA(r=cafut_perri, label="(g)", title="", xlim=c(xmin_per, xmax_per), ylim=c(ymin_per, ymax_per))
k3t <- plot_future_SDA(r=cafut_rubro, label="(k)", title="")
o4t <- plot_future_SDA(r=cafut_suare, label="(o)", title="", xlim=c(xmin_sua, xmax_sua), ylim=c(ymin_sua, ymax_sua))
# Zero-dispersal hypothesis
d1t <- plot_future_SDA(r=caZD_mada, label="(d)", title="Zero Dispersal\nRCP 8.5 2085")
h2t <- plot_future_SDA(r=caZD_perri, label="(h)", title="", xlim=c(xmin_per, xmax_per), ylim=c(ymin_per, ymax_per))
l3t <- plot_future_SDA(r=caZD_rubro, label="(l)", title="")
p4t <- plot_future_SDA(r=caZD_suare, label="(p)", title="", xlim=c(xmin_sua, xmax_sua), ylim=c(ymin_sua, ymax_sua))

# =============================================
# Text
# =============================================

## Non-threatened species
tgrob_dig <- textGrob("A. digitata",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                      hjust=0.5, vjust=2)
tgrob_gran <- textGrob("A. grandidieri",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                       hjust=0.5, vjust=2)
tgrob_za <- textGrob("A. za",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                     hjust=0.5, vjust=2)

## Threatened species
tgrob_mada <- textGrob("A. madagascariensis",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                       hjust=0.5, vjust=2)
tgrob_perri <- textGrob("A. perrieri",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                        hjust=0.5, vjust=2)
tgrob_rubro <- textGrob("A. rubrostipa",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                        hjust=0.5, vjust=2)
tgrob_suare <- textGrob("A. suarezensis",rot=90, gp=gpar(cex=1.25,fontface="italic"),
                        hjust=0.5, vjust=2)

# =============================================
# Combine plot and save
# =============================================

## Non-threatened species
lay <- rbind(c(1, 4:7), c(2, 8:11), c(3, 12:15))
## Column widths and row heights
w <- c(1, 4, 4, 4, 4)
h <- rep(1, 3)
## Arrange plots
plot_baobabs_non_threat <- grid.arrange(
  tgrob_dig, tgrob_gran, tgrob_za,
  a1,b1,c1,d1,e2,f2,g2,h2,
  i3,j3,k3,l3,
  layout_matrix=lay, widths=w, heights=h)
## Save
ggsave(file=here("outputs/plot_SDA_non_threat.png"),
       plot=plot_baobabs_non_threat, width=9, height=10, dpi="print")

## Threatened species
lay <- rbind(c(1, 5:8), c(2, 9:12),
             c(3, 13:16), c(4, 17:20))
## Column widths and row heights
w <- c(1, 4, 4, 4, 4)
h <- rep(1, 4)
plot_baobabs_threatened <- grid.arrange(
  tgrob_mada, tgrob_perri, tgrob_rubro,tgrob_suare,
  a1t,b1t,c1t,d1t,e2t,f2t,g2t,h2t,
  i3t,j3t,k3t,l3t,m4t,n4t,o4t,p4t,
  layout_matrix=lay, widths=w, heights=h)
## Save
ggsave(file=here("outputs/plot_SDA_threat.png"),
       plot=plot_baobabs_threatened, width=9, height=12, dpi="print")

# ===========
# End of file
# ===========