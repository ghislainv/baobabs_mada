# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com / mario.tagliari@posgrad.ufsc.br
# license         :GPLv3
# ==============================================================================

require(raster)

# Automatic zoom
fun.extent <- function(sp,s) {
  width <- xmax(sp)-xmin(sp)
  height <- ymax(sp)-ymin(sp)
  ## Zoom if extent less than ~1/3 of Madagascar extent (NS: 1600 km x EW: 600km)
  if (width<200000 & height<500000) {
    zoom <- TRUE
    x1 <- floor(xmin(sp)-50000)
    x2 <- ceiling(xmax(sp)+50000)
    y1 <- floor(ymin(sp)-100000)
    y2 <- ceiling(ymax(sp)+100000)
    e.map <- c(max(x1,xmin(s)),min(x2,xmax(s)),
               max(y1,ymin(s)),min(y2,ymax(s)))
    r.mar <- 4
  }  else {  ## Else return Madagascar extent 
    zoom <- FALSE
    e.map <- extent(s)
    r.mar <- 1
  }
  return(list(zoom=zoom,e.map=e.map,r.mar=r.mar))
}