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

# Accuracy_indices (see Liu 2011 and Pontius 2008)
accuracy_indices <- function(pred, obs, digits=2) {
  
  if (identical(dim(pred),as.integer(c(2,2)))) {
    df <- pred
    n00 <- df[1,1]
    n10 <- df[2,1]
    n01 <- df[1,2]
    n11 <- df[2,2]
  } else {
    
    # Create data-frame
    df <- data.frame(pred=pred, obs=obs)
    
    # Confusion matrix
    n00 <- sum((df["pred"] == 0) & (df["obs"] == 0))
    n10 <- sum((df["pred"] == 1) & (df["obs"] == 0))
    n01 <- sum((df["pred"] == 0) & (df["obs"] == 1))
    n11 <- sum((df["pred"] == 1) & (df["obs"] == 1))
  }
  
  # Accuracy indices
  N <- n11 + n10 + n00 + n01
  OA <- (n11 + n00) / N
  FOM <- n11 / (n11 + n10 + n01)
  Sensitivity <- n11 / (n11 + n01)
  Specificity <- n00 / (n00 + n10)
  TSS <- Sensitivity + Specificity - 1
  Prob_1and1 <- ((n11 + n10) / N) * ((n11 + n01) / N)
  Prob_0and0 <- ((n00 + n01) / N) * ((n00 + n10) / N)
  Expected_accuracy <- Prob_1and1 + Prob_0and0
  Kappa <- (OA - Expected_accuracy) / (1 - Expected_accuracy)
  
  # Results         
  r <- data.frame(OA=round(OA, digits), 
                  EA=round(Expected_accuracy, digits),
                  FOM=round(FOM, digits),
                  Sen=round(Sensitivity, digits),
                  Spe=round(Specificity, digits),
                  TSS=round(TSS, digits),
                  K=round(Kappa, digits))
  
  return(r)
}