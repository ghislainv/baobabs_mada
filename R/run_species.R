#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

# Libraries
require(rgdal)
require(sp)
require(raster)
require(biomod2)
require(grDevices) # for colorRampPalette
require(here)

# Source miscellaneous functions
source(here("R/miscellaneous_functions.R"))

# Colors for legend
gcolors <- colorRampPalette(c("#568203","#013220"))

# Run computation for each species
run.species <- function (i, path_to_maxent.jar, run.models=TRUE) {
  
  spdir <- sp.dir[i]
  spname <- sp.names[i]
  
  ## Message
  cat(paste(spname,"\n",sep=""))
  
  ##===============================
  ## Select data for target species
  df.tsp <- df.sp[df.sp$Species==spname,]
  
  ##===================
  ## Remove duplicates
  cell.pres <- cellFromXY(s,df.tsp)
  
  ## Spatial points at center of each raster cell with presences
  list.cell.pres <- sort(unique(cell.pres))
  cell.pres.sp <- xyFromCell(s,list.cell.pres,spatial=TRUE)
  
  ## Build data-set for presence
  d.presence <- as.data.frame(extract(s,cell.pres.sp)) 
  Coords.presence <- coordinates(cell.pres.sp)
  colnames(Coords.presence) <- c("x","y")
  data.xy <- cbind(d.presence,Coords.presence)
  
  ## Data points with complete environmental information
  wcomp <- which(complete.cases(data.xy))
  
  ## Transform as a SpatialPointsDataFrame and SpatialPoints (for presence only)
  proj_38S <- "+proj=utm +zone=38 +south +datum=WGS84 +units=m +no_defs +type=crs"
  d <- SpatialPointsDataFrame(coords=Coords.presence[wcomp,], data=data.xy[wcomp,],
                              proj4string=CRS(proj_38S))
  
  p <- SpatialPoints(d) ## This is used for presence-only data
  crs(p)  <- proj_38S
  ## Save "d" as data-frame
  dir.create(paste0(spdir,"/figures"),recursive=TRUE)
  write.table(d@data,file=paste0(spdir,"/figures/presences.txt"),sep="\t",row.names=FALSE)
  ## Save "p" as data-frame
  points_occurs <- as.data.frame(p)
  write.table(points_occurs,file=paste0(spdir,"/function_ready.txt"),sep="\t",row.names=FALSE)
  
  ## Extent for species distribution area maps
  ext <- fun.extent(p,s)
  zoom <- ext$zoom
  e.map <- ext$e.map
  r.mar <- ext$r.mar
  
  ## Number of 1km pixels with at least one presence
  npix <- nrow(d)
  write_csv(npix, here("data/baobabs/1km_pixels_one_presence.csv"))
  
  ## BIOMOD_FormatingData
  set.seed(1234) ## Reproducible pseudo-absences
  BiomodData <- BIOMOD_FormatingData(resp.var=p,
                                     expl.var=s,
                                     resp.name=spname,
                                     PA.nb.rep=1,
                                     PA.nb.absences=10000,
                                     PA.strategy="random",
                                     na.rm=TRUE)
  saveRDS(BiomodData,paste0(spdir,"/BiomodData.rds"))
  
  ## BIOMOD_ModelingOptions
  BiomodOption <- BIOMOD_ModelingOptions(GLM=list(type="quadratic", interaction.level=0, myFormula=NULL,
                                                  family=binomial(link="logit"), test="AIC"),
                                         GAM=list(algo="GAM_mgcv", type="s_smoother", k=4, interaction.level=0, 
                                                  myFormula=NULL, 
                                                  family=binomial(link="logit")),
                                         RF=list(do.classif=TRUE, ntree=500),
                                         MAXENT.Phillips=list(path_to_maxent.jar='D:/OneDrive/Cap_1_outros_papers/script_art_1/maxent', 
                                                              visible=FALSE, maximumiterations=500,
                                                              memory_allocated=512,
                                                              # To avoid overparametrization (Merow  et al.  2013)
                                                              product=FALSE, threshold=FALSE, hinge=FALSE))
  
  ## BIOMOD_Modeling
  if (run.models) {
    set.seed(1234) ## Reproducible results
    BiomodModel <- BIOMOD_Modeling(BiomodData,
                                   models=c("GLM","GAM","RF","MAXENT.Phillips"),
                                   models.options=BiomodOption,
                                   NbRunEval=5,
                                   DataSplit=70,
                                   VarImport=3,
                                   models.eval.meth=c("KAPPA","TSS","ROC"),
                                   rescal.all.models=F,
                                   do.full.models=T, # useful for rare species, or few PrPts, but it doesnt separate data for 
                                   # models evaluation, which lead to over-optimistic eval. scores
                                   modeling.id="4mod") ## 4 statistical models
  } else {
    BiomodModel <- get(load(paste0(spdir,"/",spdir,".4mod.models.out")))
  }
  
  ## Building ensemble-models
  if (run.models) {
    BiomodEM <- BIOMOD_EnsembleModeling(modeling.output = BiomodModel,
                                        chosen.models = grep("_Full_",
                                                             get_built_models(BiomodModel),
                                                             value=TRUE), #Full models only
                                        em.by = 'all', 
                                        eval.metric=c("TSS"),
                                        # In the ensemble, we keep only models for which TSS >= 0.6
                                        eval.metric.quality.threshold=0.6,
                                        models.eval.meth=c('TSS','ROC'),
                                        prob.mean=F,
                                        prob.ci=F,
                                        prob.ci.alpha=0.05,
                                        committee.averaging=TRUE,
                                        prob.mean.weight=F,
                                        prob.mean.weight.decay="proportional")
  } else {
    BiomodEM <- get(load(paste0(spdir,"/",spdir,".4modensemble.models.out")))
  }
  
  ## BIOMOD_Projection == PRESENT == ## Individual model projection
  if (run.models) {
    BiomodProj <- BIOMOD_Projection(modeling.output=BiomodModel,
                                    new.env=s,
                                    proj.name="current",
                                    selected.models= grep("_Full_",
                                                          get_built_models(BiomodModel),
                                                          value=TRUE),
                                    binary.meth=c("TSS"),
                                    filtered.meth=c("TSS"),
                                    compress=TRUE,
                                    build.clamping.mask=FALSE, ## No need for present
                                    omi.na=TRUE,
                                    on_0_1000=TRUE,
                                    output.format=".grd")
  } else {
    BiomodProj <- get(load(paste0(spdir,"/proj_current/",spdir,".current.projection.out")))
  }
  
  ## BIOMOD_EnsembleForecasting == PRESENT == ## Ensemble forecasting
  if (run.models) {
    BiomodEF <- BIOMOD_EnsembleForecasting(EM.output=BiomodEM, ## Rules for assembling
                                           projection.output=BiomodProj, ## Individual model projection
                                           binary.meth=c("TSS"),
                                           filtered.meth=c("TSS"),
                                           compress=TRUE,
                                           on_0_1000=TRUE)
  } else {
    BiomodEF <- get(load(paste0(spdir,"/proj_current/",spdir,".current.ensemble.projection.out")))
  }
  
  ## Future distribution with Future Data - MadaClim   
  mod <- c("no","he","gs") # For global climate models (GCMs): NorESM1-M, HadGEM2-ES, GISS-E2-R,  
  rcp <- c("45","85") # For representative concentration pathways (RCPs): RCP 4.5, RCP 8.5
  yr <- c("2050","2080") # For 2050, 2080
  n.mod <- length(mod)*length(rcp)*length(yr)
  
  if (run.models) {
    for (mc in 1:length(mod)) {
      for (j in 1:length(rcp)) {
        for (l in 1:length(yr)) {
          
          ## Message
          i.mod <- (mc-1)*length(rcp)*length(yr) + (j-1)*length(yr) + l
          cat(paste0("\n","Model ",i.mod,"/",n.mod,": ",mod[mc],"_",rcp[j],"_",yr[l],"\n"))
          
          ## Load climatic data
          future <- stack(paste0(dir_var_sdm,mod[mc],"_",rcp[j],"_",yr[l],".tif"))
          names(future) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),
                             paste("prec",1:12,sep=""),paste("bio",1:19,sep=""),
                             paste("pet",1:12,sep=""),"pet","cwd","ndm")
          ## Stack of explicative variables
          sf <- stack(future[[wc]]) ## See wc and we indexes before
          names(sf) <- c("tmean","tseas","prec","cwd")
          
          ## Projections by model
          BiomodProjFuture <- BIOMOD_Projection(modeling.output=BiomodModel,
                                                new.env=sf,
                                                proj.name=paste0(mod[mc],"_",rcp[j],"_",yr[l]),
                                                selected.models=grep("_Full_",
                                                                     get_built_models(BiomodModel),
                                                                     value=TRUE), ## Full models only
                                                binary.meth=c("TSS"),
                                                filtered.meth=c("TSS"),
                                                compress=TRUE,
                                                build.clamping.mask=FALSE, #
                                                omi.na=TRUE,
                                                on_0_1000=TRUE,
                                                output.format=".grd")
          
          ## BIOMOD_EnsembleForecasting == FUTURE == ## Ensemble forecasting
          BiomodEF_Future <- BIOMOD_EnsembleForecasting(EM.output=BiomodEM, ## Rules for assembling
                                                        projection.output=BiomodProjFuture, ## Individual model projection
                                                        binary.meth=c("TSS"),
                                                        filtered.meth=c("TSS"),
                                                        compress=TRUE,
                                                        on_0_1000=TRUE)
          
        }
      }
    }
  } # End of if(run.models)
  
  ##=====================
  ## Current distribution
  
  ## Presence points and altitude
  # Legend specifications
  a.arg <- list(at=seq(0,3000,length.out=4),labels=seq(0,3000,length.out=4),cex.axis=1.5)
  l.arg <- list(text="Elevation (m)",side=2, line=0.5, cex=2.5)
  # Plot
  pdf(paste0(spdir,"/figures/presence_alt.pdf"),width=6.5,height=10)
  par(mar=c(0,0,0,1),cex=1.4)
  plot(environ$alt,col=terrain.colors(255)[255:1],
       legend.width=1.5,legend.shrink=0.6,legend.mar=7,
       axis.args=a.arg,legend.arg=l.arg,
       axes=FALSE,box=FALSE,zlim=c(0,3000))
  plot(p,pch=1,add=TRUE,cex=1)
  if (zoom) {rect(e.map[1],e.map[3],e.map[2],e.map[4],border="black",lwd=1.5)}
  dev.off()
  
  ## Committee averaging
  # Legend specifications
  breakpoints <- c(-125,125,375,625,875,1125)
  colors <- c(grey(seq(0.9,0.7,-0.2)),gcolors(3))
  a.arg <- list(at=seq(0,1000,length.out=5), labels=c("0","1","2","3","4"),cex.axis=1.5)
  l.arg <- list(text="Vote",side=2, line=0.5, cex=1.5)
  # Load predictions and update extent
  pred <- stack(paste0(spdir,"/proj_current/proj_current_",spdir,"_ensemble.grd"))
  ca <- pred[[1]]
  
  # Plot
  pdf(paste0(spdir,"/figures/ca_current.pdf"),width=6.5,height=10)
  par(mar=c(0,0,0,r.mar),cex=1.4)
  plot(ca,col=colors,ext=e.map, breaks=breakpoints,
       legend.width=1.5,legend.shrink=0.6,legend.mar=7,
       axis.args=a.arg,legend.arg=l.arg,
       axes=FALSE, box=FALSE, zlim=c(0,1000))
  dev.off()
  
  # Export as tif
  writeRaster(ca,paste0(spdir,"/committee_averaging.tif"),datatype="INT2S",
              overwrite=TRUE,
              options=c("COMPRESS=LZW","PREDICTOR=2"))
  
  ## Present species distribution area (km2)
  ## values(ca)>=500 gives vote >=2
  SDA.pres <- sum(values(ca)>=500,na.rm=TRUE) # Just the sum because one pixel is 1km2.
  
  ##=================
  ## Ecological niche
  
  ## 95% quantiles for alt, temp, prec, tseas, cwd
  wC <- which(values(ca)>=500)
  niche.df <- as.data.frame(s)[wC,]
  niche.df$alt <- environ$alt[wC]
  Mean <- round(apply(niche.df,2,mean,na.rm=TRUE))
  q <- round(apply(niche.df,2,quantile,c(0.025,0.975),na.rm=TRUE))
  niche <- as.data.frame(rbind(Mean,q))
  write.table(niche,paste0(spdir,"/niche.txt"),sep="\t")
  
  ## Draw points in the SDA and extract environmental variables
  # In SDA
  nC <- length(wC)
  Samp <- if (nC>1000) {sample(wC,1000,replace=FALSE)} else {wC}
  mapmat.df <- as.data.frame(s)[Samp,]
  mapmat.df$alt <- environ$alt[Samp]
  mapmat.df$species <- rep(c(spdir))
  write.csv2(mapmat.df,paste0(spdir,"/niche_graph_species.csv"))
  
  ## Model performance of committee averaging
  
  # Individual models
  Perf.mods <- as.data.frame(as.table(get_evaluations(BiomodModel)))
  names(Perf.mods) <- c("wIndex","Index","Model","Run","PA","Value")
  write.table(Perf.mods,paste0(spdir,"/current_model_evaluation.txt"),sep="\t")
  
  # Observations and committee averaging predictions
  ObsData <- BiomodData@data.species
  ObsData[is.na(ObsData)] <- 0
  caData <- values(ca)[cellFromXY(ca,xy=BiomodData@coord)]
  PredData <- rep(0,length(caData))
  PredData[caData>=500] <- 1
  
  # Computing accuracy indices to evaluate model performance
  Perf.ca <- accuracy_indices(pred=PredData,obs=ObsData,digits=3)
  write.table(Perf.ca,paste0(spdir,"/performance_ca.txt"),sep="\t")
  
  ## Variable importance
  VarImp <- as.data.frame(get_variables_importance(BiomodModel)[,,"Full","PA1"])
  Rank <- as.data.frame(apply(-VarImp,2,rank))
  VarImp$mean.rank <- apply(Rank,1,mean)
  VarImp$rank <- rank(VarImp$mean.rank,ties.method="max")
  write.table(VarImp,paste0(spdir,"/varimp.txt"),sep="\t")
  
  ##====================
  ## Future distribution
  
  ## Considering RCP 45,85 and year 2050,2080
  ## Table for change in area
  SDA.fut <- data.frame(area.pres=SDA.pres,
                        rcp=rep(c("45","85"),each=4),yr=rep(rep(c("2050","2080"),each=2),2),
                        rcp=rep(c("45","85"),each=4),yr=rep(rep(c("2050","2080"),each=2),2),
                        disp=rep(c("full","zero"),4),area.fut=NA)
  ## Change in altitude
  Alt.fut <- data.frame(mean.pres=niche$alt[1],q1.pres=niche$alt[2],q2.pres=niche$alt[3],
                        rcp=rep(c("45","85"),each=4),yr=rep(rep(c("2050","2080"),each=2),2),
                        disp=rep(c("full","zero"),4),mean.fut=NA,q1.fut=NA,q2.fut=NA)
  
  ## Committee averaging for the three GCMs (sum)
  for (j in 1:length(rcp)) {
    for (l in 1:length(yr)) {
      # Create stack of model predictions
      caS <- stack() 
      for (mc in 1:length(mod)) {
        name.f <- paste0("/proj_",mod[mc],"_",rcp[j],"_",yr[l])
        pred.f <- stack(paste0(spdir,name.f,name.f,"_",spdir,"_ensemble.grd"))
        caS <- addLayer(caS,pred.f[[1]])
      }
      
      # Compute sum
      caFut <- sum(caS)
      # Remove data for Comoro Island
      values(caFut)[cellsCom] <- NA
      # Legend
      breakpoints <- seq(-125,3125,by=250)
      colors <- c(grey(c(0.90,seq(0.7,0.50,-0.05))),gcolors(7))
      a.arg <- list(at=seq(0,3000,length.out=13), labels=c(0:12), cex.axis=1.2)
      l.arg <- list(text="Vote",side=2, line=0.5, cex=1.2)
      
      # Plot (Committee Averaging Full Dispersal)
      writeRaster(caFut,paste0(spdir,"/caFut_",rcp[j],"_",yr[l],".tif"),datatype="INT2S",
                  overwrite=TRUE,
                  options=c("COMPRESS=LZW","PREDICTOR=2"))
      pdf(paste0(spdir,"/figures/cafd_",rcp[j],"_",yr[l],".pdf"),width=6.5,height=10)
      par(mar=c(0,0,0,r.mar),cex=1.4)
      plot(caFut,col=colors,breaks=breakpoints,ext=e.map,
           legend.width=1.5,legend.shrink=0.75,legend.mar=7,
           axis.args=a.arg,legend.arg=l.arg,
           axes=FALSE,box=FALSE,zlim=c(0,3000))
      dev.off()
      
      # Zero-dispersal
      caZD <- caFut
      values(caZD)[values(ca)<500] <- 0 
      writeRaster(caZD,paste0(spdir,"/caZD_",rcp[j],"_",yr[l],".tif"),datatype="INT2S",
                  overwrite=TRUE,
                  options=c("COMPRESS=LZW","PREDICTOR=2"))
      pdf(paste0(spdir,"/figures/cazd_",rcp[j],"_",yr[l],".pdf"),width=6.5,height=10)
      par(mar=c(0,0,0,r.mar),cex=1.4)
      plot(caZD,col=colors,breaks=breakpoints,ext=e.map,
           legend.width=1.5,legend.shrink=0.75,legend.mar=7,
           axis.args=a.arg,legend.arg=l.arg,
           axes=FALSE,box=FALSE,zlim=c(0,3000))
      dev.off()
      
      # SDA (in km2)
      SDA.fut$area.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="full"] <- sum(values(caFut)>=1500,na.rm=TRUE)
      SDA.fut$area.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="zero"] <- sum(values(caZD)>=1500,na.rm=TRUE)
      
      # Altitude
      # fd
      if (SDA.fut$area.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="full"] != 0) {
        Alt.fut$mean.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="full"] <- mean(values(environ$alt)[values(caFut)>=1500],na.rm=TRUE)
        Alt.fut$q1.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="full"] <- quantile(values(environ$alt)[values(caFut)>=1500],0.025,na.rm=TRUE)
        Alt.fut$q2.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="full"] <- quantile(values(environ$alt)[values(caFut)>=1500],0.975,na.rm=TRUE)
      }
      # zd
      if (SDA.fut$area.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="zero"] != 0) {
        Alt.fut$mean.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="zero"] <- mean(values(environ$alt)[values(caZD)>=1500],na.rm=TRUE)
        Alt.fut$q1.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="zero"] <- quantile(values(environ$alt)[values(caZD)>=1500],0.025,na.rm=TRUE)
        Alt.fut$q2.fut[SDA.fut$rcp==rcp[j] & SDA.fut$yr==yr[l] & SDA.fut$disp=="zero"] <- quantile(values(environ$alt)[values(caZD)>=1500],0.975,na.rm=TRUE)
      }
    }
  }
  # SDA
  SDA.fut$perc.change <- round(100*((SDA.fut$area.fut-SDA.fut$area.pres)/SDA.fut$area.pres))
  write.table(SDA.fut,paste0(spdir,"/sda_fut.txt"),sep="\t")
  # Alt
  Alt.fut[,7:9] <- round_data_frame(Alt.fut[,7:9])
  Alt.fut$change <- Alt.fut$mean.fut-Alt.fut$mean.pres
  write.table(Alt.fut,paste0(spdir,"/alt_fut.txt"),sep="\t")
  
  ##========================================
  ## Save main objects
  save(list=c("SDA.fut","Alt.fut","niche","Perf.mods","VarImp"),file=paste0(spdir,"/main_obj.rda"))
  
}

# ===========
# End of file
# ===========
