#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

## Loop on species
for (i in 1:length(sp.dir)) {
  ## Draw points in the SDA and extract future environmental variables in the SDA
  pred <-  stack(paste0(sp.dir[i],"/proj_current/proj_current_",sp.dir[i],"_ensemble.grd"))
  ca <- pred[[1]]  
  wC.anomalies <- which(values(ca)>=500) 
  nC.anomalies <- length(wC.anomalies)
  Samp.anomalies <- if (nC.anomalies>1000) {sample(wC.anomalies,1000,replace=FALSE)} else {wC.anomalies}
  mapmat.df.anomalies <- as.data.frame(var_85_2080)[Samp.anomalies,] 
  
  ## Table to compare current and future climate in the SDA of each species 
  names(mapmat.df.anomalies) <- c("tmeanf","tseasf","precf","cwdf")
  mapmat.df <- read.csv(file=paste0("./",sp.dir[i],"/","niche_graph_species.csv"),sep=";")
  mapmat.final <- cbind(mapmat.df.anomalies,mapmat.df)
  mapmat.f <- mapmat.final[,c(1,2,3,4,6,7,8,9,10,11)]
  write.csv2(mapmat.final,paste0("./",sp.dir[i],"/","niche_graph_species_compared_anomaly.csv"))
  write.table(mapmat.final,paste0("./",sp.dir[i],"/niche_graph_species_compared_anomaly.txt"),sep="\t")
  
  #### Future climate in current SDA
  wC.future <- which(values(ca)>=500)
  mapmat.df.future <- as.data.frame(var_85_2080)[wC.future,]
  names(mapmat.df.future) <- c("tmeanf","tseasf","precf","cwdf")
  Mean_ok_future_fut <- round(apply(mapmat.df.anomalies,2,mean,na.rm=TRUE))
  q_ok_future_fut <- round(apply(mapmat.df.anomalies,2,quantile,c(0.025,0.975),na.rm=TRUE))
  niche_ok_future <- as.data.frame(rbind(Mean_ok_future_fut,q_ok_future_fut))
  write.csv2(niche_ok_future,paste0("./",sp.dir[i],"/","mean_niche_with_future.csv"))
  write.table(niche_ok_future,paste0("./",sp.dir[i],"/mean_niche_with_future.txt"),sep="\t")
}

###################################################################################
######## Generating variables histograms graphs for all seven baobab species ######
###################################################################################

## Read species environmental data
database1_suare <- read.csv(file=paste0("Adansonia.suarezensis/","niche_graph_species.csv"), header=T,sep=";")
database2_digi <- read.csv(file=paste0("Adansonia.digitata/","niche_graph_species.csv"), header=T,sep=";")
database3_grand <- read.csv(file=paste0("Adansonia.grandidieri/","niche_graph_species.csv"), header=T,sep=";")
database4_mada <- read.csv(file=paste0("Adansonia.madagascariensis/","niche_graph_species.csv"), header=T,sep=";")
database5_perri <- read.csv(file=paste0("Adansonia.perrieri/","niche_graph_species.csv"), header=T,sep=";")
database6_rubro <- read.csv(file=paste0("Adansonia.rubrostipa/","niche_graph_species.csv"), header=T,sep=";")
database7_za <- read.csv(file=paste0("Adansonia.za/","niche_graph_species.csv"), header=T,sep=";")
data_teste <- rbind(database2_digi,database3_grand,
              database4_mada,database5_perri,database6_rubro,database1_suare,database7_za)
levels(data_teste$species)
quantiles <- data_teste[,c(2,3,4,5,6)]
quantiles <- round(apply(quantiles,2,quantile,c(0.025,0.975),na.rm=TRUE)) #use it to set new plot

data_teste$species <- factor(data_teste$species,
                            levels = c("Adansonia.digitata","Adansonia.grandidieri",
                                       "Adansonia.madagascariensis","Adansonia.perrieri",
                                       "Adansonia.rubrostipa","Adansonia.suarezensis",
                                       "Adansonia.za"),
                            labels = c("A. digitata","A. grandidieri",
                                       "A. madagascariensis","A. perrieri",
                                       "A. rubrostipa","A. suarezensis",
                                       "A. za"))

## ggplot theme parameters
theme_par <-theme(axis.title.x = element_text(size = rel(1.8),colour="black"),
                  axis.title.y = element_text(size = rel(1.8),colour="black"),
                  axis.text.x = element_text(size = rel(1.6),colour="black"),
                  axis.text.y = element_text(size = rel(1.6), colour="black"))

# Density plots ### Temp. Seasonality
range(data_teste$tseas)
# generate break positions
breaks <- c(round(seq(min(data_teste$tseas),max(data_teste$tseas),length=6)))
# and labels
labels = as.character(breaks) # labels must be the same as breaks, otherwise error
# plot the map
my_plot = ggplot(data_teste, aes(x=tseas, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + #specific palette
  geom_vline(xintercept = range(data_teste$tseas), show.legend = T, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$tseas), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))
my_plot_tseas <-  my_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 #legend.text = element_text(face= "italic",size=15),
                                 legend.title=element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Temp. Seasonality (sd x 100 ºC)", y = "Density",size=5) +
  theme(legend.position="none") + theme_par
  
# Density plots ### Annual Mean Temperature
# generate break positions
breaks <- c(round(seq(min(data_teste$tmean),max(data_teste$tmean),length=6)))
labels <- as.character(breaks)
my_plot2 <- ggplot(data_teste, aes(x=tmean, color=species)) + 
  geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = range(data_teste$tmean), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$tmean), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001))

my_plot_tmean <-  my_plot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                  legend.text = element_text(face= "italic",size=12),
                                  legend.title=element_blank(),
                                  panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Mean Annual Temperature (ºC x 10) ", y = "Density",size=5) +
  theme(legend.key.size = unit(0.5, "cm")) + theme_par
library(lemon)  
my_plot_tmean <- lemon::reposition_legend(my_plot_tmean, 'top left', offset = 0.05)

# Density plots ### Annual Mean Precipitation
# generate break positions
breaks <- c(round(seq(min(data_teste$prec),max(data_teste$prec),length=6)))
labels <- as.character(breaks)
# plot
my_plot3 <- ggplot(data_teste, aes(x=prec, color=species)) + geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = range(data_teste$prec), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$prec), breaks = breaks, labels = labels)
my_plot_prec <-  my_plot3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             #legend.text = element_text(face= "italic",size=15),
                             legend.title=element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "Density",size=5) +
  theme(legend.position="none") + theme_par
  
# Density plots ### Climatic Water Deficit
# generate break positions
breaks <- c(round(seq(min(data_teste$cwd),max(data_teste$cwd),length=6)))
labels = as.character(breaks)
my_plot4 = ggplot(data_teste, aes(x=cwd, color=species)) + geom_density(size=1.3)+
  scale_color_brewer(palette = "Set1") + 
  geom_vline(xintercept = range(data_teste$cwd), show.legend = F, colour="red", linetype="dashed") +
  scale_x_continuous(limits = range(data_teste$cwd), breaks = breaks, labels = labels)
my_plot_cwd =  my_plot4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             #legend.text = element_text(face= "italic",size=15),
                             legend.title=element_blank(),
                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Climatic Water Deficit (mm)", y = "Density",size=5) +
  theme(legend.position="none") + theme_par
  
## Combine plots 
lay_2 <- cbind(c(rep(seq(1,4,by=1),each=3)),
               c(rep(seq(2,4,by=2),each=3)))
colnames(lay_2) <-c("a","b")
lay_2 <- subset( lay_2, select = -b )

plot_densities <- grid.arrange(my_plot_tmean, my_plot_tseas, my_plot_prec, my_plot_cwd,
                               layout_matrix=lay_2)

ggsave(filename="outputs/Fig_ap_sps_niche.png", plot=plot_densities,
       width=textwidth, height=textwidth*(11/8), scale=1.2, dpi=300, unit="cm")



### Comparing bioclimatic niche of each species inside current SDA according to the 
### most important variables ######
# A_suarezensis
# Import data_set
a.suare <- read.csv(file=paste0("Adansonia.suarezensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
# 1st importance variable - Tseas
range(a.suare$tseas) 
range(a.suare$tseasf) 
mean(a.suare$tseas) 
mean(a.suare$tseasf) 
breaks <- c(round(seq(min(a.suare$tseas),max(a.suare$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_suar = ggplot(a.suare, aes(x=tseas, y=..density..)) + 
  geom_density(aes(tseas,fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(1175,1550), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_suar =  plot.seas_suar + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                               legend.text = element_blank(),
                               panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Temp. Seasonality (ºC sd x 100)", y = "A. suarezensis",size=5) +
  labs(col = "") +
  annotate("text", x  = 1190, y = 0.0105 , size=7, label = "(k)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


### second importance variable - Prec
range(a.suare$prec) 
range(a.suare$precf) 
breaks <- c(round(seq(min(a.suare$precf),max(a.suare$prec),length=8)))
labels = as.character(breaks)

plot.prec_suar = ggplot(a.suare, aes(x=prec, y=..density..)) + 
  geom_density(aes(prec,fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(800,1600), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_suar =  plot.prec_suar + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "",size=5) +
  labs(col = "") +
  annotate("text", x  = 825, y = 0.0032 , size=7, label = "(l)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")
  
   
################################
# Import data_set A perrieri
#1st importance variable - Tseas 
a.perrieri <- read.csv(file=paste0("Adansonia.perrieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.perrieri)
range(a.perrieri$tseas) 
range(a.perrieri$tseasf) 
round(mean(a.perrieri$tseas))
round(mean(a.perrieri$tseasf))
breaks <- c(round(seq(min(a.perrieri$tseas),max(a.perrieri$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_per = ggplot(a.perrieri, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(750, 2250), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_perrieri =  plot.seas_per + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        legend.text = element_blank(),
                                        panel.background = element_blank(), 
                                        axis.line = element_line(colour = "black"))+
  
  labs(x="Temp. Seasonality (ºC sd x 100)",  y = "A. perrieri",size=5) +
  labs(col = "") +
  annotate("text", x  = 805, y = 0.0025 , size=7, label = "(g)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

## 2nd Most IV A perrieri - climatic water deficit
breaks <- c(round(seq(min(a.perrieri$cwd),max(a.perrieri$cwdf),length=8)))
labels = as.character(breaks)

plot.cwd_perrieri = ggplot(a.perrieri, aes(x=cwd, y=..density..)) + 
  geom_density(aes(cwd, fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(30, 1700), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.cwd_perrieri =  plot.cwd_perrieri + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Climatic Water Deficit (mm)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 80, y = 0.0020 , size=7, label = "(h)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

############## A.rubrostipa
# Import dataset
a.rubrostipa <- read.csv(file=paste0("Adansonia.rubrostipa/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

## 1st importance variable - Climatic Water Deficit
range(a.rubrostipa$cwd)
range(a.rubrostipa$cwdf) 
mean(a.rubrostipa$cwd) 
mean(a.rubrostipa$cwdf) 
breaks <- c(round(seq(min(a.rubrostipa$cwd),max(a.rubrostipa$cwdf),length=6)))
labels = as.character(breaks)

plot.cwd_rubro = ggplot(a.rubrostipa, aes(x=cwd, y=..density..)) + 
  geom_density(aes(fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(650, 2100), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.cwd_rubro =  plot.cwd_rubro + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             legend.text = element_blank(),
                             panel.background = element_blank(), 
                             axis.line = element_line(colour = "black"))+
  
  labs(x="Climatic Water Deficit (mm)", y = "A. rubrostipa",size=5) +
  labs(col = "") +
  annotate("text", x  = 698, y = 0.0065 , size=7, label = "(i)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

### 2nd Importance Variable - Prec
range(a.rubrostipa$prec) 
round(range(a.rubrostipa$precf)) 
mean(a.rubrostipa$prec) 
mean(a.rubrostipa$precf) 
breaks <- c(round(seq(min(a.rubrostipa$prec),max(a.rubrostipa$precf),length=8)))
labels = as.character(breaks)

plot.prec_rubro = ggplot(a.rubrostipa, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(150, 1900), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_rubro =  plot.prec_rubro + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             legend.text = element_blank(),
                                             panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 200, y = 0.00125 , size=7, label = "(j)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


############## A. madagascariensis 
## Import dataset
a.madagascariensis <- read.csv(file=paste0("Adansonia.madagascariensis/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

# 1st importance variable - Tseas
range(a.madagascariensis$tseas) 
range(a.madagascariensis$tseasf) 

breaks <- c(round(seq(min(a.madagascariensis$tseas),max(a.madagascariensis$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_mada = ggplot(a.madagascariensis, aes(x=tseas, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tseas'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(800,2300), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_mada =  plot.seas_mada + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                         legend.text = element_blank(),
                                         panel.background = element_blank(), 
                                         axis.line = element_line(colour = "black"))+
  
  labs(x="Temp. Seasonality (ºC sd x 100)", y = "A. madagascariensis",size=5) +
  labs(col = "") +
  annotate("text", x  = 875, y = 0.0023 , size=7, label = "(e)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.5),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


### 2nd importance variable - Tmean
range(a.madagascariensis$tmean)
range(a.madagascariensis$tmeanf) 
mean(a.madagascariensis$tmean) 
mean(a.madagascariensis$tmeanf) 

breaks <- c(round(seq(min(a.madagascariensis$tmean),max(a.madagascariensis$tmeanf),length=8)))
labels = as.character(breaks)

plot.mean_mada = ggplot(a.madagascariensis, aes(x=tmean, y=..density..)) + 
  geom_density(aes(fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(240,320), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.mean_mada =  plot.mean_mada + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 legend.text = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Annual Mean Temperature (ºC x 10)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 241.5, y = 0.059 , size=7, label = "(f)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

############## A. grandidieri
## Import dataset
a.grandidieri <- read.csv(file=paste0("Adansonia.grandidieri/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")

# 1st Importance Variable - Prec
range(a.grandidieri$prec) 
range(a.grandidieri$precf) 
mean(a.grandidieri$prec)
mean(a.grandidieri$precf)
breaks <- c(round(seq(min(a.grandidieri$precf),max(a.grandidieri$prec),length=8)))
labels = as.character(breaks)

plot.prec_grand = ggplot(a.grandidieri, aes(x=prec, y=..density..)) + 
  geom_density(aes(fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) + 
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(300,1100), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_grand =  plot.prec_grand + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           legend.text = element_blank(),
                                           panel.background = element_blank(), 
                                           axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "A. grandidieri",size=5) +
  labs(col = "") +
  annotate("text", x  = 336, y = 0.005 , size=7, label = "(c)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


## 2nd Importance Variable - Tmean
range(a.grandidieri$tmean)
range(a.grandidieri$tmeanf)
mean(a.grandidieri$tmean)
mean(a.grandidieri$tmeanf)
breaks <- c(round(seq(min(a.grandidieri$tmean),max(a.grandidieri$tmeanf),length=8)))
labels = as.character(breaks)

plot.mean_grand = ggplot(a.grandidieri, aes(x=tmean, y=..density..)) + 
  geom_density(aes(tmean,fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(237, 305), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.mean_grand =  plot.mean_grand + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                           legend.text = element_blank(),
                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Annual Mean Temperature (ºC x 10)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 238, y = 0.08 , size=7, label = "(d)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


############## A za
## Import dataset
a.za <- read.csv(file=paste0("Adansonia.za/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.za)
# 1st Importance Variable - Prec
range(a.za$prec)
range(a.za$precf)

breaks <- c(round(seq(min(a.za$precf),max(a.za$prec),length=8)))
labels = as.character(breaks)

plot.prec_za = ggplot(a.za, aes(x=prec, y=..density..)) + 
  geom_density(aes(prec,fill=species,color= 'prec'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(prec)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(prec,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(prec,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(precf, fill=species,color='precf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(precf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(precf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(precf,0.025)),color="black",linetype="dashed",size=1) +
  
  scale_x_continuous(limits = c(200,1950), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('prec' = 'darkorange', 'precf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.prec_za =  plot.prec_za + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                     legend.text = element_blank(),
                                     panel.background = element_blank(), 
                                     axis.line = element_line(colour = "black"))+
  
  labs(x="Mean Annual Precipitation (mm.y-¹)", y = "A. za",size=5) +
  labs(col = "") +
  annotate("text", x  = 250, y = 0.0019 , size=7, label = "(m)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

## 2nd Importance Variable - Tmean
range(a.za$tmean)
range(a.za$tmeanf)
mean(a.za$tmean)
mean(a.za$tmeanf)
breaks <- c(round(seq(min(a.za$tmean),max(a.za$tmeanf),length=8)))
labels = as.character(breaks)

plot.mean_za = ggplot(a.za, aes(x=tmean, y=..density..)) + 
  geom_density(aes(tmean,fill=species,color= 'tmean'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(tmean)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmean,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmean,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(tmeanf, fill=species,color='tmeanf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(tmeanf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tmeanf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tmeanf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(200, 320), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tmean' = 'darkorange', 'tmeanf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.mean_za =  plot.mean_za + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                     legend.text = element_blank(),
                                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Annual Mean Temperature (ºC x 10)", y="",size=5) +
  labs(col = "") +
  annotate("text", x  = 204, y = 0.03 , size=7, label = "(n)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")

##############A digitata 
## Import dataset
a.digitata <- read.csv(file=paste0("Adansonia.digitata/niche_graph_species_compared_anomaly.csv"), header=T,sep=";",dec=",")
head(a.digitata)

##1st Importance Variable - Tseas

# range tseas 
range(a.digitata$tseas) 
range(a.digitata$tseasf)
mean(a.digitata$tseas)
mean(a.digitata$tseasf) 
breaks <- c(round(seq(min(a.digitata$tseas),max(a.digitata$tseasf),length=8)))
labels = as.character(breaks)

plot.seas_dig = ggplot(a.digitata, na.rm=T, aes(x=tseas, y=..density..)) + 
  
  geom_density(aes(tseas, fill=species,color= 'tseas'), alpha=.5, na.rm=T) + 
  geom_vline(aes(xintercept=mean(tseas)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseas,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseas,0.025)),color="darkorange",linetype="dashed",size=1) +
  geom_density(aes(tseasf, fill=species,color='tseasf'),alpha=.5) +
  
  geom_vline(aes(xintercept=mean(tseasf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(tseasf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(tseasf,0.025)),color="black",linetype="dashed",size=1) +
  scale_x_continuous(limits = c(800,2700), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('tseas' = 'darkorange', 'tseasf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

plot.seas_dig =  plot.seas_dig + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                       legend.text = element_blank(),
                                       panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  
  labs(x="Temp. Seasonality (ºC sd x 100)",y="A. digitata",size=5) +
  labs(col = "") +
  annotate("text", x  = 882, y = 0.004 , size=7, label = "(a)") +
  theme(axis.title.x = element_text(size = rel(2),colour="black")) +
  theme(axis.title.y = element_text(size = rel(2),colour="black",face="italic")) +
  theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
  theme(legend.position="none")


## 2nd Importance Variable - Climatic Water Deficit
range(a.digitata$cwd) 
range(a.digitata$cwdf) 
mean(a.digitata$cwd)
mean(a.digitata$cwdf) 
breaks <- c(round(seq(min(a.digitata$cwd),max(a.digitata$cwdf),length=8)))
labels = as.character(breaks)

plot.cwd_dig = ggplot(a.digitata, aes(x=cwd, y=..density..)) + 
  
  geom_density(aes(cwd,fill=species,color= 'cwd'), alpha=.5) + 
  geom_vline(aes(xintercept=mean(cwd)),color="darkorange",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwd,0.975)),color="darkorange",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwd,0.025)),color="darkorange",linetype="dashed",size=1) +
  
  geom_density(aes(cwdf, fill=species,color='cwdf'),alpha=.5) +
  geom_vline(aes(xintercept=mean(cwdf)),color="black",linetype="solid",size=1)+
  geom_vline(aes(xintercept=quantile(cwdf,0.975)),color="black",linetype="dashed",size=1) +
  geom_vline(aes(xintercept=quantile(cwdf,0.025)),color="black",linetype="dashed",size=1) +

  scale_x_continuous(limits = c(600, 2100), breaks = breaks, labels = labels) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  scale_color_manual(values = c('cwd' = 'darkorange', 'cwdf' = 'black'))+
  scale_fill_viridis(alpha=1,begin= 0.65,end=0.7, discrete=T,option="D")

  plot.cwd_dig =  plot.cwd_dig + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                        legend.text = element_blank(),
                                        panel.background = element_blank(), 
                                        axis.line = element_line(colour = "black"))+
    
  labs(x="Climatic Water Deficit (mm)", y = "",size=5) +
         labs(col = "") +
         annotate("text", x  = 628, y = 0.007 , size=7, label = "(b)") +
         theme(axis.title.x = element_text(size = rel(2),colour="black")) +
         theme(axis.title.y = element_text(size = rel(2),colour="black")) +
         theme(axis.text.x = element_text(size = rel(1.75),colour="black")) +
         theme(axis.text.y = element_text(size = rel(1.75), colour="black")) +
         theme(legend.position="none")

## Combine plots 
lay_3 <- rbind(c(rep(seq(1,2,by=1),each=3)),
               c(rep(seq(3,4,by=1),each=3)),
               c(rep(seq(5,6,by=1),each=3)),
               c(rep(seq(7,8,by=1),each=3)),
               c(rep(seq(9,10,by=1),each=3)),
               c(rep(seq(11,12,by=1),each=3)),
               c(rep(seq(13,14,by=1),each=3)))

plot_densities_curves <- grid.arrange(plot.seas_dig,plot.cwd_dig,
                                      plot.prec_grand,plot.mean_grand,
                                      plot.seas_mada,plot.mean_mada,
                                      plot.seas_perrieri,plot.cwd_perrieri,
                                      plot.cwd_rubro,plot.prec_rubro,
                                      plot.seas_suar, plot.prec_suar,
                                      plot.prec_za,plot.mean_za,
                                      layout_matrix=lay_3)

ggsave(file=paste0("./outputs/all_species_current_future_niche_comparison.png"),
       plot=plot_densities_curves,width=18,height=15,dpi=300,scale=1)

# ===========
# End of file
# ===========
