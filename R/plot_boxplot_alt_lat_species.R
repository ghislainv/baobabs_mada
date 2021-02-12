#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================

library(grid)
require(ggplot2)

####################### Box plots according to each SDAp and SDAf compared

##### Non threatened baobab species ####
### A digitata ###

pred_dig <- stack(paste0("Adansonia.digitata/proj_current/proj_current_Adansonia.digitata_ensemble.grd"))
ca_dig <- pred_dig[[1]]
ca_test <- stack(ca_dig,environ$alt)
ca_test_na <- na.omit(ca_test)

ca_test_na2 <- rasterToPoints(ca_test_na)
ca_test_na3 <- as.data.frame(ca_test_na2[complete.cases(ca_test_na2), ] )

alegria <- ca_test_na3 %>% filter(Adansonia.digitata_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

alegria2 <- sample_n(alegria, size = 1000, replace = F)
alegria2$Proj <- rep(c("Present"),1000)
names(alegria2) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_dig <- raster(paste0("Adansonia.digitata/caFut_85_2080.tif"))
ca_test_fut <- stack(caZD_dig,environ$alt)
ca_test_na_fut <- na.omit(ca_test_fut)

ca_test_na2_fut <- rasterToPoints(ca_test_na_fut)
ca_test_na3_fut <- as.data.frame(ca_test_na2_fut[complete.cases(ca_test_na2_fut), ] )
alegria_fut <- ca_test_na3_fut %>% filter(caFut_85_2080 >= 1500)

alegria2_fut <- sample_n(alegria_fut, size = 1000, replace = F)
alegria2_fut$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test <- rbind(alegria2,alegria2_fut)

### Future testing 2050 
caZD_digi_2055 <- raster(paste0("Adansonia.digitata/caFut_85_2050.tif"))
ca_test_fut_digi_2055 <- stack(caZD_digi_2055,environ$alt)
ca_test_na_fut_digi_55 <- na.omit(ca_test_fut_digi_2055)

ca_test_na2_fut_digi_55 <- rasterToPoints(ca_test_na_fut_digi_55)
ca_test_na3_fut_digi_55 <- as.data.frame(ca_test_na2_fut_digi_55[complete.cases(ca_test_na2_fut_digi_55), ] )
alegria_fut_digi_55 <- ca_test_na3_fut_digi_55 %>% filter(caFut_85_2050 >= 1500)

alegria2_fut_digi_55 <- sample_n(alegria_fut_digi_55, size = 1000, replace = F)
alegria2_fut_digi_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_digi_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_digi_finale <- rbind(alegria_test, alegria2_fut_digi_55)


alegria_test_digi_finale$Scenario = factor(alegria_test_digi_finale$Scenario,
                                           levels = c("Present", "Future_2055","Future_2085"),
                                           labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot

latitude_digi <- ggplot(alegria_test_digi_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8650000 , size=7, label = "(b)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

## Altitude plot
altitude_digi <- ggplot(alegria_test_digi_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 800 , size=7, label = "(a)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

########
######## A grandidieri ##############

pred_grand <- stack(paste0("Adansonia.grandidieri/proj_current/proj_current_Adansonia.grandidieri_ensemble.grd"))
ca_grand <- pred_grand[[1]]
ca_test_grand <- stack(ca_grand,environ$alt)
ca_test_na_grand <- na.omit(ca_test_grand)

ca_test_na2_grand <- rasterToPoints(ca_test_na_grand)
ca_test_na3_grand <- as.data.frame(ca_test_na2_grand[complete.cases(ca_test_na2_grand), ] )

alegria_grand <- ca_test_na3_grand %>% filter(Adansonia.grandidieri_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

alegria2_grand <- sample_n(alegria_grand, size = 1000, replace = F)
alegria2_grand$Proj <- rep(c("Present"),1000)
names(alegria2_grand) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_grand <- raster(paste0("Adansonia.grandidieri/caFut_85_2080.tif"))
ca_test_fut_grand <- stack(caZD_grand,environ$alt)
ca_test_na_fut_grand <- na.omit(ca_test_fut_grand)

ca_test_na2_fut_grand <- rasterToPoints(ca_test_na_fut_grand)
ca_test_na3_fut_grand <- as.data.frame(ca_test_na2_fut_grand[complete.cases(ca_test_na2_fut_grand), ] )
alegria_fut_grand <- ca_test_na3_fut_grand %>% filter(caFut_85_2080 >= 1500)

alegria2_fut_grand <- sample_n(alegria_fut_grand, size = 1000, replace = F)
alegria2_fut_grand$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_grand) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_grand <- rbind(alegria2_grand,alegria2_fut_grand)


### Future testing 2050 
caZD_grand_2055 <- raster(paste0("Adansonia.grandidieri/caFut_85_2050.tif"))
ca_test_fut_grand_2055 <- stack(caZD_grand_2055,environ$alt)
ca_test_na_fut_grand_55 <- na.omit(ca_test_fut_grand_2055)

ca_test_na2_fut_grand_55 <- rasterToPoints(ca_test_na_fut_grand_55)
ca_test_na3_fut_grand_55 <- as.data.frame(ca_test_na2_fut_grand_55[complete.cases(ca_test_na2_fut_grand_55), ] )
alegria_fut_grand_55 <- ca_test_na3_fut_grand_55 %>% filter(caFut_85_2050 >= 1500)


alegria2_fut_grand_55 <- sample_n(alegria_fut_grand_55, size = 1000, replace = F)
alegria2_fut_grand_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_grand_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_grand_finale <- rbind(alegria2_fut_grand_55,alegria_test_grand)

alegria_test_grand_finale$Scenario = factor(alegria_test_grand_finale$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_grand <- ggplot(alegria_test_grand_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8250000 , size=7, label = "(d)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

## Altitude plot
altitude_grand <- ggplot(alegria_test_grand_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1050 , size=7, label = "(c)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

######## A za ###############

pred_za <- stack(paste0("Adansonia.za/proj_current/proj_current_Adansonia.za_ensemble.grd"))
ca_za <- pred_za[[1]]
ca_test_za <- stack(ca_za,environ$alt)
ca_test_na_za <- na.omit(ca_test_za)

ca_test_na2_za <- rasterToPoints(ca_test_na_za)
ca_test_na3_za <- as.data.frame(ca_test_na2_za[complete.cases(ca_test_na2_za), ] )

alegria_za <- ca_test_na3_za %>% filter(Adansonia.za_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

alegria2_za <- sample_n(alegria_za, size = 1000, replace = F)
alegria2_za$Proj <- rep(c("Present"),1000)
names(alegria2_za) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_za <- raster(paste0("Adansonia.za/caFut_85_2080.tif"))
ca_test_fut_za <- stack(caZD_za,environ$alt)
ca_test_na_fut_za <- na.omit(ca_test_fut_za)

ca_test_na2_fut_za <- rasterToPoints(ca_test_na_fut_za)
ca_test_na3_fut_za <- as.data.frame(ca_test_na2_fut_za[complete.cases(ca_test_na2_fut_za), ] )
alegria_fut_za <- ca_test_na3_fut_za %>% filter(caFut_85_2080 >= 1500)

alegria2_fut_za <- sample_n(alegria_fut_za, size = 1000, replace = F)
alegria2_fut_za$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_za) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_za <- rbind(alegria2_za,alegria2_fut_za)

### Future testing 2050 
caZD_za_2055 <- raster(paste0("Adansonia.za/caFut_85_2050.tif"))
ca_test_fut_za_2055 <- stack(caZD_za_2055,environ$alt)
ca_test_na_fut_za_55 <- na.omit(ca_test_fut_za_2055)

ca_test_na2_fut_za_55 <- rasterToPoints(ca_test_na_fut_za_55)
ca_test_na3_fut_za_55 <- as.data.frame(ca_test_na2_fut_za_55[complete.cases(ca_test_na2_fut_za_55), ] )
alegria_fut_za_55 <- ca_test_na3_fut_za_55 %>% filter(caFut_85_2050 >= 1500)

alegria2_fut_za_55 <- sample_n(alegria_fut_za_55, size = 1000, replace = F)
alegria2_fut_za_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_za_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_za_finale <- rbind(alegria2_fut_za_55,alegria_test_za)

alegria_test_za_finale$Scenario = factor(alegria_test_za_finale$Scenario,
                                         levels = c("Present", "Future_2055","Future_2085"),
                                         labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_za <- ggplot(alegria_test_za_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8650000 , size=7, label = "(f)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

## Altitude plot
altitude_za <- ggplot(alegria_test_za_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1350 , size=7, label = "(e)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

###########Threatened baobab species ###############
######## A madagascariensis ############

pred_mada <- stack(paste0("Adansonia.madagascariensis/proj_current/proj_current_Adansonia.madagascariensis_ensemble.grd"))
ca_mada <- pred_mada[[1]]
ca_test_mada <- stack(ca_mada,environ$alt)
ca_test_na_mada <- na.omit(ca_test_mada)

ca_test_na2_mada <- rasterToPoints(ca_test_na_mada)
ca_test_na3_mada <- as.data.frame(ca_test_na2_mada[complete.cases(ca_test_na2_mada), ] )

alegria_mada <- ca_test_na3_mada %>% filter(Adansonia.madagascariensis_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

alegria2_mada <- sample_n(alegria_mada, size = 1000, replace = F)
alegria2_mada$Proj <- rep(c("Present"),1000)
names(alegria2_mada) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caFut_mada <- raster(paste0("Adansonia.madagascariensis/caFut_85_2080.tif"))
ca_test_fut_mada <- stack(caFut_mada,environ$alt)
ca_test_na_fut_mada <- na.omit(ca_test_fut_mada)

ca_test_na2_fut_mada <- rasterToPoints(ca_test_na_fut_mada)
ca_test_na3_fut_mada <- as.data.frame(ca_test_na2_fut_mada[complete.cases(ca_test_na2_fut_mada), ] )
alegria_fut_mada <- ca_test_na3_fut_mada %>% filter(caFut_85_2080 >= 1500)

alegria2_fut_mada <- sample_n(alegria_fut_mada, size = 1000, replace = F)
alegria2_fut_mada$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_mada) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_mada <- rbind(alegria2_mada,alegria2_fut_mada)
tail(alegria_test_mada)

### Future testing 2050 
caFut_mada_2055 <- raster(paste0("Adansonia.madagascariensis/caFut_85_2050.tif"))
ca_test_fut_mada_2055 <- stack(caFut_mada_2055,environ$alt)
ca_test_na_fut_mada_55 <- na.omit(ca_test_fut_mada_2055)

ca_test_na2_fut_mada_55 <- rasterToPoints(ca_test_na_fut_mada_55)
ca_test_na3_fut_mada_55 <- as.data.frame(ca_test_na2_fut_mada_55[complete.cases(ca_test_na2_fut_mada_55), ] )
alegria_fut_mada_55 <- ca_test_na3_fut_mada_55 %>% filter(caFut_85_2050 >= 1500)


alegria2_fut_mada_55 <- sample_n(alegria_fut_mada_55, size = 1000, replace = F)
alegria2_fut_mada_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_mada_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_mada_finale <- rbind(alegria2_fut_mada_55,alegria_test_mada)

alegria_test_mada_finale$Scenario = factor(alegria_test_mada_finale$Scenario,
                                           levels = c("Present", "Future_2055","Future_2085"),
                                           labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_mada <- ggplot(alegria_test_mada_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8650000 , size=7, label = "(b)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +    
  theme(legend.position = "none")


## Altitude plot

altitude_mada <- ggplot(alegria_test_mada_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 800 , size=7, label = "(a)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

######## A perrieri ###############

pred_perrieri <- stack(paste0("Adansonia.perrieri/proj_current/proj_current_Adansonia.perrieri_ensemble.grd"))
ca_perrieri <- pred_perrieri[[1]]
ca_test_perrieri <- stack(ca_perrieri,environ$alt)
ca_test_na_perrieri <- na.omit(ca_test_perrieri)

ca_test_na2_perrieri <- rasterToPoints(ca_test_na_perrieri)
ca_test_na3_perrieri <- as.data.frame(ca_test_na2_perrieri[complete.cases(ca_test_na2_perrieri), ] )

alegria_perrieri <- ca_test_na3_perrieri %>% filter(Adansonia.perrieri_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

alegria2_perrieri <- sample_n(alegria_perrieri, size = 1000, replace = F)
alegria2_perrieri$Proj <- rep(c("Present"),1000)
names(alegria2_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 (RCP 4.5!!!) Atention!
caFut_perrieri <- raster(paste0("Adansonia.perrieri/caFut_45_2080.tif"))
ca_test_fut_perrieri <- stack(caFut_perrieri,environ$alt)
ca_test_na_fut_perrieri <- na.omit(ca_test_fut_perrieri)

ca_test_na2_fut_perrieri <- rasterToPoints(ca_test_na_fut_perrieri)
ca_test_na3_fut_perrieri <- as.data.frame(ca_test_na2_fut_perrieri[complete.cases(ca_test_na2_fut_perrieri), ] )
alegria_fut_perrieri <- ca_test_na3_fut_perrieri %>% filter(caFut_45_2080 >= 1500)

alegria2_fut_perrieri <- sample_n(alegria_fut_perrieri, size = 416, replace = F)
alegria2_fut_perrieri$Proj <- rep(c("Future_2085"),416)
names(alegria2_fut_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_perrieri <- rbind(alegria2_perrieri,alegria2_fut_perrieri)
tail(alegria_test_perrieri)

### Future testing 2050 
caFut_perrieri_2055 <- raster(paste0("Adansonia.perrieri/caFut_85_2050.tif"))
ca_test_fut_perrieri_2055 <- stack(caFut_perrieri_2055,environ$alt)
ca_test_na_fut_perrieri_55 <- na.omit(ca_test_fut_perrieri_2055)

ca_test_na2_fut_perrieri_55 <- rasterToPoints(ca_test_na_fut_perrieri_55)
ca_test_na3_fut_perrieri_55 <- as.data.frame(ca_test_na2_fut_perrieri_55[complete.cases(ca_test_na2_fut_perrieri_55), ] )
alegria_fut_perrieri_55 <- ca_test_na3_fut_perrieri_55 %>% filter(caFut_85_2050 >= 1500)


alegria2_fut_perrieri_55 <- sample_n(alegria_fut_perrieri_55, size = 1000, replace = F)
alegria2_fut_perrieri_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_perrieri_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_perrieri_finale <- rbind(alegria2_fut_perrieri_55,alegria_test_perrieri)

alegria_test_perrieri_finale$Scenario = factor(alegria_test_perrieri_finale$Scenario,
                                               levels = c("Present", "Future_2055","Future_2085"),
                                               labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_perrieri <- ggplot(alegria_test_perrieri_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8600000 , size=7, label = "(d)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +    
  theme(legend.position = "none")


## Altitude plot
altitude_perrieri <- ggplot(alegria_test_perrieri_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1450 , size=7, label = "(c)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

######## A rubrostipa ###############
pred_rubro <- stack(paste0("Adansonia.rubrostipa/proj_current/proj_current_Adansonia.rubrostipa_ensemble.grd"))
ca_rubro <- pred_rubro[[1]]
ca_test_rubro <- stack(ca_rubro,environ$alt)
ca_test_na_rubro <- na.omit(ca_test_rubro)

ca_test_na2_rubro <- rasterToPoints(ca_test_na_rubro)
ca_test_na3_rubro <- as.data.frame(ca_test_na2_rubro[complete.cases(ca_test_na2_rubro), ] )

alegria_rubro <- ca_test_na3_rubro %>% filter(Adansonia.rubrostipa_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

alegria2_rubro <- sample_n(alegria_rubro, size = 1000, replace = F)
alegria2_rubro$Proj <- rep(c("Present"),1000)
names(alegria2_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caFut_rubro <- raster(paste0("Adansonia.rubrostipa/caFut_85_2080.tif"))
ca_test_fut_rubro <- stack(caFut_rubro,environ$alt)
ca_test_na_fut_rubro <- na.omit(ca_test_fut_rubro)

ca_test_na2_fut_rubro <- rasterToPoints(ca_test_na_fut_rubro)
ca_test_na3_fut_rubro <- as.data.frame(ca_test_na2_fut_rubro[complete.cases(ca_test_na2_fut_rubro), ] )
alegria_fut_rubro <- ca_test_na3_fut_rubro %>% filter(caFut_85_2080 >= 1500)

alegria2_fut_rubro <- sample_n(alegria_fut_rubro, size = 1000, replace = F)
alegria2_fut_rubro$Proj <- rep(c("Future_2085"),1000)
names(alegria2_fut_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_rubro <- rbind(alegria2_rubro,alegria2_fut_rubro)

### Future testing 2050 
caFut_rubro_2055 <- raster(paste0("Adansonia.rubrostipa/caFut_85_2050.tif"))
ca_test_fut_rubro_2055 <- stack(caFut_rubro_2055,environ$alt)
ca_test_na_fut_rubro_55 <- na.omit(ca_test_fut_rubro_2055)

ca_test_na2_fut_rubro_55 <- rasterToPoints(ca_test_na_fut_rubro_55)
ca_test_na3_fut_rubro_55 <- as.data.frame(ca_test_na2_fut_rubro_55[complete.cases(ca_test_na2_fut_rubro_55), ] )
alegria_fut_rubro_55 <- ca_test_na3_fut_rubro_55 %>% filter(caFut_85_2050 >= 1500)


alegria2_fut_rubro_55 <- sample_n(alegria_fut_rubro_55, size = 1000, replace = F)
alegria2_fut_rubro_55$Proj <- rep(c("Future_2055"),1000)
names(alegria2_fut_rubro_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_rubro_finale <- rbind(alegria2_fut_rubro_55,alegria_test_rubro)

alegria_test_rubro_finale$Scenario = factor(alegria_test_rubro_finale$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_rubro <- ggplot(alegria_test_rubro_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8600000 , size=7, label = "(f)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")


## Altitude plot
altitude_rubro <- ggplot(alegria_test_rubro_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 1000 , size=7, label = "(e)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +  
  theme(legend.position = "none")

######## A suarezensis ###############

pred_suare <- stack(paste0("Adansonia.suarezensis/proj_current/proj_current_Adansonia.suarezensis_ensemble.grd"))
ca_suare <- pred_suare[[1]]
ca_test_suare <- stack(ca_suare,environ$alt)
ca_test_na_suare <- na.omit(ca_test_suare)

ca_test_na2_suare <- rasterToPoints(ca_test_na_suare)
ca_test_na3_suare <- as.data.frame(ca_test_na2_suare[complete.cases(ca_test_na2_suare), ] )

alegria_suare <- ca_test_na3_suare %>% filter(Adansonia.suarezensis_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

alegria2_suare <- sample_n(alegria_suare, size = 1000, replace = F)
alegria2_suare$Proj <- rep(c("Present"),1000)
names(alegria2_suare) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 Atention RCP 4.5 (SAME AS A. perrieri)
caFut_suare <- raster(paste0("Adansonia.suarezensis/caFut_45_2080.tif"))
ca_test_fut_suare <- stack(caFut_suare,environ$alt)
ca_test_na_fut_suare <- na.omit(ca_test_fut_suare)

ca_test_na2_fut_suare <- rasterToPoints(ca_test_na_fut_suare)
ca_test_na3_fut_suare <- as.data.frame(ca_test_na2_fut_suare[complete.cases(ca_test_na2_fut_suare), ] )
alegria_fut_suare <- ca_test_na3_fut_suare %>% filter(caFut_45_2080 >= 1500)

alegria2_fut_suare <- sample_n(alegria_fut_suare, size = 105, replace = F)
alegria2_fut_suare$Proj <- rep(c("Future_2085"),105)
names(alegria2_fut_suare) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_suare <- rbind(alegria2_suare,alegria2_fut_suare)

### Future testing 2050 
caFut_suare_2055 <- raster(paste0("Adansonia.suarezensis/caFut_45_2050.tif"))
ca_test_fut_suare_2055 <- stack(caFut_suare_2055,environ$alt)
ca_test_na_fut_suare_55 <- na.omit(ca_test_fut_suare_2055)

ca_test_na2_fut_suare_55 <- rasterToPoints(ca_test_na_fut_suare_55)
ca_test_na3_fut_suare_55 <- as.data.frame(ca_test_na2_fut_suare_55[complete.cases(ca_test_na2_fut_suare_55), ] )
alegria_fut_suare_55 <- ca_test_na3_fut_suare_55 %>% filter(caFut_45_2050 >= 1500)

alegria2_fut_suare_55 <- sample_n(alegria_fut_suare_55, size = 15, replace = F)
alegria2_fut_suare_55$Proj <- rep(c("Future_2055"),15)
names(alegria2_fut_suare_55) <- c("Long","Lat","Prediction","Alt","Scenario")
alegria_test_suare_finale <- rbind(alegria2_fut_suare_55,alegria_test_suare)

alegria_test_suare_finale$Scenario = factor(alegria_test_suare_finale$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_suare <- ggplot(alegria_test_suare_finale) +
  aes(y= Lat, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Latitude (UTM)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 8670000 , size=7, label = "(h)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +    
  theme(legend.position = "none")


## Altitude plot
altitude_suare <- ggplot(alegria_test_suare_finale) +
  aes(y= Alt, x= Scenario) +
  geom_jitter(alpha =.5,
              height = 0,
              width = .25) +
  aes(col = Scenario) +
  geom_boxplot(alpha= .25) +
  aes(fill= Scenario) +  
  scale_colour_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  scale_fill_manual(values = c("#31688EFF","#440154FF","#6DCD59FF")) +
  xlab("") +
  ylab("Elevation (m)") +
  theme_bw() +
  labs(col = "") +
  annotate("text", x  = 0.5, y = 950 , size=7, label = "(g)") +
  theme(axis.title.x = element_text(size = rel(1.75),colour="black")) +
  theme(axis.title.y = element_text(size = rel(1.75),colour="black")) +
  theme(axis.text.x = element_text(size = rel(1.65),colour="black")) +
  theme(axis.text.y = element_text(size = rel(1.65), colour="black")) +   
  theme(legend.position = "none")

## Legend names for plot
grob_digitata <- textGrob("A. digitata",rot=90, gp=gpar(cex=2,fontface="italic"),
                          hjust=0.3, vjust=2)
grob_grand <- textGrob("A. grandidieri",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=2)
grob_mada <- textGrob("A. madagascariensis",rot=90, gp=gpar(cex=2,fontface="italic"),
                      hjust=0.4, vjust=2)
grob_perri <- textGrob("A. perrieri",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=2)
grob_rubro <- textGrob("A. rubrostipa",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=2)
grob_suare <- textGrob("A. suarezensis",rot=90, gp=gpar(cex=2,fontface="italic"),
                       hjust=0.3, vjust=2)
grob_za <- textGrob("A. za",rot=90, gp=gpar(cex=2,fontface="italic"),
                    hjust=0.3, vjust=2)
## Combine plots
### Threatened species
day_1 <- rbind(c(1,rep(seq(5,6,by=1),each=5)),
               c(2,rep(seq(7,8,by=1),each=5)),
               c(3,rep(seq(9,10,by=1),each=5)),
               c(4,rep(seq(11,12,by=1),each=5)))

### Non-threatened species
day_2 <- rbind(c(1,rep(seq(5,6,by=1),each=5)),
               c(2,rep(seq(7,8,by=1),each=5)),
               c(3,rep(seq(9,10,by=1),each=5)))             

#### Threatened species
a5 <- grid.arrange(grob_mada, grob_perri,
                   grob_rubro,grob_suare,
                   altitude_mada,latitude_mada,altitude_perrieri,latitude_perrieri,
                   altitude_rubro,latitude_rubro,altitude_suare,latitude_suare,
                   layout_matrix = day_1)

#### Non-threatened species

a6 <- grid.arrange(grob_digitata, grob_grand,
                   grob_za,
                   altitude_digi,latitude_digi,altitude_grand,latitude_grand,
                   altitude_za,latitude_za,layout_matrix = day_2)

ggsave(file=paste0("./outputs/lat_ele_all_species_within_SDAcf_threat.png"),
       plot=a5,width=18,height=15,dpi="print")

ggsave(file=paste0("./outputs/lat_ele_all_species_within_SDAcf_non_threat.png"),
       plot=a6,width=18,height=15,dpi="print")

# ===========
# End of file
# ===========