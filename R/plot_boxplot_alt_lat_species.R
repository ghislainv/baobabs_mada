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

a_dig <- ca_test_na3 %>% filter(Adansonia.digitata_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

a_dig_2 <- sample_n(a_dig, size = 1000, replace = F)
a_dig_2$Proj <- rep(c("Present"),1000)
names(a_dig_2) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_dig <- raster(paste0("Adansonia.digitata/caFut_85_2080.tif"))
ca_test_fut <- stack(caZD_dig,environ$alt)
ca_test_na_fut <- na.omit(ca_test_fut)

ca_test_na2_fut <- rasterToPoints(ca_test_na_fut)
ca_test_na3_fut <- as.data.frame(ca_test_na2_fut[complete.cases(ca_test_na2_fut), ] )
a_dig_fut <- ca_test_na3_fut %>% filter(caFut_85_2080 >= 1500)

a_dig_2_fut <- sample_n(a_dig_fut, size = 1000, replace = F)
a_dig_2_fut$Proj <- rep(c("Future_2085"),1000)
names(a_dig_2_fut) <- c("Long","Lat","Prediction","Alt","Scenario")
combine <- rbind(a_dig_2,a_dig_2_fut)

### Future testing 2050 
caZD_digi_2055 <- raster(paste0("Adansonia.digitata/caFut_85_2050.tif"))
ca_test_fut_digi_2055 <- stack(caZD_digi_2055,environ$alt)
ca_test_na_fut_digi_55 <- na.omit(ca_test_fut_digi_2055)

ca_test_na2_fut_digi_55 <- rasterToPoints(ca_test_na_fut_digi_55)
ca_test_na3_fut_digi_55 <- as.data.frame(ca_test_na2_fut_digi_55[complete.cases(ca_test_na2_fut_digi_55), ] )
a_dig_fut_55 <- ca_test_na3_fut_digi_55 %>% filter(caFut_85_2050 >= 1500)

a_dig_2_fut_55 <- sample_n(a_dig_fut_55, size = 1000, replace = F)
a_dig_2_fut_55$Proj <- rep(c("Future_2055"),1000)
names(a_dig_2_fut_55) <- c("Long","Lat","Prediction","Alt","Scenario")
final_df <- rbind(combine, a_dig_2_fut_55)


final_df$Scenario = factor(final_df$Scenario,
                                           levels = c("Present", "Future_2055","Future_2085"),
                                           labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot

latitude_digi <- ggplot(final_df) +
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
altitude_digi <- ggplot(final_df) +
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

data_grand <- ca_test_na3_grand %>% filter(Adansonia.grandidieri_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

data_2_grand <- sample_n(data_grand, size = 1000, replace = F)
data_2_grand$Proj <- rep(c("Present"),1000)
names(data_2_grand) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_grand <- raster(paste0("Adansonia.grandidieri/caFut_85_2080.tif"))
ca_test_fut_grand <- stack(caZD_grand,environ$alt)
ca_test_na_fut_grand <- na.omit(ca_test_fut_grand)

ca_test_na2_fut_grand <- rasterToPoints(ca_test_na_fut_grand)
ca_test_na3_fut_grand <- as.data.frame(ca_test_na2_fut_grand[complete.cases(ca_test_na2_fut_grand), ] )
data_fut_grand <- ca_test_na3_fut_grand %>% filter(caFut_85_2080 >= 1500)

data_2_fut_grand <- sample_n(data_fut_grand, size = 1000, replace = F)
data_2_fut_grand$Proj <- rep(c("Future_2085"),1000)
names(data_2_fut_grand) <- c("Long","Lat","Prediction","Alt","Scenario")
data_1_grand <- rbind(data_2_grand,data_2_fut_grand)


### Future testing 2050 
caZD_grand_2055 <- raster(paste0("Adansonia.grandidieri/caFut_85_2050.tif"))
ca_test_fut_grand_2055 <- stack(caZD_grand_2055,environ$alt)
ca_test_na_fut_grand_55 <- na.omit(ca_test_fut_grand_2055)

ca_test_na2_fut_grand_55 <- rasterToPoints(ca_test_na_fut_grand_55)
ca_test_na3_fut_grand_55 <- as.data.frame(ca_test_na2_fut_grand_55[complete.cases(ca_test_na2_fut_grand_55), ] )
data_fut_grand_55 <- ca_test_na3_fut_grand_55 %>% filter(caFut_85_2050 >= 1500)


data_2_fut_grand_55 <- sample_n(data_fut_grand_55, size = 1000, replace = F)
data_2_fut_grand_55$Proj <- rep(c("Future_2055"),1000)
names(data_2_fut_grand_55) <- c("Long","Lat","Prediction","Alt","Scenario")
data_final_grand <- rbind(data_2_fut_grand_55,data_1_grand)

data_final_grand$Scenario = factor(data_final_grand$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_grand <- ggplot(data_final_grand) +
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
altitude_grand <- ggplot(data_final_grand) +
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

data_za <- ca_test_na3_za %>% filter(Adansonia.za_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

data_2_za <- sample_n(data_za, size = 1000, replace = F)
data_2_za$Proj <- rep(c("Present"),1000)
names(data_2_za) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caZD_za <- raster(paste0("Adansonia.za/caFut_85_2080.tif"))
ca_test_fut_za <- stack(caZD_za,environ$alt)
ca_test_na_fut_za <- na.omit(ca_test_fut_za)

ca_test_na2_fut_za <- rasterToPoints(ca_test_na_fut_za)
ca_test_na3_fut_za <- as.data.frame(ca_test_na2_fut_za[complete.cases(ca_test_na2_fut_za), ] )
data_fut_za <- ca_test_na3_fut_za %>% filter(caFut_85_2080 >= 1500)

data_2_fut_za <- sample_n(data_fut_za, size = 1000, replace = F)
data_2_fut_za$Proj <- rep(c("Future_2085"),1000)
names(data_2_fut_za) <- c("Long","Lat","Prediction","Alt","Scenario")
data_1_za <- rbind(data_2_za,data_2_fut_za)

### Future testing 2050 
caZD_za_2055 <- raster(paste0("Adansonia.za/caFut_85_2050.tif"))
ca_test_fut_za_2055 <- stack(caZD_za_2055,environ$alt)
ca_test_na_fut_za_55 <- na.omit(ca_test_fut_za_2055)

ca_test_na2_fut_za_55 <- rasterToPoints(ca_test_na_fut_za_55)
ca_test_na3_fut_za_55 <- as.data.frame(ca_test_na2_fut_za_55[complete.cases(ca_test_na2_fut_za_55), ] )
data_fut_za_55 <- ca_test_na3_fut_za_55 %>% filter(caFut_85_2050 >= 1500)

data_2_fut_za_55 <- sample_n(data_fut_za_55, size = 1000, replace = F)
data_2_fut_za_55$Proj <- rep(c("Future_2055"),1000)
names(data_2_fut_za_55) <- c("Long","Lat","Prediction","Alt","Scenario")
data_final_za <- rbind(data_2_fut_za_55,data_1_za)

data_final_za$Scenario = factor(data_final_za$Scenario,
                                         levels = c("Present", "Future_2055","Future_2085"),
                                         labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_za <- ggplot(data_final_za) +
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
altitude_za <- ggplot(data_final_za) +
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

data_mada <- ca_test_na3_mada %>% filter(Adansonia.madagascariensis_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

data_2_mada <- sample_n(data_mada, size = 1000, replace = F)
data_2_mada$Proj <- rep(c("Present"),1000)
names(data_2_mada) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caFut_mada <- raster(paste0("Adansonia.madagascariensis/caFut_85_2080.tif"))
ca_test_fut_mada <- stack(caFut_mada,environ$alt)
ca_test_na_fut_mada <- na.omit(ca_test_fut_mada)

ca_test_na2_fut_mada <- rasterToPoints(ca_test_na_fut_mada)
ca_test_na3_fut_mada <- as.data.frame(ca_test_na2_fut_mada[complete.cases(ca_test_na2_fut_mada), ] )
data_fut_mada <- ca_test_na3_fut_mada %>% filter(caFut_85_2080 >= 1500)

data_2_fut_mada <- sample_n(data_fut_mada, size = 1000, replace = F)
data_2_fut_mada$Proj <- rep(c("Future_2085"),1000)
names(data_2_fut_mada) <- c("Long","Lat","Prediction","Alt","Scenario")
data_1_mada <- rbind(data_2_mada,data_2_fut_mada)

### Future testing 2050 
caFut_mada_2055 <- raster(paste0("Adansonia.madagascariensis/caFut_85_2050.tif"))
ca_test_fut_mada_2055 <- stack(caFut_mada_2055,environ$alt)
ca_test_na_fut_mada_55 <- na.omit(ca_test_fut_mada_2055)

ca_test_na2_fut_mada_55 <- rasterToPoints(ca_test_na_fut_mada_55)
ca_test_na3_fut_mada_55 <- as.data.frame(ca_test_na2_fut_mada_55[complete.cases(ca_test_na2_fut_mada_55), ] )
data_fut_mada_55 <- ca_test_na3_fut_mada_55 %>% filter(caFut_85_2050 >= 1500)


data_2_fut_mada_55 <- sample_n(data_fut_mada_55, size = 1000, replace = F)
data_2_fut_mada_55$Proj <- rep(c("Future_2055"),1000)
names(data_2_fut_mada_55) <- c("Long","Lat","Prediction","Alt","Scenario")
data_final_mada <- rbind(data_2_fut_mada_55,data_1_mada)

data_final_mada$Scenario = factor(data_final_mada$Scenario,
                                           levels = c("Present", "Future_2055","Future_2085"),
                                           labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_mada <- ggplot(data_final_mada) +
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

altitude_mada <- ggplot(data_final_mada) +
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

data_perrieri <- ca_test_na3_perrieri %>% filter(Adansonia.perrieri_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

data_2_perrieri <- sample_n(data_perrieri, size = 1000, replace = F)
data_2_perrieri$Proj <- rep(c("Present"),1000)
names(data_2_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 (RCP 4.5!!!) Atention!
caFut_perrieri <- raster(paste0("Adansonia.perrieri/caFut_45_2080.tif"))
ca_test_fut_perrieri <- stack(caFut_perrieri,environ$alt)
ca_test_na_fut_perrieri <- na.omit(ca_test_fut_perrieri)

ca_test_na2_fut_perrieri <- rasterToPoints(ca_test_na_fut_perrieri)
ca_test_na3_fut_perrieri <- as.data.frame(ca_test_na2_fut_perrieri[complete.cases(ca_test_na2_fut_perrieri), ] )
data_fut_perrieri <- ca_test_na3_fut_perrieri %>% filter(caFut_45_2080 >= 1500)

data_2_fut_perrieri <- sample_n(data_fut_perrieri, size = 416, replace = F)
data_2_fut_perrieri$Proj <- rep(c("Future_2085"),416)
names(data_2_fut_perrieri) <- c("Long","Lat","Prediction","Alt","Scenario")
data_1_perrieri <- rbind(data_2_perrieri,data_2_fut_perrieri)

### Future testing 2050 
caFut_perrieri_2055 <- raster(paste0("Adansonia.perrieri/caFut_85_2050.tif"))
ca_test_fut_perrieri_2055 <- stack(caFut_perrieri_2055,environ$alt)
ca_test_na_fut_perrieri_55 <- na.omit(ca_test_fut_perrieri_2055)

ca_test_na2_fut_perrieri_55 <- rasterToPoints(ca_test_na_fut_perrieri_55)
ca_test_na3_fut_perrieri_55 <- as.data.frame(ca_test_na2_fut_perrieri_55[complete.cases(ca_test_na2_fut_perrieri_55), ] )
data_fut_perrieri_55 <- ca_test_na3_fut_perrieri_55 %>% filter(caFut_85_2050 >= 1500)


data_2_fut_perrieri_55 <- sample_n(data_fut_perrieri_55, size = 1000, replace = F)
data_2_fut_perrieri_55$Proj <- rep(c("Future_2055"),1000)
names(data_2_fut_perrieri_55) <- c("Long","Lat","Prediction","Alt","Scenario")
data_final_perrieri <- rbind(data_2_fut_perrieri_55,data_1_perrieri)

data_final_perrieri$Scenario = factor(data_final_perrieri$Scenario,
                                               levels = c("Present", "Future_2055","Future_2085"),
                                               labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_perrieri <- ggplot(data_final_perrieri) +
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
altitude_perrieri <- ggplot(data_final_perrieri) +
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

data_rubro <- ca_test_na3_rubro %>% filter(Adansonia.rubrostipa_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

data_2_rubro <- sample_n(data_rubro, size = 1000, replace = F)
data_2_rubro$Proj <- rep(c("Present"),1000)
names(data_2_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080
caFut_rubro <- raster(paste0("Adansonia.rubrostipa/caFut_85_2080.tif"))
ca_test_fut_rubro <- stack(caFut_rubro,environ$alt)
ca_test_na_fut_rubro <- na.omit(ca_test_fut_rubro)

ca_test_na2_fut_rubro <- rasterToPoints(ca_test_na_fut_rubro)
ca_test_na3_fut_rubro <- as.data.frame(ca_test_na2_fut_rubro[complete.cases(ca_test_na2_fut_rubro), ] )
data_fut_rubro <- ca_test_na3_fut_rubro %>% filter(caFut_85_2080 >= 1500)

data_2_fut_rubro <- sample_n(data_fut_rubro, size = 1000, replace = F)
data_2_fut_rubro$Proj <- rep(c("Future_2085"),1000)
names(data_2_fut_rubro) <- c("Long","Lat","Prediction","Alt","Scenario")
data_1_rubro <- rbind(data_2_fut_rubro,data_2_rubro)

### Future testing 2050 
caFut_rubro_2055 <- raster(paste0("Adansonia.rubrostipa/caFut_85_2050.tif"))
ca_test_fut_rubro_2055 <- stack(caFut_rubro_2055,environ$alt)
ca_test_na_fut_rubro_55 <- na.omit(ca_test_fut_rubro_2055)

ca_test_na2_fut_rubro_55 <- rasterToPoints(ca_test_na_fut_rubro_55)
ca_test_na3_fut_rubro_55 <- as.data.frame(ca_test_na2_fut_rubro_55[complete.cases(ca_test_na2_fut_rubro_55), ] )
data_fut_rubro_55 <- ca_test_na3_fut_rubro_55 %>% filter(caFut_85_2050 >= 1500)


data_2_fut_rubro_55 <- sample_n(data_fut_rubro_55, size = 1000, replace = F)
data_2_fut_rubro_55$Proj <- rep(c("Future_2055"),1000)
names(data_2_fut_rubro_55) <- c("Long","Lat","Prediction","Alt","Scenario")
data_final_rubro <- rbind(data_2_fut_rubro_55,data_1_rubro)

data_final_rubro$Scenario = factor(data_final_rubro$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))

### Latitude plot
latitude_rubro <- ggplot(data_final_rubro) +
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
altitude_rubro <- ggplot(data_final_rubro) +
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

data_suare <- ca_test_na3_suare %>% filter(Adansonia.suarezensis_EMcaByTSS_mergedAlgo_mergedRun_mergedData >= 500)

data_2_suare <- sample_n(data_suare, size = 1000, replace = F)
data_2_suare$Proj <- rep(c("Present"),1000)
names(data_2_suare) <- c("Long","Lat","Prediction","Alt","Scenario")

### Future testing 2080 Atention RCP 4.5 (SAME AS A. perrieri)
caFut_suare <- raster(paste0("Adansonia.suarezensis/caFut_45_2080.tif"))
ca_test_fut_suare <- stack(caFut_suare,environ$alt)
ca_test_na_fut_suare <- na.omit(ca_test_fut_suare)

ca_test_na2_fut_suare <- rasterToPoints(ca_test_na_fut_suare)
ca_test_na3_fut_suare <- as.data.frame(ca_test_na2_fut_suare[complete.cases(ca_test_na2_fut_suare), ] )
data_fut_suare <- ca_test_na3_fut_suare %>% filter(caFut_45_2080 >= 1500)

data_2_fut_suare <- sample_n(data_fut_suare, size = 105, replace = F)
data_2_fut_suare$Proj <- rep(c("Future_2085"),105)
names(data_2_fut_suare) <- c("Long","Lat","Prediction","Alt","Scenario")
data_test_suare <- rbind(data_2_suare,data_2_fut_suare)

### Future testing 2050 
caFut_suare_2055 <- raster(paste0("Adansonia.suarezensis/caFut_45_2050.tif"))
ca_test_fut_suare_2055 <- stack(caFut_suare_2055,environ$alt)
ca_test_na_fut_suare_55 <- na.omit(ca_test_fut_suare_2055)

ca_test_na2_fut_suare_55 <- rasterToPoints(ca_test_na_fut_suare_55)
ca_test_na3_fut_suare_55 <- as.data.frame(ca_test_na2_fut_suare_55[complete.cases(ca_test_na2_fut_suare_55), ] )
data_fut_suare_55 <- ca_test_na3_fut_suare_55 %>% filter(caFut_45_2050 >= 1500)

data_2_fut_suare_55 <- sample_n(data_fut_suare_55, size = 15, replace = F)
data_2_fut_suare_55$Proj <- rep(c("Future_2055"),15)
names(data_2_fut_suare_55) <- c("Long","Lat","Prediction","Alt","Scenario")
data_final_suare <- rbind(alegria2_fut_suare_55,data_test_suare)

data_final_suare$Scenario = factor(data_final_suare$Scenario,
                                            levels = c("Present", "Future_2055","Future_2085"),
                                            labels = c("Present", "Future 2055","Future 2085"))
### Latitude plot
latitude_suare <- ggplot(data_final_suare) +
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
altitude_suare <- ggplot(data_final_suare) +
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
       plot=a5,width=18,height=15,dpi=300)

ggsave(file=paste0("./outputs/lat_ele_all_species_within_SDAcf_non_threat.png"),
       plot=a6,width=18,height=15,dpi=300)

# ===========
# End of file
# ===========
