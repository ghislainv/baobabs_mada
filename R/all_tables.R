#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent / Mario Tagliari
# email           :ghislain.vieilledent@cirad.fr / mario_tagliari@hotmail.com
# license         :GPLv3
# ==============================================================================
library(data.table)
require(dplyr)
require(here)
 
##### Generating study tables

## area change
digitata <- read.table(paste0("Adansonia.digitata/sda_fut.txt"), header=T,sep="\t")
grandidieri <- read.table(paste0("Adansonia.grandidieri/sda_fut.txt"), header=T,sep="\t")
mada <- read.table(paste0("Adansonia.madagascariensis/sda_fut.txt"), header=T,sep="\t")
perrieri <- read.table(paste0("Adansonia.perrieri/sda_fut.txt"), header=T,sep="\t")
rubrostipa <- read.table(paste0("Adansonia.rubrostipa/sda_fut.txt"), header=T,sep="\t")
suare <- read.table(paste0("Adansonia.suarezensis/sda_fut.txt"), header=T,sep="\t")
za <- read.table(paste0("Adansonia.za/sda_fut.txt"), header=T,sep="\t")

species <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                 "A. perrieri", "A. rubrostipa","A. suarezensis",
                 "A. za"),each=8)
all_species <- rbind(digitata,grandidieri,mada,perrieri,rubrostipa,suare,za)
all_species <- cbind(species,all_species)

## altitudinal change
digitata_alt <- read.table(paste0("Adansonia.digitata/alt_fut.txt"), header=T,sep="\t")
grandidieri_alt <- read.table(paste0("Adansonia.grandidieri/alt_fut.txt"), header=T,sep="\t")
mada_alt <- read.table(paste0("Adansonia.madagascariensis/alt_fut.txt"), header=T,sep="\t")
perrieri_alt <- read.table(paste0("Adansonia.perrieri/alt_fut.txt"), header=T,sep="\t")
rubrostipa_alt <- read.table(paste0("Adansonia.rubrostipa/alt_fut.txt"), header=T,sep="\t")
suare_alt <- read.table(paste0("Adansonia.suarezensis/alt_fut.txt"), header=T,sep="\t")
za_alt <- read.table(paste0("Adansonia.za/alt_fut.txt"), header=T,sep="\t")

all_species_alt <- rbind(digitata_alt,grandidieri_alt,mada_alt,
                         perrieri_alt,rubrostipa_alt,suare_alt,za_alt)

all_species_cool <- cbind(all_species,all_species_alt)
tail(all_species_cool)

final_table1 <- subset(all_species_cool, select = -c(5,6,11,12,13,14,15,17,18))

colnames(final_table1)<- c("Baobab species","SDAp (km²)","RCP","Year","Dispersal",
                           "SDAf (km²)","Change SDApf (%)","Current altitude mean (m)",
                           "Future altitude mean (m)","Altitude range shift (%)")
head(final_table1)
write.table(final_table1,paste0("./outputs/table1.txt"),sep="\t")

############################
### importance variable table
digitata_vi <- read.table(paste0("Adansonia.digitata/varimp.txt"), header=T,sep="\t")
grandidieri_vi <- read.table(paste0("Adansonia.grandidieri/varimp.txt"), header=T,sep="\t")
mada_vi <- read.table(paste0("Adansonia.madagascariensis/varimp.txt"), header=T,sep="\t")
perrieri_vi <- read.table(paste0("Adansonia.perrieri/varimp.txt"), header=T,sep="\t")
rubrostipa_vi <- read.table(paste0("Adansonia.rubrostipa/varimp.txt"), header=T,sep="\t")
suare_vi <- read.table(paste0("Adansonia.suarezensis/varimp.txt"), header=T,sep="\t")
za_vi <- read.table(paste0("Adansonia.za/varimp.txt"), header=T,sep="\t")

all_species_vi <- rbind(digitata_vi,grandidieri_vi,mada_vi,
                        perrieri_vi,rubrostipa_vi,suare_vi,za_vi)
all_species_vi$mean <- c("NA")
all_species_final <- subset(all_species_vi, select = -c(5,6))

setDT(all_species_final)
trying <- all_species_final[, .(Mean = rowMeans(.SD)), by = mean]
trying$mean <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                     "A. perrieri", "A. rubrostipa","A. suarezensis",
                     "A. za"),each=4)
trying <- as.data.frame(trying)
all_species_vi_done <- cbind(trying,all_species_final)
all_species_vi_done <- subset(all_species_vi_done, select = -c(7))

all_species_vi_done$clim_var <- rep(c("tmean","tseas","prec","cwd"),7)
all_species_vi_done <- all_species_vi_done[, c(1, 7, 3, 4, 5, 6, 2)]
colnames(all_species_vi_done)<- c("Baobab_species","Clim_var",
                                  "GLM","GAM","RF","Maxent","Mean")
try <- all_species_vi_done %>% group_by(Baobab_species) %>%
  mutate(Rank = row_number(max(Mean) - Mean))
View(try)

write.table(try,paste0("./outputs/table2vi.txt"),sep="\t")

#### Table 3 with performance indexes

suare_dig <- read.table(paste0("Adansonia.digitata/performance_ca.txt"), header=T,sep="\t")
suare_grandi <- read.table(paste0("Adansonia.grandidieri/performance_ca.txt"), header=T,sep="\t")
suare_mada <- read.table(paste0("Adansonia.madagascariensis/performance_ca.txt"), header=T,sep="\t")
suare_perri <- read.table(paste0("Adansonia.perrieri/performance_ca.txt"), header=T,sep="\t")
suare_rubro <- read.table(paste0("Adansonia.rubrostipa/performance_ca.txt"), header=T,sep="\t")
suare_suare <- read.table(paste0("Adansonia.suarezensis/performance_ca.txt"), header=T,sep="\t")
suare_za <- read.table(paste0("Adansonia.za/performance_ca.txt"), header=T,sep="\t")
ca_species <- rbind(suare_dig,suare_grandi,suare_mada,suare_perri,suare_rubro,
                    suare_suare,suare_za)

rownames(ca_species) <- c("A. digitata","A. grandidieri","A. madagascariensis", 
                          "A. perrieri", "A. rubrostipa","A. suarezensis","A. za")

only_sen_spe_tss <- subset(ca_species, select=c("Sen", "Spe","TSS"))

write.table(only_sen_spe_tss,paste0("./outputs/table3_performance_ready.txt"),sep="\t")


####### #### Table A2 with performance indexes

perf_dig <- read.table(paste0("Adansonia.digitata/current_model_evaluation.txt"), header=T,sep="\t")
perf_grand <- read.table(paste0("Adansonia.grandidieri/current_model_evaluation.txt"), header=T,sep="\t")
perf_mada <- read.table(paste0("Adansonia.madagascariensis/current_model_evaluation.txt"), header=T,sep="\t")
perf_perri <- read.table(paste0("Adansonia.perrieri/current_model_evaluation.txt"), header=T,sep="\t")
perf_rubro <- read.table(paste0("Adansonia.rubrostipa/current_model_evaluation.txt"), header=T,sep="\t")
perf_suare <- read.table(paste0("Adansonia.suarezensis/current_model_evaluation.txt"), header=T,sep="\t")
perf_za <- read.table(paste0("Adansonia.za/current_model_evaluation.txt"), header=T,sep="\t")

perf_species <- rbind(perf_dig,perf_grand,perf_mada,perf_perri,perf_rubro,
                      perf_suare,perf_za)

# using subset function
newdata <- subset(perf_species, Index == "Testing.data",
                  select=c(wIndex,Model, Run, Value))


newdatap <- newdata[newdata$Run != 'Full', ]
newdatap$species <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                          "A. perrieri", "A. rubrostipa","A. suarezensis","A. za")
                        ,each=60) 

partial_final <-  aggregate(Value~species + wIndex + Model, FUN=mean, data=newdatap)
## to calculate mean of the 5 runs

partial_final1 <- partial_final %>%
  filter(wIndex != 'KAPPA')#omit tempcol in output


write.table(partial_final1,paste0("./outputs/tableA2_performance_partial.txt"),sep="\t")

### To calculate mean over the full dataset
newdataf <- newdata[newdata$Run == 'Full', ]
newdataf$species <- rep(c("A. digitata","A. grandidieri","A. madagascariensis", 
                          "A. perrieri", "A. rubrostipa","A. suarezensis","A. za")
                        ,each=12)
full_final <-  aggregate(Value~species + wIndex + Model, FUN=mean, data=newdataf)
full_final1 <- full_final %>%
  filter(wIndex != 'KAPPA')#omit tempcol in output

partial_final1
full_final1
write.table(full_final1,paste0("./outputs/tableA2_performance_full.txt"),sep="\t")

############# Table 4

nichesf_dig <- read.table(paste0("Adansonia.digitata/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_dig <- read.table(paste0("Adansonia.digitata/niche.txt"), header=T,sep="\t")
niches_dig <- cbind(nichesf_dig,nichesc_dig)
niches_dig$species <- rep(c("A. digitata"),each=3)

nichesf_grand <- read.table(paste0("Adansonia.grandidieri/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_grand <- read.table(paste0("Adansonia.grandidieri/niche.txt"), header=T,sep="\t")
niches_grand <- cbind(nichesf_grand,nichesc_grand)
niches_grand$species <- rep(c("A. grandidieri"),each=3)

nichesf_mada <- read.table(paste0("Adansonia.madagascariensis/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_mada <- read.table(paste0("Adansonia.madagascariensis/niche.txt"), header=T,sep="\t")
niches_mada <- cbind(nichesf_mada,nichesc_mada)
niches_mada$species <- rep(c("A. madagascariensis"),each=3)

nichesf_perri <- read.table(paste0("Adansonia.perrieri/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_perri <- read.table(paste0("Adansonia.perrieri/niche.txt"), header=T,sep="\t")
niches_perri <- cbind(nichesf_perri,nichesc_perri)
niches_perri$species <- rep(c("A. perrieri"),each=3)

nichesf_rubro <- read.table(paste0("Adansonia.rubrostipa/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_rubro <- read.table(paste0("Adansonia.rubrostipa/niche.txt"), header=T,sep="\t")
niches_rubro <- cbind(nichesf_rubro,nichesc_rubro)
niches_rubro$species <- rep(c("A. rubrostipa"),each=3)

nichesf_suare <- read.table(paste0("Adansonia.suarezensis/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_suare <- read.table(paste0("Adansonia.suarezensis/niche.txt"), header=T,sep="\t")
niches_suare <- cbind(nichesf_suare,nichesc_suare)
niches_suare$species <- rep(c("A. suarezensis"),each=3)

nichesf_za <- read.table(paste0("Adansonia.za/mean_niche_with_future.txt"), header=T,sep="\t")
nichesc_za <- read.table(paste0("Adansonia.za/niche.txt"), header=T,sep="\t")
niches_za <- cbind(nichesf_za,nichesc_za)
niches_za$species <- rep(c("A. za"),each=3)

all_niches <- cbind(niches_dig,niches_grand,niches_mada,niches_perri,
                    niches_rubro,niches_suare,niches_za)

write.table(all_niches,paste0("./outputs/tableA3_niches_fut_cur_comparison.txt"),sep="\t")

# ===========
# End of file
# ===========