#!/usr/bin/Rscript

# ==============================================================================
# authors         :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# license         :GPLv3
# ==============================================================================

# Libraries
require(readr)
require(dplyr)

# ==============================================================================
# Import and format each dataset
# ==============================================================================

# Set URL for raw datasets

# Data-set from Jean-Michel Leong Pock Tsy sent on 20/11/2012
in1_1 <- read_csv("data/baobabs_raw_data/1_1_data_Adansonia_2004_2012.csv")
cor1_2 <- read_csv("data/baobabs_raw_data/1_2_data_corrections_A_suar_mada_confusion.csv") %>%
  dplyr::mutate(TreeID=paste(Species, TreeName, Lat, Long, sep="_"))
out1 <- in1_1 %>%
  dplyr::filter(!is.na(Lat) & !is.na(Long)) %>%
  dplyr::mutate(Year="2004-2015", Institute="CIRAD", Collector="Leong Pock Tsy Jean-Michel; Danthu Pascal",
         Source="field observation", Dataset="dataset1") %>%
  # Correct two points which are far from the coast line
  dplyr::mutate(Lat=ifelse(paste(Species, TreeName, sep="_") == "A_grandidieri_AdBefa_97", -22.069, Lat)) %>%
  dplyr::mutate(Lat=ifelse(paste(Species, TreeName, sep="_") == "A_suarezensis_ANDAC8", -12.3283833333333, Lat)) %>%
  # Correct misidentified A_suarezensis trees
  dplyr::mutate(TreeID=paste(Species, TreeName, Lat, Long, sep="_")) %>%
  dplyr::mutate(Species=ifelse(TreeID %in% cor1_2$TreeID, "A_madagascariensis", Species)) %>%
  # Select columns
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# Adansonia suarezensis data from Cyrille Cornu field observations in 2014
in2 <- read_csv("data/baobabs_raw_data/2_data_A_suarezensis_cc_2014.csv")
out2 <- in2 %>%
  dplyr::rename(Lat=Y, Long=X) %>%
  dplyr::mutate(Species="A_suarezensis", Year="2014", Institute="CIRAD",
         Source="field observation",
         Collector="Cornu Cyrille",
         TreeName=paste0("ADSU2010CCfield_",dplyr::row_number()),
         Dataset="dataset2") %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# Adansonia grandidieri data from Cyrille Cornu photo-interpretation in 2009
in3 <- read_csv("data/baobabs_raw_data/3_data_A_grandidieri_sat_2009.csv")
out3 <- in3 %>%
  dplyr::rename(Lat=Y, Long=X) %>%
  dplyr::mutate(Species="A_grandidieri", Year="2009", Institute="CIRAD",
         Source="photo-interpretation of satellite images",
         Dataset="dataset3",
         Collector="Cornu Cyrille", TreeName=paste0("ADGR2009CCsat",dplyr::row_number())) %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# Adansonia suarezensis data from Cyrille Cornu photo-interpretation in 2010
in4 <- read_csv("data/baobabs_raw_data/4_data_A_suarezensis_sat_2010.csv")
out4 <- in4 %>%
  dplyr::rename(Lat=Y, Long=X) %>%
  dplyr::mutate(Species="A_suarezensis", Year="2010", Institute="CIRAD",
         Source="photo-interpretation of satellite images",
         Dataset="dataset4",
         Collector="Cornu Cyrille", TreeName=paste0("ADSU2010CCsat_",dplyr::row_number())) %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# Comores by Jean-Michel Leong Pock Tsy
in5 <- read_csv("data/baobabs_raw_data/5_data_A_digitata_comoro_islands_jm_2011_2014.csv")
out5 <- in5 %>%
  dplyr::mutate(Institute="CIRAD",
         Year=as.character(Year),
         Source=ifelse(grepl("GoogleEarth",TreeName),
                       "photo-interpretation of satellite images",
                       "field observation"),
         Collector="Leong Pock Tsy Jean-Michel") %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source)

# Missions by Wilfried Ramahafaly
# Data sent by Jean-Michel on 09/05/2015
# Attention: In this data-set, coordinates are given in Lat/Long hddd°mm.mmm'
# Function to convert hdd°mm.mmm' in decimal degrees
ddmm2latlong <- function (x, coord="Lat") {
  format_x <- format(x, nsmall=5, drop0trailing=FALSE)
  split <- unlist(strsplit(format_x, "\\."))
  mat <- matrix(split, ncol=2, byrow=TRUE)
  d <- as.numeric(mat[,1])
  dd <- as.numeric(mat[,2])/60000
  s <- d+dd
  if (coord=="Lat") {y <- (-s)} else {y <- s}
  return (y)
}
# Input data
in6 <- read_csv("data/baobabs_raw_data/6_data_wr_missions_2013_2015.csv")
# Output data
out6 <- in6 %>%
  # Convert ddmm to lat/long format
  dplyr::mutate(Lat_ddmm=paste0(substr(Position,2,3),".",
                         substr(Position,5,6),
                         substr(Position,8,10)),
                Long_ddmm=paste0(substr(Position,13,14),".",
                          substr(Position,16,17),
                          substr(Position,19,21))) %>%
  dplyr::mutate(Lat=ddmm2latlong(Lat_ddmm, coord="Lat"),
                Long=ddmm2latlong(Long_ddmm, coord="Long")) %>%
  # Get species name when possible
  dplyr::mutate(SpInitials=substr(Name,1,3)) %>%
  dplyr::mutate(Species=ifelse(Mission %in% c("AZA_Betioky_2013","Betioky_2013"), "A_za", NA)) %>%
  dplyr::mutate(Species=ifelse((Mission == "Betioky_2013") & grepl("RUBT",Name), "A_rubrostipa", Species)) %>%
  dplyr::mutate(Species=ifelse(Mission == "ZAPER_2013", "A_za_perrieri", Species)) %>%
  dplyr::mutate(Species=ifelse(SpInitials == "AM ", "A_madagascariensis", Species)) %>%
  dplyr::mutate(Species=ifelse(SpInitials == "AR ", "A_rubrostipa", Species)) %>%
  dplyr::mutate(Species=ifelse(SpInitials == "AS ", "A_suarezensis", Species)) %>%
  dplyr::mutate(Species=ifelse(SpInitials == "AZ ", "A_za", Species)) %>%
  dplyr::mutate(Species=ifelse(grepl("CITER", Name), "A_za", Species)) %>%
  # We remove observations for which we don't have the species name
  dplyr::filter(!is.na(Species)) %>%
  # Retrieve the year
  dplyr::mutate(Year=paste0(20,substr(Description,8,9))) %>%
  # Additional info
  dplyr::mutate(Institute="CIRAD",
         Source="field observation",
         Dataset="dataset6",
         Collector="Ramahafaly Wilfried",
         TreeName=Name) %>%
  # Column selection
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# A_perrieri points including raw gps positions in dd°mm.mmm'
in7 <- read_csv("data/baobabs_raw_data/7_data_A_perrieri_gv_raw_gps_2013.csv")
out7 <- in7 %>%
  # Convert ddmm to lat/long format
  dplyr::mutate(Lat_ddmm=paste0(substr(Position,2,3),".",
                         substr(Position,5,6),
                         substr(Position,8,10)),
                Long_ddmm=paste0(substr(Position,13,14),".",
                          substr(Position,16,17),
                          substr(Position,19,21))) %>%
  dplyr::mutate(Lat_calc=ddmm2latlong(Lat_ddmm, coord="Lat"),
                Long_calc=ddmm2latlong(Long_ddmm, coord="Long")) %>%
  dplyr::mutate(Lat=ifelse(is.na(Lat), Lat_calc, Lat),
         Long=ifelse(is.na(Long), Long_calc, Long)) %>%
  dplyr::mutate(Species="A_perrieri", Year="2004-2015", Institute="CIRAD", TreeName=Name,
         Dataset="dataset7",
         Source="field observation", Collector="Leong Pock Tsy Jean-Michel; Danthu Pascal") %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# Some additional A_suarezensis points sent by email from Cyrille Cornu on 20/04/2015
in8_1 <- read_csv("data/baobabs_raw_data/8_1_data_A_suarezensis_cc_2015.csv")
out8_1 <- in8_1 %>%
  dplyr::mutate(Lat=Y, Long=X, Species="A_suarezensis", Year="2015", Institute="CIRAD",
                Source="field observation", Collector="Cornu Cyrille",
                Dataset="dataset8_1",
                TreeName=paste0("ADSU2015CCfield_",dplyr::row_number())) %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# Some additional Adansonia points sent by email from Cyrille Cornu on 22/05/2015
in8_2 <- read_csv("data/baobabs_raw_data/8_2_data_Adansonia_cc_2015.csv")
out8_2 <- in8_2 %>%
  dplyr::mutate(Lat=Y, Long=X,
         Year="2015", Institute="CIRAD",
         Source="field observation", Collector="Cornu Cyrille",
         Dataset="dataset8_2",
         TreeName=paste0("AD2015CCfield_",TreeID)) %>%
  dplyr::mutate(Species=gsub(" ", "_", Name)) %>%
  dplyr::mutate(Species=ifelse(Name=="A mada", "A_madagascariensis", Species)) %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# A_madagascariensis points close to Montagne des Français by Ghislain Vieilledent and Mario Tagliari
in9 <- read_csv("data/baobabs_raw_data/9_data_A_madagascariensis_gv_mt_2015.csv")
out9 <- in9 %>%
  dplyr::mutate(Year="2015", Institute="CIRAD", Dataset="dataset9",
         Source="field observation", Collector="Vieilledent Ghislain; Muniz-Tagliari Mario",
         TreeName=paste0("ADMA2015GVMTfield_",dplyr::row_number())) %>%
  dplyr::select(Species, TreeName, Lat, Long, Year, Institute, Collector, Source, Dataset)

# ==============================================================================
# Group data-sets together and remove duplicates
# ==============================================================================

# Bind and create TreeID
data_bind <- 
  # Bind data-sets
  dplyr::bind_rows(out1, out2, out3, out4, out5, out6, 
            out7, out8_1, out8_2, out9) %>%
  # Replace space in TreeName and define TreeID
  dplyr::mutate(TreeName=gsub(" ", "_", TreeName),
         TreeID=paste(Species, TreeName, Lat, Long, sep="_"))

# To look at duplicates (inside or between datasets), eg. A_za
data_test <- data_bind %>% 
  dplyr::filter(duplicated(TreeID) | duplicated(TreeID, fromLast=TRUE)) %>%
  dplyr::filter(Species == "A_za") %>%
  dplyr::arrange(TreeID)

# Remove duplicates and select columns
data_Adansonia <- data_bind %>%
  dplyr::distinct(TreeID, .keep_all=TRUE) %>%
  dplyr::select(-Dataset) %>%
  dplyr::filter(Species != "A_za_perrieri")
  
# Export data-set
dir.create("data/baobabs", showWarnings=FALSE)
write_csv(data_Adansonia, "data/baobabs/data_Adansonia.csv")

# Number of observations per species
obs_per_sp <- data_Adansonia %>%
  dplyr::group_by(Species) %>%
  dplyr::summarise(n())
write_csv(obs_per_sp, "data/baobabs/obs_per_sp.csv")

# ==============================================================================
# Message
# ==============================================================================

cat("Occurrence dataset available: data/baobabs/data_Adansonia.csv\n")

# End of script