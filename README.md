# Baobabs Mada

This GitHub repository includes the code and data for the following scientific article:

**MUNIZ TAGLIARI Mario, Pascal DANTHU, Jean-Michel LEONG POCK TSY, Cyrille CORNU, Jonathan LENOIR, Vítor CARVALHO-ROCHA, Ghislain VIEILLEDENT.** 2021. Not all species will migrate poleward as the climate warms: the case of the seven baobab species in Madagascar.

## Baobab occurrence data

The `data/baobabs/` folder includes the occurrence dataset for the seven baobab species present in Madagascar (file `data_Adansonia.csv`). This file is the result of years of field inventories in Madagascar by botanists and ecologists working at [CIRAD](https://www.cirad.fr).

## R script

All the results of the study can be reproduced executing the R script `R/baoabas.R`. To do so, you can either execute the `baobabs.sh` shellscript or run the R script using a R GUI such as RStudio.

## Results

<<<<<<< HEAD
All the results of the study (tables and figures) are saved in the `outputs` folder.
=======
All the results are saved in the `outputs` folder.

<img alt="Bioclimatic niche of studied species" src="outputs/Fig_ap_sps_niche.png" width="500">

Figure 1: **Bioclimatic niche of the seven baobabs species found in Madagascar**.

Check the most threatened baobab species.

<img alt="Threatened baobab species" src="outputs/plot_SDA_threat.png" width="500">

Figure 2: **The most threatened Malagasy baobab species (2085) under two dispersal hypothesis (Full-Dispersal and Zero-Dispersal) and RCP 8.5**.

Also check how latitude and elevation will change inside threatened baobab projected distribution in the present and in the future (2055 and 2085).

<img alt="Latitudinal and elevation change of  threatened baobab species" src="outputs/lat_ele_all_species_within_SDAcf_threat.png" width='500'>
>>>>>>> 8402e913a3de0d5e6f4638c2a7170aed1b2ef69a

<img alt="SDA" src="outputs/plot_SDA_threat.png" width="500">

Figure 1: **Species range contraction under climate change for the four threatened baobab species.** The four species are _A. madagascariensis_, _A. perrieri_, _A. rubrostipa_, and _A. suarezensis_ (one species per row). (a, e, i, m) Occurrence points over Madagascar elevation map (elevation in m); (b, f, j, n) Current predicted species distribution. Legend indicates the number of models (0-4) predicting the species presence; (c,g,k,o) Projected species distribution in 2085 under scenario RCP 8.5 and the full-dispersal hypothesis. Legend indicates the number of models (0-12) predicting the species presence; (d,h,l,p) Projected species distribution in 2085 under scenario RCP 8.5 and the zero-dispersal hypothesis. Legend indicates the number of models (0-12) predicting the species presence. For the distribution maps, the species is assumed to be present (green areas) when a majority of models predicts a presence (votes >= 2 in the present, and >= 6 in the future). The species is considered absent (grey areas) when no model (votes = 0), or a minority of models (votes < 2 in the present, and < 6 in the future), predicts a presence. Maps for _A. perrieri_ and _A. suarezensis_, two species distributed at the extreme North of Madagascar, have been zoomed in (black squares).

## Archive

An archive of this repository is available on the Cirad Dataverse: \[DOI: [10.18167/DVN1/KMIO8N](http://dx.doi.org/10.18167/DVN1/KMIO8N)\].
