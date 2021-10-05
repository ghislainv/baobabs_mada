# =============================================
# Tseas unit in WorlClim v1.4 or 2.1
# =============================================

# Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>

# Directory for outputs
dir.create("wc_comp")

# Version 1.4
zipfile="wc_comp/bio_10m_bil.zip"
download.file("https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip", destfile=zipfile, method = "wget")
unzip(zipfile, files=c("bio4.bil", "bio4.hdr"), exdir="wc_comp")
tseas_wc_1_4 <- raster("wc_comp/bio4.bil")
pdf("wc_comp/tseas_wc_1_4.pdf")
plot(tseas_wc_1_4, main="tseas WorldClim 1.4") # goes from 0 to 20000 (**twenty** thousand)
dev.off()

# Version 2.1
zipfile="wc_comp/wc2.1_10m_bio.zip"
download.file("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip", destfile=zipfile, method = "wget")
unzip(zipfile, files=c("wc2.1_10m_bio_4.tif"), exdir="wc_comp")
tseas_wc_2_1 <- raster("wc_comp/wc2.1_10m_bio_4.tif")
pdf("wc_comp/tseas_wc_2_1.pdf")
plot(tseas_wc_2_1, main="tseas WorldClim 2.1") # goes from 0 to 2000 (**two** thousand)
dev.off()
# crop to common extent
tseas_wc_2_1_crop <- raster::crop(tseas_wc_2_1, tseas_wc_1_4)

# Sample points and compare
tseas <- stack(tseas_wc_1_4, tseas_wc_2_1_crop)
set.seed(1234)
samp <- as.data.frame(sampleRandom(tseas, 1000))
names(samp) <- c("wc_1_4", "wc_2_1")
# Simple linear model to find correlation
mod <- lm(wc_1_4~ -1 + wc_2_1, data=samp)
print(mod)
# Call:
# lm(formula = wc_1_4 ~ -1 + wc_2_0, data = samp)
# 
# Coefficients:
# wc_2_0  
#  9.706
pdf("wc_comp/comp.pdf")
plot(samp$wc_2_1, samp$wc_1_4,
     xlab="bio4 (tseas) WorldClim 2.1",
     ylab="bio4 (tseas) WorldClim 1.4")
abline(a=0, b=1, col="red")
dev.off()

# End of file