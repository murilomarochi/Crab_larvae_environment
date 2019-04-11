library(raster)
library(ncdf4)
library(colorRamps)
library(sdmpredictors)
library(leaflet)
library(csv)

atlantic.ext <- extent(-104,-28.58, -39.31,30.69)# atlantic crop for the area I will use - long lat
atlantic.ext2 <- extent(256,331.52, -39.31,30.69)#atlantic crop for the area I will use - geographic coordinates
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#color to use in graphs

###sst.daily.mean, dissolved oxygen, ph, salinity mean and range, current velocity mean and range
####phytoplankton mean and range### data from Bio-Oracle
setwd ("C:/Users/Murilo/Desktop/Data/Download data_set")
list_datasets() # Explore datasets in the package 
list_layers()
list_layers(marine=TRUE)
list_code
test<-list_layers(marine=TRUE)
write.csv(test, file = "test.csv",row.names=FALSE)

setwd ("C:/Users/Murilo/Desktop/Data/Abiotic_data_set/NOAA/Sea surface daily mean temp/convert to tif")
dir()
sst_daily_1981_2018 <- dir(pattern = ".tif")
sst_daily_1981_2018
sst_daily_1981_2018.raster<-raster::stack(sst_daily_1981_2018)
extent(sst_daily_1981_2018.raster)
sst_daily_1981_2018.crop <- crop(sst_daily_1981_2018.raster, atlantic.ext2)
plot(sst_daily_1981_2018.crop,col=my.colors(1000),axes=TRUE, box=TRUE) 
title(cex.sub = 1.25, sub = "sst.mean")

oxygen<-load_layers("BO2_dissoxmean_ss",datadir="C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
oxygen
oxygen.crop <- crop(oxygen, atlantic.ext)
plot(oxygen.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "Dissolved Oxygen")

ph<-load_layers("BO_ph",datadir="C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
ph
ph.crop <- crop(ph, atlantic.ext)
plot(ph.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "ph")

salinity_mean<-load_layers("BO2_salinitymean_ss")
salinity_mean
salinity_mean.crop <- crop(salinity_mean, atlantic.ext)
plot(salinity_mean.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "salinity_mean")

salinity_range<-load_layers("BO2_salinityrange_ss",datadir="C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
salinity_range.crop <- crop(salinity_range, atlantic.ext)
plot(salinity_range.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "salinity_range")                   

current_velocity_mean<-load_layers("BO2_curvelmean_ss",datadir="C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
current_velocity_mean.crop <- crop(current_velocity_mean, atlantic.ext)
plot(current_velocity_mean.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "current_velocity_mean")   

current_velocity_range<-load_layers("BO2_curvelrange_ss",datadir="C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
current_velocity_range.crop <- crop(current_velocity_range, atlantic.ext)
plot(current_velocity_range.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "current_velocity_range")

phyto_mean<-load_layers("BO2_carbonphytomean_ss", datadir="C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
phyto_mean.crop <- crop(phyto_mean, atlantic.ext)
plot(phyto_mean.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "phyto_mean")

phyto_range<-load_layers("BO2_carbonphytorange_ss", datadir="C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
phyto_range.crop <- crop(phyto_range, atlantic.ext)
plot(phyto_range.crop,col=my.colors(1000),axes=FALSE, box=FALSE) 
title(cex.sub = 1.25, sub = "phyto_range")


###all environment variables together, except for stt
setwd("C:/Users/Murilo/Desktop/Data/Abiotic_data_set/Bio_Oracle")
envi<- dir(pattern = ".tif")
envi
environment.variables<-raster::stack(envi)
extent(environment.variables)
environment.variables_crop <- crop(environment.variables, atlantic.ext)
plot(environment.variables_crop[[1:4]], col = my.colors(1000),axes=FALSE, box=FALSE)
