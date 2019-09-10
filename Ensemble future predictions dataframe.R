library(cmsaf)
library(ncdf4)
library(fields)
library(colorRamps)

my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#colors to use in graphs

#-----------
#LOAD DATA

setwd("C:/Users/Murilo/Desktop/Future predictions CMIP5/RCP 4.5")

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/ACESS1.3 RCP4.5/")
Acess_2016_2025=nc_open("tos_day_ACCESS1-3_rcp45_r1i1p1_20160101-20251231.nc")
print(Acess_2016_2025)

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/MRI-CGCM3 RCP4.5/")
mri_2016_2025=nc_open("tos_day_MRI-CGCM3_rcp45_r1i1p1_20160101-20251231.nc")
print(mri_2016_2025)

setwd("/Volumes/GoogleDrive/Shared Drives/TrEnCh/Projects/Murilo/SST future predictions/CNRM-CM5 RCP4.5/")
cnrm_2016_2025=nc_open("tos_day_CNRM-CM5_rcp45_r1i1p1_20160101-20251231.nc")
print(cnrm_2016_2025)

#------------
#Acess

#extract variables tos,lon,lat and time
sea_surface=ncvar_get(Acess_2016_2025,"tos")#get infor bout SST
print(sea_surface)
dim(sea_surface)
sea_surface = sea_surface - 273.15# convert units from kelvin to celcius
lons= ncvar_get(Acess_2016_2025,"lon") #get info about long
lons= lons[,1]
dim(lons)
lats= ncvar_get(Acess_2016_2025,"lat") #get info about latitude
dim(lats)
lats= lats[1,]
times= ncvar_get(Acess_2016_2025,"time") #get info about time
dim(times)
times1= as.Date(times, origin="0001-01-01")#convert to Julian date, since 0001-01-01
times1

#estimate years and doys
years=as.numeric(format(times1,"%Y"))
doys= as.numeric(format(times1,"%j"))

#convert longituide
convert.lon= function(r0) ifelse(r0 > 180, -360 + r0, r0)
lons= convert.lon(lons)

#subset data to area of interest ##FIX
lats.sel= which(lats> (-30) &lats<(20))#specific area of interest in atlantic ocean
lons.sel= which(lons>(-60)&lons<(-25))#specific area of interest in atlantic ocean
surface_temp=sea_surface[lons.sel,lats.sel,] #using the specified area
dim(surface_temp)

#---
#for each latitude, delete values before NA
for(lat.id in 1:50){
  na1=which.max(is.na(surface_temp[,lat.id,]))
  surface_temp[1:na1,lat.id,]<-NA
}

# selecting days from January 1st to March 31th and October 1st to December 31th of yeach year
times2=c(1:91,275:457,640:821,1005:1186,1370:1553,1737:1917,2101:2282,2466:2647,2831:3013,3197:3378,3562:3653)
#select seasons, start and end dates
#inds=cbind(c(275, 640, 1005, 1370, 1737, 2101, 2466, 2831, 3197, 3562),c(457,821,1186,1553,1917,2282,2647,3013,3378,3653))
inds=cbind(c(274, 639, 1004, 1370, 1735, 2100, 2465, 2831, 3196),c(455,820,1186,1551,1916,2281,2647,3012,3377))
years.inds= years[times2]
doys.inds= doys[times2]

#sst= surface_temp[,,1:91]

times3=c(275:366,1:91)

#update latitudes and longitudes
lons= lons[lons.sel]
lats= lats[lats.sel]

#plot
image.plot(lons,lats,surface_temp[,,1700],col=my.colors(1000))

#----------
#MRI

#extract variables tos,lon,lat and time
sea_surface.mri=ncvar_get(mri_2016_2025,"tos")#get infor bout SST
print(sea_surface.mri)
dim(sea_surface.mri)
sea_surface.mri = sea_surface.mri - 273.15# convert units from kelvin to celcius
lons.mri= ncvar_get(mri_2016_2025,"rlat") #get info about long
dim(lons.mri)
lats.mri= ncvar_get(mri_2016_2025,"rlon") #get info about latitude
dim(lats.mri)
times.mri= ncvar_get(mri_2016_2025,"time") #get info about time
dim(times.mri)
times1.mri= as.Date(times.mri, origin="1850-01-01")#convert to Julian date, since 0001-01-01
times1.mri
### NEED TO ROTATE GRID

#estimate years and doys
years.mri=as.numeric(format(times1.mri,"%Y"))
doys.mri= as.numeric(format(times1.mri,"%j"))

#----------
#CNRM

#extract variables tos,lon,lat and time
sea_surface.cnrm=ncvar_get(cnrm_2016_2025,"tos")#get infor bout SST
print(sea_surface.cnrm)
dim(sea_surface.cnrm)
sea_surface.cnrm = sea_surface.cnrm - 273.15# convert units from kelvin to celcius
lons.cnrm= ncvar_get(cnrm_2016_2025,"lon") #get info about long
dim(lons.cnrm)
lons.cnrm= lons.cnrm[,1]
lats.cnrm= ncvar_get(cnrm_2016_2025,"lat") #get info about latitude
dim(lats.cnrm)
lats.cnrm= lats.cnrm[1,]
times.cnrm= ncvar_get(cnrm_2016_2025,"time") #get info about time
dim(times.cnrm)
times1.cnrm= as.Date(times.cnrm, origin="2006-01-01")#convert to Julian date, since 0001-01-01
times1.cnrm

#convert longituide
lons.cnrm= convert.lon(lons.cnrm)

#estimate years and doys
years.cnrm=as.numeric(format(times1.cnrm,"%Y"))
doys.cnrm= as.numeric(format(times1.cnrm,"%j"))

#----------
#FIND CNRM VALUES THAT CORRESPOND TO ACESS GRID
sst= surface_temp

#make arrays for mri output
sst.cnrm= array(NA, dim= dim(surface_temp))

#match times
match.times= match(times1, times1.cnrm)

#matching function
findit = function(x,vec){ 
y = abs(vec - x) 
if(all(is.na(y)))NA else which.min(y) } 

#find closest lats and lons
match.lats= sapply(lats,findit,lats.cnrm) 
match.lons= sapply(lons,findit,lons.cnrm) 

sst.cnrm= sea_surface.cnrm[match.lons, match.lats, match.times]

#plot
image.plot(lons,lats,sst.cnrm[,,1700],col=my.colors(1000))


