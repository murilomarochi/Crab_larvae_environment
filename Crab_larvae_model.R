#===================LIBRARIES===================
library(raster)
library(ncdf4)
library(fields)
library(colorRamps)
library(sdmpredictors)
library(leaflet)
library(abind)
library(csv)

#===============================================
##################atlantic crop area and nice colors for maps
atlantic.ext <- extent(-104,-28.58, -39.31,30.69)# atlantic crop for the area I will use - long lat
atlantic.ext2 <- extent(256,331.52, -39.31,30.69)#atlantic crop for the area I will use - geographic coordinates
my.colors = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))#color to use in graphs

#--------------------
#LOAD SURVIVAL DATA

#CALCULATE SURVIVAL FUNCTION
#something like: surv= lm( poly(time,3)*temp)
#extract coefficients from final model to make a survival function
#survival= function(time, temp){ "Put model here" }

##################load SST 

setwd ("C:/Users/Murilo/Desktop/Data/Abiotic_data_set/NOAA/Sea surface daily mean temp")
setwd("/Volumes/GoogleDrive/Team Drives/TrEnCh/Projects/Murilo/Sea_surface_daily_temp/Sea surface daily mean temp")

dir()

#years to analyze, correct years?
years= 1982:2018

#construct file names for first year, Can eventually loop through years
year.ind=1
file.name= paste("sst.day.mean.", years[year.ind],".nc", sep="")

#load data for 1 year
sst<-nc_open(file.name)
print(sst)

#extract temperature data
temp=ncvar_get(sst,"sst")
dim(temp) #dimensions lon, lat, time

#extract lon, lat, time
lons= ncvar_get(sst,"lon") #get info about long
dim(lons)
lats= ncvar_get(sst,"lat") #get info about latitude
dim(lats)
times= ncvar_get(sst,"time") #julian date, calendar day, since 1800-01-01
times
#change to day of year
times1= as.Date(times, origin="1800-01-01")
times1
dim(times1)

#--------------
#subset data to area of interest
lats.sel= which(lats>(-39.31)&lats<(30.69))#specific area of interest in atlantic ocean
lons.sel= which(lons>(260)&lons<(331.52))#specific area of interest in atlantic ocean

temp.sub=temp[lons.sel,lats.sel,] #using the specified area
dim(temp.sub)
#plot out
image.plot(lons[lons.sel],lats[lats.sel],temp[lons.sel,lats.sel,365],col=my.colors(1000))

#get rid of Pacific values
temp.sub2<-temp.sub #lon, lat, time
#for each latitude, delete values before NA
for(lat.id in 1:280){
  na1=which.max(is.na(temp.sub2[,lat.id,1]))
  temp.sub2[1:na1,lat.id,]<-NA
  print(lat.id)
}
temp.sub2
dim(temp.sub2)

#update latitudes and longitudes
lons= lons[lons.sel]
lats= lats[lats.sel]
#update temp data to subset
temp=temp.sub2

#image.plot format longitude, latitude, matrix of data
image.plot(lons,lats,temp[,,1],col=my.colors(1000))
#--------------

#Set species' parameters
To=17.2 ###lower developmental temperature
G=204.08 ### daily growing degree days 

#define seasonal timing in terms of days of year
times2=1:90 # selecting days from January 1st to March 31th
times2.1=257:365 #selecting days from Setempber to December 31th    

#make arrays for output
devel.out= array(NA, dim=c(length(lons),length(lats),length(times2)))

for(lon.k in 1:length(lons)){ #loop through longitude
    
  #find latitudes with data
    lats.sel2= which(!is.na(temp[lon.k,,1]))
    
    if(length(lats.sel2)>0){ #check there's cells with data
    for(lat.k in lats.sel2){ #loop through latitudes with data
    
      #estimate growing degree days
      GDDs= temp[lon.k,lat.k,]-To #how many growing degree days per day
    
      #cummulative sum of GDDs
      cumGDDs=cumsum(GDDs) 
    
      #estimate development time starting at different days
      for(time.k in times2){
      
        cumGDDs= cumsum( GDDs[times2[time.k]:365] )
        #substract off starting value
        cumGDDs= cumGDDs-cumGDDs[1]
        
        #find first date the exceed G
        devel.out[lon.k,lat.k,time.k]= which.max(cumGDDs>G)+time.k
    
        #calculate survival
        #average temperature across development period
        #use survival function calculated above to predict suvival as a function of duration of development and average temperature
        
         } #end loop timing
      
  }#end latitude loop
} #end check for data
} #end longitude loop

devel.out[lon.k,lat.k,]

#plot out doy development assuming development starts on doy 50
image.plot(lons,lats,devel.out[,,50],col=my.colors(1000))



###########I try to implement a specific period of time and a date to larvae growth on grid cells, but did not work#######
######## I tried to inform the specific dates that the larvae goes to the ocean, but I still do not know how inform a specific grid cell.
times2=times1[1:90] # selecting days from January 1st to March 31th
times2.1=times1[257:365] #selecting days from Setempber to December 31th
dates.k=temp[lons.sel,lats.sel,1:90]# selecting days from January 1st to March 31th in the grid
time.develop=25 #days to reach megalopae stage

out2= array(NA, dim=c(286,280,365))
out2
for(lon.k in 1:length(lons.sel)){#dimension in 286
  #find latitudes with data
  lats.sel2= which(!is.na(dates.k[lon.k,,]))
  if(length(lats.sel2)>0){
    for(lat.k in lats.sel2){ #dimension in 720
      GDDs2= dates.k[lon.k,lat.k,]-To #how many growing degree days per day
      cumGDDs2=cumsum(GDDs) #cummulative sum of GDDs
            out2[lon.k,lat.k,time.k]=which.max(cumGDDs>G)#find first date the exceed G
    }#end latitude loop
  } #end check for data
} #end longitude loop
GDDs2
cumGDDs2
out2

######## Attempt number 2
out3= array(NA, dim=c(286,280,365))
out3
for(lon.k in 1:length(lons.sel)){ #dimension in 286
    lats.sel2= which(!is.na(temp.sub2[lon.k,,time.k])) #find latitudes with data
  if(length(lats.sel2)>0){
    for(lat.k in lats.sel2){ #dimension in 720
      GDDs3= temp.sub2[lon.k,lat.k,time.k]-To #how many growing degree days per day
      cumGDDs3=cumsum(GDDs) #cummulative sum of GDDs
       for(dates in times1){ #dimension in 365
        times2=which.max(90(temp.sub2[,,dates]))
         out3[lon.k,lat.k,times2]=which.max(cumGDDs>G)#find first date the exceed G
     }#end latitude loop
   } #end check for data
 } #end longitude loop
}
GDDs3
cumGDDs3
out3

