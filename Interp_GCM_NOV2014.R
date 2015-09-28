######################################################################################
# R program to read in GCM data 
# and interpolate to station sites 
# and take a smoothed curve from the climate anomalies

options(warn=2)
library(ncdf)
library(akima)
#library(fUtilities)	# this allows us to use linearInterpp but this gives us NA
#library(geoR)	# will now convert the GCM data into a geodata class
library(maps)
library(RColorBrewer)
source("paddedinterp.R")

######################################################################################
# NOTES



# Kate Willett July 2014
######################################################################################
# FUNCTIONS
######################################################################################
# get_climanoms_func

#C Calculate monthly climatology over entire period of record (minimum 10 years - sorted out already by using 'LONG' file)
#C fit a loess curve representing a smooth trend through the data
#C remove the smoothed fit from the data (detrends and removes inhomogeneous features hopefully)
#C standardise by dividing by climatological month st devs

get_climanoms_func<-function(absarr,spanval) {    
  
  bunk<-0
  
  longoo<-length(absarr)
  anomsarr<-array(NA,longoo)
  stanomsarr<-array(NA,longoo)
  smoothie<-array(NA,longoo)
  climsarr<-array(NA,12)
  sdsarr<-array(NA,12)
  
  valrebin<-array(absarr,dim=c(12,longoo/12))	# reform to 12 rows by nyears
 
  sizeclims<-apply(valrebin,c(1),function(x) sum(table(x)))	# should sum number of 'TRUE' values only

  # IF minimum month count is greater than 10 (out of 26 years of data) then we're good to go, IF not then we'll have to ditch the station
  if (min(sizeclims) >= 10) {

    # efficient creation of climate anoms (scale=FALSE means do not divide by st. dev)
    anomsarr<-array(scaltro<-t(apply(valrebin,1,scale,scale=FALSE)),dim=(longoo))  
    climsarr<-apply(valrebin,1,mean,na.rm=TRUE)

    smoothie<-predict(loess(anomsarr~seq(longoo),span=spanval),array(seq(longoo),dim=(longoo)))
    stanomsarr<-array(NA,dim=(longoo))
    stanomsarr[]<-anomsarr-smoothie

    sdsarr<-apply(array(stanomsarr,dim=c(12,longoo/12)),1,sd,na.rm=TRUE)
    stanomsarr<-array(t(apply(array(stanomsarr,dim=c(12,longoo/12)),1,scale)),dim=(longoo))

    #huges<-which(abs(stanomsarr) > 4.5)
    #if (length(huges) > 0) {
    #  stanomsarr[huges]<-NA
    #  # now rerun interpolate_missing_func to cover up the isolated missing outliers ***
    #  stanomsarr=interpolate_missing_func(stanomsarr) 	
    #}
    
  } else bunk<-9 #stop("THIS STATION IS BUNK!")
  return(list(climsarr=climsarr,sdsarr=sdsarr,anomsarr=anomsarr,stanomsarr=stanomsarr,smoothie=smoothie,bunk=bunk))
}

######################################################################################
# MAIN PROGRAM
######################################################################################
# TUNEABLE PARAMETERS------------------------------------------------------------------------------------

# LOESS FIT
# New longer series (1908 months) changes spanval a little
# NEW: 763 = 0.4,600/1908=0.31 400/1908=0.21, 290/1908 = 0.15
# OLD: 650/1368=0.45 500/1368=0.36  15*12=180(0.6), 12*12=144(0.48 or 0.5), 10*12=120(0.4), 8*12=96(0.32 or 0.3), 5*12=60(0.2)
spanmonthsFAT<-763		
spanmonthsNORM<-600
spanmonthsWAVY<-400		
spanmonthsWIGGLY<-290		
		

# EDITABLE PARAMETERS------------------------------------------------------------------------------------
nstations<-32522 		# fullLONG or selected?
styr     <-1860			# Start year for station
edyr     <-2018			# end year for station
# missing data indicator
mdi	 <-(-99.99)

# HadCM3 list:************************************
hadcm3list<-as.list(rep(list(0),14))
names(hadcm3list)<-c('dir','filee','decsarr',
                     'decst','deced','ndecs','modstyr',
		     'modedyr','nlats','nlons','modlats',
		     'modlons')
hadcm3list.dir<-"/data/local/hadkw/MODEL_DATA/A1B/TMEAN_nc/aenwl/"
hadcm3list.filee<-"Tmeandaily_2HadCM3qump_"
hadcm3list.decsarr  <-c("1970-79aenwl.nc",	# running from 1970 to 2083 to match data length
             	    "1980-89aenwl.nc",
	     	    "1990-99aenwl.nc",
	     	    "2000-09aenwl.nc",
	     	    "2010-19aenwl.nc",
	     	    "2020-29aenwl.nc",
	     	    "2030-39aenwl.nc",
	     	    "2040-49aenwl.nc",
	     	    "2050-59aenwl.nc",
	     	    "2060-69aenwl.nc",
	     	    "2070-79aenwl.nc",
	     	    "2080-89aenwl.nc")		# for reading in the model data
hadcm3list.decst <-c(1,1,1,1,1,1,1,1,1,1,1,1)	# start month for each dec
hadcm3list.deced <-c(120,120,120,120,120,120,120,120,120,120,120,48) # end month for each dec
hadcm3list.ndecs    <-12			# total number of decs (files)
hadcm3list.modstyr  <-1970	# GCM start year
hadcm3list.modedyr  <-2083	# GCM end year - probably want to stop this in 2083!
hadcm3list.nlats    <-73	# 2.5 deg lats give 73 boxes -91.25 to 91.25
hadcm3list.nlons    <-96			# 3.75 deg lons give 96 boxes

hadcm3list.modlats  <-array(seq(90,-90,-2.5))	# GRIDBOX CENTRES 90to-90
# prog will add one extra box and pad, no need to shift either - 0 or -180 ok
hadcm3list.modlons  <-array(seq(0,356.25,3.75))
# shift to be -180 to 180
hadcm3list.modlons  <- hadcm3list.modlons+180.	# need to shift modgrids lons + 48

# End of HadCM3 list******************************
# HadGEM2-ES historical+historicalExt list:************************************
# atmos.Amon.r2i1p1.
hadgemeslist<-as.list(rep(list(0),14))
names(hadgemeslist)<-c('dir','filee','decsarr',
                     'decst','deced','ndecs','modstyr',
		     'modedyr','nlats','nlons','modlats',
		     'modlons')
hadgemeslist.dir<-"/data/local/hadkw/CMIP5/HadGEM2ES/ATMOSAMON/r2i1p1/tas/"
hadgemeslist.filee<-"tas_Amon_HadGEM2-ES_"
hadgemeslist.decsarr  <-c("historical_r2i1p1_185912-188411.nc",	# running from 1970 to 2083 to match data length
             	    "historical_r2i1p1_188412-190911.nc",
	     	    "historical_r2i1p1_190912-193411.nc",
	     	    "historical_r2i1p1_193412-195911.nc",
	     	    "historical_r2i1p1_195912-198411.nc",
	     	    "historical_r2i1p1_198412-200512.nc",
	     	    "historicalExt_r2i1p1_200512-201911.nc")		# for reading in the model data
hadgemeslist.decst <-c(2,1,1,1,1,1,2)			# start month for each dec	
hadgemeslist.deced <-c(300,300,300,300,300,253,157)	# end month for each dec
hadgemeslist.ndecs    <-7			# number of decs (files)	
hadgemeslist.modstyr  <-1860	# GCM start year
hadgemeslist.modedyr  <-2018 	# GCM end year - probably want to stop this in 2083!
hadgemeslist.nlats    <-145	# 1.25 deg lats give 145 boxes -90.625 to 90.625
hadgemeslist.nlons    <-192	# 1.875 deg lons give 192 boxes -0.9375 to 359.0625 
hadgemeslist.modlats  <-array(seq(-90.,90.,1.25)) # GRIDBOX CENTRES 90to-90
# prog will add one extra box and pad, no need to shift either - 0 or -180 ok
hadgemeslist.modlons  <-array(seq(0,358.125,1.875))	# need to shift modgrids lons + 96
hadgemeslist.modlons  <-hadgemeslist.modlons-180.	# need to shift modgrids lons + 96
# End of HadGEM-ES3 list******************************

decsarr	 <-hadgemeslist.decsarr		# for reading in the model data
decst    <-hadgemeslist.decst									# for reading in the model data
deced    <-hadgemeslist.deced									# for reading in the model data
ndecs	 <-hadgemeslist.ndecs													# for reading in the model data
modstyr  <-hadgemeslist.modstyr		# GCM start year
modedyr  <-hadgemeslist.modedyr	# GCM end year - probably want to stop this in 2083!
nlats    <-hadgemeslist.nlats	# 2.5 deg lats give 73 boxes -91.25 to 91.25
nlons    <-hadgemeslist.nlons			# 3.75 deg lons give 96 boxes
modlats <-hadgemeslist.modlats	# GRIDBOX CENTRES
modlons <-hadgemeslist.modlons	# need to shift modgrids lons + add one box for overlapping?

# SET IN STONE PARAMETERS--------------------------------------------------------------------------------
nyrs	 <-(edyr-styr)+1		# number of years in station (input) data
nm 	 <-((edyr+1)-styr)*12		# number of time points (months) 1987-2011 is 25 years * 12 = 300
clims	 <-c(0,nm-1)			# whole period clims
nmodyrs  <-(modedyr-modstyr)+1		# GCM years
nmodmons <-nmodyrs*12 			# GCM months
nmoddays <-nmodmons*30 			# GCM days
modstmatches<-seq(nm)+((styr-modstyr)*12)	# a pointer to the model months to map over to the station if doing a pure match
#spanval<-spanmonths/nm			# Loess smoothing parameter for the station input number of months
spanvalmodFAT<-spanmonthsFAT/nmodmons		# 0.4 Loess smoothing parameter for the model input number of months
spanvalmodNORM<-spanmonthsNORM/nmodmons		# 0.31 Loess smoothing parameter for the model input number of months
spanvalmodWAVY<-spanmonthsWAVY/nmodmons		# 0.21 Loess smoothing parameter for the model input number of months
spanvalmodWIGGLY<-spanmonthsWIGGLY/nmodmons	# 0.15 Loess smoothing parameter for the model input number of months


# FILES AND DIRECTORIES------------------------------------------------------------------------------
dirlist<-"/data/local/hadkw/ISTI/LISTS/v101_JUL2015/"
dirGCMdata<-hadgemeslist.dir
dirINTERPdata<-"/data/local/hadkw/ISTI/LISTS/BAWG/"

inlist<-"ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat"	# GOOD station list
infil<-hadgemeslist.filee	
## Smooth anomaly method 2 CLS - climatology-loess-standardise
paramtag<-paste("loess0",round(spanvalmodFAT*100),"CLS",sep="")		# output filename tag
outfilFAT<-paste("HadGEM2ESLOESS_ISTI_stage3proxyelevs_",paramtag,"_JUL2015.txt",sep="") 
## Smooth anomaly method 2 CLS - climatology-loess-standardise
paramtag<-paste("loess0",round(spanvalmodNORM*100),"CLS",sep="")		# output filename tag
outfilNORM<-paste("HadGEM2ESLOESS_ISTI_stage3proxyelevs_",paramtag,"_JUL2015.txt",sep="") 
## Smooth anomaly method 2 CLS - climatology-loess-standardise
paramtag<-paste("loess0",round(spanvalmodWAVY*100),"CLS",sep="")		# output filename tag
outfilWAVY<-paste("HadGEM2ESLOESS_ISTI_stage3proxyelevs_",paramtag,"_JUL2015.txt",sep="") 
## Smooth anomaly method 2 CLS - climatology-loess-standardise
paramtag<-paste("loess0",round(spanvalmodWIGGLY*100),"CLS",sep="")		# output filename tag
outfilWIGGLY<-paste("HadGEM2ESLOESS_ISTI_stage3proxyelevs_",paramtag,"_JUL2015.txt",sep="") 


# VARIABLES AND ARRAYS-------------------------------------------------------------------------------
ds_info     <-list(statid=array("XXXXXX",dim=(nstations)),statlats=array(mdi,dim=(nstations)),
                   statlons=array(mdi,dim=(nstations)),statelvs=array(mdi,dim=(nstations)))

#modelworld<-array(mdi,dim=c(nlats,nlons,nmodmons))	# dimensions are longitudes, latitudes and times
modelstations<-array(mdi,dim=c(nstations,nmodmons))	# model field interpolated to each station - clim anoms, smoothed loess (probably a better way of doing this)
tmpmodelloess<-array(mdi,dim=c(nmodmons))	# dimensions are longitudes, latitudes and times

######################################################################################
# ACTION

# set up modlatarr, modlonarr  ROWWISE so rows(latitudes first) VERY IMPORTANT
modlatarr<-array(rep(modlats,nlons),c(nlats,nlons))
modlonarr<-t(array(rep(modlons,nlats),c(nlons,nlats)))

nclr<-9
plotclr<-brewer.pal(nclr,"RdPu")


# 1) READ IN LIST OF STATIONS INTO mcdw_info
mush<-readLines(con=paste(dirlist,inlist,sep=""),n=-1)	# read entire file
for (linoo in 1:nstations) {
  ds_info$statid[linoo]<-substr(mush[linoo],1,12)			# station ID
  ds_info$statlats[linoo]<-type.convert(substr(mush[linoo],68,75))	# station Latitude
  ds_info$statlons[linoo]<-type.convert(substr(mush[linoo],78,86))	# station Longitude
  ds_info$statelvs[linoo]<-type.convert(substr(mush[linoo],88,95))	# station elevnation
} 
rm(mush)

# 2) Read in GCM fields into array

# READ IN MODEL DATA FROM NETCDF for Equivalent length of data
# use library(ncdf)
# read data from 1970 to 2083 (pseudo 1900-2013) - in decades:
#  1970-79, 1980-89, 1990-99, 2000-09, 2010-19, 2020-29, 2030-39, 2040-49, 2050-59, 2060-69, 2070-79, 2080-89
# this is daily and so will need averaging to monthly
ptmst<-1	        # start month for wanted model time
ptmed<-(deced[1]-decst[1])+1	# end month for wanted model time
for (dd in 1:ndecs) {						# 
  print(paste(c(dd,ptmst,ptmed,decst[dd],deced[dd]),sep=" "))  

  modelworld<-array(NA,dim=c(nlats,nlons,(deced[dd]-decst[dd])+1))	# dimensions are longitudes, latitudes and times
  
  # set up file name for correct netCDF file and open and get data
  filoo<-paste(dirGCMdata,infil,decsarr[dd],sep="")
  inny<-open.ncdf(filoo)
  tmpmodmons<-get.var.ncdf(inny)
  
  # move lons to be -180 to 180 (NOTE tmpmodmons is lon, lat, tim
  tmpmodmons<-tmpmodmons[c(((nlons/2)+1):nlons,1:(nlons/2)),,]
  
  # HadGEM2 is monthly data. HadCM3 is daily? average to months and fill in appropriate amount into modelworld (93,73,nmodmonths)
  for (lnn in 1:nlons) {
    for (ltt in 1:nlats) {
#      subarr<-array(tmpmodmons[lnn,ltt,]-273.15,dim=c(30,120))	# these are all in kelvin, reform into months and years of days - a decade
#      modelworld[ltt,lnn,]<-colMeans(subarr[,decst[dd]:deced[dd]])	# average all days within the month to get the month mean
      modelworld[ltt,lnn,]<-tmpmodmons[lnn,ltt,decst[dd]:deced[dd]]-273.15	# average all days within the month to get the month mean
    }
  }
  # close netCDF file
  close.ncdf(inny)
  
  # 3) now do the interpolation to save memory
  # Use Finn's code to interpolate: should cope with 3D
  # needs to outout a station (row) by months (columns) array
  tmp<-padded.interpp(modlonarr,modlatarr,modelworld,ds_info$statlons,ds_info$statlats,padding=2)
  modelstations[,ptmst:ptmed]<-tmp

#  for (mm in decst[dd]:deced[dd]){
#    print(paste("Month:",mm,"of",deced[dd],sep=" "))
#   
#  # check this:
#    contour(modlons,modlats,t(modelworld[,,mm]))
#    colcode<-array(0,c(nstations))
#    maxt<-max(modelstations[,ptmst+(mm-decst[dd])])
#    mint<-min(modelstations[,ptmst+(mm-decst[dd])])-0.0001
#    ranget=maxt-mint
#    for (i in 1:nstations){
#      colcode[i]<-plotclr[ceiling(((modelstations[i,ptmst+(mm-decst[dd])]-mint)/ranget)/(1./9.))]
#    }
#    map("world")
#    points(ds_info$statlons,ds_info$statlats,pch=10,col=colcode,cex=0.5)
#
#    stop("Test ds_info$statlons[1],ds_info$statlats[1],modelstations[1,ptmst+(mm-decst[dd])] against plot")
#  
#  } # end of 'for' in 2)
  
  if (dd != ndecs) {
    ptmst<-ptmed+1
    ptmed<-ptmed+((deced[dd+1]-decst[dd+1])+1)					# start month for full model time
  }
#  if (dd == ndecs) {
#    ptmst<-ptmed+1
#    ptmed<-ptmed+((deced[1]-decst[1])+1)					# start month for full model time
#  }
}

print("Done sorting model")

# 4) Calculate climatology, remove, fit smooth loess, keep the curve

restarter<-"--------"	#"--------"

for (stattoo in 1:nstations){
    
  if (restarter != "--------" & restarter != ds_info$statid[stattoo]) {
    next
  }
  restarter<-"--------"
  
  print(paste("STATION:",stattoo,"of",nstations,sep=" "))
  
  # do two different levels of smooth - the same as obs and a lower amount as models are too smooth anyway
  
  res=get_climanoms_func(modelstations[stattoo,],spanvalmodFAT)
  modCLIMANOMstation<-res$anomsarr
  modSDANOMstation<-res$stanomsarr
  tmpmodelloess<-res$smoothie
  modCLIM<-res$climsarr
  modSD<-res$sdsarr
  bunk<-res$bunk
  
  #stop("Test loess smooth of GCM station by plotting and points and lines")
  
  # 5) WRITE TO FILE: GCM loess
  write.table(paste(c(ds_info$statid[stattoo],format(round(tmpmodelloess,digits=3),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirINTERPdata,outfilFAT,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

  res=get_climanoms_func(modelstations[stattoo,],spanvalmodNORM)
  modCLIMANOMstation<-res$anomsarr
  modSDANOMstation<-res$stanomsarr
  tmpmodelloess<-res$smoothie
  modCLIM<-res$climsarr
  modSD<-res$sdsarr
  bunk<-res$bunk
  
  #stop("Test loess smooth of GCM station by plotting and points and lines")
  
  # 5) WRITE TO FILE: GCM loess
  write.table(paste(c(ds_info$statid[stattoo],format(round(tmpmodelloess,digits=3),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirINTERPdata,outfilNORM,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

  res=get_climanoms_func(modelstations[stattoo,],spanvalmodWAVY)
  modCLIMANOMstation<-res$anomsarr
  modSDANOMstation<-res$stanomsarr
  tmpmodelloess<-res$smoothie
  modCLIM<-res$climsarr
  modSD<-res$sdsarr
  bunk<-res$bunk
  
  #stop("Test loess smooth of GCM station by plotting and points and lines")
  
  # 5) WRITE TO FILE: GCM loess
  write.table(paste(c(ds_info$statid[stattoo],format(round(tmpmodelloess,digits=3),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirINTERPdata,outfilWAVY,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

  res=get_climanoms_func(modelstations[stattoo,],spanvalmodWIGGLY)
  modCLIMANOMstation<-res$anomsarr
  modSDANOMstation<-res$stanomsarr
  tmpmodelloess<-res$smoothie
  modCLIM<-res$climsarr
  modSD<-res$sdsarr
  bunk<-res$bunk
  
  #stop("Test loess smooth of GCM station by plotting and points and lines")
  
  # 5) WRITE TO FILE: GCM loess
  write.table(paste(c(ds_info$statid[stattoo],format(round(tmpmodelloess,digits=3),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirINTERPdata,outfilWIGGLY,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

} # end of 'for' in 4)

######################################################################################
# END
######################################################################################
