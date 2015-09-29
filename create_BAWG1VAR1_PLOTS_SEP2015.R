# R
#
# Author: Kate Willett
# Created: 14 September 2015
# Last update: 29 September 2015
# Location: /data/local/hadkw/ISTI/PROGS/
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This code makes analysis plots from the clean world benchmarks
# It has 8 plots to run:
#	1) Scatter plots of old vs new cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
#	2) Histograms of old cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
#	3) Histograms of new cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
#	4) Histograms of old difference series sd and autocorrelation (st anoms and clim anoms)
#	5) Histograms of new difference series sd and autocorrelation (st anoms and clim anoms)
# 	6) Individual station climate anomalies from old and new, and four nearest neighbours, may add loess smooth on raw and clean (from GCM)
#	 
# Willett, K. M., C. N. Williams, I. Jolliffe, R. Lund, L. Alexander, S. Brönniman, L. A. Vincent, S. Easterbrook, V. Venema, 
# D. Berry, R. E. Warren, G. Lopardo, R. Auchmann, E. Aguilar, M. Menne, C. Gallagher, Z. Hausfather, T. Thorarinsdottir, 
# P. W. Thorne, 2014: A framework for benchmarking of homogenisation algorithm performance on the global scale, Geoscientific 
# Instrumentation, Methods and Data Systems, 3, 187-200, doi:10.5194/gi-3-187-2014.
# 
# Lund, R. and Willett, K. M.,,in prep.: Simulation of temperature networks from data. ???, .
#
# Willett, K. M., Lund, R and Chandler, R. E., in prep.: Simulating clean monthly mean surface temperature records on the 
# global scale. ???, .
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# R packages:
#
# Kate's R modules:
# read_station_func.R - reads in an ISTI station file, written by Kate Willett
# interpolate_missing_func.R - interpolates over single missing data points to extend data use for stats and AR(1) residuals, written by Kate Willett
# get_climanoms_func.R - cleans up data and creates standardised and climate anomalies and also returns linear trend, lowess trend, monthly climatological means and standard deviations, written by Kate Willett
#
# -----------------------
# DATA
# -----------------------
# Raw data are stored here:
# dirdata<-"/data/local/hadkw/ISTI/DATA/"
# 	Real ISTI (in ISTI format) data are read in from here:
# 	infilraw<-dirdata+"ISTIv101_JUL2015/results_merged/merge_"			
# Masked clean data are stored here:
# dirdata<-"/data/local/hadkw/ISTI/DATA/"
# 	Simulated ISTI (in ISTI format) data are read in from here:
# 	infilclean<-dirdata+"CLEANWORLDS/v101_JUL2015/ISTI_TYPE/merge_"			
# Station lists for the stations are stored here:
# dirlist<-"/data/local/hadkw/ISTI/LISTS/v101_JUL2015/"
# 	Reduced ISTI stage 3 list to stations with >= 3 years of data, no ships, no location matches??:
# 	infillist<-dirlist+"ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat"	# Hoping that we just need the one list now - should all be simulatable
#	List of Neighbours:
#	infilNEIGHlist<-paste("CORRNEIGHBOURS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
#	Chosen smoothed climate anomaly curve from a GCM: (there are a choice of wigglinesses) - see READMEBAWG_SEP2015
#	infilGCMloess<-dirlist+"HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess015CLS_JUL2015.txt"	 
# dirstats<-"/data/local/hadkw/ISTI/LISTS/BAWG/SEP2015/" 
# 	OLD Covariance matrix of station with neighbours at lag 0 - CHECK - already done in python for Clim Anoms and St Anoms
#	infilOLDcovs<-paste("OLDCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	OLD Covariance matrix of station with neighbours at lag 1 - CHECK - clready done in python for Clim Anoms and St Anoms
#	infilOLDl1covs<-paste("OLDCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	NEW Covariance matrix of station with neighbours at lag 0
#	infilNEWcovs<-paste("NEWCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	NEW Covariance matrix of station with neighbours at lag 1
#	infilNEWl1covs<-paste("NEWCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
#
# 	DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - standardised anoms
# 	infilOLDDIFFSDstats<-paste("OLDDIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	infilOLDDIFFACstats<-paste("OLDDIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - climate anoms
#	infilOLDDIFFSDCMstats<-paste("OLDDIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	infilOLDDIFFACCMstats<-paste("OLDDIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - standardised anoms
#	infilNEWDIFFSDstats<-paste("DIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	infilNEWDIFFACstats<-paste("DIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - climate anoms
#	infilNEWDIFFSDCMstats<-paste("DIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	infilNEWDIFFACCMstats<-paste("DIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#		
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# This code can be run in one go or for selected parts only
# Full run:
#	Set all part_numbers to 0
#	Ensure that all sections of create_BAWG1VAR1_MAIN_SEP2015.R have been run
#	Ensure filepaths are correct
#	Ensure output goes to correct directory version
#	source("create_BAWG1VAR1_PLOTS_SEP2015.R")
# Part run:
#	As above in terms of file paths and pre-preparing files
# 	Set desired part_numbers to 1
# 	source("create_BAWG1VAR1_PLOTS_SEP2015.R")
#
# -----------------------
# OUTPUT
# Plots are stored here:
# dirplot<-"/data/local/hadkw/ISTI/IMAGES/SEP2015/STATIONS/"
#	1) Scatter plots of old vs new cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
#	outplotSCATTERStAn="Scatter_Covs_StAnoms_BNCHCAAA_SEP2015.eps"
#	outplotSCATTERClAn="Scatter_Covs_ClAnoms_BNCHCAAA_SEP2015.eps"
#	2) Histograms of new cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
#	outplotHISTNEWCOVStAn="Hist_NEW_Covs_StAnoms_BNCHCAAA_SEP2015.eps"
#	outplotHISTNEWCOVClAn="Hist_NEW_Covs_ClAnoms_BNCHCAAA_SEP2015.eps"
#	3) Histograms of old cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
#	outplotHISTOLDCOVStAn="Hist_OLD_Covs_StAnoms_BNCHCAAA_SEP2015.eps"
#	outplotHISTOLDCOVClAn="Hist_OLD_Covs_ClAnoms_BNCHCAAA_SEP2015.eps"
#	4) Histograms of new difference series sd and autocorrelation (st anoms and clim anoms)
#	outplotHISTNEWDIFFStAn="Hist_NEW_Diffs_StAnoms_BNCHCAAA_SEP2015.eps"
#	outplotHISTNEWDIFFClAn="Hist_NEW_Diffs_ClAnoms_BNCHCAAA_SEP2015.eps"
#	5) Histograms of old difference series sd and autocorrelation (st anoms and clim anoms)
#	outplotHISTOLDDIFFStAn="Hist_OLD_Diffs_StAnoms_BNCHCAAA_SEP2015.eps"
#	outplotHISTOLDDIFFClAn="Hist_OLD_Diffs_ClAnoms_BNCHCAAA_SEP2015.eps"
# 	6) Individual station climate anomalies from old and new, and four nearest neighbours, may add loess smooth on raw and clean (from GCM)
#	outplotSTATION=IDNUMBER999"_OLDvsNEW_ClAnoms_BNCHCAAA_SEP2015.eps"
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 29th September
# ---------
#  
# Enhancements
# Added output to histograms for the mean and standard deviation of distributions
# Also output to STDOUT the quantiles of the distribution
# Station time series plots now start from 1800 for the Real ISTI databank
#  
# Changes
#  
# Bug fixes
#
# Version 1 21st September
# ---------
#  
# Enhancements
#  
# Changes
#  
# Bug fixes
#  
# -----------------------
# OTHER INFORMATION
# -----------------------
#
########################################################################################
# MAIN PROGRAM 
#############################################################################################
options(warn=2)

# call functions that are needed:
source("read_station_func.R")
source("interpolate_missing_func.R")
source("get_climanoms_func.R")
#----------------------------------------------------------------------------------------

restarter<-"--------"	#"--------"

#----------------------------------------------------------------------------------------
# set up common variables and arrays

# TUNEABLE PARAMETERS

# EDITABLE PARAMETERS
styr     <-1800			# Start year for station (was 1850)
edyr     <-2015			# end year for station
nstations<-32522 		# Only stations with sufficient correlating neighbours to create a VAR model (22697,22256)
modstyr  <-1800	# GCM start year is 1860 so need to reverse 1860 to 1920 when adding GCM loess!!!
RealModStYr<-1860
ModExt=(RealModStYr-modstyr)*12
modedyr  <-2018	# GCM end year - probably want to stop this in 2083!

# missing data indicator
mdi	 <-(-99.99)

# SET IN STONE PARAMETERS
nyrs	 <-(edyr-styr)+1		# number of years in station (input) data
nm 	 <-((edyr+1)-styr)*12		# number of time points (months) 1987-2011 is 25 years * 12 = 300
clims	 <-c(0,nm-1)			# whole period clims
nmodyrs  <-(modedyr-modstyr)+1		# GCM years
nmodmons <-nmodyrs*12 			# GCM months
nmoddays <-nmodmons*30 			# GCM days
#spanval<-spanmonths/nm			# Loess smoothing parameter for the station input number of months
#spanvalmod<-spanmonths/nmodmons		# Loess smoothing parameter for the model input number of months

## Project Name
paramtag<-"BNCHCAAA"		# output filename tag

# plotting labels and ticks
xyearsnull<-c(array("",dim=(nyrs+3)))		# empty array for year labels
yrlist<-seq(styr-1,edyr+2) 			# integer array of years ffor station input
labs<-which((yrlist/5.)-floor(yrlist/5.) == 0.)	# pointer for creating year labels for those years ending in 0 or 5
miniyears<-yrlist[labs]				# built year labels for station input
xtickies<-(seq(nyrs+3)*12)-23			# tick mark locations for station input
minitickies<-xtickies[labs]			# pointer to tick marks for labelled years only
xmodyearsnull<-c(array("",dim=(nmodyrs+3)))		       # empty array for year labels for GCM data/output
yrmodlist<-seq(modstyr-1,modedyr+2) 	       # integer array of years for GCM data/output
modlabs<-which((yrmodlist/5.)-floor(yrmodlist/5.) == 0.)       # pointer for creating year labels for those years ending in 0 or 5
minimodyears<-yrmodlist[modlabs]			       # built year labels for GCM data/output
xmodtickies<-(seq(nmodyrs+3)*12)-23			       # tick mark locations for GCM data/output
minimodtickies<-xmodtickies[modlabs]			       # pointer to tick marks for labelled years only for GCM data/output

#-----------------------------------------------------------
# SET UP FILE PATHS AND FILES @ enric: I had to modify this to run in 
dirdata<-"/data/local/hadkw/ISTI/DATA/"
dirlist<-"/data/local/hadkw/ISTI/LISTS/v101_JUL2015/"
dirstats<-"/data/local/hadkw/ISTI/LISTS/BAWG/SEP2015/"
dirplot<-"/data/local/hadkw/ISTI/IMAGES/SEP2015/"

infilraw <-"ISTIv101_JUL2015/results_merged/merge_"				# station data
infillist<-"ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat"	# Hoping that we just need the one list now - should all be simulatable

# List of Neighbours (by distance - 40 nearest)
infilNEIGHlist<-paste("CORRNEIGHBOURS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 

# There are four of these: loess04, 0.31, 0.21 and loess015 - need to test wigglyness - does it lead to inhomogeneity/over correlation?.
infilGCMloess<-"HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess015CLS_JUL2015.txt"	 

# OLD Covariance matrix of station with neighbours at lag 0 - CHECK - already done in python for Clim Anoms and St Anoms
infilOLDcovs<-paste("OLDCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# OLD Covariance matrix of station with neighbours at lag 1 - CHECK - clready done in python for Clim Anoms and St Anoms
infilOLDl1covs<-paste("OLDCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# NEW Covariance matrix of station with neighbours at lag 0
infilNEWcovs<-paste("NEWCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# NEW Covariance matrix of station with neighbours at lag 1
infilNEWl1covs<-paste("NEWCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 

# DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - standardised anoms
infilOLDDIFFSDstats<-paste("OLDDIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
infilOLDDIFFACstats<-paste("OLDDIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - climate anoms
infilOLDDIFFSDCMstats<-paste("OLDDIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
infilOLDDIFFACCMstats<-paste("OLDDIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - standardised anoms
infilNEWDIFFSDstats<-paste("DIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
infilNEWDIFFACstats<-paste("DIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - climate anoms
infilNEWDIFFSDCMstats<-paste("DIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
infilNEWDIFFACCMstats<-paste("DIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")

infilclean <-paste("CLEANWORLDS/v101_JUL2015/ISTI_TYPE/merge_",sep="") 
# may later change to paramtag+'_'

#1) Scatter plots of old vs new cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
outplotSCATTERStAn="Scatter_Covs_StAnoms_BNCHCAAA_SEP2015.eps"
outplotSCATTERClAn="Scatter_Covs_ClAnoms_BNCHCAAA_SEP2015.eps"
#2) Histograms of old cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
outplotHISTOLDCOVStAn="Hist_OLD_Covs_StAnoms_BNCHCAAA_SEP2015.eps"
outplotHISTOLDCOVClAn="Hist_OLD_Covs_ClAnoms_BNCHCAAA_SEP2015.eps"
#3) Histograms of new cross-correlations at lag 0 and lag 1 (st anoms and clim anoms)
outplotHISTNEWCOVStAn="Hist_NEW_Covs_StAnoms_BNCHCAAA_SEP2015.eps"
outplotHISTNEWCOVClAn="Hist_NEW_Covs_ClAnoms_BNCHCAAA_SEP2015.eps"
#4) Histograms of old difference series sd and autocorrelation (st anoms and clim anoms)
outplotHISTOLDDIFFStAn="Hist_OLD_Diffs_StAnoms_BNCHCAAA_SEP2015.eps"
outplotHISTOLDDIFFClAn="Hist_OLD_Diffs_ClAnoms_BNCHCAAA_SEP2015.eps"
#5) Histograms of new difference series sd and autocorrelation (st anoms and clim anoms)
outplotHISTNEWDIFFStAn="Hist_NEW_Diffs_StAnoms_BNCHCAAA_SEP2015.eps"
outplotHISTNEWDIFFClAn="Hist_NEW_Diffs_ClAnoms_BNCHCAAA_SEP2015.eps"
#6) Individual station climate anomalies from old and new, and four nearest neighbours, may add loess smooth on raw and clean (from GCM)
cleanplotdir="STATIONS/"
outplotSTATION="_OLDvsNEW_ClAnoms_BNCHCAAA_SEP2015.eps"
#----------------------------------------------------------

############################################################################################
#partone<-0	# switch - if 0, scatter plots of covs
#		# if 1 - do not do anything
parttwo<-1	# switch - if 0, histograms of new covs
		# if 1 - do not do anything
partthree<-1	# switch - if 0, histograms of old covs
		# if 1 - do not do anything
partfour<-0	# switch - if 0, histograms of new diffs
		# if 1 - do not do anything
partfive<-0	# if 0, histograms of old diffs
		# if 1 - do not do anything
partsix<-1	# if 0, for each station and four nearest neighbours plot climate anomaly time series and distribution - old vs new
		# if 1 - do not do anything
###########################################################################################
if (parttwo == 0) {
  print("Part Two 0 - plots of NEW covs")
  # set up full arrays
  FullCCclimanom<-0
  FullCCsdanom<-0
  FullCCl1climanom<-0
  FullCCl1sdanom<-0
  CCmeanCA<-0
  CCsdCA<-0
  CCl1meanCA<-0
  CCl1sdCA<-0
  CCmeanSA<-0
  CCsdSA<-0
  CCl1meanSA<-0
  CCl1sdSA<-0

  # read in covs lag 0
  print('Reading in NEW COVS')
  
  mycols <- rep("NULL",81)
  mycols[2:81] <- "character"
  mooCOV<-read.table(paste(dirstats,infilNEWcovs,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullCCclimanom<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullCCclimanom<-FullCCclimanom[which(FullCCclimanom != "NA")]
  FullCCsdanom<-as.numeric(unlist(unname(mooCOV[,41:80])))
  FullCCsdanom<-FullCCsdanom[which(FullCCsdanom != "NA")]
  
  mooCOV<-read.table(paste(dirstats,infilNEWl1covs,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullCCl1climanom<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullCCl1climanom<-FullCCl1climanom[which(FullCCl1climanom != "NA")]
  FullCCl1sdanom<-as.numeric(unlist(unname(mooCOV[,41:80])))
  FullCCl1sdanom<-FullCCl1sdanom[which(FullCCl1sdanom != "NA")]

  rm(mooCOV,mycols)
  gc()
  
  # Get means and standard deviations of distributions
  CCmeanCA<-mean(FullCCclimanom)
  CCsdCA<-sd(FullCCclimanom)
  CCl1meanCA<-mean(FullCCl1climanom)
  CCl1sdCA<-sd(FullCCl1climanom)
  CCmeanSA<-mean(FullCCsdanom)
  CCsdSA<-sd(FullCCsdanom)
  CCl1meanSA<-mean(FullCCl1sdanom)
  CCl1sdSA<-sd(FullCCl1sdanom)
  
  # Print out Quantiles
  print("St Anoms: lag 0s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCsdanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("Clim Anoms: lag 0s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCclimanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)

  print("St Anoms: lag 1s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCl1sdanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("Clim Anoms: lag 1s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCl1climanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  # now make the plot and save
  setEPS()
  postscript(paste(dirplot,outplotHISTNEWCOVStAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  hist(FullCCsdanom[-1],main="Cross Correlations at lag 0 (New Std Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('c)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCmeanSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCsdSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  hist(FullCCl1sdanom[-1],main="Cross Correlations at lag 1 (New Std Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1))
  mtext('d)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCl1meanSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCl1sdSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()

  setEPS()
  postscript(paste(dirplot,outplotHISTNEWCOVClAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  hist(FullCCclimanom[-1],main="Cross Correlations at lag 0 (New Clim Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('c)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCmeanCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCsdCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  hist(FullCCl1climanom[-1],main="Cross Correlations at lag 1 (New Clim Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1))
  mtext('d)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCl1meanCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCl1sdCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()

  rm(FullCCclimanom,FullCCsdanom,FullCCl1climanom,FullCCl1sdanom)
  gc()

}
#######################################################################################
if (partthree == 0) {
  print("Part Three 0 - plots of OLD covs")
  # set up full arrays
  FullCCclimanom<-0
  FullCCsdanom<-0
  FullCCl1climanom<-0
  FullCCl1sdanom<-0
  CCmeanCA<-0
  CCsdCA<-0
  CCl1meanCA<-0
  CCl1sdCA<-0
  CCmeanSA<-0
  CCsdSA<-0
  CCl1meanSA<-0
  CCl1sdSA<-0

  # read in covs lag 0
  print('Reading in OLD COVS')
  
  mycols <- rep("NULL",81)
  mycols[2:81] <- "character"
  mooCOV<-read.table(paste(dirstats,infilOLDcovs,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullCCclimanom<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullCCclimanom<-FullCCclimanom[which(FullCCclimanom != "NA")]
  FullCCsdanom<-as.numeric(unlist(unname(mooCOV[,41:80])))
  FullCCsdanom<-FullCCsdanom[which(FullCCsdanom != "NA")]
  
  mooCOV<-read.table(paste(dirstats,infilOLDl1covs,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullCCl1climanom<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullCCl1climanom<-FullCCl1climanom[which(FullCCl1climanom != "NA")]
  FullCCl1sdanom<-as.numeric(unlist(unname(mooCOV[,41:80])))
  FullCCl1sdanom<-FullCCl1sdanom[which(FullCCl1sdanom != "NA")]

  rm(mooCOV,mycols)
  gc()

  # Get means and standard deviations of distributions
  CCmeanCA<-mean(FullCCclimanom)
  CCsdCA<-sd(FullCCclimanom)
  CCl1meanCA<-mean(FullCCl1climanom)
  CCl1sdCA<-sd(FullCCl1climanom)
  CCmeanSA<-mean(FullCCsdanom)
  CCsdSA<-sd(FullCCsdanom)
  CCl1meanSA<-mean(FullCCl1sdanom)
  CCl1sdSA<-sd(FullCCl1sdanom)
  
  # Print out Quantiles
  print("St Anoms: lag 0s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCsdanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("Clim Anoms: lag 0s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCclimanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)

  print("St Anoms: lag 1s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCl1sdanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("Clim Anoms: lag 1s")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullCCl1climanom,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)

  # now make the plot and save
  setEPS()
  postscript(paste(dirplot,outplotHISTOLDCOVStAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  hist(FullCCsdanom[-1],main="Cross Correlations at lag 0 (Old Std Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('a)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCmeanSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCsdSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  hist(FullCCl1sdanom[-1],main="Cross Correlations at lag 1 (Old Std Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1))
  mtext('b)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCl1meanSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCl1sdSA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()

  setEPS()
  postscript(paste(dirplot,outplotHISTOLDCOVClAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  hist(FullCCclimanom[-1],main="Cross Correlations at lag 0 (Old Clim Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('a)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCmeanCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCsdCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  hist(FullCCl1climanom[-1],main="Cross Correlations at lag 1 (Old Clim Anoms)",xlab="Correlation", ylab="Frequency",breaks=seq(-0.8,1.,0.1))
  mtext('b)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.1)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(CCl1meanCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(CCl1sdCA,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()

  rm(FullCCclimanom,FullCCsdanom,FullCCl1climanom,FullCCl1sdanom)
  gc()

}
#######################################################################################
if (partfour == 0) { # end of part two 
  print("Part Four:0: Plotting new diff series stats")
  # set up storage
  FullDiffSD<-0  	# growing array to save all diffs SDs to plot as hist (REMOVE 1st ELEMENT)
  FullDiffAC<-0	# growing array to save all diffs lag 1 cor to plot as hist (REMOVE 1st ELEMENT)
  SDmean<-0
  SDsd<-0
  ACmean<-0
  ACsd<-0
  
  # read in st anoms
  print('Reading in NEW DIFFS and plot')
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilNEWDIFFSDstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffSD<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffSD<-FullDiffSD[which(FullDiffSD != "NA")]
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilNEWDIFFACstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffAC<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffAC<-FullDiffAC[which(FullDiffAC != "NA")]

  rm(mooCOV,mycols)
  gc()

  # Get means and standard deviations of distributions
  SDmean<-mean(FullDiffSD)
  SDsd<-sd(FullDiffSD)
  ACmean<-mean(FullDiffAC)
  ACsd<-sd(FullDiffAC)
  
  # Print out Quantiles
  print("St Anoms: St Devs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffSD,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("St Anoms: ACs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffAC,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)

  # now make the plot and save
  setEPS()
  postscript(paste(dirplot,outplotHISTNEWDIFFStAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  #pin(8,10)
  #par(mar=c(5,3,2,2)+0.1)
  #mai(c(0.7,1.2,1,0.7)) # bottom, left, top, right
  hist(FullDiffSD[-1],main="Difference Series St Dev (New Std Anoms)",xlab="Standard Deviation", ylab="Frequency",breaks=seq(0,12,0.2)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('c)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(SDmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(SDsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  #par(mar=c(5,3,2,2)+0.1)
  #mai(c(0.7,1.2,1,0.7)) # bottom, left, top, right
  hist(FullDiffAC[-1],main="Difference Series Lag 1 Autocorrelation (New Std Anoms)",xlab="Autocorrelation (lag 1)", ylab="Frequency",breaks=seq(-0.6,1.,0.1))
  mtext('d)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(ACmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(ACsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()

  rm(FullDiffSD,FullDiffAC)
  gc()
    
  # reset storage
  FullDiffSD<-0  	# growing array to save all diffs SDs to plot as hist (REMOVE 1st ELEMENT)
  FullDiffAC<-0	# growing array to save all diffs lag 1 cor to plot as hist (REMOVE 1st ELEMENT)  
  SDmean<-0
  SDsd<-0
  ACmean<-0
  ACsd<-0
  
  # read in clim anoms
  print('Reading in NEW DIFFS and plot')
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilNEWDIFFSDCMstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffSD<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffSD<-FullDiffSD[which(FullDiffSD != "NA")]
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilNEWDIFFACCMstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffAC<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffAC<-FullDiffAC[which(FullDiffAC != "NA")]

  rm(mooCOV,mycols)
  gc()

  # Get means and standard deviations of distributions
  SDmean<-mean(FullDiffSD)
  SDsd<-sd(FullDiffSD)
  ACmean<-mean(FullDiffAC)
  ACsd<-sd(FullDiffAC)
  
  # Print out Quantiles
  print("Clim Anoms: St Devs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffSD,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("Clim Anoms: ACs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffAC,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)

  # now make the plot and save
  setEPS()
  postscript(paste(dirplot,outplotHISTNEWDIFFClAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  hist(FullDiffSD[-1],main="Difference Series St Dev (New Clim Anoms)",xlab="Standard Deviation (degrees C)", ylab="Frequency",breaks=seq(0,12,0.2)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('c)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(SDmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(SDsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  hist(FullDiffAC[-1],main="Difference Series Lag 1 Autocorrelation (New Clim Anoms)",xlab="Autocorrelation (lag 1)", ylab="Frequency",breaks=seq(-0.6,1.,0.1))
  mtext('d)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(ACmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(ACsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()
  rm(FullDiffSD,FullDiffAC)
  gc()

} 
#stop()
################################################################################################

if (partfive == 0) { # end of part two 
  print("Part Five:0: Plotting old diff series stats")
  # set up storage
  FullDiffSD<-0  	# growing array to save all diffs SDs to plot as hist (REMOVE 1st ELEMENT)
  FullDiffAC<-0	# growing array to save all diffs lag 1 cor to plot as hist (REMOVE 1st ELEMENT)
  SDmean<-0
  SDsd<-0
  ACmean<-0
  ACsd<-0
  
  # read in st anoms
  print('Reading in OLD DIFFS and plot')
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilOLDDIFFSDstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffSD<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffSD<-FullDiffSD[which(FullDiffSD != "NA")]
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilOLDDIFFACstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffAC<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffAC<-FullDiffAC[which(FullDiffAC != "NA")]

  rm(mooCOV,mycols)
  gc()

  # Get means and standard deviations of distributions
  SDmean<-mean(FullDiffSD)
  SDsd<-sd(FullDiffSD)
  ACmean<-mean(FullDiffAC)
  ACsd<-sd(FullDiffAC)
  
  # Print out Quantiles
  print("St Anoms: St Devs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffSD,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("St Anoms: ACs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffAC,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)

  # now make the plot and save
  setEPS()
  postscript(paste(dirplot,outplotHISTOLDDIFFStAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  #pin(8,10)
  #par(mar=c(5,3,2,2)+0.1)
  #mai(c(0.7,1.2,1,0.7)) # bottom, left, top, right
  hist(FullDiffSD[-1],main="Difference Series St Dev (Old Std Anoms)",xlab="Standard Deviation", ylab="Frequency",breaks=seq(0,12,0.2)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('a)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(SDmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(SDsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  #par(mar=c(5,3,2,2)+0.1)
  #mai(c(0.7,1.2,1,0.7)) # bottom, left, top, right
  hist(FullDiffAC[-1],main="Difference Series Lag 1 Autocorrelation (Old Std Anoms)",xlab="Autocorrelation (lag 1)", ylab="Frequency",breaks=seq(-0.6,1.,0.1))
  mtext('b)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(ACmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(ACsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()
  rm(FullDiffSD,FullDiffAC)
  gc()
    
  # reset storage
  FullDiffSD<-0  	# growing array to save all diffs SDs to plot as hist (REMOVE 1st ELEMENT)
  FullDiffAC<-0	# growing array to save all diffs lag 1 cor to plot as hist (REMOVE 1st ELEMENT)  
  SDmean<-0
  SDsd<-0
  ACmean<-0
  ACsd<-0
  
  # read in clim anoms
  print('Reading in OLD DIFFS and plot')
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilOLDDIFFSDCMstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffSD<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffSD<-FullDiffSD[which(FullDiffSD != "NA")]
  
  mycols <- rep("NULL",41)
  mycols[2:41] <- "character"
  mooCOV<-read.table(paste(dirstats,infilOLDDIFFACCMstats,sep=""),colClasses=mycols)	# each line is now a character string to be split
  FullDiffAC<-as.numeric(unlist(unname(mooCOV[,1:40])))
  FullDiffAC<-FullDiffAC[which(FullDiffAC != "NA")]

  rm(mooCOV,mycols)
  gc()

  # Get means and standard deviations of distributions
  SDmean<-mean(FullDiffSD)
  SDsd<-sd(FullDiffSD)
  ACmean<-mean(FullDiffAC)
  ACsd<-sd(FullDiffAC)
  
  # Print out Quantiles
  print("Clim Anoms: St Devs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffSD,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)
  
  print("Clim Anoms: ACs")
  qporbs<-seq(100)/100.
  qtls<-quantile(FullDiffAC,qporbs) # gives the % of points less than each percentile threshold
  print(qtls)

  # now make the plot and save
  setEPS()
  postscript(paste(dirplot,outplotHISTOLDDIFFClAn,sep=""),width=6, height=8)

  par(mfrow=c(2,1))
  hist(FullDiffSD[-1],main="Difference Series St Dev (Old Clim Anoms)",xlab="Standard Deviation (degrees C)", ylab="Frequency",breaks=seq(0,12,0.2)) #  xlim=c(xmin, xmax), ylim=c(ymin, ymax)
  mtext('a)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(SDmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(SDsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  hist(FullDiffAC[-1],main="Difference Series Lag 1 Autocorrelation (Old Clim Anoms)",xlab="Autocorrelation (lag 1)", ylab="Frequency",breaks=seq(-0.6,1.,0.1))
  mtext('b)',side=3,adj=0,padj=0)
  plotcoords<-par("usr") # returns (xleft,xright, ybottom, ytop)
  xpos=((plotcoords[2]-plotcoords[1])*0.7)+plotcoords[1]
  ypos=((plotcoords[4]-plotcoords[3])*0.9)+plotcoords[3]
  ypos2=((plotcoords[4]-plotcoords[3])*0.8)+plotcoords[3]
  text(xpos,ypos,paste('Mean = ',format(round(ACmean,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  text(xpos,ypos2,paste('St Dev = ',format(round(ACsd,digits=2),trim=FALSE,nsmall=2,width=6),sep=""),pos=4)
  
  dev.off()
  rm(FullDiffSD,FullDiffAC)
  gc()

} 
#stop()
################################################################################################
if (partsix == 0) {
  # make plots inidividual stations within a gridbox comparison (SLOW - INVOLVES A LOT OF READING IN AT THIS STAGE)
  # STRUCTURES FOR STORING STATION INFORMATION
  # STATIONS ON ROWS, TIME ON COLUMNS ACTUALLY - THIS SHOULD BE THE OTHER WAY AROUND!!!
  # LIST of arrays (vectors) and data frames like an IDL structure for the dataset
  ds_info     <-list(statid=array("XXXXXX",dim=(nstations)),statlats=array(mdi,dim=(nstations)),
                     statlons=array(mdi,dim=(nstations)),statelvs=array(mdi,dim=(nstations)))
  corrneighbours.id<-vector("list",nstations)
  
  Letteree1<-c('a)','b)','c)','d)','e)')
  Letteree2<-c('f)','g)','h)','i)','j)')

  tmpStation<-array(mdi,nm)
  tmpNEWStation<-array(mdi,nmodmons)
  tmpNEWSmooth<-array(mdi,nmodmons)
  tmpNeighbour<-array(mdi,c(4,nm))
  tmpNEWNeighbour<-array(mdi,c(4,nmodmons))
  tmpNEWNeighbourSmooth<-array(mdi,c(4,nmodmons))
  
  #--------------------------------------------------------------------------------
  # READ IN LIST OF STATIONS 
  mush<-readLines(con=paste(dirlist,infillist,sep=""),n=-1)	# read entire file
  for (linoo in 1:nstations) {
    ds_info$statid[linoo]<-substr(mush[linoo],2,12)			# station ID
    ds_info$statlats[linoo]<-type.convert(substr(mush[linoo],68,75))	# station Latitude
    ds_info$statlons[linoo]<-type.convert(substr(mush[linoo],78,86))	# station Longitude
    ds_info$statelvs[linoo]<-type.convert(substr(mush[linoo],88,95))	# station elevnation
  }
  nStations=length(ds_info$statid) 
  rm(mush)
  #-----------------------------------------------------------------

  # Read in GCM smooths in prep
  mooGCM<-read.table(paste(dirstats,infilGCMloess,sep=""),header=FALSE)	# each line is now a character string to be split
# NOTE: Extending GCM back beyond 1860 by flipping the loess

  # Read in neighbour list ( 40 nearest stations) in prep
  mooNEIGHS<-read.table(paste(dirstats,infilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))

  print("Finding Corrneighbours")  
  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 
  rm(moo,mooNEIGHS,founds)
  gc()  
  
  # Loop through each station
  for (nS in 1:nStations) {
    
    print(ds_info$statid[nS])
    
    # empty previously filled arrays
    tmpStation<-array(mdi,nm)
    tmpNEWStation<-array(mdi,nmodmons)
    tmpNEWSmooth<-array(mdi,nmodmons)
    tmpNeighbour<-array(mdi,c(4,nm))
    tmpNEWNeighbour<-array(mdi,c(4,nmodmons))
    tmpNEWNeighbourSmooth<-array(mdi,c(4,nmodmons))
    
    # read in station Old and get climate anomaly
    tmpStation=read_station_func(styr,paste(dirdata,infilraw,ds_info$statid[nS],"_stage3",sep=""),tmpStation)
    tmpStation[which(tmpStation == mdi)]<-NA
    res=get_climanoms_func(tmpStation,9) # any number for spanval with do
    tmpStation<-res$anomsarr
    
    # read in station New climate anomaly
    tmpNEWStation=read_station_func(modstyr,paste(dirdata,infilclean,ds_info$statid[nS],"_stage3",sep=""),tmpNEWStation)
    tmpNEWStation[which(tmpNEWStation == mdi)]<-NA
    res=get_climanoms_func(tmpNEWStation,9) # any number for spanval with do
    tmpNEWStation<-res$anomsarr
    
    # read in station GCM loess smooth
    tmpNEWSmooth<-unname(as.numeric(mooGCM[nS,-1]))
    if (RealModStYr != modstyr) {
      #extend tmploess
      tmpNEWSmooth<-c(rev(tmpNEWSmooth[1:ModExt]),tmpNEWSmooth) 
    }  
    
    # Use first four nearest neighbours in each station's corrneighbours_id
    for (nN in 1:4) {
      # Read in four nearest neighbours Old climate anomalies
      tmpNeighbour[nN,]=read_station_func(styr,paste(dirdata,infilraw,ds_info$statid[corrneighbours.id[[nS]][nN]],"_stage3",sep=""),tmpNeighbour[nN,])
      tmpNeighbour[nN,which(tmpNeighbour[nN,] == mdi)]<-NA
      res=get_climanoms_func(tmpNeighbour[nN,],9) # any number for spanval with do
      tmpNeighbour[nN,]<-res$anomsarr
    
      # Read in four nearest neighbours New climate anomlies
      tmpNEWNeighbour[nN,]=read_station_func(modstyr,paste(dirdata,infilclean,ds_info$statid[corrneighbours.id[[nS]][nN]],"_stage3",sep=""),tmpNEWNeighbour[nN,])
      tmpNEWNeighbour[nN,which(tmpNEWNeighbour[nN,] == mdi)]<-NA
      res=get_climanoms_func(tmpNEWNeighbour[nN,],9) # any number for spanval with do
      tmpNEWNeighbour[nN,]<-res$anomsarr
    
      # Read in GCM loess smooth for four nearest stations
      tmpNEWNeighbourSmooth[nN,(ModExt+1):nmodmons]<-unname(as.numeric(mooGCM[corrneighbours.id[[nS]][nN],-1]))
      if (RealModStYr != modstyr) {
        #extend tmploess
        tmpNEWNeighbourSmooth[nN,]<-c(rev(tmpNEWNeighbourSmooth[nN,(ModExt+1):(ModExt+ModExt)]),tmpNEWNeighbourSmooth[nN,(ModExt+1):nmodmons]) 
      }  
    }
    
    # Plot 5 by 2 panel plot for Old (col 1) and New (col 2) time series with GCM loess overlaying new
    setEPS()
    postscript(paste(dirplot,cleanplotdir,ds_info$statid[nS],outplotSTATION,sep=""),width=12, height=8)

    par(mfcol=c(5,2))	# fills by columns
    # par(fig=c(0,0.8,0,0.8), new=TRUE) - use to set x1, x2, y1, y2 pos of plot explicitly
    par(mar=c(4,4,2.5,1),mgp=c(2,0.7,0)) # margins set for bottom, left, top, right in lines of text and distance of (label, ticklabels, ticks)
    plot(seq(nm),tmpStation,type="b",xaxt="n",cex.lab=1,xlim=c(1,(nyrs*12)+1),main=paste("Old Clim Anoms:",ds_info$statid[nS],
           format(round(ds_info$statlats[nS],digits=3),trim=FALSE,nsmall=3,width=7),
	   format(round(ds_info$statlons[nS],digits=3),trim=FALSE,nsmall=3,width=8),
	   format(round(ds_info$statelvs[nS],digits=0),trim=FALSE,nsmall=0,width=5),sep=" "),ylab="Degrees C",xlab="Year")
    axis(1,xtickies,labels=FALSE,tck=0.015,cex.axis=1)
    axis(1,minitickies,miniyears,tck=0.03,cex.axis=1)
    axis(3,xtickies,xyearsnull,tck=0.015)
    axis(3,minitickies,xyearsnull[labs],tck=0.03)
    mtext(Letteree1[1],side=3,line=1,adj=0,padj=0)
    for (nR in 2:5) {
#      par(mar=c(2,2,2,2)) # margins set for bottom, left, top, right in lines of text
      plot(seq(nm),tmpNeighbour[nR-1,],type="b",xaxt="n",cex.lab=1,xlim=c(1,(nyrs*12)+1),main=paste("Old Clim Anoms:",ds_info$statid[corrneighbours.id[[nS]][nR-1]],
           format(round(ds_info$statlats[corrneighbours.id[[nS]][nR-1]],digits=3),trim=FALSE,nsmall=3,width=7),
	   format(round(ds_info$statlons[corrneighbours.id[[nS]][nR-1]],digits=3),trim=FALSE,nsmall=3,width=8),
	   format(round(ds_info$statelvs[corrneighbours.id[[nS]][nR-1]],digits=0),trim=FALSE,nsmall=0,width=5),sep=" "),ylab="Degrees C",xlab="Year")
      axis(1,xtickies,labels=FALSE,tck=0.015,cex.axis=1)
      axis(1,minitickies,miniyears,tck=0.03,cex.axis=1)
      axis(3,xtickies,xyearsnull,tck=0.015)
      axis(3,minitickies,xyearsnull[labs],tck=0.03)
      mtext(Letteree1[nR],side=3,line=1,adj=0,padj=0)
    }
#    par(mar=c(2,2,2,2)) # margins set for bottom, left, top, right in lines of text
    plot(seq(nmodmons),tmpNEWStation,type="b",xaxt="n",cex.lab=1,xlim=c(1,(nmodyrs*12)+1),main="New Clim Anoms",ylab="Degrees C",xlab="Year")
    lines(seq(nmodmons),tmpNEWSmooth,col='Red')
    axis(1,xmodtickies,labels=FALSE,tck=0.015,cex.axis=1)
    axis(1,minimodtickies,minimodyears,tck=0.03,cex.axis=1)
    axis(3,xmodtickies,xmodyearsnull,tck=0.015)
    axis(3,minimodtickies,xmodyearsnull[modlabs],tck=0.03)
    mtext(Letteree2[1],side=3,line=1,adj=0,padj=0)
    for (nR in 2:5) {
#      par(mar=c(2,2,2,2)) # margins set for bottom, left, top, right in lines of text
      plot(seq(nmodmons),tmpNEWNeighbour[nR-1,],type="b",xaxt="n",cex.lab=1,xlim=c(1,(nmodyrs*12)+1),main="New Clim Anoms",ylab="Degrees C",xlab="Year")
      lines(seq(nmodmons),tmpNEWNeighbourSmooth[nR-1,],col='Red')
      axis(1,xmodtickies,labels=FALSE,tck=0.015,cex.axis=1)
      axis(1,minimodtickies,minimodyears,tck=0.03,cex.axis=1)
      axis(3,xmodtickies,xmodyearsnull,tck=0.015)
      axis(3,minimodtickies,xmodyearsnull[modlabs],tck=0.03)
      mtext(Letteree2[nR],side=3,line=1,adj=0,padj=0)
    }
    dev.off()
    
    # Plot 5 by 2 panel plot for Old (col 1) and New (col 2) distributions of climate anomalies 
#    stop()
  }
}  

########################################################################################
#######################################################################################
# END
##########################################################################################
