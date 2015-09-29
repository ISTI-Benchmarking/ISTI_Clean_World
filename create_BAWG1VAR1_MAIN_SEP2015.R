# R
#
# Author: Kate Willett
# Created: 14 September 2015
# Last update: 14 September 2015
# Location: /data/local/hadkw/ISTI/PROGS/
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This code performs the ISTI Benchmark Clean World simulation
# It assumes that all of the pre-required station lists, correlation lists and distance lists have already been set up - see READMEBAWG_SEP2015
# It has 9 processes to run:
#	0) Clean up data (remove outliers) and calculate real station statistics: 
#		PRODUCES:
#		standardised anomalies (removed climatology, removed lowess trend, divided by climatological monthly standard deviation)
#		climate anomalies (removed climatology)
#		lowess trend (fitted to the anomalies using a 600 month span)
#		linear decadal trend
#		climatology (using entire station record - should be at least three years)
#		standard deviation (using entire station record, should be at least three years) 
#		autocorrelation at lag 1 (of both the standardised anomalies and climate anomalies)
#		40 nearest neighbour list
#		AR(1) residuals (or MVN random numbers if real data are not long/complete enough - 60 month minimum), and mean and standard deviation of those residuals
#	1) Use the distance elevation function (already built - see READMEBAWG_SEP2015) to produce VAR(1) parameters for each station+40 neighbour matrix:
#		PRODUCES:
#		Gamma_lag0
#		Gamma_lag1
#		VAR(1) parameter
#	2) Create MVN residual shocks for all time points/stations using factorisation and the distance elevation function, 
#	   Transform to shape of real resids by Q-Q mapping
#	   Give them real mean and standard deviation 
#		PRODUCES:
#		New residual shocks
#	3) Run VAR(1) to create the simulated standardised anomalies using the neighbour disconnect method
#		PRODUCES:
#		Simulated Standardised anomalies
#	4) Add a GCM lowess trend (already created - see READMEBAWG_SEP2015)
#	   Multiply by real climatological monthly standard deviations to create simulated climate anomalies 
#	   Calculated simulated clean world stats: linear decadal trend, climatology, standard deviation, autocorrelation at lag 1
#	   Add back real monthly climatology to create simulated monthly means
#	   Look at nearest neighbour difference series for simulated data and calculate the standard deviation and autocorrelation at lag 1 for standardised anomalies and climate anomalies
#		PRODUCES:
#		Simulated Climate anomalies
#		Simulated Monthly means
#		Simulated Clean world stats:
#		Simulated Standardised anomaly differences series standard deviations
#		Simulated Standardised anomaly autocorrelations at lag 1
#		Simulated Climate anomaly difference series standard deviations
#		Simulated Climate anomaly difference series autocorrelations at lag 1
#	5) Look at Simulated nearest neighbour correlations at lag 0 and lag 1 for standardised anomalies and climate anomalies
#		Neighbour correlations at lag 0 for Simulated Standardised anomalies and Simulated Climate anomalies 		 
#		Neighbour correlations at lag 1 for Simulated Standardised anomalies and Simulated Climate anomalies 		 
#	6) Look at nearest neighbour difference series for real data and calculate the standard deviation and autocorrelation at lag 1 for standardised anomalies and climate anomalies
#		Real Standardised anomaly differences series standard deviations
#		Real Standardised anomaly autocorrelations at lag 1
#		Real Climate anomaly difference series standard deviations
#		Real Climate anomaly difference series autocorrelations at lag 1
#	7) Look at Real nearest neighbour correlations at lag 0 and lag 1 for standardised anomalies and climate anomalies
#		Neighbour correlations at lag 0 for Real Standardised anomalies and Real Climate anomalies 		 
#		Neighbour correlations at lag 1 for Real Standardised anomalies and Real Climate anomalies 	
#	8) Read in masked simulated data and recalculate all statistics
#		NOTE*** you must have run ConvertISTItoGHCN_NOV2014.py before running this section!!! ***
#		PRODUCES:
#		linear decadal trend
#		climatology (using entire station record - should be at least three years)
#		standard deviation (using entire station record, should be at least three years) 
#		autocorrelation at lag 1 (of both the standardised anomalies and climate anomalies)
#
# See create_BAWG1VAR1_PLOTS_SEP2015.R for plotting code
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
# BigMVNgen.r - performs factorisation, written by Richard Chandler
# QMapTransform_func.R - performs Q-Q mapping, written by Kate Willett
# read_station_func.R - reads in an ISTI station file, written by Kate Willett
# interpolate_missing_func.R - interpolates over single missing data points to extend data use for stats and AR(1) residuals, written by Kate Willett
# get_climanoms_func.R - cleans up data and creates standardised and climate anomalies and also returns linear trend, lowess trend, monthly climatological means and standard deviations, written by Kate Willett
# get_ARresids_func.R - if there is a chunk of consecutative data long enough (60 months) - use AR(1) to work out real resids, or 
#			substitute with MVN proxy resids, written by Kate Willett
# get_autocovs_func.R - uses distance function to create lag 0 and lag 1 cross-correlations and VAR(1) parameter
# distance_elevation_func.R - contains three functions for getting correlations based on distance and elevation, written by Kate Willett and Richard Chandler
#	expdh.corr - calulates a matrix of correlations for a set of stations given distance (km), elevation (km) and a function (matern-ish)
#	howfar - calculates horizontal distance in kilometers of stations with each other - for expdh.corr
#	howhigh - calculates vertical distance in kilometers of stations with each other - for expdh.corr
# -----------------------
# DATA
# -----------------------
# Raw data are stored here:
# dirdata<-"/data/local/hadkw/ISTI/DATA/"
# 	Real ISTI (in ISTI format) data are read in from here:
# 	infilraw<-dirdata+"ISTIv101_JUL2015/results_merged/merge_"			
# Post-processing read in masked clean data for redoing stats: (run part eight after python code to mask)
#	infilNEWMASKED <-paste("CLEANWORLDS/v101_JUL2015/MASKEDCLEAN_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# Various lists for the stations are stored here:
# dirlist<-"/data/local/hadkw/ISTI/LISTS/v101_JUL2015/"
# 	Reduced ISTI stage 3 list to stations with >= 3 years of data, no ships, no location matches??:
# 	infillist<-dirlist+"ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat"	# Hoping that we just need the one list now - should all be simulatable
#	List of 100 nearest stations for Reduced ISTI list:
#	infildists<-dirlist+"ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat"	# Distance station list
#	List of stations for each Factorisation network (within and neighbouring) - see READMEBAWG_SEP2015
#	inNTWKstats<-dirlist+"ISTILONGNetwork_StOnly_stage3proxyelevs_JUL2015_"	# stations within each network
#	inNTWKnbs<-dirlist+"ISTILONGNetwork_NbOnly_stage3proxyelevs_JUL2015_"	# neighbours for each network
# Input GCM and then any file read in that was output during processing:
# dirbawg<-"/data/local/hadkw/ISTI/LISTS/BAWG/SEP2015/"
#	Chosen smoothed climate anomaly curve from a GCM: (there are a choice of wigglinesses) - see READMEBAWG_SEP2015
#	infilGCMloess<-dirlist+"HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess015CLS_JUL2015.txt"	 
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# This code can be run in one go or for selected parts only
# There are various versions of each part should you need to repeat at various stages
# Full run:
#	Set all part_numbers to 0
#	Ensure filepaths are correct
#	Ensure that reduced station list, distance file, factorisation networks, GCM curves have been created
#	Ensure that distance/elevation function has been defined
#	Ensure output goes to correct directory version
#	source("create_BAWG1VAR1_MAIN_SEP2015.R")
# Part run:
#	As above in terms of file paths and pre-preparing files
# 	Set desired part_numbers to 1 or 2
# 	source("create_BAWG1VAR1_MAIN_SEP2015.R")
#
# -----------------------
# OUTPUT
# Cleanworld data are stored here:
# outfildata <-paste(dirdata,"CLEANWORLDS/v101_JUL2015/CLEAN_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# All output lists and stats are stored here:
# dirbawg<-"/data/local/hadkw/ISTI/LISTS/BAWG/SEP2015/"
# 	Correlation Neighbour list for each station based on 40 nearest neighbours
#	outfilNEIGHlist<-paste("CORRNEIGHBOURS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
#
# 	VAR parameters for each station r1c1,r2c1,r3c1,r1c2,r2c2,r3c2 etc
#	outfilarps<-paste("VARPS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	GAM0 parameters for each station r1c1,r2c1,r3c1,r1c2,r2c2,r3c2 etc
#	outfilgam0s<-paste("GAM0S_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	GAM1 parameters for each station r1c1,r2c1,r3c1,r1c2,r2c2,r3c2 etc
#	outfilgam1s<-paste("GAM1S_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
#
# 	OLD Covariance matrix of station with neighbours at lag 0 - CHECK - already done in python for Clim Anoms and St Anoms
#	outfilOLDcovs<-paste("OLDCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	OLD Covariance matrix of station with neighbours at lag 1 - CHECK - clready done in python for Clim Anoms and St Anoms
#	outfilOLDl1covs<-paste("OLDCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	NEW Covariance matrix of station with neighbours at lag 0
#	outfilNEWcovs<-paste("NEWCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	NEW Covariance matrix of station with neighbours at lag 1
#	outfilNEWl1covs<-paste("NEWCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
#
# 	Smooth (loess) curve removed from each station
#	outfilloess<-paste("LOESS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# 	OLD Whole period Trend, Annual clim, month clims, annual St Dev, month st devs, Clim anom and St. Anom autocorrelation for each station - CHECK - needs code change
#	outfilOLDstats<-paste("OLDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	NEW Whole period Month climatologies and standard deviations, Clim anom and St. Anom autocorrelation for each station - CHECK - needs code change
#	outfilNEWstats<-paste("NEWSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	NEW Masked Whole period Month climatologies and standard deviations, Clim anom and St. Anom autocorrelation for each station - CHECK - needs code change
#	outfilNEWMASKstats<-paste("NEWMASKEDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#
# 	NOT IN USE: INTERIM ACs for VAR test for each station - CHECK IS THIS STILL NEEDED???
#	outfilINTstats<-paste("INTERIMSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#
# 	OLD shock (residuals) for each station stationID, mean, sd, values...
#	outfilOLDshocks<-paste("OLDSHOCKS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	NEW shock (residuals) for each station (converted to actual distribution from MVN), mean, sd, values
#	outfilNEWshocks<-paste("NEWSHOCKS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#
# 	OLD climate anomalies for all stations (rows) with time points (columns) (LARGE FILE)	
#	outfilOLDcanoms <-paste("OLDCLIMANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	
# 	OLD standardised anomalies (loess removed) for all stations (rows) with time points (columns) (LARGE FILE)	
#	outfilOLDsdanoms <-paste("OLDSDANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	
# 	NEW climate anomalies for all stations (rows) with time points (columns) (LARGE FILE)	
#	outfilNEWcanoms <-paste("NEWCLIMANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	
# 	NEW standardised anomalies (loess removed) for all stations (rows) with time points (columns) (LARGE FILE)	
#	outfilNEWsdanoms <-paste("NEWSDANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	
#
# 	DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - standardised anoms
# 	outfilOLDDIFFSDstats<-paste("OLDDIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	outfilOLDDIFFACstats<-paste("OLDDIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - climate anoms
#	outfilOLDDIFFSDCMstats<-paste("OLDCOVDIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	outfilOLDDIFFACCMstats<-paste("OLDCOVDIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - standardised anoms
#	outfilDIFFSDstats<-paste("DIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	outfilDIFFACstats<-paste("DIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# 	DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - climate anoms
#	outfilDIFFSDCMstats<-paste("DIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#	outfilDIFFACCMstats<-paste("DIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
#		
#
# Plots are stored here:
# dirplot<-"/data/local/hadkw/ISTI/IMAGES/SEP2015/"
# see create_BAWG1VAR1_PLOTS_SEP2015.R
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 14th September
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
#----------------------------------------------------------------------------------------
# call functions that are needed:
source("QMapTransform_func.R")
source("BigMVNgen.r")
source("read_station_func.R")
source("interpolate_missing_func.R")
source("get_climanoms_func.R")
source("get_ARresids_func.R")
source("get_autocovs_func.R")
source("distance_elevation_func.R") 
#----------------------------------------------------------------------------------------
# set up common variables and arrays

restarter<-"--------"	#"--------"

# TUNEABLE PARAMETERS

# LOESS FIT
# for 1850-2015 (Only use data from 1850 onwards), for 1800 to 2015 use 800
spanmonths<-600	#600/1992=0.30, 787/1968=0.4

# Exponential decay function for deriving spatial correlations modelled on real data with a fixed sigma
CC<-c(0.97,0.0006,0.07)	# with a 1:1 fix at 1.0
CCDiag=1.0
CCl1<-c(0.23,0.0006,0.07)	# with a 1:1 fix at 0.23
CCl1Diag=0.23

# EDITABLE PARAMETERS
styr     <-1850			# Start year for station (was 1850)
edyr     <-2015			# end year for station
nstations<-32522 		# Only stations with sufficient correlating neighbours to create a VAR model (22697,22256)
nnetworks<-192
modstyr  <-1800	# GCM start year is 1860 so need to reverse 1860 to 1920 when adding GCM loess!!!
RealModStYr<-1860
ModExt=(RealModStYr-modstyr)*12
modedyr  <-2018	# GCM end year - probably want to stop this in 2083!
monlabels<-array(c("01","02","03","04","05","06","07","08","09","10","11","12"))

# missing data indicator
mdi	 <-(-99.99)

# SET IN STONE PARAMETERS
nyrs	 <-(edyr-styr)+1		# number of years in station (input) data
nm 	 <-((edyr+1)-styr)*12		# number of time points (months) 1987-2011 is 25 years * 12 = 300
clims	 <-c(0,nm-1)			# whole period clims
nmodyrs  <-(modedyr-modstyr)+1		# GCM years
nmodmons <-nmodyrs*12 			# GCM months
nmoddays <-nmodmons*30 			# GCM days
spanval<-spanmonths/nm			# Loess smoothing parameter for the station input number of months
spanvalmod<-spanmonths/nmodmons		# Loess smoothing parameter for the model input number of months

## Project Name
paramtag<-"BNCHCAAA"		# output filename tag

#-----------------------------------------------------------
# SET UP FILE PATHS AND FILES @ enric: I had to modify this to run in 
dirdata<-"/data/local/hadkw/ISTI/DATA/"
dirlist<-"/data/local/hadkw/ISTI/LISTS/v101_JUL2015/"
dirbawg<-"/data/local/hadkw/ISTI/LISTS/BAWG/SEP2015/"
dirplot<-"/data/local/hadkw/ISTI/IMAGES/SEP2015/"

infilraw <-"ISTIv101_JUL2015/results_merged/merge_"				# station data
infillist<-"ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat"	# Hoping that we just need the one list now - should all be simulatable

infildists<-"ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat"	# Distance station list
inNTWKstats<-"ISTILONGNetwork_StOnly_stage3proxyelevs_JUL2015_"	# stations within each network
inNTWKnbs<-"ISTILONGNetwork_NbOnly_stage3proxyelevs_JUL2015_"	# neighbours for each network

# There are four of these: loess04, 0.31, 0.21 and loess015 - need to test wigglyness - does it lead to inhomogeneity/over correlation?.
infilGCMloess<-"HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess015CLS_JUL2015.txt"	 

# Correlation Neighbour list for each station (<=40 highest correlating neighbours CHECK - based on distance?)
outfilNEIGHlist<-paste("CORRNEIGHBOURS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# VAR parameters for each station r1c1,r2c1,r3c1,r1c2,r2c2,r3c2 etc
outfilarps<-paste("VARPS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# GAM0 parameters for each station r1c1,r2c1,r3c1,r1c2,r2c2,r3c2 etc
outfilgam0s<-paste("GAM0S_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# GAM1 parameters for each station r1c1,r2c1,r3c1,r1c2,r2c2,r3c2 etc
outfilgam1s<-paste("GAM1S_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 

# OLD Covariance matrix of station with neighbours at lag 0 - CHECK - already done in python for Clim Anoms and St Anoms
outfilOLDcovs<-paste("OLDCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# OLD Covariance matrix of station with neighbours at lag 1 - CHECK - clready done in python for Clim Anoms and St Anoms
outfilOLDl1covs<-paste("OLDCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# NEW Covariance matrix of station with neighbours at lag 0
outfilNEWcovs<-paste("NEWCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# NEW Covariance matrix of station with neighbours at lag 1
outfilNEWl1covs<-paste("NEWCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 

# Smooth (loess) curve removed from each station
outfilloess<-paste("LOESS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
# OLD Whole period Trend, Annual clim, month clims, annual St Dev, month st devs, Clim anom and St. Anom autocorrelation for each station - CHECK - needs code change
outfilOLDstats<-paste("OLDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# NEW Whole period Month climatologies and standard deviations, Clim anom and St. Anom autocorrelation for each station - CHECK - needs code change
outfilNEWstats<-paste("NEWSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# NEW MASKED DATA Whole period Month climatologies and standard deviations, Clim anom and St. Anom autocorrelation for each station - CHECK - needs code change
outfilNEWMASKstats<-paste("NEWMASKEDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")

# INTERIM ACs for VAR test for each station - NOT CURRENTLY IN USE
outfilINTstats<-paste("INTERIMSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")

# OLD shock (residuals) for each station stationID, mean, sd, values...
outfilOLDshocks<-paste("OLDSHOCKS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# NEW shock (residuals) for each station (converted to actual distribution from MVN), mean, sd, values
outfilNEWshocks<-paste("NEWSHOCKS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")

# OLD climate anomalies for all stations (rows) with time points (columns) (LARGE FILE)	
outfilOLDcanoms <-paste("OLDCLIMANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	
# OLD standardised anomalies (loess removed) for all stations (rows) with time points (columns) (LARGE FILE)	
outfilOLDsdanoms <-paste("OLDSDANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	
# NEW climate anomalies for all stations (rows) with time points (columns) (LARGE FILE)	
outfilNEWcanoms <-paste("NEWCLIMANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	
# NEW standardised anomalies (loess removed) for all stations (rows) with time points (columns) (LARGE FILE)	
outfilNEWsdanoms <-paste("NEWSDANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")   	

# DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - standardised anoms
outfilOLDDIFFSDstats<-paste("OLDDIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
outfilOLDDIFFACstats<-paste("OLDDIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# DIFF STATS OLD SD and ACs for each station's station-neighbour diff series - climate anoms
outfilOLDDIFFSDCMstats<-paste("OLDDIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
outfilOLDDIFFACCMstats<-paste("OLDDIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - standardised anoms
outfilDIFFSDstats<-paste("DIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
outfilDIFFACstats<-paste("DIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
# DIFF STATS NEW SD and ACs for each station's station-neighbour diff series - climate anoms
outfilDIFFSDCMstats<-paste("DIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")
outfilDIFFACCMstats<-paste("DIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="")

outfildata <-paste("CLEANWORLDS/v101_JUL2015/CLEAN_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 
infilNEWMASKED <-paste("CLEANWORLDS/v101_JUL2015/MASKEDCLEAN_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt",sep="") 

#----------------------------------------------------------
# STRUCTURES FOR STORING STATION INFORMATION
# STATIONS ON ROWS, TIME ON COLUMNS ACTUALLY - THIS SHOULD BE THE OTHER WAY AROUND!!!
# LIST of arrays (vectors) and data frames like an IDL structure for the dataset
ds_info     <-list(statid=array("XXXXXX",dim=(nstations)),statlats=array(mdi,dim=(nstations)),
                   statlons=array(mdi,dim=(nstations)),statelvs=array(mdi,dim=(nstations)))
DistStatID<-0	# list of IDs in the Correlations file to search through to find right line to read

#--------------------------------------------------------------------------------
# READ IN LIST OF STATIONS 
mush<-readLines(con=paste(dirlist,infillist,sep=""),n=-1)	# read entire file
for (linoo in 1:nstations) {
  ds_info$statid[linoo]<-substr(mush[linoo],2,12)			# station ID
  ds_info$statlats[linoo]<-type.convert(substr(mush[linoo],68,75))	# station Latitude
  ds_info$statlons[linoo]<-type.convert(substr(mush[linoo],78,86))	# station Longitude
  ds_info$statelvs[linoo]<-type.convert(substr(mush[linoo],88,95))	# station elevnation
} 
rm(mush)

############################################################################################
partzero<-2	# switch - if 0, run part 0 (clean up, get station stats, get ARresids - transform to MVN?
		# if 1, read in neighbour network
		# if 2 or more - do not do anything
partone<-2	# switch - if 0, run part 1 (get VAR parameters)
		# if 1, read in neighbour network station stats and VARS and resids from file
		# if 2 or more - do not do anything
parttwo<-2	# switch - if 0, get new residuals and reverse transform
		# if 1, read in station stats and VARS and resids and NEWshocks from file
		# if 2 or more - do not do anything
partthree<-2	# switch - if 0, produce new CleanWorld standardised anomalies
		# if 1, read in CleanWorld standardised anomalies and station stats from file
partfour<-3	# switch - if 0, produce new CleanWorld climate anomalies and monthly means - and diffs
		# if 1, read in new CleanWorld climate anomalies and do diff stats on clim anoms (if messed up in previous loop)
		# if 2, read in new CleanWorld climate anomalies
partfive<-1	# if 0, new covs

partsix<-1	# if 0, old covs

partseven<-0	# if 0, old diffs

parteight<-1	# if 0, New masked stats - only run if you have prepared a maaksed data version

###########################################################################################

if (partzero == 0) {	# Need to run from the beginning
  print("Part Zero:0")
  # Set up initial storage
  # maps for locating stations within the main list ds_info
  corrneighbours.id<-vector("list",nstations)	# map using full list
  
  # while all real data will be stored in temporary arrays the output data must be stored for good or printed to file/read in again.
  # OLD parameters
  tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
  tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
  tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
  tmpLOESS<-array(mdi,nm)		# temporary array for storing the loess smoothed timeseries removed from the data
  tmpCNstations<-array(mdi,dim=c(40,nm))	# temporary array for storing neighbour time series - abs, clim anoms, std anoms
  tmpACclimanom<-0		# temporary scalar for station autocorrelation (climate anomalies with loess) 
  tmpACsdanom<-0			# temporary scalar for station autocorrelation (standardised anomalies) 
  tmpresids<-0			# temporary array for storing candidate station residuals
  tmpresidsCLIM<-0.		# element for storing candidate station residuals MEAN
  tmpresidsSD<-0.			# element for storing candidate station residuals Sd

  stationTREND<-array(0,c(nstations))	# array for station climatology, annual and then months
  stationCLIM<-array(0,c(nstations,13))	# array for station climatology, annual and then months
  stationSD<-array(0,c(nstations,13))	# array for station standard deviation, all period and then months

  # get stations stats, neighbour lists, and AR resids

  # first read in list of stations in Dists file to search through to find right line
  mycols <- rep("NULL",201)
  mycols[1] <- "character"
  DistStatID <- read.table(paste(dirlist,infildists,sep=""),colClasses=mycols)

  for (stattoo in 1:nstations){

    if (restarter != "--------" & restarter != ds_info$statid[stattoo]) {
      next
    }
    restarter<-"--------"

    # reset all tmp arrays
    tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
    tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
    tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
    tmpLOESS<-array(mdi,nm)		# temporary array for storing the loess smoothed timeseries removed from the data
    tmpCNstations<-array(mdi,dim=c(40,nm))	# temporary array for storing neighbour time series - abs, clim anoms, std anoms
    tmpACclimanom<-0		# temporary scalar for station autocorrelation (climate anomalies with loess) 
    tmpACsdanom<-0			# temporary scalar for station autocorrelation (standardised anomalies) 
    tmpCCclimanom<-array(0,40)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)
    tmpCCsdanom<-array(0,40)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)
    tmpresids<-0			# temporary array for storing candidate station residuals
    tmpresidsCLIM<-0.		# element for storing candidate station residuals MEAN
    tmpresidsSD<-0.			# element for storing candidate station residuals Sd

    print(paste("Working on station ",stattoo,": ",ds_info$statid[stattoo],sep=""))
  
    # read in the candidate station and filter out data to correct time points
    tmpABSstation=read_station_func(styr,paste(dirdata,infilraw,ds_info$statid[stattoo],"_stage3",sep=""),tmpABSstation)
  
    # INTERPOLATE OVER MISSING VALUES WHERE THEY ARE ISOLATED CASES OF 1-2 MAXIMUM CONSECUTIVE MISSING DATA
    tmpABSstation[which(tmpABSstation == mdi)]<-NA   # set missing values to R standard mdi
    tmpABSstation=interpolate_missing_func(tmpABSstation) 	

    # CREATE STATION STATS: clims, st devs, anoms, residuals from a smoothed filter
    # method CLS - calculate Climatology and remove, fit Loess and remove, calculated Standard devations and divide by them.
    res=get_climanoms_func(tmpABSstation,spanval)
    tmpCLIMANOMstation<-res$anomsarr
    tmpSDANOMstation<-res$stanomsarr
    tmpLOESS<-res$smoothie
    stationTREND[stattoo]<-res$trendy
    stationCLIM[stattoo,]<-res$climsarr
    stationSD[stattoo,]<-res$sdsarr
      
    # calculate autocorrelation at lag 1 for climate anomalies and standardised (loess removed) anomalies
    moo<-acf(tmpCLIMANOMstation,lag.max=24,plot=FALSE,na.action=na.pass)
    tmpACclimanom<-moo$acf[2:25,1,1]
    moo<-acf(tmpSDANOMstation,lag.max=24,plot=FALSE,na.action=na.pass)
    tmpACsdanom<-moo$acf[2:25,1,1]
    # test for optimal AR order using aic - not perfect method
    # goo<-ar(tmpSTANOMstation[which(is.na(tmpSTANOMstation == "FALSE")],aic=TRUE)

    # FIND and read in line from distance file to get closest neighbours
    foundit<-which(DistStatID == ds_info$statid[stattoo])
    moo<-scan(paste(dirlist,infildists,sep=""),what="list",skip=(foundit-1),nlines=1)
    
    # rearrange moo to get just the stations in order of distance and take the first 40 closest neighbours
    TheNeighbours<-array(moo[2:201],c(2,100))[1,1:40]
    #corrneighbours.id[[stattoo]]<-(1:nstations)[ds_info$statid %in% TheNeighbours]		# make a list pointer to the FULL list of stations
    corrneighbours.id[[stattoo]]<-match(TheNeighbours,ds_info$statid)

    # get REAL or MVN AR resids

    res=get_ARresids_func(tmpSDANOMstation)  
    tmpresids<-res$ar_z_resids
    RealResids<-res$RealResids
    
    # WRITE TO FILE: station trend, ann and month clims, ann and month st devs, clim anom ac, stanom ac, smooth line, clim. anoms and st. anoms
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpLOESS,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
    	file=paste(dirbawg,outfilloess,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpCLIMANOMstation,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
    	file=paste(dirbawg,outfilOLDcanoms,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpSDANOMstation,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
    	file=paste(dirbawg,outfilOLDsdanoms,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(stationTREND[stattoo],digits=2),trim=FALSE,nsmall=2,width=6),
        format(round(stationCLIM[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6),
    	format(round(stationSD[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6),
    	format(round(tmpACclimanom[1],digits=3),trim=FALSE,nsmall=3,width=6),
    	format(round(tmpACsdanom[1],digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
    	file=paste(dirbawg,outfilOLDstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    # WRITE TO FILE: list of neighbours used for this station (in order) and LONGEST COUNTS and cross correlations   
    write.table(paste(c(ds_info$statid[stattoo],ds_info$statid[corrneighbours.id[[stattoo]]]),collapse=" "),
    	  file=paste(dirbawg,outfilNEIGHlist,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

# 12) get residuals stats, and store with mean and sd used to transform new ones later (not standardised here)
    tmpresidsCLIM<-mean(tmpresids)
    tmpresidsSD<-sd(tmpresids)

    # WRITE TO FILE: residuals original (+ clim and sd)
    write.table(paste(c(ds_info$statid[stattoo],RealResids,format(round(tmpresidsCLIM,digits=3),trim=FALSE,nsmall=3,width=6),
    	  format(round(tmpresidsSD,digits=3),trim=FALSE,nsmall=3,width=6),
    	  format(round(tmpresids,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
    	  file=paste(dirbawg,outfilOLDshocks,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
    #stop() # only using this with restarter
  
  } # end of 'for'
    
  # clean up
  rm(DistStatID,tmpABSstation,tmpCLIMANOMstation,tmpSDANOMstation,tmpLOESS,tmpCNstations,
   tmpACclimanom,tmpACsdanom,tmpresids,tmpresidsCLIM,tmpresidsSD,stationCLIM,stationSD)
  gc()
} else if (partzero == 1) { # end of partone read in VARPS, station stats, neighbours list
  # set up necessary variables and arrays

  # read in corrneighbours.id
  print("Part Zero:1 - Reading in station corrneighbours")
  # set up storage
  # maps for locating stations within the main list ds_info
  corrneighbours.id<-vector("list",nstations)	# map using full list

  mooNEIGHS<-read.table(paste(dirbawg,outfilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
#  founds<-unname(as.matrix(mooNEIGHS[,-c(1,42,43,44)]))
  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))
  nNeighbours<-40	# 4 elements of candidate station ID and LONGEST COUNTS: 99
  
  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 
  
  rm(mooNEIGHS,moo,founds)
  gc()
}
#stop()
#############################################################################
if (partone == 0) {	# Need to run from the beginning
  print("Part One:0")
  # create varps
  tmpstationCoords<-array(mdi,dim=c(41,2))    # temporary array for storing lats, longs of network
  tmpstationElevs<-array(mdi,dim=c(41))	    # temporary array for storing elevations of network
  tmpGAM1s<-array(NA,c(41,41))
  tmpGAM0s<-array(NA,c(41,41))

  VARphi.vals=vector("list",nstations)	# container for VAR autoregressive parameter matrices for each station

  for (stattoo in 1:nstations){

    if (restarter != "--------" & restarter != ds_info$statid[stattoo]) {
      next
    }
    restarter<-"--------"

    # reset all tmp arrays
    tmpstationCoords<-array(mdi,dim=c(41,2))	# temporary array for storing lats, longs of network
    tmpstationElevs<-array(mdi,dim=c(41))	# temporary array for storing elevations of network

    print(paste("Working on station ",stattoo,": ",ds_info$statid[stattoo],sep=""))
  
    tmpstationCoords[,1]<-c(ds_info$statlats[stattoo],ds_info$statlats[corrneighbours.id[[stattoo]]])
    tmpstationCoords[,2]<-c(ds_info$statlons[stattoo],ds_info$statlons[corrneighbours.id[[stattoo]]])
    tmpstationElevs<-c(ds_info$statelvs[stattoo],ds_info$statelvs[corrneighbours.id[[stattoo]]])

    # Fit VAR and get autoregressive parameters and candidate station rwhichesiduals
    # a) calculate sample autocovariances for lags 0 to 1 - MATRIX MULTIPLICATION OF val_for_all_stations by val+for_all_stations(t+lag)
    # upsidedownLhat(lag0)=(1/ntimes)*SUM(val#val(t+0))
    # upsidedownLhat(lag1)=(1/ntimes)*SUM(val#val(t+1))
    # b) get actual residuals from the AR 1 model estimation 
    # Zt = At- TOTAL(phi*At-k) where k=1 for AR(1)  - so TOTAL is unnecessary in an AR(1) context
    # [[At.st1r1c1],[At.st2r2c1],[At.st3r3c1]] ## ([[At.st1r1c1],[At.st2r1c2],[At.st3r1c3]]) = 

    res=get_autocovs_func(c(CC,CCDiag),c(CCl1,CCl1Diag),tmpstationCoords,tmpstationElevs)  
    VARphi.vals[[stattoo]]<-res$phi_val     
    tmpVARacs<-res$newACs
    tmpGAM0s<-res$gamma_lag0
    tmpGAM1s<-res$gamma_lag1
    
    # WRITE TO FILE: VAR matrix for station 
    write.table(paste(c(ds_info$statid[stattoo],format(round(VARphi.vals[[stattoo]],digits=4),trim=FALSE,nsmall=4,width=8)),collapse=" "),
    	  file=paste(dirbawg,outfilarps,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpGAM0s,digits=4),trim=FALSE,nsmall=4,width=8)),collapse=" "),
    	  file=paste(dirbawg,outfilgam0s,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpGAM1s,digits=4),trim=FALSE,nsmall=4,width=8)),collapse=" "),
    	  file=paste(dirbawg,outfilgam1s,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
    #stop() # only using this with restarter
  
  } # end of 'for'
  
  # clean up
  rm (DistStatID,tmpstationCoords,tmpstationElevs,tmpGAM0s,tmpGAM1s)
  gc()

} else if (partone == 1) { # end of partone read in VARPS, station stats, neighbours list
 
  # read in corrneighbours.id
  print("Part One:1 - No need to read in anything")
#  # set up storage
#  # maps for locating stations within the main list ds_info
#  corrneighbours.id<-vector("list",nstations)	# map using full list
#
#  mooNEIGHS<-read.table(paste(dirbawg,outfilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
##  founds<-unname(as.matrix(mooNEIGHS[,-c(1,42,43,44)]))
#  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))
#  nNeighbours<-40	# 4 elements of candidate station ID and LONGEST COUNTS: 99
#  
#  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
#  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 
#  
#  rm(mooNEIGHS,moo,founds)
#  gc()

  VARphi.vals=vector("list",nstations)	# container for VAR autoregressive parameter matrices for each station
  
  # read in VARS and populate VAR_phi.vals

  print("Finding VARPS")

  mooVARS<-read.table(paste(dirbawg,outfilarps,sep=""),header=FALSE)
  moo<-unname(as.matrix(mooVARS[,-1]))
  moo<-apply(moo,1,function(x) list(matrix(x,(nNeighbours+1),(nNeighbours+1)))) 
  VARphi.vals<-unlist(moo,recursive=FALSE)
  
  rm(mooVARS,moo,founds)
  gc()

 
}
#stop()
#############################################################################
if (parttwo == 0) {
  print("Part Two:0")
  
  # set up storage
  # maps for locating stations within the main list ds_info
  stationnetworks.id<-array(replicate(nstations,0))	
  neighbournetworks.id<-vector("list",nnetworks)

  # NETWORK parameters
  CORRarr<-vector("list",nnetworks)	# list of arrays containing covariance matrix for each network		# 

  # NEW parameters
  newRESarr<-array(mdi,dim=c(nstations,nmodmons))# array for new MVN station residuals (in order) for each station (can be removed after single use)
  tmpOLDRESarr<-array(mdi,nmodmons)		# temporary array for the old residuals for a station
  tmpresidsCLIM<-0
  tmpresidsSD<-0

  # set up station locs for Resids factorisation - use FULL list
  site.coords<-matrix(c(ds_info$statlats,ds_info$statlons),ncol=2,nrow=nstations)	# read in lons and lats of all stations
  
  # 11) Read in network boundaries and loop through
  for (nN in 1:nnetworks) {
    print(paste("Network:",nN,sep=" "))
    mush<-readLines(con=paste(dirlist,inNTWKstats,nN,'.dat',sep=""),n=-1)	# read entire file
    for (linoo in 1:length(mush)) {
      stationnetworks.id[type.convert(substr(mush[linoo],1,5))]<-nN			# station lat
    }
    mush<-readLines(con=paste(dirlist,inNTWKnbs,nN,'.dat',sep=""),n=-1)	# read entire file
    for (linoo in 1:length(mush)) {
      if (linoo == 1) {
        neighbournetworks.id[[nN]]<-type.convert(substr(mush[linoo],1,5))
      } else {
        neighbournetworks.id[[nN]]<-append(neighbournetworks.id[[nN]],type.convert(substr(mush[linoo],1,5)))
      }
    }

    # Get covariance matrix of the residuals
    # PROBLEM - the residuals are created using different networks by necessity - very small ones that do not garantee neighbours
    # It is HIGHLY unlikely that we will find a common period over which all stations within the network and its neighbours can be assessed
    # In fact it would be easier/quicker to just use distance as a proxy for now!!!!
    # TEST COR ON SOME RESIDUALS TO SEE HOW DEGRADED THE RELATIONSHIP IS FROM EXPONENTIAL

    wanted.sites <- ((1:nstations) %in% neighbournetworks.id[[nN]]) | (stationnetworks.id == nN)

# NOW FAFF AROUND DOING IT PROPERLY USING GAMMAs AND PHIs
# g0 using CC
# g1 using CCl1
# phi using t(g1)%*%solve(g0)
# cov(Z) using g0 - phi%*%t(g1)
    tmpgam0<-expdh.corr(CC,site.coords[wanted.sites,],ds_info$statelvs[wanted.sites]/1000.,CCDiag)
    tmpgam1<-expdh.corr(CCl1,site.coords[wanted.sites,],ds_info$statelvs[wanted.sites]/1000.,CCl1Diag)
    tmpphi<-t(tmpgam1)%*%solve(tmpgam0) 
    CORRarr[[nN]]<-round(tmpgam0 - tmpphi%*%tmpgam1,digits=3)
    # Now need some tweaking:
    # force diags to be 1.0
    # constrain to 3 significant digits to try and prevent floating point sillyness
    # check for positive definiteness using eigenvalues - if some are negative, rebuild as positive
    # SHOULD BE OK AS THESE ARE THE TINIEST EIGENVALUES
    diag(CORRarr[[nN]])<-1.0

  } # end of 'for' in 11)

  # 12) Resample the residuals to get new ones using the Factorise method
  cat("Setting up neighbourhoods and covariance structures ...\n")
  test.setup <- rbigmvn.setup(mu=1:nstations,	# must have a mean of zero so subtract later
                            coords=site.coords,		# lat, lon coordinates
                            groups=stationnetworks.id,	# stations represented by their network number id
			    neighbours=neighbournetworks.id,	# a list of neighbours for each network by station listing number
                            coord.type="geographical",
			    method="factorize",
                            Sigmalist=CORRarr)		# list of covariance matrices for each network

  # Run - returns a time (rows) by station (column) array
  newRESarr<- rbigmvnorm(nmodmons,setup=test.setup,nburnin=100,monitor=FALSE)

  print(paste(c(mean(newRESarr[,100]),newRESarr[1:5,100]),sep=" "))

  # Transform to station (rows) and time (columns) and subtract numbers 1:nstations from each row
  # row wise (1) apply needs transforming after, rowwise (2) does not
  newRESarr<-apply(t(newRESarr),2,function(x) x-c(1:nstations))	# should now be mean zero

  print(paste(c(mean(newRESarr[100,]),newRESarr[100,1:5]),sep=" "))

  # force standardise for newRESarr as they are a little skewed
  newRESarr<-t(apply(newRESarr,1,scale))	# should now be mean zero

  print(paste(c(mean(newRESarr[100,]),newRESarr[100,1:5]),sep=" "))
 
  # 13) Loop through each station and transform residuals back to what they should be in reality 
  countreals<-1
  
  mooRES<-scan(paste(dirbawg,outfilOLDshocks,sep=""),what="list",sep="\n")	# each line is now a character string to be split
  
  for (stattoo in 1:nstations) {  # NEED TO DO THIS BECAUSE FULL STATION LIST IS USED TO BUILD RESIDS NETWORK
    
    # use line from file of old residuals (ignore stat ID, , REAL/NORM, mean and sd, standardise
    # read in residsCLIM and residsSD (no locs for faux stations - 0, 1 are assumed)
    moo<-strsplit(mooRES[stattoo],split=" ")       # silly format doesn't work here so need some trickery
    goo<-moo[[1]][moo[[1]] != ""]
    tmpresidsCLIM<-type.convert(goo[3])
    tmpresidsSD<-type.convert(goo[4])
    tmpOLDRESarr<-type.convert(goo[5:length(goo)])
    tmpOLDRESarr<-(tmpOLDRESarr - tmpresidsCLIM)/tmpresidsSD # should now be mean zero/1 sd
  
    # transform to new residuals, unstandardise
    res=QMapTransform_func(newRESarr[stattoo,],tmpOLDRESarr)
    newRESarr[stattoo,]<-res$newfirstdist
    # This doesn't seem to be coming out as MVN so need to standardise CHECK THIS!!!
    newRESarr[stattoo,]<-(newRESarr[stattoo,] - mean(newRESarr[stattoo,]))/sd(newRESarr[stattoo,])
    # now adjust to mean and sd of real data 
    newRESarr[stattoo,]<-(newRESarr[stattoo,] * tmpresidsSD) + tmpresidsCLIM

    # WRITE TO FILE: residuals original (+ clim and sd)
    write.table(paste(c(ds_info$statid[stattoo],format(round(mean(newRESarr[stattoo,]),digits=3),trim=FALSE,nsmall=3,width=6),
                format(round(sd(newRESarr[stattoo,]),digits=3),trim=FALSE,nsmall=3,width=6),
		format(round(newRESarr[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
                file=paste(dirbawg,outfilNEWshocks,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
#    stop()
    
  } # end of 'for' in 13

  #stop()
  # CLEAN UP
  rm(mooRES,moo,goo,mush,stationnetworks.id,neighbournetworks.id,CORRarr,tmpOLDRESarr,
     wanted.sites,tmpresidsCLIM,tmpresidsSD,site.coords)
  gc()

} else if (parttwo == 1) { # end of part two 
  print("Part Two:1")
  
  print("Part Two:1 reading in NEIGHBOUR LISTS AND VARPS")
  # set up storage
  # maps for locating stations within the main list ds_info
  corrneighbours.id<-vector("list",nstations)	# map using full list
  VARphi.vals=vector("list",nstations)	# container for VAR autoregressive parameter matrices for each station
  newRESarr<-array(mdi,dim=c(nstations,nmodmons))# array for new MVN station residuals (in order) for each station (can be removed after single use)
  
  # read in correlating neighbours and populate corrneighbours.id
  # read in VARS and populate VAR_phi.vals
  # use actual station list here as opposed to full station list

  mooNEIGHS<-read.table(paste(dirbawg,outfilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
#  founds<-unname(as.matrix(mooNEIGHS[,-c(1,42,43,44)]))
  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))
  nNeighbours<-40	# 4 elements of candidate station ID and LONGEST COUNTS: 99

  print("Finding Corrneighbours")
  
  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 

  print("Finding VARPS")

  mooVARS<-read.table(paste(dirbawg,outfilarps,sep=""),header=FALSE)
  moo<-unname(as.matrix(mooVARS[,-1]))
  moo<-apply(moo,1,function(x) list(matrix(x,(nNeighbours+1),(nNeighbours+1)))) 
  VARphi.vals<-unlist(moo,recursive=FALSE)
  
  rm(mooNEIGHS,mooVARS,moo,founds)
  gc()

  # read in new shocks for full list
  print("Part Two:1 - reading in NEW RANDOM SHOCKS")
  mooRES<-read.table(paste(dirbawg,outfilNEWshocks,sep=""),header=FALSE)	# each line is now a character string to be split
  newRESarr<-unname(as.matrix(mooRES[,-c(1,2,3)])) # not first three elements
  rm(mooRES)
  gc()
} 
#stop()
################################################################################################

if (partthree == 0) { 
  print("Part Three:0 - Building Clean SDANOMS")
  
  # 14) Now, for each time point, for each station apply VAR, add residual - spin up for 10, begin
  # need both FullWorld and CleanWorld - use FullWorld for builing the TmpWorld and then sample over to CleanWorld for real stations
  
  #Initialise world with spin up
  FullWorldPRIOR<-array(newRESarr[,10],c(nstations))
  FullWorldPOST<-array(0,c(nstations))
  CleanWorldSDANOM<-array(NA,c(nstations,nmodmons)) # array for new data for each station - std anoms, 

  for (tt in 2:10) {	# MAY NEED A LONGER SPIN UP TO MAKE SURE THERE IS NOTHING WEIRD
    for (stattoo in 1:nstations) {
      # running residuals 10,9,,8....to 1 ready for lift off
      # Need a clever way of not overwriting the REAL stations - use tmp and fauxneighbours.id
      tmp<-(VARphi.vals[[stattoo]] %*% FullWorldPRIOR[c(stattoo,corrneighbours.id[[stattoo]])]) + 
				      newRESarr[c(stattoo,corrneighbours.id[[stattoo]]),(11-tt)]
      FullWorldPOST[stattoo]<-tmp[1]
    }
    FullWorldPRIOR<-array(FullWorldPOST,(nstations))
  } 

  # reset initial value to end of spin up - mapping across only the real stations
  CleanWorldSDANOM[,1]<-FullWorldPRIOR	# tt continues to exist outside of loop

  # Make a CleanWorld for real
  for (tt in 2:nmodmons) {
    for (stattoo in 1:nstations) {
      # Need a clever way of not overwriting the REAL stations - use tmp and fauxneighbours.id
      tmp<-(VARphi.vals[[stattoo]] %*% FullWorldPRIOR[c(stattoo,corrneighbours.id[[stattoo]])]) + 
					newRESarr[c(stattoo,corrneighbours.id[[stattoo]]),tt]
       FullWorldPOST[stattoo]<-tmp[1]
    }
    FullWorldPRIOR<-FullWorldPOST
    CleanWorldSDANOM[,tt]<-FullWorldPOST
  } # end of 'for' in 14 

  # WRITE TO FILE: NEW SDANOMS
  for (stattoo in 1:nstations) {

    # WRITE TO FILE: new SDANOMS
    write.table(paste(c(ds_info$statid[stattoo],format(round(CleanWorldSDANOM[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilNEWsdanoms,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
  }  
#  plot(CleanWorldSDANOM[stattoo,])
#  stop()
# CLEAN UP
  rm(FullWorldPOST,FullWorldPRIOR,VARphi.vals,newRESarr)
  gc()
} else if (partthree == 1) {	# end of part three
  print("Part Three:1 -  reading in NEIGHBOUR LISTS AND NEWSDS")
  # set up storage
  # maps for locating stations within the main list ds_info
  corrneighbours.id<-vector("list",nstations)	# map using full list
  CleanWorldSDANOM<-array(NA,c(nstations,nmodmons)) # array for new data for each station - std anoms, 

  # read in correlating neighbours and populate corrneighbours.id

  mooCLEAN<-read.table(paste(dirbawg,outfilNEWsdanoms,sep=""),header=FALSE)	# each line is now a character string to be split

  mooNEIGHS<-read.table(paste(dirbawg,outfilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
#  founds<-unname(as.matrix(mooNEIGHS[,-c(1,42,43,44)]))
  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))
  nNeighbours<-40	# 4 elements of candidate station ID and LONGEST COUNTS: 99

  print("Finding Corrneighbours")
  
  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 

  # read in clean world data
  print("Finding Clean SDANOMS")
  
  CleanWorldSDANOM<-unname(as.matrix(mooCLEAN[,-1]))
  
  rm(mooNEIGHS,mooCLEAN,founds,moo)
  gc()
}
#stop()
##################################################################################################
if (partfour == 0) {
  print("Part Four:0 - Building Clean CLEANANOMS, clean monthlies and diffs")

  # Set up storage
  tmpACsdanom<-0
  tmpACclimanom<-0
  tmpdiffos<-0	# 2d array of station - neighbour difference series
  tmpdiffsds<-0   # array of difference series standard deviations
  tmpdiffcorr<-0	# array of difference series lag 1 autocorrelation
  FullDiffSD<-0  	# growing array to save all diffs SDs to plot as hist (REMOVE 1st ELEMENT)
  FullDiffCORR<-0	# growing array to save all diffs lag 1 cor to plot as hist (REMOVE 1st ELEMENT)
  CleanWorldCLIMANOM<-array(0,c(nstations,nmodmons)) # array for new data for each station - clim anoms with loess, 
  CleanWorld<-array(0,c(nstations,nmodmons)) # array for new data for each station - monthly means, 
  tmpclims<-0
  tmpsds<-0
  FullSD<-0
  FullCLIM<-0
  
  # 16) Now loop through each station, comp with neighbours AND
  # calc NEW SDANOM ac, multiply by sd, add back mean, add in smoothed GCM curve calc new CLIMANOM AC
  mooSTATS<-read.table(paste(dirbawg,outfilOLDstats,sep=""),header=FALSE)	# each line is now a character string to be split
  mooGCM<-read.table(paste(dirbawg,infilGCMloess,sep=""),header=FALSE)	# each line is now a character string to be split
# NOTE: Extending GCM back beyond 1860 by flipping the loess

  for (stattoo in 1:nstations) {

    if (restarter != "--------" & restarter != ds_info$statid[stattoo]) {
      next
    }
    restarter<-"--------"

    # get set of station minus neighbours
    tmpdiffos<-t(apply(CleanWorldSDANOM[corrneighbours.id[[stattoo]],],1,function(x) CleanWorldSDANOM[stattoo,]-x))
    tmpdiffsds<-apply(tmpdiffos,1,sd)
    tmpdiffcors<-apply(tmpdiffos,1,function(x) cor(x[1:(nmodmons-1)],x[2:nmodmons]))
    
    FullDiffSD<-c(FullDiffSD,tmpdiffsds)
    FullDiffCORR<-c(FullDiffCORR,tmpdiffcors)
    
    # WRITE TO FILE DIFF SERIES SDS AND COR
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpdiffsds,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilDIFFSDstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpdiffcors,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilDIFFACstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
          
    # Get NEW SDANOM AC
    moo<-acf(CleanWorldSDANOM[stattoo,],lag.max=24,plot=FALSE,na.action=na.pass)
    tmpACsdanom<-moo$acf[2:25,1,1]
  
    # Read in clims and SD (new file format - JUL2015)
    tmpclims<-unname(as.numeric(mooSTATS[stattoo,3:14]))
    tmpsds<-unname(as.numeric(mooSTATS[stattoo,16:27]))
    
    # Multiply by SD 
    CleanWorldCLIMANOM[stattoo,]<-array(apply(array(CleanWorldSDANOM[stattoo,],dim=c(12,nmodmons/12)),2,
                                  function(x) x*tmpsds),dim=(nmodmons))

    # Recalculate SD
    tmpsds<-apply(array(CleanWorldCLIMANOM[stattoo,],dim=c(12,nmodmons/12)),1,sd,na.rm=TRUE)
    FullSD<-sd(CleanWorldCLIMANOM[stattoo,],na.rm=TRUE)
  
    # Read in the GCM LOESS for this station
    # if start year is lower than 1860 (RealModStYr) then flip to fill ***
    tmploess<-unname(as.numeric(mooGCM[stattoo,-1]))
    if (RealModStYr != modstyr) {
      #extend tmploess
      tmploess<-c(rev(tmploess[1:ModExt]),tmploess) 
    }  
    CleanWorldCLIMANOM[stattoo,]<-CleanWorldCLIMANOM[stattoo,]+tmploess
    
    # Get NEW CLIMANOM AC
    moo<-acf(CleanWorldCLIMANOM[stattoo,],lag.max=24,plot=FALSE,na.action=na.pass)
    tmpACclimanom<-moo$acf[2:25,1,1]

    # WRITE TO FILE: NEW Clim anoms
    write.table(paste(c(ds_info$statid[stattoo],format(round(CleanWorldCLIMANOM[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilNEWcanoms,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
  
    # add back clim
    CleanWorld[stattoo,]<-array(apply(array(CleanWorldCLIMANOM[stattoo,],dim=c(12,nmodmons/12)),2,
                                  function(x) x+tmpclims),dim=(nmodmons))
     
    # REcalculate CLIMS and ANOMS
    tmpclims<-apply(array(CleanWorld[stattoo,],dim=c(12,nmodmons/12)),1,mean,na.rm=TRUE)
    FullCLIM<-mean(tmpclims,na.rm=TRUE)
    
    # Calculate trends
    trendy<-120.*coef(lm(CleanWorldCLIMANOM[stattoo,]~seq(nmodmons)))[2]	# this is the decadal rate of change linear trend, should be omitting NAs - not perfect but sufficient    

    # WRITE TO FILE: NEW stations Clims, sds, clim AC, sdanom AC anoms
    write.table(paste(c(ds_info$statid[stattoo],format(round(trendy,digits=2),trim=FALSE,nsmall=2,width=6),
              format(round(FullCLIM,digits=2),trim=FALSE,nsmall=2,width=6),
              format(round(tmpclims,digits=2),trim=FALSE,nsmall=2,width=6),
              format(round(FullSD,digits=2),trim=FALSE,nsmall=2,width=6),
              format(round(tmpsds,digits=2),trim=FALSE,nsmall=2,width=6),
	      format(round(tmpACclimanom[1],digits=3),trim=FALSE,nsmall=3,width=6),
	      format(round(tmpACsdanom[1],digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
              file=paste(dirbawg,outfilNEWstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    # WRITE TO FILE: NEW stations 
    write.table(paste(c(ds_info$statid[stattoo],
                format(round(CleanWorld[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirdata,outfildata,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

  } # end of 'for' in 16)
  rm(tmpACsdanom,tmpACclimanom,tmpdiffos,tmpdiffsds,tmpdiffcorr,FullDiffSD,FullDiffCORR,
     CleanWorld,tmpclims,tmpsds,FullSD,FullCLIM,mooSTATS,mooGCM)
  gc()
  
  tmpdiffos<-0	# 2d array of station - neighbour difference series
  tmpdiffsds<-0   # array of difference series standard deviations
  tmpdiffcorr<-0	# array of difference series lag 1 autocorrelation
  FullDiffSD<-0  	# growing array to save all diffs SDs to plot as hist (REMOVE 1st ELEMENT)
  FullDiffCORR<-0	# growing array to save all diffs lag 1 cor to plot as hist (REMOVE 1st ELEMENT)

  for (stattoo in 1:nstations) {
    # get set of station minus neighbours
    tmpdiffos<-t(apply(CleanWorldCLIMANOM[corrneighbours.id[[stattoo]],],1,function(x) CleanWorldCLIMANOM[stattoo,]-x))
    tmpdiffsds<-apply(tmpdiffos,1,sd)
    tmpdiffcors<-apply(tmpdiffos,1,function(x) cor(x[1:(nmodmons-1)],x[2:nmodmons]))
    
    FullDiffSD<-c(FullDiffSD,tmpdiffsds)
    FullDiffCORR<-c(FullDiffCORR,tmpdiffcors)
    
    # WRITE TO FILE DIFF SERIES SDS AND COR
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpdiffsds,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilDIFFSDCMstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpdiffcors,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilDIFFACCMstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

  } # end of 'for' in 16)
  rm( tmpdiffos,tmpdiffsds,tmpdiffcorr,FullDiffSD,FullDiffCORR)
  gc()
} else if (partfour == 1) { # has fallen over in/just before clim anoms bit so restart here (perhaps because 'restarter' has been set?
  
  print("Part Four:1 - Reading in CleanClimAnoms and Neighbours")

  # set up storage
  # maps for locating stations within the main list ds_info
  corrneighbours.id<-vector("list",nstations)	# map using full list
  CleanWorldCLIMANOM<-array(0,c(nstations,nmodmons)) # array for new data for each station - clim anoms with loess, 
  
  mooNEIGHS<-read.table(paste(dirbawg,outfilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
#  founds<-unname(as.matrix(mooNEIGHS[,-c(1,42,43,44)]))
  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))
  nNeighbours<-40	# 4 elements of candidate station ID and LONGEST COUNTS: 99

  print("Finding Corrneighbours")
  
  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 
  
  rm(mooNEIGHS,founds,moo)
  gc()

  # read in NEW clim anoms
  mooCLEAN<-read.table(paste(dirbawg,outfilNEWcanoms,sep=""),header=FALSE)	# each line is now a character string to be split
  CleanWorldCLIMANOM<-unname(as.matrix(mooCLEAN[,-1]))
  rm(mooCLEAN)
  gc()
  
  tmpdiffos<-0	# 2d array of station - neighbour difference series
  tmpdiffsds<-0   # array of difference series standard deviations
  tmpdiffcorr<-0	# array of difference series lag 1 autocorrelation
  FullDiffSD<-0  	# growing array to save all diffs SDs to plot as hist (REMOVE 1st ELEMENT)
  FullDiffCORR<-0	# growing array to save all diffs lag 1 cor to plot as hist (REMOVE 1st ELEMENT)

  for (stattoo in 1:nstations) {
    # get set of station minus neighbours
    tmpdiffos<-t(apply(CleanWorldCLIMANOM[corrneighbours.id[[stattoo]],],1,function(x) CleanWorldCLIMANOM[stattoo,]-x))
    tmpdiffsds<-apply(tmpdiffos,1,sd)
    tmpdiffcors<-apply(tmpdiffos,1,function(x) cor(x[1:(nmodmons-1)],x[2:nmodmons]))
    
    FullDiffSD<-c(FullDiffSD,tmpdiffsds)
    FullDiffCORR<-c(FullDiffCORR,tmpdiffcors)
    
    # WRITE TO FILE DIFF SERIES SDS AND COR
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpdiffsds,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilDIFFSDCMstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpdiffcors,digits=2),trim=FALSE,nsmall=2,width=6)),collapse=" "),
              file=paste(dirbawg,outfilDIFFACCMstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
  } # end of 'for'
  rm(FullDiffSD,FullDiffCORR,tmpdiffos,tmpdiffsds,tmpdiffcors)
  gc()
  
  print("Now reading in Standardised Anoms too ready for part 5")
  # read in clean world data
  mooCLEAN<-read.table(paste(dirbawg,outfilNEWsdanoms,sep=""),header=FALSE)	# each line is now a character string to be split
  print("Finding Clean SDANOMS")
  
  CleanWorldSDANOM<-unname(as.matrix(mooCLEAN[,-1]))
  rm(mooCLEAN)
  gc()
  
}  else if (partfour == 2) { # end of part four
  print("Part Four:2 - read in NEW clim anoms, sd anoms and corrneighbours")
  
  # read in NEW clim anoms, new sd anoms and neighbour lists
  mooCLEAN<-read.table(paste(dirbawg,outfilNEWcanoms,sep=""),header=FALSE)	# each line is now a character string to be split
  CleanWorldCLIMANOM<-unname(as.matrix(mooCLEAN[,-1]))

  mooCLEAN<-read.table(paste(dirbawg,outfilNEWsdanoms,sep=""),header=FALSE)	# each line is now a character string to be split
  CleanWorldSDANOM<-unname(as.matrix(mooCLEAN[,-1]))
  
  rm(mooCLEAN)
  gc()

  mooNEIGHS<-read.table(paste(dirbawg,outfilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
#  founds<-unname(as.matrix(mooNEIGHS[,-c(1,42,43,44)]))
  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))
  nNeighbours<-40	# 4 elements of candidate station ID and LONGEST COUNTS: 99

  print("Finding Corrneighbours")
  
  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 

  rm(mooNEIGHS,moo,founds)
  gc()

}
#stop()
#######################################################################################

if (partfive == 0) {
  print("Part Five 0 - NEW covs")
  # 17) FINALLY - LOOP THROUGH EACH STATIOn, READ IN CORR NEIGHBOURS AND LOOK AT COVARIANCES OF CLIMANOMS AND SDNAOMS
   # set up full arrays
  fullCCclimanom<-0
  fullCCsdanom<-0
  fullCCl1climanom<-0
  fullCCl1sdanom<-0
  tmpCCclimanom<-array(NA,40)
  tmpCCsdanom<-array(NA,40)
  tmpCCl1climanom<-array(NA,40)
  tmpCCl1sdanom<-array(NA,40)
    
  for (stattoo in 1:nstations) {
    print(ds_info$statid[stattoo])
    nNs<-length(corrneighbours.id[[stattoo]])
    
    # set up arrays
    tmpCCclimanom<-array(NA,nNs)
    tmpCCsdanom<-array(NA,nNs)

    tmpCCl1climanom<-array(NA,nNs)
    tmpCCl1sdanom<-array(NA,nNs)
    
    # loop through corrneighbours.id (miss out fauxneighbours.id as these have not been processed)
    for (nN in 1:nNs) {
      
      # Calculate cor with candidate CLIM ANOM
      tmpCCclimanom[nN]<-cor(CleanWorldCLIMANOM[stattoo,],CleanWorldCLIMANOM[corrneighbours.id[[stattoo]][nN],])
      fullCCclimanom<-c(fullCCclimanom,tmpCCclimanom[nN])
      # Calculate cor with candidate SD ANOM
      tmpCCsdanom[nN]<-cor(CleanWorldSDANOM[stattoo,],CleanWorldSDANOM[corrneighbours.id[[stattoo]][nN],])
      fullCCsdanom<-c(fullCCsdanom,tmpCCsdanom[nN])

      # Calculate lag 1 cor with candidate CLIM ANOM
      tmpCCl1climanom[nN]<-cor(CleanWorldCLIMANOM[stattoo,1:nmodmons-1],CleanWorldCLIMANOM[corrneighbours.id[[stattoo]][nN],2:nmodmons])
      fullCCl1climanom<-c(fullCCl1climanom,tmpCCl1climanom[nN])
      # Calculate lag 1 cor with candidate SD ANOM
      tmpCCl1sdanom[nN]<-cor(CleanWorldSDANOM[stattoo,1:nmodmons-1],CleanWorldSDANOM[corrneighbours.id[[stattoo]][nN],2:nmodmons])
      fullCCl1sdanom<-c(fullCCl1sdanom,tmpCCl1sdanom[nN])
    }
    
    # WRITE TO FILE: NEW corrs
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpCCclimanom,digits=3),trim=FALSE,nsmall=3,width=6),
                format(round(tmpCCsdanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilNEWcovs,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    # WRITE TO FILE: NEW lag 1 corrs
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpCCl1climanom,digits=3),trim=FALSE,nsmall=3,width=6),
                format(round(tmpCCl1sdanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilNEWl1covs,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
  # retain all of these for making a scatter plot afterwards

  }
  rm(fullCCclimanom,fullCCsdanom,fullCCl1climanom,fullCCl1sdanom,tmpCCclimanom,tmpCCsdanom,
  tmpCCl1climanom,tmpCCl1sdanom)
  gc()
  #stop()  
} else {

} 

#######################################################################################
if (partsix == 0) {
  print("Part Six 0 - OLD covs")
  # 17) FINALLY - LOOP THROUGH EACH STATIOn, READ IN CORR NEIGHBOURS AND LOOK AT COVARIANCES OF CLIMANOMS AND SDNAOMS
  # set up full arrays
  fullCCclimanom<-0
  fullCCsdanom<-0
  fullCCl1climanom<-0
  fullCCl1sdanom<-0
  tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
  tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
  tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
  tmpCCclimanom<-array(NA,40)
  tmpCCsdanom<-array(NA,40)
  tmpCCl1climanom<-array(NA,40)
  tmpCCl1sdanom<-array(NA,40)
  tmpNABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
  tmpNCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
  tmpNSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
  
  for (stattoo in 1:nstations) {

    if (restarter != "--------" & restarter != ds_info$statid[stattoo]) {
      next
    }
    restarter<-"--------"

    # reset all tmp arrays
    tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
    tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
    tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms

    print(paste("Working on station ",stattoo,": ",ds_info$statid[stattoo],sep=""))
  
    # read in the candidate station and filter out data to correct time points
    tmpABSstation=read_station_func(styr,paste(dirdata,infilraw,ds_info$statid[stattoo],"_stage3",sep=""),tmpABSstation)
  
    # INTERPOLATE OVER MISSING VALUES WHERE THEY ARE ISOLATED CASES OF 1-2 MAXIMUM CONSECUTIVE MISSING DATA
    tmpABSstation[which(tmpABSstation == mdi)]<-NA   # set missing values to R standard mdi
    tmpABSstation=interpolate_missing_func(tmpABSstation) 	

    # CREATE STATION STATS: clims, st devs, anoms, residuals from a smoothed filter
    # method CLS - calculate Climatology and remove, fit Loess and remove, calculated Standard devations and divide by them.
    res=get_climanoms_func(tmpABSstation,spanval)
    tmpCLIMANOMstation<-res$anomsarr
    tmpSDANOMstation<-res$stanomsarr
       
    nNs<-length(corrneighbours.id[[stattoo]])
   
    # set up arrays
    tmpCCclimanom<-array(NA,nNs)
    tmpCCsdanom<-array(NA,nNs)

    tmpCCl1climanom<-array(NA,nNs)
    tmpCCl1sdanom<-array(NA,nNs)
    
    # loop through corrneighbours.id (miss out fauxneighbours.id as these have not been processed)
    for (nN in 1:nNs) {
      tmpNABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
      tmpNCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
      tmpNSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
      
      # read in the candidate station and filter out data to correct time points
      tmpNABSstation=read_station_func(styr,paste(dirdata,infilraw,ds_info$statid[corrneighbours.id[[stattoo]][nN]],"_stage3",sep=""),tmpNABSstation)
  
      # INTERPOLATE OVER MISSING VALUES WHERE THEY ARE ISOLATED CASES OF 1-2 MAXIMUM CONSECUTIVE MISSING DATA
      tmpNABSstation[which(tmpNABSstation == mdi)]<-NA   # set missing values to R standard mdi
      tmpNABSstation=interpolate_missing_func(tmpNABSstation) 	

      # CREATE STATION STATS: clims, st devs, anoms, residuals from a smoothed filter
      # method CLS - calculate Climatology and remove, fit Loess and remove, calculated Standard devations and divide by them.
      res=get_climanoms_func(tmpNABSstation,spanval)
      tmpNCLIMANOMstation<-res$anomsarr
      tmpNSDANOMstation<-res$stanomsarr

      # Match up to see if there are enough points in common
      gots<-which(!is.na(tmpCLIMANOMstation) & !is.na(tmpNCLIMANOMstation))
      if (length(gots) > 60) { # there are more than 60 months (5 years) of data in common
      
        # Calculate cor with candidate CLIM ANOM
        tmpCCclimanom[nN]<-cor(tmpCLIMANOMstation[gots],tmpNCLIMANOMstation[gots])
        fullCCclimanom<-c(fullCCclimanom,tmpCCclimanom[nN])
        # Calculate cor with candidate SD ANOM
        tmpCCsdanom[nN]<-cor(tmpSDANOMstation[gots],tmpNSDANOMstation[gots])
        fullCCsdanom<-c(fullCCsdanom,tmpCCsdanom[nN])

        # Calculate lag 1 cor with candidate CLIM ANOM
        tmpCCl1climanom[nN]<-cor(tmpCLIMANOMstation[gots[1:length(gots)-1]],tmpNCLIMANOMstation[gots[2:length(gots)]])
        fullCCl1climanom<-c(fullCCl1climanom,tmpCCl1climanom[nN])
        # Calculate lag 1 cor with candidate SD ANOM
        tmpCCl1sdanom[nN]<-cor(tmpSDANOMstation[gots[1:length(gots)-1]],tmpNSDANOMstation[gots[2:length(gots)]])
        fullCCl1sdanom<-c(fullCCl1sdanom,tmpCCl1sdanom[nN])
      }
    }
    
    # WRITE TO FILE: NEW corrs
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpCCclimanom,digits=3),trim=FALSE,nsmall=3,width=6),
               format(round(tmpCCsdanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilOLDcovs,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    # WRITE TO FILE: NEW lag 1 corrs
    write.table(paste(c(ds_info$statid[stattoo],format(round(tmpCCl1climanom,digits=3),trim=FALSE,nsmall=3,width=6),
                format(round(tmpCCl1sdanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilOLDl1covs,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
  # retain all of these for making a scatter plot afterwards

  }
  rm(fullCCclimanom,fullCCsdanom,fullCCl1climanom,fullCCl1sdanom,tmpABSstation,tmpCLIMANOMstation,
     tmpSDANOMstation,tmpCCclimanom,tmpCCsdanom,tmpCCl1climanom,tmpCCl1sdanom,tmpNABSstation,tmpNCLIMANOMstation,tmpNSDANOMstation,)
  gc()
  #stop()  
} else { # JUST READ IN AND PLOT

}
#######################################################################################
if (partseven == 0) {
  print("Part Seven 0 - OLD diffs")
  # 17) FINALLY - LOOP THROUGH EACH STATIOn, READ IN CORR NEIGHBOURS AND LOOK AT COVARIANCES OF CLIMANOMS AND SDNAOMS
  # set up full arrays
  fulldiffsdclimanom<-0
  fulldiffsdsdanom<-0
  fulldiffacclimanom<-0
  fulldiffacsdanom<-0
  tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
  tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
  tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
  tmpdiffsdclimanom<-array(NA,40)
  tmpdiffsdsdanom<-array(NA,40)
  tmpdiffacclimanom<-array(NA,40)
  tmpdiffacsdanom<-array(NA,40)
  tmpNABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs  
  tmpNCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
  tmpNSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms

  # set up storage
  # maps for locating stations within the main list ds_info
  corrneighbours.id<-vector("list",nstations)	# map using full list

  mooNEIGHS<-read.table(paste(dirbawg,outfilNEIGHlist,sep=""),header=FALSE,colClasses="character")	# each line is now a character string to be split
#  founds<-unname(as.matrix(mooNEIGHS[,-c(1,42,43,44)]))
  founds<-unname(as.matrix(mooNEIGHS[,-c(1)]))
  
  nNeighbours<-40	# 4 elements of candidate station ID and LONGEST COUNTS: 99

  print("Finding Corrneighbours")
  
  moo<-apply(founds,1,function(x) list(match(x,ds_info$statid))) 
  corrneighbours.id<-lapply(moo,function(x) unlist(x)) 
  
  rm(mooNEIGHS,founds)
  gc()
      
  for (stattoo in 1:nstations) {

    if (restarter != "--------" & restarter != ds_info$statid[stattoo]) {
      next
    }
    restarter<-"--------"

    # reset all tmp arrays
    tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
    tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
    tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms

    print(paste("Working on station ",stattoo,": ",ds_info$statid[stattoo],sep=""))
  
    # read in the candidate station and filter out data to correct time points
    tmpABSstation=read_station_func(styr,paste(dirdata,infilraw,ds_info$statid[stattoo],"_stage3",sep=""),tmpABSstation)
  
    # INTERPOLATE OVER MISSING VALUES WHERE THEY ARE ISOLATED CASES OF 1-2 MAXIMUM CONSECUTIVE MISSING DATA
    tmpABSstation[which(tmpABSstation == mdi)]<-NA   # set missing values to R standard mdi
    tmpABSstation=interpolate_missing_func(tmpABSstation) 	

    # CREATE STATION STATS: clims, st devs, anoms, residuals from a smoothed filter
    # method CLS - calculate Climatology and remove, fit Loess and remove, calculated Standard devations and divide by them.
    res=get_climanoms_func(tmpABSstation,spanval)
    tmpCLIMANOMstation<-res$anomsarr
    tmpSDANOMstation<-res$stanomsarr
       
    nNs<-length(corrneighbours.id[[stattoo]])
   
    # set up arrays
    tmpdiffsdclimanom<-array(NA,nNs)
    tmpdiffsdsdanom<-array(NA,nNs)

    tmpdiffacclimanom<-array(NA,nNs)
    tmpdiffacsdanom<-array(NA,nNs)
    
    # loop through corrneighbours.id (miss out fauxneighbours.id as these have not been processed)
    for (nN in 1:nNs) {
      tmpNABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
      tmpNCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
      tmpNSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
      
      # read in the candidate station and filter out data to correct time points
      tmpNABSstation=read_station_func(styr,paste(dirdata,infilraw,ds_info$statid[corrneighbours.id[[stattoo]][nN]],"_stage3",sep=""),tmpNABSstation)
  
      # INTERPOLATE OVER MISSING VALUES WHERE THEY ARE ISOLATED CASES OF 1-2 MAXIMUM CONSECUTIVE MISSING DATA
      tmpNABSstation[which(tmpNABSstation == mdi)]<-NA   # set missing values to R standard mdi
      tmpNABSstation=interpolate_missing_func(tmpNABSstation) 	

      # CREATE STATION STATS: clims, st devs, anoms, residuals from a smoothed filter
      # method CLS - calculate Climatology and remove, fit Loess and remove, calculated Standard devations and divide by them.
      res=get_climanoms_func(tmpNABSstation,spanval)
      tmpNCLIMANOMstation<-res$anomsarr
      tmpNSDANOMstation<-res$stanomsarr

      # Match up to see if there are enough points in common
      gots<-which(!is.na(tmpCLIMANOMstation) & !is.na(tmpNCLIMANOMstation))
      if (length(gots) > 60) { # there are more than 60 months (5 years) of data in common
      
        # Calculate cor with candidate CLIM ANOM
        if (sum(tmpCLIMANOMstation[gots[1:length(gots)-1]]-tmpNCLIMANOMstation[gots[1:length(gots)-1]]) != 0.0) {
          tmpdiffsdclimanom[nN]<-sd(tmpCLIMANOMstation[gots]-tmpNCLIMANOMstation[gots])
          fulldiffsdclimanom<-c(fulldiffsdclimanom,tmpdiffsdclimanom[nN])
          # Calculate cor with candidate SD ANOM
          tmpdiffsdsdanom[nN]<-sd(tmpSDANOMstation[gots]-tmpNSDANOMstation[gots])
          fulldiffsdsdanom<-c(fulldiffsdsdanom,tmpdiffsdsdanom[nN])

          # Calculate lag 1 cor with candidate CLIM ANOM
	  tmpdiffacclimanom[nN]<-cor(tmpCLIMANOMstation[gots[1:length(gots)-1]]-tmpNCLIMANOMstation[gots[1:length(gots)-1]],
	                           tmpCLIMANOMstation[gots[2:length(gots)]]-tmpNCLIMANOMstation[gots[2:length(gots)]])
          fulldiffacclimanom<-c(fulldiffacclimanom,tmpdiffacclimanom[nN])
          # Calculate lag 1 cor with candidate SD ANOM
          tmpdiffacsdanom[nN]<-cor(tmpSDANOMstation[gots[1:length(gots)-1]]-tmpNSDANOMstation[gots[1:length(gots)-1]],
	                         tmpSDANOMstation[gots[2:length(gots)]]-tmpNSDANOMstation[gots[2:length(gots)]])
          fulldiffacsdanom<-c(fulldiffacsdanom,tmpdiffacsdanom[nN])
        } # else do not add because Corr = 1.0
      }
    }
    
    # WRITE TO FILE: OLD diffs SDs
    write.table(paste(c(ds_info$statid[stattoo],
                format(round(tmpdiffsdclimanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilOLDDIFFSDCMstats,sep=""),
		append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    write.table(paste(c(ds_info$statid[stattoo],
                format(round(tmpdiffsdsdanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilOLDDIFFSDstats,sep=""),
		append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    # WRITE TO FILE: OLD diff ACs
    write.table(paste(c(ds_info$statid[stattoo],
                format(round(tmpdiffacclimanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilOLDDIFFACCMstats,sep=""),
		append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)

    write.table(paste(c(ds_info$statid[stattoo],
                format(round(tmpdiffacsdanom,digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
                file=paste(dirbawg,outfilOLDDIFFACstats,sep=""),
		append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
  # retain all of these for making a scatter plot afterwards

  }
  rm(fulldiffsdclimanom,fulldiffsdsdanom,fulldiffacclimanom,fulldiffacsdanom,tmpABSstation,
     tmpCLIMANOMstation,tmpSDANOMstation,tmpdiffsdclimanom,tmpdiffsdsdanom,tmpdiffacclimanom,
     tmpdiffacsdanom,tmpNABSstation,tmpNCLIMANOMstation,tmpNSDANOMstation)
  gc()
  stop()  
} 
#######################################################################################
if (parteight == 0) {
  print("Part Eight:0 - read in new masked data and recalc all stats")
  # Set up initial storage
  
  # while all real data will be stored in temporary arrays the output data must be stored for good or printed to file/read in again.
  # OLD parameters
  tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
  tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
  tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
  tmpACclimanom<-0		# temporary scalar for station autocorrelation (climate anomalies with loess) 
  tmpACsdanom<-0			# temporary scalar for station autocorrelation (standardised anomalies) 

  stationTREND<-array(0,c(nstations))	# array for station climatology, annual and then months
  stationCLIM<-array(0,c(nstations,13))	# array for station climatology, annual and then months
  stationSD<-array(0,c(nstations,13))	# array for station standard deviation, all period and then months

  # read in masked clean world
  mooCLEAN<-read.table(paste(dirdata,infilNEWMASKED,sep=""),header=FALSE)	# each line is now a character string to be split
  print("Finding Masked Clean Data")
  
  CleanWorld<-unname(as.matrix(mooCLEAN[,-1]))
  rm(mooCLEAN)
  gc()

  # get stations stats, neighbour lists, and AR resids

  for (stattoo in 1:nstations){

    if (restarter != "--------" & restarter != ds_info$statid[stattoo]) {
      next
    }
    restarter<-"--------"

    # reset all tmp arrays
    tmpABSstation<-array(mdi,nm)	# temporary array for storing candidate station time series - abs
    tmpCLIMANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - clim anoms
    tmpSDANOMstation<-array(mdi,nm)	# temporary array for storing candidate station time series - std anoms
    tmpACclimanom<-0		# temporary scalar for station autocorrelation (climate anomalies with loess) 
    tmpACsdanom<-0			# temporary scalar for station autocorrelation (standardised anomalies) 
    tmpCCclimanom<-array(0,40)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)
    tmpCCsdanom<-array(0,40)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)

    print(paste("Working on station ",stattoo,": ",ds_info$statid[stattoo],sep=""))
  
    # read in the candidate station masked clean world actual values
    tmpABSstation=CleanWorld[stattoo,]
  
    # INTERPOLATE OVER MISSING VALUES WHERE THEY ARE ISOLATED CASES OF 1-2 MAXIMUM CONSECUTIVE MISSING DATA
    tmpABSstation[which(tmpABSstation == mdi)]<-NA   # set missing values to R standard mdi
    tmpABSstation=interpolate_missing_func(tmpABSstation) 	

    # CREATE STATION STATS: clims, st devs, anoms, residuals from a smoothed filter
    # method CLS - calculate Climatology and remove, fit Loess and remove, calculated Standard devations and divide by them.
    res=get_climanoms_func(tmpABSstation,spanval)
    tmpCLIMANOMstation<-res$anomsarr
    tmpSDANOMstation<-res$stanomsarr
    stationTREND[stattoo]<-res$trendy
    stationCLIM[stattoo,]<-res$climsarr
    stationSD[stattoo,]<-res$sdsarr
      
    # calculate autocorrelation at lag 1 for climate anomalies and standardised (loess removed) anomalies
    moo<-acf(tmpCLIMANOMstation,lag.max=24,plot=FALSE,na.action=na.pass)
    tmpACclimanom<-moo$acf[2:25,1,1]
    moo<-acf(tmpSDANOMstation,lag.max=24,plot=FALSE,na.action=na.pass)
    tmpACsdanom<-moo$acf[2:25,1,1]
    # test for optimal AR order using aic - not perfect method
    # goo<-ar(tmpSTANOMstation[which(is.na(tmpSTANOMstation == "FALSE")],aic=TRUE)
    
    # WRITE TO FILE: station trend, ann and month clims, ann and month st devs, clim anom ac, stanom ac, smooth line, clim. anoms and st. anoms
    write.table(paste(c(ds_info$statid[stattoo],format(round(stationTREND[stattoo],digits=2),trim=FALSE,nsmall=2,width=6),
        format(round(stationCLIM[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6),
    	format(round(stationSD[stattoo,],digits=2),trim=FALSE,nsmall=2,width=6),
    	format(round(tmpACclimanom[1],digits=3),trim=FALSE,nsmall=3,width=6),
    	format(round(tmpACsdanom[1],digits=3),trim=FALSE,nsmall=3,width=6)),collapse=" "),
    	file=paste(dirbawg,outfilNEWMASKstats,sep=""),append=TRUE,sep=" ",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
    #stop() # only using this with restarter
  
  } # end of 'for'
    
  # clean up
  rm(tmpABSstation,tmpCLIMANOMstation,tmpSDANOMstation,tmpACclimanom,tmpACsdanom,stationCLIM,stationSD)
  gc()
}  
#######################################################################################
#######################################################################################
# END
##########################################################################################
