# R
#
# Author: Kate Willett
# Created: 2 September 2015
# Last update: 23 September 2015
# Location: /data/local/hadkw/ISTI/PROGS/
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This function calculates annual and monthly climatology over entire period of record (climsarr)
#	(minimum 3 years - sorted out already see READMEBAWG_SEP2015)
# Removes climatological means to create climate anomalies (anomsarr)
# Fits a linear trend model to the climate anomalies and stores as a decadal trend (trendy)
# Fits a lowess curve (smoothie) representing a smooth trend through the data using spanval (proportion of data to smooth)
# Removes the smoothed fit from the climate anomalies (detrends and removes inhomogeneous features hopefully)
# Calculates whole period and monthly st devs (sdsarr)
# Creates standardised anomalies (stanomsarr) by dividing by climatological month st devs
#
# If the station is too short then it returns bunk=9
# NOTE: data are first cleaned to remove outliers: USS0004B02S is a particular issue with 0 values and -12 to 35 December range 
# outlier if < quantile(data,0.25)-2.5*IQR(data)
# outlier if > quantile(data,0.75)+2.5*IQR(data)
# 
# Output is a list of all of these things that must be disected
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# R packages:
#
# Kate's R modules:
# interpolate_missing_func.R - interpolates over single missing values, written by Kate Willett
# -----------------------
# DATA
# -----------------------
# absarr - a vector or 1byN array of monthly means
# 	Missing data are NA
#	Assumes TheStation array starts in January and ends in December
# spanval - a scalar float stating the proportion of data to smooth between 0 and 1
#	This has already been set as monthspan/totalmonths (~0.3 or 600 months)
#	If spanval=9 then this switches off the loess smooth and creation of standardised anomalies
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# res=get_climanoms_func(monthlymeans,spanval)
#
# -----------------------
# OUTPUT
# -----------------------
# A list of things:
#
# climateanomalies<-res$anomsarr
# standardisedanomalies<-res$stanomsarr
# smoothlowess<-res$smoothie
# decadaltrend<-res$trendy
# climatology_means<-res$climsarr
# climatology_stdevs<-res$sdsarr
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# Version 3 (23 September 2015)
# ---------
#  
# Enhancements
# Only runs standardised anomalies if needed
# If you just want climate anomalies then set spanval to 9
# No loess smooth or standardised anomalies will be created and blank arrays will be returned.
#  
# Changes
#  
# Bug fixes
#
# Version 2 (9 September 2015)
# ---------
#  
# Enhancements
# Some stations are really messy. e.g. USS0004B02S has 0 values and a December range from -11.65 to +35
# This results in terrible 'clean' output with crazy variability and large inhomogeneities
# We can clean up the data to some degree before calculating the mean,sd and resids
# We can also make the minimum data amount requirement higher for resids calculation - see get_autocovs_func.R
# I have added an outlier test which swaps outliers for missing data:
#	outliers = data < quantile(data,0.25)-2.5*IQR(data)
#	outliers = data > quantile(data,0.75)+2.5*IQR(data)
# Hopefully more robust than mean, st dev.
#  
# Changes
#  
# Bug fixes
#
# 
# Version 1 (1 September 2015)
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
# Missing data must have been set to NA
#
########################################################################################
get_climanoms_func<-function(absarr,spanval) {    
  
  # set up functions to call
  source("interpolate_missing_func.R")
  
  # Set up variables
  
  longoo<-length(absarr)	# length of station time series
  anomsarr<-array(NA,longoo)	# empty array for climate anomalies
  stanomsarr<-array(NA,longoo)	# empty array for standardised anomalies
  smoothie<-array(NA,longoo)	# empty array for smoothed trend
  climsarr<-array(NA,13)	# empty array for annual and monthly climatological means
  sdsarr<-array(NA,13)		# empty array for full period and monthly climatological standard deviations
  trendy<-0.
  
  # reshaped the actual monthly means array to be 12 rows and Nyears columns
  valrebin<-array(absarr,dim=c(12,longoo/12))	# reform to 12 rows by nyears
  
  # remove outliers (tested and ok'd)
  valrebin<-t(apply(valrebin,1,function(x) get_outliers_func(x)))
  
  # sums all of the years present for each month (missing data must have been set to NA!)
  sizeclims<-apply(valrebin,c(1),function(x) sum(table(x)))	# should sum number of 'TRUE' values only

  # IF minimum month count is greater than 3 (must have 3 years of complete data) then we're good to go
  # IF not then we'll have to ditch the station  
  if (min(sizeclims) >= 3) { # I THINK THIS SHOULD BE 3!!!

    # efficient creation of climate anomalies (scale=FALSE means do not divide by st. dev)
    anomsarr<-array(scaltro<-t(apply(valrebin,1,scale,scale=FALSE)),dim=(longoo))  
    
    # efficient creation of climatological means
    climsarr[2:13]<-apply(valrebin,1,mean,na.rm=TRUE)
    
    # Mean of means will do for annual climatological mean
    climsarr[1]<-mean(climsarr[2:13],na.rm=TRUE)
    
    # basic linear trend fit multiplied by 120 months to be decadal
    trendy<-120.*coef(lm(anomsarr~seq(longoo)))[2]	# this is the decadal rate of change linear trend, should be omitting NAs - not perfect but sufficient
    
    if (spanval != 9) {	# then carry on and build standardised anomalies
      # a smooth lowess curve fit to the data given spanval
      smoothie<-predict(loess(anomsarr~seq(longoo),span=spanval),array(seq(longoo),dim=(longoo)))
    
      # subtract smooth trend from data
      stanomsarr[]<-anomsarr-smoothie

      # efficient calculation of climatological standard deviations for each month
      sdsarr[2:13]<-apply(array(stanomsarr,dim=c(12,longoo/12)),1,sd,na.rm=TRUE)
      # standard deviation of entire record
      sdsarr[1]<-sd(stanomsarr,na.rm=TRUE)
    
      # divide by standard deviation to create standardised anomalies
      stanomsarr<-array(t(apply(array(stanomsarr,dim=c(12,longoo/12)),1,scale)),dim=(longoo))
    
      # we would like as complete data as possible and with a bit of QC for making sure real residuals aren't too screwy
      # first remove any ridiculously large anomalies
      huges<-which(abs(stanomsarr) > 4.5)
      if (length(huges) > 0) {
        stanomsarr[huges]<-NA
        # now rerun interpolate_missing_func to cover up the isolated missing outliers ***
        stanomsarr=interpolate_missing_func(stanomsarr) 	
      }
    }
  }
  
  # return a list of all of the elements
  return(list(climsarr=climsarr,sdsarr=sdsarr,anomsarr=anomsarr,stanomsarr=stanomsarr,smoothie=smoothie,trendy=trendy))
}
###################################################################################################
get_outliers_func<-function(datavec) {
  # uses quantiles and IQR to identify outliers
  # works on vectors
  
  Low_Cutoff<-quantile(datavec,0.25,na.rm=TRUE)-(2.5*IQR(datavec,na.rm=TRUE))
  High_Cutoff<-quantile(datavec,0.75,na.rm=TRUE)+(2.5*IQR(datavec,na.rm=TRUE))
  
  #print(c(Low_Cutoff,High_Cutoff))
  
  datavec[datavec < Low_Cutoff]=NA
  datavec[datavec > High_Cutoff]=NA
  
  return(datavec)
}
###################################################################################################
