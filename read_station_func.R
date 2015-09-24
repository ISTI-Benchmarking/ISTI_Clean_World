# R
#
# Author: Kate Willett
# Created: 1 September 2015
# Last update: 1 September 2015
# Location: /data/local/hadkw/ISTI/PROGS/
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This code reads in an ISTI station in ISTI format
# It currently just pulls out only Tmean and divides by 100. to get in degrees Celsius
# It filters the data into a full time array of given dimensions (start year and provided array)
# It returns the uncompressed station data
# 
# -----------------------
# LIST OF MODULES
# -----------------------
# R packages:
#
# Kate's R modules:
#
# -----------------------
# DATA
# -----------------------
# TheStation - a vector or 1byN array, set to mdi (-99.99), to be filled with real data
# TheStYr - the first year of TheStation, which could be before or after the earliest year present
# TheFilee - a string containing the directory and filename
# 	Missing data are -99.99
#	Assumes TheStation array starts in January and ends in December
#	Does not have an end year to read to so reads in the most recent available data
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# outputarray=read_station_func(styr,paste(directoryname,filename,sep=""),outputarray)
# 
# -----------------------
# OUTPUT
# -----------------------
# TheStation - a vector or 1byN array that has been filled with Real data
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
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
#
#
########################################################################################
read_station_func<-function(TheStYr,TheFilee,TheStation) {
  
  # Set up variables
  tmpyr<-0
  tmpmn<-0
  tmpvals<-0
  
  # read in the station file data
  mush<-read.fwf(TheFilee,widths=c(61,5,2,2,6,6,6,71),comment.char="%") # each column has one of the variables and R interprets the type (usually) correctly
  
  # The Year
  tmpyr<-mush[[2]][]
  
  # The Month
  tmpmn<-mush[[3]][]
  
  # Tmean, needs dividing by 100.
  tmpvals<-mush[[7]][]/100.
  
  # Make a map of times in the compressed data array (mypoints) to times in the uncompressed data array (pointers)
  pointers<-((tmpyr[]-TheStYr)*12)+tmpmn[]	# should be an array of length of compressed data array with numbers from -N to N
  mypoints<-array(seq(length(tmpyr)))		# should be an array of length of compressed data array from 1 to N
  mypoints<-mypoints[which(pointers > 0)]	
  pointers<-pointers[which(pointers > 0)]
  
  # Fill the station data with data, in some cases this will be missing
  TheStation[pointers]<-tmpvals[mypoints]
  
  # Clean Up
  rm(mush,tmpyr,tmpmn,tmpvals,pointers,mypoints)
  gc()
  
  # Return the station data array
  return(TheStation)
}
########################################################################################
