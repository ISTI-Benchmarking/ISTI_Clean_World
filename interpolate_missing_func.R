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
# interpolate_missing_func
# This function finds missing data over single missing points and interpolates values to create longer non-missing periods.
# For any points that are not the first or last in the series, a value is given for the missing data based on interpolation between the time point before and the time point after the missing data.
#
# interpolate_missing_func_2
# This function finds missing data over two consecutive missing points and interpolates values to create longer non-missing periods.
# For any points that are not the first or last in the series, a value is given for the missing data based on interpolation between the time point before and the time point after the missing data.
# Two values are given by interpolating from the time point before and after by:
# 	dividing the difference by 3
#	adding 1/3rd to missing time point 1 
#	adding 2/3rds to missing time point 2
#
# interpolate_missing_func_3
# This function finds missing data over three consecutive missing points and interpolates values to create longer non-missing periods.
# For any points that are not the first or last in the series, a value is given for the missing data based on interpolation between the time point before and the time point after the missing data.
# Three values are given by interpolating from the time point before and after by:
# 	dividing the difference by 4
#	adding 1/4th to missing time point 1 
#	adding 2/4ths to missing time point 2
#	adding 3/4ths to missing time point 3

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
# TheData - a vector or 1byN array to be interpolated where missing
# 	Missing data are -99.99
#	Assumes TheStation array starts in January and ends in December
#	Does not have an end year to read to so reads in the most recent available data
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# outputarray=interpolate_missing_func(inputarray) 	
# 
# -----------------------
# OUTPUT
# -----------------------
# TheData - a vector or 1byN array that has been interpolated
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
# the default method (used here) is strictly monotonic. Pearhaps other options could improve the adjustment. 
# How about using harmonics for the interpolation? 
# INPUT:
# datos: a data vector; missing must be "NA"
# method: parameter methods from the spline function; set to the default monotonic interpolation. 
# OUTPUT:
# one: the interpolated vector. 
#
########################################################################################
interpolate_missing_func<-function(TheData,method='fmm'){
  # Fills in all single missing points
  
  # Make a copy that we change
  NewData<-TheData 
  
  # Length of the data array minus two points
  ne<-length(TheData)-2
  
  # One is missing
  for(i in 1:ne){
    # the potential missing point
    j<-i+1
    # the point after the potential missing point
    k<-i+2 	
    # here i is two observarions appart from k; both must exist and we interpolate  j from them
    if (!is.na(TheData[i]) & is.na(TheData[j]) & !is.na(TheData[k])){
      # calculate interpolated value
      NewData[j]<-spline(TheData[c(i,k)],xout=1.5,method=method)$y
    }
  }
  
  # Return the more filled in array (any points at ends or greater than one will not be filled)
  return(NewData)

}
########################################################################################
interpolate_missing_func_2<-function(TheData,method='fmm'){
  #Fills in all two-consecutive missing points     

  # Make a copy that we change
  NewData<-TheData 
  
  # Length of the data array minus three points
  ne<-length(TheData-3)
  
  # Two are missing
  for(i in 1:ne){
    # the potential missing points
    na1<-i+1
    na2<-i+2
    # the point after the potential missing points
    k=i+3
    # here i and k are 3 obs appart; both exist and we interpolate na1 & na2 from them
    if (!is.na(TheData[i]) & is.na(TheData[na1]) & is.na(TheData[na2]) & !is.na(TheData[k])){
      # calculate interpolated values
      sp<-spline(TheData[c(i,k)],xout=c(1+1/3,1+2/3),method=method)$y
      NewData[na1]<-sp[1]
      NewData[na2]<-sp[2]
    }
  }
  # Return the more filled in array (any points at ends or greater than one will not be filled)
  return(NewData)
}
#########################################################################################
interpolate_missing_func_3<-function(TheData,method='fmm'){
  #Fills in all three-consecutive missing points     

  # Make a copy that we change
  NewData<-TheData 
  
  # Length of the data array minus four points
  ne<-length(TheData-4)
  
  #Three are missing
  for(i in 1:ne){
    # the potential missing points
    na1<-i+1
    na2<-i+2
    na3<-i+3
    # the point after the potential missing points
    k=i+4
    # here i and k are 3 obs appart; both exist and we interpolate na1 & na2 from them
    if (!is.na(TheData[i]) & is.na(TheData[na1]) & is.na(TheData[na2]) & is.na(TheData[na3])& !is.na(TheData[k])){
      # calculate interpolated values
      sp<-spline(TheData[c(i,k)],xout=c(1+1/4,1+2/4,1+3/4),method=method)$y
      NewData[na1]<-sp[1]
      NewData[na2]<-sp[2]
      NewData[na3]<-sp[3]
    }
  }
  # Return the more filled in array (any points at ends or greater than one will not be filled)
  return(NewData)
}
#########################################################################################
