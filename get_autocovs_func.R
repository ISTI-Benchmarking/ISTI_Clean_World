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
# This function creates matrices of Vector Autoregression (VAR) parameters for networks of stations
# It uses distance and elevation difference between stations as a proxy to estimate correlations
# The distance elevation function must be provided and has already been defined - see READMEBAWG_SEP2015
# This code only deals with VAR(1) At = PHI*At-1 + Zt. 
# Given a network of stations: latitude, longitude. elevation (m)
# Given a distance elevation function
# The code calls expdh.corr to get the cross correlation matrices
# The code in turn calls howfar and howhigh to get the distance matrices
# This code returns:
#	N by N matrix of cross-correlations at lag 0 where diagonals are forced to be 1.0 (gamma_lag0)
#	N by N matrix of cross-correlations at lag 1 where diagonals are forced to be a given value (gamma_lag1)
#	N by N matrix of VAR(1) parameters (phi_val): phi_val<-t(gamma_lag1)%*%solve(gamma_lag0)
#		For stability of matrix inversion the Kronecker product method is used instead of solve()
#
# -----------------------
# LIST OF MODULES
# -----------------------
# R packages:
#
# Kate's R modules:
# distance_elevation_func.R - written by Kate Willett
# 	expdh.corr - uses vertical and horizontal distance to estimate cross-correlation, written by Kate Willett and Richard Chandler
#  	howfar - given latitude and longitudes for the network a matrix of distances is returned, written by Richard Chandler
# 	howhigh - given elevations (km), a vertical distance matrix is returned, written by Kate Willett
# -----------------------
# DATA
# -----------------------
# thetaparms0 - parameters for estimating cross correlations at lag0 - this a four element array:
#     1: sigma, (e.g, 0.97)
#     2: dist_horizontal(km) (e.g., 0.0006) 
#     3: dist_vertical(km),  (e.g., 0.07)
#     4: diagonal value (e.g. 1)
# thetaparms1 - parameters for cross correlations at lag1 - this a four element array:
#     1: sigma, (e.g., 0.23)
#     2: dist_horizontal(km) (e.g., 0.0006) 
#     3: dist_vertical(km),  (e.g., 0.07)
#     4: diagonal value (e.g. 0.23)
# TheLatLongCoords - N rows by 2 column array
#	Column 1: Latitudes
# 	Column 2: Longitudes
# TheElevs - N element array of station elevations in metres (later converted to kilometers)
#
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# NetworkCoords[,1]<-c(StationLatitude,NeighbourLatitudes)
# NetworkCoords[,2]<-c(StationLongitude,NeighbourLongitudes)
#
# res=get_autocovs_func(c(0.97,0.0006,0.07,1.0),c(0.23,0.0006,0.07,0.23),NetworkCoords,NetworkElevs)  
#
# -----------------------
# OUTPUT
# -----------------------
# A list of things:
#
# VAR_parameters<-res$phi_val - VAR(1) parameters     
# Gamma_Lag0<-res$gamma_lag0 - cross-correlation matrix at lag 0 based on distances (vertical and horizontal)
# Gamma_Lag1<-res$gamma_lag1 - cross-correlation matrix at lag 1 based on distances (vertical and horizontal)
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 2 (14 September 2015)
# ---------
#  
# Enhancements
#  
# Changes:
# I have removed the section that calculated the real or proxy residuals to get_ARresids_func.R
#  
# Bug fixes
#  
# Version 1 (2 September 2015)
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
get_autocovs_func<-function(thetaparms0,thetaparms1,TheLatLongCoords,TheElevs) {
  # Set up functions to call
  source("distance_elevation_func.R")

  # Set up variables
  phi_val<-0
  newACs<-0
  testAC<-0
  testee<-"BAD"
  
  nsts<-length(TheElevs)	# number of stations within the network (should be 41)
  
  # Estimate cross-correlations at lag 0   
  gamma_lag0<-expdh.corr(thetaparms0[1:3],TheLatLongCoords,TheElevs/1000.,thetaparms0[4])  
  # Estimate cross-correlations at lag 1   
  gamma_lag1<-expdh.corr(thetaparms1[1:3],TheLatLongCoords,TheElevs/1000.,thetaparms1[4])  

  # Kronecker product method to work out VAR(1) parameters
  # Tested this JUL2015 and it matches Lund and Willett paper for output
  I=diag(nsts)
  AA=I%x%gamma_lag0
  c=matrix(gamma_lag1,nsts*nsts,1)
  phi_val<-t(matrix(solve(AA,c),nsts,nsts)) 

#  phi_val<-t(gamma_lag1)%*%solve(gamma_lag0)
  
#  # TEST IT OUT
#  newACs<-array(0,5)
##  for (nn in 1:5) {
##    shocks<-rnorm(nsts)
##    moo<-phi_val%*%rnorm(nsts) + (((shocks - mean(shocks))/sd(shocks))*1.) # this was 1.5 which is why its still in
##    shocks<-rnorm(nsts)
##    moo<-phi_val%*%moo + (((shocks - mean(shocks))/sd(shocks))*1.)
##    shocks<-rnorm(nsts)
##    moo<-phi_val%*%moo + (((shocks - mean(shocks))/sd(shocks))*1.)
##    shocks<-rnorm(nsts)
##    moo<-phi_val%*%moo + (((shocks - mean(shocks))/sd(shocks))*1.)
##    shocks<-rnorm(nsts)
##    moo<-phi_val%*%moo + (((shocks - mean(shocks))/sd(shocks))*1.)
##    ts<-moo[1]
##    for (tt in 1:600) {
##      shocks<-rnorm(nsts)
##      moo<-phi_val%*%moo + (((shocks - mean(shocks))/sd(shocks))*1.)
##      ts<-c(ts,moo[1])
##    }
##    plot(ts)
##    print(paste(nn,"NEW AC:",cor(ts[2:600],ts[1:599]),sep=" "))
##    newACs[nn]<-cor(ts[2:600],ts[1:599])
##  }
  
  # return the things 
  return(list(phi_val=phi_val,newACs=newACs,gamma_lag0=gamma_lag0,gamma_lag1=gamma_lag1))
}
########################################################################################
