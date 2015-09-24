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
# This code contains functions to work with the get_autocovs_func.R and estimate cross-correlations based on distance and elevation
# 	expdh.corr - uses vertical and horizontal distance to estimate cross-correlation, written by Kate Willett and Richard Chandler
#   		Function to define an exponential covariance structure with horizontal and vertical distance. The covariance 
#		between sites separated by a horizontal distance of dh km and a vertical distance of dv km is:
#   			sigma^2 * exp((-phi_h*dh)+(-phi_v*dv)). 
#   		theta contains three elements: (sigma^2,phi_h,phi_v). 
#		coords is a two-column matrix of latitudes and longitudes
#		elevs is a one column matrix of station elevations in km
#		diagons is the value that should replace all diagonals to ensure positive definiteness 
#			(1.0 for lag0 cross-correlations, 0.23 for lag 1 cross-correlations).
#   		This substitutes a real covariance structure by assuming that covariances decay exponentially with distance 
#  	howfar - given latitude and longitudes for the network a matrix of distances is returned in km, written by Richard Chandler
# 
# 	howhigh - given elevations (km) for a network of stations, a vertical distance matrix is returned, written by Kate Willett
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
# expdh.corr
#   	theta: three elements: (sigma^2,phi_h,phi_v). 
#	coords: a two-column matrix of latitudes and longitudes
#	elevs: a one column matrix of station elevations in km
#	diagons: the value that should replace all diagonals to ensure positive definiteness 
# 
# howfar
# 	sites1: matrix with 2 columns, indicating the (lat, long) coordinates of each site in the first set (in degrees!)
# 	sites2:	same for the second set of sites
# 	units:  Either "km" (the default) or "nm" for distances in 
#		nautical miles (*not* nanometers!)
# 
# howhigh
#	sites1: vector with elevations in km for all stations
# 	sites2: vector with elevations in km for all stations
#
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# theta<-c(0.97,0.0006,0.07)
# coords[,1]<-c(Latitude_list[candidate_station],Latitude_list[network_station[candidate_stations]])
# coords[,2]<-c(Longitude_list[candidate_station],Longitude_list[network_station[candidate_stations]])
# elevs<-c(Elevation_list[candidate_station],Elevation_list[network_station[candidate_stations]])
# diagons<-1.0
#
# gamma_lag0<-expdh.corr(theta,coords,elevs,diagons)
#
# -----------------------
# OUTPUT
# -----------------------
# expdh.corr:
#	z: matrix of estimated correlations (using given parameters for function) of each site with every other site
#	   diagonals are forced to be the given value for the diagonals
#
# howfar:
# 	z: matrix of distances in km, with rows corresponding to the sites in sites1 and columns corresponding to 
#	   those in sites2. If sites1 and sites2 have row names, these are copied over as appropriate
#
# howhigh:
#	z: matrix of vertical distances in km, with rows corresponding to the sites in sites1 and columns corresponding to 
#	   those in sites2. If sites1 and sites2 have row names, these are copied over as appropriate
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 (14 September 2015)
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
#******************************************************
# expdh.corr
expdh.corr <- function(theta,coords,elevs,diagons) {
  # call function howfar to get matrix of horizontal distances in km
  dh <- howfar(coords,coords)
  # call function howhigh to get matrix of vertical distances in km
  dv <- howhigh(elevs,elevs)
  # work out the cross-correlations
  z <- theta[1] * exp((-theta[2]*dh) + (-theta[3]*dv))
  # force diagonals to be desired quantity
  diag(z)<-diagons
  # return the cross-correlation matrix
  return(z)
}
#-----------------------------------------------------------
# howfar
howfar <- function(sites1,sites2,units="km") {
  # work out number of rows (n1) and number of columns (n2) - should be
  # identical
  n1 <- dim(sites1)[1]; n2 <- dim(sites2)[1]
  # create matrix to hold horizontal distance matrix
  z <- matrix(nrow=n1,ncol=n2)
  # set up names of the lists
  rownames(z) <- rownames(sites1) 
  colnames(z) <- rownames(sites2) 
  # set up constant for calculating distance
  pi180 <- pi/180
  # Loop through each location and get distance
  for (i in 1:n1) {
    lat1 <- sites1[i,1]
    long1 <- sites1[i,2]
    lat2 <- sites2[,1]
    long2 <- sites2[,2]
    z[i,] <- sin(sites1[i,1]*pi180)*sin(sites2[,1]*pi180) + 
             cos(sites1[i,1]*pi180)*cos(sites2[,1]*pi180) * 
             cos((sites2[,2]-sites1[i,2])*pi180)
#
#   Next line deals with rounding errors that take the previous
#   line slightly outside the range (-1,1)
#
    z[i,] <- (60*180/pi)*acos(pmax(pmin(z[i,],1),-1))
  }
  if (units == "km") {
    z <- 1.852*z
  } else if (units != "nm") {
    stop("units must be either 'km' or 'nm'")
  }
  # return z
  z
}
#-----------------------------------------------------------
# howhigh
howhigh <- function(sites1,sites2) {
  n1 <- length(sites1); n2 <- length(sites2)
  z <- matrix(nrow=n1,ncol=n2)
  rownames(z) <- rownames(sites1) 
  colnames(z) <- rownames(sites2) 
  for (i in 1:n1) {
    z[i,] <- abs(sites2-sites1[i])	#/1000. NO lONGER NEED THIS AS ALL ELEVS ARE PASSED IN KM!!!
    
    # should be positive if sites2 station is higher and negative if sites2 station is lower
    # BUT THIS IS NOT IMPORTANT - ONLY ABSOLUTE DISTANCE
  }
  z
}

##############################################################################################
