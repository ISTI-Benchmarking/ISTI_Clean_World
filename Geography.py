#!/usr/local/sci/bin/python

#***************************************
# Code to work out useful things from spatial data:
#
#1. TwoPointDistanceKm
#   Distance between two points in space in km
#   Needs latitudes (-90 to 90) and longitudes (-180 to 180)
#   www.johndcook.com/python_longitude_latotude.html
#   6 June 2014 (D-Day!) KMW - v1
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
#************************************************************************
# Set up python imports
import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scipy.optimize import curve_fit,fsolve,leastsq
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.stats
from math import sqrt,pi
import struct
import math
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

#************************************************************************
# Subroutines
#************************************************************************
# TwoPointDistanceKm
def TwoPointDistanceKm(latC, longC, latN, longN): 
    ''' Calculates distance point to point  on unit sphere in km'''
    ''' Can work on single elements '''
    ''' For 2+ elements, a 1 by npoints array is returned '''
    ''' LatC/LongC is the candidate location '''
    ''' LatN/LongN are the other neighbours that the candidate should be compared to '''

# Ensure that inputs are numpy arrays
    latC=np.array(latC,ndmin=1)
    latN=np.array(latN,ndmin=1)
    longC=np.array(longC,ndmin=1)
    longN=np.array(longN,ndmin=1)
    
# Create output Array
    npoints=len(latN)
    TheKms=np.empty([npoints],dtype=object)

# Convert latitude and longitude to sperical coordinates in radians
    degs_to_rads=math.pi/180.
    
# phi = 90 - latitude
    phiC=(90.0-latC)*degs_to_rads
    phiN=(90.0-latN)*degs_to_rads
    
# theta=longitude
    thetaC=longC*degs_to_rads
    thetaN=longN*degs_to_rads
    
# Compute spherical distance from spherical coordinates
# For two locations in sperical coordinates
# (1,theta, phi) and (1,theta, phi)
# cosine(arc length)=sin phi sin phi' cos(theta-theta') + cose phi cos phi'
# distance = rho * arc length

    if npoints > 1:
        for loopoo in range(npoints):
	    if phiC == phiN[loopoo] and thetaC == thetaN[loopoo]: 
	        TheKms[loopoo]=0.0
	    else:
	        cos=(math.sin(phiC)*math.sin(phiN[loopoo])*math.cos(thetaC-thetaN[loopoo])+math.cos(phiC)*math.cos(phiN[loopoo]))
            #print(loopoo,cos)
	        arc=math.acos(cos)
                TheKms[loopoo]=arc*6373.
	    #print(arc*6373.)
    else:
        cos=(math.sin(phiC)*math.sin(phiN)*math.cos(thetaC-thetaN)+math.cos(phiC)*math.cos(phiN))
        arc=math.acos(cos)    
        TheKms=arc*6373.
    
    return  TheKms # TwoPointDistanceKms

#************************************************************************
