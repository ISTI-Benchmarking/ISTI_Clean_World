#!/usr/local/sci/bin/python

#***************************************
# Code to fit linear trends to time series data:
#
#1. Median of Pairwise slopes estimation of a linear trend
#   5th and 95th percentile slopes given as 5-95 percent confidence intervals
#   Copes with missing data
#   After Sen, 1961
#   Based on IDL code - Mark McCarthy
#
#
#
# 11 October 2012 KMW - v1
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

#************************************************************************
# Subroutines
#************************************************************************
# MedianPairwise
def MedianPairwise(TheData,TheMDI,TheSlope): # ,Lowee=Lowee,Highee=Highee):
    ''' Calculates slope from every point to every other point '''
    ''' Outputs the median of those slopes at a rate of unit per time step '''
    ''' Optionally outputs the 5th and 95th percentiles of those slopes as uncertainty ranges '''
    ''' If Lowee and/or Highee are set they will come out as changed values '''
    TheSlope=[0.,0.,0.]		# median, 5th adn 95th percentile rate of change per time step
    PairwiseSlopes=[]
    for i,ii in enumerate(TheData):
#        print i,ii
        if ii > TheMDI:
	    for r in range(len(TheData)-1,i,-1):
	        if TheData[r] > TheMDI:
		    PairwiseSlopes=np.append(PairwiseSlopes,((TheData[r]-ii)/(r-i)))
#		    print i,r,ii,TheData[r],(TheData[r]-ii)/(r-i)
    
#    print "%5.2f"*len(PairwiseSlopes) % tuple(PairwiseSlopes)
    if len(PairwiseSlopes) > 10:	# check there are sufficient data
	TheSlope[0]=np.median(PairwiseSlopes)
        nData=len(list(TheData[np.where(TheData > TheMDI)]))
        DegofFree=(nData * (nData-1))/2.
	weight=np.sqrt(nData * (nData-1) * ((2*nData)+5) / 18.)
#	print "No. PW, No. Data Present, Deg Free ",len(PairwiseSlopes),nData, DegofFree
#	print "WEIGHT ", weight
        rankL=int(((DegofFree-1.96*weight)/2.)-1)
        rankU=int(((DegofFree+1.96*weight)/2.)+1)
#	print "RANKS ",rankL,rankU
	PairwiseSlopes.sort()
	TheSlope[1]=PairwiseSlopes[rankL]
	TheSlope[2]=PairwiseSlopes[rankU]
         	
#        print(TheSlope)
#        num_bins = 50
#        # the histogram of the data
#        n, bins, patches = plt.hist(PairwiseSlopes, num_bins, normed=1, facecolor='green', alpha=0.5)
#        plt.show()
    
    ''' At present I think the 5th adn 95th are too sensitive to variability '''
    ''' This is a quicklook Kate version of the code and needs to be statistically verified '''
    return  TheSlope # ReadData

#************************************************************************
