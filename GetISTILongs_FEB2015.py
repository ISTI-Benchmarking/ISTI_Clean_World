#!/usr/local/sci/bin/python

#***************************************
# 06 June 2014 KMW - v1
# 
# list ISTI stations with greater than or equal to 3 years of TMean data (3 of each month)
# that all occur within a decade of data
# 
# 29 July 2015
# changed line in GetStationLength:    
#	fullyears=(2015+1)-1850 (was 2014 instead of 2015)
# changed line in GetStationLength:
#	if (TMeans[timee] > TheMDI) & (timpointer >= 0): (added the second parentheses in case ISTI ever predates 1850)
# changed line in GetStationLength for 'years' instead of 'tims' NO DIFFERENCE FOUND
#       for window in range((fullyears+1)-PeriodMax): (was looping through times in an array that was only years long)
# added a loop to look at whether a station can have a 30 year climatology calculated (must have at least 15 complete years within a 20 year window)
#       in GetStationLength
#       outputs total number of stations for which we can calculate a climatology	
#***************************************
#************************************************************************
# USE python2.7
# python2.7 GetISTILongs_FEB2015.py
#
# REQUIRES
#************************************************************************
# Set up python imports
import datetime as dt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import date2num,num2date
import sys, os
from scipy.optimize import curve_fit,fsolve,leastsq
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.stats
from math import sqrt,pi
import struct

from Geography import TwoPointDistanceKm

# RESTART VALUE
Restarter='------'	#'------'		#' RSM00022140' station ID

# Set up file locations
STATLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/INVENTORY_monthly_merged_stage3_proxyelevs'
#STATLIST='/data/local/hadkw/ISTI/LISTS/ISTIINVENTORY_corrs_stage3_JUN2014.dat'
INFILEE='/data/local/hadkw/ISTI/DATA/ISTIv101_JUL2015/results_merged/merge_'	#_stage3
OUTLONGS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGS_stage3proxyelevs_JUL2015.dat'
OUTSHORTS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTISHORTS_stage3proxyelevs_JUL2015.dat'

# Set up variables and holders
MDI=-9999
nstations=0	# defined after reading in station list
StationIDs=[]	# nstations list filled after reading in station list

MinMonthCount=3 # Minimum count for months to be present e.g., 3 Januarys
MaxDuration=10	 # Maximum duration of years over which MinMonthCount can be spread
#************************************************************************
# Subroutines
#************************************************************************
# READDATA
def ReadData(FileName,typee,delimee,ASTruth,ColumnChoice):
    ''' Use numpy genfromtxt reading to read in all rows from a complex array '''
    ''' Need to specify format as it is complex '''
    ''' outputs an array of tuples that in turn need to be subscripted by their names defaults f0...f8 '''
    ''' ISTI INVENTORY USES # WHICH GENFROMTXT ASSUMES ARE COMMENTS - IGNORES ALL TEXT AFTERWARDS '''
    ''' HENCE comments=';' - HOPEFULLY NO ; IN ISTI. '''
    ''' NOTE ASTruth needs to be either True or False '''
    ''' NOTE Column choice either needs to be a list of numbers or 'XXX' '''
    if ColumnChoice == 'XXX':
        return np.genfromtxt(FileName, dtype=typee,comments="%",delimiter=delimee) # ReadData    
    else:
        return np.genfromtxt(FileName, dtype=typee,comments="%",delimiter=delimee,autostrip=ASTruth,usecols=ColumnChoice) # ReadData

#************************************************************************
# GETSTATIONARRAY
def GetStationLength(TheFilee,TheMDI,NMonthsMin,PeriodMax,TestClim):
    ''' Read in the data from the ISTI file format '''
    ''' Set up full array to include missing months of data '''
    ''' Check how many Jans, Febs, etc '''
    ''' There must be at least NmonthsMin of each month '''
    ''' These must be within a PeriodMin year period '''
    ''' Returns a flag with 'stubby' or 'lanky' '''
    
# Read in the raw data    
    MyTypes=("|S18","|S14","float","float","float","|S1","int","int","|S2","int","int","int","|S71")
    MyDelimiters=[18,14,9,11,9,1,4,2,2,6,6,6,71]
    MyColumns=(6,7,11)
    RawData=ReadData(TheFilee,MyTypes,MyDelimiters,True,MyColumns)
    Years=np.array(RawData['f6'],ndmin=1)
    #print(Years)
    Months=np.array(RawData['f7'],ndmin=1)
    TMeans=np.array(RawData['f11'],ndmin=1)
    ntims=len(Years)
    #print(Years[0],Years[len(Years)-1],TMeans[0:10]/100.)
   
# Fill in the missing data from 1850 onwards (end of 2014 at present) (ignore time before that)
    fullyears=(2015+1)-1850
    fulltims=fullyears*12
    #print(fulltims)
    TheStationData=np.empty((fulltims),dtype=float)
    TheStationData[:]=TheMDI
    
    for timee in range(ntims):
        timpointer=((Years[timee]-1850)*12)+Months[timee]-1
	#print(ntims,timee,timpointer,Years[timee],Months[timee])
	if (TMeans[timee] > TheMDI) & (timpointer >= 0):
	    TheStationData[timpointer]=TMeans[timee]/100.
	    

# Look at month presence   
    IsItLanky='SHORT'
    TheStationData=TheStationData.reshape(fullyears,12)
    # move through the data 10 years at a time - if you find 3 years set lankyflags to 1 and continue
    for window in range((fullyears+1)-PeriodMax):
        lankyflags=np.zeros(12)
        for mm in range(12):
    	    gots=np.where(TheStationData[window:(window+PeriodMax),mm] > TheMDI)[0]
	    if len(gots) >= NMonthsMin:
	        #print(window,gots)
	        lankyflags[mm]=1
        #print(window,lankyflags)
	if sum(lankyflags) == 12:
            IsItLanky='LANKY'
	    break

# Test to see if this station would pass the 30 year climatology test
# Minimum of 15 years complete years of data within a 30 year period  
    # move through the data 30 years at a time - if you find 15 complete years set TestClim to 1 and continue
    #stop
    for window in range((fullyears+1)-30):
        climflags=0
        for yy in range(30):
    	    gots=np.where(TheStationData[window+yy,:] > TheMDI)[0]
	    if len(gots) == 12:
	        #print(window,gots)
	        climflags=climflags+1
        if climflags >= 15:
            TestClim=1
	    #stop
	    break
    
    print(lankyflags,TestClim)
#    stop()
		    
    return IsItLanky,TestClim # GetStationLength
#***********************************************************************
# MAIN PROGRAM
#***********************************************************************
MyTypes=("|S12","|S31","|S24","float","float","float","int","int","int","int","int","int","|S22","|S11")
MyDelimiters=[12,31,24,8,11,9,5,5,5,5,5,5,22,11]
MyColumns='XXX'
RawData=ReadData(STATLIST,MyTypes,MyDelimiters,False,MyColumns)
StationIDs=np.array(RawData['f0'],ndmin=1)
nstations=len(StationIDs)
print('Read in Giant file...')

My_FLhandle=file(OUTLONGS,'a')
My_FShandle=file(OUTSHORTS,'a')
CountClims=0	# Total number of stations passing climatology test

for lop in range(nstations):
    #print(Restarter,StationIDs[lop])
    if Restarter != '------' and Restarter != StationIDs[lop]:
        continue
    else: 
        Restarter='------'
    
    print(StationIDs[lop])
    # sort out RawData format
    outraw=list(RawData[lop])
    outraw[3]="{:8.4f}".format(outraw[3])
    outraw[4]="{:11.4f}".format(outraw[4])
    outraw[5]="{:9.2f}".format(outraw[5])
    outraw[6:12]=["{:5d}".format(dd) for dd in outraw[6:12]]
    outraw=np.reshape(outraw,(1,len(outraw)))     

# read in the station data for the candidate file
    MyFile=INFILEE+StationIDs[lop].strip()+'_stage3'
    TestClim=0	# turn it to 1 if it passes
    IsItLanky,TestClim=GetStationLength(MyFile,MDI,MinMonthCount,MaxDuration,TestClim)
    print(IsItLanky,TestClim)
    CountClims=CountClims+TestClim

    if IsItLanky == 'LANKY':
        np.savetxt(My_FLhandle,outraw,fmt='%s',delimiter='')
    else:
        np.savetxt(My_FShandle,outraw,fmt='%s',delimiter='')
    

My_FLhandle.close()
My_FShandle.close()

print("Total No. Stations passing climatology test: ",CountClims)

print("And, we are done!")
