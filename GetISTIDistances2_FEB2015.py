#!/usr/local/sci/bin/python

#***************************************
# 06 June 2014 KMW - v1
# 
# Rerun of distances to narrow to final station listing
# Loop through new ISTI station inventory
# calculate the distance between every station
# output:
#       Closest 100 stations in order: stA st1 dist1 st2 dist2 etc
# 	Closest 1000 stations in order: stA st1 dist1 st2 dist2 etc
# 	Complete distance matrix in 9+ 10000 by 10000 station files
#
# Find all locations that match (just double checkgin)
# Remove the matching station from the Distance Lists and 
# do not include later on
# Make new INVENTORY list and list bad stations
#
# Output 40 nearest neighbours for each station: FIXCORRNEIGHBOURS...
# 
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 GetISTIDistances2_FEB2015.py
#
# REQUIRES
# Geography.py
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
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

from Geography import TwoPointDistanceKm

# RESTART VALUE
Restarter='------'	#'------'		#'    01010000' station ID

# Set up file locations
STATLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat'
OUThundred='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat'
OUTthousand='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGDISTANCES_thousand_stage3proxyelevs_JUL2015.dat'
OUTGOODS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORYcheck_stage3proxyelevs_JUL2015.dat'
OUTBADS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGSLOCMATCHcheck_stage3proxyelevs_JUL2015.dat'
OUTneighbours='/data/local/hadkw/ISTI/LISTS/BAWG/FIXCORRNEIGHBOURS_ISTI_stage3proxyelevs_JUL2015.dat'

nstations=0	# defined after reading in station list
ngoods=0	# stations passing unique location criteria
nbads=0		# shortest stations with matching locs
StationIDs=[]	# nstations list filled after reading in station list
StationLats=[]	# nstations list filled after reading in station list
StationLons=[]	# nstations list filled after reading in station list

StatDistsAll=[]
StatDistsSorted=[]
StatIDsSorted=[]

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
    if ColumnChoice == 'XXX':
        return np.genfromtxt(FileName, dtype=typee,comments="%",delimiter=delimee) # ReadData    
    else:
        return np.genfromtxt(FileName, dtype=typee,comments="%",delimiter=delimee,autostrip=ASTruth,usecols=ColumnChoice) # ReadData

#************************************************************************
# GETDISTANCEARRAYS
def GetDistanceArrays(TheCID,TheDistsAll,TheSortedDists,TheSortedIDs,
                      TheCLat,TheCLon,TheLats,TheLons,TheIDs,TheSCount,AnyBads):
    ''' Call TwoPointDistancesKms to get a matrix of distances '''
    ''' Sort into closest 100 and 1000 '''
    
# Call TwoPointDistancesKms to work out matrix of distances for all stations
    #print(TheCLat,TheCLon,TheLats[4993:4996],TheLons[4993:4996],TheIDs[4993:4996])
    TheDistsAll=TwoPointDistanceKm(TheCLat,TheCLon,TheLats,TheLons)
    #print(TheDistsAll[4993:4996])
    
# For each station, sort distances and pull out 100, 1000 closest    
    SortIndex=np.argsort(TheDistsAll)
    #print(SortIndex)
    TheSortedDists=TheDistsAll[SortIndex]
    TheSortedIDs=TheIDs[SortIndex]
    #print(TheSortedDists[0:10])

# Remove the distance for Candidate station with Candidate station
    got=np.where(TheCID == TheSortedIDs)[0]
    TheSortedIDs=np.delete(TheSortedIDs,got)
    TheSortedDists=np.delete(TheSortedDists,got)

# Is there a 0.0 distance? A Location Match!
# If so - flag - this station will then be removed
# If there are multiple stations then all but the last will be removed
    AnyBads=np.where(TheSortedDists == 0.0)[0]
    	
    return TheDistsAll,TheSortedDists,TheSortedIDs,AnyBads # GetDistanceArrays
  
#************************************************************************
# WRITEOUTSORTED
def WriteOutSorted(TheStDists,TheStIDs,TheFile,TheCandidateID):
    ''' Output lines to text of StationID, list of stations and distances '''
    
# Convert all distances to set length strings
    #print(TheStDists[0:10])
    #TheStDistsStr=np.array(["{:9.3f}".format(dd) for dd in TheStDists.reshape(TheStDists.size)])
    #TheStDistsStr=TheStDistsStr.reshape(TheStDists.shape)
    TheStDistsStr=["{:12.3f}".format(dd) for dd in TheStDists]
    #print(TheStDistsStr[0:10])
     
# Make a nstations (rows) by 2 column array, reform such that it becomes r1c1,r2c1,r1c2,r2c2,r1c3,r2c3 etc
    TheData=np.reshape(zip(*np.vstack((TheStIDs,TheStDistsStr))),len(TheStIDs)*2)	# a one by nstations array 
    goo=np.reshape(np.append(TheCandidateID,TheData),(1,(len(TheStIDs)*2)+1))
    
    np.savetxt(TheFile,goo,fmt='%s',delimiter=' ')
    return #WriteOutSorted
    
#***********************************************************************
# WRITEOUT
def WriteOut(TheStIDs,TheFile,TheStationID):
    ''' Output a line for each station of the station ID '''
    ''' and its 40 nearest neighbours '''

    # Remove white space
    TheStationID=TheStationID.strip()
    TheStIDs=[dd.strip() for dd in TheStIDs]
    goo=np.reshape(np.append(TheStationID,TheStIDs),(1,len(TheStIDs)+1))
    
    np.savetxt(TheFile,goo,fmt='%s',delimiter=' ')
    return #WriteOut
    
#***********************************************************************
# MAIN PROGRAM
#***********************************************************************
# read in station list
MyTypes=("|S12","|S31","|S24","float","float","float","int","int","int","int","int","int","|S16","|S17")
MyDelimiters=[12,31,24,8,11,9,5,5,5,5,5,5,16,17]
MyColumns='XXX'
RawData=ReadData(STATLIST,MyTypes,MyDelimiters,False,MyColumns)
StationIDs=np.array(RawData['f0'],ndmin=1)
StationLats=np.array(RawData['f3'],ndmin=1)
StationLons=np.array(RawData['f4'],ndmin=1)
COPYStationIDs=np.array(RawData['f0'],ndmin=1)
COPYStationLats=np.array(RawData['f3'],ndmin=1)
COPYStationLons=np.array(RawData['f4'],ndmin=1)
nstations=len(StationIDs)


## Output the full matrix title of all station IDs
#if Restarter == '------':
#    My_Fhandle=file(OUTmatrix+".dat",'a')
#    goo=np.reshape(np.append("   STATIONID",StationIDs),(1,nstations+1))
#    np.savetxt(My_Fhandle,goo,fmt='%s',delimiter=' ')
#    My_Fhandle.close()

print(Restarter)
My_FGhandle=file(OUTGOODS,'a')
My_FBhandle=file(OUTBADS,'a')
# Get the distances for each file individually relative to all others
ngoods=nstations
for ss in range(nstations):
    #print(StationIDs[ss])
    if Restarter != '------' and Restarter != StationIDs[ss]:
        continue
    else: 
        Restarter='------'

    print(StationIDs[ss])
    # sort out RawData format for outputting station list
    outraw=list(RawData[ss])
    outraw[3]="{:8.4f}".format(outraw[3])
    outraw[4]="{:11.4f}".format(outraw[4])
    outraw[5]="{:9.2f}".format(outraw[5])
    outraw[6:12]=["{:5d}".format(dd) for dd in outraw[6:12]]
    outraw=np.reshape(outraw,(1,len(outraw)))     
	
    #print(StationIDs[ss])
# Create appropriate size arrays 
    StationDistsAll=np.zeros([ngoods])
    StatDistsSorted=np.zeros([ngoods])
    StatIDsSorted=np.empty([ngoods],dtype=object)	# object allows strings of any length and other types
    LocsMatch=[]
    #print(StationLats[ss],StationLons[ss])
    StationDistsAll,StatDistsSorted,StatIDsSorted,LocsMatch=GetDistanceArrays(StationIDs[ss],
                    StationDistsAll,StatDistsSorted,StatIDsSorted,StationLats[ss],
		    StationLons[ss],COPYStationLats,COPYStationLons,COPYStationIDs,ngoods,LocsMatch)
# If there is a LocsMatch value then remove this station from the list
    if len(LocsMatch) > 0:
        nbads=nbads+1
	ngoods=ngoods-1
	#findit=np.array([np.where(StationIDs == i) for i in LocsMatch]) # match multiple elements
        findit=np.where(COPYStationIDs == StationIDs[ss])[0] # match single elements
	COPYStationIDs=np.delete(COPYStationIDs,findit)
        COPYStationLats=np.delete(COPYStationLats,findit)
        COPYStationLons=np.delete(COPYStationLons,findit)	
	# output file to BAD list
        np.savetxt(My_FBhandle,outraw,fmt='%s',delimiter='')
	print("FOUND A LOCMATCH: ",ngoods,nbads,len(COPYStationIDs))
	
    else:
        # outpur file to GOOD List
        np.savetxt(My_FGhandle,outraw,fmt='%s',delimiter='')
   
	# Output the sorted arrays
	StCounts=100
	My_Fhandle=file(OUThundred,'a')
	WriteOutSorted(StatDistsSorted[0:100],StatIDsSorted[0:100],My_Fhandle,StationIDs[ss])
	My_Fhandle.close()

	StCounts=1000
	My_Fhandle=file(OUTthousand,'a')
	WriteOutSorted(StatDistsSorted[0:1000],StatIDsSorted[0:1000],My_Fhandle,StationIDs[ss])
	My_Fhandle.close()

        StCounts=40
        My_Fhandle=file(OUTneighbours,'a')
        WriteOut(StatIDsSorted[0:40],My_Fhandle,StationIDs[ss])
        My_Fhandle.close()

My_FBhandle.close()
My_FGhandle.close()

    #stop()

print("And, we are done!")
