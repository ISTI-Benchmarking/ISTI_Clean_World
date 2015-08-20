#!/usr/local/sci/bin/python

#***************************************
# 06 June 2014 KMW - v1
# 
# Loop through ISTI station distances inventory
# calculate the correlation between each station and its closest 1000 stations
# output:
#       Closest 100 stations in order: stA st1 dist1 st2 dist2 etc
# 	Closest 1000 stations in order: stA st1 dist1 st2 dist2 etc
#
# Search for stations that either correlate at 1.0
# Output lists of highcorr, sameloc but do not remove - not a problem
# 
# 30 JUL2015
# changed 2014 to 2015 in GetStationArray
# 	fullyears=(2015+1)-1850
# added (timpointer >= 0) check to future proof ISTI predating 1850
#	if (TMeans[timee] > TheMDI) & (timpointer >= 0):
#
# 3rd August 2015
# Added in standardised anomalies
# At present all corrs are calculated based on climate anomalies
# Really we should be building VAR parameters on the standardised anomalies because that is what VAR will reproduce!!!
# Shouldn't be too wildly different but best to be clear
# So, this now:
#	: Removes climatology
#	: Fits Loess and removes - detrends
# 	: Devides by the standard deviation - standardises
# Outputs corrs based on climate anomalies AND corrs based on standardised anomalies
# Will then have to retest the distelev function build. Boo!
# 
# 3rd August 2015
# I really need to know how few data points are in a station so that
# the really data sparse ones aren't used to build the VARPS distelev function
# This is because they are highly likely to result in screwy means/standard deviataions
# and thus standardised anomalies
# I now output another list of each station and the number of datapoints within (different from first to last year)
# This can then be used by fixdistelev... program to only use those stations that have > 180 data points
# which should in theory be 15 years of data - of course they could well still be screwy!


#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 GetISTICorrs_FEB2015.py
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
import statsmodels.api as sm
from Geography import TwoPointDistanceKm
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

# RESTART VALUE
Restarter='------'	#'------'		#'SWE00138858' station ID

# Set up file locations
STATLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGDISTANCES_thousand_stage3proxyelevs_JUL2015.dat'
INFILEE='/data/local/hadkw/ISTI/DATA/ISTIv101_JUL2015/results_merged/merge_'	#_stage3
OUThundred='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGCORRS_hundred_stage3proxyelevs_JUL2015.dat'
OUThundredLag1='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGCORRSLAG1_hundred_stage3proxyelevs_JUL2015.dat'
#OUThundredLag1FD='/data/local/hadkw/ISTI/LISTS/ISTILONGCORRSLAG1firstdiff_hundred_stage3proxyelevs_FEB2015.dat'
OUThundredAC='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGAUTOCORRS_hundred_stage3proxyelevs_JUL2015.dat'
OUThundredC='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGCOMMONS_hundred_stage3proxyelevs_JUL2015.dat'
OUTBADS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGBADCORRS_stage3proxyelevs_JUL2015.dat'
OUTSAhundred='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGCORRSStAnom_hundred_stage3proxyelevs_JUL2015.dat'
OUTSAhundredLag1='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGCORRSLAG1StAnom_hundred_stage3proxyelevs_JUL2015.dat'
OUTSAhundredAC='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGAUTOCORRSStAnom_hundred_stage3proxyelevs_JUL2015.dat'
OUTLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_counts_stage3proxyelevs_JUL2015.dat'

# Set up variables and holders
MDI=-9999
nstations=0	# defined after reading in station list
StationPoints=[]	# number of datapoints present for that station
StationIDs=[]	# nstations list filled after reading in station list
StationCorrs=[]	# nstations list filled after reading in station list
StatDistsSorted=[]
StatDistsIDsSorted=[]
StatCorrsSorted=[]
StatCorrsLag1=[]
#StatCorrsLag1FD=[]
StatCorrsAC=[]
StatCommonsSorted=[]
StatCorrsIDsSorted=[]
StationCorrsSA=[]	# nstations list filled after reading in station list
StatDistsSortedSA=[]
StatDistsIDsSortedSA=[]
StatCorrsSortedSA=[]
StatCorrsLag1SA=[]
#StatCorrsLag1FD=[]
StatCorrsACSA=[]
StatCommonsSortedSA=[]
StatCorrsIDsSortedSA=[]

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
        return np.genfromtxt(FileName, dtype=typee,comments="%",delimiter=delimee,autostrip=ASTruth) # ReadData    
    else:
        return np.genfromtxt(FileName, dtype=typee,comments="%",delimiter=delimee,autostrip=ASTruth,usecols=ColumnChoice) # ReadData

#************************************************************************
# GETSTATIONARRAY
def GetStationArray(TheFilee,TheMDI):
    ''' Read in the data from the ISTI file format '''
    ''' Set up full array to include missing months of data '''
    ''' Calculate anomalies relative to the entire period of record '''
    ''' ANOMALY PERIOD MAY NOT BE THE SAME FOR EVERY STATION '''
    
# Read in the raw data    
    MyTypes=("|S18","|S14","float","float","float","|S1","int","int","|S2","int","int","int","|S71")
    MyDelimiters=[18,14,9,11,9,1,4,2,2,6,6,6,71]
    MyColumns=(6,7,11)
    RawData=ReadData(TheFilee,MyTypes,MyDelimiters,True,MyColumns)
    Years=np.array(RawData['f6'],ndmin=1)
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
	#print(timee,timpointer,Years[timee],Months[timee])
	if (TMeans[timee] > TheMDI) & (timpointer >= 0):
	    TheStationData[timpointer]=TMeans[timee]/100.

# Calculate anomalies relative to the entire period (not ideal but no common overlap period)    
    TheStationData=TheStationData.reshape(fullyears,12)
    for mm in range(12):
        gots=(TheStationData[:,mm] > TheMDI)
        #print(TheStationData[gots,mm])
	#print(np.mean(TheStationData[gots,mm]))
	TheStationData[gots,mm]=TheStationData[gots,mm]-np.mean(TheStationData[gots,mm])
        #print(TheStationData[gots,mm])
	#stop()

# Create Standardised anomalies too: remove lowess trend, divide by standard deviation
# Can only really do this if there is sufficient data - although can't really compare stats for short rubbish data so may not matter?
    # Temporarily change TheMDIs to NaN - 
# Make a copy for the standardised anomalies 
    TheStationDataSA=np.copy(TheStationData.reshape(fulltims))
    #stop
    TheStationDataSA[np.where(TheStationDataSA == TheMDI)[0]]=np.nan
    # Fit the lowess - using a fraction of 0.3 (1850 to 2014 is 1992 months, 600months/1992months=0.3)
    lowess = sm.nonparametric.lowess
    fraccie=600./fulltims
    Lowie=lowess(TheStationDataSA,range(fulltims),frac=fraccie,missing='drop')
    # Remove the lowess
    TheStationDataSA[np.where(np.isfinite(TheStationDataSA))[0]]=TheStationDataSA[np.where(np.isfinite(TheStationDataSA))[0]]-Lowie[:,1]
    # Change NaN back to TheMDI and then reshape to years,months
    TheStationDataSA[np.where(np.isnan(TheStationDataSA))[0]]=TheMDI
    TheStationDataSA=TheStationDataSA.reshape(fullyears,12)
    # Find climatological standard deviations and divide by them
    for mm in range(12):
        gots=(TheStationDataSA[:,mm] > TheMDI)
        #print(TheStationDataSA[gots,mm])
	#print(mm,np.std(TheStationDataSA[gots,mm]))
	TheStationDataSA[gots,mm]=TheStationDataSA[gots,mm]/np.std(TheStationDataSA[gots,mm])
        #print(TheStationDataSA[gots,mm])
	#stop()
    #stop	
    TheStationData=TheStationData.reshape(fulltims)	
    TheStationDataSA=TheStationDataSA.reshape(fulltims)	
    
    gots=np.where(TheStationData != TheMDI)[0]
    ThePoints=len(gots)
    
    #print(TheStationData[gots])
    
    return TheStationData, TheStationDataSA,ThePoints # GetStationArray
#************************************************************************
# GETCORRARRAYS
def GetCorrArrays(TheStationIDs,TheStationData,TheStationDataSA,TheMDI,TheNFilee):  #,TheType):
    ''' For each candidate neighbour pair: '''
    ''' 		Temporally match to only common time points - store'''
    ''' 		Create first difference series '''
    ''' 		Correlate '''
    ''' Sort the array of correlations and the IDs '''
    
# Create the necessary arrays
    TheCorrsList=np.zeros(200)	#1000
    TheSortedCorrs=np.zeros(200)
    TheCorrsLag1List=np.zeros(200)	#1000
    TheSortedCorrsLag1=np.zeros(200)
#    TheCorrsLag1FDList=np.zeros(200)	#1000
#    TheSortedCorrsLag1FD=np.zeros(200)
    TheAutoCorr=0	#1000
    TheCommonsList=np.zeros(200)
    TheSortedCorrsCommons=np.zeros(200)
    TheSortedIDs=np.empty(200,object)

# Create the necessary arrays for standardised anomalies
    TheCorrsListSA=np.zeros(200)	#1000
    TheSortedCorrsSA=np.zeros(200)
    TheCorrsLag1ListSA=np.zeros(200)	#1000
    TheSortedCorrsLag1SA=np.zeros(200)
    TheAutoCorrSA=0	#1000
#    TheCommonsListSA=np.zeros(200)
    TheSortedCorrsCommonsSA=np.zeros(200)
    TheSortedIDsSA=np.empty(200,object)
    TheStationIDsSA=np.copy(TheStationIDs)	# need own version for SA

# Loop through each neighbour station
    for s,ss in enumerate(TheStationIDs):

# Read in the neighbour station and convert to anomalies
        NFilee=TheNFilee+ss+'_stage3'	# ss should be the string with no white space
	#print(s,NFilee)
        TheNeighbourData,TheNeighbourDataSA,ThePoints=GetStationArray(NFilee,TheMDI)
	 
## If TheType == 'Standardised' then use standardised anomalies
#        if (TheType == 'Standardised'):
#	    TheNeighbourData=TheNeighbourDataSA	 
	 
# Pattern match the candidate and neighbour station
        gots=(TheStationData == TheMDI)   
        TheNeighbourData[gots]=TheMDI	 
	gots=(TheNeighbourData == TheMDI)
	TheStationDataCopy=np.array(TheStationData)
	TheStationDataCopy[gots]=TheMDI

        gots=(TheStationDataSA == TheMDI)   
        TheNeighbourDataSA[gots]=TheMDI	 
	gots=(TheNeighbourDataSA == TheMDI)
	TheStationDataCopySA=np.array(TheStationDataSA)
	TheStationDataCopySA[gots]=TheMDI
	 
# Create first difference series
        gots=(TheStationDataCopy > TheMDI)
	commons=len(TheStationDataCopy[gots])
	#print('COMMONS: ',commons)
	TheCommonsList[s]=commons
	if commons > 12:	# must be at least a year to correlate with
	    TheStationComps=np.array(TheStationDataCopy[gots])
	    TheNeighbourComps=np.array(TheNeighbourData[gots])
	    TheStationDiffs=TheStationComps[1:commons]-TheStationComps[0:commons-1]
	    TheNeighbourDiffs=TheNeighbourComps[1:commons]-TheNeighbourComps[0:commons-1]

	    TheStationCompsSA=np.array(TheStationDataCopySA[gots])
	    TheNeighbourCompsSA=np.array(TheNeighbourDataSA[gots])
	    TheStationDiffsSA=TheStationCompsSA[1:commons]-TheStationCompsSA[0:commons-1]
	    TheNeighbourDiffsSA=TheNeighbourCompsSA[1:commons]-TheNeighbourCompsSA[0:commons-1]

# Correlate and populate TheCorrList
            TheCorrsList[s]=np.corrcoef(TheStationDiffs,TheNeighbourDiffs)[0, 1]
            TheCorrsListSA[s]=np.corrcoef(TheStationDiffsSA,TheNeighbourDiffsSA)[0, 1]
	    #print('CORRS diff, comp: ',np.corrcoef(TheStationDiffs,TheNeighbourDiffs)[0, 1],np.corrcoef(TheStationComps,TheNeighbourComps)[0, 1])
            # NOTE CORRS ARE LARGER FOR COMPS THAN DIFFS - NOT MUCH IN IT THOUGH
	    #stop()
	    # No check for consecutive streak of data here so missing data will
	    # make the AC lower than it should be
#            TheCorrsLag1FDList[s]=np.corrcoef(TheStationDiffs[0:commons-2],TheNeighbourDiffs[1:commons-1])[0, 1]
            TheCorrsLag1List[s]=np.corrcoef(TheStationComps[0:commons-1],TheNeighbourComps[1:commons])[0, 1]
            TheCorrsLag1ListSA[s]=np.corrcoef(TheStationCompsSA[0:commons-1],TheNeighbourCompsSA[1:commons])[0, 1]
	else:
            TheCorrsList[s]=0.0
            TheCorrsLag1List[s]=0.0
            TheCorrsListSA[s]=0.0
            TheCorrsLag1ListSA[s]=0.0
#            TheCorrsLag1FDList[s]=0.0
	     

# Sort the Corrs and station IDs
        SortIndex=np.argsort(TheCorrsList)
	SortIndex=SortIndex[::-1]
        TheSortedCorrs=TheCorrsList[SortIndex]
        TheSortedCorrsCommons=TheCommonsList[SortIndex]
        TheSortedIDs=TheStationIDs[SortIndex]

        SortIndexSA=np.argsort(TheCorrsListSA)
	SortIndexSA=SortIndexSA[::-1]
        TheSortedCorrsSA=TheCorrsListSA[SortIndexSA]
        TheSortedIDsSA=TheStationIDsSA[SortIndexSA]
	
	# Is there a 0.0 distance? A Location Match!
	# If so - flag - this station will then be listed with its duplicate(s)
	AnyBads=[]
	AnyBads=np.where(np.array(TheSortedCorrs) == 1.0)[0]
	
        TheSortedCorrsLag1=TheCorrsLag1List[SortIndex]
        TheSortedCorrsLag1SA=TheCorrsLag1ListSA[SortIndexSA]
#        TheSortedCorrsLag1FD=TheCorrsLag1FDList[SortIndex]

	gots=np.where(TheStationData != TheMDI)[0]  
	# Find the longest section of continuous data - hopefully at least 4 values
	#print(gots)
	diffsies=gots[1:len(gots)]-gots[0:len(gots)-1]
        biggs=np.where(diffsies > 1)[0]
        countbiggs=len(biggs)
        if countbiggs > 0:
            newbiggs=np.concatenate((np.array([0]),np.transpose(biggs),np.array([len(gots)-1])),0)
	    #print(newbiggs)
            biggdiffs=np.array(newbiggs[1:len(newbiggs)]-newbiggs[0:len(newbiggs)-1]) # Rememeber you've added two extra elements so countbiggs+2 is correct
            #print(biggdiffs)
	    gotcha=np.where(biggdiffs == max(biggdiffs))[0]
	    #print(gotcha)
# need to add 1 to beginning of biggs if record is continuous from the beginning?
            contgots=gots[newbiggs[gotcha[0]]+1:newbiggs[gotcha[0]+1]+1]
	    #print(contgots)
        else: 
            contgots=gots
	    #print("No Gaps: ",contgots)
    	 
	if len(contgots) >= 4:
	    TheAutoCorr=np.corrcoef(TheStationData[contgots[0:len(contgots)-1]],TheStationData[contgots[1:len(contgots)]])[0,1]
	    TheAutoCorrSA=np.corrcoef(TheStationDataSA[contgots[0:len(contgots)-1]],TheStationDataSA[contgots[1:len(contgots)]])[0,1]
	else: 
	    TheAutoCorr=0.0
	    TheAutoCorrSA=0.0
	    
    return TheSortedCorrs,TheSortedIDs,TheSortedCorrsCommons,TheSortedCorrsLag1,TheAutoCorr,AnyBads,TheSortedCorrsSA,TheSortedIDsSA,TheSortedCorrsLag1SA,TheAutoCorrSA # GetCorrArrays  TheSortedCorrsLag1FD,
  
#************************************************************************
# WRITEOUTSORTED
def WriteOutSorted(TheSTCount,TheSortedCorrs,TheSortedIDs,TheFile,TheCandidateID):
    ''' Output lines to text of StationID, list of stations and distances '''
    
# Convert all distances to set length strings
    #print(TheStDists[0:10])
    #TheStDistsStr=np.array(["{:9.3f}".format(dd) for dd in TheStDists.reshape(TheStDists.size)])
    #TheStDistsStr=TheStDistsStr.reshape(TheStDists.shape)
    if TheSTCount > 1:
        TheSortedCorrsStr=["{:12.3f}".format(dd) for dd in TheSortedCorrs]
        # Make a nstations (rows) by 2 column array, reform such that it becomes r1c1,r2c1,r1c2,r2c2,r1c3,r2c3 etc
        TheData=np.reshape(zip(*np.vstack((TheSortedIDs,TheSortedCorrsStr))),len(TheSortedIDs)*2)	# a one by nstations array 
        goo=np.reshape(np.append(TheCandidateID,TheData),(1,(len(TheSortedIDs)*2)+1))
    else:
        TheSortedCorrsStr="{:12.3f}".format(TheSortedCorrs[0])
        goo=np.reshape(np.append(TheCandidateID,TheSortedCorrsStr),(1,2))
     
    
    np.savetxt(TheFile,goo,fmt='%s',delimiter=' ')
    return #WriteOutSorted
    
#***********************************************************************
# MAIN PROGRAM
#***********************************************************************
# read in Distance lists station list
a=range(1000)
MyTypes=np.append("|S12",np.reshape(np.array([["|S13","float"] for i in a]),2000))
MyDelimiters=np.append(12,[13]*2000)
MyColumns=0
RawDists=ReadData(STATLIST,MyTypes,MyDelimiters,True,MyColumns)
StationIDs=np.array(RawDists['f0'],ndmin=1)
nstations=len(StationIDs)
#print(StationIDs[0:10])
#stop()

print('Read in giant file')
# read in the distances/IDs for each station file
f=open(STATLIST,'r')
for lop in range(nstations):
    line=f.readline()		# moo=list(RawDists[lop])

    if Restarter != '------' and Restarter != StationIDs[lop]:
        continue
    else: 
        Restarter='------'
    
    line=line.split()
    #print(line[0:10])
    MakeItStrings=np.reshape(np.array(line[1:2001],dtype="|S13"),(1000,2))
    #print(MakeItStrings[0:10,0])
    StatDistsIDsSorted=np.array(MakeItStrings[0:200,0])	# normally just : or 0:1000 for printing out 1000 corrs
    #print(StatDistsIDsSorted[0:10])
#    print('Read in epic list of distance/IDs')
    
# read in the station data for the candidate file
    MyFile=INFILEE+StationIDs[lop]+'_stage3'
    CandidateStation,CandidateStationSA,StationPoints=GetStationArray(MyFile,MDI)
#    print('Opened the station file')

# get the correlation list (sorted) from all neighbours
    DupsFound=[]
#    AnomType='Climate'
    StatCorrsSorted,StatCorrsIDsSorted,StatCommonsSorted,StatCorrsLag1,StatCorrsAC,DupsFound,StatCorrsSortedSA,StatCorrsIDsSortedSA,StatCorrsLag1SA,StatCorrsACSA=GetCorrArrays(StatDistsIDsSorted,CandidateStation,CandidateStationSA,MDI,INFILEE) #StatCorrsLag1FD, 
#    print('Done the corrs',StatCorrsAC)  
    print(StationIDs[lop],StatCorrsSorted[0])

## get the correlation list (sorted) from all neighbours using STANDARDISED ANOMALIES
#    DupsFoundSA=[]
#    AnomType='Standardisd'
#    StatCorrsSortedSA,StatCorrsIDsSortedSA,StatCommonsSortedSA,StatCorrsLag1SA,StatCorrsACSA,DupsFoundSA=GetCorrArrays(StatDistsIDsSorted,CandidateStationSA,MDI,INFILEE,AnomType) #StatCorrsLag1FD, 
##    print('Done the corrs',StatCorrsAC)  
#    print(StationIDs[lop],StatCorrsSortedSA[0])

# Write out the station and the number of time points it contains
    My_Fhandle=file(OUTLIST,'a')
    np.savetxt(My_Fhandle,np.reshape(np.append(StationIDs[lop],"{:5d}".format(StationPoints)),(1,2)),fmt='%s',delimiter=' ')
    My_Fhandle.close()
   
# write out the Standardised Anomalies sorted correlation list to files    
    StCounts=100
    My_Fhandle=file(OUTSAhundred,'a')
    WriteOutSorted(StCounts,StatCorrsSortedSA[0:100],StatCorrsIDsSortedSA[0:100],My_Fhandle,StationIDs[lop])
    My_Fhandle.close()

    StCounts=100
    My_Fhandle=file(OUTSAhundredLag1,'a')
    WriteOutSorted(StCounts,StatCorrsLag1SA[0:100],StatCorrsIDsSortedSA[0:100],My_Fhandle,StationIDs[lop])
    My_Fhandle.close()

    StCounts=1
    My_Fhandle=file(OUTSAhundredAC,'a')
    WriteOutSorted(StCounts,np.array([StatCorrsACSA]),StationIDs[lop],My_Fhandle,StationIDs[lop])
    My_Fhandle.close()
	
# write out the sorted correlation list to files    
    StCounts=100
    My_Fhandle=file(OUThundred,'a')
    WriteOutSorted(StCounts,StatCorrsSorted[0:100],StatCorrsIDsSorted[0:100],My_Fhandle,StationIDs[lop])
    My_Fhandle.close()

    StCounts=100
    My_Fhandle=file(OUThundredLag1,'a')
    WriteOutSorted(StCounts,StatCorrsLag1[0:100],StatCorrsIDsSorted[0:100],My_Fhandle,StationIDs[lop])
    My_Fhandle.close()

#    StCounts=100
#    My_Fhandle=file(OUThundredLag1FD,'a')
#    WriteOutSorted(StCounts,StatCorrsLag1FD[0:100],StatCorrsIDsSorted[0:100],My_Fhandle,StationIDs[lop])
#    My_Fhandle.close()

    StCounts=1
    My_Fhandle=file(OUThundredAC,'a')
    WriteOutSorted(StCounts,np.array([StatCorrsAC]),StationIDs[lop],My_Fhandle,StationIDs[lop])
    My_Fhandle.close()

    if len(DupsFound) > 0:
        StCounts=len(DupsFound)+1
	My_Fhandle=file(OUTBADS,'a')
        WriteOutSorted(StCounts,np.array(StatCorrsSorted[DupsFound]),StatCorrsIDsSorted[0:len(DupsFound)],My_Fhandle,StationIDs[lop])
        My_Fhandle.close()

#    StCounts=1000
#    My_Fhandle=file(OUTthousand,'a')
#    WriteOutSorted(StCounts,StatCorrsSorted[0:1000],StatCorrsIDsSorted[0:1000],My_Fhandle,StationIDs[lop])
#    My_Fhandle.close()

# write out the sorted commons list to files    
    StCounts=100
    My_Fhandle=file(OUThundredC,'a')
    WriteOutSorted(StCounts,StatCommonsSorted[0:100],StatCorrsIDsSorted[0:100],My_Fhandle,StationIDs[lop])
    My_Fhandle.close()

#    StCounts=1000
#    My_Fhandle=file(OUTthousandC,'a')
#    WriteOutSorted(StCounts,StatCommonsSorted[0:1000],StatCorrsIDsSorted[0:1000],My_Fhandle,StationIDs[lop])
#    My_Fhandle.close()

    #stop
f.close()
stop()

print("And, we are done!")
