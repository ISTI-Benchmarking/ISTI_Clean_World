#!/usr/local/sci/bin/python

#***************************************
# 06 June 2014 KMW - v1
# 
# Reads in clean world ISTI data file and outputs a file for each station
# in GHCN format
# Additionally reads in the ISTI station list file and outputs in 
# GHCN format
# Additionally outputs clean world ISTI into ISTI format

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 ConvertISTItoGHCN_NOV2014.py
#
# JUL 2015
# This has been adapted to output my single ASCII format (stations in rows, months in columns) to ISTI and GHCN format
# Could later be adapted to output netCDF too
#
# AUG 2015
# Added the missing data mask from the real data - so we now output only the missing data versions in ISTI and GHCN format
# BUT also output a missing data version in my singla ASCII format
#
# SEP 2015
# Now only prints data from the first year of data in that station to the last year of data in that station to save space, and PHA processing
# Also updated so that the GHCN list is correct. Fudged country code by using first two letters of Station ID just to fill space.
#
# REQUIRES
#************************************************************************
# Set up python imports
import datetime as dt
import numpy as np
import sys, os
from scipy.optimize import curve_fit,fsolve,leastsq
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.stats
from math import sqrt,pi
import struct
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

# RESTART VALUE
Restarter='------'	#'------'		#'    01010000' station ID

# Set up file locations
STATLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat'
INFILEE='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/CLEAN_ISTI_stage3proxyelevs_BNCHCAAA_SEP2015.txt'	#_stage3
INREALS='/data/local/hadkw/ISTI/DATA/ISTIv101_JUL2015/results_merged/merge_'
OUTSTATLIST='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/ISTILONGINVENTORYghcn_stage3proxyelevs_SEP2015.dat'
OUTFILEEGHCN='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/GHCN_TYPE/' # add ID number here
OUTFILEEISTI='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/ISTI_TYPE/' # add ID number here
OUTFILEEMASK='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/MASKEDCLEAN_ISTI_stage3proxyelevs_BNCHCAAA_SEP2015.txt'

# Set up variables and holders
MDI=-9999
#RStYr=1800	# Actual start year
#REdYr=2018	# Actual end year
#RYrCount=(REdYr-RStYr)+1
#RMonCount=RYrCount*12
#RYrArr=np.array(range(RYrCount))+RStYr
StYr=1800	# Simulated start year
EdYr=2015	# Simulated end year
YrCount=(EdYr-StYr)+1
MonCount=YrCount*12
YrArr=np.array(range(YrCount))+StYr
nstations=0	# defined after reading in station list
StatIDs=[]	# nstations list filled after reading in station list
StatLats=[]	# nstations list filled after reading in station list
StatLons=[]
StatElevs=[]
StatCOs=[]
StatNames=[]
Benchmark_ID='BCHMCAAA'	# Benchmark Clean World AAA

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
# WRITEGHCNLIST
def WriteGHCNList(TheIDs,TheLats,TheLons,TheElevs,TheCOs,TheNames,TheFilee,TheSCount):
    ''' Output station lat lon elev co name for each station in GHCN format '''

    My_Fhandle=file(TheFilee,'a')
    for ss in range(TheSCount):
        outlist=list(("{:11s}".format(TheIDs[ss].strip()),
	             "{:7.2f}".format(TheLats[ss]),
		     "{:10.2f}".format(TheLons[ss]),
		     "{:13d}".format(int(round(TheElevs[ss]))),
		     " ",
		     "{:2s}".format(TheIDs[ss][1:3]), #format(TheCOs[ss][1:3])
		     " ",
		     "{:30s}".format(TheNames[ss][1:31])))
 
        # Add this to the noelevs list
        outlist=np.reshape(outlist,(1,len(outlist)))     
        np.savetxt(My_Fhandle,outlist,fmt='%s',delimiter='')	    

    My_Fhandle.close()
    
    return #WriteGHCNList
    
#***********************************************************************
# WRITEGHCNSTATIONS
def WriteGHCNStations(TheID,TheStation,TheYears,TheFilee,TheMDI):
    ''' Output station data in GHCN format '''
    
    FindPresent=np.where(TheStation > TheMDI/100.) # a 2 array list of (years, months) for each location where there is data
    
    # Use the first location and last location
    StartYear=FindPresent[0][0] 	# switch to start printing to file when data are present
    EndYear=FindPresent[0][len(FindPresent[0])-1] # switch to stop printing to file when data are no longer present
    #print(StartYear,EndYear)
    #pdb.set_trace()
    My_Fhandle=file(TheFilee,'a')
    for ss in range(StartYear,EndYear+1):

        outlist=list(("{:11s}".format(TheID),
                   " ",
		   "{:4d}".format(TheYears[ss]),
                   "{:6d}".format(int(100*float(TheStation[ss,0]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,1]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,2]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,3]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,4]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,5]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,6]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,7]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,8]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,9]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,10]))),
		   "   ",
                   "{:6d}".format(int(100*float(TheStation[ss,11]))),
		   "   "))
    
        outlist=np.reshape(outlist,(1,len(outlist)))     
        np.savetxt(My_Fhandle,outlist,fmt='%s',delimiter='')	    

    My_Fhandle.close()

    return #WriteGHCNStations
    
#***********************************************************************
# WRITEISTISTATIONS
def WriteISTIStations(TheName,TheLat,TheLon,TheElev,TheStation,TheYears,TheBCHM,TheFilee,TheMDI):
    ''' Output station data in ISTI format '''
    ''' 30 Name, '''
    ''' 11.4 latitude, '''
    ''' 11.4 longitude, '''
    ''' 9.2.elevation, '''
    ''' x, '''
    ''' 4 year, ''' 
    ''' 2 month, '''
    ''' 'XX', '''
    ''' 6 tmax, ''' 
    ''' 6 tmin, '''
    ''' 6 tavg, ''' 
    ''' x, '''
    ''' '999/999/999/999/999/999/999/999/999/999/999 XXXXXXXX/XXXXXXXX/BCHMCAAA' '''
    
    FindPresent=np.where(TheStation > TheMDI/100.) # a 2 array list of (years, months) for each location where there is data
    
    # Use the first location and last location
    StartYear=FindPresent[0][0] 	# switch to start printing to file when data are present
    EndYear=FindPresent[0][len(FindPresent[0])-1] # switch to stop printing to file when data are no longer present
    #print(StartYear,EndYear)
   
    TheMon=['01','02','03','04','05','06','07','08','09','10','11','12']
    My_Fhandle=file(TheFilee,'a')
    for ss in range(StartYear,EndYear+1):
        for mm in range(12):
            outlist=list(("{:30s}".format(TheName),
		   "{:11.4f}".format(TheLat),
		   "{:11.4f}".format(TheLon),
		   "{:9.2f}".format(TheElev),
                   " ",
		   "{:4d}".format(TheYears[ss]),
		   "{:2s}".format(TheMon[mm]),
                   "XX -9999 -9999",
                   "{:6d}".format(int(100*float(TheStation[ss,mm]))),
		   " ",
		   "999/999/999/999/999/999/999/999/999/999/999 XXXXXXXX/XXXXXXXX/",
		   "{:8s}".format(TheBCHM)))
    
            outlist=np.reshape(outlist,(1,len(outlist)))  
	    #print('ISTI',outlist)   
            np.savetxt(My_Fhandle,outlist,fmt='%s',delimiter='')	    

    My_Fhandle.close()

    return #WriteISTIStations
    
#***********************************************************************
# MAIN PROGRAM
#***********************************************************************
# read in station list
MyTypes=("|S12","|S31","|S24","float","float","float","int","int","int","int","int","int","|S16","|S17")
MyDelimiters=[12,31,24,8,11,9,5,5,5,5,5,5,16,17]
MyColumns='XXX'
RawData=ReadData(STATLIST,MyTypes,MyDelimiters,False,MyColumns)
StatIDs=np.array(RawData['f0'])
StatLats=np.array(RawData['f3'])
StatLons=np.array(RawData['f4'])
StatElevs=np.array(RawData['f5'])
StatCOs=np.array(RawData['f2'])
StatNames=np.array(RawData['f1'])
nstations=len(np.array(RawData['f0']))
print('Read in Giant file...')

# read in the data for each station line by line
f=open(INFILEE,'r')
for lop in range(nstations):
    statty=f.readline()		# moo=list(RawDists[lop])#

    if Restarter != '------' and Restarter != StatIDs[lop]:
        continue
    else: 
        Restarter='------'
    
    statty=statty.split()
    
    print(StatIDs[lop],statty[0])
    tmpdata=statty[1:len(statty)]
    
    # read in real data and mask
    masked_data=np.empty(MonCount,dtype=float)
    masked_data.fill(MDI/100.)
    MyTypes=("|S19","|S13","float","float","float","|S1","int","int","|S2","int","int","int","|S71")
    MyDelimiters=[19,13,9,11,9,1,4,2,2,6,6,6,71]
    MyColumns=(6,7,11)
    TheFilee=INREALS+StatIDs[lop].strip()+'_stage3'
    RealRawData=ReadData(TheFilee,MyTypes,MyDelimiters,True,MyColumns)
    Years=np.array(RealRawData['f6'],ndmin=1)
    Months=np.array(RealRawData['f7'],ndmin=1)
    TMeans=np.array(RealRawData['f11'],ndmin=1)
    ntims=len(Years)
# Fill in the simulated data where real data are NOT missing from 1860 onwards (end of 2014 at present) (ignore time before that)
    #stop
    for timee in range(ntims):
        timpointer=((Years[timee]-StYr)*12)+Months[timee]-1
	if (TMeans[timee] > MDI) & (timpointer >= 0):
	    masked_data[timpointer]=np.float(tmpdata[timpointer])

    My_Fhandle=file(OUTFILEEMASK,'a')
    outraw=np.append(StatIDs[lop].strip(),["{:7.2f}".format(dd) for dd in masked_data])
    #stop
    outraw=np.reshape(outraw,(1,len(outraw)))     
    np.savetxt(My_Fhandle,outraw,fmt='%s',delimiter='')
    My_Fhandle.close()    
    
    masked_data=np.reshape(masked_data,(YrCount,12))
    #masked_data[masked_data == 'NA']=-99.99
    
    FilOuttee=OUTFILEEGHCN+StatIDs[lop].strip()+'.raw.tavg'
    WriteGHCNStations(StatIDs[lop].strip(),masked_data,YrArr,FilOuttee,MDI)
    
    FilOuttee=OUTFILEEISTI+'merge_'+StatIDs[lop].strip()+'_stage3'
    WriteISTIStations(StatNames[lop].strip(),StatLats[lop],StatLons[lop],StatElevs[lop],masked_data,YrArr,Benchmark_ID,FilOuttee,MDI)

    
WriteGHCNList(StatIDs,StatLats,StatLons,StatElevs,StatCOs,StatNames,OUTSTATLIST,nstations)
   
f.close()
stop()

print("And, we are done!")
