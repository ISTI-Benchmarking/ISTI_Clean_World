# PYTHON2.7
# 
# Author: Kate Willett
# Created: 8 September 2015
# Last update: 8 September 2015
# Location: /home/h04/hadkw/Desktop/SurfaceTemperaturesWorkshop/BENCHMARKINGASSESS_GROUP/TEAM_CREATION/VARPROGS/	# this will probably change
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Reads in homogenised monthlies from PHA
# Mask FLs.r00 filled and adjusted output to the missing data of raw and output to:
#	/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/HOMOG/GHCN_TYPE and ISTI_TYPE and MASKED... all station file
# Copy over raw versions of those stations that were not homogenised
# Establish the list of neighbours from PHA corr. file
# Read in the raw time series and raw neighbour time series
# Plot the raw vs homog vs neighbours with median_pairwise trends (NEEDS WORK ON SIGNIFICANCE!)
#
# -----------------------
# LIST OF MODULES
# -----------------------
# import datetime as dt
# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.dates import date2num,num2date
# import sys, os
# from scipy.optimize import curve_fit,fsolve,leastsq
# from scipy import pi,sqrt,exp
# from scipy.special import erf
# import scipy.stats
# from math import sqrt,pi
# import struct
# import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)
# import subprocess as sp
#
# Kate's Modules
# from LinearTrends import MedianPairwise - Linear trend method,written by Kate Willett
#
# -----------------------
# DATA
# -----------------------
# Raw benchmark data in GHCN format from /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/monthly/raw/
# Raw benchmark data in ISTI format (for direct copying) /data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/ISTI_TYPE/
# Homogenised and infilled data in GHCN format from /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/monthly/FLs.r00/
# Complete station list from /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/meta/v101JUL2015_stnlist.tavg
# Homogenised station list from /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/corr/meta...
# Non-homogenised station lists from /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/corr/meta...input_not_stnlist
# Neighbour correlation file from /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/corr/corr.
#
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# python2.7 OutputPHAASCIIPLOT_SEP2015.py
#	Ensure file paths are set up correctly
# Set up ISTI/IMAGES/PHA_ADJ/ directory for plots
# Set up ISTI/DATA/CLEANWORLD/v101_JUL2015/HOMOG/GHCN_TYPE and ISTI_TYPE
# 
# -----------------------
# OUTPUT
# -----------------------
# Outputs to: /data/local/hadkw/ISTI/IMAGES/PHA_ADJ/
#	      /data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/HOMOG/
#			MASKEDFIXCLEAN1_ISTI_stage3proxyelevs_PHAADJ_SEP2015.txt
#			GHCN_TYPE/XXXXXXXXXXX.BNCHCAAA.tavg
#			ISTI_TYPE/BNCHCAAA_XXXXXXXXXXX_stage3
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 September 8th 2015
# ---------
#  
# Enhancements
#  
# Changes
#  
# Bug fixes
#  
# Version 2 (release date)
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
#########################################################################
# modules to import
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
import subprocess as sp
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

from LinearTrends import MedianPairwise

# RESTART VALUE
Restarter='------'				#'------'		#'681040'
Spin='TRUE'	#TRUE: loop through, FALSE: perform one stations only
Plotonly='FALSE'	#TRUE or FALSE
AddLetter='---'		#'---','a)'

# Set up initial run choices
ProjectName='BNCHCAAA'

# Set up file locations
STATSUFFIXOUT='_'+ProjectName+'.txt'

CORRFIL='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/BNCHCAAA_0915/corr/corr.BNCHCAAA_0915.tavg.r00.1509161733'	
INRAWISTI='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/ISTI_TYPE/'
INRAW='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/BNCHCAAA_0915/monthly/raw/'
STATSUFFIXINRAW='.raw.tavg'
STATLISTALL='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/BNCHCAAA_0915/meta/BNCHCAAA_0915_stnlist.tavg'	
STATLISTHOM='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/BNCHCAAA_0915/corr/meta.BNCHCAAA_0915.tavg.r00.1509161733'	
STATNOTLIST1='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/BNCHCAAA_0915/corr/meta.BNCHCAAA_0915.tavg.r00.1509161733.1.input_not_stnlist'	
STATNOTLIST2='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/BNCHCAAA_0915/corr/meta.BNCHCAAA_0915.tavg.r00.1509161733.2.input_not_stnlist'	
INHOM='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/BNCHCAAA_0915/monthly/FLs.r00/'
STATSUFFIXINHOM='.FLs.r00.tavg'
OUTPLOTDIR='/data/local/hadkw/ISTI/IMAGES/SEP2015/PHA_ADJ/'
OUTHOMMSK='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/HOMOG/MASKEDCLEAN_ISTI_stage3proxyelevs_BNCHCAAA_PHAADJ_SEP2015.txt'
OUTHOMGHCN='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/HOMOG/GHCN_TYPE/'
STATSUFFIXOUTGHCN='.'+ProjectName+'.tavg'
OUTHOMISTI='/data/local/hadkw/ISTI/DATA/CLEANWORLDS/v101_JUL2015/HOMOG/ISTI_TYPE/'
STATPREFIXOUTISTI=ProjectName+'_'


# Set up variables and arrays needed
mdi=-99.99
unit='$^{o}$C'

styr=1800
edyr=2015
DATASTART=dt.datetime(styr,1,1,0,0)
DATAEND=dt.datetime(edyr,12,1,0,0)
clmst=1976
clmed=2005
clmsty=(clmst-styr)
clmedy=(clmed-styr)
clmstm=(clmst-styr)*12
clmedm=((clmed-styr)*12)+11
CLIMSTART=dt.datetime(clmst,1,1,0,0)
CLIMEND=dt.datetime(clmed,12,1,0,0)
nmons=((edyr+1)-styr)*12
monarr=range(nmons)
nyrs=(edyr-styr)+1
yrarr=np.arange(nyrs)+styr

StationListALLID=[]	# nstations list filled after reading in station list for all
StationListHOMID=[]	# nstations list filled after reading in station list for homogenised
StationListNOTID=[]	# nstations list filled after reading in station list for not homogenised
nAstations=0	# defined after reading station list all
nHstations=0	# defined after reading station list homogenised
nNstations=0	# defined after reading station list not homogenised

nNeighbours=0	# defined after reading in corr station list
NeighbourList=[] # nNstations list filled after reading in corr station list

MyStation=[]	# filled after reading in candidate station
MyRAWStation=np.empty((nyrs,12))	# full time array filled after reading in candidate station
MyRAWStation.fill(-9999)
MyHOMStation=np.empty((nyrs,12))	# full time array filled after reading in candidate station
MyHOMStation.fill(-9999)
MyClims=[]	# 12 element array of mean months 1976-2005
MyAnomalies=[]	# filled with anomalies after subtracting climatology
MyHomogAnoms=[] # filled with homogenised anomalies
MyHomogAbs=[]	# filled with climatology+homogenised anomalies
MyClimMeanShift=[] # flat value across complete climatology period that the homogenised values differ from zero by - to rezero anoms and adjust clims/abs

NeighbourStations=[]	# nNstations by nmons array filled after reading in all neighbour stations
NeighbourAnomsStations=[]	# nNstations by nmons array filled after anomalising all neighbour stations relative to climatology
NeighbourClimsStations=[]	# nNstations by nmons array filled after anomalising all neighbour stations relative to climatology
NeighbourDiffStations=[]	# nNstations by nmons array filled after creating candidate minus neighbour difference series

MyFile=' '	#string containing file name

##########################################################################
# Subroutines
##########################################################################
#************************************************************************
# READDATA
def ReadData(FileName,typee,delimee):
    ''' Use numpy genfromtxt reading to read in all rows from a complex array '''
    ''' Need to specify format as it is complex '''
    ''' outputs an array of tuples that in turn need to be subscripted by their names defaults f0...f8 '''
    return np.genfromtxt(FileName, dtype=typee,delimiter=delimee) # ReadData

#************************************************************************
# FINDNEIGHBOURS
def FindNeighbours(FileName,CandID,neighbourcount,neighbourlist):
    ''' open the corr file and find the line beginning with the candidate station '''
    ''' list all neighbouring stations up to 40'''
    ''' be sure not to count 0s'''
    ''' return neighbour count and neighbour list '''
    for line in open(FileName):
        neighbourlist=[]			# make sure its blank to start
        neighbourlist=str.split(line)		# makes a list
        if neighbourlist[0] == CandID:	# found the line
            neighbourcount=len(neighbourlist)		# this doesn't include the zeros but does include the candidate in the count.
            break				# don't waste time, exit the loop

    return neighbourcount,neighbourlist # FindNeighbours
  
#************************************************************************
# READINNETWORKS
def ReadInNetworks(TheCount,TheList,TheCStation,TheFileDir,TheFileSuffix,StYear,YCount,TheData):
    ''' Loop through all neighbour station raw files '''
    ''' IGNORE FIRST FILE AS THIS IS THE CANDIDATE STATION '''
    ''' DOUBLE CHECK ALL OTHER STATIONS ARE NOT CANDIDATE AS THIS IS A KNOWN PROBLEM '''
    ''' read in using ReadStations and add to array '''    
    TheNewCount=0	# setting up new variables to output
    TheNewList=[]
    TheData=np.array(TheData)	# was an empty list
    for n,TheNStation in enumerate(TheList[1:]):	# 1: starts at second element
	if TheNStation == TheCStation:
	    continue
	    
        TheFile=TheFileDir+TheNStation+TheFileSuffix
        
        TempStation=np.empty((nyrs,12))	# full time array filled after reading in candidate station
        TempStation.fill(-9999)
        TheTypes=np.append("|S12",["int"]*13)
        TheDelimiters=np.append([12,4,6],[9]*11)
        RawData=ReadData(TheFile,TheTypes,TheDelimiters)
        for yy in range(len(RawData['f0'])):
	    moo=list(RawData[yy])
	    ypoint=moo[1]-StYear
	    TempStation[ypoint,:]=moo[2:14] 
        TempStation=np.reshape(TempStation.astype(np.float)/100.,(1,nmons)) 
        
	if TheData.size:		# if empty array then use first element, otherwise append
	    TheData=np.append(TheData,TempStation,axis=0) #np.reshape(TempStation/100.,(1,len(TempStation))),axis=0)	# now in proper units, fill the Neighbour array
	else:
	    TheData=TempStation  #np.reshape(TempStation/100.,(1,len(TempStation)))
        if any(TheNewList):		# if empty array then use first element, otherwise append
	    TheNewList=np.append(TheNewList,TheNStation)
	else:
	    TheNewList=[TheNStation]
    
    TheNewCount=len(TheNewList)		# Now this only includes the neighbours and not the candidate, as in FingNeighbours
    return TheData,TheNewList,TheNewCount #ReadInNetworks 	

#************************************************************************
# MAKEANOMALIES
def MakeAnomalies(TheData,TheAnomalies,TheClims,TheYCount,TheStClim,TheEdClim,TheMDI):
    ''' Working on both 1D and 2D (multiple station) arrays '''
    ''' Use given climatology period to create monthly clims and anomalies '''
    
    sizoo=TheData.shape			# returns a tuple of rows,columns
    TheClims=np.empty((sizoo[0],12))	# initialise clims array for nstations (rows) by 12 months (columns)
    TheClims.fill(TheMDI)
    TheAnomalies=np.empty(sizoo)
    TheAnomalies.fill(TheMDI)
    for t,TempStation in enumerate(TheData):	# row by row so ok as long as each station is a row
        #print(t,len(TempStation))
	Mooch=np.reshape(TempStation,(TheYCount,12))	# years(rows) by months(columns)
	Mooch2=np.empty_like(Mooch)		# To make sure I don't overwrite the absolute data
	Mooch2.fill(TheMDI)
	for mm in range(12):
	    subarr=Mooch[TheStClim:TheEdClim+1,mm]
	    #print(mm,subarr)
	    gots=(subarr > TheMDI)
	    if len(subarr[gots]) >= 15:		# more sophisticated checking has been done previously 
	        TheClims[t,mm]=np.mean(subarr[gots])
		gots2=(Mooch[:,mm] > TheMDI)
	        Mooch2[gots2,mm]=Mooch[gots2,mm]-TheClims[t,mm]
		#print " %6.2f"*40 % tuple(Mooch[:,mm])
	TheAnomalies[t,]=np.reshape(Mooch2,(1,12*TheYCount))    
    return TheAnomalies,TheClims #MakeAnomalies

#************************************************************************
# PLOTHOMOGTS
def PlotHomogTS(TheFile,TheStation,TheNeighbours,TheHStation,TheNCount,TheMDI,TheStYr,TheYCount,unit,typee,Letteree):
    ''' Plot raw candidate and neighbours with homogenised candidate '''
    ''' Add medianpairwise trends - from code medianpairwise.py '''
    '''MAKE MEDIANPAIRWISE.PY and COMPLETE WHEN HOMOG SERIES IS DONE '''
 
    # create annual averages and years and titles
    TheStationAnn=np.empty(TheYCount)
    TheStationAnn.fill(TheMDI)
    TheHStationAnn=np.empty(TheYCount)
    TheHStationAnn.fill(TheMDI)
    if TheNCount > 1:
        TheNeighboursAnn=np.empty((len(TheNeighbours[:,0]),TheYCount))
        TheNeighboursAnn.fill(TheMDI)

    TheStation=np.reshape(TheStation,(TheYCount,12))
    TheHStation=np.reshape(TheHStation,(TheYCount,12))    
    
    for yy in range(TheYCount):
        if np.sum(TheStation[yy,] != TheMDI) >= 9:
	    TheStationAnn[yy]=np.mean(TheStation[yy,np.where(TheStation[yy,] != TheMDI)])
        if np.sum(TheHStation[yy,] != TheMDI) >= 9:
	    TheHStationAnn[yy]=np.mean(TheHStation[yy,np.where(TheHStation[yy,] != TheMDI)])
 
    TheStation=np.reshape(TheStation,(TheYCount*12))
    TheHStation=np.reshape(TheHStation,(TheYCount*12))    
   
    if TheNCount > 1:
        for n,Neighbour in enumerate(TheNeighbours):
            Neighbour=np.reshape(Neighbour,(TheYCount,12))
            for yy in range(TheYCount):
                if np.sum(Neighbour[yy,] != TheMDI) >= 9:
	            TheNeighboursAnn[n,yy]=np.mean(Neighbour[yy,np.where(Neighbour[yy,] != TheMDI)])
        
    
    TheYears=np.reshape(range(TheStYr,TheStYr+TheYCount),TheYCount)
    ytitlee=typee+' ('+unit+')'
    xtitlee='Years'
    
    # get decadal trends and 5th-9th conf
    rawtrend=[0.,0.,0.]
    homtrend=[0.,0.,0.]
    rawtrend=MedianPairwise(TheStationAnn,TheMDI,rawtrend)
    homtrend=MedianPairwise(TheHStationAnn,TheMDI,homtrend)
    #pdb.set_trace()    
    # set up plot
 
    plt.clf()
    plt.figure(1,figsize=(8,4))
    plt.axes([0.1,0.1,0.85,0.80])
    if TheNCount > 1:
        PileItUp=np.append(TheNeighboursAnn,np.append(np.reshape(TheStationAnn,(1,TheYCount)),
             np.reshape(TheHStationAnn,(1,TheYCount)),axis=0),axis=0)
    else:
        PileItUp=np.append(np.reshape(TheStationAnn,(1,TheYCount)),
             np.reshape(TheHStationAnn,(1,TheYCount)),axis=0)
    
    plt.ylim([np.floor(min(PileItUp[PileItUp != TheMDI]))-2,
              np.ceil(max(PileItUp[PileItUp != TheMDI]))+2])
    plt.xlim([TheStYr,TheStYr+TheYCount])
    plt.tick_params(axis='both', which='major', labelsize=16)
   
    if TheNCount > 1:
        for n,Neighbour in enumerate(TheNeighboursAnn):
            line,=plt.plot(TheYears[np.where(Neighbour > TheMDI)],Neighbour[np.where(Neighbour > TheMDI)],color='black',linewidth=0.25)
 	
    line,=plt.plot(TheYears[np.where(TheStationAnn > TheMDI)],TheStationAnn[np.where(TheStationAnn > TheMDI)],'r',linewidth=2)	
    line,=plt.plot(TheYears[np.where(TheHStationAnn > TheMDI)],TheHStationAnn[np.where(TheHStationAnn > TheMDI)],'b',linewidth=2)
    if typee=='anomalies':
        line,=plt.plot(np.append(TheYears,TheStYr+TheYCount+1),np.zeros(TheYCount+1),'black',linewidth=1)        	
    
    plt.xlabel(xtitlee,size=16)
    plt.ylabel(ytitlee,size=16)
    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    
    rawstr="%5.2f +/- %5.2f to %5.2f %s /decade " % (rawtrend[0]*10,rawtrend[1]*10,rawtrend[2]*10,unit)
    homstr="%5.2f +/- %5.2f to %5.2f %s /decade " % (homtrend[0]*10,homtrend[1]*10,homtrend[2]*10,unit)

    plt.figtext(0.1,0.84,rawstr,color='r',size=16)
    plt.figtext(0.1,0.78,homstr,color='b',size=16)
    if Letteree != '---':
       plt.figtext(0.05,0.95,Letteree,color='Black',size=18)
       
    #plt.show()
    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotHomogTS
   
#***********************************************************************
#***********************************************************************
# WRITEGHCNSTATIONS
def WriteGHCNStations(TheID,TheStation,TheYears,TheFilee,TheMDI):
    ''' Output station data in GHCN format '''
    
    FindPresent=np.where(TheStation > TheMDI) # a 2 array list of (years, months) for each location where there is data
    
    # Use the first location and last location
    #pdb.set_trace()
    StartYear=FindPresent[0][0] 	# switch to start printing to file when data are present
    EndYear=FindPresent[0][len(FindPresent[0])-1] # switch to stop printing to file when data are no longer present
    #print(StartYear,EndYear)
    #pdb.set_trace()
    My_Fhandle=file(TheFilee,'a')
    for ss in range(StartYear,EndYear+1):
        #print(ss)
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
    
    FindPresent=np.where(TheStation > TheMDI) # a 2 array list of (years, months) for each location where there is data
    
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
########################################################################
# MAIN PROGRAM
########################################################################
# read in full station list
MyTypes=("|S11","float","float","int","|S4","|S30")
MyDelimiters=[11,7,10,13,4,30]
RawData=ReadData(STATLISTALL,MyTypes,MyDelimiters)
StationListALLID=np.array(RawData['f0'])
StatLats=np.array(RawData['f1'])
StatLons=np.array(RawData['f2'])
StatElevs=np.array(RawData['f3'])
StatNames=np.array(RawData['f5'])
nAstations=len(StationListALLID)

# read in homogenised station list
MyTypes=("|S11","float","float","int","|S34")
MyDelimiters=[11,7,10,13,34]
RawData=ReadData(STATLISTHOM,MyTypes,MyDelimiters)
StationListHOMID=np.array(RawData['f0'])
nHstations=len(StationListALLID)

# read in not homogenised station list 1
MyTypes=("|S11","float","float","int","|S34")
MyDelimiters=[11,7,10,13,34]
RawData=ReadData(STATNOTLIST1,MyTypes,MyDelimiters)
StationListNOTID=np.array(RawData['f0'])

# read in not homogenised station list 2
MyTypes=("|S11","float","float","int","|S34")
MyDelimiters=[11,7,10,13,34]
RawData=ReadData(STATNOTLIST2,MyTypes,MyDelimiters)
StationListNOTID=np.append(StationListNOTID,np.array(RawData['f0']))
nNstations=len(StationListNOTID)

# loop through station by station (full list)
for st in range(nAstations):

    # check if restart necessary
    if Restarter != '------' and Restarter != StationListALLID[st]:
        continue
    Restarter='------'

    # Set up empty file and then read in raw data file and populate 
    MyRAWStation=np.empty((nyrs,12))	# full time array filled after reading in candidate station
    MyRAWStation.fill(-9999)
    MyFile=INRAW+StationListALLID[st]+STATSUFFIXINRAW  
    MyTypes=np.append("|S12",["int"]*13)
    MyDelimiters=np.append([12,4,6],[9]*11)
    RawData=ReadData(MyFile,MyTypes,MyDelimiters)
    for yy in range(len(RawData['f0'])):
	moo=list(RawData[yy])
	ypoint=moo[1]-styr
	MyRAWStation[ypoint,:]=moo[2:14] 
    print(st,MyFile)  

    MyRAWStation=np.reshape(MyRAWStation.astype(np.float)/100.,(1,nmons))    # now in proper units and an array not list
    #pdb.set_trace() # reads in ok but array is 2d so must be Array[0] to access the columns
    
    # Is this station homogenised? 
    if (len(np.where(StationListNOTID == StationListALLID[st])[0]) > 0):
        # If no then: 
	# copy GHCN raw to HOMOG/GHCN_TYPE
	sp.call(['cp '+INRAW+StationListALLID[st]+STATSUFFIXINRAW+' '+OUTHOMGHCN+StationListALLID[st]+STATSUFFIXOUTGHCN], shell=True)
	
	# copy ISTI to ISTI_TYPE
	sp.call(['cp '+INRAWISTI+'merge_'+StationListALLID[st]+'_stage3 '+OUTHOMISTI+STATPREFIXOUTISTI+StationListALLID[st]+'_stage3'], shell=True)
		
	# Write line in MASKED file from 1800 onwards
        My_Fhandle=file(OUTHOMMSK,'a')
        goo=np.reshape(np.append(StationListALLID[st].strip(),["{:7.2f}".format(dd) for dd in MyRAWStation[0]]),(1,nmons+1))
        np.savetxt(My_Fhandle,goo,fmt='%s',delimiter=' ')
        My_Fhandle.close()
	
        # Continue with next station
        continue
	
    # If yes then carry on...
    	
    # set up clean arrays and variables
    nNeighbours=0	# defined after reading corr station list
    NeighbourList=[] # nNstations list filled after reading in corr station list

    MyClims=[]	# 12 element array of mean months 1976-2005
#    MyAnomalies=[] # filled with anomalies after subtracting climatology
#    MyHomogAnoms=[] # filled with homogenised anomalies
#    MyHomogAbs=[]  # filled with climatology+homogenised anomalies
#    MyClimMeanShift=[] # flat value across complete climatology period that the homogenised values differ from zero by - to rezero anoms and adjust clims/abs

    NeighbourStations=[]	# nNstations by nmons array filled after reading in all neighbour stations
    NeighbourAnomsStations=[]	# nNstations by nmons array filled after anomalising all neighbour stations relative to climatology
#    NeighbourClimsStations=[]	# nNstations by nmons array filled after anomalising all neighbour stations relative to climatology
#    NeighbourDiffStations=[]	# nNstations by nmons array filled after creating candidate minus neighbour difference series

    # set up empty, read in the PHA HOMOGENISED station file and populate
    MyHOMStation=np.empty((nyrs,12))	# full time array filled after reading in candidate station
    MyHOMStation.fill(-9999)
    MyFile=INHOM+StationListALLID[st]+STATSUFFIXINHOM 
    MyTypes=np.append(("|S12","int"),np.tile(("int","|S1","|S1","|S1"),12))
    MyDelimiters=np.append((12,4,6,1,1,1),(6,1,1,1)*11)
    RawData=ReadData(MyFile,MyTypes,MyDelimiters)
    for yy in range(len(RawData['f0'])):
	moo=list(RawData[yy])
	ypoint=moo[1]-styr
	# get the non-str bits of moo
	noomoo=np.reshape(moo[2:50],(12,4))
	MyHOMStation[ypoint,:]=noomoo[:,0] 

    MyHOMStation=np.reshape(MyHOMStation.astype(np.float)/100.,(1,nmons))    # now in proper units and an array not list
     
    #pdb.set_trace()          
    nNeighbours,NeighbourList=FindNeighbours(CORRFIL,StationListALLID[st],nNeighbours,NeighbourList)
    print("No. of Neighbours: ",nNeighbours-1)	# not including candidate but may have duplicate
    
    if nNeighbours > 1:
    # read in the neighbour files
        NeighbourStations,NeighbourList,nNeighbours=ReadInNetworks(nNeighbours,NeighbourList,
	                                           StationListALLID[st],INRAW,
						   STATSUFFIXINRAW,styr,nyrs,NeighbourStations)
        print("Actual No. of Neighbours: ",nNeighbours)	# not including candidate but may have duplicate

## convert all to anomalies (storing station climatology)
#    MyAnomalies,MyClims=MakeAnomalies(MyStation,MyAnomalies,MyClims,nyrs,clmsty,clmedy,mdi)
#	
#    NeighbourAnomsStations,NeighbourClimsStations=MakeAnomalies(NeighbourStations,NeighbourAnomsStations,
#	                                              NeighbourClimsStations,nyrs,clmsty,clmedy,mdi)
        
    # PLOT CANDIDATE AND NEIGHBOURS UNHOMOG WITH HOMOG ON TOP - ABS, ANOMS with MedianPairwiseTrends
    # REZEROD HOMOG MAY MEAN ITS NOW OFFSET COMPARED TO ORIGINAL
    MyPlotFile=OUTPLOTDIR+StationListALLID[st]+'_pha_adj_v101JUL2015_abs_'+ProjectName
    print(MyPlotFile)
    PlotHomogTS(MyPlotFile,MyRAWStation,NeighbourStations,MyHOMStation,nNeighbours,mdi,styr,nyrs,unit,'absolutes',AddLetter)
#    MyPlotFile=OUTPLOT+StationListWMO[st]+StationListWBAN[st]+'_trendcomp_7312OCT2013anoms'
#    PlotHomogTS(MyPlotFile,MyAnomalies,NeighbourAnomsStations,MyHomogAnoms,mdi,styr,nyrs,unit,'anomalies')

    # print out homogenised data, masked to raw because its infilled in more than just the unapplied adjustment places
    if Plotonly == 'FALSE':
        # Mask the homogenised data with missing data from raw
	#pdb.set_trace() # tested to check that this works
	MyHOMStation[0,np.where(MyRAWStation[0] == mdi)[0]]=mdi
	
	# Print out line to MASKED
        My_Fhandle=file(OUTHOMMSK,'a')
        goo=np.reshape(np.append(StationListALLID[st].strip(),["{:7.2f}".format(dd) for dd in MyHOMStation[0]]),(1,nmons+1))
        np.savetxt(My_Fhandle,goo,fmt='%s',delimiter=' ')
        My_Fhandle.close()
	
	# Reshape for printing out
	MyHOMStation=np.reshape(MyHOMStation[0],(nyrs,12))
	
	# Print out file to GHCN
        FilOuttee=OUTHOMGHCN+StationListALLID[st].strip()+STATSUFFIXOUTGHCN
        WriteGHCNStations(StationListALLID[st].strip(),MyHOMStation,yrarr,FilOuttee,mdi)
	
	# Print out file to ISTI
        FilOuttee=OUTHOMISTI+STATPREFIXOUTISTI+StationListALLID[st].strip()+'_stage3'
        WriteISTIStations(StatNames[st].strip(),StatLats[st],StatLons[st],StatElevs[st],MyHOMStation,yrarr,ProjectName,FilOuttee,mdi)
    
    # If we're looking at one file in isolation then stop
    if Spin == 'FALSE':
        break
# end loop of stations

#    stop()

print("And, we are done!")
