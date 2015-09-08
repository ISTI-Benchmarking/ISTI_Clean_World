# PYTHON2.7
# 
# Author: Kate Willett
# Created: 7 September 2015
# Last update: 7 September 2015
# Location: /home/h04/hadkw/Desktop/SurfaceTemperaturesWorkshop/BENCHMARKINGASSESS_GROUP/TEAM_CREATION/VARPROGS/	# this will probably change
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# Code reads adjustment locations, magnitudes and uncertainties from the PHA log.
# It looks at both ALL adjustments, and only the 'significant' ones - where UNC < ADJ
# Various statistics and graphics are output:
#	A list of only stations with significant adjustments (unc < adj)
# 	List of adjustments, largest to smallest
# 	Average number of changepoints per station
#	Mean (median) absolute size and standard deviation of adjustment
# 	Mean (median) actual size and standard deviation of adjustment
# 	Time series of changepoint frequency for all stations combined
#	Histogram of adjustment magnitude (could add gaussian fit and missed adj uncertainty?
#
# 	Could add:
#		Gaussian fit and missed adjustment uncertainty estimate on histogram
#		
#
# -----------------------
# LIST OF MODULES
# -----------------------
# import matplotlib.pyplot as plt
# import numpy as np
# import sys, os
# import datetime as dt
# from matplotlib.dates import date2num,num2date
# from RandomsRanges import LetterRange
# import re	# regular expression stuff for character replacement within strings
# import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)
# import subprocess as sp
# 
# Kate's Modules
# from RandomsRanges import LetterRange - written by Kate Willett to provide a list of letters from A to Z 
#
# -----------------------
# DATA
# -----------------------
# PHA logfile from Adj write: onwards compress to just Adj write lines by grep ^"Adj write:" BNCHCAAA_v101JUL2015_PHA.log > BNCHCAAA_v101JUL2015_PHA_SEP2015.log
# /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/output/BNCHCAAA_v101JUL2015_PHA_SEP2015.log
# 
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# python2.7 PlotPHAAdjsStats_SEP2015.py
#	Ensure file paths are set up correctly
# Prepare input files first:
# PHA homogenised stations = 31410 (out of 32522)
# /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/corr/meta.v101JUL2015.tavg.r00.1509041626
# Reduce PHAv52j.FAST.MLY.TEST.1509041626.tavg.v101JUL2015.r00.out to only last Adj write: part
# /data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/output/BNCHCAAA_v101JUL2015_PHA.log
# 
# -----------------------
# OUTPUT
# -----------------------
# Outputs to: /data/local/hadkw/ISTI/IMAGES/
#		PHAAdj_stats_ALL_BNCHCAAA_SEP2015.eps/png
#		PHAAdjs_stats_SIG_BNCHCAAA_SEP2015.eps/png
#	    /data/local/hadkw/ISTI/LISTS/
#		Largest_PHAAdjs_ALL_BNCHCAAA_SEP2015.txt
#		Largest_PHAAdjs_SIG_BNCHCAAA_SEP2015.txt
#		PHAAdjustmentlist_ALL_BNCHCAAA_SEP2015.txt
#		PHAAdjustmentlist_SIG_BNCHCAAA_SEP2015.txt
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 September 7th 2015
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
# BNCHCAAA Notes ALL, SIG
# Adj per station (32522 stations from 1800 to 2015, v101JUL2015 missing data): 0.141, 0.073
# Mean absolute adjustment size:0.410, 0.445
# Standard deviation of absolute adjustment size: 0.302, 0.319
# Mean actual adjustment size: -0.003, 0.007
# Standard deviation of actual adjustments: 0.509, 0.548
# 
#########################################################################
# modules to import
import matplotlib.pyplot as plt
import numpy as np
import sys, os
import datetime as dt
from matplotlib.dates import date2num,num2date
from RandomsRanges import LetterRange
import re	# regular expression stuff for character replacement within strings
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)
import subprocess as sp

# Set up parameters
ProjectName='BNCHCAAA'

# directories and file names
# PHA homogenised stations = 31410 (out of 32522)
instats='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/corr/meta.v101JUL2015.tavg.r00.1509041626'
# reduced PHAv52j.FAST.MLY.TEST.1509041626.tavg.v101JUL2015.r00.out to only last Adj write: part
inlog='/data/local/hadkw/ISTI/PROGS/PHA2014v52j/PHA52j_full/pha_v52j/data/isti/v101JUL2015/output/BNCHCAAA_v101JUL2015_PHA_SEP2015.log'
outplotdir='/data/local/hadkw/ISTI/IMAGES/'
outlistdir='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/'
outplotstatALL='PHAAdj_stats_ALL_BNCHCAAA_SEP2015'
outplotstatSIG='PHAAdj_stats_SIG_BNCHCAAA_SEP2015'
outlistlargeALL='Largest_PHAAdjs_ALL_BNCHCAAA_SEP2015.txt'
#outlistlargeSIG='Largest_PHAAdjs_SIG_BNCHCAAA_SEP2015.txt'
outlistallALL='PHAAdjustmentlist_ALL_BNCHCAAA_SEP2015.txt'
#outlistallSIG='PHAAdjustmentlist_SIG_BNCHCAAA_SEP2015.txt'

# hard variables
MDI=-1e+30 # is this needed?

StYear=1800
EdYear=2015
NYrs=(EdYear-StYear)+1
NMonths=NYrs*12

# soft variables
nPHAStations=0 # filled in after read in
AllStations=0	# will be a np.array of all ISTI IDs
StationIDs=np.array((),dtype=str)	# empty np.array to grow of Station IDs where adj > abs(0.) for each station
# These adjustment sizes are the relative sizes (compared to next (more recent) HSP, log reports actual (relative to present day)
StationAdjs=np.array(())	# empty np.array to grow of Adjustment values > abs(0.) for each station
StationUncs=np.array(())	# empty np.array to grow of Uncertainty values where adj > abs(0.) for filtering
StationTotAdjs=np.array(())	# empty np.array to grow of Adjustment values > abs(0.) for each station
StationTotUncs=np.array(())	# empty np.array to grow of Uncertainty values where adj > abs(0.) for filtering
AdjLocTots=np.zeros(NMonths)	# np. array to count number of changepoints occurring each month for all stations
AdjLocTotsSigs=np.zeros(NMonths)	# np. array to count number of significant changepoints occurring each month for all stations

#########################################################################
# Functions
#########################################################################
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
# WRITEOUTLIST
def WriteOutList(TheStations,TheAdjs,TheUncs,TheFile):
    ''' Writes out list of Station IDs, Adjustment size, Adjustment uncs from largest to smallest (relative) '''
    # Get length of file
    nlines=len(TheStations)
    
    # Open file for writing
    My_Fhandle=file(TheFile,'a')
    
    # Loop through lines
    for ll in range(nlines):            
        # write to file - other lines
	goo=np.reshape(np.array(("{:11s}".format(TheStations[ll]),
	                             "{:7.2f}".format(TheAdjs[ll]),
				     "{:7.2f}".format(TheUncs[ll]))),(1,3))
        np.savetxt(My_Fhandle,goo,fmt='%s',delimiter=' ')
    My_Fhandle.close()
    
    return # WriteOutList
#************************************************************************
# PLOTHISTTIMESERIES
def PlotHistTimeSeries(TheFile,StationAdjs,AdjLocTots,StartYr,EndYr,ProjName,nStations):
    ''' Plots histogram of adjustment sizes over laid with mean and standard devations '''
    ''' Lower panels shows counts of changepoints over time '''
    
    # Set up date time thingies
    YrCount=(EndYr-StartYr)+1    
    MonCount=(YrCount*12)
    TheMonths=np.array(())
    yr=StartYr
    mon=1
    for m in range(MonCount):
        TheMonths=np.append(TheMonths,(dt.date(yr,mon,1)))
	mon=mon+1
	if mon == 13:
	    mon=1
	    yr=yr+1   
    
    # Set up figure and plot space
    f,ax=plt.subplots(2,figsize=(6,8))
    
    # Plot Histogram
    ax[0].set_position([0.12,0.55,0.8,0.40])
    ax[0].set_xlim([-5,5])
    ax[0].set_ylim([0,1000])
    ax[0].hist(StationAdjs,bins=40)
    ax[0].set_title('Adjustment size frequency for '+ProjName,fontdict=None,loc='center',size=12)
    ax[0].set_xlabel('Adjustment size ($^{o}$C)',fontsize=12)
    ax[0].set_ylabel('Change Point count',fontsize=12)
    MeanAbs="{:7.3f}".format(np.mean(abs(StationAdjs)))
    MeanAct="{:7.3f}".format(np.mean(StationAdjs))
    SdAbs="{:7.3f}".format(np.std(abs(StationAdjs)))
    SdAct="{:7.3f}".format(np.std(StationAdjs))
    ax[0].annotate('Mean: '+MeanAbs+', '+MeanAct,xy=(0.02,0.9),xycoords='axes fraction',size=12,ha='left')
    ax[0].annotate('Std Dev: '+SdAbs+', '+SdAct,xy=(0.02,0.8),xycoords='axes fraction',size=12,ha='left')
    
    # Plot time series
    ax[1].set_position([0.12,0.05,0.8,0.40])
    ax[1].set_xlim([TheMonths[0],TheMonths[MonCount-1]])
    ax[1].set_ylim([0,25])
    ax[1].plot(TheMonths[0:len(TheMonths)],AdjLocTots,linewidth=1)
    ax[1].set_title('Change point total over time for '+ProjName,fontdict=None,loc='center',size=12)
    ax[1].set_ylabel('Change point frequency',fontsize=12)
    ax[1].set_xlabel('Time',fontsize=12)
    CpFreq="{:7.3f}".format(np.float(len(StationAdjs))/np.float(nStations))
    ax[1].annotate('Mean Change Points per station: '+CpFreq,xy=(0.02,0.9),xycoords='axes fraction',size=12,ha='left')

    # Save
    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")

    return # PlotHistTimeSeries
#************************************************************************
#########################################################################
# Main program
#########################################################################
# Read in list of all stations processed by PHA
MyTypes=("|S11","float","float","int","|S34")
MyDelimiters=[11,7,10,13,34]
MyColumns='XXX'
RawData=ReadData(instats,MyTypes,MyDelimiters,False,MyColumns)
AllStations=np.array(RawData['f0'],ndmin=1)
nPHAStations=len(AllStations)

# Read in list of all logged adjustments
# read in this temporary file
MyTypes=np.append(np.append(np.append("|S10","|S11"),np.repeat("int",11)),np.repeat("float",2))
MyDelimiters=[10,11,3,5,2,5,8,5,2,5,8,5,5,7,7]
MyColumns='XXX'
RawData=ReadData(inlog,MyTypes,MyDelimiters,False,MyColumns)
AdjIDs=np.array(RawData['f1'],ndmin=1)
AdjStYr=np.array(RawData['f3'],ndmin=1)
AdjStMn=np.array(RawData['f4'],ndmin=1)
AdjEdYr=np.array(RawData['f7'],ndmin=1)
AdjEdMn=np.array(RawData['f8'],ndmin=1)
AdjLoc=np.array(RawData['f9'],ndmin=1)
# This is total size/uncertainty relative to present day
AdjSize=np.array(RawData['f13'],ndmin=1)
AdjUnc=np.array(RawData['f14'],ndmin=1)

# Loop through each station and pull out Adjustment info if there is any
for ss in range(nPHAStations):
    
    print(AllStations[ss])
    
    ## Set up string to search for
    #findline='^"Adj write:"'+AllStations[ss]
    #print(findline)
    
    # TOO SLOW!!
    ## spawn a grep command to check if there are any adjustments
    #foo=sp.check_output('grep -c '+findline+' '+inlog+' || true',shell=True).split('\n')
    
    # find candidate station in AdjIDs list
    foo=np.where(AdjIDs == AllStations[ss])[0]
    #pdb.set_trace()
    # There must be more than one line present
    if (len(foo) > 1):
        
	# TOO SLOW!!!
	## spawn linux grep command and output to file
        #sp.call(['grep '+findline+' '+inlog+' > tmp.arr'],shell=True)

	nAdjs=len(foo)
	print("Got adjustments: ",nAdjs)
	
	# Loop through adjustments and output to full file list and fill global arrays
	My_Fhandle=file(outlistdir+outlistallALL,'a')
	# write to file - first line (no adj)
	tmpreladj=0.
	tmprelunc=0.
	goo=np.reshape(np.array(("{:11s}".format(AllStations[ss]),
	                         "{:5d}".format(AdjStYr[foo[0]]),
			 	 "{:2d}".format(AdjStMn[foo[0]]),
				 "{:5d}".format(AdjEdYr[foo[0]]),
				 "{:2d}".format(AdjEdMn[foo[0]]),
				 "{:7.2f}".format(tmpreladj),
				 "{:7.2f}".format(tmprelunc),
				 "{:7.2f}".format(-(AdjSize[foo[0]])),
				 "{:7.2f}".format(AdjUnc[foo[0]]))),(1,9))
        np.savetxt(My_Fhandle,goo,fmt='%s',delimiter=' ')
        for ll in range(1,nAdjs):            
	    # push stats to np.arrays
            StationIDs=np.append(StationIDs,AllStations[ss])	# empty np.array to grow of Station IDs where adj > abs(0.) for each station
            if (ll == 1): # first line is reference period, no adj, second line is first adj
	        tmpreladj=-(AdjSize[foo[ll]])
		tmprelunc=AdjUnc[foo[ll]]
            else:
	        tmpreladj=(-(AdjSize[foo[ll]]))-(-(AdjSize[foo[ll-1]]))	# subtract more recent HSP
                tmprelunc=np.sqrt((AdjUnc[foo[ll-1]]**2) - (AdjUnc[foo[ll]]**2))	# =SQRT(more_recent_HSP^2 + this_HSP^2)
	    StationAdjs=np.append(StationAdjs,tmpreladj)	# empty np.array to grow of Adjustment values > abs(0.) for each station
            StationUncs=np.append(StationUncs,tmprelunc)	# empty np.array to grow of Uncertainty values where adj > abs(0.) for filtering
	    StationTotAdjs=np.append(StationTotAdjs,AdjSize[foo[ll]])	# empty np.array to grow of Adjustment values > abs(0.) for each station
            StationTotUncs=np.append(StationTotUncs,AdjUnc[foo[ll]])	# empty np.array to grow of Uncertainty values where adj > abs(0.) for filtering
	    AdjLocTots[AdjLoc[foo[ll]]-1]=AdjLocTots[AdjLoc[foo[ll]]-1]+1	# np. array to count number of changepoints occurring each month for all stations
	    if (AdjUnc[foo[ll]] < abs(AdjSize[foo[ll]])):
	        AdjLocTotsSigs[AdjLoc[foo[ll]]-1]=AdjLocTotsSigs[AdjLoc[foo[ll]]-1]+1	# np. array to count number of changepoints occurring each month for all stations

	    # write to file - other lines
	    goo=np.reshape(np.array(("{:11s}".format(AllStations[ss]),
	                             "{:5d}".format(AdjStYr[foo[ll]]),
				     "{:2d}".format(AdjStMn[foo[ll]]),
				     "{:5d}".format(AdjEdYr[foo[ll]]),
				     "{:2d}".format(AdjEdMn[foo[ll]]),
				     "{:7.2f}".format(tmpreladj),
				     "{:7.2f}".format(tmprelunc),
				     "{:7.2f}".format(-(AdjSize[foo[ll]])),
				     "{:7.2f}".format(AdjUnc[foo[ll]]))),(1,9))
            np.savetxt(My_Fhandle,goo,fmt='%s',delimiter=' ')
	My_Fhandle.close()

print("Plotting")        
# make plots of all adjustments
PlotHistTimeSeries(outplotdir+outplotstatALL,StationAdjs,AdjLocTots,StYear,EdYear,ProjectName,nPHAStations)

# make plots of only significant adjustments (unc < abs(adj))
PlotHistTimeSeries(outplotdir+outplotstatSIG,StationAdjs[np.where(StationTotUncs < abs(StationTotAdjs))[0]],
				  AdjLocTotsSigs,StYear,EdYear,ProjectName,nPHAStations)

print("Sorting")
# sort list from largest to smallest of absolute values
SortedAdjPointers=np.flipud(np.argsort(abs(StationAdjs)))

# output list of largest to smallest
WriteOutList(StationIDs[SortedAdjPointers],StationAdjs[SortedAdjPointers],StationUncs[SortedAdjPointers],outlistdir+outlistlargeALL)

# output stats
print("Changepoint Frequency: ",np.float(len(StationIDs))/np.float(nPHAStations),', ',
      np.float(len(np.where(StationUncs < abs(StationAdjs))[0]))/np.float(nPHAStations))
print('ABSOLUTE MEAN: ',np.mean(abs(StationAdjs)),', ',np.mean(abs(StationAdjs[np.where(StationTotUncs < abs(StationTotAdjs))[0]])))
print('ABSOLUTE MEDIAN: ',np.median(abs(StationAdjs)),', ',np.median(abs(StationAdjs[np.where(StationTotUncs < abs(StationTotAdjs))[0]])))
print('ABSOLUTE ST DEV',np.std(abs(StationAdjs)),', ',np.std(abs(StationAdjs[np.where(StationTotUncs < abs(StationTotAdjs))[0]])))
print(' ')
print('ACTUAL MEAN: ',np.mean(StationAdjs),', ',np.mean(abs(StationAdjs[np.where(StationTotUncs < StationTotAdjs)[0]])))
print('ACTUAL MEDIAN: ',np.median(StationAdjs),', ',np.median(abs(StationAdjs[np.where(StationTotUncs < StationTotAdjs)[0]])))
print('ACTUAL ST DEV',np.std(StationAdjs),', ',np.std(abs(StationAdjs[np.where(StationTotUncs < StationTotAdjs)[0]])))

pdb.set_trace()

print("And, were done!")
#########################################################################
