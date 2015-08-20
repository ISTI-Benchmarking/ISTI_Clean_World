#!/usr/local/sci/bin/python

#***************************************
# 17 November 2014 KMW - v1
# 
# Loop through all stations 
# If elevation then print out - new file with proxy elevations where necessary
# If no elevationt then find their approximate elevation from a 0.25 degree digital elevation map 
# Print out
# 
# 6 Feb 2015
# Noticed that INVENTORY includes ship data 
# Firstly - it shouldn't so these are now removed and a list output
# Secondly - these and other low elev stations can result in screwy proxyelevs
# because of the coarse resolution of the Digital Elevation Model used - most
# likely giving an ocean basin depth rather than the true height on land.
# So - we now have a minimum elevation threshold of -100m. This recognises that
# some land is below sea level but not much, and not many stations - minium dry
# land elevation is 418m (shore of Dead Sea, Jordan) according to Wikipedia
#
# 29 Jul 2015 
# Noticed that some of the ships are only pretending to be ships!
# We've identified these because their ID begins with XX but this is obviously not always true
# For example: Appleby in New Zealand, Chatham Islands, Goree off the West African coastline, The North Pole.
# Not sure how/whether to seperate these out with an extra loop BECAUSE
# ERROR FOUND:
# in the case of Appleby - the longitude is given as -173.2 when it should be 173.2!!!
# Left all of these 'ships' out for now

#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 GetStationElevation_FEB2015.py
#
# REQUIRES
#************************************************************************
# Set up python imports
import numpy as np
import sys, os
from scipy.optimize import curve_fit,fsolve,leastsq
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.stats
import struct
#from netCDF4 import Dataset
from scipy.io import netcdf
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

# RESTART VALUE
Restarter='------'	#'------'		#'    01010000' station ID

# Set up file locations
STATLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/INVENTORY_monthly_merged_stage3'
ELEV_DEL='/data/local/hadkw/ISTI/DATA/elev.0.25-deg.nc' #jisao.washington.edu/datasets/elevation/ -> 0.25-degree latitude-longitude resolution elevation
OUTLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/INVENTORY_monthly_merged_stage3_proxyelevs'
OUTNOELEVS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/INVENTORY_monthly_merged_stage3_noelevs'
OUTSHIPS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/INVENTORY_monthly_merged_stage3_ships' # ID beginning with 'XX'

# Set up variables and holders
MDI=-9999
nstations=0	# defined after reading in station list
StationIDs=[]	# nstations list filled after reading in station list
StationLats=[]	# nstations list filled after reading in station list
StationLons=[]	# nstations list filled after reading in station list
StationElevs=[]	# nstations list filled after reading in station list

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
# FINDPROXYELEVATION
def FindProxyElevation(TheLatitude,TheLongitude,TheDELFile,DELLats,DELLons):
    ''' TheLatitude: Station latitude '''
    ''' TheLongitude: Station longitude '''
    ''' TheDELFile: a 2d array of gridbox elevations '''
    ''' DELLats: vector of gridbox lat northern boundaries from north to south 90 to -90 '''
    ''' DELLongs: vector of gridbox lon western boundaries from west to east -180 to 180 '''  
    ''' Find the nearest gridbox in TheDELFile to match TheLatitude and TheLongitude for that station '''
    ''' Return DEL elevation '''
    
    TheElevation=0.0
    
    # This should work because if station is equal to a lat/lon (e.g, 90 N and 180 W) 
    # then that is the box it will be in
    # If station is 90 N and/or 180E then we have a problem!
    # In this special case then use the last Lat/Lon box
    if TheLongitude == 180.0:
        GotLons=len(DELLons)
    else:
        GotLons=np.array(np.where(DELLons > TheLongitude))[0][0] # first element of this
    
    if TheLatitude == -90.0:
        GotLats=len(DELLats)
    else:
        GotLats=np.array(np.where(DELLats < TheLatitude))[0][0] # first element of this
    
    TheElevation=TheDELFile[GotLats-1,GotLons-1]
    
    # check to ensure its not screwy (loads below sea level) - probably due to coarse res/proximity of deep ocean?
    if TheElevation < -100:
        TheElevation=-100
    
    return TheElevation # FindProxyElevation
#***********************************************************************
# MAIN PROGRAM
#***********************************************************************
# read in station list
MyTypes=("|S12","|S31","|S24","float","float","float","int","int","int","int","int","int","|S22","|S11")
MyDelimiters=[12,31,24,8,11,9,5,5,5,5,5,5,22,11]
MyColumns='XXX'
RawData=ReadData(STATLIST,MyTypes,MyDelimiters,False,MyColumns)
StationIDs=np.array(RawData['f0'],ndmin=1)
StationLats=np.array(RawData['f3'],ndmin=1)
StationLons=np.array(RawData['f4'],ndmin=1)
StationElevs=np.array(RawData['f5'],ndmin=1)
nstations=len(StationIDs)
print('Read in Giant file...')

# read in DEL file and get lats and lons 
f=netcdf.netcdf_file(ELEV_DEL,'r')
var=f.variables['data']
DELFilee=np.array(var.data[0,:,:]) # gets rid of 1 element axis
var=f.variables['lat']
DELLats=np.array(var.data)
var=f.variables['lon']
DELLons=np.array(var.data)
f.close()

# make sure its 90 to -90 
if DELLats[0] < 0:
    DELLats.reverse()
    DELFilee=DELFilee[::-1,:]	# reverses rows

# make sure its -180 to 180 
if DELLons[0] >= 0:
    DELLons=DELLons-180
    DELFilee=np.roll(DELFilee,len(DELLons)/2,1)	# shifts columns
    
# make sure these are now north and west gridbox boundaries
LonSize=DELLons[1]-DELLons[0]
LatSize=DELLats[0]-DELLats[1]

if DELLats[0] != 90.0:
    DELLats=DELLats+(LatSize/2.)

if DELLons[0] > -180.0:
    DELLons=DELLons-(LonSize/2.)

# open the new file lists
My_ALLhandle=file(OUTLIST,'a')
My_NOELhandle=file(OUTNOELEVS,'a')
My_SHIPhandle=file(OUTSHIPS,'a')
for lop in range(nstations):

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

    # Is it one of those pesky ships pretending to be land? (may miss the North Pole now but tough)
    Tester=StationIDs[lop].strip()
    if Tester[0:2] == 'XX':
        outraw=np.reshape(outraw,(1,len(outraw)))     
        np.savetxt(My_SHIPhandle,outraw,fmt='%s',delimiter='')	    
        
    
    else:# Look to see if the elevation is missing (-9999.00)
        if StationElevs[lop] <= -9999.:
            # Now look at DEL to find elevation
	    ProxyElevation=FindProxyElevation(StationLats[lop],StationLons[lop],DELFilee,DELLats,DELLons)
	
	    # reformat outraw[5]
            outraw[5]="{:9.2f}".format(ProxyElevation)
 
            # Add this to the noelevs list
            outraw=np.reshape(outraw,(1,len(outraw)))     
            np.savetxt(My_NOELhandle,outraw,fmt='%s',delimiter='')	    
        else:
            outraw=np.reshape(outraw,(1,len(outraw)))     
	
        # print station listing to file (with new elevation in some cases)		
        np.savetxt(My_ALLhandle,outraw,fmt='%s',delimiter='')	    

My_ALLhandle.close()
My_NOELhandle.close()
My_SHIPhandle.close()
stop()

print("And, we are done!")
