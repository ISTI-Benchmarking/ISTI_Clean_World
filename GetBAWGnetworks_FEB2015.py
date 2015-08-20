#!/usr/local/sci/bin/python

#***************************************
# 13 June 2014 KMW - v1
# 
# Loop through ISTI station inventory
# Find 40th corr list for each station
# Only keep stations that have a non zero fortieth correlation 
# For each network (specified by BAWGLONGNetworks3_FEB2015.txt)
#	find all stations within
#	find all 'neighbour' stations (within an extra half gridbox surrounding candidate gridbox
# output:
#       Stations within each network with neighbours appended afterwards (and their number in full list)
#	Stations within network only (and their number in full list)
# 	Neighbour stations for each network only (and their number in full list)
# 	Complete station list - that has correlating neighbours
# 	Complete station list of stations removed because they do not have any correlating neighbours
# 	Network boundaries (easy to work out on the fly but here anyway
#************************************************************************
# NOTES:
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 GetBAWGnetworks_FEB2015.py
#
# REQUIRES
#************************************************************************
# Set up python imports
import datetime as dt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from mpl_toolkits.basemap import Basemap
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

# SWITCH TO READ IN NETWORK BOXES OR CREATE REGULAR ONES
BoxSwitch=1	#1=Read from file, 0=Create regular ones

# Set up file locations
STATLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat'
DISTLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat'
NETWORKLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/BAWGLONGNetworks3_JUL2015.txt'

OUTSTonly='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGNetwork_StOnly_stage3proxyelevs_JUL2015_'
OUTNBonly='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGNetwork_NbOnly_stage3proxyelevs_JUL2015_'
OUTSTNB='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGNetwork_StNb_stage3proxyelevs_JUL2015_'
OUTNETWORKS='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGNetwork_Boundaries_stage3proxyelevs_JUL2015.dat'
OUTPLOTTOTS='/data/local/hadkw/ISTI/IMAGES/ISTILONGNetwork_Boundaries_stage3proxyelevs_JUL2015'
OUTPLOTSTATS='/data/local/hadkw/ISTI/IMAGES/ISTILONGNetwork_Stations_stage3proxyelevs_JUL2015'

nstations=0	# defined after reading in station list
StationIDs=[]	# nstations list filled after reading in station list
StationLats=[]	# nstations list filled after reading in station list
StationLons=[]	# nstations list filled after reading in station list
StationElevs=[]	# nstations list filled after reading in station list
StationNum=[]	# integer array numbering stations from 1 to n
NetworkLocs=[]

StatNbsAll=[]
StatsOnly=[]
NbsOnly=[]

# Set up any key variables
LonWd=18.
LatWd=9.
StLon=-180.
StLat=-90.
nBoxes=int((360/LonWd)*(180/LatWd))	# should be 400

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

#***********************************************************************
# MAIN PROGRAM
#***********************************************************************
# read in station list
MyTypes=("|S12","|S31","|S24","float","float","float","int","int","int","int","int","int","|S16","|S17")
MyDelimiters=[12,31,24,8,11,9,5,5,5,5,5,5,16,17]
MyColumns='XXX'
RawData=ReadData(STATLIST,MyTypes,MyDelimiters,False,MyColumns)
StationIDs=np.array(RawData['f0'],ndmin=1)
COPYStationIDs=np.array(RawData['f0'],ndmin=1)
StationLats=np.array(RawData['f3'],ndmin=1)
StationLons=np.array(RawData['f4'],ndmin=1)
StationElevs=np.array(RawData['f5'],ndmin=1)
StationElevsCopy=np.array(RawData['f5'],ndmin=1)
nstations=len(StationIDs)
print('Read in Giant file...')

# read in 100th dists and station list (just to be sure)
a=range(100)
MyTypes=np.append("|S12",np.reshape(np.array([["|S13","float"] for i in a]),200))
MyDelimiters=np.append(12,[13]*200)
MyColumns=(0,1,2)
RawDists=ReadData(DISTLIST,MyTypes,MyDelimiters,True,MyColumns)
StationDIDs=np.array(RawDists['f0'],ndmin=1)
StationCLOSEIDs=np.array(RawDists['f1'],ndmin=1)
Dist100s=np.array(RawDists['f2'],ndmin=1) # actually the first correlating station
nDstations=len(StationDIDs)
print(nstations,nDstations)
RawDists=[]
#stop()

# Set up station pointer list
StationNums=np.array(range(nstations))+1
	
# If BoxSwitch is set to 1 then read in Network boxes from file
if BoxSwitch == 1:
    MyTypes=("int","float","float","float","float")
    MyDelimiters=[2,9,9,9,9]
    MyColumns=(1,2,3,4)
    RawData=ReadData(NETWORKLIST,MyTypes,MyDelimiters,False,MyColumns)
    #BoxLon1=np.array(RawData['f1'],ndmin=1)
    #BoxLat1=np.array(RawData['f2'],ndmin=1)
    #BoxLon2=np.array(RawData['f3'],ndmin=1)
    #BoxLat2=np.array(RawData['f4'],ndmin=1)
    BoxLon1=RawData['f1']
    BoxLat1=RawData['f2']
    BoxLon2=RawData['f3']
    BoxLat2=RawData['f4']
    #print(list(BoxLon1))
    #stop()
    NetworkLocs=np.transpose(np.array([list(BoxLon1),list(BoxLat1),list(BoxLon2),list(BoxLat2)]))
    print(NetworkLocs)
    nBoxes=len(BoxLon1)
    print('Read in Network Boxes from file...',nBoxes)
 
# Set up map to plot boundaries and station totals
    #pltb.clf()
    plt.figure(1,figsize=(16,10))
    plt.axes([0.05,0.05,0.9,0.95])
    
    # plot map without continents and coastlines
    mb = Basemap(projection='cyl',lon_0=0)
    # draw map boundary, transparent
    mb.drawmapboundary()
    mb.drawcoastlines()
    # draw paralells and medians, no labels
    mb.drawparallels(np.arange(-90,90.,30.))
    mb.drawmeridians(np.arange(-180,180.,60.))

    #plts.clf()
    plt.figure(2,figsize=(16,10))
    plt.axes([0.05,0.05,0.9,0.95])
    
    # plot map without continents and coastlines
    ms = Basemap(projection='cyl',lon_0=0)
    # draw map boundary, transparent
    ms.drawmapboundary()
    ms.drawcoastlines()
    # draw paralells and medians, no labels
    ms.drawparallels(np.arange(-90,90.,30.))
    ms.drawmeridians(np.arange(-180,180.,60.))
        	
    colourtemp=('DarkRed','Red','Tomato','DarkOrange','Yellow','DarkGreen','LimeGreen',
               'DarkBlue','DodgerBlue','Indigo','DarkViolet','HotPink','Pink','LightSalmon',
               'PeachPuff','LightGreen','DarkSeaGreen','LightSkyBlue','Aqua',
               'BlueViolet','Lavender')

    	
# Set up and loop through Network boxes and boundaries
My_Fhandle=file(OUTNETWORKS,'a')	# NetworkNumb,lonW,latS,lonE,latN,TotStations,TotNeighbours
mvlat=0	# increment to 20
mvlon=0	# increment to 20
colcount=0	# maxes out at 20 so reset when 21
for bb in range(nBoxes):
    StationFinds=[]
    TotalStations=0
    TotalNeighbours=0
    NbsStationFinds=[]
    
    if BoxSwitch == 0:	# in this case LonWd and LatWd are set above
        lat1=StLat+(mvlat*LatWd)
        lat2=StLat+(mvlat*LatWd)+LatWd
        lon1=StLon+(mvlon*LonWd)
        lon2=StLon+(mvlon*LonWd)+LonWd
    
        mvlon=mvlon+1
        if mvlon == 20:
          mvlon=0
	  mvlat=mvlat+1
    
    else:
        LonWd=10
	LatWd=10
        # Set up the neighbour catchment boundary to be half the lat/lon length of the network box
	# unless that is larger than 10 deg.
        lon1,lat1,lon2,lat2=NetworkLocs[bb,:]
        if lon2 > lon1:
            testee=(lon2-lon1)
        else:
            testee=((180.-lon1)+(lon2-(-180.)))
	if testee < 10:
	    LonWd=testee
	
	if (lat2-lat1) < 10:
	    LatWd=(lat2-lat1)
       
	
    locs=np.array((lon1,lat1,lon2,lat2))
    print(bb,locs,LonWd,LatWd)
        
    # find all stations within this network
    if lat2 < 90.:
        stgotlats=np.where(np.logical_and(StationLats>=lat1,StationLats<lat2))[0]
    else:
        stgotlats=np.where(np.logical_and(StationLats>=lat1,StationLats<=lat2))[0]
    
    if lon2 > lon1:
        stgotlons=np.where(np.logical_and(StationLons>=lon1,StationLons<lon2))[0]
    else:
        stgotlons=np.where(StationLons>=lon1)[0]
	stgotlons=np.append(stgotlons,np.where(StationLons<lon2)[0])
    
    if len(stgotlats) > 0 and len(stgotlons) > 0:
        test=list(set(list(stgotlats)).intersection(list(stgotlons)))
        if len(test) > 0: 
            StationFinds=np.array(test)
            TotalStations=len(StationFinds)
            print('Station Total: ',TotalStations)
        
	    stgotlons=[]
            stgotlats=[]
    
    # find all neighbours for this network
    # look at western set (includes SW and NW
    # SPECIAL WHEN lon1=-180., AND lat1=-90, AND lat2=90
            nlon1=lon1-(LonWd/2.)
	    nlon2=lon1
            if nlon1 < -180:
	        nlon1=180-(-180 - nlon1)
            if lat1 != -90.:
                nlat1=lat1-(LatWd/2)
            else:
                nlat1=lat1
            if lat2 != 90.:	    
                nlat2=lat2+(LatWd/2)
            else:
                nlat2=lat2
            # special case at poles where boxes circumnavigate therefore no W/E or poleward neighbours
	    if lon1 != -180 or lon2 != 180:            
                if nlat2 < 90.:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<nlat2))[0]
                else:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<=nlat2))[0]
                if nlon2 > nlon1:
                    gotlons=np.where(np.logical_and(StationLons>=nlon1,StationLons<nlon2))[0]
                else:
                    gotlons=np.where(StationLons>=nlon1)[0]
	            gotlons=np.append(gotlons,np.where(StationLons<nlon2)[0])
                if len(gotlats) > 0 and len(gotlons) > 0:
	            test=list(set(list(gotlats)).intersection(list(gotlons)))
                    if len(test) > 0:
	                NbsStationFinds=np.array(test)
                    else:
	                NbStationFinds=[]    
	        else:
	            NbStationFinds=[]
	    else:
	        NbStationFinds=[]
	    print('W NBs & Total: ',nlat1,nlon1,nlat2,nlon2,len(NbsStationFinds))
		
    # look at southern set (no SW or SE)
    # SPECIAL WHEN lat1=-90
            if lat1 != -90.:
                nlat1=lat1-(LatWd/2)
                nlat2=lat1
                nlon1=lon1
                nlon2=lon2
    
                if nlat2 < 90.:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<nlat2))[0]
                else:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<=nlat2))[0]
                if nlon2 > nlon1:
                    gotlons=np.where(np.logical_and(StationLons>=nlon1,StationLons<nlon2))[0]
                else:
                    gotlons=np.where(StationLons>=nlon1)[0]
	            gotlons=np.append(gotlons,np.where(StationLons<nlon2)[0])
                if len(gotlats) > 0 and len(gotlons) > 0:
	            test=list(set(list(gotlats)).intersection(list(gotlons)))
                    if len(test) > 0:
	                NbsStationFinds=np.append(NbsStationFinds,np.array(test))
               
	        print('S NBs & Total: ',nlat1,nlon1,nlat2,nlon2,len(NbsStationFinds))
    
    # look at eastern set (includes SE adn NE)
    # SPECIAL WHEN lon2=180., AND lat1=-90, AND lat2=90
            nlon2=lon2+(LonWd/2)
	    nlon1=lon2
            if nlon2>180.:
	        nlon2=(-180)+(nlon2-180)
            if lat1 != -90.:
                nlat1=lat1-(LatWd/2)
            else:
                nlat1=lat1
            if lat2 != 90.:	    
                nlat2=lat2+(LatWd/2)
            else:
                nlat2=lat2

            # special case at poles where boxes circumnavigate therefore no W/E or poleward neighbours
	    if lon1 != -180 or lon2 != 180:
                if nlat2 < 90.:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<nlat2))[0]
                else:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<=nlat2))[0]
                if nlon2 > nlon1:
                    gotlons=np.where(np.logical_and(StationLons>=nlon1,StationLons<nlon2))[0]
                else:
                    gotlons=np.where(StationLons>=nlon1)[0]
	            gotlons=np.append(gotlons,np.where(StationLons<nlon2)[0])
                if len(gotlats) > 0 and len(gotlons) > 0:
	            test=list(set(list(gotlats)).intersection(list(gotlons)))
                    if len(test) > 0.:
	                NbsStationFinds=np.append(NbsStationFinds,np.array(test))
            
	    print('E NBs & Total: ',nlat1,nlon1,nlat2,nlon2,len(NbsStationFinds))
    
    # look at northern set (no NW or NE)
    # SPECIAL WHEN lat1=-90, AND lat2=90
            if lat2 != 90.:
                nlat2=lat2+(LatWd/2)
                nlat1=lat2
                nlon1=lon1
                nlon2=lon2
    
                if nlat2 < 90.:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<nlat2))[0]
                else:
                    gotlats=np.where(np.logical_and(StationLats>=nlat1,StationLats<=nlat2))[0]
                if nlon2 > nlon1:
                    gotlons=np.where(np.logical_and(StationLons>=nlon1,StationLons<nlon2))[0]
                else:
                    gotlons=np.where(StationLons>=nlon1)[0]
	            gotlons=np.append(gotlons,np.where(StationLons<nlon2)[0])
                if len(gotlats) > 0 and len(gotlons) > 0:
	            test=list(set(list(gotlats)).intersection(list(gotlons)))
                    if len(test) > 0.:
	                NbsStationFinds=np.append(NbsStationFinds,np.array(test))
            
	    print('N NBs & Total: ',nlat1,nlon1,nlat2,nlon2,len(NbsStationFinds))

            TotalNeighbours=len(NbsStationFinds)
            gotlons=[]
            gotlats=[]
    
    
    # output lists of stations in each network to file
            My_FOUThandle=file(OUTSTonly+str(bb+1)+'.dat','a')
            My_FOUTALLhandle=file(OUTSTNB+str(bb+1)+'.dat','a')
            for ss in range(TotalStations):
                goo=np.array("{:5d}".format(StationNums[StationFinds[ss]]))
	        goo=np.append(goo,np.array("{:8s}".format(StationIDs[StationFinds[ss]])))
	        goo=np.append(goo,np.array("{:9.3f}".format(StationLats[StationFinds[ss]])))
	        goo=np.append(goo,np.array("{:9.3f}".format(StationLons[StationFinds[ss]])))
                np.savetxt(My_FOUThandle,np.reshape(goo,(1,4)),fmt='%s',delimiter=' ')
                np.savetxt(My_FOUTALLhandle,np.reshape(goo,(1,4)),fmt='%s',delimiter=' ')
	
            My_FOUThandle.close()

            My_FOUThandle=file(OUTNBonly+str(bb+1)+'.dat','a')
            for ss in range(TotalNeighbours):
                goo=np.array("{:5d}".format(StationNums[NbsStationFinds[ss]]))
	        goo=np.append(goo,np.array("{:8s}".format(StationIDs[NbsStationFinds[ss]])))
	        goo=np.append(goo,np.array("{:9.3f}".format(StationLats[NbsStationFinds[ss]])))
	        goo=np.append(goo,np.array("{:9.3f}".format(StationLons[NbsStationFinds[ss]])))
                np.savetxt(My_FOUThandle,np.reshape(goo,(1,4)),fmt='%s',delimiter=' ')
                np.savetxt(My_FOUTALLhandle,np.reshape(goo,(1,4)),fmt='%s',delimiter=' ')
	
            My_FOUThandle.close()
            My_FOUTALLhandle.close()
    
        else:
            TotalStations=0
	    TotalNeighbours=0
    
    totals=np.array((TotalStations,TotalNeighbours))
    counter="{:3d}".format(bb+1)
    strlocs=["{:9.3f}".format(ll) for ll in locs]
    strtotals=["{:5d}".format(tt) for tt in totals]
    goo=np.reshape(np.append(np.append(counter,strlocs),strtotals),(1,7))
    np.savetxt(My_Fhandle,goo,fmt='%s',delimiter=' ')

# plot the boundary and station total on the map
    if lon2 > lon1:
        plotlons = np.array((lon1,lon2,lon2,lon1))
        plotlats = np.array((lat1,lat1,lat2,lat2))
        plt.figure(1)
	x,y = mb(plotlons,plotlats)
        poly = Polygon(zip(x,y),facecolor=colourtemp[colcount],edgecolor='none')
#        poly = Polygon(zip(x,y),edgecolor=colourtemp[colcount],facecolor='none')
        plt.gca().add_patch(poly)
        plt.figure(2)
        x,y = ms(plotlons,plotlats)
        poly = Polygon(zip(x,y),edgecolor=colourtemp[colcount],facecolor='none',linewidth=2.)
        plt.gca().add_patch(poly)
    else:
        plotlons = np.array((lon1,180,180,lon1))
        plotlats = np.array((lat1,lat1,lat2,lat2))
        plt.figure(1)
        x,y = mb(plotlons,plotlats)
        poly = Polygon(zip(x,y),facecolor=colourtemp[colcount],edgecolor='none')
#        poly = Polygon(zip(x,y),edgecolor=colourtemp[colcount],facecolor='none')
        plt.gca().add_patch(poly)    
        plt.figure(2)
        x,y = ms(plotlons,plotlats)
        poly = Polygon(zip(x,y),edgecolor=colourtemp[colcount],facecolor='none',linewidth=2.)
        plt.gca().add_patch(poly)
        plotlons = np.array((-180,lon2,lon2,-180))
        plotlats = np.array((lat1,lat1,lat2,lat2))
        plt.figure(1)
        x,y = mb(plotlons,plotlats)
        poly = Polygon(zip(x,y),facecolor=colourtemp[colcount],edgecolor='none')
#        poly = Polygon(zip(x,y),edgecolor=colourtemp[colcount],facecolor='none')
        plt.gca().add_patch(poly)    
        plt.figure(2)
        x,y = ms(plotlons,plotlats)
        poly = Polygon(zip(x,y),edgecolor=colourtemp[colcount],facecolor='none',linewidth=2.)
        plt.gca().add_patch(poly)

    if TotalStations > 0:
        for ss in range(TotalStations):
	    x,y = ms(StationLons[StationFinds[ss]],StationLats[StationFinds[ss]])
            ms.scatter(x,y,s=5,marker='o',color=colourtemp[colcount])
           
	    
    plt.figure(1)
    plt.annotate(TotalStations,xy=(lon1+0.5,lat1+0.5),xytext=None, xycoords='data',size=4,color="black")
    #plt.figure(2)
    #fig2ax.annotate(TotalStations,xy=(lon1+0.5,lat1+0.5),xytext=None, xycoords='data',size=3,color="black")
    
    colcount=colcount+1
    if colcount == 21:
        colcount=0
    #stop()

My_Fhandle.close()

    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
plt.figure(1)
plt.figtext(0.5,0.9,'Factorisation Networks',size=16,ha='center')
plt.figtext(0.5,0.05,'Total No. Boxes: '+str(nBoxes),size=16,ha='center')
plt.savefig(OUTPLOTTOTS+".eps")
plt.savefig(OUTPLOTTOTS+".png")

plt.figure(2)
plt.figtext(0.5,0.9,'Factorisation Networks',size=16,ha='center')
plt.figtext(0.5,0.05,'Total No. Boxes: '+str(nBoxes),size=16,ha='center')
plt.savefig(OUTPLOTSTATS+".eps")
plt.savefig(OUTPLOTSTATS+".png")

print("And, we are done!")
