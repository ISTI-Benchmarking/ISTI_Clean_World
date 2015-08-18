#!/usr/local/sci/bin/python

#***************************************
# 30 July 2015 KMW - v1
# Plots station locations for ISTI:
# 	a) stations that can be simulated (32522 - JUL2015)
#	b) stations that are too short to be simulated (2870 - JUL2015) and
#	   stations that have identical locations with another station (510 - JUL2015)
#		only one station for each location is kept in the the subset to simulate 
#          stations that are pretending to be ships (30 - JUL2015)
#
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotISTIStationLocs_JUL2015.py
#
# REQUIRES
# 
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import matplotlib as mpl
import numpy as np
import sys, os
import scipy.stats
import struct
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
from netCDF4 import Dataset


# RESTART VALUE
Restarter='------'				#'------'		#'681040'

# Set up initial run choices

# Set up directories and files
INGOOD='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat'
INSHORT='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTISHORTS_stage3proxyelevs_JUL2015.dat'
INMATCH='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGSLOCMATCH_stage3proxyelevs_JUL2015.dat'
INSHIP='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/INVENTORY_monthly_merged_stage3_ships'
OUTDIR='/data/local/hadkw/ISTI/IMAGES/'

OUTPLOT='BAWGISTIStationLocationsMap_JUL2015'

# Set up variables
ngoods=0	#set once file read in
nshorts=0		#set once file read in
nmatches=0		#set once file read in
nships=0

GoodWMOs=[]
GoodLats=[]
GoodLons=[]
GoodStarts=[]
GoodEnds=[]
GoodLengths=[]
ShortWMOs=[]
ShortLats=[]
ShortLons=[]
ShortStarts=[]
ShortEnds=[]
ShortLengths=[]
MatchWMOs=[]
MatchLats=[]
MatchLons=[]
MatchStarts=[]
MatchEnds=[]
MatchLengths=[]
ShipWMOs=[]
ShipLats=[]
ShipLons=[]
ShipStarts=[]
ShipEnds=[]
ShipLengths=[]

Letty=['a)','b)']
Namey=['Simulatable Stations (32522)','Non-simulatable Stations (3410)']

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
# PlotNiceDotsMap
def PlotNiceDotsMap(TheFile,TheGLats,TheGLons,TheGLengths,TheSLats,TheSLons,TheMLats,TheMLons,TheSHLats,TheSHLons,TheLetter,TheNamee):
    ''' Plot all good lats and lons TheGLats, TheGLons  - coloured by length of record'''
    ''' Plot all of the other lats and lons - red, black, blue'''
    ''' Label plot with totals '''
    ''' Save as png and eps '''
      
    # set up plot
 
    plt.clf()
    f=plt.figure(2,figsize=(8,10))
        # set up plot
    xpos=[0.05,0.05]
    ypos=[0.55,0.05]
    xfat=[0.9,0.9]
    ytall=[0.4,0.4]

    plt.axes([xpos[0],ypos[0],xfat[0],ytall[0]])
    
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))

    # set up colour scheme
    cmap=plt.get_cmap('rainbow')	#'gnuplot2' 
    # extract all of the colours    
    cmaplist=[cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey
    cmaplist[0] = (.5,.5,.5,1.0)
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    
    bounds=np.array([0,5,10,20,30,50,100,200])
    strbounds=["%3d" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    # make the scatter
    #scat = ax.scatter(x,y,c=tag,s=np.random.randint(100,500,20),cmap=cmap, norm=norm)
    x,y = m(TheGLons,TheGLats)
    maps=m.scatter(x,y,c=TheGLengths,s=3,marker='o',cmap=cmap,norm=norm,edgecolor='0.0',linewidth=0.0)

    # create a second axes for the colorbar
    cbax=f.add_axes([0.05,0.52,0.9,0.02])
    cb=plt.colorbar(maps,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=10) 
    plt.figtext(0.5,0.49,'Station Record Length (years)',size=12,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.05,0.92,TheLetter[0],size=16)
    plt.figtext(0.5,0.96,TheNamee[0],size=16,ha='center')
    
    #plt.annotate(str(len(TheGLons))+" good stations",xy=(0.52,0.265),xytext=None, xycoords='axes fraction',color="black")


    plt.axes([xpos[1],ypos[1],xfat[1],ytall[1]])
    
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))

    x,y = m(TheSLons,TheSLats)
    m.scatter(x,y,c='red',s=3,marker='o',edgecolor='0.0',linewidth=0.0,)

    x,y = m(TheMLons,TheMLats)
    m.scatter(x,y,c='grey',s=5,marker='o',edgecolor='0.0',linewidth=0.0,)

    x,y = m(TheSHLons,TheSHLats)
    m.scatter(x,y,c='dodgerblue',s=10,marker='o',edgecolor='0.0',linewidth=0.0,)

    #ax.set_title(TheNamee[1])
        
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
    plt.figtext(0.05,0.42,TheLetter[1],size=16)
    plt.figtext(0.5,0.46,TheNamee[1],size=16,ha='center')
    
    plt.figtext(0.15,0.02,str(len(TheSLons))+" short stations",size=16,ha='center',color="red")
    plt.figtext(0.5,0.02,str(len(TheMLons))+" location match stations",size=16,ha='center',color="grey")
    plt.figtext(0.85,0.02,str(len(TheSHLons))+" ship stations",size=16,ha='center',color="dodgerblue")

#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotNiceDotsMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list
MyTypes=("|S12","|S31","|S24","float","float","float","int","int","int","int","int","int","|S16","|S17")
MyDelimiters=[12,31,24,8,11,9,5,5,5,5,5,5,16,17]
MyColumns='XXX'
RawData=ReadData(INGOOD,MyTypes,MyDelimiters,False,MyColumns)
GoodWMOs=np.array(RawData['f0'],ndmin=1)
GoodLats=np.array(RawData['f3'],ndmin=1)
GoodLons=np.array(RawData['f4'],ndmin=1)
GoodStarts=np.array(RawData['f10'],ndmin=1)
GoodEnds=np.array(RawData['f11'],ndmin=1)
GoodLengths=(GoodEnds-GoodStarts)+1
ngoods=len(GoodWMOs)

RawData=ReadData(INSHORT,MyTypes,MyDelimiters,False,MyColumns)
ShortWMOs=np.array(RawData['f0'],ndmin=1)
ShortLats=np.array(RawData['f3'],ndmin=1)
ShortLons=np.array(RawData['f4'],ndmin=1)
ShortStarts=np.array(RawData['f10'],ndmin=1)
ShortEnds=np.array(RawData['f11'],ndmin=1)
ShortLengths=(ShortEnds-ShortStarts)+1
nshorts=len(ShortWMOs)

RawData=ReadData(INMATCH,MyTypes,MyDelimiters,False,MyColumns)
MatchWMOs=np.array(RawData['f0'],ndmin=1)
MatchLats=np.array(RawData['f3'],ndmin=1)
MatchLons=np.array(RawData['f4'],ndmin=1)
MatchStarts=np.array(RawData['f10'],ndmin=1)
MatchEnds=np.array(RawData['f11'],ndmin=1)
MatchLengths=(MatchEnds-MatchStarts)+1
nmatches=len(MatchWMOs)

RawData=ReadData(INSHIP,MyTypes,MyDelimiters,False,MyColumns)
ShipWMOs=np.array(RawData['f0'],ndmin=1)
ShipLats=np.array(RawData['f3'],ndmin=1)
ShipLons=np.array(RawData['f4'],ndmin=1)
ShipStarts=np.array(RawData['f10'],ndmin=1)
ShipEnds=np.array(RawData['f11'],ndmin=1)
ShipLengths=(ShipEnds-ShipStarts)+1
nships=len(ShipWMOs)

# pass to plotter
    
MyFile=OUTDIR+OUTPLOT
PlotNiceDotsMap(MyFile,GoodLats,GoodLons,GoodLengths,ShortLats,ShortLons,MatchLats,MatchLons,ShipLats,ShipLons,Letty,Namey)
		
#    stop()

print("And, we are done!")

