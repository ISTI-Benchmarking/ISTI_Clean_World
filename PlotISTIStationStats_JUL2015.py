#!/usr/local/sci/bin/python

#***************************************
# 30 July 2015 KMW - v1
# Plots station stats in station locations for ISTI:
# 	a) Station climatological average (annual) - based on entire period
#	b) Station standard deviation - based on entire period
#	c) Station linear trends - based on entire period (where there are at least 15 years of data - CHECK)
#	d) Station autocorrelation from Standardised Anomalies
#
# This can be plotted for the actual stations or the simulated clean world stations
# Should probably have some add on options to plot differences and look at these in more detail
# 	- scatter plots?
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
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)


# RESTART VALUE
Restarter='------'				#'------'		#'681040'

# Set up initial run choices

# Set up directories and files
INGOOD='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat'
INGOODCHECK='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_counts_stage3proxyelevs_JUL2015.dat'
#INSTATS='/data/local/hadkw/ISTI/LISTS/BAWG/SEP2015/OLDSTATS_ISTI_stage3proxyelevs_BNCHCAAA_SEP2015.txt'
INSTATS='/data/local/hadkw/ISTI/LISTS/BAWG/SEP2015/NEWSTATS_ISTI_stage3proxyelevs_BNCHCAAA_SEP2015.txt'
OUTDIR='/data/local/hadkw/ISTI/IMAGES/SEP2015/'

#OUTPLOT='BAWGISTIStationStatsMap_REALDATA_SEP2015'
OUTPLOT='BAWGISTIStationStatsMap_SIMULATEDDATA_SEP2015'

# Set up variables
ngoods=0	#set once file read in
nverygoods=0		#set once file read in

GoodWMOs=[]
GoodLats=[]
GoodLons=[]
GoodLengths=[]
GoodClims=[]
GoodStDevs=[]
GoodTrends=[]
GoodStAnACs=[]

Letty=['a)','b)','c)','d)']
Namey=['Climatological Mean','Standard Deviation','Long-term Decadal Trend','Autocorrelation in Standardised Anomalies']
ColBarTitle=['Temperature ($^{o}$C)','Temperature ($^{o}$C)','Temperature ($^{o}$C)','Autocorrelation at lag 1']

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
def PlotNiceDotsMap(TheFile,TheGLats,TheGLons,TheGLengths,TheGClims,TheGStDevs,TheGTrends,TheGACs,TheLetter,TheNamee,TheCBarTitle):
    ''' Plot all clims  - coloured by deg C, open circles if length < 30 years? TOO SMALL'''
    ''' Plot all st devs  - coloured by deg C, open circles if length < 30 years? TOO SMALL'''
    ''' Plot all trends  - coloured by deg C - no stations fewer than 30 years'''
    ''' Plot all St Anom ACs  - coloured by 0 to 1, open circles if length < 30 years? TOO SMALL'''
    ''' Save as png and eps '''
      
    # set up plot
 
    plt.clf()
    f=plt.figure(4,figsize=(12,8))
        # set up plot
    xpos=[0.025,0.525,0.025,0.525]
    ypos=[0.57,0.57,0.07,0.07]
    xfat=[0.45,0.45,0.45,0.45]
    ytall=[0.4,0.4,0.4,0.4]

    # set up colour scheme
    cmap=plt.get_cmap('coolwarm')	#'gnuplot2' 
    # extract all of the colours    
    cmaplist=[cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    
    MinLength=360
    
    print("Only plotting ",MinLength,"+ month stations: ",len(TheGLengths[np.where(TheGLengths > MinLength)[0]]))

#**************************
    for loo in range(4):
        plt.axes([xpos[loo],ypos[loo],xfat[loo],ytall[loo]])
    
        # plot map without continents and coastlines
        m = Basemap(projection='kav7',lon_0=0)
        # draw map boundary, transparent
        m.drawmapboundary()
        m.drawcoastlines()
        # draw paralells and medians, no labels
        m.drawparallels(np.arange(-90,90.,30.))
        m.drawmeridians(np.arange(-180,180.,60.))
    
        if (loo == 0):
	    bounds=np.array([-50,-40,-30,-20,-15,-10,-5,0,5,10,15,20,25,30])
	    ColourPoints=TheGClims
        elif (loo == 1):
	    bounds=np.array([0.,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0])
	    ColourPoints=TheGStDevs
        elif (loo == 2):
	    bounds=np.array([-2,-1.5,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1.0,1.5,2.])
	    ColourPoints=TheGTrends
        elif (loo == 3):
	    bounds=np.array([-0.8,-0.6,-0.4,-0.2,-0.0,0.2,0.4,0.6,0.8])
	    ColourPoints=TheGACs
	
	strbounds=["%3d" % i for i in bounds]
        norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

	# make the scatter - fill for long stations, just outline for short (<15 years) stations
        x,y = m(TheGLons[np.where(TheGLengths >= MinLength)[0]],TheGLats[np.where(TheGLengths >= MinLength)[0]])
        maps=m.scatter(x,y,c=ColourPoints[np.where(TheGLengths >= MinLength)[0]],s=3,marker='o',cmap=cmap,norm=norm,edgecolor='0.0',linewidth=0.0)
#        # TOO SMALL - NO POINT IN THESE
#	x,y = m(TheGLons[np.where(TheGLengths < MinLength)[0]],TheGLats[np.where(TheGLengths < MinLength)[0]])
#        maps=m.scatter(x,y,c=ColourPoints[np.where(TheGLengths < MinLength)[0]],s=3,marker='o',cmap=cmap,norm=norm,linewidth=1,facecolors='none',edgecolor='0.0')
      
        # create a second axes for the colorbar
        cbax=f.add_axes([xpos[loo],ypos[loo]-0.02,xfat[loo],0.02])
        cb=plt.colorbar(maps,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
        cb.ax.tick_params(labelsize=10) 
        plt.figtext(xpos[loo]+0.225,ypos[loo]-0.065,TheCBarTitle[loo],size=12,ha='center')

    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)
        plt.figtext(xpos[loo],ypos[loo]+0.35,TheLetter[loo],size=16)
        plt.figtext(xpos[loo]+0.225,ypos[loo]+0.405,TheNamee[loo],size=14,ha='center')

#**********************************************

#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotNiceDotsMap
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list
MyTypes=["|S12","|S31","|S24","float","float","float","int","int","int","int","int","int","|S16","|S17"]
MyDelimiters=[12,31,24,8,11,9,5,5,5,5,5,5,16,17]
MyColumns='XXX'
RawData=ReadData(INGOOD,MyTypes,MyDelimiters,False,MyColumns)
GoodWMOs=np.array(RawData['f0'],ndmin=1)
GoodLats=np.array(RawData['f3'],ndmin=1)
GoodLons=np.array(RawData['f4'],ndmin=1)
ngoods=len(GoodWMOs)

#MyTypes=["|S11","int"]
#MyDelimiters=[11,6]
#MyColumns='XXX'
#RawData=ReadData(INGOODCHECK,MyTypes,MyDelimiters,False,MyColumns)
#GoodLengths=np.array(RawData['f1'],ndmin=1,dtype=int)
GoodLengths=np.empty(ngoods,dtype=int)
GoodLengths.fill(1000)

MyTypes=np.append("|S11",np.repeat("float",29))
MyDelimiters=np.append(11,np.repeat(7,29))
MyColumns=[1,2,15,29]
RawData=ReadData(INSTATS,MyTypes,MyDelimiters,False,MyColumns)
GoodClims=np.array(RawData['f2'],ndmin=1)
GoodStDevs=np.array(RawData['f15'],ndmin=1)
GoodTrends=np.array(RawData['f1'],ndmin=1)
GoodStAnACs=np.array(RawData['f29'],ndmin=1)

## TEMPORARY CROP TO 1819 WHILE TESTING
#GoodLats=GoodLats[0:25000]
#GoodLons=GoodLons[0:25000]
#GoodLengths=GoodLengths[0:25000]
#GoodClims=GoodClims[0:25000]
#GoodStDevs=GoodStDevs[0:25000]
#GoodTrends=GoodTrends[0:25000]
#GoodStAnACs=GoodStAnACs[0:25000]

# pass to plotter
    
MyFile=OUTDIR+OUTPLOT
PlotNiceDotsMap(MyFile,GoodLats,GoodLons,GoodLengths,GoodClims,GoodStDevs,GoodTrends,GoodStAnACs,Letty,Namey,ColBarTitle)
		
#    stop()

print("And, we are done!")

