#!/usr/local/sci/bin/python

#***************************************
# 28th April 2015
# This version reads in from /data/local/hadkw/HADCRUH2

# Reads in station file for CRUTEM4 (downloaded APR 2015) and for HadISDH.landT
# Plots a map for each with gridboxes coloured by number of stations contributing
# Plots a three panel HadiSDH, CRUTEM-HadISDH, GHCNM-HadISDH difference maps
# Outputs a list of these station groups for each dataset
 
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 PlotStationDensityGrids_JUN2015.py
#
# REQUIRES
# 
#************************************************************************
# Set up python imports
import matplotlib.pyplot as plt
import matplotlib.colors as mc
import matplotlib.cm as mpl_cm
import numpy.ma as ma
import numpy as np
import sys, os
import scipy.stats
import struct
import os.path
import math
from mpl_toolkits.basemap import Basemap
import datetime as dt
from matplotlib.dates import date2num,num2date
#from netCDF4 import Dataset
from scipy.io import netcdf
from scipy.stats.stats import pearsonr   
from RandomsRanges import LetterRange
from LinearTrends import MedianPairwise
import re	# regular expression stuff for character replacement within strings
from scipy.stats import itemfreq # a way of looking for unique values and their frequency of occurence
import pdb	# for stopping and restarting with editability (stop is pdb.set_trace(),restart is c)

# Input Files
INFILEH='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/PosthomogIDPHAt_goodsHadISDH.2.0.1.2014p_JAN2015.txt'
INFILEC='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/CRUTEM4_StationList_APR2015.txt'
INFILEG='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/ghcnm.tavg.v3.2.2.20150430.qca.inv'

# Output Files
OUTPLOTH='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/StationDensity_HadISDH_JUN2015'
OUTPLOTC='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/StationDensity_CRUTEM4_JUN2015'
OUTPLOTG='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/StationDensity_GHCNM3_JUN2015'
OUTPLOTA='/data/local/hadkw/HADCRUH2/UPDATE2014/IMAGES/OTHER/StationDensity_ALL_JUN2015'

OUTLIST='/data/local/hadkw/HADCRUH2/UPDATE2014/LISTS_DOCS/StationDensityAllGBs_HadISDHCRUTEMGHCNM_JUN2015.txt'

# Variables
CRUTEM_count=0
HadISDH_count=0
GHCNM_count=0

CRUTEM_stations=list()	# list of arrays to store CRUTEM station info
HadISDH_stations=list()	# list of arrays to store HadISDH station info
GHCNM_stations=list()

DataLabels=['HadISDH','CRUTEM4','GHCNM3']

nlats=36
nlons=72
LatList=np.flipud((np.arange(nlats)*5.)-87.5)
LonList=(np.arange(nlons)*5)-177.5

HadISDH_grids=np.zeros((nlats,nlons))
CRUTEM_grids=np.zeros((nlats,nlons))
GHCNM_grids=np.zeros((nlats,nlons))
CRUHAD_grids=np.zeros((nlats,nlons))
GHCNMHAD_grids=np.zeros((nlats,nlons))
CRUHAD_grids.fill(-999)
GHCNMHAD_grids.fill(-999)

#************************************************************************
# Subroutines
#************************************************************************
# READDATA
def ReadData(FileName,typee,delimee):
    ''' Use numpy genfromtxt reading to read in all rows from a complex array '''
    ''' Need to specify format as it is complex '''
    ''' outputs an array of tuples that in turn need to be subscripted by their names defaults f0...f8 '''
    return np.genfromtxt(FileName, dtype=typee,delimiter=delimee) # ReadData

#************************************************************************
# FINDSTATIONS
def FindStationsList(StationList,Slat,Nlat,Wlon,Elon,FoundStationsList):
    ''' Loop through the given list of station info '''
    ''' Append to FoundStationList all info when a station is found within the box '''
    ''' Station lat/lon can be within or equal to gridbox boundaries '''
    # find lats and lons that are within the box
    # now find where these match (may be able to quadruple loop the where
    Gots=np.where((StationList[0] >= Wlon) & (StationList[0] <= Elon) & (StationList[1] >= Slat) & (StationList[1] <= Nlat))[0]
    
    if (len(Gots) > 0):
        for loo in range (5):
	    FoundStationsList.append(StationList[loo][Gots])
    
    return FoundStationsList # FindStations
#************************************************************************
# WriteOut
def WriteOut(TheFile,Lat,Lon,TheLabels,FoundH,FoundC,FoundG):
    ''' Append station info out to file given '''

    filee=open(TheFile,'a+')
    countH=0
    countC=0
    countG=0
    if (len(FoundH) > 0):
        countH=len(FoundH[0])
    if (len(FoundC) > 0):
        countC=len(FoundC[0])
    if (len(FoundG) > 0):
        countG=len(FoundG[0])
    
    filee.write(str('Long: '+'{:8.3f}'.format(Lon)+' Lat: '+'{:7.3f}'.format(Lat)+' HadISDH: '+'{:3d}'.format(countH)+' CRUTEM: '+'{:3d}'.format(countC)+' GHCNM: '+'{:3d}'.format(countG)+'\n'))
    
    if (countH > 0):
        for outt in range(len(FoundH[0])):
            filee.write(str('{:11s}'.format(FoundH[4][outt])+' '+'{:8.3f}'.format(FoundH[0][outt])+' '+'{:7.3f}'.format(FoundH[1][outt])+' '+'{:7.1f}'.format(FoundH[2][outt])+' '+'{:30s}'.format(FoundH[3][outt])+' '+'{:7s}'.format(TheLabels[0])+'\n'))

    if (countC > 0):
        for outt in range(len(FoundC[0])):
            filee.write(str('{:11s}'.format(FoundC[4][outt])+' '+'{:8.3f}'.format(FoundC[0][outt])+' '+'{:7.3f}'.format(FoundC[1][outt])+' '+'{:7.1f}'.format(FoundC[2][outt])+' '+'{:30s}'.format(FoundC[3][outt])+' '+'{:7s}'.format(TheLabels[1])+'\n'))

    if (countG > 0):
        for outt in range(len(FoundG[0])):
            filee.write(str('{:11s}'.format(FoundG[4][outt])+' '+'{:8.3f}'.format(FoundG[0][outt])+' '+'{:7.3f}'.format(FoundG[1][outt])+' '+'{:7.1f}'.format(FoundG[2][outt])+' '+'{:30s}'.format(FoundG[3][outt])+' '+'{:7s}'.format(TheLabels[2])+'\n'))

    filee.close()

    return #WriteOut
#************************************************************************
# PlotStationDensity
def PlotStationDensity(TheFile,DataGrids,TheLats,TheLons,Lablees):
    ''' Plot a map with colour showing number of stations '''
    ''' Label plot with totals for 1, 2, 5, 10+'''
    ''' Save as png and eps '''

    # set up dimensions and plot - this is a 2 column nvar rows plot
    f=plt.figure(1,figsize=(6,4))	#6,18
    plt.clf()

    plt.axes([0.1,0.12,0.8,0.8])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
      
    # set up colour scheme
    cmap=plt.get_cmap('YlGn') # PiYG
        
    cmaplist=[cmap(i) for i in range(cmap.N)]
#    for loo in range((cmap.N/2)-30,(cmap.N/2)+30):
#        cmaplist.remove(cmaplist[(cmap.N/2)-30]) # remove the very pale colours in the middle
    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle

# remove the darkest and lightest (white and black) - and reverse
    for loo in range(40):
        cmaplist.remove(cmaplist[0])
#    cmaplist.reverse()
#    for loo in range(10):
#        cmaplist.remove(cmaplist[0])
#    cmaplist.reverse()

    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    bounds=np.array([1,2,3,4,5,10,20,100])
    strbounds=["%3d" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)
    
    # make a grid of lats and lons
    #print(TheLons,TheLats)
    LngArrLons,LngArrLats=np.meshgrid(np.append(TheLons-2.5,180.),np.append(90.,TheLats-2.5))
        
    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    MSKData=0
    MSKData=ma.masked_where(DataGrids == 0,DataGrids)
    grids=m.pcolor(LngArrLons,LngArrLats,MSKData,cmap=cmap,norm=norm,latlon='TRUE')
    
#    plt.annotate(str(len(H_only))+" HadISDH stations",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)
#    plt.annotate(str(len(H_only))+" CRUTEM4 stations",xy=(0.07,0.5),xytext=None, xycoords='axes fraction',color="Black",size=12)

    cbax=f.add_axes([0.1,0.11,0.8,0.03])
    cb=plt.colorbar(grids,cax=cbax,orientation='horizontal',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=12) 
    cb.ax.set_xticklabels(strbounds)
    plt.figtext(0.5,0.02,'No. Stations',size=14,ha='center')    

    print(Lablees)
    print('1 station: ',len(np.where(DataGrids == 1)[0]), len(np.where(DataGrids == 1)[0])/float(len(np.where(DataGrids > 0)[0])))
    print('2 stations: ',len(np.where(DataGrids == 2)[0]), len(np.where(DataGrids == 2)[0])/float(len(np.where(DataGrids > 0)[0])))
    print('3 to 4 stations: ',len(np.where((DataGrids >= 3) & (DataGrids <= 4))[0]), len(np.where((DataGrids >= 3) & (DataGrids <= 4))[0])/float(len(np.where(DataGrids > 0)[0])))
    print('5 to 9 stations: ',len(np.where((DataGrids >= 5) & (DataGrids <= 9))[0]), len(np.where((DataGrids >= 5) & (DataGrids <= 9))[0])/float(len(np.where(DataGrids > 0)[0])))
    print('10+ stations: ',len(np.where(DataGrids >= 10)[0]), len(np.where(DataGrids >= 10)[0])/float(len(np.where(DataGrids > 0)[0])))

    plt.figtext(0.5,0.91,Lablees,size=18,ha='center')    
    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)

    #stop
    #plt.show()
    #stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotStationDensity
#************************************************************************
# PlotStationDensityTriple
def PlotStationDensityTriple(TheFile,DataGridsH,DataGridsCH,DataGridsGH,TheLats,TheLons,Lablees):
    ''' Plot three maps with colour showing number of stations and then differences '''
    ''' Label plot with totals for 1, 2, 5, 10+'''
    ''' Save as png and eps '''

    # make a grid of lats and lons
    LngArrLons,LngArrLats=np.meshgrid(np.append(TheLons-2.5,180.),np.append(90.,TheLats-2.5))

    # set up dimensions and plot - this is a 2 column nvar rows plot
    nplots=3
    xpos=[0.05,0.05,0.05]
    ypos=[0.67,0.34,0.01]
    xfat=[0.8,0.8,0.8]
    ytall=[0.29,0.29,0.29]
      
    f=plt.figure(3,figsize=(8,12))	#6,18
    plt.clf()

    # set up colour scheme
    cmap=plt.get_cmap('YlGn') # PiYG
        
    cmaplist=[cmap(i) for i in range(cmap.N)]
#    for loo in range((cmap.N/2)-30,(cmap.N/2)+30):
#        cmaplist.remove(cmaplist[(cmap.N/2)-30]) # remove the very pale colours in the middle
    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle

# remove the darkest and lightest (white and black) - and reverse
    for loo in range(40):
        cmaplist.remove(cmaplist[0])
#    cmaplist.reverse()
#    for loo in range(10):
#        cmaplist.remove(cmaplist[0])
#    cmaplist.reverse()

    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    bounds=np.array([1,2,3,4,5,10,20,100])
    strbounds=["%3d" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    plt.axes([xpos[0],ypos[0],xfat[0],ytall[0]])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
                  
    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    MSKData=ma.masked_where(DataGridsH == 0,DataGridsH)
    grids=m.pcolor(LngArrLons,LngArrLats,MSKData,cmap=cmap,norm=norm,latlon='TRUE')
    
    plt.figtext(0.45,0.97,Lablees[0],size=18,ha='center')    

    cbax=f.add_axes([0.87,0.67,0.03,0.29])
    cb=plt.colorbar(grids,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=12) 
    cb.ax.set_xticklabels(strbounds)
    plt.figtext(0.9,0.97,'No. Stations',size=14,ha='center')    


    # set up colour scheme
    cmap=plt.get_cmap('PiYG') # PiYG
        
    cmaplist=[cmap(i) for i in range(cmap.N)]
#    for loo in range((cmap.N/2)-30,(cmap.N/2)+30):
#        cmaplist.remove(cmaplist[(cmap.N/2)-30]) # remove the very pale colours in the middle
    #cmaplist.remove(cmaplist[(cmap.N/2)-10:(cmap.N/2)+10]) # remove the very pale colours in the middle

# remove the darkest and lightest (white and black) - and reverse
##    for loo in range(30):
##        cmaplist.remove(cmaplist[0])
#    cmaplist.reverse()
#    for loo in range(10):
#        cmaplist.remove(cmaplist[0])
#    cmaplist.reverse()

    cmap=cmap.from_list('this_cmap',cmaplist,cmap.N)
    
    bounds=np.array([-10,-5,-3,-2,-1,0,1,2,3,5,10])
    strbounds=["%3d" % i for i in bounds]
    norm=mpl_cm.colors.BoundaryNorm(bounds,cmap.N)

    plt.axes([xpos[1],ypos[1],xfat[1],ytall[1]])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
                  
    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    MSKData=ma.masked_where(DataGridsCH == -999,DataGridsCH)
    grids=m.pcolor(LngArrLons,LngArrLats,MSKData,cmap=cmap,norm=norm,latlon='TRUE')
    
    plt.figtext(0.45,0.64,Lablees[1]+'-'+Lablees[0],size=18,ha='center')    

    cbax=f.add_axes([0.87,0.34,0.03,0.29])
    cb=plt.colorbar(grids,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=12) 
    cb.ax.set_xticklabels(strbounds)
    plt.figtext(0.9,0.64,'No. Stations',size=14,ha='center')    

#    print(Lablees[1],' minus ',Lablees[0])
#    print('<=-10 station: ',len(np.where(MSKData <= -10)[0]), len(np.where(MSKData <= -10)[0])/float(len(np.where(MSKData)[0])))
#    print('<=-5 stations: ',len(np.where((MSKData > -10) & (MSKData <= -5))[0]), len(np.where((MSKData > -10) & (MSKData <= -5))[0])/float(len(np.where(MSKData)[0])))
#    print('<0 stations: ',len(np.where((MSKData > -5) & (MSKData < 0))[0]), len(np.where((MSKData > -5) & (MSKData < 0))[0])/float(len(np.where(MSKData)[0])))
#    print('0 stations: ',len(np.where((MSKData == 0) & (DataGridsH != 0))[0]), len(np.where((MSKData == 0) & (DataGridsH != 0))[0])/float(len(np.where(MSKData)[0])))
#    print('>0 stations: ',len(np.where((MSKData > 0) & (MSKData < 5))[0]), len(np.where((MSKData > 0) & (MSKData < 5))[0])/float(len(np.where(MSKData)[0])))
#    print('>=5 stations: ',len(np.where((MSKData >= 5) & (MSKData < 10))[0]), len(np.where((MSKData >= 5) & (MSKData < 10))[0])/float(len(np.where(MSKData)[0])))
#    print('>=10 stations: ',len(np.where(MSKData >= 10)[0]), len(np.where(MSKData >= 10)[0])/float(len(np.where(MSKData)[0])))

    plt.axes([xpos[2],ypos[2],xfat[2],ytall[2]])
     
    # plot map without continents and coastlines
    m = Basemap(projection='kav7',lon_0=0)
    # draw map boundary, transparent
    m.drawmapboundary()
    m.drawcoastlines()
    # draw paralells and medians, no labels
    m.drawparallels(np.arange(-90,90.,30.))
    m.drawmeridians(np.arange(-180,180.,60.))
                  
    # Only pcolor can do boxing on masked arrays, not really sure why we use pcolormesh at all
    MSKData=ma.masked_where(DataGridsGH == -999,DataGridsGH)

    grids=m.pcolor(LngArrLons,LngArrLats,MSKData,cmap=cmap,norm=norm,latlon='TRUE')
    
    plt.figtext(0.45,0.31,Lablees[2]+'-'+Lablees[0],size=18,ha='center')    

    cbax=f.add_axes([0.87,0.01,0.03,0.29])
    cb=plt.colorbar(grids,cax=cbax,orientation='vertical',ticks=bounds) #, extend=extend
    cb.ax.tick_params(labelsize=12) 
    cb.ax.set_xticklabels(strbounds)
    plt.figtext(0.9,0.31,'No. Stations',size=14,ha='center')    

#    print(Lablees[2],' minus ',Lablees[0])
#    print('<=-10 station: ',len(np.where(MSKData <= -10)[0]), len(np.where(MSKData <= -10)[0])/float(len(np.where(MSKData)[0])))
#    print('<=-5 stations: ',len(np.where((MSKData > -10) & (MSKData <= -5))[0]), len(np.where((MSKData > -10) & (MSKData <= -5))[0])/float(len(np.where(MSKData)[0])))
#    print('<0 stations: ',len(np.where((MSKData > -5) & (MSKData < 0))[0]), len(np.where((MSKData > -5) & (MSKData < 0))[0])/float(len(np.where(MSKData)[0])))
#    print('0 stations: ',len(np.where((MSKData == 0) & (DataGridsH != 0))[0]), len(np.where((MSKData == 0) & (DataGridsH != 0))[0])/float(len(np.where(MSKData)[0])))
#    print('>0 stations: ',len(np.where((MSKData > 0) & (MSKData < 5))[0]), len(np.where((MSKData > 0) & (MSKData < 5))[0])/float(len(np.where(MSKData)[0])))
#    print('>=5 stations: ',len(np.where((MSKData >= 5) & (MSKData < 10))[0]), len(np.where((MSKData >= 5) & (MSKData < 10))[0])/float(len(np.where(MSKData)[0])))
#    print('>=10 stations: ',len(np.where(MSKData >= 10)[0]), len(np.where(MSKData >= 10)[0])/float(len(np.where(MSKData)[0])))


    
    # add labals and watermark    
#    watermarkstring="/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename( __file__ )+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
#    plt.figtext(0.01,0.01,watermarkstring,size=6)


#    plt.show()
#    stop()

    plt.savefig(TheFile+".eps")
    plt.savefig(TheFile+".png")
     
    return #PlotStationDensityTriple    
    
#************************************************************************
# MAIN PROGRAM
#************************************************************************
# read in station list for HadISDH
MyTypes=("|S11","float","float","float","|S4","|S30","|S7")
MyDelimiters=[11,9,10,7,4,30,7]
MyFile=INFILEH
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
HadISDH_stations.append(np.array(RawData['f2']))	# longitudes
HadISDH_stations.append(np.array(RawData['f1']))	# latitudes
HadISDH_stations.append(np.array(RawData['f3']))	# elevations
HadISDH_stations.append(np.array(RawData['f5']))	# names
HadISDH_stations.append(np.array(RawData['f0']))	# ids
HadISDH_count=len(HadISDH_stations[0])

# read in station list for CRUTEM4 and append to list of arrays
MyTypes=("|S6","|S3","|S20","|S15","float","float","int","int","int") # need to remove '--' characters from station name
MyDelimiters=[6,3,20,15,7,7,6,5,5]
MyFile=INFILEC
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
CRUTEM_stations.append(-(np.array(RawData['f5'])))	# longitudes WHY ARE THESE IN REVERSE???
CRUTEM_stations.append(np.array(RawData['f4']))	# latitudes
CRUTEM_stations.append(np.array(RawData['f6']))	# elevations
CRUTEM_stations.append(np.array(RawData['f2']))	# names
CRUTEM_stations.append(np.array(RawData['f0']))	# ids
CRUTEM_count=len(CRUTEM_stations[0])

# read in station list for GHCNM3 and append to list of arrays
MyTypes=("|S11","float","float","float","|S1","|S20","|S49") # need to remove '--' characters from station name
MyDelimiters=[11,9,10,7,1,20,49]
MyFile=INFILEG
RawData=ReadData(MyFile,MyTypes,MyDelimiters)
GHCNM_stations.append(np.array(RawData['f2']))	# longitudes WHY ARE THESE IN REVERSE???
GHCNM_stations.append(np.array(RawData['f1']))	# latitudes
GHCNM_stations.append(np.array(RawData['f3']))	# elevations
GHCNM_stations.append(np.array(RawData['f5']))	# names
GHCNM_stations.append(np.array(RawData['f0']))	# ids
GHCNM_count=len(GHCNM_stations[0])
print('Read in data...')

# Loop through each gridbox and then through each station list to find stations within
for ltt in range(nlats):
    for lnn in range(nlons):
        latS=LatList[ltt]-2.5
        latN=LatList[ltt]+2.5
        lonW=LonList[lnn]-2.5
        lonE=LonList[lnn]+2.5
    
        Found_Hs=list() # list of arrays for ID, Lat, Lon, Elev, Name
        Found_Cs=list() # list of arrays for ID, Lat, Lon, Elev, Name
        Found_Gs=list() # list of arrays for ID, Lat, Lon, Elev, Name

        # Loop through HadISDH
        Found_Hs=FindStationsList(HadISDH_stations,latS,latN,lonW,lonE,Found_Hs)
        # Loop through CRUTEM
        Found_Cs=FindStationsList(CRUTEM_stations,latS,latN,lonW,lonE,Found_Cs)    
        # Loop through GHCNM
        Found_Gs=FindStationsList(GHCNM_stations,latS,latN,lonW,lonE,Found_Gs)
        
	count=0
	if (len(Found_Hs) > 0):
            count=len(Found_Hs[0])
	    HadISDH_grids[ltt,lnn]=count
	    if (len(Found_Cs) > 0):
	        CRUHAD_grids[ltt,lnn]=len(Found_Cs[0])-count
	    if (len(Found_Gs) > 0):
	        GHCNMHAD_grids[ltt,lnn]=len(Found_Gs[0])-count
        print('HadISDH: ',count)
    
        count=0
        if (len(Found_Cs) > 0):
            count=len(Found_Cs[0])
	    CRUTEM_grids[ltt,lnn]=count
	print('CRUTEM: ',count)
    
        count=0
        if (len(Found_Gs) > 0):
            count=len(Found_Gs[0])
	    GHCNM_grids[ltt,lnn]=count
        print('GHCNM: ',count)

        # Write to file by appending
        print('Writing to file....',lonW,lonE,latS,latN,ltt,lnn)
        WriteOut(OUTLIST,LatList[ltt],LonList[lnn],DataLabels,Found_Hs,Found_Cs,Found_Gs)


PlotStationDensity(OUTPLOTH,HadISDH_grids,LatList,LonList,DataLabels[0])
PlotStationDensity(OUTPLOTC,CRUTEM_grids,LatList,LonList,DataLabels[1])
PlotStationDensity(OUTPLOTG,GHCNM_grids,LatList,LonList,DataLabels[2])

PlotStationDensityTriple(OUTPLOTA,HadISDH_grids,CRUHAD_grids,GHCNMHAD_grids,LatList,LonList,DataLabels)

print(DataLabels[1],' minus ',DataLabels[0])
print('<=-10 station: ',len(np.where((CRUHAD_grids > -999) & (CRUHAD_grids <= -10))[0]),len(np.where((CRUHAD_grids > -999) & (CRUHAD_grids <= -10))[0])/float(len(np.where(CRUHAD_grids > -999)[0])))
print('<=-5 stations: ',len(np.where((CRUHAD_grids > -10) & (CRUHAD_grids <= -5))[0]), len(np.where((CRUHAD_grids > -10) & (CRUHAD_grids <= -5))[0])/float(len(np.where(CRUHAD_grids > -999)[0])))
print('<0 stations: ',len(np.where((CRUHAD_grids > -5) & (CRUHAD_grids < 0))[0]), len(np.where((CRUHAD_grids > -5) & (CRUHAD_grids < 0))[0])/float(len(np.where(CRUHAD_grids > -999)[0])))
print('0 stations: ',len(np.where(CRUHAD_grids == 0)[0]), len(np.where(CRUHAD_grids == 0)[0])/float(len(np.where(CRUHAD_grids > -999)[0])))
print('>0 stations: ',len(np.where((CRUHAD_grids > 0) & (CRUHAD_grids < 5))[0]), len(np.where((CRUHAD_grids > 0) & (CRUHAD_grids < 5))[0])/float(len(np.where(CRUHAD_grids > -999)[0])))
print('>=5 stations: ',len(np.where((CRUHAD_grids >= 5) & (CRUHAD_grids < 10))[0]), len(np.where((CRUHAD_grids >= 5) & (CRUHAD_grids < 10))[0])/float(len(np.where(CRUHAD_grids > -999)[0])))
print('>=10 stations: ',len(np.where(CRUHAD_grids >= 10)[0]), len(np.where(CRUHAD_grids >= 10)[0])/float(len(np.where(CRUHAD_grids > -999)[0])))

print(DataLabels[2],' minus ',DataLabels[0])
print('<=-10 station: ',len(np.where((GHCNMHAD_grids > -999) & (GHCNMHAD_grids <= -10))[0]),len(np.where((GHCNMHAD_grids > -999) & (GHCNMHAD_grids <= -10))[0])/float(len(np.where(GHCNMHAD_grids > -999)[0])))
print('<=-5 stations: ',len(np.where((GHCNMHAD_grids > -10) & (GHCNMHAD_grids <= -5))[0]), len(np.where((GHCNMHAD_grids > -10) & (GHCNMHAD_grids <= -5))[0])/float(len(np.where(GHCNMHAD_grids > -999)[0])))
print('<0 stations: ',len(np.where((GHCNMHAD_grids > -5) & (GHCNMHAD_grids < 0))[0]), len(np.where((GHCNMHAD_grids > -5) & (GHCNMHAD_grids < 0))[0])/float(len(np.where(GHCNMHAD_grids > -999)[0])))
print('0 stations: ',len(np.where(GHCNMHAD_grids == 0)[0]), len(np.where(GHCNMHAD_grids == 0)[0])/float(len(np.where(GHCNMHAD_grids > -999)[0])))
print('>0 stations: ',len(np.where((GHCNMHAD_grids > 0) & (GHCNMHAD_grids < 5))[0]), len(np.where((GHCNMHAD_grids > 0) & (GHCNMHAD_grids < 5))[0])/float(len(np.where(GHCNMHAD_grids > -999)[0])))
print('>=5 stations: ',len(np.where((GHCNMHAD_grids >= 5) & (GHCNMHAD_grids < 10))[0]), len(np.where((GHCNMHAD_grids >= 5) & (GHCNMHAD_grids < 10))[0])/float(len(np.where(GHCNMHAD_grids > -999)[0])))
print('>=10 stations: ',len(np.where(GHCNMHAD_grids >= 10)[0]), len(np.where(GHCNMHAD_grids >= 10)[0])/float(len(np.where(GHCNMHAD_grids > -999)[0])))
		
#    stop()

print("And, we are done!")

