#!/usr/local/sci/bin/python

#***************************************
# 12 August 2015 KMW - v1
# 
# Loop through ISTI station inventory used for benchmarks 
# If it doesn't have a country name then:
#	- find that station in the distance file and pull out neighbours
#       - find the nearest neighbour with the same first two letters and within 1000km that has a Country listed
# 	- OR else find the nearest neighbour within 500km that has a country listed
#	- or else - fill as NONE
# Also test that the country name is not one of the 32 duplicate names (SortingISTICountries.xlsx)
# If it is then reconcile these with one unique country name - results in 219 countries/networks plus NONE
# output:
#       Number of NONE stations
# 	Number of pseudo-country stations
# 	New list of stations with countries
# 
#************************************************************************
#                                 START
#************************************************************************
# USE python2.7
# python2.7 GetISTICountries_AUG2015.py
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

from Geography import TwoPointDistanceKm

# RESTART VALUE
Restarter='------'	#'------'		#'    01010000' station ID

# Set up file locations
STATLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat'
INDIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat'
OUTLIST='/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015_countries.dat'

nstations=0	# defined after reading in station list
ngoods=0	# stations passing unique location criteria
nbads=0		# shortest stations with matching locs
StationIDs=[]	# nstations list filled after reading in station list
StationLats=[]	# nstations list filled after reading in station list
StationLons=[]	# nstations list filled after reading in station list

StatDistsAll=[]
StatDistsSorted=[]
StatIDsSorted=[]

# tuples of desired country names for duplicates, and duplicate country names
DuplicateCountries=np.array(("AFGHANISTAN__ISLAMIC",  
 			"ALGERIA_/_ALGERIE",     
 			"BERMUDA_[UNITED_KING",  
 			"BOLIVIA_/_BOLIVIE",     
 			"BOSNIA_AND_HERZEGOVI",  
 			"BRAZIL_/_BRESIL",      
 			"CHILE_/_CHILI",         
 			"COTE_D'IVOIRE",         
 			"DOMINICAN_REPUBLIC_/",  
 			"EGYPT_/_EGYPTE",        
 			"EL_SALVADOR",
 			"ESTONIA_/_ESTONIE",     
 			"FRENCH POLYNESIA (FR",  
 			"GERMANY_/_ALLEMAGNE",   
 			"GREECE_/_GRECE",        
 			"INDONESIA_/_INDONESI",  
 			"JAMAICA_/_JAMAIQUE",    
 			"JAPAN_/_JAPON",         
 			"LIBYA_/_LIBYE",         
 			"MEXICO_/_MEXIQUE",      
 			"MOROCCO_/_MAROC",       
 			"NORWAY_/_NORVEGE",      
 			"RUSSIAN_FEDERATION_(",  
 			"SIERRA LEONE",
 			"SOLOMON ISLANDS",       
 			"SRI LANKA",
 			"THE BAHAMAS",
 			"THE GAMBIA",
 			"TRINIDAD AND TOBAGO",   
 			"TURKEY_/_TURQUIE",      
 			"UNITED_STATES_OF_AME",  
 			"WALLIS AND FUTUNA (F")) 

DesiredCountries=np.array(("AFGHANISTAN",
 			"ALGERIA",
 			"BERMUDA",
 			"BOLIVIA",
 			"BOSNIA",
 			"BRAZIL", 
 			"CHILE",  
 			"COTE D'IVOIRE",         
 			"DOMINICAN_REPUBLIC",    
 			"EGYPT",  
 			"EL SALVADOR",
 			"ESTONIA",
 			"FRENCH_POLYNESIA",      
 			"GERMANY",
 			"GREECE", 
 			"INDONESIA",
 			"JAMAICA",
 			"JAPAN",  
 			"LIBYA",  
 			"MEXICO", 
 			"MOROCCO",
 			"NORWAY", 
 			"RUSSIAN_FEDERATION",    
 			"SIERRA_LEONE",
 			"SOLOMON_ISLANDS",       
 			"SRI_LANKA",
 			"BAHAMAS,_THE",
 			"GAMBIA_/_GAMBIE",       
 			"TRINIDAD_AND_TOBAGO",   
 			"TURKEY", 
 			"UNITED_STATES",         
 			"WALLIS_AND_FUTUNA_[F"))  


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
# MAIN PROGRAM
#***********************************************************************
# read in station list
MyTypes=("|S12","|S32","|S23","float","float","float","int","int","int","int","int","int","|S16","|S17")
MyDelimiters=[12,32,23,8,11,9,5,5,5,5,5,5,16,17]
MyColumns='XXX'
RawData=ReadData(STATLIST,MyTypes,MyDelimiters,False,MyColumns)
StationIDs=np.array(RawData['f0'],ndmin=1)
StationLats=np.array(RawData['f3'],ndmin=1)
StationLons=np.array(RawData['f4'],ndmin=1)
StationCC=np.array(RawData['f2'],ndmin=1)
nstations=len(StationIDs)

CountNones=0	# Howmany NONES
CountPseudos=0 	# Howmany Pseudo countries?

print('Read in giant file',nstations)
# read in the distances/IDs for each station file
My_Fhandle=file(OUTLIST,'a')
f=open(INDIST,'r')
lines=f.readlines()
f.close()
for lop in range(nstations):
    		# moo=list(RawDists[lop])
    if Restarter != '------' and Restarter != StationIDs[lop]:
        continue
    else: 
        Restarter='------'

   # Is there a country missing?
    print(StationIDs[lop],StationCC[lop])
    if (StationCC[lop][0] == ' '):
        print('A NONE!',lop)
	#if (StationIDs[lop].strip() == 'AOM00066152'):
	#    stop
        # go through each neighbour within 500km to see if it has a country  
	# careful because there are some oddities - e.g. Congo has 7 stations that are listed as USA and Australian - very strange, no evidence of military bases there
	# So, add a check to first look for another 'close' station with first three letters of ID matching
	# If that doesn't work then allow any station within 500km   
        linoo=lines[lop].split()
        MakeItStrings=np.reshape(np.array(linoo[1:201],dtype="|S13"),(100,2))
        #print(MakeItStrings[0:10,0])
        StatDistsIDsSorted=np.array(MakeItStrings[:,0])	
        StatDistsKMsSorted=np.array(MakeItStrings[:,1])	
	StationCC[lop]='NONE'    
	CountNones=CountNones+1               
        Switch_point=0
	#if (StationIDs[lop].strip() == 'AOM00066152'):
	#    stop
	for dd in range(100):
            #print(StatDistsIDsSorted[dd],np.float(StatDistsKMsSorted[dd]))  
	    #stop          
	    if (StationIDs[lop][1:3] == StatDistsIDsSorted[dd][0:2]) & (StationCC[np.where(StationIDs == ' '+StatDistsIDsSorted[dd])[0]][0][0] != ' ') & (np.float(StatDistsKMsSorted[dd]) <= 1000.):
		print("Matched a neighbour with ID")
		StationCC[lop]=StationCC[np.where(StationIDs == ' '+StatDistsIDsSorted[dd])[0]][0]
	        CountPseudos=CountPseudos+1
		CountNones=CountNones-1
		Switch_point=1
		break
        if (Switch_point == 0):
	    for dd in range(100):
                #print(StatDistsIDsSorted[dd],np.float(StatDistsKMsSorted[dd]))  
	        #stop          
	        if (StationCC[np.where(StationIDs == ' '+StatDistsIDsSorted[dd])[0]][0][0] != ' ') & (np.float(StatDistsKMsSorted[dd]) <= 500.):
		    print("Matched nearest neighbour regardless")
		    StationCC[lop]=StationCC[np.where(StationIDs == ' '+StatDistsIDsSorted[dd])[0]][0]
	            CountPseudos=CountPseudos+1
		    CountNones=CountNones-1
		    Switch_point=1
		    break
		
    # Check for non-desireable country name (a duplicate)
    if (len(np.where(DuplicateCountries == StationCC[lop].strip())[0]) > 0):
        #stop
        StationCC[lop]=DesiredCountries[DuplicateCountries == StationCC[lop].strip()][0]		
     
    # sort out RawData format for outputting station list
    outraw=list(RawData[lop])
    outraw[2]="{:23s}".format(StationCC[lop])
    outraw[3]="{:8.4f}".format(outraw[3])
    outraw[4]="{:11.4f}".format(outraw[4])
    outraw[5]="{:9.2f}".format(outraw[5])
    outraw[6:12]=["{:5d}".format(dd) for dd in outraw[6:12]]
    outraw=np.reshape(outraw,(1,len(outraw)))     
    np.savetxt(My_Fhandle,outraw,fmt='%s',delimiter='')

My_Fhandle.close()

print('Number of NONES: ',CountNones)
print('Number of Pseudo-countries: ',CountPseudos)
	
    #stop()

print("And, we are done!")
