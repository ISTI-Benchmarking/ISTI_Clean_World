# -*- coding: iso-8859-1 -*-
# Finding land use at an observing station in ISTI located at (stlat, stlong) using the ESA-CCI land use data set.  N/S are +/-;   W/E are -/+.
# LAND USE FLAG MEANINGS. VALUES ARE IN BYTE FORMAT, RANGE -128 to +128
#0 no_data 
#10 cropland_rainfed 
#11 herbaceous_cover 
#12 tree_or_shrub_cover 
#20 cropland_irrigated 
#30 mosaic_cropland 
#40 mosaic_natural_vegetation 
#50 tree_broadleaved_evergreen_closed_to_open
#60 tree_broadleaved_deciduous_closed_to_open 
#61 tree_broadleaved_deciduous_closed 
#62 tree_broadleaved_deciduous_open 
#70 tree_needleleaved_evergreen_closed_to_open 
#71 tree_needleleaved_evergreen_closed 
#72 tree_needleleaved_evergreen_open 
#80 tree_needleleaved_deciduous_closed_to_open 
#81 tree_needleleaved_deciduous_closed 
#82 tree_needleleaved_deciduous_open 
#90 tree_mixed 
#100 mosaic_tree_and_shrub 
#110 mosaic_herbaceous 
#120 shrubland 
#-126 grassland 
#-116 lichens_and_mosses 
#-106 sparse_vegetation 
#-96 tree_cover_flooded_fresh_or_brakish_water 
#-86 tree_cover_flooded_saline_water 
#-76 shrub_or_herbaceous_cover_flooded 
#-66 urban 
#-56 bare_areas 
#-46 inland_water 
#-45 ocean_water 
#-36 snow_and_ice

#  
# Run using cd PROGRAMS/Python/Urbanisation; python2.7; import station_9pixels_landuse_ISTI as stlu_isti then linestart = 0; nlines = 100; stlu_isti.station_landtype (linestart, nlines)
# This processes nlines stations beginning at linestart in the input station file which has one line per station
def station_landtype (linestart, nlines):
  import iris
  import numpy as np 
  cubeg = iris.load('/project/LandCoverCCI/V1/ESACCI-LC-Map-nc/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.1-inlandwater.nc') # This loads the global file 
  cubeg0 = cubeg[0] 
  latvalues = cubeg0.coord('latitude').points
  nlats = len(cubeg0.coord('latitude').points) # nlats is 64800
  longvalues = cubeg0.coord('longitude').points 
# Latitude decreases (moves south) but longitude increases (moves east) with increasing array element. The increments are approximately 0.002777 degrees of latitude or longitude but are not absolutely fixed.
  nlongs = len(cubeg0.coord('longitude').points) # nlongs is 129600
  dy_max = 0.00278473 # These are the maximum and minimum latitude increments in the land use dataset 
  dy_min = 0.0027771
  dx_max = 0.00279236 # These are the maximum and minimum longitude increments in the land use dataset 
  dx_min = 0.0027771

#  station_in1 = open('/home/h06/hadux/DATA/LandAirTemperatures/UrbanStudiesStationLists/ISTIinventory_monthly_merged_stage3.txt', 'r+b') # Station locations to 4 decimal places of latitude and longitude
#  station_out = open('/home/h06/hadux/DATA/LandAirTemperatures/UrbanStudiesStationLists/ISTI_monthlymergedstage4_Statlocs_landuse_and_urbansignal.txt', 'w')
  station_in1 = open('/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015_countries.dat', 'r+b')
#  station_in1.readline()  # scroll through 2 initial documentation lines in the input file of station locations
#  station_in1.readline()
  if linestart > 0:
    for skip in range(linestart):
       station_in1.readline()    # skip the first linestart lines of data.
  station_out = open('/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015_landuse3.dat', 'a')
  station_out.write('International Surface Temperature Initiative (ISTI) monthly merged stage 4 station locations (degrees to 4 decimal places; south and west are negative); land use from ESA-CCI  \n')    
  station_out.write(' Land use code: 0 no_data; 10 cropland_rainfed; 11 herbaceous_cover; 12 tree_or_shrub_cover; 20 cropland_irrigated; 30 mosaic_cropland; 40 mosaic_natural_vegetation \n')
  station_out.write(' 50 tree_broadleaved_evergreen_closed_to_open; 60 tree_broadleaved_deciduous_closed_to_open; 61 tree_broadleaved_deciduous_closed; 62 tree_broadleaved_deciduous_open \n')
  station_out.write(' 70 tree_needleleaved_evergreen_closed_to_open; 71 tree_needleleaved_evergreen_closed; 72 tree_needleleaved_evergreen_open \n')
  station_out.write(' 80 tree_needleleaved_deciduous_closed_to_open; 81 tree_needleleaved_deciduous_closed; 82 tree_needleleaved_deciduous_open; 90 tree_mixed; 100 mosaic_tree_and_shrub; 110 mosaic_herbaceous; 120 shrubland \n')
  station_out.write(' -126 grassland; -116 lichens_and_mosses; -106 sparse_vegetation; -96 tree_cover_flooded_fresh_or_brakish_water; -86 tree_cover_flooded_saline_water; -76 shrub_or_herbaceous_cover_flooded  \n')
  station_out.write(' -66 urban; -56 bare_areas; -46 inland_water; -45 ocean_water; -36 snow_and_ice \n')
  station_out.write('  Station   Latitude   Longitude             Station name, country                         Land use (-66 = urban)    a. 1.0 = urban; 0.0 = non-urban pixel  b. Simple average of 9 pixels centred on station \n')
  station_out.write('            c. As b but weight 0.6, 0.15/(1+cos(lat)), 0.15*cos(lat)/(1+cos(lat)), 0.025 for centre, E & W side pixels,  N & S side pixels, & corner pixels     d. Signal = c*function of neighbouring 8 pixels not water, barren or sparse veg.    NW etc: land use code in 9 pixels centred on the station  \n')
  station_out.write('                                                                                                        a       b        c        d         NW       N      NE       W    target     E      SW       S      SE\n')
  station_out.close()

  blank = ' '
  blank4 = '    '
  blank5 = '     '
  blank9 = '         '

  for j in range(linestart, linestart+nlines):
#    blank4 = station_in1.readline(4)  # 4 leading blanks
    stid = station_in1.readline(12)  # station number  
    txt1 = station_in1.readline(55) # descriptive text 
    lat = station_in1.readline(8)  # latitude e.g.   1.2345 or -32.1834 with 0, 1, or 2 leading blanks
    blank2 = station_in1.readline(2)  # blanks 
    lon = station_in1.readline(9)  # longitude e.g.  1.2345 or -175.1504 with 0, 1, 2 or 3 leading blanks
    txt2 = station_in1.readline(60) # elevation, length of record and other information not to be used here
    station_in1.readline() # read line-throw

    stlat = float (lat) 
    rlat = np.radians(stlat) 
    stlong = float (lon)
    
    if stlat > -89.5: # treat the South pole separately later in this script

      imx = int(1.0 +(90.0 - stlat)/dy_min) # maximum possible latitude array-element corresponding to the station given the range of latitude increments
      imn = int(-1.0 +(90.0 - stlat)/dy_max) # minimum possible latitude array-element corresponding to the station given the range of latitude increments
      jmx = int(1.0 +(stlong + 180)/dx_min) # maximum possible longitude array-element corresponding to the station given the range of latitude increments
      jmn = int(-1.0 +(stlong + 180)/dx_max) # minimum possible longitude array-element corresponding to the station given the range of latitude increments
      print stid, blank, stlat, blank, stlong, blank, imn, imx, jmn, jmx 
# Now iterate through the land use tiles, finding the latitude and longitude ranges of each one until a tile matches the station location. 
# To reduce this loop, run from imn = (90-stlat)/dy_max to imx = (90-stlat)/dy_min and from jmn = (stlong+180)/dx_max to jmx = (stlong+180)dx_min, each one to the nearest integer mins(maxs) rounded down(up)

      cubeg1 = cubeg0[imn:imx+1, jmn:jmx+1]
      cubegd = cubeg1.data   # download only the data for the range of array-elements imn through imx, jmn through jmx
      stlanduse = -256 # this value does not occur in the input land-use dataset; if it is output, the station has not been found!
      pxlu9 = np.zeros((3,3), dtype = int)  # array of land use in station pixel and in surrounding pixels
      pxur9 = np.zeros((3,3), dtype = float)  # array of urban indicators (1.0 urban, 0.0 not urban) in station pixel and in surrounding pixels
      for i in range(imn, imx+1):
        for j in range(jmn, jmx+1): 
          latmax = 90.0 # North Pole
          latmin = -89.0 # most southerly latitude considered here
          longmax = 180.0 # Dateline at eastern edge
          longmin = -180.0 # Dateline at western edge
          if i != 0: latmax = latvalues[i] + 0.5*(latvalues[i-1]  - latvalues[i])
          if i != nlats-1: latmin = latvalues[i] - 0.5*(latvalues[i]  - latvalues[i+1])  
          if j != nlongs-1: longmax = longvalues[j] + 0.5*(longvalues[j+1]  - longvalues[j])
          if j != 0: longmin = longvalues[j] - 0.5*(longvalues[j]  - longvalues[j-1]) 
          if stlat <= latmax and stlat >= latmin and stlong <= longmax and stlong >= longmin: 
            for ii in range(3):
              for jj in range(3):
                pxlu9[ii,jj] = cubegd[i-imn+ii-1, j-jmn+jj-1]
                if pxlu9[ii,jj] == -66: pxur9[ii,jj] = 1.0  # pxur9 gives the urban status of each of the 9 pixels centred on the station
            break
# stlanduse is the desired land use at the station sited at (stlat, stlong) 

    if stlat <= -89.5:          # Fix for South pole: -36 denotes snow and ice.
      for ii in range(3):
        for jj in range(3):
          pxlu9[ii,jj] = -36
          pxur9[ii,jj] = 0.0
#    print pxlu9[1,1]    # This is the land use in the station pixel

    stlup = str (pxlu9[1,1]) # station pixel land use
    if len(stlup) == 1: stlup = '        ' + stlup
    if len(stlup) == 2: stlup = '       ' + stlup
    if len(stlup) == 3: stlup = '      ' + stlup
    if len(stlup) == 4: stlup = '     ' + stlup
    stur = str (pxur9[1,1]) # Station pixel urban indicator: (1.0 urban, 0.0 not urban)
    coslat = np.cos(rlat)  # cosine of latitude
    nswt = 0.15*coslat/(1.0+coslat)
    ewwt = 0.15/(1.0+coslat) # e.g. at 60N coslat=0.5: nswt=0.5, ewwt=0.1 are the weights of the north & south, east & west side boxes
    pxavur = (pxur9[0,0] + pxur9[0,1] + pxur9[0,2] + pxur9[1,0] + pxur9[1,1] + pxur9[1,2] + pxur9[2,0] + pxur9[2,1] + pxur9[2,2])/9.0
    stpxavur = str(round(pxavur, 2)) # Average urbanization in the 9-pixel block centred on the station.
    pxwtur = (0.025*pxur9[0,0] + nswt*pxur9[0,1] + 0.025*pxur9[0,2] + ewwt*pxur9[1,0] + 0.6*pxur9[1,1] + ewwt*pxur9[1,2] + 0.025*pxur9[2,0] + nswt*pxur9[2,1] + 0.025*pxur9[2,2]) 
    stpxwt = str(round(pxwtur, 2)) # Weighted average urbanization (centre 0.6; sides 0.3; corners 0.1) in the 9-pixel block centred on the station. 
                                   # East & west sides weighted more than north & south sides by factor 1/coslat
    watdes = 0.0 # indicator for summing the prevalence of water or desert in surrounding pixels
                 # then count neighbouring pixels with either water or desert are counted, with relative weighting the same as just above, thus having a maximum total weight of 0.4. 
    for ii in range(3):
      for jj in range(3):
        if pxlu9[ii,jj] == -45 or pxlu9[ii,jj] == -45 or pxlu9[ii,jj] == -56 or pxlu9[ii,jj] == -106: # water or desert or sparse veg pixels
          if (ii == 0 or ii == 2) and jj == 1: watdes = watdes + nswt
          if ii == 1 and (jj == 0 or jj == 2): watdes = watdes + ewwt
          if (ii == 0 and jj == 0) or (ii == 0 and jj == 2) or (ii == 2 and jj == 0) or (ii == 2 and jj == 2): watdes = watdes + 0.025
            # a proportion (total weight of water plus desert plus sparse veg pixels * 2.5) is subtracted from the urban influence coming to zero if city is fully surrounded by water or desert
    watdesur = 0.0
    if pxwtur > 0.0: watdesur = pxwtur*(1 - 2.5*watdes)
    stwdur = str(round(watdesur, 2)) 
    nxtspace = blank
    if len(stwdur) == 3: nxtspace = '  '           # formatting of output into tidy columns
    if len(stpxwt) == 3: stwdur = blank + stwdur   # formatting of output into tidy columns
    if len(stpxavur) == 3: stpxwt = blank + stpxwt # formatting of output into tidy columns
    textout = stid + blank4 + lat + blank + lon + txt1 + stlup + blank9 + stur + blank5 + stpxavur + blank5 + stpxwt + blank5 + stwdur + nxtspace
    
    sp = blank
    for ii in range(3):
      for jj in range(3):
        sp  = str(pxlu9[ii,jj])
        lenp = len(sp)
        if lenp == 1: sp = '       ' + sp 
        if lenp == 2: sp = '      ' + sp 
        if lenp == 3: sp = '     ' + sp
        if lenp == 4: sp = '    ' + sp # Land use code is 1 to 4 digits - see above documentation
        textout = textout + sp

    textout = textout + ' \n'
    station_out = open('/data/local/hadkw/ISTI/LISTS/v101_JUL2015/ISTILONGINVENTORY_stage3proxyelevs_JUL2015_landuse3.dat', 'a')
    station_out.write(textout)
    station_out.close()

