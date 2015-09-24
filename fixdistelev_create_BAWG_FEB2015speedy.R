# R program to redo VARs using distance and elevation as a proxy

source("BigMVNgen.r")


#####################################################
# This should create much higher correlations (more stable?)
# and allow more stations to be brought in to do the VAR
# because co-existant data is no longer essential

# Only use GHCN stations as these are likely higher in quality?
# Not really sure that's a fair assumption
# At least means its fewer stations

# Read in each station's 100 highest correlating neighbours:
# lag 0 (first diff clim anoms), lag 1 (clim anoms) and the station autocorrelation
# Read in the neighbour distances and sort to match with 100 highest corrs
# find those station and neighbour elevations to work out vertical distance

# plot ALL stations distance from candidate verses cross-corr
# fit exponential (best fit) and retain equation

# plot ALL stations distance from candidate verses
# cross-corr at lag 1
# fit exponential (best fit) and retain equation

# for each station, use 40 nearest neighbours (save)
# using their distance (vertical and horizontal) make up g0 and g1 using the equations
# above (save) and then the VARPS

# With any luck (doubtful) this should improve the stability
# and prevent craziness

# so from this, you can just swap the CORR, COVS, VARPS files and run as before.

# AUG 2015
# This now uses Standardised anomalies rather than climate anomalies for the cors at lag 0 and lag1
# This now uses a third of GHCN stations
# This now only uses candidate stations, and neighbours where each station's total counts of months are greater than or equal to 360 (30 years)
# This now only considers points where corrLag0 > 0.0 and corLag1 > 0.1
# This now includes 100 neighbours to try and stretch the coverage
# Use the 3 element fits of lag0 and lag1 to chose the intercept for the 2 element fits
# Use the mean of the first ~8000 points to set the diagonals for the autocorrelation - test by plotting OCC[1;9000] and cut off where they drop < 1.0 (~4000)
# Base the diagonals on mean of ACs (not lag1 cross correlations) - e.g., 0.22 (mean=0.214 so rounded up)
# Base the intercept on mean of all + 1 sd - e.g., 0.2 (0.1415 + 0.0689) rounded down
######################################################
# Functions
#----------------------------------------------------
# exp.corr
#   Function to define an exponential covariance structure. The 
#   covariance between sites separated by a distance of d km is
#   sigma^2 * exp(-phi*d). The vector theta below contains 
#   two elements: (sigma^2,phi). coords is a two-column
#   matrix of latitudes and longitudes
# KW: This substitutes a real covariance structure by assuming that 
# covariances decay exponentially with distance from 1 at 0 km to 
# 1/e(ish 0.37) at 1000km to 0.0 at 5000km 
exp.corr <- function(thetaplus,coords,elevations) {
  d <- howfar(coords,coords)
  h <- howhigh(elevations,elevations)
  z <- theta[1] * ( ( exp(-theta[2]*d)^2 + exp(-theta[3]*e)^2 ) / 2 )
  #print(z)
  #z[(z < 0.4)] <- 0.4
  #print(z)
  #stop()
  # force all non-diags to be <0.99
  # force diags to be 1.00
#  z[(z >= 0.989)]<-0.989
#  diag(z)<-1.00
  return(z)	#list(z=z,d=d))
}

#-----------------------------------------------------------
# howfar
# Function to calculate the distance between sets of points on the 
# earth's surface, with co-ordinates given in latitude and longitude.
# Arguments:
#
# sites1  A matrix with 2 columns, indicating the (lat, long)
#		coordinates of each site in the first set (in 
#		degrees!)
# sites2	The same for the second set of sites
# units		Either "km" (the default) or "nm" for distances in 
#		nautical miles (*not* nanometers!)
#
# Value:	A matrix of distances, with rows corresponding to 
#		the sites in sites1 and columns corresponding to 
#		those in sites2. If sites1 and sites2 have row names,
#		these are copied over as appropriate
howfar <- function(sites1,sites2,units="km") {
  n1 <- dim(sites1)[1]; n2 <- dim(sites2)[1]
  z <- matrix(nrow=n1,ncol=n2)
  rownames(z) <- rownames(sites1) 
  colnames(z) <- rownames(sites2) 
  pi180 <- pi/180
  for (i in 1:n1) {
    lat1 <- sites1[i,1]
    long1 <- sites1[i,2]
    lat2 <- sites2[,1]
    long2 <- sites2[,2]
    z[i,] <- sin(sites1[i,1]*pi180)*sin(sites2[,1]*pi180) + 
             cos(sites1[i,1]*pi180)*cos(sites2[,1]*pi180) * 
             cos((sites2[,2]-sites1[i,2])*pi180)
#
#   Next line deals with rounding errors that take the previous
#   line slightly outside the range (-1,1)
#
    z[i,] <- (60*180/pi)*acos(pmax(pmin(z[i,],1),-1))
  }
  if (units == "km") {
    z <- 1.852*z
  } else if (units != "nm") {
    stop("units must be either 'km' or 'nm'")
  }
  z
}
#-----------------------------------------------------------
# howhigh
# Function to calculate the vertical distance between sets of points on the 
# earth's surface, with co-ordinates given in height in km
# Arguments:
#
# sites1  A vector with elevations in km for all stations
# sites2  A vector with elevations in km for all stations
#
# Value:	A matrix of vertical distances, with rows corresponding to 
#		the sites in sites1 and columns corresponding to 
#		those in sites2. If sites1 and sites2 have row names,
#		these are copied over as appropriate
howhigh <- function(sites1,sites2) {
  n1 <- dim(sites1)[1]; n2 <- dim(sites2)[1]
  z <- matrix(nrow=n1,ncol=n2)
  rownames(z) <- rownames(sites1) 
  colnames(z) <- rownames(sites2) 
  for (i in 1:n1) {
    z[i,] <- abs(sites2-sites1[i])
    
    # should be positive if sites2 station is higher and negative if sites2 station is lower
    # BUT THIS IS NOT IMPORTANT - ONLY ABSOLUTE DISTANCE
  }
  z
}

#----------------------------------------------------
# read_station_func

#C read in an ISTI station file and pull out only the TMean (div by 100.)
#C filter into a fulltime array
#C return the uncompressed temperature data
#C clean up to save memory

read_station_func<-function(TheStYr,TheFilee,TheStation) {
  tmpyr<-0
  tmpmn<-0
  tmpvals<-0
  mush<-read.fwf(TheFilee,widths=c(61,5,2,2,6,6,6,71),comment.char=";") # each columns has one of the variables and R interprets the type (usually) correctly
  tmpyr<-mush[[2]][]
  tmpmn<-mush[[3]][]
  tmpvals<-mush[[7]][]/100.
  pointers<-((tmpyr[]-TheStYr)*12)+tmpmn[]	# should be an array from 1 to 1368
  mypoints<-array(seq(length(tmpyr)))
  mypoints<-mypoints[which(pointers > 0)]
  pointers<-pointers[which(pointers > 0)]
  TheStation[pointers]<-tmpvals[mypoints]
  rm(mush,tmpyr,tmpmn,tmpvals,pointers,mypoints)

  return(TheStation)
}
#--------------------------------------------------------------------------------------------

######################################################
# PARAMETERS
# TUNEABLE PARAMETERS

# EDITABLE PARAMETERS
nstations<-32522 	# Only stations with sufficient correlating neighbours to create a VAR model (22697,22256)

# missing data indicator
mdi	 <-(-99.99)

# SET IN STONE PARAMETERS

######################################################
# Files and Directories

#-----------------------------------------------------------
# SET UP FILE PATHS AND FILES @ enric: I had to modify this to run in 
dirlist<-"/data/local/hadkw/ISTI/LISTS/v101_JUL2015/"
outdir<-"/data/local/hadkw/ISTI/IMAGES/"
infillist<-"ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat"	# Actual VAR compatible station list (if parttone == 1)
infildists<-"ISTILONGDISTANCES_thousand_stage3proxyelevs_JUL2015.dat"	# Correlation station list
infilcounts<-"ISTILONGINVENTORY_counts_stage3proxyelevs_JUL2015.dat"	# Correlation station list

# Lag 0 Correlation Neighbour list for each station (<=40 highest correlating neighbours)
#inCORRlag0<-"ISTILONGCORRS_hundred_stage3proxyelevs_JUL2015.dat" 
#inCORRlag1<-"ISTILONGCORRSLAG1_hundred_stage3proxyelevs_JUL2015.dat" 
#inCORRAC<-"ISTILONGAUTOCORRS_hundred_stage3proxyelevs_JUL2015.dat" 
inCORRlag0<-"ISTILONGCORRSStAnom_hundred_stage3proxyelevs_JUL2015.dat" 
inCORRlag1<-"ISTILONGCORRSLAG1StAnom_hundred_stage3proxyelevs_JUL2015.dat" 
inCORRAC<-"ISTILONGAUTOCORRSStAnom_hundred_stage3proxyelevs_JUL2015.dat" 

# output plots
# 3 element fit CCl0 - fit to real data
outplotFit3CCl0=paste(outdir,"DistElevDecayFunc3CCl0_AUG2015stanom100N00AC.eps",sep="")
# 2 element fit CCl0 - fit to real data
outplotFit2CCl0=paste(outdir,"DistElevDecayFunc2CCl0_AUG2015stanom100N00AC.eps",sep="")
# 3 element fit CCl1 - fit to real data
outplotFit3CCl1=paste(outdir,"DistElevDecayFunc3CCl1_AUG2015stanom100N00AC.eps",sep="")
# 2 element fit CCl1 - fit to real data
outplotFit2CCl1=paste(outdir,"DistElevDecayFunc2CCl1_AUG2015stanom100N00AC.eps",sep="")
# 3 element fit CCl1 using ccl0 - fit to real data
outplotFit3CCl1l0=paste(outdir,"DistElevDecayFunc3CCl1_usel0_AUG2015stanom100N00AC.eps",sep="")
# 2 element fit CCl1 using ccl0 - fit to real data
outplotFit2CCl1l0=paste(outdir,"DistElevDecayFunc2CCl1_usel0_AUG2015stanom100N00AC.eps",sep="")

# 3 element fit CCl0 - scatter real vs estimated data
outplotScat3CCl0=paste(outdir,"ScatterDistElevDecayFunc3CCl0_AUG2015stanom100N00AC.eps",sep="")
# 2 element fit CCl0 - scatter real vs estimated data
outplotScat2CCl0=paste(outdir,"ScatterDistElevDecayFunc2CCl0_AUG2015stanom100N00AC.eps",sep="")
# 3 element fit CCl1 - scatter real vs estimated data
outplotScat3CCl1=paste(outdir,"ScatterDistElevDecayFunc3CCl1_AUG2015stanom100N00AC.eps",sep="")
# 2 element fit CCl1 - scatter real vs estimated data
outplotScat2CCl1=paste(outdir,"ScatterDistElevDecayFunc2CCl1_AUG2015stanom100N00AC.eps",sep="")
# 3 element fit CCl1 using ccl0 - scatter real vs estimated data
outplotScat3CCl1l0=paste(outdir,"ScatterDistElevDecayFunc3CCl1_usel0_AUG2015stanom100N00AC.eps",sep="")
# 2 element fit CCl1 using ccl0 - scatter real vs estimated data
outplotScat2CCl1l0=paste(outdir,"ScatterDistElevDecayFunc2CCl1_usel0_AUG2015stanom100N00AC.eps",sep="")

######################################################
# Variables and arrays
ds_info     <-list(statid=array("XXXXXX",dim=(nstations)),statlats=array(mdi,dim=(nstations)),
                   statlons=array(mdi,dim=(nstations)),statelvs=array(mdi,dim=(nstations)),
		   statsource=array("XXXXXX",dim=(nstations)),monthcounts=array(0,dim=(nstations)))

DistStatID<-0	# list of IDs in the Correlations file to search through to find right line to read

# maps for locating stations within the main list ds_info
corrneighbours.id<-vector("list",nstations)	# map using full list
ghcnstationsmap<-0 				# map for GHCN stations within full list

# while all real data will be stored in temporary arrays the output data must be stored for good or printed to file/read in again.
tmpACs<-0		# temporary scalar for station autocorrelation (climate anomalies with loess) 
tmpCClag1<-array(0,100)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)
tmpCClag0<-array(0,100)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)
tmpNEIGHS<-array(0,100)	# temporary array for a list of the neighbours
tmpHORIZDISTS<-array(0,100) # temporary array for a list of horizontal dists
tmpVERTDISTS<-array(0,100)  # temporary array for a list of vertical dists

AllCrossCorrs<-0		# This array is grown - first nGstations (<29377) are 1.0
AllCrossCorrsLag1<-0		# This array is grown - first nGstations (<29377) are ACs
AllDistances<-0			# This array is grown - first nGstations (<29377) are 0.0
AllVertDistances<-0		# This array is grown - first nGstations (<29377) are 0.0
AllCounts<-0			# This array is grown using actual count of each station - later used to filter (no check on 'commons' between stations though

######################################################
######################################################
########             MAIN PROGRAM              #######
######################################################
######################################################

# Read in stations and pull out GHCN ones
# Read in those neighbour correlations and distances 
# Match up distances
# Work out vertical distances

# plot ALL stations distance from candidate verses cross-corr
# fit exponential (best fit) and retain equation

# plot ALL stations distance from candidate verses
# cross-corr at lag 1
# fit exponential (best fit) and retain equation

# 1) READ IN LIST OF STATIONS
mush<-readLines(con=paste(dirlist,infillist,sep=""),n=-1)	# read entire file
bush<-readLines(con=paste(dirlist,infilcounts,sep=""),n=-1)	# read entire file
for (linoo in 1:nstations) {
  ds_info$statid[linoo]<-substr(mush[linoo],2,12)			# station ID
  ds_info$statlats[linoo]<-type.convert(substr(mush[linoo],68,75))	# 68:75station Latitude
  ds_info$statlons[linoo]<-type.convert(substr(mush[linoo],78,86))	# 78:86station Longitude
  ds_info$statelvs[linoo]<-type.convert(substr(mush[linoo],88,95))	# 88:95station elevnation
  ds_info$statsource[linoo]<-type.convert(substr(mush[linoo],127,128))	# 88:95station elevnation
  ds_info$monthcounts[linoo]<-type.convert(substr(bush[linoo],12,17))	# month counts
} 
rm(mush)
rm(bush)

# map reals to full
ghcnstationsmap<-(1:nstations)[ds_info$statsource == 1] # should now point to all GHCN stations
## too many so halve these
#ghcnstationsmap<-array(ghcnstationsmap[1:24810],c(2,12405))
#ghcnstationsmap<-ghcnstationsmap[1,]
# too many so third these (lose last station to make divisible by three)
ghcnstationsmap<-array(ghcnstationsmap[1:27087],c(3,9029))
ghcnstationsmap<-ghcnstationsmap[1,]
#-----------------------------------------------------------------
mooIN<-read.table(paste(dirlist,inCORRAC,sep=""),header=FALSE)
moo<-array(unlist(mooIN[2]))
GHCNACS<-moo[ghcnstationsmap]
moo<-0
nGstations<-length(GHCNACS)	# this should now be ~9000 not 27000

AllCrossCorrs<-array(replicate(nGstations,1.0)) 	# first 9029 points are correlation of the candidate station with itself (1.0)
AllCrossCorrsLag1<-array(GHCNACS)			# first 9029 points are Autocorrelation of the candidate station
AllDistances<-array(replicate(nGstations,0.0))		# first 9029 points are horizontal distance from candidate station to itself (0.0)
AllVertDistances<-array(replicate(nGstations,0.0))	# first 9029 points are horizontal distance from candidate station to itself (0.0)
AllCounts<-array(ds_info$monthcounts[ghcnstationsmap])	# all counts, later filtered to only keep those >= 30 years/360 months

# Read in the other files for later
mooD<-scan(paste(dirlist,infildists,sep=""),what="list",sep="\n")
mooCCl0<-scan(paste(dirlist,inCORRlag0,sep=""),what="list",sep="\n")
mooCCl1<-scan(paste(dirlist,inCORRlag1,sep=""),what="list",sep="\n")

for (stattoo in 1:nGstations){
  
  # is this station a GHCN station?
  if ((ds_info$statsource[ghcnstationsmap[stattoo]] == 1) & (ds_info$monthcounts[ghcnstationsmap[stattoo]] >= 360)) { # must have thirty years of data in candidate
    
    print(stattoo)
    
    tmpCClag1<-array(0,100)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)
    tmpCClag0<-array(0,100)	# temporary array for station+neighbours cross correlation (climate anomalies with loess)
    tmpNEIGHS<-array(0,100)	# temporary array for a list of the neighbours
    tmpHORIZDISTS<-array(0,100) # temporary array for a list of horizontal dists
    tmpVERTDISTS<-array(0,100)  # temporary array for a list of vertical dists
    tmpCOUNTS<-array(0,100)	# temporary array for a list of month counts

    # Find the list of distances for this station
    foundit<-which(ds_info$statid == ds_info$statid[ghcnstationsmap[stattoo]])
#    mooD<-scan(paste(dirlist,infildists,sep=""),what="list",skip=(foundit-1),nlines=1)
#    mooCCl0<-scan(paste(dirlist,inCORRlag0,sep=""),what="list",skip=(foundit-1),nlines=1)
#    mooCCl1<-scan(paste(dirlist,inCORRlag1,sep=""),what="list",skip=(foundit-1),nlines=1)

    # pull out 40 highest corrs and their ID
    tmpCClag0<-array(unlist(strsplit(mooCCl0[foundit],split=" ")))
    tmpCClag0<-tmpCClag0[tmpCClag0 != ""]
    tmpCClag0<-array(tmpCClag0[2:201],c(2,100)) # This should be r1=statid, r2=corr
    tmpNEIGHS<-tmpCClag0[1,]
    tmpCClag0<-tmpCClag0[2,]
    
    tmpCClag1<-array(unlist(strsplit(mooCCl1[foundit],split=" ")))
    tmpCClag1<-tmpCClag1[tmpCClag1 != ""]
    tmpCClag1<-array(tmpCClag1[2:201],c(2,100)) # This should be r1=statid, r2=corr
    tmpCClag1<-tmpCClag1[2,]

    tmpHORIZDISTS<-array(unlist(strsplit(mooD[foundit],split=" ")))
    tmpHORIZDISTS<-tmpHORIZDISTS[tmpHORIZDISTS != ""]
    tmpHORIZDISTS<-array(tmpHORIZDISTS[2:2001],c(2,1000)) # This should be r1=statid, r2=corr
    
    # search through distances and match up
    founddists<-apply(t(tmpNEIGHS),1,function(x) (1:1000)[tmpHORIZDISTS[1,] %in% x])
    tmpHORIZDISTS<-tmpHORIZDISTS[2,founddists]
    
    # calculate vertical distances in km
    foundneighs<-apply(t(tmpNEIGHS),1,function(x) (1:nstations)[ds_info$statid %in% x])
    tmpVERTDISTS<-abs((replicate(100,ds_info$statelvs[ghcnstationsmap[stattoo]])-ds_info$statelvs[foundneighs])/1000.)
    
    # get the counts of all of these neighbour stations
    tmpCOUNTS<-array(ds_info$monthcounts[foundneighs])
    
    # populate mother arrays
    AllCrossCorrs<-c(AllCrossCorrs,array(unlist(tmpCClag0)))
    AllCrossCorrsLag1<-c(AllCrossCorrsLag1,array(unlist(tmpCClag1)))
    AllDistances<-c(AllDistances,array(unlist(tmpHORIZDISTS)))
    AllVertDistances<-c(AllVertDistances,array(unlist(tmpVERTDISTS)))
    AllCounts<-c(AllCounts,array(unlist(tmpCOUNTS)))
#    stop()
 
  }
}  
 
# AllDistances<-type.convert(AllDistances)

noo<-order(type.convert(AllDistances))									      
OAD<-type.convert(AllDistances[noo])
# THE QUANTILES OCC>0.0, OCCl1>0.0, >=360 months
#	 1%	    2%         3%	  4%	     5% 	6%	   7% 
#   0.00000    9.88276   20.67612   27.43052   32.51300   36.78484   40.68066 
#	 8%	    9%        10%	 11%	    12%        13%	  14% 
#  44.41100   47.92784   51.19100   54.36436   57.46356   60.48700   63.31932 
#	15%	   16%        17%	 18%	    19%        20%	  21% 
#  66.07170   68.73600   71.30130   73.83984   76.40322   78.74660   81.05896 
#	22%	   23%        24%	 25%	    26%        27%	  28% 
#  83.47236   85.79474   88.05636   90.30200   92.54000   94.71200   96.92864 
#	29%	   30%        31%	 32%	    33%        34%	  35% 
#  99.12002  101.28180  103.42200  105.46416  107.62300  109.70892  111.74230 
#	36%	   37%        38%	 39%	    40%        41%	  42% 
# 113.82200  115.89600  117.92844  120.02282  122.11120  124.19158  126.25900 
#	43%	   44%        45%	 46%	    47%        48%	  49% 
# 128.36934  130.51000  132.75000  134.92248  137.10686  139.32248  141.53762 
#	50%	   51%        52%	 53%	    54%        55%	  56% 
# 143.85500  146.12300  148.39152  150.78728  153.13004  155.57300  158.03428 
#	57%	   58%        59%	 60%	    61%        62%	  63% 
# 160.60900  163.25704  165.96200  168.78360  171.65136  174.65212  177.82994 
#	64%	   65%        66%	 67%	    68%        69%	  70% 
# 180.99100  184.42600  187.98700  191.71030  195.77700  200.16966  204.83960 
#	71%	   72%        73%	 74%	    75%        76%	  77% 
# 210.15298  216.29780  223.50374  231.87572  241.75500  252.70364  265.61530 
#	78%	   79%        80%	 81%	    82%        83%	  84% 
# 280.04944  296.58210  315.10660  334.74524  355.38736  376.97972  399.64704 
#	85%	   86%        87%	 88%	    89%        90%	  91% 
# 423.09970  447.45876  473.16730  502.36244  534.43102  568.64200  605.39896 
# 	92%	   93%        94%	 95%	    96%        97%	  98% 
#  645.94484  693.82370  748.24460  811.13290  887.48400  993.93834 1152.30448 
# 	99%	  100% 
# 1405.27364 5226.19300 
# 29 % of neighbours are within 100km
# 97 % of 'neighbours' are within 1000km
# DO WE NEED TO SAMPLE OUTSIDE OF THIS RANGE OR CAN WE ASSUME CURVE CONTINUES AT SAME DECAY RATE?

OAVD<-AllVertDistances[noo]
# THE QUANTILES OCC>0.0, OCCl1>0.0, >=360 months
#	1%	  2%	    3%        4%	5%	  6%	    7%        8% 
#0.0000000 0.0000000 0.0010000 0.0020000 0.0030000 0.0040000 0.0050000 0.0061000 
#	9%	 10%	   11%       12%       13%	 14%	   15%       16% 
#0.0072000 0.0086000 0.0097000 0.0110000 0.0122000 0.0137000 0.0150000 0.0162000 
#      17%	 18%	   19%       20%       21%	 22%	   23%       24% 
#0.0180000 0.0190000 0.0210000 0.0220000 0.0240000 0.0253000 0.0270000 0.0290000 
#      25%	 26%	   27%       28%       29%	 30%	   31%       32% 
#0.0305000 0.0320000 0.0339000 0.0359000 0.0375000 0.0395400 0.0411000 0.0430000 
#      33%	 34%	   35%       36%       37%	 38%	   39%       40% 
#0.0451000 0.0472000 0.0491000 0.0518000 0.0540000 0.0560000 0.0580000 0.0609000 
#      41%	 42%	   43%       44%       45%	 46%	   47%       48% 
#0.0631000 0.0659000 0.0683000 0.0710000 0.0735000 0.0762000 0.0792000 0.0823000 
#      49%	 50%	   51%       52%       53%	 54%	   55%       56% 
#0.0851000 0.0884000 0.0914000 0.0945000 0.0980000 0.1015000 0.1051000 0.1091000 
#      57%	 58%	   59%       60%       61%	 62%	   63%       64% 
#0.1130000 0.1170000 0.1211000 0.1253000 0.1300000 0.1350000 0.1400000 0.1450000 
#      65%	 66%	   67%       68%       69%	 70%	   71%       72% 
#0.1503000 0.1560000 0.1620000 0.1682000 0.1750000 0.1820000 0.1900000 0.1980000 
#      73%	 74%	   75%       76%       77%	 78%	   79%       80% 
#0.2063000 0.2152000 0.2250000 0.2359880 0.2470000 0.2590000 0.2719000 0.2860000 
#      81%	 82%	   83%       84%       85%	 86%	   87%       88% 
#0.3010000 0.3170000 0.3342724 0.3532000 0.3734000 0.3956000 0.4204000 0.4480000 
#      89%	 90%	   91%       92%       93%	 94%	   95%       96% 
#0.4780000 0.5130000 0.5540000 0.6020000 0.6526340 0.7130000 0.7890000 0.8830000 
#      97%	 98%	   99%      100% 
#1.0100000 1.1991720 1.5823340 6.6042000 
# 53% within 100m (0.1km)
# 96% within 1km
# 6.6 km seems a little far fetched

OCC<-type.convert(AllCrossCorrs[noo])
#qporbs<-seq(100)/100.
#qtls<-quantile(OCC,qporbs) # gives the % of points less than each percentile threshold
# THE QUANTILES OCC>0.0, OCCl1>0.0, >=360 months
#   1%    2%	3%    4%    5%    6%	7%    8%    9%   10%   11%   12%   13% 
#0.190 0.267 0.329 0.384 0.437 0.488 0.536 0.579 0.616 0.646 0.671 0.693 0.712 
#  14%   15%   16%   17%   18%   19%   20%   21%   22%   23%   24%   25%   26% 
#0.729 0.744 0.758 0.770 0.780 0.790 0.798 0.807 0.814 0.820 0.827 0.832 0.837 
#  27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39% 
#0.842 0.847 0.852 0.856 0.860 0.864 0.867 0.871 0.874 0.877 0.880 0.883 0.886 
#  40%   41%   42%   43%   44%   45%   46%   47%   48%   49%   50%   51%   52% 
#0.888 0.891 0.893 0.896 0.898 0.901 0.903 0.905 0.907 0.909 0.910 0.912 0.914 
#  53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65% 
#0.916 0.917 0.919 0.920 0.922 0.923 0.924 0.926 0.927 0.928 0.930 0.931 0.932 
#  66%   67%   68%   69%   70%   71%   72%   73%   74%   75%   76%   77%   78% 
#0.933 0.934 0.935 0.937 0.938 0.939 0.940 0.941 0.942 0.943 0.944 0.945 0.946 
#  79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
#0.947 0.948 0.949 0.950 0.951 0.952 0.953 0.954 0.955 0.957 0.958 0.959 0.961 
#  92%   93%   94%   95%   96%   97%   98%   99%  100% 
#0.962 0.964 0.966 0.968 0.970 0.974 0.980 1.000 1.000 
# mean = 0.8479, sd = 0.1658 
# 92% correlated > 0.6, 80% > 0.8, 56% > 0.9

OCCl1<-type.convert(AllCrossCorrsLag1[noo])
# THE QUANTILES OCC>0.0, OCCl1>0.0, >=360 months
#   1%    2%	3%    4%    5%    6%	7%    8%    9%   10%   11%   12%   13% 
#0.010 0.018 0.025 0.031 0.037 0.042 0.046 0.051 0.054 0.058 0.062 0.065 0.068 
#  14%   15%   16%   17%   18%   19%   20%   21%   22%   23%   24%   25%   26% 
#0.071 0.074 0.076 0.079 0.081 0.084 0.086 0.088 0.090 0.092 0.094 0.096 0.098 
#  27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39% 
#0.100 0.101 0.103 0.105 0.107 0.108 0.110 0.112 0.113 0.115 0.117 0.118 0.120 
#  40%   41%   42%   43%   44%   45%   46%   47%   48%   49%   50%   51%   52% 
#0.121 0.123 0.124 0.126 0.127 0.129 0.130 0.132 0.133 0.135 0.136 0.138 0.140 
#  53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65% 
#0.141 0.143 0.144 0.146 0.147 0.149 0.151 0.152 0.154 0.156 0.157 0.159 0.161 
#  66%   67%   68%   69%   70%   71%   72%   73%   74%   75%   76%   77%   78% 
#0.162 0.164 0.166 0.168 0.170 0.172 0.174 0.176 0.178 0.180 0.182 0.185 0.187 
#  79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
#0.190 0.192 0.195 0.198 0.201 0.204 0.207 0.210 0.214 0.218 0.222 0.227 0.232 
#  92%   93%   94%   95%   96%   97%   98%   99%  100% 
#0.238 0.244 0.252 0.261 0.271 0.285 0.303 0.334 0.903 
# mean = 0.1415, sd = 0.0689
# 26% < 0.1, 82% < 0.2, 97% < 0.3, 99% <0.4
# 73% > 0.1, 18% > 0.2, 2% > 0.3, 1% > 0.4 

# first nGstations should be lag0 and lag1 autocorrelation - for plotting later

# remove really bad corrs and autocorrs to try and improve the fit (start with 370189 points)
# could exclude very large autocorrelation but there are only a few > 0.8 (37) and in a 
# very stable location (e.g., tropics) AC could be high.
print(paste("No. of points: ",length(OAD),sep=""))
# 442529 (could have been 9029*100 if no exclusion from monthcounts)

# reduce to only those points where the station has > 30 years (360 months) of data
AllCounts<-AllCounts[noo]
OAD<-OAD[which(AllCounts >= 360)]
OAVD<-OAVD[which(AllCounts >= 360)]
OCCl1<-OCCl1[which(AllCounts >= 360)]
OCC<-OCC[which(AllCounts >= 360)]
print(paste("No. of points after removal of counts<360: ",length(OAD),sep=""))
# 303043

#OAD<-OAD[which(OCC > 0.1)]
#OAVD<-OAVD[which(OCC > 0.1)]
#OCCl1<-OCCl1[which(OCC > 0.1)]
#OCC<-OCC[which(OCC > 0.1)]
#print(paste("No. of points after removal of OCC<0.1: ",length(OAD),sep=""))
## 367818
OAD<-OAD[which(OCC > 0.0)]
OAVD<-OAVD[which(OCC > 0.0)]
OCCl1<-OCCl1[which(OCC > 0.0)]
OCC<-OCC[which(OCC > 0.0)]
print(paste("No. of points after removal of OCC<0.0: ",length(OAD),sep=""))
# 301283

OAD<-OAD[which(OCCl1 > 0.0)]
OAVD<-OAVD[which(OCCl1 > 0.0)]
OCC<-OCC[which(OCCl1 > 0.0)]
OCCl1<-OCCl1[which(OCCl1 > 0.0)]
print(paste("No. of points after removal of OCCl1>0.0: ",length(OAD),sep=""))
# 334932 (OCC>0.1), 284939 (OCC>0.0)

print(paste("OCC mean and sd:",mean(OCC), sd(OCC),sep=" "))
print(paste("OCCl1 mean and sd:",mean(OCCl1), sd(OCCl1),sep=" "))

# COMPUTE RELATIONSHIP OF CC and CClag1 with DISTANCE
datccl0<-data.frame(OAD,OAVD,OCC)   
datccl1<-data.frame(OAD,OAVD,OCCl1)   

# DISTANCE AND ELEVATION
# allow all elements to be fitted
f3<-function(x,y,a,b,c){ a * exp((-b*x) + (-c*y))}
#plot(y~x)	# y as a function of x
fm<-nls(OCC~f3(OAD,OAVD,a,b,c),data=datccl0,start=c(a=0.999,b=0.001,c=0.01))
co3ccl0<-coef(fm)
print(paste("lag 0 3 element: ",co3ccl0))
#co3ccl0[1]=1.008 co3ccl0[2]=0.0007 co3ccl0[3]=0.1

fm<-nls(OCCl1~f3(OAD,OAVD,a,b,c),data=datccl1,start=c(a=0.499,b=0.001,c=0.01)) # makes v little diff is a=0.999 < 0.000001
co3ccl1<-coef(fm)
print(paste("lag 1 3 element: ",co3ccl1))
#co3ccl1[1]=0.146 co3ccl1[2]=-0.0001, co3ccl1[3]=0.1

# restrict fit to decay elements only
f2ccl0<-function(x,y,b,c){ 0.97 * exp((-b*x) + (-c*y))}
fm<-nls(OCC~f2ccl0(OAD,OAVD,b,c),data=datccl0,start=c(b=0.001,c=0.01))
co2ccl0<-coef(fm)
print(paste("lag 0 2 element: ",co2ccl0))
#co2ccl0[1]=0.0006 co2ccl0[2]=0.07 

f2ccl1<-function(x,y,b,c){ 0.23 * exp((-b*x) + (-c*y))}	# DIAGS 0.23
fm<-nls(OCCl1~f2ccl1(OAD,OAVD,b,c),data=datccl1,start=c(b=0.001,c=0.01))
co2ccl1<-coef(fm)
print(paste("lag 1 2 element: ",co2ccl1))
#co2ccl1[1]=0.001 co2ccl1[2]=0.82

# BEST FIT PLOTS
library(scatterplot3d)

setEPS()
# 3 element fit CCl0 - fit to real data
postscript(outplotFit3CCl0)
# object with lines
plotobj<-scatterplot3d(OAD,OAVD,OCC,pch='.',highlight.3d=TRUE,
         xlim=c(0,ceiling(max(OAD)/1000.)*1000),
	 ylim=c(0,ceiling(max(OAVD))),
	 zlim=c(0,1),
	 xlab='Horizontal Distance (km)',
	 ylab='Vertical Distance (km)',
	 zlab='Correlation at lag 0')
# function curve
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,ceiling(max(OAVD)),length.out=100),
        f3(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,ceiling(max(OAVD)),length.out=100),
	   a=co3ccl0[1],b=co3ccl0[2],c=co3ccl0[3]),
	   col='blue',type='l',lwd=5)
#function curve if no height difference
plotobj$points3d(seq(0,0,length.out=100),seq(0,ceiling(max(OAVD)),length.out=100),
        f3(seq(0,0,length.out=100),seq(0,ceiling(max(OAVD)),length.out=100),
	   a=co3ccl0[1],b=co3ccl0[2],c=co3ccl0[3]),
	   col='cyan3',type='l',lwd=5)
#function curve if no horizontal distance away
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
        f3(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
	   a=co3ccl0[1],b=co3ccl0[2],c=co3ccl0[3]),
	   col='cornflowerblue',type='l',lwd=5)
# show function as an equation on the plot
#text(0.15,0.1,paste("R = ",round(co3ccl0[1],3)," * exp((-",
#              round(co3ccl0[2],4),"* dh) + (-",
#	      round(co3ccl0[3],2),"* dv))"),pos=4,cex=1.2)
mtext(paste("R = ",round(co3ccl0[1],3)," * exp((-",round(co3ccl0[2],4),"* dh) + (-",round(co3ccl0[3],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()


setEPS()
# 2 element fit CCl0 - fit to real data
postscript(outplotFit2CCl0)
# object with lines
plotobj<-scatterplot3d(OAD,OAVD,OCC,pch='.',highlight.3d=TRUE,
         xlim=c(0,ceiling(max(OAD)/1000.)*1000),
	 ylim=c(0,ceiling(max(OAVD))),
	 zlim=c(0,1),
	 xlab='Horizontal Distance (km)',
	 ylab='Vertical Distance (km)',
	 zlab='Correlation at lag 0')
# function curve
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f2ccl0(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   b=co2ccl0[1],c=co2ccl0[2]),
	   col='blue',type='l',lwd=5)
#function curve if no height difference
plotobj$points3d(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f2ccl0(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   b=co2ccl0[1],c=co2ccl0[2]),
	   col='cyan3',type='l',lwd=5)
#function curve if no horizontal distance away
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
        f2ccl0(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
	   b=co2ccl0[1],c=co2ccl0[2]),
	   col='cornflowerblue',type='l',lwd=5)
mtext(paste("R = 0.97 * exp((-",round(co2ccl0[1],4),"* dh) + (-",round(co2ccl0[2],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()
 

setEPS()
# 3 element fit CCl1 - fit to real data
postscript(outplotFit3CCl1)
# object with lines
plotobj<-scatterplot3d(OAD,OAVD,OCCl1,pch='.',highlight.3d=TRUE,
         xlim=c(0,ceiling(max(OAD)/1000.)*1000),
	 ylim=c(0,ceiling(max(OAVD))),
	 zlim=c(0,1),
	 xlab='Horizontal Distance (km)',
	 ylab='Vertical Distance (km)',
	 zlab='Correlation at lag 1')
# function curve
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f3(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   a=co3ccl1[1],b=co3ccl1[2],c=co3ccl1[3]),
	   col='blue',type='l',lwd=5)
#function curve if no height difference
plotobj$points3d(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f3(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   a=co3ccl1[1],b=co3ccl1[2],c=co3ccl1[3]),
	   col='cyan3',type='l',lwd=5)
#function curve if no horizontal distance away
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
        f3(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
	   a=co3ccl1[1],b=co3ccl1[2],c=co3ccl1[3]),
	   col='cornflowerblue',type='l',lwd=5)
mtext(paste("R = ",round(co3ccl1[1],3)," * exp((-",round(co3ccl1[2],4),"* dh) + (-",round(co3ccl1[3],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()


setEPS()
# 2 element fit CCl1 - fit to real data
postscript(outplotFit2CCl1)
# object with lines
plotobj<-scatterplot3d(OAD,OAVD,OCCl1,pch='.',highlight.3d=TRUE,
         xlim=c(0,ceiling(max(OAD)/1000.)*1000),
	 ylim=c(0,ceiling(max(OAVD))),
	 zlim=c(0,1),
	 xlab='Horizontal Distance (km)',
	 ylab='Vertical Distance (km)',
	 zlab='Correlation at lag 1')
# function curve
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f2ccl1(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   b=co2ccl1[1],c=co2ccl1[2]),
	   col='blue',type='l',lwd=5)
#function curve if no height difference
plotobj$points3d(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f2ccl1(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   b=co2ccl1[1],c=co2ccl1[2]),
	   col='cyan3',type='l',lwd=5)
#function curve if no horizontal distance away
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
        f2ccl1(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
	   b=co2ccl1[1],c=co2ccl1[2]),
	   col='cornflowerblue',type='l',lwd=5)
mtext(paste("R = 0.23 * exp((-",round(co2ccl1[1],4),"* dh) + (-",round(co2ccl1[2],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()

setEPS()
# 3 element fit CCl1 using ccl0 - fit to real data
postscript(outplotFit3CCl1l0)
# object with lines
plotobj<-scatterplot3d(OAD,OAVD,OCCl1,pch='.',highlight.3d=TRUE,
         xlim=c(0,ceiling(max(OAD)/1000.)*1000),
	 ylim=c(0,ceiling(max(OAVD))),
	 zlim=c(0,1),
	 xlab='Horizontal Distance (km)',
	 ylab='Vertical Distance (km)',
	 zlab='Correlation at lag 1')
# function curve
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f3(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   a=0.3,b=co3ccl0[2],c=co3ccl0[3]),
	   col='blue',type='l',lwd=5)
#function curve if no height difference
plotobj$points3d(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f3(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   a=0.3,b=co3ccl0[2],c=co3ccl0[3]),
	   col='cyan3',type='l',lwd=5)
#function curve if no horizontal distance away
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
        f3(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
	   a=0.3,b=co3ccl0[2],c=co3ccl0[3]),
	   col='cornflowerblue',type='l',lwd=5)
mtext(paste("R = 0.23 * exp((-",round(co3ccl0[2],4),"* dh) + (-",round(co3ccl0[3],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()

setEPS()
# 2 element fit CCl1 using ccl0 - fit to real data
postscript(outplotFit2CCl1l0)
# object with lines
plotobj<-scatterplot3d(OAD,OAVD,OCCl1,pch='.',highlight.3d=TRUE,
         xlim=c(0,ceiling(max(OAD)/1000.)*1000),
	 ylim=c(0,ceiling(max(OAVD))),
	 zlim=c(0,1),
	 xlab='Horizontal Distance (km)',
	 ylab='Vertical Distance (km)',
	 zlab='Correlation at lag 1')
# function curve
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f2ccl1(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   b=co2ccl0[1],c=co2ccl0[2]),
	   col='blue',type='l',lwd=5)
#function curve if no height difference
plotobj$points3d(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
        f2ccl1(seq(0,0,length.out=100),seq(0,(ceiling(max(OAVD))/1000.)*1000,length.out=100),
	   b=co2ccl0[1],c=co2ccl0[2]),
	   col='cyan3',type='l',lwd=5)
#function curve if no horizontal distance away
plotobj$points3d(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
        f2ccl1(seq(0,ceiling(max(OAD)/1000.)*1000,length.out=100),seq(0,0,length.out=100),
	   b=co2ccl0[1],c=co2ccl0[2]),
	   col='cornflowerblue',type='l',lwd=5)
mtext(paste("R = 0.23 * exp((-",round(co2ccl0[1],4),"* dh) + (-",round(co2ccl0[2],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()


# SCATTER PLOTS - NEED TO SIMULATE


# 3 element fit CCl0 - scatter real vs estimated data
setEPS()
postscript(outplotScat3CCl0)
simsies<-f3(OAD,OAVD,co3ccl0[1],co3ccl0[2],co3ccl0[3])
plot(OCC,simsies,pch='.',
         xlim=c(0,1),ylim=c(0,1),
	 xlab='GHCN Cross Correlations at lag 0',
	 ylab='Estimated Cross Correlations at lag 0')
# plot 1:1 line
lines(seq(0,1,0.1),seq(0,1,0.1),col='red',lw=5)
# show self correlation in red 
# (shows that we need to force diagonals - but that we could in theory cope with same locs)
points(OCC[which(OAD == 0.0)],simsies[which(OAD == 0.0)],col='red',pch=20)
# show function as an equation on the plot
#text(0.15,0.1,paste("R = ",round(co3ccl0[1],3)," * exp((-",
#              round(co3ccl0[2],4),"* dh) + (-",
#	      round(co3ccl0[3],2),"* dv))"),pos=4,cex=1.2)
mtext(paste("R = ",round(co3ccl0[1],3)," * exp((-",round(co3ccl0[2],4),"* dh) + (-",round(co3ccl0[3],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()
qporbs<-seq(100)/100.
print("QUANTILES: Simulations for 3CCl0")
qtls<-quantile(simsies,qporbs)
print(qtls)
print(paste("Simulated OCC mean and sd:",mean(simsies), sd(simsies),sep=" "))
#[1] "QUANTILES: Simulations for 3CCl0"
#	1%	  2%	    3%        4%	5%	  6%	    7%        8% 
#0.3566581 0.4250671 0.4754246 0.5134567 0.5436813 0.5687574 0.5917287 0.6119314 
#	9%	 10%	   11%       12%       13%	 14%	   15%       16% 
#0.6305278 0.6471922 0.6629742 0.6786301 0.6939629 0.7078895 0.7212883 0.7339511 
#      17%	 18%	   19%       20%       21%	 22%	   23%       24% 
#0.7458028 0.7576031 0.7687444 0.7796165 0.7899080 0.7992802 0.8083314 0.8164456 
#      25%	 26%	   27%       28%       29%	 30%	   31%       32% 
#0.8238851 0.8303563 0.8360698 0.8410313 0.8456402 0.8497762 0.8535718 0.8569831 
#      33%	 34%	   35%       36%       37%	 38%	   39%       40% 
#0.8601062 0.8630233 0.8656806 0.8681598 0.8704772 0.8727018 0.8748612 0.8768945 
#      41%	 42%	   43%       44%       45%	 46%	   47%       48% 
#0.8788624 0.8807404 0.8825618 0.8843821 0.8861753 0.8878563 0.8894880 0.8911011 
#      49%	 50%	   51%       52%       53%	 54%	   55%       56% 
#0.8926974 0.8942842 0.8958643 0.8973796 0.8988874 0.9003700 0.9018674 0.9033750 
#      57%	 58%	   59%       60%       61%	 62%	   63%       64% 
#0.9048545 0.9063780 0.9077840 0.9092016 0.9106781 0.9120891 0.9135634 0.9150246 
#      65%	 66%	   67%       68%       69%	 70%	   71%       72% 
#0.9165041 0.9180096 0.9194828 0.9209593 0.9224461 0.9239613 0.9255307 0.9270323 
#      73%	 74%	   75%       76%       77%	 78%	   79%       80% 
#0.9285530 0.9301372 0.9316963 0.9333227 0.9349589 0.9366725 0.9383874 0.9400926 
#      81%	 82%	   83%       84%       85%	 86%	   87%       88% 
#0.9418373 0.9435718 0.9454356 0.9473183 0.9493157 0.9513086 0.9533877 0.9555781 
#      89%	 90%	   91%       92%       93%	 94%	   95%       96% 
#0.9577931 0.9602368 0.9627850 0.9655549 0.9684655 0.9717302 0.9754175 0.9797602 
#      97%	 98%	   99%      100% 
#0.9849499 0.9937186 1.0083089 1.0083089 
# 93% > 0.6, 78% > 0.8, 47% > 0.9
#[1] "Simulated OCC mean and sd: 0.848466093060066 0.136506543109533"
# ACTUAL DATA:
# mean = 0.8479, sd = 0.1658 
# 92% correlated > 0.6, 80% > 0.8, 56% > 0.9

# 2 element fit CCl0 - scatter real vs estimated data
setEPS()
postscript(outplotScat2CCl0)
simsies<-f2ccl0(OAD,OAVD,co2ccl0[1],co2ccl0[2])
plot(OCC,simsies,pch='.',
         xlim=c(0,1),ylim=c(0,1),
	 xlab='GHCN Cross Correlations at lag 0',
	 ylab='Estimated Cross Correlations at lag 0')
# plot 1:1 line
lines(seq(0,1,0.1),seq(0,1,0.1),col='red',lw=5)
# show self correlation in red 
# (shows that we need to force diagonals - but that we could in theory cope with same locs)
points(OCC[which(OAD == 0.0)],simsies[which(OAD == 0.0)],col='red',pch=20)
# show function as an equation on the plot
#text(0.15,0.0,paste("R = 0.97 * exp((-",
#             round(co2ccl0[1],4),"* dh) + (-",
#	      round(co2ccl0[2],2),"* dv))"),pos=4,cex=1.2)
mtext(paste("R = 0.97 * exp((-",round(co2ccl0[1],4),"* dh) + (-",round(co2ccl0[2],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()
print("QUANTILES: Simulations for 2CCl0")
qtls<-quantile(simsies,qporbs)
print(qtls)
print(paste("Simulated OCC mean and sd:",mean(simsies), sd(simsies),sep=" "))
#[1] "QUANTILES: Simulations for 2CCl0"
#	1%	  2%	    3%        4%	5%	  6%	    7%        8% 
#0.3925084 0.4587840 0.5059150 0.5417550 0.5688781 0.5917272 0.6125328 0.6311512 
#	9%	 10%	   11%       12%       13%	 14%	   15%       16% 
#0.6477000 0.6626707 0.6767643 0.6903715 0.7040465 0.7163057 0.7280244 0.7390100 
#      17%	 18%	   19%       20%       21%	 22%	   23%       24% 
#0.7495419 0.7598586 0.7696900 0.7793682 0.7885986 0.7972049 0.8048245 0.8118882 
#      25%	 26%	   27%       28%       29%	 30%	   31%       32% 
#0.8180711 0.8235856 0.8282302 0.8326145 0.8364451 0.8397785 0.8427787 0.8454994 
#      33%	 34%	   35%       36%       37%	 38%	   39%       40% 
#0.8480027 0.8503332 0.8524541 0.8544361 0.8563266 0.8581644 0.8599102 0.8615990 
#      41%	 42%	   43%       44%       45%	 46%	   47%       48% 
#0.8631877 0.8647730 0.8662815 0.8677821 0.8692151 0.8706433 0.8719899 0.8733229 
#      49%	 50%	   51%       52%       53%	 54%	   55%       56% 
#0.8746304 0.8759550 0.8772779 0.8785491 0.8798314 0.8810695 0.8823122 0.8835323 
#      57%	 58%	   59%       60%       61%	 62%	   63%       64% 
#0.8847563 0.8860263 0.8872390 0.8884336 0.8896362 0.8908805 0.8920732 0.8932884 
#      65%	 66%	   67%       68%       69%	 70%	   71%       72% 
#0.8944956 0.8957237 0.8969920 0.8982176 0.8994533 0.9007432 0.9019832 0.9032298 
#      73%	 74%	   75%       76%       77%	 78%	   79%       80% 
#0.9045446 0.9058339 0.9071800 0.9085394 0.9099754 0.9113698 0.9128182 0.9142068 
#      81%	 82%	   83%       84%       85%	 86%	   87%       88% 
#0.9156240 0.9171294 0.9186693 0.9202641 0.9219104 0.9236205 0.9253705 0.9271257 
#      89%	 90%	   91%       92%       93%	 94%	   95%       96% 
#0.9290084 0.9309767 0.9330974 0.9353765 0.9377786 0.9404440 0.9434677 0.9469851 
#      97%	 98%	   99%      100% 
#0.9513493 0.9585115 0.9700000 0.9700000 
# 94% > 0.6, 78% > 0.8, 30% > 0.9
#[1] "Simulated OCC mean and sd: 0.83560643445337 0.119733928168211"
# ACTUAL DATA:
# mean = 0.8479, sd = 0.1658 
# 92% correlated > 0.6, 80% > 0.8, 56% > 0.9

# 3 element fit CCl1 - scatter real vs estimated data
setEPS()
postscript(outplotScat3CCl1)
simsies<-f3(OAD,OAVD,co3ccl1[1],co3ccl1[2],co3ccl1[3])
plot(OCCl1,simsies,pch='.',
         xlim=c(0,1),ylim=c(0,1),
	 xlab='GHCN Cross Correlations at lag 1',
	 ylab='Estimated Cross Correlations at lag 1')
# plot 1:1 line
lines(seq(0,1,0.1),seq(0,1,0.1),col='red',lw=5)
# show self correlation in red 
# (shows that we need to force diagonals - but that we could in theory cope with same locs)
points(OCCl1[which(OAD == 0.0)],simsies[which(OAD == 0.0)],col='red',pch=20)
# show function as an equation on the plot
#text(0.15,0.1,paste("R = ",round(co3ccl1[1],3)," * exp((-",
#              round(co3ccl1[2],4),"* dh) + (-",
#	      round(co3ccl1[3],2),"* dv))"),pos=4,cex=1.2)
mtext(paste("R = ",round(co3ccl1[1],3)," * exp((-",round(co3ccl1[2],4),"* dh) + (-",round(co3ccl1[3],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()
print("QUANTILES: Simulations for 3CCl1")
qtls<-quantile(simsies,qporbs)
print(qtls)
print(paste("Simulated OCCl1 mean and sd:",mean(simsies), sd(simsies),sep=" "))
#[1] "QUANTILES: Simulations for 3CCl1"
#	1%	  2%	    3%        4%	5%	  6%	    7%        8% 
#0.1216912 0.1265639 0.1291397 0.1307771 0.1320407 0.1330249 0.1338840 0.1345919 
#	9%	 10%	   11%       12%       13%	 14%	   15%       16% 
#0.1352259 0.1357757 0.1362856 0.1367468 0.1371645 0.1375449 0.1378849 0.1382119 
#      17%	 18%	   19%       20%       21%	 22%	   23%       24% 
#0.1385158 0.1387926 0.1390588 0.1393116 0.1395471 0.1397756 0.1399880 0.1401869 
#      25%	 26%	   27%       28%       29%	 30%	   31%       32% 
#0.1403774 0.1405580 0.1407285 0.1408955 0.1410545 0.1412005 0.1413436 0.1414847 
#      33%	 34%	   35%       36%       37%	 38%	   39%       40% 
#0.1416179 0.1417450 0.1418654 0.1419778 0.1420910 0.1422018 0.1423051 0.1424072 
#      41%	 42%	   43%       44%       45%	 46%	   47%       48% 
#0.1425026 0.1425945 0.1426817 0.1427659 0.1428498 0.1429324 0.1430114 0.1430873 
#      49%	 50%	   51%       52%       53%	 54%	   55%       56% 
#0.1431560 0.1432252 0.1432925 0.1433607 0.1434257 0.1434856 0.1435447 0.1436030 
#      57%	 58%	   59%       60%       61%	 62%	   63%       64% 
#0.1436602 0.1437125 0.1437642 0.1438140 0.1438647 0.1439131 0.1439610 0.1440072 
#      65%	 66%	   67%       68%       69%	 70%	   71%       72% 
#0.1440507 0.1440942 0.1441362 0.1441788 0.1442190 0.1442566 0.1442953 0.1443346 
#      73%	 74%	   75%       76%       77%	 78%	   79%       80% 
#0.1443725 0.1444106 0.1444474 0.1444847 0.1445220 0.1445587 0.1445970 0.1446347 
#      81%	 82%	   83%       84%       85%	 86%	   87%       88% 
#0.1446719 0.1447098 0.1447489 0.1447871 0.1448258 0.1448662 0.1449082 0.1449500 
#      89%	 90%	   91%       92%       93%	 94%	   95%       96% 
#0.1449955 0.1450423 0.1450902 0.1451420 0.1451986 0.1452628 0.1453351 0.1454179 
#      97%	 98%	   99%      100% 
#0.1455248 0.1456927 0.1460153 0.1460153 
# 0% < 0.1, 100% < 0.2
# 100% > 0.1, 0% > 0.2 
#[1] "Simulated OCCl1 mean and sd: 0.141493294810322 0.0049505327365448"
# ACTUAL DATA:
# mean = 0.1415, sd = 0.0689 
# 26% < 0.1, 82% < 0.2, 97% < 0.3, 99% <0.4
# 73% > 0.1, 18% > 0.2, 2% > 0.3, 1% > 0.4 

# 2 element fit CCl1 - scatter real vs estimated data
setEPS()
postscript(outplotScat2CCl1)
simsies<-f2ccl1(OAD,OAVD,co2ccl1[1],co2ccl1[2])
plot(OCCl1,simsies,pch='.',
         xlim=c(0,1),ylim=c(0,1),
	 xlab='GHCN Cross Correlations at lag 1',
	 ylab='Estimated Cross Correlations at lag 1')
# plot 1:1 line
lines(seq(0,1,0.1),seq(0,1,0.1),col='red',lw=5)
# show self correlation in red 
# (shows that we need to force diagonals - but that we could in theory cope with same locs)
points(OCCl1[which(OAD == 0.0)],simsies[which(OAD == 0.0)],col='red',pch=20)
# show function as an equation on the plot
#text(0.15,0.1,paste("R = 0.2 * exp((-",
#              round(co2ccl1[1],4),"* dh) + (-",
#	      round(co2ccl1[2],2),"* dv))"),pos=4,cex=1.2)
mtext(paste("R = 0.23 * exp((-",round(co2ccl1[1],4),"* dh) + (-",round(co2ccl1[2],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()
print("QUANTILES: Simulations for 2CCl1")
qtls<-quantile(simsies,qporbs)
print(qtls)
print(paste("Simulated OCCl1 mean and sd:",mean(simsies), sd(simsies),sep=" "))
#[1] "QUANTILES: Simulations for 2CCl1"
#	 1%	    2%         3%	  4%	     5% 	6%	   7% 
#0.02344470 0.03393317 0.04214756 0.04917638 0.05547151 0.06141558 0.06675826 
#	 8%	    9%        10%	 11%	    12%        13%	  14% 
#0.07182202 0.07642571 0.08077857 0.08505576 0.08900604 0.09294162 0.09659322 
#	15%	   16%        17%	 18%	    19%        20%	  21% 
#0.10026035 0.10379575 0.10727371 0.11077174 0.11413899 0.11727709 0.12025206 
#	22%	   23%        24%	 25%	    26%        27%	  28% 
#0.12314063 0.12589269 0.12853499 0.13114150 0.13368293 0.13610189 0.13852635 
#	29%	   30%        31%	 32%	    33%        34%	  35% 
#0.14082694 0.14299557 0.14510606 0.14712705 0.14908684 0.15103126 0.15284205 
#	36%	   37%        38%	 39%	    40%        41%	  42% 
#0.15449655 0.15616083 0.15769168 0.15922492 0.16063935 0.16196938 0.16322541 
#	43%	   44%        45%	 46%	    47%        48%	  49% 
#0.16445520 0.16564566 0.16676424 0.16784684 0.16889462 0.16990027 0.17087658 
#	50%	   51%        52%	 53%	    54%        55%	  56% 
#0.17182242 0.17270299 0.17359270 0.17440844 0.17524157 0.17603725 0.17682966 
#	57%	   58%        59%	 60%	    61%        62%	  63% 
#0.17757997 0.17834784 0.17909370 0.17985103 0.18060100 0.18134994 0.18206550 
#	64%	   65%        66%	 67%	    68%        69%	  70% 
#0.18277121 0.18346712 0.18419808 0.18489974 0.18561175 0.18631764 0.18703666 
#	71%	   72%        73%	 74%	    75%        76%	  77% 
#0.18773272 0.18846139 0.18917459 0.18988028 0.19063377 0.19139277 0.19215479 
#	78%	   79%        80%	 81%	    82%        83%	  84% 
#0.19290834 0.19370165 0.19450340 0.19530149 0.19615867 0.19706072 0.19794302 
#	85%	   86%        87%	 88%	    89%        90%	  91% 
#0.19885326 0.19978482 0.20077713 0.20180245 0.20290172 0.20403346 0.20529937 
#	92%	   93%        94%	 95%	    96%        97%	  98% 
#0.20662181 0.20808818 0.20972405 0.21161559 0.21378899 0.21659612 0.22106958 
#	99%	  100% 
#0.23000000 0.23000000 
# 14% < 0.1, 86% < 0.2, 100% < 0.3
# 85% > 0.1, 13% >= 0.2
#[1] "Simulated OCCl1 mean and sd: 0.156537888390555 0.0478468728777009"
# ACTUAL DATA:
# mean = 0.1415, sd = 0.0689 
# 26% < 0.1, 82% < 0.2, 97% < 0.3, 99% <0.4
# 73% > 0.1, 18% > 0.2, 2% > 0.3, 1% > 0.4 
# THIS IS THE BEST!!!!

# 3 element fit CCl1 using ccl0 - scatter real vs estimated data
setEPS()
postscript(outplotScat3CCl1l0)
simsies<-f3(OAD,OAVD,0.23,co3ccl0[2],co3ccl0[3])
plot(OCCl1,simsies,pch='.',
         xlim=c(0,1),ylim=c(0,1),
	 xlab='GHCN Cross Correlations at lag 1',
	 ylab='Estimated Cross Correlations at lag 1')
# plot 1:1 line
lines(seq(0,1,0.1),seq(0,1,0.1),col='red',lw=5)
# show self correlation in red 
# (shows that we need to force diagonals - but that we could in theory cope with same locs)
points(OCCl1[which(OAD == 0.0)],simsies[which(OAD == 0.0)],col='red',pch=20)
# show function as an equation on the plot
#text(0.15,0.1,paste("R = 0.3 * exp((-",
#              round(co3ccl0[2],4),"* dh) + (-",
#	      round(co3ccl0[3],2),"* dv))"),pos=4,cex=1.2)
mtext(paste("R = 0.23 * exp((-",round(co3ccl0[2],4),"* dh) + (-",round(co3ccl0[3],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()
print("QUANTILES: Simulations for 3CCl1 use l0")
qtls<-quantile(simsies,qporbs)
print(qtls)
print(paste("Simulated OCCl1 mean and sd:",mean(simsies), sd(simsies),sep=" "))
#[1] "QUANTILES: Simulations for 3CCl1 use l0"
#	 1%	    2%         3%	  4%	     5% 	6%	   7% 
#0.08135539 0.09695980 0.10844657 0.11712189 0.12401625 0.12973623 0.13497610 
#	 8%	    9%        10%	 11%	    12%        13%	  14% 
#0.13958442 0.14382634 0.14762757 0.15122752 0.15479871 0.15829618 0.16147291 
#	15%	   16%        17%	 18%	    19%        20%	  21% 
#0.16452925 0.16741769 0.17012113 0.17281282 0.17535421 0.17783418 0.18018173 
#	22%	   23%        24%	 25%	    26%        27%	  28% 
#0.18231956 0.18438418 0.18623506 0.18793205 0.18940816 0.19071145 0.19184318 
#	29%	   30%        31%	 32%	    33%        34%	  35% 
#0.19289451 0.19383795 0.19470373 0.19548186 0.19619426 0.19685966 0.19746582 
#	36%	   37%        38%	 39%	    40%        41%	  42% 
#0.19803133 0.19855993 0.19906738 0.19955995 0.20002375 0.20047264 0.20090101 
#	43%	   44%        45%	 46%	    47%        48%	  49% 
#0.20131648 0.20173170 0.20214075 0.20252418 0.20289639 0.20326434 0.20362847 
#	50%	   51%        52%	 53%	    54%        55%	  56% 
#0.20399043 0.20435086 0.20469650 0.20504044 0.20537862 0.20572018 0.20606408 
#	57%	   58%        59%	 60%	    61%        62%	  63% 
#0.20640155 0.20674908 0.20706980 0.20739316 0.20772994 0.20805181 0.20838810 
#	64%	   65%        66%	 67%	    68%        69%	  70% 
#0.20872140 0.20905889 0.20940231 0.20973833 0.21007513 0.21041427 0.21075992 
#	71%	   72%        73%	 74%	    75%        76%	  77% 
#0.21111790 0.21146043 0.21180730 0.21216865 0.21252431 0.21289530 0.21326851 
#	78%	   79%        80%	 81%	    82%        83%	  84% 
#0.21365940 0.21405057 0.21443953 0.21483752 0.21523315 0.21565829 0.21608775 
#	85%	   86%        87%	 88%	    89%        90%	  91% 
#0.21654337 0.21699796 0.21747220 0.21797184 0.21847710 0.21903451 0.21961577 
#	92%	   93%        94%	 95%	    96%        97%	  98% 
#0.22024761 0.22091152 0.22165622 0.22249730 0.22348789 0.22467169 0.22667188 
#	99%	  100% 
#0.23000000 0.23000000 
# 2% < 0.1, 39% < 0.2, 100% < 0.3
# 92% > 0.1, 60% >= 0.2
#[1] "Simulated OCCl1 mean and sd: 0.193539097839627 0.0311377831344345"
# ACTUAL DATA:
# mean = 0.1415, sd = 0.0689 
# 26% < 0.1, 82% < 0.2, 97% < 0.3, 99% <0.4
# 73% > 0.1, 18% > 0.2, 2% > 0.3, 1% > 0.4 


# 2 element fit CCl1 using ccl0 - scatter real vs estimated data
setEPS()
postscript(outplotScat2CCl1l0)
simsies<-f2ccl1(OAD,OAVD,co2ccl0[1],co2ccl0[2])
plot(OCCl1,simsies,pch='.',
         xlim=c(0,1),ylim=c(0,1),
	 xlab='GHCN Cross Correlations at lag 1',
	 ylab='Estimated Cross Correlations at lag 1')
# plot 1:1 line
lines(seq(0,1,0.1),seq(0,1,0.1),col='red',lw=5)
# show self correlation in red 
# (shows that we need to force diagonals - but that we could in theory cope with same locs)
points(OCCl1[which(OAD == 0.0)],simsies[which(OAD == 0.0)],col='red',pch=20)
# show function as an equation on the plot
#text(0.15,0.1,paste("R = 0.3 * exp((-",
#              round(co2ccl0[1],4),"* dh) + (-",
#	      round(co2ccl0[2],2),"* dv))"),pos=4,cex=1.2)
mtext(paste("R = 0.23 * exp((-",round(co2ccl0[1],4),"* dh) + (-",round(co2ccl0[2],2),"* dv))"),
      side = 3, line = 1, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 1.2, col = NA, font = NA)
dev.off()
print("QUANTILES: Simulations for 2CCl1 use l0")
qtls<-quantile(simsies,qporbs)
print(qtls)
print(paste("Simulated OCCl1 mean and sd:",mean(simsies), sd(simsies),sep=" "))
#[1] "QUANTILES: Simulations for 2CCl1 use l0"
#	 1%	    2%         3%	  4%	     5% 	6%	   7% 
#0.09306901 0.10878384 0.11995922 0.12845738 0.13488862 0.14030646 0.14523972 
#	 8%	    9%        10%	 11%	    12%        13%	  14% 
#0.14965440 0.15357834 0.15712811 0.16046989 0.16369635 0.16693887 0.16984568 
#	15%	   16%        17%	 18%	    19%        20%	  21% 
#0.17262435 0.17522918 0.17772642 0.18017266 0.18250382 0.18479864 0.18698731 
#	22%	   23%        24%	 25%	    26%        27%	  28% 
#0.18902796 0.19083468 0.19250956 0.19397561 0.19528317 0.19638447 0.19742406 
#	29%	   30%        31%	 32%	    33%        34%	  35% 
#0.19833235 0.19912275 0.19983412 0.20047924 0.20107280 0.20162539 0.20212829 
#	36%	   37%        38%	 39%	    40%        41%	  42% 
#0.20259824 0.20304652 0.20348228 0.20389623 0.20429667 0.20467336 0.20504926 
#	43%	   44%        45%	 46%	    47%        48%	  49% 
#0.20540696 0.20576276 0.20610254 0.20644119 0.20676050 0.20707657 0.20738658 
#	50%	   51%        52%	 53%	    54%        55%	  56% 
#0.20770068 0.20801434 0.20831576 0.20861981 0.20891338 0.20920806 0.20949734 
#	57%	   58%        59%	 60%	    61%        62%	  63% 
#0.20978757 0.21008870 0.21037625 0.21065951 0.21094466 0.21123972 0.21152252 
#	64%	   65%        66%	 67%	    68%        69%	  70% 
#0.21181064 0.21209689 0.21238810 0.21268884 0.21297942 0.21327243 0.21357829 
#	71%	   72%        73%	 74%	    75%        76%	  77% 
#0.21387231 0.21416789 0.21447965 0.21478536 0.21510453 0.21542688 0.21576736 
#	78%	   79%        80%	 81%	    82%        83%	  84% 
#0.21609800 0.21644144 0.21677069 0.21710672 0.21746367 0.21782879 0.21820695 
#	85%	   86%        87%	 88%	    89%        90%	  91% 
#0.21859731 0.21900280 0.21941775 0.21983393 0.22028033 0.22074704 0.22124991 
#	92%	   93%        94%	 95%	    96%        97%	  98% 
#0.22179031 0.22235987 0.22299189 0.22370885 0.22454286 0.22557766 0.22727592  
#	99%	  100% 
#0.23000000 0.23000000
# 1% < 0.1, 31% < 0.2, 100% < 0.3
# 99%> 0.1, 69% >= 0.2
#[1] "Simulated OCCl1 mean and sd: 0.198133484458016 0.0283905190501944"
# ACTUAL DATA:
# mean = 0.1415, sd = 0.0689 
# 26% < 0.1, 82% < 0.2, 97% < 0.3, 99% <0.4
# 73% > 0.1, 18% > 0.2, 2% > 0.3, 1% > 0.4 

# STILL ERRING TOWARDS THIS ONE BUT ITS CLEAR THAT AC DOES NOT DECAY EXPONENTIALLY WITH DISTANCE AND ELEVATION!
# Probably have to use function based on 2 elements only as need to ensure that two close stations are sufficiently
# different.
# Looks like the decay has to be the same for the VAR to work so quite stuck really
# Could still be that autocorrelation is very sensitive to missing data/poor quality data. 
# We're using 0.0 and above. Could narrow down more?
# Did try lag1 with different diagonals to the intercept but this resulted in very large variability - best to keep it the same
# Now going with 1, 0.97, 0.0007, 0.08 for lag 0 and 0.23, 0.23, 0.0007 and 0.08 for lag1 - could potentially increase the diagonal and intercept



########################################################
#  END
########################################################
