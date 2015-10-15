# ISTI_Clean_Worlds
Code for generating the homogeneous synthetic data for the International Surface Temperature Initiative Benchmarking

R and Python code set up to work with the ISTI databank and create analog synthetic station data

#######################################################################################################################
# LATEST VERSION OF MAIN CODE: SEP 2015~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################################################################################################################

1) Which stations to keep:
   Remove the stations that are pretending to be land and are actually ships
   This first gave me 30 but on inspection, Appleby is actually a place in New Zealand. AND - its longitude is
   wrong - should be 173.2, not -173.2!!!
   Goree and Chatham Island are also NOT ships - left them all out for now.
   Fill in missing elevations using a digital elevation model but do not believe anything below -100m
   The lowest point on land (wikipedia) is -418m (shore of Dead Sea in Jordan)
   We were getting screwy estimated elevations of -6000m
	-> GetStationElevations_FEB2015.py
	
	INPUTS:
	INVENTORY_monthly_merged_stage3 (v101 as of Jul 2015) 35932 stations
	elev.0.25-deg.nc digital elevation map from JISAO, Uni Washington
	
	OUTPUTS:
	INVENTORY_monthly_merged_stage3_ships 30 stations beginning with 'XX'
	INVENTORY_monthly_merged_stage3_noelevs 657 stations
	INVENTORY_monthly_merged_stage3_proxyelevs 35902 stations with missing elevs estimated
	

2) List ISTI stations that have at least 3 years of data in a 15 year period - in fact, need a 
   minimum of four years of data - 2yrs NO INHOM time at beginning and end of record.
   This now also counts the number of stations which have enough data to calculate a climatology:
   at least 15 years of complete data within any 30 year period.
	-> GetISTILongs_FEB2015.py
	
	INPUTS:
	INVENTORY_monthly_merged_stage3_proxyelevs (35902 stations)
	ISTIv100/results_merged/merge_XXX00000000_stage3 etc.
	
	OUTPUT:
	ISTILONGS_stage3proxyelevs_JUL2015.dat (33032)
	ISTISHORTS_stage3proxyelevs_JUL2015.dat (2870)
	19614 stations have enough data to calculate climatology = 59.38%
	
3) Using only longs, get ISTI Distance matrix (100 and 1000) for the long enough stations
   Check for stations with the same locations - these will screw up the VAR matrices
   Remove the matching station (keep the last version of the duplicate - no test for longest version) and output lists
   Distance is based on longitude and latidue (not elevation)
   After this stage you can also make up the figure 2 map using: PlotISTIStationLocs_JUL2015.py
	-> GetISTIDistances_FEB2015.py
	
	INPUTS:
	ISTILONGS_stage3proxyelevs_JUL2015.dat (33032)
	
	OUTPUTS:
	#ISTILONGDISTANCESprelim_hundred_stage3proxyelevs_JUL2015.dat (32522)
	#ISTILONGDISTANCESprelim_thousand_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGSLOCMATCH_stage3proxyelevs_JUL2015.dat (510)
	
   Run this a second time so that the distance lists only contain unique stations
   This also outputs the 40 nearest neighbours for each station - for later
        -> GetISTIDistances2_Feb2015.py
	INPUTS:
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat (32522)
	
	OUTPUTS:
	ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGDISTANCES_thousand_stage3proxyelevs_JUL2015.dat (32522)
	FIXCORRNEIGHBOURS_ISTI_stage3proxyelevs_JUL2015.dat (32522)
	
   At this stage - also fill in the missing countries with either the nearest neighbour within 500km with a country listed, or NONE
   	-> GetISTICountries_AUG2015.py	
	
	INPUTS:
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat (32522)

        OUTPUTS:
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015_countries.dat (32522)
	
	NONE: 5
	Pseudo-countries: 3076

   At this stage - also a good time to get the landuse classifications for each station based on the latest ESACCI snapshot
   	-> station_9pixels_landuse_ISTI.py	(written by David Parker)
	
	INPUTS:
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015_countries.dat (32522)

        OUTPUTS:
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015_landuse.dat (32522)
			
4) Using only longs, based on distances, get 100 highest correlating stations for each long station
   Also print out the corrs at lag 1 for the same 100 stations
   List the stations that corrlate at 1.0 with other stations 
   *** On investigation this is most likely stations correlating
   *** at >0.999 but not actually 1.000. I have put a catch in to 
   *** list stations with 1.000. Still,
   *** its a little odd that two stations would correlate so highly.
   Tried this with first diffs and it isn't sensible - all negative
   Also print out the autocorrelation for each station
   Check for duplicates (same location, correlation of 1.0)
   As long as these are not the same location they can still be simulated
   The distance function means that they will not have a corr=1.0
   Output list
   Corrs (lag 0) are calculated from the first difference series of climate anomalies
   This is to remove effect of inhomogeneities Menne and Williams 2009?
   Can we use first difference series for correlation at lag 1? Not sure. Think I will
   just use climate anomalies for this. Should test...  NO WE CANNOT 
   Modified this to output all of the correlation files also for Standardised Anoms (CLS method, 0.3 LOWESS frac or 600 months)
   Made an extra output list of number of present obs - for later cutting down the stations we use to build distelev function to only longer
   stations
   TAKES AGES!
	-> GetISTICorrs_FEB2015.py
	
	INPUTS:
	ISTILONGDISTANCES_thousand_stage3proxyelevs_JUL2015.dat (32522)
	ISTIv100/results_merged/merge_XXX00000000_stage3 etc.
	
	OUTPUTS:
	ISTILONGCORRS_hundred_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGCORRSLAG1_hundred_stage3proxyelevs_JUL2015.dat (32522)
	#ISTILONGCORRSLAG1firstdiff_hundred_stage3proxyelevs_JUL2015.dat
	ISTILONGAUTOCORRS_hundred_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGCOMMONS_hundred_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGBADCORRS_stage3proxyelevs_JUL2015.dat (6 - but many more > 0.999)
	ISTILONGCORRSStAnom_hundred_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGCORRSLAG1StAnom_hundred_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGAUTOCORRSStAnom_hundred_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGINVENTORY_counts_stage3proxyelevs_JUL2015.dat (32522)

6) Look at corrs of GHCN stations only at lag 0 and lag 1 to fit distance function


   Look at lag1 autocorrelation and cross-correlation (using climate anomalies or clim anom first diffs?)

   	->fixdistelev_create_BAWG_FEB2015speedy.R

   See notes below for more details.
   Lag 0: 1=0.97, b=0.0006, c=0.07 - diagonals forced to be 1 (does need to be higher on the diagonals)
   Lag 1: 1=0.23, b=0.0006, c=0.07 - diagonals forced to be 0.23 (does not need to be higher on the diagonals)
	
7) Given a set of long/lat box boundaries (specified by me based on making sure there are at least 2
stations and a neighbour present but not too many) sort out lists of networks and networkneighbours for
the factorisation
Will need to test these to make sure they work using: FactoriseISTI_JUN2014.R (tested ok for 2CCL0 and 2CCL1(useL0), not ok for 2CCL0 and 2CCL1
       	-> GetBAWGnetworks_FEB2015.py
	
	INPUTS:
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat (32522)
	ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat (32522)
	BAWGLONGNetworks3_JUL2015.txt 192 sub-network boxes (192)
	
	OUTPUTS:
	ISTILONGNetwork_StOnly_stage3proxyelevs_JUL2015_1.dat etc     	only based on _corrs_
	ISTILONGNetwork_NbOnly_stage3proxyelevs_JUL2015_1.dat etc		only based on _corrs_
	ISTILONGNetwork_StNb_stage3proxyelevs_JUL2015_1.dat etc		only based on _corrs_
	ISTILONGNetwork_Boundaries_stage3proxyelevs_JUL2015.dat
	ISTILONGNetwork_Boundaries_stage3proxyelevs_JUL2015.eps .png
	ISTILONGNetwork_Stations_stage3proxyelevs_JUL2015.eps .png
		
8) Get your choice of GCM loess for all LONG stations
        -> Interp_GCM_NOV2014.R
	
	INPUTS:
	ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat (32522)
	tas_Amon_HadGEM2-ES_historical_r2i1p1_YYYYMM-YYYYMM.nc	
	
	OUTPUTS:
        HadGEM2RSLOESS_ISTI_stage3proxyelevs_loess04CLS_JUL2015.txt FAT 0.4
        HadGEM2RSLOESS_ISTI_stage3proxyelevs_loess03CLS_JUL2015.txt NORMAL 0.31
        HadGEM2RSLOESS_ISTI_stage3proxyelevs_loess02CLS_JUL2015.txt WAVY 0.21
        HadGEM2RSLOESS_ISTI_stage3proxyelevs_loess02CLS_JUL2015.txt WIGGLY 0.15

And some plots to decide which in IMAGES/
	
9) Run the VAR
	-> create_BAWG1VAR1_MAIN_SEP2015.R
	-> ConvertISTItoGHCN_NOV.py (must be run between section 7 and 8 noted below!!!
	
	INPUTS:
        merge_00000000000_stage3 etc.		       
        ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat
        ISTILONGDISTANCES_hundred_stage3proxyelevs_JUL2015.dat
        ISTILONGNetwork_StOnly_stage3proxyelevs_JUL2015_
        ISTILONGNetwork_NbOnly_stage3proxyelevs_JUL2015_
        HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess015CLS_JUL2015.txt # choice of 4!	
        INPUT and OUTPUT:
        MASKEDCLEAN_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt # paramtag=BNCHCAAA

	OUTPUTS: # paramtag='BNCHCAAA'
	CLEAN_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	CORRNEIGHBOURS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	VARPS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	GAM0S_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt 
	GAM1S_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	NEWCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	NEWCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	LOESS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	NEWSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	NEWMASKEDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	INTERIMSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt ???? NOT IN CURRENT USE
	OLDSHOCKS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	NEWSHOCKS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDCLIMANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDSDANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt 	
	NEWCLIMANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt  	
	NEWSDANOMS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt  	
	OLDDIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDDIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDCOVDIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	OLDCOVDIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	DIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	DIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	DIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	DIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
	
	
It has 9 processes to run:
      0) Clean up data (remove outliers) and calculate real station statistics: 
              PRODUCES:
              standardised anomalies (removed climatology, removed lowess trend, divided by climatological monthly standard deviation)
              climate anomalies (removed climatology)
              lowess trend (fitted to the anomalies using a 600 month span)
              linear decadal trend
              climatology (using entire station record - should be at least three years)
              standard deviation (using entire station record, should be at least three years) 
              autocorrelation at lag 1 (of both the standardised anomalies and climate anomalies)
              40 nearest neighbour list
              AR(1) residuals (or MVN random numbers if real data are not long/complete enough - 60 month minimum), and mean and standard deviation of those residuals
      1) Use the distance elevation function (already built - see READMEBAWG_SEP2015) to produce VAR(1) parameters for each station+40 neighbour matrix:
              PRODUCES:
              Gamma_lag0
              Gamma_lag1
              VAR(1) parameter
      2) Create MVN residual shocks for all time points/stations using factorisation and the distance elevation function, 
         Transform to shape of real resids by Q-Q mapping
         Give them real mean and standard deviation 
              PRODUCES:
              New residual shocks
      3) Run VAR(1) to create the simulated standardised anomalies using the neighbour disconnect method
              PRODUCES:
              Simulated Standardised anomalies
      4) Add a GCM lowess trend (already created - see READMEBAWG_SEP2015)
         Multiply by real climatological monthly standard deviations to create simulated climate anomalies 
         Calculated simulated clean world stats: linear decadal trend, climatology, standard deviation, autocorrelation at lag 1
         Add back real monthly climatology to create simulated monthly means
         Look at nearest neighbour difference series for simulated data and calculate the standard deviation and autocorrelation at lag 1 for standardised anomalies and climate anomalies
              PRODUCES:
              Simulated Climate anomalies
              Simulated Monthly means
              Simulated Clean world stats:
              Simulated Standardised anomaly differences series standard deviations
              Simulated Standardised anomaly autocorrelations at lag 1
              Simulated Climate anomaly difference series standard deviations
              Simulated Climate anomaly difference series autocorrelations at lag 1
      5) Look at Simulated nearest neighbour correlations at lag 0 and lag 1 for standardised anomalies and climate anomalies
              Neighbour correlations at lag 0 for Simulated Standardised anomalies and Simulated Climate anomalies	       
              Neighbour correlations at lag 1 for Simulated Standardised anomalies and Simulated Climate anomalies	       
      6) Look at nearest neighbour difference series for real data and calculate the standard deviation and autocorrelation at lag 1 for standardised anomalies and climate anomalies
              Real Standardised anomaly differences series standard deviations
              Real Standardised anomaly autocorrelations at lag 1
              Real Climate anomaly difference series standard deviations
              Real Climate anomaly difference series autocorrelations at lag 1
      7) Look at Real nearest neighbour correlations at lag 0 and lag 1 for standardised anomalies and climate anomalies
              Neighbour correlations at lag 0 for Real Standardised anomalies and Real Climate anomalies	       
              Neighbour correlations at lag 1 for Real Standardised anomalies and Real Climate anomalies      
      8) Read in masked simulated data and recalculate all statistics
              NOTE*** you must have run ConvertISTItoGHCN_NOV2014.py before running this section!!! ***
              PRODUCES:
              linear decadal trend
              climatology (using entire station record - should be at least three years)
              standard deviation (using entire station record, should be at least three years) 
              autocorrelation at lag 1 (of both the standardised anomalies and climate anomalies)

9a) Convert the output simulated data into ISTI and GHCN format in addition to masking with missing data from the real data, and create a GHCN-like station list
	-> create_BAWG1VAR1_MAIN_SEP2015.R
		INPUTS:
		ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat'
		CLEAN_ISTI_stage3proxyelevs_BNCHCAAA_SEP2015.txt'	#_stage3
		merge_00000000000_stage etc.
		
		OUTPUTS:
		ISTILONGINVENTORYghcn_stage3proxyelevs_SEP2015.dat'
		GHCN_TYPE/00000000000_raw.tavg
		ISTI_TYPE/merge_00000000000_stage3 etc
		MASKEDCLEAN_ISTI_stage3proxyelevs_BNCHCAAA_SEP2015.txt		
		
9) Run the VAR Plotter to get validation plots from simulated data
	-> create_BAWG1VAR1_PLOTS_SEP2015.R
		INPUTS: # paramtag='BNCHCAAA'
		merge_00000000000_stage3		
		ISTI_TYPE/merge_00000000000_stage3			
		ISTILONGINVENTORY_stage3proxyelevs_JUL2015.dat
		CORRNEIGHBOURS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		HadGEM2ESLOESS_ISTI_stage3proxyelevs_loess015CLS_JUL2015.txt	 # choice of 4 
		OLDCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt 
		OLDCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		NEWCOVS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		NEWCOVSlag1_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		OLDDIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		OLDDIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		OLDDIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		OLDDIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		DIFFSDSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		DIFFACSTATS_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		DIFFSDSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		DIFFACSTATSCM_ISTI_stage3proxyelevs_",paramtag,"_SEP2015.txt
		
		OUTPUTS:
		Scatter_Covs_StAnoms_BNCHCAAA_SEP2015.eps
		Scatter_Covs_ClAnoms_BNCHCAAA_SEP2015.eps
		Hist_NEW_Covs_StAnoms_BNCHCAAA_SEP2015.eps
		Hist_NEW_Covs_ClAnoms_BNCHCAAA_SEP2015.eps
		Hist_OLD_Covs_StAnoms_BNCHCAAA_SEP2015.eps
		Hist_OLD_Covs_ClAnoms_BNCHCAAA_SEP2015.eps
		Hist_NEW_Diffs_StAnoms_BNCHCAAA_SEP2015.eps
		Hist_NEW_Diffs_ClAnoms_BNCHCAAA_SEP2015.eps
		Hist_OLD_Diffs_StAnoms_BNCHCAAA_SEP2015.eps
		Hist_OLD_Diffs_ClAnoms_BNCHCAAA_SEP2015.eps
		IDNUMBER999"_OLDvsNEW_ClAnoms_BNCHCAAA_SEP2015.eps

See NOTES: for quantiles of output
		
10) Run the PHA pairwise homogenisation algorithm code on the masked simulated data to see how many changepoints are detectable


11) Look at the output from PHA in terms of frequency and size of adjustments


12) Look at the raw vs homogenised output on a station by station basis


13) Plot the summary of climate statistics for old and new (masked) data


######################################################################################
EXTRA NOTES

****************************
Choice of Climate Model:

The cmip5 model data is stored at:

/project/champ/data/cmip5/output1/

with a whole forest of sub-directories for each institution, model, experiment, frequency etc etc.

Managecmip is the tool that downloads it all - and will also search through the currently held local archives to see what is available already:

http://fcm1.metoffice.com/projects/ACEM/wiki/UserGuide

At the bottom of Gareths' D&A Twiki page (http://www-twiki/Main/DandALinks) there are links to the data repositories,

HadGEM2-ES
located in: /data/local/hadkw/CMIP5/HadGEM2ES/ATMOSAMON/:
r2i1p1,r3i1p1,r4i1p1

taken from: 
/project/champ/data/cmip5/output1/MOHC/HadGEM2-ES/historical/mon/atmos/Amon/
r1i1p1,r2i1p1,r3i1p1,r4i1p1
12/1859 to 11/2005
/project/champ/data/cmip5/output1/MOHC/HadGEM2-ES/historicalExt/mon/atmos/Amon/
r2i1p1, r3i1p1, r4i1p1
12/2005 to 11/2019

Can then also get other scenarios e.g., rcp85:
/project/champ/data/cmip5/output1/MOHC/HadGEM2-ES/rcp85/mon/atmos/Amon/:
r1i1p1,r2i1p1,r3i1p1,r4i1p1
12/2005 to 11/2100

These need to be interpolated to station point for each station

*****************************
Distance Elevation Function:

At some point I had to decide upon my distance function. I wanted the stations correlations to decay
with both horizontal distance and vertical distance. I believe that vertical distance has a stronger influence
than horizontal distance although I haven't tested this explicitly. I could look at radiosondes or reanalyses
pressure levels and correlation of series at different heights to look at this. It could be that elevation 
makes very little difference but it is one way of introducing topography into the VAR.

I read in the stations and their nearest 100 neighbours. For all neighbour pairs where there is a data overlap of 
> 60 months (5 years). I compute standardised anomalies. (detrend, remove clim, divide by standard deviation). I
compute the correlation at lag 0 and then at lag 1 for all pairs. I also save the horizonal and vertical differences
of those pairs.

This takes ages to run so I only ran it on 1 in 3 GHCN stations (~9000) is ISTIv101 (Jul 2015 download). Ths gave ~300000
datapoints.

After sorting all data based on horizontal distance I tried two functions:
lag 0:
ffix<-function(x,y,b,c){ 0.97 * exp((-b*x) + (-c*y))} # forces 0 distance correlations to be high but not identical
f<-function(x,y,a,b,c){ a * exp((-b*x) + (-c*y))} # allows 0 distance to be anything

lag 1:
ffix<-function(x,y,b,c){ 0.23 * exp((-b*x) + (-c*y))} # forces 0 distance correlations to be high but not identical
f<-function(x,y,a,b,c){ a * exp((-b*x) + (-c*y))} # allows 0 distance to be anything

I limit the correlations included to be only those stations where:
- there are at least 360 months / 30 years of data
- the correlation at lag 0 and lag 1 is greater than 0
= 284939 points

The distributions are as follows:
Horizontal distance:
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

Vertical Distance:
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

Lag 0 correlations:
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

lag 1 correlations:
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

RESULTS:
Allow all parameters to be modelled:
Lag 0:
a=1.008 b=0.0007 c=0.1
Lag 1:
a=0.146 b=-0.0001, c=0.1

Allow only b and c parameters to be modelled:
Lag 0 (a = 0.97):
b=0.0006 c=0.07 
Lag 1 (a = 0.23):
b=0.001 c=0.82

I am choosing ffix (2 parameters) over f (3 parameters) because this enables stations to correlate much more strongly when they are very close while
still preventing too high a correlation to screw up the covariance matrices. Overall, predicted correlations are
higher than the observed. We would expect this - real data contain errors and other local topographical feaures that
we have not attempted to account for here. The horizonal and vertical distance function has a faster decay than my 
original horizontal only choice which I think is better and prevents too much overcorrelation.

I am forcing lag 1 correlations to decay at the same rate as lag 0 because otherwise I end up with matrix inversion issues due to negative 
eigenvalues. Something to do with non-positive-definiteness?
Lag 0: 1=0.97, b=0.0006, c=0.07 - diagonals forced to be 1 (does need to be higher on the diagonals)
Lag 1: 1=0.23, b=0.0006, c=0.07 - diagonals forced to be 0.23 (does not need to be higher on the diagonals)

See figures:
DistElevDecayFunc2CCl0_AUG2015stanom100N00AC.eps
DistElevDecayFunc2CCl1_usel0_AUG2015stanom100N00AC.eps
ScatterDistElevDecayFunc2CCl0_AUG2015stanom100N00AC.eps
ScatterDistElevDecayFunc2CCl1_AUG2015stanom100N00AC.eps

This distributions of estimated correlations are as follows:
Lag 0 2 element fit:
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

Lag 1 2 element fit using lag 0 decay:
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

These could definitely be improved upon but they only need to be good enough to produce some spatio-temporal 
correlations that are ball-park realistic. Fingers crossed...

##################################
Output of New vs Old Covariances and Difference Series
NEW COVS:
[1] "St Anoms: lag 0s"
   1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13% 
0.483 0.566 0.606 0.632 0.653 0.669 0.684 0.696 0.707 0.716 0.725 0.734 0.742 
  14%   15%   16%   17%   18%   19%   20%   21%   22%   23%   24%   25%   26% 
0.750 0.757 0.764 0.770 0.776 0.781 0.786 0.791 0.795 0.799 0.803 0.806 0.810 
  27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39% 
0.812 0.815 0.818 0.820 0.822 0.824 0.826 0.828 0.830 0.831 0.833 0.834 0.836 
  40%   41%   42%   43%   44%   45%   46%   47%   48%   49%   50%   51%   52% 
0.837 0.839 0.840 0.841 0.842 0.844 0.845 0.846 0.847 0.848 0.849 0.850 0.851 
  53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65% 
0.852 0.853 0.854 0.855 0.856 0.857 0.858 0.859 0.860 0.861 0.862 0.862 0.863 
  66%   67%   68%   69%   70%   71%   72%   73%   74%   75%   76%   77%   78% 
0.864 0.865 0.866 0.867 0.868 0.868 0.869 0.870 0.871 0.872 0.873 0.874 0.874 
  79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
0.875 0.876 0.877 0.878 0.879 0.880 0.881 0.882 0.883 0.884 0.885 0.887 0.888 
  92%   93%   94%   95%   96%   97%   98%   99%  100% 
0.889 0.891 0.892 0.894 0.896 0.898 0.901 0.905 0.926 
mean=0.82
sd = 0.09
[1] "Clim Anoms: lag 0s"
   1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13% 
0.522 0.592 0.626 0.650 0.667 0.682 0.694 0.704 0.714 0.722 0.729 0.736 0.743 
  14%   15%   16%   17%   18%   19%   20%   21%   22%   23%   24%   25%   26% 
0.749 0.754 0.759 0.764 0.768 0.772 0.776 0.779 0.783 0.786 0.789 0.792 0.795 
  27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39% 
0.797 0.800 0.802 0.804 0.807 0.809 0.811 0.813 0.815 0.817 0.818 0.820 0.822 
  40%   41%   42%   43%   44%   45%   46%   47%   48%   49%   50%   51%   52% 
0.823 0.825 0.827 0.828 0.830 0.831 0.833 0.834 0.835 0.837 0.838 0.839 0.841 
  53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65% 
0.842 0.843 0.844 0.845 0.847 0.848 0.849 0.850 0.851 0.852 0.853 0.855 0.856 
  66%   67%   68%   69%   70%   71%   72%   73%   74%   75%   76%   77%   78% 
0.857 0.858 0.859 0.860 0.861 0.862 0.863 0.864 0.865 0.866 0.867 0.869 0.870 
  79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
0.871 0.872 0.873 0.874 0.875 0.876 0.878 0.879 0.880 0.881 0.883 0.884 0.886 
  92%   93%   94%   95%   96%   97%   98%   99%  100% 
0.887 0.889 0.891 0.893 0.896 0.899 0.903 0.909 0.965 
mean=0.82
sd=0.08
~6% between 0.6 and 0.7
~20% between 0.7 and 0.8
~70% between 0.8 and 0.9
~3% between 0.9 and 1.0
99% < 0.91
[1] "St Anoms: lag 1s"
   1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13% 
0.103 0.123 0.134 0.140 0.145 0.150 0.153 0.156 0.158 0.161 0.163 0.165 0.166 
  14%   15%   16%   17%   18%   19%   20%   21%   22%   23%   24%   25%   26% 
0.168 0.170 0.171 0.172 0.174 0.175 0.176 0.177 0.179 0.180 0.181 0.182 0.183 
  27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39% 
0.184 0.185 0.186 0.186 0.187 0.188 0.189 0.190 0.191 0.191 0.192 0.193 0.194 
  40%   41%   42%   43%   44%   45%   46%   47%   48%   49%   50%   51%   52% 
0.195 0.195 0.196 0.197 0.197 0.198 0.199 0.200 0.200 0.201 0.202 0.202 0.203 
  53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65% 
0.204 0.204 0.205 0.206 0.206 0.207 0.208 0.208 0.209 0.210 0.210 0.211 0.212 
  66%   67%   68%   69%   70%   71%   72%   73%   74%   75%   76%   77%   78% 
0.212 0.213 0.214 0.214 0.215 0.216 0.217 0.217 0.218 0.219 0.220 0.220 0.221 
  79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
0.222 0.223 0.223 0.224 0.225 0.226 0.227 0.228 0.229 0.230 0.231 0.232 0.233 
  92%   93%   94%   95%   96%   97%   98%   99%  100% 
0.235 0.236 0.238 0.239 0.241 0.244 0.247 0.252 0.313 
mean=0.20
sd=0.03
[1] "Clim Anoms: lag 1s"
   1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13% 
0.149 0.161 0.169 0.175 0.181 0.186 0.190 0.195 0.199 0.203 0.207 0.210 0.213 
  14%   15%   16%   17%   18%   19%   20%   21%   22%   23%   24%   25%   26% 
0.217 0.219 0.222 0.224 0.227 0.229 0.231 0.233 0.235 0.236 0.238 0.239 0.241 
  27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39% 
0.242 0.244 0.245 0.246 0.247 0.249 0.250 0.251 0.252 0.253 0.254 0.255 0.256 
  40%   41%   42%   43%   44%   45%   46%   47%   48%   49%   50%   51%   52% 
0.257 0.258 0.259 0.260 0.261 0.262 0.263 0.264 0.265 0.266 0.266 0.267 0.268 
  53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65% 
0.269 0.270 0.271 0.272 0.273 0.274 0.274 0.275 0.276 0.277 0.278 0.279 0.280 
  66%   67%   68%   69%   70%   71%   72%   73%   74%   75%   76%   77%   78% 
0.281 0.282 0.283 0.284 0.285 0.286 0.288 0.289 0.290 0.291 0.292 0.294 0.295 
  79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
0.296 0.298 0.299 0.301 0.303 0.305 0.307 0.309 0.311 0.314 0.317 0.320 0.324 
  92%   93%   94%   95%   96%   97%   98%   99%  100% 
0.329 0.334 0.342 0.353 0.369 0.393 0.431 0.497 0.787 
mean = 0.27
sd = 0.06
~9% between 0.1 and 0.2
~62% between 0.2 and 0.3
OLD COVS:
[1] "St Anoms: lag 0s"
    1%     2%     3%     4%     5%     6%     7%     8%     9%    10%    11% 
0.0950 0.1840 0.2450 0.2930 0.3340 0.3720 0.4060 0.4370 0.4660 0.4920 0.5160 
   12%    13%    14%    15%    16%    17%    18%    19%    20%    21%    22% 
0.5370 0.5560 0.5730 0.5880 0.6020 0.6150 0.6260 0.6370 0.6460 0.6550 0.6640 
   23%    24%    25%    26%    27%    28%    29%    30%    31%    32%    33% 
0.6720 0.6800 0.6870 0.6940 0.7000 0.7060 0.7120 0.7180 0.7230 0.7280 0.7330 
   34%    35%    36%    37%    38%    39%    40%    41%    42%    43%    44% 
0.7380 0.7430 0.7480 0.7520 0.7560 0.7610 0.7650 0.7690 0.7730 0.7770 0.7810 
   45%    46%    47%    48%    49%    50%    51%    52%    53%    54%    55% 
0.7840 0.7880 0.7920 0.7950 0.7990 0.8020 0.8050 0.8090 0.8120 0.8150 0.8184 
   56%    57%    58%    59%    60%    61%    62%    63%    64%    65%    66% 
0.8220 0.8250 0.8280 0.8310 0.8340 0.8370 0.8400 0.8430 0.8460 0.8490 0.8520 
   67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
0.8550 0.8580 0.8610 0.8630 0.8660 0.8690 0.8720 0.8750 0.8780 0.8800 0.8830 
   78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88% 
0.8860 0.8890 0.8920 0.8940 0.8970 0.9000 0.9030 0.9060 0.9090 0.9120 0.9150 
   89%    90%    91%    92%    93%    94%    95%    96%    97%    98%    99% 
0.9190 0.9220 0.9250 0.9290 0.9330 0.9370 0.9410 0.9470 0.9530 0.9600 0.9700 
  100% 
1.0000 
mean=0.75
sd=0.19
[1] "Clim Anoms: lag 0s"
   1%    2%    3%    4%    5%    6%    7%    8%    9%   10%   11%   12%   13% 
0.062 0.174 0.248 0.306 0.354 0.396 0.433 0.467 0.498 0.527 0.552 0.575 0.595 
  14%   15%   16%   17%   18%   19%   20%   21%   22%   23%   24%   25%   26% 
0.614 0.630 0.645 0.659 0.671 0.683 0.693 0.702 0.711 0.719 0.727 0.735 0.741 
  27%   28%   29%   30%   31%   32%   33%   34%   35%   36%   37%   38%   39% 
0.748 0.754 0.760 0.766 0.772 0.777 0.782 0.787 0.791 0.796 0.800 0.804 0.809 
  40%   41%   42%   43%   44%   45%   46%   47%   48%   49%   50%   51%   52% 
0.812 0.816 0.820 0.824 0.827 0.831 0.834 0.837 0.841 0.844 0.847 0.850 0.853 
  53%   54%   55%   56%   57%   58%   59%   60%   61%   62%   63%   64%   65% 
0.856 0.859 0.861 0.864 0.867 0.869 0.872 0.875 0.877 0.880 0.882 0.885 0.887 
  66%   67%   68%   69%   70%   71%   72%   73%   74%   75%   76%   77%   78% 
0.889 0.892 0.894 0.896 0.899 0.901 0.903 0.906 0.908 0.910 0.912 0.915 0.917 
  79%   80%   81%   82%   83%   84%   85%   86%   87%   88%   89%   90%   91% 
0.919 0.921 0.924 0.926 0.928 0.930 0.933 0.935 0.938 0.940 0.943 0.946 0.948 
  92%   93%   94%   95%   96%   97%   98%   99%  100% 
0.951 0.954 0.957 0.961 0.964 0.968 0.973 0.979 1.000 
mean=0.78
sd=0.19
~7% between 0.6 and 0.7
~17% between 0.7 and 0.8
~34% between 0.8 and 0.9
~30% between 0.9 and 1.0
74% < 0.91
[1] "St Anoms: lag 1s"
    1%     2%     3%     4%     5%     6%     7%     8%     9%    10%    11% 
-0.172 -0.133 -0.108 -0.090 -0.076 -0.063 -0.053 -0.043 -0.034 -0.027 -0.019 
   12%    13%    14%    15%    16%    17%    18%    19%    20%    21%    22% 
-0.013 -0.006  0.000  0.005  0.010  0.015  0.020  0.024  0.028  0.032  0.036 
   23%    24%    25%    26%    27%    28%    29%    30%    31%    32%    33% 
 0.040  0.043  0.047  0.050  0.054  0.057  0.060  0.063  0.066  0.069  0.072 
   34%    35%    36%    37%    38%    39%    40%    41%    42%    43%    44% 
 0.074  0.077  0.080  0.082  0.085  0.088  0.090  0.092  0.095  0.097  0.100 
   45%    46%    47%    48%    49%    50%    51%    52%    53%    54%    55% 
 0.102  0.104  0.106  0.109  0.111  0.113  0.115  0.117  0.120  0.122  0.124 
   56%    57%    58%    59%    60%    61%    62%    63%    64%    65%    66% 
 0.126  0.128  0.130  0.133  0.135  0.137  0.139  0.141  0.144  0.146  0.148 
   67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
 0.150  0.153  0.155  0.158  0.160  0.163  0.165  0.168  0.170  0.173  0.176 
   78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88% 
 0.179  0.182  0.185  0.188  0.191  0.195  0.199  0.202  0.206  0.211  0.215 
   89%    90%    91%    92%    93%    94%    95%    96%    97%    98%    99% 
 0.220  0.225  0.231  0.237  0.244  0.252  0.260  0.271  0.284  0.301  0.329 
  100% 
 0.799 
 mean=0.11
 sd=0.10
[1] "Clim Anoms: lag 1s"
    1%     2%     3%     4%     5%     6%     7%     8%     9%    10%    11% 
-0.108 -0.060 -0.033 -0.013  0.002  0.015  0.025  0.034  0.043  0.050  0.057 
   12%    13%    14%    15%    16%    17%    18%    19%    20%    21%    22% 
 0.063  0.069  0.075  0.080  0.084  0.089  0.094  0.098  0.102  0.106  0.110 
   23%    24%    25%    26%    27%    28%    29%    30%    31%    32%    33% 
 0.113  0.117  0.120  0.123  0.127  0.130  0.133  0.136  0.139  0.141  0.144 
   34%    35%    36%    37%    38%    39%    40%    41%    42%    43%    44% 
 0.147  0.150  0.152  0.155  0.158  0.160  0.163  0.165  0.168  0.171  0.173 
   45%    46%    47%    48%    49%    50%    51%    52%    53%    54%    55% 
 0.176  0.178  0.181  0.183  0.186  0.189  0.191  0.194  0.196  0.199  0.201 
   56%    57%    58%    59%    60%    61%    62%    63%    64%    65%    66% 
 0.204  0.207  0.209  0.212  0.215  0.217  0.220  0.223  0.226  0.229  0.232 
   67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
 0.235  0.238  0.241  0.245  0.248  0.251  0.255  0.258  0.262  0.266  0.269 
   78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88% 
 0.273  0.278  0.282  0.286  0.291  0.296  0.301  0.307  0.312  0.318  0.325 
   89%    90%    91%    92%    93%    94%    95%    96%    97%    98%    99% 
 0.332  0.339  0.347  0.356  0.365  0.375  0.387  0.400  0.417  0.440  0.482 
  100% 
 0.911 
mean=0.19
sd=0.12
~39% between 0.1 and 0.2
~33% between 0.2 and 0.3
NEW DIFFS:
[1] "Clim Anoms: St Devs"
  1%   2%   3%   4%   5%   6%   7%   8%   9%  10%  11%  12%  13%  14%  15%  16% 
0.29 0.35 0.39 0.43 0.47 0.50 0.52 0.55 0.57 0.59 0.61 0.62 0.64 0.66 0.67 0.69 
 17%  18%  19%  20%  21%  22%  23%  24%  25%  26%  27%  28%  29%  30%  31%  32% 
0.70 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83 0.83 0.84 0.85 
 33%  34%  35%  36%  37%  38%  39%  40%  41%  42%  43%  44%  45%  46%  47%  48% 
0.86 0.87 0.87 0.88 0.89 0.89 0.90 0.91 0.92 0.92 0.93 0.94 0.94 0.95 0.96 0.96 
 49%  50%  51%  52%  53%  54%  55%  56%  57%  58%  59%  60%  61%  62%  63%  64% 
0.97 0.98 0.98 0.99 1.00 1.00 1.01 1.02 1.03 1.03 1.04 1.05 1.06 1.07 1.08 1.08 
 65%  66%  67%  68%  69%  70%  71%  72%  73%  74%  75%  76%  77%  78%  79%  80% 
1.09 1.10 1.11 1.12 1.14 1.15 1.16 1.17 1.19 1.20 1.21 1.23 1.25 1.26 1.28 1.30 
 81%  82%  83%  84%  85%  86%  87%  88%  89%  90%  91%  92%  93%  94%  95%  96% 
1.32 1.34 1.36 1.39 1.42 1.44 1.47 1.50 1.53 1.57 1.60 1.64 1.68 1.74 1.80 1.88 
 97%  98%  99% 100% 
1.99 2.14 2.42 6.73 
mean=1.04
sd=0.42
~10% between 0.7 and 0.8
~14% between 0.8 and 0.9
~15% between 0.9 and 1.0
~45% between 1.0 and 1.1
~3 % > 2
[1] "Clim Anoms: ACs"
  1%   2%   3%   4%   5%   6%   7%   8%   9%  10%  11%  12%  13%  14%  15%  16% 
0.01 0.03 0.03 0.04 0.05 0.05 0.05 0.06 0.06 0.06 0.06 0.07 0.07 0.07 0.07 0.08 
 17%  18%  19%  20%  21%  22%  23%  24%  25%  26%  27%  28%  29%  30%  31%  32% 
0.08 0.08 0.08 0.08 0.08 0.09 0.09 0.09 0.09 0.09 0.09 0.09 0.10 0.10 0.10 0.10 
 33%  34%  35%  36%  37%  38%  39%  40%  41%  42%  43%  44%  45%  46%  47%  48% 
0.10 0.10 0.10 0.10 0.11 0.11 0.11 0.11 0.11 0.11 0.11 0.11 0.12 0.12 0.12 0.12 
 49%  50%  51%  52%  53%  54%  55%  56%  57%  58%  59%  60%  61%  62%  63%  64% 
0.12 0.12 0.12 0.12 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.14 0.14 0.14 0.14 
 65%  66%  67%  68%  69%  70%  71%  72%  73%  74%  75%  76%  77%  78%  79%  80% 
0.14 0.14 0.15 0.15 0.15 0.15 0.15 0.16 0.16 0.16 0.16 0.16 0.17 0.17 0.17 0.17 
 81%  82%  83%  84%  85%  86%  87%  88%  89%  90%  91%  92%  93%  94%  95%  96% 
0.18 0.18 0.18 0.18 0.19 0.19 0.19 0.20 0.20 0.21 0.21 0.21 0.22 0.22 0.23 0.24 
 97%  98%  99% 100% 
0.25 0.26 0.29 0.78 
mean=0.13
sd=0.06
~28% 0.0 and 0.1
~59% 0.1 and 0.2
~12% 0.2 and 0.3
~1% > 0.7
OLD DIFFS:
[1] "Clim Anoms: St Devs"
     1%      2%      3%      4%      5%      6%      7%      8%      9%     10% 
 0.3630  0.4050  0.4340  0.4560  0.4750  0.4920  0.5080  0.5220  0.5360  0.5490 
    11%     12%     13%     14%     15%     16%     17%     18%     19%     20% 
 0.5610  0.5730  0.5850  0.5950  0.6060  0.6160  0.6260  0.6350  0.6450  0.6540 
    21%     22%     23%     24%     25%     26%     27%     28%     29%     30% 
 0.6640  0.6730  0.6820  0.6910  0.7000  0.7090  0.7170  0.7260  0.7340  0.7426 
    31%     32%     33%     34%     35%     36%     37%     38%     39%     40% 
 0.7500  0.7590  0.7670  0.7750  0.7830  0.7900  0.7990  0.8070  0.8150  0.8230 
    41%     42%     43%     44%     45%     46%     47%     48%     49%     50% 
 0.8300  0.8380  0.8460  0.8540  0.8620  0.8710  0.8790  0.8870  0.8950  0.9030 
    51%     52%     53%     54%     55%     56%     57%     58%     59%     60% 
 0.9120  0.9200  0.9290  0.9370  0.9460  0.9550  0.9640  0.9730  0.9830  0.9920 
    61%     62%     63%     64%     65%     66%     67%     68%     69%     70% 
 1.0020  1.0120  1.0220  1.0320  1.0430  1.0540  1.0640  1.0760  1.0870  1.0990 
    71%     72%     73%     74%     75%     76%     77%     78%     79%     80% 
 1.1120  1.1240  1.1370  1.1500  1.1650  1.1790  1.1940  1.2100  1.2260  1.2440 
    81%     82%     83%     84%     85%     86%     87%     88%     89%     90% 
 1.2620  1.2820  1.3030  1.3250  1.3490  1.3740  1.4010  1.4310  1.4630  1.4980 
    91%     92%     93%     94%     95%     96%     97%     98%     99%    100% 
 1.5390  1.5850  1.6370  1.7000  1.7780  1.8770  2.0130  2.2140  2.5860 11.2740 
mean=0.99
sd=0.46
~13% between 0.7 and 0.8
~12% between 0.8 and 0.9
~11% between 0.9 and 1.0
~26% between 1.0 and 1.1
~4 % > 2
[1] "Clim Anoms: ACs"
    1%     2%     3%     4%     5%     6%     7%     8%     9%    10%    11% 
-0.122 -0.065 -0.031 -0.006  0.013  0.029  0.043  0.055  0.066  0.076  0.086 
   12%    13%    14%    15%    16%    17%    18%    19%    20%    21%    22% 
 0.094  0.102  0.110  0.117  0.124  0.131  0.137  0.144  0.150  0.156  0.162 
   23%    24%    25%    26%    27%    28%    29%    30%    31%    32%    33% 
 0.167  0.173  0.178  0.184  0.189  0.194  0.199  0.204  0.209  0.214  0.219 
   34%    35%    36%    37%    38%    39%    40%    41%    42%    43%    44% 
 0.224  0.229  0.234  0.239  0.244  0.248  0.253  0.258  0.263  0.267  0.272 
   45%    46%    47%    48%    49%    50%    51%    52%    53%    54%    55% 
 0.277  0.282  0.287  0.292  0.297  0.301  0.306  0.311  0.317  0.322  0.327 
   56%    57%    58%    59%    60%    61%    62%    63%    64%    65%    66% 
 0.332  0.337  0.342  0.348  0.353  0.359  0.364  0.370  0.375  0.381  0.386 
   67%    68%    69%    70%    71%    72%    73%    74%    75%    76%    77% 
 0.392  0.398  0.404  0.410  0.417  0.423  0.430  0.437  0.443  0.451  0.458 
   78%    79%    80%    81%    82%    83%    84%    85%    86%    87%    88% 
 0.465  0.473  0.481  0.489  0.497  0.506  0.515  0.524  0.534  0.544  0.555 
   89%    90%    91%    92%    93%    94%    95%    96%    97%    98%    99% 
 0.566  0.579  0.592  0.606  0.621  0.638  0.657  0.679  0.705  0.738  0.787 
  100% 
 0.969
mean=0.31
sd=0.20
~8% 0.0 and 0.1
~17% 0.1 and 0.2
~20% 0.2 and 0.3
~19% 0.3 and 0.4
~14% between 0.4 and 0.5
~9% between 0.5 and 0.6
~5% between 0.6 and 0.7
~3% > 0.7
