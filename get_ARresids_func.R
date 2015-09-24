# R
#
# Author: Kate Willett
# Created: 14 September 2015
# Last update: 14 September 2015
# Location: /data/local/hadkw/ISTI/PROGS/
# GitHub: https://github.com/SurfaceTemp/ISTI_Clean_Worlds/
# -----------------------
# CODE PURPOSE AND OUTPUT
# -----------------------
# This code searches through the real standardised anomalies to find the longest period of consecutive (non-missing) data.
# If it is 60 months or longer then it applies AR(1) in reverse to work out real shock residuals
#	phi<-cov(st_anoms[1:countall-1],st_anoms[2:countall])
#	shock_t=st_anoms[t]-(phi*st_anoms[t-1])
# If there are fewer than 60 months then random numbers are generated from an MVN distribution to proxy the real residuals
# An array of real or proxy shock residuals are returned, the same length as the provide array (with NAs compressed out in case of
# REAL)
#
# -----------------------
# LIST OF MODULES
# -----------------------
# R packages:
#
# Kate's R modules:
#
# -----------------------
# DATA
# -----------------------
# anomsarr - standardised anomalies for a single station
#
# -----------------------
# HOW TO RUN THE CODE
# -----------------------
# res=get_ARresids_func(StandardisedAnomalies)  
#
# -----------------------
# OUTPUT
# -----------------------
# A list of things:
# actualresids<-res$ar_z_resids - actual shock residuals based on AR(1) if station is long enough (5 years consecutive)
# RealResids<-res$RealResids - label where REAL = actual shock residuals are available and NORM means that are assumed MVN
# 
# -----------------------
# VERSION/RELEASE NOTES
# -----------------------
# 
# Version 1 (14 September 2015)
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
# 
########################################################################################
get_ARresids_func<-function(anomsarr) {

  # Set up variables
  phi_val<-0
  ar_z_resids<-0
  RealResids<-'NORM'
  
  # Compute the actual AR (not VAR) coeficients for the candidate station and get the residuals
  # each station may be a different length of series that has different missing data
  # find longest consecutive period (station will already have had interpolation over missing single elements) 
  # if less than 5 years (60 months) then just provide MVN ar_z_resids
  
  # get pointers for all data present
  Thegotsall<-which(!is.na(anomsarr))
  if (length(Thegotsall) >= 60) {  # 5 years - IF not then assume MVN resids
    # total number of months present
    countall<-length(Thegotsall)    
    # compressed array of only those data present   
    myanoms<-anomsarr[Thegotsall]
    # find the length of the gaps in the data
    findbreaks<-Thegotsall[2:countall]-Thegotsall[1:(countall-1)] # array of hopefully mostly 1s with some larger numbers reflecting the breaks
    # find all gaps - where number is > 1
    biggs<-which(findbreaks > 1)
    # find how many gaps there are
    countbiggs<-length(biggs)
    # if there are some gaps then find the longest consecutive segment
    if (countbiggs > 0) { # there are some gaps present
      # make an array of pointers to gaps with 0 and last data point too 
      biggs<-c(0,biggs,countall)
      # find the length of consecutive sections
      biggdiffs<-c(biggs[1],biggs[2:(countbiggs+2)]-biggs[1:(countbiggs+1)]) # Rememeber you've added two extra elements so countbiggs+2 is correct
      # identify the longest consecutive section
      gotcha<-which(biggdiffs == max(biggdiffs))
      # need to add 1 to beginning of biggs if record is continuous from the beginning?
      # reduce the pointers to data present to only one consecutive segment that is the longest
      Thegotsall<-Thegotsall[(biggs[(gotcha[1]-1)]+1):(biggs[gotcha[1]])]
      # get the length of the segment
      countall<-length(Thegotsall)
      # get the anomalies in that segment
      myanoms<-anomsarr[Thegotsall]
    } 
    # if the segment is at least 60 months then go ahead and get resids
    if (countall >= 60) { # now carry on to make resids or assume MVN
      # restandardise these as they may be a subset of the period and so not have a zero mean
      myanoms<-(myanoms-mean(myanoms))/sd(myanoms)
      
      # in AR(1) case the cross-correlation at lag 0 with self is 1.0
      gam0tmp<-1.0
      # get the lag 1 autocorrelation
      gam1tmp<-cov(myanoms[1:countall-1],myanoms[2:countall])

      # OR
      # gam1tmp<-sum(myanoms[1:countall-1]*myanoms[2:countall])*(1./(countall-1))

      # the lag 1 autocorrelation is essentially the AR(1) parameter
      tmpphi<-gam1tmp			# equal to the lag one AR coefficient Wilks page 412
      
      # now reverse At = At-1 + Zt to get ar_z_resids
      ar_z_resids<-array(NA,dim=c(countall-1))    # one element short as this assesses t and t-1 together
      for (tt in 2:countall) {
  	ar_z_resids[tt-1]<-myanoms[tt]-(tmpphi*myanoms[tt-1])     
      }
      # label these as REALS
      RealResids<-'REAL'	# flag that real resids have been calculated
      # browser()
    }
  }
  # if we haven't been able to establish real resids then make some random numbers that are MVN
  if (RealResids == 'NORM') {
    ar_z_resids<-t(apply(array(rnorm(length(anomsarr)),dim=c(1,length(anomsarr))),1,scale)) 
  }
  
  # return the things 
  return(list(ar_z_resids=ar_z_resids,RealResids=RealResids))
}
########################################################################################
