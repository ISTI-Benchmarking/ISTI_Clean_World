# R function to convert one distribution to another
# If only one set of data are supplied then it assumes you want to convert to MVN (N(0,1))
# If two sets are supplied the first is given the same distribution as the second
# Please supply standardised anomalies in the first instance


# Kate JUN2014
#-----------------------------------------------------------------------
# QMapTransform_func

QMapTransform_func <- function(firstdist,seconddist=0){

# firstdist: this is a vector of standardised anomalies
# seconddist: this is an optional vector of the standardised anomalies
#	with distribution that we want firstdist to have

# CAN COPE WITH MISSING DATA - JUST REPLACE AFTERWARDS
# order lists NAs at the end


if (missing(seconddist)) {
# creating an MVN distribution to map to
    seconddist<-rnorm(10000)
    seconddist<-(seconddist-mean(seconddist))/sd(seconddist)
} 
nfirsts<-length(firstdist)	# no. of elements in original distribution
nseconds<-length(seconddist)	# no. of elements in distribution to map to

## TEST PLOTS MAY WANT TO SWITCH OFF
## set up test plotting window:
#if (dev.cur()==1) x11(width=8,height=6)
#par(lwd=2,ask=TRUE,mar=c(3,3,2,2),mgp=c(2,0.75,0),oma=c(1,1,1,1),mfrow=c(1,1))
#par(mfrow=c(4,2))
#plot(density(firstdist),main="Candidate Timseries Distribution")
#plot(density(seconddist),main="Target Timseries Distribution")
#
#plot(firstdist,main="Candidate Timseries")
#plot(seconddist,main="Target Timseries")

# sort both vectors smallest to largest, keeping the index for later
indF<-order(firstdist)
indS<-order(seconddist)

# compress to data present only and get actual counts
sortfirstdist<-firstdist[indF]
sortfirstdist<-sortfirstdist[which(sortfirstdist!="NA")]
anfirsts<-length(sortfirstdist)
totFnas<-nfirsts-anfirsts	# use this later to add NAs at end of resampled vector

sortseconddist<-seconddist[indS]
sortseconddist<-sortseconddist[which(sortseconddist!="NA")]
anseconds<-length(sortseconddist)
totSnas<-nseconds-anseconds	# use this later to add NAs at end of resampled vector

# replace the elements of newfirstdist with newseconddist elements at corresponding quantiles
newsortfirstdist<-array(0.,dim=(anfirsts))

# easy peasy if first is smaller than second 
for (i in 1:anfirsts) {
    if (anfirsts < anseconds) {
        newsortfirstdist[i]<-sortseconddist[round(i*(anseconds/anfirsts))]
    } else {
        if (((i*(anseconds/anfirsts)) - floor(i*(anseconds/anfirsts))) == 0) {
            newsortfirstdist[i]<-sortseconddist[round(i*(anseconds/anfirsts))]
        } else if (floor(i*(anseconds/anfirsts))==0) {
            newsortfirstdist[i]<-sortseconddist[1]
        } else if (floor(i*(anseconds/anfirsts))==anseconds) {    
            newsortfirstdist[i]<-sortseconddist[anseconds]
        } else {
            multi<-(i*(anseconds/anfirsts))-(floor(i*(anseconds/anfirsts)))
            newsortfirstdist[i]<-sortseconddist[floor(i*(anseconds/anfirsts))]
                                 +multi*(sortseconddist[ceiling(i*(anseconds/anfirsts))]
                                 -sortseconddist[floor(i*(anseconds/anfirsts))])
	}
    }
}

# Now add back the NAs and unsort the new array using the old indices
if (totFnas > 0) {
    newsortfirstdist<-append(newsortfirstdist,replicate(totFnas,"NA"))
}
newfirstdist<-array(0.,dim=(nfirsts))
newfirstdist[indF]<-newsortfirstdist

## TEST PLOTS MAY WANT TO SWITCH OFF
#plot(density(newfirstdist),main="New Candidate Timseries Distribution")
#plot(density(seconddist),main="Target Timseries Distribution")
#
#plot(newfirstdist,main="New Candidate Timseries")
#plot(seconddist,main="Target Timseries")

  return(list(newfirstdist=newfirstdist,sortfirstdist=sortfirstdist))
}
#-----------------------------------------------------------------------
