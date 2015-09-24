require(akima)


## gridlon must include westernmost longitude and be a nr,nc matrix
## gridlon can optionally include easternmost longitude
## gridlat must include 90 and -90, (90,...,-90) or (-90,...,90) nr,nc matrix

## NB.. THIS IS ONLY TESTED USING A MATRIX 
## KATE HAS MODIFIED TO DEAL WITH A 3D ARRAY OF VALS NOV 2014

padded.grid <- function(gridlon, gridlat, gridvals, padding=2)
{
    min.lon <- min(gridlon)
    gridlon <- gridlon - min.lon ## Shift to simplify construction (if -180to180) UNDONE LATER
    step.lon <- diff(gridlon[1,1:2])
    step.lat <- diff(gridlat[1:2,1])
    n.lon <- ncol(gridvals)
    n.lat <- nrow(gridvals)
    if (!isTRUE(all.equal(gridlon[1,n.lon], 360))) { # check that last is 360
        ## Add duplicate meridian
        gridlon <- cbind(gridlon, gridlon[,1]+360)
        gridlat <- cbind(gridlat, gridlat[,1])
        if (length(dim(gridvals)) == 3) {
	  tmp<-gridvals
	  mons<-dim(tmp)[3]
	  gridvals<-array(NA,c(n.lat,n.lon+1,mons))
	  for (mon in 1:mons) {
	    gridvals[,,mon] <- cbind(tmp[,,mon], tmp[,1,mon])
          }
	  n.lon <- ncol(gridvals)
        } else {
	  gridvals <- cbind(gridvals, gridvals[,1])
          n.lon <- ncol(gridvals)
	}
    }
    ## Aproximate index of lon+180
    mid.lon <- which.min(abs(gridlon[1,]-180)) #
    shifted.lon <- c(mid.lon:n.lon, 1:(mid.lon-1)) # creates a sequence

    if (padding >= 1) {
        ## Add latitude padding
	# Note drop=FALSE means that if one axis is 1 then it is still a dimension
        gridlon <- rbind(gridlon[(1:padding)+1,,drop=FALSE],
                         gridlon,
                         gridlon[n.lat-(1:padding),,drop=FALSE])
        gridlat <- rbind(gridlat[(1:padding),,drop=FALSE]-padding*step.lat,
                         gridlat,
                         (gridlat[n.lat-(padding:1)+1,,drop=FALSE]+
                          padding*step.lat))
        if (length(dim(gridvals)) == 3) {
	  tmp<-gridvals
	  mons<-dim(tmp)[3]
	  gridvals<-array(NA,c(n.lat+(padding*2),n.lon,mons))
	  for (mon in 1:mons) {
	    #browser()
	    gridvals[,,mon] <- rbind(tmp[(padding:1)+1,shifted.lon,mon],
                          tmp[,,mon],
                          tmp[n.lat-(1:padding),shifted.lon,mon])	
	  }
	} else {
	  gridvals <- rbind(gridvals[(padding:1)+1,shifted.lon,drop=FALSE],
                          gridvals,
                          gridvals[n.lat-(1:padding),shifted.lon,drop=FALSE])
        }
	n.lat <- nrow(gridvals)
        ## Add longitude padding
        gridlon <- cbind(gridlon[, n.lon+(1:padding)-padding-1, drop=FALSE]-360,
                         gridlon,
                         gridlon[, (1:padding)+1, drop=FALSE]+360)
        gridlat <- cbind(gridlat[, n.lon+(1:padding)-padding-1, drop=FALSE],
                         gridlat,
                         gridlat[, (1:padding)+1, drop=FALSE])
        if (length(dim(gridvals)) == 3) {
	  tmp<-gridvals
	  mons<-dim(tmp)[3]
	  gridvals<-array(NA,c(n.lat,n.lon+(padding*2),mons))
	  for (mon in 1:mons) {
	    gridvals[,,mon] <- cbind(tmp[, n.lon+(1:padding)-padding-1,mon],
                          tmp[,,mon],
                          tmp[, (1:padding)+1,mon])	
	  }
	} else {
	  gridvals <- cbind(gridvals[, n.lon+(1:padding)-padding-1, drop=FALSE],
                          gridvals,
                          gridvals[, (1:padding)+1, drop=FALSE])
        }
	n.lon <- ncol(gridvals)
    }

    if (step.lat > 0) {
        return(list(lon=gridlon + min.lon, ## Undo shift
                    lat=gridlat,
                    vals=gridvals))
    } else {
        return(list(lon=gridlon[n.lat:1,,,drop=FALSE] + min.lon, ## Undo shift
                    lat=gridlat[n.lat:1,,,drop=FALSE],
                    vals=gridvals[n.lat:1,,,drop=FALSE]))
    }
}

padded.interpp <- function(gridlon, gridlat, gridvals, lon, lat,
                           padding=2, ...)
{
    require(akima)
    grid <- padded.grid(gridlon, gridlat, gridvals, padding=padding)

# ADDED A BIT TO INTERPP OVER EACH FIELD SINGULARLY THEN PATCH TOGETHER
    if (length(dim(gridvals)) == 3) {
      mons<-dim(gridvals)[3]
      pointees<-array(NA,c(length(lon),mons))
      for (mon in 1:mons) {
        moo<-interpp(grid$lon, grid$lat, grid$vals[,,mon], lon, lat, ...)
        pointees[,mon]<-moo$z
      }
      return(pointees)
    } else {
      moo<-interpp(grid$lon, grid$lat, grid$vals, lon, lat, ...)
      pointees<-moo$z
      return(pointees)
    }
}
