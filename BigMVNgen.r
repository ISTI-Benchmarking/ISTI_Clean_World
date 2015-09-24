#######################################################################
#######################################################################
#######################################################################
######                                                           ######
######      Routines for simulation of high-dimensionsal         ######
######                   multivariate normals                    ######
######                                                           ######
#######################################################################
#######################################################################
rbigmvn.setup <- function(mu,covfunc,covpars,Sigmalist,coords,
                          groups,neighbours,nNN,MaxDist,coord.type,
                          method="Gibbs",vars.to.monitor) {
#######################################################################
#
#   To set up the quantities required to simulate from a 
#   high-dimensional multivariate normal distribution using
#   routine rbigmvnorm. Arguments:
#
#   mu        mean vector of the distribution
#   covfunc   Name of a function to compute covariances between
#             pairs of points. This function itself should have
#             two arguments: the first a vector of parameters,
#             and the second a matrix of "co-ordinates" or 
#             covariates. Denoting the parameter vector by theta,
#             and the co-ordinates for variable i by x[i], the 
#             covariance between variables i and j can be
#             computed as f(x[i],x[j],theta). It should return a
#             matrix of covariances. Note that the argument
#             "Sigmalist" (see below) offers an alternative
#             way to define the covariance structure. 
#   covpars   Vector of parameters for the covariance function
#   coords    Matrix giving the "coordinates" or covariate values
#             for each variable, for use when evaluating the 
#             covariance function. Row i corresponds to variable 
#             i; the number of rows should match the length of mu.
#   groups    A numeric vector, the same length as mu, containing group
#             numbers. Variables in the same group will be simulated
#             as a block within each Gibbs sampler iteration. If 
#             there are M groups then they must be numbered from 
#             1 to M.
#   neighbours  A list of M vectors, where M is the number of groups.
#             The jth element of this list should be a vector of 
#             indices for the variables that are not in group j but
#             are in its neighbourhood. This argument should not
#             be supplied if either nNN or MaxDist is provided;
#             but it must be supplied if neither of these is provided.
#   nNN       Optional scalar, giving the number of nearest neighbours
#             to include in the neighbourhood of each group. If supplied,
#             then distfunc must also be supplied; and neither 
#             neighbours nor MaxDist should be supplied. Conversely,
#             if this is not supplied then either neighbours or 
#             MaxDist must be supplied
#   MaxDist   Optional scalar, or vector of length M (number of groups)
#             giving the maximum distance to consider as being within 
#             the neighbourhood of each group. A variable will be 
#             considered as in the neighbourhood of a group if it is
#             within MaxDist of at least one site in the group, but 
#             not otherwise. If this argument is supplied, then
#             distfunc must also be supplied; and neither neighbours 
#             nor nNN should be supplied. Conversely, if this is not 
#             supplied then either neighbours or nNN must be supplied
#   Sigmalist An optional list of covariance matrices, with length M
#             where M is the total number of groups. This should not
#             be used if cov.func is supplied. If provided, its jth 
#             element should be the covariance matrix of all pairs of 
#             variables within both group j and its neighbours (as defined 
#             via the "neighbours" argument, which must therefore be 
#             provided in this case). The ordering of variables in 
#             this matrix must correspond with that in the full joint
#             distribution (i.e. it consists of the rows and columns
#             of the full covariance matrix corresponding to the 
#             current group and its neighbours).
#             Note that the routine does *not* check that the elements
#             of the list are in mutual agreement (for example, that
#             the supplied covariance between any pair of sites is the
#             same for all groups where that pair appears). It is the
#             user's responsibility to check this, therefore.
#   coord.type  Either "geographical" (if coords is a two-column
#             matrix containing latitudes and longitudes - NB in 
#             this case the first column should be latitude) or
#             "euclidean". This controls how distances between 
#             variables are calculated.
#   method    Either "Gibbs" or "factorize". If "Gibbs", the routine
#             will set things up for a Gibbs sampler; if "factorize"
#             then it will set things up so that simulation is done
#             in a single pass (no iteration) using factorisation of
#             the full joint distribution into products of conditionals.
#             The latter is potentially quicker, but this increase in
#             speed may be offset by the need to specify considerably
#             larger neighbourhoods in order to achieve the same 
#             level of accuracy.
#   vars.to.monitor   An optional vector of indices in the range
#             1:length(mu). If provided, the routine will  
#             calculate the quantities needed to produce a trace 
#             plot of the logarithm of the joint density of the 
#             corresponding variables during random variate 
#             generation when method="Gibbs". This allows 
#             convergence monitoring. 
#
# Note that if neighbourhoods are defined using nNN or MaxDist, it
# becomes necessary to calculate all of the pairwise intersite distances,
# which has a computational cost quadratic in the number of sites ...
#
#######################################################################
#
#   Checks on inputs and set up storage for results
#
  p <- length(mu)
  if (!is.matrix(coords)) stop("coords must be a matrix")
  if (nrow(coords) != p) {
    stop("coords has wrong number of rows - should be length(mu)")
  }
  if (length(groups) != p) stop("groups should have same length as mu")
  M <- max(groups)
  if (!all(sort(unique(groups))==1:M)) {
    stop("groups should be numbered from 1 to M where M is total number of groups")
  }
  if (missing(covfunc)) {
    if (missing(Sigmalist)) stop("Either covfunc or Sigmalist must be supplied")
  } else {
    if (!is.function(covfunc)) stop("covfunc should be the name of a function")
  }
  if (as.numeric(missing(neighbours))+as.numeric(missing(nNN))+
        as.numeric(missing(MaxDist)) != 2) {
    stop("You must provide exactly one of neighbours, nNN and MaxDist")
  }
  if (!missing(neighbours)) {
    if (!is.list(neighbours)) stop("neighbours must be a list")
    if (length(neighbours) != M) {
      stop("neighbours must have a list element for each group")
    }
  }
  if (!missing(nNN)) {
    if (!(length(nNN) %in% c(1,M))) {
      stop("nNN must be either a scalar, or a vector with an element for each group")
    }
  }
  if (!missing(MaxDist)) {
    if (!(length(MaxDist) %in% c(1,M))) {
      stop("MaxDist must be either a scalar, or a vector with an element for each group")
    }
  }
  if (!missing(Sigmalist)) {
    if (!missing(covfunc)) stop("Only one of covfunc and Sigmalist can be specified")
    if (missing(neighbours)) stop("'neighbours' must be specified when using 'Sigmalist'")
    if (length(Sigmalist) != M) stop("Length of Sigmalist doesn't match number of groups")
  }
  if (coord.type=="geographical") {
    distfunc <- howfar
  } else if (coord.type=="euclidean") {
    distfunc <- function(coords1,coords2) {
      sqrt(outer(coords1[,1],coords2[,1],"-")^2 +
             outer(coords1[,2],coords2[,2],"-")^2)
    }
  } else {
    stop("coord.type must be either 'geographical' or 'euclidean'")
  }
  cat(paste("Tested all inputs and these are ok\n"))
  
  group.details <- matrix(0,nrow=M,ncol=p)
  Sig12.Sig22Inv <- CondCovChol <- vector("list",M)
#
#   Set up neighbourhood structures for each group
#
  for (i in 1:M) {
    cat(paste("Working on block",i,"of",M,"...                \n"))
    group.details[i,groups==i] <- 1 # in network stations=1
    if (!missing(neighbours)) {
      tmp <- neighbours[[i]]
      tmp <- tmp[!(tmp %in% (1:p)[groups==i])] # double check no neighbours are actual stations?
      group.details[i,tmp] <- 2 # network neighbours =2? 
    }
    if (!missing(nNN)) {
      nNN.cur <- nNN
      if (length(nNN) > 1) nNN.cur <- nNN[i]
      dist.cur <- distfunc(coords[groups==i,],coords[groups!=i,])
      min.dist <- apply(dist.cur,MARGIN=2,FUN=min)
      wanted.sites <- order(min.dist)[1:nNN.cur]
      tmp <- (1:p)[groups!=i]
      group.details[i,tmp[wanted.sites]] <- 2
    }
    if (!missing(MaxDist)) {
      MD.cur <- MaxDist
      if (length(MaxDist) > 1) MD.cur <- MaxDist[i]
      dist.cur <- distfunc(coords[groups==i,],coords[groups!=i,])
      wanted.sites <- (dist.cur <= MD.cur)
      tmp <- (1:p)[groups!=i]
      group.details[i,tmp[wanted.sites]] <- 2
    }
#
#   Now set up the components of the partitioned covariance matrix, and 
#   required Choleski factors. The variables to condition on are the ones
#   that are in the neighbourhood of the current group and, if method
#   method="factorize", that are also in lower-numbered groups (because 
#   these will already have been allocated by the time we get to the 
#   current group).
#

# Use this whereever you want to stop with access to variables
#    browser()

    if (method=="Gibbs") {
      wanted.vars <- (group.details[i,]>0)
    } else if (method=="factorize") {
      wanted.vars <- (group.details[i,]>0 & groups <= i)
    } else {
      stop("method must be either 'Gibbs' or 'factorize'")
    }
    if (missing(Sigmalist)) {
      Sigma <- covfunc(covpars,coords[wanted.vars,])
    } else {
      Sigma <- Sigmalist[[i]]
      
      if (nrow(Sigma) != sum(group.details[i,]>0)) {
        stop(paste("Dimensions of Sigmalist[[",i,
                   "]] don't match required number of variables",sep=""))
      }
      if (method=="factorize") {                     # Get rid of elements for 
        tmp <- (groups[group.details[i,]>0] <= i)    # groups that haven't been 
        Sigma <- matrix(Sigma[tmp,tmp],
	         nrow=length(which(tmp==TRUE)),
		 ncol=length(which(tmp==TRUE)))                      # done yet
      }
    }
#    cat(paste("Now in the block and working on Sigma",i,"\n",sep=""))

    tmp <- group.details[i,wanted.vars]
#    Sig11 <- matrix(Sigma[tmp==1,tmp==1],nrow=length(which(tmp==1)),ncol=length(which(tmp==1)))
#    Sig12 <- matrix(Sigma[tmp==1,tmp==2],nrow=length(which(tmp==1)),ncol=length(which(tmp==1)))
#    Sig22 <- matrix(Sigma[tmp==2,tmp==2],nrow=length(which(tmp==1)),ncol=length(which(tmp==1)))
    Sig11 <- Sigma[tmp==1,tmp==1]
    Sig12 <- Sigma[tmp==1,tmp==2]
    Sig22 <- Sigma[tmp==2,tmp==2]
    	# for first few there may be no neighbours
#    print(tmp)
    #print(length(tmp))
#    if (nrow(Sig22)>0) {
    if (length(which(tmp==2))>0) {
      Sig12.Sig22Inv[[i]] <- Sig12 %*% solve(Sig22)
      CondCovChol[[i]] <- chol(Sig11 - tcrossprod(Sig12.Sig22Inv[[i]],Sig12))
    } else {
      CondCovChol[[i]] <- chol(Sig11)
    }
  }
#
# Choleski factor of inverse of covariance matrix for monitoring 
# locations
#
  if (!missing(vars.to.monitor)) {
    if (!(method == "Gibbs")) {
      warning("argument vars.to.monitor is only used when method='Gibbs'")
      SigInvChol <- NULL
      vtm <- NULL
    } else {
      Sigma <- covfunc(covpars,coords[vars.to.monitor,])
      SigInvChol <- solve(chol(Sigma))
      vtm <- vars.to.monitor
    }
  } else {
    SigInvChol <- NULL
    vtm <- NULL
  }
  invisible(list(mu=mu,groups=groups,group.details=group.details,
                 Sig12.Sig22Inv=Sig12.Sig22Inv,CondCovChol=CondCovChol,
                 SigInvChol=SigInvChol,method=method,vars.to.monitor=vtm))
}
#######################################################################
#######################################################################
#######################################################################
rbigmvnorm <- function(n,setup,nburnin=0,nthin=1,init,monitor=FALSE) {
#######################################################################
#
#   To simulate from a high-dimensional multivariate normal 
#   distribution. Arguments:
#
#   n         Number of realisations required
#   setup     A list object containing all of the information
#             required to carry out the simulation. This 
#             should be produced via a call to rbigmvn.setup
#   nburnin   Number of Gibbs Sampler iterations to use as "burnin".
#             Only used if setup$method="Gibbs".
#   nthin     Thinning rate for the sampler after the burnin period:
#             the n returned values will correspong to iterations
#             nburnin+nthin, nburnin+(2*nthin), ... , nburnin+(n*nthin)
#             Only used if setup$method="Gibbs".
#   init      Starting value for the Gibbs sampler. If supplied, 
#             this can be used (for example) to generate further 
#             draws from a chain that is already known to have 
#             converged to its limiting distribution, thereby
#             eliinating the need for any further burnin.
#             Only used if setup$method="Gibbs".
#   monitor   A logical scale indicated whether to produce a 
#             trace plot of the logarithm of the joint density 
#             for the variables defined by setup$vars.to.monitor.
#             The purpose is to assess convergence of the Gibbs 
#             sampler. To produce the plot *does* slow down the 
#             iterations, however. Only used if setup$method="Gibbs".
#
#######################################################################
  p <- length(setup$mu)
  M <- nrow(setup$group.details)
  if (!missing(init)) {
    if (length(init) != p) stop("Lengths of 'init' and 'mu' don't match")
  }
  if (monitor & is.null(setup$vars.to.monitor)) {
    warning("No monitoring variables defined - trace plot will not be produced")
    monitor <- FALSE
  }
  sim <- matrix(nrow=n,ncol=p)
#
# Initialisation if we're using a Gibbs sampler
#
  if (setup$method=="Gibbs") {
    logL <- rep(0,nburnin+(n*nthin))
    ntt <- nthin
    if (missing(init))  {
      cat("Initialising Gibbs sampler ...                           \r")
      z <- rep(NA,p)
      for (i in 1:M) {
        wanted <- (setup$groups==i)
        e <- rnorm(sum(wanted))
        z[wanted] <- setup$mu[wanted] + (setup$CondCovChol[[i]] %*% e)
      }
    } else {
      z <- init
    }
#
# Burnin for Gibbs sampler
#
    for (iter in (if (nburnin < 1) NULL else 1:nburnin)) {
      cat(paste("Gibbs sampler burnin period: iteration",iter,"of",nburnin,"...\r"))
      for (i in 1:M) {
        wanted <- (setup$groups==i)
        mu1 <- setup$mu[wanted]
        y2 <- z[setup$group.details[i,]==2]
        mu2 <- setup$mu[setup$group.details[i,]==2]
        e <- rnorm(length(mu1))
        z[wanted] <- mu1 + (setup$Sig12.Sig22Inv[[i]] %*% (y2-mu2)) + 
          t(setup$CondCovChol[[i]]) %*% e
      }    
      if (monitor) {
        logL[iter] <- -sum((t(setup$SigInvChol) %*% 
                              (z[setup$vars.to.monitor]-
                                 setup$mu[setup$vars.to.monitor]))^2)/2
      }
    }
    if (monitor) {
      pp <- 1/(nburnin+(nthin*n))
      ylim <- -0.5*qchisq(c(pp,1-pp),df=length(setup$vars.to.monitor))
      if (nburnin > 0) ylim <- c(ylim,logL[1:nburnin])
      ylim <- extendrange(ylim)
      plot(1:(nburnin+(nthin*n)),type="n",xlab="Iteration",ylab="Score",
           main="Gibbs sampler iterations: log-density score for monitoring variables",
           ylim=ylim)
      abline(v=nburnin,col="red")
      if (nburnin > 0) lines(1:nburnin,logL[1:nburnin])
    } else {
      ntt <- 1
    }
  }

  for (cur.sim in 1:n) {
    cat(paste("Generating sample",cur.sim,"of",n,"...               \r"))
    if (setup$method == "factorize") z <- rep(NA,p)
    for (iter in 1:nthin) {
      for (i in 1:M) {
        wanted <- (setup$groups==i)
        mu1 <- setup$mu[wanted]
        e <- rnorm(length(mu1))
        if (setup$method=="Gibbs") {
          cond.vars <- (setup$group.details[i,]==2)
        } else if (setup$method=="factorize") {
          cond.vars <- (setup$group.details[i,]==2) & (setup$groups < i)
        }
        if (sum(cond.vars) > 0) {
          y2 <- z[cond.vars]
          mu2 <- setup$mu[cond.vars]
          z[wanted] <- mu1 + (setup$Sig12.Sig22Inv[[i]] %*% (y2-mu2)) + 
            t(setup$CondCovChol[[i]]) %*% e
        } else {
            z[wanted] <- mu1 + t(setup$CondCovChol[[i]]) %*% e
        }
      }    
      if (monitor & setup$method=="Gibbs") {
        cur.idx <- nburnin + ((cur.sim-1)*ntt) + iter
        logL[cur.idx] <- -sum((t(setup$SigInvChol) %*% 
                                 (z[setup$vars.to.monitor]-
                                    setup$mu[setup$vars.to.monitor]))^2)/2
        if (cur.idx > 1) lines(c(cur.idx-1,cur.idx),logL[c(cur.idx-1,cur.idx)])
      }
    }
    sim[cur.sim,] <- z
  }
  sim
}
######################################################################
######################################################################
######################################################################
howfar <- function(sites1,sites2,units="km") {
######################################################################
#
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
######################################################################
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
