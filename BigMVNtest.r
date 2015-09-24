#######################################################################
#######################################################################
#######################################################################
######                                                           ######
######          Testing routines for simulating high-            ######
######      dimensional multivariate normal distributions        ######
######                                                           ######
#######################################################################
#######################################################################
#######################################################################
source("BigMVNgen.r")
#
# Generate some sampling locations uniformly distributed on the sphere
# NB to achieve this, latitudes must be sampled from a density 
# proportional to cos(lat) with latitude in radians, otherwise points
# will be denser at the poles. Not that it really matters for 
# illustrative purposes. 
#
set.seed(2000)
Nsites <- 1000
site.coords <- matrix(ncol=2,nrow=Nsites)
site.coords[,1] <- 180*(asin(2*runif(Nsites)-1)/pi)
site.coords[,2] <- 360*runif(Nsites)-180
#
#   Split latitude and longitude ranges into 5 bins each - so 25 groups in
#   total, which means about 400 sites in each group on average if there 
#   are 10000 sites in total (although the equatorial grid cells have more
#   more sites, obviously, because they're bigger)
#
N.grid <- 5
group.id <- 1 + (N.grid*floor(N.grid*(site.coords[,1]+90)/180)) +
                        floor(N.grid*(site.coords[,2]+180)/360)
cat("NUMBERS OF SITES PER GROUP:\n")
print(table(group.id))
#
#   Now define the "neighbourhoods" of each grid cell. This could be 
#   done on a distance basis; but this would require that we compute
#   all pairwise distances (which kind of defeats the purpose of 
#   trying to get everything done in something like linear 
#   computational time). The cunning plan is to calculate the 
#   distances fromthe grid cell centres to each site, and to 
#   define the neighbourhoods of each site to be all sites that
#   are within an appropriate distance of the corresponding grid
#   cell centre. For "appropriate distance", it seems reasonable
#   to take, say, the maximum distance between adjacent 
#   grid cell centres (this means we're guaranteed to include 
#   anything half-way into the next cell)
#
cell.centres <- cbind(rep(-90 + (((1:N.grid)-0.5)*180/N.grid),each=N.grid),
                     rep(-180 + (((1:N.grid)-0.5)*360/N.grid),N.grid))
sites.to.cells <- howfar(cell.centres,site.coords)
cells.to.cells <- howfar(cell.centres,cell.centres)
max.dist <- max(cells.to.cells[row(sites.to.cells)-col(sites.to.cells) == 1])
neighbours <- vector("list",N.grid^2)
for (i in 1:(N.grid^2)) {
  wanted.sites <- (group.id != i) & (sites.to.cells[i,] < max.dist)
  neighbours[[i]] <- (1:Nsites)[wanted.sites]
}
#
#   Function to define an exponential covariance structure. The 
#   covariance between sites separated by a distance of d km is
#   sigma^2 * exp(-phi*d). The vector theta below contains 
#   two elements: (sigma^2,phi). coords is a two-column
#   matrix of latitudes and longitudes
# KW: This substitutes a real covariance structure by assuming that 
# covariances decay exponentially with distance from 1 at 0 km to 
# 1/e(ish 0.37) at 1000km to 0.0 at 5000km 
#
exp.corr <- function(theta,coords) {
  d <- howfar(coords,coords)
  z <- theta[1] * exp(-theta[2]*d)
}
theta <- c(1,0.001)
#
#   For this test script, could just pass the correlation function 
#   through to rbigmvn.setup(). However, for illustrative 
#   purposes the next fewlines demonstrate how to prepare a
#   list of covariance matrices that can be passed over instead. 
#   In practice, the call to exp.corr() below would be replaced 
#   by code that calculates the appropriate empirical covariance
#   matrices, but hopefully the general idea is clear. As it stands,
#   the code is a bit inefficient because some pairs of sites 
#   will appear in more than one group and the code doesn't 
#   "remember" the pairs that it has already done. Worse things
#   happen at sea. 
#
#   NB the "wanted.sites" command below is critical: as supplied,
#   it ensures that the sites are extracted in *exactly* the same
#   order that they appear in the full site list, which is 
#   necessary to avoid generating nonsense!
#
M <- N.grid^2
Sigmalist <- vector("list",M)
for (i in 1:M) {
  cat(paste("Calculating covariances: group",i,"of",M,"...\r"))
  wanted.sites <- ((1:Nsites) %in% neighbours[[i]]) | (group.id == i)
  Sigmalist[[i]] <- exp.corr(theta,site.coords[wanted.sites,])
}
cat("\n")
#stop()
#
#   Now do all the preliminary calculations for the sampling
#
cat("Setting up neighbourhoods and covariance structures ...\n")
test.setup <- rbigmvn.setup(mu=1:Nsites,coords=site.coords,
                            groups=group.id,neighbours=neighbours,
                            coord.type="geographical",method="factorize",
                            Sigmalist=Sigmalist)
#
#   And generate 10000 realisations at those 1000 sites (the graphics
#   window is for the trace plot)
#
if (dev.cur()==1) x11(width=8,height=6)
par(lwd=2,ask=TRUE,mar=c(3,3,2,2),mgp=c(2,0.75,0),oma=c(1,1,1,1),mfrow=c(1,1))
cat("Simulating ...\n")
test <- rbigmvnorm(10000,setup=test.setup,nburnin=100,monitor=FALSE)
#
# Check the properties for, say, the first 10 sites
#
cat("PROPERTIES OF FIRST 12 SITES:\n")
cat("Sample means: "); print(round(colMeans(test[,1:10]),3))
cat("Theoretical means: "); print(round(1:10,3))
cat("Sample covariance matrix:\n")
print(round(cov(test[,1:10]),2))
cat("Theoretical covariance matrix:\n")
print(round(exp.corr(theta,site.coords[1:10,]),2))
par(mfrow=c(3,4),oma=c(1,1,3,1))
for (i in 1:12) {
  qqnorm(test[,i],main=paste("Site",i))
}
mtext("Normal Q-Q plots for first 12 sites",outer=TRUE,line=1)
#
# That looks OK, but the sites are obviously well separated and hence
# effectively independent. What about sites that are close together -
# and in different groups? Take sites within 10 degrees (lat,long) of
# (18,36) which is a boundary between two grid cells. 
#
wanted.sites <- (abs(site.coords[,1]-18) < 10) & (abs(site.coords[,2]-36) < 10)
cat("PROPERTIES OF SITES WITHIN 10 DEGREES OF (18,36):\n")
cat("Site numbers and groups:\n")
print(rbind((1:Nsites)[wanted.sites],group.id[wanted.sites]))
cat("Sample means: "); print(round(colMeans(test[,wanted.sites]),3))
cat("Theoretical means: "); print(round((1:Nsites)[wanted.sites],3))
cat("Sample covariance matrix:\n")
print(round(cov(test[,wanted.sites]),2))
cat("Theoretical covariance matrix:\n")
print(round(exp.corr(theta,site.coords[wanted.sites,]),2))
par(mfrow=c(2,4))
for (i in (1:Nsites)[wanted.sites]) {
  qqnorm(test[,i],main=paste("Site",i))
}
mtext("Normal Q-Q plots for sites within 10 degrees of (18,36)",outer=TRUE,line=1)
#
# Finally, check results for sites around +/-180 longitude, to ensure
# that wrapping is handled correctly (NB the latitude range from -13.5
# to +13.5 is merely to ensure that there are 12 sites in total, which 
# allows a 3*4 array of plots)
#
wanted.sites <- (abs(site.coords[,1]) < 13.5) & (abs(site.coords[,2]) > 170)
cat("PROPERTIES OF SITES AROUND (0,180):\n")
cat("Site numbers and groups:\n")
print(rbind((1:Nsites)[wanted.sites],group.id[wanted.sites]))
cat("Sample means: "); print(round(colMeans(test[,wanted.sites]),3))
cat("Theoretical means: "); print(round((1:Nsites)[wanted.sites],3))
cat("Sample covariance matrix:\n")
print(round(cov(test[,wanted.sites]),2))
cat("Theoretical covariance matrix:\n")
print(round(exp.corr(theta,site.coords[wanted.sites,]),2))
par(mfrow=c(3,4))
for (i in (1:Nsites)[wanted.sites]) {
  qqnorm(test[,i],main=paste("Site",i))
}
mtext("Normal Q-Q plots for sites around (0,180)",outer=TRUE,line=1)

#
#   Just out of interest, here's a demonstration of how one might
#   fit a VAR model to those data using least squares. Consider
#   the rows of the matrix we've just generated to be "time points",
#   and the columns to be "variables" (which they are). Exploit
#   the fact that if the response variable in a call to lm() is 
#   a matrix, a separate regression model will be fitted for 
#   each column; and the coefficients themselves will be stored
#   in columns (this doubtless exploits the fact that the design
#   matrix is the same for all of the regressions, so that 
#   things like (X'X)^-1 have to be calculated only once). The
#   estimate of Phi is then just the transpose of the returned
#   coefficient matrix. Neat!
#
#   NB the "-1" in the model formula is to exclude the intercept
#   in order to fit a VAR model with zero mean. 
#
cat("\nFitting VAR(1) model to output ...\n")
VAR.fit <- lm(test[-1,] ~ test[-nrow(test),] - 1)
Phi <- t(coef(VAR.fit))
par(mfrow=c(1,1))
plot(density(diag(Phi)),
     main="Diagonal coefficients of Phi matrix in VAR(1) model for simulated sequence")
rug(diag(Phi))


stop()

#
#   All looking good. Now: just out of interest, how long will it take to 
#   do the setup for 100,000 sites (which is the order of magnitude for
#   the ISTI project)? Probably split into something like a 30*30 grid
#   here, so that there are ~100 sites per block and 900 blocks. Storage
#   could be an issue given that i've just stored a 10000*1000 matrix, so
#   get rid of it.
#
rm(test.setup); rm(test); gc()
cat("Generating 100000 new site locations ...\n")
Nsites <- 100000
site.coords <- matrix(ncol=2,nrow=Nsites)
site.coords[,1] <- 180*(asin(2*runif(Nsites)-1)/pi)
site.coords[,2] <- 360*runif(Nsites)-180
N.grid <- 30
group.id <- 1 + (N.grid*floor(N.grid*(site.coords[,1]+90)/180)) +
  floor(N.grid*(site.coords[,2]+180)/360)
cell.centres <- cbind(rep(-90 + (((1:N.grid)-0.5)*180/N.grid),each=N.grid),
                      rep(-180 + (((1:N.grid)-0.5)*360/N.grid),N.grid))
sites.to.cells <- howfar(cell.centres,site.coords)
cells.to.cells <- howfar(cell.centres,cell.centres)
max.dist <- max(cells.to.cells[row(sites.to.cells)-col(sites.to.cells) == 1])
neighbours <- vector("list",N.grid^2)
for (i in 1:(N.grid^2)) {
  wanted.sites <- (group.id != i) & (sites.to.cells[i,] < max.dist)
  neighbours[[i]] <- (1:Nsites)[wanted.sites]
}
#
#   sites.to.cells is VERY big, so remove it and do garbage collection
#
rm(sites.to.cells); gc()
theta <- c(1,0.001)
cat("Setting up neighbourhoods and covariance structures ...\n")
test.setup <- rbigmvn.setup(mu=1:Nsites,covfunc=exp.corr,covpars=theta,
                            coords=site.coords,groups=group.id,neighbours=neighbours,
                            coord.type="geographical")
test <- rbigmvnorm(1,setup=test.setup)
