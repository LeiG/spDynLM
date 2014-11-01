#-----------------------------------------------------------------------------
# This program uses the R package "spBayes" to generate MCMC samples from the 
# Bayesian dynamic space-time models (Gelfand et al., 2005).
#
# The aim is to precisely estimate the "truth" of the parameters through many
# parallel runs with sufficient length in order use the "truth" in other 
# programs to evaluate the performances of stopping criteria utilized.
#-----------------------------------------------------------------------------

require(spBayes)

#### spatial Bayesian model setup ####
## manipulate raw data
data("NETemp.dat")
ne.temp <- NETemp.dat
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]
y.t <- ne.temp[,4:15]
N.t <- ncol(y.t) ##number of months
n <- nrow(y.t) ##number of observation per months
## scale to km
coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))
## set starting and priors
p <- 2 #number of regression parameters in each month
starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
                 "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                 "sigma.eta"=diag(rep(0.01, p)))
tuning <- list("phi"=rep(5, N.t))
priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
               "phi.Unif"=list(rep(3/(0.9*max.d), N.t), 
                               rep(3/(0.05*max.d), N.t)),
               "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
               "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
               "sigma.eta.IW"=list(2, diag(0.001,p)))
## symbolic model formula for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)

#### generate MCMC samples ####
n.samples <- 1000000
n.parallel<- 100
sample.means<- list()
sample.vars<- list()
for(i in 1:n.parallel){
  ## get MCMC samples
  m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), 
                 coords=coords,starting=starting, tuning=tuning, 
                 priors=priors, get.fitted =TRUE,cov.model="exponential",
                 n.samples=n.samples, n.report=0.1*n.samples)
  samples<- cbind(m.1$p.beta.0.samples, m.1$p.beta.samples, 
                  m.1$p.sigma.eta.samples, m.1$p.theta.samples, 
                  t(m.1$p.u.samples))
  rm(m.1)
  
  ## summary samples
  sample.means[length(sample.means)+1]<- list(apply(samples, 2, mean))
  sample.vars[length(sample.vars)+1]<- list(apply(samples, 2, var))
  rm(samples)
  write.table(i, "counter.txt")
}

## the truth
sample.means<- simplify2array(sample.means)
sample.vars<- simplify2array(sample.vars)
param.mean<- apply(sample.means, 1, mean)
param.var<- (n.samples-1)/(n.parallel*n.samples-1)*apply(sample.vars, 1, sum)+
  n.samples/(n.parallel*n.samples-1)*apply((sample.means-pop.mean)^2, 1, sum)
sample.truth<- cbind(param.mean, sqrt(param.var))
write.table(sample.truth, "truth.txt", col.names=c("mean", "sd"),
            row.names=FALSE)

