#-----------------------------------------------------------------------------
# This program use Geweke diagnostic as a stopping criterion to terminate the
# MCMC simulation from the spDynLM() from the R package "spBayes". The 
# evaluation of the performance is carried out by comparing to the "truth" 
# from truth.txt.
#
# input
# ------
# args: a unique integer (e.g. i) indicates that it is the i-th run. Check 
# the spBayes_geweke.sh file for more details.
#
# output
# -------
# args[1]_output_geweke.txt: matrix of parameter-wise total simulation 
#effort, standard effective sample size and proposed effective sample size.
#
# args[1]_mcse_geweke.txt: sequence of parameter-wise markov chain standard 
# error.
#
# args[1]_sd_geweke.txt: sequence of parameter-wise posterior standard 
# deivation.
#
# args[1]_mean_geweke.txt: sequence of parameter-wise posterior mean.
#
# args[1]_time_geweke.txt: record the actual run time.
#
# args[1]_probcover_geweke.txt: 0/1 valued sequence used to estimate the
# coverage probabilities of each parameter.
#-----------------------------------------------------------------------------

require(mcmcse)
require(boa)
require(coda)
require(spBayes)

#### set random seed ####
args<-commandArgs(TRUE)
set.seed(args[1])

#### spatial Bayesian model setup ####
## manipulate raw data
data("NETemp.dat")
ne.temp <- NETemp.dat
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]
y.t <- ne.temp[,4:15]
N.t <- ncol(y.t) ## number of months
n <- nrow(y.t) ## number of observation per months
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
##make symbolic model formula statement for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)

#### generate MCMC samples ####
n.samples<- 50000
m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), 
               coords=coords,starting=starting, tuning=tuning, priors=priors,
               get.fitted =TRUE,cov.model="exponential", n.samples=n.samples,
               n.report=0.1*n.samples)
samples<- cbind(m.1$p.beta.0.samples, m.1$p.beta.samples, 
                m.1$p.sigma.eta.samples, m.1$p.theta.samples,
                t(m.1$p.u.samples))
rm(m.1)

## define functions
bm<- function(x){
  n= length(x)
  b= floor(sqrt(n))
  a= floor(n/b)
  
  y = sapply(1:a, function(k) return(mean(x[((k - 1) * b + 1):(k * b)])))
  mu.hat = mean(y)
  var.hat = b * sum((y - mu.hat)^2)/(a - 1)
  se = sqrt(var.hat/n)
  return(se)
}

#### record present time ####
ptm<- proc.time()

#### check diagnostic ####
burnin<- 1000
z<- 1.96
n<- 15000
X<- samples[(burnin+1):(burnin+n),]
# boa.geweke(X, 0.1, 0.5)
diag<- unlist(geweke.diag(X))
diag<- diag[1:(length(diag)-2)] #remove frac1 frac2
# check if P(p-value < 0.05)
print(length(which(2*pnorm(-abs(diag))<0.05))/dim(X)[2])
## calculate mcse
diag.mcse<- apply(X, 2, bm)
diag.sd<- apply(X, 2, sd)
diag.mean<- apply(X, 2, mean)
diag.ess.old<- apply(X, 2, ess)
diag.ess.new<- diag.sd^2/diag.mcse^2
write.table(diag.mcse, paste(args[1], "mcse_geweke.txt", sep="_"), 
            row.names=FALSE)
write.table(diag.sd, paste(args[1], "sd_geweke.txt", sep="_"), 
            row.names=FALSE)
write.table(diag.mean, paste(args[1], "mean_geweke.txt", sep="_"), 
            row.names=FALSE)
write.table(list(n=n, ess.old=diag.ess.old, ess.new=diag.ess.new),
            paste(args[1], "output_geweke.txt", sep="_"), row.names=FALSE)

#### run time ####
write((proc.time()-ptm)[1:3], 
      paste(args[1], "time_geweke.txt", sep="_"))

#### check coverage probability ####
truth<- read.table("truth.txt", header = TRUE)
upper<- diag.mean+z*diag.mcse
lower<- diag.mean-z*diag.mcse
# coverage of the resulting confidence interval (0/1)
prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper), ncol=1)*1
write.table(prob.coverage, 
            paste(args[1], "probcover_geweke.txt", sep="_"),
            row.names=FALSE)


