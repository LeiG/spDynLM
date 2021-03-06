## apply stopping rule with ORIGINAL bm
require(mcmcse)

## set threshold for stopping rule
# args<-commandArgs(TRUE)
# eps<- args[1]

## mcmc samples
# beta.0.samples<- read.table("beta.0.samples.txt")
# beta.samples<- read.table("beta.samples.txt")
# sigma.eta.samples<- read.table("sigma.eta.samples.txt")
# theta.samples<- read.table("theta.samples.txt")
# samples<- cbind(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples)
# rm(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples)
# samples<- as.matrix(read.table("samples.txt"))

## spatial Bayesian model
require(spBayes)
## manipulate raw data
data("NETemp.dat")
ne.temp <- NETemp.dat
set.seed(3)
##take a chunk of New England
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]
##subset first 1 years (Jan 2000 - Dec. 2001)
y.t <- ne.temp[,4:15]
N.t <- ncol(y.t) ##number of months
n <- nrow(y.t) ##number of observation per months
##add some missing observations to illistrate prediction
# miss <- sample(1:N.t, 10)
# holdout.station.id <- 5
# y.t.holdout <- y.t[holdout.station.id, miss]
# y.t[holdout.station.id, miss] <- NA
##scale to km
coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))
## set starting and priors
p <- 2 #number of regression parameters in each month
starting <- list("beta"=rep(0,N.t*p), "phi"=rep(3/(0.5*max.d), N.t),
                 "sigma.sq"=rep(2,N.t), "tau.sq"=rep(1, N.t),
                 "sigma.eta"=diag(rep(0.01, p)))
tuning <- list("phi"=rep(5, N.t))
priors <- list("beta.0.Norm"=list(rep(0,p), diag(1000,p)),
               "phi.Unif"=list(rep(3/(0.9*max.d), N.t), rep(3/(0.05*max.d), N.t)),
               "sigma.sq.IG"=list(rep(2,N.t), rep(10,N.t)),
               "tau.sq.IG"=list(rep(2,N.t), rep(5,N.t)),
               "sigma.eta.IW"=list(2, diag(0.001,p)))
##make symbolic model formula statement for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)
n.samples <- 3000000
## get MCMC samples
m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
               starting=starting, tuning=tuning, priors=priors, get.fitted =TRUE,
               cov.model="exponential", n.samples=n.samples, n.report=0.1*n.samples)
# beta.0.samples<- m.1$p.beta.0.samples
# beta.samples<- m.1$p.beta.samples
# sigma.eta.samples<- m.1$p.sigma.eta.samples
# theta.samples<- m.1$p.theta.samples
# u.samples<- m.1$p.u.samples
# y.samples<- m.1$p.y.samples
samples<- cbind(m.1$p.beta.0.samples, m.1$p.beta.samples, m.1$p.sigma.eta.samples, m.1$p.theta.samples, t(m.1$p.u.samples))
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

## stopping rule
eps<- 0.05
jump<- 20
check<- 2^14
z<- 1.96
while(1){
#   X<- samples[0:check,]
  print(check)
  
#   mcse<- apply(X, 2, bm)
#   std<- apply(X, 2, sd)
  mcse<- apply(samples[0:check,], 2, bm)
  std<- apply(samples[0:check,], 2, sd)
  cond<- 2*z*mcse-eps*std #95% confidence
  write.table(cond, paste(eps, "cond_origin.txt", sep="_"))
  if(all(cond<=0)){
#     ess.old<- apply(X, 2, ess)
    ess.old<- apply(samples[0:check,], 2, ess)
    ess.new<- std^2/mcse^2
    ess.app<- 4*(z/eps)^2
#     X.mean<- apply(X, 2, mean)
    X.mean<- apply(samples[0:check,], 2, mean)
    out<- list(n= check, app= ess.app, old= ess.old, new= ess.new)
    write.table(out, paste(eps, "output_origin.txt", sep="_"))
    write.table(mcse, paste(eps, "mcse_origin.txt", sep="_"))
    write.table(std, paste(eps, "sd_origin.txt", sep="_"))
    write.table(X.mean, paste(eps, "mean_origin.txt", sep="_"))
    break
  }
  
  b.size<- floor(sqrt(check))
  check<- check+jump*b.size
  nbatch<- floor(check/b.size)
  if(!(nbatch%%2)){check<- check+b.size}
}