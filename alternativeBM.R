## apply stopping rule with ALTERNATIVE bm
require(mcmcse)

## set threshold for stopping rule
args<-commandArgs(TRUE)
eps<- args[1]

## mcmc samples
# beta.0.samples<- read.table("beta.0.samples.txt")
# beta.samples<- read.table("beta.samples.txt")
# sigma.eta.samples<- read.table("sigma.eta.samples.txt")
# theta.samples<- read.table("theta.samples.txt")
# samples<- as.matrix(cbind(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples))
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
n.samples <- 2000000
## get MCMC samples
m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
               starting=starting, tuning=tuning, priors=priors, get.fitted =TRUE,
               cov.model="exponential", n.samples=n.samples, n.report=0.1*n.samples)
samples<- cbind(m.1$p.beta.0.samples, m.1$p.beta.samples, m.1$p.sigma.eta.samples, m.1$p.theta.samples, t(m.1$p.u.samples))
rm(m.1)

## stopping rule
eps<- 0.05
z<- 1.96
d<- dim(samples)[2]
b.size<- c(2^7, 2^7)
b.range<- 2^seq(7, 15)
batch<- matrix(, b.size[1], d)
b.count<- 0 # number of iterations within batch
batch.count<- 0 # number of batches
n<- 0 # number of iterations
tank.std<- matrix(0, 3, d)

# n = 1
b.count= b.count+1
n= n+1
batch[b.count,]<- samples[n,]
tank.std[1,]<- samples[1,]

while(1){
  b.count= b.count+1
  n= n+1
  batch[b.count,]<- samples[n,]
  tank.std[2,]<- ((n-1)*tank.std[1,]+batch[b.count,])/n
  tank.std[3,]<- tank.std[3,]+n*(tank.std[2,]-batch[b.count,])^2/(n-1)
  tank.std[1,]<- tank.std[2,]
  
  # check if batch is full
  if(b.count == b.size[1]){
    b.count<- 0
    if(batch.count == 0){
      tank<- matrix(apply(batch, 2, mean), nrow=1)
      batch.count<- 1
    }
    else{
      tank<- rbind(tank, matrix(apply(batch, 2, mean), nrow=1))
      batch.count= batch.count+1
    }
  }
  
  # check every 20/21 batches
  if((n>=2^14) && (batch.count>20) && (!(dim(tank)[1]%%2))){
    print(n)
    batch.count<- 1
    b.size[2]<- b.size[1]
    b.size[1]<- 2^max(which(sqrt(n)>=b.range)+6)
    nbatch<- floor(n/b.size[1])
    
    #check if batch size changes
    if(b.size[1] != b.size[2]){
      batch<- matrix(, b.size[1], d)
      tank<- t(sapply(split.data.frame(tank, rep(1:(dim(tank)[1]/2), times=rep(2, (dim(tank)[1]/2)))), apply, 2, mean))
    }
    
    tank.mcse<- apply(tank, 2, sd)*sqrt(b.size[1]/n)
    cond<- 2*z*tank.mcse-eps*sqrt(tank.std[3,]/(n-1))
    write.table(cond, "cond_alternative.txt")
    if(all(cond<=0)){
      ess.new<- (tank.std[3,]/(n-1)^2)/tank.mcse^2*n
      ess.app<- 4*(z/eps)^2
      out<- list(n= n, app= ess.app, new= ess.new)
      write.table(out, paste(eps, "output_alternative.txt", sep="_"))
      write.table(tank.mcse, paste(eps, "mcse_alternative.txt", sep="_"))
      break
    }
  }
}