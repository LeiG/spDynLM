#-----------------------------------------------------------------------------
# This package use the new batch means (BM) calculation to estimate the markov
# standard error. The new method does not require the storage of the entire 
# chain for recalculation purposes. It significantly reduces the memory usage 
# from O(n) to O(n^(1-\tau)), where a common choice of \tau is 1/2.
#
# The effective sample size(ESS) as stopping criterion with tunable parameters
# (n.ess) is implemented to terminate the MCMC simulation from the spDynLM()
# from the R package "spBayes". The evaluation of the performance is carried 
# out by comparing to the "truth" from truth.txt.
# 
# input
# ------
# args: a unique integer (e.g. i) indicates that it is the i-th run. Check 
# the spBayes_newBM.sh file for more details.
#
# output
# -------
# "$n.ess"_args[1]_cond_ess_new.txt: 
#
# "$n.ess"_args[1]_output_ess_new.txt: matrix of parameter-wise total 
# simulation effort, standard effective sample size and proposed effective 
# sample size
#
# "$n.ess"_args[1]_mcse_ess_new.txt: sequence of parameter-wise markov chain 
# standard error.
#
# "$n.ess"_args[1]_sd_ess_new.txt: sequence of parameter-wise posterior 
# standard deivation.
#
# "$n.ess"_args[1]_mean_ess_new.txt: sequence of parameter-wise posterior mean
#
# “$n.ess”_args[1]_time_ess_new.txt: record the actual run time.
#
# "$n.ess"_args[1]_probcover_ess_new.txt: 0/1 valued sequence used to estimate 
# the coverage probabilities of each parameter.
#
# "$n.ess"_args[1]_probdist_ess_new.txt: 0/1 valued sequence used to estimate
# the prob. of whether the distance of estimates from truth is larger than the
# threshold.
#-----------------------------------------------------------------------------

require(mcmcse)
require(spBayes)

#### set random seed ####
args<-commandArgs(TRUE)
set.seed(args[1])

#### record present time ####
ptm<- proc.time()

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
## make symbolic model formula statement for each month
mods <- lapply(paste(colnames(y.t),'elev',sep='~'), as.formula)

#### generate MCMC samples ####
n.samples<- 500000
#n.samples<- 2000000
m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), 
               coords=coords,starting=starting, tuning=tuning, priors=priors,
               get.fitted =TRUE,cov.model="exponential", n.samples=n.samples,
               n.report=0.1*n.samples)
samples<- cbind(m.1$p.beta.0.samples, m.1$p.beta.samples, 
                m.1$p.sigma.eta.samples, m.1$p.theta.samples,
                t(m.1$p.u.samples))
rm(m.1)

#### apply stopping rule (ESS) retrospectively ####
n.ess<- 1000
#n.ess<- 5000
z<- 1.96
eps<- sqrt(4*1.96^2/n.ess)
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
      tank<- t(sapply(
        split.data.frame(tank, rep(1:(dim(tank)[1]/2), 
                                   times=rep(2, (dim(tank)[1]/2)))),
        apply, 2, mean))
    }
    
    tank.mcse<- apply(tank, 2, sd)*sqrt(b.size[1]/n)
    ess.new<- (tank.std[3,]/(n-1))/(tank.mcse^2)
    write.table(ess.new, paste(n.ess, args[1], "cond_ess_new", sep="_"))
    if(all(ess.new>=n.ess)){
      tank.mean<- apply(tank, 2, mean)
      ess.app<- 4*(z/eps)^2
      out<- list(n= n, app= ess.app, new= ess.new)
      write.table(out, paste(n.ess, args[1], "output_ess_new.txt", sep="_"),
                  row.names=FALSE)
      write.table(tank.mcse, paste(n.ess, args[1], "mcse_ess_new.txt", sep="_"),
                  row.names=FALSE)
      write.table(sqrt(tank.std[3,]/(n-1)), 
                  paste(n.ess, args[1], "sd_ess_new.txt", sep="_"),
                  row.names=FALSE)
      write.table(tank.mean, 
                  paste(n.ess, args[1], "mean_ess_new.txt", sep="_"),
                  row.names=FALSE)
      break
    }
  }
}

#### run time ####
write((proc.time()-ptm)[1:3], paste(n.ess, args[1], 
                                    "time_ess_new.txt", sep="_"))

#### check coverage probability ####
truth<- read.table("truth.txt", header = TRUE)
upper<- tank.mean+z*tank.mcse
lower<- tank.mean-z*tank.mcse
# coverage of the resulting confidence interval (0/1)
prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper), ncol=1)*1
write.table(prob.coverage, 
            paste(n.ess, args[1], "probcover_ess_new.txt", sep="_"), 
            row.names=FALSE)
# distance between estimates and truth
prob.distance<- matrix(abs(tank.mean-truth[,1])>=eps*truth[,2], ncol=1)*1
write.table(prob.distance, 
            paste(n.ess, args[1], "probdist_ess_new.txt", sep="_"), 
            row.names=FALSE)

