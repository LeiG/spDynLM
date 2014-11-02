#-----------------------------------------------------------------------------
# Several stopping criteria, i.e.
# 1. effective sample size (ESS)
# 2. fixed-width stopping rule (FWSR)
# 3. geweke diagnostic
# with different variance estimation techniques, i.e. 
# 1. standard batch means estimator
# 2. proposed lower bound batch means estimator
# 3. proposed upper bound batch means estimator
# and tuning parameter settings are implemented to terminate the MCMC 
# simulation from the spDynLM() of the R package "spBayes". Evaluation of the
# performances is carried out by comparing the estimates to the pre-estimated
# "truth" from truth.txt.
# 
# input
# ------
# args: a unique integer (e.g. i) indicates that it is the i-th run. Check 
# the spBayes_newBM.sh file for more details.
#
# output
# -------
# *_args[1]_cond_new.txt: 
#
# *_args[1]_output_new.txt: matrix of parameter-wise total simulation
# effort, standard effective sample size and proposed effective sample size
#
# *_args[1]_mcse_new.txt: sequence of parameter-wise markov chain 
# standard error.
#
# *_args[1]_sd_new.txt: sequence of parameter-wise posterior standard
# deivation.
#
# *_args[1]_mean_new.txt: sequence of parameter-wise posterior mean.
#
# *_args[1]_time_new.txt: record the actual run time.
#
# *_args[1]_probcover_new.txt: 0/1 valued sequence used to estimate the 
# coverage probabilities of each parameter.
#
# *_args[1]_probdist_new.txt: 0/1 valued sequence used to estimate the
# prob. of whether the distance of estimates from truth is larger than the
# threshold.
#-----------------------------------------------------------------------------

require(mcmcse)
require(spBayes)
require(coda)

#### set random seed ####
args<-commandArgs(TRUE)
set.seed(args[1])

#### mk directories to store results ####
if(!file.exists("./geweke")){dir.create("./geweke")}
if(!file.exists("./sBM")){dir.create("./sBM")}
if(!file.exists("./sBM_ess")){dir.create("./sBM_ess")}
if(!file.exists("./nBM_lower")){dir.create("./nBM_lower")}
if(!file.exists("./nBM_lower_ess")){dir.create("./nBM_lower_ess")}
if(!file.exists("./nBM_upper")){dir.create("./nBM_upper")}
if(!file.exists("./nBM_upper_ess")){dir.create("./nBM_upper_ess")}

#### pre-estimated truth ####
truth<- read.table("truth.txt", header = TRUE)

#### spatial Bayesian model setup ####
## manipulate raw data
data("NETemp.dat")
ne.temp <- NETemp.dat
ne.temp <- ne.temp[ne.temp[,"UTMX"] > 5500000 & ne.temp[,"UTMY"] > 3250000,]
y.t <- ne.temp[,4:15]
N.t <- ncol(y.t) # number of months
n <- nrow(y.t) # number of observation per months
coords <- as.matrix(ne.temp[,c("UTMX", "UTMY")]/1000)
max.d <- max(iDist(coords))
## set starting and priors
p <- 2 # number of regression parameters in each month
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
n.samples<- 3000000
m.1<- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), 
              coords=coords, starting=starting, tuning=tuning, priors=priors,
              get.fitted =TRUE, cov.model="exponential", n.samples=n.samples,
              n.report=0.1*n.samples)
samples<- cbind(m.1$p.beta.0.samples, m.1$p.beta.samples, 
                m.1$p.sigma.eta.samples, m.1$p.theta.samples,
                t(m.1$p.u.samples))
rm(m.1)

#### implement stopping criteria retrospectively ####
## standard batch means
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
## parameter settings
z<- 1.96
n.ess<- c(1000, 2000, 4000, 1000000)
eps<- sqrt(4*z^2/n.ess) # equivalence between ess and eps

## geweke diagnostic #########################################################
## record present time
ptm<- proc.time()

n.geweke<- 20000
burn_in<- 5000
X<- samples[(burn_in+1):(burn_in+n.geweke),]
diag<- unlist(geweke.diag(X))
diag<- diag[1:(length(diag)-2)] #remove frac1 frac2
# check if P(p-value < 0.05)
#print(length(which(2*pnorm(-abs(diag))<0.05))/dim(X)[2])
## calculate mcse
diag.mcse<- apply(X, 2, bm)
diag.sd<- apply(X, 2, sd)
diag.mean<- apply(X, 2, mean)
diag.ess.old<- apply(X, 2, ess)
diag.ess.new<- diag.sd^2/diag.mcse^2
write.table(diag.mcse, 
            paste("./geweke/", args[1], "_geweke_mcse.txt", sep=""), 
            row.names=FALSE)
write.table(diag.sd, 
            paste("./geweke/", args[1], "_geweke_sd.txt", sep=""), 
            row.names=FALSE)
write.table(diag.mean, 
            paste("./geweke/", args[1], "_geweke_mean.txt", sep=""), 
            row.names=FALSE)
write.table(list(n=n, ess.old=diag.ess.old, ess.new=diag.ess.new),
            paste("./geweke/", args[1], "_geweke_output.txt", sep=""),
            row.names=FALSE)
write((proc.time()-ptm)[1:3], 
      paste("./geweke/", args[1], "_geweke_time.txt", sep=""))
## coverage probabilities
upper<- diag.mean+z*diag.mcse
lower<- diag.mean-z*diag.mcse
prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper), ncol=1)*1
write.table(prob.coverage, 
            paste("./geweke/", args[1], "_geweke_probcover.txt", sep=""),
            row.names=FALSE)
rm(n.geweke, burn_in, X, diag, diag.mcse, diag.sd, diag.mean, 
   diag.ess.new, diag.ess.old, upper, lower, prob.coverage)
print("==============Geweke finished================")

## standard batch means ######################################################
i<- matrix(1, ncol = 2) # indicator of the tuning parameter value
colnames(i)<- c("sBM", "sBMess")
check<- 2^14

## record present time
ptm<- proc.time()

while(1){
  sBM.mcse<- apply(samples[0:check,], 2, bm)
  sBM.sd<- apply(samples[0:check,], 2, sd)
  
  ## FWSR
  cond<- 2*z*sBM.mcse-eps[i[1]]*sBM.sd #95% confidence
  if(all(cond<=0)){
    ess.old<- apply(samples[0:check,], 2, ess)
    ess.new<- sBM.sd^2/sBM.mcse^2
    ess.app<- 4*(z/eps[i[1]])^2
    sBM.mean<- apply(samples[0:check,], 2, mean)
    sBM.out<- list(n= check, app= ess.app, old= ess.old, new= ess.new)
    write.table(sBM.out, 
                paste("./sBM/", args[1], "_", n.ess[i[1]], 
                      "_standard_output.txt", sep=""),
                row.names=FALSE)
    write.table(sBM.mcse, 
                paste("./sBM/", args[1], "_", n.ess[i[1]],
                      "_standard_mcse.txt", sep=""),
                row.names=FALSE)
    write.table(sBM.sd, 
                paste("./sBM/", args[1], "_", n.ess[i[1]], 
                      "_standard_sd.txt", sep=""), 
                row.names=FALSE)
    write.table(sBM.mean, 
                paste("./sBM/", args[1], "_", n.ess[i[1]], 
                      "_standard_mean.txt", sep=""), 
                row.names=FALSE)
    upper<- sBM.mean+z*sBM.mcse
    lower<- sBM.mean-z*sBM.mcse
    # coverage of the resulting confidence interval (0/1)
    prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper), ncol=1)*1
    write.table(prob.coverage, 
                paste("./sBM/", args[1], "_", n.ess[i[1]],
                      "_standard_probcover.txt", sep=""),
                row.names=FALSE)
    # distance between estimates and truth
    prob.distance<- matrix(abs(sBM.mean-truth[,1])>=eps[i[1]]*truth[,2], 
                           ncol=1)*1
    write.table(prob.distance, 
                paste("./sBM/", args[1], "_", n.ess[i[1]],
                      "_standard_probdist.txt", sep=""), 
                row.names=FALSE)
    write((proc.time()-ptm)[1:3], 
          paste("./sBM/", args[1], "_", n.ess[i[1]], 
                "_standard_time.txt", sep=""))
    
    i[1]= i[1]+1
  }
  
  ## ESS
  cond<- (sBM.sd/sBM.mcse)^2 #ess
  if(all(cond>=n.ess[i[2]])){
    ess.old<- apply(samples[0:check,], 2, ess)
    ess.app<- 4*(z/eps[i[2]])^2
    sBM.mean<- apply(samples[0:check,], 2, mean)
    sBM.out<- list(n= check, app= ess.app, old= ess.old)
    write.table(sBM.out, 
                paste("./sBM_ess/", args[1], "_", n.ess[i[2]],
                      "_ess_standard_output.txt",sep=""),
                row.names=FALSE)
    write.table(sBM.mcse, 
                paste("./sBM_ess/", args[1], "_", n.ess[i[2]], 
                      "_ess_standard_mcse.txt", sep=""),
                row.names=FALSE)
    write.table(sBM.sd, 
                paste("./sBM_ess/", args[1], "_", n.ess[i[2]], 
                      "_ess_standard_sd.txt", sep=""),
                row.names=FALSE)
    write.table(sBM.mean, 
                paste("./sBM_ess/", args[1], "_", n.ess[i[2]], 
                      "_ess_standard_mean.txt",sep=""),
                row.names=FALSE)
    upper<- sBM.mean+z*sBM.mcse
    lower<- sBM.mean-z*sBM.mcse
    # coverage of the resulting confidence interval (0/1)
    prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper), ncol=1)*1
    write.table(prob.coverage, 
                paste("./sBM_ess/", args[1], "_", n.ess[i[2]],
                      "_ess_standard_probcover.txt", sep=""),
                row.names=FALSE)
    # distance between estimates and truth
    prob.distance<- matrix(abs(sBM.mean-truth[,1])>=eps[i[2]]*truth[,2], 
                           ncol=1)*1
    write.table(prob.distance, 
                paste("./sBM_ess/", args[1], "_", n.ess[i[2]],
                      "_ess_standard_probdist.txt", sep=""), 
                row.names=FALSE)
    write((proc.time()-ptm)[1:3], 
          paste("./sBM_ess/", args[1], "_", n.ess[i[2]], 
                "_ess_standard_time.txt", sep=""))
    
    i[2]= i[2]+1
  }
  
  b.size<- floor(sqrt(check))
  check<- check+20*b.size
  nbatch<- floor(check/b.size)
  if(!(nbatch%%2)){check<- check+b.size}
  
  if(all(i == 4)){break}
}
rm(i, sBM.mcse, sBM.sd, cond, ess.old, ess.new, sBM.mean, sBM.out, prob.distance,
   prob.coverage, upper, lower)
print("==============sBM finished================")

## proposed lower bound batch means ##########################################
i<- matrix(1, ncol = 2) # indicator of the tuning parameter value
colnames(i)<- c("nBM_lower", "nBMess_lower")

d<- dim(samples)[2]
b.size<- c(2^7, 2^7)
b.range<- 2^seq(7, 15)
batch<- matrix(, b.size[1], d)
b.count<- 0 # number of iterations within batch
batch.count<- 0 # number of batches
n<- 0 # number of iterations
tank.std<- matrix(0, 3, d)

b.count= b.count+1
n= n+1
batch[b.count,]<- samples[n,]
tank.std[1,]<- samples[1,]

#### record present time ####
ptm<- proc.time()

while(1){
  b.count= b.count+1
  n= n+1
  batch[b.count,]<- samples[n,]
  ## recursive standard deviation estimation
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
  if((n>=2^14) && (batch.count>=20) && (!(dim(tank)[1]%%2))){
    batch.count<- 1
    b.size[2]<- b.size[1]
    b.size[1]<- 2^max(which(sqrt(n)>=b.range)+6)  # lower bound
    nbatch<- floor(n/b.size[1])
    
    #check if batch size changes
    if(b.size[1] != b.size[2]){
      batch<- matrix(, b.size[1], d)
      tank<- t(sapply(split.data.frame(
        tank, rep(1:(dim(tank)[1]/2), 
                  times=rep(2, (dim(tank)[1]/2)))), apply, 2, mean))
    }

    tank.mcse<- apply(tank, 2, sd)*sqrt(b.size[1]/n)

    ## FWSR
    cond<- 2*z*tank.mcse-eps[i[1]]*sqrt(tank.std[3,]/(n-1))
    if(all(cond<=0)){
      tank.mean<- apply(tank, 2, mean)
      ess.new<- (tank.std[3,]/(n-1))/(tank.mcse^2)
      ess.app<- 4*(z/eps[i[1]])^2
      out<- list(n= n, app= ess.app, new= ess.new)
      write.table(out, 
                  paste("./nBM_lower/", args[1], "_", n.ess[i[1]],
                        "_new_lower_output.txt", sep=""), 
                  row.names=FALSE)
      write.table(tank.mcse, 
                  paste("./nBM_lower/", args[1], "_", n.ess[i[1]],
                        "_new_lower_mcse.txt", sep=""), 
                  row.names=FALSE)
      write.table(sqrt(tank.std[3,]/(n-1)), 
                  paste("./nBM_lower/", args[1], "_", n.ess[i[1]],
                        "_new_lower_sd.txt",sep=""), row.names=FALSE)
      write.table(tank.mean, 
                  paste("./nBM_lower/", args[1], "_", n.ess[i[1]], 
                        "_new_lower_mean.txt", sep=""), 
                  row.names=FALSE)
      upper<- tank.mean+z*tank.mcse
      lower<- tank.mean-z*tank.mcse
      # coverage of the resulting confidence interval (0/1)
      prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper),ncol=1)*1
      write.table(prob.coverage, 
                  paste("./nBM_lower/", args[1], "_", n.ess[i[1]], 
                        "_new_lower_probcover.txt", sep=""), 
                  row.names=FALSE)
      # distance between estimates and truth
      prob.distance<- matrix(abs(tank.mean-truth[,1])>=eps[i[1]]*truth[,2],
                             ncol= 1)*1
      write.table(prob.distance, 
                  paste("./nBM_lower/", args[1], "_", n.ess[i[1]], 
                        "_new_lower_probdist.txt", sep=""), 
                  row.names=FALSE)
      write((proc.time()-ptm)[1:3], 
            paste("./nBM_lower/", args[1], "_", n.ess[i[1]], 
                  "_new_time.txt", sep=""))
      
      i[1]= i[1] + 1
    }

    ## ESS
    ess.new<- (tank.std[3,]/(n-1))/(tank.mcse^2)
    if(all(ess.new>=n.ess[i[2]])){
      tank.mean<- apply(tank, 2, mean)
      ess.app<- 4*(z/eps[i[2]])^2
      out<- list(n= n, app= ess.app, new= ess.new)
      write.table(out, paste("./nBM_lower_ess/", args[1], "_", n.ess[i[2]], 
                             "_ess_new_lower_output.txt", sep=""),
                  row.names=FALSE)
      write.table(tank.mcse, 
                  paste("./nBM_lower_ess/", args[1], "_", n.ess[i[2]], 
                        "_ess_new_lower_mcse.txt", sep=""),
                row.names=FALSE)
      write.table(sqrt(tank.std[3,]/(n-1)), 
                  paste("./nBM_lower_ess/", args[1], "_", n.ess[i[2]],
                        "_ess_new_lower_sd.txt", sep=""),
                row.names=FALSE)
      write.table(tank.mean, 
                paste("./nBM_lower_ess/", args[1], "_", n.ess[i[2]],
                      "_ess_new_lower_mean.txt", sep=""),
                row.names=FALSE)
      upper<- tank.mean+z*tank.mcse
      lower<- tank.mean-z*tank.mcse
      # coverage of the resulting confidence interval (0/1)
      prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper),ncol=1)*1
      write.table(prob.coverage, 
                  paste("./nBM_lower_ess/", args[1], "_", n.ess[i[2]], 
                        "_ess_new_lower_probcover.txt", sep=""), 
                  row.names=FALSE)
      # distance between estimates and truth
      prob.distance<- matrix(abs(tank.mean-truth[,1])>=eps[i[2]]*truth[,2],
                             ncol= 1)*1
      write.table(prob.distance, 
                  paste("./nBM_lower_ess/", args[1], "_", n.ess[i[2]], 
                        "_ess_new_lower_probdist.txt", sep=""), 
                  row.names=FALSE)
      write((proc.time()-ptm)[1:3], 
            paste("./nBM_lower_ess/", args[1], "_", n.ess[i[2]], 
                  "_ess_new_time.txt", sep=""))
      
      i[2]= i[2] + 1
    }
  }

  if(all(i == 4)){break}
}
rm(i, d, b.size, b.range, batch, b.count, batch.count, tank.std,
   tank, nbatch, tank.mcse, cond, tank.mean, ess.new, ess.app, 
   out, prob.coverage, prob.distance)
print("==============nBM_lower finished================")

## proposed upper bound batch means ##########################################
i<- matrix(1, ncol = 2) # indicator of the tuning parameter value
colnames(i)<- c("nBM_upper", "nBMess_upper")

d<- dim(samples)[2]
b.size<- c(2^7, 2^7)
b.range<- 2^seq(7, 15)
batch<- matrix(, b.size[1], d)
b.count<- 0 # number of iterations within batch
batch.count<- 0 # number of batches
n<- 0 # number of iterations
tank.std<- matrix(0, 3, d)

b.count= b.count+1
n= n+1
batch[b.count,]<- samples[n,]
tank.std[1,]<- samples[1,]

#### record present time ####
ptm<- proc.time()

while(1){
  b.count= b.count+1
  n= n+1
  batch[b.count,]<- samples[n,]
  ## recursive standard deviation estimation
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
  if((n>=2^14) && (batch.count>=20) && (!(dim(tank)[1]%%2))){
    batch.count<- 1
    b.size[2]<- b.size[1]
    b.size[1]<- 2^min(which(sqrt(n)<=b.range)+6) # upper bound
    nbatch<- floor(n/b.size[1])
    
    #check if batch size changes
    if(b.size[1] != b.size[2]){
      batch<- matrix(, b.size[1], d)
      tank<- t(sapply(split.data.frame(
        tank, rep(1:(dim(tank)[1]/2), 
                  times=rep(2, (dim(tank)[1]/2)))), apply, 2, mean))
    }
    
    tank.mcse<- apply(tank, 2, sd)*sqrt(b.size[1]/n)
    
    ## FWSR
    cond<- 2*z*tank.mcse-eps[i[1]]*sqrt(tank.std[3,]/(n-1))
    if(all(cond<=0)){
      tank.mean<- apply(tank, 2, mean)
      ess.new<- (tank.std[3,]/(n-1))/(tank.mcse^2)
      ess.app<- 4*(z/eps[i[1]])^2
      out<- list(n= n, app= ess.app, new= ess.new)
      write.table(out, 
                  paste("./nBM_upper/", args[1], "_", n.ess[i[1]],
                        "_new_upper_output.txt", sep=""), 
                  row.names=FALSE)
      write.table(tank.mcse, 
                  paste("./nBM_upper/", args[1], "_", n.ess[i[1]],
                        "_new_upper_mcse.txt", sep=""), 
                  row.names=FALSE)
      write.table(sqrt(tank.std[3,]/(n-1)), 
                  paste("./nBM_upper/", args[1], "_", n.ess[i[1]],
                        "_new_upper_sd.txt",sep=""), row.names=FALSE)
      write.table(tank.mean, paste("./nBM_upper/", args[1], "_", n.ess[i[1]], 
                                   "_new_upper_mean.txt", sep=""), 
                  row.names=FALSE)
      upper<- tank.mean+z*tank.mcse
      lower<- tank.mean-z*tank.mcse
      # coverage of the resulting confidence interval (0/1)
      prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper),ncol=1)*1
      write.table(prob.coverage, 
                  paste("./nBM_upper/", args[1], "_", n.ess[i[1]], 
                        "_new_upper_probcover.txt", sep=""), 
                  row.names=FALSE)
      # distance between estimates and truth
      prob.distance<- matrix(abs(tank.mean-truth[,1])>=eps[i[1]]*truth[,2],
                             ncol= 1)*1
      write.table(prob.distance, 
                  paste("./nBM_upper/", args[1], "_", n.ess[i[1]], 
                        "_new_upper_probdist.txt", sep=""), 
                  row.names=FALSE)
      write((proc.time()-ptm)[1:3], 
            paste("./nBM_upper/", args[1], "_", n.ess[i[1]], 
                  "_new_time.txt", sep=""))
      
      i[1]= i[1] + 1
    }
    
    ## ESS
    ess.new<- (tank.std[3,]/(n-1))/(tank.mcse^2)
    if(all(ess.new>=n.ess[i[2]])){
      tank.mean<- apply(tank, 2, mean)
      ess.app<- 4*(z/eps[i[2]])^2
      out<- list(n= n, app= ess.app, new= ess.new)
      write.table(out, paste("./nBM_upper_ess/", args[1], "_", n.ess[i[2]], 
                             "_ess_new_upper_output.txt", sep=""),
                  row.names=FALSE)
      write.table(tank.mcse, 
                  paste("./nBM_upper_ess/", args[1], "_", n.ess[i[2]], 
                        "_ess_new_upper_mcse.txt", sep=""),
                  row.names=FALSE)
      write.table(sqrt(tank.std[3,]/(n-1)), 
                  paste("./nBM_upper_ess/", args[1], "_", n.ess[i[2]],
                        "_ess_new_upper_sd.txt", sep=""),
                  row.names=FALSE)
      write.table(tank.mean, 
                  paste("./nBM_upper_ess/", args[1], "_", n.ess[i[2]],
                        "_ess_new_upper_mean.txt", sep=""),
                  row.names=FALSE)
      upper<- tank.mean+z*tank.mcse
      lower<- tank.mean-z*tank.mcse
      # coverage of the resulting confidence interval (0/1)
      prob.coverage<- matrix((truth[,1]>=lower)&(truth[,1]<=upper),ncol=1)*1
      write.table(prob.coverage, 
                  paste("./nBM_upper_ess/", args[1], "_", n.ess[i[2]], 
                        "_ess_new_upper_probcover.txt", sep=""), 
                  row.names=FALSE)
      # distance between estimates and truth
      prob.distance<- matrix(abs(tank.mean-truth[,1])>=eps[i[2]]*truth[,2],
                             ncol= 1)*1
      write.table(prob.distance, 
                  paste("./nBM_upper_ess/", args[1], "_", n.ess[i[2]], 
                        "_ess_new_upper_probdist.txt", sep=""), 
                  row.names=FALSE)
      write((proc.time()-ptm)[1:3], 
            paste("./nBM_upper_ess/", args[1], n.ess[i[2]], 
                  "_ess_new_time.txt", sep=""))
      
      i[2]= i[2] + 1
    }
  }
  
  if(all(i == 4)){break}
}
rm(i, d, b.size, b.range, batch, b.count, batch.count, tank.std,
   tank, nbatch, tank.mcse, cond, tank.mean, ess.new, ess.app, 
   out, prob.coverage, prob.distance)
print("==============nBM_upper finished================")


