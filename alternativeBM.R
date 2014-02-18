## apply stopping rule with ALTERNATIVE bm
require(mcmcse)

## mcmc samples
beta.0.samples<- read.table("beta.0.samples.txt")
beta.samples<- read.table("beta.samples.txt")
sigma.eta.samples<- read.table("sigma.eta.samples.txt")
theta.samples<- read.table("theta.samples.txt")
samples<- as.matrix(cbind(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples))
rm(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples)

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
    
    tank.mcse<- sqrt(apply((tank-apply(tank, 2, mean))^2, 2, sum)*b.size[1]/((nbatch-1)*n))
    cond<- 2*z*tank.mcse-eps*sqrt(tank.std[3,]/(n-1))
    write.table(cond, "cond_alternative.txt")
    if(all(cond<=0)){
      ess.new<- (tank.std[3,]/(n-1))/tank.mcse^2*n
      ess.app<- 4*(z/eps)^2
      out<- list(n= check, app= ess.app, new= ess.new)
      write.table(out, "output_alternative.txt")
      write.table(tank.mcse, "mcse_alternative.txt")
      break
    }
  }
}