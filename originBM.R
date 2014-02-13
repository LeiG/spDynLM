## apply stopping rule with ORIGINAL bm
library(mcmcse)

## mcmc samples
beta.0.samples<- read.table("beta.0.samples.txt")
beta.samples<- read.table("beta.samples.txt")
sigma.eta.samples<- read.table("sigma.eta.samples.txt")
theta.samples<- read.table("theta.samples.txt")
samples<- cbind(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples)

## define functions
bm<- function(x){
  n= length(x)
  b= floor(sqrt(n))
  a= floor(n/b)
  
  y = sapply(1:a, function(k) return(mean(x[k:(k + b - 1)])))
  mu.hat = mean(x)
  var.hat = n * b * sum((y - mu.hat)^2)/(a - 1)/a
  se = sqrt(var.hat/n)
  return(c(mu.hat, se))
}

## stopping rule
eps<- 0.05
jump<- 20
check<- 2^12
z<- 1.96
while(1){
  X<- samples[0:check,]
  print(check)
  
  mcse<- apply(X, 2, bm)
  std<- apply(X, 2, sd)
  cond<- 2*z*mcse[2,]/sqrt(check)-eps*std #95% confidence
  if(all(cond<=0)){
    ess.old<- apply(X, 2, ess)
    ess.new<- std^2/mcse^2
    ess.app<- 4*(z/eps)^2
    out<- list(n= check, app= ess.app, old= ess.old, new= ess.new)
    write.table(out, "output_origin.txt")
    write.table(mcse, "mcse_origin.txt")
    break
  }
  
  b.size<- sqrt(check)
  check<- check+jump*b.size
  nbatch<- check/b.size
  if(!(nbatch%%2)){check<- check+b.size}
}