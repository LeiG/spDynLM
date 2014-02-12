## apply stopping rule with ORIGINAL bm

## mcmc samples
beta.0.samples<- read.table("beta.0.samples.txt")
beta.samples<- read.table("beta.samples.txt")
sigma.eta.samples<- read.table("sigma.eta.samples.txt")
theta.samples<- read.table("theta.samples.txt")

## define functions
bm<- function(x){
  n= length(x)
  b= floor(sqrt(n))
  a= floor(n/b)
  
  y = sapply(1:a, function(k) return(mean(x[k:(k + b - 1)])))
  mu.hat = mean(x)
  var.hat = n * b * sum((y - mu.hat)^2)/(a - 1)/a
  se = sqrt(var.hat/n)
  list(est= mu.hat, se= se)
}

## stopping rule
epsilon<- 0.05
while(1){
  check<- 2^12
  b<- sqrt(check)
}
mcse.beta.0<- apply(beta.0.samples, 2, bm)
mcse.beta<- apply(beta.samples, 2, bm)
mcse.sigma.eta<- apply(sigma.eta.samples, 2, bm)
mcse.theta<- apply(theta.samples, 2, bm)