## apply stopping rule with ALTERNATIVE bm
library(mcmcse)

## mcmc samples
beta.0.samples<- read.table("beta.0.samples.txt")
beta.samples<- read.table("beta.samples.txt")
sigma.eta.samples<- read.table("sigma.eta.samples.txt")
theta.samples<- read.table("theta.samples.txt")
samples<- cbind(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples)

## stopping rule
eps<- 0.05
jump<- 20
check<- 2^12
z<- 1.96
d<- dim(samples)[2]
b<- matrix(, 20, d)
X<- y = sapply(1:a, function(k) return(mean(g(x[((k - 1) * 
                                                   b + 1):(k * b)]))))
X<- matrix(, check/sqrt(check), d)
sapply(0:b, function(j) return(mean((g(x[1:(n - j)]) - 
                                       mu.hat) * (g(x[(j + 1):n]) - mu.hat))))
while(1){
  b.size<- sqrt(check)
  nbatch<- check/b.size
  
  
}