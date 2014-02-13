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
while(1){
  temp<- samples[]
  
}