## apply Geweke convergence diagnostics
require(mcmcse)
require(boa)

## mcmc samples
beta.0.samples<- read.table("beta.0.samples.txt")
beta.samples<- read.table("beta.samples.txt")
sigma.eta.samples<- read.table("sigma.eta.samples.txt")
theta.samples<- read.table("theta.samples.txt")
samples<- as.matrix(cbind(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples))
rm(beta.0.samples, beta.samples, sigma.eta.samples, theta.samples)

## diagnostic
burnin<- 
n<- 
X<- samples[burnin:(burnin+n),]