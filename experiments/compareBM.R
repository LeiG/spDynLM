#-----------------------------------------------------------------------------
# An experiment to compare the performance of standard non-overlap batch means
# and proposed non-overlap batch means estimation.
#-----------------------------------------------------------------------------

#### define batch means estimator ####
newBM<- function(x){
  n= length(x)
  b.range= 2^seq(1, 15)
#   b= 2^max(which(sqrt(n)>=b.range)) # lower bound
  b= 2^(max(which(sqrt(n)>=b.range))+1) # upper bound
  a= floor(n/b)
  
  y = sapply(1:a, function(k) return(mean(x[((k - 1) * b + 1):(k * b)])))
  mu.hat = mean(y)
  var.hat = b * sum((y - mu.hat)^2)/(a - 1)
  se = sqrt(var.hat/n)
  return(se)
}

standardBM<- function(x){
  n= length(x)
  b= floor(sqrt(n))
  a= floor(n/b)
  
  y = sapply(1:a, function(k) return(mean(x[((k - 1) * b + 1):(k * b)])))
  mu.hat = mean(y)
  var.hat = b * sum((y - mu.hat)^2)/(a - 1)
  se = sqrt(var.hat/n)
  return(se)
}

#### compare mcse estiamtes ####
plot.mcse<- function(x){
  x.new<- matrix(, ncol=length(x)-100)
  x.standard<- matrix(, ncol=length(x)-100)
  for(i in 1:(length(x)-100)){
    x.new[i]<- newBM(x[1:(i+100)])
    x.standard[i]<- standardBM(x[1:(i+100)])
  }
  plot(101:length(x), x.new, "l")
  lines(101:length(x), x.standard, col=2)
}

#### generate MCMC samples ####
x.indep<- rnorm(10000)
x.ar<- arima.sim(n = 10000, list(order =c(1,0,0), ar=0.5))

#### plots ####
plot.mcse(x.indep)
plot.mcse(x.ar)
