## spatial Bayesian model
require(spBayes)

## manipulate raw data
data("NETemp.dat")
ne.temp <- NETemp.dat
set.seed(3)
#take a chunk of New England
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
n.samples <- 2000

## get MCMC samples
m.1 <- spDynLM(mods, data=cbind(y.t,ne.temp[,"elev",drop=FALSE]), coords=coords,
               starting=starting, tuning=tuning, priors=priors, get.fitted =TRUE,
               cov.model="exponential", n.samples=n.samples, n.report=0.1*n.samples)
burn.in <- floor(0.75*n.samples)
quant <- function(x){quantile(x, prob=c(0.5, 0.025, 0.975))}
beta <- apply(m.1$p.beta.samples[burn.in:n.samples,], 2, quant)
beta.0 <- beta[,grep("Intercept", colnames(beta))]
beta.1 <- beta[,grep("elev", colnames(beta))]
plot(m.1$p.beta.0.samples)
par(mfrow=c(2,1))
plot(1:N.t, beta.0[1,], pch=19, cex=0.5, xlab="months", ylab="beta.0", ylim=range(beta.0))
arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[3,], length=0.02, angle=90)
arrows(1:N.t, beta.0[1,], 1:N.t, beta.0[2,], length=0.02, angle=90)
plot(1:N.t, beta.1[1,], pch=19, cex=0.5, xlab="months", ylab="beta.1", ylim=range(beta.1))
arrows(1:N.t, beta.1[1,], 1:N.t, beta.1[3,], length=0.02, angle=90)
arrows(1:N.t, beta.1[1,], 1:N.t, beta.1[2,], length=0.02, angle=90)
theta <- apply(m.1$p.theta.samples[burn.in:n.samples,], 2, quant)
sigma.sq <- theta[,grep("sigma.sq", colnames(theta))]
tau.sq <- theta[,grep("tau.sq", colnames(theta))]
phi <- theta[,grep("phi", colnames(theta))]
par(mfrow=c(3,1))
plot(1:N.t, sigma.sq[1,], pch=19, cex=0.5, xlab="months", ylab="sigma.sq", ylim=range(sigma.sq))
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[3,], length=0.02, angle=90)
arrows(1:N.t, sigma.sq[1,], 1:N.t, sigma.sq[2,], length=0.02, angle=90)
plot(1:N.t, tau.sq[1,], pch=19, cex=0.5, xlab="months", ylab="tau.sq", ylim=range(tau.sq))
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[3,], length=0.02, angle=90)
arrows(1:N.t, tau.sq[1,], 1:N.t, tau.sq[2,], length=0.02, angle=90)
plot(1:N.t, 3/phi[1,], pch=19, cex=0.5, xlab="months", ylab="eff. range (km)", ylim=range(3/phi))
arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[3,], length=0.02, angle=90)
arrows(1:N.t, 3/phi[1,], 1:N.t, 3/phi[2,], length=0.02, angle=90)
y.hat <- apply(m.1$p.y.samples[,burn.in:n.samples], 1, quant)
y.hat.med <- matrix(y.hat[1,], ncol=N.t)
y.hat.up <- matrix(y.hat[3,], ncol=N.t)
y.hat.low <- matrix(y.hat[2,], ncol=N.t)
y.obs <- as.vector(as.matrix(y.t[-holdout.station.id, -miss]))
y.obs.hat.med <- as.vector(y.hat.med[-holdout.station.id, -miss])
y.obs.hat.up <- as.vector(y.hat.up[-holdout.station.id, -miss])
y.obs.hat.low <- as.vector(y.hat.low[-holdout.station.id, -miss])
y.ho <- as.matrix(y.t.holdout)
y.ho.hat.med <- as.vector(y.hat.med[holdout.station.id, miss])
y.ho.hat.up <- as.vector(y.hat.up[holdout.station.id, miss])
y.ho.hat.low <- as.vector(y.hat.low[holdout.station.id, miss])
par(mfrow=c(2,1))
plot(y.obs, y.obs.hat.med, pch=19, cex=0.5, xlab="observed",
     ylab="fitted", main="Observed vs. fitted")
arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.up, length=0.02, angle=90)32 spGLM
arrows(y.obs, y.obs.hat.med, y.obs, y.obs.hat.low, length=0.02, angle=90)
lines(-50:50, -50:50, col="blue")
plot(y.ho, y.ho.hat.med, pch=19, cex=0.5, xlab="observed",
     ylab="predicted", main="Observed vs. predicted")
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.up, length=0.02, angle=90)
arrows(y.ho, y.ho.hat.med, y.ho, y.ho.hat.low, length=0.02, angle=90)
lines(-50:50, -50:50, col="blue")

