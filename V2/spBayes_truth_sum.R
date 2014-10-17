#-----------------------------------------------------------------------------
# This program takes the outputs from spBayes_truth_para.R (i.e. ) as inputs 
# and then summarize them using a modified variance calculation. It estimates
# the mean and standard deviation of the unknown parameters.
#
# input
# ------
# args: the totol number of parallel runs and the length of each run.
# 
# output
# -------
# truth.txt: the overall estimation of mean and standard deviation of the 
# parameters from the parallel runs.
#-----------------------------------------------------------------------------

#### total number of parallel runs ####
args<-commandArgs(TRUE)
n.parallel<- args[1]
n.samples<- args[2]

#### read files ####
sample.means<- list()
sample.vars<- list()
for(i in 1:n.parallel){
  sample<- read.table(i, "truth.txt", sep="_"), header = TRUE)
  sample.means[length(sample.means)+1]<- list(sample[,1])
  sample.vars[length(sample.vars)+1]<- list(sample[,2])
  rm(sample)
}

#### estimate parameters ####
sample.means<- simplify2array(sample.means)
sample.vars<- simplify2array(sample.vars)
param.mean<- apply(sample.means, 1, mean)
param.var<- (n.samples-1)/(n.parallel*n.samples-1)*apply(sample.vars, 1, sum) +
  n.samples/(n.parallel*n.samples-1)*apply((sample.means-pop.mean)^2, 1, sum)
sample.truth<- cbind(param.mean, sqrt(param.var))
write.table(sample.truth, "truth.txt", col.names=c("mean", "sd"),
            row.names=FALSE)


