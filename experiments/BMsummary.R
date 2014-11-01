#-----------------------------------------------------------------------------
# This program summarize the outputs from multiple parallel runs of the 
# spBayes simulations with stopping rule implemented. 
#
# input
# ------
# args[1]: the totol number of parallel runs 
#
# args[2]: a string ["new" or "standard"] indicating the source of results
#
# args[3]: effective sample size, i.e. n.ess
#
# output
# ------
# args[1]_args[2]_args[3]_prob.txt: sequence of parameter-wise coverage 
# probabilities and probability of whether the distance of estimates from 
# truth is larger than the threshold
#
# args[1]_args[2]_args[3]_output.txt: summarized output includes average chain
# length, average run time, average coverage probability and max Probability.
#-----------------------------------------------------------------------------

#### inputs ####
args<-commandArgs(TRUE)
n.parallel<- as.numeric(args[1])
met<- args[2]
n.ess<- as.numeric(args[3])

#### read files ####
prob.cover<- list()
prob.dist<- list()
avg.time<- list()
avg.length<- list()
for(i in 1:n.parallel){
  prob.cover[length(prob.cover)+1]<- read.table(
    paste(n.ess, "_", i, "_", "probcover_", met, ".txt", sep=""), 
    header = TRUE, row.names = NULL)
  prob.dist[length(prob.dist)+1]<- read.table(
    paste(n.ess, "_", i, "_", "probdist_", met, ".txt", sep=""), 
    header = TRUE, row.names = NULL)
  avg.length[length(avg.length)+1]<- read.table(
    paste(n.ess, "_", i, "_", "output_", met, ".txt", sep=""), 
    header = TRUE, row.names = NULL)[1]
  avg.time[length(avg.time)+1]<- read.table(
    paste(n.ess, "_", i, "_", "time_", met, ".txt", sep=""), 
    row.names = NULL)[1]
}

#### estimate parameters ####
prob.cover<- simplify2array(prob.cover)
prob.dist<- simplify2array(prob.dist)
prob.cover<- apply(prob.cover, 1, mean)
prob.dist<- apply(prob.dist, 1, mean)
write.table(cbind(prob.cover, prob.dist), 
            paste(n.ess, met, "prob.txt", sep="_"), 
            col.names=c("coverage", "distance"),
            row.names=FALSE)
avg.length<- simplify2array(avg.length)
avg.time<- simplify2array(avg.time)
avg.length<- mean(avg.length)
avg.time<- mean(avg.time)
avg.cover<- mean(prob.cover)
max.dist<- max(prob.dist)
write.table(matrix(c(avg.length, avg.time, avg.cover, max.dist), nrow=1), 
            paste(n.ess, met, "output.txt", sep="_"),
            col.names=c("length", "time", "coverage", "maxProb"),
            row.names=FALSE)
