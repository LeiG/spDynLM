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
# args[1]_args[2]_args[3]_cover.txt: sequence of parameter-wise coverage 
# probabilities
#
# args[1]_args[2]_args[3]_dist.txt: sequence of parameter-wise probability of
# whether the distance of estimates from truth is larger than the threshold
#-----------------------------------------------------------------------------

#### inputs ####
args<-commandArgs(TRUE)
n.parallel<- as.numeric(args[1])
met<- args[2]
n.ess<- as.numeric(args[3])

#### read files ####
prob.cover<- list()
prob.dist<- list()
for(i in 1:n.parallel){
  prob.cover[length(prob.cover)+1]<- read.table(
    paste(n.ess, "_", i, "_", "probcover_", met, ".txt", sep=""), 
    header = TRUE, row.names = NULL)
  prob.dist[length(prob.dist)+1]<- read.table(
    paste(n.ess, "_", i, "_", "probdist_", met, ".txt", sep=""), 
    header = TRUE, row.names = NULL)
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
