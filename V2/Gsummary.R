#-----------------------------------------------------------------------------
# This program summarize the outputs from multiple parallel runs of the 
# spBayes simulations with Geweke diagnostic implemented. 
#
# input
# ------
# args[1]: the totol number of parallel runs 
#
# args[2]: a string ["geweke"] indicating the source of results
#
# output
# ------
# args[1]_args[2]_prob.txt: sequence of parameter-wise coverage probabilities.
#
# args[1]_args[2]_output.txt: summarized output includes average chain length,
# average run time, average coverage probability.
#-----------------------------------------------------------------------------

#### inputs ####
args<-commandArgs(TRUE)
n.parallel<- as.numeric(args[1])
met<- args[2]

#### read files ####
prob.cover<- list()
avg.time<- list()
avg.length<- list()
for(i in 1:n.parallel){
  prob.cover[length(prob.cover)+1]<- read.table(
    paste(i, "_", "probcover_", met, ".txt", sep=""), 
    header = TRUE, row.names = NULL)
  avg.length[length(avg.length)+1]<- read.table(
    paste(i, "_", "output_", met, ".txt", sep=""), 
    header = TRUE, row.names = NULL)[1]
  avg.time[length(avg.time)+1]<- read.table(
    paste(i, "_", "time_", met, ".txt", sep=""), 
    row.names = NULL)[1]
}

#### estimate parameters ####
prob.cover<- simplify2array(prob.cover)
prob.cover<- apply(prob.cover, 1, mean)
write.table(prob.cover, 
            paste(met, "prob.txt", sep="_"), 
            col.names=c("coverage"),
            row.names=FALSE)
avg.length<- simplify2array(avg.length)
avg.time<- simplify2array(avg.time)
avg.length<- mean(avg.length)
avg.time<- mean(avg.time)
avg.cover<- mean(prob.cover)
write.table(matrix(c(avg.length, avg.time, avg.cover), nrow=1), 
            paste(met, "output.txt", sep="_"),
            col.names=c("length", "time", "coverage"),
            row.names=FALSE)