#-----------------------------------------------------------------------------
# This program summarize the outputs from multiple parallel runs of the 
# spBayes_stoppingrules simulations. 
#
# input
# ------
# args[1]: the totol number of parallel runs 
#
# args[2]: a string indicating the method of interest
#
# output
# ------
# args[1]_args[2]_prob.txt: sequence of parameter-wise coverage probabilities
# and probability of whether the distance of estimates from truth is larger 
# than the threshold
#
# args[1]_args[2]_output.txt: summarized output includes average chain length,
# maximum memory usage, average coverage probability and max Probability.
#-----------------------------------------------------------------------------

#### inputs ####
args<-commandArgs(TRUE)
n.parallel<- as.numeric(args[1])
method<- args[2]

#### summarize information ####
if(length(grep(method, 'geweke'))){
  prob.cover<- list()
  run.length<- list()
  for(i in 1:n.parallel){
    prob.cover[length(prob.cover)+1]<- read.table(
      paste(i, method, "probcover.txt", sep="_"), 
      header = TRUE, row.names = NULL)
    run.length[length(run.length)+1]<- read.table(
      paste(i, method, "output.txt", sep="_"), 
      header = TRUE, row.names = NULL)[1]
  }
  
  ## estimate parameters
  prob.cover<- simplify2array(prob.cover)
  prob.cover<- apply(prob.cover, 1, mean)
  write.table(prob.cover, 
              paste(method, "prob.txt", sep="_"), 
              col.names=c("coverage"),
              row.names=FALSE)
  run.length<- simplify2array(run.length)
  sd.length<- sd(run.length)
  avg.length<- mean(run.length)
#   run.memory<- apply(run.length, 1, 
#                      function(x){object.size(matrix(runif(186*x), 186, x))})
# due to a bug
  run.memory<- apply(array(15000, n.parallel), 1, 
                     function(x){object.size(matrix(runif(186*x), 186, x))})
  avg.memory<- mean(run.memory)
  sd.memory<- sd(run.memory)
  avg.cover<- mean(prob.cover)
  write.table(matrix(c(avg.length, sd.length, avg.memory, sd.memory, avg.cover, 0), nrow=2), 
              paste(method, "output.txt", sep="_"),
              col.names=c("length", "memory", "coverage"),
              row.names=FALSE)
}else{
  ## memory usage calculation
  b.size<- 2^seq(7, 15)
  cal_memory_lower<- function(x){
    b<-2^max(which(sqrt(x)>=b.size)+6)
    a<- floor(x/b)
    return(object.size(matrix(runif(186*a),186,a)))
  }
  cal_memory_upper<- function(x){
    b<-2^min(which(sqrt(x)<=b.size)+6)
    a<- floor(x/b)
    return(object.size(matrix(runif(186*a),186,a)))
  }
  
  n.ess<- c(1000, 2000, 4000)
  for(n in n.ess){
    prob.cover<- list()
    prob.dist<- list()
    run.length<- list()
    for(i in 1:n.parallel){
      prob.cover[length(prob.cover)+1]<- read.table(
        paste(i, n, method, "probcover.txt", sep="_"), 
        header = TRUE, row.names = NULL)
      prob.dist[length(prob.dist)+1]<- read.table(
        paste(i, n, method, "probdist.txt", sep="_"), 
        header = TRUE, row.names = NULL)
      run.length[length(run.length)+1]<- read.table(
        paste(i, n, method, "output.txt", sep="_"), 
        header = TRUE, row.names = NULL)[1]
    }
    
    ## estimate parameters
    prob.cover<- simplify2array(prob.cover)
    prob.cover<- apply(prob.cover, 1, mean)
    prob.dist<- simplify2array(prob.dist)
    prob.dist<- apply(prob.dist, 1, mean)
    write.table(cbind(prob.cover, prob.dist), 
                paste(method, n, "prob.txt", sep="_"), 
                col.names=c("coverage", "distance"),row.names=FALSE)
    run.length<- simplify2array(run.length)
    sd.length<- sd(run.length)
    avg.length<- mean(run.length)
    # calculate max memory
    if(length(grep(method, c('standard', 'ess_standard')))){
      run.memory<- apply(run.length, 2, 
                         function(x){object.size(matrix(runif(186*x),186,x))})
    }else if(length(grep(method, c('new_lower', 'ess_new_lower')))){
      run.memory<- apply(run.length, 2, cal_memory_lower)      
    }else if(length(grep(method, c('new_upper', 'ess_new_upper')))){
      run.memory<- apply(run.length, 2, cal_memory_upper)      
    }else{print("WRONG method!!")}
    avg.memory<- mean(run.memory)
    sd.memory<- sd(run.memory)
    avg.cover<- mean(prob.cover)
    max.dist<- max(prob.dist)
    write.table(matrix(c(avg.length, sd.length, avg.memory, sd.memory, 
                         avg.cover, 0, max.dist, 0), nrow=2), 
                paste(method, n, "output.txt", sep="_"),
                col.names=c("length", "memory", "coverage", "distance"),
                row.names=FALSE)
  }
}