## comparative study

## input files
mcse.alter.05<- read.table("0.05_mcse_alternative.txt")
mcse.orig.05<- read.table("0.05_mcse_origin.txt")
mcse.orig.1<- read.table("0.1_mcse_origin.txt")
mean.alter.05<- read.table("0.05_mean_alternative.txt")
mean.orig.05<- read.table("0.05_mean_origin.txt")
mean.orig.1<- read.table("0.1_mean_origin.txt")
sd.alter.05<- read.table("0.05_sd_alternative.txt")
sd.orig.05<- read.table("0.05_sd_origin.txt")
sd.orig.1<- read.table("0.1_sd_origin.txt")
out.alter.05<- read.table("0.05_output_alternative.txt")
out.orig.05<- read.table("0.05_output_origin.txt")
out.orig.1<- read.table("0.1_output_origin.txt")

## paramters
n.alter.05<- out.alter.05[1,1]
n.orig.05<- out.orig.05[1,1]
n.orig.1<- out.orig.1[1,1]

## ratios
ratio1<- unlist((mcse.alter.05*sqrt(n.alter.05))/(mcse.orig.05*sqrt(n.orig.05)))
ratio2<- unlist((mcse.orig.1*sqrt(n.orig.1))/(mcse.orig.05*sqrt(n.orig.05)))
hist(ratio1)
hist(ratio2)

ratio3<- unlist(mean.alter.05/mean.orig.05)
ratio4<- unlist(mean.orig.1/mean.orig.05)
hist(ratio3)
hist(ratio4)

ratio5<- unlist(sd.alter.05/sd.orig.05)
ratio6<- unlist(sd.orig.1/sd.orig.05)
hist(ratio5)
hist(ratio6)

ratio7<- unlist(2*1.96*mcse.alter.05/sd.alter.05)
ratio8<- unlist(2*1.96*mcse.orig.05/sd.orig.05)
ratio9<- unlist(2*1.96*mcse.orig.1/sd.orig.1)
hist(ratio7)
hist(ratio8)
hist(ratio9)