####################### R input file ##################################
respath = "/mnt/gaiagpfs/users/workdirs/fbaldini/PD/allreads/Final" #path to result folder containing csv and bam/sam
#####################  ATable creator from bam table/ files #######################

library(parallel)

setwd(respath)
bams = list.files(respath)
bams = bams[grep('filter.sort.bam',bams)]
bams = bams[-grep('bai',bams)]

cl <- makeCluster(length(bams)/2, type="PSOCK") # PSOCK works with win/mac/lin

print(system.time(reps <- parLapply(cl, bams, function(x){
  library(Rsamtools)
  bam <- scanBam(x)[[1]] # the result comes in nested lists
  ind <- ! is.na(bam$pos)
  bam <- lapply(bam, function(x) x[ind])
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
  return(list(name=x,coverage=mean(coverage(ranges))))    
})))
#names(reps)=bams
save(reps,file="covlist.RData")
