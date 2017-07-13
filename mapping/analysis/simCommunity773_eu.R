library(parallel)

setwd("/work/users/ebauer/CommunityModels773/FINAL/filtered")
rdats = list.files()
rdats = rdats[grep(".RData",rdats)]
rmdat <- grep("Result",rdats)
if(length(rmdat)!=0){rdats <- rdats[-rmdat]}

cl <- makeCluster(length(rdats), type="PSOCK") # PSOCK works with win/mac/lin

print(system.time(reps <- parLapply(cl, rdats, function(x){
  library(sybil)
  library(Rcpp)
  library(RcppArmadillo)
  library(ReacTran)
  source(file="/work/users/ebauer/BacArena/R/Arena.R")
  source(file="/work/users/ebauer/BacArena/R/Stuff.R")
  source(file="/work/users/ebauer/BacArena/R/Substance.R")
  source(file="/work/users/ebauer/BacArena/R/Organism.R")
  source(file="/work/users/ebauer/CommunityModels773/reInitProb.R")
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
  
  load(x)
  result = list()
  for(i in 1:3){
    result[[i]] = simEnv(reInitProb(patient[[i]]), time=12)
  }
  save(result,file=paste("Result_",x,sep="_"))
  return("")
}) ))
