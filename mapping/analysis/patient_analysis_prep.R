library(vegan)
library(ggplot2)
library(igraph)
setwd('P:/GitRep/BacArena')
source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

# cats = read.csv("P:/MAPPING/Model_data/class773.csv",row.names=1)
# all = as.matrix(read.csv("P:/MAPPING/Model_data/alldat773.csv", header=T, row.names=1))
# alldiet = read.csv("P:/MAPPING/Model_data/diets.csv", header=T, row.names=1)
# pclass = read.csv("P:/MAPPING/Model_data/class773.csv", header=T, row.names=1)
# gradmet = c(rownames(alldiet)[which(alldiet$type=="Mucin")],"EX_o2(e)","EX_no(e)")

setwd("P:/MAPPING/Pediatric_crohns")
alldiet = read.csv("P:/MAPPING/Pediatric_crohns/diet_unhealthy.csv", header=T, row.names=1)
all_filter = as.matrix(read.csv("P:/MAPPING/Pediatric_crohns/alldat773_pediatric_filtered.csv", header=T, row.names=1))
cats = read.csv("P:/MAPPING/Pediatric_crohns/sra_result_just_nutrition.csv",row.names=2)
catsel = data.frame(Reads=gsub(" S","S",rownames(cats)),Sample.Name=gsub(" S","S",rownames(cats)),category=cats$Disease)
catsel = cbind(catsel,cats)
catsel$category = as.character(catsel$category)
catsel$category[which(catsel$category == "Control")] = "Healthy"
catsel$category[which(catsel$category == "Crohn")] = "Crohn's disease"
catsel$category = as.factor(catsel$category)
pclass = data.frame(V2=catsel$category,V3=catsel$Reads)
rownames(pclass) = as.character(pclass$V3)


#setwd("P:/MAPPING/CommunityModels773/FINAL/filtered/mix/")
setwd("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/FINAL/")
tit = 25

######################### Are there difference in abundance and concentration between healthy and disease
########################################### generating the data

res = list.files()
res = res[grep(".RData",res)]
res = res[grep("Result",res)]

sublist = list()
poplist = list()
for(i in res){
  print(i)
  load(i)
  if(class(result)=="Eval"){result = list(result)}
  patnam = gsub("Result_","",i)
  patnam = gsub(".RData","",patnam)
  sublist[[patnam]] = list()
  for(j in 1:length(result)){
    dat = plotCurves(result[[j]], retdata=T, graph=F)
    sublist[[patnam]][[j]] = dat$Substances
    poplist[[patnam]][[j]] = dat$Population
  }
  rm(result)
  gc()
}
names(sublist) = gsub("_","",names(sublist))
names(poplist) = gsub("_","",names(poplist))
save(sublist,file="sublist_eu.RData")
save(poplist,file="poplist_eu.RData")

######################### Are there difference in abundance and concentration between healthy and disease
########################################### generating the data

# res = list.files()
# res = res[grep(".RData",res)]
# res = res[grep("Result",res)]
# 
# poplist = list()
# sublist = list()
# gralist = list()
# pops = matrix(0,nrow=nrow(all),ncol=length(res),
#               dimnames=list(rownames(all),gsub("Result_","",res)))
# colnames(pops) = gsub(".RData","",colnames(pops))
# subs = matrix(0,nrow=nrow(alldiet),ncol=length(res),
#               dimnames=list(rownames(alldiet),colnames(pops)))
# for(i in res){
#   print(i)
#   load(i)
#   if(class(result)=="Eval"){result = list(result)}
#   patnam = gsub("Result_","",i)
#   patnam = gsub(".RData","",patnam)
#   pop = matrix(0,nrow=length(result[[1]]@specs),ncol=length(result),
#                dimnames=list(names(result[[1]]@specs),1:3))
#   sub = matrix(0,nrow=length(result[[1]]@mediac),ncol=length(result),
#                dimnames=list(unname(result[[1]]@mediac),1:3))
#   gradk = list()
#   for(k in 1:length(result[[1]]@simlist)){
#     gradj = list()
#     for(j in 1:length(result)){
#       dat = plotCurves(result[[j]], retdata=T, graph=F)
#       pop[rownames(dat$Population),j] = dat$Population[,k]
#       sub[rownames(dat$Substances),j] = dat$Substances[,k]
#       gradj[[j]] = result[[j]]@medlist[[k]][gradmet]
#     }
#     gradk[[k]] = gradj
#     pops[names(apply(pop,1,mean)),patnam] = apply(pop,1,mean)
#     subs[names(apply(sub,1,mean)),patnam] = apply(sub,1,mean)
#     poplist[[k]] = pops
#     sublist[[k]] = subs
#   }
#   gralist[[i]] = gradk
#   rm(result)
#   gc()
# }
# names(gralist) = gsub(".RData","",gsub("Result_","",res))
# 
# save(gralist,file="gralist_eu.RData")
# names(gralist) = gsub("_","",names(gralist))
# #load("result_lists_eu.RData")

#########################
###########################################

res = list.files()
res = res[grep(".RData",res)]
res = res[grep("Result",res)]

fluxlist = list()
for(i in res){
  print(i)
  load(i)
  if(class(result)=="Eval"){result = list(result)}
  patnam = gsub("Result_","",i)
  patnam = gsub(".RData","",patnam)
  fluxlist[[patnam]] = list()
  for(j in 1:length(result)){
    fluxlist[[patnam]][[j]] = result[[j]]@mfluxlist
  }
  rm(result)
  gc()
}
save(fluxlist,file="fluxlist_eu.RData")

fluxlist_red = list()
for(i in 1:length(fluxlist)){
  print(i)
  flj = list()
  for(j in names(fluxlist[[i]][[1]][[tit]])){
    bacflux = fluxlist[[i]][[1]][[tit]][[j]]
    # for(k in 2:length(fluxlist[[i]])){
    #   bacflux = bacflux + fluxlist[[i]][[k]][[tit]][[j]]
    # }
    # bacflux = bacflux/length(fluxlist[[i]])
    flj[[j]] = bacflux
  }
  fluxlist_red[[i]] = flj
}
names(fluxlist_red) = names(fluxlist)
names(fluxlist_red) = gsub("_","",names(fluxlist_red))
save(fluxlist_red,file="fluxlist_red_eu.RData")


######################### get mean value of all substances for all positions
##########################################

# res = list.files()
# res = res[grep(".RData",res)]
# res = res[grep("Result",res)]
# 
# subgrid = matrix(0, nrow=length(res), ncol=(nrow(alldiet)),
#                  dimnames=list(gsub(".RData","",gsub("Result_","",res)),rownames(alldiet)))
# fluxlist = list()
# for(i in res){
#   print(i)
#   load(i)
#   if(class(result)=="Eval"){result = list(result)}
#   patnam = gsub("Result_","",i)
#   patnam = gsub(".RData","",patnam)
#   subgridit = matrix(0, nrow=length(result),ncol=ncol(subgrid))
#   colnames(subgridit) = colnames(subgrid)
#   for(j in 1:length(result)){
#     arenait = getArena(result[[j]],tit-1) 
#     subgridit[j,names(lapply(arenait@media,name))] = unlist(lapply(arenait@media,function(x){mean(as.matrix(x@diffmat))}))
#   }
#   subgrid[patnam,] = apply(subgridit,2,mean)
#   rm(result)
#   gc()
# }
# 
# write.csv(subgrid,file="subgrid.csv")

######################### reducing the size of the gradient lists
##################################################################

# gralist_reduce = list()
# tit = 13
# for(i in 1:length(gralist)){
#   itgra = gralist[[i]][[tit]][[1]]
#   for(j in 2:length(gralist[[i]][[tit]])){
#     for(k in 1:length(gralist[[i]][[tit]][[j]])){
#       itgra[[k]] = itgra[[k]] + gralist[[i]][[tit]][[j]][[k]]
#     }
#   }
#   for(j in 1:length(itgra)){
#     itgra[[j]] = itgra[[j]]/length(gralist[[i]][[tit]])
#   }
#   gralist_reduce[[i]] = itgra
# }
# names(gralist_reduce) = names(gralist)
# gralist_init = list()
# tit = 1
# for(i in 1:length(gralist)){
#   itgra = gralist[[i]][[tit]][[1]]
#   for(j in 2:length(gralist[[i]][[tit]])){
#     for(k in 1:length(gralist[[i]][[tit]][[j]])){
#       itgra[[k]] = itgra[[k]] + gralist[[i]][[tit]][[j]][[k]]
#     }
#   }
#   for(j in 1:length(itgra)){
#     itgra[[j]] = itgra[[j]]/length(gralist[[i]][[tit]])
#   }
#   gralist_init[[i]] = itgra
# }
# names(gralist_init) = names(gralist)
# save(gralist_init,file="gralist_init_eu.RData")
# save(gralist_reduce,file="gralist_reduce_eu.RData")
