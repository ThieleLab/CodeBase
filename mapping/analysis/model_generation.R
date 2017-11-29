##Loading the required packages## 
library(ReacTran)
library(sybil)
library(sybilSBML)
library(Rcpp)
library(RcppArmadillo)
library(sybil)
library(igraph)
library(parallel)
library(gridExtra)
library(ggplot2)
setwd('P:/GitRep/BacArena')

source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
source(file="R/Stuff.R")
Rcpp::sourceCpp("src/diff.cpp")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

# setwd("P:/Mapping/sbml/28_04")
# files = list.files(".", pattern = ".xml")
# specs = list()
# for(i in 1:length(files)){
#   print(i)
#   specs[[i]] = readSBMLmod(files[i])
# }
# names(specs) = unlist(lapply(specs,mod_desc))
# for(i in 1:length(specs)){
#   specs[[i]] = changeBounds(specs[[i]],specs[[i]]@react_id[grep('EX_',specs[[i]]@react_id)],lb=-1000)
#   objnam = specs[[i]]@react_id[grep("biomass",specs[[i]]@react_id)]
#   objnam = objnam[-which(objnam=="EX_biomass(e)")]
#   specs[[i]] = changeObjFunc(specs[[i]],objnam)
# }
# names(specs)[which(names(specs)=="Gemella_moribillum_M424")] = "Gemella_morbillorum_M424"
# save(specs,file="P:/MAPPING/Models773/agora773_constrain.RData")

alldiet = read.csv("P:/MAPPING/Model_data/diets_eu.csv", header=T, row.names=1)
newflux = read.csv("P:/MAPPING/Model_data/fluxes.csv",header=T,row.names = 6)
alldiet[intersect(rownames(alldiet),rownames(newflux)),"EU"] = newflux[intersect(rownames(alldiet),rownames(newflux)),"fluxValue"]

pclass = read.csv("P:/MAPPING/Model_data/class773.csv", header=F, row.names=1)
all_filter = as.matrix(read.csv("P:/MAPPING/Results773/bwa_human_filtered/filter1/alldat773_sel_filter1.csv", header=T, row.names=1))
all_unfilter = as.matrix(read.csv("P:/MAPPING/Results773/bwa_human_filtered/unfiltered/alldat773_sel_unfiltered.csv", header=T, row.names=1))
load("P:/MAPPING/Models773/agora773_constrain.RData")

#alldiet = read.csv("P:/MAPPING/Pediatric_crohns/diet_unhealthy.csv", header=T, row.names=1)
all_filter = as.matrix(read.csv("P:/MAPPING/Pediatric_crohns/alldat773_pediatric_dysbiotic_filtered.csv", header=T, row.names=1))
load("P:/MAPPING/Models773/agora773_constrain.RData")

rownames(all_filter)[which(rownames(all_filter)=="Gemella_moribillum_M424")] = "Gemella_morbillorum_M424"
#rownames(all_unfilter)[which(rownames(all_filter)=="Gemella_moribillum_M424")] = "Gemella_morbillorum_M424"
#all_filter = all_filter[,-(grep("R.",colnames(all_filter)))]
#all_unfilter = all_filter[,-(grep("R.",colnames(all_unfilter)))]


length(intersect(names(specs),rownames(all_filter)))/nrow(all_filter)

#############################################################################################################
#################################### Scripts to create defined models #######################################
#############################################################################################################

creIndMod = function(specs, abundance, diet, dsel, disease, maxind=500, replicates, econc, gconc, n=100, speed=5){
  vol = 2
  reps = list()
  for(i in 1:replicates){
    print(paste("Building replicate",i,sep=" "))
    arena <- Arena(n, n, seed=i, stir=F, Lx=0.025*(n/100), Ly=0.025*(n/100), tstep=1)
    for(j in names(abundance)){
      model = specs[[j]]
      if(ceiling(abundance[j]*maxind)>0){
        bac = Bac(model=model, growtype="exponential", speed=speed)
        arena = addOrg(arena, bac, amount=ceiling(abundance[j]*maxind))
      }
    }
    dietsel = diet[arena@mediac,]
    muc = rownames(dietsel)[which(dietsel$type == "Mucin")]
    arena = addSubs(arena,smax=dietsel[,dsel]/vol,difspeed=dietsel[,"diffconstant"],unit="mM",mediac=rownames(dietsel))
    arena = addSubs(arena,smax=econc,unit="mM",mediac=rownames(dietsel)[which(dietsel$class=="essential")],add=T)
    arena = createGradient(arena,smax=gconc,mediac=muc,position='left',steep=0.5,add=F,unit="mM")
    if(disease){
      if("EX_o2(e)" %in% arena@mediac){
        arena = createGradient(arena,smax=gconc,mediac="EX_o2(e)",position='left',steep=0.5,add=F,unit="mM")
      }
      if("EX_no(e)" %in% arena@mediac){
        arena = createGradient(arena,smax=gconc,mediac="EX_no(e)",position='left',steep=0.5,add=F,unit="mM")
      }
    }
    reps[[i]] = arena
  }
  return(reps)
}

creIndModMix = function(specs, abundance, diet, dsel, maxind=500, replicates, econc, n=100, speed=5, tstep=1){
  vol = 2
  reps = list()
  for(i in 1:replicates){
    print(paste("Building replicate",i,sep=" "))
    arena <- Arena(n, n, seed=i, stir=T, Lx=0.025*(n/100), Ly=0.025*(n/100), tstep=1)
    for(j in names(abundance)){
      model = specs[[j]]
      if(ceiling(abundance[j]*maxind)>0){
        bac = Bac(model=model, growtype="exponential")
        arena = addOrg(arena, bac, amount=ceiling(abundance[j]*maxind))
      }
    }
    dietsel = diet[arena@mediac,]
    arena = addSubs(arena,smax=dietsel[,dsel]/vol,unit="mM",mediac=rownames(dietsel))
    arena = addSubs(arena,smax=econc,unit="mM",mediac=rownames(dietsel)[which(dietsel$unclass=="essential")],add=T)
    reps[[i]] = arena
  }
  if(replicates == 1){reps = reps[[1]]}
  return(reps)
}

creIndModMixRich = function(specs, abundance, maxind=500, replicates, econc, n=100, speed=5, tstep=1){
  reps = list()
  for(i in 1:replicates){
    print(paste("Building replicate",i,sep=" "))
    arena <- Arena(n, n, seed=i, stir=T, Lx=0.025*(n/100), Ly=0.025*(n/100), tstep=1)
    for(j in names(abundance)){
      model = specs[[j]]
      if(ceiling(abundance[j]*maxind)>0){
        bac = Bac(model=model, growtype="exponential")
        arena = addOrg(arena, bac, amount=ceiling(abundance[j]*maxind))
      }
    }
    arena = addSubs(arena,smax=econc,unit="mM")
    reps[[i]] = arena
  }
  if(replicates == 1){reps = reps[[1]]}
  return(reps)
}

creIndModMixRichTreat = function(specs, abundance, maxind=500, replicates, econc, n=100, speed=5, tstep=1, treatmet,tconc){
  reps = list()
  for(i in 1:replicates){
    print(paste("Building replicate",i,sep=" "))
    arena <- Arena(n, n, seed=i, stir=T, Lx=0.025*(n/100), Ly=0.025*(n/100), tstep=1)
    for(j in names(abundance)){
      model = specs[[j]]
      if(ceiling(abundance[j]*maxind)>0){
        bac = Bac(model=model, growtype="exponential")
        arena = addOrg(arena, bac, amount=ceiling(abundance[j]*maxind))
      }
    }
    arena = addSubs(arena,smax=econc,unit="mM")
    arena = addSubs(arena,intersect(arena@mediac,treatmet),smax=tconc,unit="mM",add=F)
    reps[[i]] = arena
  }
  if(replicates == 1){reps = reps[[1]]}
  return(reps)
}

#############################################################################################################
#################################### create defined models ##################################################
#############################################################################################################

conc = 0.0001
mconc = 1
setwd("P:/MAPPING/CommunityModels773/FINAL/filtered")
all = all_filter
for(i in 1:ncol(all)){
  pid = colnames(all)[i]
  print(pid)
  if(pclass[pid,1]=="CD"){dis=T}else{dis=F}
  patient = creIndMod(specs, all[,i], alldiet, dsel="EU", disease=dis, maxind=300, 
                      replicates=3, econc=conc, gconc=mconc, n=200, speed=10)
  save(patient,file=paste(pid,"RData",sep="."))
  if(dis){
    patient = creIndMod(specs, all[,i], alldiet, dsel="EU", disease=F, maxind=300, 
                        replicates=3, econc=conc, gconc=mconc, n=200, speed=10)
    save(patient,file=paste(pid, "control","RData",sep="."))
  }
}

conc = 0.0001
mconc = 1
setwd("P:/MAPPING/CommunityModels773/FINAL/unfiltered")
all = all_unfilter
for(i in 1:ncol(all)){
  pid = colnames(all)[i]
  print(pid)
  if(pclass[pid,1]=="CD"){dis=T}else{dis=F}
  patient = creIndMod(specs, all[,i], alldiet, dsel="EU", disease=dis, maxind=500, 
                      replicates=3, econc=conc, gconc=mconc, n=200, speed=10)
  save(patient,file=paste(pid,"RData",sep="."))
  if(dis){
    patient = creIndMod(specs, all[,i], alldiet, dsel="EU", disease=F, maxind=500, 
                        replicates=3, econc=conc, gconc=mconc, n=200, speed=10)
    save(patient,file=paste(pid, "control","RData",sep="."))
  }
}

conc = 0.0001
mconc = 1
setwd("P:/MAPPING/CommunityModels773/FINAL/filtered/mix")
all = all_filter
for(i in 1:ncol(all)){
  pid = colnames(all)[i]
  print(pid)
  patient = creIndModMix(specs, all[,i], alldiet, dsel="EU", maxind=300, 
                      replicates=3, econc=conc, n=200, speed=10)
  save(patient,file=paste(pid,"mix","RData",sep="."))
  rm(patient)
  gc()
}


conc = 0.0001
mconc = 1
setwd("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/newest")
all = all_filter
for(i in 1:ncol(all)){
  pid = colnames(all)[i]
  print(pid)
  patient = creIndModMix(specs, all[,i], alldiet, dsel="EU", maxind=500, 
                         replicates=3, econc=conc, n=200, speed=10,tstep=2)
  save(patient,file=paste(pid,"RData",sep="."))
  rm(patient)
  gc()
}

##############FINAL CODE!
conc = 0.0002
setwd("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST")
all = all_filter
for(i in 1:ncol(all)){
  pid = colnames(all)[i]
  print(pid)
  patient = creIndModMixRich(specs, all[,i], maxind=500, 
                         replicates=1, econc=conc, n=100, speed=10,tstep=1)
  save(patient,file=paste(pid,"RData",sep="."))
  rm(patient)
  gc()
}
patient = creIndModMixRich(specs, all[,i], maxind=500, 
                           replicates=1, econc=conc, n=100, speed=10,tstep=1)
res = simEnv(patient,time=24,reduce=T)


conc = 0.0002
tconc = 0.02
setwd("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/FINAL")
load("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/mtreatment.RData")
all = all_filter[,names(tlist)]
for(i in 1:ncol(all)){
  pid = colnames(all)[i]
  print(pid)
  if(length(tlist[[pid]])==0){
    patient = creIndModMixRich(specs, all[,i], maxind=500, 
                               replicates=1, econc=conc, n=100, speed=10,tstep=1)
  }else{
    patient = creIndModMixRichTreat(specs, all[,i], maxind=500, 
                                    replicates=1, econc=conc, n=100, speed=10,tstep=1,treatmet=tlist[[pid]],tconc)
  }
  save(patient,file=paste(pid,"RData",sep="."))
  rm(patient)
  gc()
}

#############################################################################################################
########################################## model tests
#############################################################################################################

load("P:/MAPPING/CommunityModels773/FINAL/filtered/V1.CD.1.RData")
source("P:/GitRep/mapping/analysis/reInitProb.R")
result = simEnv(reInitProb(patient[[1]]), time=12)


#############################################################################################################
########################################## create modified disease models
#############################################################################################################
alldietmod = read.csv("P:/MAPPING/Model_data/diets_eu_mod.csv", header=T, row.names=1)

vol = 2
gconc = 1
load("P:/MAPPING/CommunityModels773/FINAL/filtered/300/Result_V1.CD.1.RData")
larena = getArena(result[[1]],12)
dietsel = alldietmod[larena@mediac,]
larena = addSubs(larena,dietsel[,"EU_mod"]/vol,difspeed=dietsel[,"diffconstant"],unit="mM",mediac=rownames(dietsel),add=F)
larena = addSubs(larena,smax=0.0001,unit="mM",mediac=rownames(dietsel)[which(dietsel$class=="essential")],add=T)
muc = rownames(dietsel)[which(dietsel$type == "Mucin")]
larena = createGradient(larena,smax=gconc,mediac=muc,position='left',steep=0.5,add=F,unit="mM")
if("EX_o2(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_o2(e)",position='left',steep=0.5,add=F,unit="mM")
}
if("EX_no(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_no(e)",position='left',steep=0.5,add=F,unit="mM")
}
larena@orgdat = larena@orgdat[-sample(1:nrow(larena@orgdat),nrow(larena@orgdat)/2),]
save(larena,file="P:/MAPPING/CommunityModels773/FINAL/filtered/300/Mod_V1.CD.1.RData")

load("P:/MAPPING/CommunityModels773/FINAL/filtered/300/Result_V1.CD.12.RData")
larena = getArena(result[[1]],12)
dietsel = alldietmod[larena@mediac,]
larena = addSubs(larena,dietsel[,"EU_mod"]/vol,difspeed=dietsel[,"diffconstant"],unit="mM",mediac=rownames(dietsel),add=F)
larena = addSubs(larena,smax=0.0001,unit="mM",mediac=rownames(dietsel)[which(dietsel$class=="essential")],add=T)
muc = rownames(dietsel)[which(dietsel$type == "Mucin")]
larena = createGradient(larena,smax=gconc,mediac=muc,position='left',steep=0.5,add=F,unit="mM")
if("EX_o2(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_o2(e)",position='left',steep=0.5,add=F,unit="mM")
}
if("EX_no(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_no(e)",position='left',steep=0.5,add=F,unit="mM")
}
larena@orgdat = larena@orgdat[-sample(1:nrow(larena@orgdat),nrow(larena@orgdat)/2),]
save(larena,file="P:/MAPPING/CommunityModels773/FINAL/filtered/300/Mod_V1.CD.12.RData")


load("P:/MAPPING/CommunityModels773/FINAL/filtered/300/Result_V1.CD.15.RData")
larena = getArena(result[[1]],12)
dietsel = alldietmod[larena@mediac,]
larena = addSubs(larena,dietsel[,"EU_mod"]/vol,difspeed=dietsel[,"diffconstant"],unit="mM",mediac=rownames(dietsel),add=F)
larena = addSubs(larena,smax=0.0001,unit="mM",mediac=rownames(dietsel)[which(dietsel$class=="essential")],add=T)
muc = rownames(dietsel)[which(dietsel$type == "Mucin")]
larena = createGradient(larena,smax=gconc,mediac=muc,position='left',steep=0.5,add=F,unit="mM")
if("EX_o2(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_o2(e)",position='left',steep=0.5,add=F,unit="mM")
}
if("EX_no(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_no(e)",position='left',steep=0.5,add=F,unit="mM")
}
larena@orgdat = larena@orgdat[-sample(1:nrow(larena@orgdat),nrow(larena@orgdat)/2),]
save(larena,file="P:/MAPPING/CommunityModels773/FINAL/filtered/300/Mod_V1.CD.15.RData")


load("P:/MAPPING/CommunityModels773/FINAL/filtered/300/Result_V1.CD.6.RData")
larena = getArena(result[[1]],12)
dietsel = alldietmod[larena@mediac,]
larena = addSubs(larena,dietsel[,"EU_mod"]/vol,difspeed=dietsel[,"diffconstant"],unit="mM",mediac=rownames(dietsel),add=F)
larena = addSubs(larena,smax=0.0001,unit="mM",mediac=rownames(dietsel)[which(dietsel$class=="essential")],add=T)
muc = rownames(dietsel)[which(dietsel$type == "Mucin")]
larena = createGradient(larena,smax=gconc,mediac=muc,position='left',steep=0.5,add=F,unit="mM")
if("EX_o2(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_o2(e)",position='left',steep=0.5,add=F,unit="mM")
}
if("EX_no(e)" %in% larena@mediac){
  larena = createGradient(larena,smax=gconc,mediac="EX_no(e)",position='left',steep=0.5,add=F,unit="mM")
}
larena@orgdat = larena@orgdat[-sample(1:nrow(larena@orgdat),nrow(larena@orgdat)/2),]
save(larena,file="P:/MAPPING/CommunityModels773/FINAL/filtered/300/Mod_V1.CD.6.RData")
