#This script applies a threshold to normCoverage and splits it in several sub files to execute with k fold batch strategy

#coverage=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/normCoverage.csv")
load("E:/BackupInesLeave/ND/HQ/tablesBS/2018-12-07workspace.RData")
coverage=agspeciesMat
load("E:/BackupInesLeave/ND/HQ/otherScripts/11Mars8pm.RData")


numcoverage=coverage[,3:length(coverage[1,])]

#Here we will replace values < 0.0001 with 0
for (i in 1:length(numcoverage[1,])){
  wildPesce=as.numeric(as.character(numcoverage[,i]))
  wildPesce[(which(wildPesce<0.0001))]=0
  numcoverage[,i]=wildPesce
}

#As last we need a normalization to 1 again
for (i in 1:length(numcoverage[1,])){
  numcoverage[,i]=as.numeric(numcoverage[,i])/sum(as.numeric(numcoverage[,i]))
}

bigCoverage=cbind(coverage[,1:2],numcoverage)

#Exporting a new normCoverage which has no empty lines
pesce=vector()
for (i in 1:length(numcoverage[,1])){
  if (sum(numcoverage[i,])==0){ 
    pesce[i]=i
  }
}
rmC=pesce[-which(is.na(pesce)==1)]
dim(numcoverage[-rmC,]);dim(numcoverage) #Dimension comparison

bigCoverage=bigCoverage[-rmC,]

#here we select only patients we are interested in
rep=colnames(bigCoverage)[3:length(bigCoverage[1,])]
rep=gsub(".STLOMN.*","",rep)
rep=gsub(".N","-N",rep)
rep=gsub("D.","D-",rep)
rep=gsub("X","",rep)
bigCoverage2=bigCoverage
colnames(bigCoverage2)[3:length(bigCoverage[1,])]=rep
bigCoverage2=bigCoverage2[,rownames(patstat)]
bigCoverage2=cbind(bigCoverage[,1:2],bigCoverage2)

write.csv(file = "E:/BackupInesLeave/ND/HQ/tablesBS/treshold19/normCoverage.csv",bigCoverage2)

#We decided to go for a k fold batch strategy...let's do th efolders and the splits
storeCoverage=bigCoverage
bigCoverage=bigCoverage2
dir="E:/BackupInesLeave/ND/HQ/tablesBS/treshold19"
setwd(dir)

kfolds=6 #number of equal or almost equal parts of whoch deivide the dataset
totObs=length(bigCoverage[1,])-2

nfolds=round(totObs/kfolds,digits = 0)

end=0
init=3
#split1
for (i in 1:kfolds){
if (i<kfolds){
  split=bigCoverage[,init:(nfolds+end+2)]
  end=end+nfolds
  init=init+nfolds
  subDir=paste("batch",as.character(i),sep="")
  dir.create(file.path(dir, subDir))
  split=cbind(bigCoverage[,1:2],split)
  write.csv(file=paste(file.path(dir, subDir),"normCoverage.csv",sep="/"), split)
  dir.create(file.path(dir,subDir, "results"))
}else{
  split=bigCoverage[,init:length(bigCoverage[1,])]
  init=init+nfolds
  subDir=paste("batch",as.character(i),sep="")
  dir.create(file.path(dir, subDir))
  split=cbind(bigCoverage[,1:2],split)
  write.csv(file=paste(file.path(dir, subDir),"normCoverage.csv",sep="/"), split)
}
}





