
coverage = do.call(cbind,lapply(reps,function(x){return(x$coverage)}))
colnames(coverage)=gsub(".filter.sort.bam","",unlist(lapply(reps,function(x){return(x$name)})))

# Code to replace NA with 0 in coverage file: NOTE THAT INDEX STARTS FROM 1 AS HEADER IS T
for (i in 1:length(coverage[,1])){
  for (j in 1:length(coverage[1,])){   #for (j in 2:length(coverage[1,])){
    if (is.na(coverage[i,j])){
      coverage[i,j]=0 #if the coverage is NA convert it to 0 for the sum function
    }
  }
}

tcoverage=coverage

for (i in 1:length(tcoverage[,1])){
  for (j in 1:length(tcoverage[1,])){
    if (tcoverage[i,j]<0.1){
      tcoverage[i,j]=0 #if the coverage is NA convert it to 0 for the sum function
    }
  }
}

for (i in 1:length(tcoverage[1,])){
  tcoverage[,i]=tcoverage[,i]/sum(tcoverage[,i]) 
}

rm=vector()
for (i in 1:length(tcoverage[,1])){
  if(sum(tcoverage[i,])==0){
    rm[i]=i
  }
}
rm=rm[-which(is.na(rm)==1)]

tcoverage=tcoverage[-rm,]
write.csv(file="Y:/Microbiome/PD/DeNoPa_scripts/normCoverage.csv",tcoverage)