##Finding documents [not to execute if inputs already set]
#path="P:/Documenti/HMP_cov/Mapped/Ear"
#out="P:/Documenti/HMP_cov"
btThr=0.80
############# Dependencies

if (!exists("btThr")){
}else{
  warning("bootstrap value found")
}


require(stringr)
#Function to move files 
my.file.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}

############part 1: exctracting abundances
files=list.files(path=path)
matList=list() #list of species abundances 
gnList=list() #list of genus abundances
readsCnt=vector() #vector containing reads number 


files=list.files(path=path)


##Loading all tables with results and create matrix with aboundances
for (f in 1:length(files)){
  it=as.character(f)
  point=paste(it,"/",as.character(length(files)))
  print(point)
  jump=F
  spingo_Tab=read.table(file=paste(path,"/",files[f],sep=""),fill=T) #changed fill 1904 -> spingo bug
  #Detecting column position -> AMBIGUOUS by default 
  index=0
  for (p in 1:length(spingo_Tab[1,])){
    catch0=which(spingo_Tab[,p]=="AMBIGUOUS")
    
    if (length(catch0)!=0){
      print("Species column detected")
      index=p
    }
    
  }
  
  if (index==0){
    print("Sorry, info not located.")
    print("Using previous tab idex")
  }
  
  tot=sum(table(spingo_Tab[, index]))
  ind=sum(table(spingo_Tab[, index])[grep("_",names(table(spingo_Tab[, index])))])+table(spingo_Tab[, index])[which(names(table(spingo_Tab[, index]))=="AMBIGUOUS")]
  container=list()
  
  if ((tot/ind)< 1.7 & index!=0){
    abb=spingo_Tab[,index] #15
    indL=index 
  }
  
  if((tot/ind)> 1.7 & index!=0) {
    warning("It seems that the table is not properly formatted, I will take care jumping one line each iter, this will take long.")
    warning("Wrong parsing might cause unwanted rarefaction. Please manually check tables parsing.")
    k=1
    container=list()
    for (i in 1:length(spingo_Tab[,1])){
      k=k+1 
      catch=as.character(spingo_Tab[k,index]) #6
      container[[i]]=catch
      k=k+1
      jump=T
    }
    
    contvec=unlist(container)
    abb=contvec
    indS=index
  }
  
  
  if ((tot/ind)< 1.7 & index==0){
    abb=spingo_Tab[,indL] #15
  }
  
  
  if((tot/ind)> 1.7 & index==0){
    print("fixing issue")
    k=1
    container=list()
    for (i in 1:length(spingo_Tab[,1])){
      k=k+1
      catch=as.character(spingo_Tab[k,indS]) #6
      container[[i]]=catch
      k=k+1
      jump=T
    }
    
    contvec=unlist(container)
    abb=contvec
  }
  
  if (length(which(is.na(container)))>length(container)/3){
    stop("probable parsing error")
  }
  
  spingMat=as.matrix(table(abb))
  if (length(abb)!=length(spingo_Tab[,1])){
    stop("dim mismatch")
  }
  
  if(sum(as.matrix(table(abb))[,1])!=length(spingo_Tab[,1])){
    stop("dim mismatch")
  }
  
  
  #################
  #ok, we introduce bootstrap score  
  if (exists("btThr")){
  print("taking into account bootstrap score")
  sel=which(spingo_Tab[,index+1] >= btThr)
  abbT=spingo_Tab[sel,]
  sel2=which(abbT[,index-5] >= 0.5)
  abbT=abbT[sel2,index]
  spingMat=as.matrix(table(abbT))
  }
  ################## 
  
  matList[[f]]=spingMat
  readsCnt[f]=sum(as.matrix(table(abb))[,1])
  
  
  #########Now the same but for the genus
  index=0
  for (p in 1:length(spingo_Tab[1,])){
    catch0=which(spingo_Tab[,p]=="Bacteroides")
    
    if (length(catch0)!=0){
      print("Genus column detected")
      index=p
    }
    
  }
  
  if (index==0){
    print("Sorry, info not located.")
    print("Using previous tab idex")
  }
  
  
  if (jump==F & index!=0){
    Gabb=spingo_Tab[,index] #15
    indL=index 
  }
  
  if(jump==T & index!=0){
    print("fixing issue")
    container=list()
    for (i in 1:length(spingo_Tab[,1])){
      k=k+1 
      catch=as.character(spingo_Tab[k,index]) #6
      container[[i]]=catch
      k=k+1
    }
    
    Gcontvec=unlist(container)
    Gabb=Gcontvec
    indS=index
  }
  
  if (jump==F & index==0){
    Gabb=spingo_Tab[,indL] #15
  }
  
  
  if(jump==T & index==0){
    print("fixing issue")
    k=1
    container=list()
    for (i in 1:length(spingo_Tab[,1])){
      k=k+1
      catch=as.character(spingo_Tab[k,indS]) #6
      container[[i]]=catch
      k=k+1
    }
    
    Gcontvec=unlist(container)
    Gabb=Gcontvec
  }
  
  if (length(which(is.na(container)))>length(container)/3){
    stop("probable parsing error")
  }
  
  spingMat=as.matrix(table(Gabb))
  if (length(Gabb)!=length(spingo_Tab[,1])){
    stop("dim mismatch")
  }
  
  if(sum(as.matrix(table(Gabb))[,1])!=length(spingo_Tab[,1])){
    stop("dim mismatch")
  }
  
  ###############
  #ok, what we introduce bootstrap score 
  if (exists("btThr")){
    print("taking into account bootstrap")
    sel=which(spingo_Tab[,index+1] >= btThr)
    abbT=spingo_Tab[sel,]
    sel2=which(abbT[,index-3] >= 0.5)
    abbT=abbT[sel2,index]
    spingMat=as.matrix(table(abbT))
  }
  
  #################
  
  gnList[[f]]=spingMat
}


#Creating and exporting reads count table
rCount=cbind(files,readsCnt)
write.csv(rCount,file=paste(out,"/","readsCount.csv",sep = ""))
#Saving checkpoint
day=as.character(Sys.Date())
save.image(paste(out,paste(day,"checkpoint.RData",sep = ""),sep="/"))

#######Part2 ->creating tables #########################################################
##Finding all the species present the different samples
finSpec=c()
for (i in 1:length(files)){
finSpec=unique(union(finSpec,unique(rownames(matList[[i]]))))
}

##Attributing abundances to all species 
mapTab=matrix(nrow=length(finSpec),ncol = length(files)+1)
mapTab[,1]=finSpec

abcol=as.numeric(c()) #vector with abundances ordered in finSpec order
for (j in 1:length(files)){
 for (i in 1:length(finSpec)){
 index=which(rownames(matList[[j]])==finSpec[i])
  if (length(index)!=0){
  abcol[i]=matList[[j]][index,1]
  }else{
   abcol[i]=0
  }
 }
mapTab[,j+1]=abcol
}

##Attributing names on samples 
tabNam=c()
tabNam[1]="Samples"
for (i in 2:length(files)+1){
sampNam= files[(i-1)]
sampNam=gsub(".out","",sampNam)
tabNam[i]=sampNam
}
mapTab=rbind(tabNam,mapTab)


sampNam= files[1]
sampNam=gsub(".out","",sampNam)
mapTab[1,2]=sampNam

#Now the same for genus level
##Finding all the genus present the different samples
finGn=c()
for (i in 1:length(files)){
  finGn=unique(union(finGn,unique(rownames(gnList[[i]]))))
}

genTab=matrix(nrow=length(finGn),ncol = length(files)+1)
genTab[,1]=finGn

abcol=as.numeric(c()) #vector with abundances ordered in finSpec order
for (j in 1:length(files)){
  for (i in 1:length(finGn)){
    index=which(rownames(gnList[[j]])==finGn[i])
    if (length(index)!=0){
      abcol[i]=gnList[[j]][index,1]
    }else{
      abcol[i]=0
    }
  }
  genTab[,j+1]=abcol
}

##Attributing names on samples 
tabNam=c()
tabNam[1]="Samples"
for (i in 2:length(files)+1){
  sampNam= files[(i-1)]
  sampNam=gsub(".out","",sampNam)
  tabNam[i]=sampNam
}
genTab=rbind(tabNam,genTab)

##Adding missing name
sampNam= files[1]
sampNam=gsub(".out","",sampNam)
genTab[1,2]=sampNam

#############here new part removing ambiguous before normalization of species
remove=which(is.na(str_extract(as.character(mapTab[,1]), "_")))
removeP=remove[-1]

mapTab=mapTab[-removeP,]
mapSt=mapTab
##Normalizing table 
for (i in 2:length(files)+1){
  sum=sum(as.numeric(mapTab[2:length(mapTab[,1]),i]))
  prova=(as.numeric(mapTab[2:length(mapTab[,1]),i])*1)/sum
  mapTab[2:length(mapTab[,1  ]),i]=prova
}

sum=sum(as.numeric(mapTab[2:length(mapTab[,1]),2]))
prova=(as.numeric(mapTab[2:length(mapTab[,1]),2])*1)/sum
mapTab[2:length(mapTab[,1  ]),2]=prova

write.csv(mapTab,file=paste(out,"/","speciesMap.csv",sep = ""),col.names = F)

#############here new part removing ambiguous before normalization of genus
genTab[genTab==""] <- "AMBIGUOUS"
remove=c(which(!is.na(str_extract(as.character(genTab[,1]), "0.00"))),which(!is.na(str_extract(as.character(genTab[,1]), "AMBIGUOUS"))))


genTab=genTab[-remove,]
##Normalizing table 
for (i in 2:length(files)+1){
  sum=sum(as.numeric(genTab[2:length(genTab[,1]),i]))
  prova=(as.numeric(genTab[2:length(genTab[,1]),i])*1)/sum
  genTab[2:length(genTab[,1  ]),i]=prova
}

sum=sum(as.numeric(genTab[2:length(genTab[,1]),2]))
prova=(as.numeric(genTab[2:length(genTab[,1]),2])*1)/sum
genTab[2:length(genTab[,1  ]),2]=prova

write.csv(genTab,file=paste(out,"/","genusMap.csv",sep = ""),col.names = F)


##Uploading list of models names
models=list.files(path="Y:/Microbiome/Eldermet/panModels")

##Formatting models names
species=c()
for ( i in 1:length(models)){
  Nam= models[i]
  Nam=gsub(".mat","",Nam)
  Nam=gsub("pan","",Nam)
  species[i]=Nam
}
length(intersect(species,mapTab[,1]))

agspeciesMat=mapSt
##Creating table shaped on AGORA species
delspecies=setdiff(mapSt[,1],species)

for (i in 1:length(delspecies)){
  index=which(agspeciesMat[,1]==delspecies[i])
  if (length(index)!=0){
    agspeciesMat=agspeciesMat[-index,]
  }
}

agspeciesMat=rbind(tabNam,agspeciesMat)
#Adding name
agspeciesMat[1,2]=sampNam
write.csv(agspeciesMat,file=paste(out,"/","Coverage.csv",sep = ""),col.names = F)

##normalising to 1 after erasing ambiguous empty lines and species not present in agora 
agspeciesMat[1,2]=mapTab[1,2]

panspec=as.character(agspeciesMat[,1])
panspec=paste("pan",panspec,sep="")
agspeciesMat[,1]=panspec

for (i in 2:length(files)+1){
  sum=sum(as.numeric(agspeciesMat[2:length(agspeciesMat[,1]),i]))
  prova=(as.numeric(agspeciesMat[2:length(agspeciesMat[,1]),i])*1)/sum
  agspeciesMat[2:length(agspeciesMat[,1  ]),i]=prova
}


sum=sum(as.numeric(agspeciesMat[2:length(agspeciesMat[,1]),2]))
prova=(as.numeric(agspeciesMat[2:length(agspeciesMat[,1]),2])*1)/sum
agspeciesMat[2:length(agspeciesMat[,1  ]),2]=prova


for (i in 2:length(files)+1){
  print(sum(as.numeric(agspeciesMat[2:length(agspeciesMat[,1]),i])))
}

colnames(agspeciesMat)=as.character(agspeciesMat[1,])
agspeciesMat=agspeciesMat[-1,]
numbers=1:length(agspeciesMat[,1])
agspeciesMat=cbind(numbers,agspeciesMat)
write.csv(agspeciesMat,file=paste(out,"/","normCoverage.csv",sep = ""),col.names = F)
day=as.character(Sys.Date())
save.image(paste(out,paste(day,"workspace.RData",sep = ""),sep="/"))


