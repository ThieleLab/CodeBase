#Script used to retreive number of unclassified/low quality classified and total reads

load("Y:/Microbiome/NDcollect/HQ/tablesBS/2018-12-07workspace.RData")


#Number of sequences for each sample classified and retained for relative abundances computation
abbAmpSp=vector()
abbAmpGn=vector()

for (i in 1:length(matList)){
  abbAmpSp[i]=sum(as.numeric(matList[[i]][,1])) #genus
  abbAmpGn[i]=sum(as.numeric(gnList[[i]][,1])) #species
}


#Now we remember that the total number of reads in each sample are called readsCnt and rCount

NoabbAmpSp=as.numeric(rCount[,2])-abbAmpSp
NoabbAmpGn=as.numeric(rCount[,2])-abbAmpGn

# On average the percentuage of excluded reads 
summary(NoabbAmpSp/as.numeric(rCount[,2]))
summary(NoabbAmpGn/as.numeric(rCount[,2]))

#Compiling a table with classification information on the reads 

sumTab=cbind(rCount,abbAmpSp,abbAmpGn,NoabbAmpSp,NoabbAmpGn)

#Formattting names
sumTab[,1]=gsub(".out","",sumTab[,1])

rep=sumTab[,1]
rep=gsub(".STLOMN.*","",rep)
rep=gsub(".N","-N",rep)
rep=gsub("D.","D-",rep)

sumTab[,1]=rep
rownames(sumTab)=rep

#Importing the final list of observations to select data of interest

finObs=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/finFiles/id_list_final.csv",header=T)
finObs=as.character(finObs[,1])

sumTab=sumTab[finObs,]
write.csv(sumTab, file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/resubmission/seqNumb.csv")



