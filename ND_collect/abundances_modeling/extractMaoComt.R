medfile=read.csv(file="C:/Users/federico.baldini/Downloads/lux_park.clinical (49).csv",header=T,fill = T)

idsel=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/finFiles/id_list_final.csv",header=T)
idsel=as.character(idsel[,1])

ind=vector()
for (i in 1:length(idsel)){
  ind[i]=which(as.character(medfile[,1])==idsel[i])
}

medsel=medfile[ind,]
write.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/resubmission/medications.csv",medsel)
