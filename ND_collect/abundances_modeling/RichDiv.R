#Computing species richness and diversity. 

# Definitions from http://www.flutterbys.com.au/stats/tut/tut13.2.html and http://www.metagenomics.wiki/pdf/definition/alpha-beta-diversity
#and http://www2.uaem.mx/r-mirror/web/packages/vegan/vignettes/diversity-vegan.pdf and https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html and
# http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html
 
                          
                                     ###Species level
#Importing species table

speciesTab=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/finFiles/speciesMap.csv")
rownames(speciesTab)=as.character(speciesTab[,1])
speciesTab=speciesTab[,-1]

#Importing the final list of observations to select data of interest
finObs=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/finFiles/id_list_final.csv",header=T)
finObs=as.character(finObs[,1])

#Filtering species Tab

speciesTab=speciesTab[finObs,]

rm=vector()

f=1
for (i in 1:length(speciesTab[1,])){
  if (sum(as.numeric(speciesTab[,i]))==0){
    rm[f]=i
    f=f+1
  }
}

speciesTab=speciesTab[,-rm]

#Species richness

rich=vector()
for (i in 1:length(speciesTab[,1])){
  rich[i]=length(which(as.numeric(speciesTab[i,])!=0))
}

speciesInfo=cbind(rownames(speciesTab),rich) #table with all information 


#Alpha diversity (shannon)
library(vegan)
adiv=diversity(speciesTab, index="shannon")
speciesInfo=cbind(speciesInfo,adiv)

#Pilou evenness

S <- apply(speciesTab>0,1,sum)
PilouEven=diversity(speciesTab, index="simpson")/log(S)

speciesInfo=cbind(speciesInfo,PilouEven)

#Beta diversity (Bray-Curtis dissimilarity)

beta_dist <- vegdist(t(speciesTab),
                     index = "bray")

beta_dist <- vegdist((speciesTab),
                     index = "bray")

mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)

#Now we load samples info

met=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/finFiles/finalMetadata.csv",header=T)
rownames(met)=as.character(met[,1])
met=met[rownames(speciesTab),]

mds_data$SampleID <- rownames(mds_data)

mds_data2=cbind(mds_data,met$patstat)

library(ggplot2)
ggplot(mds_data, aes(x = MDS1, y = MDS2, color = met$patstat)) +
  geom_point()

speciesInfo=cbind(speciesInfo,mds_data2$MDS1,mds_data2$MDS2)
colnames(speciesInfo)=c("ID","Richness (species)","Alpha diversity (Shannon)","Eveness (Pilou)", "Beta diversity (Bray MDS1)", "Beta diversity (Bray MDS2)" )
write.csv(speciesInfo, file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/resubmission/spDiversity.csv")

#Detecting if any feature is different for the main predictor (PD Vs. Control)
tags=c("PD","Control")
speciesInfo=cbind(speciesInfo,as.character(met$patstat))

fmat=as.numeric(c())

for (i in 2:(length(speciesInfo[1,])-1)){
f=wilcox.test(as.numeric(as.character(speciesInfo[which(speciesInfo[,7]==tags[1]),i])),as.numeric(as.character(speciesInfo[which(speciesInfo[,7]==tags[2]),i])),exact=F)
fmat[i]=f$p.value
}

names(fmat)=colnames(speciesInfo)[-7]

p.adjust(fmat, method = "bonferroni", n = 5)
p.adjust(fmat, method = "fdr", n = 5)

#Plotting

ggplot(mds_data2, aes(x = met$patstat, y = MDS1)) +
  geom_text(data = data.frame(),
             aes(x = (met$patstat), y = max(mds_data2$MDS1) , label = met$patstat),
             col = 'black',
             size = 10) +
  geom_boxplot() +
  ggtitle("Beta diversity MDS1") +
  xlab("Site") +
  ylab("Beta diversity index")

#Adonis test

spBdist=as.matrix(beta_dist)
write.csv(spBdist,file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/resubmission/spBdist.csv")
resSp=adonis(formula = beta_dist ~ patstat, data = met) 


#Anosim


attach(met)
sp.ano <- anosim(beta_dist, patstat)
summary(sp.ano)
plot(sp.ano)



########Now the same but for the genus level#########################


#Importing genus table

genusTab=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/finFiles/genusMap.csv")
rownames(genusTab)=as.character(genusTab[,1])
genusTab=genusTab[,-1]




#Filtering species Tab

genusTab=genusTab[finObs,]

rm=vector()

f=1
for (i in 1:length(genusTab[1,])){
  if (sum(as.numeric(genusTab[,i]))==0){
    rm[f]=i
    f=f+1
  }
}

genusTab=genusTab[,-rm]

#Genus richness

Grich=vector()
for (i in 1:length(genusTab[,1])){
  Grich[i]=length(which(as.numeric(genusTab[i,])!=0))
}

genusInfo=cbind(rownames(genusTab),Grich) #table with all information 


#Alpha diversity (shannon)

Gadiv=diversity(genusTab, index="shannon")
genusInfo=cbind(genusInfo,Gadiv)

#Pilou evenness

Sg <- apply(genusTab>0,1,sum)
GPilouEven=diversity(genusTab, index="simpson")/log(Sg)

genusInfo=cbind(genusInfo,GPilouEven)

#Beta diversity (Bray-Curtis dissimilarity)

# beta_dist <- vegdist(t(speciesTab),
#                      index = "bray")

Gbeta_dist <- vegdist((genusTab),
                     index = "bray")

Gmds <- metaMDS(Gbeta_dist)
Gmds_data <- as.data.frame(Gmds$points)

#Now we load samples info


Gmds_data$SampleID <- rownames(Gmds_data)

Gmds_data2=cbind(Gmds_data,met$patstat)

ggplot(Gmds_data, aes(x = MDS1, y = MDS2, color = met$patstat)) +
  geom_point()

genusInfo=cbind(genusInfo,Gmds_data2$MDS1,Gmds_data2$MDS2)
colnames(genusInfo)=c("ID","Richness (genus)","Alpha diversity (Shannon)","Eveness (Pilou)", "Beta diversity (Bray MDS1)", "Beta diversity (Bray MDS2)" )
write.csv(genusInfo, file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/resubmission/gnDiversity.csv")

#Detecting if any feature is different for the main predictor (PD Vs. Control)
tags=c("PD","Control")
genusInfo=cbind(genusInfo,as.character(met$patstat))

Gfmat=as.numeric(c())

for (i in 2:(length(genusInfo[1,])-1)){
  f=wilcox.test(as.numeric(as.character(genusInfo[which(genusInfo[,7]==tags[1]),i])),as.numeric(as.character(genusInfo[which(genusInfo[,7]==tags[2]),i])),exact=F)
  Gfmat[i]=f$p.value
}

names(Gfmat)=colnames(genusInfo)[-7]

p.adjust(Gfmat, method = "bonferroni", n = 5) #!errror
p.adjust(Gfmat, method = "fdr", n = 5) #!errror

#Plotting

ggplot(Gmds_data2, aes(x = met$patstat, y = MDS2)) +
  geom_text(data = data.frame(),
            aes(x = (met$patstat), y = max(Gmds_data2$MDS2) , label = met$patstat),
            col = 'black',
            size = 10) +
  geom_boxplot() +
  ggtitle("Beta diversity MDS2") +
  xlab("Site") +
  ylab("Beta diversity index")

#Adonis test

gnBdist=as.matrix(Gbeta_dist)
write.csv(gnBdist,file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/resubmission/gnBdist.csv")
resGn=adonis(formula = Gbeta_dist ~ patstat, data = met) 


#Anosim

gn.ano <- anosim(Gbeta_dist, patstat)
summary(gn.ano)
plot(gn.ano)




################ Analysis on number of sequences per sample ##############

seqNumb=read.csv(file="Y:/Microbiome/NDcollect/HQ/tablesBS/treshold19/resubmission/seqNumb.csv",header = T)
rownames(seqNumb)=as.character(seqNumb[,1])
seqNumb=seqNumb[,-1]


seqNumb=cbind(seqNumb,as.character(met$patstat))

seqMat=as.numeric(c())

for (i in 2:(length(seqNumb[1,])-1)){
  f=wilcox.test(as.numeric(as.character(seqNumb[which(seqNumb[,7]==tags[1]),i])),as.numeric(as.character(seqNumb[which(seqNumb[,7]==tags[2]),i])),exact=F)
  seqMat[i]=f$p.value
}

names(seqMat)=colnames(seqNumb)[-7]

seqMat



  

