#make sure all these packages are installed and loaded
library(NMF)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(Hmisc)
library(spatstat)
library(scales)
library(RColorBrewer)
library(sybil)
library(sybilSBML)

# createSpecies <- function(sbml){
#   model <- readSBMLmod(as.character(sbml), bndCond = FALSE)
#   return(model)
# }
# setwd("Y:/Studies/Microbiome/Stefania/Microbiota_models/SBML_Files_AGORA")
# spec <- list.files(pattern=".xml")
# specs <- list()
# for(i in seq_along(spec)){
#   model <- createSpecies(spec[i])
#   specs[[i]] <- model
#   print(i)
# }
# save(specs,file="P:/AGORA/Almut/specs_agora700.RData")
# load("P:/AGORA/Almut/specs_agora700.RData")
# 
# allrxns = vector()
# for(i in 1:length(specs)){
#   allrxns = union(allrxns,react_id(specs[[i]]))
# }
# rpa = matrix(0,ncol=length(specs),nrow=length(allrxns))
# rownames(rpa) = allrxns
# colnames(rpa) = unlist(lapply(specs,mod_desc))
# for(i in 1:ncol(rpa)){
#   rpa[react_id(specs[[i]]),i] = 1
# }
# rpa = t(rpa)

subsys = read.csv("P:/AGORA/Revision3/TranslatedReactionDatabase_AGORA2.csv",header=T,row.names=2)
subsys$Subsystem = as.factor(ifelse(as.character(subsys$Subsystem)=="",'Unassigned',as.character(subsys$Subsystem)))
#rpass <- rpa[,intersect(colnames(rpa),rownames(subsys))]
#setdiff(colnames(rpa),rownames(subsys))
modstats = read.csv("P:/AGORA/Revision3/ModelInformation.csv",header=T,row.names=1)
rownames(modstats) = modstats$ModelAGORA
# 
# setwd("P:/AGORA/FINAL_DATA") #change this part to reference to your directory

# modstats = read.csv("mod_info.csv",header=T,row.names=1)
# rownames(modstats) = modstats$ModelTrans
# subsys = read.csv("TranslatedReactionFinal.csv",header=T,row.names=1)
# subsys$Subsystem = as.factor(ifelse(as.character(subsys$Subsystem)=="",'Unassigned',as.character(subsys$Subsystem)))
rpa = t(read.csv("P:/AGORA/Revision3/rxn_pa.csv",header=T,row.names=1))
rpass <- rpa[,intersect(colnames(rpa),rownames(subsys))]

############################################################################################################
##################### Heatmaps showing the similarity between the microbes
############################################################################################################

rdist = as.matrix(vegdist(rpass,'jaccard'))
length(which(apply(rpass,2,sum)==nrow(rpass)))
#modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Order, modstats$Family, modstats$Genus),]
modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Genus, modstats$Species),]
rdist = rdist[rownames(modtax),rownames(modtax)]
annotation = data.frame(modtax[,c('Phylum','Class')])
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
ann_colors = list(Class = getPalette(length(levels(modtax$Class))))
aheatmap(rdist,scale="none", Rowv=NA, Colv=NA, labRo=NA, labCol=NA,
         annCol = annotation, annRow=annotation, annColors=ann_colors,
         #color=colorRampPalette(c("navy", "white", "firebrick3"))(100))#9x7
         color=colorRampPalette(c("white", "darkred"))(100))#9x7 #11x9

sspa <- levels(subsys$Subsystem)
ssmat <- matrix(0,nrow=nrow(rpass),ncol=length(sspa))
rownames(ssmat) = rownames(rpass)
colnames(ssmat) = sspa
for(i in 1:nrow(rpass)){
  print(i)
  for(j in 1:ncol(rpass)){
    if(rpass[i,j]==1){
      ssmat[rownames(rpass)[i],as.character(subsys[colnames(rpass)[j],'Subsystem'])] = 
        ssmat[rownames(rpass)[i],as.character(subsys[colnames(rpass)[j],'Subsystem'])] + 1
    }
  }
}
#remove some subsystems based on visual inspections
rmss <- c('Glycerophospholipid metabolism','Exchange/demand reaction','Transport, extracellular', 
          'Fatty acid synthesis','Nucleotide interconversion','Cell wall biosynthesis','Glycolysis/gluconeogenesis', 
          'Peptide metabolism','Tropane, piperidine and pyridine alkaloid biosynthesis', 
          'Polycyclic aromatic hydrocarbon degradation','Hyaluronan metabolism','Chloroalkane and chloroalkene degradation', 
          'Others', 
          'N-glycan synthesis','Miscellaneous','Aminophosphonate Metabolism','Phosphonate and phosphinate metabolism', 
          'Selenoamino acid metabolism','Lipoate metabolism','tRNA Charging','Oxidative phosphorylation', 
          'Biosynthesis of siderophore group nonribosomal peptides', 
          'Aminobenzoate degradation','Aminosugar metabolism','beta-Alanine metabolism','Bile acid synthesis','Cholesterol metabolism', 
          'Chondroitin sulfate degradation','Citric acid cycle','CoA synthesis','Glutathione metabolism','Heparan sulfate degradation', 
          'NAD metabolism','Glycine, serine, alanine and threonine metabolism','Naphthalene degradation','Nucleotide salvage pathway', 
          'Stickland reaction','Tetrahydrobiopterin metabolism','Unassigned','Wood-Ljungdahl Pathway','Urea cycle','ROS detoxification', 
          'Squalene and cholesterol synthesis','Nucleotide sugar metabolism','Purine synthesis','Alanine and aspartate metabolism', 
          'Benzoate degradation',         
          'Ascorbate and aldarate metabolism','Glutamate metabolism','Nitrogen metabolism','Propanoate metabolism',
          'Sulfur metabolism','Thiamine metabolism','Tyrosine metabolism')
modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Order, modstats$Family, modstats$Genus),]
#splot = t(ssmat[rownames(modtax),-which(colnames(ssmat) %in% rmss)])
splot = t(ssmat)
#for(i in 1:ncol(splot)){splot[,i] = splot[,i]/sum(splot[,i])}
annotation = data.frame('Phylum'=modtax[,c('Phylum')],'Class'=modtax[,c('Class')])
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
ann_colors = list(Class = getPalette(length(levels(modtax$Class))))
aheatmap(splot,scale="none", Colv=NA, labCol=NA,
         annCol = annotation, annColors=ann_colors,
         color=colorRampPalette(c("seashell", "darkred"))(100)) #9x7#8x8

############################################################################################################
##################### doing the Pan plot analysis
############################################################################################################

##################### Now doing the same with the HMP data

labs = read.csv("P:/AGORA/MicrobeLabels.csv",header=F)

load("P:/AGORA/Revision3/pan_sample1000.RData")

hmp = read.csv("P:/AGORA/Revision3/HMP_abundance_table_agora_modelid.csv",header=T,row.names=1)
trans = modstats
#rownames(trans) = trans$Organism
#rownames(hmp) = trans[rownames(hmp),'ModelTrans']
#hmpo = hmp
hmp = ifelse(hmp==0,0,1)
#write.csv(hmp,file="P:/AGORA/Supplement/Hmp_microbe_presence_absence.csv")
length(intersect(rownames(trans),rownames(hmp)))
length(which(apply(hmp,1,sum)==ncol(hmp))) #core microbes

pansub = data.frame(X2=1:ncol(pan_sample1000),value=apply(pan_sample1000,2,mean))
maxr = pansub[nrow(pansub),2]
means = pansub[,2]
which(round(means)==maxr*0.95)
#maxr*0.95 #223 which(round(means)== 3032)
#maxr*0.9 #90 which(round(means)== 2873)
#maxr*0.75 #13 which(round(means)== 2397)

mrpass <- rpass

modab_pa <- hmp
sampr <- vector()
sampm <- vector()
for(i in 1:ncol(modab_pa)){
  rnam <- mrpass[rownames(modab_pa)[as.logical(modab_pa[,i])],]
  sampr[i] <- length(which(apply(rnam,2,sum)!=0))
  sampm[i] <- nrow(rnam)
}
length(intersect(rownames(modab_pa),rownames(mrpass)))
sampshmp = data.frame('mic'=sampm, 'rxn'=sampr)

library(RColorBrewer)

ggplot(pansub, aes(X2, value)) + 
  geom_ribbon(aes(x = X2, ymin = value-2*apply(pan_sample1000,2,sd), ymax = value+2*apply(pan_sample1000,2,sd)),fill = "darkgrey", alpha = 0.2) +
  geom_ribbon(aes(x = X2, ymin = value-apply(pan_sample1000,2,sd), ymax = value+apply(pan_sample1000,2,sd)),fill = "darkgrey", alpha = 0.6) +
  geom_line(color='black',size=1.4) +
  ylab("Number of unique reactions") +
  xlab("Microbiota size (number of strains)") +
  scale_y_continuous(limits = c(min(pansub$value), max(pansub$value)),labels = comma) +
  scale_x_continuous(labels = comma) +
  geom_hline(yintercept=maxr,color='red',linetype=2,size=1) +
  #geom_text(aes(730,maxr,label='3,192 reactions',vjust=2),color='red',show_guide = FALSE) +
  geom_segment(aes(x=0,y=maxr*0.75,xend=13,yend=maxr*0.75)) +
  geom_segment(aes(x=13,y=min(pansub$value),xend=13,yend=maxr*0.75)) + #y=753
  #geom_text(aes(1,maxr*0.75,label='75%',vjust=-1)) +
  geom_segment(aes(x=0,y=maxr*0.9,xend=90,yend=maxr*0.9)) +
  geom_segment(aes(x=90,y=min(pansub$value),xend=90,yend=maxr*0.9)) +
  #geom_text(aes(1,maxr*0.9,label='90%',vjust=-1)) +
  geom_segment(aes(x=0,y=maxr*0.95,xend=223,yend=maxr*0.95)) +
  geom_segment(aes(x=223,y=min(pansub$value),xend=223,yend=maxr*0.95)) +
  #geom_text(aes(1,maxr*0.95,label='95%',vjust=-1)) +
  geom_point(data=sampshmp, aes(x=mic, y=rxn), size=3) +
  scale_colour_manual(values=c('red','grey40')) + 
  theme_bw(base_size = 30) +
  theme(legend.position='none',
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=20)) #10x8 #12x8

#######################################################

#function to create one random microbiota pan plot
getRandPan <- function(rdat){
  rdat = rdat[sample(nrow(rdat),nrow(rdat)),]
  rxnall = colnames(rdat)[which(rdat[1,]==1)]
  rxnew = vector()
  rxnew[1] = length(rxnall)
  for(i in 2:nrow(rdat)){
    rxspec = colnames(rdat)[as.logical(rdat[i,])]
    rxnall = union(rxnall,rxspec)
    rxnew[i] = length(rxnall)
  }
  return(rxnew)
}

#doing the 10000 random microbiotas -> takes a long time to compute thats why I saved the computation and commented this part of the code
pan_sample1000 = matrix(0,ncol=nrow(rpass),nrow=1000)
for(i in 1:1000){
  print(i)
  pan_sample1000[i,] = getRandPan(rpass)
}

save(pan_sample1000, file="P:/AGORA/Revision3/pan_sample1000.RData")

#doing the same for the pan models
prpa = t(read.csv("P:/AGORA/Revision3/pan_rxn_pa.csv",header=T,row.names=1))
pan_sample_pan10000 = matrix(0,ncol=nrow(prpa),nrow=1000)
for(i in 1:1000){
  print(i)
  pan_sample_pan10000[i,] = getRandPan(prpa)
}
save(pan_sample_pan10000, file="P:/AGORA/Revision3/pan_sample_pan1000.RData")

##################### plotting now the pan plots with patient data

load("P:/AGORA/Revision3/pan_sample_pan1000.RData")

prpa = t(read.csv("P:/AGORA/Revision3/pan_rxn_pa.csv",header=T,row.names=1))
cleclass = read.csv("P:/AGORA/FINAL_DATA/ClaessonNewData.csv",header=F,row.names=1) #load the claessons patient classes
cle = read.csv("P:/AGORA/Revision3/ClaessonAbundance.csv",header=T,row.names=1) #load the claessons
rownames(cleclass) = gsub('-T[0-9]','',rownames(cleclass))
clepa = ifelse(cle!=0,1,0) #binarize the claessons
#write.csv(clepa,file="P:/AGORA/Supplement/Cleasson_microbe_presence_absence.csv")
length(intersect(rownames(cleclass),colnames(cle)))
#convert the cleasson data into reactions per individual
mrpass <- prpa
modab_pa <- clepa
sampr <- vector()
sampm <- vector()
for(i in 1:ncol(modab_pa)){
  rnam <- mrpass[rownames(modab_pa)[as.logical(modab_pa[,i])],]
  if(is.na(sum(rnam[,1]))){
    rnam <- rnam[-which(is.na(rnam[,1])),]
  }
  sampr[i] <- length(which(apply(rnam,2,sum)!=0))
  sampm[i] <- nrow(rnam)
}

#calculate the mean of the pan microbiotas
pansub = data.frame(X2=1:ncol(pan_sample_pan10000),value=apply(pan_sample_pan10000,2,mean))
maxr = pansub[nrow(pansub),2]
means = pansub[,2]
which(round(means)==round(maxr*0.95))
#maxr*0.95 #123 which(round(means)==2942)
#maxr*0.9 #60 which(round(means)==2787)
#maxr*0.75 #12 which(round(means)==2310)
samps = data.frame(mic=sampm,rxn=sampr,type=cleclass[colnames(modab_pa),6],meth=as.factor(clepa['Methanobrevibacter_smithii',colnames(modab_pa)]))
levels(samps$meth) = c('No Methanogens','Methanogens')

#Making the actual plot
getPalette = colorRampPalette(brewer.pal(8, "Set3"))
cols = getPalette(length(levels(samps$type)))

pan <- ggplot(pansub, aes(X2, value)) + 
  geom_ribbon(aes(x = X2, ymin = value-2*apply(pan_sample_pan10000,2,sd), ymax = value+2*apply(pan_sample_pan10000,2,sd)),fill = "darkgrey", alpha = 0.2) +
  geom_ribbon(aes(x = X2, ymin = value-apply(pan_sample_pan10000,2,sd), ymax = value+apply(pan_sample_pan10000,2,sd)),fill = "darkgrey", alpha = 0.6) +
  geom_line(color='black',size=1.4) +
  ylab("Number of unique reactions") +
  xlab("Microbiota size (number of species)") +
  scale_y_continuous(limits = c(min(pansub$value), max(pansub$value)),labels = comma) +
  scale_x_continuous(limits = c(0, nrow(pansub)),labels = comma) +
  geom_hline(yintercept=maxr,color='red',linetype=2,size=1) +
  #geom_text(aes(300,maxr,label='3,096 reactions',vjust=2),color='red',show_guide = FALSE) +
  geom_segment(aes(x=0,y=maxr*0.75,xend=12,yend=maxr*0.75)) +
  geom_segment(aes(x=12,y=min(pansub$value),xend=12,yend=maxr*0.75)) + #y=753
  #geom_text(aes(1,maxr*0.75,label='75%',vjust=-1)) +
  geom_segment(aes(x=0,y=maxr*0.9,xend=60,yend=maxr*0.9)) +
  geom_segment(aes(x=60,y=min(pansub$value),xend=60,yend=maxr*0.9)) +
  #geom_text(aes(1,maxr*0.9,label='90%',vjust=-1)) +
  geom_segment(aes(x=0,y=maxr*0.95,xend=123,yend=maxr*0.95)) +
  geom_segment(aes(x=123,y=min(pansub$value),xend=123,yend=maxr*0.95)) +
  #geom_text(aes(1,maxr*0.95,label='95%',vjust=-1)) +
  geom_point(data=samps, aes(x=mic, y=rxn, col=type), size=3) +
  #geom_point(data=samps, aes(x=mic, y=rxn), size=3) +
  theme_bw(base_size = 30) +
  theme(legend.position = 'none',
    axis.text.x = element_text(size=20),
    axis.text.y = element_text(size=20),
    axis.title.y = element_text(size=30,vjust=0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour='black',size=2),
    axis.ticks = element_line(size=1,color='black'),
    plot.title = element_text(size=20)) #12x8

subpan <- ggplot(pansub, aes(X2, value)) + 
  scale_x_continuous(limits = c(min(samps$mic)-2,max(samps$mic)+2)) +
  scale_y_continuous(limits = c(min(samps$rxn)-2,max(samps$rxn)+2)) +
  geom_point(data=samps, aes(x=mic, y=rxn, col=type), size=3) +
  theme_bw(base_size = 30) +
  theme(legend.position=c(0.5, 1.3),
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=1.5),
        axis.ticks = element_blank(),
        plot.title = element_blank()) #12x8

pan + geom_rect(xmin=min(samps$mic)-2, ymin=min(samps$rxn)-2, xmax=max(samps$mic)+2, ymax=max(samps$rxn)+2, 
                fill=NA, colour='black',size=0.5) +
  annotation_custom(ggplotGrob(subpan), xmin=50+100, xmax=200+100, ymin=1200, ymax=2200)

############################################################################################################
##################### Calculating some statistics on the microbe selection
############################################################################################################

########################################### now doing the PCoA with the pan models

mrpass <- prpa
modab_pa <- clepa
patrxn = matrix(0,nrow=ncol(modab_pa),ncol=ncol(mrpass))
colnames(patrxn) = colnames(mrpass)
rownames(patrxn) = colnames(modab_pa)
for(i in 1:ncol(modab_pa)){
  rnam <- mrpass[rownames(modab_pa)[as.logical(modab_pa[,i])],]
  if(is.na(sum(rnam[,1]))){
    rnam <- rnam[-which(is.na(rnam[,1])),]
  }
  patrxn[i,names(which(apply(rnam,2,sum)!=0))] = 1
}
length(which(apply(patrxn,2,sum)==nrow(patrxn))) #number of reactions that occur in all individuals
hpat = patrxn[,-which(apply(patrxn,2,sum)==nrow(patrxn))] #remove those reactions that occur in all individuals
annotation = data.frame("Groups"=cleclass[rownames(hpat),6],'Individuals'=cleclass[rownames(hpat),2]) #cle$V7,gam=cle$V
aheatmap(hpat, scale="none", labRo=NA, labCol=NA,
         annRow = annotation, #annColors=ann_colors,
         color=colorRampPalette(c("seashell", "brown"))(100))#9x7.5

#reactions within patient groups
old = patrxn[rownames(cleclass[which(cleclass$V6=='Old'),]),]
young = patrxn[rownames(cleclass[which(cleclass$V3=='Young control'),]),]
rehab = patrxn[rownames(cleclass[which(cleclass$V3=='Rehabilitation'),]),]
mean(apply(old,1,sum))
sd(apply(old,1,sum))
mean(apply(young,1,sum))
sd(apply(young,1,sum))
mean(apply(rehab,1,sum))
sd(apply(rehab,1,sum))

patreac = as.matrix(prpa[,names(which(apply(patrxn,2,sum)==nrow(patrxn)))])
ncol(patreac) #number of reactions in pan reactom
length(which(apply(patreac,2,sum)==nrow(patreac))) #number of reactions present in all individuals
length(which(apply(patreac,2,sum)==1)) #number of reactions present in all individuals
plot(sort(apply(patreac,2,sum)),xlab='Number of reactions',ylab='Number of pan microbes',type='l')

dires <- capscale(patrxn~1, distance="jaccard")
dires_sum <- summary(dires)
imp = 150
abseigen <- abs(dires_sum$species)
absort <- sort(apply(abseigen, 1, sum), decreasing=T)
sel <- names(head(absort, imp))

selrxn <- dires_sum$species[,1:2]
selrxn <- selrxn[names(which(absort!=0)),]
write.csv(data.frame("PCmode1"=selrxn[,1],"PCmode2"=selrxn[,2],'Subsystem'=subsys[rownames(selrxn),'Subsystem'],'Summed_eigenvalues'=absort[rownames(selrxn)]),
          file='P:/AGORA/Revision3/pcoa_eigenvalues.csv')

gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], "Groups"=cleclass[rownames(dires_sum$sites),6], 
                   'meth'=as.factor(clepa['Methanobrevibacter_smithii',rownames(hpat)]),
                   'R1'=c(dires_sum$species[sel,1],rep(NA,(nrow(dires_sum$sites)-imp))),
                   'R2'=c(dires_sum$species[sel,2],rep(NA,(nrow(dires_sum$sites)-imp))),
                   'labs'=c(sel,rep(NA,(nrow(dires_sum$sites)-imp))))
levels(gdat$meth) = c('No Methanogens','Methanogens')
levels(gdat$Groups) = c('G-proteobacteria','No G-proteobacteria')
pcoap <- ggplot(data=gdat, aes(x=PC1, y=PC2)) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=0.5) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=0.5) +   
  #geom_segment(data=gdat, aes(x=0, y=0, xend=R1*5, yend=R2*5), arrow=arrow(length=unit(0.2,"cm")), alpha=0.5, color="darkgrey", size=0.4)+
  #geom_text(aes(x=R1*5, y=R2*5, label=labs), color="black", size=2)+
  geom_point(aes(x=PC1, y=PC2, colour=Groups, shape=meth), size=3)+#, colour="black") +
  xlab(paste("PC1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
  ylab(paste("PC2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
  theme_bw(base_size = 15) +
  theme(legend.position='top',
        legend.text=element_text(size=14),
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=20))
pcoap

pan + geom_rect(xmin=min(samps$mic)-2, ymin=min(samps$rxn)-2, xmax=max(samps$mic)+2, ymax=max(samps$rxn)+2, 
                fill=NA, colour='black',size=0.5) +
  annotation_custom(ggplotGrob(pcoap), xmin=50+100, xmax=200+100, ymin=1200, ymax=2550)
