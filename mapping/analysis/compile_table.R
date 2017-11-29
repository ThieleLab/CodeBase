library(vegan)
library(ggplot2)
library(Rsamtools)
library(tsne)
library(devtools)
library(digest)
library(gridExtra)
#source_url("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R")    

setwd('P:/MAPPING/Pediatric_crohns')
#setwd('P:/MAPPING/Results773/bwa_human_filtered/filter1')
coverage = as.matrix(read.csv("coverage773.csv",row.names=1,header=T))
coverage = as.matrix(read.csv("coverage773_pediatric_dysbiotic.csv",row.names=1,header=T))

modstats = read.csv("P:/AGORA/Revision3/ModelInformation.csv",header=T,row.names=3)
cats = read.csv("P:/MAPPING/MetaHit/spain_metahit.csv")
pclass = read.csv("P:/MAPPING/Model_data/class773.csv", header=F, row.names=1)
cats = read.csv("P:/MAPPING/Model_data/dataset2_class.csv")

cats = read.csv("P:/MAPPING/Pediatric_crohns/dysbiotic_cluster.csv",row.names=2)
catsel = data.frame(Reads=gsub(" S","S",rownames(cats)),Sample.Name=gsub(" S","S",rownames(cats)),category=cats$Disease)
catsel = cbind(catsel,cats)
catsel$FCP = ifelse(is.na(cats$FCP),0,cats$FCP)/max(ifelse(is.na(cats$FCP),0,cats$FCP))
catsel$Time = ifelse(is.na(catsel$Time),1,catsel$Time)
rownames(catsel) = as.character(catsel$Reads)

test0001 = cov2relab(coverage, modstats, catsel, mincoverage=0.01)
if(length(which(is.na(apply(test0001,2,sum))))!=0){test0001 = test0001[,-which(is.na(apply(test0001,2,sum)))]}
if(length(which(apply(test0001,1,sum)==0))!=0){test0001 = test0001[-which(apply(test0001,1,sum)==0),]}
catsel$category = as.character(catsel$category)
catsel$category[which(catsel$category == "Control")] = "Healthy"
catsel$category[which(catsel$category == "Crohn")] = "Crohn's disease"
catsel$category = as.factor(catsel$category)
pclass = data.frame(V2=catsel$category,V3=catsel$Reads)
rownames(pclass) = as.character(catsel$Reads)

apply(ifelse(test0001>0,1,0),2,sum)
min(apply(ifelse(test0001>0,1,0),2,sum))
hist(apply(ifelse(test0001>0,1,0),2,sum))

#write.csv(test0001,file="alldat773_pediatric_dysbiotic_filtered.csv")
#write.csv(test0001,file="alldat773_allibd.csv")

#catsel$FCP=ifelse(is.na(cats$FCP),0,cats$FCP)/max(ifelse(is.na(cats$FCP),0,cats$FCP))*10
fcpcoa <- plotPCoA(test0001, imp=15, catsel, gtit="FCP", secvar="FCP")
fcpcoa #9x9
catsel$PUCAI=(ifelse(is.na(cats$PUCAI),0,cats$PUCAI)/max(ifelse(is.na(cats$PUCAI),0,cats$PUCAI)))*10
pucaipcoa <- plotPCoA(test0001, imp=15, catsel, gtit="PUCAI", secvar="PUCAI")
catsel$PCDAI=(ifelse(is.na(cats$PCDAI),0,cats$PCDAI)/max(ifelse(is.na(cats$PCDAI),0,cats$PCDAI)))*10
pcdaipcoa <- plotPCoA(test0001, imp=15, catsel, gtit="PCDAI", secvar="PCDAI")
catsel$BristolScore=ifelse(is.na(cats$BristolScore),0,cats$BristolScore)
brispcoa <- plotPCoA(test0001, imp=15, catsel, gtit="Bristol Score", secvar="BristolScore")
catsel$Cluster=ifelse(is.na(cats$Cluster),0,cats$Cluster)
cluspcoa <- plotPCoA(test0001[,which(catsel[colnames(test0001),]$Time==1)], imp=15, catsel, gtit="Bristol Cluster", secvar="Cluster")

grid.arrange(fcpcoa,pucaipcoa,pcdaipcoa,brispcoa,ncol=2)#9x17
#dotsne(test0001, catsel)

timesel = test0001[,which(catsel[colnames(test0001),]$Time==1)]
plotPCoA(timesel, imp=15, catsel, gtit="Bristol Cluster", secvar="Cluster", dismeth="bray")#manhattan, jaccard, cranberra
#dotsne(timesel, catsel)
barphyla(timesel,modstats,pclass)
catsel$Treatment = ifelse(is.na(catsel$Treatment),'healthy',as.character(catsel$Treatment))
plotPCoA(timesel[,which(catsel[colnames(timesel),]$Treatment%in%c("EEN","healthy"))], 
         imp=15, catsel, gtit="Bristol Cluster", secvar="Cluster", dismeth="manhattan")

#write.csv(timesel,file="alldat773_pediatric_filtered.csv")

cov2relab <- function(coverage, modstats, cats, mincoverage, rm=NULL){
  rownames(modstats) = modstats$ModelAGORA
  rownames(modstats)[which(rownames(modstats)=="Ruminococcus_albus_7_DSM_20455")] = "Ruminococcus_albus_7"
  rownames(modstats)[which(rownames(modstats)=="Bifidobacterium_thermacidophilum_subsp_thermacidophilum_DSM_158")] = "Bifidobacterium_thermacidophilum_subsp_thermacidophilum_DSM_15837"
  dis = as.character(cats$category)
  names(dis) = gsub("-",".",as.character(cats$Sample.Name))
  coverage = ifelse(is.na(coverage),0,coverage)
  cover = ifelse(coverage<mincoverage,0,coverage)
  cover = apply(cover,2,function(x){x/sum(x)})
  colnames(cover) = gsub("-",".",colnames(cover))
  if(!is.null(rm)){cover = cover[,-which(dis[colnames(cover)] == rm)]}
  return(cover)
}
barphyla <- function(cover, modstats, pclass){
  tcover = t(cover)
  pdat = data.frame()
  for(i in c("Actinobacteria","Bacteroidetes","Firmicutes","Proteobacteria","Verrucomicrobia")){
    ab = tcover[,intersect(colnames(tcover),rownames(modstats)[which(modstats$Phylum==i)])]
    if(sum(ab)!=0){
      if(class(ab)!="numeric"){
        ab=apply(ab,1,sum)
      }
      pdat = rbind(pdat,data.frame(abundance=ab,patient=rownames(tcover),
                                   type=as.character(pclass[rownames(tcover),1]),phylum=rep(i,nrow(tcover))))
    }
  }
  pdat$patient = factor(as.character(pdat$patient),levels=rev(names(sort(apply(tcover[,rownames(modstats)[which(modstats$Phylum=="Bacteroidetes")]],1,sum)))))
  ggplot(pdat, aes(x=patient, y=abundance, fill=phylum)) +
    geom_bar(stat="identity")
  # ggplot(pdat[which(pdat$phylum=="Firmicutes"),], aes(x=type, y=abundance)) +
  #   geom_boxplot(aes(fill = factor(type)))
  # 
  # pdat = data.frame()
  # for(i in c("Faecalibacterium","Bacteroides","Enterococcus","Eubacterium","Escherichia")){
  #   ab = tcover[,intersect(colnames(tcover),rownames(modstats)[which(modstats$Genus==i)])]
  #   if(sum(ab)!=0){
  #     if(class(ab)!="numeric"){
  #       ab=apply(ab,1,sum)
  #     }
  #     pdat = rbind(pdat,data.frame(abundance=ab,patient=rownames(tcover),
  #                                  type=as.character(pclass[rownames(tcover),1]),phylum=rep(i,nrow(tcover))))
  #   }
  # }
  # ggplot(pdat[which(pdat$phylum=="Escherichia"),], aes(x=type, y=abundance)) +
  #   geom_boxplot(aes(fill = factor(type)))
}
plotPCoA <- function(cover, imp=15, cats, gtit="", secvar=NA, dismeth="bray"){
  rownames(cats) = as.character(cats$Reads)
  dis = as.character(cats$category)
  names(dis) = gsub("-",".",as.character(cats$Sample.Name))
  dires <- capscale(t(cover)~1, distance=dismeth)
  #dires <- capscale(t(ifelse(cover==0,0,1))~1, distance="jaccard")
  dires_sum <- summary(dires)
  abseigen <- abs(dires_sum$species)
  absort <- sort(abs(abseigen[,2]), decreasing=T)
  sel <- names(head(absort, imp))
  gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                     "Phenotype"=as.factor(dis[rownames(dires_sum$sites)]))
  gdat2 <- data.frame('R1'=dires_sum$species[sel,1],
                      'R2'=dires_sum$species[sel,2],
                      'labs'=sel)
  levels(gdat$Phenotype) = c("Crohn's disease", "Healthy")
  gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC1"])
  gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
  gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC2"])
  gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
  gdat$secvar = cats[rownames(dires_sum$sites),secvar]
  fac = 2
  #ggplot(data=gdat, aes(x=PC1, y=PC2, color=Phenotype, size=secvar)) +
  ggplot(data=gdat, aes(x=PC1, y=PC2, color=Phenotype)) +
    geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=1) +
    geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=1) +
    stat_ellipse(level = 0.95, size=1) +
    #geom_segment(data=gdat2, aes(x=0, y=0, xend=R1*fac, yend=R2*fac), arrow=arrow(length=unit(0.2,"cm")), alpha=0.5, color="darkgrey", size=0.4)+
    #geom_text(data=gdat2, aes(x=R1*fac, y=R2*fac, label=labs), color="black", size=3) +
    #scale_colour_manual(values=taxcols) +
    #scale_shape_manual(values=c(1,16,18,17,2)) +#c(1,13,4,3,2)
    ggtitle(gtit) +
    scale_colour_manual(values = c("firebrick2","royalblue2")) +
    #stat_ellipse(level = 0.95,type="t", size=1) +
    geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=0.8) +
    geom_point(aes(x=PC1, y=PC2),size=3) +
    #geom_point(aes(x=mean.x, y=mean.y), shape=4, size=5, stroke=2) +
    xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
    ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
    theme_bw(base_size = 30) +
    theme(legend.text=element_text(size=20),
          legend.position="top",
          legend.key=element_blank(),
          legend.title=element_blank(),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.border = element_rect(colour='black',size=2),
          axis.ticks = element_line(size=1,color='black'),
          axis.ticks.length=unit(0.3,'cm'),
          plot.title = element_text(size=25,hjust = 0.5)) +
    guides(colour = guide_legend(override.aes = list(shape = 15))) #width=12, height=8)
}
dotsne = function(cover, cats, dismeth="bray"){
  dis = as.character(cats$category)
  names(dis) = gsub("-",".",as.character(cats$Sample.Name))
  bdist = vegdist(t(cover), method=dismeth)
  ecb = function(x,y){plot(x,t='n', axes=FALSE, frame.plot=F, xlab = NA, ylab = NA);
    text(x,labels=colnames(cover),
         col=as.numeric(as.factor(dis[colnames(cover)])),cex=0.6)}
  rtsne = tsne(bdist, epoch_callback = ecb, perplexity=10, max_iter = 3000) #8x8
  gdat <- data.frame("axis1"=rtsne[,1], "axis2"=rtsne[,2],
                     "Phenotype"=as.factor(dis[colnames(cover)]))
  levels(gdat$Phenotype) = c("Crohn's disease", "Healthy")
  ggplot(data=gdat, aes(x=axis1, y=axis2)) +
    geom_point(aes(x=axis1, y=axis2, colour=Phenotype, shape=Phenotype), size=4)+#, colour="black") +
    xlab("axis1") +
    ylab("axis2") +
    theme_bw(base_size = 30) +
    theme(legend.text=element_text(size=20),
          legend.key=element_blank(),
          legend.title=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour='black',size=2),
          axis.ticks = element_blank(),
          plot.title = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank())#15x7
}

#checking how different filtered and unfiltered data is
setwd('P:/MAPPING/Results773/bwa_human_filtered/filter1')
coverage = as.matrix(read.csv("coverage773.csv",row.names=1,header=T))
coverage = ifelse(is.na(coverage),0,coverage)
length(which(apply(coverage,1,sum)>0))
setwd('P:/MAPPING/Results773/bwa_human_filtered/unfiltered/')
coverage = as.matrix(read.csv("coverage773.csv",row.names=1,header=T))
coverage = ifelse(is.na(coverage),0,coverage)
length(which(apply(coverage,1,sum)>0))


setwd('P:/MAPPING/Results773/bwa_human_filtered/filter1')
coverage = as.matrix(read.csv("coverage773.csv",row.names=1,header=T))
coverage = ifelse(is.na(coverage),0,coverage)
length(which(apply(coverage,1,sum)>0))


setwd('P:/MAPPING/Results773/bwa_human_filtered/unfiltered')
coverage = as.matrix(read.csv("coverage773.csv",row.names=1,header=T))
test0001 = cov2relab(coverage, modstats, cats, pclass, mincoverage=0.01,rm="UC")
plotPCoA(test0001, imp=15, cats, gtit="Coverage > 0.01")
plotPCoA(test0001[grep("Faecalibacterium",rownames(test0001)),], imp=15, cats, gtit="Coverage > 0.01")
dotsne(test0001, cats)

test001 = cov2relab(coverage, modstats, cats, pclass, mincoverage=0.01)
test01 = cov2relab(coverage, modstats, cats, pclass, mincoverage=0.1)
test1 = cov2relab(coverage, modstats, cats, pclass, mincoverage=1)
test2 = cov2relab(coverage, modstats, cats, pclass, mincoverage=2)
test5 = cov2relab(coverage, modstats, cats, pclass, mincoverage=5)
test10 = cov2relab(coverage, modstats, cats, pclass, mincoverage=10)
plotPCoA(test, imp=15, cats)

grid.arrange(plotPCoA(test001, imp=15, cats, gtit="Coverage > 0.01"), #8x8
             plotPCoA(test01, imp=15, cats, gtit="Coverage > 0.1"),
             plotPCoA(test1, imp=15, cats, gtit="Coverage > 1"),
             plotPCoA(test2, imp=15, cats, gtit="Coverage > 2"),
             plotPCoA(test5, imp=15, cats, gtit="Coverage > 5"),
             plotPCoA(test10, imp=15, cats, gtit="Coverage > 10"),
             ncol=3)#24x16

#sort(coverage[which(coverage[,"V1.CD.1"]>1),"V1.CD.1"],decreasing=T)
#mean(apply(coverage,2,function(x){length(x[which(x>1)])}))

#coverage["Abiotrophia_defectiva_ATCC_49176","V1.CD.1"]

##############################################################
########################### Export Data
##############################################################

setwd('P:/MAPPING/Results773/bwa_human_filtered/unfiltered')
coverage = as.matrix(read.csv("coverage773.csv",row.names=1,header=T))
cover = cov2relab(coverage, modstats, cats, pclass, mincoverage=1)

setwd('P:/MAPPING/Results773/bwa_human_filtered/filter1')
coverage = as.matrix(read.csv("coverage773.csv",row.names=1,header=T))
cover = cov2relab(coverage, modstats, cats, pclass, mincoverage=1)

#create random microbiotas
randab = matrix(0, nrow=nrow(cover), ncol=5) 
for(i in 1:ncol(randab)){
  set.seed(i)
  rands <- runif(nrow(cover),min=0,max=max(ifelse(is.na(coverage),1,0)))
  randab[,i] = rands/sum(rands)
}
colnames(randab) = paste(rep("R",ncol(randab)),1:ncol(randab),sep=".")
newdat = cbind(cover, randab)

write.csv(newdat,file="alldat773_sel_unfiltered.csv")

write.csv(newdat,file="alldat773_sel_filter1.csv")

##############################################################
########################### Plot Data
##############################################################

classes = dis[colnames(alldat)]
classes[which(is.na(classes))] = "random"
names(classes)[(ncol(cover)+1):ncol(alldat)] = colnames(alldat)[(ncol(cover)+1):ncol(alldat)]
write.csv(classes,file="P:/MAPPING/Model_data/class773_sel.csv")

dires <- capscale(t(newdat)~1, distance="bray")
dires_sum <- summary(dires)

imp=15
abseigen <- abs(dires_sum$species)
absort <- sort(abs(abseigen[,2]), decreasing=T)
sel <- names(head(absort, imp))

gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                   "Phenotype"=as.factor(c(dis[rownames(dires_sum$sites)])))
levels(gdat$Phenotype) = c("Crohn's disease", "Healthy", "Random")
gdat$Phenotype[which(is.na(gdat$Phenotype))] = factor("Random")

gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC1"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC2"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
gdat[which(gdat$Phenotype=="Random"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Random"),"PC1"])
gdat[which(gdat$Phenotype=="Random"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Random"),"PC2"])

fac = 2
ggplot(data=gdat, aes(x=PC1, y=PC2, colour=Phenotype)) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=1) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=1) +   
  stat_ellipse(level = 0.95) +
  
  geom_point(aes(x=PC1, y=PC2, shape=Phenotype), size=4) +
  geom_point(aes(x=mean.x, y=mean.y), shape=4, size=2, stroke=2) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=0.7) +
  xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
  ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.key=element_blank(),
        legend.title=element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=20)) +
  guides(colour = guide_legend(override.aes = list(shape = 15))) #width=12, height=8)

scatter3D(dires_sum$sites[,1],dires_sum$sites[,2],dires_sum$sites[,3],phi=10,pch=19,
          colvar=as.numeric(as.factor(classes[rownames(dires_sum$sites)])))
#text3D(dires_sum$sites[,1],dires_sum$sites[,2],dires_sum$sites[,3],
#       labels=rownames(dires_sum$sites),phi=10,add=T)

scatter3Drgl(dires_sum$sites[,1],dires_sum$sites[,2],dires_sum$sites[,3],phi=10,pch=19,
             colvar=as.numeric(as.factor(classes[rownames(dires_sum$sites)])))

library(RColorBrewer)
piedat = alldat
piedat = piedat[,1:(ncol(piedat)-3)]
piedat = piedat[which(apply(ifelse(piedat==0,0,1),1,sum)==ncol(piedat)),]
piedat = piedat[names(tail(sort(apply(piedat,1,sum)),10)),]
cols <- colorRampPalette(brewer.pal(8,"Set3"))(nrow(piedat))
win.metafile("P:/MAPPING/patient_pie.wmf", width = 100, height = 100)
op <- par(mfrow = c(4, 5))
for(i in 1:ncol(piedat)){
  par(lwd=25,cex=8)
  pie(piedat[,i],main=colnames(piedat)[i],labels=NA,col=cols)
}
dev.off()
