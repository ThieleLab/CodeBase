library(vegan)
library(ggplot2)
library(reshape)
library(gridExtra)
library(NMF)
library(scales)
library(RColorBrewer)
library(gridBase)
library(grid)
library(sybil)
require(e1071) # for svm()     
require(rgl) # for 3d graphics. 
library(plot3D)
library(Hmisc)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(cowplot)

setwd('P:/GitRep/BacArena')
source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

##################################################################################################################
################################################## Preparing the mapping results
##################################################################################################################

plotPCoA <- function(cover, pclass, dismeth="bray", tit="", legend=F){
  dis = as.character(pclass[,1])
  names(dis) = gsub("-",".",pclass[,2])
  dires <- capscale(t(cover)~1, distance=dismeth)
  dires_sum <- summary(dires)

  gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                     "Phenotype"=as.factor(dis[rownames(dires_sum$sites)]))
  levels(gdat$Phenotype) = c("Crohn's disease", "Healthy")
  gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC1"])
  gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
  gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC2"])
  gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
  gp <- ggplot(data=gdat, aes(x=PC1, y=PC2, color=Phenotype)) +
    geom_hline(aes(yintercept=0),color="grey", size=1.5) +
    geom_vline(aes(xintercept=0),color="grey", size=1.5) +
    stat_ellipse(level = 0.95, size=2) +
    scale_colour_manual(values = c("firebrick2","royalblue2")) +
    geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=1.8) +
    geom_point(aes(x=PC1, y=PC2),size=5) +
    xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
    ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
    theme_bw(base_size = 30) 
  if(tit==""){
    gp = gp
  }else{
    gp = gp + ggtitle(tit)
  }
  if(legend){
    gp <- gp + theme(legend.text=element_text(size=20),
          legend.justification = c(1, 0), 
          legend.position = c(0.9, 0.1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.key=element_blank(),
          legend.title=element_blank(),
          panel.border = element_rect(colour='black',size=3),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          #axis.ticks = element_line(size=1,color='black'),
          #axis.ticks.length=unit(0.3,'cm'),
          plot.title = element_text(hjust = 0.5)) +
      guides(colour = guide_legend(override.aes = list(shape = 15))) #width=12, height=8)
  }else{
    gp <- gp + theme(legend.text=element_text(size=20),
          legend.position="none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.key=element_blank(),
          legend.title=element_blank(),
          panel.border = element_rect(colour='black',size=3),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          #axis.ticks = element_line(size=1,color='black'),
          #axis.ticks.length=unit(0.3,'cm'),
          plot.title = element_text(hjust = 0.5)) +
      guides(colour = guide_legend(override.aes = list(shape = 15))) #width=12, height=8)
  }
  return(list("plot"=gp,"rtab"=dires_sum$species))
}
getFluxMet <- function(fluxlist_red,meti,meanmat_pop,modtax,subsys,ylab=F,tax){
  fluxlist_red = fluxlist_red[rownames(meanmat_pop)]
  fmat = matrix(0, nrow=nrow(meanmat_pop),ncol=ncol(meanmat_pop),
                dimnames=list(rownames(meanmat_pop),colnames(meanmat_pop)))
  for(i in 1:length(fluxlist_red)){
    for(j in 1:length(fluxlist_red[[i]])){
      if(meti %in% names(fluxlist_red[[i]][[j]])){
        fmat[names(fluxlist_red)[i],names(fluxlist_red[[i]])[j]] = fluxlist_red[[i]][[j]][meti]
      }
    }
  }
  pfmat = matrix(0, nrow=nlevels(modtax[,"Class"]), ncol=nrow(fmat),
                 dimnames=list(levels(modtax[,"Class"]), rownames(fmat)))
  for(i in rownames(pfmat)){
    if(class(fmat[,which(modtax[colnames(fmat),"Class"]==i)])!="numeric"){
      msum = apply(fmat[,which(modtax[colnames(fmat),"Class"]==i)],1,sum)
    }else{
      msum = fmat[,which(modtax[colnames(fmat),"Class"]==i)]
    }
    pfmat[i,names(msum)] = msum
  }
  if(nlevels(pclass$V2)>2){
    pfsm = matrix(0,ncol=3,nrow=nrow(pfmat))
    rownames(pfsm) = rownames(pfmat)
    for(i in 1:nrow(pfmat)){
      pfsm[i,1] = median(pfmat[i,which(pclass[colnames(pfmat),1] == "Crohn's disease")])
      #pfsm[i,2] = sd(pfmat[i,which(pclass[colnames(pfmat),1] == "Crohn's disease")])
      pfsm[i,2] = median(pfmat[i,which(pclass[colnames(pfmat),1] == "Healthy")])
      pfsm[i,2] = median(pfmat[i,which(pclass[colnames(pfmat),1] == "Treatment")])
    }
    colnames(pfsm) = c("Crohn's disease","Healthy","Treatment")
    pfsm = pfsm[tax,]
    pfdat = melt(pfsm)
    colnames(pfdat) = c("Class","Patient","Flux")
    #pfdat$Type = pclass[as.character(pfdat$Patient),1]
    pfdat$Met = meti
  }else{
    pfsm = matrix(0,ncol=2,nrow=nrow(pfmat))
    rownames(pfsm) = rownames(pfmat)
    for(i in 1:nrow(pfmat)){
      pfsm[i,1] = median(pfmat[i,which(pclass[colnames(pfmat),1] == "Crohn's disease")])
      #pfsm[i,2] = sd(pfmat[i,which(pclass[colnames(pfmat),1] == "Crohn's disease")])
      pfsm[i,2] = median(pfmat[i,which(pclass[colnames(pfmat),1] == "Healthy")])
      #pfsm[i,4] = sd(pfmat[i,which(pclass[colnames(pfmat),1] == "Healthy")])
    }
    colnames(pfsm) = c("Crohn's disease","Healthy")
    pfsm = pfsm[tax,]
    pfdat = melt(pfsm)
    colnames(pfdat) = c("Class","Patient","Flux")
    #pfdat$Type = pclass[as.character(pfdat$Patient),1]
    pfdat$Met = meti
  }
  return(pfdat)
}

setwd('P:/MAPPING/Pediatric_crohns')
coverage = as.matrix(read.csv("coverage773_pediatric_dysbiotic.csv",row.names=1,header=T))
cats = read.csv("P:/MAPPING/Model_data/dataset2_class.csv")

mincoverage=0.01
cats=catsel
rm=NULL
dis = as.character(cats$category)
names(dis) = gsub("-",".",as.character(cats$Sample.Name))
coverage = ifelse(is.na(coverage),0,coverage)
cover = ifelse(coverage<mincoverage,0,coverage)
cover = apply(cover,2,function(x){x/sum(x)})
colnames(cover) = gsub("-",".",colnames(cover))
if(!is.null(rm)){cover = cover[,-which(dis[colnames(cover)] == rm)]}

test0001 = cover
if(length(which(is.na(apply(test0001,2,sum))))!=0){test0001 = test0001[,-which(is.na(apply(test0001,2,sum)))]}
if(length(which(apply(test0001,1,sum)==0))!=0){test0001 = test0001[-which(apply(test0001,1,sum)==0),]}

##################################################################################################################
################################################## Preparing the simulation results
##################################################################################################################

cats = read.csv("P:/MAPPING/Model_data/class773.csv",row.names=1)
mrpass = t(read.csv("P:/AGORA/Revision3/rxn_pa.csv",header=T,row.names=1))
all = as.matrix(read.csv("P:/MAPPING/Model_data/alldat773.csv", header=T, row.names=1))
rownames(all)[which(rownames(all)=="Gemella_moribillum_M424")] = "Gemella_morbillorum_M424"
alldiet = read.csv("P:/MAPPING/Model_data/diets_eu.csv", header=T, row.names=1)
pclass = read.csv("P:/MAPPING/Model_data/class773.csv", header=F, row.names=1)
rpa = t(read.csv("P:/AGORA/Revision3/rxn_pa.csv",header=T,row.names=1))
meta = read.csv("P:/MAPPING/Pediatric_crohns/dysbiotic_cluster.csv",row.names=2)
gradmet = c(rownames(alldiet)[which(alldiet$type=="Mucin")],"EX_o2(e)","EX_no(e)")
modstats = read.csv("P:/AGORA/Revision3/ModelInformation.csv",header=T,row.names=1)
rownames(modstats) = modstats$ModelAGORA
rownames(modstats)[which(rownames(modstats)=="Ruminococcus_albus_7_DSM_20455")] = "Ruminococcus_albus_7"
rownames(modstats)[which(rownames(modstats)=="Bifidobacterium_thermacidophilum_subsp_thermacidophilum_DSM_158")] = "Bifidobacterium_thermacidophilum_subsp_thermacidophilum_DSM_15837"
rownames(modstats)[which(rownames(modstats)=="Gemella_moribillum_M424")] = "Gemella_morbillorum_M424"
subsys = read.csv("P:/AGORA/Revision3/TranslatedReactionDatabase_AGORA2.csv",header=T,row.names=2)
subsys$Subsystem = as.factor(ifelse(as.character(subsys$Subsystem)=="",'Unassigned',as.character(subsys$Subsystem)))

pclass[,1] = as.character(pclass[,1])

#for the pediatric chron's disease
cats = read.csv("P:/MAPPING/Pediatric_crohns/dysbiotic_cluster.csv",row.names=2)
catsel = data.frame(Reads=gsub(" S","S",rownames(cats)),Sample.Name=gsub(" S","S",rownames(cats)),category=cats$Disease)
catsel$category = as.character(catsel$category)
catsel$category[which(catsel$category == "Control")] = "Healthy"
catsel$category[which(catsel$category == "Crohn")] = "Crohn's disease"
catsel$category = as.factor(catsel$category)
pclass = data.frame(V2=catsel$category,V3=catsel$Reads)
rownames(pclass) = as.character(catsel$Reads)

setwd("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/")
load("sublist_eu.RData")
load("poplist_eu.RData")
#load("gralist_reduce_eu.RData")
#load("gralist_init_eu.RData")
load("fluxlist_red_eu.RData")
load("P:/MAPPING/Models773/agora773_constrain.RData")
names(poplist) = gsub(".mix","",names(poplist))
names(sublist) = gsub(".mix","",names(sublist))

##########################################################################
############### Figure 1
##########################################################################

# Figure S1
numic = apply(ifelse(cover!=0,1,0),2,sum)
nmdat = data.frame("patient"=names(numic),"numic"=numic,"type"=pclass[names(numic),"V2"])
nmdat$patient = factor(as.character(nmdat$patient),levels=as.character(nmdat$patient[order(nmdat$numic)]))
ggplot(nmdat, aes(x=patient, y=numic, fill=type)) +
  geom_bar(position="dodge", stat="identity") +
  #scale_y_continuous(limits = c(min(pfdat$Flux),max(pfdat$Flux)),breaks=seq(min(pfdat$Flux),max(pfdat$Flux),by=round(max(pfdat$Flux)/10000000,0)*10000)) +
  coord_flip() +
  scale_fill_manual(values = c("firebrick2","royalblue2")) +
  #ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
  xlab("") +
  ylab("Number of microbes") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        #legend.position="none",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        #axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        #axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_blank())#16x4
min(nmdat$numic)
max(nmdat$numic)

meanmat_pop = matrix(0,nrow=length(poplist),ncol=nrow(all),
                     dimnames=list(names(poplist),rownames(all)))
sdmat = matrix(0,nrow=length(poplist),ncol=nrow(all),
               dimnames=list(names(poplist),rownames(all)))
for(i in 1:length(poplist)){
  repmat = matrix(0,nrow=length(poplist[[i]]),ncol=nrow(poplist[[i]][[1]]))
  for(j in 1:length(poplist[[i]])){
    repmat[j,] = poplist[[i]][[j]][,ncol(poplist[[i]][[j]])]
  }
  meanmat_pop[i,rownames(poplist[[i]][[1]])] = apply(repmat,2,mean)
  sdmat[i,rownames(poplist[[i]][[1]])] = apply(repmat,2,sd)
}
mean(apply(meanmat_pop,1,sum))
sd(apply(meanmat_pop,1,sum))
mean(apply(meanmat_pop,1,sum)/10000)
sd(apply(meanmat_pop,1,sum)/10000)
meanmat_pop = meanmat_pop/apply(meanmat_pop,1,sum)

meanmat_sub = matrix(0,nrow=length(sublist),ncol=nrow(alldiet),
                     dimnames=list(names(sublist),rownames(alldiet)))
sdmat = matrix(0,nrow=length(sublist),ncol=nrow(alldiet),
               dimnames=list(names(sublist),rownames(alldiet)))
for(i in 1:length(sublist)){
  repmat = matrix(0,nrow=length(sublist[[i]]),ncol=nrow(sublist[[i]][[1]]))
  for(j in 1:length(sublist[[i]])){
    repmat[j,] = sublist[[i]][[j]][,ncol(sublist[[i]][[j]])]
  }
  meanmat_sub[i,rownames(sublist[[i]][[1]])] = apply(repmat,2,mean)
  sdmat[i,rownames(sublist[[i]][[1]])] = apply(repmat,2,sd)
}
meanmat_sub = ifelse(meanmat_sub<0,0,meanmat_sub)
#meanmat_sub = meanmat_sub/apply(meanmat_sub,1,sum)

psim = plotPCoA(t(meanmat_pop), pclass, tit="Simulated abundance")[[1]]
#psim = plotPCoA(t(meanmat_sub), pclass)
pmap = plotPCoA(cover, pclass, legend=T, tit="Mapped abundance")[[1]]

mean(apply(ifelse(cover!=0,1,0),2,sum))
sd(apply(ifelse(cover!=0,1,0),2,sum))

meanmat_pa = ifelse(meanmat_pop!=0,1,0)
allreact = unique(unlist(lapply(specs,react_id)))
reactmat = matrix(0, nrow=nrow(meanmat_pa), ncol=length(allreact), dimnames=list(rownames(meanmat_pa),allreact))
allreact = vector()
for(i in 1:nrow(meanmat_pa)){
  allreact[i] = length(unlist(lapply(specs[names(which(meanmat_pa[i,]==1))],react_id)))
  patreact = unique(unlist(lapply(specs[names(which(meanmat_pa[i,]==1))],react_id)))
  reactmat[i, patreact] = 1
}
prxn = plotPCoA(t(reactmat), pclass, dismeth="jaccard", tit="Reaction difference")[[1]]
rxntab = as.data.frame(plotPCoA(t(reactmat), pclass, dismeth="jaccard", tit="Reaction difference")[[2]][,1:2])
rxntab$Subsystem = subsys[rownames(rxntab),"Subsystem"]
rxntab = rxntab[-grep("biomass",rownames(rxntab)),]
rxntab = rxntab[order(abs(rxntab[,1]),decreasing=T),]
mean(allreact)
sd(allreact)
mean(apply(t(reactmat),2,sum))
sd(apply(t(reactmat),2,sum))

write.csv(rxntab,file="P:/MAPPING/manuscript/Supplementary_Table_S1.csv")

psm1 = grid.arrange(pmap,psim,prxn,ncol=3)#20x6


grid.arrange(plotPCoA(t(meanmat_pop), pclass, tit="Manhattan distance", dismeth="manhattan")[[1]],
             plotPCoA(t(meanmat_pop), pclass, tit="Euclidean distance", dismeth="euclidean")[[1]],ncol=2)#20x6

##########################################################################
############### Figure 1.2
##########################################################################

modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Genus, modstats$Species),]
relab = meanmat_pop/apply(meanmat_pop,1,sum)
pdat = data.frame()
pvals = vector()
#c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")
for(i in c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")){
  ab = relab[,intersect(rownames(modstats)[which(modstats$Class==i)],colnames(relab))]
  if(sum(ab)!=0){
    if(class(ab)!="numeric"){
      ab=apply(ab,1,sum)
    }
    pdat = rbind(pdat,data.frame(abundance=ab,patient=rownames(relab),
                                 type=as.character(pclass[rownames(relab),1]),phylum=rep(i,nrow(relab))))
    #pval = t.test(ab[which(pclass[names(ab),1]=="Healthy")],ab[which(pclass[names(ab),1]=="Crohn's disease")])$p.value
    pval = wilcox.test(ab[which(pclass[names(ab),1]=="Healthy")],ab[which(pclass[names(ab),1]=="Crohn's disease")],exact=F)$p.value
    #pval = ks.test(ab[which(pclass[names(ab),1]=="Healthy")],ab[which(pclass[names(ab),1]=="Crohn's disease")])$p.value
    pvals[i] = round(pval,5)
  }
}
pdat$phylum = gsub("Gammaproteobacteria","GProteobacteria",pdat$phylum)
pdat$phylum = factor(as.character(pdat$phylum),levels=rev(c("Clostridia","Bacteroidia","GProteobacteria","Bacilli")))
sp <- ggplot(pdat, aes(factor(phylum), abundance)) +
  geom_boxplot(aes(color=type),lwd=2, position=position_dodge(0.9),outlier.color=NA) +
  #coord_flip() +
  geom_line(data=data.frame(a=c(0.75,0.75,1.25,1.25),b=c(0.95,1,1,0.95)),aes(x=a,y=b)) +
  annotate("text", x=1, y=1.05, label=paste("p=",round(pvals["Bacilli"],3),sep=""), size=9) +
  geom_line(data=data.frame(a=c(1.75,1.75,2.25,2.25),b=c(0.95,1,1,0.95)),aes(x=a,y=b)) +
  annotate("text", x=2, y=1.05, label=paste("p=",round(pvals["Gammaproteobacteria"],3),sep=""), size=9) +
  geom_line(data=data.frame(a=c(2.75,2.75,3.25,3.25),b=c(0.55,0.6,0.6,0.55)),aes(x=a,y=b)) +
  annotate("text", x=3, y=0.65, label=paste("p<","0.001",sep=""), size=9) +
  geom_line(data=data.frame(a=c(3.75,3.75,4.25,4.25),b=c(0.55,0.6,0.6,0.55)),aes(x=a,y=b)) +
  annotate("text", x=4, y=0.65, label=paste("p<","0.001",sep=""), size=9) +
  xlab("") +
  ylab("Relative abundance") +
  scale_color_manual(values = c("firebrick2","royalblue2","darkorchid1")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="none",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))#7x9.5/14x8
sp

pdat = data.frame()
pvals = vector()
#selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
#selsel=c("EX_ala_L(e)","EX_arg_L(e)","EX_asp_L(e)","EX_cys_L(e)","EX_glu_L(e)","EX_gly(e)","EX_his_L(e)","EX_ile_L(e)","EX_leu_L(e)","EX_lys_L(e)","EX_met_L(e)","EX_phe_L(e)","EX_pro_L(e)","EX_ser_L(e)","EX_thr_L(e)","EX_trp_L(e)","EX_tyr_L(e)","EX_val_L(e)")
selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
for(i in rev(selsel)){
  #conc = meanmat_sub[,i]-initsub[i,]#/max(meanmat_sub[,i])
  conc = meanmat_sub[,i]
  metnam = gsub("exchange reaction for ","",gsub(" exchange","",subsys[i,"Reaction.Name"]))
  if(metnam == "Butyrate (n-C4:0)"){metnam = "Butyrate"}
  if(metnam == "H2"){metnam = "Hydrogen"}
  pdat = rbind(pdat,data.frame(conc=conc,patient=names(conc),type=as.character(pclass[names(conc),1]),met=rep(metnam,length(conc))))
  #pval = t.test(conc[which(pclass[names(conc),1]=="Healthy")],conc[which(pclass[names(conc),1]=="Crohn's disease")])$p.value
  pval = wilcox.test(conc[which(pclass[names(conc),1]=="Healthy")],conc[which(pclass[names(conc),1]=="Crohn's disease")],exact=F)$p.value
  #pval = ks.test(conc[which(pclass[names(conc),1]=="Healthy")],conc[which(pclass[names(conc),1]=="Crohn's disease")])$p.value
  #print(shapiro.test(conc[which(pclass[names(conc),1]=="Healthy")]))
  #print(shapiro.test(conc[which(pclass[names(conc),1]=="Crohn's disease")]))
  pvals[i] = pval             
}
#pdat$conc = pdat$conc/(10^12* 0.01 * 6.25e-08)
#pdat$conc = pdat$conc/max(pdat$conc)
met <- ggplot(pdat, aes(factor(met), conc)) +
  geom_boxplot(aes(color=type),lwd=2, position=position_dodge(0.9),outlier.color=NA) +
  scale_y_continuous(limits = c(0, 180)) +
  #scale_y_continuous(limits = c(0, 5)) +
  #coord_flip() +
  annotate("text", x=5, y=78, label=paste("p=",round(pvals["EX_ac(e)"],3),sep=""), size=9) +
  geom_line(data=data.frame(a=c(4.75,4.75,5.25,5.25),b=c(65,70,70,65)),aes(x=a,y=b)) +
  annotate("text", x=4, y=38, label=paste("p=",round(pvals["EX_lac_L(e)"],3),sep=""), size=9) +
  geom_line(data=data.frame(a=c(3.75,3.75,4.25,4.25),b=c(25,30,30,25)),aes(x=a,y=b)) +
  annotate("text", x=3, y=38, label=paste("p<","0.001",sep=""), size=9) +
  geom_line(data=data.frame(a=c(2.75,2.75,3.25,3.25),b=c(25,30,30,25)),aes(x=a,y=b)) +
  annotate("text", x=2, y=168, label=paste("p<","0.001",sep=""), size=9) +
  geom_line(data=data.frame(a=c(1.75,1.75,2.25,2.25),b=c(155,160,160,155)),aes(x=a,y=b)) +
  annotate("text", x=1, y=78, label=paste("p<","0.001",sep=""), size=9) +
  geom_line(data=data.frame(a=c(0.75,0.75,1.25,1.25),b=c(65,70,70,65)),aes(x=a,y=b)) +
xlab("") +
  ylab("Concentration in mM") +
  scale_color_manual(values = c("firebrick2","royalblue2","darkorchid1")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="none",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))#8x10/14x8
met


popflx = rbind(getFluxMet(fluxlist_red,"EX_but(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
               getFluxMet(fluxlist_red,"EX_ppa(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
               getFluxMet(fluxlist_red,"EX_isobut(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
               getFluxMet(fluxlist_red,"EX_lac_L(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
               getFluxMet(fluxlist_red,"EX_ac(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")))
popflx$Met[which(popflx$Met=="EX_but(e)")] = "Butyrate"
popflx$Met[which(popflx$Met=="EX_ppa(e)")] = "Propionate"
popflx$Met[which(popflx$Met=="EX_isobut(e)")] = "Isobutyrate"
popflx$Met[which(popflx$Met=="EX_lac_L(e)")] = "L-Lactate"
popflx$Met[which(popflx$Met=="EX_ac(e)")] = "Acetate"
popflx$Met = factor(popflx$Met,levels=c("Butyrate","Propionate","Isobutyrate","L-Lactate","Acetate"))
popflx$Class = gsub("Gammaproteobacteria","GProteobacteria",popflx$Class)
popflx$Class = factor(as.character(popflx$Class),levels=c("Clostridia","Bacteroidia","GProteobacteria","Bacilli"))
gflx = ggplot(popflx, aes(x=Class, y=(Flux/max(Flux))*100, fill=Patient)) +
  geom_bar(position="dodge", stat="identity") +
  #scale_y_continuous(limits = c(min(pfdat$Flux),max(pfdat$Flux)),breaks=seq(min(pfdat$Flux),max(pfdat$Flux),by=round(max(pfdat$Flux)/10000000,0)*10000)) +
  coord_flip() +
  facet_grid(~Met) +
  scale_fill_manual(values = c("firebrick2","royalblue2")) +
  #ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
  xlab("") +
  ylab("Relative population flux percentage") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="none",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        #axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=25),
        #axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_blank())#16x4

psm = grid.arrange(sp,met,ncol=2)#20x6
grid.arrange(psm1,psm,gflx,ncol=1,heights=c(5,5,3.5))#25x20

##########################################################################
############### Figure 2
##########################################################################

pdat = data.frame()

ratios = vector()
selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
for(i in rev(selsel)){
  #conc = meanmat_sub[,i]-initsub[i,]#/max(meanmat_sub[,i])
  conc = meanmat_sub[,i]
  metnam = gsub("exchange reaction for ","",gsub(" exchange","",subsys[i,"Reaction.Name"]))
  if(metnam == "Butyrate (n-C4:0)"){metnam = "Butyrate"}
  if(metnam == "H2"){metnam = "Hydrogen"}
  ratios[metnam] = mean(conc[which(pclass[names(conc),1]=="Healthy")])/mean(conc[which(pclass[names(conc),1]=="Crohn's disease")])
}

exp = read.csv("P:/MAPPING/literature/Validation/hove1995_data.csv",row.names=1)
comp = cbind(ratios,exp[names(ratios),"Mean"])
#comp = cbind(1/ratios,1/exp[names(ratios),"Mean"])
colnames(comp) = c("Simulation","Experiment")
compdat = melt(comp)
compdat$X1 = factor(as.character(compdat$X1),levels=c("Butyrate","Propionate","Isobutyrate","L-Lactate","Acetate"))
metcomp <- ggplot(compdat, aes(x=X1, y=value, fill=X2)) +
  geom_bar(width=0.6, stat="identity", position=position_dodge(), color="black") +
  scale_fill_manual(values = c("darkorange2","seashell3")) + #grey40
  geom_hline(aes(yintercept=1), linetype="dashed", color="black", size=2) +
  #ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
  ylab("Controls to CD patients ratio") +
  xlab("") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=25),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        #axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        #axis.title.x = element_text(size=25),
        #axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_blank())#12x7

exbac = read.csv("P:/MAPPING/literature/Validation/Lewis_bac.csv",row.names=1,header=T)
for(i in rownames(exbac)){
  cdm = t(cover)[rownames(pclass)[which(pclass[,1]=="Crohn's disease")],rownames(modstats)[which(as.character(modstats$Genus)==i)]]
  healthm = t(cover)[rownames(pclass)[which(pclass[,1]=="Healthy")],rownames(modstats)[which(as.character(modstats$Genus)==i)]]
  if(class(cdm)!="numeric"){cdm=apply(cdm,1,sum)}
  if(class(healthm)!="numeric"){healthm=apply(healthm,1,sum)}
  cd = meanmat_pop[rownames(pclass)[which(pclass[,1]=="Crohn's disease")],rownames(modstats)[which(as.character(modstats$Genus)==i)]]
  health = meanmat_pop[rownames(pclass)[which(pclass[,1]=="Healthy")],rownames(modstats)[which(as.character(modstats$Genus)==i)]]
  if(class(cd)!="numeric"){cd=apply(cd,1,sum)}
  if(class(health)!="numeric"){health=apply(health,1,sum)}
  exbac[i,"Health.map"] = median(healthm*100)
  exbac[i,"CD.map"] = median(cdm*100)
  exbac[i,"Health.sim"] = median(health*100)
  exbac[i,"CD.sim"] = median(cd*100)
}
mapcor = data.frame("Experiment"=c(exbac$Control,exbac$CD.far.cluster),"Mapped"=c(exbac$Health.map,exbac$CD.map))
ggplot(mapcor,aes(y=Experiment,x=Mapped)) +
  #scale_fill_manual(low = "grey",high = "black")+
  geom_point(size=3)  +
  #scale_y_log10() +
  #scale_x_log10() +
  xlab("Relative abundance based on reference based mapping") +
  ylab("Relative abundance based on reference independent analysis") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=25),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        #axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        #axis.title.x = element_text(size=25),
        #axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_blank())#12x7

rel = ifelse(cbind(exbac[,1]-exbac[,2],exbac[,3]-exbac[,4],exbac[,5]-exbac[,6])<0,"Higher in CD","Higher in control")
colnames(rel) = c("Experiment","Mapping","Simulation")
rownames(rel) = rownames(exbac)
rdat = melt(rel)
rdat$X2 = factor(as.character(rdat$X2),levels=c("Experiment","Mapping","Simulation"))
rdat$match = as.factor(rep(rel[,1]==rel[,3],3))
exbac[c("Bifidobacterium","Enterococcus","Lactobacillus","Veillonella"),]

spcomp <- ggplot(melt(rel),aes(y=factor(X2),x=factor(X1))) +
  geom_tile(aes(fill=rdat$match),color = "black",show.legend=F) +
  #scale_fill_manual(low = "grey",high = "black")+
  #geom_point(aes(colour=Growth, size=Growth))  +
  scale_fill_manual(values = c("grey90","white")) +
  scale_colour_manual(values = c("firebrick2","royalblue2")) +
  geom_point(aes(color=value),size=10) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(),
        legend.title=element_blank(),
        axis.title = element_blank(),
        panel.border = element_rect(colour='black',size=3),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_blank(),
        #axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=25,hjust = 0.5))

spcomp #20x6
metcomp #18x7
grid.arrange(spcomp,metcomp,ncol=1,heights=c(3.5,5))#25x20

##########################################################################
############### Figure 3
##########################################################################

selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
pconc = meanmat_sub[,selsel]#/(10^12*0.01*6.25e-08)

barplot(apply(meanmat_pop,1,diversity))

#mdat_red = meanmat_pop[as.character(pclass[which(pclass=="Crohn's disease"),2]),]
mdat_red = meanmat_pop
mdat_red = mdat_red[,-which(apply(mdat_red,2,sum)==0)]
mdat_red = ifelse(mdat_red==0,NA,mdat_red)
modstats2 = modstats[colnames(mdat_red),]
modtax <- modstats2[order(modstats2$Phylum, modstats2$Class, modstats2$Genus, modstats2$Species),]
mdat_red = mdat_red[,rownames(modtax)]

pclass$names=NA
pclass[which(pclass$V2=="Crohn's disease"),"names"] = paste("CD",1:length(which(pclass$V2=="Crohn's disease")),sep="")
pclass[which(pclass$V2=="Healthy"),"names"] = paste("HC",1:length(which(pclass$V2=="Healthy")),sep="")
catex = cats
catex$ID = pclass[rownames(catex),3]
write.csv(catex,file="P:/MAPPING/manuscript/TableS1.csv")

andf = data.frame("class"=modtax[colnames(mdat_red),c('Class')])
andf$class = as.character(andf$class)
rownames(andf) = colnames(mdat_red)
#andf[which(!(andf[,1] %in% c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli"))),] = "Other"
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
andf$class = as.factor(andf$class)
colpal = getPalette(nlevels(andf$class))
names(colpal) = levels(andf$class)
hatop = HeatmapAnnotation(df = andf, col = list(class=colpal))
haright = rowAnnotation(df = data.frame(type=pclass[rownames(mdat_red),c(1)]), 
                        col = list(type=c("Crohn's disease"="firebrick2","Healthy"="royalblue2")))
splitval = as.numeric(pclass[rownames(mdat_red),1])
mdat_red_plot = mdat_red
rownames(mdat_red_plot) = pclass[rownames(mdat_red_plot),"names"]
ha = Heatmap(mdat_red_plot,top_annotation=hatop, cluster_columns=F, show_row_names=T, show_column_names=F,
        col=colorRampPalette(c("seashell3","blue1","blue2","blue3","blue4"))(100), na_col="white",row_dend_width=unit(1, "cm"),
        top_annotation_height = unit(0.5, "cm"), row_title="Patients",row_title_side="left", split=splitval) + haright
getPalette = colorRampPalette(brewer.pal(5, "Dark2"))
metcol = getPalette(5)
rowdat = pconc[rownames(mdat_red),]

mcfluxsel = list()
for(j in selsel){
  fluxmat = matrix(0,nrow=length(fluxlist_red),ncol=ncol(mdat_red_plot),dimnames=list(names(fluxlist_red),colnames(mdat_red_plot)))
  for(i in rownames(fluxmat)){
    fluxmat[i,names(fluxlist_red[[i]])] = unlist(lapply(fluxlist_red[[i]],function(x){sum(x[intersect(names(x),j)])}))
  }
  mainc = c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli", "Other")
  mcflux = matrix(0,nrow=length(fluxlist_red),ncol=length(mainc),dimnames=list(names(fluxlist_red),mainc))
  for(i in mainc){
    if(i!="Other"){
      patflx = apply(fluxmat[,which(modstats[colnames(fluxmat),"Class"]==i)],1,sum)
      mcflux[names(patflx),i] = patflx
    }else{
      patflx = apply(fluxmat[,which(!(modstats[colnames(fluxmat),"Class"] %in% mainc))],1,sum)
      mcflux[names(patflx),i] = patflx
    }
  }
  mcfluxsel[[j]] = ifelse(mcflux<0,0,mcflux)
  mcfluxsel[[j]] = mcfluxsel[[j]]/apply(mcfluxsel[[j]],1,sum)
  mcfluxsel[[j]] = ifelse(is.na(mcfluxsel[[j]]),0,mcfluxsel[[j]])
}
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
colpal = getPalette(nlevels(andf$class))
names(colpal) = levels(andf$class)
Acetate = rowdat[,1]*mcfluxsel[["EX_ac(e)"]] 
hac = rowAnnotation(Acetate=row_anno_barplot(Acetate,
                                             gp=gpar(fill=c(colpal[c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")],"black"))),
                    annotation_width=unit(2.5,"cm"), show_annotation_name=T, annotation_name_rot=0)
Lactate = rowdat[,2]*mcfluxsel[["EX_lac_L(e)"]] 
hla = rowAnnotation(Lactate=row_anno_barplot(Lactate,
                                             gp=gpar(fill=c(colpal[c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")],"black"))),
                    annotation_width=unit(2.5,"cm"), show_annotation_name=T, annotation_name_rot=0)
Isobutyrate = rowdat[,3]*mcfluxsel[["EX_isobut(e)"]] 
his = rowAnnotation(Isobutyrate=row_anno_barplot(Isobutyrate,
                                                 gp=gpar(fill=c(colpal[c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")],"black"))),
                    annotation_width=unit(2.5,"cm"), show_annotation_name=T, annotation_name_rot=0)
Propionate = rowdat[,4]*mcfluxsel[["EX_ppa(e)"]] 
hpp = rowAnnotation(Propionate=row_anno_barplot(Propionate,
                                                gp=gpar(fill=c(colpal[c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")],"black"))),
                    annotation_width=unit(2.5,"cm"), show_annotation_name=T, annotation_name_rot=0)
Butyrate = rowdat[,5]*mcfluxsel[["EX_but(e)"]] 
hbu = rowAnnotation(Butyrate=row_anno_barplot(Butyrate,
                                              gp=gpar(fill=c(colpal[c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")],"black"))),
                    annotation_width=unit(2.5,"cm"), show_annotation_name=T, annotation_name_rot=0)
ht = ha+hbu+hpp+his+hla+hac
draw(ht, padding = unit(c(8, 2, 2, 2), "mm")) #16x9


#some stats
Acetate = rowdat[,1]
Lacatate = rowdat[,2]
Isobutyrate = rowdat[,3]
Propionate = rowdat[,4]
Butyrate = rowdat[,5]
scfas = Acetate+Isobutyrate+Propionate+Butyrate
mean(scfas[rownames(pclass)[which(pclass$V2=="Crohn's disease")]])
mean(scfas[rownames(pclass)[which(pclass$V2=="Healthy")]])

mean(Butyrate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]])
mean(Butyrate[rownames(pclass)[which(pclass$V2=="Healthy")]])
mean(Isobutyrate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]])
mean(Isobutyrate[rownames(pclass)[which(pclass$V2=="Healthy")]])
mean(Propionate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]])
mean(Propionate[rownames(pclass)[which(pclass$V2=="Healthy")]])
mean(Lactate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]])
mean(Lactate[rownames(pclass)[which(pclass$V2=="Healthy")]])
mean(Acetate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]])
mean(Acetate[rownames(pclass)[which(pclass$V2=="Healthy")]])
length(which(Butyrate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Butyrate[rownames(pclass)[which(pclass$V2=="Healthy")]])))
pclass[which(Butyrate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Butyrate[rownames(pclass)[which(pclass$V2=="Healthy")]])),3]
length(which(Propionate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Propionate[rownames(pclass)[which(pclass$V2=="Healthy")]])))
pclass[which(Propionate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Propionate[rownames(pclass)[which(pclass$V2=="Healthy")]])),3]
length(which(Isobutyrate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Isobutyrate[rownames(pclass)[which(pclass$V2=="Healthy")]])))
pclass[which(Isobutyrate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Isobutyrate[rownames(pclass)[which(pclass$V2=="Healthy")]])),3]
length(which(Lactate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Lactate[rownames(pclass)[which(pclass$V2=="Healthy")]])))

length(which(Acetate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Acetate[rownames(pclass)[which(pclass$V2=="Healthy")]])))
pclass[which(Acetate[rownames(pclass)[which(pclass$V2=="Crohn's disease")]]>=mean(Acetate[rownames(pclass)[which(pclass$V2=="Healthy")]])),3]

fluxmat = matrix(NA,nrow=length(fluxlist_red),ncol=ncol(mdat_red_plot),dimnames=list(names(fluxlist_red),colnames(mdat_red_plot)))
for(i in rownames(fluxmat)){
  fluxmat[i,names(fluxlist_red[[i]])] = unlist(lapply(fluxlist_red[[i]],function(x){sum(x[intersect(names(x),selsel)])}))
}
fluxmat = ifelse(fluxmat<0,0,fluxmat)
ha = Heatmap(fluxmat,top_annotation=hatop, cluster_columns=F, show_row_names=T, show_column_names=F,
             col=colorRampPalette(c("seashell3","red1","red2","red3","red4"))(100), na_col="white",row_dend_width=unit(1, "cm"),
             top_annotation_height = unit(0.5, "cm"), row_title="Patients",row_title_side="left", split=splitval) + haright

plot(apply(ifelse(is.na(fluxmat),0,fluxmat),2,sum),
     apply(ifelse(is.na(mdat_red_plot),0,mdat_red_plot),2,sum))

##########################################################################
############### Figure 4
##########################################################################

pconc = meanmat_sub[,selsel]
pcd = pconc[which(pclass[rownames(pconc),1]=="Crohn's disease"),]
phe = pconc[which(pclass[rownames(pconc),1]=="Healthy"),]
hconc = apply(phe,2,mean)
probs = matrix(T,nrow=nrow(pcd),ncol=ncol(pcd),dimnames=list(rownames(pcd),colnames(pcd)))
for(i in rownames(pcd)){
  if(pcd[i,"EX_ac(e)"] >= hconc["EX_ac(e)"]){probs[i,"EX_ac(e)"]=F}
  if(pcd[i,"EX_lac_L(e)"] <= hconc["EX_lac_L(e)"]){probs[i,"EX_lac_L(e)"]=F}
  if(pcd[i,"EX_isobut(e)"] >= hconc["EX_isobut(e)"]){probs[i,"EX_isobut(e)"]=F}
  if(pcd[i,"EX_ppa(e)"] >= hconc["EX_ppa(e)"]){probs[i,"EX_ppa(e)"]=F}
  if(pcd[i,"EX_but(e)"] >= hconc["EX_but(e)"]){probs[i,"EX_but(e)"]=F}
}

sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

rxnmat = matrix(0,nrow=length(specs),ncol=ncol(meanmat_sub),dimnames=list(names(specs),colnames(meanmat_sub)))
obs = vector()
for(i in 1:length(specs)){
  med = findExchReact(specs[[i]])@react_id[grep("EX",findExchReact(specs[[i]])@react_id)]
  mod = changeBounds(specs[[i]], med, lb=-1)
  rxnmat[names(specs)[i],med]=1
  if("EX_isobut(e)" %in% med){
    mod = changeObjFunc(mod,"EX_isobut(e)")
    obs[i] = optimizeProb(mod,retOptSol=F)$obj
  }
}
hist(obs)

patient = rownames(pclass)[which(pclass$V2=="Crohn's disease")]
patlist = list()
#selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
for(i in patient){
  patmic = rownames(poplist[[i]][[1]])
  patmat = matrix(0, length(patmic), ncol=length(selsel), dimnames=list(patmic,selsel))
  for(j in patmic){
    med = findExchReact(specs[[j]])@react_id[grep("EX",findExchReact(specs[[j]])@react_id)]
    mod = changeBounds(specs[[j]], med, lb=-1)
    #sol = optimizeProb(mod,retOptSol=F)
    #patmat[j,"bio"] = sol$obj
    for(k in selsel){
      if(k %in% med){
        modex = changeObjFunc(mod,k)
        if(k == "EX_lac_L(e)"){
          patmat[j,k] = optimizeProb(modex,retOptSol=F,lpdir="min")$obj
        }else{
          patmat[j,k] = optimizeProb(modex,retOptSol=F)$obj
        }
      }
    }
  }
  patlist[[i]] = patmat
}
patlist_origin = patlist

patient = rownames(pclass)[which(pclass$V2=="Crohn's disease")]
patlist = list()
selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
for(i in patient){
  patmic = rownames(poplist[[i]][[1]])
  gbtest = which(modstats[patmic,"Class"] %in% c("Bacilli","Gammaproteobacteria")) 
  bctest = which(modstats[patmic,"Class"] %in% c("Bacteroidia","Clostridia")) 
  if(length(gbtest)!=0){if(length(gbtest)!=1){gb=apply(rxnmat[patmic[gbtest],],2,sum)}else{gb=patmic[gbtest]}}
  if(length(bctest)!=0){if(length(bctest)!=1){bc=apply(rxnmat[patmic[bctest],],2,sum)}else{bc=patmic[bctest]}}
  if(length(gbtest)*length(bctest)!=0){metest = setdiff(names(bc)[which(bc!=0)],names(gb)[which(gb!=0)])}else{metest = NA}
  if(is.null(metest) || length(metest)==0){metest = NA}
  if(!is.na(metest)){metest=setdiff(metest,selsel)}
  if(is.null(metest) || length(metest)==0){metest = NA}
  #patmat = matrix(0, length(patmic), ncol=length(selsel)+1, dimnames=list(patmic,c("bio",selsel)))
  miclist = list()
  for(j in patmic){
    med = findExchReact(specs[[j]])@react_id[grep("EX",findExchReact(specs[[j]])@react_id)]
    mod = changeBounds(specs[[j]], med, lb=-1)
    patmat = matrix(0, length(intersect(metest,med)), ncol=length(selsel), dimnames=list(intersect(metest,med),selsel))
    if(!is.na(metest)){
      for(l in intersect(metest,med)){
        modi = changeBounds(mod, l, lb=-1000)
        for(k in selsel){
          if(k %in% med && probs[i,k]){
            modex = changeObjFunc(modi,k)
            if(k == "EX_lac_L(e)"){
              patmat[l,k] = optimizeProb(modex,retOptSol=F,lpdir="min")$obj
            }else{
              patmat[l,k] = optimizeProb(modex,retOptSol=F)$obj
            }
          }
        }
      }
    }
    miclist[[j]] = patmat
  }
  patlist[[i]] = miclist
}

tol=0.1
tlist = list()
for(i in names(patlist)){
  tlist[[i]] = vector()
  for(j in names(patlist[[i]])){
    if(nrow(patlist[[i]][[j]])!=0){
      for(k in rownames(patlist[[i]][[j]])){
        if(patlist[[i]][[j]][k,"EX_ac(e)"]-patlist_origin[[i]][j,"EX_ac(e)"] >= tol){
          tlist[[i]] = unique(c(tlist[[i]],k))
        }
        if(patlist[[i]][[j]][k,"EX_ppa(e)"]-patlist_origin[[i]][j,"EX_ppa(e)"] >= tol){
          tlist[[i]] = unique(c(tlist[[i]],k))
        }
        if(patlist[[i]][[j]][k,"EX_but(e)"]-patlist_origin[[i]][j,"EX_but(e)"] >= tol){
          tlist[[i]] = unique(c(tlist[[i]],k))
        }
        if(patlist[[i]][[j]][k,"EX_isobut(e)"]-patlist_origin[[i]][j,"EX_isobut(e)"] >= tol){
          tlist[[i]] = unique(c(tlist[[i]],k))
        }
        if(patlist[[i]][[j]][k,"EX_lac_L(e)"]-patlist_origin[[i]][j,"EX_lac_L(e)"] <= -tol){
          tlist[[i]] = unique(c(tlist[[i]],k))
        }
      }
    }
  }
}
length(which(unlist(lapply(tlist,length))==0))

save(tlist,file="P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/mtreatment.RData")
#############################################################
###############################################################################################################################################
#############################################################
load("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/mtreatment.RData")

tmet = sort(table(unlist(tlist)))
subs = gsub("exchange reaction for ","",gsub(" exchange","",subsys[names(tmet),"Reaction.Name"]))
tdat = melt(tmet)
tdat$indices = factor(as.character(tdat$Var.1),levels=tdat$Var.1)
tdat$Var.1 = factor(as.character(subs),levels=subs)

barall <- ggplot(tdat, aes(x=Var.1, y=value)) +
  geom_bar(width=0.6, stat="identity", position=position_dodge(), color="black") +
  #ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
  ylab("Number of patients") +
  xlab("") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=25),
        legend.position="right",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        #axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=15,angle = 45, hjust = 1),
        #axis.title.x = element_text(size=25),
        #axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        axis.ticks.x=element_blank(),
        plot.title = element_blank())#22x11

#write.csv(tdat,file="P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/tdat.csv")

tcats = read.csv("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/tdat_cats.csv")
rownames(tcats) = tcats$indices
catlist = lapply(tlist,function(x){as.character(tcats[x,"Atype"])})
tmet = sort(table(unlist(catlist)))
tdat = melt(tmet)
tdat$indices = factor(as.character(tdat$Var.1),levels=tdat$Var.1)
tdat$value = tdat$value/sum(tdat$value)
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
gpie <- ggplot(tdat, aes(x="", y=value, fill=indices)) + #do a pie chart here
  geom_bar(width=1, stat="identity", color="white", size=2) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=getPalette(10)) +
  ggtitle("Metabolite distribution") +
  ylab("") +
  xlab("") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=25),
        legend.position="none",
        legend.key=element_blank(),
        #legend.title=element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        #axis.text.y = element_text(size=25),
        axis.text.x = element_blank(),
        #axis.title.x = element_text(size=25),
        #axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5))#22x9

catmat = matrix(0,nrow=length(tlist),ncol=length(unique(unlist(catlist))),dimnames=list(names(tlist),unique(unlist(catlist))))
for(i in names(tlist)){
  cats = table(as.character(tcats[tlist[[i]],"Atype"]))
  if(length(cats)!=0){
    catmat[i,names(cats)] = cats
  }
}
catdato = melt(catmat)
rownames(catmat) = pclass[rownames(catmat),"names"]
catdat = melt(catmat)
catdat$X1 = factor(as.character(catdat$X1),levels=names(sort(apply(catmat,1,sum))))
catdat$X2 = factor(as.character(catdat$X2),levels=c("Carbohydrate metabolism","Central metabolism","Fermentation pathways","Host polysaccharides","Mucins",                 
                                                    "Phenylalanine metabolism","Plant polysaccharides","Respiration","TCA cycle","Vitamin metabolism"))
gpat <- ggplot(catdat, aes(x=X1, y=value, fill=X2)) +
  geom_bar(width=0.8, stat="identity", position="stack", color="black") +
  coord_flip() +
  scale_fill_manual(values=getPalette(10)) +
  scale_x_discrete(position = "top", expand = c(0.01, 0.01)) +
  scale_y_continuous(expand = c(0.01, 0.01)) +
  #ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
  ylab("Number of metabolites") +
  #xlab("Patients") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=15),
        #legend.position="none",
        legend.justification = c(1, 0), 
        legend.position = c(0.9, 0.1),
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(size=25),
        #axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        #axis.ticks.y=element_blank(),
        plot.title = element_blank())#22x9
gpat
grid.arrange(gpat,gpie,ncol=2)#22x8
pclass2 = pclass

apply(ifelse(meanmat_pop!=0,1,0)[names(which(apply(catmat,1,sum)==0)),],1,sum)
nonresponder = ifelse(meanmat_pop!=0,1,0)[names(which(apply(catmat,1,sum)==0)),]
sort(table(modstats[names(which(nonresponder[1,]==1)),"Class"]))
sort(table(modstats[names(which(nonresponder[2,]==1)),"Class"]))
sort(table(modstats[names(which(nonresponder[3,]==1)),"Class"]))
sort(table(modstats[names(which(nonresponder[4,]==1)),"Class"]))


# sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")
# 
# obs = vector()
# for(i in 1:length(specs)){
#   med = findExchReact(specs[[i]])@react_id[grep("EX",findExchReact(specs[[i]])@react_id)]
#   mod = changeBounds(specs[[i]], med, lb=-10)
#   obs[i] = optimizeProb(mod,retOptSol=F)$obj
# }
# hist(obs)
# 
# patient = rownames(pclass)[which(pclass$V2=="Crohn's disease")]
# patlist = list()
# #selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
# selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)","EX_ala_L(e)","EX_asp_L(e)","EX_arg_L(e)","EX_gly(e)","EX_met_L(e)","EX_trp_L(e)")
# for(i in patient){
#   patmic = rownames(poplist[[i]][[1]])
#   patmat = matrix(0, length(patmic), ncol=length(selsel)+1, dimnames=list(patmic,c("bio",selsel)))
#   for(j in patmic){
#     med = findExchReact(specs[[j]])@react_id[grep("EX",findExchReact(specs[[j]])@react_id)]
#     mod = changeBounds(specs[[j]], med, lb=-10)
#     sol = optimizeProb(mod,retOptSol=F)
#     patmat[j,"bio"] = sol$obj
#     for(k in selsel){
#       if(k %in% mod@react_id){patmat[j,k] = sol$fluxes[which(mod@react_id==k)]}
#     }
#   }
#   patlist[[i]] = patmat
# }
# patlist_origin = patlist
# 
# patient = rownames(pclass)[which(pclass$V2=="Crohn's disease")]
# patlist = list()
# selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)","EX_ala_L(e)","EX_asp_L(e)","EX_arg_L(e)","EX_gly(e)","EX_met_L(e)","EX_trp_L(e)")
# for(i in patient){
#   patmic = rownames(poplist[[i]][[1]])
#   #patmat = matrix(0, length(patmic), ncol=length(selsel)+1, dimnames=list(patmic,c("bio",selsel)))
#   miclist = list()
#   for(j in patmic){
#     med = findExchReact(specs[[j]])@react_id[grep("EX",findExchReact(specs[[j]])@react_id)]
#     mod = changeBounds(specs[[j]], med, lb=-10)
#     patmat = matrix(0, length(med), ncol=length(selsel)+1, dimnames=list(med,c("bio",selsel)))
#     for(l in med){
#       modi = changeBounds(mod, l, lb=-1000)
#       sol = optimizeProb(modi,retOptSol=F)
#       patmat[l,"bio"] = sol$obj
#       for(k in selsel){
#         if(k %in% modi@react_id){patmat[l,k] = sol$fluxes[which(modi@react_id==k)]}
#       }
#     }
#     miclist[[j]] = patmat
#   }
#   patlist[[i]] = miclist
# }
# save(patlist,file="P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/newest/patients_mtreatment.RData")
# save(patlist_origin,file="P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/newest/patients_original.RData")
# 
# ################
# 
# load("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/newest/patients_mtreatment.RData")
# load("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/newest/patients_original.RData")
# 
# orgmat = do.call(rbind,lapply(patlist_origin,function(x){apply(x,2,sum)}))
# 
# btreatment <- lapply(rownames(orgmat),function(x){
#   matrix(0,nrow=ncol(meanmat_sub),ncol=ncol(meanmat_pop),dimnames=list(colnames(meanmat_sub),colnames(meanmat_pop)))
# })
# names(btreatment) = names(patlist)
# bctreat = list()
# for(i in names(patlist)){
#   for(j in names(patlist[[i]])){
#     btreatment[[i]][names(patlist[[i]][[j]][,1]),j] = patlist[[i]][[j]][,1] - patlist_origin[[i]][j,1]
#     bctreat[[i]] = cbind(apply(btreatment[[i]][,which(modstats[colnames(btreatment[[i]]),"Class"]=="Bacilli")],1,sum),
#           apply(btreatment[[i]][,which(modstats[colnames(btreatment[[i]]),"Class"]=="Gammaproteobacteria")],1,sum),
#           apply(btreatment[[i]][,which(modstats[colnames(btreatment[[i]]),"Class"]=="Bacteroidia")],1,sum),
#           apply(btreatment[[i]][,which(modstats[colnames(btreatment[[i]]),"Class"]=="Clostridia")],1,sum))
#     colnames(bctreat[[i]]) = c("Bacilli","Gammaproteobacteria","Bacteroidia","Clostridia")
#   }
# }
# 
# findMetBio = function(pdat,tol,tabu){
#   target = intersect(intersect(intersect(which(pdat[,"Bacilli"]<=-tol),
#                                                    which(pdat[,"Gammaproteobacteria"]<=-tol)),
#                                          which(pdat[,"Bacteroidia"]>=tol)),
#                                which(pdat[,"Clostridia"]>=tol))
#   target = rownames(pdat)[target]
#   test = which(target %in% tabu)
#   if(length(test)!=0){target=target[-test]}
#   test2 = which(apply(abs(pdat[target,]),1,sum)==0)
#   if(length(test2)!=0){target=target[-test2]}
#   return(target)
# }
# 
# tmetbio = list()
# for(i in names(bctreat)){
#   pdat = as.data.frame(bctreat[[i]])
#   target = findMetBio(pdat,0,c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)","EX_biomass(e)"))
#   tmetbio[[i]] = target
# }
# #pdat[tmetbio[[i]],]
# 
# #compare simple model with predictions of big model
# #plot(pconc[rownames(orgmat),1],orgmat[,2])
# #plot(pconc[rownames(orgmat),2],orgmat[,3])
# #plot(pconc[rownames(orgmat),3],orgmat[,4])
# #plot(pconc[rownames(orgmat),4],orgmat[,5])
# #plot(pconc[rownames(orgmat),5],orgmat[,6])
# 
# treatmet <- lapply(colnames(meanmat_sub),function(x){
#   matrix(0,nrow=nrow(orgmat),ncol=ncol(orgmat),dimnames=list(rownames(orgmat),colnames(orgmat)))
#   })
# names(treatmet) = colnames(meanmat_sub)
# for(i in names(patlist)){
#   for(j in colnames(orgmat)){
#     for(k in names(patlist[[i]])){
#       for(l in rownames(patlist[[i]][[k]])){
#         treatmet[[l]][i,j] = treatmet[[l]][i,j] + patlist[[i]][[k]][l,j]
#       }
#     }
#   }
# }
# 
# ptreat <- lapply(rownames(treatmet[[1]]),function(x){
#   matrix(0,nrow=ncol(treatmet[[1]]),ncol=length(treatmet),dimnames=list(colnames(treatmet[[1]]),names(treatmet)))
# })
# names(ptreat) = rownames(orgmat)
# for(i in names(treatmet)){
#   dif = treatmet[[i]] - orgmat
#   for(j in rownames(dif)){
#     ptreat[[j]][,i] = dif[j,]
#   }
# }
# ##########################
# 
# findMet = function(pdat,tol,tabu){
#   target = intersect(intersect(intersect(which(pdat[,"EX_lac_L(e)"]<=tol),
#                                 which(pdat[,"EX_isobut(e)"]>=-tol)),
#                       which(pdat[,"EX_ppa(e)"]>=-tol)),
#             which(pdat[,"EX_but(e)"]>=-tol))
#   target = rownames(pdat)[target]
#   test = which(target %in% tabu)
#   if(length(test)!=0){target=target[-test]}
#   return(target)
# }
# 
# tmet = list()
# for(i in names(ptreat)){
#   pdat = as.data.frame(t(ptreat[[i]]))
#   target = findMet(pdat,0,c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)","EX_biomass(e)"))
#   if(length(target)==0){
#     target = findMet(pdat,1,c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)","EX_biomass(e)"))
#     if(length(target)==0){
#       target = findMet(pdat,10,c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)","EX_biomass(e)"))
#       if(length(target)==0){
#         target = findMet(pdat,100,c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)","EX_biomass(e)"))
#       }
#     }
#   }
#   tmet[[i]] = target
# }
# 
# treatall = list()
# for(i in names(tmet)){
#   treatall[[i]] = union(tmet[[i]],tmetbio[[i]])
# }
# 
# save(treatall,file="P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/treatment/mtreatment.RData")

##########################################################################
############### Figure 5
##########################################################################

#load in treatment
setwd("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/rich/high_fin/FINAL_TEST/treat/FINAL")
load("sublist_eu.RData")
load("poplist_eu.RData")

meanmat_popb = meanmat_pop
meanmat_subb = meanmat_sub

meanmat_pop = matrix(0,nrow=length(poplist),ncol=nrow(all),
                     dimnames=list(names(poplist),rownames(all)))
sdmat = matrix(0,nrow=length(poplist),ncol=nrow(all),
               dimnames=list(names(poplist),rownames(all)))
for(i in 1:length(poplist)){
  repmat = matrix(0,nrow=length(poplist[[i]]),ncol=nrow(poplist[[i]][[1]]))
  for(j in 1:length(poplist[[i]])){
    repmat[j,] = poplist[[i]][[j]][,ncol(poplist[[i]][[j]])]
  }
  meanmat_pop[i,rownames(poplist[[i]][[1]])] = apply(repmat,2,mean)
  sdmat[i,rownames(poplist[[i]][[1]])] = apply(repmat,2,sd)
}
meanmat_pop = meanmat_pop/apply(meanmat_pop,1,sum)

meanmat_sub = matrix(0,nrow=length(sublist),ncol=nrow(alldiet),
                     dimnames=list(names(sublist),rownames(alldiet)))
sdmat = matrix(0,nrow=length(sublist),ncol=nrow(alldiet),
               dimnames=list(names(sublist),rownames(alldiet)))
for(i in 1:length(sublist)){
  repmat = matrix(0,nrow=length(sublist[[i]]),ncol=nrow(sublist[[i]][[1]]))
  for(j in 1:length(sublist[[i]])){
    repmat[j,] = sublist[[i]][[j]][,ncol(sublist[[i]][[j]])]
  }
  meanmat_sub[i,rownames(sublist[[i]][[1]])] = apply(repmat,2,mean)
  sdmat[i,rownames(sublist[[i]][[1]])] = apply(repmat,2,sd)
}
meanmat_sub = ifelse(meanmat_sub<0,0,meanmat_sub)
#meanmat_sub = meanmat_sub/apply(meanmat_sub,1,sum)

rownames(meanmat_sub) = paste(rownames(meanmat_sub),"treat",sep="_")
rownames(meanmat_pop) = paste(rownames(meanmat_pop),"treat",sep="_")

meanmat_pop_all = rbind(meanmat_popb,meanmat_pop)
meanmat_sub_all = rbind(meanmat_subb,meanmat_sub)

pclass = rbind(pclass[,c(1,2)],data.frame(V2="Treatment",V3=rownames(meanmat_pop),row.names=rownames(meanmat_pop)))

cdmet = meanmat_sub_all[rownames(pclass)[which(pclass$V2 == "Crohn's disease")],selsel]
treatmet = meanmat_sub_all[rownames(pclass)[which(pclass$V2 == "Treatment")],selsel][paste("treat",rownames(cdmet),"_treat",sep=""),]

response = matrix("No response",nrow=nrow(cdmet),ncol=ncol(cdmet),dimnames=list(rownames(cdmet),colnames(cdmet)))
for(i in 1:ncol(treatmet)){
  if(colnames(treatmet)[i]!="EX_lac_L(e)"){
    ratio = cdmet[,i]/treatmet[,i]
  }else{
    ratio = treatmet[,i]/cdmet[,i]
  }
  response[which(ratio<0.7),i] = "Response"
}
colnames(response) = c("Acetate","L-Lactate","Isobutyrate","Propionate","Butyrate")
rownames(response) = pclass2[rownames(response),"names"]
rdat = melt(response)
rdat$X1 = factor(as.character(rdat$X1),levels=levels(catdat$X1))
rdat$X2 = factor(as.character(rdat$X2),levels=rev(c("Acetate","L-Lactate","Isobutyrate","Propionate","Butyrate")))
rpat = ggplot(rdat,aes(y=factor(X1),x=factor(X2))) +
  geom_tile(aes(fill=value),color =NA,show.legend=F) +
  #scale_fill_manual(low = "grey",high = "black")+
  #geom_point(aes(colour=Growth, size=Growth))  +
  xlab("") +
  scale_fill_manual(values = c("white","darkorchid1")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key=element_blank(),
        legend.title=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour='black',size=3),
        axis.text.x = element_text(size=25),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        #axis.ticks.y = element_blank(),
        #axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=25,hjust = 0.5))

pan1 = grid.arrange(gpat,rpat,ncol=2, widths=c(5,3.5))#22x8



modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Genus, modstats$Species),]
relab = meanmat_pop_all
pdat = data.frame()
pvals = vector()
pvals2 = vector()
#c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")
for(i in c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")){
  ab = relab[,intersect(rownames(modstats)[which(modstats$Class==i)],colnames(relab))]
  if(sum(ab)!=0){
    if(class(ab)!="numeric"){
      ab=apply(ab,1,sum)
    }
    pdat = rbind(pdat,data.frame(abundance=ab,patient=rownames(relab),
                                 type=as.character(pclass[rownames(relab),1]),phylum=rep(i,nrow(relab))))
    pval = wilcox.test(ab[which(pclass[names(ab),1]=="Treatment")],ab[which(pclass[names(ab),1]=="Crohn's disease")],exact=F)$p.value
    pvals[i] = pval  
    pval = wilcox.test(ab[which(pclass[names(ab),1]=="Treatment")],ab[which(pclass[names(ab),1]=="Healthy")],exact=F)$p.value
    pvals2[i] = pval  
  }
}
pdat$phylum = gsub("Gammaproteobacteria","GProteobacteria",pdat$phylum)
pdat$phylum = factor(as.character(pdat$phylum),levels=rev(c("Clostridia","Bacteroidia","GProteobacteria","Bacilli")))
pdat$type = factor(as.character(pdat$type),levels=c("Crohn's disease","Treatment","Healthy"))
sp <- ggplot(pdat, aes(factor(phylum), abundance)) +
  geom_boxplot(aes(color=type),lwd=2,position=position_dodge(0.9),outlier.color=NA) +
  #coord_flip() +
  geom_line(data=data.frame(a=c(0.7,0.7,1,1),b=c(0.96,0.99,0.99,0.96)),aes(x=a,y=b)) + #0.1
  annotate("text", x=0.85, y=1.02, label=paste("p=",round(pvals["Bacilli"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(1,1,1.3,1.3),b=c(0.96,0.99,0.99,0.96)),aes(x=a,y=b)) + #0.1
  annotate("text", x=1.15, y=1.02, label=paste("p=",round(pvals2["Bacilli"],3),sep=""), size=5) +
  
  geom_line(data=data.frame(a=c(1.7,1.7,2,2),b=c(0.96,0.99,0.99,0.96)),aes(x=a,y=b)) + #0.1
  annotate("text", x=1.85, y=1.02, label=paste("p=",round(pvals["Gammaproteobacteria"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(2,2,2.3,2.3),b=c(0.96,0.99,0.99,0.96)),aes(x=a,y=b)) + #0.1
  annotate("text", x=2.15, y=1.02, label=paste("p=",round(pvals2["Gammaproteobacteria"],3),sep=""), size=5) +
  
  geom_line(data=data.frame(a=c(2.7,2.7,3,3),b=c(0.56,0.59,0.59,0.56)),aes(x=a,y=b)) + #0.1
  annotate("text", x=2.85, y=0.62, label=paste("p=",round(pvals["Bacteroidia"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(3,3,3.3,3.3),b=c(0.56,0.59,0.59,0.56)),aes(x=a,y=b)) + #0.1
  annotate("text", x=3.15, y=0.62, label="p<0.001", size=5) +
  
  geom_line(data=data.frame(a=c(3.7,3.7,4,4),b=c(0.56,0.59,0.59,0.56)),aes(x=a,y=b)) + #0.1
  annotate("text", x=3.85, y=0.62, label=paste("p=",round(pvals["Clostridia"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(4,4,4.3,4.3),b=c(0.56,0.59,0.59,0.56)),aes(x=a,y=b)) + #0.1
  annotate("text", x=4.15, y=0.62, label="p<0.001", size=5) +

  #geom_line(data=data.frame(a=c(1.75,1.75,2.25,2.25),b=c(0.15,0.19,0.19,0.15)),aes(x=a,y=b)) +
  #annotate("text", x=2, y=0.24, label=paste("p=",round(pvals["Proteobacteria"],3),sep=""), size=5,angle=270) +
  #geom_line(data=data.frame(a=c(2.75,2.75,3.25,3.25),b=c(0.9,0.94,0.94,0.9)),aes(x=a,y=b)) +
  #annotate("text", x=3, y=0.99, label=paste("p=",pvals["Firmicutes"],"*",sep=""), size=5,angle=270) +
  #geom_line(data=data.frame(a=c(3.75,3.75,4.25,4.25),b=c(0.9,0.94,0.94,0.9)),aes(x=a,y=b)) +
  #annotate("text", x=4, y=0.99, label=paste("p=",pvals["Bacteroidetes"],"*",sep=""), size=5,angle=270) +
  xlab("") +
  ylab("Relative abundance") +
  scale_color_manual(values = c("firebrick2","darkorchid1","royalblue2")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="none",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=3),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))#7x9.5/14x8
sp

pdat = data.frame()
pvals = vector()
pvals2 = vector()
#selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
#selsel=c("EX_ala_L(e)","EX_arg_L(e)","EX_asp_L(e)","EX_cys_L(e)","EX_glu_L(e)","EX_gly(e)","EX_his_L(e)","EX_ile_L(e)","EX_leu_L(e)","EX_lys_L(e)","EX_met_L(e)","EX_phe_L(e)","EX_pro_L(e)","EX_ser_L(e)","EX_thr_L(e)","EX_trp_L(e)","EX_tyr_L(e)","EX_val_L(e)")
selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
for(i in rev(selsel)){
  #conc = meanmat_sub[,i]-initsub[i,]#/max(meanmat_sub[,i])
  conc = meanmat_sub_all[,i]
  metnam = gsub("exchange reaction for ","",gsub(" exchange","",subsys[i,"Reaction.Name"]))
  if(metnam == "Butyrate (n-C4:0)"){metnam = "Butyrate"}
  if(metnam == "H2"){metnam = "Hydrogen"}
  pdat = rbind(pdat,data.frame(conc=conc,patient=names(conc),type=as.character(pclass[names(conc),1]),met=rep(metnam,length(conc))))
  pval = wilcox.test(conc[which(pclass[names(conc),1]=="Treatment")],conc[which(pclass[names(conc),1]=="Crohn's disease")],exact=F)$p.value
  pvals[i] = pval   
  pval = wilcox.test(conc[which(pclass[names(conc),1]=="Treatment")],conc[which(pclass[names(conc),1]=="Healthy")],exact=F)$p.value
  pvals2[i] = pval   
}
#pdat$conc = pdat$conc/(10^12* 0.01 * 6.25e-08)
#pdat$conc = pdat$conc/max(pdat$conc)
pdat$type = factor(as.character(pdat$type),levels=c("Crohn's disease","Treatment","Healthy"))
met <- ggplot(pdat, aes(factor(met), conc)) +
  geom_boxplot(aes(color=type),lwd=2, position=position_dodge(0.9),outlier.color=NA) +
  scale_y_continuous(limits = c(0, 200)) +
  #coord_flip() +
  geom_line(data=data.frame(a=c(0.7,0.7,1,1),b=c(145,150,150,145)),aes(x=a,y=b)) + #0.1
  annotate("text", x=0.85, y=155, label=paste("p=",round(pvals["EX_but(e)"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(1,1,1.3,1.3),b=c(145,150,150,145)),aes(x=a,y=b)) + #0.1
  annotate("text", x=1.15, y=155, label=paste("p=",round(pvals2["EX_but(e)"],3),sep=""), size=5) +
  
  geom_line(data=data.frame(a=c(1.7,1.7,2,2),b=c(145,150,150,145)),aes(x=a,y=b)) + #0.1
  annotate("text", x=1.85, y=155, label=paste("p=",round(pvals["EX_ppa(e)"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(2,2,2.3,2.3),b=c(145,150,150,145)),aes(x=a,y=b)) + #0.1
  annotate("text", x=2.15, y=155, label="p<0.001", size=5) +
  
  geom_line(data=data.frame(a=c(2.7,2.7,3,3),b=c(45,50,50,45)),aes(x=a,y=b)) + #0.1
  annotate("text", x=2.85, y=55, label=paste("p=",round(pvals["EX_isobut(e)"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(3,3,3.3,3.3),b=c(45,50,50,45)),aes(x=a,y=b)) + #0.1
  annotate("text", x=3.15, y=55, label="p<0.001", size=5) +
  
  geom_line(data=data.frame(a=c(3.7,3.7,4,4),b=c(45,50,50,45)),aes(x=a,y=b)) + #0.1
  annotate("text", x=3.85, y=55, label=paste("p=",round(pvals["EX_lac_L(e)"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(4,4,4.3,4.3),b=c(45,50,50,45)),aes(x=a,y=b)) + #0.1
  annotate("text", x=4.15, y=55, label=paste("p=",round(pvals2["EX_lac_L(e)"],3),sep=""), size=5) +
  
  geom_line(data=data.frame(a=c(4.7,4.7,5,5),b=c(145,150,150,145)),aes(x=a,y=b)) + #0.1
  annotate("text", x=4.85, y=155, label=paste("p=",round(pvals["EX_ac(e)"],3),sep=""), size=5) +
  geom_line(data=data.frame(a=c(5,5,5.3,5.3),b=c(145,150,150,145)),aes(x=a,y=b)) + #0.1
  annotate("text", x=5.15, y=155, label=paste("p=",round(pvals2["EX_ac(e)"],3),sep=""), size=5) +
  # annotate("text", x=5, y=1.05, label=paste("p=",round(pvals["EX_ac(e)"],4),"*",sep=""), size=5,angle=270) +
  # geom_line(data=data.frame(a=c(4.75,4.75,5.25,5.25),b=c(0.95,1,1,0.95)),aes(x=a,y=b)) +
  # annotate("text", x=4, y=0.9, label=paste("p=",round(pvals["EX_for(e)"],4),"*",sep=""), size=5,angle=270) +
  # geom_line(data=data.frame(a=c(3.75,3.75,4.25,4.25),b=c(0.8,0.85,0.85,0.8)),aes(x=a,y=b)) +
  # annotate("text", x=3, y=0.6, label=paste("p=",round(pvals["EX_but(e)"],4),"*",sep=""), size=5,angle=270) +
  # geom_line(data=data.frame(a=c(2.75,2.75,3.25,3.25),b=c(0.5,0.55,0.55,0.5)),aes(x=a,y=b)) +
  # annotate("text", x=2, y=0.6, label=paste("p=",round(pvals["EX_h2(e)"],4),"*",sep=""), size=5,angle=270) +
  # geom_line(data=data.frame(a=c(1.75,1.75,2.25,2.25),b=c(0.5,0.55,0.55,0.5)),aes(x=a,y=b)) +
  # annotate("text", x=1, y=0.6, label=paste("p=",round(pvals["EX_but(e)"],4),"*",sep=""), size=5,angle=270) +
  # geom_line(data=data.frame(a=c(0.75,0.75,1.25,1.25),b=c(0.5,0.55,0.55,0.5)),aes(x=a,y=b)) +
xlab("") +
  ylab("Concentration in mM") +
  scale_color_manual(values = c("firebrick2","darkorchid1","royalblue2")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.justification = c(1, 0), 
        legend.position = c(0.9, 0.8),
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))#8x10/14x8
met

pan2 = grid.arrange(sp,met,ncol=2)#12x12

pan3 = grid.arrange(pan1,pan2,ncol=1)#25x15

grid.arrange(pan3,ncol=1)#30x22

