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

setwd('P:/GitRep/BacArena')
source(file="R/Arena.R")
source(file="R/Stuff.R")
source(file="R/Substance.R")
source(file="R/Organism.R")
Rcpp::sourceCpp("src/diff.cpp")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI")

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

setwd("P:/MAPPING/Pediatric_crohns/Models/dysbiotic/EU/newest")
load("sublist_eu.RData")
load("poplist_eu.RData")
load("gralist_reduce_eu.RData")
load("gralist_init_eu.RData")
load("fluxlist_red_eu.RData")
load("P:/MAPPING/Models773/agora773_constrain.RData")
names(poplist) = gsub(".mix","",names(poplist))
names(sublist) = gsub(".mix","",names(sublist))

##########################################################################
############### Principle coordinate analysis of Population and Susbtances
##########################################################################

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
dires <- capscale(meanmat_pop~1, distance="bray")
dires_sum <- summary(dires)
gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                   "Phenotype"=as.character(pclass[rownames(dires_sum$sites),1]),
                   "var"=(apply(sdmat,1,mean)/sum(apply(sdmat,1,mean)))[rownames(dires_sum$sites)])
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC1"])
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC2"])
gdat[which(gdat$Phenotype=="Disease control"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Disease control"),"PC1"])
gdat[which(gdat$Phenotype=="Disease control"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Disease control"),"PC2"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
gdat[which(gdat$Phenotype=="Random"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Random"),"PC1"])
gdat[which(gdat$Phenotype=="Random"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Random"),"PC2"])
ppcoa <- ggplot(data=gdat, aes(x=PC1, y=PC2, colour=Phenotype)) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=1) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=1) +   
  stat_ellipse(level = 0.95,size=1) +
  scale_colour_manual(values = c("firebrick2","royalblue2")) +
  geom_point(aes(x=PC1, y=PC2), size=3) +
  geom_point(aes(x=mean.x, y=mean.y), shape=4, size=2, stroke=2) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=0.8) +
  xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
  ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
  ggtitle("Population") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))

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
dires <- capscale(meanmat_sub[,-which(apply(meanmat_sub,2,sum)==0)]~1, distance="bray")
dires_sum <- summary(dires)
gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                   "Phenotype"=as.character(pclass[rownames(dires_sum$sites),1]),
                   "var"=(apply(sdmat,1,mean)/sum(apply(sdmat,1,mean)))[rownames(dires_sum$sites)])
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC1"])
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC2"])
gdat[which(gdat$Phenotype=="Disease control"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Disease control"),"PC1"])
gdat[which(gdat$Phenotype=="Disease control"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Disease control"),"PC2"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
gdat[which(gdat$Phenotype=="Random"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Random"),"PC1"])
gdat[which(gdat$Phenotype=="Random"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Random"),"PC2"])
spcoa <- ggplot(data=gdat, aes(x=PC1, y=PC2, colour=Phenotype)) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=1) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=1) +   
  stat_ellipse(level = 0.85) +
  scale_colour_manual(values = c("firebrick2","darkorchid1","royalblue2","darkgrey")) +
  geom_point(aes(x=PC1, y=PC2), size=4) +
  geom_point(aes(x=mean.x, y=mean.y), shape=4, size=2, stroke=2) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=0.7) +
  xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
  ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
  ggtitle("Metabolites") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))

grid.arrange(ppcoa,spcoa,ncol=1)#9x17


##########################################################################
############### #dimensional plot with SVM
##########################################################################

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
dires <- capscale(meanmat_pop~1, distance="bray")
dires_sum <- summary(dires)

t <- data.frame(x=dires_sum$sites[,1], y=dires_sum$sites[,2], z=dires_sum$sites[,3], cl=NA)
t$cl <- pclass[rownames(dires_sum$sites),1]
#t$cl[which(t$cl=="Crohn's disease")] = -1
#t$cl[which(t$cl=="Healthy")] = 1
t$cl = as.factor(t$cl)
#t = t[-which(rownames(t) %in% c("V1.CD.2","V1.UC.7","V1.CD.11")),]

cweigth = c(table(t$cl)[2]/table(t$cl)[1],1)
names(cweigth)=levels(t$cl)

svm_model <- svm(cl~x+y+z, t, type='C-classification', kernel='linear', class.weights=cweigth)
w <- t(svm_model$coefs) %*% svm_model$SV

detalization <- 100  
data <- list(
  list(
    z = (svm_model$rho- w[1,1]*grid[,1] - w[1,2]*grid[,2]) / w[1,3], 
    x = seq(from=min(t$x),to=max(t$x),length.out=detalization), 
    y = seq(from=min(t$y),to=max(t$y),length.out=detalization), 
    type = "surface"
  )
)
plotly(data)


scatter3D(t$x,t$y,t$z,type="n",bty="u",colvar=NA,shade=3,pch=16,cex=2,phi=20,theta=-120,
          col=c('firebrick2','royalblue2')[as.numeric(t$cl)],
          xlim=c(min(grid[,1]),max(grid[,1])),ylim=c(min(grid[,2]),max(grid[,2])),zlim=c(min(z),max(z)),
          xlab=paste("PCo1",round(dires_sum$cont$importance[2,1]*100, 2),"%","explained variance"),
          ylab=paste("PCo2",round(dires_sum$cont$importance[2,2]*100, 2),"%","explained variance"), 
          zlab=paste("PCo3",round(dires_sum$cont$importance[2,3]*100, 2),"%","explained variance"),
          surf=list(x=matrix(grid[,1],50,50)[1:50,],y=matrix(grid[,2],50,50)[1:50,],z=matrix(z,50,50),
                    shade=0.2, col="black", alpha=0.5, facets = T))


##########################################################################
################################# Investigating the alpha diversity
##########################################################################

poplist=poplist[rownames(pclass)[which(pclass[,1] %in% c("Crohn's disease","Healthy","Random"))]]
meanmat = matrix(0,nrow=length(poplist),ncol=ncol(poplist[[1]][[1]]))
sdmat = matrix(0,nrow=length(poplist),ncol=ncol(poplist[[1]][[1]]))
rownames(meanmat) = names(poplist)
rownames(sdmat) = names(poplist)
for(i in 1:length(poplist)){
  repmat = matrix(0,nrow=length(poplist[[i]]),ncol=ncol(poplist[[i]][[1]]))
  for(j in 1:length(poplist[[i]])){
    repmat[j,] = apply(poplist[[i]][[j]],2,vegan::diversity)
  }
  meanmat[i,] = apply(repmat,2,mean)
  sdmat[i,] = apply(repmat,2,sd)
}
divdat = cbind(melt(meanmat),melt(sdmat)[,3])
colnames(divdat) = c("patient","time","mean","sd")
divdat$type = pclass[as.character(divdat$patient),1]
gpdiv <- ggplot(divdat, aes(x=time, y=mean)) +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,group=patient,color=type,fill=type),alpha=0.5) +
  geom_line(data=data.frame(mean=apply(meanmat[rownames(pclass)[which(pclass[,1]=="Healthy")],],2,median),
                            time=1:ncol(meanmat)),aes(x=time,y=mean),lwd=1.2) +
  scale_fill_manual(values = c("firebrick2","royalblue2","gray50")) +
  scale_color_manual(values = c("firebrick2","royalblue2","gray50")) +
  ggtitle("Population") +
  xlab("Time in h") +
  ylab("Diversity (Shannon index)") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.key=element_blank(),
        legend.position="none",
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=30))
psdat = data.frame(patient=rownames(meanmat),pop=meanmat[,ncol(meanmat)],type=pclass[rownames(meanmat),1])
sublist=sublist[names(poplist)]
meanmat = matrix(0,nrow=length(sublist),ncol=ncol(sublist[[1]][[1]]))
sdmat = matrix(0,nrow=length(sublist),ncol=ncol(sublist[[1]][[1]]))
rownames(meanmat) = names(sublist)
rownames(sdmat) = names(sublist)
for(i in 1:length(sublist)){
  repmat = matrix(0,nrow=length(sublist[[i]]),ncol=ncol(sublist[[i]][[1]]))
  for(j in 1:length(sublist[[i]])){
    repmat[j,] = apply(sublist[[i]][[j]],2,vegan::diversity)
  }
  meanmat[i,] = apply(repmat,2,mean)
  sdmat[i,] = apply(repmat,2,sd)
}
divdat = cbind(melt(meanmat),melt(sdmat)[,3])
colnames(divdat) = c("patient","time","mean","sd")
divdat$type = pclass[as.character(divdat$patient),1]
gsdiv <- ggplot(divdat, aes(x=time, y=mean)) +
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd,group=patient,color=type,fill=type),alpha=0.5) +
  geom_line(data=data.frame(mean=apply(meanmat[rownames(pclass)[which(pclass[,1]=="healthy")],],2,median),
                            time=1:ncol(meanmat)),aes(x=time,y=mean),lwd=1.2) +
  scale_fill_manual(values = c("firebrick2","royalblue2","gray50")) +
  scale_color_manual(values = c("firebrick2","royalblue2","gray50")) +
  ggtitle("Metabolites") +
  xlab("Time in h") +
  ylab("Diversity (Shannon index)") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.key=element_blank(),
        legend.position="none",
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=30))

psdat$sub = meanmat[,ncol(meanmat)]
gspdiv <- ggplot(psdat, aes(x=pop, y=sub)) +
  geom_point(aes(x=pop,y=sub,group=patient,color=type,fill=type),size=4) +
  scale_fill_manual(values = c("firebrick2","royalblue2","gray50")) +
  scale_color_manual(values = c("firebrick2","royalblue2","gray50")) +
  xlab("Population diversity") +
  ylab("Metabolite diversity") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=30))#10x6

modab_pa = t(ifelse(meanmat_pop==0,0,1))
rownames(modab_pa)[which(rownames(modab_pa)=="Ruminococcus_albus_7")] = "Ruminococcus_albus_7_DSM_20455"
rownames(modab_pa)[which(rownames(modab_pa)=="Bifidobacterium_thermacidophilum_subsp_thermacidophilum_DSM_15837")] = "Bifidobacterium_thermacidophilum_subsp_thermacidophilum_DSM_158"
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
psdat$rxn = apply(patrxn,1,sum)[as.character(psdat$patient)]
gspdivrm <- ggplot(psdat, aes(x=rxn, y=sub)) +
  geom_point(aes(x=rxn,y=sub,group=patient,color=type,fill=type),size=4) +
  scale_fill_manual(values = c("firebrick2","royalblue2","gray50")) +
  scale_color_manual(values = c("firebrick2","royalblue2","gray50")) +
  xlab("Number of reactions") +
  ylab("Metabolite diversity") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.key=element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=20),
        axis.title.y = element_text(size=30,vjust=0.5),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        plot.title = element_text(size=30))#10x6

grid.arrange(gpdiv,gspdiv,gsdiv,gspdivrm,ncol=2)#,width=c(5,6,5,6))#19x12

#pclass = data.frame(pclass[rownames(pclass)[which(pclass[,1] %in% c("Crohn's disease","Healthy"))],])
##########################################################################
################################# What species/metabolites are depleted
##########################################################################

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
meanmat_pop = meanmat_pop[intersect(rownames(pclass)[which(pclass[,1] %in% c("Crohn's disease","Healthy"))],rownames(meanmat_pop)),]
meanmat_sub = meanmat_sub[intersect(rownames(pclass)[which(pclass[,1] %in% c("Crohn's disease","Healthy"))],rownames(meanmat_sub)),]

#meanmat_pop = rbind(meanmat_pop,t(popmod))
#meanmat_sub = rbind(meanmat_sub,t(submod))

relab = meanmat_pop/apply(meanmat_pop,1,sum)
melt(apply(relab[,rownames(modstats)[which(modstats$Phylum=="Firmicutes")]],1,sum))
melt(apply(relab[,rownames(modstats)[which(modstats$Class=="Clostridia")]],1,sum))
melt(apply(relab[,rownames(modstats)[which(modstats$Family=="Ruminococcaceae")]],1,sum))

melt(apply(relab[,rownames(modstats)[which(modstats$Genus=="Bifidobacterium")]],1,sum))
melt(apply(relab[,rownames(modstats)[which(modstats$Genus=="Lactobacillus")]],1,sum))
melt(apply(relab[,rownames(modstats)[which(modstats$Species=="Faecalibacterium prausnitzii")]],1,sum))
pdat = data.frame()
pvals = vector()
for(i in rev(c("Bacteroidetes","Firmicutes","Proteobacteria","Actinobacteria"))){
  ab = relab[,intersect(rownames(modstats)[which(modstats$Phylum==i)],colnames(relab))]
  if(sum(ab)!=0){
    if(class(ab)!="numeric"){
      ab=apply(ab,1,sum)
    }
    pdat = rbind(pdat,data.frame(abundance=ab,patient=rownames(relab),
                      type=as.character(pclass[rownames(relab),1]),phylum=rep(i,nrow(relab))))
    pval = t.test(ab[which(pclass[names(ab),1]=="Healthy")],ab[which(pclass[names(ab),1]=="Crohn's disease")])$p.value
    pvals[i] = round(pval,5)
  }
}
sp <- ggplot(pdat, aes(factor(phylum), abundance)) +
  geom_boxplot(aes(color=type),lwd=1.1,position=position_dodge(0.9)) +
  #coord_flip() +
  #geom_line(data=data.frame(a=c(0.75,0.75,1.25,1.25),b=c(0.15,0.19,0.19,0.15)),aes(x=a,y=b)) +
  #annotate("text", x=1, y=0.24, label=paste("p=",round(pvals["Actinobacteria"],3),"*",sep=""), size=5,angle=270) +
  #geom_line(data=data.frame(a=c(1.75,1.75,2.25,2.25),b=c(0.15,0.19,0.19,0.15)),aes(x=a,y=b)) +
  #annotate("text", x=2, y=0.24, label=paste("p=",round(pvals["Proteobacteria"],3),sep=""), size=5,angle=270) +
  #geom_line(data=data.frame(a=c(2.75,2.75,3.25,3.25),b=c(0.9,0.94,0.94,0.9)),aes(x=a,y=b)) +
  #annotate("text", x=3, y=0.99, label=paste("p=",pvals["Firmicutes"],"*",sep=""), size=5,angle=270) +
  #geom_line(data=data.frame(a=c(3.75,3.75,4.25,4.25),b=c(0.9,0.94,0.94,0.9)),aes(x=a,y=b)) +
  #annotate("text", x=4, y=0.99, label=paste("p=",pvals["Bacteroidetes"],"*",sep=""), size=5,angle=270) +
  xlab("") +
  ylab("Relative abundance") +
  scale_color_manual(values = c("firebrick2","royalblue2","darkorchid1")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))#7x9.5/14x8

pdat = data.frame()
pvals = vector()
#c("Faecalibacterium","Roseburia","Eubacterium")
for(i in c("Faecalibacterium","Roseburia","Eubacterium")){
  ab = relab[,intersect(rownames(modstats)[which(modstats$Genus==i)],colnames(relab))]
  if(sum(ab)!=0){
    if(class(ab)!="numeric"){
      ab=apply(ab,1,sum)
    }
    pdat = rbind(pdat,data.frame(abundance=ab,patient=rownames(relab),
                                 type=as.character(pclass[rownames(relab),1]),phylum=rep(i,nrow(relab))))
    pval = t.test(ab[which(pclass[names(ab),1]=="Healthy")],ab[which(pclass[names(ab),1]=="Crohn's disease")])$p.value
    pvals[i] = round(pval,5)
  }
}

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
    pval = t.test(ab[which(pclass[names(ab),1]=="Healthy")],ab[which(pclass[names(ab),1]=="Crohn's disease")])$p.value
    pvals[i] = round(pval,5)
  }
}



pdat = data.frame()
pvals = vector()
cors = vector()
cors2 = vector()
lookfor = "PUCAI" 
for(i in levels(modstats$ModelAGORA)){
  ab = relab[,intersect(rownames(modstats)[which(modstats$ModelAGORA==i)],colnames(relab))]
  if(sum(ab)!=0){
    if(class(ab)!="numeric"){ab=apply(ab,1,sum)}
    pdat = rbind(pdat,data.frame(abundance=ab,patient=rownames(relab),
                                 type=as.character(pclass[rownames(relab),1]),phylum=rep(i,nrow(relab))))
    pval = t.test(ab[which(pclass[names(ab),1]=="Healthy")],ab[which(pclass[names(ab),1]=="Crohn's disease")])$p.value
    pvals[i] = round(pval,5)
    cors[i] = rcorr(cbind(ab,meta[names(ab),lookfor]),type="spearman")$P[1,2]
    cors2[i] = cor(cbind(ab,meta[names(ab),lookfor]),method="spearman",use="complete.obs")[1,2]
  }
}
ab = relab[,names(which.min(na.omit(cors)))]
if(class(ab)!="numeric"){ab=apply(ab,1,sum)}
plot(cbind(ab,meta[names(ab),lookfor]),main=names(which.min(na.omit(cors))))
#Distance


modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Genus, modstats$Species),]
mdat = t(relab[,rownames(modtax)])
mdat = mdat[,rownames(pclass[colnames(mdat),])[order(pclass[colnames(mdat),1])]]
rowannot = data.frame(modtax[,c('Phylum','Class')])
colannot = data.frame(Phenotype=as.character(pclass[colnames(mdat),1]))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
ann_colors = list(Class=getPalette(length(levels(modtax$Class))),
                  Phenotype=c("firebrick2","royalblue2"))
mdat_red = mdat[-which(apply(mdat,1,sum)==0),]
#mdat_red = mdat[rownames(modstats)[which(modstats$Genus %in% names(which(pvals<0.0005)))],]#names(which(pvals<0.05))
aheatmap(t(mdat_red),scale="column", Colv=NA, Rowv=NA, labRo=NA, labCol=NA,
         annCol=rowannot[rownames(mdat_red),], annColors=ann_colors, annRow=colannot,
         color=colorRampPalette(c("seashell","darkorchid2","darkorchid4"))(100))#30x9


modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Genus, modstats$Species),]
mdat = t(relab[,rownames(modtax)])
mdat = mdat[,rownames(pclass[colnames(mdat),])[order(pclass[colnames(mdat),1])]]
rowannot = data.frame(modtax[,c('Phylum')])
rownames(rowannot) = rownames(modtax); colnames(rowannot) = "Phylum"
colannot = data.frame(Phenotype=as.character(pclass[colnames(mdat),1]))
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
ann_colors = list(Phylum=getPalette(length(levels(modtax$Phylum))),
                  Phenotype=c("firebrick2","royalblue2"))
mdat_red = mdat[-which(apply(mdat,1,sum)==0),]
#mdat_red = mdat[rownames(modstats)[which(modstats$Genus %in% names(which(pvals<0.0005)))],]#names(which(pvals<0.05))
aheatmap(t(mdat_red),scale="column", Colv=NA, Rowv=NA, labRo=NA, labCol=NA,
         annCol=data.frame(Phylum=rowannot[rownames(mdat_red),]), annColors=ann_colors, annRow=colannot,
         color=colorRampPalette(c("seashell","lightsalmon3","lightsalmon4"))(100))#12x3

dires <- capscale(meanmat_pop~1, distance="bray")
dires_sum <- summary(dires)
imp=20
abseigen <- abs(dires_sum$species)
absort <- sort(abs(abseigen[,1]), decreasing=T)
sel <- names(head(absort, imp))
beigen = data.frame(eigen=absort/sum(absort),bac=factor(names(absort),levels=names(absort)))
speigen = ggplot(beigen,aes(x=bac,y=eigen)) +
  geom_bar(stat="identity") +
  geom_vline(xintercept=20, col="red", lwd=1.1) +
  xlab("Microbes") +
  ylab("Relative contribution") +
  ggtitle("Population") +
  theme_bw(base_size = 30) +
  theme(legend.position="none",
        panel.border = element_rect(colour='black',size=2),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30))
selnam = unlist(lapply(strsplit(sel,"_"),function(x){paste(substr(x[1],1,1),".",x[2]," ",paste(x[3:length(x)],collapse="_"),sep="")}))
selnam = gsub("F.cf prausnitzii_KLE1255","F.prausnitzii KLE1255",selnam)
eigsel = data.frame(eigen=dires_sum$species[sel,1],bac=factor(selnam,levels=selnam),order=as.character(modstats[sel,"Order"]))
speigen_sel = ggplot(eigsel,aes(x=bac,y=eigen)) +
  geom_bar(stat="identity",aes(fill=order)) +
  scale_y_continuous(limits = c(-max(abs(eigsel$eigen)), max(abs(eigsel$eigen)))) +
  #geom_text(aes(y=0,hjust=ifelse(eigen>0,1.02,-0.02),label=bac),colour="black",size=4.2) + 
  coord_flip() +
  ylab("Relative contribution") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_blank(),
        plot.title = element_text(size=30))
gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                   "Phenotype"=as.character(pclass[rownames(dires_sum$sites),1]),
                   "var"=(apply(sdmat,1,mean)/sum(apply(sdmat,1,mean)))[rownames(dires_sum$sites)])
gdat2 <- data.frame('R1'=dires_sum$species[sel,1],
                    'R2'=dires_sum$species[sel,2],
                    'labs'=sel,"order"=as.character(modstats[sel,"Order"]))
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC1"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC2"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
ppcoa <- ggplot(data=gdat, aes(x=PC1, y=PC2, colour=Phenotype)) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=1) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=1) +   
  #geom_segment(data=gdat2, aes(x=0, y=0, xend=R1, yend=R2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.5, color="black", size=1) +
  #geom_text(data=gdat2, aes(x=R1, y=R2, label=labs), color="black", size=3) +
  stat_ellipse(level = 0.85) +
  scale_colour_manual(values = c("firebrick2","royalblue2")) +
  geom_point(aes(x=PC1, y=PC2), size=5) +
  geom_point(aes(x=mean.x, y=mean.y), shape=4, size=4, stroke=2) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=1) +
  xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
  ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
  ggtitle("Population") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))

dires <- capscale(meanmat_sub[,-which(apply(meanmat_sub,2,sum)==0)]~1, distance="bray")
dires_sum <- summary(dires)
imp=20
abseigen <- abs(dires_sum$species)
absort <- sort(abs(abseigen[,1]), decreasing=T)
sel <- names(head(absort, imp))
beigen = data.frame(eigen=absort/sum(absort),bac=factor(names(absort),levels=names(absort)))
smeigen = ggplot(beigen,aes(x=bac,y=eigen)) +
  geom_bar(stat="identity") +
  geom_vline(xintercept=20, col="red", lwd=1.1) +
  xlab("Substances") +
  ylab("Relative contribution") +
  #ggtitle("Metabolites") +
  theme_bw(base_size = 30) +
  theme(legend.position="none",
        panel.border = element_rect(colour='black',size=2),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30))
eigsel = data.frame(eigen=dires_sum$species[sel,1],bac=factor(gsub("exchange reaction for ","",gsub(" exchange","",subsys[sel,"Reaction.Name"])),
                                                              levels=gsub("exchange reaction for ","",gsub(" exchange","",subsys[sel,"Reaction.Name"]))))
smeigen_sel = ggplot(eigsel,aes(x=bac,y=eigen)) +
  geom_bar(stat="identity",aes(fill="black")) +
  scale_y_continuous(limits = c(-max(abs(eigsel$eigen)), max(abs(eigsel$eigen)))) +
  #geom_text(aes(y=0,hjust=ifelse(eigen>0,1.05,-0.05),label=bac),colour="black",size=6) + 
  #geom_text(aes(y=0,hjust=1,label=bac),colour="black",size=6) + 
  scale_fill_manual(values=c("grey30")) +
  coord_flip() +
  ylab("Relative contribution") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=17),
        axis.title.y = element_blank(),
        plot.title = element_text(size=30))

gdat <- data.frame("PC1"=dires_sum$sites[,1], "PC2"=dires_sum$sites[,2], 
                   "Phenotype"=as.character(pclass[rownames(dires_sum$sites),1]),
                   "var"=(apply(sdmat,1,mean)/sum(apply(sdmat,1,mean)))[rownames(dires_sum$sites)])
gdat2 <- data.frame('R1'=dires_sum$species[sel,1],
                    'R2'=dires_sum$species[sel,2],
                    'labs'=sel)
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC1"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.x"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC1"])
gdat[which(gdat$Phenotype=="Crohn's disease"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Crohn's disease"),"PC2"])
gdat[which(gdat$Phenotype=="Healthy"),"mean.y"] = mean(gdat[which(gdat$Phenotype=="Healthy"),"PC2"])
spcoa <- ggplot(data=gdat, aes(x=PC1, y=PC2, colour=Phenotype)) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=1) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=1) +   
  #geom_segment(data=gdat2, aes(x=0, y=0, xend=R1, yend=R2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.5, color="black", size=1) +
  stat_ellipse(level = 0.85) +
  scale_colour_manual(values = c("firebrick2","royalblue2")) +
  geom_point(aes(x=PC1, y=PC2), size=5) +
  geom_point(aes(x=mean.x, y=mean.y), shape=4, size=4, stroke=2) +
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2), lwd=1) +
  xlab(paste("PCo1 ", round(dires_sum$cont$importance[2,1]*100, 2), "%", " explained variance", sep="")) +
  ylab(paste("PCo2 ", round(dires_sum$cont$importance[2,2]*100, 2), "%", " explained variance", sep="")) +
  ggtitle("Metabolites") +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))

grid.arrange(ppcoa,speigen_sel,spcoa,smeigen_sel,ncol=2,widths=c(8,6.5))#14x18

grid.arrange(ppcoa,sp,spcoa,met,ncol=2)#,width=c(5,6,5,6))#15x15
##########################################################################
################################# Analyses of metabolite selections
##########################################################################

initsub = matrix(0, nrow=nrow(alldiet),ncol=length(sublist),
                 dimnames = list(rownames(alldiet),names(sublist)))
for(i in 1:length(sublist)){
  initsub[rownames(sublist[[i]][[1]]),i] = sublist[[i]][[1]][,1]
}
initsub = initsub[,rownames(meanmat_sub)]
intersect(rownames(alldiet)[which(alldiet$EU!=0)],sel)
pdat = data.frame()
pvals = vector()
#selsel = c("EX_for(e)","EX_ac(e)","EX_lac_L(e)","EX_but(e)","EX_acald(e)","EX_succ(e)","EX_arg_L(e)","EX_isobut(e)","EX_h2(e)","EX_acac(e)",
#           "EX_acgam(e)","EX_ch4(e)")
#selsel=c("EX_acgam(e)","EX_T_antigen(e)","EX_hspg_rest(e)")
#selsel=c("EX_lac_D(e)","EX_succ(e)","EX_for(e)","EX_ac(e)","EX_lac_L(e)","EX_but(e)","EX_ppa(e)","EX_acald(e)","EX_isobut(e)","EX_h2(e)","EX_acac(e)","EX_fum(e)")

selsel=c("EX_for(e)","EX_ac(e)","EX_lac_L(e)","EX_but(e)","EX_acald(e)","EX_isobut(e)","EX_h2(e)","EX_ppa(e)","EX_fum(e)","EX_succ(e)")
#selsel=c("EX_ac(e)","EX_for(e)","EX_but(e)","EX_h2(e)","EX_isobut(e)")
#names(head(absort, length(absort)))
selsel=c("EX_ac(e)","EX_lac_L(e)","EX_isobut(e)","EX_ppa(e)","EX_but(e)")
#selsel=c("EX_ac(e)","EX_lac_L(e)","EX_h2(e)","EX_ppa(e)","EX_but(e)","EX_isobut(e)","EX_isoval(e)")


for(i in rev(selsel)){
  #conc = meanmat_sub[,i]-initsub[i,]#/max(meanmat_sub[,i])
  conc = meanmat_sub[,i]
  metnam = gsub("exchange reaction for ","",gsub(" exchange","",subsys[i,"Reaction.Name"]))
  if(metnam == "Butyrate (n-C4:0)"){metnam = "Butyrate"}
  if(metnam == "H2"){metnam = "Hydrogen"}
  pdat = rbind(pdat,data.frame(conc=conc,patient=names(conc),type=as.character(pclass[names(conc),1]),met=rep(metnam,length(conc))))
  pval = t.test(conc[which(pclass[names(conc),1]=="Healthy")],conc[which(pclass[names(conc),1]=="Crohn's disease")])$p.value
  pvals[i] = pval             
}
pdat$conc = pdat$conc/(10^12* 0.01 * 6.25e-08)
#pdat$conc = pdat$conc/max(pdat$conc)
met <- ggplot(pdat, aes(factor(met), conc)) +
  geom_boxplot(aes(color=type),lwd=1.1,position=position_dodge(0.9)) +
  #coord_flip() +
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
  scale_color_manual(values = c("firebrick2","royalblue2","darkorchid1")) +
  theme_bw(base_size = 30) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=30)) +
  guides(colour = guide_legend(override.aes = list(shape = 15)))#8x10/14x8
met

test = names(which(pvals<0.001))
View(cbind(names(head(absort, length(absort))),
      gsub("exchange reaction for ","",gsub(" exchange","",subsys[names(head(absort, length(absort))),"Reaction.Name"]))))

##########################################################################
################################# Analyses of glycan gradients
##########################################################################

initsub = matrix(0, nrow=nrow(alldiet),ncol=length(sublist),
                 dimnames = list(rownames(alldiet),names(sublist)))
for(i in 1:length(sublist)){
  initsub[rownames(sublist[[i]][[1]]),i] = sublist[[i]][[1]][,1]
}
initsub = initsub[,rownames(meanmat_sub)]
pdat = data.frame()
pvals = vector()
ratio = vector()
for(i in rownames(alldiet)[which(alldiet$type=="Mucin")]){
  conc = meanmat_sub[,i]-initsub[i,]
  metnam = gsub("exchange reaction for ","",gsub(" exchange","",subsys[i,"Reaction.Name"]))
  pdat = rbind(pdat,data.frame(conc=conc,patient=names(conc),type=as.character(pclass[names(conc),1]),met=rep(metnam,length(conc))))
  pval = t.test(conc[which(pclass[names(conc),1]=="Healthy")],conc[which(pclass[names(conc),1]=="Crohn's disease")])$p.value
  pvals[i] = pval
  if(sum(conc)<0){
    ratio[i] = mean(meanmat_sub[which(pclass[names(conc),1]=="Crohn's disease"),i])/mean(meanmat_sub[which(pclass[names(conc),1]=="Healthy"),i])
  }
}
length(which(sort(ratio)<1))/length(ratio)
#mucs = intersect(names(which(sort(ratio)<1)),names(which(pvals<0.1)))
# mucs = c("EX_f1a(e)","EX_gncore1(e)","EX_core5(e)","EX_core4(e)","EX_core7(e)")
# par(mfrow = c(5,2))
# for(i in mucs){
#   fuclist = lapply(gralist_reduce[rownames(pclass)[which(pclass[,1]=="Healthy")]],function(x){x[[i]]})
#   image(matrix(apply(do.call(rbind,fuclist),2,mean),nrow=200,ncol=200),axes=FALSE,
#         col=colorRampPalette(c("seashell", "indianred4"))(10))
#   title(main=i)
#   fuclist = lapply(gralist_reduce[rownames(pclass)[which(pclass[,1]=="Crohn's disease")]],function(x){x[[i]]})
#   image(matrix(apply(do.call(rbind,fuclist),2,mean),nrow=200,ncol=200),axes=FALSE,
#         col=colorRampPalette(c("seashell", "indianred4"))(10))
#   title(main=i)
# }

par(mfrow = c(2,2))
fuclist = lapply(gralist_reduce[rownames(pclass)[which(pclass[,1]=="Healthy")]],function(x){x[["EX_f1a(e)"]]})
image(matrix(apply(do.call(rbind,fuclist),2,mean),nrow=200,ncol=200),axes=FALSE,
      col=colorRampPalette(c("seashell", "darkblue"))(50))
title(main="Healthy individuals")
fuclist = lapply(gralist_reduce[rownames(pclass)[which(pclass[,1]=="Crohn's disease")]],function(x){x[["EX_f1a(e)"]]})
image(matrix(apply(do.call(rbind,fuclist),2,mean),nrow=200,ncol=200),axes=FALSE,
      col=colorRampPalette(c("seashell", "darkblue"))(50))
title(main="Crohn's disease patients")
fuclist = lapply(gralist_reduce[rownames(pclass)[which(pclass[,1]=="Healthy")]],function(x){x[["EX_gncore1(e)"]]})
image(matrix(apply(do.call(rbind,fuclist),2,mean),nrow=200,ncol=200),axes=FALSE,
      col=colorRampPalette(c("seashell", "darkblue"))(50))
title(main="Healthy individuals")
fuclist = lapply(gralist_reduce[rownames(pclass)[which(pclass[,1]=="Crohn's disease")]],function(x){x[["EX_gncore1(e)"]]})
image(matrix(apply(do.call(rbind,fuclist),2,mean),nrow=200,ncol=200),axes=FALSE,
      col=colorRampPalette(c("seashell", "darkblue"))(50))
title(main="Crohn's disease patients")
# image(matrix(gralist_reduce[["V1.CD.1"]][["EX_f1a(e)"]],nrow=200,ncol=200),axes=FALSE,
#       col=colorRampPalette(c("seashell", "darkblue"))(100))

par(mfrow = c(4,5))
for(i in as.character(pclass[which(pclass$V2 %in% c("Crohn's disease","Healthy")),2])){
  dat <- datlist[[i]][[13]][[1]]
  dat$type = modstats[dat$type,"Phylum"]
  dat$type = as.character(dat$type)
  dat$type[which(!(dat$type %in% c("Firmicutes","Bacteroidetes")))] = "Xother"
  dat$type = as.factor(dat$type)
  plot(dat$x,dat$y,xlim=c(0,200),ylim=c(0,200),col=c(rainbow(nlevels(dat$type)-1),NA)[dat$type],axes=FALSE,cex=0.1,xlab='',ylab='',pch=16)
}

##########################################################################
################################# Analyses of secretion/consumption paterns
##########################################################################

modtax <- modstats[order(modstats$Phylum, modstats$Class, modstats$Genus, modstats$Species),]

getFluxMet(fluxlist_red,"EX_but(e)",meanmat_pop,modtax,subsys,ylab=T)

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
  return(pfdat)
}
#Looking at fermentation pattern
plotFluxMet <- function(fluxlist_red,meti,meanmat_pop,modtax,subsys,ylab=F){
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
  pfmat = matrix(0, nrow=nlevels(modtax[,"Phylum"]), ncol=nrow(fmat),
                 dimnames=list(levels(modtax[,"Phylum"]), rownames(fmat)))
  for(i in rownames(pfmat)){
    if(class(fmat[,which(modtax[colnames(fmat),"Phylum"]==i)])!="numeric"){
      msum = apply(fmat[,which(modtax[colnames(fmat),"Phylum"]==i)],1,sum)
    }else{
      msum = fmat[,which(modtax[colnames(fmat),"Phylum"]==i)]
    }
    pfmat[i,names(msum)] = msum
  }
  pfdat = melt(pfmat)
  colnames(pfdat) = c("Phylum","Patient","Flux")
  pfdat$Type = pclass[as.character(pfdat$Patient),1]
  pfdat$Met = meti
  gp <- ggplot(pfdat, aes(factor(Phylum), Flux)) +
    geom_boxplot(aes(color=Type),lwd=1.1,position=position_dodge(0.9)) +
    #scale_y_continuous(limits = c(min(pfdat$Flux),max(pfdat$Flux)),breaks=seq(min(pfdat$Flux),max(pfdat$Flux),by=round(max(pfdat$Flux)/10000000,0)*10000)) +
    coord_flip() +
    scale_color_manual(values = c("firebrick2","royalblue2")) +
    ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
    xlab("") +
    ylab("Total population flux") +
    theme_bw(base_size = 20) +
    theme(legend.text=element_text(size=20),
          legend.position="top",
          legend.key=element_blank(),
          legend.title=element_blank(),
          panel.border = element_rect(colour='black',size=2),
          axis.text.y = element_text(size=25),
          axis.text.x = element_text(size=15),
          axis.title.x = element_text(size=25),
          axis.title.y = element_blank(),
          axis.ticks = element_line(size=1,color='black'),
          axis.ticks.length=unit(0.3,'cm'),
          plot.title = element_text(size=25,hjust = 0.5))
  if(!ylab){
    gp <- gp + theme(axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
  }
  return(gp)
}

gbu <- getFluxMet(fluxlist_red,"EX_but(e)",meanmat_pop,modtax,subsys,ylab=T)
gpa <- getFluxMet(fluxlist_red,"EX_ppa(e)",meanmat_pop,modtax,subsys,ylab=F)
gh2 <- getFluxMet(fluxlist_red,"EX_h2(e)",meanmat_pop,modtax,subsys,ylab=F)
gla <- getFluxMet(fluxlist_red,"EX_lac_L(e)",meanmat_pop,modtax,subsys,ylab=F)

grid.arrange(gbu,gpa,gh2,gla,ncol=4,widths=c(8,5,5,5))#20x6

popflx = rbind(getFluxMet(fluxlist_red,"EX_but(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
      getFluxMet(fluxlist_red,"EX_ppa(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
      getFluxMet(fluxlist_red,"EX_h2(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
      getFluxMet(fluxlist_red,"EX_lac_L(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")),
      getFluxMet(fluxlist_red,"EX_ac(e)",meanmat_pop,modtax,subsys,ylab=T,tax=c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")))
popflx$Met[which(popflx$Met=="EX_but(e)")] = "Butyrate"
popflx$Met[which(popflx$Met=="EX_ppa(e)")] = "Propionate"
popflx$Met[which(popflx$Met=="EX_h2(e)")] = "Hydrogen"
popflx$Met[which(popflx$Met=="EX_lac_L(e)")] = "L-Lactate"
popflx$Met[which(popflx$Met=="EX_ac(e)")] = "Acetate"
popflx$Met = factor(popflx$Met,levels=c("Butyrate","Propionate","Hydrogen","L-Lactate","Acetate"))
popflx$Class = factor(as.character(popflx$Class),levels=rev(c("Clostridia","Bacteroidia","Gammaproteobacteria","Bacilli")))
ggplot(popflx, aes(x=Class, y=Flux/max(Flux), fill=Patient)) +
  geom_bar(position="dodge", stat="identity") +
  #scale_y_continuous(limits = c(min(pfdat$Flux),max(pfdat$Flux)),breaks=seq(min(pfdat$Flux),max(pfdat$Flux),by=round(max(pfdat$Flux)/10000000,0)*10000)) +
  coord_flip() +
  facet_grid(~Met) +
  scale_fill_manual(values = c("firebrick2","royalblue2")) +
  ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
  xlab("") +
  ylab("Relative population flux") +
  theme_bw(base_size = 20) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_blank())#16x4

##########################################################################
#################### Trace back the most active butyrate producer
##########################################################################

meti="EX_but(e)"
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

fmat_red = fmat[,rownames(modstats)[which(modstats$Phylum=="Firmicutes")]]
fmat_red = fmat_red[,-which(apply(fmat_red,2,sum)==0)]
fmat_red = fmat_red + 1
fchange = vector()
for(i in 1:ncol(fmat_red)){
  fchange[i] = median(fmat_red[which(pclass[rownames(fmat_red),1]=="Healthy"),i])/median(fmat_red[which(pclass[rownames(fmat_red),1]=="Crohn's disease"),i])
}
topbut = colnames(fmat_red[,which(fchange>1)])

fdat = melt(fmat[,topbut])
colnames(fdat) = c("Patient","Microbe","Flux")
selnam = unlist(lapply(strsplit(as.character(fdat$Microbe),"_"),function(x){paste(substr(x[1],1,1),".",x[2]," ",paste(x[3:length(x)],collapse="_"),sep="")}))
fdat$Microbe = factor(selnam)
fdat$Type = pclass[as.character(fdat$Patient),1]
gpop <- ggplot(fdat, aes(factor(Microbe), Flux)) +
  geom_boxplot(aes(color=Type),lwd=1.1,position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(min(fdat$Flux),1300000)) +
  coord_flip() +
  scale_color_manual(values = c("firebrick2","royalblue2")) +
  ylab("Total population flux") +
  ggtitle(gsub(" exchange","",subsys[meti,"Reaction.Name"])) +
  theme_bw(base_size = 20) +
  theme(legend.text=element_text(size=20),
        legend.position="top",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_text(size=25,hjust = 0.5))
gpop #7x8

##########################################################################
#################### Reduction to individual models
##########################################################################

butspecs = specs[topbut]
# oldbut = vector()
# for(i in 1:length(butspecs)){
#   exsrxn = butspecs[[i]]@react_id[grep('EX_',butspecs[[i]]@react_id)]
#   spmod = changeObjFunc(butspecs[[i]],"EX_but(e)")
#   spmod = changeBounds(spmod,exsrxn,lb=-(alldiet[exsrxn,"EU"]+0.0001),ub=1000)
#   lpsol = optimizeProb(spmod,retOptSol=F)
#   oldbut[i] = lpsol$obj
# }
fac=1
oldbio = vector()
oldbut = vector()
for(i in 1:length(butspecs)){
  exsrxn = butspecs[[i]]@react_id[grep('EX_',butspecs[[i]]@react_id)]
  spmod = changeBounds(butspecs[[i]],exsrxn,lb=-(alldiet[exsrxn,"EU"]+0.0001),ub=1000)
  lpsol = optimizeProb(spmod,retOptSol=F)
  oldbio[i] = lpsol$obj
  oldbut[i] = lpsol$fluxes[which(spmod@react_id=="EX_but(e)")]
}
oldbut = ifelse(oldbut<=0,1,oldbut)
exlist = vector()
for(i in 1:length(butspecs)){
  exsrxn = butspecs[[i]]@react_id[grep('EX_',butspecs[[i]]@react_id)]
  spmod = changeBounds(butspecs[[i]],exsrxn,lb=-(alldiet[exsrxn,"EU"]+0.0001),ub=1000)
  for(j in exsrxn){
    spmod = changeBounds(spmod,j,lb=-1000,ub=1000)
    lpsol = optimizeProb(spmod,retOptSol=F)
    if(lpsol$obj>oldbio[i]*fac && lpsol$fluxes[which(spmod@react_id=="EX_but(e)")]>oldbut[i]*fac){
      exlist = c(exlist,j)
    }
  }
}
busubs = unique(exlist)

pbac = meanmat_pop[rownames(pclass)[which(pclass[,1]=="Crohn's disease")],rownames(modstats[which(modstats$Phylum=="Bacteroidetes"),])]
pbac = pbac[,-which(apply(pbac,2,sum)==0)]
rbac = rpa[colnames(pbac),]
rbac = rbac[,grep("EX_",colnames(rbac))]
rbac = rbac[,-which(apply(rbac,2,sum)==0)]

#design the diet
alldiet[,"EU_mod"] = alldiet[,"EU"]
alldiet[setdiff(busubs,colnames(rbac)),"EU_mod"] = 10000
write.csv(alldiet,file="P:/MAPPING/Model_data/diets_eu_mod.csv")

##########################################################################
#################### Find out who the most prominent butyrate producer are
##########################################################################

meti="EX_but(e)"
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
fdat = melt(fmat[,names(tail(sort(apply(fmat,2,mean)),10))])
colnames(fdat) = c("Patient","Microbe","Flux")
selnam = unlist(lapply(strsplit(as.character(fdat$Microbe),"_"),function(x){paste(substr(x[1],1,1),".",x[2]," ",paste(x[3:length(x)],collapse="_"),sep="")}))
selnam = gsub("F.cf prausnitzii_KLE1255","F.prausnitzii KLE1255",selnam)
fdat$Microbe = factor(selnam)
fdat$Type = pclass[as.character(fdat$Patient),1]
gpop <- ggplot(fdat, aes(factor(Microbe), Flux)) +
  geom_boxplot(aes(color=Type),lwd=1.1,position=position_dodge(0.9)) +
  scale_y_continuous(limits = c(min(fdat$Flux),1300000)) +
  coord_flip() +
  scale_color_manual(values = c("firebrick2","royalblue2")) +
  ylab("Total population butyrate flux") +
  theme_bw(base_size = 20) +
  theme(legend.text=element_text(size=20),
        legend.position="none",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.text.y = element_text(size=25),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_blank())
fmatrel = fmat/meanmat_pop
fmatrel = ifelse(is.na(fmatrel),0,fmatrel)
fdatrel = melt(fmatrel[,names(tail(sort(apply(fmat,2,mean)),10))])
colnames(fdatrel) = c("Patient","Microbe","Flux")
selnam = unlist(lapply(strsplit(as.character(fdatrel$Microbe),"_"),function(x){paste(substr(x[1],1,1),".",x[2]," ",paste(x[3:length(x)],collapse="_"),sep="")}))
selnam = gsub("F.cf prausnitzii_KLE1255","F.prausnitzii KLE1255",selnam)
fdatrel$Microbe = factor(selnam)
fdatrel$Type = pclass[as.character(fdatrel$Patient),1]
grel <- ggplot(fdatrel, aes(factor(Microbe), Flux)) +
  geom_boxplot(aes(color=Type),lwd=1.1,position=position_dodge(0.9)) +
  coord_flip() +
  scale_color_manual(values = c("firebrick2","royalblue2")) +
  ylab("Mean individual butyrate flux") +
  theme_bw(base_size = 20) +
  theme(legend.text=element_text(size=20),
        legend.position="none",
        legend.key=element_blank(),
        legend.title=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.text.x = element_text(size=20),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        plot.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

msel = names(tail(sort(apply(fmat,2,mean)),10))
grid.arrange(gpop,grel,ncol=2,widths=c(8,4.5))#15x8

##########################################################################
#################### Reduction to individual models
##########################################################################

msel = topbut
exflux = lapply(fluxlist_red[[1]][msel],function(x,exs){x[intersect(names(x),exs)]},exs=colnames(meanmat_sub)[grep("EX_",colnames(meanmat_sub))])
specmlist = lapply(exflux,function(x,pats){
  matrix(0, nrow=length(x), ncol=length(pats),dimnames=list(names(x),pats))
},pats=rownames(meanmat_sub))
for(j in rownames(meanmat_sub)){
  exflux = lapply(fluxlist_red[[j]][msel],function(x,exs){x[intersect(names(x),exs)]},exs=colnames(meanmat_sub)[grep("EX_",colnames(meanmat_sub))])
  for(i in names(exflux)){
    iflux = exflux[[i]]/meanmat_pop[j,i]
    specmlist[[i]][names(iflux[which(iflux<0)]),j] = iflux[which(iflux<0)]
  }
}
specmlistr = specmlist
for(i in names(specmlistr)){rownames(specmlistr[[i]]) = paste(i,rownames(specmlistr[[i]]),sep=".")}
specm = do.call(rbind,specmlistr)
specm = specm[-which(apply(specm,1,sum)==0),]

gmat = matrix(0,nrow=length(specmlist),ncol=nrow(meanmat_sub),dimnames=list(names(specmlist),rownames(meanmat_sub)))
bmat = matrix(0,nrow=length(specmlist),ncol=nrow(meanmat_sub),dimnames=list(names(specmlist),rownames(meanmat_sub)))
for(j in colnames(gmat)){
  for(i in names(specmlist)){
    mod = specs[[i]]
    medium = react_id(findExchReact(mod))[grep("EX_",react_id(findExchReact(mod)))]
    testmod = changeBounds(mod, medium, lb=0, ub=1000)
    testmod = changeBounds(testmod, names(specmlist[[i]][,j]), lb=ifelse(specmlist[[i]][,j]!=0,-0.1,0), ub=1000)
    lpobj = optimizeProb(testmod)
    gmat[i,j] = lpobj@lp_obj
    testmod = changeBounds(testmod, "biomass0", lb=lpobj@lp_obj, ub=lpobj@lp_obj)
    testmod = changeObjFunc(testmod,"EX_but(e)")
    bmat[i,j] = optimizeProb(testmod)@lp_obj
    #bmat[i,j] = getFluxDist(lpobj)[which(testmod@react_id == "EX_but(e)")]
  }
}

dmat = melt(gmat)
colnames(dmat) = c("Microbe","Patient","Growth")
dmat[,"Butyrate"] = melt(bmat)$value
selnam = unlist(lapply(strsplit(as.character(dmat$Microbe),"_"),function(x){paste(substr(x[1],1,1),".",x[2]," ",paste(x[3:length(x)],collapse="_"),sep="")}))
selnam = gsub("F.cf prausnitzii_KLE1255","F.prausnitzii KLE1255",selnam)
dmat$Microbe = factor(selnam)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
dmat$Type = pclass[as.character(dmat$Patient),1]
dmat$Patient = factor(as.character(dmat$Patient), levels=as.character(pclass[unique(as.character(dmat$Patient)),2]))
ggplot(dmat,aes(y=factor(Microbe),x=factor(Patient))) +
  geom_tile(aes(fill=Butyrate), color = "black") +
  geom_rect(xmin=0,xmax=4.5,ymin=0,ymax =18,color="red",fill=NA,lwd=3) +
  #scale_fill_manual(low = "grey",high = "black")+
  #geom_point(aes(colour=Growth, size=Growth))  +
  geom_point(aes(size=Growth))  +
  scale_color_gradient(low = "grey",high = "black")+
  scale_fill_gradientn(colours=YlOrBr) +
  scale_size(range = c(1, 15)) +
  xlab("Patients") +
  theme_bw(base_size = 20) +
  theme(legend.title=element_text(size=15),
        legend.text=element_text(size=10),
        legend.key=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        axis.ticks.x = element_blank(),
        plot.title = element_blank())#15x8
# ggplot(dmat, aes(factor(Microbe), log10(Butyrate))) +
#   geom_boxplot(aes(color=Type),lwd=1.1,position=position_dodge(0.9)) +
#   #scale_y_continuous(limits = c(min(fdat$Flux),1300000)) +
#   coord_flip() +
#   scale_color_manual(values = c("firebrick2","royalblue2")) +
#   ylab("Total population butyrate flux") +
#   theme_bw(base_size = 20) +
#   theme(legend.text=element_text(size=20),
#         legend.position="none",
#         legend.key=element_blank(),
#         legend.title=element_blank(),
#         panel.border = element_rect(colour='black',size=2),
#         axis.text.y = element_text(size=25),
#         axis.text.x = element_text(size=15),
#         axis.title.x = element_text(size=25),
#         axis.title.y = element_blank(),
#         axis.ticks = element_line(size=1,color='black'),
#         axis.ticks.length=unit(0.3,'cm'),
#         plot.title = element_blank())
pbmat = bmat[,rownames(pclass)[which(pclass$V2=="Crohn's disease")]]
pgmat = gmat[,rownames(pclass)[which(pclass$V2=="Crohn's disease")]]
gsptestl = list()
bsptestl = list()
for(j in colnames(pbmat)){
  gmratp = matrix(0,nrow=length(specmlist),ncol=ncol(meanmat_sub),dimnames=list(names(specmlist),colnames(meanmat_sub)))
  bmratp = matrix(0,nrow=length(specmlist),ncol=ncol(meanmat_sub),dimnames=list(names(specmlist),colnames(meanmat_sub)))
  for(i in names(specmlist)){
    bnds = ifelse(specmlist[[i]][,j]!=0,-0.1,0)
    for(k in names(which(bnds==0))){
      bndsk = bnds
      mod = specs[[i]]
      medium = react_id(findExchReact(mod))[grep("EX_",react_id(findExchReact(mod)))]
      testmod = changeBounds(mod, medium, lb=0, ub=1000)
      bndsk[k] = -0.1
      testmod = changeBounds(testmod, names(bndsk), lb=bndsk, ub=1000)
      lpobj = optimizeProb(testmod)
      gmratp[i,k] = lpobj@lp_obj/pgmat[i,j]
      testmod = changeBounds(testmod, "biomass0", lb=lpobj@lp_obj, ub=lpobj@lp_obj)
      testmod = changeObjFunc(testmod,"EX_but(e)")
      bmratp[i,k] = optimizeProb(testmod)@lp_obj/pbmat[i,j]
    }
  }
  bmratp = bmratp[,-which(apply(bmratp,2,sum)==0)]
  gmratp = gmratp[,-which(apply(gmratp,2,sum)==0)]
  bsptestl[[j]] = bmratp
  gsptestl[[j]] = gmratp
}
pat1 = names(tail(sort(apply(bsptestl[[1]],2,sum)),5))
pat2 = names(tail(sort(apply(bsptestl[[2]],2,sum)),5))
pat3 = names(tail(sort(apply(bsptestl[[3]],2,sum)),5))
pat4 = names(tail(sort(apply(bsptestl[[4]],2,sum)),5))
bgen = unique(c(pat1,pat2,pat3,pat4))
pat1 = names(tail(sort(apply(gsptestl[[1]],2,sum)),5))
pat2 = names(tail(sort(apply(gsptestl[[2]],2,sum)),5))
pat3 = names(tail(sort(apply(gsptestl[[3]],2,sum)),5))
pat4 = names(tail(sort(apply(gsptestl[[4]],2,sum)),5))
ggen = unique(c(pat1,pat2,pat3,pat4))

ggmat = pgmat
gbmat = pbmat
for(j in colnames(ggmat)){
  for(i in names(specmlist)){
    mod = specs[[i]]
    medium = react_id(findExchReact(mod))[grep("EX_",react_id(findExchReact(mod)))]
    testmod = changeBounds(mod, medium, lb=0, ub=1000)
    bnds = ifelse(specmlist[[i]][,j]!=0,-0.1,0)
    bnds[intersect(names(bnds),bgen)] = -0.1
    testmod = changeBounds(testmod, names(bnds), lb=bnds, ub=1000)
    lpobj = optimizeProb(testmod)
    ggmat[i,j] = lpobj@lp_obj
    testmod = changeBounds(testmod, "biomass0", lb=lpobj@lp_obj, ub=lpobj@lp_obj)
    testmod = changeObjFunc(testmod,"EX_but(e)")
    gbmat[i,j] = optimizeProb(testmod)@lp_obj
    #bmat[i,j] = getFluxDist(lpobj)[which(testmod@react_id == "EX_but(e)")]
  }
}
sgmat = pgmat
sbmat = pbmat
for(j in colnames(ggmat)){
  for(i in names(specmlist)){
    mod = specs[[i]]
    medium = react_id(findExchReact(mod))[grep("EX_",react_id(findExchReact(mod)))]
    testmod = changeBounds(mod, medium, lb=0, ub=1000)
    bnds = ifelse(specmlist[[i]][,j]!=0,-0.1,0)
    bnds[intersect(names(bnds),c(bgen,ggen))] = -0.1
    testmod = changeBounds(testmod, names(bnds), lb=bnds, ub=1000)
    lpobj = optimizeProb(testmod)
    sgmat[i,j] = lpobj@lp_obj
    testmod = changeBounds(testmod, "biomass0", lb=lpobj@lp_obj, ub=lpobj@lp_obj)
    testmod = changeObjFunc(testmod,"EX_but(e)")
    sbmat[i,j] = optimizeProb(testmod)@lp_obj
    #bmat[i,j] = getFluxDist(lpobj)[which(testmod@react_id == "EX_but(e)")]
  }
}

gmatall = cbind(pgmat,ggmat,sgmat)
colnames(gmatall) = paste("P",1:ncol(gmatall),sep=".")
bmatall = cbind(pbmat,gbmat,sbmat)
colnames(bmatall) = paste("P",1:ncol(bmatall),sep=".")
dmat = melt(gmatall)
colnames(dmat) = c("Microbe","Patient","Growth")
dmat$Butyrate = melt(bmatall)$value
selnam = unlist(lapply(strsplit(as.character(dmat$Microbe),"_"),function(x){paste(substr(x[1],1,1),".",x[2]," ",paste(x[3:length(x)],collapse="_"),sep="")}))
selnam = gsub("F.cf prausnitzii_KLE1255","F.prausnitzii KLE1255",selnam)
dmat$Microbe = factor(selnam)
YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
#dmat$Type = pclass[as.character(dmat$Patient),1]
dmat$Patient = factor(as.character(dmat$Patient), levels=paste("P",1:ncol(bmatall),sep="."))
ggplot(dmat,aes(y=factor(Microbe),x=factor(Patient))) +
  geom_tile(aes(fill=Butyrate), color = "black") +
  geom_vline(xintercept = 4.5, color="black",lwd=3) +
  geom_vline(xintercept = 8.5, color="black",lwd=3) +
  #geom_rect(xmin=0,xmax=4.5,ymin=0,ymax =18,color="red",fill=NA,lwd=3) +
  #scale_fill_manual(low = "grey",high = "black")+
  #geom_point(aes(colour=Growth, size=Growth))  +
  geom_point(aes(size=Growth))  +
  scale_color_gradient(low = "grey",high = "black")+
  scale_fill_gradientn(colours=YlOrBr) +
  scale_size(range = c(1, 15)) +
  xlab("Patients") +
  theme_bw(base_size = 20) +
  theme(legend.title=element_text(size=15),
        legend.text=element_text(size=10),
        legend.key=element_blank(),
        panel.border = element_rect(colour='black',size=2),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=25),
        axis.title.y = element_blank(),
        axis.ticks = element_line(size=1,color='black'),
        axis.ticks.length=unit(0.3,'cm'),
        axis.ticks.x = element_blank(),
        plot.title = element_blank())#15x8

gsub("exchange reaction for ","",gsub(" exchange","",subsys[bgen,"Reaction.Name"]))
gsub("exchange reaction for ","",gsub(" exchange","",subsys[ggen,"Reaction.Name"]))

######################################################################################
mod = specs[[i]]
medium = react_id(findExchReact(mod))[grep("EX_",react_id(findExchReact(mod)))]
testmod = changeBounds(mod, medium, lb=0, ub=1000)
testmod = changeBounds(testmod, names(bnds), lb=bnds, ub=1000)

testmod = changeBounds(testmod, "biomass0", lb=lpobj@lp_obj, ub=lpobj@lp_obj)
testmod = changeObjFunc(testmod,"EX_but(e)")
lpobj = optimizeProb(testmod)
meanmat_pop[1,names(exflux)]

mod = specs[[i]]
testmod = changeBounds(mod, names(bnds), lb=bnds, ub=1000)
lpob = optimizeProb(testmod)
getFluxDist(lpob)[which(testmod@react_id == "EX_but(e)")]

load("P:/MAPPING/Models773/agora773_constrain.RData")
specsel = specs[msel]

dietmod = ifelse(subgrid[rownames(meanmat_sub),]<0,0,subgrid[rownames(meanmat_sub),])*10^6
#dietmod = meanmat_sub
#dietmod = subgrid[rownames(meanmat_sub),]

#dietmod = rbind(alldiet$EU,meanmat_sub)
#rownames(dietmod)[1] = "EU"
#dietmod = dietmod*10^15
#dietmod = (dietmod/apply(dietmod,1,sum))*10^14

#dietmod = ifelse(meanmat_sub<0,0,meanmat_sub)
dietmod = (dietmod/apply(dietmod,1,sum))*10^6

gmat = matrix(0,nrow=length(specsel),ncol=nrow(dietmod),dimnames=list(unlist(lapply(specsel,mod_desc)),rownames(dietmod)))
bmat = matrix(0,nrow=length(specsel),ncol=nrow(dietmod),dimnames=list(unlist(lapply(specsel,mod_desc)),rownames(dietmod)))
for(j in 1:length(specsel)){
  modelsel = specsel[[j]]
  medium = react_id(findExchReact(modelsel))[grep("EX_",react_id(findExchReact(modelsel)))]
  growths = vector()
  for(i in 1:nrow(dietmod)){
    #testmod = changeBounds(modelsel, medium, lb=0, ub=1000)
    bnds = ifelse(dietmod[i,medium]!=0,1,0)
    testmod = changeBounds(modelsel, medium, lb=-(bnds), ub=1000)
    #testmod = changeBounds(testmod, "EX_but(e)", lb=0, ub=1000)
    lpobj = optimizeProb(testmod)
    gmat[j,i] = lpobj@lp_obj
    testmod = changeBounds(testmod, "biomass0", lb=lpobj@lp_obj, ub=lpobj@lp_obj)
    testmod = changeObjFunc(testmod,"EX_but(e)")
    lpobj = optimizeProb(testmod)
    bmat[j,i] = lpobj@lp_obj
    #bmat[j,i] = getFluxDist(lpobj)[which(testmod@react_id == "EX_but(e)")]
  }
}

