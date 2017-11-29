library(vegan)
library(ggplot2)
library(Rsamtools)
library(tsne)
library(devtools)
library(digest)
source_url("https://raw.github.com/low-decarie/FAAV/master/r/stat-ellipse.R")    

bamcoverage <- function (bamfile) {
  # read in the bam file
  bam <- scanBam(bamfile)[[1]] # the result comes in nested lists
  # filter reads without match position
  ind <- ! is.na(bam$pos)
  ## remove non-matches, they are not relevant to us
  bam <- lapply(bam, function(x) x[ind])
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  ## names of the bam data frame:
  ## "qname"  "flag"   "rname"  "strand" "pos"    "qwidth"
  ## "mapq"   "cigar"  "mrnm"   "mpos"   "isize"  "seq"    "qual"
  ## construc: genomic ranges object containing all reads
  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname )
  ## returns a coverage for each reference sequence (aka. chromosome) in the bam file
  return (mean(coverage(ranges)))      
}
setwd("P:/MAPPING/Results/BAM/")
bams = list.files("P:/MAPPING/Results/BAM/")
bams = bams[grep('sorted.bam',bams)]
bams = bams[-grep('bai',bams)]

test = bamcoverage(bams[1])
coverage = matrix(0, nrow=length(test), ncol=length(bams))
colnames(coverage) = gsub("_sorted.bam","",bams)
rownames(coverage) = names(test)
coverage[names(test),1] = test
for(i in 2:length(bams)){
  print(i)
  test = bamcoverage(bams[i])
  coverage[names(test),i] = test
}
write.csv(coverage,file="P:/MAPPING/Results/coverage.csv")


coverage = as.matrix(read.csv("P:/MAPPING/Results/coverage.csv",row.names=1,header=T))
cats = read.csv("P:/MAPPING/MetaHit/spain_metahit.csv")
dis = as.character(cats$category)
names(dis) = gsub("-",".",as.character(cats$Sample.Name))
# setwd("P:/MAPPING/Results/BAM/abtables")
# org1 = read.csv(tabs[1], header=F, sep='\t',fileEncoding="UCS-2LE", row.names=1)
# abtable = org1[-nrow(org1),-3]
# for(i in 2:length(tabs)){
#   taborg = read.csv(tabs[i], header=F, sep='\t',fileEncoding="UCS-2LE", row.names=1)
#   abtable[rownames(taborg),i+1] = taborg[,2]
# }
# colnames(abtable) = c("Genome.size",gsub(".csv","",tabs))
# abtable = abtable[-nrow(abtable),]
# abtable = abtable[,sort(colnames(abtable))]
# abred = matrix(0,nrow=nrow(abtable),ncol=length(unique(unlist(lapply(strsplit(colnames(abtable),"_"),function(x){return(x[1])})))))
# rownames(abred) = rownames(abtable)
# colnames(abred) = unique(unlist(lapply(strsplit(colnames(abtable),"_"),function(x){return(x[1])})))
# abred[,1] = abtable[,1]
# j=2
# for(i in seq(2,ncol(abtable), 2)){
#   abred[,j] = abtable[,i] + abtable[,i+1]
#   j = j+1
# }
# abtable = abred
# abtable = abtable[,-which(dis[colnames(abtable)] == "UC")]
# #unique(unlist(lapply(strsplit(colnames(abtable),"_"),function(x){return(x[2])})))
# #cover <- abtable[,2:ncol(abtable)]
# cover <- abtable[,2:ncol(abtable)]*75
# for(i in 1:ncol(cover)){
#   cover[,i] = cover[,i]/abtable[,1]
# }

coverage = ifelse(is.na(coverage),0,coverage)
cover = matrix(0, nrow=nrow(coverage), ncol=ncol(coverage)/2)
rownames(cover) = rownames(coverage)
colnames(cover) = unique(unlist(lapply(strsplit(colnames(coverage),"_"),function(x){return(x[1])})))
j=1
for(i in seq(2,ncol(coverage), 2)){
  cover[,j] = (coverage[,i-1] + coverage[,i])
  j = j+1
}
cover = ifelse(cover<0.01,0,cover)
cover = apply(cover,2,function(x){x/sum(x)})
cover = cover[,-which(dis[colnames(cover)] == "UC")]
#cover = ifelse(cover<0.1,0,cover)
#cover = ifelse(is.na(cover),0,cover)
#if(length(which(apply(cover,1,sum)==0)) != 0){
#  cover = cover[-which(apply(cover,1,sum)==0),]
#}
#cover = cover[-which(rownames(cover)=="Prevotella_copri_CB7_DSM_18205"),]

# bdist = vegdist(t(cover), method="bray")
# pcoa <- cmdscale(bdist)
# plot(pcoa[,1],pcoa[,2],pch=19,col=
#        as.numeric(as.factor(dis[unlist(lapply(strsplit(rownames(pcoa),"_"),function(x){return(x[1])}))]))+1)

##############################################################
########################### plotting PCoA
##############################################################
#cover = cover[which(apply(ifelse(cover>0,1,0),1,sum)==ncol(cover)),]
#dires <- capscale(t(cover[sel,])~1, distance="bray")
dires <- capscale(t(cover)~1, distance="bray")
dires_sum <- summary(dires)

imp=15
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

fac = 2
ggplot(data=gdat, aes(x=PC1, y=PC2, colour=Phenotype)) +
  geom_hline(aes(yintercept=0), linetype="dashed", color="grey", size=1) +
  geom_vline(aes(xintercept=0), linetype="dashed", color="grey", size=1) +   
  #geom_segment(data=gdat2, aes(x=0, y=0, xend=R1*fac, yend=R2*fac), arrow=arrow(length=unit(0.2,"cm")), alpha=0.5, color="darkgrey", size=0.4)+
  #geom_text(data=gdat2, aes(x=R1*fac, y=R2*fac, label=labs), color="black", size=3) +
  
  #stat_ellipse(level = 0.95) +
  #scale_colour_manual(values=taxcols) +
  #scale_shape_manual(values=c(1,16,18,17,2)) +#c(1,13,4,3,2)
  scale_colour_manual(values = c("firebrick2","royalblue2")) +
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

bdist = vegdist(t(cover), method="bray")
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

#heatmap(cover[sel,]*100)
#heatmap(cover)

#cover[,c(rownames(gdat[which(gdat$Phenotype=="Crohn's disease"),]),rownames(gdat[which(gdat$Phenotype=="Healthy"),]))]
#View(cover[sel,c(rownames(gdat[which(gdat$Phenotype=="Crohn's disease"),]),rownames(gdat[which(gdat$Phenotype=="Healthy"),]))]*100)

##############################################################
########################### Export Data
##############################################################

#create random microbiotas
randab = matrix(0, nrow=321, ncol=3) 
for(i in 1:ncol(randab)){
  set.seed(i)
  rands <- runif(321,min=0,max=1)
  #rands <- rnorm(321,100,50)
  randab[,i] = rands/sum(rands)
}
lowab = which(apply(randab,1,sum) == 0)
if(length(lowab)!=0){randab = randab[-lowab,]}
colnames(randab) = paste(rep("R",ncol(randab)),1:ncol(randab),sep=".")
newdat = cbind(cover, randab)

alldat = newdat
seldat = newdat[sel,]
specdat = newdat[which(apply(ifelse(cover==0,0,1),1,sum)==ncol(cover)),]

write.csv(alldat,file="P:/MAPPING/Model_data/alldat.csv")
write.csv(seldat,file="P:/MAPPING/Model_data/seldat.csv")
write.csv(specdat,file="P:/MAPPING/Model_data/specdat.csv")

classes = dis[colnames(alldat)]
classes[which(is.na(classes))] = "random"
names(classes)[(ncol(cover)+1):ncol(alldat)] = colnames(alldat)[(ncol(cover)+1):ncol(alldat)]
write.csv(classes,file="P:/MAPPING/Model_data/class.csv")

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
