library(ggplot2)
library(reshape)

setwd('P:/MAPPING/Pediatric_crohns')
filter5 <- as.character(read.table("mappstats_original_CD.txt", header=F, sep="\t")[,1])
original = filter5
filter3 = filter5
filter1 = filter5
unfilter = filter5

original <- as.character(read.table("mappstats_original.txt", header=F, sep="\t")[,1])
filter5 <- as.character(read.table("filter5/mappstats_filter5.txt", header=F, sep="\t")[,1])
filter3 <- as.character(read.table("filter3/mappstats_filter3.txt", header=F, sep="\t")[,1])
filter1 <- as.character(read.table("filter1/mappstats_filter1.txt", header=F, sep="\t")[,1])
unfilter <- as.character(read.table("unfiltered/mappstats_unfiltered.txt", header=F, sep="\t")[,1])

numpat = length(original)/12

patnames = vector()
for(i in 0:(numpat-1)){patnames = c(patnames,original[i*12+1])}

mapdat = data.frame(row.names = patnames)
for(i in 0:(numpat-1)){
  mapdat[original[i*12+1],'All reads'] = as.numeric(unlist(lapply(strsplit(original[i*12+2]," "),function(x){x[1]})))
  mapdat[original[i*12+1],'Mapped'] = as.numeric(unlist(lapply(strsplit(original[i*12+4]," "),function(x){x[1]})))
  mapdat[unfilter[i*12+1],'Unfiltered'] = as.numeric(unlist(lapply(strsplit(unfilter[i*12+4]," "),function(x){x[1]})))
  mapdat[unfilter[i*12+1],'MAPQ>0'] = as.numeric(unlist(lapply(strsplit(filter1[i*12+4]," "),function(x){x[1]})))
  mapdat[unfilter[i*12+1],'MAPQ>2'] = as.numeric(unlist(lapply(strsplit(filter3[i*12+4]," "),function(x){x[1]})))
  mapdat[unfilter[i*12+1],'MAPQ>4'] = as.numeric(unlist(lapply(strsplit(filter5[i*12+4]," "),function(x){x[1]})))
}
boxplot(as.matrix(mapdat),ylab="Number of mapped reads")
maprel = as.matrix(mapdat[,-1])
maprel = apply(maprel,2,function(x,org){x/org},org=mapdat[,1])
boxplot(maprel*100,ylab="Percentage of mapped reads")

dat = melt(maprel)
dat$X2 = factor(as.character(dat$X2), levels = c("Mapped","Unfiltered","MAPQ>0","MAPQ>2","MAPQ>4"))
gmrel = ggplot(dat,aes(factor(X2),value*100)) + 
  geom_boxplot() +
  xlab("") +
  ylab("Percentage of mapped reads") +
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
        plot.title = element_text(size=20))

dat = melt(mapdat)
dat$X2 = factor(as.character(dat$variable), levels = c("All reads","Mapped","Unfiltered","MAPQ>0","MAPQ>2","MAPQ>4"))
gmreads = ggplot(dat,aes(factor(variable),value)) + 
  geom_boxplot() +
  xlab("") +
  ylab("Number of mapped reads") +
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
        plot.title = element_text(size=20))

grid.arrange(gmreads,gmrel,ncol=1)#16x12

barplot(maprel[,'filter1']*100,horiz=T,las=2)


library(ggplot2)
library(reshape)

setwd('P:/MAPPING/Pediatric_crohns')
cdreads <- as.character(read.table("mappstats_original_CD.txt", header=F, sep="\t")[,1])

numpat = length(cdreads)/12

patnames = vector()
for(i in 0:(numpat-1)){patnames = c(patnames,cdreads[i*12+1])}

mapdat = data.frame(row.names = patnames)
for(i in 0:(numpat-1)){
  mapdat[cdreads[i*12+1],'All reads'] = as.numeric(unlist(lapply(strsplit(cdreads[i*12+2]," "),function(x){x[1]})))
  mapdat[cdreads[i*12+1],'Mapped'] = as.numeric(unlist(lapply(strsplit(cdreads[i*12+4]," "),function(x){x[1]})))
}

setwd('P:/MAPPING/Pediatric_crohns')
hreads <- as.character(read.table("mappstats_original_healthy.txt", header=F, sep="\t")[,1])

numpat = length(hreads)/12

patnames = vector()
for(i in 0:(numpat-1)){patnames = c(patnames,hreads[i*12+1])}

hmapdat = data.frame(row.names = patnames)
for(i in 0:(numpat-1)){
  hmapdat[hreads[i*12+1],'All reads'] = as.numeric(unlist(lapply(strsplit(hreads[i*12+2]," "),function(x){x[1]})))
  hmapdat[hreads[i*12+1],'Mapped'] = as.numeric(unlist(lapply(strsplit(hreads[i*12+4]," "),function(x){x[1]})))
}

alldat = rbind(mapdat,hmapdat)

cats = read.csv("P:/MAPPING/Pediatric_crohns/dysbiotic_cluster.csv",row.names=2)
catsel = data.frame(Reads=gsub(" S","S",rownames(cats)),Sample.Name=gsub(" S","S",rownames(cats)),category=cats$Disease)
catsel$category = as.character(catsel$category)
catsel$category[which(catsel$category == "Control")] = "Healthy"
catsel$category[which(catsel$category == "Crohn")] = "Crohn's disease"
catsel$category = as.factor(catsel$category)
pclass = data.frame(V2=catsel$category,V3=catsel$Reads)
rownames(pclass) = as.character(catsel$Reads)

alldat = alldat[rownames(pclass),]

mean(alldat[,2]/alldat[,1])*100
sd(alldat[,2]/alldat[,1])*100
