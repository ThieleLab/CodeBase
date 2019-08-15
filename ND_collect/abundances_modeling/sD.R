
tags=c("PD","Control")

c1=as.character(patstat[which(patstat[,2]==tags[2]),1])
c2=as.character(patstat[which(patstat[,2]==tags[1]),1])

sip=diversity(data[c1,])
sic=diversity(data[c2,])

sipm=cbind(sip,tags[2])
sicm=cbind(sic,tags[1])

divstat=as.data.frame(rbind(sipm,sicm))
colnames(divstat)=c("val","id")