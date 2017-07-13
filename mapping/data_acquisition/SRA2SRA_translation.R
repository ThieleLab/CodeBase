setwd("/mnt/nfs/users/workdirs/ebauer/Mapping/pediatric/dysbiotic/ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX113")

dirs = list.files()
trans = dirs
names(trans) = dirs
for(i in dirs){
  setwd(i)
  trans[i] = list.files()[1]
  setwd("/mnt/nfs/users/workdirs/ebauer/Mapping/pediatric/dysbiotic/ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX113")
}
write.csv(trans,file="SRX2SRR.csv")
