#! /bin/bash
source /etc/profile

module use $RESIF_ROOTINSTALL/lcsb/modules/all
module load bio/BEDTools
module load bio/BWA

for file in fix_*.1.fastq; do
	name="${file%.1.*}"
	reads1="$name.1.fastq"
	reads2="$name.2.fastq"
	outsam="$name.sam"
	echo $name
	bwa mem -t 11 -M GCA_000001405.24_GRCh38.p9_genomic.fna $reads1 $reads2 > $outsam
done
