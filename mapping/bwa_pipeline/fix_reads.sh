#! /bin/bash
source /etc/profile

module use $RESIF_ROOTINSTALL/lcsb/modules/all
module load bio/BEDTools
module load bio/BWA

for file in *.1.fastq; do
	name="${file%.1.*}"
	reads1="$name.1.fastq"
	reads2="$name.2.fastq"
	reads1out="fix_$name.1.fastq"
	reads2out="fix_$name.2.fastq"
	echo $name
	/work/users/ebauer/bbmap/repair.sh in1=$reads1 in2=$reads2 out1=$reads1out out2=$reads2out
done
