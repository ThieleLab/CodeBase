#! /bin/bash
source /etc/profile

module use $RESIF_ROOTINSTALL/lcsb/modules/all
module load bio/BEDTools
module load bio/BWA

for file in *.sam; do
name="${file%.*}"
bamfile="$name.bam"
bamsort="$name.sort"
echo $name

samtools view -bS -f 4 $file > $bamfile
samtools sort $bamfile $bamsort
bamToFastq -i "$bamsort.bam" -fq "$name.fixed.1.fastq" -fq2 "$name.fixed.2.fastq"

#rm $file
done