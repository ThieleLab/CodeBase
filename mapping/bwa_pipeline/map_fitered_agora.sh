#! /bin/bash
source /etc/profile

module use $RESIF_ROOTINSTALL/lcsb/modules/all
module load bio/BEDTools
module load bio/BWA

for file in *.fixed.1.fastq; do
name="${file%.fixed.1.*}"
reads1="$name.fixed.1.fastq"
reads2="$name.fixed.2.fastq"
echo $name

bwa mem -t 11 -M joined_genomes.fasta $reads1 $reads2 > "$name.sam"
samtools view -bS "$name.sam" > "$name.bam"
samtools sort "$name.bam" "$name.sort"
samtools index "$name.sort.bam"
samtools view -q 1 -b "$name.bam" > "$name.filter.bam"
samtools sort "$name.filter.bam" "$name.filter.sort"
samtools index "$name.filter.sort.bam"

rm "$name.sam"
done

