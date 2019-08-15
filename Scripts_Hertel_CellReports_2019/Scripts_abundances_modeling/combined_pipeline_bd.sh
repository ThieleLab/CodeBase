#!/bin/bash
source /etc/profile

module use $RESIF_ROOTINSTALL/lcsb/modules/all
module load bio/BEDTools
module load bio/BWA
module load bio/Trimmomatic


arr=(*_1.fastq)

#B=("${arr[@]:0:5}")
#B=("${arr[@]:6:10}")
#B=("${arr[@]:11:15}")
#B=("${arr[@]:16:20}")
#B=("${arr[@]:21:25}")
#B=("${arr[@]:26:30}")
#B=("${arr[@]:31:35}")
#B=("${arr[@]:36:40}")
#B=("${arr[@]:41:45}")
#B=("${arr[@]:46:50}")
#B=("${arr[@]:51:55}")
#B=("${arr[@]:56:60}")


for i in `seq 0 60`; do
file="${arr[$i]}"
name="${file%_1.*}"
reads1=$name"_1.fastq"
reads2=$name"_2.fastq"
#outpair1="$name.trimmed.1.fastq"
#outpair2="$name.trimmed.2.fastq"
outupair1="$name.unpaired.1.fastq"
outupair2="$name.unpaired.2.fastq"
reads1out="$name.fixed.1.fastq"
reads2out="$name.fixed.2.fastq"
echo $name
#java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.33.jar PE -phred33 $reads1 $reads2 $outpair1 $outupair1 $outpair2 $outupair2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#java -jar /usr/share/java/trimmomatic-0.32.jar PE -phred33 $reads1 $reads2 $outpair1 $outupair1 $outpair2 $outupair2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
#rm $outupair1
#rm $outupair2
#Fixing read orientation
$WORK/bbmap/repair.sh in1=$reads1 in2=$reads2 out1=$reads1out out2=$reads2out
#rm $outpair1
#rm $outpair2
#Map to human genome and filter reads
bwa mem -t 5 -M GCA_000001405.24_GRCh38.p9_genomic.fna $reads1out $reads2out > "$name.sam"
rm $reads1out
rm $reads2out
samtools view -bS -f 4 "$name.sam" > "$name.bam" #taking reads that did not map to human (-f 4 flag) 
samtools sort "$name.bam" "$name.sort"
bamToFastq -i "$name.sort.bam" -fq "$name.filter.1.fastq" -fq2 "$name.filter.2.fastq"
rm "$name.sam"
rm "$name.bam"
rm "$name.sort.bam"
#Map to microbe references
bwa mem -t 5 -M joined_genomes.fasta "$name.filter.1.fastq" "$name.filter.2.fastq" > "$name.sam"
samtools view -bS "$name.sam" > "$name.bam"
samtools sort "$name.bam" "$name.sort"
samtools index "$name.sort.bam"
samtools view -q 5 -b "$name.bam" > "$name.filter.bam" #q5 reads lower then quality 5 are removed
samtools sort "$name.filter.bam" "$name.filter.sort"
samtools index "$name.filter.sort.bam"
rm "$name.sam"
done
