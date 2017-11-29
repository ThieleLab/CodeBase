#! /bin/bash
source /etc/profile

module use $RESIF_ROOTINSTALL/lcsb/modules/all
module load bio/BEDTools
module load bio/BWA
module load bio/Trimmomatic

for file in *1.fastq; do
name="${file%_1.*}"
reads1=$name"_1.fastq"
reads2=$name"_2.fastq"
outpair1="$name.trimmed.1.fastq"
outpair2="$name.trimmed.2.fastq"
outupair1="$name.unpaired.1.fastq"
outupair2="$name.unpaired.2.fastq"
reads1out="$name.fixed.1.fastq"
reads2out="$name.fixed.2.fastq"
echo $name
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.33.jar PE -phred33 $reads1 $reads2 $outpair1 $outupair1 $outpair2 $outupair2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
rm $outupair1
rm $outupair2
#Fixing read orientation
$WORK/bbmap/repair.sh in1=$outpair1 in2=$outpair2 out1=$reads1out out2=$reads2out
rm $outpair1
rm $outpair2
#Map to human genome and filter reads
bwa mem -t 11 -M GCA_000001405.24_GRCh38.p9_genomic.fna $reads1out $reads2out > "$name.sam"
rm $reads1out
rm $reads2out
samtools view -bS -f 4 "$name.sam" > "$name.bam"
samtools sort "$name.bam" "$name.sort"
bamToFastq -i "$name.sort.bam" -fq "$name.filter.1.fastq" -fq2 "$name.filter.2.fastq"
rm "$name.sam"
rm "$name.bam"
rm "$name.sort.bam"
#Map to microbe references
bwa mem -t 11 -M joined_genomes.fasta "$name.filter.1.fastq" "$name.filter.2.fastq" > "$name.sam"
samtools view -bS "$name.sam" > "$name.bam"
samtools sort "$name.bam" "$name.sort"
samtools index "$name.sort.bam"
samtools view -q 1 -b "$name.bam" > "$name.filter.bam"
samtools sort "$name.filter.bam" "$name.filter.sort"
samtools index "$name.filter.sort.bam"
rm "$name.sam"
done
