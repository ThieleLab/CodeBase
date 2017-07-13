#! /bin/bash
source /etc/profile

module use $RESIF_ROOTINSTALL/lcsb/modules/all
module load bio/BEDTools
module load bio/BWA

for file in /work/users/ebauer/Mapping/reads/Mapped_bam/*.bam; do
name="${file%.bam}"; echo $name; samtools sort $file "$name.unfiltered.sort" &
done

for file in /work/users/ebauer/Mapping/reads/Mapped_bam/*.bam; do
name="${file%.bam}"; echo $name; samtools view -q 1 -b $file > "$name.filter1.bam"; samtools sort "$name.filter1.bam" "$name.filter1.sort" &
done

for file in /work/users/ebauer/Mapping/reads/Mapped_bam/*.bam; do
name="${file%.bam}"; echo $name; samtools view -q 3 -b $file > "$name.filter3.bam"; samtools sort "$name.filter3.bam" "$name.filter3.sort" &
done