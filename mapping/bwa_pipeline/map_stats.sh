
#mapping statistics for original reads
for i in *.bam; do
	name="${i%.bam}"
	echo $name >> mappstats_original.txt
	samtools flagstat $i >> mappstats_original.txt
done

#mapping statistics for MAPQ=5 filtered reads
for i in *.filter.sort.bam; do
	name="${i%.filter.*}"
	echo $name >> mappstats_filter5.txt
	samtools flagstat $i >> mappstats_filter5.txt
done

#mapping statistics for MAPQ=1 filtered reads
for i in *.filter1.sort.bam; do
	name="${i%.filter.*}"
	echo $name >> mappstats_filter1.txt
	samtools flagstat $i >> mappstats_filter1.txt
done

#mapping statistics for MAPQ=3 filtered reads
for i in *.filter3.sort.bam; do
	name="${i%.filter.*}"
	echo $name >> mappstats_filter3.txt
	samtools flagstat $i >> mappstats_filter3.txt
done

#mapping statistics for unfiltered reads
for i in *.unfiltered.sort.bam; do
	name="${i%.unfiltered.*}"
	echo $name >> mappstats_unfiltered.txt
	samtools flagstat $i >> mappstats_unfiltered.txt
done

