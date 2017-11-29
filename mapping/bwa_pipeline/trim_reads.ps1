$array = @(
"SRR579278_1.fastq",
"SRR579290_2.fastq",
"SRR579277_2.fastq",
"SRR579276_2.fastq",
"SRR579290_1.fastq",
"SRR579274_2.fastq",
"SRR579291_1.fastq",
"SRR579274_1.fastq",
"SRR579280_1.fastq",
"SRR579280_2.fastq",
"SRR579275_1.fastq",
"SRR579281_2.fastq",
"SRR579277_1.fastq",
"SRR579292_2.fastq",
"SRR579275_2.fastq",
"SRR579279_2.fastq",
"SRR579292_1.fastq",
"SRR579281_1.fastq",
"SRR579276_1.fastq",
"SRR579291_2.fastq",
"SRR579278_2.fastq",
"SRR579279_1.fastq")

foreach ($element in $array) {
	$element
	$read1 = 'P:\MAPPING\MetaHit\raw_reads\' + $element + '_1.fastq'
	$read2 = 'P:\MAPPING\MetaHit\raw_reads\' + $element + '_2.fastq'
	$outpair1 = 'P:\MAPPING\MetaHit\raw_reads\trimmed\paired\output_paired_' + $element + '_1.fastq'
	$outupair1 = 'P:\MAPPING\MetaHit\raw_reads\trimmed\unpaired\output_unpaired_' + $element + '_1.fastq'
	$outpair2 = 'P:\MAPPING\MetaHit\raw_reads\trimmed\paired\output_paired_' + $element + '_2.fastq'
	$outupair2 = 'P:\MAPPING\MetaHit\raw_reads\trimmed\unpaired\output_unpaired_' + $element + '_2.fastq'

	java -jar C:\Trimmomatic-0.35\trimmomatic-0.35.jar PE -phred33 $read1 $read2 $outpair1 $outupair1 $outpair2 $outupair2 ILLUMINACLIP:TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
}