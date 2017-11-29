#! /bin/bash
source /etc/profile

while IFS=, read col1 col2 col3
do	
$WORK/sratoolkit.2.8.1-2-ubuntu64/bin/fastq-dump --outdir $WORK/Mapping/pediatric --split-files 'ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX113/'$col1'/'$col2'/'$col2'.sra'
done < sra_result_just_nutrition.csv
