#! /bin/bash
source /etc/profile

while IFS=, read col1 col2 col3
do	
wget -r 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX113/'$col1'/*' #&
done < sra_result_just_nutrition.csv