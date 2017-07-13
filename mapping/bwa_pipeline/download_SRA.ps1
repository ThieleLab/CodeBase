$sraid = import-csv  P:\MAPPING\Pediatric_crohns\sra_result_just_nutrition.csv | % {$_.ExperimentAccession}

foreach ($element in $sraid) {
	$link = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX/SRX113/SRX1133410/' + $element
	C:\wget64.exe -r 
}