# Quality contol and preparations
First the reads need to checked for quality and trimmed, then the reads are mapped to the target sequences.

## Trimming
Trim the raw reads according to quality information in fastq files using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). The commands and their parameters are given in *trim_reads.ps1*.

## Concatenating genomes
To map the reads to a set of target microbe genomes, the genomes need to be concatenated to one sequence.
The script *parse2one.py* takes each genome in the directory, concatenates the contigs to one continous sequence, and combines all genomes into one multi-fasta file.
In this file each genome is one independent sequence and interpreted as a chromosome in the mapping steps, thus making it possible to distinguish the mapped reads per chromosome/genome. The commands and their parameters are given in *parse2one.py*.

## Fixing and combining read files
To combine reads from multiple libraries into single files, the files are concatenated using the commands in *catfiles.py*. In some cases the paired reads need to be sorted to ensure that each pair is in the same position in the different files. This can be done using the *repair.sh* script in the [BBMap package](https://sourceforge.net/projects/bbmap/). The corresponding code can be found in *fix_reads.sh*.

# General pipeline
After the quality control, human reads need to be filtered out and reads neew to be mapped to the concatenated microbe genomes.

## Filtering of human reads
The most recent sequence of the human genome needs to be downloaded [here](https://www.ncbi.nlm.nih.gov/grc). The genome fasta file needs to be converted into index files using *bwa index* of the [bwa package](https://sourceforge.net/projects/bio-bwa/files). Then the trimmed reads are mapped to the index file using the [bwa package](https://sourceforge.net/projects/bio-bwa/files) exporting in the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf). The generated SAM files are then processed into BAM files with [SAMtools](http://samtools.sourceforge.net/) while fitering out unmapped reads. Using *bamToFastq* from the [BEDtools package](http://bedtools.readthedocs.io/en/latest/), the resulting BAM files are then converted into paired fastq files. The commands and their parameters are given in *bwa_human_filter.sh* and *filter_human_reads.sh*.

## Mapping to concatenated microbe genomes
The concatenated genome is first converted into index files using the *bwa index* of the [bwa package](https://sourceforge.net/projects/bio-bwa/files). Then the trimmed and filtered reads are mapped to the index files using the [bwa package](https://sourceforge.net/projects/bio-bwa/files) exporting in the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf). The generated SAM files are then processed to BAM files and afterwards indexed using [SAMtools](http://samtools.sourceforge.net/). The commands and their parameters are given in *map_fitered_agora.sh*.

## Result pre-processing and filtering
To filter cross-mapped reads as well as reads with low mapping quality, the BAM files are filtered by the MAPQ quality score with [SAMtools](http://samtools.sourceforge.net/). The commands and their parameters are given in *samtools_filter.sh*.

## Creating abundance information
To create a table of an abundance measure for each microbe and each patient, the number of mapped reads per genome (or coverage) per sample can be used as a proxy.
This table can be generated using *samtools idxstats* of [SAMtools](http://samtools.sourceforge.net/) or using the custom *R* script *compile_table.R* using [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html).

## Downstream analysis
Different cutoffs can be used to determine at what point the microbe is really present (number of mapped reads).
The abundance information can be integrated in community models giving an estimation on how much of one particular microbe is present in a human individual.