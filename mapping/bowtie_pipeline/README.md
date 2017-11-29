# General Pipeline
First the reads need to checked for quality and trimmed, then the reads are mapped to the target sequences.

## Trimming
Trim the raw reads according to quality information in fastq files using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).

## Concatenating genomes
To map the reads to a set of target microbe genomes, the genomes need to be concatenated to one sequence.
The script *parse2one.py* takes each genome in the directory, concatenates the contigs to one continous sequence, and combines all genomes into one multi-fasta file.
In this file each genome is one independent sequence and interpreted as a chromosome in the mapping steps, thus making it possible to distinguish the mapped reads per chromosome/genome.

## Mapping
The concatenated genome is first converted into index files using the *bowtie2-build* command of the [bowtie2 package](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).
Then the trimmed reads are mapped to the index file using the [bowtie2 package](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) exporting in the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf).

## Result pre-processing
The generated SAM files are then processed to BAM files and afterwards indexed using [SAMtools](http://samtools.sourceforge.net/).

## Creating abundance information
To create a table of an abundance measure for each microbe and each patient, the number of mapped reads per genome (or coverage) per sample can be used as a proxy.
This table can be generated using *samtools idxstats* of [SAMtools](http://samtools.sourceforge.net/) or using the custom *R* script *compile_table.R* using [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html).

## Downstream analysis
Different cutoffs can be used to determine at what point the microbe is really present (number of mapped reads).
The abundance information can be integrated in community models giving an estimation on how much of one particular microbe is present in a human individual.