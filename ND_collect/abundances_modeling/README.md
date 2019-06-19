# Computing relative abundances and modeling 

## Introduction

The following files represent the scripts created for processing 16S rRNA gene sequencing data obtaining relative abundance tables, 
personalized microbiome metabolic models and relative metabolite resolved secretion profiles (NMPCs).

**WARNING:** Please take into consideration only the files listed in this
document. Everything present in the folder but not listed and explained in this
document is to be considered not relevant.

## Requirements

The following study requires `MeFiT`, `SPINGO classifier`, `R programming`, `Matrix Laboratory` the `Parallel Computing Toolbox`,`Microbiome modeling Toolbox`,
as well as, the COBRA Toolbox (commit:de88b809b11fbaf894432066a5996fa5555273d5 )to be installed. Please refer to
the [installation instructions](https://opencobra.github.io/cobratoolbox/stable/installation.html).
The usage of `ILOG CPLEX` solver is strongly advised to obtain the best speed performance
(required for the fastFVA.m function).

## Main Folder Structure and Files

The following scripts werer used to do the different computation steps of this work.

| Script                                         | Purpose                                                                |
| -----------------------------------------------|------------------------------------------------------------------------|
| United_workflow.sh                             | *to prepare fastq reads and do classification* |
| Spingo.csv.R                                   | *to parse classification results, create relative abundance tables* |
| tresholdBatchSplitter.m                        | *applay presence treshold and create metabolic modeling inputs*|
| Starter.m                                      | *script to initiate the mgPipe module of the Microbiome modeling toolbox*|


## Funding

This study was funded by grants from the Luxembourg National Research Fund (FNR) within the National Centre of Excellence in Research (NCER) 
on Parkinson’s disease (FNR/ NCER13/BM/11264123) and the PEARL programme (FNR/P13/6682797 to RK), and by the European Research Council (ERC) 
under the European Union’s Horizon 2020 research and innovation programme (grant agreement No 757922) and (grant agreement No 692320).

## Author & Documentation Date

*Federico Baldini, 28.05.19*

*Luxembourg Centre for Systems Biomedicine, University of Luxembourg, Campus Belval, Esch-sur-Alzette, Luxembourg*

*[federico.baldini@uni.lu](federico.baldini@uni.lu)*


