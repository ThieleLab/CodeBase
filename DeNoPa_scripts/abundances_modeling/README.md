# Computing relative abundances and modeling 

## Introduction

The following files represent the scripts and the related files created for obtaining 
sample information and reference mapping of the samples deposited in the European Bioinformatics 
Institute-Sequence Read Archive database, under accession number ERP019674. The related study can 
be found at https://doi.org/10.1186/s13073-017-0428-y.

**WARNING:** Please take into consideration only the files listed in this
document. Everything present in the folder but not listed and explained in this
document is to be considered not relevant.

## Requirements

The following study requires `R programming`, `Matrix Laboratory` the `Parallel Computing Toolbox`,`Microbiome modeling Toolbox`,
as well as, the COBRA Toolbox to be installed. Please refer to
the [installation instructions](https://opencobra.github.io/cobratoolbox/stable/installation.html).
The usage of `ILOG CPLEX` solver is strongly advised to obtain the best speed performance
(required for the fastFVA.m function).

The steps and necessary files for the reference mapping were the same used for DOI:10.1038/s41540-018-0063-2 
and deposited in [bwa_pipeline](https://github.com/ThieleLab/CodeBase/tree/master/mapping/bwa_pipeline.html).

## Main Folder Structure and Files

The following scripts were used to do the different computation steps of this work.

| Script                                         | Purpose                                                                |
| -----------------------------------------------|------------------------------------------------------------------------|
| transITableCreator.R                           | *to extract samples information (parkinson or control) from downloaded .xml file* |
| conc.sh                                        | *to concatenate different fastq files together, rename, unzip and move to different folder* |
| combined_pipeline_bd.sh                        | *to do reference mapping (human filtering and removal of low quality aligned sequeces)* |
| compute_coverage.R                             | *to coverage from mapping* |
| compute_relAbb.R                               | *from coverage to relative abundances* |
| adaptVMHDietToAGORA.m                          | *convert a specific diet from VMH into an AGORA and add essential elements to microbiota models*|


## Steps done

The fastq files were downloaded from EBI, and multiple fastq files correspondent to the same individual were joint together. 
Sample information was also downloaded from EBI. Reference genome assembly and mapping was done following the steps and the scripts in [bwa_pipeline](https://github.com/ThieleLab/CodeBase/tree/master/mapping/bwa_pipeline.html) and as described in . 
Relative abundances were computed using the scrips provided here which represents a customized version of the scripts in [bwa_pipeline](https://github.com/ThieleLab/CodeBase/tree/master/mapping/bwa_pipeline.html). 
Relative abundances were the input for the modeling and simulation part done using 
the [mgPipe] (https://github.com/opencobra/cobratoolbox/tree/develop/src/analysis/multiSpecies/microbiomeModelingToolbox/mgPipe) module of the Microbiome Modeling Toolbox. 



## Funding

This study received funding from the Luxembourg National Research Fund(FNR), through the ATTRACT programme (FNR/A12/01), and the OPEN
grant (FNR/O16/11402054), as well as the European Research Council(ERC) under the European Union’s Horizon 2020 research and innovation
programme (grant agreement No 757922).

## Author & Documentation Date

*Federico Baldini, 24.09.18*

*Luxembourg Centre for Systems Biomedicine, University of Luxembourg, Campus Belval, Esch-sur-Alzette, Luxembourg*

*[federico.baldini@uni.lu](federico.baldini@uni.lu)*


