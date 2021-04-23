### Readme for scripts in the folder Simulations_Heinken_Microbiome_2019 ###
### Almut Heinken, 02/2018 ###

This folder contains scripts needed to execute the data anlysis performed in Heinken et al., Systematic assessment of secondary bile acid metabolism in gut microbes reveals distinct metabolic capabilities in inflammatory bowel disease (Microbiome. 2019 May 15;7(1):75.).
Please cite the above reference if you use these scripts.

Includes the following steps:
Calculation of reaction abundances
Creation of output tables summarizing the results shown in the paper
Statistical analysis
Analysis of shadow prices

Please note that pairwise modeling results, strain-level flux contributions,and shadow prices are provided already extracted from the computed flux solutions. The source flux solutions are several GB of data and thus could not be provided.

To execute the analysis, please first download AGORA version 1.02 (unconstrained) from https://vmh.uni.lu/#downloadview and place the reconstructions in mat file format into a subfolder named "AGORA_models" in the folder Simulations_Heinken_Microbiome_2019.
Afterwards, execute the script "executeDataAnalysis".
