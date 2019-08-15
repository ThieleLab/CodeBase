*Master-File for DeNoPa metabolome paper
*Performing this scripts performs the whole statistical analyses for the paper "Integrated Analyses of Microbiome and Longitudinal Metabolome Data Reveal Microbial-Host Interactions on Sulfur Metabolism in Parkinson’s Disease"
*Johannes Hertel, 2018


*1.Data preparation: Merging of data-sets, restructuring of data-sets for longitudinal data-analyses
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\Scripts"
do denopa_preparation_finalization.do

*2.Longitudinal Analyses
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\Scripts"
do denopa_standard_longitudinal_finalization.do

*3.Prediction of disease progression
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\Scripts"
do denopa_progression_finalization.do

*4.Generating figures
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\Scripts"
do script_figures_raw_final.do

*5.Bivariate analyses. Be careful: needs around a day to complete
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\Scripts"
*do denopa_twostage_finalization.do

*6 Bedarf et al. Metagenomic analyses + flux analyses
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\Scripts"
do bedarf_metagenomics_finalization.do

 