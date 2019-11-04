#########################################################
#### look-up MWAS paper for Parkinsons Disease       ####
#### Maik Pietzner                        11/04/2019 ####
#########################################################

rm(list=ls())
setwd("V:/Programme1_DiabetesAetiology/People/Maik/side projects/replication PD and metabolomics/data/")
options(stringsAsFactors = F)
load(".RData")

##########################################
####          load MWAS results       ####
##########################################

pd.res = read.table("V:/Programme1_DiabetesAetiology/People/Maik/MWAS/regression/output/all/inc_parkinsn.txt", sep="\t", header= T)
head(pd.res)

## select requested subset of metabolites
tmp = subset(pd.res, (SUB_PATHWAY == "Methionine, Cysteine, SAM and Taurine Metabolism") | BIOCHEMICAL %in% grep("taur", pd.res$BIOCHEMICAL, value = T) | BIOCHEMICAL %in% c("3-indoxyl sulfate"))
write.table(tmp, "look.up.parkinsons.MWAS.20191104.txt", sep="\t", row.names = F)

##########################################
####     run more extensive models    ####
##########################################

## updated data files: 06/03/2019
load("../../../MWAS/data processing/data/epic_traits.RData")
load("../../../MWAS/data processing/data/epic_mwas.RData")

## create list of available metabolites
label = read.delim("../../../epic metabolite data/data/label_epic_metabolon_20181101.txt", sep="\t", header=T)
mets  = label$ID

## rename epic data set
epic <- tmp
rm(tmp)

## run the models
source("../scripts/meta_cox_reg.R")
pd.res.aug <- meta_cox_reg(epic, "ndate", "inc_parkinsn_date", "inc_parkinsn", tmp$exposure, 
                           c('age + sex + bmi + csmoke + log_crp_n'), label)
write.table(pd.res.aug, "look.up.parkinsons.MWAS.20191704.txt", sep="\t", row.names = F)

## bile acids only
res.bile <- meta_cox_reg(epic, "ndate", "inc_parkinsn_date", "inc_parkinsn", subset(label, SUB_PATHWAY %in% grep("Bile", label$SUB_PATHWAY, value=T))$ID, 
                         c('age + sex + bmi + csmoke + log_crp_n'), label)
