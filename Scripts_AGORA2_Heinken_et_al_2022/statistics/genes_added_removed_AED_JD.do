*Revision AGORA2 analyses - gene added/removed - comparison JD vs AED - accuracy drug prediction
*Johannes Hertel

*generates figures S1, S7

clear
clear mata
clear matrix
set more off
set maxvar 32000
capture log close

cd A:\AGORA_2_New\Files_for_Johannes_revision\processed_data

import excel "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\Reactions_Removed_Added.xlsx", sheet("Sheet1") firstrow
rename Model_ID model_id
gen reactions_ref=Reactions_only_in_curated+Reactions_in_both
gen reactions_draft=Reactions_only_in_draft+Reactions_in_both
rename Reactions_only_in_curated Reactions_added
rename Reactions_only_in_draft Reactions_removed
gen rxn_fold_change=reactions_ref/reactions_draft
save rxn_removed_added.dta, replace
clear

import excel "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\Genes_Removed_Added.xlsx", sheet("Sheet1") firstrow
rename Model_ID model_id
gen genes_ref=Genes_only_in_curated+Genes_in_both
gen genes_draft=Genes_only_in_draft+Genes_in_both
rename Genes_only_in_curated Genes_added
rename Genes_only_in_draft Genes_removed
gen genes_fold_change=genes_ref/genes_draft
save genes_removed_added.dta, replace
clear
/*
import delimited A:\AGORA_2_New\Final_data_for_analysis\Revision\Draft_Refined_Comparison\ReconstructionFeatures_Draft_AGORA2.txt, varnames(1) 
drop growth* atp*
rename reactions_supported rxn_supported
foreach j of varlist _all{
	rename `j' `j'_draft
	}
rename model_id_draft modelid
 egen test=seq(), by(modelid)
 drop if test!=1
 drop test
save reconstruction_features_draft.dta, replace
clear
import delimited A:\AGORA_2_New\Final_data_for_analysis\Revision\Draft_Refined_Comparison\ReconstructionFeatures_Refined_AGORA2.txt, varnames(1) 
drop growth* atp*
rename reactions_supported rxn_supported
foreach j of varlist _all{
	rename `j' `j'_ref
	}
rename model_id_ref modelid
save reconstruction_features_refined.dta, replace
clear */

use rxn_removed_added.dta
merge 1:1 model_id using genes_removed_added.dta
drop _merge
merge 1:1 model_id using Agora2_infofile.dta
drop _merge

sum Reactions_removed Reactions_added genes_fold rxn_fold, det

cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
histogram Reactions_added, fraction fcolor(navy) fintensity(50) lcolor(navy) kdensity kdenopts(lcolor(dknavy) lwidth(thick) gaussian) xtitle(Number of reactions added in refinement) graphregion(fcolor(white) lcolor(white)) saving(histogram_added_raw, replace) xsize(7) ysize(5)
histogram Reactions_removed, fraction fcolor(navy) fintensity(50) lcolor(navy) kdensity kdenopts(lcolor(dknavy) lwidth(thick) gaussian) xtitle(Reactions removed in refinement) graphregion(fcolor(white) lcolor(white)) saving(histogram_removed_raw, replace) xsize(7) ysize(5)
histogram rxn_fold_change, fraction fcolor(navy) fintensity(50) lcolor(navy) kdensity kdenopts(lcolor(dknavy) lwidth(thick) gaussian) xtitle(Net fold change in reaction numbers after refinement) graphregion(fcolor(white) lcolor(white)) saving(histogram_fold_change, replace) xsize(7) ysize(5)
graph combine histogram_added_raw.gph histogram_removed_raw.gph histogram_fold_change.gph, col(1) graphregion(fcolor(white) lcolor(white)) xsize(4) ysize(7) saving(histograms_combined.gph, replace)
graph box genes_fold_change rxn_fold_change if Reactions_added!=., ytitle(Fold change due to refinement) by(, graphregion(fcolor(white) lcolor(white))) ylabel(#3) by(Phylum, col(4)) graphregion(fcolor(white) lcolor(white)) saving(box_plots_phylum_raw.gph, replace)
graph combine histograms_combined.gph box_plots_phylum_raw.gph, col(2) graphregion(fcolor(white) lcolor(white)) xsize(10) ysize(7) saving(Figure_S1_raw.gph, replace)


clear
cd A:\AGORA_2_New\Files_for_Johannes_revision\processed_data
import delimited A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\AGORA2_CRC_Objectives_AED.txt, varnames(1) 
rename objective v1
egen test=seq(), by(v1)
drop if test==2
drop test
replace v1=substr(v1,23,.)
rename ex_dh5furafe ex_dh5fura_fcsn_fe
rename v6 ex_dh5fura_5fura_fe
rename ex_ac5asafe ex_ac5asa_5asa_fe
rename v11 ex_ac5asa_bzd_fe

foreach j of varlist ex*{
	rename `j' `j'_AED
	}
drop in 1
save AGORA2_CRC_Objectives_AED.dta, replace
clear

import delimited A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\AGORA2_CRC_Objectives_JD.txt, varnames(1) 
rename objective v1 
rename ex_dh5furafe ex_dh5fura_fcsn_fe
rename v6 ex_dh5fura_5fura_fe
rename ex_ac5asafe ex_ac5asa_5asa_fe
rename v11 ex_ac5asa_bzd_fe

foreach j of varlist ex*{
	rename `j' `j'_JD
	}
drop in 1
replace v1=substr(v1,23,. )
save AGORA2_CRC_Objectives_JD.dta, replace

merge 1:1 v1 using AGORA2_CRC_Objectives_AED.dta
 
drop if strlen(v1)>11
destring(ex_sn38fe_JD-ex_cholatefe_AED), replace
local list=""
local t=1
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
foreach j in ex_sn38fe_ ex_r406fe_ ex_5furafe_ ex_dh5fura_fcsn_fe_ ex_dh5fura_5fura_fe_ ex_dfdurife_ ex_dihydro_digoxinfe_ ex_ac5asa_5asa_fe_ ex_5asafe_ ex_ac5asa_bzd_fe_ ex_nchlphnclfe_ ex_bvufe_ ex_dopafe_ ex_mtymfe_ ex_pcresolfe_ ex_cholatefe_{
	reg `j'AED `j'JD
	display "`j'"
	display sqrt(e(r2))
	quietly test `j'JD
	display r(p)
	local ylab=subinstr(substr("`j'",4,.),"_fe__","(AED)",1)
	local xlab=subinstr(substr("`j'",4,.),"_fe__","(JD)",1)
	scatter `j'AED `j'JD, ytitle("`ylab'", size(medlarge)) xtitle("`xlab'", size(medlarge)) mcolor(navy) msymbol(smcircle_hollow) graphregion(fcolor(white) lcolor(white)) saving(g`t'.gph, replace) xsize(5) ysize(5)
	local list="`list'"+" "+"g`t'.gph"
	local t=`t'+1
	}

graph combine `list', graphregion(fcolor(white) lcolor(white)) saving(Figure_S7_raw, replace)

clear

cd A:\AGORA_2_New\Files_for_Johannes_revision\processed_data\ 

import excel "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\Overview_drug_validation.xlsx", sheet("Sheet1") firstrow

egen group_strain_species=group(Strainspecies)
gen in_silico_num=0 if Insilico=="No"
replace in_silico_num=1 if Insilico=="Yes"
gen in_vitro_num=0 if Invitro=="No"
replace in_vitro_num=1 if Invitro=="Yes"
save Overview_drug_validation.dta, replace

log using A:\AGORA_2_New\Files_for_Johannes_revision\results\logs\prediction_drug_capacities.log, replace

tab in_vitro in_silico, ex
display r(p_exact)
kap in_vitro in_silico
xtset group_strain_species
xtlogit in_vitro in_silico , or
test in_silico_num
display r(p)
log close

clear