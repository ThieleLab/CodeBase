*Data preparation DENOPA

clear
clear mata
clear matrix
capture log close
set maxvar 32000
set more off
cd A:\Metabolomics_PD\Paper\Finalization\data\

log using A:\Metabolomics_PD\Paper\Finalization\log_files\data_preparation.log, replace

import delimited "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\Copy of 160221_DIAG_sample_shipment_Leiden.csv", delimiter(";") varnames(3) 
rename v1 ID
replace ID="" in 1
replace ID="" in 34
replace ID="00"+ID if strlen(ID)==1
replace ID="0"+ID if strlen(ID)==2
gen group="DKK" in 2/31
replace group="DKP" in 35/64
gen group_num=0 if group=="DKK"
replace group_num=1 if group=="DKP"
gen Probanden_Nr=group+ID
drop if Probanden_Nr==""
rename boxid box_1
rename v9 box_2
rename v14 box_3
foreach j of varlist box_*{
	replace `j'=substr(`j',8,4)
	destring(`j'), replace force
	}

rename baselineibblidkit01 labor_id_1	
rename fu1ibblidkit02 labor_id_2
rename fu2ibblidkit03 labor_id_3
foreach j of varlist labor_*{
	local l=strlen(`j') in 1
	local l=`l'-1
	replace `j'=substr(`j',1,9)+substr(`j',`l',2)
	}
	
keep Probanden_Nr group group_num box_* labor_*
save denopa_labor_keys_wide.dta, replace
keep Probanden_Nr labor_id_1-box_3
reshape long box_ labor_id_, i(Probanden_Nr) j(wave)
rename labor_id_ labor_id
save denopa_labor_keys_long.dta, replace

clear

import excel "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DeNoPa Disease Duration.xlsx", sheet("denopa.clinical_baseline") firstrow
gen sex=0 if Geschlecht=="Male"
replace sex=1 if Geschlecht=="Female"
rename a_Erkrankungsdauer_Monate months_of_disease
rename a_Alter age
save denopa_disease_duration.dta, replace
merge 1:1 Probanden_Nr using A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\denopa_labor_keys.dta
drop if _merge!=3
drop _merge
save denopa_merged.dta, replace

clear

local s1="Visit 1"
local s2="Visit 2"
local s3="Visit 3"

foreach j of numlist 1 2 3{
	import excel "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\Denopa.xlsx", sheet("`s`j''") firstrow
	save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_data_`j'.dta", replace
	clear
	import excel "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\Kopie von Creatinine Levels DeNoPa.xlsx", sheet("`s`j''") firstrow
	rename a_Labor crea_`j'
	save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_crea_`j'.dta", replace
	clear
	}

use "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_data_1.dta"
merge 1:1 Probanden_Nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_data_2.dta"
drop _merge
merge 1:1 Probanden_Nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_data_3.dta"
drop _merge
merge 1:1 Probanden_Nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_crea_1.dta"
drop _merge
merge 1:1 Probanden_Nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_crea_2.dta"
drop _merge
merge 1:1 Probanden_Nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\denopa_crea_3.dta"
drop _merge
destring(crea_*), replace force
foreach j of varlist a_* b_* c_*{
	local l=strlen("`j'")-2
	local nam=substr("`j'", 3,`l')
	local t=substr("`j'",1,1)
	if "`t'"=="a"{
		rename `j' `nam'_1
		}
	if "`t'"=="b"{
		rename `j' `nam'_2
		}
	if "`t'"=="c"{
		rename `j' `nam'_3
		}
	}
	
reshape long crea_ UPDRS_I_sum_ UPDRS_II_sum_ UPDRS_III_sum_ UPDRS_IV_sum_ UPDRS_sum_, i(Probanden_Nr) j(wave)
save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\denopa_scores.dta", replace

clear

local s1="Bile Acids data"
local s2="Bile Acids data caution"
local s3="Oxidative Stress High pH data"
local s4="Ox. Stress High pH data_Caution"
local s5="Oxidative Stress Low pH data"
local s6="Ox. Stress Low pH data_Caution"
local s7="Acylcarnitine data"
local s8="Acylcarnitine data with caution"
local s9="Amines data"
local s10="Amines data with caution"
local s11="Organic acids data"
local s12="Organic acids data with caution"
local s13="PosLip_NonTGs " 
local s14="PosLip_NonTGs Caution"
local s15="PosLip_TGs"
local s16="PosLip_TGs with caution"
local s17="Neg Lipids"
local s18="Neg Lipids Caution"

local r1="L"
local l1=182
local r2="G"
local l2=182
local r3="S"
local l3=177
local r4="H"
local l4=177
local r5="D"
local l5=177
local r6="I"
local l6=177
local r7="X"
local l7=182
local r8="E"
local l8=182
local r9="AU"
local l9=182
local r10="H"
local l10=182
local r11="J"
local l11=180
local r12="D"
local l12=180
local r13="BY"
local l13=181
local r14="F"
local l14=181
local r15="AI"
local l15=181
local r16="I"
local l16=181
local r17="R"
local l17=182
local r18="H"
local l18=182


foreach j of numlist 1 2 to 18{
	import excel "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\NMC-17-004_DeNoPa_Data - V03_2.xlsx", sheet("`s`j''") cellrange(A3:`r`j''`l`j'') firstrow
	rename Original labor_id
	save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\met_data_`j'.dta", replace
	clear
	}
	
use "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\met_data_1.dta", clear

foreach j of numlist 2 3 to 18{
	merge 1:1 labor_id using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\temp\met_data_`j'.dta"
	drop _merge
	}
gen wave=substr(labor_id,10,2)
destring(wave), replace force
save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\denopa_metabolomics.dta", replace

clear

local k=1
foreach j in baseline first_visit second_visit{
	import delimited "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa.clinical_`j'.csv", delimiter("|") varnames(1)
	gen wave=`k'
	capture rename id Probanden_Nr
	save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_`j'.dta", replace
	local k=`k'+1
	clear
	}

use "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\denopa_metabolomics.dta"
merge 1:1 labor_id using "denopa_labor_keys_long.dta"
drop _merge
merge m:1 Probanden_Nr using "denopa_merged.dta"
drop _merge
merge 1:1 Probanden_Nr wave using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\denopa_scores.dta"
drop _merge
gen pr_id=substr(Probanden_Nr, 4,3)
destring(pr_id), replace force
gen pr_id_n=pr_id
replace pr_id=pr_id+1000 if group_num==1
drop if labor_id==""
save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\denopa_final_data_long.dta", replace
drop Id NMCname a_Gruppe labor_id_1- box_3 pr_id_n
reshape wide labor_id box_ CA-FA_240 crea_ a_MDS_UPDRS_H_Y a_MDS_UPDRS_sum UPDRS_I_sum_-UPDRS_sum_, i(pr_id) j(wave)
save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\denopa_final_data_wide.dta", replace

clear

local k=1
foreach j in baseline first_visit second_visit{
	import delimited "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa.clinical_`j'.csv", delimiter("|") varnames(1)
	gen wave=`k'
	capture rename id Probanden_Nr
	save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_`j'.dta", replace
	local k=`k'+1
	clear
	}
use "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_baseline.dta"
append using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_first_visit.dta",force
append using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_second_visit.dta",force
replace Probanden_Nr=probanden_nr if wave==3
keep Probanden_Nr wave a_h_y b_h_y c_h_y a_mds_updrs_h_y b_mds_updrs_h_y c_mds_updrs_h_y *piribedile *pramipexole *ropinirole *rotigotine *azilect* b_l_dopa_agonisten_einnahme c_agonisten_einnahme a_labor_kreatinin_creas_krea b_labor_kreatinin_creas_krea c_labor_kreatinin_creas_krea a_labor_gamma_gt_ggt b_labor_gamma_gt_ggt c_labor_gamma_gt_ggt a_labor_cholesterin_chol b_labor_cholesterin_chol c_labor_cholesterin_chol *labor_triglyzeride_trig b_l_dopa_dosis_mg *crf_gewicht a_crf_groesse a_anderediagnosen b_anderediagnosen c_anderediagnosen a_labor_crp_mg_dl b_labor_crp_mg_dl c_labor_crp_mg_dl a_bdi_sum b_bdi_sum c_bdi_sum b_levodopa_aequivalenzdosis_dopa b_levodopa_aequivalenzdosis_ande b_levodopa_dosis_mg b_levodopa_einnahme c_levodopa_aequivalenzdosis_agon c_levodopa_aequivalenzdosis_ande c_levodopa_dosis_mg c_levodopa_einnahme c_levodopa_gesamt_aequivalenzdos

*renaming drugs

foreach j in piribedile pramipexole ropinirole rotigotine{
	gen `j'=0 if wave==1
	replace `j'=1 if b_dopa_agonisten_`j'=="true" & wave==2
	replace `j'=0 if b_dopa_agonisten_`j'=="false" & wave==2
	replace `j'=1 if c_agonisten_`j'==1 & wave==3
	replace `j'=0 if c_agonisten_`j'==0 & wave==3
	}
gen azilect=0 if wave==1
replace azilect=1 if b_azilect=="true" & wave==2
replace azilect=0 if b_azilect=="false" & wave==2
replace azilect=1 if c_azilect=="true" & wave==3
replace azilect=0 if c_azilect=="false" & wave==3

gen dopa_agonisten_intake=0 if wave==1
replace dopa_agonisten_intake=1 if b_l_dopa_agonisten_einnahme=="true" & wave==2
replace dopa_agonisten_intake=0 if b_l_dopa_agonisten_einnahme=="false" & wave==2
replace dopa_agonisten_intake=1 if c_agonisten_einnahme=="true" & wave==3
replace dopa_agonisten_intake=0 if c_agonisten_einnahme=="false" & wave==3


foreach j in anderediagnosen h_y mds_updrs_h_y labor_crp_mg_dl bdi_sum labor_triglyzeride_trig crf_gewicht labor_kreatinin_creas_krea labor_gamma_gt_ggt labor_cholesterin_chol{
	gen `j'=a_`j' if wave==1
	replace `j'=b_`j' if wave==2
	replace `j'=c_`j' if wave==3
	}
foreach j in levodopa_dosis_mg levodopa_einnahme{
	gen `j'=b_`j' if wave==2
	replace `j'=c_`j' if wave==3
	}
gen levo_equivalent=b_l_dopa_dosis_mg
replace levo_equivalent=c_levodopa_gesamt_aequivalenzdos if wave==3

save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\clinic_extract_data.dta", replace
merge m:m Probanden_Nr wave using "A:\Metabolomics_PD\denopa_final_data_long.dta", update
keep if _merge==3 | _merge==4
drop _merge
merge m:1 Probanden_Nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_baseline.dta", keepusing(a_crf_groesse) update
keep if _merge==3 | _merge==4
drop _merge
*save "A:\Metabolomics_PD\denopa_final_data_long.dta", replace
save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\denopa_final_data_long.dta", replace

clear

use "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\denopa_final_data_wide.dta", replace
merge 1:1 Probanden_Nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_first_visit.dta", keepusing(*h_y *piribedile *pramipexole *ropinirole *rotigotine *azilect* b_l_dopa_agonisten_einnahme b_labor_triglyzeride_trig b_l_dopa_dosis_mg b_crf_gewicht b_anderediagnosen b_labor_crp_mg_dl b_bdi_sum b_levodopa_aequivalenzdosis_dopa b_levodopa_aequivalenzdosis_ande b_levodopa_dosis_mg b_levodopa_einnahme)
keep if _merge==3
drop _merge
rename Probanden_Nr probanden_nr
merge 1:1 probanden_nr using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_denopa_trend\DENOPA_all\denopa_clinical_second_visit.dta", keepusing(*h_y *piribedile *pramipexole *ropinirole *rotigotine *azilect* c_agonisten_einnahme c_labor_triglyzeride_trig c_crf_gewicht c_anderediagnosen c_labor_crp_mg_dl c_bdi_sum c_levodopa_gesamt_aequivalenzdos c_levodopa_aequivalenzdosis_ande c_levodopa_dosis_mg c_levodopa_einnahme)
keep if _merge==3
drop _merge
rename probanden_nr Probanden_Nr
*save "A:\Metabolomics_PD\denopa_final_data_wide.dta", replace
save "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\denopa_final_data_wide.dta", replace

clear

log close 

clear


