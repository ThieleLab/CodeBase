*Performance comparison across resources
*Johannes Hertel


clear
clear mata
clear matrix
set more off
capture log close

cd A:\AGORA_2_New\Files_for_Johannes_revision\results\logs

*NJC19
log using comparison_NJC19.log, replace
import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\NJC19_summary.csv", varnames(1) 

rename organism Organism
egen group_organism=group(Organism)
rename agora2 AGORA2
rename kbase KBase
rename bigg BiGG
rename carveme CarveMe
rename magma MAGMA
rename direction Direction

gen true_value=1 if substr(AGORA2,1,2)=="TP" | substr(KBase,1,2)=="TP" | substr(BiGG,1,2)=="TP" | substr(CarveMe,1,2)=="TP" | substr(MAGMA,1,2)=="TP" | substr(gapseq,1,2)=="TP"
replace true_value=1 if (substr(AGORA2,1,2)=="FN" | substr(KBase,1,2)=="FN" | substr(BiGG,1,2)=="FN" | substr(CarveMe,1,2)=="FN" | substr(MAGMA,1,2)=="FN" | substr(gapseq,1,2)=="FN") & true_value!=1 
replace true_value=0 if substr(AGORA2,1,2)=="TN" | substr(KBase,1,2)=="TN" | substr(BiGG,1,2)=="TN" | substr(CarveMe,1,2)=="TN" | substr(MAGMA,1,2)=="TN" | substr(gapseq,1,2)=="TN" 
replace true_value=0 if (substr(AGORA2,1,2)=="FP" | substr(KBase,1,2)=="FP" | substr(BiGG,1,2)=="FP" | substr(CarveMe,1,2)=="FP" | substr(MAGMA,1,2)=="FP" | substr(gapseq,1,2)=="FP") & true_value==.
egen group_direction=group(Direction)

recode group_dir (1=1) (2=-1)
gen gold_standard= group_direction
replace gold_standard=0 if true_value==0

foreach j of varlist AGORA2 KBase BiGG CarveMe MAGMA gapseq{
	gen `j'_recoded=1 if substr(`j',2,1)=="P"
	replace `j'_recoded=0 if substr(`j',2,1)=="N"
	}
	
xtset group_organism
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	display "`j'"
	display "NUmber of models Secretion"
	egen count_models=group(Organism) if `j'_rec!=. & group_dir==1
	sum count_models if `j'_rec!=. & group_dir==1
	drop count_models
	display "NUmber of models uptake"
	egen count_models=group(Organism) if `j'_rec!=. & group_dir==-1
	sum count_models if `j'_rec!=. & group_dir==-1
	drop count_models
	}

foreach j of varlist AGORA2 KBase CarveMe MAGMA gapseq{
	*Secretion
	xtlogit true `j'_rec if group_dir==1, or 
	test `j'
	display r(p)
	*Uptake
	xtlogit true `j'_rec if group_dir==-1, or
	test `j'
	display r(p)
	}
	
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	gen `j'_right=1 if true==1 & `j'_recoded==1
	replace `j'_right=1 if true==0 & `j'_recoded==0
	replace `j'_right=0 if true==1 & `j'_recoded==0
	replace `j'_right=0 if true==0 & `j'_recoded==1
	egen acc_all_`j'=mean(`j'_right), by(group_organism)
	egen acc_uptake_`j'=mean(`j'_right) if group_dir==-1, by(group_organism)
	egen acc_secretion_`j'=mean(`j'_right) if group_dir==1, by(group_organism)
	}
	
foreach j of varlist KBase CarveMe BiGG MAGMA gapseq{
	display "AGORA2 vs `j'"
	egen sq=seq() if acc_uptake_AGORA2!=. & acc_uptake_`j'!=., by(group_organism)
	signrank acc_uptake_AGORA2=acc_uptake_`j' if sq==1
	drop sq 
	egen sq=seq() if acc_secretion_AGORA2!=. & acc_secretion_`j'!=., by(group_organism)
	signrank acc_secretion_AGORA2=acc_secretion_`j' if sq==1
	drop sq
	}

foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	display "*********"
	display "secretion"
	display "*********"
	tab true `j'_recoded if group_direction==1, ex
	display "*********"
	display "uptake"
	display "*********"
	tab true `j'_recoded if group_direction==-1, ex
	}	
	
local i=1
local m=1
gen method=.
gen accuracy=.
gen lb_accuracy=.
gen ub_accuracy=.
gen sec_up=.
label define repo 1 "AGORA2" 2 "KBase" 3 "CarveMe" 4 "BiGG" 5 "MAGMA" 6 "GapSeq"
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	foreach k in -1 1{
		replace method=`m' in `i'
		replace sec_up=`k' in `i'
		ci proportion `j'_right if group_dir==`k'
		replace accuracy=r(mean) in `i'
		replace lb_accuracy=r(lb) in `i'
		replace ub_accuracy=r(ub) in `i'
		local i=`i'+1
		}
	local m=`m'+1
	}
local i=1
local m=1

label values method repo
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
twoway (bar accuracy method if sec_up==-1, fcolor(dknavy) lcolor(dknavy) barwidth(0.7)) (rcap lb_accuracy ub_accuracy method if sec_up==-1), ytitle(Prediction accuracy) ylabel(0(0.25)1, nogrid) xlabel(1(1)6, valuelabel) title("Overall accuracy" "uptake (NCJ19)", size(medium) color(black)) legend(off) xsize(4) ysize(6) graphregion(fcolor(white) lcolor(white)) saving(accuracy_uptake_NCJ19.gph, replace)	
twoway (bar accuracy method if sec_up==1, fcolor(dknavy) lcolor(dknavy) barwidth(0.7)) (rcap lb_accuracy ub_accuracy method if sec_up==1), ytitle(Prediction accuracy) ylabel(0(0.25)1, nogrid) xlabel(1(1)6, valuelabel) title("Overall accuracy" "secretion (NCJ19)", size(medium) color(black)) legend(off) xsize(4) ysize(6) graphregion(fcolor(white) lcolor(white)) saving(accuracy_secretion_NCJ19.gph, replace)		
	
clear
log close

*Madin
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\logs
log using comparison_Madin.log, replace

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\Madin_summary.csv", varnames(1) 
rename organism Organism
egen group_organism=group(Organism)
rename agora2 AGORA2
rename kbase KBase
rename bigg BiGG
rename carveme CarveMe
rename magma MAGMA
rename direction Direction
drop BiGG

gen true_value=1 if substr(AGORA2,1,2)=="TP" | substr(KBase,1,2)=="TP" | substr(CarveMe,1,2)=="TP" | substr(MAGMA,1,2)=="TP" | substr(gapseq,1,2)=="TP" 
replace true_value=1 if (substr(AGORA2,1,2)=="FN" | substr(KBase,1,2)=="FN" | substr(CarveMe,1,2)=="FN" | substr(MAGMA,1,2)=="FN" | substr(gapseq,1,2)=="FN") & true_value!=1 
replace true_value=0 if substr(AGORA2,1,2)=="TN" | substr(KBase,1,2)=="TN" | substr(CarveMe,1,2)=="TN" | substr(MAGMA,1,2)=="TN" | substr(gapseq,1,2)=="TN" 
replace true_value=0 if (substr(AGORA2,1,2)=="FP" | substr(KBase,1,2)=="FP" | substr(CarveMe,1,2)=="FP" | substr(MAGMA,1,2)=="FP" | substr(gapseq,1,2)=="FP") & true_value==.

egen group_direction=group(Direction)
recode group_dir (1=-1)
gen gold_standard= group_direction
replace gold_standard=0 if true_value==0

foreach j of varlist AGORA2 KBase CarveMe MAGMA gapseq{
	gen `j'_recoded=1 if substr(`j',2,1)=="P"
	replace `j'_recoded=0 if substr(`j',2,1)=="N"
	}
	
xtset group_organism
foreach j of varlist AGORA2 KBase CarveMe MAGMA gapseq{
	display "`j'"
	display "NUmber of models Secretion"
	egen count_models=group(Organism) if `j'_rec!=. & group_dir==1
	sum count_models if `j'_rec!=. & group_dir==1
	drop count_models
	display "NUmber of models uptake"
	egen count_models=group(Organism) if `j'_rec!=. & group_dir==-1
	sum count_models if `j'_rec!=. & group_dir==-1
	drop count_models
	}
	
foreach j of varlist AGORA2 KBase CarveMe MAGMA gapseq{
	gen `j'_right=1 if true==1 & `j'_recoded==1
	replace `j'_right=1 if true==0 & `j'_recoded==0
	replace `j'_right=0 if true==1 & `j'_recoded==0
	replace `j'_right=0 if true==0 & `j'_recoded==1
	egen acc_all_`j'=mean(`j'_right), by(group_organism)
	egen acc_uptake_`j'=mean(`j'_right) if group_dir==-1, by(group_organism)
	egen acc_secretion_`j'=mean(`j'_right) if group_dir==1, by(group_organism)
	}
	
foreach j of varlist KBase CarveMe MAGMA gapseq{
	display "AGORA2 vs `j'"
	egen sq=seq() if acc_uptake_AGORA2!=. & acc_uptake_`j'!=., by(group_organism)
	signrank acc_uptake_AGORA2=acc_uptake_`j' if sq==1
	drop sq
	}

foreach j of varlist AGORA2 KBase CarveMe MAGMA gapseq{
	tab true `j'_recoded
	}	

local i=1
local m=1
gen method=.
gen accuracy=.
gen lb_accuracy=.
gen ub_accuracy=.
gen sec_up=.
gen BiGG_right=0
gen BiGG=0
label define repo 1 "AGORA2" 2 "KBase" 3 "CarveMe" 4 "BiGG" 5 "MAGMA" 6 "GapSeq"
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	foreach k in -1 1{
		replace method=`m' in `i'
		replace sec_up=`k' in `i'
		ci proportion `j'_right if group_dir==`k'
		replace accuracy=r(mean) in `i'
		replace lb_accuracy=r(lb) in `i'
		replace ub_accuracy=r(ub) in `i'
		local i=`i'+1
		}
	local m=`m'+1
	}
local i=1
local m=1

label values method repo
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
twoway (bar accuracy method if sec_up==-1, fcolor(dknavy) lcolor(dknavy) barwidth(0.7)) (rcap lb_accuracy ub_accuracy method if sec_up==-1), ytitle(Prediction accuracy) ylabel(0(0.25)1, nogrid) xlabel(1(1)6, valuelabel) title("Overall accuracy" "uptake (Madin)", size(medium) color(black)) legend(off) xsize(4) ysize(6) graphregion(fcolor(white) lcolor(white)) saving(accuracy_uptake_Madin.gph, replace)	
		
log close
clear
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\logs

*BacDive enzymes
log using comparison_BacDive_enzymes.log, replace

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\BacDive_Enzymes_summary.csv", varnames(1) 

rename organism Organism
egen group_organism=group(Organism)
rename agora2 AGORA2
rename kbase KBase
rename bigg BiGG
rename carveme CarveMe
rename magma MAGMA
rename direction Direction

gen true_value=1 if substr(AGORA2,1,2)=="TP" | substr(KBase,1,2)=="TP" | substr(BiGG,1,2)=="TP" | substr(CarveMe,1,2)=="TP" | substr(MAGMA,1,2)=="TP" | substr(gapseq,1,2)=="TP"
replace true_value=1 if (substr(AGORA2,1,2)=="FN" | substr(KBase,1,2)=="FN" | substr(BiGG,1,2)=="FN" | substr(CarveMe,1,2)=="FN" | substr(MAGMA,1,2)=="FN" | substr(gapseq,1,2)=="FN") & true_value!=1 
replace true_value=0 if substr(AGORA2,1,2)=="TN" | substr(KBase,1,2)=="TN" | substr(BiGG,1,2)=="TN" | substr(CarveMe,1,2)=="TN" | substr(MAGMA,1,2)=="TN" | substr(gapseq,1,2)=="TN" 
replace true_value=0 if (substr(AGORA2,1,2)=="FP" | substr(KBase,1,2)=="FP" | substr(BiGG,1,2)=="FP" | substr(CarveMe,1,2)=="FP" | substr(MAGMA,1,2)=="FP" | substr(gapseq,1,2)=="FP") & true_value==.
egen group_direction=group(Direction)

gen gold_standard= group_direction
replace gold_standard=0 if true_value==0

foreach j of varlist AGORA2 KBase BiGG CarveMe MAGMA gapseq{
	gen `j'_recoded=1 if substr(`j',2,1)=="P"
	replace `j'_recoded=0 if substr(`j',2,1)=="N"
	}
	
xtset group_organism
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	display "`j'"
	display "NUmber of models Presence"
	egen count_models=group(Organism) if `j'_rec!=. & group_dir==1
	sum count_models if `j'_rec!=. & group_dir==1
	drop count_models
	}

foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	*Presence
	xtlogit true `j'_rec if group_dir==1, or 
	test `j'
	display r(p)
	}
	
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	gen `j'_right=1 if true==1 & `j'_recoded==1
	replace `j'_right=1 if true==0 & `j'_recoded==0
	replace `j'_right=0 if true==1 & `j'_recoded==0
	replace `j'_right=0 if true==0 & `j'_recoded==1
	egen acc_all_`j'=mean(`j'_right), by(group_organism)
	egen acc_secretion_`j'=mean(`j'_right) if group_dir==1, by(group_organism)
	}
	
foreach j of varlist KBase CarveMe BiGG MAGMA gapseq{
	display "AGORA2 vs `j'"
	egen sq=seq() if acc_secretion_AGORA2!=. & acc_secretion_`j'!=., by(group_organism)
	signrank acc_secretion_AGORA2=acc_secretion_`j' if sq==1
	drop sq
	}

foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	display "*********"
	display "presence"
	display "*********"
	tab true `j'_recoded if group_direction==1, ex
	}
	
local i=1
local m=1
gen method=.
gen accuracy=.
gen lb_accuracy=.
gen ub_accuracy=.
gen sec_up=.
label define repo 1 "AGORA2" 2 "KBase" 3 "CarveMe" 4 "BiGG" 5 "MAGMA" 6 "GapSeq"
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	foreach k in 1{
		replace method=`m' in `i'
		replace sec_up=`k' in `i'
		ci proportion `j'_right if group_dir==`k'
		replace accuracy=r(mean) in `i'
		replace lb_accuracy=r(lb) in `i'
		replace ub_accuracy=r(ub) in `i'
		local i=`i'+1
		}
	local m=`m'+1
	}
local i=1
local m=1

label values method repo
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
twoway (bar accuracy method if sec_up==1, fcolor(dknavy) lcolor(dknavy) barwidth(0.7)) (rcap lb_accuracy ub_accuracy method if sec_up==1), ytitle(Prediction accuracy) ylabel(0(0.25)1, nogrid) xlabel(1(1)6, valuelabel) title("Overall accuracy" "enzymes (BacDive)", size(medium) color(black)) legend(off) xsize(4) ysize(6) graphregion(fcolor(white) lcolor(white)) saving(accuracy_secretion_raw_BacDive_enzymes.gph, replace)	

log close
clear

cd A:\AGORA_2_New\Files_for_Johannes_revision\results\logs

*Backdive metabolites


log using comparison_BacDive_metabolites.log, replace
import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\BacDive_Metabolites_summary.csv", varnames(1) 

rename organism Organism
egen group_organism=group(Organism)
rename agora2 AGORA2
rename kbase KBase
rename bigg BiGG
rename carveme CarveMe
rename magma MAGMA
rename direction Direction

gen true_value=1 if substr(AGORA2,1,2)=="TP" | substr(KBase,1,2)=="TP" | substr(BiGG,1,2)=="TP" | substr(CarveMe,1,2)=="TP" | substr(MAGMA,1,2)=="TP" | substr(gapseq,1,2)=="TP"
replace true_value=1 if (substr(AGORA2,1,2)=="FN" | substr(KBase,1,2)=="FN" | substr(BiGG,1,2)=="FN" | substr(CarveMe,1,2)=="FN" | substr(MAGMA,1,2)=="FN" | substr(gapseq,1,2)=="FN") & true_value!=1 
replace true_value=0 if substr(AGORA2,1,2)=="TN" | substr(KBase,1,2)=="TN" | substr(BiGG,1,2)=="TN" | substr(CarveMe,1,2)=="TN" | substr(MAGMA,1,2)=="TN" | substr(gapseq,1,2)=="TN" 
replace true_value=0 if (substr(AGORA2,1,2)=="FP" | substr(KBase,1,2)=="FP" | substr(BiGG,1,2)=="FP" | substr(CarveMe,1,2)=="FP" | substr(MAGMA,1,2)=="FP" | substr(gapseq,1,2)=="FP") & true_value==.
egen group_direction=group(Direction)

recode group_dir (1=1) (2=-1)
gen gold_standard= group_direction
replace gold_standard=0 if true_value==0

foreach j of varlist AGORA2 KBase BiGG CarveMe MAGMA gapseq{
	gen `j'_recoded=1 if substr(`j',2,1)=="P"
	replace `j'_recoded=0 if substr(`j',2,1)=="N"
	}
	
xtset group_organism
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	display "`j'"
	display "NUmber of models Secretion"
	egen count_models=group(Organism) if `j'_rec!=. & group_dir==1
	sum count_models if `j'_rec!=. & group_dir==1
	drop count_models
	display "NUmber of models uptake"
	egen count_models=group(Organism) if `j'_rec!=. & group_dir==-1
	sum count_models if `j'_rec!=. & group_dir==-1
	drop count_models
	}

foreach j of varlist AGORA2 KBase CarveMe MAGMA gapseq{
	*Secretion
	xtlogit true `j'_rec if group_dir==1, or 
	test `j'
	display r(p)
	*Uptake
	xtlogit true `j'_rec if group_dir==-1, or
	test `j'
	display r(p)
	}
	
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	gen `j'_right=1 if true==1 & `j'_recoded==1
	replace `j'_right=1 if true==0 & `j'_recoded==0
	replace `j'_right=0 if true==1 & `j'_recoded==0
	replace `j'_right=0 if true==0 & `j'_recoded==1
	egen acc_all_`j'=mean(`j'_right), by(group_organism)
	egen acc_uptake_`j'=mean(`j'_right) if group_dir==-1, by(group_organism)
	egen acc_secretion_`j'=mean(`j'_right) if group_dir==1, by(group_organism)
	}
	
foreach j of varlist KBase CarveMe BiGG MAGMA gapseq{
	display "AGORA2 vs `j'"
	egen sq=seq() if acc_uptake_AGORA2!=. & acc_uptake_`j'!=., by(group_organism)
	signrank acc_uptake_AGORA2=acc_uptake_`j' if sq==1
	drop sq 
	egen sq=seq() if acc_secretion_AGORA2!=. & acc_secretion_`j'!=., by(group_organism)
	signrank acc_secretion_AGORA2=acc_secretion_`j' if sq==1
	drop sq
	}

foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	display "*********"
	display "secretion"
	display "*********"
	tab true `j'_recoded if group_direction==1, ex
	display "*********"
	display "uptake"
	display "*********"
	tab true `j'_recoded if group_direction==-1, ex
	}	
	
local i=1
local m=1
gen method=.
gen accuracy=.
gen lb_accuracy=.
gen ub_accuracy=.
gen sec_up=.
label define repo 1 "AGORA2" 2 "KBase" 3 "CarveMe" 4 "BiGG" 5 "MAGMA" 6 "GapSeq"
foreach j of varlist AGORA2 KBase CarveMe BiGG MAGMA gapseq{
	foreach k in -1 1{
		replace method=`m' in `i'
		replace sec_up=`k' in `i'
		ci proportion `j'_right if group_dir==`k'
		replace accuracy=r(mean) in `i'
		replace lb_accuracy=r(lb) in `i'
		replace ub_accuracy=r(ub) in `i'
		local i=`i'+1
		}
	local m=`m'+1
	}
local i=1
local m=1

label values method repo
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
twoway (bar accuracy method if sec_up==-1, fcolor(dknavy) lcolor(dknavy) barwidth(0.7)) (rcap lb_accuracy ub_accuracy method if sec_up==-1), ytitle(Prediction accuracy) ylabel(0(0.25)1, nogrid) xlabel(1(1)6, valuelabel) title("Overall accuracy" "uptake (BacDive)", size(medium) color(black)) legend(off) xsize(4) ysize(6) graphregion(fcolor(white) lcolor(white)) saving(accuracy_uptake_BacDive.gph, replace)	
twoway (bar accuracy method if sec_up==1, fcolor(dknavy) lcolor(dknavy) barwidth(0.7)) (rcap lb_accuracy ub_accuracy method if sec_up==1), ytitle(Prediction accuracy) ylabel(0(0.25)1, nogrid) xlabel(1(1)6, valuelabel) title("Overall accuracy" "secretion (BacDive)", size(medium) color(black)) legend(off) xsize(4) ysize(6) graphregion(fcolor(white) lcolor(white)) saving(accuracy_secretion_BacDive.gph, replace)		
	
clear
log close

*graph raw
graph combine accuracy_uptake_NCJ19.gph accuracy_secretion_NCJ19.gph  accuracy_uptake_BacDive.gph accuracy_secretion_BacDive.gph accuracy_secretion_raw_BacDive_enzymes.gph accuracy_uptake_Madin.gph, xsize(25) ysize(10) row(1) graphregion(fcolor(white) lcolor(white)) saving(accuracy_combined_raw.gph, replace)		
	