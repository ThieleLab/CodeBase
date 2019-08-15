clear
clear mata
clear matrix
capture log close
set maxvar 32000
set more off

log using "A:\Metabolomics_PD\Paper\Finalization\log_files\figure_finalization.log", replace


********************************************
***Figures from metabolomic data************
********************************************


use "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\denopa_final_data_long.dta", replace
xtset pr_id wave
drop NMCname
set obs 4000

*Outlier filtering	
foreach j of varlist CA-FA_240{
	local l=strlen("`j'")
	if `l'>10{
		local k=substr("`j'",1,12)
		rename `j' `k'
		}
	}
*log-trafo and outlier definition

foreach j of varlist CA-FA_240{
	gen ln_`j'=ln(`j')
	egen std_`j'=std(ln_`j')
	gen filter_`j'=0
	replace filter_`j'=1 if abs(std_`j')>4
	}

gen Methoxy_bin=0 if ln_Methoxytyros<-4
replace Methoxy_bin=1 if ln_Methoxytyros>-4 & ln_Methoxytyros<999
replace levodopa_dosis_mg=0 if Methoxy_bin==0
gen levodopa_dos2=levodopa_dosis_mg
replace levodopa_dosis_mg=0 if Methoxy_bin==0 
replace levo_equi=0 if (wave==1 & Methoxy_bin==0) | group_num==0
replace azilect=0 if (wave==1) | group_num==0
replace dopa_agonisten_intake=0 if (wave==1) | group_num==0
replace pramipexole=0 if (wave==1) | group_num==0 | pramipexole==.
gen levo_diff=levo_equi-levodopa_dosis
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\figure\raw"
gen var_=""


***********************************************
***Figure 1B***********************************
***********************************************
gen group_2=group_num
replace group_2=2 if levodopa_dosis_mg>0 & levodopa_dosis_mg<9999
replace group_2=. if wave>1 & levodopa_dosis_mg>9999
graph box ln_Methoxytyros, over(group_2)
graph save figure_1B_raw.gph, replace
graph export figure_1B_raw.tif, replace

***********************************************
***Figure 2C***********************************
***********************************************

capture gen ind_w=.
capture gen ind_g=.
capture gen b_coeff=.
capture gen CI_l=.
capture gen CI_h=.


local i=1

foreach k in 1 2 3 4 5{ 
	foreach l in 0 1{
		if `k'==4{
			quietly xtreg ln_LMethionine age sex i.wave##group_num group_num##c.ln_LHomoserine if filter_LMethionine==0 & filter_LHomoserine==0, vce(robust)
			if `l'==0{
				lincom c.ln_LHomoserine
				}
			if `l'==1{
				lincom 1.group_num#c.ln_LHomoserine+c.ln_LHomoserine
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==5{
			quietly xtreg ln_LMethionine age sex i.wave##group_num group_num##c.ln_LHomoserine c.levo_equi##i.wave if filter_LMethionine==0 & filter_LHomoserine==0, vce(robust)
			if `l'==0{
				lincom c.ln_LHomoserine
				}
			if `l'==1{
				lincom 1.group_num#c.ln_LHomoserine+c.ln_LHomoserine
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==1{
			quietly xtreg ln_LMethionine age sex c.ln_LHomoserine##group_num##i.wave if filter_LMethionine==0 & filter_LHomoserine==0, vce(robust)
			if `l'==0{
				lincom c.ln_LHomoserine
				}
			if `l'==1{
				lincom 1.group_num#c.ln_LHomoserine+c.ln_LHomoserine
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `k'==2 | `k'==3{
			quietly xtreg ln_LMethionine age sex c.ln_LHomoserine##group_num##i.wave if filter_LMethionine==0 & filter_LHomoserine==0, vce(robust)
			if `l'==0{
				lincom c.ln_LHomoserine+`k'.wave#c.ln_LHomoserine
			}
			if `l'==1{
				lincom 1.group_num#c.ln_LHomoserine+c.ln_LHomoserine+`k'.wave#c.ln_LHomoserine+`k'.wave#c.ln_LHomoserine#1.group
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `k'==1 | `k'==2 | `k'==3{
			if `l'==0{
				replace ind_w=`k'+0.05 in `i'
				}
			else{
				replace ind_w=`k'-0.05 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		}
	}


graph twoway connected b_coeff ind_w if ind_g==0 & ind_w<3.5, legend(off) xscale(range(0 6)) xlabel(1 "BL" 2 "FU1" 3 "FU2" 4 "Overall" 5 "adj. for med.") ytitle("Regression coefficients of homoserine" "regarding methionine") xtitle("")|| rcap CI_h CI_l ind_w if ind_g==0 & ind_w<3.5|| connected b_coeff ind_w if ind_g==1 & ind_w<3.5|| rcap CI_h CI_l ind_w if ind_g==1 & ind_w<3.5 ||scatter b_coeff ind_w if ind_g==0 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==0 & ind_w>3.5 || scatter b_coeff ind_w if ind_g==1 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==1 & ind_w>3.5,saving(g1.gph, replace)


drop ind_w ind_g b_coeff CI_l CI_h


capture gen ind_w=.
capture gen ind_g=.
capture gen b_coeff=.
capture gen CI_l=.
capture gen CI_h=.


foreach k in 1 2 3 4 5{ 
	foreach l in 0 1{
		if `k'==4{
			quietly xtreg ln_LAlpha age sex i.wave##group_num group_num##c.ln_Cystathion if filter_Cystathion==0 & filter_LAlpha==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathion
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathion+c.ln_Cystathion
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==5{
			quietly xtreg ln_LAlpha age sex i.wave##group_num group_num##c.ln_Cystathion c.levo_equi##i.wave if filter_Cystathion==0 & filter_LAlpha==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathion
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathion+c.ln_Cystathion
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==1{
			quietly xtreg ln_LAlpha age sex c.ln_Cystathionin##group_num##i.wave if filter_Cystathion==0 & filter_LAlpha==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathionin
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathionin+c.ln_Cystathionin
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `k'==2 | `k'==3{
			quietly xtreg ln_LAlpha age sex c.ln_Cystathionin##group_num##i.wave if filter_Cystathion==0 & filter_LAlpha==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin
			}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathionin+c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin#1.group
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `l'==0{
			replace ind_w=`k'+0.05 in `i'
			}
		else{
			replace ind_w=`k'-0.05 in `i'
			}
		replace ind_g=`l' in `i'
		local i=`i'+1
		}
	}

graph twoway connected b_coeff ind_w if ind_g==0 & ind_w<3.5, legend(off) xscale(range(0 6)) xlabel(1 "BL" 2 "FU1" 3 "FU2" 4 "Overall" 5 "adj. for med.") ytitle("Regression coefficients of Cystathionine" "regarding alpha-aminobutyrate") xtitle("")|| rcap CI_h CI_l ind_w if ind_g==0 & ind_w<3.5|| connected b_coeff ind_w if ind_g==1 & ind_w<3.5|| rcap CI_h CI_l ind_w if ind_g==1 & ind_w<3.5 ||scatter b_coeff ind_w if ind_g==0 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==0 & ind_w>3.5 || scatter b_coeff ind_w if ind_g==1 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==1 & ind_w>3.5,saving(g2.gph, replace)


drop ind_w ind_g b_coeff CI_l CI_h


capture gen ind_w=.
capture gen ind_g=.
capture gen b_coeff=.
capture gen CI_l=.
capture gen CI_h=.


local i=1

foreach k in 1 2 3 4 5{ 
	foreach l in 0 1{
		if `k'==4{
			quietly xtreg ln_OA01 age sex i.wave##group_num group_num##c.ln_Cystathion if filter_Cystathion==0 & filter_OA01==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathion
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathion+c.ln_Cystathion
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==5{
			quietly xtreg ln_OA01 age sex i.wave##group_num group_num##c.ln_Cystathion c.levo_equi##i.wave if filter_Cystathion==0 & filter_OA01==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathion
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathion+c.ln_Cystathion
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}		
		if `k'==1{
			quietly xtreg ln_OA01 age sex c.ln_Cystathionin##group_num##i.wave if filter_Cystathion==0 & filter_OA01==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathionin
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathionin+c.ln_Cystathionin
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `k'==2 | `k'==3{
			quietly xtreg ln_OA01 age sex c.ln_Cystathionin##group_num##i.wave if filter_Cystathion==0 & filter_OA01==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin
			}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathionin+c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin#1.group
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `l'==0{
			replace ind_w=`k'+0.05 in `i'
			}
		else{
			replace ind_w=`k'-0.05 in `i'
			}
		replace ind_g=`l' in `i'
		local i=`i'+1
		}
	}

graph twoway connected b_coeff ind_w if ind_g==0 & ind_w<3.5, legend(off) xscale(range(0 6)) xlabel(1 "BL" 2 "FU1" 3 "FU2" 4 "Overall" 5 "adj. for med.") ytitle("Regression coefficients of Cystathionine" "regarding 2-Hydroxybutyrate") xtitle("")|| rcap CI_h CI_l ind_w if ind_g==0 & ind_w<3.5|| connected b_coeff ind_w if ind_g==1 & ind_w<3.5|| rcap CI_h CI_l ind_w if ind_g==1 & ind_w<3.5 ||scatter b_coeff ind_w if ind_g==0 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==0 & ind_w>3.5 || scatter b_coeff ind_w if ind_g==1 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==1 & ind_w>3.5,saving(g3.gph, replace)


drop ind_w ind_g b_coeff CI_l CI_h


capture gen ind_w=.
capture gen ind_g=.
capture gen b_coeff=.
capture gen CI_l=.
capture gen CI_h=.


local i=1
foreach k in 1 2 3 4 5{ 
	foreach l in 0 1{
		if `k'==4{
			quietly xtreg ln_Cystein age sex i.wave##group_num group_num##c.ln_Cystathion if filter_Cystathion==0 & filter_Cystein==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathion
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathion+c.ln_Cystathion
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==5{
			quietly xtreg ln_Cystein age sex i.wave##group_num group_num##c.ln_Cystathion Methoxy_bin c.levo_equi##i.wave if filter_Cystathion==0 & filter_Cystein==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathion
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathion+c.ln_Cystathion
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==1{
			quietly xtreg ln_Cystein age sex c.ln_Cystathionin##group_num##i.wave if filter_Cystathion==0 & filter_Cystein==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathionin
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathionin+c.ln_Cystathionin
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `k'==2 | `k'==3{
			quietly xtreg ln_Cystein age sex c.ln_Cystathionin##group_num##i.wave if filter_Cystathion==0 & filter_Cystein==0, vce(robust)
			if `l'==0{
				lincom c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin
			}
			if `l'==1{
				lincom 1.group_num#c.ln_Cystathionin+c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin+`k'.wave#c.ln_Cystathionin#1.group
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `l'==0{
			replace ind_w=`k'+0.05 in `i'
			}
		else{
			replace ind_w=`k'-0.05 in `i'
			}
		replace ind_g=`l' in `i'
		local i=`i'+1
		}
	}

graph twoway connected b_coeff ind_w if ind_g==0 & ind_w<3.5, legend(off) xscale(range(0 6)) xlabel(1 "BL" 2 "FU1" 3 "FU2" 4 "Overall" 5 "adj. for med.") ytitle("Regression coefficients of Cystathionine" "regarding cystein") xtitle("")|| rcap CI_h CI_l ind_w if ind_g==0 & ind_w<3.5|| connected b_coeff ind_w if ind_g==1 & ind_w<3.5|| rcap CI_h CI_l ind_w if ind_g==1 & ind_w<3.5 ||scatter b_coeff ind_w if ind_g==0 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==0 & ind_w>3.5 || scatter b_coeff ind_w if ind_g==1 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==1 & ind_w>3.5,saving(g4.gph, replace)

graph combine g1.gph g2.gph g3.gph g4.gph, saving(figure_2C_raw.gph, replace)
graph export figure_2C_raw.tif, replace
***********************************************
***Figure 2C***********************************
***********************************************
foreach j of varlist CA-FA_240{
	egen mean_`j'=mean(ln_`j') if group_num==1, by(pr_id)
	gen `j'_cen=ln_`j'-mean_`j'
	egen std_cen=std(`j'_cen)
	replace `j'_cen=. if abs(std_cen)>3
	drop std_cen  mean_`j'
	}

local k=1
foreach j of varlist LHomoserine LMethionine Homocysteine Methoxytyros  LSerine Cystathionin LAlphaaminob OA263Hydroxy Cysteine Taurine Glutathione{
	graph box `j'_cen if filter_`j'==0, over(wave) saving(g`k'.gph, replace)
	local k=`k'+1
	}
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph g6.gph g7.gph g8.gph g9.gph g10.gph g11.gph, saving(Figure_2D_raw.gph, replace)
graph export figure_2D_raw.tif, replace
***********************************************
***Figure 3B***********************************
***********************************************
	
drop *_cen

foreach j of varlist CA-FA_240{
	egen mean_`j'=mean(ln_`j') if wave!=3, by(pr_id)
	egen m2_`j'=mean(ln_`j'), by(group_num)
	gen `j'_cen1=ln_`j'-mean_`j'
	egen std_cen=std(`j'_cen1)
	replace `j'_cen1=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}
	
	
foreach j of varlist CA-FA_240{
	egen mean_`j'=mean(ln_`j') if wave!=1, by(pr_id)
	egen m2_`j'=mean(ln_`j'), by(group_num)
	gen `j'_cen2=ln_`j'-mean_`j'
	egen std_cen=std(`j'_cen2)
	replace `j'_cen2=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}
	
foreach j of varlist levo_equivalent levodopa_dosis_mg levo_diff UPDRS_*{
	egen mean_`j'=mean(`j') if wave!=3, by(pr_id)
	egen m2_`j'=mean(`j'), by(group_num)
	gen `j'_cen1=`j'-mean_`j'
	egen std_cen=std(`j'_cen1)
	replace `j'_cen1=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}
	
	
foreach j of varlist levo_equivalent levodopa_dosis_mg  levo_diff UPDRS_I_sum_-UPDRS_sum_{
	egen mean_`j'=mean(`j') if wave!=1, by(pr_id)
	egen m2_`j'=mean(`j'), by(group_num)
	gen `j'_cen2=`j'-mean_`j'
	egen std_cen=std(`j'_cen2)
	replace `j'_cen2=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}	
	
foreach j of varlist CA-FA_240 levo_equivalent levodopa_dosis_mg levo_diff UPDRS_I_sum_-UPDRS_sum_{
	replace `j'_cen1=`j'_cen2 if wave==3
	replace `j'_cen1=. if wave==1
	}
gen wave2=wave if wave!=1

scatter Cystathionin_cen1 levodopa_dosis_mg_cen1 if wave!=1 & group_num==1 & filter_Cystathionin==0|| lfit Cystathionin_cen1 levodopa_dosis_mg_cen1 if wave!=1 & group_num==1 & filter_Cystathionin==0, range(-250 350) by(wave2) saving(g1.gph, replace)
scatter Methoxytyros_cen1 levodopa_dosis_mg_cen1 if wave!=1 & group_num==1 & filter_Cystathionin==0|| lfit Methoxytyros_cen1 levodopa_dosis_mg_cen1 if wave!=1 & group_num==1 & filter_Cystathionin==0, range(-250 350) by(wave2) saving(g2.gph, replace)

graph twoway (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP006") (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP014", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP016", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP018", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP044", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP050", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP051", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP060", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP062", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP063", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP081", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP092", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP096", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP097", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP100", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP110", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP111", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP112", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP114", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP116", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP119", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP123", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP128", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP131", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP132", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP144", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP145", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP151", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP161", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP162", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (scatter ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1, msymbol(circle_hollow)) (lfit ln_Cystathionin levodopa_dosis_mg if group_num==1 & wave!=1, lwidth(thick)), saving(g3.gph, replace) legend(off)
graph twoway (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP006") (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP014", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP016", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP018", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP044", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP050", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP051", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP060", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP062", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP063", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP081", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP092", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP096", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP097", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP100", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP110", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP111", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP112", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP114", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP116", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP119", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP123", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP128", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP131", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP132", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP144", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP145", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP151", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP161", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (connected ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1 & Probanden_Nr=="DKP162", msymbol(circle_hollow) mcolor(navy) lwidth(thin) lcolor(navy)) (scatter ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1, msymbol(circle_hollow)) (fpfit ln_Methoxytyros levodopa_dosis_mg if group_num==1 & wave!=1, lwidth(thick)), saving(g4.gph, replace) legend(off)

graph combine g1.gph g2.gph g3.gph g4.gph, saving(figure_3B_raw.gph, replace)
graph export figure_3B_raw.tif, replace
***********************************************
***Figure 3C***********************************
***********************************************

scatter FA_180_cen1 levodopa_dosis_mg_cen1 if wave!=1 & group_num==1 & filter_Cystathionin==0|| lfit FA_180_cen1 levodopa_dosis_mg_cen1 if wave!=1 & group_num==1 & filter_FA_180==0, range(-200 250) by(wave2) saving(g1.gph, replace)
scatter TG589_cen1 levo_equivalent_cen1 if wave!=1 & group_num==1 & filter_Cystathionin==0|| lfit TG589_cen1 levo_equivalent_cen1 if wave!=1 & group_num==1 & filter_TG589==0, range(-250 500) by(wave2) saving(g2.gph, replace)
scatter PCO341_cen1 levo_diff_cen1 if wave!=1 & group_num==1 & filter_Cystathionin==0|| lfit PCO341_cen1 levo_diff_cen1 if wave!=1 & group_num==1 & filter_PCO341==0, range(-250 250) by(wave2) saving(g3.gph, replace)

graph combine g2.gph g3.gph, row(1) saving(figure_3C_raw.gph, replace)
graph export figure_3C_raw.tif, replace

**********************************************
***Supplementary Figure S1A*******************
**********************************************

drop *cen1 *cen2

foreach j of varlist CA-FA_240{
	egen mean_`j'=mean(ln_`j'), by(pr_id)
	egen m2_`j'=mean(ln_`j'), by(group_num)
	gen `j'_cen=ln_`j'-mean_`j'+m2_`j'
	egen std_cen=std(`j'_cen)
	replace `j'_cen=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}

local k=1
foreach j of varlist Methoxytyros LMethionine LSerine DL3aminoisob LHomoserine LAlphaaminob LLeucine LAsparagine Cystathionin OA263Hydroxy LPhenylalani FA_140 FA_171 FA_201{
	graph box `j'_cen if filter_`j'==0, over(wave) by(group_num) saving(g`k'.gph, replace)
	local k=`k'+1
	}
graph combine g1.gph g2.gph g3.gph g4.gph g5.gph g6.gph g7.gph g8.gph g9.gph g10.gph g11.gph g12.gph g13.gph g14.gph, saving(figure_S1A_raw.gph, replace)
graph export figure_S1A_raw.tif, replace
**************************************************
***Supplementary Figure S1D/S1B*******************
**************************************************

pca ln_LMethionine ln_LSerine ln_LPhenylala ln_LAsparag ln_LLeucine ln_DL3 ln_LAlphaa ln_OA263 ln_FA_140 ln_FA_171 ln_FA_201 if filter_LPhenylala==0, components(2)
*Table for figure S1B
predict pca1 pca2

foreach j of varlist pca1 pca2{
	egen mean_`j'=mean(`j'), by(pr_id)
	egen m2_`j'=mean(`j'), by(group_num)
	gen `j'_cen=`j'-mean_`j'
	egen std_cen=std(`j'_cen)
	replace `j'_cen=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}

graph box pca1_cen, over(wave) by(group_num) saving(g1.gph, replace)
graph box pca2_cen, over(wave) by(group_num) saving(g2.gph, replace)

graph combine g1.gph g2.gph, saving(figure_S1D_raw.gph, replace)
graph export figure_S1D_raw.tif, replace
 
***********************************************
***Supplementary Figure S2A********************
***********************************************

drop ind_w ind_g b_coeff CI_l CI_h

capture gen ind_w=.
capture gen ind_g=.
capture gen b_coeff=.
capture gen CI_l=.
capture gen CI_h=.


local i=1
foreach k in 1 2 3 4 5{ 
	foreach l in 0 1{
		if `k'==4{
			quietly xtreg ln_SMd181_250 age sex i.wave##group_num group_num##c.ln_SMd181_251 if filter_SMd181_250==0 & filter_SMd181_251==0, vce(robust)
			if `l'==0{
				lincom c.ln_SMd181_251
				}
			if `l'==1{
				lincom 1.group_num#c.ln_SMd181_251+c.ln_SMd181_251
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==5{
			quietly xtreg ln_SMd181_250 age sex i.wave##group_num group_num##c.ln_SMd181_251 c.levo_equi##i.wave if filter_SMd181_250==0 & filter_SMd181_251==0, vce(robust)
			if `l'==0{
				lincom c.ln_SMd181_251
				}
			if `l'==1{
				lincom 1.group_num#c.ln_SMd181_251+c.ln_SMd181_251
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==1{
			quietly xtreg ln_SMd181_250 age sex c.ln_SMd181_251##group_num##i.wave if filter_SMd181_250==0 & filter_SMd181_251==0, vce(robust)
			if `l'==0{
				lincom c.ln_SMd181_251
				}
			if `l'==1{
				lincom 1.group_num#c.ln_SMd181_251+c.ln_SMd181_251
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `k'==2 | `k'==3{
			quietly xtreg ln_SMd181_250 age sex c.ln_SMd181_251##group_num##i.wave if filter_SMd181_250==0 & filter_SMd181_251==0, vce(robust)
			if `l'==0{
				lincom c.ln_SMd181_251+`k'.wave#c.ln_SMd181_251
			}
			if `l'==1{
				lincom 1.group_num#c.ln_SMd181_251+c.ln_SMd181_251+`k'.wave#c.ln_SMd181_251+`k'.wave#c.ln_SMd181_251#1.group
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `l'==0{
			replace ind_w=`k'+0.05 in `i'
			}
		else{
			replace ind_w=`k'-0.05 in `i'
			}
		replace ind_g=`l' in `i'
		local i=`i'+1
		}
	}

graph twoway connected b_coeff ind_w if ind_g==0 & ind_w<3.5, legend(off) xscale(range(0 6)) xlabel(1 "BL" 2 "FU1" 3 "FU2" 4 "Overall" 5 "adj. for med.") ytitle("Regression coefficients of SM(d18:18/25:1)" "regarding SM(d18:18/25:0)") xtitle("")|| rcap CI_h CI_l ind_w if ind_g==0 & ind_w<3.5|| connected b_coeff ind_w if ind_g==1 & ind_w<3.5|| rcap CI_h CI_l ind_w if ind_g==1 & ind_w<3.5 ||scatter b_coeff ind_w if ind_g==0 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==0 & ind_w>3.5 || scatter b_coeff ind_w if ind_g==1 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==1 & ind_w>3.5,saving(g1.gph, replace)

drop ind_w ind_g b_coeff CI_l CI_h

capture gen ind_w=.
capture gen ind_g=.
capture gen b_coeff=.
capture gen CI_l=.
capture gen CI_h=.


local i=1
foreach k in 1 2 3 4 5{ 
	foreach l in 0 1{
		if `k'==4{
			quietly xtreg ln_LPC160 age sex i.wave##group_num group_num##c.ln_Homocitrull if filter_LPC160==0 & filter_Homocitrull==0, vce(robust)
			if `l'==0{
				lincom c.ln_Homocitrull
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Homocitrull+c.ln_Homocitrull
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==5{
			quietly xtreg ln_LPC160 age sex i.wave##group_num group_num##c.ln_Homocitrull c.levo_equi##i.wave if filter_LPC160==0 & filter_Homocitrull==0, vce(robust)
			if `l'==0{
				lincom c.ln_Homocitrull
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Homocitrull+c.ln_Homocitrull
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			if `l'==0{
				replace ind_w=`k'+0.2 in `i'
				}
			else{
				replace ind_w=`k'-0.2 in `i'
				}
			replace ind_g=`l' in `i'
			local i=`i'+1
			}
		if `k'==1{
			quietly xtreg ln_LPC160 age sex c.ln_Homocitrull##group_num##i.wave if filter_LPC160==0 & filter_Homocitrull==0, vce(robust)
			if `l'==0{
				lincom c.ln_Homocitrull
				}
			if `l'==1{
				lincom 1.group_num#c.ln_Homocitrull+c.ln_Homocitrull
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `k'==2 | `k'==3{
			quietly xtreg ln_LPC160 age sex c.ln_Homocitrull##group_num##i.wave if filter_LPC160==0 & filter_Homocitrull==0, vce(robust)
			if `l'==0{
				lincom c.ln_Homocitrull+`k'.wave#c.ln_Homocitrull
			}
			if `l'==1{
				lincom 1.group_num#c.ln_Homocitrull+c.ln_Homocitrull+`k'.wave#c.ln_Homocitrull+`k'.wave#c.ln_Homocitrull#1.group
				}
			replace b_coeff=r(estimate) in `i'
			local se=r(se)
			replace CI_l=r(estimate)-invnormal(0.975)*r(se) in `i' 
			replace CI_h=r(estimate)+invnormal(0.975)*r(se)  in `i'
			}
		if `l'==0{
			replace ind_w=`k'+0.05 in `i'
			}
		else{
			replace ind_w=`k'-0.05 in `i'
			}
		replace ind_g=`l' in `i'
		local i=`i'+1
		}
	}

graph twoway connected b_coeff ind_w if ind_g==0 & ind_w<3.5, legend(off) xscale(range(0 6)) xlabel(1 "BL" 2 "FU1" 3 "FU2" 4 "Overall" 5 "adj. for med.") ytitle("Regression coefficients of homocitrulline" "regarding LPC 16:0") xtitle("")|| rcap CI_h CI_l ind_w if ind_g==0 & ind_w<3.5|| connected b_coeff ind_w if ind_g==1 & ind_w<3.5|| rcap CI_h CI_l ind_w if ind_g==1 & ind_w<3.5 ||scatter b_coeff ind_w if ind_g==0 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==0 & ind_w>3.5 || scatter b_coeff ind_w if ind_g==1 & ind_w>3.5 || rcap CI_h CI_l ind_w if ind_g==1 & ind_w>3.5,saving(g2.gph, replace)

graph combine g1.gph g2.gph, saving(figure_S2A_raw.gph, replace)
graph export figure_S2A_raw.tif, replace

***********************************************
***Supplementary Figure S2B********************
***********************************************

scatter ln_SMd181_250 ln_SMd181_251 if group_num==1 || scatter ln_SMd181_250 ln_SMd181_251 if group_num==0 || lfit ln_SMd181_250 ln_SMd181_251 if group_num==1 || lfit ln_SMd181_250 ln_SMd181_251 if group_num==0, by(wave, row(1)) saving(g1.gph, replace)
scatter ln_LPC160 ln_Homocitrulli if group_num==1 || scatter ln_LPC160 ln_Homocitrulli if group_num==0 || lfit ln_LPC160 ln_Homocitrulli if group_num==1 || lfit ln_LPC160 ln_Homocitrulli if group_num==0, by(wave, row(1)) saving(g2.gph, replace)
graph combine g1.gph g2.gph, row(2) saving(figure_S2B_raw.gph, replace)
graph export figure_S2B_raw.tif, replace

***********************************************
***Supplementary Figure S2C********************
***********************************************
gen prami_azi=pramipexole+azilect
graph box ln_LAlphaa if wave!=1, over(prami_azi) saving(figure_S2C_raw.gph, replace)
graph export figure_S2C_raw.tif, replace

***********************************************
***Supplementary Figure S2D********************
***********************************************
local k=1
foreach j of varlist TG542 SMd181_161 LPAC22_5 PCO343{
	scatter ln_`j' levodopa_dos2 if wave!=1 & group_num==1 & Methoxy_bin==1 ||lfit ln_`j' levodopa_dos2 if wave!=1 & group_num==1 & Methoxy_bin==1, saving(g`k'.gph, replace)
	local k=`k'+1
	}
graph combine g1.gph g2.gph g3.gph g4.gph, saving(figure_S2D_raw.gph, replace)
graph export figure_S2D_raw.tif, replace

***********************************************
***Supplementary Figure S3B********************
***********************************************

foreach j of varlist LHistidine TLCA TDCA TCDCA{
	scatter UPDRS_III_sum ln_`j' if group_num==1 || lfit UPDRS_III_ ln_`j' if group_num==1 , by(wave, row(1))
	graph save  UPDRS_III_`j'.gph, replace
	}
graph combine UPDRS_III_LHistidine.gph UPDRS_III_TCDCA.gph UPDRS_III_TLCA.gph UPDRS_III_TDCA.gph, row(2) saving(figure_S3B_raw.gph, replace)
graph export figure_S3B_raw.tif, replace

clear


********************************************************
***Figures from metagenomic data/simulations************
********************************************************

********************************************
***Figure 3C********************************
********************************************

import delimited A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_bedarf\met_strain_flux.csv
gen group=0 if v3=="control"
replace group=1 if v3=="parkinson"
drop if group==.

foreach j of varlist ex_12dgr180fe-ex_zn2fe{
	gen ln_`j'=ln(`j')
	}
local k=1
foreach j of varlist ex_met_lfe ex_asn_lfe ex_so3fe ex_h2sfe{
	graph box ln_`j', over(group) saving(g`k'.gph, replace)
	local k=`k'+1
	}
graph combine g1.gph g2.gph g3.gph g4.gph, row(2) saving(figure_4C_raw.gph, replace)
graph export figure_4C_raw.tif, replace

********************************************
***Supplementary Figure S4B*****************
********************************************

gen ln_akker=ln(akkermansia_muciniphila_atcc_baa)
local k=1
foreach j of varlist ex_met_lfe ex_asn_lfe ex_h2sfe ex_so3fe {
	quietly mfp reg ln_`j' ln_akker
	predict `j'_pred
	scatter ln_`j' ln_akker || fpfit `j'_pred ln_akker, saving(g`k'.gph, replace) legend(off)
	local k=`k'+1
	}
graph combine g1.gph g2.gph g3.gph g4.gph, row(2) saving(figure_S4B_raw.gph, replace)	
graph export figure_S4B_raw.tif, replace

clear
********************************************
***Figure 3D********************************
********************************************


import excel "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_bedarf\TaxonAbundances.xlsx", firstrow

gen group=0 if B=="control"
replace group=1 if B=="parkinson"
drop if group==.

local k=1
local s=1
foreach j of varlist Actinobacteria-Yokenellaregensburgei{
	local s=`s'+1
	gen _`j'=0 if `j'==0
	replace _`j'=1 if `j'>0
	quietly sum _`j'
	if r(mean)<0.1{
		drop `j'
		local k=`k'+1
		}
	}
display `k'	
display `s'
drop _*

set obs 4000 

foreach j of varlist Actinobacteria-Veillonellaratti{
	gen _`j'=ln(`j')
	egen `j'1=std(`j') if group==1
	egen `j'0=std(`j') if group==0
	drop _`j'
	replace `j'=. if (abs(`j'0)>4 & group==0) |(abs(`j'1)>4 & group==1)
	drop `j'0 `j'1
	}

graph box Bilophila, over(group) saving(figure_4D_raw.gph, replace)
graph export figure_4D_raw.tif, replace



log close

clear
