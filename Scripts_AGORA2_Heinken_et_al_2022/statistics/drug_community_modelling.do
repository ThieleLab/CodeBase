*Analyses CRC simulated drug metabolites

clear
clear mata
clear matrix
set more off
set maxvar 32000
capture log close 


cd A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision


import delimited A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\AGORA2_CRC_Objectives_JD.txt, varnames(1) clear 
drop in 1
drop if substr(objective,-2,.)=="_1"
gen model_id=substr(objective,23,.)
replace model_id=subinstr(model_id, "e", "e_",1)
merge 1:1 model_id using meta_data.dta
drop if _merge==2
drop _merge
destring(ex_sn38fe-ex_cholatefe), replace
rename ex_dh5furafe ex_dh5furafe_fcsn
rename v6 ex_dh5furafe_5fura
rename ex_ac5asafe ex_ac5asafe_5asa
rename v11 ex_ac5asafe_bzd

gen pattern_EX=""
foreach j of varlist ex_*{
	gen b_`j'=`j'
	replace b_`j'=1 if `j'>0.00001 & `j'!=.
	replace b_`j'=0 if `j'<0.00001
	quietly sum b_`j'
	if r(mean)>0 & r(mean)<1{
		tostring(b_`j'), gen(s_`j')
		replace pattern_EX=pattern_EX+s_`j'
		}
	}

egen gr_pattern_EX=group(pattern_EX)
tab gr_pattern_EX, gen(pattern)
egen sum_secreted=rowtotal(b_ex*), missing
drop if sum_secreted==0 | sum_secreted==.

mkspline spline_age=age, cubic nknots(4)
mkspline spline_bmi=bmi, cubic nknots(4)

cd A:\AGORA_2_New\Files_for_Johannes_revision\results\logs
log using community_modelling_drugs.log, replace

*pat_stat 
foreach j of varlist ex_*{
	quietly sum b_`j'
	if r(mean)>0.05 & r(mean)<0.95{
		logit b_`j' spline_age* sex bmi pat_stat, or
		test pat_stat
		test bmi
		}
	quietly sum b_`j'
	if r(mean)>0.5{
		reg `j' spline_age* sex bmi pat_stat, vce(robust)
		test pat_stat
		test bmi
		}
	}
*age
foreach j of varlist ex_*{
	quietly sum b_`j'
	if r(mean)>0.05 & r(mean)<0.95{
		logit b_`j' spline_age* sex pat_stat, or
		test spline_age1 spline_age2 spline_age3
		display r(p)
		test spline_age2 spline_age3
		display r(p)
		}
	quietly sum `j'
	if r(mean)>0.5{
		reg `j' spline_age* sex pat_stat, vce(robust)
		test spline_age1 spline_age2 spline_age3
		display r(p)
		test spline_age2 spline_age3
		display r(p)
		}
	}

	
*Descriptive Table Drug metabolites

cd A:\AGORA_2_New\Files_for_Johannes_revision\results\logs
*log using log_analyses_CRC.log, replace
set obs 1000
local i=1
foreach k in all CRC control{
	gen mean_`k'=.
	gen sd_`k'=.
	gen det_rate_`k'=.
	}
gen _var=""
gen p_detect=.
gen p_mean_diff=.
gen p_val=.
gen CI_l=.
gen CI_h=.
gen b_coeff=.
gen CI_95=""
gen p_val_age=.
gen p_val_sex=.
gen p_val_bmi=.
gen p_val_age_nl=.
gen p_val_bmi_nl=.

foreach j of varlist ex_*fe*{
	replace _var="`j'" in `i'
	foreach k in all CRC control{
		if "`k'"=="all"{
			quietly sum `j'
			local n=r(N)
			replace mean_`k'=r(mean) in `i'
			replace sd_`k'=r(sd) in `i'
			quietly sum `j' if `j'>0
			local n2=r(N)
			local drate=`n2'/`n'
			replace det_rate_`k'=`n2'/`n' in `i'
			}
		if "`k'"=="CRC"{
			quietly sum `j' if pat_stat==1
			local n=r(N)
			replace mean_`k'=r(mean) in `i'
			replace sd_`k'=r(sd) in `i'
			quietly sum `j' if `j'>0 & pat_stat==1
			local n2=r(N)
			replace det_rate_`k'=`n2'/`n' in `i'
			}
		if "`k'"=="control"{
			quietly sum `j' if pat_stat==2
			local n=r(N)
			replace mean_`k'=r(mean) in `i'
			replace sd_`k'=r(sd) in `i'
			quietly sum `j' if `j'>0 & pat_stat==2
			local n2=r(N)
			replace det_rate_`k'=`n2'/`n' in `i'
			}
		}
	gen `j'_bin=`j'
	replace `j'_bin=1 if `j'>0 & `j'!=.
	quietly sum `j'_bin
	if r(mean)<0.5 & r(mean)>0.05{
		logit `j'_bin spline_age* sex spline_bmi* pat_stat, or
		quietly test pat_stat
		replace p_detect=r(p) in `i'
		quietly test spline_bmi1 spline_bmi2 spline_bmi3
		replace p_val_bmi=r(p) in `i'
		quietly test spline_bmi2 spline_bmi3
		replace p_val_bmi_nl=r(p) in `i'
		logit `j'_bin spline_age* sex pat_stat, or
		quietly test spline_age1 spline_age2 spline_age3
		replace p_val_age=r(p) in `i'
		quietly test spline_age2 spline_age3
		replace p_val_age_nl=r(p) in `i'
		quietly test sex
		replace p_val_sex=r(p) in `i'
		}
	if `drate'>0.5{
		quietly ttest `j', by(pat_stat) welch
		replace p_mean_diff=r(p) in `i'
		reg `j' spline_age* sex spline_bmi* pat_stat, vce(robust)
		local df=e(df_r)
		matrix V=e(V)
		local indV=colsof(V)-1
		matrix B=e(b)
		local indB=colsof(B)-1
		replace b_coeff=B[1,`indB'] in `i'
		replace CI_l=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
		replace CI_h=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
		tostring(CI_l), replace force
		tostring(CI_h), replace force
		tostring(b_coeff),replace force 
		replace CI_95=substr(b_coeff, 1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
		if B[1,`indB']>0{
			replace p_val=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
			}
		else{
			replace p_val=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
			}
		destring(CI_l), replace force
		destring(CI_h), replace force
		destring(b_coeff), replace force
		reg `j' spline_age* sex spline_bmi* pat_stat, vce(robust)
		quietly test spline_bmi1 spline_bmi2 spline_bmi3
		replace p_val_bmi=r(p) in `i'
		quietly test spline_bmi2 spline_bmi3
		replace p_val_bmi_nl=r(p) in `i'
		reg `j' spline_age* sex pat_stat, vce(robust)
		quietly test spline_age1 spline_age2 spline_age3
		replace p_val_age=r(p) in `i'
		quietly test spline_age2 spline_age3
		replace p_val_age_nl=r(p) in `i'
		quietly test sex
		replace p_val_sex=r(p) in `i'
		}
	local i=`i'+1
	drop `j'_bin
	}

export excel _var mean_all sd_all det_rate_all mean_CRC sd_CRC det_rate_CRC mean_control sd_control det_rate_control p_detect b_coeff CI_95 p_val p_val_age p_val_age_nl p_val_sex p_val_bmi p_val_bmi_nl using descr_table_drug_metabolites.xlsx, replace firstrow(variables)

local i=1
foreach k in all CRC control{
	replace mean_`k'=.
	replace sd_`k'=.
	replace det_rate_`k'=.
	}
replace _var=""
replace p_detect=.
replace p_mean_diff=.
replace p_val=.
replace CI_l=.
replace CI_h=.
replace b_coeff=.
replace CI_95=""
replace p_val_age=.
replace p_val_sex=.
replace p_val_bmi=.
replace p_val_age_nl=.
replace p_val_bmi_nl=.

*Graphics


*age
cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
foreach j of varlist ex_ac5asafe_bzd ex_5furafe ex_ac5asafe_5asa ex_bvufe{
	quietly reg `j' spline_age* pat_stat
	predict xb_`j'
	twoway (scatter `j' age if pat_stat==1, mcolor(navy) msymbol(circle_hollow)) (scatter `j' age if pat_stat==2, mcolor(maroon) msymbol(circle_hollow)) (mband xb_`j' age if pat_stat==1, lcolor(dknavy) lwidth(thick)) (mband xb_`j' age if pat_stat==2, lcolor(cranberry) lwidth(thick)), ytitle(Secretion potential) ylabel(, nogrid) xtitle(Age in years) title("`j'") graphregion(fcolor(white) lcolor(white)) saving(`j'_age, replace)
	}

graph combine ex_ac5asafe_bzd_age.gph ex_5furafe_age.gph ex_ac5asafe_5asa_age.gph ex_bvufe_age.gph, graphregion(fcolor(white) lcolor(white)) xsize(5) ysize(5) saving(Figure_5B_raw, replace)
 drop xb_*

*BMI

foreach j of varlist ex_r406fe{
	quietly reg `j' spline_bmi1 pat_stat
	predict xb_`j'
	twoway (scatter `j' bmi if pat_stat==1, mcolor(navy) msymbol(circle_hollow)) (scatter `j' bmi if pat_stat==2, mcolor(maroon) msymbol(circle_hollow)) (mband xb_`j' bmi if pat_stat==1, lcolor(dknavy) lwidth(thick)) (mband xb_`j' bmi if pat_stat==2, lcolor(cranberry) lwidth(thick)), ytitle(Secretion potential) ylabel(, nogrid) xtitle(BMI [kg/m^2]) title("`j'") graphregion(fcolor(white) lcolor(white)) saving(`j'_bmi, replace) xsize(5) ysize(5)
	}
	
*Sex

foreach j of varlist ex_cholatefe ex_dfdurife{
	graph box `j', over(sex) graphregion(fcolor(white) lcolor(white)) ylabel(,nogrid) ytitle("Secretion potential") title("`j'") saving(`j'_sex, replace)
	}
	
graph combine ex_cholatefe_sex.gph ex_dfdurife_sex.gph, graphregion(fcolor(white) lcolor(white)) xsize(10) ysize(5) saving(sex_differences_raw, replace) row(1)
graph combine ex_r406fe_bmi.gph sex_differences_raw.gph, graphregion(fcolor(white) lcolor(white)) xsize(10) ysize(5) saving(Figure_S9_raw, replace) row(1)

*ex_nchlphnclfe
*graph box ex_nchlphnclfe, over(pat_stat) graphregion(fcolor(white) lcolor(white)) ylabel(,nogrid) ytitle("Secretion potential") title("`j'") saving(`j'_sex, replace)
	
log close 

 




	