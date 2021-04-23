* IBD Sulfur species diversity script

clear
clear mata
clear matrix
set more off
set maxvar 32000
capture log close

cd A:\IBD_project\IBD_simulation_data

log using manuscript_analyses.log, replace

use data_merge_all.dta
merge 1:1 id using reactions.dta
drop _merge

tab stratification, gen(group_)

*General classification

foreach j of varlist __12dhchol-__xtsn{
	gen bin`j'=0 if `j'!=.
	replace bin`j'=1 if `j'!=. & `j'>0
	}

gen class_sec=""
foreach j of varlist bin*{
	tostring(`j'), gen(s_`j')
	replace class_sec=class_sec+s_`j'
	drop s_`j'
	}
tab class_sec
egen gr_class_sec=group(class_sec)

egen total_number_secr=rowtotal(bin__*)

gen hspg_cspg_bin=bin__hspg_degr_1
replace hspg_cspg_bin=1 if bin__cspg_a_degr==1

*Classification Sulfur metabolites


local list="gam26s Lcystin thm hspg_degr_1 cspg_a_degr coa ch4s acgalidour2s tsul so3 so4 gthrd gthox ahcys 5mta h2s met_L cys_L taur hspg_degr_10 hspg_degr_11 hspg_degr_12 hspg_degr_13 hspg_degr_1 hspg_degr_2 hspg_degr_3 hspg_degr_4 hspg_degr_5 hspg_degr_6 hspg_degr_7 hspg_degr_8 hspg_degr_9 cspg_a_degr cspg_ab_rest cspg_b_degr cspg_c_degr cspg_c_rest" 

gen class_sulfur_sec=""
gen total_number_sulfur=0
foreach k in `list'{
	replace total_number_sulfur=total_number_sulfur+bin__`k'
	tostring(bin__`k'), gen(s_EX_`k'_e__bin)
	replace class_sulfur_sec=class_sulfur_sec+s_EX_`k'_e__bin
	drop s_EX_`k'_e__bin
	}
tab class_sulfur_sec
egen gr_class_sulfur_sec=group(class_sulfur_sec)

gen total_number_minusS=total_number_secr-total_number_sulfur
graph box total_number_sulfur, over(stratification) saving(temp1, replace)
graph box total_number_secr, over(stratification)  saving(temp2, replace)
graph box total_number_minusS, over(stratification)  saving(temp3, replace)

graph combine temp1.gph temp2.gph temp3.gph, row(1) 

graph box total_number_sulfur, over(bin__cspg_a_degr) saving(temp1, replace)
graph box total_number_secr, over(bin__cspg_a_degr)  saving(temp2, replace)
graph box total_number_minusS, over(bin__cspg_a_degr)  saving(temp3, replace)

graph combine temp1.gph temp2.gph temp3.gph, row(1)

gen total_number_strains=0
foreach j of varlist _ATCC_49176-_ATCC_43003{
	gen bin`j'=0 if `j'==0
	replace bin`j'=1 if `j'>0 & `j'<1.5
	replace total_number_strains=total_number_strains+bin`j'
	}

scatter total_number_secr total_number_strains
scatter total_number_sulfur total_number_strains

gen ln_total_number_strains=ln(total_number_strains)

reg total_number_secr ln_total_number_strains 2.group
	
local list_ind="thm hspg_degr_1 cspg_a_degr coa ch4s tsul so3 so4 gthrd ahcys 5mta met_L"

gen class_sulfur_sec_ind=""
gen total_number_sulfur_ind=0
foreach k in `list_ind'{
	replace total_number_sulfur_ind=total_number_sulfur_ind+bin__`k'
	tostring(bin__`k'), gen(s_EX_`k'_e__bin)
	replace class_sulfur_sec_ind=class_sulfur_sec_ind+s_EX_`k'_e__bin
	drop s_EX_`k'_e__bin
	}
tab class_sulfur_sec
egen gr_class_sulfur_sec_ind=group(class_sulfur_sec)

local list_end="tsul h2s so3 so4 ch4s"


gen class_sulfur_sec_end=""
gen total_number_sulfur_end=0
foreach k in `list_end'{
	replace total_number_sulfur_end=total_number_sulfur_end+bin__`k'
	tostring(bin__`k'), gen(s_EX_`k'_e__bin)
	replace class_sulfur_sec_end=class_sulfur_sec_end+s_EX_`k'_e__bin
	drop s_EX_`k'_e__bin
	}
tab class_sulfur_sec_end
egen gr_class_sulfur_sec_end=group(class_sulfur_sec_end)

*Table: Screen secretions on dependency on group-status

set obs 500
gen name=""
gen perc_1=.
gen perc_2=.
gen perc_3=.
gen p_value=.
local i=1

foreach j of varlist __12dhchol-__xtsn{
	local nam=substr("`j'", 3,.)
	replace name="`nam'" in `i'
	quietly tab bin`j' group, ex
	replace p_value=r(p_exact) in `i'
	foreach k in 1 2 3{
		quietly sum bin`j' if group==`k'
		replace perc_`k'=r(mean) in `i'
		}
	local i=`i'+1
	}

export excel name perc_1 perc_2 perc_3 p_value using "secretions_binary_group.xlsx", replace firstrow(variables)

local i=1
replace name=""
replace perc_1=.
replace perc_2=.
replace perc_3=.
replace p_value=.
	

local list_all="__phenol __dtmp __tsul __acmana __pydxn __dttp __ha_deg1 __ha_pre1 __chor __ncam __hspg_degr_1 __cspg_c_degr __Ser_Thr __acgal __idon_L __isobut __12dhchol __btd_RR __pac __pydx __uchol __so4 __gthrd __so3 __13ppd __icdchol __ch4s __phpyr __g6p __3dhcdchol __26dap_M __coprost __2mbut __ctbt __gbbtn __12ppd_S __oaa __2hyoxplac __glyclt __4hphac __4hbz __ahcys __xtsn __csn __xan __HC02191 __4hoxpacd __5mthf __thf __fol __fald __n2o __n2 __oxa __actn_R __ind3ac __dopa __hista __coa __diact __5mta __ch4 __adocbl __M03134 __but __34dhpha __mnl __pyr __dgsn __4mcat __h2o2 __ispre __dhpppn __tmao __ura __Lkynr __glutar __bglc __btoh __thm __15dap __2ddglcn __3mop __glcn __glyb __pi __rib_D __tma __met_L"
	
gen total_number_ind=0
foreach k in `list_all'{
	replace total_number_ind=total_number_ind+bin`k'
	}

reg total_number_sulfur_ind group_2 if group_3==0
bootstrap r(ind_eff) r(dir_eff), reps(1000): sgmediation total_number_sulfur_ind, mv(hspg_csp) iv(group_2)

tab class_sulfur_sec_end, gen(dummy)
foreach j of varlist dummy*{
	tab `j' group, ex
	}

*Mediation analyses
bootstrap r(ind_eff) r(dir_eff), reps(1000): sgmediation total_number_sulfur_ind if group_3==0, mv(hspg_csp) iv(group_2)
	
	
*Figures

*A

reg total_number_sec ln_total_number_strains 2.group
matrix B=e(b)
scatter total_number_sec total_number_strains if group==1 ||scatter total_number_sec total_number_strains if group==2 ||scatter total_number_sec total_number_strains if group==3 || function y=B[1,1]*ln(x)+B[1,3], range(15 300) || function y=B[1,1]*ln(x)+B[1,3]+B[1,2], range(15 300) saving(fig_A_raw, replace)

*B
graph box total_number_secr, over(stratification) saving(fig_B_raw, replace)

*C
graph box total_number_sulfur_ind, over(stratification) by(hspg_csp) saving(fig_C_raw, replace)

*D	

tab class_sulfur_sec_end, gen(sulfur_pattern)
foreach j of varlist *pattern*{
	tab `j' stratification, ex
	display r(p_exact)
	}

clear 
log close
