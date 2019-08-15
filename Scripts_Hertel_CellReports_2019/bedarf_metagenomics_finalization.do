*Analyses fluxes and metagenomics from Bedarf et al. 2017

clear
clear mata
clear matrix
set more off
capture log close
set maxvar 32000
import delimited A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_bedarf\Bedarf_FluxSpans.csv

gen _varname=substr(v1,1,8)
gen _varname2=substr(v1, -16,.)
replace _varname=_varname+"_"+_varname2
drop  v1 strain reaction _varname2
xpose, clear varname
gen id="ERS"+substr(_varname,4,.)
save A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_bedarf\Bedarf_FluxSpans.dta,replace
clear
import delimited A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_bedarf\met_strain_flux.csv

rename v1 id 
merge 1:1 id using A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_bedarf\Bedarf_FluxSpans.dta

*Histograms for secretions
local t=1
local s=1
foreach j of varlist ex_12dgr180fe-ex_zn2fe{
	gen `j'_bin=0 if `j'==0
	replace `j'_bin=1 if `j'>0 & `j'<999999999
	quietly sum `j'_bin
	gen `j'_det=r(mean)
	if r(mean)>0.5{
		local t=`t'+1
		}
	local s=`s'+1
	}
display `t'
display `s'
drop *_bin 

cd A:\Hertel\Finalized\Hertel_2018_DeNoPa\figure\descriptive_bedarf\histograms
foreach j of varlist ex_12dgr180fe-ex_zn2fe{
	gen ln_`j'=ln(`j')
	histogram ln_`j'
	graph save ln_`j'.gph, replace
	graph export ln_`j'.png, replace
	histogram `j'
	graph save z_`j'.gph, replace
	graph export z_`j'.png, replace
	}

*Simple tests

gen group=0 if v3=="control"
replace group=1 if v3=="parkinson"

*box plots

cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\figure\descriptive_bedarf\boxplots"
foreach j of varlist ex_12dgr180fe-ex_zn2fe{
	quietly sum `j'_det
	if r(mean)>0.5{
		graph box ln_`j', by(group)
		graph save box_ln_`j'.gph, replace
		graph export box_ln_`j'.png, replace
		}
	}


*skewness after ln

foreach j of varlist ex_12dgr180fe-ex_zn2fe{
	egen std_ln_`j'=std(ln_`j')
	}

foreach j of varlist ex_12dgr180fe-ex_zn2fe{
	quietly sum `j'_det
	if r(mean)>0.5{
		quietly sktest ln_`j'
		if r(P_chi2)<0.05{
			display "`j'"
			display r(P_chi2)
			}
		}
	}

set obs 4000 

gen var_=""
gen p_val=.
gen CI_l=.
gen CI_h=.
gen b_coeff=.
gen CI_95=""
local i=1
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\tables"
foreach j of varlist ex_12dgr180fe-ex_zn2fe{
	quietly sum `j'_det
	if r(mean)>0.3{
		replace var_="`j'" in `i'
		quietly reg ln_`j' ln_ex_so3 group if abs(std_ln_`j')<4, vce(robust)
		matrix V=e(V)
		local indV=colsof(V)-1
		matrix B=e(b)
		local indB=colsof(B)-1
		replace b_coeff=B[1,`indB'] in `i'
		replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
		replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV']) in `i'
		tostring(CI_l), replace force
		tostring(CI_h), replace force
		tostring(b_coeff),replace force 
		replace CI_95=substr(b_coeff, 1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
		if B[1,`indB']>0{
			replace p_val=2*normal(-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
			}
		else{
			replace p_val=2*normal(B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
			}
		destring(CI_l), replace force
		destring(CI_h), replace force
		destring(b_coeff), replace force
		local i=`i'+1
		}
	}

export excel var_ CI_95 p_val using fluxes_bedarf_PD_cond_Akker.xlsx, replace firstrow(variables)

replace var_=""
replace p_val=.
replace CI_l=.
replace CI_h=.
replace b_coeff=.
replace CI_95=""

local i=1
	
*variance contribution of b.wadsworthia/a.muciniphila
capture log close
log using A:\Hertel\Finalized\Hertel_2018_DeNoPa\tables\variance_contribution.log, replace
egen std_bilo=std(bilophila)
foreach j of varlist ex_met_lfe ex_asn_lfe ex_h2sfe ex_so3fe {
	mfp reg `j' akker 
	}
foreach j of varlist ex_met_lfe ex_asn_lfe ex_h2sfe ex_so3fe {
	mfp reg `j' bilo if abs(std_bilo<4)
	}
clear
log close

import excel "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\raw_bedarf\Abundances_Studies_combined_update.xlsx", sheet("Bedarf_paper_abundances") firstrow


cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\tables"

gen group=0 if B=="control"
replace group=1 if B=="parkinson"
drop if group==.

foreach j of varlist DGR180ti-sink_s{
	gen `j'_bin0=0 if `j'==0
	replace `j'_bin0=1 if `j'>0
	gen `j'_bin1=0 if `j'==1
	replace `j'_bin1=1 if `j'<1
	quietly sum `j'_bin0
	if r(mean)<0.1{
		drop `j'
		}
	quietly sum `j'_bin1
	if r(mean)<0.1{
		capture drop `j'
		}
	drop `j'_bin0 `j'_bin1
	}
foreach j of varlist DGR180ti-sink_s{
	quietly sum `j'
	if r(mean)>0.99{
		capture drop `j'
		}
	if r(mean)<0.01{
		capture drop `j'
		}
	}
set obs 4000 

gen var_=""
gen p_val=.
gen CI_l=.
gen CI_h=.
gen OR=.
gen CI_95=""
local i=1

foreach j of varlist DGR180ti-sink_s{
	local lab: variable label `j'
	replace var_="`lab'" in `i'
	fracreg logit `j' group
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace OR=exp(B[1,`indB']) in `i'
	replace CI_l=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
	replace CI_h=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(OR), replace force
	replace CI_95=substr(OR, 1,5)+"("+substr(CI_l,1,5)+","+substr(CI_h,1,5)+")" in `i'
	if B[1,`indB']>0{
		replace p_val=2*normal(-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
		}
	else{
		replace p_val=2*normal(B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
		}
	destring(CI_l), replace force
	destring(CI_h), replace force
	destring(OR), replace force
	local i=`i'+1
	}

export excel var_ CI_95 p_val using reaction_abundancies_PD.xlsx, replace firstrow(variables)

replace var_=""
replace p_val=.
replace CI_l=.
replace CI_h=.
replace OR=.
replace CI_95=""

local i=1

clear

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


gen var_=""
gen p_val=.
gen CI_l=.
gen CI_h=.
gen OR=.
gen CI_95=""
local i=1

foreach j of varlist Actinobacteria-Veillonellaratti{
	local lab: variable label `j'
	replace var_="`lab'" in `i'
	quietly fracreg logit `j' group
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace OR=exp(B[1,`indB']) in `i'
	replace CI_l=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
	replace CI_h=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(OR), replace force
	replace CI_95=substr(OR, 1,5)+"("+substr(CI_l,1,5)+","+substr(CI_h,1,5)+")" in `i'
	if B[1,`indB']>0{
		replace p_val=2*normal(-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
		}
	else{
		replace p_val=2*normal(B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
		}
	quietly fracreg logit `j' group
	quietly test group
	destring(CI_l), replace force
	destring(CI_h), replace force
	destring(OR), replace force
	local i=`i'+1
	}
cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\tables\"

export excel var_ CI_95 p_val using taxon_abundancies_PD.xlsx, replace firstrow(variables)

replace var_=""
replace p_val=.
replace CI_l=.
replace CI_h=.
replace OR=.
replace CI_95=""

local i=1
