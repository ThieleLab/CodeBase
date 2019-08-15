//Two-stage MWAS, Medication on pairs of metabolites in PD-cases


clear
clear mata
clear matrix
capture log close
set maxvar 32000
set more off
cd A:\Metabolomics_PD\


use "A:\Metabolomics_PD\denopa_final_data_long.dta", replace

xtset pr_id wave
drop NMCname
set obs 4000

*Outlier filtering	
foreach j of varlist CA-TG5810{
	local l=strlen("`j'")
	if `l'>10{
		local k=substr("`j'",1,12)
		rename `j' `k'
		}
	}
*log-trafo and outlier definition

foreach j of varlist CA-TG5810{
	gen ln_`j'=ln(`j')
	egen std_`j'=std(ln_`j')
	gen filter_`j'=0
	replace filter_`j'=1 if abs(std_`j')>4
	drop std_`j'
	}

cd A:\Metabolomics_PD\Analyses_Denopa\Tables\

gen Methoxy_bin=0 if ln_Methoxytyros<-4
replace Methoxy_bin=1 if ln_Methoxytyros>-4 & ln_Methoxytyros<999
gen ln_crp=ln(labor_crp_mg_dl)
replace levodopa_dosis_mg=0 if Methoxy_bin==0

log using "A:\Metabolomics_PD\Analyses_Denopa\Tables\full_twostage_ratios.log", replace
 
foreach j of varlist CA-TG5810{ 
	local a="`j'"
	foreach k of varlist CA-TG5810{ 
		local b="`k'"
		quietly xtreg ln_`j' c.ln_`k' c.levo_equivalent i.wave age sex group_num##c.ln_`k' if filter_`j'==0 & filter_`k'==0, vce(robust) 
		quietly test c.ln_`k'#1.group_num 
		if r(p)<0.0001 & "`a'"!="`b'"{ 
			display "`j'" 
			display "`k'"
			display r(p)
			}
		} 
	}

log close

log using "A:\Metabolomics_PD\Analyses_Denopa\Tables\full_twostage_ratios_update.log", replace

foreach j of varlist FA_180-FA_240{ 
	local a="`j'"
	foreach k of varlist CA-FA_240{ 
		local b="`k'"
		quietly xtreg ln_`j' c.ln_`k' c.levo_equivalent i.wave age sex group_num##c.ln_`k' if filter_`j'==0 & filter_`k'==0, vce(robust) 
		quietly test c.ln_`k'#1.group_num 
		if r(p)<0.0001 & "`a'"!="`b'"{ 
			display "`j'" 
			display "`k'"
			display r(p)
			}
		} 
	}
log close 

log using "A:\Metabolomics_PD\Analyses_Denopa\Tables\full_twostage_ratios_update_2.log", replace
	
foreach j of varlist CA-FA_240{ 
	local a="`j'"
	foreach k of varlist FA_180-FA_240{ 
		local b="`k'"
		quietly xtreg ln_`j' c.ln_`k' c.levo_equivalent i.wave age sex group_num##c.ln_`k' if filter_`j'==0 & filter_`k'==0, vce(robust) 
		quietly test c.ln_`k'#1.group_num 
		if r(p)<0.0001 & "`a'"!="`b'"{ 
			display "`j'" 
			display "`k'"
			display r(p)
			}
		} 
	}
log close
	
/*
log using "screening_Methoxy_medication_ratios.log", replace

foreach j of varlist CA-Valerylcarni Citrulline-TG5810{
	foreach k of varlist CA-Valerylcarni Citrulline-TG5810{
		quietly xtreg ln_`j' c.ln_`k'##c.levodopa_dosis_mg c.ln_`k'##Methoxy_bin i.wave age sex months_of if filter_`j'==0 & filter_`k'==0 & group_num==1
		quietly test c.ln_`k'#1.Methoxy_bin c.ln_`k'#c.levodopa_dosis_mg 
		if r(p)<0.0001{
			display "`j'"
			display "`k'"
			display r(p)
			}
		}
	}

log using PD_specific_bivariate_trajectories.log, replace
	
foreach j of varlist CA-FA_240{ 
	local a="`j'"
	foreach k of varlist CA-FA_240{ 
		local b="`k'"
		quietly xtreg ln_`j' c.ln_`k'  age sex i.wave##group_num##c.ln_`k' if filter_`j'==0 & filter_`k'==0, vce(robust) 
		quietly test c.ln_`k'#1.group_num 2.wave#1.group_num#c.ln_`k' 3.wave#1.group_num#c.ln_`k'
		if r(p)<0.0001 & "`a'"!="`b'"{ 
			display "`j'" 
			display "`k'"
			display r(p)
			}
		} 
	}	
	
log close



log using "screening_dosis_only_ratios.log", replace

foreach j of varlist CA-Valerylcarni Citrulline-TG5810{
	foreach k of varlist CA-Valerylcarni Citrulline-TG5810{
		quietly xtreg ln_`j' c.ln_`k'##c.levodopa_dosis_mg i.wave age sex months_of if filter_`j'==0 & filter_`k'==0 & group_num==1 & Methoxy_bin==1
		quietly test c.ln_`k'#c.levodopa_dosis_mg 
		if r(p)<0.0001{
			display "`j'"
			display "`k'"
			display r(p)
			}
		}
	}


log using "screening_dosis_ratios.log", replace

foreach j of varlist CA-Valerylcarni Citrulline-TG5810{
	foreach k of varlist CA-Valerylcarni Citrulline-TG5810{
		quietly xtreg ln_`j' c.ln_`k'##c.levodopa_dosis_mg i.wave age sex months_of if filter_`j'==0 & filter_`k'==0 & group_num==1 
		quietly test c.ln_`k'#c.levodopa_dosis_mg 
		if r(p)<0.0001{
			display "`j'"
			display "`k'"
			display r(p)
			}
		}
	}
		
log close*/


clear
