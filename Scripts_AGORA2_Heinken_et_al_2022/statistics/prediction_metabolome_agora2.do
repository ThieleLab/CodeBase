*data preparation AGROA2 Yachida integration for in silico in vivo pattern analyses

clear
clear mata
clear matrix 
set more off
capture log close
set maxvar 32000

cd A:\AGORA_2_New\Files_for_Johannes_revision\processed_data
use "CRC_AGORA2_merged.dta"

local n=0
local list_met=""
foreach j of varlist _C00024-_C02242{
	gen b_`j'=`j'
	replace b_`j'=1 if b_`j'>0 & b_`j'!=.
	quietly sum b_`j'
	if r(mean)>0.5{
		local list_met="`list_met'"+" "+"`j'"
		local n=`n'+1
		}
	}
display `n'
display "`list_met'"

local n=0
local list_sec=""
foreach j of varlist EX_12dgr180_fe_-EX_zn2_fe_{
	gen b_`j'=`j'
	replace b_`j'=1 if b_`j'>0 & b_`j'!=.
	quietly sum b_`j'
	if r(mean)>0.5{
		local list_sec="`list_sec'"+" "+"`j'"
		local n=`n'+1
		}
	}
display `n'
display "`list_sec'"


clear
/*
import delimited "A:\AGORA_2_New\Final_data_for_analysis\Revision\VMHIds.csv", varnames(1)
replace keggid="_"+keggid
forvalues k=1(1)517{
	local i=`k'
	local s=keggid in `i'
	local t=0
	foreach j in `list_met'{
		if "`s'"=="`j'"{
			local t=1
			}
		}
	if `t'==0{
		replace keggid="" in `i'
		}
	}
keep if keggid!=""

*filling in missing vmh ids manually

drop if vmhid==""
save "A:\AGORA_2_New\Final_data_for_analysis\Revision\list_AGORA2_crc_metabolome.dta", replace

replace vmhid="EX_"+vmhid+"_fe_"
forvalues k=1(1)136{
	local i=`k'
	local s=vmhid in `i'
	local t=0
	foreach j in `list_sec'{
		if "`s'"=="`j'"{
			local t=1
			}
		}
	if `t'==0{
		replace vmhid="" in `i'
		}
	}
keep if vmhid!=""
save "A:\AGORA_2_New\Final_data_for_analysis\Revision\list_AGORA2_crc_metabolome_matched.dta", replace
gen vmhid2=strreverse(substr(strreverse(substr(vmhid,4,.)), 5,.))
drop vmhid
rename vmhid2 vmhid
clear
*/
use "CRC_AGORA2_merged.dta"


rename _C00147 _ade
rename _C00212 _adn
rename _C00041 _ala_L
rename _C00062 _arg_L
rename _C00152 _asn_L
rename _C00049 _asp_L
rename _C00719 _glyb
rename _C00246 _but
rename _C01672 _15dap
rename _C00695 _cholate
rename _C00475 _cytd
rename _C02679 _ddca
rename _C00364 _dtmp
rename _C00334 _4abut
rename _C00064 _gln_L
rename _C00025 _glu_L
rename _C00329 _gam
rename _C00037 _gly
rename _C00093 _glyc3p
rename _C00387 _gsn
rename _C00135 _his_L
rename _C00262 _hxan
rename _C00407 _ile_L
rename _C00294 _ins
rename _C00186 _lac_L
rename _C00123 _leu_L
rename _C00047 _lys_L
rename _C00073 _met_L
rename _C00140 _acgam
rename _C00270 _acnam
rename _C00153 _ncam
rename _C00077 _orn
rename _C00864 _pnto_R
rename _C00079 _phe_L
rename _C00148 _pro_L
rename _C00163 _ppa
rename _C00134 _ptrc
rename _C00255 _ribflv
rename _C00065 _ser_L
rename _C00315 _spmd
rename _C00042 _succ
rename _C00245 _taur
rename _C00378 _thm
rename _C00188 _thr_L
rename _C00214 _thymd
rename _C00423 _cinnm
rename _C00078 _trp_L
rename _C00082 _tyr_L
rename _C00483 _tym
rename _C00086 _urea
rename _C00299 _uri
rename _C00183 _val_L
rename _C08262 _isoval
rename _C00170 _5mta

local list_compare="ade adn ala_L arg_L asn_L asp_L glyb but 15dap cytd ddca dtmp 4abut gln_L glu_L gam gly glyc3p gsn his_L hxan ile_L ins lac_L leu_L lys_L met_L acgam acnam ncam orn pnto_R phe_L pro_L ppa ptrc ribflv ser_L spmd succ thm thr_L thymd cinnm trp_L tyr_L tym urea uri val_L isoval 5mta"


drop if _but==.

gen pattern_EX=""
foreach j of varlist EX_12dgr180_fe_-EX_zn2_fe_{
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
tab gr_pattern_EX
egen sum_secreted=rowtotal(b_EX*), missing
drop if sum_secreted==0 | sum_secreted==.

foreach j of varlist Abiotrophia_defectiva-Weissella_paramesenteroides{
	gen b_`j'=`j'
	replace b_`j'=1 if `j'>0 & `j'!=.
	}
egen sum_species=rowtotal(b_Abiotrophia_defectiva-b_Weissella_paramesenteroides)


******************
*pattern analyses*
******************

foreach j in `list_compare'{
	gen ln_`j'=ln(_`j')
	gen b_`j'=0 if _`j'==0
	replace b_`j'=1 if _`j'!=. & _`j'>0
	}
	
cd "A:\AGORA_2_New\Files_for_Johannes_revision\results\Invivo_insilico_tables"

set obs 1000
local i=1
gen _var=""
foreach j in conc flux{
	gen p_val_detect_`j'=.
	gen p_val_global_`j'=.
	gen p_val_linear_`j'=.
	gen r2_`j'=.
	gen CI_l_`j'=.
	gen CI_h_`j'=.
	gen CI_95_`j'=""
	gen b_coeff_`j'=.
	gen CI_l_`j'_detect=.
	gen CI_h_`j'_detect=.
	gen b_coeff_`j'_detect=.
	gen CI_95_`j'_detect=""
	gen d_rate_`j'=.
	}


foreach j of varlist Abiotrophia_defectiva-Weissella_paramesenteroides{
	quietly sum b_`j'
	if r(mean)>0.1{
		local file="`j'"+".xlsx"	
		foreach k in `list_compare'{
			replace _var="`k'" in `i'
			*Fluxes linear
			quietly sum b_`j'
			if r(mean)>0.5{
				quietly reg EX_`k'_fe_ age sex bmi pat_stat `j', vce(robust)
				local df=e(df_r)
				matrix V=e(V)
				local indV=colsof(V)-1
				matrix B=e(b)
				local indB=colsof(B)-1
				replace b_coeff_flux=B[1,`indB'] in `i'
				replace CI_l_flux=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
				replace CI_h_flux=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
				tostring(CI_l_flux), replace force
				tostring(CI_h_flux), replace force
				tostring(b_coeff_flux),replace force 
				replace CI_95_flux=substr(b_coeff_flux, 1,6)+"("+substr(CI_l_flux,1,6)+","+substr(CI_h_flux,1,6)+")" in `i'
				destring(CI_l_flux), replace force
				destring(CI_h_flux), replace force
				destring(b_coeff_flux), replace force
				quietly reg EX_`k'_fe_ age sex bmi pat_stat `j', vce(robust)
				quietly test `j'
				replace p_val_linear_flux=r(p) in `i'
				}
			*Fluxes binary 
			quietly sum b_`j'
			if r(mean)>0.1 & r(mean)<0.9{
				quietly reg EX_`k'_fe_ age sex bmi pat_stat b_`j', vce(robust)
				local df=e(df_r)
				matrix V=e(V)
				local indV=colsof(V)-1
				matrix B=e(b)
				local indB=colsof(B)-1
				replace b_coeff_flux_detect=B[1,`indB'] in `i'
				replace CI_l_flux_detect=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
				replace CI_h_flux_detect=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
				tostring(CI_l_flux_detect), replace force
				tostring(CI_h_flux_detect), replace force
				tostring(b_coeff_flux_detect),replace force 
				replace CI_95_flux_detect=substr(b_coeff_flux_detect, 1,6)+"("+substr(CI_l_flux_detect,1,6)+","+substr(CI_h_flux_detect,1,6)+")" in `i'
				destring(CI_l_flux_detect), replace force
				destring(CI_h_flux_detect), replace force
				destring(b_coeff_flux_detect), replace force
				quietly reg EX_`k'_fe_ age sex bmi b_`j', vce(robust)
				quietly test b_`j'
				replace p_val_detect_flux=r(p) in `i'
				quietly reg EX_`k'_fe_ age sex bmi pat_stat b_`j' `j', vce(robust)
				quietly test b_`j' `j'
				replace p_val_global_flux=r(p) in `i'
				quietly sum b_EX_`k'_fe_
				replace d_rate_flux=r(mean) in `i'
				quietly reg EX_`k'_fe_ b_`j' `j'
				replace r2_flux=e(r2) in `i'
				}
			*Concentrations linear
			quietly sum b_`j'
			if r(mean)>0.5{
				quietly reg ln_`k' age sex bmi pat_stat `j', vce(robust)
				local df=e(df_r)
				matrix V=e(V)
				local indV=colsof(V)-1
				matrix B=e(b)
				local indB=colsof(B)-1
				replace b_coeff_conc=B[1,`indB'] in `i'
				replace CI_l_conc=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
				replace CI_h_conc=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
				tostring(CI_l_conc), replace force
				tostring(CI_h_conc), replace force
				tostring(b_coeff_conc),replace force 
				replace CI_95_conc=substr(b_coeff_conc, 1,6)+"("+substr(CI_l_conc,1,6)+","+substr(CI_h_conc,1,6)+")" in `i'
				destring(CI_l_conc), replace force
				destring(CI_h_conc), replace force
				destring(b_coeff_conc), replace force
				quietly reg ln_`k' age sex bmi pat_stat `j', vce(robust)
				quietly test `j'
				replace p_val_linear_conc=r(p) in `i'
				}
			*Concentrations binary
			quietly sum b_`j'
			if r(mean)>0.1 & r(mean)<0.9{
				quietly reg ln_`k' age sex bmi pat_stat b_`j', vce(robust)
				local df=e(df_r)
				matrix V=e(V)
				local indV=colsof(V)-1
				matrix B=e(b)
				local indB=colsof(B)-1
				replace b_coeff_conc_detect=B[1,`indB'] in `i'
				replace CI_l_conc_detect=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
				replace CI_h_conc_detect=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
				tostring(CI_l_conc_detect), replace force
				tostring(CI_h_conc_detect), replace force
				tostring(b_coeff_conc_detect),replace force 
				replace CI_95_conc_detect=substr(b_coeff_conc_detect, 1,6)+"("+substr(CI_l_conc_detect,1,6)+","+substr(CI_h_conc_detect,1,6)+")" in `i'	
				destring(CI_l_conc_detect), replace force
				destring(CI_h_conc_detect), replace force
				destring(b_coeff_conc_detect), replace force
				quietly reg ln_`k' age sex bmi pat_stat b_`j', vce(robust)
				quietly test b_`j'
				replace p_val_detect_conc=r(p) in `i'
				quietly reg ln_`k' age sex bmi  pat_stat b_`j' `j', vce(robust)
				quietly test b_`j' `j'
				replace p_val_global_conc=r(p) in `i'
				quietly sum b_`j'
				replace d_rate_conc=r(mean) in `i'
				quietly reg ln_`k' age sex bmi pat_stat
				local r=e(r2)
				quietly reg ln_`k' b_`j' `j' age sex bmi pat_stat
				replace r2_conc=e(r2)-`r' in `i'
				}
			local i=`i'+1
			}
		export excel _var d_rate_flux b_coeff_flux CI_95_flux b_coeff_flux_detect CI_95_flux_detect r2_flux p_val_linear_flux p_val_detect_flux p_val_global_flux d_rate_conc b_coeff_conc CI_95_conc b_coeff_conc_detect CI_95_conc_detect r2_conc p_val_linear_conc p_val_detect_conc p_val_global_conc using "`file'", replace firstrow(variables)
		local i=1
		replace _var=""
		foreach j in conc flux{
			replace p_val_detect_`j'=.
			replace p_val_global_`j'=.
			replace p_val_linear_`j'=.
			replace r2_`j'=.
			replace CI_l_`j'=.
			replace CI_h_`j'=.
			replace b_coeff_`j'=.
			replace CI_95_`j'=""
			replace d_rate_`j'=.
			}
		}
	}
	
clear

*reorganisation of results 

cd A:\AGORA_2_New\Files_for_Johannes_revision\processed_data
use "CRC_AGORA2_merged.dta"

foreach j of varlist EX_12dgr180_fe_-EX_zn2_fe_{
	gen b_`j'=`j'
	replace b_`j'=1 if `j'>0.000001 & `j'!=.
	replace b_`j'=0 if `j'<0.000001
	egen std_`j'=std(`j')
	}
	
egen sum_secreted=rowtotal(b_EX*), missing
drop if sum_secreted==0 | sum_secreted==.
local list=""
drop if _C00246==.
local i=1
gen species=""
gen d_rate_species=.
foreach j of varlist Abiotrophia_defectiva-Weissella_paramesenteroides{
	gen b_`j'=`j'
	replace b_`j'=1 if `j'>0 & `j'!=.
	quietly sum b_`j'
	if r(mean)>0.1{
		local list="`list'"+" "+ "`j'"
		replace species="`j'" in `i'
		replace d_rate_species=r(mean) in `i'
		local i=`i'+1
		}
	}
cd "A:\AGORA_2_New\Files_for_Johannes_revision\results\Invivo_insilico_tables"
export excel species d_rate_species using d_rate_spec.xlsx, replace firstrow(variables)
clear
import excel d_rate_spec.xlsx, sheet("Sheet1") firstrow
drop if species==""
save d_rate_species.dta, replace
clear

cd "A:\AGORA_2_New\Files_for_Johannes_revision\results\Invivo_insilico_tables"

foreach j in `list'{
	local file="`j'"+".xlsx"
	import excel "`file'", sheet("Sheet1") firstrow
	drop if _var==""
	gen species="`j'"
	tostring(CI_95_flux CI_95_conc), replace
	local file="`j'"+".dta"
	save "`file'", replace
	clear
	}
	
use Abiotrophia_defectiva.dta, replace
foreach j in `list'{
	if "`j'"!="Abiotrophia_defectiva"{
		local file="`j'"+".dta"
		append using "`file'"
		}
	}
merge m:1 species using d_rate_species.dta
save results_final.dta, replace
drop if _var=="cholate"

*FDR correction

sort p_val_linear_conc
gen t=0 if p_val_linear_conc!=.
egen rank=seq() if p_val_linear_conc!=., by(t)
quietly sum rank
local n=r(max)
gen FDR_linear_conc=p_val_linear_conc*`n'/rank
forvalues j=2(1)`n'{
    local k=`j'-1
	local l=`j'
	local num1=FDR_linear_conc in `k'
	local num2=FDR_linear_conc in `l'
	if `num1'>`num2'{
	    replace FDR_linear_conc=`num1' in `l'
		}
	}
drop t rank
sort p_val_detect_conc
gen t=0 if p_val_detect_conc!=.
egen rank=seq() if p_val_detect_conc!=., by(t)
quietly sum rank
local n=r(max)
gen FDR_detect_conc=p_val_detect_conc*`n'/rank
forvalues j=2(1)`n'{
    local k=`j'-1
	local l=`j'
	local num1=FDR_detect_conc in `k'
	local num2=FDR_detect_conc in `l'
	if `num1'>`num2'{
	    replace FDR_detect_conc=`num1' in `l'
		}
	}
drop t rank


*In silico in vivo pattern analyses
replace b_coeff_conc=b_coeff_conc/1000
replace b_coeff_flux=b_coeff_flux/1000

egen group_met=group(_var)

cd A:\AGORA_2_New\Files_for_Johannes_revision\results\insilico_invivo_figures
			
local i=1
gen metabolite=""
gen sign_detect_flux=1 if b_coeff_flux_detect>0
gen sign_detect_conc=1 if b_coeff_conc_detect>0
replace sign_detect_flux=-1 if b_coeff_flux_detect<0
replace sign_detect_conc=-1 if b_coeff_conc_detect<0
gen sign_linear_flux=1 if b_coeff_flux>0
gen sign_linear_conc=1 if b_coeff_conc>0
replace sign_linear_flux=-1 if b_coeff_flux<0
replace sign_linear_conc=-1 if b_coeff_conc<0

foreach j in 05 1 FDR{
	foreach k in linear detect{
		gen accuracy_`j'_`k'=.
		gen exp_accuracy_`j'_`k'=.
		gen p_val_fisher_`j'_`k'=.
		gen p_val_`j'_`k'=.
		gen CI_l_`j'_`k'=.
		gen CI_h_`j'_`k'=.
		gen b_coeff_`j'_`k'=.
		gen CI_95_`j'_`k'=""
		gen r2_`j'_`k'=.
		}
	}
  
sort _var
forvalues j=1(1)52{
	local t=(`j'-1)*174+1
	local name=_var in `t'
	replace metabolite="`name'" in `i'
	foreach l in 05 1 FDR{
		foreach k in linear detect{
			if "`l'"=="05"{
				quietly tab sign_`k'_flux sign_`k'_conc if p_val_`k'_conc<0.05 & group_met==`j', ex
				replace p_val_fisher_`l'_`k'=r(p_exact) in `i'
				quietly kap sign_`k'_flux sign_`k'_conc if p_val_`k'_conc<0.05 & group_met==`j'
				replace accuracy_`l'_`k'=r(prop_o) in `i'
				replace exp_accuracy_`l'_`k'=r(prop_e) in `i'
				if "`k'"=="detect"{
					quietly reg b_coeff_conc_detect b_coeff_flux_detect if group_met==`j' & p_val_detect_conc<0.05
					replace r2_05_detect=e(r2) in `i'
					local df=e(df_r)
					matrix V=e(V)
					local indV=colsof(V)-1
					matrix B=e(b)
					local indB=colsof(B)-1
					replace b_coeff_05_detect=B[1,`indB'] in `i'
					replace CI_l_05_detect=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
					replace CI_h_05_detect=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
					tostring(CI_l_05_detect), replace force
					tostring(CI_h_05_detect), replace force
					tostring(b_coeff_05_detect),replace force 
					replace CI_95_05_detect=substr(b_coeff_05_detect, 1,6)+"("+substr(CI_l_05_detect,1,6)+","+substr(CI_h_05_detect,1,6)+")" in `i'
					if B[1,`indB']>0{
						replace p_val_05_detect=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					else{
						replace p_val_05_detect=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					destring(CI_l_05_detect), replace force
					destring(CI_h_05_detect), replace force
					destring(b_coeff_05_detect), replace force
					}
				if "`k'"=="linear"{
					quietly reg b_coeff_conc b_coeff_flux if group_met==`j' & p_val_linear_conc<0.05
					replace r2_05_linear=e(r2) in `i'
					local df=e(df_r)
					matrix V=e(V)
					local indV=colsof(V)-1
					matrix B=e(b)
					local indB=colsof(B)-1
					replace b_coeff_05_linear=B[1,`indB'] in `i'
					replace CI_l_05_linear=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
					replace CI_h_05_linear=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
					tostring(CI_l_05_linear), replace force
					tostring(CI_h_05_linear), replace force
					tostring(b_coeff_05_linear),replace force 
					replace CI_95_05_linear=substr(b_coeff_05_linear, 1,6)+"("+substr(CI_l_05_linear,1,6)+","+substr(CI_h_05_linear,1,6)+")" in `i'
					if B[1,`indB']>0{
						replace p_val_05_linear=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					else{
						replace p_val_05_linear=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					destring(CI_l_05_linear), replace force
					destring(CI_h_05_linear), replace force
					destring(b_coeff_05_linear), replace force
					}
				}
			if "`l'"=="1"{
				quietly tab sign_`k'_flux sign_`k'_conc if p_val_`k'_conc<1 & group_met==`j', ex
				replace p_val_fisher_`l'_`k'=r(p_exact) in `i'
				quietly kap sign_`k'_flux sign_`k'_conc if p_val_`k'_conc<1 & group_met==`j'
				replace accuracy_`l'_`k'=r(prop_o) in `i'
				replace exp_accuracy_`l'_`k'=r(prop_e) in `i'
				if "`k'"=="detect"{
					quietly reg b_coeff_conc_detect b_coeff_flux_detect if group_met==`j' & p_val_detect_conc<1
					replace r2_1_detect=e(r2) in `i'
					local df=e(df_r)
					matrix V=e(V)
					local indV=colsof(V)-1
					matrix B=e(b)
					local indB=colsof(B)-1
					replace b_coeff_1_detect=B[1,`indB'] in `i'
					replace CI_l_1_detect=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
					replace CI_h_1_detect=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
					tostring(CI_l_1_detect), replace force
					tostring(CI_h_1_detect), replace force
					tostring(b_coeff_1_detect),replace force 
					replace CI_95_1_detect=substr(b_coeff_1_detect, 1,6)+"("+substr(CI_l_1_detect,1,6)+","+substr(CI_h_1_detect,1,6)+")" in `i'
					if B[1,`indB']>0{
						replace p_val_1_detect=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					else{
						replace p_val_1_detect=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					destring(CI_l_1_detect), replace force
					destring(CI_h_1_detect), replace force
					destring(b_coeff_1_detect), replace force
					}
				if "`k'"=="linear"{
					quietly reg b_coeff_conc b_coeff_flux if group_met==`j' & p_val_linear_conc<1
					replace r2_1_linear=e(r2) in `i'
					local df=e(df_r)
					matrix V=e(V)
					local indV=colsof(V)-1
					matrix B=e(b)
					local indB=colsof(B)-1
					replace b_coeff_1_linear=B[1,`indB'] in `i'
					replace CI_l_1_linear=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
					replace CI_h_1_linear=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
					tostring(CI_l_1_linear), replace force
					tostring(CI_h_1_linear), replace force
					tostring(b_coeff_1_linear),replace force 
					replace CI_95_1_linear=substr(b_coeff_1_linear, 1,6)+"("+substr(CI_l_1_linear,1,6)+","+substr(CI_h_1_linear,1,6)+")" in `i'
					if B[1,`indB']>0{
						replace p_val_1_linear=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					else{
						replace p_val_1_linear=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
						}
					destring(CI_l_1_linear), replace force
					destring(CI_h_1_linear), replace force
					destring(b_coeff_1_linear), replace force
					}
				}
			if "`l'"=="FDR"{
				quietly sum sign_`k'_flux if FDR_`k'_conc<0.05 & group_met==`j'
				if r(N)>9{
					quietly tab sign_`k'_flux sign_`k'_conc if FDR_`k'_conc<0.05 & group_met==`j', ex
					replace p_val_fisher_`l'_`k'=r(p_exact) in `i'
					capture quietly kap sign_`k'_flux sign_`k'_conc if FDR_`k'_conc<0.05 & group_met==`j'
					replace accuracy_`l'_`k'=r(prop_o) in `i'
					replace exp_accuracy_`l'_`k'=r(prop_e) in `i'
					if "`k'"=="detect"{
						quietly reg b_coeff_conc_detect b_coeff_flux_detect if group_met==`j' & FDR_`k'_conc<0.05
						replace r2_FDR_detect=e(r2) in `i'
						local df=e(df_r)
						matrix V=e(V)
						local indV=colsof(V)-1
						matrix B=e(b)
						local indB=colsof(B)-1
						replace b_coeff_FDR_detect=B[1,`indB'] in `i'
						replace CI_l_FDR_detect=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
						replace CI_h_FDR_detect=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
						tostring(CI_l_FDR_detect), replace force
						tostring(CI_h_FDR_detect), replace force
						tostring(b_coeff_FDR_detect),replace force 
						replace CI_95_1_detect=substr(b_coeff_FDR_detect, 1,6)+"("+substr(CI_l_FDR_detect,1,6)+","+substr(CI_h_FDR_detect,1,6)+")" in `i'
						if B[1,`indB']>0{
							replace p_val_FDR_detect=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'	
						}
						else{
							replace p_val_FDR_detect=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
							}
						destring(CI_l_FDR_detect), replace force
						destring(CI_h_FDR_detect), replace force
						destring(b_coeff_FDR_detect), replace force
						}
					if "`k'"=="linear"{
						quietly reg b_coeff_conc b_coeff_flux if group_met==`j' & FDR_`k'_conc<0.05
						replace r2_FDR_linear=e(r2) in `i'
						local df=e(df_r)
						matrix V=e(V)
						local indV=colsof(V)-1
						matrix B=e(b)
						local indB=colsof(B)-1
						replace b_coeff_1_linear=B[1,`indB'] in `i'
						replace CI_l_FDR_linear=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
						replace CI_h_FDR_linear=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
						tostring(CI_l_FDR_linear), replace force
						tostring(CI_h_FDR_linear), replace force
						tostring(b_coeff_FDR_linear),replace force 
						replace CI_95_FDR_linear=substr(b_coeff_FDR_linear, 1,6)+"("+substr(CI_l_FDR_linear,1,6)+","+substr(CI_h_FDR_linear,1,6)+")" in `i'
						if B[1,`indB']>0{
							replace p_val_FDR_linear=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
							}
						else{
							replace p_val_FDR_linear=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
							}
						destring(CI_l_FDR_linear), replace force
						destring(CI_h_FDR_linear), replace force
						destring(b_coeff_FDR_linear), replace force
						}
					}
				}
			}
		}
	local file="`name'"+"binary"+".tif"
	*twoway (lfitci b_coeff_conc_detect b_coeff_flux_detect if _var=="`name'" & d_rate_species>0.1 & p_val_detect_conc<0.05, lcolor(red) clwidth(medthick) clpattern(dash) ciplot(rline) lcolor(red) blwidth(medthin) blpattern(longdash_shortdash) estopts()) (scatter b_coeff_conc_detect b_coeff_flux_detect, mcolor(gs7) msymbol(smdiamond_hollow) mlabel() mlabsize(tiny) mlabcolor(gs10) mlabposition(3)) (scatter b_coeff_conc_detect b_coeff_flux_detect if _var=="`name'" & d_rate_species>0.1 & p_val_detect_conc<0.05, mcolor(maroon) msymbol(smdiamond_hollow)) (scatter b_coeff_conc_detect b_coeff_flux_detect if _var=="`name'" & d_rate_species>0.1 & FDR_detect_conc<0.05, mcolor(navy) msymbol(smdiamond_hollow)) if _var=="`name'" & d_rate_species>0.1, ytitle("In vivo association statistic" "Δlog c[nmol/g]/Δspecies presence[yes/no]") ytitle(, size(medsmall)) yline(0, lwidth(medthick) lcolor(black)) ylabel(, nogrid) xtitle("In silico association statistic" "Δj[mmol/d]/Δspecies presence[yes/no]") xtitle(, size(medium)) xline(0, lwidth(medthick) lcolor(black)) title("`name'", size(vlarge) color(black)) legend(order(1 "95% Confidence interval"  2 "Linear fit (species p<0.05)" 3 "species with p>0.05" 4 "species with p<0.05" 5 "species with FDR<0.05") region(lcolor(white))) graphregion(fcolor(white) lcolor(white)) saving("`name'", replace)  xsize(5.5) ysize(5.5)
	capture twoway (lfitci b_coeff_conc_detect b_coeff_flux_detect if _var=="`name'" & d_rate_species>0.1 & p_val_detect_conc<0.05, lcolor(red) clwidth(medthick) clpattern(dash) ciplot(rline) lcolor(red) blwidth(medthin) blpattern(longdash_shortdash) estopts()) (scatter b_coeff_conc_detect b_coeff_flux_detect if _var=="`name'" & d_rate_species>0.1 & p_val_detect_conc<0.05, mcolor(maroon) msymbol(smdiamond_hollow)) (scatter b_coeff_conc_detect b_coeff_flux_detect if _var=="`name'" & d_rate_species>0.1 & FDR_detect_conc<0.05, mcolor(navy) msymbol(smdiamond_hollow)) if _var=="`name'" & d_rate_species>0.1, ytitle("In vivo association statistic" "Δlog c[nmol/g]/Δspecies presence[yes/no]") ytitle(, size(medsmall)) yline(0, lwidth(medthick) lcolor(black)) ylabel(, nogrid) xtitle("In silico association statistic" "Δj[mmol/d]/Δspecies presence[yes/no]") xtitle(, size(medium)) xline(0, lwidth(medthick) lcolor(black)) title("`name'", size(vlarge) color(black)) legend(order(1 "95% Confidence interval"  2 "Linear fit (species p<0.05)" 3 "species with p<0.05" 4 "species with FDR<0.05") region(lcolor(white))) graphregion(fcolor(white) lcolor(white)) saving("`name'", replace)  xsize(5.5) ysize(5.5)
	graph export "`file'", replace as(tif) height(975) width(1338)
	local file="`name'"+"_linear_all"+".tif"
	local name2="`name'"+"_linear_all"
	*twoway (lfitci b_coeff_conc b_coeff_flux if _var=="`name'" & d_rate_species>0.1 & p_val_linear_conc<0.05, lcolor(red) clwidth(medthick) clpattern(dash) ciplot(rline) lcolor(red) blwidth(medthin) blpattern(longdash_shortdash) estopts()) (scatter b_coeff_conc b_coeff_flux, mcolor(gs7) msymbol(smdiamond_hollow) mlabel() mlabsize(tiny) mlabcolor(gs10) mlabposition(3)) (scatter b_coeff_conc b_coeff_flux if _var=="`name'" & d_rate_species>0.1 & p_val_linear_conc<0.05, mcolor(maroon) msymbol(smdiamond_hollow)) (scatter b_coeff_conc b_coeff_flux if _var=="`name'" & d_rate_species>0.1 & FDR_linear_conc<0.05, mcolor(navy) msymbol(smdiamond_hollow)) if _var=="`name'" & d_rate_species>0.1, ytitle( "In vivo association statistic" "Δlog c[nmol/g]/Δspecies abundance[‰]") ytitle(, size(medsmall)) yline(0, lwidth(medthick) lcolor(black)) ylabel(, nogrid) xtitle("In silico association statistic" "Δj[mmol/d]/Δspecies abundance[‰]") xtitle(, size(medium)) xline(0, lwidth(medthick) lcolor(black)) title("`name'", size(vlarge) color(black)) legend(order(1 "95% Confidence interval"  2 "Linear fit (species p<0.05)" 3 "species with p>0.05" 4 "species with p<0.05" 5 "species with FDR<0.05") region(lcolor(white))) graphregion(fcolor(white) lcolor(white)) saving("`name2'", replace) xsize(5.5) ysize(5.5)
	capture twoway (lfitci b_coeff_conc b_coeff_flux if _var=="`name'" & d_rate_species>0.1 & p_val_linear_conc<0.05, lcolor(red) clwidth(medthick) clpattern(dash) ciplot(rline) lcolor(red) blwidth(medthin) blpattern(longdash_shortdash) estopts()) (scatter b_coeff_conc b_coeff_flux if _var=="`name'" & d_rate_species>0.1 & p_val_linear_conc<0.05, mcolor(maroon) msymbol(smdiamond_hollow)) (scatter b_coeff_conc b_coeff_flux if _var=="`name'" & d_rate_species>0.1 & FDR_linear_conc<0.05, mcolor(navy) msymbol(smdiamond_hollow)) if _var=="`name'" & d_rate_species>0.1, ytitle( "In vivo association statistic" "Δlog c[nmol/g]/Δspecies abundance[‰]") ytitle(, size(medsmall)) yline(0, lwidth(medthick) lcolor(black)) ylabel(, nogrid) xtitle("In silico association statistic" "Δj[mmol/d]/Δspecies abundance[‰]") xtitle(, size(medium)) xline(0, lwidth(medthick) lcolor(black)) title("`name'", size(vlarge) color(black)) legend(order(1 "95% Confidence interval"  2 "Linear fit (species p<0.05)" 3 "species with p<0.05" 4 "species with FDR<0.05") region(lcolor(white))) graphregion(fcolor(white) lcolor(white)) saving("`name2'", replace) xsize(5.5) ysize(5.5)
	graph export "`file'", replace as(tif) height(975) width(1338)
	local i=`i'+1
	}
local i=1

*Export summary statistics for Supplement
cd "A:\AGORA_2_New\Files_for_Johannes_revision\results\Invivo_insilico_tables"
save in_silico_in_vivo_pattern.dta, replace
drop if metabolite==""
keep metabolite accuracy_05_detect exp_accuracy_05_detect p_val_fisher_05_detect

sort p_val_fisher_05_detect
gen t=0 if p_val_fisher_05_detect!=.
egen rank=seq() if p_val_fisher_05_detect!=., by(t)
quietly sum rank
local n=r(max)
gen FDR=p_val_fisher_05_detect*`n'/rank
forvalues j=2(1)`n'{
    local k=`j'-1
	local l=`j'
	local num1=FDR in `k'
	local num2=FDR in `l'
	if `num1'>`num2'{
	    replace FDR=`num1' in `l'
		}
	}
drop t rank
replace FDR=1 if FDR>1
replace p_val=1 if p_val==.

cd "A:\AGORA_2_New\Files_for_Johannes_revision\results\Invivo_insilico_tables"
export excel metabolite accuracy_05_detect exp_accuracy_05_detect p_val_fisher_05_detect FDR using TableS11.xlsx, replace firstrow(variables)

clear


