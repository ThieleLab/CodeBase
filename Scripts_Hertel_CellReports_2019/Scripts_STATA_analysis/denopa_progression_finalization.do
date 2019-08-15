*Analyses DENOPA
*Longitudinal: Predicting disease progression operationalized by the UPDRS-III scale

clear
clear mata
clear matrix
set maxvar 32000
set more off
cd A:\Metabolomics_PD
capture log close

log using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\log_files\analyses_longitudinal.log", replace

use "A:\Hertel\Finalized\Hertel_2018_DeNoPa\data\denopa_final_data_wide.dta", replace

set obs 4000

*Outlier filtering	

foreach j of varlist CA1-FA_2401 CA2-FA_2402 CA3-FA_2403{
	gen ln_`j'=ln(`j')
	egen std_`j'=std(ln_`j')
	replace ln_`j'=. if abs(std_`j')>4
	replace std_`j'=. if abs(std_`j')>4
	}


gen var_=""

gen Methoxy_bin=0 if ln_Methoxytyrosine2<-4
replace Methoxy_bin=1 if ln_Methoxytyrosine2>-4 & ln_Methoxytyrosine2<999999
egen azilect=group(b_azilect_einnahme)
egen pramipexole=group(b_dopa_agonisten_pramipexole)
replace pramipexole=1 if b_dopa_agonisten_pramipexole==""
egen agonisten_einnahme=group(b_l_dopa_agonisten_einnahme)

cd "A:\Hertel\Finalized\Hertel_2018_DeNoPa\tables"
*****************************************************************************
***Predicting progression with baseline-values: Follow-up 1/2 UPDRS-total****
*****************************************************************************

foreach k in _I_sum _II_sum _III_sum _IV_sum _sum{
	gen p_val`k'=.
	gen CI_l`k'=.
	gen CI_h`k'=.
	gen b_coeff`k'=.
	gen OR`k'=.
	gen CI_95`k'=""
	}


foreach k in _I_sum _II_sum _III_sum _IV_sum _sum {
	local i=1
	local nam="`k'"
	foreach j of varlist GCDCA1 TDCA1 TLCA1 TCDCA1 LHistidine1{
		if "`nam'"=="_I_sum"{
			local len=strlen("`j'")-1
			replace var_=substr("`j'",1,`len') in `i'
			}
		if "`nam'"=="_III_sum" | "`nam'"=="_sum"{	
			quietly reg UPDRS`k'_2 UPDRS`k'_1 c.age sex months_of_disease std_`j' if group_num==1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'),replace force 
			replace CI_95`k'=substr(b_coeff`k', 1,6)+"("+substr(CI_l`k',1,6)+","+substr(CI_h`k',1,6)+")" in `i'
			if B[1,`indB']>0{
				replace p_val`k'=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			else{
				replace p_val`k'=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'),replace force
			}
		else{
			quietly ologit UPDRS`k'_2 UPDRS`k'_1 c.age sex months_of_disease std_`j' if group_num==1
			local aux=e(k_aux)
			matrix V=e(V)
			local indV=colsof(V)-`aux'
			matrix B=e(b)
			local indB=colsof(B)-`aux'
			replace OR`k'=exp(B[1,`indB']) in `i'
			replace CI_l`k'=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
			replace CI_h`k'=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(OR`k'), replace force
			replace CI_95`k'=substr(OR`k', 1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			if B[1,`indB']>0{
				replace p_val`k'=2*normal(-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			else{
				replace p_val`k'=2*normal(B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(OR`k'), replace force
			}
		local i=`i'+1
		}
	}
export excel var_ CI_95_I_sum p_val_I_sum CI_95_II_sum p_val_II_sum CI_95_III_sum p_val_III_sum CI_95_IV_sum p_val_IV_sum CI_95_sum p_val_sum using predict_UPDRS_all_21_baseline.xlsx, replace firstrow(variables)

foreach k in _I_sum _II_sum _III_sum _IV_sum _sum{
	replace p_val`k'=.
	replace CI_l`k'=.
	replace CI_h`k'=.
	replace b_coeff`k'=.
	replace OR`k'=.
	replace CI_95`k'=""
	}
replace var_=""

local i=1


*****************************************************************************
***Predicting progression with baseline-values: Follow-up 1/3 UPDRS-total****
*****************************************************************************
/*
foreach k in _I_sum _II_sum _III_sum _IV_sum _sum{
	gen p_val`k'=.
	gen CI_l`k'=.
	gen CI_h`k'=.
	gen b_coeff`k'=.
	gen OR`k'=.
	gen CI_95`k'=""
	}*/


foreach k in _I_sum _II_sum _III_sum _IV_sum _sum {
	local i=1
	local nam="`k'"
	foreach j of varlist GCDCA1 TDCA1 TLCA1 TCDCA1 LHistidine1{
		if "`nam'"=="_I_sum"{
			local len=strlen("`j'")-1
			replace var_=substr("`j'",1,`len') in `i'
			}
		if "`nam'"=="_III_sum" | "`nam'"=="_sum"{	
			quietly reg UPDRS`k'_3 UPDRS`k'_1 c.age sex Methoxy_bin b_levodopa_aequivalenzdosis_dopa months_of_disease std_`j' if group_num==1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'),replace force 
			replace CI_95`k'=substr(b_coeff`k', 1,6)+"("+substr(CI_l`k',1,6)+","+substr(CI_h`k',1,6)+")" in `i'
			if B[1,`indB']>0{
				replace p_val`k'=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			else{
				replace p_val`k'=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		else{
			quietly ologit UPDRS`k'_3 UPDRS`k'_1 c.age sex Methoxy_bin b_levodopa_aequivalenzdosis_dopa months_of_disease std_`j' if group_num==1
			local aux=e(k_aux)
			matrix V=e(V)
			local indV=colsof(V)-`aux'
			matrix B=e(b)
			local indB=colsof(B)-`aux'
			replace OR`k'=exp(B[1,`indB']) in `i'
			replace CI_l`k'=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
			replace CI_h`k'=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(OR`k'), replace force
			replace CI_95`k'=substr(OR`k', 1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			if B[1,`indB']>0{
				replace p_val`k'=2*normal(-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			else{
				replace p_val`k'=2*normal(B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(OR`k'), replace force
			}
		local i=`i'+1
		}
	}
export excel var_ CI_95_I_sum p_val_I_sum CI_95_II_sum p_val_II_sum CI_95_III_sum p_val_III_sum CI_95_IV_sum p_val_IV_sum CI_95_sum p_val_sum using predict_UPDRS_all_31_baseline.xlsx, replace firstrow(variables)

foreach k in _I_sum _II_sum _III_sum _IV_sum _sum{
	replace p_val`k'=.
	replace CI_l`k'=.
	replace CI_h`k'=.
	replace b_coeff`k'=.
	replace OR`k'=.
	replace CI_95`k'=""
	}
replace var_=""

local i=1



**************************************************************************
***Predicting progression with baseline-values: Follow-up 2/3 UPDRS-total****
**************************************************************************
/*
foreach k in _I_sum _II_sum _III_sum _IV_sum _sum{
	gen p_val`k'=.
	gen CI_l`k'=.
	gen CI_h`k'=.
	gen b_coeff`k'=.
	gen OR`k'=.
	gen CI_95`k'=""
	}*/


foreach k in _I_sum _II_sum _III_sum _IV_sum _sum {
	local i=1
	local nam="`k'"
	foreach j of varlist GCDCA2 TDCA2 TLCA2 TCDCA2 LHistidine2{
		if "`nam'"=="_I_sum"{
			local len=strlen("`j'")-1
			replace var_=substr("`j'",1,`len') in `i'
			}
		if "`nam'"=="_III_sum" | "`nam'"=="_sum"{	
			quietly reg UPDRS`k'_3 UPDRS`k'_2 c.age sex Methoxy_bin b_levodopa_aequivalenzdosis_dopa months_of std_`j' if group_num==1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invt(`df',0.975)*sqrt(V[`indV',`indV']) in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'),replace force 
			replace CI_95`k'=substr(b_coeff`k', 1,6)+"("+substr(CI_l`k',1,6)+","+substr(CI_h`k',1,6)+")" in `i'
			if B[1,`indB']>0{
				replace p_val`k'=2*t(`df',-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			else{
				replace p_val`k'=2*t(`df',B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		else{
			quietly ologit UPDRS`k'_3 UPDRS`k'_2 c.age sex Methoxy_bin b_levodopa_aequivalenzdosis_dopa months_of std_`j' if group_num==1
			local aux=e(k_aux)
			matrix V=e(V)
			local indV=colsof(V)-`aux'
			matrix B=e(b)
			local indB=colsof(B)-`aux'
			replace OR`k'=exp(B[1,`indB']) in `i'
			replace CI_l`k'=exp(B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i' 
			replace CI_h`k'=exp(B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])) in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(OR`k'), replace force
			replace CI_95`k'=substr(OR`k', 1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			if B[1,`indB']>0{
				replace p_val`k'=2*normal(-B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			else{
				replace p_val`k'=2*normal(B[1,`indB']/sqrt(V[`indV',`indV'])) in `i'
				}
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(OR`k'), replace force
			}
		local i=`i'+1
		}
	}
export excel var_ CI_95_I_sum p_val_I_sum CI_95_II_sum p_val_II_sum CI_95_III_sum p_val_III_sum CI_95_IV_sum p_val_IV_sum CI_95_sum p_val_sum using predict_UPDRS_all_32_baseline.xlsx, replace firstrow(variables)

foreach k in _I_sum _II_sum _III_sum _IV_sum _sum{
	replace p_val`k'=.
	replace CI_l`k'=.
	replace CI_h`k'=.
	replace b_coeff`k'=.
	replace OR`k'=.
	replace CI_95`k'=""
	}
replace var_=""

local i=1


clear

log close
