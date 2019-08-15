*Analyses DENOPA
*Matching considered to be done on group level

clear
clear mata
clear matrix
capture log close
set maxvar 32000
set more off

log using "A:\Hertel\Finalized\Hertel_2018_DeNoPa\log_files\analyses_longitudinal_final.log", replace

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

/*
*Box plots after outlier filtering
cd "A:\Metabolomics_PD\Analyses_Denopa\Graphs\BoxPlots"

foreach j of varlist CA-FA_240{
	graph box ln_`j' if filter_`j'==0, over(wave) by(group_num)
	graph export `j'.png, replace
	}
*Histograms after outlier filtering
cd "A:\Metabolomics_PD\Analyses_Denopa\Graphs\Histograms"
foreach j of varlist CA-FA_240{
	histogram ln_`j' if filter_`j'==0
	graph export `j'_hist.png, replace
	}
*/	
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
gen ln_crp=ln(labor_crp_mg_dl)
gen ln_ggt=ln(labor_gamma_gt_ggt)


*Descriptive statistics, Table 1
ttest age if wave==1, by(group_num) welch
sum months_of if wave==1, det
tab sex group_num if wave==1, ex col	
foreach j in 1 2 3{
	tab Methoxy_bin group_num if wave==`j', col ex
	ttest Methoxytyros if wave==`j', welch by(group_num)
	ttest UPDRS_sum_ if wave==`j', welch by(group_num)
	ttest bdi_sum if wave==`j', welch by(group_num)
	ttest labor_triglyzeride_trig if wave==`j', welch by(group_num)
	ttest labor_gamma_gt_ggt if wave==`j', welch by(group_num)
	ttest ln_ggt if wave==`j', welch by(group_num)
	}
	
********************************************************************	
****Descriptive Table: Concentrations over waves and study group****
********************************************************************

*Summary Statistics S1
*Supplementary table S1

cd A:\Hertel\Finalized\Hertel_2018_DeNoPa\tables
gen var_=""

foreach l in control PD{
	gen p_val_`l'=.
	gen ICC_`l'=.
	foreach t in 1 2 3{
		gen mean_`l'_`t'=.
		gen sd_`l'_`t'=.
		gen m_sd_`l'_`t'=""
		}
	}
local i=1
foreach j in . 1{
	foreach k of varlist CA-FA_240 UPDRS_I_sum_ UPDRS_II_sum_ UPDRS_III_sum_ UPDRS_IV_sum_ UPDRS_sum_ Methoxy_bin{
		local name="`k'"
		replace var_="`k'" in `i'
		foreach l in control PD{
			foreach t in 1 2 3{
				if substr("`name'",1,5)!="UPDRS" & "`name'"!="Methoxy_bin"{
					if "`l'"=="PD"{
						quietly sum `k' if filter_`k'!=`j' & group_num==1 & wave==`t'
						replace mean_`l'_`t'=round(1000*r(mean))/1000 in `i'
						replace sd_`l'_`t'=round(1000*r(sd))/1000 in `i'
						tostring(mean_`l'_`t' sd_`l'_`t'), replace force
						replace m_sd_`l'_`t'=substr(mean_`l'_`t',1,5)+"("+substr(sd_`l'_`t',1,5)+")"
						destring(mean_`l'_`t' sd_`l'_`t'), replace force
						}
					else{
						quietly sum `k' if filter_`k'!=`j' & group_num==0 & wave==`t'
						replace mean_`l'_`t'=round(1000*r(mean))/1000 in `i'
						replace sd_`l'_`t'=round(1000*r(sd))/1000 in `i'
						tostring(mean_`l'_`t' sd_`l'_`t'), replace force
						replace m_sd_`l'_`t'=substr(mean_`l'_`t',1,5)+"("+substr(sd_`l'_`t',1,5)+")"
						destring(mean_`l'_`t' sd_`l'_`t'), replace force
						}
					}
				else{
					if "`l'"=="PD"{
						quietly sum `k' if group_num==1 & wave==`t'
						replace mean_`l'_`t'=round(1000*r(mean))/1000 in `i'
						replace sd_`l'_`t'=round(1000*r(sd))/1000 in `i'
						tostring(mean_`l'_`t' sd_`l'_`t'), replace force
						replace m_sd_`l'_`t'=substr(mean_`l'_`t',1,5)+"("+substr(sd_`l'_`t',1,5)+")"
						destring(mean_`l'_`t' sd_`l'_`t'), replace force
						}
					else{
						quietly sum `k' if  group_num==0 & wave==`t'
						replace mean_`l'_`t'=round(1000*r(mean))/1000 in `i'
						replace sd_`l'_`t'=round(1000*r(sd))/1000 in `i'
						tostring(mean_`l'_`t' sd_`l'_`t'), replace force
						replace m_sd_`l'_`t'=substr(mean_`l'_`t',1,5)+"("+substr(sd_`l'_`t',1,5)+")"
						destring(mean_`l'_`t' sd_`l'_`t'), replace force
						}
					}
				}
			if "`name'"=="Methoxy_bin"{	
				if "`l'"=="PD"{
					quietly xtlogit Methoxy_bin i.wave if group_num==1 , re
					replace ICC_`l'=round(100*e(rho))/100 in `i'
					quietly test 2.wave 3.wave
					replace p_val_`l'=r(p) in `i'
					}
				}
			if "`name'"=="UPDRS_I_sum_" | "`name'"=="UPDRS_II_sum_" | "`name'"=="UPDRS_IV_sum_" | "`name'"=="UPDRS_sum_"| "`name'"=="UPDRS_III_sum_"{
				if "`l'"=="PD"{
					quietly xtreg `k' i.wave if group_num==1,re
					replace ICC_`l'=round(100*e(rho))/100 in `i'
					quietly test 2.wave 3.wave
					replace p_val_`l'=r(p) in `i'
					}
				else{
					quietly xtreg `k' i.wave if group_num==0,re
					replace ICC_`l'=round(100*e(rho))/100 in `i'
					quietly test 2.wave 3.wave
					replace p_val_`l'=r(p) in `i'
					}
				}
			if "`name'"=="UPDRS_I_sum_" | "`name'"=="UPDRS_II_sum_" | "`name'"=="UPDRS_IV_sum_" {
				if "`l'"=="PD"{
					quietly xtologit `k' i.wave if group_num==1
					quietly test 2.wave 3.wave
					replace p_val_`l'=r(p) in `i'
					}
				else{
					quietly xtologit `k' i.wave if group_num==0
					quietly test 2.wave 3.wave
					replace p_val_`l'=r(p) in `i'
					}
				}
			if substr("`name'",1,5)!="UPDRS" & "`name'"!="Methoxy_bin"{
				if "`l'"=="PD"{
					quietly xtreg ln_`k' i.wave if filter_`k'!=`j' & group_num==1 , re
					replace ICC_`l'=round(100*e(rho))/100 in `i'
					quietly test 2.wave 3.wave
					replace p_val_`l'=r(p) in `i'
					}
				else{
					quietly xtreg ln_`k' i.wave if filter_`k'!=`j' & group_num==0 , re
					replace ICC_`l'=round(100*e(rho))/100 in `i'
					quietly test 2.wave 3.wave
					replace p_val_`l'=r(p) in `i'
					}
				}
			}
		local i=`i'+1
		}
	if `j'==.{
		export excel var_ m_sd_control_1 m_sd_control_2 m_sd_control_3 ICC_control p_val_control m_sd_PD_1 m_sd_PD_2 m_sd_PD_3 ICC_PD p_val_PD using descr_table_with_out.xlsx, replace firstrow(variables)
		}
	else{
		export excel var_ m_sd_control_1 m_sd_control_2 m_sd_control_3 ICC_control p_val_control m_sd_PD_1 m_sd_PD_2 m_sd_PD_3 ICC_PD p_val_PD using descr_table.xlsx, replace firstrow(variables)
		}
	foreach s in control PD{
		replace p_val_`s'=.
		replace ICC_`s'=.
		foreach q in 1 2 3{
			replace mean_`s'_`q'=.
			replace sd_`s'_`q'=.
			replace m_sd_`s'_`q'=""
			}
		}
	replace var_=""
	local i=1
	}
drop p_val_control-m_sd_PD_3 			

***************************************************************
****Simple longitudinal association studies using RE and FE****
***************************************************************

*Summary Statistics S2

foreach k in nad ad{ 
	foreach j in re fe{
		gen p_val_main_`j'_`k'=.
		gen p_val_inter_`j'_`k'=.
		gen p_val_global_`j'_`k'=.
		}
	gen CI_l`k'=.
	gen CI_h`k'=.
	gen b_coeff`k'=.
	gen CI_95`k'=""
	}
gen p_val_haus=.
local i=1		
replace levodopa_dosis_mg=0 if Methoxy_bin==0
*categorical parametrization
foreach l in nad ad{
	local i=1
	local ad="`l'"
	foreach j of varlist CA-FA_240{
		replace var_="`j'" in `i'
		foreach k in fe re{
			local name="`k'"
			if "`name'"=="re"{
				if "`ad'"=="nad"{
					quietly xtreg ln_`j' age sex i.wave group_num if filter_`j'==0, vce(robust)
					}
				if "`ad'"=="ad"{
					quietly xtreg ln_`j' age sex i.wave levo_equivalent group_num if filter_`j'==0, vce(robust)
					}
				local df=e(df_r)
				matrix V=e(V)
				local indV=colsof(V)-1
				matrix B=e(b)
				local indB=colsof(B)-1
				replace b_coeff`l'=B[1,`indB'] in `i'
				replace CI_l`l'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
				replace CI_h`l'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
				tostring(CI_l`l'), replace force
				tostring(CI_h`l'), replace force
				tostring(b_coeff`l'), replace force
				replace CI_95`l'=substr(b_coeff`l',1,6)+"("+substr(CI_l`l',1,6)+","+substr(CI_h`l',1,6)+")" in `i'
				destring(b_coeff`l'), replace force
				destring(CI_l`l'), replace force
				destring(CI_h`l'), replace force
				if "`ad'"=="nad"{
					quietly xtreg ln_`j' age sex i.wave group_num if filter_`j'==0, `k' vce(robust)
					quietly test group_num
					replace p_val_main_`k'_`l'=r(p) in `i'
					quietly xtreg ln_`j' age sex i.wave##group_num if filter_`j'==0, `k' vce(robust)
					quietly test 2.wave#1.group_num 3.wave#1.group_num
					replace p_val_inter_`k'_`l'=r(p) in `i'
					quietly test 2.wave#1.group_num 3.wave#1.group_num 1.group_num
					replace p_val_global_`k'_`l'=r(p) in `i'
					}
				if "`ad'"=="ad"{
					quietly xtreg ln_`j' age sex i.wave levo_equivalent group_num if filter_`j'==0, `k' vce(robust)
					quietly test group_num
					replace p_val_main_`k'_`l'=r(p) in `i'
					quietly xtreg ln_`j' age sex c.levo_equivalent##i.wave##group_num if filter_`j'==0, `k' vce(robust)
					quietly test 2.wave#1.group_num 3.wave#1.group_num
					replace p_val_inter_`k'_`l'=r(p) in `i'
					quietly test 2.wave#1.group_num 3.wave#1.group_num 1.group_num
					replace p_val_global_`k'_`l'=r(p) in `i'
					}
				}
			else{
				if "`ad'"=="nad"{
					quietly xtreg ln_`j' age sex i.wave##group_num  if filter_`j'==0, `k' vce(robust)
					}
				if "`ad'"=="ad"{
					quietly xtreg ln_`j' age sex i.wave##group_num##c.levo_equivalent if filter_`j'==0, `k' vce(robust)
					}
				quietly test 2.wave#1.group_num 3.wave#1.group_num
				replace p_val_inter_`k'_`l'=r(p) in `i'
				}
			}
		quietly xtreg `j' i.wave if filter_`j'==0, fe
		est store FE
		quietly xtreg `j' i.wave if filter_`j'==0, re
		est store RE
		quietly hausman FE RE
		replace p_val_haus=r(p) in `i'
		local i=`i'+1
		}
	}
export excel var_ p_val_inter_fe_nad CI_95nad p_val_main_re_nad p_val_inter_re_nad p_val_global_re_nad p_val_inter_fe_ad CI_95ad p_val_main_re_ad p_val_inter_re_ad p_val_global_re_ad p_val_haus using "summary_statistics_S2.xlsx", replace firstrow(variables)
	
foreach k in nad ad{ 
	foreach j in re fe{
		replace p_val_main_`j'_`k'=.
		replace p_val_inter_`j'_`k'=.
		replace p_val_global_`j'_`k'=.
		}
	replace CI_l`k'=.
	replace CI_h`k'=.
	replace b_coeff`k'=.
	replace CI_95`k'=""
	}
replace p_val_haus=.
local i=1		
replace var_=""

***********************************************
***Associations with 3MOT**********************
***********************************************

gen p_val_main=.
gen p_val_inter=.
gen p_val_global=.
gen CI_l=.
gen CI_h=.
gen b_coeff=.
gen CI_95=""
local i=1	

foreach j of varlist CA-Valerylcarni Citrulline-FA_240{
	replace var_="`j'" in `i'
	quietly xtreg ln_`j' age sex i.wave ln_Methoxytyros if filter_`j'==0 & group_num==1, re vce(robust)
	quietly test ln_Methoxy
	replace p_val_main=r(p) in `i'
	quietly xtreg ln_`j' age sex i.wave ln_Methoxytyros if filter_`j'==0 & group_num==1, re vce(robust)
	local df=e(df_r)
	matrix V=e(V)
	local indV=colsof(V)-1
	matrix B=e(b)
	local indB=colsof(B)-1
	replace b_coeff=B[1,`indB'] in `i'
	replace CI_l=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
	replace CI_h=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
	tostring(CI_l), replace force
	tostring(CI_h), replace force
	tostring(b_coeff), replace force
	replace CI_95=substr(b_coeff,1,6)+"("+substr(CI_l,1,6)+","+substr(CI_h,1,6)+")" in `i'
	destring(b_coeff), replace force
	destring(CI_l), replace force
	destring(CI_h), replace force
	local i=`i'+1
	}
export excel var_ CI_95 p_val_main using "metabolome_methoxy_PD.xlsx", replace firstrow(variables)

replace p_val_main=.
replace p_val_inter=.
replace p_val_global=.
replace CI_l=.
replace CI_h=.
replace b_coeff=.
replace CI_95=""
replace var_=""

local i=1

***********************************************
****association symptoms in PD cases***********
***********************************************
*replace levo_equivalent=0 if Methoxy_bin==0
foreach k in _I_sum _II_sum _III_sum _IV_sum _sum {
	gen p_val_main`k'=.
	gen p_val_inter`k'=.
	gen p_val_global`k'=.
	gen CI_l`k'=.
	gen CI_h`k'=.
	gen b_coeff`k'=.
	gen OR`k'=.
	gen CI_95`k'=""
	}

foreach k in  _I_sum _II_sum _III_sum _IV_sum _sum {
	local i=1
	local nam="`k'"
	foreach j of varlist CA-FA_240{
		if "`nam'"=="_I_sum" | "`nam'"=="_II_sum" | "`nam'"=="_IV_sum"{
			replace var_="`j'" in `i'
			quietly xtologit c.UPDRS`k' age sex i.wave months_of levo_equi std_`j' if filter_`j'==0 & group_num==1
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
			replace CI_95`k'=substr(OR`k',1,4)+"("+substr(CI_l`k',1,4)+","+substr(CI_h`k',1,4)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(OR`k'), replace force
			quietly xtologit c.UPDRS`k' std_`j' age sex i.wave months_of levo_equi if filter_`j'==0 & group_num==1
			quietly test c.std_`j'
			replace p_val_main`k'=r(p) in `i'
			quietly xtologit c.UPDRS`k' age sex i.wave##c.std_`j' months_of i.wave##c.levo_equi if filter_`j'==0 & group_num==1
			quietly test 2.wave#c.std_`j' 3.wave#c.std_`j'
			replace p_val_inter`k'=r(p) in `i'
			quietly test 2.wave#c.std_`j' 3.wave#c.std_`j' c.std_`j'
			replace p_val_global`k'=r(p) in `i'
			local i=`i'+1
			}
		else{
			replace var_="`j'" in `i'
			quietly xtreg c.UPDRS`k' age sex i.wave months_of levo_equi std_`j' if filter_`j'==0 & group_num==1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,6)+"("+substr(CI_l`k',1,6)+","+substr(CI_h`k',1,6)+")" in `i'
			destring(b_coeff`k'), replace force
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			quietly xtreg c.UPDRS`k' std_`j' age sex i.wave months_of  levo_equi if filter_`j'==0 & group_num==1, vce(robust)
			quietly test c.std_`j'
			replace p_val_main`k'=r(p) in `i'
			quietly xtreg c.UPDRS`k' age sex i.wave##c.std_`j' months_of i.wave##c.levo_equi if filter_`j'==0 & group_num==1, vce(robust)
			quietly test 2.wave#c.std_`j' 3.wave#c.std_`j'
			replace p_val_inter`k'=r(p) in `i'
			quietly test 2.wave#c.std_`j' 3.wave#c.std_`j' c.std_`j'
			replace p_val_global`k'=r(p) in `i'
			local i=`i'+1
			}
		}
	}


export excel var_ CI_95_I_sum p_val_main_I_sum p_val_inter_I_sum p_val_global_I_sum CI_95_II_sum p_val_main_II_sum p_val_inter_II_sum p_val_global_II_sum CI_95_III_sum p_val_main_III_sum p_val_inter_III_sum p_val_global_III_sum CI_95_IV_sum p_val_main_IV_sum p_val_inter_IV_sum p_val_global_IV_sum CI_95_sum p_val_main_sum p_val_inter_sum p_val_global_sum using "summary_statistics_S3.xlsx", replace firstrow(variables)
	
foreach k in _I_sum _II_sum _III_sum _IV_sum _sum {
	replace p_val_main`k'=.
	replace p_val_inter`k'=.
	replace p_val_global`k'=.
	replace CI_l`k'=.
	replace CI_h`k'=.
	replace b_coeff`k'=.
	replace OR`k'=.
	replace CI_95`k'=""
	}

replace var_=""
local i=1

**********************************************************************************************
****Simple longitudinal association studies on Medication status and levodopa dosis***********
**********************************************************************************************
foreach k in bin azi prami dos ago equi{
	gen p_val`k'=.
	gen CI_l`k'=.
	gen CI_h`k'=.
	gen b_coeff`k'=.
	gen CI_95`k'=""
	}
gen p_int_global=.
gen p_dos_global=.

foreach j of varlist prami ropi roti piri{
	replace `j'=0 if b_l_dopa_agonisten_einnahme=="false"
	}
replace levodopa_dos2=levodopa_dos2/100	
replace levo_equi=levo_equi/100
replace levo_diff=levo_equi-levodopa_dos2
foreach k in bin dos azi prami ago equi{	
	local i=1
	local nam="`k'"
	foreach j of varlist CA-FA_240{
		replace var_="`j'" in `i'
		if "`nam'"=="bin"{
			quietly xtreg ln_`j' age sex i.wave months_of Methoxy_bin pramipexole azilect if filter_`j'==0 & group_num==1 & wave!=1, vce(robust)
			quietly test pramipexole azilect Methoxy_bin
			replace p_int_global=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of Methoxy_bin if filter_`j'==0 & group_num==1 , vce(robust)
			quietly test Methoxy_bin
			replace p_val`k'=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of Methoxy_bin if filter_`j'==0 & group_num==1 , vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		if "`nam'"=="dos"{
			quietly xtreg ln_`j' age sex i.wave months_of levodopa_dos2 levo_equi if filter_`j'==0 & group_num==1 & wave!=1, vce(robust)
			quietly test levodopa_dos2 levo_equi
			replace p_dos_global=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of levodopa_dos2 if filter_`j'==0 & group_num==1 & Methoxy_bin==1 & wave!=1, vce(robust)
			quietly test levodopa_dos2
			replace p_val`k'=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of levodopa_dos2 if filter_`j'==0 & group_num==1 & Methoxy_bin==1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		if "`nam'"=="equi"{
			quietly xtreg ln_`j' age sex i.wave months_of levo_equi if filter_`j'==0 & group_num==1 & wave!=1 & levo_equi>0, vce(robust)
			quietly test levo_equi
			replace p_val`k'=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of levo_equi if filter_`j'==0 & group_num==1 & wave!=1 & levo_equi>0, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		if "`nam'"=="prami"{
			quietly xtreg ln_`j' age sex i.wave months_of prami if filter_`j'==0 & group_num==1 & wave!=1 , vce(robust)
			quietly test prami
			replace p_val`k'=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of prami if filter_`j'==0 & group_num==1 & wave!=1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		if "`nam'"=="azi"{
			quietly xtreg ln_`j' age sex i.wave months_of azilect if filter_`j'==0 & group_num==1 & wave!=1 , vce(robust)
			quietly test azilect
			replace p_val`k'=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of azilect if filter_`j'==0 & group_num==1 & wave!=1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		if "`nam'"=="ago"{
			quietly xtreg ln_`j' age sex i.wave months_of dopa_agonisten_intake if filter_`j'==0 & group_num==1 & wave!=1 , vce(robust)
			quietly test dopa_agonisten_intake
			replace p_val`k'=r(p) in `i'
			quietly xtreg ln_`j' age sex i.wave months_of dopa_agonisten_intake if filter_`j'==0 & group_num==1 & wave!=1, vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}	
		local i=`i'+1	
		}
	}
export excel var_ CI_95bin p_valbin CI_95azi p_valazi CI_95prami p_valprami p_int_global CI_95dos p_valdos CI_95equi p_valequi CI_95ago p_valago p_dos_global using "summary_statistics_S5.xlsx", replace firstrow(variables)
	
foreach k in bin dos azi prami ago equi{
	replace p_val`k'=.
	replace CI_l`k'=.
	replace CI_h`k'=.
	replace b_coeff`k'=.
	replace CI_95`k'=""
	}
replace p_int_global=.
replace p_dos_global=.	
replace var_=""
local i=1


***********************************************************************
***Associations changes of medication with changes in concentrations***
***********************************************************************


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
	
foreach j of varlist levo_equivalent levodopa_dos2 levo_diff UPDRS_*{
	egen mean_`j'=mean(`j') if wave!=3, by(pr_id)
	egen m2_`j'=mean(`j'), by(group_num)
	gen `j'_cen1=`j'-mean_`j'
	egen std_cen=std(`j'_cen1)
	replace `j'_cen1=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}
	
	
foreach j of varlist levo_equivalent levodopa_dos2  levo_diff UPDRS_I_sum_-UPDRS_sum_{
	egen mean_`j'=mean(`j') if wave!=1, by(pr_id)
	egen m2_`j'=mean(`j'), by(group_num)
	gen `j'_cen2=`j'-mean_`j'
	egen std_cen=std(`j'_cen2)
	replace `j'_cen2=. if abs(std_cen)>3
	drop std_cen m2_`j' mean_`j'
	}	
	
foreach j of varlist CA-FA_240 levo_equivalent levodopa_dos2 levo_diff UPDRS_I_sum_-UPDRS_sum_{
	replace `j'_cen1=`j'_cen2 if wave==3
	replace `j'_cen1=. if wave==1
	}
gen wave2=wave if wave!=1

foreach k in dos equi diff{
	capture gen p_val`k'=.
	capture gen CI_l`k'=.
	capture gen CI_h`k'=.
	capture gen b_coeff`k'=.
	capture gen CI_95`k'=""
	}
foreach k in dos equi diff{	
	local i=1
	local nam="`k'"
	foreach j of varlist CA-FA_240{
		replace var_="`j'" in `i'
		if "`nam'"=="dos"{
			quietly xtreg `j'_cen1 age sex i.wave2 months_of levodopa_dos2_cen1 if filter_`j'==0 & group_num==1, vce(robust)
			quietly test levodopa_dos2_cen1
			replace p_val`k'=r(p) in `i'
			quietly xtreg `j'_cen1 age sex i.wave2 months_of levodopa_dos2_cen1 if filter_`j'==0 & group_num==1 , vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		if "`nam'"=="equi"{
			quietly xtreg `j'_cen1 age sex i.wave2 months_of levo_equivalent_cen1 if filter_`j'==0 & group_num==1, vce(robust)
			quietly test levo_equivalent_cen1
			replace p_val`k'=r(p) in `i'
			quietly xtreg `j'_cen1 age sex i.wave2 months_of levo_equivalent_cen1 if filter_`j'==0 & group_num==1 , vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		if "`nam'"=="diff"{
			quietly xtreg `j'_cen1 age sex i.wave2 months_of levo_diff_cen1 if filter_`j'==0 & group_num==1, vce(robust)
			quietly test levo_diff_cen1
			replace p_val`k'=r(p) in `i'
			quietly xtreg `j'_cen1 age sex i.wave2 months_of levo_diff_cen1 if filter_`j'==0 & group_num==1 , vce(robust)
			local df=e(df_r)
			matrix V=e(V)
			local indV=colsof(V)-1
			matrix B=e(b)
			local indB=colsof(B)-1
			replace b_coeff`k'=B[1,`indB'] in `i'
			replace CI_l`k'=B[1,`indB']-invnormal(0.975)*sqrt(V[`indV',`indV']) in `i' 
			replace CI_h`k'=B[1,`indB']+invnormal(0.975)*sqrt(V[`indV',`indV'])  in `i'
			tostring(CI_l`k'), replace force
			tostring(CI_h`k'), replace force
			tostring(b_coeff`k'), replace force
			replace CI_95`k'=substr(b_coeff`k',1,5)+"("+substr(CI_l`k',1,5)+","+substr(CI_h`k',1,5)+")" in `i'
			destring(CI_l`k'), replace force
			destring(CI_h`k'), replace force
			destring(b_coeff`k'), replace force
			}
		local i=`i'+1	
		}
	}
export excel var_ CI_95equi p_valequi CI_95dos p_valdos CI_95diff p_valdiff using "summary_statistics_S4.xlsx", replace firstrow(variables)
	
foreach k in  dos equi diff{
	replace p_val`k'=.
	replace CI_l`k'=.
	replace CI_h`k'=.
	replace b_coeff`k'=.
	replace CI_95`k'=""
	}
replace var_=""
local i=1
log close

log using bivariate_transsulfuration.log, replace

*Bivariate analyses: Transsulfuration

xtreg ln_LMethionine age sex i.wave##group_num group_num##c.ln_LHomoserine if filter_LMethionine==0 & filter_LHomoserine==0, vce(robust)
xtreg ln_LMethionine age sex i.wave##group_num group_num##c.ln_LHomoserine c.levodopa_dosis_mg##i.wave if filter_LMethionine==0 & filter_LHomoserine==0, vce(robust)
xtreg ln_LMethionine age sex i.wave##group_num group_num##c.ln_LHomoserine c.levo_equivalent##i.wave if filter_LMethionine==0 & filter_LHomoserine==0, vce(robust)

*Cystein-Cystathionine
xtreg ln_Cystein age sex i.wave##group_num group_num##c.ln_Cystathion if filter_Cystathion==0 & filter_Cystein==0, vce(robust)
xtreg ln_Cystein age sex i.wave##group_num group_num##c.ln_Cystathion c.levodopa_dosis_mg##i.wave if filter_Cystathion==0 & filter_Cystein==0, vce(robust)
xtreg ln_Cystein age sex i.wave##group_num group_num##c.ln_Cystathion c.levo_equivalent##i.wave if filter_Cystathion==0 & filter_Cystein==0, vce(robust)

*2obut-Cystathionine
xtreg ln_OA01 age sex i.wave##group_num group_num##c.ln_Cystathion if filter_Cystathion==0 & filter_OA01==0, vce(robust)
xtreg ln_OA01 age sex i.wave##group_num group_num##c.ln_Cystathion c.levodopa_dosis_mg##i.wave if filter_Cystathion==0 & filter_OA01==0, vce(robust)
xtreg ln_OA01 age sex i.wave##group_num group_num##c.ln_Cystathion c.levo_equivalent##i.wave if filter_Cystathion==0 & filter_OA01==0, vce(robust)

*Aalpha-Cystathionine
xtreg ln_LAlpha age sex i.wave##group_num group_num##c.ln_Cystathion if filter_Cystathion==0 & filter_LAlpha==0, vce(robust)
xtreg ln_LAlpha age sex i.wave##group_num group_num##c.ln_Cystathion c.levodopa_dosis_mg##i.wave if filter_Cystathion==0 & filter_LAlpha==0, vce(robust)
xtreg ln_LAlpha age sex i.wave##group_num group_num##c.ln_Cystathion c.levo_equivalent##i.wave if filter_Cystathion==0 & filter_LAlpha==0, vce(robust)

*Bivariate: Top-hits from screening

*SM(181/250)-SM(181/250)
xtreg ln_SMd181_250 age sex i.wave##group_num group_num##c.ln_SMd181_251 if filter_SMd181_251==0 & filter_SMd181_250==0, vce(robust)
test 1.group_num#c.ln_SMd181_251
display r(p)
xtreg ln_SMd181_250 age sex i.wave##group_num group_num##c.ln_SMd181_251 c.levodopa_dosis_mg##i.wave if filter_SMd181_251==0 & filter_SMd181_250==0, vce(robust)
test 1.group_num#c.ln_SMd181_251
display r(p)
xtreg ln_SMd181_250 age sex i.wave##group_num group_num##c.ln_SMd181_251 c.levo_equivalent##i.wave if filter_SMd181_251==0 & filter_SMd181_250==0, vce(robust)
test 1.group_num#c.ln_SMd181_251
display r(p)

*LPC160-Homocitrulline
xtreg ln_LPC160 age sex i.wave##group_num group_num##c.ln_Homocitrulli if filter_Homocitrulli==0 & filter_LPC160==0, vce(robust)
test 1.group_num#c.ln_Homocitrulli
display r(p)
xtreg ln_LPC160 age sex i.wave##group_num group_num##c.ln_Homocitrulli c.levodopa_dosis_mg##i.wave if filter_Homocitrulli==0 & filter_LPC160==0, vce(robust)
test 1.group_num#c.ln_Homocitrulli
display r(p)
xtreg ln_LPC160 age sex i.wave##group_num group_num##c.ln_Homocitrulli c.levo_equivalent##i.wave if filter_Homocitrulli==0 & filter_LPC160==0, vce(robust)
test 1.group_num#c.ln_Homocitrulli
display r(p)

	
log close

clear

