* Data preparation for integrating metabolome data with AGORA2 community models using Yachida et al. 2019 data
* Johannes Hertel


clear
clear mata
clear matrix
set more off
capture log close
set maxvar 50000

cd A:\AGORA_2_New\Files_for_Johannes_revision\processed_data

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\inputDiet_net_secretion_fluxes.csv"
rename netsecretion v1
sort v1
gen len=strlen(v1)
replace v1=substr(v1,1,25) if len>28
egen t=seq(), by(v1)
replace v1=subinstr(v1,"[","_",2)
replace v1=subinstr(v1,"]","_",2)
replace v1=subinstr(v1,"(","_",2)
replace v1=subinstr(v1,")","_",2)
replace v1=subinstr(v1,"-","_",2)
drop len t
rename v1 _varname
xpose, clear varname
order _varname
rename _varname id
drop if substr(strreverse(id),1,2)=="1_"
save netsecretion_CRC_AGORA2.dta,replace
clear

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\inputDiet_net_uptake_fluxes.csv"
rename netuptake v1
sort v1
gen len=strlen(v1)
replace v1=substr(v1,1,25) if len>28
egen t=seq(), by(v1)
replace v1=subinstr(v1,"[","_",2)
replace v1=subinstr(v1,"]","_",2)
replace v1=subinstr(v1,"(","_",2)
replace v1=subinstr(v1,")","_",2)
replace v1=subinstr(v1,"-","_",2)
drop len t
replace v1=subinstr(v1,"EX","UP",1)
rename v1 _varname
xpose, clear varname
order _varname
rename _varname id
drop if substr(strreverse(id),1,2)=="1_"
save netuptake_CRC_AGORA2.dta,replace
clear

import delimited "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\normCoverage_CRC.csv"
replace id=substr(id,4,.)
gen len=strlen(id)
replace id=substr(id,1,29) if len>29
egen t=seq(), by(id)
tostring(t), gen(ts)
replace id=id+ts if t>1
drop t ts len
rename id _varname 
xpose, clear varname
order _varname
rename _var id

merge 1:1 id using netuptake_CRC_AGORA2.dta
keep if _merge==3
drop _merge
merge 1:1 id using netsecretion_CRC_AGORA2.dta
keep if _merge==3
drop _merge


rename id model_id
replace model_id=subinstr(model_id, "le","le_",1)
merge 1:1 model_id using "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\metabolites_kegg.dta" 
drop if _merge==2
drop _merge
merge 1:1 model_id using "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\meta_data.dta" 
keep if _merge==3
drop _merge

save CRC_AGORA2_merged.dta, replace
clear

