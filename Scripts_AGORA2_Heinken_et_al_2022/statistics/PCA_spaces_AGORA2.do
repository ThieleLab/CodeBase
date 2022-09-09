*Secretion spaces/ Uptake Spaces via PCAs
* Johannes Hertel
*generates raw figures S3A, S3B


clear
clear mata
clear matrix
set more off
set maxvar 32000
capture log close 

cd A:\AGORA_2_New\Files_for_Johannes_revision\processed_data

import excel "A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\AGORA2_infoFile.xlsx", sheet("AGORA2_Reconstructions") firstrow
drop AQ-CO
rename MicrobeID model_id
*drop empty rows
egen test=seq(), by(model_id)
drop if test>1
drop test
save Agora2_infofile.dta, replace
clear

import delimited A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\uptakeFluxes_AGORA2_refined.txt
foreach j of varlist ex_12dgr180e-ex_zn2e{
	rename `j' UP_`j'
	}
save AGORA2_uptakefluxes.dta, replace
clear

import delimited A:\AGORA_2_New\Files_for_Johannes_revision\raw_data_revision\secretionFluxes_AGORA2_refined.txt
save AGORA2_secretionfluxes.dta, replace
merge 1:1 model_id using AGORA2_uptakefluxes.dta
drop _merge
merge 1:1 model_id using AGORA2_infofile.dta
drop _merge

foreach j of varlist Genus-Phylum{
	egen gr_`j'=group(`j')
	tab gr_`j'
	}

gen big_phyla=.
replace big_phyla=1 if Phylum=="Firmicutes" 
replace big_phyla=2 if Phylum=="Proteobacteria"
replace big_phyla=3 if Phylum=="Actinobacteria"
replace big_phyla=4 if Phylum=="Bacteroidetes"
replace big_phyla=5 if Phylum=="Fusobacteria"
replace big_phyla=6 if Phylum=="Tenericutes"
replace big_phyla=7 if Phylum=="Spirochaetes"

*foreach j of varlist ex_12dgr180e-ex_zn2e{
*	gen diff_`j'=UP_`j'-`j'
*	}

foreach j of varlist ex_12dgr180e-ex_zn2e UP_*{
	gen b_`j'=`j'
	replace b_`j'=1 if `j'>0 & `j'!=.
	}
*dropping variables with secretion=0 for all reconstructions
local n=0
foreach j of varlist ex_12dgr180e-ex_zn2e UP_*{
	quietly sum b_`j'
	if r(mean)>0{
		local n=`n'+1
		}
	if r(mean)==0{
		drop `j' b_`j'
		}
	}
display `n'

*foreach j of varlist ex_12dhchole-ex_zn2e{
*	quietly loneway `j' big_phyla
*	display "`j'"
*	display r(rho)
*	}
	

local list_all=""
foreach j of varlist ex_12dgr180e-ex_zn2e{
	local list_all="`list_all'"+ " "+ "`j'"
	}

local list_ind=""	
foreach j of varlist `list_all'{
	if strlen("`list_ind'")==0{
		local list_ind="`j'"
		}
	else{
		local l=0
		quietly corr `j' `list_ind' 
		matrix M=r(C)
		local col=colsof(M)
		forvalues k=2(1)`col'{
			gen test=M[1,`k']
			quietly sum test
			if r(mean)==1{
				local l=1
				}
			*display `l'
			drop test
			}
		if `l'==0{
			local list_ind="`list_ind'"+" "+"`j'"
			}
		}
	}
log using A:\AGORA_2_New\Files_for_Johannes_revision\results\logs\PCA_AGORA2_secretion_uptake_spaces.log, replace
pca `list_ind', component(10)
predict pca1 pca2 pca3 pca4 pca5

cd A:\AGORA_2_New\Files_for_Johannes_revision\results\raw_figures
scatter pca1 pca2 if big_phyla==1 || scatter pca1 pca2 if big_phyla==2 || scatter pca1 pca2 if big_phyla==3 || scatter pca1 pca2 if big_phyla==4 || scatter pca1 pca2 if big_phyla==5 || scatter pca1 pca2 if big_phyla==6 || || scatter pca1 pca2 if big_phyla==7, saving(Figure_S3A_raw, replace)


local list_all_up=""
foreach j of varlist UP_*{
	local list_all_up="`list_all_up'"+ " "+ "`j'"
	}

local list_ind_up=""	
foreach j of varlist `list_all_up'{
	if strlen("`list_ind_up'")==0{
		local list_ind_up="`j'"
		}
	else{
		local l=0
		quietly corr `j' `list_ind_up' 
		matrix M=r(C)
		local col=colsof(M)
		forvalues k=2(1)`col'{
			gen test=M[1,`k']
			quietly sum test
			if r(mean)==1{
				local l=1
				}
			*display `l'
			drop test
			}
		if `l'==0{
			local list_ind_up="`list_ind_up'"+" "+"`j'"
			}
		}
	}

pca `list_ind_up', component(10)
predict pca1_up pca2_up pca3_up pca4_up pca5_up

scatter pca1_up pca2_up if big_phyla==1 || scatter pca1_up pca2_up if big_phyla==2 || scatter pca1_up pca2_up if big_phyla==3 || scatter pca1_up pca2_up if big_phyla==4 || scatter pca1_up pca2_up if big_phyla==5 || scatter pca1 pca2 if big_phyla==6 || || scatter pca1 pca2 if big_phyla==7, saving(Figure_S3B_raw, replace)


log close
clear
