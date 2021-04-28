*IBD analyses

clear
clear mata
clear matrix
set more off
set maxvar 32000
capture log close

cd A:\IBD_project\IBD_simulation_data\

import delimited A:\IBD_project\IBD_simulation_data\netProductionFluxes_IBD.csv, clear
gen _varname="__"+metabolite_vmh_id
xpose, clear varname
drop in 1/2
rename _varname id
save netProductionFluxes_IBD.dta, replace
clear


import delimited A:\IBD_project\IBD_simulation_data\MicrobeContributionsToMetabolites_IBD.csv
gen _varname="_"+substr(strain,-10,.)+"_"+metabolite_vmh_id
drop v1-species
xpose, clear varname
rename _varname id
save MicrobeContributionsToMetabolites_IBD.dta, replace
clear

import delimited A:\IBD_project\IBD_simulation_data\metadata_IBD.csv, varnames(1) 
gen id="srr"+substr(ind,4,.)
egen group=group(stratification)
save metadata_IBD.dta, replace
clear

import delimited A:\IBD_project\IBD_simulation_data\IBD_species_data.csv, varnames(1) 
gen _varname="_"+substr(x,-10,.)
drop v1 x
xpose, clear varname
rename _varname id
save IBD_species_data.dta, replace
clear

use metadata_IBD.dta
merge 1:1 id using netProductionFluxes_IBD.dta
drop _merge
merge 1:1 id using MicrobeContributionsToMetabolites_IBD.dta
drop _merge
merge 1:1 id using IBD_species_data.dta
drop _merge

save data_merge_all.dta,replace
clear
