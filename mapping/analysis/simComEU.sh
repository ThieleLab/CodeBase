#! /bin/bash
source /etc/profile

module use $PROJECTWORK/cplex/soft/modules
module load CPLEX/12.6.3
module load lang/R/3.2.0-ictce-7.3.5-bare

R CMD BATCH simCommunity773_eu.R