#! /bin/bash
source /etc/profile

module load lang/R/3.2.0-ictce-7.3.5-bare

R CMD BATCH extract_coverage.R