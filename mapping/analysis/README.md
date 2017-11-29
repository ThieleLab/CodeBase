# Analysis scripts
This repository contains scripts to simulate individual microbiotas and analyze data

## Pre-processing
The script *compile_diet.R* compiles the information of different diets from [vmh](http://vmh.uni.lu/#nutrition). 
The script *compile_table.R* compiles the abundance data from mapped coverages. 
The script *model_generation.R* generates microbial community models using [BacArena](https://cran.r-project.org/web/packages/BacArena/index.html).

## Simulation
The script *reInitProb.R* reinitiates LP problem objects for warmstart.
The script *simCommunity773_eu.R* simulates microbial community models.
The script *simComEU.sh* is a launcher script to start the simulation on a Linux machine or HPC.

## Analysis
The script *patient_analysis_prep.R* prepares and extracts simulation objects for subsequent analysis.
The script *patient_analysis_figures.R* analyzes prepared simulation results and plots figures.