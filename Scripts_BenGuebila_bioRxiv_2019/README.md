The set of scripts describe the dWBM/dHarvey T1D model that combines WBM/Harvey, the organ-resolved whole-body metabolic constraint-based model and GIM, the gluocse insulin ODE model.

Reference: Pan-organ model integration of metabolic and regulatory processes in type 1 diabetes. https://www.biorxiv.org/content/10.1101/859876v1.abstract

Authors: Marouen Ben Guebila and Ines Thiele.

Corresponding author: Ines Thiele ines.thiele@gmail.com

## System requirements

- Ubuntu 16.04+ / Windows 10

- MATLAB 2013b+ (https://www.mathworks.com/products/matlab.html)

- the COBRA Toolbox v3.0 (https://github.com/opencobra/cobratoolbox/)

## Installation guide

You can add the files in this folder to your MATLAB path.

Note: running the dWBM whole-body metabolic model requires significant computer memory and may requires up to 48h.

## Demonstration

Since the actual simulation of the Whole-body model that covers more than 10k ODE and 80k linear equations, may take several days and consequent resources, you can use the `common\ToyModel.m` script to test the an implementation of the CRONICS framework for coupling ODE and FBA models using Ecoli core model. The script should run in less than a minute. You need to add both the folder and the COBRA Toolbox to the path and simply run the script.

After running the demo script, you should get a two figures printed to the screen. The first figure illustrates the growth rates of E.coli over time using an indirectly coupled ODE-FBA model. The seoncd figure represents the time-course of the metabolite g6p in the cyotplasm using a directly coupled ODE-FBA model.

## Instructions

You can generate the data for the study from scratch by going through the instructions in `common\main.m`.

The repo has the following folders:

### TrialSimulation

 - IVGTT: code and results of the Intravenous Glucose Tolerance Test

 - IVITT: code and results of the Intravenous Insulin Tolerance Test

 - Noinf: code and results of the steady state glucose levels

 - SCIB: code and results of the Subcutaneous Insulin Bolus

 - SCII: code and results of the Subcutaneous Insulin Infusion

 - WBLiquid: code and results of the liquid glucose tolerance tests

 - WBSolid: code and results of the solid meal tolerance test


### common

 - lib: common library functions used for the analysis.

 - FVA: Flux Variability Analysis

 - diffGenes: conversion of T1D differential gene expression into constraints in the metabolic model. The T1D gene expression was taken from (https://www.ncbi.nlm.nih.gov/pubmed/20708640)
