%This script creates the variables through which the required parameters 
%and files are inputted to the metagenomic pipeline (MgPipe). Input 
%variables should be changed by the user according to what specified in the 
%documentation. Running this script will automatically launch the pipeline. 

% Federico Baldini, 2017-2018
%HPC adding files and testing 
parpool(7)
addpath(genpath('/mnt/gaiagpfs/users/workdirs/fbaldini/cobratoolbox'))
addpath(genpath('/mnt/gaiagpfs/apps/resif/data/production/v0.3-20170713/default/software/math/CPLEX/12.7.1-foss-2017a/cplex/matlab/x86-64_linux'))
initCobraToolbox()
changeCobraSolver('ibm_cplex')
%Testing HPC simulation functions
load ecoli_core_model.mat
val0=optimizeCbModel(model)
val=solveCobraLP(model)
[min,max]=fastFVA(model)
%Clearing workspace
clear('ans','max','min','model','val')
changeCobraSolver('ibm_cplex')

%REQUIRED INPUT VARIABLES
% path to microbe models
modPath='/mnt/gaiagpfs/users/workdirs/fbaldini/Models0518/';
% path where to save results
resPath='/mnt/gaiagpfs/users/workdirs/fbaldini/PDbatch/batch1/results';
% path to where the COBRA Toolbox is located
global CBTDIR
toolboxPath=CBTDIR;
% path to and name of the file with dietary information.
dietFilePath=[CBTDIR filesep 'papers' filesep '2018_microbiomeModelingToolbox' filesep 'resources' filesep 'AverageEuropeanDiet'];
% path to and name of the file with abundance information.
abunFilePath='/mnt/gaiagpfs/users/workdirs/fbaldini/PDbatch/batch1/normCoverage.csv';
% name of objective function of organisms 
objre={'EX_biomass(e)'};
%the output is vectorized picture, change to '-dpng' for .png
figForm = '-depsc'
% number of cores dedicated for parallelization 
numWorkers = 7;
% autofix for names mismatch
autoFix = 1; 
% if outputs in open formats should be produced for each section (1=T)
compMod = 0; 
% if documentation (.csv) on stratification criteria is available
indInfoFilePath='none';
% to enable also rich diet simulations 
rDiet = 0; 
% if if to use an external solver and save models with diet
extSolve = 0; 
% the type of FVA function to use to solve
fvaType = 1; 
% To tourn off the autorun to be able to manually execute each part of the pipeline.
autorun=1; 

%END OF REQUIRED INPUT VARIABLES

%%
%PIPELINE LAUNCHER 
[init,modPath,toolboxPath,resPath,dietFilePath,abunFilePath,indInfoFilePath,objre,figForm,numWorkers,autoFix,compMod,rDiet,extSolve,fvaType,autorun]= initMgPipe(modPath, toolboxPath, resPath, dietFilePath, abunFilePath, indInfoFilePath, objre, figForm, numWorkers, autoFix, compMod, rDiet,extSolve,fvaType,autorun);

