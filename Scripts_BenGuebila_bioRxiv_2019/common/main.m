%%
%The following main code allows to generate the data in dHarvey/dWBM
%simulations. 
%%
%Augment Harvey model with organ, geneRxn matrix and Recon metabolite
%vectors for easier analysis
% Requires the COBRA Toolbox https://github.com/opencobra/cobratoolbox/
% Please download WBM/Harvey from vmh.com
load 2017_05_18_HarveyJoint_06_30_constraintHMDB_EUDiet_d;
m = length(modelOrganAllCoupled.mets);
n = length(modelOrganAllCoupled.rxns);
%manually adding a missing gene in the model
modelOrganAllCoupled.genes = [modelOrganAllCoupled.genes;'9481.1'];
g = length(modelOrganAllCoupled.genes);
modelOrganAllCoupled.rxnGeneMat = zeros(n,g);
nCores = 4;
modelOrganAllCoupled=augmentHarvey(modelOrganAllCoupled,nCores);
cd('../data')
save('augHarvey','modelOrganAllCoupled');
%%
%Generate the data
%Simulate the WB ODE model alone in different settings : Healthy, T1D
%and different administration routes and doses : IVGTT,IVITT,SCIB,SCII,
%WBLiquid,WBSolid, Noinf(steady-state)
simODE;
%%
%simulate Healthy and T1D setting in the different adminsitration routes
%(the simulation time is within ~48h, use the saved results instead)
%The simulation is sequential because of time dependency in CRONICS. In FBA
%and pFBA MATLAB reportedly mishandles local global varibales in parfor
%loops, which is why the simulation is set to sequential.
%Healthy
%FBA
system('./batchSimHealthyFBANoOff.sh');
%pFBA
system('./batchSimHealthypFBANoOff.sh');
%CRONICS
system('./batchSimHealthyCRONICSNoOff.sh');
%T1D - without gene expression (disease off target)
%FBA
system('./batchSimT1DFBANoOff.sh');
%pFBA
system('./batchSimT1DpFBANoOff.sh');
%CRONICS
system('./batchSimT1DCRONICSNoOff.sh');
%%
%FVA on Healthy model (Reversible and Irreversible models)
%Irreversible means that reversible reactions are decomposed into two
%reactions with lower bound=0. This is useful to get all fluxes positives
%and perform differential analysis between T1D and healthy.
healthy;
%%
%FVA on unconstrained model (Harvey without constraints from dynamical model)
unc;
%%
%T1D - with gene expression (disease off target)
%First generate the gene expression constraints
geneConversion;
%%
%FBA
system('./batchSimT1DFBAOffGenEx.sh');
%pFBA
system('./batchSimT1DpFBAOffGenEx.sh');
%CRONICS
system('./batchSimT1DCRONICSOffGenEx.sh');
%%
%FVA on T1D model
t1d;