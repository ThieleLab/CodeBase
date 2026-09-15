function personalizedModel = personalizeWholeBodyModel(model,Diet,Physiological_param)
% This function computes personalized whole body model with physiolgical parameters based on the
% provided input data.
%
% function personalizedModel = personalizeWholeBodyModel(model,Diet,inputData)
%
% INPUT
% model                     model structure, whole-body metabolic model
% IndividualParameters      structure with Individual parameters to be
%                           personalized or updated based on the input
%                           data. This structure, with default parameters, can be obtained using the
%                           script standardPhysiolDefaultParameters.m
% InputData                 InputData
% optionCardiacOutput       Different ways of calculatibg the cardiac
%                           output based on physiological data have been implemented.
%                           - Based on heart rate and stroke volume: CardiacOutput = HeartRate * StrokeVolume. (optionCardiacOutput = 1, default)
%                           - Assume that cardiac output = blood volume.(optionCardiacOutput = 2)
%                           - Based on the polynomial suggested by Youndg et al [http://www.ams.sunysb.edu/~hahn/psfile/pap_obesity.pdf]: CardiacOutput = 9119-exp(9.164-2.91e-2*Wt+3.91e-4*Wt^2-1.91e-6*Wt^3); Wt = weight in kg (optionCardiacOutput = 0)
%                           - Based on Fick's principle, which requires that the oxygen consumption rate is known: VO_2 = (CO* C_a) - (CO *C_v);  where CO = Cardiac Output, Ca = Oxygen concentration of arterial blood and Cv = Oxygen concentration of mixed venous blood. We assume that C_a is 200 ml O2/L and C_v is 150 ml O2/L. (optionCardiacOutput = 3)
%                           - Based on Fick's principle, while estimating the VO_2 based on the body surface area (BSA):   We use the Du Bois formula,[4][5]: BSA=0.007184 * Wt^0.425* Ht^0.725; Wt = weight in kg; Ht in cm; (optionCardiacOutput = 4)
%
% OUTPUT
% modelPersonalized         Updated model structure
% IndividualParametersNew   Updated, personalized individual parameters
%
% Faiz Khan Mohammad, 2023

    % Step 1: Load or create your whole-body model (e.g., biomechanical model)
    global useSolveCobraLPCPLEX
    if (Physiological_param.gender == 0)
        sex ='male'; 
    else 
        sex = 'female'; 
    end
    
    standardPhysiolDefaultParameters;
    IndividualParameters_model = IndividualParameters;
    % Step 2: Preprocess and align the input data if necessary
    InputData =  {'sex' ''  sex
                  'age' '' num2str(Physiological_param.age)
                  'Height' '' num2str(Physiological_param.height)
                  'Weight' '' num2str(Physiological_param.bodyMass)
                  'Fat-free' '' num2str((100-Physiological_param.per_bodyFat)*Physiological_param.bodyMass/100)
                  'BMI' '' num2str(Physiological_param.BMI) 
                  'BMR' '' num2str(Physiological_param.BMR(1))
                  };
    model.sex = sex;
      
    % Step 3: Personalize the model using the input data
    model = physiologicalConstraintsHMDBbased(model,IndividualParameters_model);
    model = setDietConstraints(model, Diet); % Diet
%     model = setSimulationConstraints(model);
    model.osense = -1;
    
    minInf = -1000000;
    maxInf = 1000000;
    % count how many reactions have a non-infinity bounds
    minConstraints = length(intersect(find(model.lb>minInf),find(model.lb)));
    maxConstraints =length(intersect(find(model.ub<maxInf),find(model.ub)));
    
    model.lb(strmatch('BBB_KYNATE[CSF]upt',model.rxns)) = -1000000; %constrained uptake
    model.lb(strmatch('BBB_LKYNR[CSF]upt',model.rxns)) = -1000000; %constrained uptake
    model.lb(strmatch('BBB_TRP_L[CSF]upt',model.rxns)) = -1000000; %constrained uptake
    model.ub(strmatch('Brain_EX_glc_D(',model.rxns)) = -100; % currently -400 rendering many of the models to be infeasible in germfree state
    
    
    
    % load standard parameters
    
    % adjust Muscle and fat atp in biomass - compute and compare BMF's with
    % BMRs
    if useSolveCobraLPCPLEX
        model.A = model.S;
    end
    optionCardiacOutput = 0; %0 4 2
    [model,IndividualParametersPersonalized] = individualizedLabReport(model,IndividualParameters, InputData,optionCardiacOutput);
    % 3. set personalized constraints
    model.IndividualParametersPersonalized=IndividualParametersPersonalized;
    % 3. set personalized constraints

    % calculate new organ weight fraction for a given personalized weight
    [listOrgan,OrganWeight,OrganWeightFract,IndividualParametersPersonalized] = calcOrganFract(model,IndividualParametersPersonalized);

    % I repeat this step, as I calculate in some options the CO based on
    % bloodVol which is only calculated in calcOrganFract
    [model,IndividualParametersPersonalized] = individualizedLabReport(model,IndividualParametersPersonalized, InputData,optionCardiacOutput);
    model.IndividualParametersPersonalized = IndividualParametersPersonalized;
    % adjust whole body maintenance reaction based on new organ weight
    % fractions
    % adjust adipocyte and muscle weight to measured one
    Fat_fraction = 1 - str2num(InputData{5,3})/str2num(InputData{4,3});
    FatDiff = Fat_fraction-OrganWeightFract(11);
    OrganWeightFract(11) = OrganWeightFract(11)+FatDiff;
    OrganWeightFract(12) = OrganWeightFract(12)-FatDiff;

    Fat_weight =  (str2num(InputData{4,3})-str2num(InputData{5,3}))*1000;
    Fat_weightOri=str2num(IndividualParametersPersonalized.OrgansWeights{11,2});
    Muscle_weightOri=str2num(IndividualParametersPersonalized.OrgansWeights{12,2});
    % substract higher fat_weight from Muscle_weight
    Muscle_weight = Muscle_weightOri-(Fat_weight - Fat_weightOri);

    IndividualParametersPersonalized.OrgansWeights{11,2} = num2str(Fat_weight);
    IndividualParametersPersonalized.OrgansWeights{11,3} = num2str(OrganWeightFract(11)*100);
    IndividualParametersPersonalized.OrgansWeights{12,2} = num2str(Muscle_weight);
    IndividualParametersPersonalized.OrgansWeights{12,3} = num2str(OrganWeightFract(12)*100);

    [model] = adjustWholeBodyRxnCoeff(model, listOrgan, OrganWeightFract);
    model.listOrgan = listOrgan;
    model.OrganWeightFract = OrganWeightFract;
    model = setDietConstraints(model);
    % set some more constraints
    model = physiologicalConstraintsHMDBbased(model,IndividualParametersPersonalized);
    BM = (find(~cellfun(@isempty,strfind(model.rxns,'Heart_DM_atp'))));
    if strcmp(sex,'male')
        model.lb(BM) = 6000;
    elseif strcmp(sex,'female')
        model.lb(BM) = 1000;
    end
    % Step 4: Return the personalized model
    personalizedModel = model;
end