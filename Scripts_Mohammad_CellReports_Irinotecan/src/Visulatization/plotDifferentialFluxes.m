function plotDifferentialFluxes(Flux_t, t, solution, model, gender)
% Generates and plots differential fluxes over time for various organs
%
% Computes differential fluxes between a drug-treated and a reference state,
% calculates the area under the curve (AUC) for these fluxes, and visualizes
% the results across different organs.
%
% USAGE:
%
%   differentialFluxesPlot(Flux_t, t, solution, model, gender)
%
% INPUT:
%    Flux_t:             (n x m matrix) Flux profile over time
%    t:                  (1 x m vector) Time vector
%    solution:          (struct) Reference state (No Drug state)
%                         * v - (n x 1 vector) Fluxes in the reference state
%    model:              (struct) Metabolic model structure with the following fields
%                         * rxns - (n x 1 cell array) Reaction names
%                         * mets - (p x 1 cell array) Metabolite names
%                         * S    - (p x n matrix) Stoichiometric matrix
%    gender:             (char) Gender of the subject ('male' or 'female')
%
% OUTPUT:
%    The function generates a figure displaying the differential fluxes for various organs.
%
% EXAMPLE:
%    differentialFluxesPlot(Flux_t, t, solution, model, 'male')
%
% INTERNAL FUNCTIONS:
%    findRxnIDs:         Finds reaction IDs for specified reactions
%    findRxnsFromMets:   Finds reactions involving specified metabolites
%
% NOTE:
%    Ensure that the provided model structure contains the required fields.
%
% Authors: % August 2024 - Faiz Khan Mohammad & Ines Thiele

    cum_flux_t = cumtrapz(t/24, Flux_t, 2);

    % Find differential fluxes (Delta_vi) of all reactions
    Delta_vi = (Flux_t - solution.v);
    % Set small differential flux values to zero for numerical stability
    Delta_vi(abs(Delta_vi) < 1e-5) = 0;

    % Extract compartments from reaction names
    compartments = strtok(model.rxns, '_');

    % Organ abbreviations
    if strcmpi(gender, 'male')
        Organ_abbrv = {'Adipocytes';'Agland';'Bcells';'Brain';'CD4Tcells';'Colon';'Gall';'Heart';'Kidney';'Liver';'Lung';
            'Monocyte';'Muscle';'Nkcells';'Pancreas';'Platelet';'Prostate';'Pthyroidgland';'RBC';'Retina';'Scord';'sIEC';'Skin';'Spleen';
            'Stomach';'Testis';'Thyroidgland';'Urinarybladder'};
    elseif strcmpi(gender, 'female')
        Organ_abbrv = {'Adipocytes';'Agland';'Bcells';'Brain';'Breast';'CD4Tcells';'Cervix';'Colon';'Gall';'Heart';'Kidney';'Liver';'Lung';
            'Monocyte';'Muscle';'Nkcells';'Ovary';'Pancreas';'Platelet';'Pthyroidgland';'RBC';'Retina';'Scord';'sIEC';'Skin';'Spleen';
            'Stomach';'Thyroidgland';'Urinarybladder';'Uterus'};
    end

    % Calculate area under the curve (AUC) for differential fluxes
    AUC = cumtrapz(t, Delta_vi, 2);

    % Identify reaction IDs for specific metabolites
    R = findRxnIDs(model, findRxnsFromMets(model, model.mets(contains(model.mets, {'cpt11[','sn38[', 'sn38g[', 'apc[', ...
        'npc[', 'apcald[', 'apca[', 'sn38cb[', 'sn38gcb[', 'cpt11cb[', 'bpca[', 'etiri[', 'tetrapeg[', 'lvlnar['}))));
    
    % Create a table of reaction names and differential fluxes
    T = table(model.rxns(R), Delta_vi(R, :));
    % Filter out reactions with zero differential flux
    T = T(any(T{:, 2:end}, 2), :);
    Tab = T;

    % Create a table of reaction names, AUC, and differential fluxes
    T = table(model.rxns(R), AUC(R, end), Delta_vi(R, :));
    Tab = T(any(T{:, 2:end}, 2), :);
    % Sort the table rows by reaction name
    Tab = sortrows(Tab, 'Var1', 'ascend');

    % Plot differential fluxes for each organ
    fig = figure('Position', [10, 10, 2400, 1200]);
    for i = 1:length(Organ_abbrv)
        subplot(5, 6, i);
        % Find indices of reactions belonging to the current organ
        IDs = find(ismember(compartments, Organ_abbrv{i}));
        % Plot differential fluxes over time for the current organ
        plot(t, Delta_vi(IDs, :));
        % Find the maximum absolute value in the y data
        maxAbsY = max(max(abs(Delta_vi(IDs, :)), [], 2));
        % Set y-axis limits symmetrically around 0
        ylim([-maxAbsY, maxAbsY]);
        title(Organ_abbrv{i});
        set(gca, 'FontSize', 15);
        set(gca, 'LineWidth', 1.25);
        set(get(gca, 'XAxis'), 'TickLength', [0 0]);
        set(get(gca, 'YAxis'), 'TickLength', [0 0]);
        % Set specific y-axis limits for certain organs
        if ismember(Organ_abbrv{i}, {'Agland','Bcells','Breast','Cervix','Monocyte','Nkcells','Ovary','Retina','Platelet', ...
                'Pthyroidgland','Thyroidgland','Urinarybladder','Gall','Prostate','Testis'})
            ylim([-0.5 0.5]);
        end
        xlim([0 40]);
    end

    % Add common labels and title to the figure
    han = axes(fig, 'Visible', 'off'); 
    han.Title.Visible = 'on';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabelText = ylabel(han, '\Delta flux (\muM per day)', 'FontSize', 20);
    ylabelPosition = get(ylabelText, 'Position');
    ylabelPosition(1) = ylabelPosition(1) - 0.02; % Adjust the ylabel position left
    set(ylabelText, 'Position', ylabelPosition);
    
    xlabelText = xlabel(han, 'time (hrs)', 'FontSize', 20);
    xlabelPosition = get(xlabelText, 'Position');
    xlabelPosition(2) = xlabelPosition(2) - 0; % Adjust the xlabel position
    set(xlabelText, 'Position', xlabelPosition);
    
    titleText = title(han, 'Differential Fluxes', 'FontSize', 25);
    titlePosition = get(titleText, 'Position');
    titlePosition(2) = titlePosition(2) + 0.02; % Adjust the title position up
    set(titleText, 'Position', titlePosition);
end