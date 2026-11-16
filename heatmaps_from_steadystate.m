
%% Generate heatmaps from steady-state parameter sweeps

clearvars;
close all;
clc;

folder = "temp/";
files = {
    'L=0.1L_V=0.1V_Reta=0_SteadyStateSweep_2026-11-01_18-47-03.mat'
    'L=0.2L_Reta=10_SteadyStateSweep_2026-11-01_18-49-54.mat'
    'L=0.1L_V=0.1V_Reta=10_SteadyStateSweep_2026-11-01_18-51-17.mat'
    'V=0.2L_Reta=0_SteadyStateSweep_2026-10-31_17-39-49.mat'
    'L=0.2L_Reta=0_SteadyStateSweep_2026-10-26_16-16-23.mat'
    'V=0.2V_Reta=10_SteadyStateSweep_2026-11-01_18-50-23.mat'
};

% Variables to plot
variableNames = {'R','S','E','I','L','V'};

% Columns corresponding to each variable in SteadyState_temp
variableColumns = 3:8;


% Loop through each file
for f = 1:length(files)

    % Load data
    data = load(folder + files{f});
    SteadyState_temp = data.SteadyState_temp;

    % Extract p and gamma values
    P = unique(SteadyState_temp(:,1));
    Gamma = unique(SteadyState_temp(:,2));

    % Loop through each variable
    for k = 1:length(variableNames)

        % Extract steady-state concentration
        values = SteadyState_temp(:,variableColumns(k));

        % Reshape into p x gamma matrix
        Z = reshape(values, length(P), length(Gamma));
        Z = Z';

        % Logarithmic color scale
        Z(Z <= 0) = NaN;

        % Create a separate figure
        figure('Name', files{f}, 'Color', 'w');

        h = heatmap(P, Gamma, Z);

        % Plot title = variable name
        h.Title = variableNames{k};

        h.XLabel = 'p';
        h.YLabel = '\gamma';

        % Logarithmic color scale
        h.ColorScaling = 'log';
    end
end

%% as pannels
% Loop through each file
for f = 1:length(files)

    % Load data
    data = load(folder + files{f});

    SteadyState_temp = data.SteadyState_temp;

    % Extract p and gamma values
    P = unique(SteadyState_temp(:,1));
    Gamma = unique(SteadyState_temp(:,2));

    % Create figure for this parameter set
    figure('Name', files{f}, 'Color', 'w');

    % Generate one heatmap per variable
    for k = 1:length(variableNames)

        % Extract steady-state concentration
        values = SteadyState_temp(:,variableColumns(k));

        % Reshape into p x gamma matrix
        Z = reshape(values, length(P), length(Gamma));

        % Transpose so rows = gamma and columns = p
        Z = Z';

        % Logarithmic color scale requires positive values
        Z(Z <= 0) = NaN;

        % Create heatmap
        subplot(2,3,k);

        h = heatmap(P, Gamma, Z);

        h.Title = variableNames{k};
        h.XLabel = 'p';
        h.YLabel = '\gamma';

        % Logarithmic color scale
        h.ColorScaling = 'log';

        % Colorbar label
        h.ColorbarVisible = 'on';
    end

    sgtitle(strrep(files{f}, '_', '\_'), ...
            'Interpreter', 'tex');
end