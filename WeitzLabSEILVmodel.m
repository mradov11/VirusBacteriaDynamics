%% SEILV Model
params;

function dxdt = growth(t, x, pars)
    R = x(1);   % resource
    S = x(2);   % susceptible bacteria
    E = x(3);   % exposed cells
    I = x(4);   % lytically infected cells
    L = x(5);   % lysogens
    V = x(6);   % free virus

    % Growth function ψ(R)
    psi = pars.umax * R / (pars.Rin + R);

    % Lysis function
    eta = pars.eta0 * R / (pars.Reta + R);
    %eta = pars.eta0;

    % Differential equations:
    % Resource
    dRdt = -pars.e * psi * (S + E + I + L) + pars.J + eta*I*(pars.e)*(0.5);
    % Susceptible
    dSdt = psi*S - pars.phi*S*V - pars.ds*S;
    % Exposed
    dEdt = pars.phi*S*V - pars.lambda*E - pars.de*E;
    % Lytically infected
    dIdt = (1 - pars.p)*pars.lambda*E + pars.gamma*L - eta*I - pars.di*I;
    % Lysogens
    dLdt = pars.p*pars.lambda*E + psi*L - pars.gamma*L - pars.dl*L;
    % Virus
    dVdt = pars.beta*eta*I - pars.phi*V*(S + E + I + L) - pars.m*V;

    % Output
    dxdt = [dRdt; dSdt; dEdt; dIdt; dLdt; dVdt];
end

% Make a figure that displays the change in the variables over the course
% of one cycle
opts = odeset('NonNegative', 1:6);
[t, x] = ode45(@(t,x) growth(t,x,pars), [t0 tf], x0, opts);

figure
hold on

colororder([
    0.45 0.45 0.45      % R - slightly darker grey
    0 0 1               % S - blue
    0.95 0.75 0.2       % E - orange/yellow
    0.75 0.3 1.0        % I - brighter purple
    0.4 0.8 1.0         % L - light blue
    1 0 0               % V - red
])

plot(t, x, 'LineWidth', 7)   % thicker lines

xlabel('Time (hr)', 'FontSize', 20)
ylabel('Concentration / Number', 'FontSize', 20)
legend({'R','S','E','I','L','V'}, 'FontSize', 14)
set(gca, 'YScale', 'log', 'FontSize', 24)
ylim([1e-4 1e12])

yline(1e-3,'r:','LineWidth',2)

%% Nest it in 3 for-loops to interate over values of R0, S0, V0

params;

% Strategies: [p, gamma]
strategies = [
    0      0;        % lytic
    0.5    0.05;     % temperate
    1      0      % lysogenic
];

% Initialize results table
datatable = table();

row = 1;  % table row counter

for strat = 1:3
    
    pars.p = strategies(strat, 1);
    pars.gamma = strategies(strat, 2);

    for R0 = 10.^(1:1:5)
        for S0 = 10.^(3:1:8)
            for V0 = 10.^(0:1:7)

                % Initial conditions
                x0 = pars.x0
                x0(1) = R0;       % update R
                x0(2) = S0;       % update S
                x0(6) = V0;       % update V

                % Solve ODE
                opts = odeset('NonNegative', 1:6);
                [t, x] = ode45(@(t,x) growth(t,x,pars), [t0 tf], x0, opts);

                R = x(:,1);

                % Find when R hits 1e-2
                idx = find(R <= 1e-2, 1);

                if isempty(idx)
                    t_hit = NaN;  % never reached threshold
                else
                    t_hit = t(idx);
                end

                % Save results to table
                datatable(row,:) = table( ...
                    strat, R0, S0, V0, ...
                    x(end,1), x(end,2), x(end,3), ...
                    x(end,4), x(end,5), x(end,6), ...
                    t_hit, ...
                    'VariableNames', ...
                    {'Strategy','R0','S0','V0', ...
                     'R_final','S_final','E_final', ...
                     'I_final','L_final','V_final','t_R_1e2'} ...
                );

                row = row + 1;

                % Plot
                f = figure('Visible','off'); % don't spam screen
                semilogy(t, x, 'LineWidth', 2)
                legend({'R','S','E','I','L','V'}, 'Location','best')
                hold on
                yline(1e-2,'r--','LineWidth',2)

                xlabel('Time (hr)')
                ylabel('Concentration / Number')
                title(sprintf('Strat=%d, R0=%.1e, S0=%.1e, V0=%.1e', ...
                    strat, R0, S0, V0))

                % Filename (unique names)
                filename = sprintf( ...
                    'Strat%d_R0_%1.0e_S0_%1.0e_V0_%1.0e.png', ...
                    strat, R0, S0, V0);

                saveas(f, filename)
                close(f)

            end
        end
    end
end
%
disp(datatable)

% sort the table by time. look for cases where, for all three strategies, the
% resources last for a long time (14 to 15 hours at least). Look at the plots

datatable.isnan_flag = isnan(datatable.t_R_1e2);
datatable = sortrows(datatable, {'isnan_flag','t_R_1e2'}, {'descend','descend'});
disp(datatable)


% count how many of each strat lasted >14h
%mask = datatable.t_R_1e2 > 14 | isnan(datatable.t_R_1e2);
%counts = groupcounts(datatable(mask, :), 'Strategy');
%disp(counts)

% Show the strat 3 rows with t>14h
%strat3_table = datatable(datatable.Strategy == 3 & (datatable.t_R_1e2 > 14 | isnan(datatable.t_R_1e2)), :);
%disp(strat3_table)

% Mark all conditions where all 3 strategies had t>14
mask_all = datatable.t_R_1e2 > 14 | isnan(datatable.t_R_1e2);
good_all = datatable(mask_all, :);
% Count how many strategies per (R0,S0,V0)
counts = groupcounts(good_all, {'R0','S0','V0'});
% Keep only those where all 3 strategies passed
valid_conditions = counts(counts.GroupCount == 3, {'R0','S0','V0'});
valid_conditions = innerjoin(datatable, valid_conditions, 'Keys', {'R0','S0','V0'});
disp(valid_conditions)

% final viral density for strategy 1 has to be 10 times the initial viral density
valid_conditions_strat1_subset = valid_conditions(valid_conditions.Strategy == 1 & ...
    valid_conditions.V_final >= 10 .* valid_conditions.V0, :);
disp(valid_conditions_strat1_subset)

%% Nest it inside of loops to vary p and gamma for a set R0, S0, V0, R.in, and 4 R.etas
clear;

params;

% Create folder to save graphs to
outdir = fullfile('~/Desktop/results');
mkdir(outdir)

% Initialize results table
datatable2 = table();


% Set the 3 R.eta values
R_eta_vals = [0, 10, 100, 1000];

datatable_all = table();

row = 1;  % table row counter

for Reta = R_eta_vals
    pars.Reta = Reta;

    for p = 0:0.1:1
        for gamma = 10.^(-6:0.5:-1)

                % Initial conditions
                % Initial conditions
                x0 = pars.x0
                pars.p = p;                 % update p
                pars.gamma = gamma;       % update gamma

                % Solve ODE
                opts = odeset('NonNegative', 1:6);
                [t, x] = ode45(@(t,x) growth(t,x,pars), [t0 tf], x0, opts);

                % Save results to table
                datatable2(row,:) = table(R0, S0, V0, ...
                x(end,1), x(end,2), x(end,3), x(end,4), x(end,5), x(end,6), ...
                p, gamma, Reta, 'VariableNames', ...
                {'R0','S0','V0', 'R_final','S_final','E_final', ...
                'I_final','L_final','V_final', 'p','gamma','Reta'});

                row = row + 1;

                % Plot
                f = figure('Visible','off'); % don't spam screen
                semilogy(t, x, 'LineWidth', 2)
                legend({'R','S','E','I','L','V'}, 'Location','best')
                hold on
                yline(1e-2,'r--','LineWidth',2)

                xlabel('Time (hr)')
                ylabel('Concentration / Number')
                title(sprintf('p=%.2f, gamma=%.1e, R0=%.1e, S0=%.1e, V0=%.1e, R\\eta=%.1e', ...
                p, gamma, R0, S0, V0, pars.Reta))

                filename = fullfile(outdir, sprintf('p=%.2f, gamma=%.1e, R0=%.1e, S0=%.1e, V0=%.1e, Reta=%.1e.png', ...
                p, gamma, R0, S0, V0, pars.Reta));

                saveas(f, filename)
                close(f)
        end
    end
end

disp(datatable2)


%% Making heatmaps:
% Axes values
p_vals = unique(datatable2.p);
gamma_vals = unique(datatable2.gamma);

% Store matrices
Z_all = cell(length(R_eta_vals),1);
for k = 1:length(R_eta_vals)
    Z = zeros(length(gamma_vals), length(p_vals));
    for i = 1:length(p_vals)
        for j = 1:length(gamma_vals)

            idx = datatable2.p == p_vals(i) & ...
                  datatable2.gamma == gamma_vals(j) & ...
                  datatable2.Reta == R_eta_vals(k);

            Z(j,i) = datatable2.L_final(idx);  % can change variable here
        end
    end
    Z_all{k} = Z;
end

% Global color scale
all_vals = cell2mat(Z_all(:));
cmin = min(all_vals(:));
cmax = max(all_vals(:));

% Plot all 3 in one figure, shared scale
figure;
for k = 1:length(R_eta_vals)
    subplot(1,3,k)

    imagesc(p_vals, log10(gamma_vals), Z_all{k});
    set(gca,'YDir','normal', 'ColorScale','log');

    caxis([cmin cmax]);   % shared scale

    xlabel('p');
    ylabel('log_{10}(\gamma)');
    title(sprintf('R\\eta = %.0f', R_eta_vals(k)));
end
h = colorbar;
sgtitle('Comparison across Retas');

% Plot all 3 separately
for k = 1:length(R_eta_vals)
    figure;

    imagesc(p_vals, log10(gamma_vals), Z_all{k});
    set(gca,'YDir','normal');
    colorbar;
    xlabel('p');
    ylabel('log_{10}(\gamma)');
    title(sprintf('R\\eta = %.0f', R_eta_vals(k)));
end


%% Make figure 2, which will display the state of the variables over the course of multiple cycles:
params;

nPassages = 4;     % Number of passages

% Initialize storage
Tall = [];          % Matrix for times
Xall = [];          % Matrix for values at various times

% Initial conditions
x0_curr = x0;
tShift = 0;

% Loop through cycles
for k = 1:nPassages

    [t, x] = ode45(@(t,x) growth(t,x,pars), [t0 tf], x0, opts);

    t = t + tShift;      % Shift time so passages are sequential

    % Add results to the matrices
    Tall = [Tall; t];
    Xall = [Xall; x]; 

    x0(1) = 1e4;       % reset resources R
    x0(2) = 1e3        % reset susceptibles S
    x0(3) = 0
    x0(4) = 0
    x0(5) = 0.1*x(end, 5)     % 10% lysogens pass to next round
    x0(6) = 0


    tShift = t(end);   % Update time shift
end

figure
plot(Tall, Xall, 'LineWidth', 2)
xlabel('Time (hr)', 'FontSize', 20)
ylabel('Concentration / Number', 'FontSize', 20)
legend({'R','S','E','I','L','V'}, 'FontSize', 14)
set(gca, 'YScale', 'log', 'FontSize', 20)
hold on
yline(1e-3,'r:','LineWidth',2)


%% Nest Figure 2 to vary p, gamma, Reta
params;

% Parameter values
pVals     = 0:0.1:1;
gammaVals = 10.^(-6:0.5:-1);
RetaVals  = [0 10 100 1000];

nPassages = 4;

% Create folder to save graphs to
outdir = fullfile('~/Desktop/results multicycle');
mkdir(outdir);

% Loop through parameter combinations
parfor ip = 1:length(pVals)
    for ig = 1:length(gammaVals)
        for ir = 1:length(RetaVals)

            % Set parameters
            pars_local = pars;
            pars_local.p     = pVals(ip);
            pars_local.gamma = gammaVals(ig);
            pars_local.Reta  = RetaVals(ir);

            % Initialize storage
            Tall = [];
            Xall = [];

            % Reset ICs
            x0_curr = x0;
            tShift = 0;

            % Passage loop
            for k = 1:nPassages

                [t, x] = ode45(@(t,x) growth(t,x,pars_local), [t0 tf], x0_curr, opts);

                % Shift time
                t = t + tShift;

                % Store results
                Tall = [Tall; t];
                Xall = [Xall; x];

                % Reset conditions for next passage
                x0_curr(1) = 1e4;
                x0_curr(2) = 1e3;
                x0_curr(3) = 0;
                x0_curr(4) = 0;
                x0_curr(5) = 0.1*x(end,5);
                x0_curr(6) = 0;

                % Update time shift
                tShift = t(end);

            end

            % Plot
            fig = figure('Visible','off');
            plot(Tall, Xall, 'LineWidth', 2)
            xlabel('Time (hr)', 'FontSize', 20)
            ylabel('Concentration / Number', 'FontSize', 20)
            legend({'R','S','E','I','L','V'}, 'FontSize', 14)
            set(gca, 'YScale', 'log', 'FontSize', 20)
            hold on
            yline(1e-3,'r:','LineWidth',2)
            title(sprintf(['R0=%g, S0=%g, V0=%g, p=%.2f, gamma=%.0e, Reta=%g'], ...
                x0(1), x0(2), x0(6), ...
                pars_local.p, pars_local.gamma, pars_local.Reta))

            % Create filename
            filename = sprintf( ...
                'multicycle_R0_%g_S0_%g_V0_%g_p_%0.2f_gamma_%0.0e_Reta_%g.png', ...
                x0(1), x0(2), x0(6), ...
                pars_local.p, pars_local.gamma, pars_local.Reta);
            fullFile = fullfile(outdir, filename);

            % Save figure
            saveas(fig, fullFile);
            close(fig);
        end
    end
end

%% graph of eta for presentation
clear; clc; close all;

eta0 = 1;
R = linspace(0, 50, 1000);
Rin_values = [0 4 20];

% Custom colors (no yellow)
colors = [
    0 0.4470 0.7410;   % blue
    0.8500 0.3250 0.0980; % red/orange
    0.1 0.6 0.1        % dark green
];

figure;
hold on;

for i = 1:length(Rin_values)
    Rin = Rin_values(i);
    eta = (eta0 .* R) ./ (R + Rin);
    plot(R, eta, 'LineWidth', 3, 'Color', colors(i,:));
end

xlabel('R', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('\eta', 'FontSize', 20, 'FontWeight', 'bold');

legend('R_{in} = 0', 'R_{in} = 4', 'R_{in} = 20', ...
       'FontSize', 16, 'Location', 'southeast');

set(gca, 'FontSize', 16, 'LineWidth', 2);
box on;
grid on;

ylim([0 1.1*eta0]);   % Headroom at top
set(gcf, 'Position', [100 100 800 600]);

hold off;