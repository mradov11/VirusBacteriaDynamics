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
    %eta = pars.eta0 * R / (pars.Reta + R);
    eta = pars.eta0;

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
ylim([1e-4 1e10])

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

disp(datatable)

% sort the table by time. look for cases where, for all three strategies, the
% resources last for a long time (14 to 15 hours at least). Look at the plots

datatable.isnan_flag = isnan(datatable.t_R_1e2);
datatable = sortrows(datatable, {'isnan_flag','t_R_1e2'}, {'descend','descend'});
disp(datatable)

%% Make figure 2, which will display the state of the variables over the course of multiple cycles:

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

    x0 = 0.01*x(end, :)     % filter what will go to the next round
    x0(1) = 100;       % reset resources R

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