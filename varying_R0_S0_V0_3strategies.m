%% Nest 3 for-loops to interate over values of R0, S0, V0

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
                [t, x, R_final, S_final, E_final, I_final, L_final, V_final] = ...
                    one_cycle(pars.p, pars.gamma, pars.Reta, pars, t0, tf, x0);

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
                    R_final, S_final, E_final, I_final, L_final, V_final, ...
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
