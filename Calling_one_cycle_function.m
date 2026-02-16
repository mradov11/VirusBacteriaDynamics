params;

% create a table to store the data
p_vals     = 0:0.01:1;
gamma_vals = 10.^(-6:-1);
nRuns = length(p_vals) * length(gamma_vals);
% preallocate storage
p_store     = zeros(nRuns,1);
gamma_store = zeros(nRuns,1);
R_store     = zeros(nRuns,1);
S_store     = zeros(nRuns,1);
E_store     = zeros(nRuns,1);
I_store     = zeros(nRuns,1);
L_store     = zeros(nRuns,1);
V_store     = zeros(nRuns,1);
row = 1;

% call the function one_cycle for various values of p and gamma
for index = 0:0.01:1
    p = index;

    for number = -6:1:-1
        gamma = 10^number;

        % run one cycle and capture outputs
        [R,S,E,I,L,V] = one_cycle(p, gamma, pars, t0, tf, x0);

        % store results
        p_store(row)     = p;
        gamma_store(row) = gamma;
        R_store(row)     = R;
        S_store(row)     = S;
        E_store(row)     = E;
        I_store(row)     = I;
        L_store(row)     = L;
        V_store(row)     = V;

        row = row + 1;
    end
end

% save the table to a .mat file
results = table( ...
    p_store, gamma_store, ...
    R_store, S_store, E_store, I_store, L_store, V_store, ...
    'VariableNames', {'p','gamma','R','S','E','I','L','V'} );

save('parameter_sweep_results.mat', 'results');

%%
% plot two ways:
% 3D scatter:
figure
scatter3(results.p, results.gamma, results.V, 40, results.V, 'filled')
set(gca,'YScale','log', 'ZScale', 'log')
xlabel('p')
ylabel('\gamma')
zlabel('V')
title('3D scatter of V')
colorbar

%%
% heatmap:
R_grid = reshape(results.R, length(gamma_vals), length(p_vals));
figure
imagesc(p_vals, gamma_vals, R_grid)
axis xy
set(gca,'YScale','log')
xlabel('p')
ylabel('\gamma')
title('Heatmap of R')
colorbar

%% 10 percent of viruses passed only, 10 cycles
nCycles = 10;

p_vals     = 0:0.01:1;
gamma_vals = 10.^(-6:-1);
nRuns = length(p_vals)*length(gamma_vals);

p_store = zeros(nRuns,1);
gamma_store = zeros(nRuns,1);
V_store = zeros(nRuns,1);
L_store = zeros(nRuns,1);

row = 1;

for p = p_vals
    for gamma = gamma_vals

        x_current = x0;

        for c = 1:nCycles

            [R,S,E,I,L,V] = one_cycle(p, gamma, pars, t0, tf, x_current);

            % 10% viruses passed only
            x_current = [
                100;
                S;
                0;
                0;
                0;
                0.1*V
            ];
        end

        % store final values after 10 cycles
        p_store(row) = p;
        gamma_store(row) = gamma;
        V_store(row) = V;
        L_store(row) = L;

        row = row + 1;

    end
end

results_virusOnly = table(p_store, gamma_store, V_store, L_store);

% Find optimum
[maxV, idxV] = max(V_store, [], 'omitnan');
opt_p_Vvirus = p_store(idxV);
opt_gamma_Vvirus = gamma_store(idxV);

disp('10% Virus Only Passed:')
disp([opt_p_Vvirus opt_gamma_Vvirus maxV])

%% %% 10 percent of lysogens passed only, 10 cycles
nCycles = 10;

p_vals     = 0:0.01:1;
gamma_vals = 10.^(-6:-1);
nRuns = length(p_vals)*length(gamma_vals);

p_store = zeros(nRuns,1);
gamma_store = zeros(nRuns,1);
V_store = zeros(nRuns,1);
L_store = zeros(nRuns,1);

row = 1;

for p = p_vals
    for gamma = gamma_vals

        x_current = x0;

        for c = 1:nCycles

            [R,S,E,I,L,V] = one_cycle(p, gamma, pars, t0, tf, x_current);

            % 10% lysogens passed only
            x_current = [
                100;
                S;
                0;
                0;
                0.1*L;
                0
            ];
        end

        % store final values after 10 cycles
        p_store(row) = p;
        gamma_store(row) = gamma;
        V_store(row) = V;
        L_store(row) = L;

        row = row + 1;

    end
end

results_lysogensOnly = table(p_store, gamma_store, V_store, L_store);

% Find optimum
[maxV, idxV] = max(V_store, [], 'omitnan');
opt_p_Vvirus = p_store(idxV);
opt_gamma_Vvirus = gamma_store(idxV);

disp('10% Lysogens Only Passed:')
disp([opt_p_Vvirus opt_gamma_Vvirus maxV])