%% 10 percent of viruses passed only, 10 cycles
clear all
params;

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

            [~,~,R,S,E,I,L,V] = one_cycle(p, gamma, pars.Reta, pars, t0, tf, x_current);

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

%% 10 percent of lysogens passed only, 10 cycles
clear all
params;

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

            [~,~,R,S,E,I,L,V] = one_cycle(p, gamma, pars.Reta, pars, t0, tf, x_current);

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
[maxL, idxL] = max(L_store, [], 'omitnan');
opt_p_Lvirus = p_store(idxL);
opt_gamma_Lvirus = gamma_store(idxL);

disp('10% Lysogens Only Passed:')
disp([opt_p_Lvirus opt_gamma_Lvirus maxL])