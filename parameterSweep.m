function results = parameterSweep(pars, t0, tf, x0)

% create a table to store the data
p_vals     = 0:0.01:1;
gamma_vals = 10.^(-6:1:-1);
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
for p = p_vals

    for gamma = gamma_vals;

        % run one cycle and capture outputs
        [~,~,R,S,E,I,L,V] = one_cycle(p, gamma, pars, t0, tf, x0);

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

end