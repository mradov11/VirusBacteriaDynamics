clear all;
params;

Rin_values = [0, 4, 40, 400];

for r = 1:length(Rin_values)

    pars.Rin = Rin_values(r);

    p_vals     = 0:0.01:1;
    gamma_vals = 10.^(-6:1:-1);

    [results, p_vals, gamma_vals] = p_gamma_parameterSweep(pars, t0, tf, x0, p_vals, gamma_vals);

    % Find max L
    [maxL, idxL] = max(results.L, [], 'omitnan');
    pL     = results.p(idxL);
    gammaL = results.gamma(idxL);

    % Find max V
    [maxV, idxV] = max(results.V, [], 'omitnan');
    pV     = results.p(idxV);
    gammaV = results.gamma(idxV);

    fprintf('\nRin = %g\n', pars.Rin)
    fprintf('Max L at p = %.2f, gamma = %.2e\n', pL, gammaL)
    fprintf('Max V at p = %.2f, gamma = %.2e\n', pV, gammaV)

end