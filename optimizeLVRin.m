params;

Rin_values = [0, 4, 40, 400];

for r = 1:length(Rin_values)

    pars.Rin = Rin_values(r);

    results = parameterSweep(pars, t0, tf, x0);

    % Find max L
    [maxL, idxL] = max(L_store, [], 'omitnan');
    pL     = results.p(idxL);
    gammaL = results.gamma(idxL);

    % Find max V
    [maxV, idxV] = max(V_store, [], 'omitnan');
    pV     = results.p(idxV);
    gammaV = results.gamma(idxV);

    fprintf('\nRin = %g\n', pars.Rin)
    fprintf('Max L at p = %.2f, gamma = %.2e\n', pL, gammaL)
    fprintf('Max V at p = %.2f, gamma = %.2e\n', pV, gammaV)

end