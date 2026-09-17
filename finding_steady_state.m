function result = find_steady_state(p, gamma, Reta, filter)

    params;

    % covergence rules
    tol = 1e-6;
    maxCycles = 10000;

    % Set parameters
    pars_local = pars;
    pars_local.p = p;
    pars_local.gamma = gamma;
    pars_local.Reta = Reta;

    % Initial starting condition
    x0_curr = x0;

    for k = 1:maxCycles

        % Starting condition of this cycle
        x0_start = x0_curr;

        % Run one cycle
        [~, x] = one_cycle( ...
            p, gamma, Reta, pars_local, t0, tf, x0_start);

        % Final condition of this cycle
        x_final = x(end,:)';

        % Apply between-cycle filter
        switch filter

            case 1
                % Reset R,S; pass 10% V
                x0_curr = [
                    1e4;
                    1e3;
                    0;
                    0;
                    0;
                    0.1*x_final(6)];

            case 2
                % Reset R,S; pass 10% L and 10% V
                x0_curr = [
                    1e4;
                    1e3;
                    0;
                    0;
                    0.1*x_final(5);
                    0.1*x_final(6)];

            case 3
                % Reset R,S; pass 10% L
                x0_curr = [
                    1e4;
                    1e3;
                    0;
                    0;
                    0.1*x_final(5);
                    0];

            otherwise
                error('filter must be 1, 2, or 3');
        end

        % Compare current cycle start to next cycle start
        relChange = max( ...
            abs(x0_curr - x0_start) ./ max(abs(x0_start), 1e-12));

        % Stop if converged
        if relChange < tol
            result = table( ...
                p, gamma, Reta, filter, k, ...
                x0_curr(1), x0_curr(2), x0_curr(3), ...
                x0_curr(4), x0_curr(5), x0_curr(6), ...
                'VariableNames', { ...
                'p', 'gamma', 'Reta', 'filter', 'cycles', ...
                'R0', 'S0', 'E0', 'I0', 'L0', 'V0'});
            return
        end
    end

    % If we reach here, it never converged
    result = table( ...
        p, gamma, Reta, filter, k, ...
        NaN, NaN, NaN, NaN, NaN, NaN, ...
        'VariableNames', { ...
        'p', 'gamma', 'Reta', 'filter', 'cycles', ...
        'R0', 'S0', 'E0', 'I0', 'L0', 'V0'});

end
