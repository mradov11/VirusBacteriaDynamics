params;

% Parameter values
pVals     = [0.25 0.5 0.75 1];
gammaVals = [1e-6 1e-4 1e-2 1e-1];
RetaVals  = [0 10 100 1000];

nPassages = 1;

poolobj = parpool(10)

% Loop through parameter combinations
parfor ip = 1:length(pVals)
    for ig = 1:length(gammaVals)
        for ir = 1:length(RetaVals)

            % Set parameters
            Params = pars;
            Params.p     = pVals(ip);
            Params.gamma = gammaVals(ig);
            Params.Reta  = RetaVals(ir);

            [t, x] = one_cycle(Params.p, Params.gamma, Params, t0, tf, x0);

            % Create filename
            filename = sprintf( ...
                'R0_%g_S0_%g_V0_%g_p_%0.2f_gamma_%0.0e_Reta_%g.mat', ...
                x0(1), x0(2), x0(6), Params.p, Params.gamma, Params.Reta);

            % Save figure
            parsave(filename, t, x)
        end
    end
end
delete(poolobj)
%%


            % Reset ICs
            x0_curr = x0;
            tShift = 0;

            % Passage loop
            for k = 1:nPassages

                [t, x] = one_cycle(Params.p, Params.gamma, Params, t0, tf, x0);

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
