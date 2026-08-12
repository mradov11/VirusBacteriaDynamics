params;

addpath(pwd)

% Create folder to save graphs to
outdir = fullfile('~/Desktop/parfor_results');
mkdir(outdir)

% Parameter values
pVals     = 0:0.1:1;
gammaVals = 10.^(-6:0.5:-1);
RetaVals  = [0, 10, 100, 1000];

%nPassages = 1;

poolobj = parpool(10);

% Loop through parameter combinations
parfor ip = 1:length(pVals)
    for ig = 1:length(gammaVals)
        for ir = 1:length(RetaVals)

            % Set parameters
            Params = pars;
            Params.p     = pVals(ip);
            Params.gamma = gammaVals(ig);
            Params.Reta  = RetaVals(ir);

            [t, x] = one_cycle(Params.p, Params.gamma, Params.Reta, Params, t0, tf, x0);

            % Create filename
            filename = sprintf( ...
                'R0_%g_S0_%g_V0_%g_p_%0.2f_gamma_%0.0e_Reta_%g.mat', ...
                x0(1), x0(2), x0(6), Params.p, Params.gamma, Params.Reta);

            % Save figure
            fullFile = fullfile(outdir, filename);
            parsave(fullFile, t, x)

        end
    end
end
delete(poolobj)