%% Make figure 2, which will display the state of the variables over the course of multiple cycles:
params;

nPassages = 4;     % Number of passages

% Initialize storage
Tall = [];          % Matrix for times
Xall = [];          % Matrix for values at various times

% Initial conditions
x0_curr = x0;
tShift = 0;

% Loop through cycles
for k = 1:nPassages

    [t, x] = one_cycle(pars.p, pars.gamma, pars.Reta, pars, t0, tf, x0_curr);

    t = t + tShift;      % Shift time so passages are sequential

    % Add results to the matrices
    Tall = [Tall; t];
    Xall = [Xall; x]; 

    x0_curr = [
        1e4;        % reset resources R
        1e3;        % reset susceptibles S
        0;
        0;
        0.1*x(end,5);     % 10% lysogens pass to next round
        0.1*x(end,6)      % 10% viruses pass to next round
    ];

    fprintf('Passage %d: final V = %.3e, next V = %.3e\n', ...
    k, x(end,6), x0_curr(6));

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


%% Nest Figure 2 to vary p, gamma, Reta
params;

% Parameter values
pVals     = 0:0.1:1;
gammaVals = 10.^(-6:0.5:-1);
RetaVals  = [0 10 100 1000];

nPassages = 4;

% Create folder to save graphs to
outdir = fullfile('~/Desktop/results multicycle');
mkdir(outdir);

% Loop through parameter combinations
parfor ip = 1:length(pVals)
    for ig = 1:length(gammaVals)
        for ir = 1:length(RetaVals)

            % Set parameters
            pars_local = pars;
            pars_local.p     = pVals(ip);
            pars_local.gamma = gammaVals(ig);
            pars_local.Reta  = RetaVals(ir);

            % Initialize storage
            Tall = [];
            Xall = [];

            % Reset ICs
            x0_curr = x0;
            tShift = 0;

                % Loop through cycles
                for k = 1:nPassages

                    [t, x] = one_cycle(pars_local.p, pars_local.gamma, pars_local.Reta, pars_local, t0, tf, x0_curr);

                    t = t + tShift;      % Shift time so passages are sequential

                    % Add results to the matrices
                    Tall = [Tall; t];
                    Xall = [Xall; x]; 

                    x0_curr = [
                        1e4;        % reset resources R
                        1e3;        % reset susceptibles S
                        0;
                        0;
                        0.1*x(end,5);     % 10% lysogens pass to next round
                        0.1*x(end,6)      % 10% viruses pass to next round
                    ];

                    fprintf('Passage %d: final V = %.3e, next V = %.3e\n', ...
                    k, x(end,6), x0_curr(6));

                    tShift = t(end);   % Update time shift
                end


            % Plot
            fig = figure('Visible','off');
            plot(Tall, Xall, 'LineWidth', 2)
            xlabel('Time (hr)', 'FontSize', 20)
            ylabel('Concentration / Number', 'FontSize', 20)
            legend({'R','S','E','I','L','V'}, 'FontSize', 14)
            set(gca, 'YScale', 'log', 'FontSize', 20)
            hold on
            yline(1e-3,'r:','LineWidth',2)
            title(sprintf(['.1 lysogens, .1 viruses']))

            % Create filename
            filename = sprintf( ...
                'multicycle_R0_%g_S0_%g_V0_%g_p_%0.2f_gamma_%0.0e_Reta_%g.png', ...
                x0(1), x0(2), x0(6), ...
                pars_local.p, pars_local.gamma, pars_local.Reta);
            fullFile = fullfile(outdir, filename);

            % Save figure
            saveas(fig, fullFile);
            close(fig);
        end
    end
end