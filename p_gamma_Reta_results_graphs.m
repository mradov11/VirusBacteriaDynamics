%% Nest loops to vary p and gamma for a set R0, S0, V0, R.in, and 4 R.etas
params;

% Create folder to save graphs to
outdir = fullfile('~/Desktop/results');
mkdir(outdir)

% Initialize results table
datatable2 = table();

% Set the parameter values
R_eta_vals = [0, 10, 100, 1000];
p_vals = 0:0.1:1;
gamma_vals = 10.^(-6:0.5:-1);

row = 1;  % table row counter

for Reta = R_eta_vals
    pars.Reta = Reta;

    for p = p_vals
        for gamma = gamma_vals

                % Initial conditions
                % Initial conditions
                x0 = pars.x0;
                pars.p = p;                 % update p
                pars.gamma = gamma;       % update gamma

                % Solve ODE
                [t, x, R_final, S_final, E_final, I_final, L_final, V_final] = ...
                    one_cycle(p, gamma, Reta, pars, t0, tf, x0);

                % Save results to table
                datatable2(row,:) = table(R0, S0, V0, ...
                R_final, S_final, E_final, I_final, L_final, V_final, ...
                p, gamma, Reta, 'VariableNames', ...
                {'R0','S0','V0', 'R_final','S_final','E_final', ...
                'I_final','L_final','V_final', 'p','gamma','Reta'});

                row = row + 1;

                % Plot
                f = figure('Visible','off'); % don't spam screen
                semilogy(t, x, 'LineWidth', 2)
                legend({'R','S','E','I','L','V'}, 'Location','best')
                hold on
                yline(1e-2,'r--','LineWidth',2)

                xlabel('Time (hr)')
                ylabel('Concentration / Number')
                title(sprintf('p=%.2f, gamma=%.1e, R0=%.1e, S0=%.1e, V0=%.1e, R\\eta=%.1e', ...
                p, gamma, R0, S0, V0, pars.Reta))

                filename = fullfile(outdir, sprintf('p=%.2f, gamma=%.1e, R0=%.1e, S0=%.1e, V0=%.1e, Reta=%.1e.png', ...
                p, gamma, R0, S0, V0, pars.Reta));

                saveas(f, filename)
                close(f)
        end
    end
end

disp(datatable2)


%% Making heatmaps:
% Axes values
p_vals = unique(datatable2.p);
gamma_vals = unique(datatable2.gamma);

% Store matrices
Z_all = cell(length(R_eta_vals),1);
for k = 1:length(R_eta_vals)
    Z = zeros(length(gamma_vals), length(p_vals));
    for i = 1:length(p_vals)
        for j = 1:length(gamma_vals)

            idx = datatable2.p == p_vals(i) & ...
                  datatable2.gamma == gamma_vals(j) & ...
                  datatable2.Reta == R_eta_vals(k);

            Z(j,i) = (datatable2.V_final(idx) + datatable2.L_final(idx));  % can change variable here
        end
    end
    Z_all{k} = Z;
end

% Global color scale
all_vals = cell2mat(Z_all(:));
cmin = min(all_vals(:));
cmax = max(all_vals(:));

% Plot all 3 in one figure, shared scale
figure;
for k = 1:length(R_eta_vals)
    subplot(1,length(R_eta_vals),k)

    imagesc(p_vals, log10(gamma_vals), Z_all{k});
    set(gca,'YDir','normal', 'ColorScale','log');

    caxis([cmin cmax]);   % shared scale

    xlabel('p');
    ylabel('log_{10}(\gamma)');
    title(sprintf('R\\eta = %.0f', R_eta_vals(k)));
end
h = colorbar;
sgtitle('Comparison across Retas');

% Plot all 3 separately
for k = 1:length(R_eta_vals)
    figure;

    imagesc(p_vals, log10(gamma_vals), Z_all{k});
    set(gca,'YDir','normal');
    colorbar;
    xlabel('p');
    ylabel('log_{10}(\gamma)');
    title(sprintf('R\\eta = %.0f', R_eta_vals(k)));
end