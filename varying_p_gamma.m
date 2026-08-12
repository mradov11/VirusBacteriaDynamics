params;

% Run the parameter sweep using the function
p_vals     = 0:0.01:1;
gamma_vals = 10.^(-6:1:-1);

[results, p_vals, gamma_vals] = p_gamma_parameterSweep(pars, t0, tf, x0, p_vals, gamma_vals);

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