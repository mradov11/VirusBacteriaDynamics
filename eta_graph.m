%% graph of eta for presentation
clear; clc; close all;

eta0 = 1;
R = linspace(0, 50, 1000);
Rin_values = [0 4 20];

% Custom colors (no yellow)
colors = [
    0 0.4470 0.7410;   % blue
    0.8500 0.3250 0.0980; % red/orange
    0.1 0.6 0.1        % dark green
];

figure;
hold on;

for i = 1:length(Rin_values)
    Rin = Rin_values(i);
    eta = (eta0 .* R) ./ (R + Rin);
    plot(R, eta, 'LineWidth', 3, 'Color', colors(i,:));
end

xlabel('R', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('\eta', 'FontSize', 20, 'FontWeight', 'bold');

legend('R_{in} = 0', 'R_{in} = 4', 'R_{in} = 20', ...
       'FontSize', 16, 'Location', 'southeast');

set(gca, 'FontSize', 16, 'LineWidth', 2);
box on;
grid on;

ylim([0 1.1*eta0]);   % Headroom at top
set(gcf, 'Position', [100 100 800 600]);

hold off;