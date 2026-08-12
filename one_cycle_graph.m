function graph = one_cycle_graph(p, gamma, Reta, pars, t0, tf, x0)

[t, x] = one_cycle(p, gamma, Reta, pars, t0, tf, x0);

figure
hold on

colororder([
    0.45 0.45 0.45
    0 0 1
    0.95 0.75 0.2
    0.75 0.3 1.0
    0.4 0.8 1.0
    1 0 0
])

plot(t, x, 'LineWidth', 7)
xlabel('Time (hr)', 'FontSize', 20)
ylabel('Concentration / Number', 'FontSize', 20)
legend({'R','S','E','I','L','V'}, 'FontSize', 14)
set(gca, 'YScale', 'log', 'FontSize', 24)
hold on
yline(1e-3,'r:','LineWidth',2)

end