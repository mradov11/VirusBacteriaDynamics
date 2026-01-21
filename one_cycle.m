function graph = one_cycle(p, gamma, pars, t0, tf, x0)

pars.p = p
pars.gamma = gamma

opts = odeset('NonNegative', 1:6);
[t, x] = ode45(@(t,x) growth(t,x,pars), [t0 tf], x0, opts);

figure
plot(t, x, 'LineWidth', 2)
xlabel('Time (hr)', 'FontSize', 20)
ylabel('Concentration / Number', 'FontSize', 20)
legend({'R','S','E','I','L','V'}, 'FontSize', 14)
set(gca, 'YScale', 'log', 'FontSize', 20)
hold on
yline(1e-3,'r:','LineWidth',2)

end