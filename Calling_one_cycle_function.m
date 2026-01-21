params;

% call the function one_cycle for various values of p and gamma
for index = 0:0.25:1
    p = index;
    gamma = 0.05;    % can be set to any function of index
    one_cycle(p, gamma, pars, t0, tf, x0)
end