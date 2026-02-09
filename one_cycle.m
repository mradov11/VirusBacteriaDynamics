function [R, S, E, I, L, V] = one_cycle(p, gamma, pars, t0, tf, x0)

    pars.p = p;
    pars.gamma = gamma;

    opts = odeset('NonNegative', 1:6);
    [t, x] = ode45(@(t,x) growth(t,x,pars), [t0 tf], x0, opts);

    % extract final values
    R = x(end,1);
    S = x(end,2);
    E = x(end,3);
    I = x(end,4);
    L = x(end,5);
    V = x(end,6);

end