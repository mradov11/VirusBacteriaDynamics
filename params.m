% Parameters to use

pars.e     = 5e-7;        % conversion efficiency
pars.umax  = 1.2;         % maximal growth rate
pars.Rin   = 400;         % half velocity constant TRY 0.4, 4, 40, 400

pars.ds = 0.04;
pars.de = 0.04;
pars.di = 0.04;
pars.dl = 0.04;

pars.lambda = 2;          % transition rate
pars.beta   = 50;         % burst size
pars.phi    = 3.4e-10;    % adsorption rate
pars.eta0   = 1           % max lysis rate
pars.m      = 1/24;       % viral decay rate
pars.J      = 0           % flow rate (1, 10, or 100)

pars.p      = 1;        % integration probability (changed depending on lytic/lysogenic/temprate)
pars.gamma  = 0;        % induction rate (changed depending on lytic/lysogenic/temprate)

x0 = [
    100;        % R(0)
    1e7;       % S(0)
    1e0;         % E(0)
    1e0;         % I(0)
    0;         % L(0)
    1e5        % V(0)
]
t0 = 0
tf = 24