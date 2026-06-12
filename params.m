% Parameters to use

pars.e     = 5e-7;        % conversion efficiency
pars.umax  = 1.2;         % maximal growth rate
pars.Rin   = 1000;         % half velocity constant for host growth
pars.Reta  = 100;         % half velocity constant for virus lysis (try 1, 10, 100, 1000, 10000)

pars.ds = 0.04;
pars.de = 0.04;
pars.di = 0.04;
pars.dl = 0.04;

pars.lambda = 2;          % transition rate
pars.beta   = 50;         % burst size
pars.phi    = 3.4e-10;    % adsorption rate
pars.eta0   = 1;           % max lysis rate
pars.m      = 1/24;       % viral decay rate
pars.J      = 0;           % flow rate (1, 10, or 100)

%Lytic
pars.p      = 0;        % integration probability (changed depending on lytic/lysogenic/temprate)
pars.gamma  = 0;        % induction rate (changed depending on lytic/lysogenic/temprate)

%Temprate
%pars.p      = 0.5;        % integration probability (changed depending on lytic/lysogenic/temprate)
%pars.gamma  = 0.05;        % induction rate (changed depending on lytic/lysogenic/temprate)

R0 = 1e4
S0 = 1e3
V0 = 1

pars.x0 = [
    R0;       % R(0)
    S0;       % S(0)
    0;       % E(0)
    0;       % I(0)
    0;         % L(0)
    V0;        % V(0)
]

x0 = pars.x0;

t0 = 0;
tf = 24;