%% SEILV Model
pars.e     = 5e-7;        % conversion efficiency
pars.umax  = 1.2;         % maximal growth rate
pars.Rin   = 4.0;         % half velocity constant

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

pars.p      = 0;        % integration probability (changed depending on lytic/lysogenic/temprate)
pars.gamma  = 0;        % induction rate (changed depending on lytic/lysogenic/temprate)

x0 = [
    100;        % R(0)
    1e7;       % S(0)
    1e0;         % E(0)
    1e0;         % I(0)
    0;         % L(0)
    1e5        % V(0)
];
t0 = 0
tf = 24


function dxdt = growth(t, x, pars)
    R = x(1);   % resource
    S = x(2);   % susceptible bacteria
    E = x(3);   % exposed cells
    I = x(4);   % lytically infected cells
    L = x(5);   % lysogens
    V = x(6);   % free virus

    % Growth function ψ(R)
    psi = pars.umax * R / (pars.Rin + R);

    % Lysis function
    eta = pars.eta0 * R / (pars.Rin + R);

    % Differential equations:
    % Resource
    dRdt = -pars.e * psi * (S + E + I + L) + pars.J;
    % Susceptible
    dSdt = psi*S - pars.phi*S*V - pars.ds*S;
    % Exposed
    dEdt = pars.phi*S*V - pars.lambda*E - pars.de*E;
    % Lytically infected
    dIdt = (1 - pars.p)*pars.lambda*E + pars.gamma*L - eta*I - pars.di*I;
    % Lysogens
    dLdt = pars.p*pars.lambda*E + psi*L - pars.gamma*L - pars.dl*L;
    % Virus
    dVdt = pars.beta*eta*I - pars.phi*V*(S + E + I + L) - pars.m*V;

    % Output
    dxdt = [dRdt; dSdt; dEdt; dIdt; dLdt; dVdt];
end

% Make a figure that displays the change in the variables over the course
% of one cycle
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

%%
% Make figure 2, which will display the state of the variables over the 
% course of multiple cycles:

nPassages = 4;     % Number of passages

% Initialize storage
Tall = [];          % Matrix for times
Xall = [];          % Matrix for values at various times

% Initial conditions
x0_curr = x0;
tShift = 0;

% Loop through cycles
for k = 1:nPassages

    [t, x] = ode45(@(t,x) growth(t,x,pars), [t0 tf], x0, opts);

    t = t + tShift;      % Shift time so passages are sequential

    % Add results to the matrices
    Tall = [Tall; t];
    Xall = [Xall; x]; 

    x0 = 0.01*x(end, :)     % filter what will go to the next round
    x0(1) = 100;       % reset resources R

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
