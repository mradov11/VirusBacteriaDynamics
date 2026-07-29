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
    eta = pars.eta0 * R / (pars.Reta + R);
    %eta = pars.eta0;

    % Differential equations:
    % Resource
    dRdt = -pars.e * psi * (S + E + I + L) + pars.J + eta*I*(pars.e)*(0.5);
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
