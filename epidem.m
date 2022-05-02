%function describing the dynamics of the controlled SIDARE model
function [y,dy] = epidem(dt, x, beta, u, v, gamma_i, gamma_d, gamma_a, ksi_i, ksi_d, mu, zeta, psi, psi_hat)
    
    %Controlled SIDARE model
    dy(1,1) = -beta*(1 - u)*x(1,1)*x(2,1) + psi*x(7,1) + psi_hat*x(5,1) - zeta*x(1,1); %Susceptible State
    dy(2,1) = beta*(1 - u)*x(1,1)*x(2,1) - gamma_i*x(2,1) - ksi_i*x(2,1) - v*x(2,1); %Infected undetected State
    dy(3,1) = v*x(2,1) - gamma_d*x(3,1) - ksi_d*x(3,1); %Detected infected State
    dy(4,1) = ksi_i*x(2,1) +  ksi_d*x(3,1) - gamma_a*x(4,1) - mu*x(4,1); %Acutely symptomatic State
    dy(5,1) = gamma_i*x(2,1) +  gamma_d*x(3,1) + gamma_a*x(4,1) + psi_hat*x(5,1); %Recovered State
    dy(6,1) = mu*x(4,1); %Extinct (Deceased) State
    dy(7,1) = zeta*x(1,1) - psi*x(7,1); %Vaccinated State
    y = x + dt*dy; %State update



