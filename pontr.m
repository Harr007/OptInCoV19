%function of costate variables update based on pontryagin's minimum principle

function [y,dy] = pontr(dt, l, x, u, v, zeta, beta, gamma_i, gamma_d, gamma_a, ksi_i, ksi_d,mu, Q, psi, psi_hat)
    
    %Equations based on Pontryagin's minimum principle
    dy = -[beta*x(2,1)*(l(2,1) - l(1,1))*(1-u) + zeta*(l(7,1) - l(1,1));
    beta*x(1,1)*(l(2,1) - l(1,1))*(1-u) + v*(l(3,1)-l(2,1)) + gamma_i*(l(5,1) - l(2,1)) + ksi_i*(l(4,1) - l(2,1));
    gamma_d*(l(5,1) - l(3,1))  + ksi_d*(l(4,1) - l(3,1));
    Q(4,4)*x(4,1) + gamma_a*(l(5,1) - l(4,1)) + mu*(l(6,1) - l(4,1));
    psi_hat*(l(1,1) - l(5,1));
    0;
    psi*(l(1,1) - l(7,1))];
    y = l - dy*dt; %backwards in time, since boundary condition for costate variables is at t = T


    
    
    
    
    
    
    
    
    
    
    
    
    
