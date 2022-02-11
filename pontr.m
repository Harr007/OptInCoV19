%function of costate variables update based on pontryagin's minimum principle

function [y,dy] = pontr(dt, l, x, u, v, zeta, beta, gamma_i, gamma_d, gamma_a, ksi_i, ksi_d,mu,mu_h, H_th,  Q, psi, psi_hat)

    %Effect of healthcare capacity
    if (x(4,1)) <= H_th
        dmu_1 = mu; %partial mu_1/partial i
    else
        dmu_1 = mu_h; %partial mu_1/partial i
    end
    
    %Equations based on Pontryagin's minimum principle
%     dy = -[beta*x(2,1)*(l(2,1) - l(1,1))*(1-u);
%          beta*x(1,1)*(l(2,1) - l(1,1))*(1-u) + v*(l(3,1)-l(2,1)) + gamma_i*(l(5,1) - l(2,1)) + ksi_i*(l(4,1) - l(2,1));
%          gamma_d*(l(5,1) - l(3,1))  + ksi_d*(l(4,1) - l(3,1));
%          Q(4,4)*x(4,1) + gamma_a*(l(5,1) - l(4,1)) + dmu_1*(l(6,1) - l(4,1));
%          0;
%          0];
%         dy = [beta*x(2,1)*(l(1,1) - l(2,1))*(1-u) + 2*zeta*(l(1,1) - l(7,1));
%          beta*x(1,1)*(l(1,1) - l(2,1))*(1-u) + v*(l(2,1)-l(3,1)) + gamma_i*(l(2,1) - l(5,1)) + ksi_i*(l(2,1) - l(4,1));
%          gamma_d*(l(3,1) - l(5,1))  + ksi_d*(l(3,1) - l(4,1));
%          Q(4,4)*x(4,1) + gamma_a*(l(4,1) - l(5,1)) + dmu_1*(l(4,1) - l(6,1));
%          psi_hat*(l(5,1) - l(1,1));
%          0;
%          psi*(l(7,1) - l(1,1))];
        dy = -[beta*x(2,1)*(l(2,1) - l(1,1))*(1-u) + zeta*(l(7,1) - l(1,1));
         beta*x(1,1)*(l(2,1) - l(1,1))*(1-u) + v*(l(3,1)-l(2,1)) + gamma_i*(l(5,1) - l(2,1)) + ksi_i*(l(4,1) - l(2,1));
         gamma_d*(l(5,1) - l(3,1))  + ksi_d*(l(4,1) - l(3,1));
         Q(4,4)*x(4,1) + gamma_a*(l(5,1) - l(4,1)) + dmu_1*(l(6,1) - l(4,1));
         psi_hat*(l(1,1) - l(5,1));
         0;
         psi*(l(1,1) - l(7,1))];
    y = l - dy*dt; %backwards in time, since boundary condition for costate variables is at t = T


    
    
    
    
    
    
    
    
    
    
    
    
    