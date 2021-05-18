function [f,J]  = Morris_Lecar_1D_kappa_alpha_beta(u,L1,L2,par,numPar,phase_cond)
% Solves for a wave train solution in the 1D Morris-Lecar equation
% Free parameter is wave number kappa
% Written in terms of potassium opening and closing rates (alpha and beta)

% Rename parameters for convience
nx = numPar.nx;
Gca = par.Gca;  Eca = par.Eca; 
Gk = par.Gk;    Ek = par.Ek;
Gl = par.Gl;    El = par.El;
ep = par.ep;    
delta1 = par.delta1;
delta2 = par.delta2;
omega = par.omega;

% Input initial condition
U = u(1:nx);
V = u(nx+1:2*nx);
kappa = u(end); % Free parameter

% Derivative of phase condition
u_phase_th = L1*phase_cond.u_old(1:numPar.nx); 

% Reaction terms
[m_inf,m_inf_prime]  = m_infty(U,par);
[alpha, beta, alpha_prime, beta_prime] = ML_alpha_beta(U,par);

fU = -Gca.*m_inf.*(U - Eca) - Gk.*V.*(U - Ek) - Gl.*(U - El) + par.I;
fV = ep.*(alpha.*(1 - V) - beta.*V);  % written in terms of alpha and beta

line1 = delta1.*kappa^2.*(L2*U) + omega.*(L1*U) + fU;
line2 = delta2.*kappa^2.*(L2*V) + omega.*(L1*V) + fV;
line3 = (u_phase_th)'*(U - phase_cond.u_old);

f = [line1;
    line2;
    line3];
   
% Jacobian              
if (nargout > 1)
    
    % Partial derivatives
    fU_U = -Gca.*m_inf_prime.*(U - Eca) - Gca.*m_inf - Gk.*V - Gl;
    fU_V = -Gk.*(U - Ek);
    fV_U = ep.*( alpha_prime.*(1 - V) - beta_prime.*V);
    fV_V = ep.*(-alpha - beta);
    
    dU = [delta1.*kappa^2.*L2 + omega.*L1 + spdiags(fU_U,0,nx,nx);
        spdiags(fV_U,0,nx,nx);
        u_phase_th'];
    
    dV = [spdiags(fU_V,0,nx,nx);
        delta2.*kappa^2.*L2 + omega.*L1 + spdiags(fV_V,0,nx,nx);
        sparse(1,nx)];
   
    
    dkappa = [(2*kappa*delta1)*L2*U;
        (2*kappa*delta2)*L2*V;
        0];
    
    J = [dU, dV, dkappa];
    J = sparse(J);
             
end


    
    
    
    
    
    
    