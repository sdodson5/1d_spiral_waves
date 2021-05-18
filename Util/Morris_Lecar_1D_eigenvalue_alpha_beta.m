function [f,J]  = Morris_Lecar_1D_eigenvalue_alpha_beta(u,u_infty,L1,L2,par,numPar)
% Compute eigenvalue-eigenvector pair
% Written in terms of opening and closeing rates of potassium channels
% (alpha and beta)

% Rename parameters for convience
nx = numPar.nx;

Gca = par.Gca;  Eca = par.Eca; 
Gk = par.Gk;    Ek = par.Ek;
Gl = par.Gl;    El = par.El;
ep = par.ep;    
delta1 = par.delta1;
delta2 = par.delta2;
kappa = par.kappa;

% Solution to linearize about
u2 = u_infty(1:nx);
v2 = u_infty(nx+1:2*nx);
omega = u_infty(end);

U = u(1:nx);
V = u(nx+1:2*nx);
lambda = u(end); 

% Linearized reaction terms
fU_U = @(u,v) -Gca.*(0.5.*(sech((u - par.u1)./par.u2).^2)./par.u2).*(u - Eca) - Gca.*m_infty(u,par) - Gk.*v - Gl;
fU_V = @(u,v) -Gk.*(u - Ek);

[alpha, beta, alpha_prime, beta_prime] = ML_alpha_beta(u2,par);
fV_U =  ep.*( alpha_prime.*(1 - v2) - beta_prime.*v2 );
fV_V =  ep.*( - alpha - beta);

    
line1 = delta1.*kappa^2.*(L2*U) + omega.*(L1*U) + fU_U(u2,v2).*U + fU_V(u2,v2).*V - lambda.*U;
line2 = delta2.*kappa^2.*(L2*V) + omega.*(L1*V) + fV_U.*U + fV_V.*V - lambda.*V;
line3 = [U;V]'*[U;V] - 1;

f = [line1;
    line2;
    line3];
   
% Jacobian              
if (nargout > 1)
    
    dU = [delta1.*kappa^2.*L2 + omega.*L1 + spdiags(fU_U(u2,v2) - lambda,0,nx,nx);
        spdiags(fV_U,0,nx,nx);
        2*U'];
    
    dV = [spdiags(fU_V(u2,v2),0,nx,nx);
        delta2.*kappa^2.*L2 + omega.*L1 + spdiags(fV_V - lambda,0,nx,nx);
        2*V'];
   
    
    dlambda = [-U;
        -V;
        0];
    
        
    J = [dU, dV, dlambda];
    J = sparse(J);
   
          
end


    
    
    
    
    
    
    