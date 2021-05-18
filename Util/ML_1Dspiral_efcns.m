function [F,J] = ML_1Dspiral_efcns(w,u_infty,par,mesh_params)
% Solve for eigenfunction 
% Free variables: eigenfunction w and temporal eigenvalue lambda

% Rename parameters for convience
L1y = mesh_params.DT;
L2x = mesh_params.D2Z;
nt = mesh_params.nt;
nz = mesh_params.nz;

Gca = par.Gca;  Eca = par.Eca; 
Gk = par.Gk;    Ek = par.Ek;
Gl = par.Gl;    
ep = par.ep;    
delta1 = par.delta1;
delta2 = par.delta2;

U_infty = u_infty(1:nt*nz); % Source defect/1D spiral solution (on short grid)
V_infty = u_infty(nt*nz+1:2*nz*nt);

U = w(1:nt*nz);  % Eigenfunctions
V = w(nt*nz+1:2*nz*nt);
lambda = w(end);

dt = -par.omega;

% Reaction terms (linearized)
fU_U = @(u,v) -Gca.*(0.5.*(sech((u - par.u1)./par.u2).^2)./par.u2).*(u - Eca) - Gca.*m_infty(u,par) - Gk.*v - Gl;
fU_V = @(u,v) -Gk.*(u - Ek);

% Write n-variable in terms of potassium channels opening and closing rates
[alpha, beta, alpha_prime, beta_prime] = ML_alpha_beta(U_infty,par);
fV_U = ep.*( alpha_prime.*(1 - V_infty) - beta_prime.*V_infty);
fV_V = ep.*(-alpha - beta);

line1 =  delta1.* (L2x * U) + dt.*L1y*U + fU_U(U_infty,V_infty).*U + fU_V(U_infty,V_infty).*V - lambda.*U;
line2 =  delta2.* (L2x * V) + dt.*L1y*V + fV_U.*U + fV_V.*V - lambda.*V;
line3 = [U;V]'*[U;V]-1;  % Phase condition: norm of U-eigenfunction = 1
 

F = [line1; line2; line3]; 

% Jacobain
if nargout > 1
        
     dwU = [delta1.*L2x  + dt.*L1y + spdiags(fU_U(U_infty,V_infty) - lambda,0,nz*nt,nz*nt);
         spdiags(fV_U,0,nz*nt,nz*nt);
         2.*U'];

     dwV = [spdiags(fU_V(U_infty,V_infty),0,nz*nt,nz*nt);
          delta2.*L2x + dt.*L1y + spdiags( fV_V- lambda,0,nz*nt,nz*nt);
         2.*V'];   
     
     dlambda = [-U;-V; 0];
   
     J = [dwU, dwV, dlambda]; 
     
     J = sparse(J);
     
end







