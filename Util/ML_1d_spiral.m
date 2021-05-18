function [F,J] = ML_1d_spiral(w,par,numPar,mesh_params)
% Solve for 1D spiral (time-periodic source defects) in the ML system as a
% 2D spatiotemporal equilibrium solution
% Written in terms of n opening and closing rates alpha and beta

global U0; global global_kappa;

% Rename parameters for convience
L1y = mesh_params.DT;
L2x = mesh_params.D2Z;
nt = mesh_params.nt;
nz = mesh_params.nz;
L2X_long = mesh_params.L2X_long;

Gca = par.Gca;  Eca = par.Eca; 
Gk = par.Gk;    Ek = par.Ek;
Gl = par.Gl;    El = par.El;
ep = par.ep;    
delta1 = par.delta1;
delta2 = par.delta2;

% Input solution
wU = w(1:nt*nz);
wV = w(nt*nz+1:2*nz*nt);
omega = w(end);

% Find the far-field solution
par.omega = omega;
[u_ff,v_ff,Du_ff,kappa] = get_ML_ff(U0,par,numPar,mesh_params,mesh_params.Dt,mesh_params.Dt2);  
u_ff = u_ff(:); v_ff = v_ff(:); Du_ff = Du_ff(:);

par.kappa = kappa;
global_kappa = kappa;
dt = -omega;

% Reaction terms
fU = @(u,v) -Gca.* m_infty(u,par).*(u - Eca) - Gk.*v.*(u - Ek) - Gl.*(u - El) + par.I;
fV = @(u,v) ep.*(ML_alpha1(u,par).*(1 - v) - ML_beta1(u,par).*v); % Written in terms of alpha and beta

D2u = L2X_long * (mesh_params.chi_ffLong.*u_ff); D2u = D2u(mesh_params.bc_idx); % Have a long & short grid to remove boundary effects - take spatial derivatives on long grid, solve for core solution on short grid
D2v = L2X_long * (mesh_params.chi_ffLong.*v_ff); D2v = D2v(mesh_params.bc_idx);

u_ff_short = mesh_params.chi_ff.*u_ff(mesh_params.bc_idx);
v_ff_short = mesh_params.chi_ff.*v_ff(mesh_params.bc_idx);
Du_ff      = mesh_params.chi_ff.*Du_ff(mesh_params.bc_idx);

line1 = dt.*L1y*(u_ff_short + wU) + delta1.*D2u + delta1.* (L2x * wU) + fU(u_ff_short + wU, v_ff_short + wV);
line2 = dt.*L1y*(v_ff_short + wV) + delta2.*D2v + delta2.* (L2x * wV) + fV(u_ff_short + wU, v_ff_short + wV);

% These set the core solution to 0 at the end of the domain
 line1(mesh_params.iend) = wU(mesh_params.iend); % dirchlet bcs at the LHS/RHS of the domain for w
 line2(mesh_params.iend) = wV(mesh_params.iend);
 
% Phase condition: one sided - far-right of the domain
L_cut = mesh_params.Lz - 2*pi/(par.kappa); 
iBm = find((mesh_params.xx(:)) >= L_cut ); % find indices near the end of the domain

z = linspace(-mesh_params.Lz,mesh_params.Lz,mesh_params.nz);  hz= z(2) - z(1);
iim= find((z) >= L_cut); nnzm = length(iim);

wz = [1, 2*ones(1,nnzm-2)+2*mod([1:nnzm-2],2),1]*hz/3;   % Simpson weights for intergration int = w*u
wt = 2*mesh_params.Lt*ones(mesh_params.nt,1)/mesh_params.nt;
wwm= kron(wz,wt');
wwm= wwm(:)';
  
wBm = wU(iBm);                       % find w on the domain x = 0:2*pi/kappa and y = 0..Ly
u_prime = Du_ff(iBm);
line3 = wwm*(u_prime .* wBm);

F = [line1; line2; line3]; 

% Jacobain
if nargout > 1
        
    fU_U = @(u,v) -Gca.*(0.5.*(sech((u - par.u1)./par.u2).^2)./par.u2).*(u - Eca) - Gca.*m_infty(u,par) - Gk.*v - Gl;
    fU_V = @(u,v) -Gk.*(u - Ek);
    alpha_prime = @(u) ( sinh((u - par.u3a)./(2.*par.u4a)).*( 0.25 - 0.25.*tanh((par.u3a - u)./(par.u4a))) + 0.5.*cosh((par.u3a - u)./(2.*par.u4a)).*(sech((par.u3a - u)./(par.u4a))).^2 )./par.u4a;
    beta_prime = @(u) ( sinh((u - par.u3b)./(2.*par.u4b)).*( 0.25.*tanh((par.u3b - u)./par.u4b) + 0.25) - 0.5.*cosh((par.u3b - u)./(2.*par.u4b)).*(sech((par.u3b - u)./par.u4b)).^2 )./par.u4b;
    fV_U = @(u,v) ep.*( alpha_prime(u).*(1 - v) - beta_prime(u).*v);
    fV_V = @(u,v) ep.*(-ML_alpha1(u,par) - ML_beta1(u,par));
    
    phase_jacob = zeros(1,nz*nt);
    phase_jacob(iBm) = wwm*spdiags(u_prime,0,nnzm*mesh_params.nt,nnzm*mesh_params.nt);
    
     dwU = [dt.*L1y + delta1.*L2x + spdiags(fU_U(u_ff_short + wU,v_ff_short+wV),0,nz*nt,nz*nt);
         spdiags(fV_U(u_ff_short+wU, v_ff_short + wV),0,nz*nt,nz*nt);
         phase_jacob];

     dwV = [spdiags(fU_V(u_ff_short + wU,v_ff_short+wV),0,nz*nt,nz*nt);
         dt.*L1y + delta2.*L2x + spdiags(fV_V(u_ff_short + wU,v_ff_short+wV),0,nz*nt,nz*nt);
         sparse(1,nz*nt)];   
     
     epsiF = 1e-8;
     dF = ML_1d_spiral([w(1:end-1); w(end)+epsiF],par,numPar,mesh_params);
     dT = (dF - F)./epsiF;  
   
     J = [dwU, dwV, dT];
   
     % Boundary conditions in Jacobian
      J(mesh_params.iend,:) = 0; J(nz*nt+mesh_params.iend, :) = 0; 
      J(mesh_params.iend,mesh_params.iend)=speye(length(mesh_params.iend));
      J(nz*nt+mesh_params.iend, nz*nt+mesh_params.iend)=speye(length(mesh_params.iend));    
     
     J = sparse(J);
     
end







