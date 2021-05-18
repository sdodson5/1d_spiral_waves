function [F,J] = ML_1Dspiral_heteroclinic_bifur(w,par,numPar,mesh_params,phase_cond_full)
% Solve for the heteroclinic bifurcation point that creates 1D spiral
% Slow pulse splits into oppositely propagating slow-fast pulse pair
% Written in terms of opening and closing rates of potassium channel (alpha
% and beta)


global Us; global Uf; 
global Uefcn; global Uadjefcn;

options = optimset('Display','off','Jacobian','on', 'DerivativeCheck','off',...
                    'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',200);
                
% Rename parameters for convience
L1x = mesh_params.L1x;  % 1D differentiation matrices
L2x = mesh_params.L2x;
L1t = mesh_params.L1t;
D1X = mesh_params.D1X;  % 2D differentiation matrices
D2X = mesh_params.D2X;
DT = mesh_params.DT;
nx = mesh_params.nx;
nt = mesh_params.nt;

Gca = par.Gca;  Eca = par.Eca; 
Gk = par.Gk;    Ek = par.Ek;
Gl = par.Gl;    El = par.El;    
delta1 = par.delta1;
delta2 = par.delta2;

% slow to slow-fast wave pattern (pattern solved for)
U = w(1:nt*nx);
V = w(nt*nx+1:2*nx*nt);

ep = w(end-2); par.ep = ep; % Free parameter 1: epsilon
cs = w(end-1); par.cs = cs; % Free parameter 2: cs - speed of slow pulse
T = w(end);                 % Free parameter 3: time domain is [-T,T]

% Solve for the slow pulse with updated cs & epsilon
phase_cond.u_old = Us(1:mesh_params.nx);  % Phase condition: previous slow pulse
uout = fsolve(@(y) Morris_Lecar_1D_alpha_beta(y,L1x,L2x,par,numPar,phase_cond),[Us;par.cs],options); % Solve for slow pulse
Us = uout(1:end-1); % New slow pulse

% Solve for the fast pulse with updated epsilon
phase_cond.u_old = Uf(1:mesh_params.nx);
uout = fsolve(@(y) Morris_Lecar_1D_alpha_beta(y,L1x,L2x,par,numPar,phase_cond),[Uf;par.cf],options); % Solve fast wave
Uf = uout(1:end-1);  % New fast pulse
par.cf = uout(end);  % Speed of fast wave

% Compute unstable eigenvector of slow wave: for use in time = -T boundary
% condition
uout = fsolve(@(y) Morris_Lecar_1D_eigenvalue_alpha_beta(y,Us,L1x,L2x,par,numPar),[Uefcn;par.lambda],options); 
Uefcn = uout(1:end-1);  % Updated unstable eigenvector of slow wave
par.lambda = uout(end); % Value not important, but is used in adjoint equation

% Reaction terms
[alpha, beta, alpha_prime, beta_prime] = ML_alpha_beta(U,par);
fU = @(u,v) -Gca.* m_infty(u,par).*(u - Eca) - Gk.*v.*(u - Ek) - Gl.*(u - El) + par.I;
fV = ep.*(alpha.*(1 - V) - beta.*V); 


% Dirichlet boundary conditions at t = -T: U,V are slow wave plus a little bit in direction of unstable eigenvector
A = (Us(1:nx)      + par.efcn_w.*abs(Uefcn(1:nx)) ); % V-component
B = (Us(nx+1:2*nx) + par.efcn_w.*abs(Uefcn(nx+1:2*nx)) ); % n-component

tmpU = -DT*[A; U]; tmpU = tmpU(nx+1:end);
tmpV = -DT*[B;V]; tmpV = tmpV(nx+1:end);

% Equations for slow wave to slow-fast wave splits: temporally scaled to [-1,1]
line1 = T.*(delta1.*D2X*U + cs.*D1X*U + fU(U,V)) + tmpU;
line2 = T.*(delta2.*D2X*V + cs.*D1X*V + fV) + tmpV;
 
u_phase_t = L1t*phase_cond_full.u_old_t./T; 
line3 = (u_phase_t)'*(U(mesh_params.phase_idx_time) - phase_cond_full.u_old_t);  % Temporal phase condition for ep

u_phase_x = L1x * phase_cond_full.u_old_x; 
line4 = (u_phase_x)'*(U(mesh_params.phase_idx_space) - phase_cond_full.u_old_x); % Spatial phase condition for cs

% Time = T: Need to be close to linear combination of slow and fast waves
M = construct_M_pulse(par.qs,par.D,par,numPar,mesh_params);
line5 = ([U(mesh_params.bc_idx1); V(mesh_params.bc_idx1)] - M)'*abs(Uadjefcn);  % Phase condition for T: Heteroclinic boundary condition at t = T

F = [line1; line2; line3; line4; line5]; 

% Jacobain
if nargout > 1
        
    fU_U = @(u,v) -Gca.*(0.5.*(sech((u - par.u1)./par.u2).^2)./par.u2).*(u - Eca) - Gca.*m_infty(u,par) - Gk.*v - Gl;
    fU_V = @(u,v) -Gk.*(u - Ek);
    fV_U =  ep.*( alpha_prime.*(1 - V) - beta_prime.*V );
    fV_V =  ep.*( - alpha - beta);

    phase_jacob_x = zeros(1,nx*nt);
    phase_jacob_x(mesh_params.phase_idx_space) = (u_phase_x)';

    phase_jacob_t = zeros(1,nx*nt);
    phase_jacob_t(mesh_params.phase_idx_time) = (u_phase_t)';
    
    phase_line5_jacob_u = zeros(1,nx*nt);
    phase_line5_jacob_u(mesh_params.bc_idx1) = abs(Uadjefcn(1:nx))';
    
    phase_line5_jacob_v = zeros(1,nx*nt);
    phase_line5_jacob_v(mesh_params.bc_idx1) = abs(Uadjefcn(nx+1:2*nx))';
    
    dwU = [-DT(numPar.nx+1:end,numPar.nx+1:end) + T.*(delta1.*D2X + cs.*D1X + spdiags(fU_U(U,V),0,nx*nt,nx*nt));
         T.*spdiags(fV_U,0,nx*nt,nx*nt);
         phase_jacob_t;
         phase_jacob_x;
         phase_line5_jacob_u];

     dwV = [T.*spdiags(fU_V(U,V),0,nx*nt,nx*nt);
       -DT(numPar.nx+1:end,numPar.nx+1:end) + T.*(delta2.*D2X + cs.*D1X + spdiags(fV_V,0,nx*nt,nx*nt));
        sparse(1,nx*nt);
        sparse(1,nx*nt);
        phase_line5_jacob_v];   
          
     epsiF = 1e-8;
     
     dF = ML_1Dspiral_heteroclinic_bifur([w(1:end-3); w(end-2)+epsiF; w(end-1:end)],par,numPar,mesh_params,phase_cond_full);
     dep = (dF - F)./epsiF;
    
     dF = ML_1Dspiral_heteroclinic_bifur([w(1:end-2); w(end-1)+ epsiF; w(end)],par,numPar,mesh_params,phase_cond_full);
     dcs = (dF - F)./epsiF;      
                  
     dF = ML_1Dspiral_heteroclinic_bifur([w(1:end-1); w(end)+ epsiF],par,numPar,mesh_params,phase_cond_full);
     dT = (dF - F)./epsiF;    
     
     J = [dwU, dwV, dep, dcs, dT];
     J = sparse(J);
     
end







