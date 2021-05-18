function [u_ff0, v_ff0,Du_ff,kappa] = get_ML_ff(u,par,numPar,mesh_params,L1,L2)
% Solves the 1D wave train problem with
% Input: u - [U;V] periodic wave train
%        L1, L2: 1D periodic differentiation matrices (for wave train
%        problem)
% Output: u_ff0, v_ff0: wave train expanded to the (x,t)-plane long grid
%         U: solved wave train
%         Du_ff: derivative of wave train
% Written in terms of potassium channel opening and closing rates (alpha
% and beta)
%%

phase_cond.u_old = u(1:numPar.nx);
u1 = [u; par.kappa];


options = optimset('Display','off','Jacobian','on', 'DerivativeCheck','off',...
                    'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',200);
                
               
uout = fsolve(@(y) Morris_Lecar_1D_kappa_alpha_beta(y,L1,L2,par,numPar,phase_cond),u1,options);
U = uout(1:numPar.nx);
V = uout(numPar.nx+1:2*numPar.nx);
kappa = uout(end);
par.kappa = uout(end);

Du = L1*U;
U2 = flip(U); U2 = [U2(end/2+1:end); U2(1:end/2)];
V2 = flip(V); V2 = [V2(end/2+1:end); V2(1:end/2)];
Du2 = L1*U2;

 % Interpolating on the eta mesh
eta  = par.kappa.*mesh_params.xx2_pos(:) - mesh_params.tt2_pos(:);  % find skewed rolls and interpolate
eta2 = par.kappa.*mesh_params.xx2_neg(:) + mesh_params.tt2_neg(:);  % find skewed rolls and interpolate
  
u_ff01 = zeros(size(eta));
v_ff01 = zeros(size(eta));
Du_ff1= zeros(size(eta));
u_ff02 = zeros(size(eta2));
v_ff02 = zeros(size(eta2));
Du_ff2= zeros(size(eta2));

nt = mesh_params.nt;
for i = 1:nt % Over the times
    zz = eta-mesh_params.tt2_pos(i,1);
    periodic_sinc = diric(zz,nt).*cos(zz/2);
    u_ff01 = u_ff01 + U(i)*periodic_sinc;
    v_ff01 = v_ff01 + V(i)*periodic_sinc;
    Du_ff1 = Du_ff1 + Du(i)*periodic_sinc;
    
    zz2 = eta2-mesh_params.tt2_neg(i,1);
    periodic_sinc2 = diric(zz2,nt).*cos(zz2/2);
    u_ff02 = u_ff02 + U2(i)*periodic_sinc2;
    v_ff02 = v_ff02 + V2(i)*periodic_sinc2;
    Du_ff2 = Du_ff2 + Du2(i)*periodic_sinc2;
   
end

u_ff0 = [u_ff02; u_ff01];
v_ff0 = [v_ff02; v_ff01];
Du_ff = [Du_ff2; Du_ff1];




