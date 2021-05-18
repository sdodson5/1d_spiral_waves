function [F, J]  = FixedPointSecantPredictorCorrector(sol,sol1,sol0,par,numPar,mesh_params,contPar,fh_solver,phase_cond)
% Secant predoctor corrector to find fixed point
% Pass in general function handle

u = sol(1:end-1); % Predictor 
p = sol(end);  % Continuation parameter

u1 = sol1(1:end-1); % second point
p1 = sol1(end); 

u0 = sol0(1:end-1); % first point
p0 = sol0(end);

d = norm( sol1 - sol0 , 2 );

alpha = (u1 - u0) / d;
beta  = (p1 - p0) / d;

par.(contPar.Name) = p; % update the continuation parameter to value of predictor

F = zeros(size(sol));
if nargin > 8
    [Fu,Ju] = fh_solver(u,par,numPar,mesh_params,phase_cond);
else
    [Fu,Ju] = fh_solver(u,par,numPar,mesh_params); % evaluate system at predictor u and p
end
F(1:end-1) = Fu;
F(end)    = alpha' * (u - u1) + beta * (p - p1) - contPar.ds; % corrector condition

if nargout > 1 % Redefine the Jacobian - need to find the Jacobian wrt the continuation parameter
    h = 1e-8;
    parpe = par;
    parpe.(contPar.Name) = p+h; % Add a small amount to the continuation parameter
    
    if nargin > 8
        Fu2 = fh_solver(u,parpe,numPar,mesh_params,phase_cond);
    else
        Fu2 = fh_solver(u,parpe,numPar,mesh_params); % evaluate system at predictor u and p
    end

    
    Fp = (Fu2 - Fu)/h; % "derivative" of the system wrt the continuation parameter    
    J = [Ju, Fp; % d_u d_v d_omega d_par 
         alpha', beta];
     
    J = sparse(J);
end

