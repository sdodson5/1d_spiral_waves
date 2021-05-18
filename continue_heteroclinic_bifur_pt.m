%% Continue the heteroclinic bifurcation point (epsilon*) in system parameters 
% Given an solution, use secant continuation methods to continue the 
% heteroclinic bifurcation point as a spatiotemporal equilibrium solution
%
% Stephanie Dodson
% May 2021
% UC Davis, sadodson@ucdavis.edu
%  Input file: contains (1) slow wave, (2) fast wave, 
%   (3) unstable eigenfunction of slow wave, 
%   (4) unstable adjoint eigenfunction of slow wave


close all; clear all;

%% File names and solver code

file_names.out_file = 'ML_heteroclinic_connection_continuation';     % Name of output file
file_names.initial_data = 'ML_SsF_soln.mat';  % Name of input file
file_names.solver = 'ML_1Dspiral_heteroclinic_bifur'; % Function to solve for the S to sF heteroclinic connection as a spatiotemporal equilibrium

%% Set up

global Us; global Uf; 
global Uefcn; global Uadjefcn;

ic = load(file_names.initial_data);
Us = ic.Us;             % Slow pulse
Uf = ic.Uf;             % Fast pulse
Uefcn = ic.Uefcn;       % Unstable eigenfunction of slow pulse
Uadjefcn = ic.Uadjefcn; % Adjoint eigenfunction of slow pulse
par = ic.par;           % System parameters
numPar = ic.numPar;     % Numerical parameters for 1D problems
mesh_params = ic.mesh_params;   % Numerical parameters and differentiation matrices for heteroclinic bifurcation
w0 = ic.SsF;            % input/guess for S to sF solution

% Ensure free parameters are up to date from previous continuations
par.ep = ic.SsF(end-2);
par.cs = ic.SsF(end-1);
par.Lt = ic.SsF(end);

% Continuation parameters
contPar.Free1 = 'ep';       % Free parameter 1: epsilon - gives epsilon *
contPar.Free2 = 'cs';       % Free parameter 2: cs - speed of slow pulse
contPar.Free3 = 'Lt';       % Free parameter 3: Time domain length
contPar.Name  = 'u4a';      % Continuation parameter
contPar.final = 5;          % Final continuation parameter
contPar.ds = 10;            % Continuation step size
contPar.initial_ds = 0.001; % Initial step of continuation: Sign gives direction of initial continuation 
contPar.numContSteps = 10;  % Max number of continuation step
contPar.plot_iter = 5;      % Plot continuation diagram every plot_iter steps

% Function set-up
addpath Util/
options = optimset('Jacobian','on','Display','iter','MaxIter',50,'TolFun',1e-6); % options for fsolve
fh_solver = str2func(file_names.solver); % Function to solve for heteroclinic bifurcation point

%% Compute initial points for secant continuation: Solve the initial condition

% Vectors to store free and continuation parameters
free_vec = zeros(3,contPar.numContSteps+2); % Stores values of free parameters
cont_vec = zeros(contPar.numContSteps+2,1); % Stores value of continuation parameter

% Solve first solution (input solution)
%  Phase conditions:
U0 = w0(1:numPar.nt*numPar.nx);
phase_cond.u_old_x = U0(mesh_params.phase_idx_space);
phase_cond.u_old_t = U0(mesh_params.phase_idx_time);

uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params,phase_cond),w0,options);
wU0 = uout(1:2*mesh_params.nx*mesh_params.nt);
par.(contPar.Free1) = uout(end-2); free_vec(1,1) = par.(contPar.Free1);
par.(contPar.Free2) = uout(end-1); free_vec(2,1) = par.(contPar.Free2);
par.(contPar.Free3) = uout(end);   free_vec(3,1) = par.(contPar.Free3);
cont_vec(1) = par.(contPar.Name);

sol0 = [uout; par.(contPar.Name)];  % First point

% Solve second solution
% Update parameter
par.(contPar.Name) = par.(contPar.Name) + contPar.initial_ds;   

% Update phase condition
U0 = wU0(1:numPar.nt*numPar.nx);
phase_cond.u_old_x = U0(mesh_params.phase_idx_space);
phase_cond.u_old_t = U0(mesh_params.phase_idx_time);

uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params,phase_cond),uout,options);
par.(contPar.Free1) = uout(end-2); free_vec(1,2) = par.(contPar.Free1);
par.(contPar.Free2) = uout(end-1); free_vec(2,2) = par.(contPar.Free2);
par.(contPar.Free3) = uout(end);   free_vec(3,2) = par.(contPar.Free3);
cont_vec(2) = par.(contPar.Name);

sol1 = [uout; par.(contPar.Name)]; % Second point

%% Continuation
for k = 1:contPar.numContSteps
    
    disp(k)
    disp(['Continiation Parameter: ' num2str(par.(contPar.Name))])
    disp(['Free 1: ' num2str(par.(contPar.Free1))])
    tmp_start = tic;
    
    % Predictor
    pred = sol1 + (sol1 - sol0)/norm(sol1 - sol0, 2) * contPar.ds;
    
    % Phase condition
    U0 = sol1(1:numPar.nt*numPar.nx);
    phase_cond.u_old_x = U0(mesh_params.phase_idx_space);
    phase_cond.u_old_t = U0(mesh_params.phase_idx_time);
     
    % Corrector
    uout = fsolve( @(sol) FixedPointSecantPredictorCorrector(sol,sol1,sol0,par,numPar,mesh_params,contPar,fh_solver,phase_cond),pred,options);
    
    par.(contPar.Free1) = uout(end-3); free_vec(1,k+2) = par.(contPar.Free1);
    par.(contPar.Free2) = uout(end-2); free_vec(2,k+2) = par.(contPar.Free2);
    par.(contPar.Free3) = uout(end-1); free_vec(3,k+2) = par.(contPar.Free3);
    par.(contPar.Name)  = uout(end);   cont_vec(k+2)  = par.(contPar.Name);
    
    % Set up for next iteration
    sol0 = sol1; 
    sol1 = uout; 
    
    tmp_end = toc(tmp_start);
    disp(['fsolve Time: ' num2str(tmp_end)]);  % fsolve time
    
    if mod(k,contPar.plot_iter)== 0 % Optional plotting and saving of results every plot_iter steps - can comment out
        
        figure(1); pcolor(mesh_params.xx',mesh_params.tt',reshape(uout(1:mesh_params.nx*mesh_params.nt),numPar.nx,numPar.nt));
        shading interp; title([contPar.Name num2str(par.(contPar.Name)) ', \epsilon_* = ' num2str(par.ep)]);
        drawnow;
        
        SsF = uout(1:end-1); % Remove continuation parameter
        %save([file_names.out_file '_' num2str(k) '.mat'],'SsF','par','numPar','mesh_params','Us','Uf','Uadjefcn','Uefcn','free_vec','cont_vec')
        
        
    end
    
     % Stopping criteria to end the loop  
    if (par.(contPar.Name)*sign(contPar.initial_ds)) >= (contPar.final*sign(contPar.initial_ds))
        disp(['Continuation parameter reached before specified number of iterations. Final parameter value: ' num2str(par.(contPar.Name))])
        
        % Solve the source defect again with final parameter value
        par.(contPar.Name) = contPar.final;
        uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params,phase_cond),uout(1:end-1),options);
        par.(contPar.Free1) = uout(end-3); free_vec(1,k+2) = par.(contPar.Free1);
        par.(contPar.Free2) = uout(end-2); free_vec(2,k+2) = par.(contPar.Free2);
        par.(contPar.Free3) = uout(end-1); free_vec(3,k+2) = par.(contPar.Free3);
        
        SsF = uout; 
        break;
    end
    
    
end

if k < contPar.numContSteps
    
    free_vec = free_vec(:,1:k+1); % remove extra zeros at end of vectors
    cont_vec = cont_vec(1:k+1);
    SsF = uout(1:end-1); % remove continuation parameter
    
end

U = reshape(uout(1:mesh_params.nx*mesh_params.nt),numPar.nx,numPar.nt);
V = reshape(uout(1+mesh_params.nx*mesh_params.nt:2*mesh_params.nx*mesh_params.nt),numPar.nx,numPar.nt);

figure(1);  pcolor(mesh_params.xx',mesh_params.tt',U); shading interp; 
figure(2);  pcolor(mesh_params.xx',mesh_params.tt',V); shading interp;

figure(3); subplot(3,1,1); 
plot(cont_vec,free_vec(1,:),'o-','linewidth',3);
xlabel(contPar.Name); ylabel('\epsilon_{*}');
set(gca,'fontsize',20,'linewidth',2); box on;

subplot(3,1,2)
plot(cont_vec,free_vec(2,:),'o-');
xlabel(contPar.Name); ylabel('c_s');
set(gca,'fontsize',20,'linewidth',2); box on;

subplot(3,1,3)
plot(cont_vec,free_vec(3,:),'o-');
xlabel(contPar.Name); ylabel('T');
set(gca,'fontsize',20,'linewidth',2); box on;


save([file_names.out_file '_final.mat'],'SsF','par','numPar','mesh_params','Us','Uf','Uadjefcn','Uefcn','free_vec','cont_vec')

