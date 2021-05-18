%% Continue the 1D spiral & it's unstable eigenvalue in system parameters 
% Given an solution, use secant continuation methods to continue the 1D
% spiral as a spatiotemporal equilibrium solutions using the core and
% far-field decomposition
% Then, with fixed parameters & 1D spiral solution, unstable eigenvalue is
% solved for
%
% Stephanie Dodson
% May 2021
% UC Davis, sadodson@ucdavis.edu


close all; clear all;

global U0; 
global global_kappa;
%% System to solve

file_names.out_file = 'ML_eigenvalue_continuation';            % Name of output file
file_names.initial_condition = 'ML_1dspiral_efcn_soln.mat';    % Initial condition

file_names.contruct_initial_solution = 'get_ML_ff';               % Construct the far-field solution and initial condition
file_names.pattern_solver = 'ML_1d_spiral';                       % Function to solve for the pattern
file_names.eval_solver = 'ML_1Dspiral_efcns';                % Function to solve the eigenvalue problem 
%% Set ups
               
ic = load(file_names.initial_condition);
mesh_params = ic.mesh_params;                    % Numerical parameters for the 1D spiral pattern
wU0 = ic.wU(1:2*mesh_params.nt*mesh_params.nz);  % 1D spiral pattern input
wW0 = ic.wW(1:2*mesh_params.nt*mesh_params.nz);  % eigenvector input
U0 = ic.U0;         % Far-field wave train
par = ic.par;       % System parameters
numPar = ic.numPar; % Some numerical parameters
global_kappa = par.kappa;

% Continuation parameters
contPar.Free = 'omega';         % Free parameter for 1D spiral solution
contPar.evalFree = 'lambda';    % Free parameter for eigenvalue computation
contPar.Name = 'Gca';           % Continuation parameter
contPar.final = 4.35;             % Final continuation parameter value
contPar.ds = 30;                % Secant continuation step size
contPar.initial_ds = -0.0001;   % Initial step of continuation: Sign gives direction of initial continuation 
contPar.numContSteps = 200;     % Maximum number of continuation steps
contPar.plot_iter = 10;         % Plot continuation diagram every plot_iter steps

% Additional parameters: Set up differentiation matrices and things
addpath Util/
options = optimset('Jacobian','on','Display','iter','MaxIter',100,'DerivativeCheck','off','TolFun',1e-6);

fh_rolls = str2func(file_names.contruct_initial_solution);  % Function to construct far-field pattern
fh_solver = str2func(file_names.pattern_solver);            % Function to solve for core solution
fh_eval = str2func(file_names.eval_solver);                 % Function to solve for unstable eigenvalue & eigenvector

%% Solve for initial continuation points

% Compute initial point (pattern)
w0 = [wU0;  par.(contPar.Free)]; % Input solution: core pattern appended with free parameter
uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params),w0,options); % Solve for core pattern
wU0 = uout(1:2*mesh_params.nz*mesh_params.nt); % 1D spiral solution  
par.(contPar.Free) = uout(end);     % Update free parameter
sol0 = [uout; par.(contPar.Name)];  % Solution 1: core solution appended with free parameter and continuation parameter

% With fixed parameters, solve eigenvalue problem: 
u_infty = form_full_source_defect_pattern(wU0,par,numPar,mesh_params,fh_rolls); % compute the full 1D spiral to linearize about
% Solve for the eigenfunction
v0 = [wW0; par.(contPar.evalFree)]; % Eigenfunction guess appended with free parameter
uout_eval = fsolve(@(y)  fh_eval(y,u_infty,par,mesh_params),v0,options); % Solve for the unstable eigenfunction & eigenvalue
par.(contPar.evalFree) = uout_eval(end);  % Update free parameter

 
%% Compute the second point (pattern)
par.(contPar.Name) = par.(contPar.Name) + contPar.initial_ds;   % Add initial step to continuation parameter
uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params),uout,options);
par.(contPar.Free) = uout(end);
sol1 = [uout; par.(contPar.Name)];

% With fixed parameters, solve eigenvalue problem
u_infty = form_full_source_defect_pattern(uout(1:end-1),par,numPar,mesh_params,fh_rolls); % compute the full 1D spiral to linearize about
% Solve for the eigenfunction
uout_eval = fsolve(@(y)  fh_eval(y,u_infty,par,mesh_params),uout_eval,options);
par.(contPar.evalFree) = uout_eval(end);

%% Continuation diagram
free_vec = zeros(contPar.numContSteps,1); % Vectors of free and continuation parameters
cont_vec = zeros(contPar.numContSteps,1);
kappa_vec = zeros(contPar.numContSteps,1);
eval_vec = zeros(contPar.numContSteps,1);

free_vec(1) = par.(contPar.Free);
cont_vec(1) = par.(contPar.Name);
kappa_vec(1) = global_kappa;
eval_vec(1) = par.(contPar.evalFree);

%% Continuation
con_start = tic;

for k = 1:contPar.numContSteps
    disp(k)
    disp(['Continiation Parameter: ' num2str(par.(contPar.Name))])
    disp(['Free: ' num2str(par.(contPar.Free))])
    tmp_start = tic;
    
    %% Continue the 1D spiral with secant method
    % Predictor
    pred = sol1 + (sol1 - sol0)/norm(sol1 - sol0, 2) * contPar.ds;
    
    % Corrector - compute linear operator in the corrector
    uout = fsolve( @(sol) FixedPointSecantPredictorCorrector(sol,sol1,sol0,par,numPar,mesh_params,contPar,fh_solver),pred,options);
    par.(contPar.Free) = uout(end-1); 
    par.(contPar.Name) = uout(end);
    
    %% With fixed parameters, solve eigenvalue problem
    u_infty = form_full_source_defect_pattern(uout(1:end-2),par,numPar,mesh_params,fh_rolls); % Construct full pattern to linearize about

    % Solve for the eigenfunction
    uout_eval = fsolve(@(y)  fh_eval(y,u_infty,par,mesh_params),uout_eval,options);
    par.(contPar.evalFree) = uout_eval(end);

    free_vec(k+1) = par.(contPar.Free);
    cont_vec(k+1) = par.(contPar.Name);
    kappa_vec(k+1) = global_kappa;
    eval_vec(k+1) = par.(contPar.evalFree);
    
    sol0 = sol1; 
    sol1 = uout; 
    
    tmp_end = toc(tmp_start);
    disp(['fsolve Time: ' num2str(tmp_end)]);  % fsolve time
    %%
    if mod(k,contPar.plot_iter) == 0 
         
        figure(1); plot(cont_vec(1:k+1),free_vec(1:k+1),'-o','LineWidth',2); % continuation diagram
        xlabel(contPar.Name); ylabel(contPar.Free); title('Continuation Parameter'); set(gca,'fontsize',14'); drawnow;
        
        figure(2);  plot(cont_vec(1:k+1),eval_vec(1:k+1),'-o','LineWidth',2); hold on; % show the update in eigenvalue
        xlabel(contPar.Name); ylabel('\lambda'); set(gca,'fontsize',14'); hold off; drawnow;

       plot_oneD_spiral_core(uout,mesh_params); drawnow;
       plot_full_1Dspiral(uout_eval,mesh_params); drawnow;
      
    end
    
     % Stopping criteria:  
    if (par.(contPar.Name)*sign(contPar.initial_ds)) >= (contPar.final*sign(contPar.initial_ds))
        disp(['Continuation parameter reached before specified number of iterations. Final parameter value: ' num2str(par.(contPar.Name))])
        
        % Solve for the 1D spiral again with final parameter value
        par.(contPar.Name) = contPar.final;
        uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params),uout(1:end-1),options);
        par.(contPar.Free) = uout(end-1); 
        par.(contPar.Name) = uout(end);
        
        % Solve for eigenvalue problem
        u_infty = form_full_source_defect_pattern(uout(1:end-1),par,numPar,mesh_params,fh_rolls); % Construct full pattern to linearize about
        uout_eval = fsolve(@(y)  fh_eval(y,u_infty,par,mesh_params),uout_eval,options);
        par.(contPar.evalFree) = uout_eval(end);
        
        free_vec(k+1) = par.(contPar.Free);
        cont_vec(k+1) = par.(contPar.Name);
        kappa_vec(k+1) = global_kappa;
        eval_vec(k+1) = par.(contPar.evalFree);
        
        wU = uout(1:end-1);
        wW = uout_eval(1:end-1);
        break;
    end
        
           
end
par.kappa = global_kappa;

if k < contPar.numContSteps
    cont_vec = cont_vec(1:nnz(cont_vec));
    free_vec = free_vec(1:nnz(cont_vec));
    kappa_vec = kappa_vec(1:nnz(cont_vec));
    eval_vec = eval_vec(1:nnz(cont_vec));    
else
    wU = uout(1:end-2);
    wW = uout_eval(1:end-1);
end

% Make a final plot of the continuation
figure(1); plot(cont_vec,free_vec,'-o','LineWidth',2); % continuation diagram for 1D spiral wave
xlabel(contPar.Name); ylabel(contPar.Free); title('Continuation Parameter'); set(gca,'fontsize',14'); drawnow;

figure(2);  plot(cont_vec,eval_vec,'-o','LineWidth',2); hold on; % Continuation diagram for unstable eigenvalue
xlabel(contPar.Name); ylabel('\lambda'); set(gca,'fontsize',14'); hold off; drawnow;


%save(file_names.out_file,'wU','wW','par','mesh_params','numPar','cont_vec','free_vec','kappa_vec','eval_vec','U0');


