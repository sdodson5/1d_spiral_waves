%% Continue the 1D spiral in system parameters 
% Given an solution, use secant continuation methods to continue the 1D
% spiral as a spatiotemporal equilibrium solutions using the core and
% far-field decomposition
%
% Stephanie Dodson
% May 2021
% UC Davis, sadodson@ucdavis.edu

close all; clear all;

%% File names and solver code
file_names.oneD_spiral = 'ML_oneDspiral_soln.mat';  % Initial condition for 1D spiral
file_names.out_file = 'oneD_spiral_continuation';    % Name of output/saved file

file_names.contruct_initial_solution = 'get_ML_ff';           % Function to construct the far-field solution 
file_names.solver = 'ML_1d_spiral';        % Function to solve for the core of the pattern

%% Set up    
global U0;  global global_kappa;

% Load initial condition for the 1D spiral and rename
ic = load(file_names.oneD_spiral);
mesh_params = ic.mesh_params;                       % Numerical parameters and differentiation matrices to solve for 1D spiral
wU0 = ic.wU(1:2*mesh_params.nt*mesh_params.nz);     % 1D spiral core solution
U0 = ic.U0;         % 1D far-field wave train
par = ic.par;       % System parameters
numPar = ic.numPar; % Numerical parameters for 1D cases

global_kappa = par.kappa;

% Continuation parameters
contPar.Free = 'omega';         % Free parameter - set up for omega
contPar.Name = 'u4b';           % Continuation parameter: can use any parameter in structure 'par'
contPar.final = 20;             % Final parameter value
contPar.ds = 50;                % Secant continuation step size
contPar.initial_ds = 0.001;     % Initial step of continuation: Sign gives direction of initial continuation 
contPar.numContSteps = 20;      % Max number of continuation steps
contPar.plot_iter = 10;         % Plot continuation diagram every plot_iter steps

% Function set-up
addpath Util/ 
options = optimset('Jacobian','on','Display','iter','MaxIter',100,'DerivativeCheck','off','TolFun',1e-6); % Options for fsolve
fh_rolls = str2func(file_names.contruct_initial_solution);  % Function to construct far-field wave train pattern
fh_solver = str2func(file_names.solver);          % Function to solve for core pattern

%% Solve input pattern for two initial points for secant method

% Compute initial point
w0 = [wU0;  par.(contPar.Free)]; % Input solution: core pattern appended with free parameter
uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params),w0,options); % Solve for core pattern
par.(contPar.Free) = uout(end);     % Update free parameter
sol0 = [uout; par.(contPar.Name)];  % Solution 1: core solution appended with free parameter and continuation parameter
%plot_boundary_sink(uout(1:end-1),mesh_params); % Plot the solved core pattern


% Compute the second point
par.(contPar.Name) = par.(contPar.Name) + contPar.initial_ds;           % Add initial step to continuation parameter
uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params),uout,options);  % Solve for core pattern
par.(contPar.Free) = uout(end);     % Update free parameter
sol1 = [uout; par.(contPar.Name)];  % solution 2

%% Continuation diagram
free_vec = zeros(contPar.numContSteps,1); % Vector of free parameter
cont_vec = zeros(contPar.numContSteps,1);
kappa_vec = zeros(contPar.numContSteps,1);
free_vec(1) = par.(contPar.Free);
cont_vec(1) = par.(contPar.Name);
kappa_vec(1) = global_kappa;

%% Continuation
con_start = tic;  % Monitor time

for k = 1:contPar.numContSteps
    
    disp(k)
    disp(['Continiation Parameter: ' num2str(par.(contPar.Name))])
    disp(['Free: ' num2str(par.(contPar.Free))])
    tmp_start = tic;
    
    % Predictor
    pred = sol1 + (sol1 - sol0)/norm(sol1 - sol0, 2) * contPar.ds;
    
    % Corrector 
    uout = fsolve( @(sol) FixedPointSecantPredictorCorrector(sol,sol1,sol0,par,numPar,mesh_params,contPar,fh_solver),pred,options);
    par.(contPar.Free) = uout(end-1); 
    par.(contPar.Name) = uout(end);
    
    % Save continuation information
    free_vec(k+1) = par.(contPar.Free);
    cont_vec(k+1) = par.(contPar.Name);
    kappa_vec(k+1) = global_kappa;
    
    % Prepare for next continuation step
    sol0 = sol1; 
    sol1 = uout; 
    
    tmp_end = toc(tmp_start);
    disp(['fsolve Time: ' num2str(tmp_end)]);  % fsolve time
    
    if mod(k,contPar.plot_iter) == 0 % Plot and save updates - can comment out to run faster 
         
        figure(1); plot(cont_vec(1:k+1),free_vec(1:k+1),'-o','LineWidth',2); % continuation diagram
        xlabel(contPar.Name,'FontSize',14); ylabel(contPar.Free,'FontSize',14); title('Continuation Parameter','FontSize',14); drawnow;
        
        % Make plot
        plot_oneD_spiral_core(uout,mesh_params); drawnow;
        plot_full_1Dspiral(uout,U0,par,numPar,mesh_params,fh_rolls); drawnow;
 
        % Save temporary file
        %save([file_names.out_file '_' num2str(k) '.mat'],'uout','par','mesh_params','numPar','cont_vec','free_vec','kappa_vec','U0')

    end
    
    % Stopping criteria to end the loop  
    if (par.(contPar.Name)*sign(contPar.initial_ds)) >= (contPar.final*sign(contPar.initial_ds))
        disp(['Continuation parameter reached before specified number of iterations. Final parameter value: ' num2str(par.(contPar.Name))])
        
        % Solve for the 1D spiral again with final parameter value
        par.(contPar.Name) = contPar.final;
        uout = fsolve(@(y)  fh_solver(y,par,numPar,mesh_params),uout(1:end-1),options);
        free_vec(k+2) = uout(end);
        cont_vec(k+2) = par.(contPar.Name);
        kappa_vec(k+2) = global_kappa;
        
        wU = uout(1:end-1); 
        break;
    end
        
           
end

if k < contPar.numContSteps % If the continuation finished early 
    
    cont_vec = cont_vec(1:nnz(cont_vec));
    free_vec = free_vec(1:nnz(cont_vec));
    kappa_vec = kappa_vec(1:nnz(cont_vec));  
    
else
    wU = uout(1:end-1);  
end

par.kappa = global_kappa;


% Plot continuation diagram and final pattern
    figure; plot(cont_vec,free_vec,'-','LineWidth',3); % continuation diagram
        xlabel(contPar.Name); ylabel(contPar.Free); 
        set(gca,'fontsize',18,'linewidth',2); box on;
        
plot_full_1Dspiral(wU0,U0,par,numPar,mesh_params,fh_rolls);
 
% Save output
save([file_names.out_file '_final.mat'],'wU','par','mesh_params','numPar','cont_vec','free_vec','kappa_vec','U0')


