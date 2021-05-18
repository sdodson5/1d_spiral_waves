function  plot_full_1Dspiral(u,u_wt,par,numPar,mesh_params,fh_rolls)
% Plots the full 1D spiral pattern: adds core and computed far-field regions

numVars = floor(length(u)/(mesh_params.nt*mesh_params.nz)); % Detect how many parameters

% Calculate far-field solution (core solution is loaded in)
[u_ff0,v_ff0] = fh_rolls(u_wt,par,numPar,mesh_params,mesh_params.Dt,mesh_params.Dt2);
 u_ff0 = [u_ff0(:); v_ff0(:)];

figure; hold on;

for j = 1:numVars
    
    subplot(numVars,1,j);
    U = u((j-1)*mesh_params.nt*mesh_params.nz+1:j*mesh_params.nt*mesh_params.nz);
    
    u_ff = u_ff0((j-1)*mesh_params.nt*mesh_params.nz_long+1:j*mesh_params.nt*mesh_params.nz_long);
    u_ff = mesh_params.chi_ff.*u_ff(mesh_params.bc_idx);
    
    wU = U + u_ff;

    pcolor(mesh_params.xx,mesh_params.tt,reshape(wU,mesh_params.nt,mesh_params.nz));shading interp; colorbar; 
   set(gca,'fontsize',18,'linewidth',2); box on; 
    
    
end


