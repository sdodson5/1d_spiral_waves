function plot_oneD_spiral_core(u,mesh_params)
% Plots the "core"-region of the 1D spiral
% Plots all components in subfigures

numVars = floor(length(u)/(mesh_params.nt*mesh_params.nz)); % Detect how many parameters

figure; hold on;

for j = 1:numVars
    
    subplot(numVars,1,j);
    U = u((j-1)*mesh_params.nt*mesh_params.nz+1:j*mesh_params.nt*mesh_params.nz);

    pcolor(mesh_params.xx,mesh_params.tt,reshape(U,mesh_params.nt,mesh_params.nz));shading interp; colorbar; 
   set(gca,'fontsize',18,'linewidth',2); box on; 
    
    
end


