function u_infty = form_full_source_defect_pattern(wU0,par,numPar,mesh_params,fh_rolls)
% Input the core-solution of the source defect/1D spiral
% Output the full solution

global U0;

% Form the far-field rolls
[u_ff0,v_ff0] = fh_rolls(U0,par,numPar,mesh_params,mesh_params.Dt,mesh_params.Dt2); % On the long grid - need to put on short grid
u_ff = mesh_params.chi_ff.*u_ff0(mesh_params.bc_idx);  % put on short grid.
v_ff = mesh_params.chi_ff.*v_ff0(mesh_params.bc_idx); 

wU = wU0(1:mesh_params.nz*mesh_params.nt) + u_ff;
wV = wU0(mesh_params.nz*mesh_params.nt + 1: 2*mesh_params.nz*mesh_params.nt) + v_ff;
u_infty = [wU;wV]; 
