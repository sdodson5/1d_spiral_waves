function M = construct_M_pulse(qs,delta,par,numPar,mesh_params)

global Us; global Uf;
global Uadjefcn; 

options = optimset('Display','off','Jacobian','on', 'DerivativeCheck','off',...
                     'TolX',1.e-6,'TolFun',1.e-6,'MaxIter',200);

tmp1 = flip(Us(1:mesh_params.nx));  % Need new slow wave to propagate in opposite direction
tmp2 = flip(Us(mesh_params.nx+1:2*mesh_params.nx));

tmp3 = Uf(1:mesh_params.nx);  % Fast wave
tmp4 = Uf(mesh_params.nx+1:2*mesh_params.nx);

%% Interpolate onto x-t mesh
x = mesh_params.xx(1,:)';

% Locations of fast and slow pulse
xs = x(end)/2 - qs; 
xf = x(end)/2 - (qs + delta);

% Fit slow wave to a function and shift
myFit_1 = fit(x,tmp1,'smoothingspline');
myFit_2 = fit(x,tmp2,'smoothingspline');

U1 = [myFit_1(x + xs); myFit_2(x + xs)];

uout = fsolve(@(y) Morris_Lecar_1D_adj_eigenvalue(y,[U1;-2*par.cs],mesh_params.L1x,mesh_params.L2x,par,numPar),[Uadjefcn;par.lambda],options); % Adjoint eigenfunction is centered for slow wave of M
Uadjefcn = uout(1:end-1);

% Fit fast wave to a function and shift
myFit_1 = fit(x,tmp3,'smoothingspline');
myFit_2 = fit(x,tmp4,'smoothingspline');

U2 = [myFit_1(x + xf); myFit_2(x + xf)];

M = [U1(1:mesh_params.nx) + U2(1:mesh_params.nx) - U2(1); U1(mesh_params.nx+1:2*mesh_params.nx) + U2(mesh_params.nx+1:2*mesh_params.nx) - U1(end)];
