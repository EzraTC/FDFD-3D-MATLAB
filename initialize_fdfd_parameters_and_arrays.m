disp('initializing FDFD parameters and arrays');

% constant parameters
eps_0 = 8.854187817e-12; % permittivity of free space               
mu_0  = 4*pi*1e-7; % permeability of free space                    
c = 1/sqrt(mu_0*eps_0); % speed of light in free space
eta0 = sqrt(mu_0/eps_0); %intrinsic impedance of free space
% Create and initialize field and current arrays
disp('creating field arrays');

Hx = zeros(nxp1,ny,nz);   
Hy = zeros(nx,nyp1,nz);   
Hz = zeros(nx,ny,nzp1);
Ex = zeros(nx,nyp1,nzp1);
Ey = zeros(nxp1,ny,nzp1);
Ez = zeros(nxp1,nyp1,nz);
