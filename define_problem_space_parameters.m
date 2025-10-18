disp('defining the problem space parameters');

% A factor determining the accuracy limit of FDFD results
number_of_cells_per_wavelength = 20;    %Do not know if this number will work

%The frequency in which this is running at 
frequency_domain.f = 5e9;

% Dimensions of a unit cell in x, y, and z directions (meters)
dx=1e-3;    
dy=1e-3;   
dz=1e-3;

% ==<boundary conditions>========
% Here we define the boundary conditions parameters 
% 'pec' : perfect electric conductor
% 'cpml': conlvolutional PML
% 'pml' : PML for FDFD
% if cpml_number_of_cells is less than zero
% CPML extends inside of the domain rather than outwards

% boundary.type_xn = 'cpml';
% boundary.air_buffer_number_of_cells_xn = 8;
% boundary.cpml_number_of_cells_xn = 8;
% 
% boundary.type_xp = 'cpml';
% boundary.air_buffer_number_of_cells_xp = 8;
% boundary.cpml_number_of_cells_xp = 8;
% 
% boundary.type_yn = 'cpml';
% boundary.air_buffer_number_of_cells_yn = 8;
% boundary.cpml_number_of_cells_yn = 8;
% 
% boundary.type_yp = 'cpml';
% boundary.air_buffer_number_of_cells_yp = 8;
% boundary.cpml_number_of_cells_yp = 8;
% 
% boundary.type_zn = 'cpml';
% boundary.air_buffer_number_of_cells_zn = 8;
% boundary.cpml_number_of_cells_zn = 8;
% 
% boundary.type_zp = 'cpml';
% boundary.air_buffer_number_of_cells_zp = 8;
% boundary.cpml_number_of_cells_zp = 8;
% 
% boundary.cpml_order = 3; 
% boundary.cpml_sigma_factor = 1.3;
% boundary.cpml_kappa_max = 7;
% boundary.cpml_alpha_min = 0;
% boundary.cpml_alpha_max = 0.05;

boundary.type_xn = 'pml';
boundary.air_buffer_number_of_cells_xn = 2;
boundary.pml_number_of_cells_xn = 4;

boundary.type_xp = 'pml';
boundary.air_buffer_number_of_cells_xp = 2;
boundary.pml_number_of_cells_xp = 4;

boundary.type_yn = 'pml';
boundary.air_buffer_number_of_cells_yn = 2;
boundary.pml_number_of_cells_yn = 4;

boundary.type_yp = 'pml';
boundary.air_buffer_number_of_cells_yp = 2;
boundary.pml_number_of_cells_yp = 4;

boundary.type_zn = 'pml';
boundary.air_buffer_number_of_cells_zn = 2;
boundary.pml_number_of_cells_zn = 4;

boundary.type_zp = 'pml';
boundary.air_buffer_number_of_cells_zp = 2;
boundary.pml_number_of_cells_zp = 4;

boundary.pml_order = 2;
boundary.pml_R_0 = 1e-8;

% ===<material types>============
% Here we define and initialize the arrays of material types
% eps_r   : relative permittivity
% mu_r    : relative permeability
% sigma_e : electric conductivity
% sigma_m : magnetic conductivity

% air
material_types(1).eps_r   = 1;
material_types(1).mu_r    = 1;
material_types(1).sigma_e = 0;
material_types(1).sigma_m = 0; 
material_types(1).color   = [1 1 1];

% PEC : perfect electric conductor
material_types(2).eps_r   = 1;
material_types(2).mu_r    = 1;
material_types(2).sigma_e = 1e10;
material_types(2).sigma_m = 0; 
material_types(2).color   = [1 0 0]; % R G B

% PMC : perfect magnetic conductor
material_types(3).eps_r   = 1;
material_types(3).mu_r    = 1;
material_types(3).sigma_e = 0;
material_types(3).sigma_m = 1e10;
material_types(3).color   = [0 1 0]; % R G B

% a dielectric 
material_types(4).eps_r   = 4;
material_types(4).mu_r    = 1;
material_types(4).sigma_e = 0;
material_types(4).sigma_m = 0; 
material_types(4).color   = [0 0 1]; %R G B

% index of material types defining air, PEC, and PMC 
material_type_index_air = 1;
material_type_index_pec = 2;
material_type_index_pmc = 3;
