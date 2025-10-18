disp('defining output parameters');

%Need to add the FDFD specfic items to this code
sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];

% figure refresh rate
plotting_step = 100;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% define animation
% field_type shall be 'e' or 'h'
% plane cut shall be 'xy', yz, or zx
% component shall be 'x', 'y', 'z', or 'm;
animation(1).field_type = 'e';
animation(1).component = 'x';
animation(1).plane_cut(1).type = 'xy';
animation(1).plane_cut(1).position  = 1e-3;
animation(1).enable = false;
animation(1).display_grid = false;
animation(1).display_objects = true;

% display problem space parameters
problem_space_display.labels = false;
problem_space_display.axis_at_origin = false;
problem_space_display.axis_outside_domain = true;
problem_space_display.grid_xn = false;
problem_space_display.grid_xp = false;
problem_space_display.grid_yn = false;
problem_space_display.grid_yp = false;
problem_space_display.grid_zn = false;
problem_space_display.grid_zp = false;
problem_space_display.outer_boundaries = true;
problem_space_display.pml_boundaries = true;

