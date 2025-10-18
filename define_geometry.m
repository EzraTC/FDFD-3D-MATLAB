disp('defining the problem geometry');

bricks = [];
spheres = [];

% Initialize Substrate
bricks(1).min_x = -1e-3;
bricks(1).min_y = -1e-3;
bricks(1).min_z = -1e-3;
bricks(1).max_x = 1e-3;
bricks(1).max_y = 1e-3;
bricks(1).max_z = 1e-3;
bricks(1).material_type = 4;

