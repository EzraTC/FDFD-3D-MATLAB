%% solve_matrix_equation.m
% Assemble the frequency–domain FDFD matrix system and solve A * x = y, where
% the unknown vector x contains the scattered electric-field components on the
% interior Yee cells and the right-hand side y holds the incident fields. The
% script expects all geometric, material, and coefficient arrays to be prepared
% by the initialization routines that precede it in ``fdfd_solve``.

disp('assembling linear system and solving for scattered electric fields');

%% ------------------------------ Dimensions ---------------------------------
N  = nx * ny * nz;      % Number of unknowns per Cartesian component
NT = 3 * N;             % Total size of the linear system (Ex,Ey,Ez blocks)

% Global ids for each block (flattened with i running fastest)
ex_id = reshape(          1:N, nx, ny, nz);
ey_id = reshape(    N + (1:N), nx, ny, nz);
ez_id = reshape(2 * N + (1:N), nx, ny, nz);

% ------------------- Cross-component neighbor ids (replicate) ----------------
ey_id_jm = cat(2, ey_id(:,1,:),     ey_id(:,1:end-1,:));   % Ey(j-1)
ez_id_km = cat(3, ez_id(:,:,1),     ez_id(:,:,1:end-1));   % Ez(k-1)
ex_id_im = cat(1, ex_id(1,:,:),     ex_id(1:end-1,:,:));   % Ex(i-1)

ez_id_km_forEy = ez_id_km;   % aliases for clarity when building the matrix
ex_id_im_forEy = ex_id_im;
ex_id_im_forEz = ex_id_im;
ey_id_jm_forEz = ey_id_jm;

% ------------------- Same-component neighbor ids (replicate) -----------------
% Ex neighbors along y and z:
ex_id_jm_sc = cat(2, ex_id(:,1,:),     ex_id(:,1:end-1,:));  % j-1
ex_id_jp    = cat(2, ex_id(:,2:end,:), ex_id(:,end,:));      % j+1
ex_id_km_sc = cat(3, ex_id(:,:,1),     ex_id(:,:,1:end-1));  % k-1
ex_id_kp    = cat(3, ex_id(:,:,2:end), ex_id(:,:,end));      % k+1

% Ey neighbors along x and z:
ey_id_im_sc = cat(1, ey_id(1,:,:),     ey_id(1:end-1,:,:));  % i-1
ey_id_ip    = cat(1, ey_id(2:end,:,:), ey_id(end,:,:));      % i+1
ey_id_km_sc = cat(3, ey_id(:,:,1),     ey_id(:,:,1:end-1));  % k-1
ey_id_kp    = cat(3, ey_id(:,:,2:end), ey_id(:,:,end));      % k+1

% Ez neighbors along x and y:
ez_id_im_sc = cat(1, ez_id(1,:,:),     ez_id(1:end-1,:,:));  % i-1
ez_id_ip    = cat(1, ez_id(2:end,:,:), ez_id(end,:,:));      % i+1
ez_id_jm_sc = cat(2, ez_id(:,1,:),     ez_id(:,1:end-1,:));  % j-1
ez_id_jp    = cat(2, ez_id(:,2:end,:), ez_id(:,end,:));      % j+1

%% ---------------------------- Sparse Assembly -------------------------------
% Original Codex terms (diag + cross couplings)  = 15N entries
% Added same-component neighbors (4 per block)   = 12N entries
% Total = 27N nonzeros (before duplicate collation by sparse()).

rows = [
    % ------ Ex rows ------
    ex_id(:);                      % diag
    ex_id(:); ex_id(:);            % to Ey(j), Ey(j-1)
    ex_id(:); ex_id(:);            % to Ez(k), Ez(k-1)
    % ------ Ey rows ------
    ey_id(:);                      % diag
    ey_id(:); ey_id(:);            % to Ex(i), Ex(i-1)
    ey_id(:); ey_id(:);            % to Ez(k), Ez(k-1)
    % ------ Ez rows ------
    ez_id(:);                      % diag
    ez_id(:); ez_id(:);            % to Ex(i), Ex(i-1)
    ez_id(:); ez_id(:);            % to Ey(j), Ey(j-1)
    % ------ Ex same-component neighbors ------
    ex_id(:); ex_id(:); ex_id(:); ex_id(:);    % Ex↔Ex(j-1), Ex(j+1), Ex(k+1), Ex(k-1)
    % ------ Ey same-component neighbors ------
    ey_id(:); ey_id(:); ey_id(:); ey_id(:);    % Ey↔Ey(k-1), Ey(k+1), Ey(i+1), Ey(i-1)
    % ------ Ez same-component neighbors ------
    ez_id(:); ez_id(:); ez_id(:); ez_id(:)     % Ez↔Ez(i-1), Ez(i+1), Ez(j+1), Ez(j-1)
];

cols = [
    % ------ Ex block couplings ------
    ex_id(:);
    ey_id(:);        ey_id_jm(:);
    ez_id(:);        ez_id_km(:);
    % ------ Ey block couplings ------
    ey_id(:);
    ex_id(:);        ex_id_im_forEy(:);
    ez_id(:);        ez_id_km_forEy(:);
    % ------ Ez block couplings ------
    ez_id(:);
    ex_id(:);        ex_id_im_forEz(:);
    ey_id(:);        ey_id_jm_forEz(:);
    % ------ Ex same-component neighbors ------
    ex_id_jm_sc(:);  ex_id_jp(:);   ex_id_kp(:);  ex_id_km_sc(:);
    % ------ Ey same-component neighbors ------
    ey_id_km_sc(:);  ey_id_kp(:);   ey_id_ip(:);  ey_id_im_sc(:);
    % ------ Ez same-component neighbors ------
    ez_id_im_sc(:);  ez_id_ip(:);   ez_id_jp(:);  ez_id_jm_sc(:)
];

vals = [
    % ------ Ex block couplings ------
    C_sxx5(:);
   -C_sxy1(:);      -C_sxy2(:);
   -C_sxz1(:);      -C_sxz2(:);
    % ------ Ey block couplings ------
    C_syy5(:);
   -C_syx1(:);      -C_syx2(:);
   -C_syz1(:);      -C_syz2(:);
    % ------ Ez block couplings ------
    C_szz5(:);
   -C_szx1(:);      -C_szx2(:);
   -C_szy1(:);      -C_szy2(:);
    % ------ Ex same-component neighbors ------
   -C_sxx1(:);      -C_sxx2(:);    -C_sxx3(:);    -C_sxx4(:);
    % ------ Ey same-component neighbors ------
   -C_syy1(:);      -C_syy2(:);    -C_syy3(:);    -C_syy4(:);
    % ------ Ez same-component neighbors ------
   -C_szz1(:);      -C_szz2(:);    -C_szz3(:);    -C_szz4(:)
];

A = sparse(rows, cols, vals, NT, NT);

%% --------------------- Incident-Field Right-Hand Side -----------------------
% Build the incident fields on the Yee electric grids if they are not already
% present (plane-wave excitation defined in initialize_waveforms).
if ~exist('waveforms', 'var') || ~isfield(waveforms, 'plane_wave')
    error('Plane-wave definition (waveforms.plane_wave) is required before solving.');
end

theta_inc = waveforms.plane_wave.theta_inc;
phi_inc   = waveforms.plane_wave.phi_inc;
E_theta   = waveforms.plane_wave.E_theta;
E_phi     = waveforms.plane_wave.E_phi;
k0        = waveforms.plane_wave.k0;
khat      = waveforms.plane_wave.khat(:).';

Ex_coef =  E_theta*cos(theta_inc)*cos(phi_inc) - E_phi*sin(phi_inc);
Ey_coef =  E_theta*cos(theta_inc)*sin(phi_inc) + E_phi*cos(phi_inc);
Ez_coef = -E_theta*sin(theta_inc);

% Coordinate vectors for the three electric-field staggered grids
x_ex = fdfd_domain.min_x + ((0:nx-1) + 0.5) * dx;
y_ex = fdfd_domain.min_y + (0:nyp1-1)       * dy;
z_ex = fdfd_domain.min_z + (0:nzp1-1)       * dz;
[xEx, yEx, zEx] = ndgrid(x_ex, y_ex, z_ex);
E_inc_x = Ex_coef * exp(-1i * k0 * (khat(1) * xEx + khat(2) * yEx + khat(3) * zEx));

x_ey = fdfd_domain.min_x + (0:nxp1-1)       * dx;
y_ey = fdfd_domain.min_y + ((0:ny-1) + 0.5) * dy;
z_ey = fdfd_domain.min_z + (0:nzp1-1)       * dz;
[xEy, yEy, zEy] = ndgrid(x_ey, y_ey, z_ey);
E_inc_y = Ey_coef * exp(-1i * k0 * (khat(1) * xEy + khat(2) * yEy + khat(3) * zEy));

x_ez = fdfd_domain.min_x + (0:nxp1-1)       * dx;
y_ez = fdfd_domain.min_y + (0:nyp1-1)       * dy;
z_ez = fdfd_domain.min_z + ((0:nz-1) + 0.5) * dz;
[xEz, yEz, zEz] = ndgrid(x_ez, y_ez, z_ez);
E_inc_z = Ez_coef * exp(-1i * k0 * (khat(1) * xEz + khat(2) * yEz + khat(3) * zEz));

% Extract the interior slices that correspond to unknowns (nx×ny×nz)
Ex_inc_int = E_inc_x(:, 2:nyp1, 2:nzp1);
Ey_inc_int = E_inc_y(2:nxp1, :, 2:nzp1);
Ez_inc_int = E_inc_z(2:nxp1, 2:nyp1, :);

y = [Ex_inc_int(:); Ey_inc_int(:); Ez_inc_int(:)];

%% ------------------------------ Solve System --------------------------------
x = A \ y;

% Reshape the scattered-field solution back onto Yee-cell interiors
Ex_scat = reshape(x(1:N),       nx, ny, nz);
Ey_scat = reshape(x(N+1:2*N),   nx, ny, nz);
Ez_scat = reshape(x(2*N+1:end), nx, ny, nz);

% Store scattered, incident, and total fields for later post-processing
E_scattered.x = Ex_scat;   E_incident.x = Ex_inc_int;   E_total.x = Ex_inc_int + Ex_scat;
E_scattered.y = Ey_scat;   E_incident.y = Ey_inc_int;   E_total.y = Ey_inc_int + Ey_scat;
E_scattered.z = Ez_scat;   E_incident.z = Ez_inc_int;   E_total.z = Ez_inc_int + Ez_scat;

% Insert total fields back into the Yee arrays (interior indices only)
Ex(:, 2:nyp1, 2:nzp1) = E_total.x;
Ey(2:nxp1, :, 2:nzp1) = E_total.y;
Ez(2:nxp1, 2:nyp1, :) = E_total.z;

clear xEx yEx zEx xEy yEy zEy xEz yEz zEz;
