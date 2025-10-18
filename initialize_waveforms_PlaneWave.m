disp('initializing source waveforms');

% === FDFD: Incident Plane Wave ============================

% % Frequency (Hz) - user-input: %Moved to define _problem space parameters
% frequency = 2.4e9;

w = frequency_domain.w;

% Direction (radians) — user-input:
theta_inc = 0;           % polar angle from +z
phi_inc   = pi/4;        % azimuth from +x toward +y

% Polarization in spherical basis:
E_theta = 1;             % amplitude along \hat{theta}
E_phi   = 0;             % amplitude along \hat{phi}

k0 = 2*pi*w / c;

% Unit propagation vector k-hat
kx_hat = sin(theta_inc)*cos(phi_inc);
ky_hat = sin(theta_inc)*sin(phi_inc);
kz_hat = cos(theta_inc);

% Coordinates at cell centers
Xc = fdfd_domain.cell_center_coordinates_x;  % [nx, ny, nz]
Yc = fdfd_domain.cell_center_coordinates_y;  % [nx, ny, nz]
Zc = fdfd_domain.cell_center_coordinates_z;  % [nx, ny, nz]

% Spatial phase factor 
phase = exp(-1j*k0*(kx_hat.*Xc + ky_hat.*Yc + kz_hat.*Zc));

% Cartesian E-field components from (11.17) with f_wf = exp(-jk0 k̂·r)
Ex_coef =  E_theta*cos(theta_inc)*cos(phi_inc) - E_phi*sin(phi_inc);
Ey_coef =  E_theta*cos(theta_inc)*sin(phi_inc) + E_phi*cos(phi_inc);
Ez_coef = -E_theta*sin(theta_inc);

E_inc_center.x = Ex_coef .* phase;   % [nx, ny, nz]
E_inc_center.y = Ey_coef .* phase;   % [nx, ny, nz]
E_inc_center.z = Ez_coef .* phase;   % [nx, ny, nz]

% -------- Add magnetic field: H = (1/eta0) * (khat x E) ----------
H_inc_center.x = (ky_hat .* E_inc_center.z - kz_hat .* E_inc_center.y) / eta0;
H_inc_center.y = (kz_hat .* E_inc_center.x - kx_hat .* E_inc_center.z) / eta0;
H_inc_center.z = (kx_hat .* E_inc_center.y - ky_hat .* E_inc_center.x) / eta0;
% --------------------------------------------------------------------------

% Stash into 'waveforms' 
waveforms.plane_wave_fdfd.theta_inc = theta_inc;
waveforms.plane_wave_fdfd.phi_inc   = phi_inc;
waveforms.plane_wave_fdfd.k0        = k0;
waveforms.plane_wave_fdfd.E_theta   = E_theta;
waveforms.plane_wave_fdfd.E_phi     = E_phi;
waveforms.plane_wave_fdfd.khat      = [kx_hat, ky_hat, kz_hat];
waveforms.plane_wave_fdfd.E_center  = E_inc_center;
waveforms.plane_wave_fdfd.H_center  = H_inc_center;
% =========================================================================

clear w;
