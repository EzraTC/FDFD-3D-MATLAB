disp('defining output parameters');

% FDFD: clear sample requests (add any you need later)
sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];

% figure refresh rate
plotting_step = 100;

% mode of operation
run_simulation      = true;
show_material_mesh  = true;
show_problem_space  = true;

% ------------------ Define animation (FDFD phasor) ------------------
% field_type shall be 'e' or 'h'
% plane_cut shall be 'xy', 'yz', or 'zx'
% component shall be 'x', 'y', or 'z'  (use 'm' only if your plotter supports magnitude)
animation(1).field_type        = 'e';
animation(1).component         = 'x';           % 'x' | 'y' | 'z'
animation(1).plane_cut(1).type = 'xy';
animation(1).enable            = true;
animation(1).display_grid      = false;
animation(1).display_objects   = true;

% FDFD phasor animation settings (synthesize frames from a single frequency)
animation(1).mode        = 'phasor';            % tell the driver this is FDFD
animation(1).num_frames  = 40;                  % phase steps over 0..2π
animation(1).phase_list  = linspace(0, 2*pi, animation(1).num_frames+1);
animation(1).phase_list(end) = [];              % drop duplicate 2π endpoint
animation(1).frequency   = frequency_domain.f;  % labeling only
animation(1).omega       = 2*pi*frequency_domain.f;

% NOTE: The plotter should generate each frame as:
%   frame_field = real(E_component .* exp(1i * animation(1).phase_list(m)));
% Do NOT put that line here; this file only defines parameters.

% ------------------ Choose & align the plane position ------------------
desired_z = 1e-3;   % meters (pick your preferred cut location)

% Safely read PML thicknesses if present (avoid referencing undefined vars)
npml_xn = 0; if exist('n_pml_xn','var'), npml_xn = n_pml_xn; end
npml_xp = 0; if exist('n_pml_xp','var'), npml_xp = n_pml_xp; end
npml_yn = 0; if exist('n_pml_yn','var'), npml_yn = n_pml_yn; end
npml_yp = 0; if exist('n_pml_yp','var'), npml_yp = n_pml_yp; end
npml_zn = 0; if exist('n_pml_zn','var'), npml_zn = n_pml_zn; end
npml_zp = 0; if exist('n_pml_zp','var'), npml_zp = n_pml_zp; end

% Valid k ranges (faces vs. centers) excluding PML
k_min_face = 1    + npml_zn;          % z-faces index: 1..nzp1
k_max_face = nzp1 - npml_zp;
k_min_ctr  = 1    + npml_zn;          % z-centers index: 1..nz
k_max_ctr  = nz   - npml_zp;

switch lower(animation(1).component)
    case 'x'  % Ex lies on z-faces: z = min_z + (k-1)*dz, k = 1..nzp1
        k = round((desired_z - fdfd_domain.min_z)/dz) + 1;
        k = max(k_min_face, min(k, k_max_face));
        animation(1).plane_cut(1).position = fdfd_domain.min_z + (k-1)*dz;

    case 'y'  % Ey lies on z-faces: z = min_z + (k-1)*dz, k = 1..nzp1
        k = round((desired_z - fdfd_domain.min_z)/dz) + 1;
        k = max(k_min_face, min(k, k_max_face));
        animation(1).plane_cut(1).position = fdfd_domain.min_z + (k-1)*dz;

    case 'z'  % Ez lies at z-centers: z = min_z + (k-0.5)*dz, k = 1..nz
        k = round((desired_z - fdfd_domain.min_z)/dz + 0.5);
        k = max(k_min_ctr, min(k, k_max_ctr));
        animation(1).plane_cut(1).position = fdfd_domain.min_z + (k-0.5)*dz;

    otherwise
        error('animation: unknown component "%s". Use x, y, or z.', animation(1).component);
end

% ------------------ Problem-space display options ------------------
problem_space_display.labels               = false;
problem_space_display.axis_at_origin       = false;
problem_space_display.axis_outside_domain  = true;
problem_space_display.grid_xn              = false;
problem_space_display.grid_xp              = false;
problem_space_display.grid_yn              = false;
problem_space_display.grid_yp              = false;
problem_space_display.grid_zn              = false;
problem_space_display.grid_zp              = false;
problem_space_display.outer_boundaries     = true;
problem_space_display.pml_boundaries       = true;
