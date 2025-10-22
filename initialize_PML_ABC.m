
% Initialize PML boundary condition
%
% This version fixes indexing/dimension mismatches and uses a consistent
% angular frequency (w). It keeps the same variable names used elsewhere
% in the codebase (eps_r_xy, eps_r_xz, ..., mu_r_yx, mu_r_zx, ...).
%
% Yee-grid sizes (from initialize_fdfd_material_grid.m):
%   eps_r_x : [nx,   nyp1, nzp1]
%   eps_r_y : [nxp1, ny,   nzp1]
%   eps_r_z : [nxp1, nyp1, nz  ]
%   mu_r_x  : [nxp1, ny,   nz  ]
%   mu_r_y  : [nx,   nyp1, nz  ]
%   mu_r_z  : [nx,   ny,   nzp1]
%
% PML layers:
%   x- (xn): first indices along x;  x+ (xp): last indices along x
%   y- (yn): first indices along y;  y+ (yp): last indices along y
%   z- (zn): first indices along z;  z+ (zp): last indices along z

disp('initializing PML boundary conditions');

% -------------------------------------------------------------------------
% Map base material arrays into the per-orientation aliases used downstream
% (No copies; MATLAB will reference the same data.)
% -------------------------------------------------------------------------
% epsilon aliases
eps_r_xy = eps_r_x;   % [nx,   nyp1, nzp1]
eps_r_xz = eps_r_x;   % [nx,   nyp1, nzp1]  (z-slices will index dim 3)
eps_r_xi = eps_r_x;
eps_r_yx = eps_r_y;   % [nxp1, ny,   nzp1]
eps_r_yz = eps_r_y;   % [nxp1, ny,   nzp1]
eps_r_yi = eps_r_y;
eps_r_zx = eps_r_z;   % [nxp1, nyp1, nz  ]
eps_r_zy = eps_r_z;   % [nxp1, nyp1, nz  ]
eps_r_zi = eps_r_z;
% mu aliases
mu_r_xy  = mu_r_x;    % [nxp1, ny,   nz  ]
mu_r_xz  = mu_r_x;    % [nxp1, ny,   nz  ]
mu_r_xi  = mu_r_x;
mu_r_yx  = mu_r_y;    % [nx,   nyp1, nz  ]
mu_r_yz  = mu_r_y;    % [nx,   nyp1, nz  ]
mu_r_yi  = mu_r_y;
mu_r_zx  = mu_r_z;    % [nx,   ny,   nzp1]
mu_r_zy  = mu_r_z;    % [nx,   ny,   nzp1]
mu_r_zi  = mu_r_z;

% -------------------------------------------------------------------------
% Angular frequency (use frequency_domain if provided)
% -------------------------------------------------------------------------
if isfield(frequency_domain, 'w') && ~isempty(frequency_domain.w)
    w = frequency_domain.w;
else
    % frequency_domain.f should be in Hz
    w = 2*pi*frequency_domain.f;
end

pml_order = boundary.pml_order;
R_0       = boundary.pml_R_0;

% Helper for graded profile (coordinate-normalized distance exponent)
grade = @(rho) rho.^pml_order;

% ============================  X- SIDE  ==================================
if is_pml_xn
    % Max sigma for x-PML (use dx)
    sigma_max = -(pml_order+1)*eps_0*c*log(R_0)/(2*dx*n_pml_xn);
    % Layer-centered distances into the PML (electric & magnetic)
    rho = ((n_pml_xn:-1:1) - 0.5) / n_pml_xn;
    sig_e_vec = sigma_max * grade(rho);          % 1 x n_pml_xn
    sig_m_vec = (mu_0/eps_0) * sig_e_vec;        % impedance-scaled

    for ii = 1:n_pml_xn
        i_eps   = ii;   % eps_r_* arrays with size nxp1 on dim-1
        i_mu_nx = ii;   % mu_r_yx/mu_r_zx have first dim nx

        % eps components touched by x-normal PML (MULTIPLICATIVE UPDATE)
        eps_r_yx(i_eps,:,:) = eps_r_yx(i_eps,:,:) .* (1 + sig_e_vec(ii)/(1i*w*eps_0)); % [ny, nzp1]
        eps_r_zx(i_eps,:,:) = eps_r_zx(i_eps,:,:) .* (1 + sig_e_vec(ii)/(1i*w*eps_0)); % [nyp1, nz]
        % mu components (MULTIPLICATIVE UPDATE)
        mu_r_yx(i_mu_nx,:,:) = mu_r_yx(i_mu_nx,:,:) .* (1 + sig_m_vec(ii)/(1i*w*mu_0)); % [nyp1, nz]
        mu_r_zx(i_mu_nx,:,:) = mu_r_zx(i_mu_nx,:,:) .* (1 + sig_m_vec(ii)/(1i*w*mu_0)); % [ny,   nzp1]
    end
end

% ============================  X+ SIDE  ==================================
if is_pml_xp
    sigma_max = -(pml_order+1)*eps_0*c*log(R_0)/(2*dx*n_pml_xp);
    rho = ((1:n_pml_xp) - 0.5) / n_pml_xp;
    sig_e_vec = sigma_max * grade(rho);
    sig_m_vec = (mu_0/eps_0) * sig_e_vec;

    for ii = 1:n_pml_xp
        i_eps   = nxp1 - n_pml_xp + ii;  % eps_r_yx/eps_r_zx: first dim is nxp1
        i_mu_nx = nx   - n_pml_xp + ii;  % mu_r_yx/mu_r_zx:  first dim is nx

        eps_r_yx(i_eps,:,:) = eps_r_yx(i_eps,:,:) .* (1 + sig_e_vec(ii)/(1i*w*eps_0));
        eps_r_zx(i_eps,:,:) = eps_r_zx(i_eps,:,:) .* (1 + sig_e_vec(ii)/(1i*w*eps_0));
        mu_r_yx(i_mu_nx,:,:) = mu_r_yx(i_mu_nx,:,:) .* (1 + sig_m_vec(ii)/(1i*w*mu_0));
        mu_r_zx(i_mu_nx,:,:) = mu_r_zx(i_mu_nx,:,:) .* (1 + sig_m_vec(ii)/(1i*w*mu_0));
    end
end

% ============================  Y- SIDE  ==================================
if is_pml_yn
    sigma_max = -(pml_order+1)*eps_0*c*log(R_0)/(2*dy*n_pml_yn);
    rho = ((n_pml_yn:-1:1) - 0.5) / n_pml_yn;
    sig_e_vec = sigma_max * grade(rho);
    sig_m_vec = (mu_0/eps_0) * sig_e_vec;

    for jj = 1:n_pml_yn
        j_nyp1 = jj;                   % for arrays with nyp1 on dim-2
        j_ny   = jj;                   % for arrays with ny on dim-2

        eps_r_xy(:,j_nyp1,:) = eps_r_xy(:,j_nyp1,:) .* (1 + sig_e_vec(jj)/(1i*w*eps_0)); % [nx, nzp1]
        eps_r_zy(:,j_nyp1,:) = eps_r_zy(:,j_nyp1,:) .* (1 + sig_e_vec(jj)/(1i*w*eps_0)); % [nxp1, nz]
        mu_r_xy(:,j_ny,:)    = mu_r_xy(:,j_ny,:)    .* (1 + sig_m_vec(jj)/(1i*w*mu_0));  % [nxp1, nz]
        mu_r_zy(:,j_ny,:)    = mu_r_zy(:,j_ny,:)    .* (1 + sig_m_vec(jj)/(1i*w*mu_0));  % [nx,   nzp1]
    end
end

% ============================  Y+ SIDE  ==================================
if is_pml_yp
    sigma_max = -(pml_order+1)*eps_0*c*log(R_0)/(2*dy*n_pml_yp);
    rho = ((1:n_pml_yp) - 0.5) / n_pml_yp;
    sig_e_vec = sigma_max * grade(rho);
    sig_m_vec = (mu_0/eps_0) * sig_e_vec;

    for jj = 1:n_pml_yp
        j_nyp1 = nyp1 - n_pml_yp + jj;  % for nyp1-sized dim-2
        j_ny   = ny   - n_pml_yp + jj;  % for ny-sized dim-2

        eps_r_xy(:,j_nyp1,:) = eps_r_xy(:,j_nyp1,:) .* (1 + sig_e_vec(jj)/(1i*w*eps_0));
        eps_r_zy(:,j_nyp1,:) = eps_r_zy(:,j_nyp1,:) .* (1 + sig_e_vec(jj)/(1i*w*eps_0));
        mu_r_xy(:,j_ny,:)    = mu_r_xy(:,j_ny,:)    .* (1 + sig_m_vec(jj)/(1i*w*mu_0));
        mu_r_zy(:,j_ny,:)    = mu_r_zy(:,j_ny,:)    .* (1 + sig_m_vec(jj)/(1i*w*mu_0));
    end
end

% ============================  Z- SIDE  ==================================
if is_pml_zn
    sigma_max = -(pml_order+1)*eps_0*c*log(R_0)/(2*dz*n_pml_zn);
    rho = ((n_pml_zn:-1:1) - 0.5) / n_pml_zn;
    sig_e_vec = sigma_max * grade(rho);
    sig_m_vec = (mu_0/eps_0) * sig_e_vec;

    for kk = 1:n_pml_zn
        k_nzp1 = kk;                   % eps arrays have nzp1 on dim-3
        k_nz   = kk;                   % mu arrays  have nz   on dim-3

        eps_r_xz(:,:,k_nzp1) = eps_r_xz(:,:,k_nzp1) .* (1 + sig_e_vec(kk)/(1i*w*eps_0)); % [nx,   nyp1]
        eps_r_yz(:,:,k_nzp1) = eps_r_yz(:,:,k_nzp1) .* (1 + sig_e_vec(kk)/(1i*w*eps_0)); % [nxp1, ny  ]
        mu_r_xz(:,:,k_nz)    = mu_r_xz(:,:,k_nz)    .* (1 + sig_m_vec(kk)/(1i*w*mu_0));  % [nxp1, ny  ]
        mu_r_yz(:,:,k_nz)    = mu_r_yz(:,:,k_nz)    .* (1 + sig_m_vec(kk)/(1i*w*mu_0));  % [nx,   nyp1]
    end
end

% ============================  Z+ SIDE  ==================================
if is_pml_zp
    sigma_max = -(pml_order+1)*eps_0*c*log(R_0)/(2*dz*n_pml_zp);
    rho = ((1:n_pml_zp) - 0.5) / n_pml_zp;
    sig_e_vec = sigma_max * grade(rho);
    sig_m_vec = (mu_0/eps_0) * sig_e_vec;

    for kk = 1:n_pml_zp
        k_nzp1 = nzp1 - n_pml_zp + kk;  % eps arrays dim-3
        k_nz   = nz   - n_pml_zp + kk;  % mu arrays  dim-3

        eps_r_xz(:,:,k_nzp1) = eps_r_xz(:,:,k_nzp1) .* (1 + sig_e_vec(kk)/(1i*w*eps_0));
        eps_r_yz(:,:,k_nzp1) = eps_r_yz(:,:,k_nzp1) .* (1 + sig_e_vec(kk)/(1i*w*eps_0));
        mu_r_xz(:,:,k_nz)    = mu_r_xz(:,:,k_nz)    .* (1 + sig_m_vec(kk)/(1i*w*mu_0));
        mu_r_yz(:,:,k_nz)    = mu_r_yz(:,:,k_nz)    .* (1 + sig_m_vec(kk)/(1i*w*mu_0));
    end
end

% No temporary sigma arrays retained; nothing to clear.