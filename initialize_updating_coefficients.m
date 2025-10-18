disp('initializing general updating coefficients');

% ============================ Angular frequency ===========================
if isfield(frequency_domain,'w') && ~isempty(frequency_domain.w)
    w = frequency_domain.w;                 % rad/s
else
    w = 2*pi*frequency_domain.f;            % rad/s
end

% ============================ Interior windows ============================
% Unknowns live on the interior Yee indices:
% Ex: (i=1..nx,   j=2..nyp1 -> ny, k=2..nzp1 -> nz)
% Ey: (i=2..nxp1 -> nx,     j=1..ny,          k=2..nzp1 -> nz)
% Ez: (i=2..nxp1 -> nx,     j=2..nyp1 -> ny,  k=1..nz)
Jc = 2:nyp1;       % maps to ny points along y
Kc = 2:nzp1;       % maps to nz points along z
Ic = 2:nxp1;       % maps to nx points along x

% =========================================================================
%                    FULL-SIZE, PML-AWARE μ (by orientation)
% =========================================================================
% Base alias sizes (after initialize_PML_ABC):
%   mu_r_x  : [nxp1, ny,   nz  ]
%   mu_r_y  : [nx,   nyp1, nz  ]
%   mu_r_z  : [nx,   ny,   nzp1]
%
% Directional aliases (these are the ones PML modified):
%   mu_r_zy : [nx,   ny,   nzp1]   % μ_z modified along y
%   mu_r_yz : [nx,   nyp1, nz  ]   % μ_y modified along z
%   mu_r_xz : [nxp1, ny,   nz  ]   % μ_x modified along z
%   mu_r_zx : [nx,   ny,   nzp1]   % μ_z modified along x
%   mu_r_yx : [nx,   nyp1, nz  ]   % μ_y modified along x
%   mu_r_xy : [nxp1, ny,   nz  ]   % μ_x modified along y
%
% We expand/shift these to match each E-grid with replication on the outer
% layer so arrays remain full-size.

%% ----------------- For Ex grid: size [nx, nyp1, nzp1] ---------------------

% ZY initialization
mu_zy__Ex_j  = zeros(nx, nyp1, nzp1);  
mu_zy__Ex_jm = zeros(nx, nyp1, nzp1); 

% ZX initialization
mu_zx__Ex_j  = zeros(nx, nyp1, nzp1);  
mu_zx__Ex_jm = zeros(nx, nyp1, nzp1);  

% YZ initialization
mu_yz__Ex_k  = zeros(nx, nyp1, nzp1); 
mu_yz__Ex_km = zeros(nx, nyp1, nzp1);  

% YX initialization
mu_yx__Ex_k  = zeros(nx, nyp1, nzp1);  
mu_yx__Ex_km = zeros(nx, nyp1, nzp1);  

% Zi initialization
mu_zi__Ex_j  = zeros(nx, nyp1, nzp1);  
mu_zi__Ex_jm = zeros(nx, nyp1, nzp1); 

% Yi initialization
mu_yi__Ex_k  = zeros(nx, nyp1, nzp1);  
mu_yi__Ex_km = zeros(nx, nyp1, nzp1);   

% ZY
mu_zy__Ex_j(:,1:ny,:)    = mu_r_zy(:,1:ny,:);
mu_zy__Ex_j(:,nyp1,:)    = mu_r_zy(:,ny,:);
mu_zy__Ex_jm(:,2:nyp1,:) = mu_r_zy(:,1:ny,:);
mu_zy__Ex_jm(:,1,:)      = mu_r_zy(:,1,:);

% ZX
mu_zx__Ex_j(:,1:ny,:)    = mu_r_zx(:,1:ny,:);
mu_zx__Ex_j(:,nyp1,:)    = mu_r_zx(:,ny,:);
mu_zx__Ex_jm(:,2:nyp1,:) = mu_r_zx(:,1:ny,:);
mu_zx__Ex_jm(:,1,:)      = mu_r_zx(:,1,:);

% YZ
mu_yz__Ex_k(:,:,1:nz)    = mu_r_yz(:,:,1:nz);
mu_yz__Ex_k(:,:,nzp1)    = mu_r_yz(:,:,nz);
mu_yz__Ex_km(:,:,2:nzp1) = mu_r_yz(:,:,1:nz);
mu_yz__Ex_km(:,:,1)      = mu_r_yz(:,:,1);

% YX
mu_yx__Ex_k(:,:,1:nz)    = mu_r_yx(:,:,1:nz);
mu_yx__Ex_k(:,:,nzp1)    = mu_r_yx(:,:,nz);
mu_yx__Ex_km(:,:,2:nzp1) = mu_r_yx(:,:,1:nz);
mu_yx__Ex_km(:,:,1)      = mu_r_yx(:,:,1);

%Zi
mu_zi__Ex_j(:,1:ny,:)    = mu_r_zi(:,1:ny,:);
mu_zi__Ex_j(:,nyp1,:)    = mu_r_zi(:,ny,:);
mu_zi__Ex_jm(:,2:nyp1,:) = mu_r_zi(:,1:ny,:);
mu_zi__Ex_jm(:,1,:)      = mu_r_zi(:,1,:);

% Yi
mu_yi__Ex_k(:,:,1:nz)    = mu_r_yi(:,:,1:nz);
mu_yi__Ex_k(:,:,nzp1)    = mu_r_yi(:,:,nz);
mu_yi__Ex_km(:,:,2:nzp1) = mu_r_yi(:,:,1:nz);
mu_yi__Ex_km(:,:,1)      = mu_r_yi(:,:,1);

%% ----------------- For Ey grid: size [nxp1, ny, nzp1] ---------------------
% XZ initialize
mu_xz__Ey_k  = zeros(nxp1, ny, nzp1);  
mu_xz__Ey_km = zeros(nxp1, ny, nzp1);  

% XY initialize
mu_xy__Ey_k  = zeros(nxp1, ny, nzp1);  
mu_xy__Ey_km = zeros(nxp1, ny, nzp1);  

% ZX initialize
mu_zx__Ey_i  = zeros(nxp1, ny, nzp1);  
mu_zx__Ey_im = zeros(nxp1, ny, nzp1);  

% ZY initialize
mu_zy__Ey_i  = zeros(nxp1, ny, nzp1);  
mu_zy__Ey_im = zeros(nxp1, ny, nzp1);  

% Xi initialize
mu_xi__Ey_k  = zeros(nxp1, ny, nzp1);  
mu_xi__Ey_km = zeros(nxp1, ny, nzp1);  

% Zi initialize
mu_zi__Ey_k  = zeros(nxp1, ny, nzp1);  
mu_zi__Ey_km = zeros(nxp1, ny, nzp1); 

% XZ
mu_xz__Ey_k(:,:,1:nz)    = mu_r_xz(:,:,1:nz);
mu_xz__Ey_k(:,:,nzp1)    = mu_r_xz(:,:,nz);
mu_xz__Ey_km(:,:,2:nzp1) = mu_r_xz(:,:,1:nz);
mu_xz__Ey_km(:,:,1)      = mu_r_xz(:,:,1);

% XY
mu_xy__Ey_k(:,:,1:nz)    = mu_r_xy(:,:,1:nz);
mu_xy__Ey_k(:,:,nzp1)    = mu_r_xy(:,:,nz);
mu_xy__Ey_km(:,:,2:nzp1) = mu_r_xy(:,:,1:nz);
mu_xy__Ey_km(:,:,1)      = mu_r_xy(:,:,1);

% ZX
mu_zx__Ey_i(1:nx,:,:)    = mu_r_zx(1:nx,:,:);
mu_zx__Ey_i(nxp1,:,:)    = mu_r_zx(nx,:,:);
mu_zx__Ey_im(2:nxp1,:,:) = mu_r_zx(1:nx,:,:);
mu_zx__Ey_im(1,:,:)      = mu_r_zx(1,:,:);

% ZY
mu_zy__Ey_i(1:nx,:,:)    = mu_r_zy(1:nx,:,:);
mu_zy__Ey_i(nxp1,:,:)    = mu_r_zy(nx,:,:);
mu_zy__Ey_im(2:nxp1,:,:) = mu_r_zy(1:nx,:,:);
mu_zy__Ey_im(1,:,:)      = mu_r_zy(1,:,:);

% Xi
mu_xi__Ey_k(:,:,1:nz)    = mu_r_xi(:,:,1:nz);
mu_xi__Ey_k(:,:,nzp1)    = mu_r_xi(:,:,nz);
mu_xi__Ey_km(:,:,2:nzp1) = mu_r_xi(:,:,1:nz);
mu_xi__Ey_km(:,:,1)      = mu_r_xi(:,:,1);

% Zi
mu_zi__Ey_i(1:nx,:,:)    = mu_r_zi(1:nx,:,:);
mu_zi__Ey_i(nxp1,:,:)    = mu_r_zi(nx,:,:);
mu_zi__Ey_im(2:nxp1,:,:) = mu_r_zi(1:nx,:,:);
mu_zi__Ey_im(1,:,:)      = mu_r_zi(1,:,:);

%% ----------------- For Ez grid: size [nxp1, nyp1, nz] ---------------------

% YX initialize
mu_yx__Ez_i  = zeros(nxp1, nyp1, nz);  
mu_yx__Ez_im = zeros(nxp1, nyp1, nz); 

% XY initialize 
mu_xy__Ez_j  = zeros(nxp1, nyp1, nz); 
mu_xy__Ez_jm = zeros(nxp1, nyp1, nz);  

% YZ initialize
mu_yz__Ez_i  = zeros(nxp1, nyp1, nz);  
mu_yz__Ez_im = zeros(nxp1, nyp1, nz);

% XZ initialize 
mu_xz__Ez_j  = zeros(nxp1, nyp1, nz); 
mu_xz__Ez_jm = zeros(nxp1, nyp1, nz);  

% Yi initialize
mu_yi__Ez_i  = zeros(nxp1, nyp1, nz);  
mu_yi__Ez_im = zeros(nxp1, nyp1, nz);

% Xi initialize 
mu_xi__Ez_j  = zeros(nxp1, nyp1, nz); 
mu_xi__Ez_jm = zeros(nxp1, nyp1, nz); 

% YX
mu_yx__Ez_i(1:nx,:,:)    = mu_r_yx(1:nx,:,:);
mu_yx__Ez_i(nxp1,:,:)    = mu_r_yx(nx,:,:);
mu_yx__Ez_im(2:nxp1,:,:) = mu_r_yx(1:nx,:,:);
mu_yx__Ez_im(1,:,:)      = mu_r_yx(1,:,:);

% XY
mu_xy__Ez_j(:,1:ny,:)    = mu_r_xy(:,1:ny,:);
mu_xy__Ez_j(:,nyp1,:)    = mu_r_xy(:,ny,:);
mu_xy__Ez_jm(:,2:nyp1,:) = mu_r_xy(:,1:ny,:);
mu_xy__Ez_jm(:,1,:)      = mu_r_xy(:,1,:);

% YZ
mu_yz__Ez_i(1:nx,:,:)    = mu_r_yz(1:nx,:,:);
mu_yz__Ez_i(nxp1,:,:)    = mu_r_yz(nx,:,:);
mu_yz__Ez_im(2:nxp1,:,:) = mu_r_yz(1:nx,:,:);
mu_yz__Ez_im(1,:,:)      = mu_r_yz(1,:,:);

% XZ
mu_xz__Ez_j(:,1:ny,:)    = mu_r_xz(:,1:ny,:);
mu_xz__Ez_j(:,nyp1,:)    = mu_r_xz(:,ny,:);
mu_xz__Ez_jm(:,2:nyp1,:) = mu_r_xz(:,1:ny,:);
mu_xz__Ez_jm(:,1,:)      = mu_r_xz(:,1,:);

% Yi
mu_yi__Ez_i(1:nx,:,:)    = mu_r_yi(1:nx,:,:);
mu_yi__Ez_i(nxp1,:,:)    = mu_r_yi(nx,:,:);
mu_yi__Ez_im(2:nxp1,:,:) = mu_r_yi(1:nx,:,:);
mu_yi__Ez_im(1,:,:)      = mu_r_yi(1,:,:);

% Xi
mu_xi__Ez_j(:,1:ny,:)    = mu_r_xi(:,1:ny,:);
mu_xi__Ez_j(:,nyp1,:)    = mu_r_xi(:,ny,:);
mu_xi__Ez_jm(:,2:nyp1,:) = mu_r_xi(:,1:ny,:);
mu_xi__Ez_jm(:,1,:)      = mu_r_xi(:,1,:);


% =========================================================================
%                       PML-AWARE ε (by orientation)
% =========================================================================
% Ex uses: eps_r_xy (y-derivative terms), eps_r_xz (z-derivative terms)
% Ey uses: eps_r_yz (z-derivative terms), eps_r_yx (x-derivative terms)
% Ez uses: eps_r_zx (x-derivative terms), eps_r_zy (y-derivative terms)

%% Epselon's and Mu's used in Coefficients

% Ex:
Eps_x_xy = eps_r_xy(1:nx, Jc,   Kc);   % εx for y-coupled Ex terms
Eps_x_xz = eps_r_xz(1:nx, Jc,   Kc);   % εx for z-coupled Ex terms
Mu_zx_j   = mu_zx__Ex_j(1:nx,     Jc,   Kc);
Mu_zx_jm  = mu_zx__Ex_jm(1:nx,    Jc,   Kc);
Mu_zy_j   = mu_zy__Ex_j(1:nx,     Jc,   Kc);
Mu_zy_jm  = mu_zy__Ex_jm(1:nx,    Jc,   Kc);
Mu_yz_k   = mu_yz__Ex_k(1:nx,     Jc,   Kc);
Mu_yz_km  = mu_yz__Ex_km(1:nx,    Jc,   Kc);
Mu_yx_k   = mu_yx__Ex_k(1:nx,     Jc,   Kc);
Mu_yx_km  = mu_yx__Ex_km(1:nx,    Jc,   Kc);
Mu_zi_j   = mu_zi__Ex_j(1:nx,     Jc,   Kc);
Mu_zi_jm  = mu_zi__Ex_jm(1:nx,    Jc,   Kc);
Mu_yi_k   = mu_yi__Ex_k(1:nx,     Jc,   Kc);
Mu_yi_km  = mu_yi__Ex_km(1:nx,    Jc,   Kc);

% Ey:
Eps_y_yz = eps_r_yz(Ic,  1:ny,  Kc);   % εy for z-coupled Ey terms
Eps_y_yx = eps_r_yx(Ic,  1:ny,  Kc);   % εy for x-coupled Ey terms
Mu_xz_k   = mu_xz__Ey_k(Ic,        1:ny,  Kc);
Mu_xz_km  = mu_xz__Ey_km(Ic,       1:ny,  Kc);
Mu_zx_i   = mu_zx__Ey_i(Ic,        1:ny,  Kc);
Mu_zx_im  = mu_zx__Ey_im(Ic,       1:ny,  Kc);
Mu_xy_k   = mu_xy__Ey_k(Ic,        1:ny,  Kc);
Mu_xy_km  = mu_xy__Ey_km(Ic,       1:ny,  Kc);
Mu_zy_i   = mu_zy__Ey_i(Ic,        1:ny,  Kc);
Mu_zy_im  = mu_zy__Ey_im(Ic,       1:ny,  Kc);
Mu_xi_k   = mu_xi__Ey_k(Ic,        1:ny,  Kc);
Mu_xi_km  = mu_xi__Ey_km(Ic,       1:ny,  Kc);
Mu_zi_i   = mu_zi__Ey_i(Ic,        1:ny,  Kc);
Mu_zi_im  = mu_zi__Ey_im(Ic,       1:ny,  Kc);

% Ez:
Eps_z_zx = eps_r_zx(Ic,  Jc,     1:nz); % εz for x-coupled Ez terms
Eps_z_zy = eps_r_zy(Ic,  Jc,     1:nz); % εz for y-coupled Ez terms
Mu_yx_i   = mu_yx__Ez_i(Ic,        Jc,    1:nz);
Mu_yx_im  = mu_yx__Ez_im(Ic,       Jc,    1:nz);
Mu_xy_j   = mu_xy__Ez_j(Ic,        Jc,    1:nz);
Mu_xy_jm  = mu_xy__Ez_jm(Ic,       Jc,    1:nz);
Mu_yz_i   = mu_yz__Ez_i(Ic,        Jc,    1:nz);
Mu_yz_im  = mu_yz__Ez_im(Ic,       Jc,    1:nz);
Mu_xz_j   = mu_xz__Ez_j(Ic,        Jc,    1:nz);
Mu_xz_jm  = mu_xz__Ez_jm(Ic,       Jc,    1:nz);
Mu_yi_i   = mu_yi__Ez_i(Ic,        Jc,    1:nz);
Mu_yi_im  = mu_yi__Ez_im(Ic,       Jc,    1:nz);
Mu_xi_j   = mu_xi__Ez_j(Ic,        Jc,    1:nz);
Mu_xi_jm  = mu_xi__Ez_jm(Ic,       Jc,    1:nz);


%% =========================================================================
%                          COEFFICIENTS (nx×ny×nz)
% ==========================================================================

% --------------------------- EX coefficients ------------------------------
C_sxy1 = 1 ./ ( w.^2 .* Eps_x_xy .* eps_0 .* Mu_zx_j  .* mu_0 .* dx .* dy );
C_sxy2 = 1 ./ ( w.^2 .* Eps_x_xy .* eps_0 .* Mu_zx_jm .* mu_0 .* dx .* dy );

C_sxz1 = 1 ./ ( w.^2 .* Eps_x_xz .* eps_0 .* Mu_yx_k  .* mu_0 .* dx .* dz );
C_sxz2 = 1 ./ ( w.^2 .* Eps_x_xz .* eps_0 .* Mu_yx_km .* mu_0 .* dx .* dz );

C_sxx1 = 1 ./ ( w.^2 .* Eps_x_xy .* eps_0 .* Mu_zy_jm .* mu_0 .* (dy.^2) );
C_sxx2 = 1 ./ ( w.^2 .* Eps_x_xy .* eps_0 .* Mu_zy_j  .* mu_0 .* (dy.^2) );
C_sxx3 = 1 ./ ( w.^2 .* Eps_x_xz .* eps_0 .* Mu_yz_k  .* mu_0 .* (dz.^2) );
C_sxx4 = 1 ./ ( w.^2 .* Eps_x_xz .* eps_0 .* Mu_yz_km .* mu_0 .* (dz.^2) );
C_sxx5 = 1 - C_sxx1 - C_sxx2 - C_sxx3 - C_sxx4;

% Contrast (for RHS in your formulation)
C_eix  = (eps_0 - (eps_r_xi .* eps_0)) ./ (eps_r_xi .* eps_0);

% Mixed “hi*” (RHS-related; kept for completeness)
C_hixz1 = (mu_0 - (Mu_zi_j  .* mu_0)) ./ ( 1i .* w .* Eps_x_xy .* eps_0 .* Mu_zi_j  .* mu_0 .* dy );
C_hixz2 = (mu_0 - (Mu_zi_jm .* mu_0)) ./ ( 1i .* w .* Eps_x_xy .* eps_0 .* Mu_zi_jm .* mu_0 .* dy );
C_hixy1 = (mu_0 - (Mu_yi_k  .* mu_0)) ./ ( 1i .* w .* Eps_x_xz .* eps_0 .* Mu_yi_k  .* mu_0 .* dz );
C_hixy2 = (mu_0 - (Mu_yi_km .* mu_0)) ./ ( 1i .* w .* Eps_x_xz .* eps_0 .* Mu_yi_km .* mu_0 .* dz );

% --------------------------- EY coefficients ------------------------------
C_syz1 = 1 ./ ( w.^2 .* Eps_y_yz .* eps_0 .* Mu_xy_k  .* mu_0 .* dy .* dz );
C_syz2 = 1 ./ ( w.^2 .* Eps_y_yz .* eps_0 .* Mu_xy_km .* mu_0 .* dy .* dz ); 

C_syx1 = 1 ./ ( w.^2 .* Eps_y_yx .* eps_0 .* Mu_zy_i  .* mu_0 .* dy .* dx );
C_syx2 = 1 ./ ( w.^2 .* Eps_y_yx .* eps_0 .* Mu_zy_im .* mu_0 .* dy .* dx );

C_syy1 = 1 ./ ( w.^2 .* Eps_y_yz .* eps_0 .* Mu_xz_km .* mu_0 .* (dz.^2) ); 
C_syy2 = 1 ./ ( w.^2 .* Eps_y_yz .* eps_0 .* Mu_xz_k  .* mu_0 .* (dz.^2) ); 
C_syy3 = 1 ./ ( w.^2 .* Eps_y_yx .* eps_0 .* Mu_zx_i  .* mu_0 .* (dx.^2) ); 
C_syy4 = 1 ./ ( w.^2 .* Eps_y_yx .* eps_0 .* Mu_zx_im .* mu_0 .* (dx.^2) );
C_syy5 = 1 - C_syy1 - C_syy2 - C_syy3 - C_syy4;

C_eiy  = (eps_0 - (eps_r_yi .* eps_0)) ./ (eps_r_yi .* eps_0);% done

C_hiyx1 = (mu_0 - (Mu_xi_k  .* mu_0)) ./ ( 1i .* w .* Eps_y_yz .* eps_0 .* Mu_xi_k  .* mu_0 .* dz );
C_hiyx2 = (mu_0 - (Mu_xi_km .* mu_0)) ./ ( 1i .* w .* Eps_y_yz .* eps_0 .* Mu_xi_km .* mu_0 .* dz );
C_hiyz1 = (mu_0 - (Mu_zi_i  .* mu_0)) ./ ( 1i .* w .* Eps_y_yx .* eps_0 .* Mu_zi_i  .* mu_0 .* dx );
C_hiyz2 = (mu_0 - (Mu_zi_im .* mu_0)) ./ ( 1i .* w .* Eps_y_yx .* eps_0 .* Mu_zi_im .* mu_0 .* dx );

% --------------------------- EZ coefficients ------------------------------
C_szx1 = 1 ./ ( w.^2 .* Eps_z_zx .* eps_0 .* Mu_yz_i  .* mu_0 .* dz .* dx );
C_szx2 = 1 ./ ( w.^2 .* Eps_z_zx .* eps_0 .* Mu_yz_im .* mu_0 .* dz .* dx );

C_szy1 = 1 ./ ( w.^2 .* Eps_z_zy .* eps_0 .* Mu_xz_j  .* mu_0 .* dz .* dy );
C_szy2 = 1 ./ ( w.^2 .* Eps_z_zy .* eps_0 .* Mu_xz_jm .* mu_0 .* dz .* dy );

C_szz1 = 1 ./ ( w.^2 .* Eps_z_zx .* eps_0 .* Mu_yx_im .* mu_0 .* (dx.^2) ); 
C_szz2 = 1 ./ ( w.^2 .* Eps_z_zx .* eps_0 .* Mu_yx_i  .* mu_0 .* (dx.^2) );
C_szz3 = 1 ./ ( w.^2 .* Eps_z_zy .* eps_0 .* Mu_xy_j  .* mu_0 .* (dy.^2) );
C_szz4 = 1 ./ ( w.^2 .* Eps_z_zy .* eps_0 .* Mu_xy_jm .* mu_0 .* (dy.^2) );
C_szz5 = 1 - C_szz1 - C_szz2 - C_szz3 - C_szz4;

C_eiz  = (eps_0 - (eps_r_zi .* eps_0)) ./ (eps_r_zi .* eps_0);

C_hizy1 = (mu_0 - (Mu_yi_k  .* mu_0)) ./ ( 1i .* w .* Eps_z_zx .* eps_0 .* Mu_yi_k  .* mu_0 .* dx );
C_hizy2 = (mu_0 - (Mu_yi_km .* mu_0)) ./ ( 1i .* w .* Eps_z_zx .* eps_0 .* Mu_yi_km .* mu_0 .* dx );
C_hizx1 = (mu_0 - (Mu_xi_j  .* mu_0)) ./ ( 1i .* w .* Eps_z_zy .* eps_0 .* Mu_xi_j  .* mu_0 .* dy );
C_hizx2 = (mu_0 - (Mu_xi_jm .* mu_0)) ./ ( 1i .* w .* Eps_z_zy .* eps_0 .* Mu_xi_jm .* mu_0 .* dy );

clear w;
