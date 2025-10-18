 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a/solve_matrix_equation.m b/solve_matrix_equation.m
index a0832cc0cda3880c61f331cbe3b52246f7b8e8f9..806f53bf54802502458c0c1b68c07cd87afb82c0 100644
--- a/solve_matrix_equation.m
+++ b/solve_matrix_equation.m
@@ -1,104 +1,155 @@
-function [A, y, x] = assemble_fdfd_Ay_solve( ...
-    nx, ny, nz, nxp1, nyp1, nzp1, ...
-    C_sxy1, C_sxy2, C_sxz1, C_sxz2, C_sxx5, ...
-    C_syz1, C_syz2, C_syx1, C_syx2, C_syy5, ...
-    C_szx1, C_szx2, C_szy1, C_szy2, C_szz5, ...
-    E_inc_x, E_inc_y, E_inc_z)
-
-% ---------- sizes ----------
-N  = nx * ny * nz;           % unknowns per component
-NT = 3 * N;                  % total unknowns
-
-% ---------- global ids for each component block (i-fastest) ----------
-ex_id = reshape(           1:N, nx, ny, nz);
-ey_id = reshape(N +        1:N, nx, ny, nz);
-ez_id = reshape(2*N +      1:N, nx, ny, nz);
-
-% ---------- neighbor id arrays (replicate at boundary to keep full size) ----------
-% j-1 shift (for coupling to j-1)
-ey_id_jm = cat(2, ey_id(:,1,:), ey_id(:,1:end-1,:));    % size [nx,ny,nz]
-% k-1 shift
-ez_id_km = cat(3, ez_id(:,:,1), ez_id(:,:,1:end-1));    % size [nx,ny,nz]
-% i-1 shift
-ex_id_im = cat(1, ex_id(1,:,:), ex_id(1:end-1,:,:));    % size [nx,ny,nz]
-
-% For Ey row (needs ex(i-1,j,k) and ez(i,j,k-1))
-ez_id_km_forEy = ez_id_km;                 % alias for clarity
-ex_id_im_forEy = ex_id_im;
-
-% For Ez row (needs ex(i-1,j,k) and ey(i,j-1,k))
-ex_id_im_forEz = ex_id_im;
-ey_id_jm_forEz = ey_id_jm;
-
-% ---------- assemble sparse triplets ----------
-% Estimate ~5 nnz per row (diag + 4 cross terms) => ~15N total
-est_nnz = 15 * N;
-rows = zeros(est_nnz,1);
-cols = zeros(est_nnz,1);
-vals = zeros(est_nnz,1);
-ptr  = 1;
-
-% --- helper to append a block of triplets ---
-    function append_triplets(R, C, V)
-        m = numel(V);
-        rows(ptr:ptr+m-1) = R(:);
-        cols(ptr:ptr+m-1) = C(:);
-        vals(ptr:ptr+m-1) = V(:);
-        ptr = ptr + m;
-    end
-
-% ===================== EX rows =====================
-% diagonal
-append_triplets(ex_id, ex_id, C_sxx5);
-
-% cross to Ey (j and j-1)
-append_triplets(ex_id, ey_id,     -C_sxy1);
-append_triplets(ex_id, ey_id_jm,  -C_sxy2);
-
-% cross to Ez (k and k-1)
-append_triplets(ex_id, ez_id,     -C_sxz1);
-append_triplets(ex_id, ez_id_km,  -C_sxz2);
-
-% ===================== EY rows =====================
-% diagonal
-append_triplets(ey_id, ey_id, C_syy5);
-
-% cross to Ex (i and i-1)
-append_triplets(ey_id, ex_id,          -C_syx1);
-append_triplets(ey_id, ex_id_im_forEy, -C_syx2);
-
-% cross to Ez (k and k-1)
-append_triplets(ey_id, ez_id,          -C_syz1);
-append_triplets(ey_id, ez_id_km_forEy, -C_syz2);
-
-% ===================== EZ rows =====================
-% diagonal
-append_triplets(ez_id, ez_id, C_szz5);
-
-% cross to Ex (i and i-1)
-append_triplets(ez_id, ex_id,          -C_szx1);
-append_triplets(ez_id, ex_id_im_forEz, -C_szx2);
-
-% cross to Ey (j and j-1)
-append_triplets(ez_id, ey_id,          -C_szy1);
-append_triplets(ez_id, ey_id_jm_forEz, -C_szy2);
-
-% Trim to actual length and build sparse
-rows = rows(1:ptr-1); cols = cols(1:ptr-1); vals = vals(1:ptr-1);
-A = sparse(rows, cols, vals, NT, NT);
-
-% ---------- build RHS y from incident fields ----------
-% Take interior windows aligned with unknowns:
-% Ex: (i=1..nx, j=2..nyp1 -> 1..ny, k=2..nzp1 -> 1..nz)
-Ex_inc_int = E_inc_x(1:nx, 2:nyp1, 2:nzp1);
-% Ey: (i=2..nxp1 -> 1..nx, j=1..ny, k=2..nzp1 -> 1..nz)
-Ey_inc_int = E_inc_y(2:nxp1, 1:ny, 2:nzp1);
-% Ez: (i=2..nxp1 -> 1..nx, j=2..nyp1 -> 1..ny, k=1..nz)
-Ez_inc_int = E_inc_z(2:nxp1, 2:nyp1, 1:nz);
-
-y = [Ex_inc_int(:); Ey_inc_int(:); Ez_inc_int(:)];
-
-% ---------- solve ----------
-x = A \ y;    % scattered-field unknowns stacked as [Ex; Ey; Ez]
-
-end
+%% solve_matrix_equation.m
+% Assemble the frequency–domain FDFD matrix system and solve A * x = y, where
+% the unknown vector x contains the scattered electric-field components on the
+% interior Yee cells and the right-hand side y holds the incident fields.  The
+% script expects all geometric, material, and coefficient arrays to be prepared
+% by the initialization routines that precede it in ``fdfd_solve``.
+
+disp('assembling linear system and solving for scattered electric fields');
+
+%% ------------------------------ Dimensions ---------------------------------
+N  = nx * ny * nz;      % Number of unknowns per Cartesian component
+NT = 3 * N;             % Total size of the linear system (Ex,Ey,Ez blocks)
+
+% Global ids for each block (flattened with i running fastest)
+ex_id = reshape(          1:N, nx, ny, nz);
+ey_id = reshape(    N + (1:N), nx, ny, nz);
+ez_id = reshape(2 * N + (1:N), nx, ny, nz);
+
+% Neighbor ids (replicated at the boundaries so every array stays full size)
+ey_id_jm = cat(2, ey_id(:,1,:),     ey_id(:,1:end-1,:));
+ez_id_km = cat(3, ez_id(:,:,1),     ez_id(:,:,1:end-1));
+ex_id_im = cat(1, ex_id(1,:,:),     ex_id(1:end-1,:,:));
+
+ez_id_km_forEy = ez_id_km;   % aliases for clarity when building the matrix
+ex_id_im_forEy = ex_id_im;
+ex_id_im_forEz = ex_id_im;
+ey_id_jm_forEz = ey_id_jm;
+
+%% ---------------------------- Sparse Assembly -------------------------------
+% Diagonal + four couplings per row -> 5 contributions per component block.
+rows = [
+    ex_id(:);
+    ex_id(:);
+    ex_id(:);
+    ex_id(:);
+    ex_id(:);
+    ey_id(:);
+    ey_id(:);
+    ey_id(:);
+    ey_id(:);
+    ey_id(:);
+    ez_id(:);
+    ez_id(:);
+    ez_id(:);
+    ez_id(:);
+    ez_id(:)
+];
+cols = [
+    ex_id(:);
+    ey_id(:);
+    ey_id_jm(:);
+    ez_id(:);
+    ez_id_km(:);
+    ey_id(:);
+    ex_id(:);
+    ex_id_im_forEy(:);
+    ez_id(:);
+    ez_id_km_forEy(:);
+    ez_id(:);
+    ex_id(:);
+    ex_id_im_forEz(:);
+    ey_id(:);
+    ey_id_jm_forEz(:)
+];
+vals = [
+    C_sxx5(:);
+    -C_sxy1(:);
+    -C_sxy2(:);
+    -C_sxz1(:);
+    -C_sxz2(:);
+    C_syy5(:);
+    -C_syx1(:);
+    -C_syx2(:);
+    -C_syz1(:);
+    -C_syz2(:);
+    C_szz5(:);
+    -C_szx1(:);
+    -C_szx2(:);
+    -C_szy1(:);
+    -C_szy2(:)
+];
+
+A = sparse(rows, cols, vals, NT, NT);
+
+%% --------------------- Incident-Field Right-Hand Side -----------------------
+% Build the incident fields on the Yee electric grids if they are not already
+% present (plane-wave excitation defined in initialize_waveforms).
+if ~exist('waveforms', 'var') || ~isfield(waveforms, 'plane_wave')
+    error('Plane-wave definition (waveforms.plane_wave) is required before solving.');
+end
+
+theta_inc = waveforms.plane_wave.theta_inc;
+phi_inc   = waveforms.plane_wave.phi_inc;
+E_theta   = waveforms.plane_wave.E_theta;
+E_phi     = waveforms.plane_wave.E_phi;
+k0        = waveforms.plane_wave.k0;
+khat      = waveforms.plane_wave.khat(:).';
+
+Ex_coef =  E_theta*cos(theta_inc)*cos(phi_inc) - E_phi*sin(phi_inc);
+Ey_coef =  E_theta*cos(theta_inc)*sin(phi_inc) + E_phi*cos(phi_inc);
+Ez_coef = -E_theta*sin(theta_inc);
+
+% Coordinate vectors for the three electric-field staggered grids
+x_ex = fdfd_domain.min_x + ((0:nx-1) + 0.5) * dx;
+y_ex = fdfd_domain.min_y + (0:nyp1-1)       * dy;
+z_ex = fdfd_domain.min_z + (0:nzp1-1)       * dz;
+[xEx, yEx, zEx] = ndgrid(x_ex, y_ex, z_ex);
+E_inc_x = Ex_coef * exp(-1i * k0 * (khat(1) * xEx + khat(2) * yEx + khat(3) * zEx));
+
+x_ey = fdfd_domain.min_x + (0:nxp1-1)       * dx;
+y_ey = fdfd_domain.min_y + ((0:ny-1) + 0.5) * dy;
+z_ey = fdfd_domain.min_z + (0:nzp1-1)       * dz;
+[xEy, yEy, zEy] = ndgrid(x_ey, y_ey, z_ey);
+E_inc_y = Ey_coef * exp(-1i * k0 * (khat(1) * xEy + khat(2) * yEy + khat(3) * zEy));
+
+x_ez = fdfd_domain.min_x + (0:nxp1-1)       * dx;
+y_ez = fdfd_domain.min_y + (0:nyp1-1)       * dy;
+z_ez = fdfd_domain.min_z + ((0:nz-1) + 0.5) * dz;
+[xEz, yEz, zEz] = ndgrid(x_ez, y_ez, z_ez);
+E_inc_z = Ez_coef * exp(-1i * k0 * (khat(1) * xEz + khat(2) * yEz + khat(3) * zEz));
+
+% Extract the interior slices that correspond to unknowns
+Ex_inc_int = E_inc_x(:, 2:nyp1, 2:nzp1);
+Ey_inc_int = E_inc_y(2:nxp1, :, 2:nzp1);
+Ez_inc_int = E_inc_z(2:nxp1, 2:nyp1, :);
+
+y = [Ex_inc_int(:); Ey_inc_int(:); Ez_inc_int(:)];
+
+%% ------------------------------ Solve System --------------------------------
+x = A \ y;
+
+% Reshape the scattered-field solution back onto Yee-cell interiors
+Ex_scat = reshape(x(1:N),       nx, ny, nz);
+Ey_scat = reshape(x(N+1:2*N),   nx, ny, nz);
+Ez_scat = reshape(x(2*N+1:end), nx, ny, nz);
+
+% Store scattered, incident, and total fields for later post-processing
+E_scattered.x = Ex_scat;
+E_scattered.y = Ey_scat;
+E_scattered.z = Ez_scat;
+
+E_incident.x = Ex_inc_int;
+E_incident.y = Ey_inc_int;
+E_incident.z = Ez_inc_int;
+
+E_total.x = Ex_inc_int + Ex_scat;
+E_total.y = Ey_inc_int + Ey_scat;
+E_total.z = Ez_inc_int + Ez_scat;
+
+% Insert total fields back into the Yee arrays (interior indices only)
+Ex(:, 2:nyp1, 2:nzp1) = E_total.x;
+Ey(2:nxp1, :, 2:nzp1) = E_total.y;
+Ez(2:nxp1, 2:nyp1, :) = E_total.z;
+
+clear xEx yEx zEx xEy yEy zEy xEz yEz zEz; 
EOF
)

