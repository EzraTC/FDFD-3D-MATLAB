function [A, y, x] = assemble_fdfd_Ay_solve( ...
    nx, ny, nz, nxp1, nyp1, nzp1, ...
    C_sxy1, C_sxy2, C_sxz1, C_sxz2, C_sxx5, ...
    C_syz1, C_syz2, C_syx1, C_syx2, C_syy5, ...
    C_szx1, C_szx2, C_szy1, C_szy2, C_szz5, ...
    E_inc_x, E_inc_y, E_inc_z)

% ---------- sizes ----------
N  = nx * ny * nz;           % unknowns per component
NT = 3 * N;                  % total unknowns

% ---------- global ids for each component block (i-fastest) ----------
ex_id = reshape(           1:N, nx, ny, nz);
ey_id = reshape(N +        1:N, nx, ny, nz);
ez_id = reshape(2*N +      1:N, nx, ny, nz);

% ---------- neighbor id arrays (replicate at boundary to keep full size) ----------
% j-1 shift (for coupling to j-1)
ey_id_jm = cat(2, ey_id(:,1,:), ey_id(:,1:end-1,:));    % size [nx,ny,nz]
% k-1 shift
ez_id_km = cat(3, ez_id(:,:,1), ez_id(:,:,1:end-1));    % size [nx,ny,nz]
% i-1 shift
ex_id_im = cat(1, ex_id(1,:,:), ex_id(1:end-1,:,:));    % size [nx,ny,nz]

% For Ey row (needs ex(i-1,j,k) and ez(i,j,k-1))
ez_id_km_forEy = ez_id_km;                 % alias for clarity
ex_id_im_forEy = ex_id_im;

% For Ez row (needs ex(i-1,j,k) and ey(i,j-1,k))
ex_id_im_forEz = ex_id_im;
ey_id_jm_forEz = ey_id_jm;

% ---------- assemble sparse triplets ----------
% Estimate ~5 nnz per row (diag + 4 cross terms) => ~15N total
est_nnz = 15 * N;
rows = zeros(est_nnz,1);
cols = zeros(est_nnz,1);
vals = zeros(est_nnz,1);
ptr  = 1;

% --- helper to append a block of triplets ---
    function append_triplets(R, C, V)
        m = numel(V);
        rows(ptr:ptr+m-1) = R(:);
        cols(ptr:ptr+m-1) = C(:);
        vals(ptr:ptr+m-1) = V(:);
        ptr = ptr + m;
    end

% ===================== EX rows =====================
% diagonal
append_triplets(ex_id, ex_id, C_sxx5);

% cross to Ey (j and j-1)
append_triplets(ex_id, ey_id,     -C_sxy1);
append_triplets(ex_id, ey_id_jm,  -C_sxy2);

% cross to Ez (k and k-1)
append_triplets(ex_id, ez_id,     -C_sxz1);
append_triplets(ex_id, ez_id_km,  -C_sxz2);

% ===================== EY rows =====================
% diagonal
append_triplets(ey_id, ey_id, C_syy5);

% cross to Ex (i and i-1)
append_triplets(ey_id, ex_id,          -C_syx1);
append_triplets(ey_id, ex_id_im_forEy, -C_syx2);

% cross to Ez (k and k-1)
append_triplets(ey_id, ez_id,          -C_syz1);
append_triplets(ey_id, ez_id_km_forEy, -C_syz2);

% ===================== EZ rows =====================
% diagonal
append_triplets(ez_id, ez_id, C_szz5);

% cross to Ex (i and i-1)
append_triplets(ez_id, ex_id,          -C_szx1);
append_triplets(ez_id, ex_id_im_forEz, -C_szx2);

% cross to Ey (j and j-1)
append_triplets(ez_id, ey_id,          -C_szy1);
append_triplets(ez_id, ey_id_jm_forEz, -C_szy2);

% Trim to actual length and build sparse
rows = rows(1:ptr-1); cols = cols(1:ptr-1); vals = vals(1:ptr-1);
A = sparse(rows, cols, vals, NT, NT);

% ---------- build RHS y from incident fields ----------
% Take interior windows aligned with unknowns:
% Ex: (i=1..nx, j=2..nyp1 -> 1..ny, k=2..nzp1 -> 1..nz)
Ex_inc_int = E_inc_x(1:nx, 2:nyp1, 2:nzp1);
% Ey: (i=2..nxp1 -> 1..nx, j=1..ny, k=2..nzp1 -> 1..nz)
Ey_inc_int = E_inc_y(2:nxp1, 1:ny, 2:nzp1);
% Ez: (i=2..nxp1 -> 1..nx, j=2..nyp1 -> 1..ny, k=1..nz)
Ez_inc_int = E_inc_z(2:nxp1, 2:nyp1, 1:nz);

y = [Ex_inc_int(:); Ey_inc_int(:); Ez_inc_int(:)];

% ---------- solve ----------
x = A \ y;    % scattered-field unknowns stacked as [Ex; Ey; Ez]

end
