% initialize boundary parameters

% define logical parameters for the conditions that  will be used often for
% cpml
is_cpml_xn = false; is_cpml_xp = false; is_cpml_yn = false; 
is_cpml_yp = false; is_cpml_zn = false; is_cpml_zp = false; 
is_any_side_cpml = false;
if strcmp(boundary.type_xn, 'cpml')
    is_cpml_xn = true;
    n_cpml_xn = abs(boundary.cpml_number_of_cells_xn);
end
if strcmp(boundary.type_xp, 'cpml')
    is_cpml_xp = true;
    n_cpml_xp = abs(boundary.cpml_number_of_cells_xp);
end
if strcmp(boundary.type_yn, 'cpml')
    is_cpml_yn = true;
    n_cpml_yn = abs(boundary.cpml_number_of_cells_yn);
end
if strcmp(boundary.type_yp, 'cpml')
    is_cpml_yp = true;
    n_cpml_yp = abs(boundary.cpml_number_of_cells_yp);
end
if strcmp(boundary.type_zn, 'cpml')
    is_cpml_zn = true;
    n_cpml_zn = abs(boundary.cpml_number_of_cells_zn);
end
if strcmp(boundary.type_zp, 'cpml')
    is_cpml_zp = true;
    n_cpml_zp = abs(boundary.cpml_number_of_cells_zp);
end

if (is_cpml_xn || is_cpml_xp || is_cpml_yn ...
        || is_cpml_yp || is_cpml_zn || is_cpml_zp)
    is_any_side_cpml = true;
end

% Call CPML initialization routine if any side is CPML
if is_any_side_cpml
    initialize_CPML_ABC;
end



% define logical parameters for the conditions that  will be used often for
% pml
is_pml_xn = false; is_pml_xp = false; is_pml_yn = false; 
is_pml_yp = false; is_pml_zn = false; is_pml_zp = false; 
is_any_side_pml = false;
if strcmp(boundary.type_xn, 'pml')
    is_pml_xn = true;
    n_pml_xn = abs(boundary.pml_number_of_cells_xn);
end
if strcmp(boundary.type_xp, 'pml')
    is_pml_xp = true;
    n_pml_xp = abs(boundary.pml_number_of_cells_xp);
end
if strcmp(boundary.type_yn, 'pml')
    is_pml_yn = true;
    n_pml_yn = abs(boundary.pml_number_of_cells_yn);
end
if strcmp(boundary.type_yp, 'pml')
    is_pml_yp = true;
    n_pml_yp = abs(boundary.pml_number_of_cells_yp);
end
if strcmp(boundary.type_zn, 'pml')
    is_pml_zn = true;
    n_pml_zn = abs(boundary.pml_number_of_cells_zn);
end
if strcmp(boundary.type_zp, 'pml')
    is_pml_zp = true;
    n_pml_zp = abs(boundary.pml_number_of_cells_zp);
end

if (is_pml_xn || is_pml_xp || is_pml_yn ...
        || is_pml_yp || is_pml_zn || is_pml_zp)
    is_any_side_pml = true;
end

% Call PML initialization routine if any side is PML
if is_any_side_pml
    initialize_PML_ABC;
end