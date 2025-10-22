disp('initializing the output parameters');
% initialize_output_parameters.m  (FDFD: single-frequency)

number_of_sampled_electric_fields  = size(sampled_electric_fields, 2);
number_of_sampled_magnetic_fields  = size(sampled_magnetic_fields, 2);
number_of_sampled_voltages         = size(sampled_voltages, 2);
number_of_sampled_currents         = size(sampled_currents, 2);
number_of_ports                    = size(ports, 2);

% ---- Frequency-domain 
frequency_domain.frequencies = frequency_domain.f;   
frequency_domain.number_of_frequencies = 1;
if ~isfield(frequency_domain, 'w') || isempty(frequency_domain.w)
    frequency_domain.w = 2*pi*frequency_domain.f;    % rad/s
end

% ========================= Electric-field samples =========================
% Map each requested sample to the correct Yee grid indices based on component.
% Ex grid: [nx,   nyp1, nzp1] at (x+1/2, y,   z)
% Ey grid: [nxp1, ny,   nzp1] at (x,     y+1/2,z)
% Ez grid: [nxp1, nyp1, nz  ] at (x,     y,    z+1/2)
for ind = 1:number_of_sampled_electric_fields
    s = sampled_electric_fields(ind);

    % Default to 'z' if component not provided
    if ~isfield(s,'component') || isempty(s.component)
        s.component = 'z';
    end

    switch lower(s.component)
        case 'x'  % Ex grid
            is = round((s.x - fdfd_domain.min_x)/dx) + 1;   % 1..nx
            js = round((s.y - fdfd_domain.min_y)/dy) + 1;   % 1..nyp1
            ks = round((s.z - fdfd_domain.min_z)/dz) + 1;   % 1..nzp1
            % clamp
            is = min(max(is,1), nx);
            js = min(max(js,1), nyp1);
            ks = min(max(ks,1), nzp1);

        case 'y'  % Ey grid
            is = round((s.x - fdfd_domain.min_x)/dx) + 1;   % 1..nxp1
            js = round((s.y - fdfd_domain.min_y)/dy) + 1;   % 1..ny
            ks = round((s.z - fdfd_domain.min_z)/dz) + 1;   % 1..nzp1
            is = min(max(is,1), nxp1);
            js = min(max(js,1), ny);
            ks = min(max(ks,1), nzp1);

        case 'z'  % Ez grid
            is = round((s.x - fdfd_domain.min_x)/dx) + 1;   % 1..nxp1
            js = round((s.y - fdfd_domain.min_y)/dy) + 1;   % 1..nyp1
            ks = round((s.z - fdfd_domain.min_z)/dz) + 1;   % 1..nz
            is = min(max(is,1), nxp1);
            js = min(max(js,1), nyp1);
            ks = min(max(ks,1), nz);

        otherwise
            error('initialize_output_parameters: unknown E-field component "%s". Use ''x'',''y'',''z''.', s.component);
    end

    sampled_electric_fields(ind).component = s.component;  % persist
    sampled_electric_fields(ind).is = is;
    sampled_electric_fields(ind).js = js;
    sampled_electric_fields(ind).ks = ks;

    % Single-frequency phasor placeholder (complex scalar)
    sampled_electric_fields(ind).sampled_value = complex(0);
    % Keep a record of the evaluation frequency
    sampled_electric_fields(ind).frequency = frequency_domain.f;
    % No time axis in FDFD
    sampled_electric_fields(ind).time = [];
end

% ========================= Magnetic-field samples =========================
for ind = 1:number_of_sampled_magnetic_fields
    is = round((sampled_magnetic_fields(ind).x - fdfd_domain.min_x)/dx) + 1;
    js = round((sampled_magnetic_fields(ind).y - fdfd_domain.min_y)/dy) + 1;
    ks = round((sampled_magnetic_fields(ind).z - fdfd_domain.min_z)/dz) + 1;

    sampled_magnetic_fields(ind).is = is;
    sampled_magnetic_fields(ind).js = js;
    sampled_magnetic_fields(ind).ks = ks;

    sampled_magnetic_fields(ind).sampled_value = complex(0);
    sampled_magnetic_fields(ind).frequency = frequency_domain.f;
    sampled_magnetic_fields(ind).time = [];
end

% =============================== Voltages ================================
for ind = 1:number_of_sampled_voltages
    is = round((sampled_voltages(ind).min_x - fdfd_domain.min_x)/dx) + 1;
    js = round((sampled_voltages(ind).min_y - fdfd_domain.min_y)/dy) + 1;
    ks = round((sampled_voltages(ind).min_z - fdfd_domain.min_z)/dz) + 1;

    ie = round((sampled_voltages(ind).max_x - fdfd_domain.min_x)/dx) + 1;
    je = round((sampled_voltages(ind).max_y - fdfd_domain.min_y)/dy) + 1;
    ke = round((sampled_voltages(ind).max_z - fdfd_domain.min_z)/dz) + 1;

    sampled_voltages(ind).is = is;  sampled_voltages(ind).ie = ie;
    sampled_voltages(ind).js = js;  sampled_voltages(ind).je = je;
    sampled_voltages(ind).ks = ks;  sampled_voltages(ind).ke = ke;

    % Field index list on the appropriate Yee faces
    switch (sampled_voltages(ind).direction(1))
        case 'x'
            fi = create_linear_index_list(Ex, is:ie-1, js:je,   ks:ke);
            sampled_voltages(ind).Csvf = -dx/((je-js+1)*(ke-ks+1));
        case 'y'
            fi = create_linear_index_list(Ey, is:ie,   js:je-1, ks:ke);
            sampled_voltages(ind).Csvf = -dy/((ke-ks+1)*(ie-is+1));
        case 'z'
            fi = create_linear_index_list(Ez, is:ie,   js:je,   ks:ke-1);
            sampled_voltages(ind).Csvf = -dz/((ie-is+1)*(je-js+1));
    end
    if strcmp(sampled_voltages(ind).direction(2), 'n')
        sampled_voltages(ind).Csvf = -sampled_voltages(ind).Csvf;
    end
    sampled_voltages(ind).field_indices = fi;

    % Single-frequency phasor placeholder
    sampled_voltages(ind).sampled_value = complex(0);
    sampled_voltages(ind).frequency = frequency_domain.f;
    sampled_voltages(ind).time = [];
end

% =============================== Currents ================================
for ind = 1:number_of_sampled_currents
    is = round((sampled_currents(ind).min_x - fdfd_domain.min_x)/dx) + 1;
    js = round((sampled_currents(ind).min_y - fdfd_domain.min_y)/dy) + 1;
    ks = round((sampled_currents(ind).min_z - fdfd_domain.min_z)/dz) + 1;

    ie = round((sampled_currents(ind).max_x - fdfd_domain.min_x)/dx) + 1;
    je = round((sampled_currents(ind).max_y - fdfd_domain.min_y)/dy) + 1;
    ke = round((sampled_currents(ind).max_z - fdfd_domain.min_z)/dz) + 1;

    sampled_currents(ind).is = is;  sampled_currents(ind).ie = ie;
    sampled_currents(ind).js = js;  sampled_currents(ind).je = je;
    sampled_currents(ind).ks = ks;  sampled_currents(ind).ke = ke;

    sampled_currents(ind).sampled_value = complex(0);
    sampled_currents(ind).frequency = frequency_domain.f;
    sampled_currents(ind).time = [];
end
