% initialize the matlab workspace
clear all; close all; clc;

% define the problem 
define_problem_space_parameters; %Need - Done *mentions  cpml* *fixed for pml*
define_geometry; %Need - Done
define_sources_and_lumped_elements; %Need - Done

define_output_parameters;%Need - Done but needs to have some additions before useage

% initialize the problem space and parameters 
initialize_fdfd_material_grid; %Need - Done
display_problem_space; 
display_material_mesh;

if run_simulation
    initialize_fdfd_parameters_and_arrays;
    initialize_sources_and_lumped_elements;
   
    initialize_boundary_conditions;
    initialize_updating_coefficients;

    initialize_output_parameters;
    initialize_display_parameters;

    solve_matrix_equation;

    % display simulation results
    % post_process_and_display_results;
    
end
