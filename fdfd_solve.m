% initialize the matlab workspace
clear all; close all; clc;

% define the problem 
define_problem_space_parameters; 
define_geometry; 
define_sources_and_lumped_elements; 

define_output_parameters;

% initialize the problem space and parameters 
initialize_fdfd_material_grid; 
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
    post_process_and_display_results;
    
end
