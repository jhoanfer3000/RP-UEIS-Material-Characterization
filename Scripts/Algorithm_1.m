% =========================================================================
% ALGORITHM 1: Sensitivity Analysis of PZT Disk Material Parameters
% =========================================================================
% Description: This script performs a local sensitivity analysis on 16 
% independent electromechanical parameters of a PZT-4 transducer. By 
% perturbing each parameter and tracking the variation in the simulated 
% impedance/admittance spectra, it quantifies their influence. This allows 
% the dimensionality of the subsequent inverse problem (RP-UEIS) to be 
% reduced to only the most influential parameters.
% Dependencies: COMSOL Multiphysics 6.1a with LiveLink for MATLAB, 
%               'evaluate_FEM' custom function.
% =========================================================================

clear all; clc;

%% 1. INITIALIZATION & SETUP
disp('--- Initializing Sensitivity Analysis ---');

% Load the axisymmetric FEM model of the isolated PZT disc
model_path = fullfile('..', 'models', 'piezoelectric_Disc');
model = mphload(model_path);

% Define the frequency window around the fundamental thickness-mode resonance
f0 = 1.03e6; % Nominal resonance frequency [Hz]
freq_set = linspace(f0*0.9, f0*1.1, 100); % +/- 10% bandwidth, 100 points

% Parameter names for PZT-4 (16 independent parameters total)
param_names = {'C11_real', 'C12_real', 'C13_real', 'C33_real', 'C44_real', ...
               'e15_real','e31_real', 'e33_real', 'eps11_real', 'eps33_real', ...
               'C11_imag', 'C12_imag', 'C13_imag', 'C33_imag', 'C44_imag', 'eps33_imag'};

% Known macroscopic properties of the PZT-4 sample
h = 2.02; % Thickness [mm]
rho_PZT = 7862; % Density [kg/m^3]

% Nominal baseline values (theta_0) derived from manufacturer datasheets.
% Note: Imaginary loss components are initialized as a fraction (e.g., 0.1%) 
% of their real counterparts to establish a baseline for sensitivity.
theta_0 = [1.38999e11, 7.78366e10, 7.42836e10, 1.15412e11, 2.5641e10, ... % Real Stiffness (Pa)
           12.718, -5.200, 15.080, 762.5, 663.2, ...                      % Piezo & Dielectric (Real)
           1.38999e8, 7.78366e7, -7.42836e7, 1.15412e8, 0.025, -0.663];   % Imaginary Losses 
       
% Logical mask to separate parameters: true for Real (Storage), false for Imaginary (Loss)
% This dictates which objective function metric (Impedance vs Admittance) is used later.
is_real = [true(1,10), false(1,6)]; 

Np = length(theta_0);  % Total number of material parameters (16)
Ns = 25;               % Number of discrete perturbation steps

%% 2. DEFINE SCALING FACTORS (alpha)
% Create an array of scaling factors to perturb the parameters over a range 
% of +/- 30% [0.7 to 1.3] relative to their nominal values. This wide range 
% captures the local convexity of the cost function.
alpha = zeros(1, Ns);
for j = 1:Ns
    alpha(j) = 0.7 + 0.6 * (j - 1) / (Ns - 1);
end

%% 3. COMPUTE REFERENCE SPECTRA
% Calculate the baseline response of the model using the nominal parameter set (theta_0).
% This serves as the "pseudo-experimental" target to evaluate the cost of perturbations.
fprintf('Computing reference FEM spectra...\n');
[Z_ref, Y_ref] = evaluate_FEM(theta_0, freq_set, h, rho_PZT, model);

mag_Z_ref = abs(Z_ref);        % Baseline Impedance Magnitude
mag_Y_ref = abs(Y_ref);        % Baseline Admittance Magnitude

%% 4. SENSITIVITY ANALYSIS LOOP
s_index = zeros(Np, 1); % Preallocate array to store the maximum sensitivity slope for each parameter

fprintf('Starting parameter perturbation loop...\n');

% --- Start Timer and Initialize Progress Bar ---
total_steps = Np * Ns; % Total number of FEM evaluations
step_count = 0;        % Step counter for the waitbar
h_wait = waitbar(0, 'Initializing sensitivity analysis...');
tic; % Start global execution timer
% -----------------------------------------------

for i = 1:Np
    fprintf('Evaluating parameter %d/%d: %s\n', i, Np, param_names{i});
    c_ij = zeros(1, Ns); % Cost array for the current parameter's perturbations
    
    for j = 1:Ns
        % --- Update Progress Bar ---
        step_count = step_count + 1;
        waitbar(step_count / total_steps, h_wait, ...
            sprintf('Evaluating %s: Step %d of %d', param_names{i}, step_count, total_steps));
        
        % Perturb the specific i-th parameter by the j-th scaling factor
        theta_pert = theta_0;
        theta_pert(i) = alpha(j) * theta_0(i);
        
        % Compute the FEM spectral response for the perturbed parameter set
        [Z_num, Y_num] = evaluate_FEM(theta_pert, freq_set, h, rho_PZT, model);
        
        % Calculate Cost based on the physical nature of the parameter
        if is_real(i)
            % REAL Parameters (Storage): Use the Logarithmic Impedance metric.
            % Highly sensitive to resonance frequency shifts.
            res_mag = log10(mag_Z_ref) - log10(abs(Z_num));
        else
            % IMAGINARY Parameters (Losses): Use the Linear Admittance metric.
            % Highly sensitive to peak broadening and damping at series resonance.
            res_mag = mag_Y_ref - abs(Y_num);
        end
        
        % Total Cost: Squared Euclidean norm of the residual vector
        R_vector = res_mag(:); 
        c_ij(j) = sum(R_vector.^2);
    end
    
    % Compute the local rate of change (derivative approximation) of the cost function
    slopes = zeros(1, Ns-1);
    for j = 1:(Ns-1)
        slopes(j) = abs((c_ij(j+1) - c_ij(j)) / (alpha(j+1) - alpha(j)));
    end
    
    % The Sensitivity Index (s_i) is defined as the maximum local rate of change
    s_index(i) = max(slopes);
end

% --- Stop Timer and Close Progress Bar ---
close(h_wait);
elapsed_time_sec = toc; 
elapsed_time_min = elapsed_time_sec / 60; 

fprintf('\nTotal execution time: %.2f minutes\n', elapsed_time_min);

%% 5. BUILD RANKING TABLE (Separated and Normalized)
% Because real and imaginary parameters utilize different residual metrics 
% (Impedance vs Admittance), their raw cost values are on different scales. 
% Therefore, they must be separated and normalized independently.

% Logical indexing to filter Real and Imaginary parameters
idx_real = (is_real == true)';
idx_imag = (is_real == false)';

% --- PROCESS REAL PARAMETERS (STORAGE) ---
names_real = param_names(idx_real)';
val_real   = theta_0(idx_real)';
s_real     = s_index(idx_real);

% Normalize against the maximum slope within the Real group (e.g., C33_real)
s_real_norm = s_real / max(s_real);

% Create and sort the Real sensitivity table
T_real = table(names_real, val_real, s_real_norm, ...
    'VariableNames', {'Parameter', 'Nominal_Value', 'Norm_Sensitivity'});
T_real_sorted = sortrows(T_real, 'Norm_Sensitivity', 'descend');

% --- PROCESS IMAGINARY PARAMETERS (LOSSES) ---
names_imag = param_names(idx_imag)';
val_imag   = theta_0(idx_imag)';
s_imag     = s_index(idx_imag);

% Normalize against the maximum slope within the Imaginary group (e.g., C33_imag)
s_imag_norm = s_imag / max(s_imag);

% Create and sort the Imaginary sensitivity table
T_imag = table(names_imag, val_imag, s_imag_norm, ...
    'VariableNames', {'Parameter', 'Nominal_Value', 'Norm_Sensitivity'});
T_imag_sorted = sortrows(T_imag, 'Norm_Sensitivity', 'descend');

% --- DISPLAY FINAL RESULTS ---
fprintf('\n======================================================\n');
disp('--- REAL PARAMETERS SENSITIVITY RANKING (Normalized) ---');
disp(T_real_sorted);

fprintf('\n======================================================\n');
disp('--- IMAGINARY PARAMETERS SENSITIVITY RANKING (Normalized) ---');
disp(T_imag_sorted);
fprintf('======================================================\n');

%% 6. EXPORT RESULTS

disp('Saving sensitivity ranking results...');

% 1. Define the relative path to the 'data' folder using '..'
data_folder = fullfile('..', 'data');

% 2. Create the 'results' folder if it doesn't exist yet
if ~exist(data_folder, 'dir')
    mkdir(data_folder);
end

% 3. Define the relative paths for the output files
mat_file = fullfile(data_folder, 'sensitivity_ranking.mat');
csv_real = fullfile(data_folder, 'Sensitivity_Real_Ranking.csv');
csv_imag = fullfile(data_folder, 'Sensitivity_Imag_Ranking.csv');

% 4. Save MATLAB workspace data
save(mat_file, 'T_real_sorted', 'T_imag_sorted', 's_index');

% 5. Export to CSV for universal access (great for GitHub visualization)
writetable(T_real_sorted, csv_real);
writetable(T_imag_sorted, csv_imag);

disp('Results saved successfully as .mat and .csv files in the ../results folder.');
disp('Done!');