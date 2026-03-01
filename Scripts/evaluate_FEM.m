function [Z_ref, Y_ref] = evaluate_FEM(theta, freq_set, h, rho_PZT, model)
    % =====================================================================
    % EVALUATE_FEM: COMSOL LiveLink Interface for PZT Simulation
    % =====================================================================
    % Description: Updates the finite element model with a given set of 
    % 16 electromechanical parameters (theta), executes the frequency 
    % sweep study, and extracts the complex electrical response.
    %
    % Inputs:
    %   - theta: 1x16 vector of PZT parameters (Real and Imaginary parts)
    %   - freq_set: Array of frequency points for the sweep [Hz]
    %   - h: Thickness of the PZT disk [mm]
    %   - rho_PZT: Density of the PZT material [kg/m^3]
    %   - model: Active COMSOL model object
    %
    % Outputs:
    %   - Z_ref: Complex electrical impedance [Ohms] (Column vector)
    %   - Y_ref: Complex electrical admittance [S] (Column vector)
    % =====================================================================

    %% 1. Extract Frequency Sweep Parameters
    % COMSOL requires the step size, start, and stop frequencies for the study
    freq_step = freq_set(2) - freq_set(1);
    freq_ini  = freq_set(1);
    freq_end  = freq_set(end);
    
    % Send frequency sweep boundaries to the COMSOL parameter node
    model.param.set('Fre_I', sprintf('%g [Hz]', freq_ini)); 
    model.param.set('Step_Fre', sprintf('%g [Hz]', freq_step)); 
    model.param.set('Fre_F', sprintf('%g [Hz]', freq_end + freq_step)); 

    %% 2. Map 'theta' to COMSOL Material Parameters
    % Note: sprintf('%g + %g*i') is used to safely format complex numbers, 
    % preventing COMSOL parsing errors during automated optimization loops.
    
    % --- Mechanical Stiffness (Elastic) Matrix [c^E] ---
    % theta(1:5) are the real (storage) parts; theta(11:15) are the imaginary (loss) parts
    model.param.set('CE11', sprintf('%g + %g*i [Pa]', theta(1), theta(11)));
    model.param.set('CE12', sprintf('%g + %g*i [Pa]', theta(2), theta(12)));
    model.param.set('CE13', sprintf('%g + %g*i [Pa]', theta(3), theta(13)));
    model.param.set('CE33', sprintf('%g + %g*i [Pa]', theta(4), theta(14)));
    model.param.set('CE44', sprintf('%g + %g*i [Pa]', theta(5), theta(15)));
    
    % --- Piezoelectric Coupling Matrix [e] ---
    % Assumed lossless (real values only) for this model formulation
    model.param.set('eE15', sprintf('%g [C/m^2]', theta(6)));
    model.param.set('eE31', sprintf('%g [C/m^2]', theta(7)));
    model.param.set('eE33', sprintf('%g [C/m^2]', theta(8)));
    
    % --- Relative Permittivity (Dielectric) Matrix [eps^S] ---
    % eps33r includes dielectric dissipation captured by theta(16)
    model.param.set('epsilon11r', sprintf('%g', theta(9)));
    model.param.set('epsilon33r', sprintf('%g + %g*i', theta(10), theta(16)));
    
    %% 3. Set Macroscopic Properties
    model.param.set('th', sprintf('%g [mm]', h)); 
    model.param.set('rho_m', sprintf('%g [kg/m^3]', rho_PZT));
    
    %% 4. Execute FEM Study
    % The 'progress', 'off' flag is highly recommended when calling this 
    % function from a loop (like fminsearch) to save CPU overhead.
    mphrun(model, 'study', 'progress', 'off');
   
    %% 5. Extract Results
    % Evaluate the global admittance at the terminal/electrodes.
    % The transpose operator (') converts the output into a column vector.
    Y_ref = mphglobal(model, 'es.Y11')';
    
    % Derive Impedance directly from Admittance
    % This is mathematically equivalent to mphglobal(model, '1/es.Y11')'
    % but avoids calling the COMSOL server twice, saving execution time.
    Z_ref = 1 ./ Y_ref; 
end

