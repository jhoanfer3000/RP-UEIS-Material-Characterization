function [f_window, Z_window, Y_window, phi_window] = window_frequency(f_exp, Z_exp_mag, Y_exp_mag, phi_exp,f0,x)
    % =====================================================================
    % WINDOW_FREQUENCY: Isolates the resonance region for optimization
    % =====================================================================
    % Description: Identifies the series resonance frequency (minimum |Z|) 
    % and crops the experimental datasets to a specific fractional bandwidth 
    % around it. It also applies signal processing techniques (unwrapping, 
    % detrending, and smoothing) to the phase data to ensure stable and 
    % noise-free cost function evaluations.
    %
    % Inputs:
    %   - f_exp: Full experimental frequency vector [Hz]
    %   - Z_exp_mag: Full experimental impedance magnitude |Z| [Ohms]
    %   - Y_exp_mag: Full experimental admittance magnitude |Y| [S]
    %   - phi_exp: Full experimental phase vector [Rad]
    %   - x: Fractional bandwidth for the window (e.g., 0.03 for +/- 3%)
    %
    % Outputs:
    %   - f_window: Cropped frequency vector
    %   - Z_window: Cropped impedance magnitude
    %   - Y_window: Cropped admittance magnitude
    %   - phi_window: Cropped, detrended, and smoothed phase vector
    % =====================================================================
    %% 1. Define Frequency Boundaries
    % Calculate the lower and upper bounds based on the fractional parameter 'x'
    f_min = f0 * (1 - x); % Lower bound (-x% of f0)
    f_max = f0 * (1 + x); % Upper bound (+x% of f0)
    
    %% 2. Apply Logical Masking
    % Create a logical array (mask) to filter the data points that fall
    % strictly within the defined frequency window.
    mask = (f_exp >= f_min) & (f_exp <= f_max);
    
    % Crop the magnitude vectors
    f_window = f_exp(mask);
    Z_window = Z_exp_mag(mask);
    Y_window = Y_exp_mag(mask);
    %% 3. Phase Signal Processing
    % Extract the phase data within the window and process it:
    %   1. unwrap(): Removes artificial 2*pi jumps caused by phase wrapping.
    %   2. detrend(): Removes any underlying linear trend from the phase 
    %                 background, isolating the dynamic resonance shift.
    phi_window = detrend(unwrap(phi_exp(mask)));
    
    % Apply a moving median filter (window size = 3) to smooth the phase.
    % A moving median is highly robust against local outliers or high-frequency 
    % noise spikes common in experimental phase measurements.
    phi_window = smoothdata(phi_window, 'movmedian', 3);   
end
