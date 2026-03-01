function C = cost_function_Z_PZT(pA, pB, pC, pD, model, f_exp, Z_exp, phi_exp, Fre_step, th, rho_m)
    % =====================================================================
    % COST_FUNCTION_Z: Objective function for PZT parameters optimization
    % =====================================================================
    % Description: Updates COMSOL parameters, evaluates the FEM model, and
    % calculates the logarithmic sum of squared residuals for impedance.
    % It supports multiple optimization stages:
    %   - Stage A: Adjusts pA (C11', C13', C33', C44') while pB is fixed.
    %   - Stage B: Adjusts pB (e15', e33', eps33r') while pA is fixed.
    % Features a real-time tracking dashboard with 3 subplots.
    % =====================================================================
    
    % --- Persistent variables to track cost history across iterations ---
    persistent eval_count cost_history
    if isempty(eval_count)
        eval_count = 0;
        cost_history = [];
    end
    
    %% 1. UPDATE COMSOL PARAMETERS
    Fre_I_Exp = f_exp(1);
    Fre_F_Exp = f_exp(end);
    
    % Set frequency sweep limits in COMSOL
    model.param.set('Fre_I', sprintf('%g [Hz]', Fre_I_Exp - Fre_step)); 
    model.param.set('Step_Fre', sprintf('%g [Hz]', Fre_step)); 
    model.param.set('Fre_F', sprintf('%g [Hz]', Fre_F_Exp + Fre_step)); 
    
    % Set physical material parameters
    % Using sprintf ensures complex numbers are formatted safely for COMSOL
    model.param.set('CE11', sprintf('%g + %g*i [Pa]', pA(1), pC(1)));
    model.param.set('CE12', sprintf('%g [Pa]', pD(1)));
    model.param.set('CE13', sprintf('%g + %g*i [Pa]', pA(2), pC(2)));
    model.param.set('CE33', sprintf('%g + %g*i [Pa]', pA(3), pC(3)));
    model.param.set('CE44', sprintf('%g [Pa]', pA(4)));
    
    model.param.set('eE15', sprintf('%g [C/m^2]', pB(1)));
    model.param.set('eE31', sprintf('%g [C/m^2]', pD(2)));
    model.param.set('eE33', sprintf('%g [C/m^2]', pB(2)));
    
    model.param.set('epsilon11r', sprintf('%g', pD(3)));
    model.param.set('epsilon33r', sprintf('%g', pB(3)));
    
    model.param.set('th', sprintf('%g [mm]', th)); 
    model.param.set('rho_m', sprintf('%g [kg/m^3]', rho_m));
    
    %% 2. RUN FEM MODEL & EXTRACT DATA
    mphrun(model, 'study', 'progress', 'off');
    
    % Extract complex admittance and convert to impedance
    Y_num = mphglobal(model, 'es.Y11'); % Get Admittance directly
    Z_num = 1 ./ Y_num;                 % Convert to Impedance
    fre_num = (Fre_I_Exp - Fre_step : Fre_step : Fre_F_Exp + Fre_step)'; % Column vector
    
    %% 3. INTERPOLATION & COST CALCULATION
    % Interpolate numerical data to match the experimental frequency points
    Z_num_interp = interp1(fre_num, Z_num(:), f_exp, 'linear', 'extrap');
    phi_num = detrend(unwrap(angle(Z_num_interp)));
    
    % Calculate magnitude residuals (Logarithmic scale)
    res_mag = log10(abs(Z_exp)) - log10(abs(Z_num_interp));
    
    % Sum of squared residuals (Cost function)
    R_vector =res_mag(:);
    C = sum(R_vector.^2);
    
    % Penalize if COMSOL fails (returns NaN)
    if isnan(C) || isinf(C)
        C = 1e6;
    end
    
    % --- Update History Tracker ---
    eval_count = eval_count + 1;
    cost_history(eval_count) = C;
    
    %% 4. REAL-TIME DASHBOARD (PLOTTING)
    line_w = 1.2;
    size_font = 12;
    
    % Setup figure and clear it to prevent text overlapping during loops
    f1 = figure(1);
    clf(f1); 
    f1.Position = [100, 100, 1300, 400]; % Wide figure for 3 plots
    
    % --- Plot 1: Impedance Magnitude ---
    subplot(1,3,1)
    plot(f_exp./1000, abs(Z_exp), '-r', 'LineWidth', line_w); hold on;
    plot(f_exp./1000, abs(Z_num_interp), '--b', 'LineWidth', line_w);
    set(gca, 'YScale', 'log')
    xlabel('Frequency (kHz)'); ylabel('Impedance (\Omega)');
    title('|Z| Spectrum');
    legend('Exp', 'Num', 'Location', 'best');
    ax = gca; ax.FontSize = size_font; grid on; grid minor; axis tight;
    
    % --- Plot 2: Impedance Phase ---
    subplot(1,3,2)
    plot(f_exp./1000, phi_exp, '-r', 'LineWidth', line_w); hold on;
    plot(f_exp./1000, phi_num, '--b', 'LineWidth', line_w);
    xlabel('Frequency (kHz)'); ylabel('Phase (Rad)');
    title('Phase Spectrum');
    legend('Exp', 'Num', 'Location', 'best');
    ax = gca; ax.FontSize = size_font; grid on; grid minor; axis tight;
    
    % --- Plot 3: Cost Function Evolution ---
    subplot(1,3,3)
    plot(1:eval_count, cost_history, '*b', 'LineWidth', 1.5);
    set(gca, 'YScale', 'log') % Log scale highlights convergence plateaus
    xlabel('Function Evaluations'); ylabel('Cost Value (C)');
    title(sprintf('Convergence (C = %.4f)', C));
    ax = gca; ax.FontSize = size_font; grid on; grid minor; axis tight;
    
    % --- Global Super Title (Dynamic Parameter Display) ---
    % Displays the parameters from BOTH stages, now correctly including C13' (pA(2)).
    param_str = sprintf(['Iteration: %d | Cost: %.4f\n' ...
                         'Stage A (Mech) -> C11'': %.2e | C13'': %.2e | C33'': %.2e | C44'': %.2e\n' ...
                         'Stage B (Piezo)-> e15'': %.2f | e33'': %.2f | eps33'': %.1f'], ...
                         eval_count, C, pA(1), pA(2), pA(3), pA(4), pB(1), pB(2), pB(3));
                     
    sgtitle(param_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0 0.2 0.5]);
    
    drawnow;
end

