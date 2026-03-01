function C = cost_function_Y(theta_s, theta, model, f_exp, Z_exp, phi_exp, ...
                             Fre_step, th, rho_m, d_ma, th_ma, rho_ma, ...
                             th_Gel, rho_Gel, c_Gel)
    % =====================================================================
    % COST_FUNCTION_Y: Objective function for Tested Material Parameters
    % =====================================================================
    % Description: Updates COMSOL parameters for the RP-UEIS assembly, 
    % evaluates the FEM model, and calculates the linear sum of squared 
    % residuals for ADMITTANCE (|Y|). 
    % It optimizes ALL complex bulk (K) and shear (G) moduli of the coupled 
    % material simultaneously.
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
    
    % Set BASE PZT physical material parameters (Fixed during this stage)
    model.param.set('CE11', sprintf('%g + %g*i [Pa]', theta(1), theta(11)));
    model.param.set('CE12', sprintf('%g [Pa]', theta(2)));
    model.param.set('CE13', sprintf('%g + %g*i [Pa]', theta(3), theta(12)));
    model.param.set('CE33', sprintf('%g + %g*i [Pa]', theta(4), theta(13)));
    model.param.set('CE44', sprintf('%g [Pa]', theta(5)));
    
    model.param.set('eE15', sprintf('%g [C/m^2]', theta(6)));
    model.param.set('eE31', sprintf('%g [C/m^2]', theta(7)));
    model.param.set('eE33', sprintf('%g [C/m^2]', theta(8)));
    
    model.param.set('epsilon11r', sprintf('%g', theta(9)));
    model.param.set('epsilon33r', sprintf('%g', theta(10)));
    
    model.param.set('th', sprintf('%g [mm]', th)); 
    model.param.set('rho_m', sprintf('%g [kg/m^3]', rho_m));
    
    % Set TESTED MATERIAL and COUPLING GEL parameters (Node: par2)
    % theta_s = [K_real, G_real, K_imag, G_imag]
    model.param("par2").set("Km", sprintf('%g + %g*i [Pa]', theta_s(1), theta_s(3)));
    model.param("par2").set("Gm", sprintf('%g + %g*i [Pa]', theta_s(2), theta_s(4)));
    
    % Geometry and Density for the assembly
    model.param("par2").set("th_Gel", sprintf('%g [um]', th_Gel));
    model.param("par2").set("rho_so", sprintf('%g [kg/m^3]', rho_ma));
    model.param("par2").set("th_so", sprintf('%g [mm]', th_ma));
    model.param("par2").set("d_so", sprintf('%g [mm]', d_ma));
    model.param("par2").set("rho_gel", sprintf('%g [kg/m^3]', rho_Gel));
    model.param("par2").set("c_gel", sprintf('%g [m/s]', c_Gel));
    
    %% 2. RUN FEM MODEL & EXTRACT DATA
    mphrun(model, 'study', 'progress', 'off');
    
    % Extract complex Admittance directly from COMSOL
    Y_num = mphglobal(model, 'es.Y11'); 
    fre_num = (Fre_I_Exp - Fre_step : Fre_step : Fre_F_Exp + Fre_step)'; % Column vector
    
    %% 3. INTERPOLATION & COST CALCULATION
    % Interpolate numerical admittance to match experimental frequency points
    Y_num_interp = interp1(fre_num, Y_num(:), f_exp, 'linear', 'extrap');
    
    % Extract phase of Admittance (Y). 
    phi_num = detrend(unwrap(angle(Y_num_interp)));
    
    % Derive experimental Admittance from Impedance
    Y_exp = 1 ./ Z_exp;
    
    % Calculate magnitude residuals (Linear scale for Admittance)
    res_mag = abs(Y_exp) - abs(Y_num_interp);
    
    % Sum of squared residuals (Cost function)
    R_vector = res_mag(:);
    C = sum(R_vector.^2);
    
    % Penalize if COMSOL fails (returns NaN or Inf) to prevent optimizer crash
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
    plot(f_exp./1000, abs(1 ./ Y_num_interp), '--b', 'LineWidth', line_w);
    set(gca, 'YScale', 'log')
    xlabel('Frequency (kHz)'); ylabel('Impedance (\Omega)');
    title('|Z| Spectrum');
    legend('Exp', 'Num', 'Location', 'best');
    ax = gca; ax.FontSize = size_font; grid on; grid minor; axis tight;
    
    % --- Plot 2: Impedance Phase ---
    subplot(1,3,2)
    plot(f_exp./1000, phi_exp, '-r', 'LineWidth', line_w); hold on;
    plot(f_exp./1000, -phi_num, '--b', 'LineWidth', line_w);
    xlabel('Frequency (kHz)'); ylabel('Phase (Rad)');
    title('Phase Spectrum');
    legend('Exp', 'Num', 'Location', 'best');
    ax = gca; ax.FontSize = size_font; grid on; grid minor; axis tight;
    
    % --- Plot 3: Cost Function Evolution ---
    subplot(1,3,3)
    plot(1:eval_count, cost_history, '*b', 'LineWidth', 1.5); % Cambiado a línea sólida negra
    set(gca, 'YScale', 'log') 
    xlabel('Function Evaluations'); ylabel('Cost Value (C)');
    title(sprintf('Convergence (C = %.4e)', C));
    ax = gca; ax.FontSize = size_font; grid on; grid minor; axis tight;
    
    % --- Global Super Title (Dynamic Parameter Display) ---
    % Displays the parameters of the Tested Material (K', K'', G', G'').
    param_str = sprintf(['Iteration: %d | Cost: %.4e\n' ...
                         'Tested Material -> K'': %.2e | G'': %.2e | K'''': %.2e | G'''': %.2e'], ...
                         eval_count, C, theta_s(1), theta_s(2), theta_s(3), theta_s(4));
                     
    sgtitle(param_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.1 0.4 0.1]); 
    
    drawnow;
end