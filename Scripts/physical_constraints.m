function [c, ceq] = physical_constraints(x, rho_s,poisson_min, poisson_max)
    % x is the parameter vector: [K', G',K'', G'']
    Kp  = x(1);
    Gp  = x(2);
    Kpp = x(3);
    Gpp = x(4);

    % 1. Base calculations for Moduli and Sound Speeds
    Lp  = Kp + (4/3)*Gp;
    Lpp = Kpp + (4/3)*Gpp;

    cl = sqrt(Lp / rho_s);
    cs = sqrt(Gp / rho_s);

    % 2. Loss Tangents
    tan_delta_l = Lpp / Lp;
    tan_delta_s = Gpp / Gp;

    % 3. CONSTRAINT 1: Shear attenuation must be greater than longitudinal attenuation
    % We require alpha_l < alpha_s, which translates to alpha_l - alpha_s <= 0 for fmincon
    c_attenuation = (tan_delta_l / cl) - (tan_delta_s / cs);

    % 4. CONSTRAINT 2: Physical Poisson's ratio (e.g., restricted between 0.30 and 0.45)
    nu = (3*Kp - 2*Gp) / (2*(3*Kp + Gp));
    c_poisson_max = nu - poisson_max; % Forces nu <= 0.49
    c_poisson_min = poisson_min - nu; % Forces nu >= 0.30

    % Pack all inequality constraints (fmincon requires c <= 0)
    c = [c_attenuation; c_poisson_max; c_poisson_min];
    
    % No equality constraints
    ceq = []; 
end

