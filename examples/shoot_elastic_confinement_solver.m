function shoot_elastic_confinement_solver
    % PARAMETERS
    % R_c = 1.0;      % Reference cavity radius
    % R_s = 2.0;      % Reference outer radius
    % R = 1.2;        % Deformed cavity radius (input parameter)
    % mu = 1.0;       % Shear modulus
    % kappa = 10.0;   % Bulk moduluss
    format long;
    R_c = 1e-3;        % 1 mm
    R_s = 2e-3;        % 2 mm
    R   = 1.2e-3;      % 1.2 mm
    mu  = 1e3;         % 1 kPa
    kappa = 1e5;       % 100 kPa
    sigma = 0.05;      % N/m (typical for water-air interface)
    R_ref = R_c;
    p_ref = 101325;    % 1 atm
    Ng_kBT = p_ref * (4/3) * pi * R_ref^3;  % Joules

    % Solve using shooting method to satisfy r(R_s) = R_s
    dr0_guess = fzero(@(dr0_init) shoot(dr0_init, R_c, R_s, R, mu, kappa), [0.1, 10]);

    % Integrate using correct initial slope
    [r0_vals, y_vals] = ode45(@(r0, y) odesystem(r0, y, mu, kappa), [R_c R_s], [R; dr0_guess]);

    % Compute stress at inner wall
    r0 = R_c;
    r = R;
    dr_dr0 = dr0_guess;
    lambda_r = dr_dr0;
    lambda_theta = r / r0;
    J = lambda_r * lambda_theta^2;

    sigma_rr = mu * (lambda_r^2 - 1) + kappa * log(J);

    fprintf('Radial stress at cavity wall: %.4f\n', sigma_rr);

    % Extract solution
    r_vals = y_vals(:,1);
    dr_vals = y_vals(:,2);

    % Compute stretches
    lambda_r = dr_vals;
    lambda_theta = r_vals ./ r0_vals;
    J = lambda_r .* lambda_theta.^2;

    % Compute strain energy density at each point
    W = 0.5 * mu .* (lambda_r.^2 + 2 .* lambda_theta.^2 - 3) ...
        - mu .* log(J) + 0.5 * kappa .* (log(J)).^2;
   
    % Compute total elastic energy by numerical integration
    F_elastic = trapz(r0_vals, 4 * pi .* r0_vals.^2 .* W);
    fprintf('Total elastic energy: %.6f\n', F_elastic);
  
    % Compute additional energies
    F_sigma = 4 * pi * R^2 * sigma;
    F_gas = -3 * Ng_kBT * log(R / R_ref);
  
    % Compute total energy
    F_total = F_elastic + F_sigma + F_gas;

    fprintf('Surface energy:       %.6f\n', F_sigma);
    fprintf('Gas energy:           %.6f\n', F_gas);
    fprintf('Total free energy:    %.6f\n', F_total);

end

function F = shoot(dr0_init, R_c, R_s, R, mu, kappa)
    % Integrate and return mismatch at r0 = R_s
    [~, y] = ode45(@(r0, y) odesystem(r0, y, mu, kappa), [R_c R_s], [R; dr0_init]);
    r_s_deformed = y(end,1);
    F = r_s_deformed - R_s;
end

function dydr0 = odesystem(r0, y, mu, kappa)
    % y(1) = r(r0), y(2) = dr/dr0
    r = y(1);
    dr_dr0 = y(2);

    lambda_r = dr_dr0;
    lambda_theta = r / r0;
    J = lambda_r * lambda_theta^2;

    % Cauchy stresses
    sigma_rr = mu * (lambda_r^2 - 1) + kappa * log(J);
    sigma_tt = mu * (lambda_theta^2 - 1) + kappa * log(J);

    % Convert to derivative in r0-space (dσ_rr/dr0)
    d_sigma_rr_dr0 = ( ...
        -2 / r * (sigma_rr - sigma_tt) ) * dr_dr0;

    % Return derivatives: [dr/dr0; d²r/dr0²]
    d2r_dr0_2 = d_sigma_rr_dr0 / ...
        (2 * mu * lambda_r + 2 * kappa * (1 / lambda_r));

    dydr0 = [dr_dr0; d2r_dr0_2];
end

