function bvp_elastic_confinement_solver
    % Parameters
    mu = 1.0;           % Shear modulus
    lambda_lame = 10.0; % First Lamé parameter
    Rc = 1.0;           % Inner reference radius
    Rs = 2.0;           % Outer reference radius
    Pin = 0.5;          % Applied internal pressure (positive means compression)

    % Mesh
    r0_mesh = linspace(Rc, Rs, 100);

    % Initial guess: no deformation
    solinit = bvpinit(r0_mesh, @guess);

    % Solve BVP
    sol = bvp4c(@(r0, y) odefun(r0, y, mu, lambda_lame), ...
                @(ya, yb) bcfun(ya, yb, mu, lambda_lame, Rc, Rs, Pin), ...
                solinit);

    % Extract solution
    r0_vals = sol.x;
    r_vals = sol.y(1, :);

    % Plot
    figure;
    plot(r0_vals, r_vals, 'LineWidth', 2);
    xlabel('Reference Radius r_0');
    ylabel('Deformed Radius r(r_0)');
    title('Radial Deformation of Compressible Neo-Hookean Sphere');
    grid on;
end

% ------------------------------
% ODE function: y = [r; dr/dr0]
% ------------------------------
function dydr0 = odefun(r0, y, mu, lambda_lame)
    r = y(1);
    drdr0 = y(2);

    lambda_r = drdr0;
    lambda_theta = r ./ r0;
    J = lambda_r .* lambda_theta.^2;

    % Cauchy stresses
    sigma_rr = mu .* (lambda_r.^2 - 1) + lambda_lame .* log(J);
    sigma_tt = mu .* (lambda_theta.^2 - 1) + lambda_lame .* log(J);

    % Equilibrium equation (rewritten as 1st-order system)
    dsigma_dr0 = -(2 ./ r) .* (sigma_rr - sigma_tt) .* lambda_r;

    dydr0 = [drdr0; dsigma_dr0];
end

% ------------------------------------
% Boundary conditions
% ya and yb are y = [r; dr/dr0]
% ------------------------------------
function res = bcfun(ya, yb, mu, lambda_lame, Rc, Rs, Pin)
    % At inner boundary (r0 = Rc)
    r_in = ya(1);
    dr_in = ya(2);
    lambda_r_in = dr_in;
    lambda_theta_in = r_in / Rc;
    J_in = lambda_r_in * lambda_theta_in^2;
    sigma_rr_in = mu * (lambda_r_in^2 - 1) + lambda_lame * log(J_in);

    % At outer boundary (r0 = Rs)
    r_out = yb(1);
    dr_out = yb(2);
    lambda_r_out = dr_out;
    lambda_theta_out = r_out / Rs;
    J_out = lambda_r_out * lambda_theta_out^2;
    sigma_rr_out = mu * (lambda_r_out^2 - 1) + lambda_lame * log(J_out);

    % Boundary conditions: σ_rr(Rc) = -P_in, σ_rr(Rs) = 0
    res = [sigma_rr_in + Pin; sigma_rr_out];
end

% ------------------------------------
% Initial guess for [r; dr/dr0]
% ------------------------------------
function g = guess(r0)
    g = [r0; 1];  % undeformed state
end

