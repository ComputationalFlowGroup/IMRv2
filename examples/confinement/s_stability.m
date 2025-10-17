% bubble_stability.m
% Stability analysis of bubble radius in hyperelastic confinement

clear; clc;

%% Physical and Material Constants
sigma = 0.072;              % Surface tension [N/m]
T = 300;                    % Temperature [K]
kB = 1.38e-23;              % Boltzmann constant [J/K]
G = 1e5;                    % Shear modulus of solid [Pa]

% Geometry
R0 = 1e-6;                  % Initial bubble radius [m]
Rc = 20e-6;                 % Outer confinement radius [m]

% Reference radius R* and typical energy scale
deltaP0 = -1e5;             % Driving pressure [Pa]
R_star = -2 * sigma / deltaP0;
F_star = (4/3) * pi * R_star^2 * sigma;

% Range of Ng (gas amount) for bifurcation diagram
Ng_array = linspace(0.5e9, 5e9, 100); % Number of gas molecules
R_guess = linspace(R0*1.01, Rc*0.99, 500); % Range of R to test

Req_stable = [];
Req_unstable = [];
Ng_stable = [];
Ng_unstable = [];

%% Loop over Ng values
for Ng = Ng_array
    for R = R_guess
        % First derivative of Free Energy
        dF1 = dF_dR(R, Ng, sigma, kB, T, R0, Rc, G);
        
        % Detect sign change in dF/dR (equilibrium point)
        persistent dF_prev R_prev
        if ~isempty(dF_prev) && dF1 * dF_prev < 0
            % Equilibrium found between R_prev and R
            Req = fzero(@(r) dF_dR(r, Ng, sigma, kB, T, R0, Rc, G), [R_prev, R]);
            d2F = d2F_dR2(Req, Ng, sigma, kB, T, R0, Rc, G);
            
            % Classify stability
            if d2F > 0
                Req_stable(end+1) = Req;
                Ng_stable(end+1) = Ng;
            else
                Req_unstable(end+1) = Req;
                Ng_unstable(end+1) = Ng;
            end
        end
        dF_prev = dF1;
        R_prev = R;
    end
    clear dF_prev R_prev
end

%% Plot bifurcation diagram
figure; hold on;
plot(Ng_stable, Req_stable * 1e6, 'b-', 'LineWidth', 2);     % Stable branch
plot(Ng_unstable, Req_unstable * 1e6, 'r--', 'LineWidth', 2); % Unstable branch
xlabel('Gas Amount \(N_g\)', 'Interpreter', 'latex');
ylabel('Equilibrium Radius \(R_{\text{eq}}\) [\mu m]', 'Interpreter', 'latex');
title('Bifurcation Diagram of Bubble Radius', 'Interpreter', 'latex');
legend('Stable', 'Unstable');
grid on;

function dF = dF_dR(R, Ng, sigma, kB, T, R0, Rc, G)
    % Compute dF/dR = surface + gas + solid
    surface_term = 8 * pi * sigma * R;
    gas_term = -3 * Ng * kB * T / R;
    
    integrand = @(r0) (1 + (R^3 - R0^3) ./ r0.^3).^(-1/3) ...
                    - (1 + (R^3 - R0^3) ./ r0.^3).^(-7/3);
    solid_term = 8 * pi * G * R^2 * integral(@(r0) integrand(r0) ./ r0, R0, Rc);
    
    dF = surface_term + gas_term + solid_term;
end

function d2F = d2F_dR2(R, Ng, sigma, kB, T, R0, Rc, G)
    % Compute second derivative of F w.r.t. R
    gas_term = 3 * Ng * kB * T / R^2;
    I = integral(@(r0) ((1 + (R^3 - R0^3)./r0.^3).^(-1/3) - ...
                        (1 + (R^3 - R0^3)./r0.^3).^(-7/3)) ./ r0, R0, Rc);
    du = @(r0) 3 * R^2 ./ r0.^3;
    u = @(r0) 1 + (R^3 - R0^3)./r0.^3;
    dI = integral(@(r0) ( -u(r0).^(-4/3) + 7*u(r0).^(-10/3) ) .* ...
                            du(r0) ./ r0.^4, R0, Rc);

    d2F = 8 * pi * sigma + gas_term + 8 * pi * G * (2 * R * I + R^4 * dI);
end




