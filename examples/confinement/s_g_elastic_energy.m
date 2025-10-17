% Parameters
G = 1e4;       % Shear modulus in Pa
R0 = 1e-4;     % Reference bubble radius in meters
R = 1.2e-3;    % Current bubble radius (expanded)

% Compute stored elastic energy
F = neoHookeanElasticEnergy(R, R0, G);

fprintf('Stored elastic energy: %.6f J\n', F);

% Define params struct with physical constants

params.sigma = 0.072;   % N/m example for water at room temp
params.Ng = 1e16;       % example gas particle number
params.kB = 1.38e-23;   % J/K
params.T = 298;         % K
params.Rref = R0/8;     % m, example ref radius
params.R0 = 1e-6;       % m, initial bubble radius
params.Rc = 1e-4;       % m, cavity radius
params.G = G;           % Pa, shear modulus

% Objective function for fzero or fsolve (numerical derivative)
F_prime = @(R) (totalFreeEnergy(R + 1e-9, params) - totalFreeEnergy(R - 1e-9, params)) / (2e-9);

% Initial guess near R0
R_guess = 1.1 * params.R0;

% Find root of dF/dR = 0
R_eq = fzero(F_prime, R_guess);

fprintf('Equilibrium radius R_eq = %.4e m\n', R_eq);

function [F_solid, lambda_c, lambda_0] = neoHookeanElasticEnergyConfined(R, R0, Rc, mu)
    % Compute stored elastic energy for a bubble expanding inside a
    % spherical cavity of radius Rc (reference).
    %
    % Inputs:
    %   R  - current bubble radius (R > R0)
    %   R0 - reference bubble radius
    %   Rc - reference cavity radius
    %   mu - shear modulus
    %
    % Outputs:
    %   F_solid - stored elastic energy
    %   lambda_c - hoop stretch at cavity wall
    %   lambda_0 - hoop stretch at bubble wall
    
    % Validate inputs
    if R <= R0
        error('Current bubble radius R must be greater than reference radius R0.');
    end
    if Rc <= R0
        error('Cavity radius Rc must be greater than bubble radius R0.');
    end
    
    % Compute delta R^3
    deltaR3 = R^3 - R0^3;
    
    % Compute hoop stretches
    lambda_0 = R / R0;
    lambda_c = ((deltaR3 + Rc^3)^(1/3)) / Rc;
    
    % Define integrand
    integrand = @(lambda) (lambda.^(-4) + 2*lambda.^2 - 3) .* (lambda.^2 ./ (1 - lambda.^3).^2);
    
    % Check integration limits order
    if lambda_c < lambda_0
        % Integrate from lambda_c to lambda_0
        integral_val = integral(integrand, lambda_c, lambda_0, 'RelTol',1e-9,'AbsTol',1e-12);
    else
        % If lambda_c >= lambda_0, no energy stored or negative? Just zero.
        integral_val = 0;
        warning('lambda_c >= lambda_0; no elastic energy stored.');
    end
    
    % Compute elastic energy
    F_solid = (2 * pi * mu / 3) * deltaR3 * integral_val;
end
function F_total = totalFreeEnergy(R, params)
    % Compute total free energy of the system at bubble radius R
    %
    % params: structure with fields
    %   sigma   : surface tension
    %   Ng      : number of gas particles
    %   kB      : Boltzmann constant
    %   T       : temperature
    %   Rref    : reference radius for gas energy
    %   R0      : reference bubble radius
    %   Rc      : reference cavity radius
    %   mu      : shear modulus of solid
    
    sigma = params.sigma;
    Ng = params.Ng;
    kB = params.kB;
    T = params.T;
    Rref = params.Rref;
    R0 = params.R0;
    Rc = params.Rc;
    mu = params.mu;
    
    % Surface energy
    F_sigma = 4 * pi * R^2 * sigma;
    
    % Gas energy
    F_g = -3 * Ng * kB * T * log(R / Rref);
    
    % Solid elastic energy
    [F_solid, ~, ~] = neoHookeanElasticEnergyConfined(R, R0, Rc, mu);
    
    % Total energy
    F_total = F_sigma + F_g + F_solid;
end

%%
function F_solid = neoHookeanElasticEnergy(R, R0, mu)
    % Compute the stored elastic energy of an incompressible neo-Hookean solid
    % around a spherical bubble expanding from radius R0 to R.
    %
    % Inputs:
    %   R  - current bubble radius (scalar)
    %   R0 - reference bubble radius (scalar)
    %   mu - shear modulus of the solid (scalar)
    %
    % Output:
    %   F_solid - stored elastic energy (scalar)
    
    % Check inputs
    if R <= R0
        error('Current radius R must be greater than reference radius R0.');
    end
    
    % Compute delta R cubed
    deltaR3 = R^3 - R0^3;
    
    % Define the lower integration limit lambda0
    lambda0 = R / R0;  % hoop stretch at bubble wall
    
    % Define the integrand function
    integrand = @(lambda) (lambda.^(-4) + 2*lambda.^2 - 3) .* (lambda.^2 ./ (1 - lambda.^3).^2);
    
    % Numerically evaluate the integral from lambda0 to 1
    % Note: Since lambda0 > 1 in expansion, integration bounds are reversed
    % So integrate from 1 to lambda0 and then negate
    integral_val = integral(integrand, lambdac, lambda0, 'RelTol',1e-9,'AbsTol',1e-12);
    integral_val = -integral_val;  % because limits reversed
    
    % Compute the stored elastic energy
    F_solid = (2 * pi * mu / 3) * deltaR3 * integral_val;
end

