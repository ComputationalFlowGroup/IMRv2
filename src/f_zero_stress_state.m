% file f_zero_stress_gm.m
% brief contains function f_zero_stress_gm

% brief This function uses a Newton solver to find and plot the zero stress
% state for a neo-Hookean Kelvin-Voigt graded material.
function [] = zero_stress_state()

% Parameters
dGdhs_values = linspace(0, 10, 100);
hs_values = zeros(size(dGdhs_values));
tol = 1e-6;
max_iter = 100;

% Loop through each value 
for i = 1:length(dGdhs_values)
    dGdhs = dGdhs_values(i);
    hs_values(i) = newton_solver(dGdhs, 1e-3, 1e3);
end

% Plot hs as a function of dGdhs
figure;
plot(dGdhs_values, hs_values, 'LineWidth',2);
xlabel('dGdhs');
ylabel('hs');
title('Variation of dGdhs on hs');
grid on;
end

function hs = newton_solver(dGdhs, tol, max_iter)
% neo-Hookean Kelvin-Voigt gm with linear profile definitions
f = @(hs) -2.5 - (1/3)*dGdhs + 0.5*hs^4 + 2*hs + (1/3)*dGdhs*hs^3 + dGdhs*log(hs);
f_prime = @(hs) -2*hs^3 + 4 + dGdhs*hs^2 + hs/dGdhs;

% Initial guess for hs
hs = 0.01;

% Newton's method iterations
for iter = 1:max_iter
   
    % Compute the function value and its derivative at the current hs
    f_val = f(hs);
    f_prime_val = f_prime(hs);

    % Update the value of hs using Newton's method
    hs_new = hs - f_val / f_prime_val;

    % Check for convergence
    if abs(hs_new - hs) < tol
        hs = hs_new;
        return;
    end

    % update hs for the next iteration
    hs = hs_new;
end

% Display warning if exit the loop without convergence
warning('Newton method did not converge within the max iterations given.');
end