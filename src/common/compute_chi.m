function I = compute_chi(l, m)
    % High-accuracy quadrature points
    N_theta = 500;
    N_phi = 500;
    
    theta = linspace(0, pi, N_theta);
    phi = linspace(0, 2*pi, N_phi);
    [Theta, Phi] = meshgrid(theta, phi);
    
    x = cos(Theta);  % x = cos(theta)

    % Compute spherical harmonic and second derivative w.r.t theta
    Y = spherical_harmonic(x, Phi, l, m);
    
    % Numerical second derivative with respect to theta
    dtheta = theta(2) - theta(1);
    Y_d2 = diff(Y, 2, 2) / dtheta^2;
    Y = Y(:, 2:end-1);  % Match shape after diff
    
    Theta_mid = Theta(:, 2:end-1);  % updated theta for weight
    integrand = sin(Theta_mid) .* Y .* Y_d2;
    
    % Integral using trapezoidal rule
    I_theta_phi = trapz(phi, trapz(theta(2:end-1), integrand, 2));
    I = real(I_theta_phi);  % Should be real
end
