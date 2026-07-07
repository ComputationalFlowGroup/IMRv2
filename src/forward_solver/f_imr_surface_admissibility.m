% file f_imr_surface_admissibility.m
% brief minimum shape factor for a real spherical-harmonic perturbation
function geom = f_imr_surface_admissibility(n, m, ep, varargin)
%F_IMR_SURFACE_ADMISSIBILITY Evaluate min(1 + sum eps_i Y_i^m).
%
%   geom = f_imr_surface_admissibility(n, m, ep)
%
%   n, m, and ep are vectors. For a matrix ep, each column is evaluated as
%   one perturbation state.

    [ntheta, nphi] = parse_options(varargin);

    n = n(:);
    m = m(:);
    if isempty(m)
        m = zeros(size(n));
    end
    if size(ep,1) ~= length(n) && size(ep,2) == length(n)
        ep = ep.';
    end
    if size(ep,1) ~= length(n)
        error('ep must have one row per mode.');
    end

    theta = linspace(0, pi, ntheta);
    if all(m == 0)
        phi = 0;
    else
        phi = linspace(0, 2*pi, nphi);
    end
    [TH, PH] = ndgrid(theta, phi);

    Y = zeros(numel(TH), length(n));
    for k = 1:length(n)
        Y(:,k) = real_spherical_harmonic(n(k), m(k), TH, PH);
    end

    shape = 1 + Y*ep;
    min_shape = min(shape, [], 1);

    geom = struct();
    geom.min_shape = min_shape;
    geom.admissible = min_shape > 0;
    geom.min_over_states = min(min_shape);
end

function [ntheta, nphi] = parse_options(args)
    ntheta = 721;
    nphi = 361;
    if mod(length(args),2) == 1
        error('Optional inputs must be name/value pairs.');
    end
    for k = 1:2:length(args)
        switch lower(args{k})
            case 'theta_points'
                ntheta = args{k+1};
            case 'phi_points'
                nphi = args{k+1};
        end
    end
end

function Y = real_spherical_harmonic(n, m, theta, phi)
    x = cos(theta(:)).';
    abs_m = abs(m);
    P_all = legendre(n, x);
    Pnm = reshape(P_all(abs_m+1,:), size(theta));
    norm_factor = sqrt((2*n + 1)/(4*pi) * ...
        exp(gammaln(n - abs_m + 1) - gammaln(n + abs_m + 1)));

    if m == 0
        Y = norm_factor.*Pnm;
    elseif m > 0
        Y = sqrt(2).*norm_factor.*Pnm.*cos(abs_m.*phi);
    else
        Y = sqrt(2).*norm_factor.*Pnm.*sin(abs_m.*phi);
    end
    Y = Y(:);
end
