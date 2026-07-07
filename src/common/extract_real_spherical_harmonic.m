function Y = extract_real_spherical_harmonic(x, phi, l, m)
%EXTRACT_REAL_SPHERICAL_HARMONIC Real spherical harmonic normalization.
%   This matches Common_functions/spherical_harmonic.m for the m = 0 and
%   positive-m cases used by the extraction routines.

    validateattributes(l, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
        mfilename, 'l');
    validateattributes(m, {'numeric'}, {'scalar', 'integer'}, mfilename, 'm');

    if abs(m) > l
        error('%s:InvalidMode', mfilename, 'Require abs(m) <= l.');
    end

    if m == 0
        prefactor = sqrt((2 .* l + 1) ./ (4 .* pi));
        Y = prefactor .* associated_legendre_value(x, l, 0);
    elseif m > 0
        prefactor = sqrt((2 .* l + 1) ./ (4 .* pi) .* ...
            factorial(l - m) ./ factorial(l + m));
        Y = (-1) .^ m .* sqrt(2) .* prefactor .* ...
            associated_legendre_value(x, l, m) .* cos(m .* phi);
    else
        m_abs = abs(m);
        prefactor = sqrt((2 .* l + 1) ./ (4 .* pi) .* ...
            factorial(l - m_abs) ./ factorial(l + m_abs));
        Y = (-1) .^ m .* sqrt(2) .* prefactor .* ...
            associated_legendre_value(x, l, m_abs) .* sin(m_abs .* phi);
    end
end

function P = associated_legendre_value(x, l, m)
    if m == 0 && l == 0
        P = ones(size(x));
    elseif l == 1 && m == 0
        P = x;
    elseif l == 1 && m == 1
        P = -(1 - x .^ 2) .^ (1 / 2);
    elseif m == l
        P = (-1) ^ l .* double_factorial_value(2 .* l - 1) .* ...
            (1 - x .^ 2) .^ (l / 2);
    elseif m == l - 1
        P = x .* (2 .* l - 1) .* associated_legendre_value(x, l - 1, l - 1);
    else
        P = ((2 .* l - 1) .* x .* associated_legendre_value(x, l - 1, m) - ...
            (l + m - 1) .* associated_legendre_value(x, l - 2, m)) ./ (l - m);
    end
end

function result = double_factorial_value(n)
    if n <= 0
        result = 1;
        return
    end

    result = 1;
    for k = n:-2:1
        result = result .* k;
    end
end
