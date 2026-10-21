% file f_mercury_pvsat.m
% brief interpolates measured mercury saturation pressure

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% T is in K, Pv is in Pa. Log pressure is linear in reciprocal temperature
% between the Fokin data nodes. No extrapolation or water fit is used.
% This is an ambient-state interpolation, not a dynamic mass-transfer law.
function Pv = f_mercury_pvsat(T)
    validateattributes(T, {'numeric'}, {'real','finite','nonempty'}, ...
        mfilename, 'T');
    p = f_mercury_parameters();
    nodes = p.vapor.saturation_T;
    if any(T(:) < nodes(1) | T(:) > nodes(end))
        error('IMR:mercurySaturation', ...
            'Mercury saturation data support 298.15-550 K only.');
    end
    x = 1./fliplr(nodes);
    y = log(fliplr(p.vapor.saturation_P));
    Pv = reshape(exp(interp1(x, y, 1./double(T(:)), 'linear')), size(T));
end
