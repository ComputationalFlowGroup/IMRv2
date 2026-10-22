% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Exponential inverse-temperature law and its exact temperature derivative.
% Use SI or dimensionless T consistently with c1 and c2; no conversion is
% performed. Scalar coefficients and arbitrary temperature arrays are valid.
function [P,dPdT] = f_saturation_pressure(T,c1,c2)
    P = c1*exp(c2./T);
    if nargout > 1
        dPdT = -c2*P./T.^2;
    end
end
