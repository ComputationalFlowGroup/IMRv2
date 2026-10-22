% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Inverse Kirchhoff map for K(T)=a*T+b, theta(1)=0 and positive a.
% Inputs use compatible scalar/array shapes; outputs follow implicit
% expansion. K is the conductivity, dT/dtheta=1/K. Optional derivatives
% with respect to a and b hold theta fixed (needed by mixture Jacobians).
function [T,K,dTda,dTdb] = f_thermal_inverse(theta,a,b)
    K = sqrt((a+b).^2+2*a.*theta);
    T = (K-b)./a;
    if nargout > 2
        dTda = ((a+b+theta)./K-T)./a;
        dTdb = ((a+b)./K-1)./a;
    end
end
