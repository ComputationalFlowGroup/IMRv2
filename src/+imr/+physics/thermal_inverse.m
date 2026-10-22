% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Inverse Kirchhoff map for K(T)=a*T+b, theta(1)=0, a>=0, a+b>0.
% Inputs use compatible scalar/array shapes; outputs follow implicit
% expansion. K is the conductivity, dT/dtheta=1/K. Optional derivatives
% with respect to a and b hold theta fixed (needed by mixture Jacobians).
function [T,K,dTda,dTdb] = thermal_inverse(theta,a,b)
    K = sqrt((a+b).^2+2*a.*theta);
    % Rationalizing K-(a+b) avoids cancellation as a -> 0 or theta -> 0
    % and includes the constant-conductivity limit without division by a.
    deltaT = 2*theta./(K+a+b);
    T = 1+deltaT;
    if nargout > 2
        % Implicit derivatives of theta=a*(T^2-1)/2+b*(T-1).
        dTda = -deltaT.*(2+deltaT)./(2*K);
        dTdb = -deltaT./K;
    end
end
