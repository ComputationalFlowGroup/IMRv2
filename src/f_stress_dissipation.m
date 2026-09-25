% file f_stress_dissipation.m
% brief contains function f_stress_dissipation

% Copyright (C) Mauro Rodriguez
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the stress dissipation term, \tau : \nabla u.
% The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [taugradu] = f_stress_dissipation(stress,spectral,Req,R,Rdot,Ca,Br, ...
    Re8,alphax,yT2,yT3,iyT3,iyT4,iyT6,X,ZZT,ivisco1,ivisco2,fnu,DRe)

% effective viscosity 1/Re8 + f/DRe, with the non-Newtonian factor f at
% the medium grid points (fnu is a scalar 0 for Newtonian runs)
if DRe == 0
    visc = 1/Re8;
else
    visc = 1/Re8 + fnu/DRe;
end
% spectral solution approach to compute the mechanical dissipation from
% the stress field
if spectral == 1
    taugradu = -2*Br*Rdot./(R*yT3).* ...
        (ZZT*(X(ivisco1)-X(ivisco2)));
    
    % no stress
elseif stress == 0
    taugradu = 0;
    
else
    % current amplification factor
    Rst = Req/R;
    % incompressible condition
    x2 = (yT3-1+Rst.^3).^(2/3);
    % inversion square of the reference coordinate
    ix2 = x2.^-1;
    % fourth power of the reference coordinate
    x4 = x2.^2;
    
    % neo-Hookean elasticity: Kelvin-Voigt (1), Maxwell, Jeffreys, Zener
    % (3), and UCM and Oldroyd-B in the Kelvin-Voigt closed form (5)
    if stress == 1 || stress == 3 || stress == 5
        taugradu = 12*Br*visc.*(Rdot/R)^2.*iyT6 + ...
            2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4);
        
        % quadratic neo-Hookean elasticity: Kelvin-Voigt (2) and Maxwell,
        % Jeffreys, Zener (4)
    elseif stress == 2 || stress == 4
        taugradu = 12*Br*visc.*(Rdot/R)^2.*iyT6 + ...
            2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4) .* ...
            (1+alphax*(x4.*iyT4 + 2*yT2.*ix2 - 3));
    end
end

end
