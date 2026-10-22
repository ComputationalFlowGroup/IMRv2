% file stress_power.m
% brief contains function imr.physics.stress_power

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the medium heating by the stress,
% taugradu = Br*Phi with Phi = tau:D, the full stress power, at the medium
% grid points. With a = R^2*Rdot/r^3 = (Rdot/R)*yT^-3 the local strain rate
% and s = tau_rr - tau_thetatheta, tau:D = -2*a*s. The heating is
%   stress 0     no stress, 0
%   stress 1, 2  Kelvin-Voigt: 12*mu*a^2 - 2*a*se, with the viscosity
%                mu = 1/Re8 + fnu*iDRe acting instantly and the neo-Hookean
%                (1) or quadratic (2) elastic stress difference se
%   stress 3-5   12*mu_s*a^2 - 2*a*sp, with the solvent viscosity
%                mu_s = LAM/Re8 + fnu*iDRe and the local polymer stress
%                difference sp, which carries the stress history: the Zener
%                field of imr.physics.zener_field (stress 3, 4) or the spectral stress
%                field tau_rr - tau_thetatheta (stress 5)
% fnu is the non-Newtonian viscosity factor at the medium points (0 for a
% Newtonian run) and iDRe the factor of its stress term, 1/DRe for stress
% 1 and 2 and LAM/DRe for stress 3 to 5, as in the radial equation. sp is
% the polymer stress difference at the medium points, needed for stress 3
% to 5 and ignored otherwise. R and Rdot may be row vectors with sp holding
% one state per column, giving one column of taugradu per state.
function [taugradu] = stress_power(stress,Req,R,Rdot,Ca,Br,Re8, ...
    LAM,alphax,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu,sp)

if stress == 0
    % no stress
    taugradu = 0;

elseif stress == 1 || stress == 2
    % Kelvin-Voigt: the viscous and elastic stresses act instantly
    visc = 1/Re8 + fnu.*iDRe;
    % current amplification factor
    Rst = Req./R;
    % incompressible condition, x2 = (r0/R)^2 for the reference position r0
    x2 = (yT3-1+Rst.^3).^(2/3);
    % inversion square of the reference coordinate
    ix2 = x2.^-1;
    % fourth power of the reference coordinate
    x4 = x2.^2;
    % with the hoop stretch lam = r/r0, yT2.*ix2 = lam^2 and
    % iyT4.*x4 = lam^-4, so the second term is -2*a*se
    if stress == 1
        % neo-Hookean
        taugradu = 12*Br*visc.*(Rdot./R).^2.*iyT6 + ...
            2*Br/Ca*iyT3.*(Rdot./R).*(yT2.*ix2 - iyT4.*x4);
    else
        % quadratic neo-Hookean, factor 1 + alphax*(I1 - 3)
        taugradu = 12*Br*visc.*(Rdot./R).^2.*iyT6 + ...
            2*Br/Ca*iyT3.*(Rdot./R).*(yT2.*ix2 - iyT4.*x4) .* ...
            (1+alphax*(x4.*iyT4 + 2*yT2.*ix2 - 3));
    end

else
    % stress 3-5: the solvent part acts instantly, the polymer part through
    % its local stress difference sp, which remembers the history
    if isempty(sp)
        error(['imr.physics.stress_power: stress = %d needs the local polymer ' ...
            'stress field for the medium heating'], stress);
    end
    visc = LAM/Re8 + fnu.*iDRe;
    taugradu = 12*Br*visc.*(Rdot./R).^2.*iyT6 - 2*Br*(Rdot./R).*iyT3.*sp;
end

end
