% file stress_power.m
% brief contains function imr.physics.stress_power

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the medium heating by the stress,
% taugradu = Br*Phi with Phi = tau:D, the full stress power, at the medium
% grid points. With a = R^2*Rdot/r^3 = (Rdot/R)*yT^-3 the local strain rate
% and s = tau_rr - tau_thetatheta, tau:D = -2*a*s. The heating is
%   stress 0        no stress, 0
%   stress 1, 2     Kelvin-Voigt: 12*mu*a^2 - 2*a*se, with the viscosity
%                   mu = 1/Re8 + fnu*iDRe acting instantly and the
%                   neo-Hookean (1) or quadratic (2) elastic stress
%                   difference se
%   stress 3, 4     same instantaneous Kelvin-Voigt dashpot, with the
%                   2nd-order Fung (3) or Mooney-Rivlin (4) elastic
%                   stress difference se in place of the neo-Hookean one
%   stress 5-11     12*mu_s*a^2 - 2*a*sp, with the solvent viscosity
%                   mu_s = LAM/Re8 + fnu*iDRe and the local polymer stress
%                   difference sp, which carries the stress history: the
%                   Zener field of imr.physics.zener_field (stress 5, 6, 7, 8) or the
%                   spectral stress field tau_rr - tau_thetatheta
%                   (stress 9, 10, 11)
% fnu is the non-Newtonian viscosity factor at the medium points (0 for a
% Newtonian run) and iDRe the factor of its stress term, 1/DRe for the
% Kelvin-Voigt family (1-4) and LAM/DRe otherwise, as in the radial
% equation. C2 is the Mooney-Rivlin second-modulus ratio (unused outside
% stress 4, 8). sp is the polymer stress difference at the medium
% points, needed for stress 5-11 and ignored otherwise. R and Rdot
% may be row vectors with sp holding one state per column, giving one
% column of taugradu per state.
function [taugradu] = stress_power(stress,Req,R,Rdot,Ca,Br,Re8, ...
    LAM,alphax,C2,iDRe,yT2,yT3,iyT3,iyT4,iyT6,fnu,sp)

if stress == 0
    % no stress
    taugradu = 0;

elseif stress >= 1 && stress <= 4
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
    elseif stress == 2 || stress == 3
        % quadratic neo-Hookean (2), factor 1 + alphax*(I1 - 3); 2nd-order
        % Fung (3) adds the next term of the same series,
        % (alphax^2/2)*(I1 - 3)^2
        x = x4.*iyT4 + 2*yT2.*ix2 - 3;
        factor = 1 + alphax*x;
        if stress == 3
            factor = factor + 0.5*alphax.^2.*x.^2;
        end
        taugradu = 12*Br*visc.*(Rdot./R).^2.*iyT6 + ...
            2*Br/Ca*iyT3.*(Rdot./R).*(yT2.*ix2 - iyT4.*x4).*factor;
    else
        % Mooney-Rivlin (4): neo-Hookean term plus the C2-scaled term of
        % the second invariant, lam^4 - lam^-2
        lam4 = (yT2.*ix2).^2;
        lam2i = x2./yT2;
        taugradu = 12*Br*visc.*(Rdot./R).^2.*iyT6 + ...
            2*Br/Ca*iyT3.*(Rdot./R).*(yT2.*ix2 - iyT4.*x4) + ...
            2*Br*C2/Ca*iyT3.*(Rdot./R).*(lam4 - lam2i);
    end

else
    % stress 5-11: the solvent part acts instantly, the polymer part
    % through its local stress difference sp, which remembers the history
    if isempty(sp)
        error(['imr.physics.stress_power: stress = %d needs the local polymer ' ...
            'stress field for the medium heating'], stress);
    end
    visc = LAM/Re8 + fnu.*iDRe;
    taugradu = 12*Br*visc.*(Rdot./R).^2.*iyT6 - 2*Br*(Rdot./R).*iyT3.*sp;
end

end
