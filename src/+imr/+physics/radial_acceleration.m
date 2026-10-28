% file radial_acceleration.m
% brief contains function imr.physics.radial_acceleration

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features the spherical radial bubble dynamics
% equations used to compute the bubble wall acceleration. The radial models
% Rayleigh-Plesset, Keller-Miksis in pressure and enthalpy, and Gilmore are
% available in the solver.
function [Rddot] = radial_acceleration(radial, P, Pdot, Pf8, Pf8dot, iWe, ...
    R, Rdot, S, Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
    ddintfnu, iDRe )

% ddintfnu is the Rddot coefficient of d/dt(intfnu); the stress enters as
% -6*intfnu*iDRe, so its Rddot part moves to the denominator as
% +6*R*ddintfnu*iDRe/C, mirroring JdotA

% Rayleigh-Plesset
if radial == 1
    Rddot = (P - 1 - Pf8 - iWe./R + S - 1.5*Rdot.^2)./R;

    % Keller-Miksis in pressure
elseif radial == 2
    Rddot = ((1+Rdot./Cstar).*(P - 1 - Pf8 - iWe./R + S) ...
        + R/Cstar.*(Pdot + iWe*Rdot./R.^2 +  Sdot - Pf8dot) ...
        - 1.5*(1-Rdot/(3*Cstar)).*Rdot.^2)./((1-Rdot/Cstar).*R + JdotA/Cstar ...
        + 6*R.*ddintfnu*iDRe/Cstar);

    % enthalpy forms: Keller-Miksis (3, 5) with the far-field sound speed
    % and Gilmore (4, 6) with the local one, for the Tait (3, 4) or
    % Mie-Gruneisen (5, 6) equation of state
else
    if radial <= 4
        % Tait EoS
        Pb = P - iWe./R + GAMa + S;
        logB = log1p((P - iWe./R + S - 1)/sam);
        logF = log1p(Pf8/sam);
        % The nstate = 1 limit is logarithmic. expm1 also avoids loss of
        % significance for nstate close to one or pressure near ambient.
        if no == 0
            hB = sam*logB;
            hF = sam*logF;
        else
            hB = (sam/no)*expm1(no*logB);
            hF = (sam/no)*expm1(no*logF);
        end
        hH = (sam./Pb).^(1/nstate);
        % The driving pressure belongs to the same barotrope as the wall.
        % Subtracting Pf8 itself assumes unit far-field density at all P.
        Pfar = 1 + Pf8 + GAMa;
        hFprime = (sam./Pfar).^(1/nstate);
        if radial == 4
            rho = (Pb/sam).^(1/nstate);
            C = sqrt(nstate*Pb./rho);
        else
            C = Cstar;
        end
    else
        % Mie-Gruneisen EoS (imr.physics.mie_gruneisen_eos, the general
        % dimensionless p(rho,e)-style evaluator; see that file)
        Pb = P - iWe./R + S;
        [C, hB, hH] = imr.physics.mie_gruneisen_eos(Pb, Cstar, hugoniot_s, nog);
        [~, hF, hFprime] = imr.physics.mie_gruneisen_eos(1 + Pf8, Cstar, hugoniot_s, nog);
        if radial == 5
            C = Cstar;
        end
    end
    Rddot = ((1 + Rdot./C).*(hB - hF) - R./C.*hFprime.*Pf8dot ...
        + R./C.*hH.*(Pdot + iWe*Rdot./R.^2 + Sdot) ...
        - 1.5*(1 - Rdot./(3*C)).*Rdot.^2) ./ ((1 - Rdot./C).*R + ...
        JdotA*hH./C + 6*hH.*R.*ddintfnu*iDRe./C);

end

end
