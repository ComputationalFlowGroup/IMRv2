% file mie_gruneisen_eos.m
% brief contains function imr.physics.mie_gruneisen_eos

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief Dimensionless Mie-Gruneisen equation of state used by the Gilmore
% (radial 6) and far-field-sound-speed Keller-Miksis enthalpy form
% (radial 5), given a dimensionless pressure (scaled by P8) and the three
% parameters any profile needs to use it: the far-field sound speed Cstar
% (c8/Uc), the Hugoniot slope hugoniot_s (the S of the linear shock
% relation Us = C0 + S*up), and nog = (nstate - 1)/2, where nstate relates
% to the material's Gruneisen parameter. These three are plain solver
% inputs (schema.m), independent of any fluid profile: a profile supplies
% the ambient state (rho8, c8) that fixes Cstar, never hugoniot_s or
% nstate themselves (docs/model/radial-dynamics.md, docs/model/fluids.md).
% This file is the general p(rho,e)-style evaluator behind that split: it
% was previously private to imr.physics.radial_acceleration, pulled out
% here so it can be tested and reused on its own.
%
% At p = 1 (the ambient reference state) mu = 0, hB = 0 and C = Cstar by
% construction, for any hugoniot_s and nog: coupling a new profile's
% ambient (rho8, c8) to this equation of state is always thermodynamically
% consistent at the reference point. hugoniot_s away from 0 is the part
% that needs a material-specific fit; hugoniot_s = 0 is the defensible
% fallback (a purely linear/acoustic equation of state, no shock
% stiffening) when no such fit exists.
%
% Inputs:
%   pval        dimensionless pressure (scalar or array)
%   Cstar       dimensionless far-field sound speed, c8/Uc
%   hugoniot_s  Hugoniot slope S (0 for a purely linear/acoustic EoS)
%   nog         (nstate - 1)/2
% Outputs:
%   C   local sound speed at pval
%   hB  enthalpy h(pval) = integral from 1 to pval of dp/rho, closed form
%   hH  1/rho(pval), so that dot{h} = hH*dot{p}
function [C, hB, hH] = mie_gruneisen_eos(pval, Cstar, hugoniot_s, nog)

    mu = f_mie_mu(pval);
    C = Cstar*sqrt(((1+2*nog*mu).*(1-hugoniot_s*mu).^2+ ...
        2*hugoniot_s*mu.*(1+nog*mu).*(1-hugoniot_s*mu))./(1-hugoniot_s*mu).^4);

    %   hH  - 1 / ρ(P) (for dot{h} = hH * dot{P})
    hH = 1./(1 + mu);

    % enthalpy h(P) = ∫_1^P dp/ρ in closed form. With p = Cstar^2*g(mu),
    % g = mu*(1 + nog*mu)/(1 - s*mu)^2, integration by parts gives
    % h = [p/(1 + mu)] + Cstar^2*∫ g/(1 + mu)^2 dmu, and the partial
    % fractions of g/(1 + mu)^2 have poles at mu = -1 and mu = 1/s. The
    % differences are written so that they hold for s = 0 and near P = 1.
    s = hugoniot_s;
    mu1 = f_mie_mu(1);
    dmu = mu - mu1;
    A1 = -(2*nog + s - 1)/(1 + s)^3;
    A2 = (nog - 1)/(1 + s)^2;
    B2 = (s + nog)/(1 + s)^2;
    G = A1*(log1p(dmu/(1 + mu1)) - log1p(-s*dmu/(1 - s*mu1))) + ...
        A2*dmu./((1 + mu)*(1 + mu1)) + B2*dmu./((1 - s*mu)*(1 - s*mu1));
    hB = pval./(1 + mu) - 1/(1 + mu1) + Cstar^2*G;

    % density strain mu = rho - 1 of the Mie-Gruneisen EoS at pressure p
    function mu = f_mie_mu(p)
        A = p/Cstar^2;
        As2 = A*hugoniot_s^2;
        a = As2 - nog;
        b = -2*A*hugoniot_s - 1;
        d = b.^2 - 4*a.*A;
        % small (physical) root, in numerically stable form: mu -> A as
        % A -> 0, i.e. drho/rho ~ p/(rho*c^2). The other root is the
        % spurious -1/nog branch since a < 0 near ambient
        mu = 2*A ./ (-b + sqrt(d));
    end

end
