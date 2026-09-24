% file f_radial_eq.m
% brief contains function f_radial_eq

% brief This function features the spherical radial bubble dynamics
% equations used to compute the bubble wall acceleration. The radial models
% Rayleigh-Plesset, Keller-Miksis in pressure and enthalpy, and Gilmore are
% available in the solver.
function [Rddot] = f_radial_eq(radial, P, Pdot, Pf8, Pf8dot, iWe, ...
    R, Rdot, S, Sdot, Cstar, sam, no, GAMa, nstate, nog, hugoniot_s, JdotA, ...
    ddintfnu, iDRe )

% ddintfnu is the Rddot coefficient of d/dt(intfnu); the stress enters as
% -6*intfnu*iDRe, so its Rddot part moves to the denominator as
% +6*R*ddintfnu*iDRe/C, mirroring JdotA

% Rayleigh-Plesset
if radial == 1
    Rddot = (P - 1 - Pf8 - iWe/R + S - 1.5*Rdot^2)/R;
    
    % Keller-Miksis in pressure
elseif radial == 2
    Rddot = ((1+Rdot./Cstar)*(P - 1 - Pf8 - iWe/R + S) ...
        + R/Cstar*(Pdot + iWe*Rdot/R^2 +  Sdot - Pf8dot) ...
        - 1.5*(1-Rdot/(3*Cstar))*Rdot^2)/((1-Rdot/Cstar)*R + JdotA/Cstar ...
        + 6*R*ddintfnu*iDRe/Cstar);
    
    % Keller-Miksis in enthalpy with Tait EoS
elseif radial == 3
    Pb = P - iWe/R + GAMa + S;
    hB = (sam/no)*((Pb/sam)^no - 1);
    hH = (sam/Pb)^(1/nstate);
    Rddot = ((1 + Rdot/Cstar)*(hB - Pf8) - R/Cstar*Pf8dot ...
        + R/Cstar*hH*(Pdot + iWe*Rdot/R^2 + Sdot) ...
        - 1.5*(1 - Rdot/(3*Cstar))*Rdot^2) / ((1 - Rdot/Cstar)*R + ...
        JdotA*hH/Cstar + 6*hH*R*ddintfnu*iDRe/Cstar);
    
    % Gilmore equation with Tait EoS
elseif radial == 4
    Pb = P - iWe/R + GAMa + S;
    rho = (Pb/sam)^(1/nstate);
    C = sqrt(nstate*Pb/rho);
    hB = sam/no*((Pb/sam)^no - 1);
    hH = (sam/Pb)^(1/nstate);
    Rddot = ((1 + Rdot/C)*(hB - Pf8) - R/C*Pf8dot ...
        + R/C*hH*(Pdot + iWe*Rdot/R^2 + Sdot) ...
        - 1.5*(1 - Rdot/(3*C))*Rdot^2) / ((1 - Rdot/C)*R + ...
        JdotA*hH/C + 6*hH*R*ddintfnu*iDRe/C);
    
    % Keller-Miksis in enthalpy with Mie-Gruneisen EoS
elseif radial == 5
    Pb = P - iWe/R + S;
    [~, hB, hH] = f_mie_gruneisen_eos_scalar(Pb);
    Rddot = ((1 + Rdot/Cstar)*(hB - Pf8) - R/Cstar*Pf8dot ...
        + R/Cstar*hH*(Pdot + iWe*Rdot/R^2 + Sdot) ...
        - 1.5*(1 - Rdot/(3*Cstar))*Rdot^2) / ((1 - Rdot/Cstar)*R + ...
        JdotA*hH/Cstar + 6*hH*R*ddintfnu*iDRe/Cstar);
    
    % Gilmore with Mie-Gruneisen EoS
elseif radial == 6
    Pb = P - iWe/R + S;
    [C, hB, hH] = f_mie_gruneisen_eos_scalar(Pb);
    Rddot = ((1 + Rdot/C)*(hB - Pf8) - R/C*Pf8dot ...
        + R/C*hH*(Pdot + iWe*Rdot/R^2 + Sdot) ...
        - 1.5*(1 - Rdot/(3*C))*Rdot^2) / ((1 - Rdot/C)*R + ...
        JdotA*hH/C + 6*hH*R*ddintfnu*iDRe/C);
    
end

% computes mu, rho, enthalpy h, and hH = 1/rho at a scalar pressure P
% using the Mie–Grüneisen EoS
function [C, hB, hH] = f_mie_gruneisen_eos_scalar(pval)
    
    % normalize pressure
    A = pval / Cstar^2;
    As = A*hugoniot_s;
    As2 = A*hugoniot_s^2;
    
    % quadratic coefficients for solving mu
    a = As2 - nog;
    b = -2*As - 1;
    % discriminant calculation
    d = b^2 - 4*a*A;
    
    % density strain, small (physical) root in numerically stable form:
    % mu -> A as A -> 0, i.e. drho/rho ~ p/(rho*c^2). The other root is the
    % spurious -1/nog branch since a < 0 near ambient
    mu = 2*A / (-b + sqrt(d));
    C = Cstar*sqrt(((1+2*nog*mu)*(1-hugoniot_s*mu)^2+ ...
        2*hugoniot_s*mu*(1+nog*mu)*(1-hugoniot_s*mu))/(1-hugoniot_s*mu)^4);
    
    %   hH  - 1 / ρ(P) (for dot{h} = hH * dot{P})
    hH = 1/(1 + mu);
    
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
        A2*dmu/((1 + mu)*(1 + mu1)) + B2*dmu/((1 - s*mu)*(1 - s*mu1));
    hB = pval/(1 + mu) - 1/(1 + mu1) + Cstar^2*G;
    
end

% density strain mu = rho - 1 of the Mie-Gruneisen EoS at pressure p
function mu = f_mie_mu(p)
    A = p/Cstar^2;
    a = A*hugoniot_s^2 - nog;
    b = -2*A*hugoniot_s - 1;
    d = b^2 - 4*a*A;
    mu = 2*A/(-b + sqrt(d));
end

end
