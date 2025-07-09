% file f_stress_dissipation_graded.m
% brief contains function f_stress_dissipation_graded

% brief This function computes the stress dissipation term for a graded material, 
% \tau : \nabla u. The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [taugradu] = f_stress_dissipation_graded(stress,spectral,Req,R,Rdot,Ca,Ca1, ...
    Br,Re8,yT,yT3,iyT3,iyT6,X,ZZT,ivisco1,ivisco2,fnu,DRe,l1,l2,v_a,v_nc)

% current ammplification factor
Rst = Req/R;
% incompressible condition
x2 = (yT3-1+Rst.^3).^(2/3);
% reference coordinate
x = x2.^0.5;

% no stress
if stress == 0
    taugradu = 0;
    
    % Kelvin-Voigt neo-Hookean solid
elseif stress == 1
    taugradu_viscous = 12*(Br/(Re8+DRe*fnu))*(Rdot/R)^2*iyT6;
    %taugradu_elastic = 2*Br/Ca*iyT3.*(Rdot/R).*(yT2.*ix2 - iyT4.*x4);
    taugradu_elastic = 2*Br*iyT3.*(Rdot/R);
    
    % function for graded region
    f_cy = (l2 - x) ./ (x - l1);
    % near field
    r0near = x;
    % zero out regions larger than the near field
    r0near(r0near > l1) = 0;
    taurr1 = (2/(3*Ca))*((r0near.^4 ./ yT.^4) - (yT.^2 ./ r0near.^2));
    % zero out the stress above and below the near field
    taurr1(isinf(taurr1) | isnan(taurr1)) = 0;
    % far field
    r0far = x;
    % zero out regions less than the far field
    r0far(r0far < l2) = 0;
    taurr3 = (2/(3*Ca1))*((r0far.^4 ./ yT.^4) - (yT.^2 ./ r0far.^2));
    % zero out the stress below the far field
    taurr3(isinf(taurr3) | isnan(taurr3)) = 0;
    % compute the graded region coordinate
    r0mid = x - r0far - r0near;
    % graded stress
    taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).* ...
    ((r0mid.^4 ./ yT.^4) - (yT.^2 ./ r0mid.^2));
    % removing the negative infinities
    taurr2(isinf(taurr2) | isnan(taurr2)) = 0;

    % debugging
    if any(imag(taurr1)) || any(isnan(taurr1)) || any(isinf(taurr1))
        error('taurr1 imaginary')
    elseif any(imag(taurr2)) || any(isnan(taurr2)) || any(isinf(taurr2))
        error('taurr2 imaginary')
    elseif any(imag(taurr3)) || any(isnan(taurr3)) || any(isinf(taurr3))
        error('taurr3 imaginary')
    end
    
    taugradu = taugradu_viscous + taugradu_elastic.*(taurr1+taurr2+taurr3);
    
else
    error('stress setting is not available'); 
end

% spectral solution approach to compute the mechanical dissipation
if spectral == 1
    taugradu = -2*Br*Rdot./(R*yT3).* ...
        (ZZT*(X(ivisco1)-X(ivisco2)));
end

end
