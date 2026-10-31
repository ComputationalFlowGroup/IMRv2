% file material_stress.m
% brief contains function imr.physics.material_stress

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features the stress integral and its time derivative
% solver. The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B, 2nd-order Fung Kelvin-Voigt
% and Zener, and Mooney-Rivlin Kelvin-Voigt and Zener. R, Rdot, intfnu and
% dintfnu may be row vectors with X holding one state per column (a batch
% of states), every operation is elementwise across the columns
function [S,Sdot,Z1dot,Z2dot] = material_stress(stress,X,Req,R,Ca,De,Re8,...
    Rdot,alphax,ivisco1,ivisco2,LAM,intfnu,dintfnu,iDRe,C2)

Z1dot = [];
Z2dot = [];

% radial stretch
Rst = Req./R;

% no stress
if stress == 0
    S = 0;
    Sdot = 0;

    % Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 1
    S = -(5 - 4*Rst - Rst.^4)/(2*Ca) - 4/Re8*Rdot./R - 6*intfnu*iDRe;
    Sdot = -2*Rdot./R.*(Rst + Rst.^4)/Ca + 4/Re8*(Rdot./R).^2 - 6*dintfnu*iDRe;

    % quadratic Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 2
    S = (3*alphax-1)*(5 - Rst.^4 - 4*Rst)/(2*Ca) - 4/Re8*Rdot./R - 6*intfnu*iDRe + ...
        (2*alphax/Ca)*(27/40 + (1/8)*Rst.^8 + (1/5)*Rst.^5 + ...
        Rst.^2 - 2./Rst);
    Sdot = (Rdot./R)*((3*alphax - 1)/(2*Ca)).*(4*Rst.^4+4*Rst) + ...
        4*(Rdot./R).^2/Re8 - 6*dintfnu*iDRe -...
        2*alphax/Ca*Rdot./R.*(Rst.^8 + Rst.^5 + 2*Rst.^2 + 2*Rst.^(-1));

    % linear Maxwell, Jeffreys, Zener -- neo-Hookean
elseif stress == 3
    % extract stress auxiliary variable
    Z1 = X(ivisco1,:);
    S = Z1./R.^3 - 4*LAM/Re8*Rdot./R - 6*intfnu*iDRe;
    % elastic shift Ze
    Ze = -0.5*(R.^3/Ca).*(5 - Rst.^4 - 4*Rst);
    % stress auxiliary variable integral derivative
    Z1dot = -(Z1 - Ze)/De + 4*(LAM - 1)/(Re8*De)*R.^2.*Rdot;
    % stress integral derivative
    Sdot = Z1dot./R.^3 - 3*Rdot./R.^4.*Z1 + 4*LAM/Re8*(Rdot./R).^2 - 6*dintfnu*iDRe;

    % linear Maxwell, Jeffreys, Zener -- quadratic neo-Hookean
elseif stress == 4
    % extract stress auxiliary variable
    Z1 = X(ivisco1,:);
    S = Z1./R.^3 - 4*LAM/Re8*Rdot./R - 6*intfnu*iDRe;
    strainhard = (3*alphax - 1) / (2*Ca);
    % simplified Ze equation with decimal fractions
    Ze = R.^3 .* (strainhard * (5 - Rst.^4 - 4*Rst) + ...
        (2 * alphax / Ca) * (0.675 + 0.125 * Rst.^8 + ...
        0.2 * Rst.^5 + Rst.^2 - 2 ./ Rst));
    % stress auxiliary derivative
    Z1dot = -(Z1-Ze)/De + 4*(LAM-1)/(Re8*De)*R.^2.*Rdot ;
    % stress integral derivative
    Sdot = Z1dot./R.^3 - 3*Rdot./R.^4.*Z1 + 4*LAM/Re8*Rdot.^2./R.^2 - 6*dintfnu*iDRe;

    % upper-convected Maxwell, Oldroyd-B
elseif stress == 5
    % extract stress sub-integrals
    Z1 = X(ivisco1,:);
    Z2 = X(ivisco2,:);
    % compute new derivatives
    Z1dot = -(1/De - 2*Rdot./R).*Z1 + 2*(LAM-1)/(Re8*De)*R.^2.*Rdot;
    Z2dot = -(1/De +   Rdot./R).*Z2 + 2*(LAM-1)/(Re8*De)*R.^2.*Rdot;
    % stress integral and derivative
    S = (Z1 + Z2)./R.^3 - 4*LAM/Re8*Rdot./R - 6*intfnu*iDRe;
    Sdot = (Z1dot+Z2dot)./R.^3 - 3*Rdot./R.^4.*(Z1+Z2) + 4*LAM/Re8*Rdot.^2./R.^2 ...
        - 6*dintfnu*iDRe;

    % 2nd-order Fung Kelvin-Voigt elasticity, factor 1 + alphax*(I1-3) +
    % (alphax^2/2)*(I1-3)^2, the next term of the same series stress 2
    % truncates after its first order
elseif stress == 8
    S = (3*alphax-1)*(5 - Rst.^4 - 4*Rst)/(2*Ca) - 4/Re8*Rdot./R - 6*intfnu*iDRe + ...
        (2*alphax/Ca)*(27/40 + (1/8)*Rst.^8 + (1/5)*Rst.^5 + Rst.^2 - 2./Rst) + ...
        (alphax.^2/Ca).*(-2909/180 + 12./Rst - (4/3)./Rst.^3 - 6*Rst.^2 + 9*Rst + ...
        (4/3)*Rst.^3 + (9/4)*Rst.^4 - (6/5)*Rst.^5 + (2/3)*Rst.^6 - (3/4)*Rst.^8 + ...
        (1/9)*Rst.^9 + (1/12)*Rst.^12 + 4*log(Rst));
    Sdot = (Rdot./R).*((-2*Rst.^4 - 2*Rst)/Ca + ...
        alphax*(-2*Rst.^8 - 2*Rst.^5 + 6*Rst.^4 - 4*Rst.^2 + 6*Rst - 4./Rst)/Ca + ...
        alphax.^2.*(-Rst.^12 - Rst.^9 + 6*Rst.^8 - 4*Rst.^6 + 6*Rst.^5 - 9*Rst.^4 - ...
        4*Rst.^3 + 12*Rst.^2 - 9*Rst - 4 + 12./Rst - 4./Rst.^3)/Ca) + ...
        4/Re8*(Rdot./R).^2 - 6*dintfnu*iDRe;

    % linear Maxwell, Jeffreys, Zener -- 2nd-order Fung spring
elseif stress == 9
    Z1 = X(ivisco1,:);
    S = Z1./R.^3 - 4*LAM/Re8*Rdot./R - 6*intfnu*iDRe;
    % elastic shift Ze = R^3 times the stress-8 elastic-only S above
    Ze = R.^3.*((3*alphax-1)*(5 - Rst.^4 - 4*Rst)/(2*Ca) + ...
        (2*alphax/Ca)*(27/40 + (1/8)*Rst.^8 + (1/5)*Rst.^5 + Rst.^2 - 2./Rst) + ...
        (alphax.^2/Ca).*(-2909/180 + 12./Rst - (4/3)./Rst.^3 - 6*Rst.^2 + 9*Rst + ...
        (4/3)*Rst.^3 + (9/4)*Rst.^4 - (6/5)*Rst.^5 + (2/3)*Rst.^6 - (3/4)*Rst.^8 + ...
        (1/9)*Rst.^9 + (1/12)*Rst.^12 + 4*log(Rst)));
    Z1dot = -(Z1 - Ze)/De + 4*(LAM - 1)/(Re8*De)*R.^2.*Rdot;
    Sdot = Z1dot./R.^3 - 3*Rdot./R.^4.*Z1 + 4*LAM/Re8*(Rdot./R).^2 - 6*dintfnu*iDRe;

    % Mooney-Rivlin Kelvin-Voigt elasticity: the neo-Hookean term (scaled
    % by C1, absorbed into Ca as usual) plus one new term scaled by the
    % second modulus ratio C2 (0 reduces exactly to stress 1)
elseif stress == 10
    S = -(5 - 4*Rst - Rst.^4)/(2*Ca) - (C2/Ca)*(2./Rst - Rst.^2 - 1) - ...
        4/Re8*Rdot./R - 6*intfnu*iDRe;
    Sdot = (Rdot./R).*(-2*(Rst + Rst.^4)/Ca - (2*C2/Ca)*(1./Rst + Rst.^2)) + ...
        4/Re8*(Rdot./R).^2 - 6*dintfnu*iDRe;

    % linear Maxwell, Jeffreys, Zener -- Mooney-Rivlin spring
elseif stress == 11
    Z1 = X(ivisco1,:);
    S = Z1./R.^3 - 4*LAM/Re8*Rdot./R - 6*intfnu*iDRe;
    % elastic shift Ze = R^3 times the stress-10 elastic-only S above
    Ze = R.^3.*(-(5 - 4*Rst - Rst.^4)/(2*Ca) - (C2/Ca)*(2./Rst - Rst.^2 - 1));
    Z1dot = -(Z1 - Ze)/De + 4*(LAM - 1)/(Re8*De)*R.^2.*Rdot;
    Sdot = Z1dot./R.^3 - 3*Rdot./R.^4.*Z1 + 4*LAM/Re8*(Rdot./R).^2 - 6*dintfnu*iDRe;

else
    error('stress setting is not available');
end

end
