% file f_tcol_calc_graded.m
% brief contains function f_tcol_calc_graded

% brief This function features the collapse time solver for a graded
% material. The solver currently assumes the Kelvin-Voigt with neo-Hookean
% elasticity.
function [S,Sdot] = f_tcol_calc_graded(stress,Req,R,Ca,Ca1,Re8,...
    Rdot,alphax,intfnu,dintfnu,iDRe,l1,l2,v_a,v_nc)

addpath ../src;

% radial stretch
Rst = R/Req;
reltol = 1e-8;
abstol = 1e-8;
x1 = 1 + ((Rst^3 -1)/((1+l1)^3))^(1/3);
x2 = 1 + ((Rst^3 -1)/((1+l2)^3))^(1/3);
%m = (1 + f^v_a)^((v_nc-1)/v_a)
ee2integ = @(x) m * ((1./x.^2) - 2 - 3.*x.^2)./((x.^3 - 1).^2);

if stress == 1
    Ee1 = (1/Ca)*(((1+ 2*x1 +2*x1^2)/(x1 + x1^2 + x1^3)) - ((1+ 2*Rst +2*Rst^2)/(Rst+ Rst^2 + Rst^3)));
    Ee2 = integral(@(x) ee2integ(x),x2,x1,'RelTol',reltol,'AbsTol',abstol);
    Ee3 = (1/Ca1)*((5/3)- (1/x2) -((1+x2)/(1+ x2 +x2^2)));
    Ee = (Ee1 + Ee2 + Ee3).*(4*pi*(R^3 - Req^2));
    
    Eke = 2*pi*R^3*Rdot^2;
    Epe = (4/3)*pi*R^3*Pref;
    %how to address Pref and surften??
    Esf = 4*pi*R^2*surften;
    Etot = Ee + Eke + Epe + Esf;
    
    dtg = (2/3)*(Pref/rho)*((Rmax/R).^3 -1) + ((2*surften)/(rho*R))*((Rmax/R).^2 - 1) + (1/Ca)*(1 - (Req/Rmax).^3)*(Ee/(4*pi*(Rmax^3 - Req^3))) + (1/Ca)*(1-(1/Rst.^3))*(Ee/(4*pi*(Rmax^3 - Req^3)));
    %need rho in 1/Ca?
    
end

end
