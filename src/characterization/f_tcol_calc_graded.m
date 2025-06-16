% file f_tcol_calc_graded.m
% brief contains function f_tcol_calc_graded

% brief This function features the collapse time solver for a graded
% material. The solver currently assumes the Kelvin-Voigt with neo-Hookean
% elasticity.
function [tg] = f_tcol_calc_graded(stress,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc)

% radial stretch
Rs = Req/R; %R_0/R in math
Rm = R0/R; %Rmax/R in math
Rst = 1/Rs; %Lambda
Rmt = R0/Req; %Lambda_m

x1 = 1 + ((Rst^3 -1)/((1+l1)^3))^(1/3); %Lambda_1
x2 = 1 + ((Rst^3 -1)/((1+l2)^3))^(1/3); %Lambda_2
xm1 = 1 + ((Rmt^3 -1)/((1+l1)^3))^(1/3); %Lambda_m1
xm2 = 1 + ((Rmt^3 -1)/((1+l2)^3))^(1/3); %Lambda_m2

f_cy = @(x) ( l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1 ) ./ ( 1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) );
m = (1 + f_cy^v_a)^((v_nc-1)/v_a);
ee2integ = @(x) (1/Ca + (1/Ca1 - 1/Ca)* m) * ((1./x.^2) + 2 - 3.*x.^2)./((x.^3 - 1).^2);

reltol = 1e-8;
abstol = 1e-8;

if stress == 1
    Ee1 = (1/Ca)*(((1+ 2*x1 +2*x1^2)/(x1 + x1^2 + x1^3)) - ((1+ 2*Rst +2*Rst^2)/(Rst+ Rst^2 + Rst^3)));
    Eem1 = (1/Ca)*(((1+ 2*xm1 +2*xm1^2)/(xm1 + xm1^2 + xm1^3)) - ((1+ 2*Rmt +2*Rmt^2)/(Rmt+ Rmt^2 + Rmt^3)));
    Ee2 = integral(@(x) ee2integ(x),x2,x1,'RelTol',reltol,'AbsTol',abstol);
    Eem2 = integral(@(x) ee2integ(x),xm2,xm1,'RelTol',reltol,'AbsTol',abstol);
    Ee3 = (1/Ca1)*((5/3)- (1/x2) -((1+x2)/(1+ x2 +x2^2)));
    Eem3 = (1/Ca1)*((5/3)- (1/xm2) -((1+xm2)/(1+ xm2 +xm2^2)));
    Ee = Ee1 + Ee2 + Ee3;
    Eem = Eem1 + Eem2 + Eem3;
    
    %need rho in 1/Ca?
    dtg = -( (2/3)*(Pref/rho)*(Rm.^3 -1) + 2*(Rm.^3 - Rs.^3).*Eem - 2*(1 - Rs.^3).*Ee ).^0.5;
    tg = integral(@(x) dtg,0,R0,'RelTol',reltol,'AbsTol',abstol);
end

end
