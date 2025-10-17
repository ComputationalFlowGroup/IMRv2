% file f_graded.m
% brief contains function f_gstressint used in ../src/forward_solver/f_stress_calc_graded
% brief contains function f_gstress used in ../examples/

% brief This function calls the stress integral and its time derivative
% function, also located here, in case the user wants to obtain all
% quantities related to graded materials.
function [Se,Sedot,taurr1,taurr,r1,r2] = f_graded(stress,Req,R,Ca,Ca1,Rdot,l1,l2,v_a,v_nc)
% addpath('src/forward_solver/')
if stress == 1
   [Se,Sedot] = f_gstressint(stress,Req,R,Ca,Ca1,Rdot,l1,l2,v_a,v_nc);
   %actually wouldn't want to couple since need r_coord for taus and NOT in S calc
   [taurr1,taurr,r1,r2] = f_gradedstress(r_coord,R,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a);
end

end

% brief This function features the stress integral and its time derivative
% solver. The solver accounts for Kelvin-Voigt with neo-Hookean elasticity.
function [Se,Sedot] = f_gstressint(stress,Req,R,Ca,Ca1,Rdot,l1,l2,v_a,v_nc)
%addpath('../src/forward_solver/')

reltol = 1e-8;
abstol = 1e-8;

reltoldtycy = 1e-4;
abstoldtycy = 1e-7;

% radial stretch
Rst = R/Req;
Rstdot = Rdot/Req;

x1 = (1 + (Rst.^3 - 1)./(l1^3))^(1/3);
x2 = (1 + (Rst.^3 - 1)./(l2^3))^(1/3);
x1dot = Rstdot*Rst.^2 ./ (l1*x1^2);
x2dot = Rstdot*Rst.^2 ./ (l2*x2^2);

%f_cy = @(x) (l2*((x.^3 - 1)/(Rst^3 - 1)).^(1/3) - 1)/(1-l1*((x.^3 - 1)/(Rst^3 - 1)).^(1/3));
fnum_cy = @(x) l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1;
fden_cy = @(x) 1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3);
f_cy = @(x) fnum_cy(x)./fden_cy(x);
fdot_cy = @(x) (((x.^3-1).^(1/3))./(Rst.^3-1)^(4/3)).*Rstdot.*(Rst.^2).*(l1-l2)./(fden_cy(x).^2);

g = @(x) (1/Ca + (1/Ca1 - 1/Ca)*(1+f_cy(x).^v_a).^((v_nc-1)/v_a)).*((1./x.^5) + (1./x.^2));
gdot = @(x) (1/Ca1 - 1/Ca).*((1./x.^5) + (1./x.^2)).*(v_nc-1).*...
    ((1+f_cy(x).^v_a).^((v_nc-1-v_a)/v_a)).*f_cy(x).^(v_a-1).*fdot_cy(x);

if stress == 1
    Se1 = (1/(2*Ca))*(1/(Rst.^4) + 4/Rst - (1./x1.^4 + 4./x1));
    Se2 = 2*integral(@(x) g(x),x1,x2,'RelTol',reltol,'AbsTol',abstol);
    Se3 = -(1/(2*Ca1))*(5 - 4./x2 - 1./x2.^4);
    Se = Se1 + Se2 + Se3;

    Se1dot = (2/Ca)*(x1dot./x1.^2 + x1dot./x1.^5 - Rstdot./Rst.^2 - Rstdot./Rst^.5);
    Se2dot = 2*(g(x2)*x2dot - g(x1)*x1dot + integral(@(x) gdot(x),x1,x2,'RelTol',reltoldtycy,'AbsTol',abstoldtycy));
    Se3dot = -(2/Ca1)*(x2dot./x2.^5 + x2dot./x2.^2);
    Sedot = Se1dot + Se2dot + Se3dot;
end
end

function [taurr1,taurr,r1,r2] = f_gradedstress(r_coord,R,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a)
   aa = r_coord.^3 - R.^3 + Req^3;
   aa = (1./(1-(aa<0))).*aa;
   r0_coord = real((aa).^(1/3));
 % r0_coord = (r_coord.^3 - R.^3 + Req^3).^(1/3);
   % valid = r0_coord > 0 & isreal(r0_coord); %where r0_coord values are negative (inside bubble)
   % r0coord = r0_coord.*valid;
    
   f_cy = (l2 - r0_coord) ./ (r0_coord - l1);
        
   taurr1 = (2/(3*Ca))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr3 = (2/(3*Ca1))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));

   taurr = taurr1 + taurr2 + taurr3;

   r1 = ((l1/R0)^3 + R.^3 - Req^3).^(1/3);
   r2 = ((l2/R0)^3 + R.^3 - Req^3).^(1/3);
end