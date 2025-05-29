% file f_stress_calc_graded.m
% brief contains function f_stress_calc_graded

% brief This function features the stress integral and its time derivative
% solver. The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [S,Sdot,Z1dot,Z2dot] = f_stress_calc_graded(stress,Req,R,Ca,Ca1,Re8,...
    Rdot,alphax,intfnu,dintfnu,iDRe,l1,l2,v_a,v_nc)

Z1dot = [];
Z2dot = [];

reltol = 1e-8;
abstol = 1e-8;

reltoldtycy = 1e-4;
abstoldtycy = 1e-7;

% radial stretch
Rst = R/Req;
Rstdot = Rdot/Req;

x1 = (1 + (Rst^3 - 1)/(l1^3))^(1/3);
x2 = (1 + (Rst^3 - 1)/(l2^3))^(1/3);
x1dot = Rstdot*Rst^2 / (l1*(x1^3 - 1))^(2/3);
x2dot = Rstdot*Rst^2 / (l2*(x1^3 - 1))^(2/3);

ycy = @(x) (1/Ca+(1/Ca-1/Ca1)*(1+( 1 - l1.*((x.^3 - 1)./(Rst.^3-1)).^(1/3)...
    ./ (l2.*((x.^3 - 1)./(Rst.^3-1)).^(1/3) - 1) ).^v_a).^((v_nc-1)/v_a))...
    .*(1./x.^5+1./x.^2);

dtycy = @(x) (1/Ca1 - 1/Ca).*(1./x.^5+1./x.^2).*(v_nc-1).* ...
    (1+( (Rst^3-1 - l1.*(x.^3 - 1).^(1/3)) ./ (l2.*(x.^3 - 1).^(1/3) - ...
    (Rst^3-1)^(1/3)) ).^v_a).^((v_nc-1-v_a)/v_a).*...
    ( (Rst^3-1 - l1.*(x.^3 - 1).^(1/3)) ./ (l2.*(x.^3 - 1).^(1/3) -...
    (Rst^3-1)^(1/3)) ).^(v_a - 1).*((Rstdot.*Rst.^2).*(l2-l1).*(x.^3 - 1).^(1/3)) ...
    ./ ( ((Rst.^3 - 1).^(2/3)).*(l2.*(x.^3 -1).^(1/3) - (Rst.^3 -1).^(1/3)).^2);

% no stress
if stress == 0
    S = 0;
    Sdot = 0;
    
    % Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 1
    Sv = - 4/Re8*Rdot/R - 6*intfnu*iDRe;
    
    Se1 = (1/(2*Ca))*(1/(Rst.^4) + 4/Rst - (1./x1.^4 + 4./x1));
    Se2 = 2*integral(@(x) ycy(x),x1,x2,'RelTol',reltol,'AbsTol',abstol);
    Se3 = -(1/(2*Ca1))*(5 - 4./x2 - 1./x2.^4);
    
    S = Sv + Se1 + Se2 + Se3;
    
    Svdot = 4/Re8*(Rdot/R)^2 - 6*dintfnu*iDRe;
    
    Se1dot = (2/Ca)*(x1dot./x1.^2 + x1dot./x1.^5 - Rstdot./Rst.^2 - Rstdot./Rst^.5);
    %Se2dot = 2*integral(@(x) dtycy(x),x1,x2,'RelTol',reltoldtycy,'AbsTol',abstoldtycy);
    Se2dot = 2*quadgk(dtycy,x1,x2,'RelTol',reltoldtycy,'AbsTol',abstoldtycy);
    % eps = 1e-6;
    % if x1 < 1 && x2 > 1
    %     Se2dot = 2*integral(@(x) dtycy(x),x1,1-eps,'RelTol',reltoldtycy,'AbsTol',abstoldtycy) + ...
    %              + 2*integral(@(x) dtycy(x),1+eps,x2,'RelTol',reltoldtycy,'AbsTol',abstoldtycy);
    % else
    %     Se2dot = 2*integral(@(x) dtycy(x),x1,x2,'RelTol',reltoldtycy,'AbsTol',abstoldtycy);
    % end
    Se3dot = -(2/Ca1)*(x2dot./x2.^5 + x2dot./x2.^2);
    Sedot = Se1dot + Se2dot + Se3dot;
    Sdot = Sedot + Svdot;
    
else
    error('stress setting is not available');
end

end

% other models
% ype = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*asinh((x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ympe = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*log((x-x1(Rst))./(x2(Rst)-x)+1)./((x-x1(Rst))./(x2(Rst)-x)).^a).*(1./x.^5+1./x.^2);
% ycr = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).*(1./x.^5+1./x.^2);
% yscr = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*1./(1+(x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ymcr = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*(1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).^a).*(1./x.^5+1./x.^2);
