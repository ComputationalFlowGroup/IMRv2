% file f_stress_calc.m
% brief contains function f_stress_calc

% brief This function features the stress integral and its time derivative
% solver. The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [S,Sdot,Z1dot,Z2dot] = f_stress_calc_graded(stress,Req,R,Ca,Ca1,Re8,...
    Rdot,alphax,intfnu,dintfnu,iDRe,l1,l2,v_a,v_nc)

Z1dot = [];
Z2dot = [];

% radial stretch
Rst = Req/R;
reltol = 1e-8;
abstol = 1e-8;
x1 = (1+(Rst.^-3-1)./(1+l1).^3).^(1/3);
x2 = (1+(Rst.^-3-1)./(1+l2).^3).^(1/3);
ycy = @(x) (1/Ca1+(1/Ca-1/Ca1)*(1+((x-x1)./(x2-x)).^v_a).^((v_nc-1)/v_a)).*(1./x.^5+1./x.^2);

% other models
% ype = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*asinh((x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ympe = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*log((x-x1(Rst))./(x2(Rst)-x)+1)./((x-x1(Rst))./(x2(Rst)-x)).^a).*(1./x.^5+1./x.^2);
% ycr = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).*(1./x.^5+1./x.^2);
% yscr = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*1./(1+(x-x1(Rst))./(x2(Rst)-x))).*(1./x.^5+1./x.^2);
% ymcr = @(x,Rst) (1/Ca1+(1/Ca-1/Ca1)*(1./(1+(x-x1(Rst))./(x2(Rst)-x)).^n).^a).*(1./x.^5+1./x.^2);

% no stress
if stress == 0
    S = 0;
    Sdot = 0;
    
    % Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 1
    Sv = - 4/Re8*Rdot/R - 6*intfnu*iDRe;
    
    Se1 = (1/(2*Ca))*(Rst.^4 + Rst - (1./x1.^4 + 1./x1));
    Se2 = 2*integral(@(x) ycy(x),x1,x2,'RelTol',reltol,'AbsTol',abstol);
    Se3 = -(1/(2*Ca1))*(5 - 4./x2 - 1./x2.^4);
    
    S = Sv + Se1 + Se2 + Se3;
    
    Sdotv = 4/Re8*(Rdot/R)^2 - 6*dintfnu*iDRe;
    Sdote = 0;% more on this later
    Sdot = Sdote + Sdotv;
    
    % quadratic Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 2
    S = (3*alphax-1)*(5 - Rst^4 - 4*Rst)/(2*Ca) - 4/Re8*Rdot/R - 6*intfnu*iDRe + ...
        (2*alphax/Ca)*(27/40 + (1/8)*Rst^8 + (1/5)*Rst^5 + ...
        Rst^2 - 2/Rst);
    Sdot = (Rdot/R)*((3*alphax - 1)/(2*Ca))*(4*Rst^4+4*Rst) + ...
        4*(Rdot/R)^2/Re8 - 6*dintfnu*iDRe -...
        2*alphax/Ca*Rdot/R*(Rst^8 + Rst^5 + 2*Rst^2 + 2*Rst^(-1));
else
    error('stress setting is not available');
end

end
