% file f_tg_calc.m
% brief contains function f_tg_calc

% brief This function features the collapse time solver for a graded
% material. The solver currently assumes the Kelvin-Voigt with neo-Hookean
% elasticity.
function [tg] = f_tg_calc(stress,Req,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8)
% let x be R
% radial stretch
Rs = @(x) Req./x; %R_0/R in math
Rm = @(x) R0./x; %Rmax/R in math
Rst = @(x) 1./Rs(x); %Lambda
Rmt = R0/Req; %Lambda_m

x1 = @(x) (1 + ((Rst(x).^3 -1)./l1^3)).^(1/3); %Lambda_1
x2 = @(x) (1 + ((Rst(x).^3 -1)./l2^3)).^(1/3); %Lambda_2
xm1 = (1 + ((Rmt^3 -1)/l1^3))^(1/3); %Lambda_m1
xm2 = (1 + ((Rmt^3 -1)/l2^3))^(1/3); %Lambda_m2

% let y be lambda
f_cy = @(y,x) ( l2.*(((y.^3 - 1)./(Rst(x).^3 - 1)).^(1/3)) - 1 ) ./ ( 1 - l1.*((y.^3 - 1)./(Rst(x).^3 - 1)).^(1/3));
m = @(y,x) (1 + f_cy(y,x).^v_a).^((v_nc-1)/v_a);
ee2integ = @(y,x) (1/Ca + (1/Ca1 - 1/Ca).* m(y,x)) .* ((1./y.^2) + 2.*y.^4 - 3.*y.^2)./((y.^3 - 1).^2);

reltol = 1e-8;
abstol = 1e-8;

if stress == 1
    Ee1 = @(x) (1/Ca).*(((1+ 2.*x1(x) +2.*x1(x).^2)./(x1(x) + x1(x).^2 + x1(x).^3)) - ((1+ 2*Rst(x) +2*Rst(x).^2)./(Rst(x)+ Rst(x).^2 + Rst(x).^3)));
    Eem1 = (1/Ca)*(((1+ 2*xm1 +2*xm1^2)/(xm1 + xm1^2 + xm1^3)) - ((1+ 2*Rmt +2*Rmt^2)/(Rmt+ Rmt^2 + Rmt^3)));
    Ee2 = @(x) arrayfun(@(xi) integral(@(y) ee2integ(y,xi),x1(xi),x2(xi),'RelTol',reltol,'AbsTol',abstol),x); %integrating wrt lambda
    Eem2 = @(x) arrayfun(@(xi) integral(@(y) ee2integ(y,xi),xm1,xm2,'RelTol',reltol,'AbsTol',abstol),x); %integrating wrt lambda
    Ee3 = @(x) (1/Ca1)*((5/3)- (1./x2(x)) -((1+x2(x))./(1+ x2(x) +x2(x).^2)));
    Eem3 = (1/Ca1)*((5/3)- (1/xm2) -((1+xm2)/(1+ xm2 +xm2^2)));
    Ee = @(x) Ee1(x) + Ee2(x) + Ee3(x);
    Eem = @(x) Eem1 + Eem2(x) + Eem3;

    dtg_sq = @(x) (2/3)*(Pref/rho8).*(Rm(x).^3 -1) + 2*(Rm(x).^3 - Rs(x).^3).*Eem(x)./rho8 - 2*(1 - Rs(x).^3).*Ee(x)./rho8;
    %fprintf('i=%d, R =%.5e, dtg_sq=%.5e\n', i , Rnow, dtg_sq);
    dtg_int = @(x) 1./sqrt(dtg_sq(x));
    %eps = 1e-7;
    %tg = integral(@(x) dtg_int(x),eps,1-eps,'RelTol',reltol,'AbsTol',abstol);
    tg = -integral(@(x) dtg_int(x),0,1,'RelTol',reltol,'AbsTol',abstol);
end


% % remove initial R when Rnow = Rmax leading to Inf at dtg_vals
% % identify where this happens (in case not at initial R)
% inf_idx = isinf(dtg_vals);
% % remove from R and dtg_vals
% R_clean = R(~inf_idx);
% dtg_clean = dtg_vals(~inf_idx);
% % at min R (first collapse), dtg_val has small imag value
% dtg_cleaner = real(dtg_clean);
% perform integration
%tg = trapz(R_clean,dtg_cleaner);
end