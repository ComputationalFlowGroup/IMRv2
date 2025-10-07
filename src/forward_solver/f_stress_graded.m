% file f_stress_graded.m
% brief contains function f_stress_graded

% brief This function features the stress integral and its time derivative
% solver. The solver accounts for the Kelvin-Voigt with neo-Hookean
% elasticity, quadratic K-V neo-Hookean elasticity, linear Maxwell, linear
% Jeffreys, linear Zener, UCM and Oldroyd-B
function [S,Sdot,Z1dot,Z2dot] = f_stress_graded(radial,stress,Req,R,...
    Ca,Ca1,Re8,Rdot,alphax,intfnu,dintfnu,iDRe,l1,l2,v_a,v_nc,yT3)

Z1dot = [];
Z2dot = [];

reltol = 1e-4; %e-8
abstol = 1e-4; %e-8
% current ammplification factor
Rstinv = Req/R;
% incompressible condition (reference coord)
r0 = (yT3-1+Rstinv.^3).^(1/3);
chi = r0/Req;

% radial stretch
Rst = R/Req;
Rstdot = Rdot/Req;

x1 = (1 + (Rst.^3 - 1)./(l1^3))^(1/3);
x2 = (1 + (Rst.^3 - 1)./(l2^3))^(1/3);
x1dot = Rstdot*Rst.^2 ./ (l1*x1^2);
x2dot = Rstdot*Rst.^2 ./ (l2*x2^2);

% try
  % if mex file exist, compile it
    % if exist('f_g_mex','file') ~= 3
    %     mex('f_g_mex.c');
    % end
    % g_handle = @(x) f_g_mex(x,Ca,Ca1,Rst,l1,l2,v_a,v_nc,Rstdot,1); % mode 1 = g
    % gdot_handle = @(x) f_g_mex(x,Ca,Ca1,Rst,l1,l2,v_a,v_nc,Rstdot,2); % mode 2 = gdot
% catch
    % fnum_cy = @(x) l2*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3) - 1;
    % fden_cy = @(x) 1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3);
    fnum_cy = @(x) l2*(x.^3 -1).^(1/3) - (Rst.^3 - 1).^(1/3);
    fden_cy = @(x) (Rst.^3 - 1).^(1/3) - l1*(x.^3 -1).^(1/3);
    f_cy = @(x) fnum_cy(x)./fden_cy(x);

   % fdot_val = fdot_cy(x,l1,l2,Rst,Rstdot);
   fdot_cy = @(x) (((x.^3-1).^(1/3))./(Rst.^3-1)^(4/3)).*Rstdot.*(Rst.^2).*(l1-l2)./(fden_cy(x).^2);
    % dtf_1 = @(x) 1 ./ ((Rst.^3 - 1).^(1/3) .* (x.^3 -1).^(2/3));
    % dtf_2 = @(x) ((x.^3 -1).^(1/3)) ./ ((Rst.^3 - 1).^(4/3) .* chi.^3);
    % fdot_cy_scalar = @(x) (l2 - l1)./(1- l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) .* (Rst.^2 .* Rstdot) .* (dtf_1(x) - dtf_2(x));
    % fdot_cy = @(x) arrayfun(fdot_cy_scalar, x);

    g_handle = @(x) (1/Ca + (1/Ca1 - 1/Ca)*(1+f_cy(x).^v_a).^((v_nc-1)/v_a)).*((1./x.^5) + (1./x.^2));
    gdot_handle = @(x) (1/Ca1 - 1/Ca).*((1./x.^5) + (1./x.^2)).*(v_nc-1).*...
        ((1+f_cy(x).^v_a).^((v_nc-1-v_a)/v_a)).*f_cy(x).^(v_a-1).*fdot_cy(x);

    f_tanh = @(x) ((Rst.^3 - 1)/(x.^3 - 1)).^(1/3);
    g_tanh = @(x) (1/Ca + (1/Ca1 - 1/Ca)*(1/2)*(1+tanh(v_a*(2*f_tanh(x,Rst)-1)))).*(1./x.^5+1./x.^2);
    %fdot_tanh = @(x) Rst.^2 .* Rstdot .* (dtf_1(x) - dtf_2(x));
    %gdot_tanh = @(x) (1/Ca1 - 1/Ca).* .5 .* sech(a*(2*(f_tanh(x) -1)))^2; % did not account for ((1./x.^5) + (1./x.^2))
    gdot_tanh = @(x) (1/Ca + (1/Ca1 - 1/Ca)*(1/2)*(1+tanh(v_a*(2*f_tanh(x,Rst)-1)))).*(-5./x.^4 - 2./x).*(Rst.^2 .* Rstdot .* 1./x.^2 .* chi.^3);
%end

% no stress
if stress == 0
    S = 0;
    Sdot = 0*alphax;
    
    % Kelvin-Voigt with neo-Hookean elasticity
elseif stress == 1
    Sv = - 4/Re8*Rdot/R - 6*intfnu*iDRe;
    
    Se1 = (1/(2*Ca))*(1/(Rst.^4) + 4/Rst - (1./x1.^4 + 4./x1));
    if Rst == 1 % fix OR do not start at R=Req
        Se2 = 0;
    else
    % fnum_cy(x1)
    % fnum_cy(x2)
    % fden_cy(x2)
    % fden_cy(x1)
    % g_handle(x1)
    % g_handle(x2)
         Se2 = 2*integral(@(x) g_handle(x),x1,x2,'RelTol',reltol,'AbsTol',abstol);
    end
    Se3 = -(1/(2*Ca1))*(5 - 4./x2 - 1./x2.^4);
    
    S = Sv + Se1 + Se2 + Se3;
    Sdot = 0;
    
    if radial > 1
        Svdot = 4/Re8*(Rdot/R)^2 - 6*dintfnu*iDRe;
        Se1dot = (2/Ca)*(x1dot./x1.^2 + x1dot./x1.^5 - ...
            Rstdot./Rst.^2 - Rstdot./Rst^.5);
        % bu_cond = @(Rst) (((1/l1)^3 .* (Rst.^3 -1)) + 1).^(1/3); % is dynamically changing
        % if (Rst == 1) || (x == bu_cond)
        %     Se2dot = 0;
        % elseif (Rstdot == 0) && (Rst == 1) || (x == bu_cond) 
        %     Se2dot = 0;
        % else
            Se2dot =  2*(g_handle(x2)*x2dot - g_handle(x1)*x1dot + ...
                integral(@(x) gdot_handle(x),x1,x2,'RelTol',reltol,'AbsTol',abstol));
                % (2/Ca1)*(x2dot./x2.^5 + x2dot./x2.^2) -...
                % (2/Ca)*(x1dot./x1.^5 + x1dot./x1.^2) + ...
                % 2*(g_handle(x2)*x2dot - g_handle(x1)*x1dot + ...
                % integral(@(x) gdot_handle(x),x1,x2,'RelTol',reltol,'AbsTol',abstol));
        %end
        Se3dot = -(2/Ca1)*(x2dot./x2.^5 + x2dot./x2.^2);
        Sedot = Se1dot + Se2dot + Se3dot;
        Sdot = Sedot + Svdot;
    end
    
else
    error('stress setting is not available');
end
% 
% function out = fdot_cy(x,l1,l2,Rst,Rstdot)
%     base = ((x.^3 - 1).^(1/3)) ./ ((Rst.^3 - 1).^(4/3));
%     den =  1-l1*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3);
%     if isnan(base / den^2)
%         out = Rstdot * Rst^2 * (l1-l2);
%     else
%         out = base * Rstdot * Rst^2 * (l1-l2) / den^2;
%     end
% end

end
