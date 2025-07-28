% file f_tcol_calc_graded.m
% brief contains function f_tcol_calc_graded

% brief This function features the collapse time solver for a graded
% material. The solver currently assumes the Kelvin-Voigt with neo-Hookean
% elasticity.
function [tg] = f_tcol_calc_graded(stress,Req,R,R0,Ca,Ca1,Pref,l1,l2,v_a,v_nc,rho8)
    dtg_vals = zeros(1,length(R));
    Eem_vals = zeros(1,length(R));
    Ee_vals = zeros(1,length(R));
    for i = 1:length(R)
        Rnow = R(i);
        % radial stretch
        Rs = Req/Rnow; %R_0/R in math
        Rm = R0/Rnow; %Rmax/R in math
        Rst = 1./Rs; %Lambda
        Rmt = R0/Req; %Lambda_m
        
        x1 = (1 + ((Rst.^3 -1)./l1^3)).^(1/3); %Lambda_1
        x2 = (1 + ((Rst.^3 -1)./l2^3)).^(1/3); %Lambda_2
        xm1 = (1 + ((Rmt^3 -1)/l1^3))^(1/3); %Lambda_m1
        xm2 = (1 + ((Rmt^3 -1)/l2^3))^(1/3); %Lambda_m2

        % el1 = l1*R0;
        % el2 = l2*R0;
        % x1 = (1 + (Rnow^3 - Req^3)/el1^3)^(1/3); 
        % x2 = (1 + (Rnow^3 - Req^3)/el2^3)^(1/3);
        % xm1 = (1 + (R0^3 - Req^3)/el1^3)^(1/3);
        % xm2 = (1 + (R0^3 - Req^3)/el2^3)^(1/3);


        f_cy = @(x) ( l2.*(((x.^3 - 1)./(Rst.^3 - 1)).^(1/3)) - 1 ) ./ ( 1 - l1.*((x.^3 - 1)./(Rst.^3 - 1)).^(1/3));
        m = @(x) (1 + f_cy(x).^v_a).^((v_nc-1)/v_a);
        ee2integ = @(x) (1/Ca + (1/Ca1 - 1/Ca).* m(x)) .* ((1./x.^2) + 2.*x.^4 - 3.*x.^2)./((x.^3 - 1).^2);
        
        reltol = 1e-8;
        abstol = 1e-8;
        
        if stress == 1
            Ee1 = (1/Ca).*(((1+ 2.*x1 +2.*x1.^2)./(x1 + x1.^2 + x1.^3)) - ((1+ 2*Rst +2*Rst.^2)./(Rst+ Rst.^2 + Rst.^3)));
            Eem1 = (1/Ca)*(((1+ 2*xm1 +2*xm1^2)/(xm1 + xm1^2 + xm1^3)) - ((1+ 2*Rmt +2*Rmt^2)/(Rmt+ Rmt^2 + Rmt^3)));
            Ee2 = integral(@(x) ee2integ(x),x1,x2,'RelTol',reltol,'AbsTol',abstol); %integrating wrt lambda
            Eem2 = integral(@(x) ee2integ(x),xm1,xm2,'RelTol',reltol,'AbsTol',abstol); %integrating wrt lambda
            Ee3 = (1/Ca1)*((5/3)- (1./x2) -((1+x2)./(1+ x2 +x2.^2)));
            Eem3 = (1/Ca1)*((5/3)- (1/xm2) -((1+xm2)/(1+ xm2 +xm2^2)));
            Ee = Ee1 + Ee2 + Ee3;
            Eem = Eem1 + Eem2 + Eem3;

            Eem_vals(i) = Eem;
            Ee_vals(i) = Ee;
            
            %need rho in 1/Ca?
            %dtg = @(R) -1./sqrt( (2/3)*(Pref/rho)*(Rm.^3 -1) + 2*(Rm.^3 - Rs.^3).*Eem - 2*(1 - Rs.^3).*Ee );
            %dtg_vals = -1./sqrt( (2/3)*(Pref/rho8)*(Rm.^3 -1) + 2*(Rm.^3 - Rs.^3).*Eem - 2*(1 - Rs.^3).*Ee );
            %tg = integral(@(R) dtg(R),0,R0,'Reltol',reltol,'AbsTol',abstol);
            first = 2/3.*(Pref/rho8).*(Rm.^3 -1);
            second = 2*(Rm.^3 - Rs.^3).*Eem./rho8;
            third = 2*(1 - Rs.^3).*Ee./rho8;
            dtg_sq = first + second - third;
            % if isreal(dtg_sq) && dtg_sq > 0
            %     dtg_vals(i) = -1/sqrt(dtg_sq);
            % else
            %     dtg_vals(i) = NaN;
            % end
            if (first+second < third)
                warning('Ee is too strong')
            end
            if dtg_sq < 0
                warning('Imaginary collapse time at R = %.3e', Rnow);
                break;
            end
            %fprintf('i=%d, R =%.5e, dtg_sq=%.5e\n', i , Rnow, dtg_sq);
            
            % dtg_vals(i) = -1/sqrt(abs(dtg_sq));
            dtg_vals(i) = -1/sqrt(dtg_sq);
        end
    end
    % figure
    % hold on
    % plot(R, Eem_vals, 'k')
    % plot(R,Ee_vals,'b--')
   % tg = trapz(R,dtg_vals);

   % remove the initial time step where Rnow = Rmax leading to Inf at dtg_vals
   % identify where this happens
   inf_idx = isinf(dtg_vals);
   % remove from R and dtg_vals
   R_clean = R(~inf_idx);
   dtg_clean = dtg_vals(~inf_idx);
   % perform integration
   %tgc = trapz(R_clean,dtg_clean);

   % at min R (first collapse), dtg_val has small imag value
   dtg_cleaner = real(dtg_clean);
   % perform integration
   tg = trapz(R_clean,dtg_cleaner);
end
