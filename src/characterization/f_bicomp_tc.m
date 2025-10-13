%% modified collapse time for bicomposite 
%function [tg] = f_bicomp_tc(Req,R,R0,G0,G1,l1,Pref)
function [tg] = f_bicomp_tc(Req,R,G0,G1,l1,Pref,rho8)
    Ca = G0/Pref; Ca1 = G1/Pref;
    dtg_vals = zeros(1,length(R));

    for i = 1:length(R)
        %R = fliplr(R);
        Rnow = R(i); %normalized from 1 to Req/Rmax
        % radial stretch
        Rs = Req/Rnow; %R_0/R in math
        %Rm = R0/Rnow; %Rmax/R in math
        Rm = 1/Rnow;
        Rst = 1./Rs; %Lambda
        %Rmt = R0/Req; %Lambda_m
        Rmt = 1/Req;

        % volumes needed in integrals for energy density
        % vol_diff_now = Rnow^3 - Req^3;
        % if vol_diff_now < 0; stop; end
        % vol_diff_max = R0^3 - Req^3;
        % if vol_diff_max < 0; stop; end
        
        x1 = nthroot(1 + ((Rst.^3 -1)./l1^3),3); %Lambda_1
        xm1 = nthroot(1 + ((Rmt.^3 -1)/l1^3),3); %Lambda

        % nondim of parameters
        % usq = Pref/rho8;
        % tchar = (R0*Req)/sqrt(usq);
        
        Ee1 = (1/Ca).*(((1+ 2.*x1 +2.*x1.^2)./(x1 + x1.^2 + x1.^3)) - ((1+ 2*Rst +2*Rst.^2)./(Rst+ Rst.^2 + Rst.^3)));
        Eem1 = (1/Ca).*(((1+ 2.*xm1 +2.*xm1.^2)./(xm1 + xm1.^2 + xm1.^3)) - ((1+ 2.*Rmt +2.*Rmt.^2)./(Rmt+ Rmt.^2 + Rmt.^3)));
        Ee3 = (1/Ca1).*((5/3)- (1./x1) -((1+x1)./(1+ x1 +x1.^2)));
        Eem3 = (1/Ca1).*((5/3)- (1./xm1) -((1+xm1)./(1+ xm1 +xm1.^2)));
        Ee = Ee1 + Ee3;
        Eem = Eem1 + Eem3;


        %dtg_sq = 2/3.*(Pref/rho8).*(Rm.^3 -1) + 2*(Rm.^3 - Rs.^3).*Eem./rho8 - 2*(1 - Rs.^3).*Ee./rho8;
        %fprintf('i=%d, R =%.5e, dtg_sq=%.5e\n', i , Rnow, dtg_sq);
        term1 = 2/3.*(Rm.^3 -1); %.*(Pref/rho8)
        term2 = 2*(Rm.^3 - Rs.^3).*(Eem./rho8);%./vol_diff_max;
        term3 = 2*(1 - Rs.^3).*(Ee./rho8);%./vol_diff_now;
        terms = term1 + term2;
        dtg_sq = terms - term3;
        dtg_vals(i) = -1/sqrt(dtg_sq);
    end
    
    % remove initial R when Rnow = Rmax leading to Inf at dtg_vals
    % identify where this happens (in case not at initial R)
    inf_idx = isinf(dtg_vals);
    % remove from R and dtg_vals
    R_clean = R(~inf_idx);
    dtg_clean = dtg_vals(~inf_idx);
    % at min R (first collapse), dtg_val has small imag value
    dtg_cleaner = real(dtg_clean);
    % perform integration
    tg = trapz(R_clean,dtg_cleaner);

    % nondimensionalizing
    % tg = tg_dim / (R0/sqrt(Pref/rho8)
end