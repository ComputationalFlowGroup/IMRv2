function [tg] = f_bicomp_predict_tg(Req, Rmax, G0, G1, l1, Pref, rho8, PG0_k, kappa,pvsat)
    % Collapse time approximation using energy balance
    % Inputs: nondimensionalized Req and Rmax (should be ~1 and >1)
    
    % Material parameters (non-dimensionalized)
    Ca = Pref / G0;
    Ca1 = Pref / G1;

    % Radius vector (from Rmax to Req)
    R = linspace(Rmax, 0, 100000);
    
    dtg_vals = zeros(1, length(R));

    for i = 1:length(R)
        Rnow = R(i);
        
        Rs = Req / Rnow;       % R0 / R (R0 = Req)
        Rm = Rmax / Rnow;      % Rmax / R
        Rst = 1 / Rs;          % Lambda
        Rmt = Rmax / Req;      % Lambda_m

        % Lambda_1 and Lambda_m1
        %x1 = nthroot(1 + ((Rst^3 - 1) / l1^3), 3);
        %xm1 = nthroot(1 + ((Rmt^3 - 1) / l1^3), 3);
        x1 = (l1^3 + Rst^3 - 1).^(1/3);
        xm1 = (l1^3 + Rmt^3 - 1).^(1/3);

        % Elastic energy terms
        Ee1 = (1 / Ca)  * (((1 + 2 * x1 + 2 * x1^2) / (x1 + x1^2 + x1^3)) ...
                         - ((1 + 2 * Rst + 2 * Rst^2) / (Rst + Rst^2 + Rst^3)));
        Eem1 = (1 / Ca) * (((1 + 2 * xm1 + 2 * xm1^2) / (xm1 + xm1^2 + xm1^3)) ...
                         - ((1 + 2 * Rmt + 2 * Rmt^2) / (Rmt + Rmt^2 + Rmt^3)));

        Ee3 = (1 / Ca1) * ((5/3) - (1 / x1) - ((1 + x1) / (1 + x1 + x1^2)));
        Eem3 = (1 / Ca1) * ((5/3) - (1 / xm1) - ((1 + xm1) / (1 + xm1 + xm1^2)));

        Ee = Ee1 + Ee3;
        Eem = Eem1 + Eem3;

        % Energy balance terms (nondimensionalized)
        term1 = (2 / 3) * (Rm^3 - 1);
        term2 = 2 * (Rm^3 - Rs^3) * (Eem / rho8);
        term3 = 2 * (1 - Rs^3)   * (Ee  / rho8);

        pb = PG0_k.*Rs.^(3*kappa) + pvsat/Pref;
        E_bie = (pb.*(4/3).*pi.*(1./Rm).^3) ./ (kappa-1).*1./(2*pi).*Rm.^3;
        pb_m = PG0_k.*Rmt.^(-3*kappa) + pvsat/Pref;
        E_bie_m = pb_m.*(4/3).*pi ./ (kappa-1)./(2*pi);
        
        dtg_sq = term1 + term2 - term3 + E_bie_m - E_bie;

        if dtg_sq <= 0
            dtg_vals(i) = NaN;  % Skip unphysical
        else
            dtg_vals(i) = -1 / sqrt(dtg_sq);
        end
    end

    % Remove NaNs or Inf
    valid_idx = isfinite(dtg_vals);
    R_clean = R(valid_idx);
    dtg_clean = dtg_vals(valid_idx);

    % Integrate
    tg = trapz(R_clean, dtg_clean);
end

