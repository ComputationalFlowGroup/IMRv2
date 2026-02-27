% file f_f_gNH.m
% brief contains function f_fit_gNH

% brief This function estimates graded collapse time for a given set of graded input data
% describing a Neo-Hookean Newtonian material ongoing inertial cavitation.
function t1approx_vals = f_fit_gNH(G0_guess,G1_guess,l1_guess,l2_guess,va_guess,vnc_guess,data_fit)
    
    % Inputs:
    %   G_guess - base elastic shear modulus (Pa)
    %   G1_guess - upper elastic shear modulus (Pa)
    %   l1_guess - l1 wrt R0 (m)
    %   l2_guess - l1 wrt R0 (m)
    %   va_guess - Carreau-Yasuda parameter
    %   vnc_guess - Carreau-Yasuda parameter
    
    % Unwrap input
    nX = length(data_fit(:,1));
    RX = data_fit(:,1);
    LX = data_fit(:,2);
    f_Ma = data_fit(:,3);
    f_We = data_fit(:,4);
    f_gas = data_fit(:,5);
    Ca_scale = data_fit(:,6);
    
    % Other constants to use:
    trc = sqrt(pi/6)*gamma(5/6)/gamma(4/3); % ~= 0.9147
    
    % fitting parameters (intended to be the same value for all experiments)
    Ca_guess = Ca_scale(1)/G0_guess(1);
    Ca1_guess = Ca_scale(1)/G1_guess(1);
    
    reltol = 1e-1;
    abstol = 1e-1;
    t1approx_vals = zeros(1:nX);
    for i = 1:nX
        Req_i = RX(i)./LX(i);
        el1 = l1_guess ./ Req_i;
        el2 = l2_guess ./ Req_i;
        % let x be R
        Rst = @(x) x./Req_i;
        x1 = @(x) (1 + (Rst(x).^3 - 1)./(el1.^3)).^(1/3);
        x2 = @(x) (1 + (Rst(x).^3 - 1)./(el2.^3)).^(1/3);
        
        % let y be lambda
        fnum_cy = @(y,x) el2.*((y.^3 - 1)./(Rst(x).^3 - 1)).^(1/3) - 1;
        fden_cy = @(y,x) 1-el1.*((y.^3 - 1)./(Rst(x).^3 - 1)).^(1/3);
        f_cy = @(y,x) fnum_cy(y,x)./fden_cy(y,x);
        g_handle = @(y,x) (1/Ca_guess + (1/Ca1_guess - 1/Ca_guess).*(1+f_cy(y,x).^va_guess).^((vnc_guess-1)/va_guess)).*...
            ((1./y.^5) + (1./y.^2));
        
        %neo-Hookean elasticity
        Se1 = @(x) (1/(2*Ca_guess)).*(Req_i.^4 ./ x.^4  +  4.*Req_i ./ x  -  1./x1(x).^4  -  4./x1(x));
        Se2 = @(x) arrayfun(@(xi) 2*integral(@(y) g_handle(y,xi),x1(xi),x2(xi),'RelTol',reltol,'AbsTol',abstol),x);
        Se3 = @(x) -(1/(2*Ca1_guess)).*(5 - 4./x2(x) - 1./x2(x).^4);
        Se = @(x) Se1(x) + Se2(x) + Se3(x);
        
        % compute \bar{f}_gradedelastic^*
        feg_int = @(x) -1./((2/3)*(1./x.^3 - 1)).^0.5; % SINGULARITY AT X=1
        eps = 1e-7;
        f_eg = -1/trc * integral(@(x) Se(x).*feg_int(x),eps,1-eps,'RelTol',reltol,'AbsTol',abstol);
        % f_eg = 1/trc * integral(@(x) Se(x).*feg_int(x),1,0,'RelTol',reltol,'AbsTol',abstol);
        
        fsum = 1 - f_eg - f_We(i) - f_Ma(i) - f_gas(i);
        % find approximate collapse time
        if min(fsum) < 0
            t1_approx = nan;
        else
            t1_approx = (fsum).^(-1/2);
        end
        
        % store values
        t1approx_vals(i) = t1_approx;
    end
end
