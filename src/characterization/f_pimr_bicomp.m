%% bicomposite pimr estimate    
    % fitting parameters (intended to be the same value for all experiments)
    function [tc_pimr] = f_pimr_bicomp(G0,G1,l1,Pref,data_fit)
    % Unwrap input
    nX = length(data_fit(:,1));
    RX = data_fit(:,1);     % Not used, but let's keep for debugging
    LX = data_fit(:,2);
    f_Ma = data_fit(:,3);
    f_We = data_fit(:,4);
    f_gas = data_fit(:,5);
    trc = sqrt(pi/6)*gamma(5/6)/gamma(4/3); % ~= 0.9147

    reltol = 1e-8;
    abstol = 1e-8;
    Ca = G0/Pref;
    Ca1 = G1/Pref;
    tc_pimr = zeros(1,nX);
    for i = 1:nX
        Req_i = RX(i)./LX(i);
        %el1 = l1 ./ Req_i;
        el1 = l1;
        % let x be R
        Rst = @(x) x./Req_i;
        %Rst = @(x) x * LX(i); %./ (Req_i / RX(i));
        x1 = @(x) (1 + (Rst(x).^3 - 1)./(el1.^3)).^(1/3);

        %neo-Hookean elasticity
        Se0 = @(x) (1/(2*Ca)).*(Req_i.^4 ./ x.^4  +  4.*Req_i ./ x  -  1./x1(x).^4  -  4./x1(x));
        Se1 = @(x) -(1/(2*Ca1)).*(5 - 4./x1(x) - 1./x1(x).^4);
        Se = @(x) Se0(x) + Se1(x);

        % compute \bar{f}_gradedelastic^*
        feg_int = @(x) -1./((2/3)*(1./x.^3 - 1)).^0.5; % SINGULARITY AT X=1
        eps = 1e-7;
        f_eg = -1/trc * integral(@(x) Se(x).*feg_int(x),eps,1-eps,'RelTol',reltol,'AbsTol',abstol);
        % f_eg = 1/trc * integral(@(x) Se(x).*feg_int(x),1,0,'RelTol',reltol,'AbsTol',abstol);
        
        fsum = 1 - f_eg - f_We(i) - f_Ma(i) - f_gas(i);
        %find approximate collapse time
        if min(fsum) < 0
            t1_approx = nan;
        else
            t1_approx = (fsum).^(-1/2);
        end
        
        %store values
        tc_pimr(i) = t1_approx;
    end
end