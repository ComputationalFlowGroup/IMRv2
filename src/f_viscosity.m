% file f_viscosity.m
% brief contains function f_viscosity

% brief This function computes the non-Newtonian viscosity integrals for
% the Carreau, Carreau-Yasuda, Powell-Eyring, modified Powell-Eyring,
% Cross, simplified Cross, and modified Cross models. The viscosity is
% mu = mu8 + (mu0 - mu8)*f(x) with x = lambda*|gammadot|, and each model
% supplies f(x) and x*f'(x). With gammadot = gammadot_R*(R/r)^3 and the
% substitution u = R/r, the stress integral and its time derivative become
% smooth integrals over u in [0, 1]:
%   I1 = int_R^Inf gammadot/r*f dr     = gammadot_R*int_0^1 u^2*f du
%   I2 = int_R^Inf (f + x*f')/r^4 dr   = R^-3*int_0^1 u^2*(f + x*f') du
% Outputs: f at the wall, intf (enters the stress integral as
% -6*intf/DRe), and the Rdot part (dintf) and Rddot coefficient (ddintf)
% of d(intf)/dt.
function [f,intf,dintf,ddintf] = f_viscosity(nu_model,Rdot,R,a,nc,lambda)
    
    gammadot_R   = -2*Rdot/R;
    % d(gammadot_R*R^3)/dt = dgammadot + ddgammadot*Rddot
    dgammadot    = -4*Rdot*Rdot*R;
    ddgammadot   = -2*R*R;
    abstol = 1E-8;
    reltol = 1E-8;
    
    % shear-rate parameter at the bubble wall
    xR = lambda*abs(gammadot_R);
    
    % viscosity at the wall
    f = f_model(xR);
    % Leibniz integration rule term at the lower limit r = R
    h = f*gammadot_R*(Rdot/R);
    % stress integral, goes into the stress term of the radial equation
    I1 = gammadot_R*integral(@(u) u.^2.*f_model(xR*u.^3), 0, 1, ...
        'RelTol',reltol,'AbsTol',abstol);
    % time derivative of the stress integral
    I2 = integral(@(u) u.^2.*f_sum(xR*u.^3), 0, 1, ...
        'RelTol',reltol,'AbsTol',abstol)/R^3;
    
    % S_def = (6/DRe)*I1 enters f_stress as -6*intf*iDRe  =>  intf = -I1
    intf = -I1;
    % d/dt intf = -(dg/dt)*I2 + Leibniz boundary term h; dg/dt splits
    % into the Rdot part (dintf) and the Rddot coefficient (ddintf)
    dintf = -dgammadot*I2 + h;
    ddintf = -ddgammadot*I2;
    
    % f(x) for the selected model
    function fx = f_model(x)
        [fx,~] = f_eval(x);
    end
    
    % f(x) + x*f'(x) for the selected model
    function s = f_sum(x)
        [fx,xdf] = f_eval(x);
        s = fx + xdf;
    end
    
    % non-Newtonian model evaluations, x = lambda*|gammadot| >= 0;
    % returns f(x) and x*df/dx
    function [fx,xdf] = f_eval(x)
        switch nu_model
            case 1
            % Carreau
            p = 1 + x.^2;
            fx = p.^((nc-1)/2);
            xdf = (nc-1)*x.^2.*p.^((nc-3)/2);
            case 2
            % Carreau-Yasuda
            p = 1 + x.^a;
            fx = p.^((nc-1)/a);
            xdf = (nc-1)*x.^a.*p.^((nc-1-a)/a);
            case 3
            % Powell-Eyring, asinh(x)/x -> 1 as x -> 0
            fx = ones(size(x));
            xdf = zeros(size(x));
            nz = x > 0;
            fx(nz) = asinh(x(nz))./x(nz);
            xdf(nz) = 1./sqrt(1 + x(nz).^2) - fx(nz);
            case 4
            % modified Powell-Eyring, log(1+x)/x^nc (nc <= 1)
            fx = zeros(size(x)) + (nc == 1);
            nz = x > 0;
            fx(nz) = log1p(x(nz))./x(nz).^nc;
            xdf = x.^(1-nc)./(1 + x) - nc*fx;
            case 5
            % Cross
            p = 1 + x.^nc;
            fx = 1./p;
            xdf = -nc*x.^nc./p.^2;
            case 6
            % simplified Cross
            fx = 1./(1 + x);
            xdf = -x./(1 + x).^2;
            case 7
            % modified Cross
            p = 1 + x.^nc;
            fx = p.^(-a);
            xdf = -a*nc*x.^nc.*p.^(-a-1);
            otherwise
            error('nu_model %d is not available', nu_model);
        end
    end
    
end
