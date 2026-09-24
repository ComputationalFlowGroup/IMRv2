% file f_viscosity.m
% brief contains function f_viscosity

% brief This function computes the non-Newtonian viscosity integrals for
% the Carreau, Carreau-Yasuda, Powell-Eyring, modified Powell-Eyring,
% Cross, simplified Cross, and modified Cross models. The viscosity is
% mu = mu8 + (mu0 - mu8)*f(x) with x = lambda*|gammadot|, and each model
% supplies f(x). With gammadot = gammadot_R*(R/r)^3 and the substitution
% u = R/r, the stress integral and its time derivative become
%   I1 = int_R^Inf gammadot/r*f dr     = gammadot_R*int_0^1 u^2*f(xR*u^3) du
%   I2 = int_R^Inf (f + x*f')/r^4 dr   = R^-3*int_0^1 u^2*(f + x*f') du
%      = f(xR)/(3*R^3),
% since u^2*du = dx/(3*xR) and f + x*f' = d(x*f)/dx. I1 is evaluated with a
% fixed Gauss-Legendre rule on panels that halve toward u = 0, where the
% integrand varies on the scale u ~ xR^(-1/3) at high shear rates.
% Outputs: f at the wall, intf (enters the stress integral as
% -6*intf/DRe), the Rdot part (dintf) and Rddot coefficient (ddintf) of
% d(intf)/dt, and, when the radial positions yq = r/R are given, f at those
% positions from the local shear rate gammadot_R*(R/r)^3 (fq).
function [f,intf,dintf,ddintf,fq] = f_viscosity(nu_model,Rdot,R,a,nc,lambda,yq)
    
    persistent uq wq
    if isempty(uq)
        [uq, wq] = f_graded_gauss(14, 10);
    end
    
    gammadot_R   = -2*Rdot/R;
    % d(gammadot_R*R^3)/dt = dgammadot + ddgammadot*Rddot
    dgammadot    = -4*Rdot*Rdot*R;
    ddgammadot   = -2*R*R;
    
    % shear-rate parameter at the bubble wall
    xR = lambda*abs(gammadot_R);
    
    % viscosity at the wall
    f = f_model(xR);
    % Leibniz integration rule term at the lower limit r = R
    h = f*gammadot_R*(Rdot/R);
    % stress integral, goes into the stress term of the radial equation
    I1 = gammadot_R*(wq'*(uq.^2.*f_model(xR*uq.^3)));
    % time derivative of the stress integral, in closed form
    I2 = f/(3*R^3);
    
    % S_def = (6/DRe)*I1 enters f_stress as -6*intf*iDRe  =>  intf = -I1
    intf = -I1;
    % d/dt intf = -(dg/dt)*I2 + Leibniz boundary term h; dg/dt splits
    % into the Rdot part (dintf) and the Rddot coefficient (ddintf)
    dintf = -dgammadot*I2 + h;
    ddintf = -ddgammadot*I2;
    
    % viscosity factor at given radial positions
    if nargin > 6
        fq = f_model(xR*yq.^-3);
    end
    
    % f(x) for the selected model, x = lambda*|gammadot| >= 0
    function fx = f_model(x)
        switch nu_model
            case 1
            % Carreau
            fx = (1 + x.^2).^((nc-1)/2);
            case 2
            % Carreau-Yasuda
            fx = (1 + x.^a).^((nc-1)/a);
            case 3
            % Powell-Eyring, asinh(x)/x -> 1 as x -> 0
            fx = ones(size(x));
            nz = x > 0;
            fx(nz) = asinh(x(nz))./x(nz);
            case 4
            % modified Powell-Eyring, log(1+x)/x^nc (nc <= 1)
            fx = zeros(size(x)) + (nc == 1);
            nz = x > 0;
            fx(nz) = log1p(x(nz))./x(nz).^nc;
            case 5
            % Cross
            fx = 1./(1 + x.^nc);
            case 6
            % simplified Cross
            fx = 1./(1 + x);
            case 7
            % modified Cross
            fx = (1 + x.^nc).^(-a);
            otherwise
            error('nu_model %d is not available', nu_model);
        end
    end
    
end

% Gauss-Legendre rule with np points on each of the panels [0, 2^-L],
% [2^-L, 2^-(L-1)], ..., [1/2, 1]; nodes from the Golub-Welsch eigenvalue
% problem
function [u, w] = f_graded_gauss(L, np)
    k = 1:np-1;
    beta = k./sqrt(4*k.^2 - 1);
    [V, D] = eig(diag(beta,1) + diag(beta,-1));
    x = diag(D);
    wx = 2*V(1,:)'.^2;
    edges = [0, 2.^(-L:0)];
    npan = numel(edges) - 1;
    u = zeros(np*npan, 1);
    w = zeros(np*npan, 1);
    for i = 1:npan
        lo = edges(i);
        hi = edges(i+1);
        idx = (i-1)*np + (1:np);
        u(idx) = (hi - lo)/2*x + (hi + lo)/2;
        w(idx) = (hi - lo)/2*wx;
    end
end
