% file f_max_pre_stress.m
% brief contains function f_max_pre_stress

% brief This function computes the stress value at the bubble maximum
% radius. The implementation is the approximate/analytical model as opposed
% to the numerical model that performs a shooting method to determine the
% initial stress (see f_initial_stress_calc.m). This model should be
% treated as approximate and of lower fidelity relative to the numerical
% model.
function [Smaxpred] = f_approx_init_stress(Ro, kappa, al_nd, pwv_nd, We, Re, De, Ca, alpha)
    
    % trc constant
    trc = sqrt(6*pi)*gamma(11/6)/(5*gamma(4/3));
    
    % stress-related terms
    fbarst = pi/(sqrt(6)*We*trc);
    
    % constants
    B = 2.1844;
    
    % barotropic correction term
    fbarbc = -(1 - pwv_nd + 1/(We*Ro))*B*Ro^(3*kappa) + pwv_nd;
    
    % compressibility correction
    Mc = 1 / al_nd;
    fbarc = -2*Mc/(Mc + sqrt(Mc^2 + 4*trc^2));
    
    % viscous correction
    C = 0.46379 + 0.56391 / Re + 5.74916 / Re^2;
    fbarv = -4*C^2/(2*C^2 + sqrt(4*C^4 + C^2*Re^2*trc^2));
    
    % max stress correction
    fbarmax = fbarv + (De/trc)*(fbarv*exp(-trc/De) - fbarv);
    
    % elastic correction
    fbare = (1 / (60 * Ca * gamma(5/6))) * gamma(1/3) * ...
        ((40 * sqrt(pi) * Ro * (1 - 3 * alpha)) + ...
        (120 * (-1 + 2 * Ro^3) * alpha * gamma(7/6)) / (Ro * gamma(2/3)) + ...
        (-50 + 177 * alpha) * gamma(5/6) / gamma(4/3));
    
    % stress loss correction
    fbarsls = fbarmax - fbare;
    
    % total stress sum
    fsum = fbarbc + fbarst + fbarc + fbarsls;
    
    if fsum > 1
        error('INPUT ERROR: The Cauchy number lower than possible for prestress');
    end
    
    % time constant tg
    tg = (5*sqrt(pi)*gamma(5/6) - ...
        6*Ro^(5/2)*gamma(4/3)*hypergeom([1/2, 5/6], 11/6, Ro^3))/...
        (5*sqrt(6 - 6*fsum)*gamma(4/3));
    
    % predicted maximum stress
    Smaxpred = fbarv*(1 - exp(-tg / De));
    
end
