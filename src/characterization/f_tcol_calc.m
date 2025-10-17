% file f_tcol_calc.m
% brief contains function f_tcol_calc

% brief This function features the collapse time solver for a nongraded
% elastic material. The solver currently assumes the Kelvin-Voigt with neo-Hookean
% elasticity, using Yang 2024.
function [tcol] = f_tcol_calc(stress,R0,Ca,Pref,rho8)
    if stress == 1
        tcol = ((2/3 + 5/(3*Ca))^(-0.5)) * 0.747 * R0 * sqrt(rho8/Pref);
        % DON'T YOU WANT TO COMPARE THIS AGAINST RAYLEIGH COLLAPSE??
    end
end
