% file f_display.m
% brief contains function f_display

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features the display output on the command window
function [] = f_display(radial, bubtherm, medtherm, masstrans, stress, ...
    spectral, nu_model, eps3, Pv_star, Re8, De, Ca, LAM, mode)

% Command window display

if radial == 1
    eqn = 'Rayleigh Plesset equation';
elseif radial == 2
    eqn = 'Keller-Miksis pressure';
elseif radial == 3
    eqn = 'Keller-Miksis enthalpy';
elseif radial == 4
    eqn = 'Gilmore';
elseif radial == 5
    eqn = 'Keller-Miksis in enthalpy with Mie-Gruneisen EoS';
elseif radial == 6
    eqn = 'Gilmore with Mie-Gruneisen EoS';
else
    error('INPUT ERROR: radial must be 1-6');
end
const = 'none';

if bubtherm == 1
    if medtherm == 1
        therm = 'bubble and medium heat transfer';
    else
        therm = 'bubble heat transfer, medium at the far-field temperature';
    end
else
    therm = 'polytropic approximation';
end

if masstrans == 1
    mass = 'mass transfer in bubble';
else
    mass = 'no mass transfer in bubble';
end

if Pv_star ~= 0
    vap = 'vapor in bubble';
else
    vap = 'no vapor in bubble';
end

% the elastic law of stress 2 and 4, which acts only with an elastic
% modulus (finite Ca)
if stress == 2 || stress == 4
    elastic = 'quadratic neo-Hookean';
else
    elastic = 'neo-Hookean';
end
if stress == 0
    const = 'no stress applied';
elseif stress == 1 || stress == 2
    if Ca == Inf && Re8 == Inf
        const = 'inviscid fluid';
    elseif Ca == Inf
        const = 'Newtonian fluid';
    elseif Re8 == Inf
        const = [elastic ' elastic solid (inviscid)'];
    else
        const = [elastic ' Kelvin-Voigt'];
    end
elseif stress == 3 || stress == 4
    % a spring (finite Ca) with a Maxwell element makes the Zener model; a
    % solvent viscosity (LAM > 0) adds a dashpot in parallel
    if Ca == Inf && LAM == 0
        const = 'linear Maxwell';
    elseif Ca == Inf
        const = 'linear Jeffreys';
    elseif LAM == 0
        const = ['linear Zener, ' elastic ' spring'];
    else
        const = ['linear Zener with solvent viscosity, ' elastic ' spring'];
    end
elseif stress == 5
    if LAM == 0
        const = 'upper-convected Maxwell';
    else
        const = 'Oldroyd-B';
    end
    if spectral == 1
        const = [const ' (spectral stress field)'];
    end
end
if eps3 ~= 0
    const = [const ', eps3 = ' num2str(eps3)];
end

if stress == 0
    visco = 'none (no stress applied)';
elseif nu_model == 0
    visco = 'Newtonian';
elseif nu_model == 1
    visco = 'Carreau';
elseif nu_model == 2
    visco = 'Carreau-Yasuda';
elseif nu_model == 3
    visco = 'Powell-Eyring';
elseif nu_model == 4
    visco = 'Modified Powell-Eyring';
elseif nu_model == 5
    visco = 'Cross';
elseif nu_model == 6
    visco = 'Simplified Cross';
elseif nu_model == 7
    visco = 'Modified Cross';
end

% how the stress of the history-dependent models is carried: the stress
% field in spectral modes (stress 5, spectral = 1, f_imr_spectral), or
% ordinary differential equations for its integrals
if stress < 3
    solut = 'closed form (no stress history)';
elseif spectral == 1
    solut = 'stress field, spectral modes';
else
    solut = 'stress integrals (ODEs)';
end

% display run settings
disp(['Mode: ' mode]);
disp(['Radial dynamics: ' eqn]);
disp(['Constitutive model: ' const]);
disp(['Viscosity rheology: ' visco]);
disp(['Thermal effects: ' therm]);
disp(['Mass effects: ' mass]);
disp(['Vapor effects: ' vap]);
disp(['Stress representation: ' solut]);
disp('--- Dimensionless numbers ---');
disp(['Re8 = ' num2str(Re8,'%10.10f')]);
disp(['De = ' num2str(De,'%10.10f')]);
disp(['Ca = ' num2str(Ca,'%10.10f')]);
disp(['LAM = ' num2str(LAM,'%10.10f')]);

end
