% file f_display.m
% brief contains function f_display

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function features the display output on the command window.
% The inputs after mode describe the non-Newtonian viscosity: DRe, iDRe
% (the factor of its stress term, 0 when it is inactive), v_lambda, v_nc
% and v_a.
function [] = f_display(radial, bubtherm, medtherm, masstrans, stress, ...
    spectral, nu_model, eps3, Pv_star, Re8, De, Ca, LAM, mode, DRe, iDRe, ...
    v_lambda, v_nc, v_a)

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
else
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

% non-Newtonian viscosity: the model name, and whether it acts. Its stress
% term carries the factor iDRe = 1/DRe (stress 1-2) or LAM/DRe (stress 3-5,
% where it is the solvent share), which vanishes for DRe = 0 or, for stress
% 3-5, LAM = 0
nu_names = {'Carreau', 'Carreau-Yasuda', 'Powell-Eyring', ...
    'Modified Powell-Eyring', 'Cross', 'Simplified Cross', 'Modified Cross'};
nu_params = '';
if stress == 0
    visco = 'none (no stress applied)';
elseif nu_model == 0
    visco = 'Newtonian';
else
    visco = nu_names{nu_model};
    if DRe == 0
        visco = [visco ' selected but inactive (DRe = 0, mu0 = mu8): ' ...
            'Newtonian'];
    elseif iDRe == 0
        visco = [visco ' selected but inactive (LAM = 0, no solvent ' ...
            'share for stress 3-5): Newtonian'];
    else
        if stress >= 3
            visco = [visco ', solvent share (LAM)'];
        end
        if v_lambda == 0
            visco = [visco ', v_lambda = 0 so constant viscosity mu0'];
        end
        % the parameters each model uses
        nu_params = ['DRe = ' num2str(DRe,'%10.10f') ...
            ', v_lambda = ' num2str(v_lambda,'%10.10f')];
        if any(nu_model == [1 2 4 5 7])
            nu_params = [nu_params ', v_nc = ' num2str(v_nc,'%10.10f')];
        end
        if any(nu_model == [2 7])
            nu_params = [nu_params ', v_a = ' num2str(v_a,'%10.10f')];
        end
    end
end

% how the stress of the history-dependent models is carried: the stress
% field in spectral modes (stress 5, spectral = 1, f_imr_spectral), or
% ordinary differential equations for its integrals, which stress 3 and 4
% with medium heat transfer complement with the local polymer stress field
% for the medium heating (f_zener_field)
if stress < 3
    solut = 'closed form (no stress history)';
elseif spectral == 1
    solut = 'stress field, spectral modes';
elseif stress <= 4 && medtherm == 1
    solut = 'stress integral (ODE), and its local field for the medium heating';
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
if ~isempty(nu_params)
    disp(nu_params);
end

end
