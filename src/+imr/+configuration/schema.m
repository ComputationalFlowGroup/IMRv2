% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Internal parameter names, public names and unconditional discrete rules.
% imr.configuration.defaults owns values, imr.configuration.validate_case owns coupled/physical constraints.
% Used by the public struct, parser, validation and documentation checks.
function schema = schema()
    data = {
    'radial', 'radial', 'choice', 1:6
    'bubtherm', 'bubtherm', 'choice', 0:1
    'medtherm', 'medtherm', 'choice', 0:1
    'stress', 'stress', 'choice', 0:7
    'eps3', 'eps3', '', []
    'vapor', 'vapor', 'choice', 0:1
    'masstrans', 'masstrans', 'choice', 0:1
    'tfin', 'TFin', '', []
    'tvector', 'TVector', '', []
    'method', 'method', 'choice', [3 15 23]
    'spectral', 'spectral', 'choice', 0:1
    'divisions', 'divisions', 'count', []
    'nt', 'Nt', 'count', []
    'mt', 'Mt', 'count', []
    'nv', 'Nv', 'count', []
    'lv', 'Lv', '', []
    'lt', 'Lt', '', []
    'reltol', 'rtol', '', []
    'abstol', 'atol', '', []
    'jacobian', 'jacobian', 'choice', 0:3
    'scheme', 'scheme', 'choice', 0:1
    'collapse', 'collapse', 'choice', 0:1
    'rdotzero', 'Rdotzero', '', []
    'req', 'Req', '', []
    'progdisplay', 'progdisplay', 'choice', 0:1
    'cstar', 'Cstar', '', []
    'gama', 'GAMa', '', []
    'nstate', 'nstate', '', []
    'hugoniot_s', 'hugoniot_s', '', []
    'ee', 'ee', '', []
    'om', 'om', '', []
    'tw', 'tw', '', []
    'dt', 'dt', '', []
    'mn', 'mn', '', []
    'wave_type', 'wave_type', 'choice', -3:3
    'we', 'We', '', []
    're', 'Re8', '', []
    'ca', 'Ca', '', []
    'de', 'De', '', []
    'lam', 'LAM', '', []
    'alphax', 'alphax', '', []
    'dre', 'DRe', '', []
    'v_a', 'v_a', '', []
    'v_nc', 'v_nc', '', []
    'v_lambda', 'v_lambda_star', '', []
    'nu_model', 'nu_model', 'choice', 0:7
    'kappa', 'kappa', '', []
    'foh', 'Foh', '', []
    'br', 'Br', '', []
    'alpha_g', 'alpha_g', '', []
    'beta_g', 'beta_g', '', []
    'alpha_v', 'alpha_v', '', []
    'beta_v', 'beta_v', '', []
    'chi', 'chi', '', []
    'kmstar', 'Km_star', '', []
    'fom', 'Fom', '', []
    'lheat', 'L_heat_star', '', []
    'rv', 'Rv_star', '', []
    'ra', 'Ra_star', '', []
    'pv', 'Pv_sat', '', []
    'thetav', 'thetav', '', []
    };
    schema = struct('name',data(:,1), 'field',data(:,2), ...
        'rule',data(:,3), 'allowed',data(:,4));
end
