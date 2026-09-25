% file f_imr_inputs.m
% brief contains function f_imr_inputs

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function returns every solver input with its default value, as
% a struct whose field names are the input names of the name/value pairs
% (p.re is the input 're'). Change fields and pass the struct to a solver,
% optionally followed by name/value pairs that override it:
%   p = f_imr_inputs();
%   p.radial = 2;
%   p.re = 300;
%   [t, R] = f_imr_nodal(p);
% Only the fields that differ from these defaults are passed on, so the
% struct and the name/value pairs describe a run the same way. The
% optional inputs without a default, p0star and stress0, can be added as
% fields. The values come from default_case.
function p = f_imr_inputs()
    
    d = default_case();
    % names in default_case that differ from the input names (the input
    % names are those of the name/value parser in f_call_params)
    % (one pair per row, rows separated by line breaks)
    renamed = {
    'Nv', 'nv'
    'Nt', 'nt'
    'Mt', 'mt'
    'Lv', 'lv'
    'Lt', 'lt'
    'rtol', 'reltol'
    'atol', 'abstol'
    'TFin', 'tfin'
    'TVector', 'tvector'
    'Req', 'req'
    'Rdotzero', 'rdotzero'
    'Cstar', 'cstar'
    'GAMa', 'gama'
    'We', 'we'
    'Re8', 're'
    'Ca', 'ca'
    'De', 'de'
    'LAM', 'lam'
    'DRe', 'dre'
    'v_lambda_star', 'v_lambda'
    'Foh', 'foh'
    'Br', 'br'
    'Km_star', 'kmstar'
    'Fom', 'fom'
    'L_heat_star', 'lheat'
    'Rv_star', 'rv'
    'Ra_star', 'ra'
    'Pv_sat', 'pv'
    };
    p = struct();
    f = fieldnames(d);
    for k = 1:numel(f)
        j = find(strcmp(renamed(:,1), f{k}), 1);
        if isempty(j)
            name = f{k};
        else
            name = renamed{j,2};
        end
        p.(name) = d.(f{k});
    end
    
end
