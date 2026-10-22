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
    defaults = default_case();
    schema = f_parameter_schema();
    p = struct();
    for k = 1:numel(schema)
        p.(schema(k).name) = defaults.(schema(k).field);
    end
end
