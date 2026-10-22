% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Internal state contract, shared by assembly, RHS, Jacobians and outputs.
% Nodal fields store node values; spectral fields store order+1 coefficients.
% Each index block is a row vector; state values are columns (or batches).
function layout = state_layout(c)
    m = c.model;
    n = c.numerics;
    spectral = strcmp(c.solver,'spectral');
    counts = [m.bubtherm*(n.Nt+spectral), m.medtherm*(n.Mt+spectral), ...
        m.masstrans*n.Nt, 0, 0, 0];
    if m.stress == 3 || m.stress == 4
        counts(4) = 1;
        if m.medtherm
            counts(6) = n.Mt;
            if spectral
                counts(6) = n.Nv;
            end
        end
    elseif m.stress == 5
        counts(4:5) = n.Nv;
    end
    names = {'bubble','medium','vapor','stress1','stress2','polymer'};
    layout.R = 1;
    layout.Rdot = 2;
    layout.P = 3;
    offset = 3;
    for k = 1:numel(names)
        layout.(names{k}) = offset+(1:counts(k));
        offset = offset+counts(k);
    end
    layout.stress = [layout.stress1 layout.stress2];
    layout.allStress = [layout.stress layout.polymer];
    layout.count = offset;
end
