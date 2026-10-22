% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Internal numerical system for a checked imr.configuration.prepare_case result.
% initial: column state, rhs(t,X): column or nodal batch, jacobian: function
% or sparsity matrix, outputs(t,X): the eight public outputs as a cell row.
% No integration is performed. This contract also permits independent RHS
% and Jacobian checks without adding diagnostic arguments to public solvers.
function system = build(c)
    if strcmp(c.solver,'nodal')
        system = imr.solvers.build_nodal(c);
    elseif strcmp(c.solver,'spectral')
        system = imr.solvers.build_spectral(c);
    else
        error('IMR:solver','Solver must be nodal or spectral.');
    end
    system.layout = c.layout;
end
