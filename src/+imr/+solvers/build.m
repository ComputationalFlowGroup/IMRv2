% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Internal numerical system for a checked imr.configuration.prepare_case result.
% initial: column state, rhs(t,X): column or nodal batch, jacobian: function
% or sparsity matrix, outputs(t,X): the eight public outputs as a cell row.
% output_rhs(t,X): the same derivatives with progress display suppressed,
% used only when the shared runner is asked to return acceleration.
% reset_output_rhs(X): optional wall-cache reset before a new output history.
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
    % A valid list of scalar inputs can still assemble an inadmissible
    % initial state: e.g. a sonic denominator or pressure outside an EoS
    % domain. Evaluate once here, rather than adding checks to each RHS
    % evaluation during integration. Underlying programming errors propagate.
    initial_rate = system.rhs(c.time(1),system.initial);
    if ~isreal(initial_rate) || any(~isfinite(initial_rate(:)))
        error('IMR:initialState', ...
            ['INPUT ERROR: the assembled initial state has a nonfinite or complex derivative. ' ...
            'Check initial pressure, velocity, stress, liquid EoS domain and forcing.']);
    end
end
