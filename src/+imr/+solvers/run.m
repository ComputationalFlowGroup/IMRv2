% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Shared public orchestration. Runs once per simulation, outside the RHS.
function varargout = run(solver,varargin)
    prepared = prepare_case(solver,varargin{:});
    system = imr.solvers.build(prepared);
    n = prepared.numerics;
    [t,X] = imr.numerics.integrate(system.rhs,system.initial,n.method, ...
        n.divisions,prepared.time,system.jacobian,n.rtol,n.atol);
    varargout = system.outputs(t,X);
end
