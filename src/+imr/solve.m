% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Solve a dimensionless case with 'solver','nodal' (default) or 'spectral'.
% Accept name/value pairs, or an input struct followed by overriding pairs.
% Outputs: [t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot]. Rddot, S (the stress
% integral) and Sdot (its time derivative) are optional and are evaluated
% from the governing equations only when requested (outputs 9-11).
% The selected solver is recorded in inputs.solver for replay.
function varargout = solve(varargin)
    [varargout{1:nargout}] = imr.solvers.run('',varargin{:});
end
