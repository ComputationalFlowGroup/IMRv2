% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Solve a dimensionless case with 'solver','nodal' (default) or 'spectral'.
% Accept name/value pairs, or an input struct followed by overriding pairs.
% Outputs: [t,R,Rdot,P,T,Tm,kv,inputs,Rddot]. Rddot is optional and is
% evaluated from the governing equations only when output 9 is requested.
% The selected solver is recorded in inputs.solver for replay.
function varargout = solve(varargin)
    [varargout{1:nargout}] = imr.solvers.run('',varargin{:});
end
