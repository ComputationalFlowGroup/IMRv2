% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% The time-integrator code (method = 23, ode23tb) that every test and
% example not itself comparing integrators should pin to. method = 3
% (CVODES, the default) needs the ode object's cvodesstiff solver, which
% requires MATLAB R2024a or later and is not available in CI; using this
% one value, instead of the bare literal, keeps that constraint in a
% single place if it ever changes.
function m = ci_safe_method()
    m = 23;
end
