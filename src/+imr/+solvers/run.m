% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Shared public orchestration. Runs once per simulation, outside the RHS.
function varargout = run(solver,varargin)
    if nargout > 9
        error('IMR:outputs','IMR returns at most 9 outputs: [t,R,Rdot,P,T,Tm,kv,inputs,Rddot].');
    end
    prepared = imr.configuration.prepare_case(solver,varargin{:});
    system = imr.solvers.build(prepared);
    n = prepared.numerics;
    [t,X] = imr.numerics.integrate(system.rhs,system.initial,n.method, ...
        n.divisions,prepared.time,system.jacobian,n.rtol,n.atol);
    varargout = system.outputs(t,X);
    if nargout == 9
        % Reuse the coupled RHS at the actual returned states and times.
        % No differentiation of sampled velocity, new ODE state or extra
        % RHS evaluations for calls requesting only the original outputs.
        if isfield(system,'reset_output_rhs')
            system.reset_output_rhs(X(1,:).');
        end
        acceleration = zeros(numel(t),1);
        for k = 1:numel(t)
            rate = system.output_rhs(t(k),X(k,:).');
            acceleration(k) = rate(system.layout.Rdot);
        end
        varargout{9} = acceleration;
    end
end
