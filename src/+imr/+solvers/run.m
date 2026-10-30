% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% Shared public orchestration. Runs once per simulation, outside the RHS.
function varargout = run(solver,varargin)
    if nargout > 11
        error('IMR:outputs',['IMR returns at most 11 outputs: ' ...
            '[t,R,Rdot,P,T,Tm,kv,inputs,Rddot,S,Sdot].']);
    end
    prepared = imr.configuration.prepare_case(solver,varargin{:});
    system = imr.solvers.build(prepared);
    n = prepared.numerics;
    [t,X] = imr.numerics.integrate(system.rhs,system.initial,n.method, ...
        n.divisions,prepared.time,system.jacobian,n.rtol,n.atol);
    % Pass nargout to outputs so expensive reconstructions (wall resolution,
    % history computation) are skipped if their outputs are not requested
    varargout = system.outputs(t,X,nargout);
    if nargout >= 9
        % Reuse the coupled RHS at the actual returned states and times.
        % No differentiation of sampled velocity, new ODE state or extra
        % RHS evaluations for calls requesting only the original outputs.
        % S and Sdot (the stress integral and its time derivative) are
        % already computed on every right-hand-side call to feed the
        % radial equation, so outputs 10 and 11 cost nothing beyond what
        % output 9 already evaluates.
        if isfield(system,'reset_output_rhs')
            system.reset_output_rhs(X(1,:).');
        end
        acceleration = zeros(numel(t),1);
        S = zeros(numel(t),1);
        Sdot = zeros(numel(t),1);
        for k = 1:numel(t)
            [rate,S(k),Sdot(k)] = system.output_rhs(t(k),X(k,:).');
            acceleration(k) = rate(system.layout.Rdot);
        end
        varargout{9} = acceleration;
        if nargout >= 10
            varargout{10} = S;
        end
        if nargout >= 11
            varargout{11} = Sdot;
        end
    end
end
