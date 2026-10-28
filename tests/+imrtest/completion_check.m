function [is_complete, reason, warnings] = completion_check(t, t_request, tol)
% Check if a simulation completed to the requested final time.
%
% Validates that the integration reached the requested endpoint within
% floating-point tolerance. Catches early exits due to NaN, event triggers,
% or integration failures that would otherwise pass silently in test suites.
%
% Syntax:
%   [is_complete, reason] = imrtest.completion_check(t, t_request)
%   [is_complete, reason, warnings] = imrtest.completion_check(t, t_request, tol)
%
% Inputs:
%   t          - Time vector (column or row, in SI units)
%   t_request  - Requested final integration time (SI units)
%   tol        - Relative tolerance (default: 1e-9)
%                Completion criterion: t(end) >= t_request - tol*t_request
%
% Outputs:
%   is_complete - Logical: true if t(end) >= t_request - tol*t_request
%   reason      - String describing completion status or failure reason
%   warnings    - Cell array of warning strings (empty if none)
%
% Example:
%   [t, R] = imr.solve_nodal(...);
%   [ok, msg, warns] = imrtest.completion_check(t, tfin, 1e-10);
%   if ~ok
%       error(['Integration failed: ', msg]);
%   end

    arguments
        t (:,:) double {mustBeVector}
        t_request (1,1) double {mustBePositive}
        tol (1,1) double = 1e-9
    end

    % Ensure t is a column vector
    if isrow(t)
        t = t(:);
    end

    % Validate input time vector
    if isempty(t)
        is_complete = false;
        reason = 'Time vector is empty';
        return;
    end

    if ~ismonotonic(t)
        is_complete = false;
        reason = 'Time vector is not monotonically increasing';
        return;
    end

    t_final = t(end);
    threshold = t_request - tol * t_request;

    if t_final >= threshold
        is_complete = true;
        margin = t_final - t_request;
        if margin > tol * t_request
            reason = sprintf('Completed with %.3e overshoot', margin);
        else
            reason = sprintf('Completed at t=%.16e (requested %.16e)', t_final, t_request);
        end
    else
        is_complete = false;
        shortfall = t_request - t_final;
        reason = sprintf('Incomplete: reached t=%.16e, requested %.16e (shortfall %.3e)', ...
                        t_final, t_request, shortfall);
    end

end

function tf = ismonotonic(v)
    tf = all(diff(v) >= 0);
end
