% file far_field_pressure.m
% brief contains function imr.physics.far_field_pressure

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the time-dependent external pressure in the
% surrounding material that drives the bubble oscillations. The function
% features wave types: experimentally-derived, impulse, Gaussian,
% histotripsy, and Heaviside impulse. The time t may be a scalar or an
% array. vararg is the cell {[om ee tw dt mn wave_type], poly, dpoly}. The
% solvers call it on every right-hand-side evaluation, so each waveform is
% evaluated inline, without nested function calls.
function [p8,p8dot] = far_field_pressure(t,vararg)

    % waveform parameters [om ee tw dt mn wave_type]
    wp = vararg{1};
    wave_type = wp(6);

    if wave_type < 0
        % experimentally-derived incoming waveform
        p8 = ppval(vararg{2},t);
        p8dot = ppval(vararg{3},t);
    elseif wave_type == 0
        % impulse wave
        p8 = wp(2)*ones(size(t));
        p8dot = zeros(size(t));
    elseif wave_type == 1
        % Gaussian wave
        ee = wp(2);
        tw = wp(3);
        dt = wp(4);
        g = exp(-((t-dt).^2) / tw^2);
        p8 = - ee * g;
        p8dot = ee * (2*(t-dt)/tw^2) .* g;
    elseif wave_type == 2
        % histotripsy wave
        om = wp(1);
        ee = wp(2);
        dt = wp(4);
        mn = wp(5);
        p8 = zeros(size(t));
        p8dot = zeros(size(t));
        if mn == 0
            % Rectangular pulse: use the same zero pointwise derivative
            % convention as the Heaviside waveform, including its edges.
            on = (t >= dt - pi/om) & (t <= dt + pi/om);
            p8(on) = ee;
        elseif ee ~= 0
            % The half-angle form avoids cancellation near the pulse edges.
            % Set the edge derivative to zero. For mn <= 1/2 the derivative
            % has a jump or a singular one-sided limit at those edges.
            on = (t > dt - pi/om) & (t < dt + pi/om);
            phase = om*(t(on)-dt)/2;
            c = cos(phase);
            p8(on) = ee*c.^(2*mn);
            p8dot(on) = -ee*mn*om*c.^(2*mn-1).*sin(phase);
        end
    elseif wave_type == 3
        % heaviside impulse, the step at t = tw is measure-zero in time
        ee = wp(2);
        tw = wp(3);
        p8 = -ee*(t < tw);
        p8dot = zeros(size(t));
    end

end
