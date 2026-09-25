% file f_pinfinity.m
% brief contains function f_pinfinity

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function computes the time-dependent external pressure in the
% surrounding material that drives the bubble oscillations. The function
% features wave types: experimentally-derived, impulse, Gaussian,
% histotripsy, and Heaviside impulse. The time t may be a scalar or an
% array. vararg is the cell {[om ee tw dt mn wave_type], poly, dpoly}. The
% solvers call it on every right-hand-side evaluation, so each waveform is
% evaluated inline, without nested function calls.
function [p8,p8dot] = f_pinfinity(t,vararg)
    
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
        on = (t >= dt - pi/om) & (t <= dt + pi/om);
        p8 = zeros(size(t));
        p8dot = zeros(size(t));
        p8(on) = ee*(0.5 + 0.5*cos(om*(t(on) - dt))).^mn;
        p8dot(on) = -ee*mn*(0.5+0.5*cos(om*(t(on)-dt))).^(mn-1)*0.5*om.* ...
            sin(om*(t(on)-dt));
    elseif wave_type == 3
        % heaviside impulse, the step at t = tw is measure-zero in time
        ee = wp(2);
        tw = wp(3);
        p8 = -ee*(t < tw);
        p8dot = zeros(size(t));
    end
    
end
