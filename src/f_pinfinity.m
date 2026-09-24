% file f_pinfinity.m
% brief contains function f_pinfinity

% brief This function computes the time-dependent external pressure in the
% surrounding material that drives the bubble oscillations. The function
% features wave types: experimentally-derived, impulse, Gaussian,
% histotripsy, and Heaviside impulse. The time t may be a scalar or an
% array. vararg is the cell {[om ee tw dt mn wave_type], poly, dpoly}.
function [p8,p8dot] = f_pinfinity(t,vararg)
    
    % extracting input variables
    om = vararg{1}(1);
    ee = vararg{1}(2);
    tw = vararg{1}(3);
    dt = vararg{1}(4);
    mn = vararg{1}(5);
    wave_type = vararg{1}(6);
    
    if wave_type < 0
        % experimentally-derived incoming waveform
        p8 = ppval(vararg{2},t);
        p8dot = ppval(vararg{3},t);
    elseif wave_type == 0
        % impulse wave
        [p8, p8dot] = impulse(t);
    elseif wave_type == 1
        % Gaussian wave
        [p8, p8dot] = gaussian(t);
    elseif wave_type == 2
        % histotripsy wave
        [p8, p8dot] = histo(t);
    elseif wave_type == 3
        % heaviside impulse
        [p8,p8dot] = heaviside_impulse(t);
    end
    
    % histotripsy waveform
    function [p,pdot] = histo(t)
        on = (t >= dt - pi/om) & (t <= dt + pi/om);
        p = zeros(size(t));
        pdot = zeros(size(t));
        p(on) = ee*(0.5 + 0.5*cos(om*(t(on) - dt))).^mn;
        pdot(on) = -ee*mn*(0.5+0.5*cos(om*(t(on)-dt))).^(mn-1)*0.5*om.* ...
            sin(om*(t(on)-dt));
    end
    
    % impulse waveform
    function [p,pdot] = impulse(t)
        p = ee*ones(size(t));
        pdot = zeros(size(t));
    end
    
    % Gaussian waveform
    function [p,pdot] = gaussian(t)
        p = - ee * exp(-((t-dt).^2) / tw^2);
        pdot = ee * (2*(t-dt)/tw^2) .* exp(-((t-dt).^2) / tw^2);
    end
    
    % heaviside impulse, the step at t = tw is measure-zero in time
    function [p,pdot] = heaviside_impulse(t)
        p = -ee*(t < tw);
        pdot = zeros(size(t));
    end
    
end
