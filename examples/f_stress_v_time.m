% file f_stress_v_time.m
% brief contains function f_stress_v_time

% brief This function interpolates stress of a graded or nongraded material
% at dynamic tracer(s) position(s)
function [gstress_all] = f_stress_v_time(isgraded,nt,nloc,r_coord,R,R0,Req,r_far,t,Ca,Ca1,l1,l2,v_nc,v_a)

% for first two columns of stress_all
gstress_all = zeros(nt,nloc);
for time_idx = 1:nt
    t_now = t(time_idx);
    % current time
    Rnow = R(time_idx);
    
    % if not graded, compute homogeneous stress
    if isgraded
        % get graded stress at tracer location

        % get stress profile at current time for graded material
        [taurr_now,r1_now,r2_now] = f_g_stress(r_coord,Rnow,Req,Ca,Ca1,l1,l2,v_nc,v_a);
    else
        % get stress profile at current time for nongraded material
        [taurr_now,r1_now,r2_now] = f_ng_stress(r_coord,Rnow,Req,Ca,l1,l2);
    end

    % if any(~isreal(taurr_now),'all')
    %     fprintf('warning: complex stress detected at this timestep')
    % end

    % location of each tracer (at each time)
    r_eval = [R0,r1_now,r2_now,r_far];

    % interpolate stress at each location of each tracer (at each time)
    for loc=1:nloc
        r_current = r_eval(loc);
        %warning('off','MATLAB:interp1:NaNstrip')
        % if r_current < min(r_coord) || r_current > max(r_coord)
        %     fprintf('warning: r_eval(%d) = %.4f is outside range [%.4f, %.4f]\n',...
        %         loc, r_current, min(r_coord), max(r_coord));
        % end
        gstress_all(time_idx,loc) = interp1(r_coord,taurr_now,r_current,'pchip');
    end
end
end