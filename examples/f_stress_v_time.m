% file f_stress_v_time.m
% brief contains function f_stress_v_time

% brief This function interpolates stress of a graded or nongraded material
% at dynamic tracer(s) position(s)
function [stress_all] = f_stress_v_time(isgraded,nt,nloc,R,R0,Req,r_far,Ca,Ca1,l1,l2,v_nc,v_a)

stress_all = zeros(nt,nloc);
for time_idx = 1:nt
    % current time
    Rnow = R(time_idx);
    %rcoord_now = r_coord(time_idx);
    r1_now = (l1^3 + Rnow.^3 - Req^3).^(1/3);
    r2_now = (l2^3 + Rnow.^3 - Req^3).^(1/3);
    % location of each tracer (at each time)
    r_eval = [R0,r1_now,r2_now,r_far];
    
    for loc = 1:nloc
        r_current = r_eval(loc);
        %r0coord = (r_coord.^3 - R.^3 + Req.^3).^(1/3);
        r0_eval = (r_current.^3 - Rnow.^3 + Req.^3).^(1/3);
        % if not graded, compute homogeneous stress
        if isgraded
            % get stress profile at current time for graded material
            taurr_now = f_g_stress(r_current,r0_eval,Ca,Ca1,l1,l2,v_nc,v_a);
        else
            % get stress profile at current time for nongraded material
            taurr_now = (2/(3*Ca))*((r0_eval.^4 ./ r_current^4) - (r_current.^2 ./ r0_eval.^2));
        end
        stress_all(time_idx,loc) = taurr_now;
    end
end
end