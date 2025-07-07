% file f_g_stress.m
% brief contains function f_g_stress

% brief This function sets computes the graded stress field 
% at the desired r_eval location as a function of time
function [taurr] = f_g_stress(r_eval,r0_eval,Ca,Ca1,l1,l2,v_nc,v_a)

% % choose subdomain that r_eval lives in
% mask = rcoord <= r_eval;
% % truncate rcoord to appropriate region
% r_subset = rcoord(mask);
% compute r0coord for the subdomain
%r0_eval = (r_eval.^3 - R.^3 + Req.^3).^(1/3);
% slight offset near bubble wall
% r0shift = (rcoord.^3 + Req^3 - (R+0.001).^3).^(1/3);
% removing the data within the bubble and slightly away from the bubble wall
% r0_eval(r0_eval < r0shift) = NaN;

% function for graded region
f_cy = @(x) (l2 - x) ./ (x - l1);

% initialize taurr
taurr = 0;

% compute stress according to the region
if r0_eval < l1
    % near field
    taurr = (2/(3*Ca))*((r0_eval.^4 ./ r_eval.^4) - (r_eval.^2 ./ r0_eval.^2));

elseif (r0_eval <= l2 && r0_eval >= l1)
     % compute graded stress
     taurr = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy(r0_eval).^v_a).^((v_nc-1)/v_a)).* ...
            ((r0_eval.^4 ./ r_eval.^4) - (r_eval.^2 ./ r0_eval.^2));

elseif r0_eval > l2
    % far field
    taurr = (2/(3*Ca1))*((r0_eval.^4 ./ r_eval.^4) - (r_eval.^2 ./ r0_eval.^2));
end

% removing any infs or nans
taurr(isnan(taurr) | isinf(taurr)) = 0;
end
