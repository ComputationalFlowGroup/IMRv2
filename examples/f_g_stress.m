% file f_g_stress.m
% brief contains function f_g_stress

% brief This function sets computes the graded stress field 
% at the desired r_eval location as a function of time
function [taurr,r1,r2] = f_g_stress(r_eval,rcoord,R,Req,Ca,Ca1,l1,l2,v_nc,v_a)

% choose subdomain that r_eval lives in
mask = rcoord <= r_eval;
% truncate rcoord to appropriate region
r_subset = rcoord(mask);
% compute r0coord for the subdomain
r0coord = (r_subset.^3 - R.^3 + Req.^3).^(1/3);
% slight offset near bubble wall
r0shift = (r_subset.^3 + Req^3 - (R+0.001).^3).^(1/3);
% removing the data within the bubble and slightly away from the bubble wall
r0coord(r0coord < r0shift) = NaN;

% function for graded region
f_cy = @(r0coord) (l2 - r0coord) ./ (r0coord - l1);

% initialize taurr
taurr = zeros(size(r_subset));

% compute stress according to the region
if r_eval <= l1
    % near field
    taurr = (2/(3*Ca))*((r0coord.^4 ./ r_subset.^4) - (r_subset.^2 ./ r0coord.^2));
    % zero out the stress above and below the near field
    taurr(isinf(taurr)) = 0;

elseif r_eval <= l2
    % near field
     near_mask = r0coord <= l1;
     r0near = r0coord(near_mask);
     r_near = r_subset(near_mask);
     % compute near stress
     taurr1 = (2/(3*Ca))*((r0near.^4 ./ r_near.^4) - (r_near.^2 ./ r0near.^2));
     % zero out the stress above and below the near field
     taurr1(isinf(taurr1)) = 0;

     % graded region
     mid_mask = r0coord > l1 & r0coord <= l2;
     r0mid = r0coord(mid_mask);
     r_mid = r_subset(mid_mask);
     % compute graded stress
     taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy(r0mid).^v_a).^((v_nc-1)/v_a)).* ...
            ((r0mid.^4 ./ r_mid.^4) - (r_mid.^2 ./ r0mid.^2));
     % removing the negative infinities
     taurr2(isinf(taurr2)) = 0;
     taurr = taurr1+taurr2;

elseif r_eval > l2
    % near field
    near_mask = r0coord <= l1;
    r0near = r0coord(near_mask);
    r_near = r_subset(near_mask);
    % compute near stress
    taurr1 = (2/(3*Ca))*((r0near.^4 ./ r_near.^4) - (r_near.^2 ./ r0near.^2));
    % zero out the stress above and below the near field
    taurr1(isinf(taurr1)) = 0;

    % compute the graded region coordinate
    mid_mask = r0coord > l1 & r0coord <= l2;
    r0mid = r0coord(mid_mask);
    r_mid = r_subset(mid_mask);
    % compute graded stress
    taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy(r0mid).^v_a).^((v_nc-1)/v_a)).* ...
        ((r0mid.^4 ./ r_mid.^4) - (r_mid.^2 ./ r0mid.^2));
    % removing the negative infinities
    taurr2(isinf(taurr2)) = 0;

    % far field
    far_mask = r0coord > l2;
    r0far = r0coord(far_mask);
    r_far = r_subset(far_mask);
    % compute far stress
    taurr3 = (2/(3*Ca1))*((r0far.^4 ./ r_far.^4) - (r_far.^2 ./ r0far.^2));
    % zero out the stress below the far field
    taurr3(isinf(taurr3)) = 0;
    taurr = taurr1 + taurr2 + taurr3;
end

r1 = (l1^3 + R.^3 - Req^3).^(1/3);
r2 = (l2^3 + R.^3 - Req^3).^(1/3);

end
