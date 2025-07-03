% file f_g_stress.m
% brief contains function f_g_stress

% brief This function sets computes the graded stress field 
% at the desired r location as a function of time
function [taurr,r1,r2] = f_g_stress(r_eval,rcoord,R,Req,Ca,Ca1,l1,l2,v_nc,v_a)

% determine the region r_eval lives in
if r_eval <= l1
    % near field
    region = 'near';
    mask = rcoord <= l1;
elseif r_eval <= l2
    % graded region
    region = 'mid';
    mask = rcoord >= l1 & rcoord <= l2;
else
    % far field
    region = 'far';
    mask = rcoord >= l2 & rcoord <= reval;
end

% truncate rcoord to appropriate region
r_subset = rcoord(mask);
% compute r0coord for the subdomain
r0coord = (r_subset.^3 - R.^3 + Req.^3).^(1/3);
% slight offset near bubble wall
r0shift = (rsubset.^3 + Req^3 - (R+0.001).^3).^(1/3);
% removing the data within the bubble and slightly away from the bubble wall
r0coord(r0coord < r0shift) = NaN;

% function for graded region
f_cy = (l2 - r0coord) ./ (r0coord - l1);

% compute stress according to the region
switch region
    case 'near'
        % near field
        r0near = r0coord;
        % zero out regions larger than the near field
        r0near(r0near > l1) = 0;
        taurr1 = (2/(3*Ca))*((r0near.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0near.^2));
        taurr2 = 0; taurr3 = 0;
        % zero out the stress above and below the near field
        taurr1(isinf(taurr1)) = 0;
    case 'mid'
        % compute the graded region coordinate
        r0mid = r0coord - r0far - r0near;
        % graded stress
        taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).* ...
            ((r0mid.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0mid.^2));
        % removing the negative infinities
        taurr2(isinf(taurr2)) = 0;
    case 'far'
        % far field
        r0far = r0coord;
        % zero out regions less than the far field
        r0far(r0far < l2) = 0;
        taurr3 = (2/(3*Ca1))*((r0far.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0far.^2));
        % zero out the stress below the far field
        taurr3(isinf(taurr3)) = 0;

end

% sum of the near, far, and graded stress fields
taurr = taurr1 + taurr2 + taurr3;

r1 = (l1^3 + R.^3 - Req^3).^(1/3);
r2 = (l2^3 + R.^3 - Req^3).^(1/3);

end
