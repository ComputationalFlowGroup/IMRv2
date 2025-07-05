% file f_ng_stress.m
% brief contains function f_ng_stress

% brief This function computes the stress field for a homogeneous, 
% isotropic material as a function of time and space
function [taurr] = f_ng_stress(r_eval,rcoord,R,Req,Ca)

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

% compute stress
taurr = (2/(3*Ca))*((r0coord.^4 ./ r_subset.^4) - (r_subset.^2 ./ r0coord.^2));
% removing infinities in stress calculation
taurr(isinf(taurr)) = 0;
end