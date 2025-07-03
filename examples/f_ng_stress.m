% file f_ng_stress.m
% brief contains function f_ng_stress

% brief This function computes the stress field for a homogeneous, 
% isotropic material as a function of time and space
function [taurr,r1,r2] = f_ng_stress(rcoord,R,Req,Ca,l1,l2)

r1 = (l1^3 + R.^3 - Req^3).^(1/3);
r2 = (l2^3 + R.^3 - Req^3).^(1/3);

% reference coordinate calculation, with slight shift away from bubble wall
% to be consistent with f_graded_stress
r0coord = rcoord.^3 - (R+0.001).^3 + Req.^3;
% identifying where r0coord is inside the bubble
%valid = r0coord > 0 & isreal(r0coord);
%r0coord = (r0coord.*valid).^(1/3);
% removing negative r0coord values
r0coord(r0coord < 0) = NaN;
r0coord = r0coord.^(1/3);

% compute stress
taurr = (2/(3*Ca))*((r0coord.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0coord.^2));
% removing infinities in stress calculation
taurr(isinf(taurr)) = 0;
end