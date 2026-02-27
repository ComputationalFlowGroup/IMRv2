% file f_graded_stress.m
% brief contains function f_graded_stress

% brief This function sets computes the graded stress field as a function
% of time and space
function [taurr_near,taurr_mid,taurr_far,mask_near,mask_grad,mask_far,r1,r2] = f_graded_stress(rcoord,R,Req,Ca,Ca1,l1,l2,v_nc,v_a)

r1 = (l1^3 + R.^3 - Req^3).^(1/3);
r2 = (l2^3 + R.^3 - Req^3).^(1/3);

% reference coordinate calculation
r0coord = (rcoord.^3 - R.^3 + Req.^3).^(1/3);
r0shift = (rcoord.^3 + Req^3 - (R+0.001).^3).^(1/3);
% removing the data within the bubble and slightly away from the bubble wall
r0coord(r0coord < r0shift) = NaN;

% masks for each region
mask_near = r0coord <= l1 - 1e-6;
mask_far = r0coord >= l2 + 1e-6;
mask_grad = ~(mask_near | mask_far);

% function for graded region
f_cy = @(r0coord) (l2 - r0coord) ./ (r0coord - l1);
tau = @(r,r0) (2/3)*(r0.^4 ./ r.^4 - r.^2 ./ r0.^2);
gtau = @(r,r0) Ca+(Ca1-Ca)*((1+ f_cy(r0).^v_a).^((v_nc-1)/v_a));

% total stress
%taurr = NaN(size(rcoord));

% near field
r0near = r0coord(mask_near);
rnear = rcoord(mask_near); %????
diff_rnear = abs(rnear - r0near);
fprintf('Min |r - r0|: %g, Max |r - r0|: %g\n', min(diff_rnear(:)), max(diff_rnear(:)))
% taurr(mask_near) = 1/Ca *tau(rnear,r0near);
% taurr_near = taurr(mask_near);
taurr_near = 1/Ca *tau(rnear,r0near);

% far field
r0far = r0coord(mask_far);
rfar = rcoord(mask_far);
diff_rfar = abs(rfar - r0far);
fprintf('Min |r - r0|: %g, Max |r - r0|: %g\n', min(diff_rfar(:)), max(diff_rfar(:)))
% taurr(mask_far) = 1/Ca1 * tau(rfar,r0far);
% taurr_far = taurr(mask_far);
taurr_far = 1/Ca1 * tau(rfar,r0far);

% mid field
r0mid = r0coord(mask_grad);
rmid = rcoord(mask_grad);
diff_rmid = abs(rmid - r0mid);
fprintf('Min |r - r0|: %g, Max |r - r0|: %g\n', min(diff_rmid(:)), max(diff_rmid(:)))
% taurr(mask_grad) = gtau(rmid,r0mid).*tau(rmid,r0mid);
% taurr_mid = taurr(mask_grad);
taurr_mid = gtau(rmid,r0mid).*tau(rmid,r0mid);

% % zero out regions larger than the near field
% r0near(r0near > l1) = 0;
% taurr1 = (2/(3*Ca))*((r0near.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0near.^2));
% % zero out the stress above and below the near field
% taurr1(isinf(taurr1)) = 0;
% % far field
% r0far = r0coord;
% % zero out regions less than the far field
% r0far(r0far < l2) = 0;
% taurr3 = (2/(3*Ca1))*((r0far.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0far.^2));
% % zero out the stress below the far field
% taurr3(isinf(taurr3)) = 0;
% % compute the graded region coordinate
% r0mid = r0coord - r0far - r0near;
% % graded stress
% taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).* ...
%     ((r0mid.^4 ./ rcoord.^4) - (rcoord.^2 ./ r0mid.^2));
% % removing the negative infinities
% taurr2(isinf(taurr2)) = 0;
% % sum of the near, far, and graded stress fields
% %taurr = taurr1 + taurr2 + taurr3;

% figure
% strain = abs(rcoord - r0coord) ./ r0coord;
% imagesc(strain); colorbar; title("Relative radial strain");

end
