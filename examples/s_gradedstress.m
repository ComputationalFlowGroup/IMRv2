

function [taurr,r1,r2] = f_gradedstress(r_coord,R,Req,Ca,Ca1,l1,l2,v_nc,v_a)

   r0coord = (r_coord.^3 + Req^3 - R.^3);
   r0coord(r0coord < 5E-2) = NaN;
   r0_coord = nthroot(r0coord,3);
    
   f_cy = (l2 - r0_coord) ./ (r0_coord - l1);
   r0near = r0_coord;
   r0near(r0near > l1) = 0;
   taurr1 = (2/(3*Ca))*((r0near.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0near.^2));
   taurr1(isinf(taurr1)) = 0;
   r0far = r0_coord;
   r0far(r0far < l2) = 0;
   taurr3 = (2/(3*Ca1))*((r0far.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0far.^2));
   taurr3(isinf(taurr3)) = 0;
   r0mid = r0_coord - r0far - r0near;
   taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).*((r0mid.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0mid.^2));
   taurr2(isinf(taurr2)) = 0;
   taurr = taurr1 + taurr2 + taurr3;

   r1 = (l1^3 + R.^3 - Req^3).^(1/3);
   r2 = (l2^3 + R.^3 - Req^3).^(1/3);
end