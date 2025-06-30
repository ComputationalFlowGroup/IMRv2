function [taurr1,taurr,r1,r2] = f_gradedstress(r_coord,R,Req,R0,Ca,Ca1,l1,l2,v_nc,v_a)
   % aa = r_coord.^3 - R.^3 + Req^3;
   % aa = (1./(1-(aa<0))).*aa;
   % r0_coord = real((aa).^(1/3));
   r0coord = (r_coord.^3 - R.^3 + Req^3);
   %valid = r0coord > 0 & isreal(r0coord); %where r0_coord values are negative (inside bubble)
   %r0_coord = (r0coord.*valid).^(1/3);
   r0coord(r0coord < 0) = NaN;
   r0_coord = nthroot(r0coord,3);
    
   f_cy = (l2 - r0_coord) ./ (r0_coord - l1);
        
   taurr1 = (2/(3*Ca))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr2 = ((1/Ca) + (1/Ca1 - 1/Ca)*(1+f_cy.^v_a).^((v_nc-1)/v_a)).*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));
   taurr3 = (2/(3*Ca1))*((r0_coord.^4 ./ r_coord.^4) - (r_coord.^2 ./ r0_coord.^2));

   taurr = taurr1 + taurr2 + taurr3;

   r1 = ((l1/R0)^3 + R.^3 - Req^3).^(1/3);
   r2 = ((l2/R0)^3 + R.^3 - Req^3).^(1/3);
end