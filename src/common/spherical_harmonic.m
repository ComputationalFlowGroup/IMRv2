function [Y] = spherical_harmonic(x, phi, l, m)
    if m == 0
        prefactor = sqrt((2.*l+1)./(4*pi));
        Y = prefactor*associated_legendre(x, l, m);
   elseif m > 0
       prefactor = sqrt((2.*l+1)./(4*pi).*factorial(l-m)./factorial(l+m));
       Y = (-1).^m.*sqrt(2).*prefactor.*associated_legendre(x, l, m).*cos(m.*phi);
    else
       prefactor = sqrt((2*l+1)/(4*pi)*factorial(l-abs(m))/factorial(l+abs(m)));
       Y = (-1).^m*sqrt(2).*prefactor.*associated_legendre(x, l, abs(m)).*sin(abs(m).*phi);
    end
end