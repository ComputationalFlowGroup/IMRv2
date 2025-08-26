function [P] = associated_legendre(x, l, m)
    if m <= 0 && l <= 0
        P = 1;
    elseif l == 1 && m <= 0
        P = x;
    elseif l == 1 && m == 1
        P = -(1-x.^2).^(1/2);
    elseif m == l
        P = (-1)^l.*double_factorial(2.*l-1).*(1-x.^2).^(l/2);
    elseif m == l-1
        P = x.*(2*l-1).*associated_legendre(x, l-1, l-1);
    else
        P = ((2*l-1).*x.*associated_legendre(x, l-1, m)-(l+m-1).*associated_legendre(x, l-2, m))./(l-m);
    end
end