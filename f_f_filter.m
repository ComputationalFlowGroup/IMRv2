function [taurr] = f_f_filter(taurr,r,R,N,M)
%F_MUFILTER Summary of this function goes here
%   Detailed explanation goes here
for k = 1:N
    for j = 1:M
        if r(k,j) < R(k,j)
            taurr(k,j) = NaN;
        end
    end
end

end