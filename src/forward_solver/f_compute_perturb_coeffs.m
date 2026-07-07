% file f_compute_perturb_coeffs.m
% brief contains function f_compute_perturb_coeffs

% brief This function computes the coefficients for the linearized bubble
% surface perturbation evolution equation
function [epddot, terms] = f_compute_perturb_coeffs(epnm, epnmd, epnmeq, epmod, R, U, Udot, n, Ro, We, Re, Ca, alpha, pertmod)

    lam = R/Ro;
    nmodes = length(n);
    xi = zeros(1,nmodes);
    eta = zeros(1,nmodes);
    elastns = zeros(1,nmodes);
    sselastns = zeros(1,nmodes);
    viscns = zeros(1,nmodes);

    if isequal(Ca, Inf)
        eta = 3.*U./R+2.*(n+1).*(n+2)./(Re.*R^2);
        xi = -(n-1).*Udot./R+4.*(n-1).*(n+1).*1./Re.*U./R^3+(n-1).*(n+ ...
              1).*(n+2)./(2.*We.*R^3);
        xi = xi'; eta = eta';
    elseif pertmod == 0

        % pre allocattion for only inertial case
        eta = 4.*(5+2.*n).*U./((4+n).*R);
        xi = 6.*(n+2)./(n+4).*U^2./R^2+(8+n-n.^2).*Udot./((n+4).*R);


        for i = 1:length(n)
            % --- First-oUer elastic term (coefficient of G) ---
            elast = 1/(lam^11.*(6+n(i)).*(9+n(i)).*(12+n(i)).*Ro^2).*(1+n(i)).*(2.*lam^5.*(1+ ...
                2.*lam^3-lam^6+n(i)).*(648+234.*n(i)+27.*n(i)^2+n(i)^3)+lam^3.*(-1+lam^3).*(13+ ...
                n(i)).*(-n(i).*(14+9.*n(i)+n(i)^2)+2.*lam^6.*(54+15.*n(i)+n(i)^2)+lam^3.*n(i).*(86+ ...
                35.*n(i)+3.*n(i)^2)).*hyp2f1(-(1/3),4+n(i)/3,5+n(i)/3,1-1/lam^3)+(-1+ ...
                lam^3).*(-2.*lam^9.*(54+15.*n(i)+n(i)^2)+lam^6.*n(i).*(54+69.*n(i)+16.*n(i)^2+ ...
                n(i)^3)+n(i).*(140+104.*n(i)+19.*n(i)^2+n(i)^3)-lam^3.*n(i).*(266+199.*n(i)+ ...
                37.*n(i)^2+2.*n(i)^3)).*hyp2f1(2/3,4+n(i)/3,5+n(i)/3,1-1/lam^3));
            
            elastns(i) = 1/(Ca.*lam^14.*Ro^2).*(1+n(i)).*(-2.*epnmeq(i).*lam^5.*(2.*lam^3+n(i))+1/((6+n(i)).*(9+n(i)).* ...
                (12+n(i))).*epnmeq(i).*(1-lam^3).*(lam^3.*(13+n(i)).*(-n(i).*(2+n(i)).*(7+n(i))+2.*lam^6.*(6+n(i)).*(9+n(i))+lam^3.*...
                n(i).*(86+n(i).*(35+3.*n(i)))).*hyp2f1(-(1/3),4+n(i)/3,5+n(i)/3,1-1/lam^3)+(-2.*lam^9.*(6+n(i)).*(9+n(i))+...
                lam^6.*n(i).*(1+n(i)).*(6+n(i)).*(9+n(i))+n(i).*(2+n(i)).*(7+n(i)).*(10+n(i))-lam^3.*n(i).*(2+n(i)).*(7+n(i)).*(19+2.*n(i))).*hyp2f1(2/3,4+n(i)/3,5+n(i)/3,1-1/lam^3)));

            % --- Second-oUer elastic term (coefficient of α.*G) ---
            sselast = (1/(3.*lam^10.*(6+n(i)).*(9+n(i)).*Ro^2)).*(1+n(i)).*(6.*(1-lam^2).*(6+n(i)).*(9+n(i)).*(1+lam^2- ...
                lam^3.*(-4+lam.*(2+lam.*(-4+lam+2.*lam^2+lam^3-2.*lam^5)))-(-2-2.*lam^2+lam^4).*n(i))-(1/(lam.*(12+ ...
                n(i)))).*3.*(1-lam^3).*(lam^3.*((14+n(i)).*(4.*lam^9.*(6+n(i)).*(9+n(i))+2.*lam^6.*(1+n(i)).*(2+n(i)).*(6+ ...
                n(i)).*(9+n(i))+n(i).*(318+n(i).*(175+n(i).*(33+2.*n(i))))-lam^3.*(-324+n(i).*(120+n(i).*(229+n(i).*(59+ ...
                4.*n(i)))))).*hyp2f1(-(2/3),4+n(i)/3,5+n(i)/3,1-1/lam^3)+3.*(13+n(i)).*(n(i).*(2+n(i)).*(7+n(i))- ...
                2.*lam^6.*(6+n(i)).*(9+n(i))-lam^3.*n(i).*(86+n(i).*(35+3.*n(i)))).*hyp2f1(-(1/3),4+n(i)/3,5+n(i)/3,1-1/lam^3))- ...
                2.*(4.*lam^12.*(6+n(i)).*(9+n(i))+2.*lam^9.*(1+n(i)).*(2+n(i)).*(6+n(i)).*(9+n(i))+lam^6.*(1+n(i)).*(6+n(i)).*(9+ ...
                n(i)).*(6+n(i).*(6+n(i)))+n(i).*(2+n(i)).*(11+n(i)).*(39+n(i).*(13+n(i)))-lam^3.*n(i).*(1398+n(i).*(1411+n(i).*(427+ ...
                2.*n(i).*(25+n(i)))))).*hyp2f1(1/3,4+n(i)/3,5+n(i)/3,1-1/lam^3)-3.*(-2.*lam^9.*(6+n(i)).*(9+n(i))+lam^6.*n(i).*(1+ ...
                n(i)).*(6+n(i)).*(9+n(i))+n(i).*(2+n(i)).*(7+n(i)).*(10+n(i))-lam^3.*n(i).*(2+n(i)).*(7+n(i)).*(19+2.*n(i))).*hyp2f1(2/3,4+n(i)/3,5+n(i)/3,1-1/lam^3)));

            sselastns(i) = 1/(Ca^2*lam^13*Ro^2)*epnmeq(i)*(1+n(i))*alpha*(-8*lam^3+12*lam^7-4*lam^9-4*n(i)+6*lam^4*n(i)-2*lam^6*n(i)+1/(lam*(6+n(i))*(9+n(i))*(12+n(i)))*(1-...
                lam^3)*(lam^3*((14+n(i))*(4*lam^9*(6+n(i))*(9+n(i))+2*lam^6*(1+n(i))*(2+n(i))*(6+n(i))*(9+n(i))+n(i)*(318+n(i)*(175+n(i)*(33+2*n(i))))-lam^3*(-324+n(i)*...
                (120+n(i)*(229+n(i)*(59+4*n(i))))))*hyp2f1(-(2/3),4+n(i)/3,5+n(i)/3,1-1/lam^3)+3*(13+n(i))*(n(i)*(2+n(i))*(7+n(i))-2*lam^6*(6+n(i))*(9+n(i))-lam^3*n(i)*...
                (86+n(i)*(35+3*n(i))))*hyp2f1(-(1/3),4+n(i)/3,5+n(i)/3,1-1/lam^3))-2*(4*lam^12*(6+n(i))*(9+n(i))+2*lam^9*(1+n(i))*(2+n(i))*(6+n(i))*(9+n(i))+lam^6*(1+n(i))*...
                (6+n(i))*(9+n(i))*(6+n(i)*(6+n(i)))+n(i)*(2+n(i))*(11+n(i))*(39+n(i)*(13+n(i)))-lam^3*n(i)*(1398+n(i)*(1411+n(i)*(427+2*n(i)*(25+n(i))))))*...
                hyp2f1(1/3,4+n(i)/3,5+n(i)/3,1-1/lam^3)-3*(-2*lam^9*(6+n(i))*(9+n(i))+lam^6*n(i)*(1+n(i))*(6+n(i))*(9+n(i))+n(i)*(2+n(i))*(7+n(i))*(10+n(i))-...
                lam^3*n(i)*(2+n(i))*(7+n(i))*(19+2*n(i)))*hyp2f1(2/3,4+n(i)/3,5+n(i)/3,1-1/lam^3)));
         
            % --- Viscous term (coefficient of μ) for eta ---
            visceta = 2.*(n(i)+1).*(n(i)+2)./R.^2;

            % --- Viscous term (coefficient of μ) for xi ---
            viscxi = (n(i).*(n(i)+1).*(20+7.*n(i))).*U./((6+n(i)).*R^3);

            viscns(i) = (2*epnmeq(i)*n(i)*(1+n(i))*U*Ro^3)/((6+n(i))*R^6);

            % --- Surface tension term (coefficient of γ) ---
            surften = ((n(i)+1).*(n(i)-1).*(n(i)+2)./R^3);

            % Correcting coefficients
            xi(i) = xi(i) + 1/Ca*elast + alpha/Ca*sselast + 1/(2*We)*surften + 1/Re*viscxi;
            eta(i) = eta(i) + 1/Re*visceta;
        end
        xi = xi'; eta = eta';
    elseif pertmod == 1 %Kazuya
        for i = 1:length(n)
            eta(i) = 3.*U./R+2.*(n(i)+1).*(n(i)+2)./(Re.*R^2);
            xi(i) = -(n(i)-1).*Udot./R+4.*(n(i)-1).*(n(i)+1).*1./Re.*U./R^3+(n(i)-1).*(n(i)+ ...
                1).*(n(i)+2)./(2.*We.*R^3)+ (n(i)+1).*(2.*Ro./(Ca.*R^3).*(1+Ro^3./R^3)+n(i).*(n(i)+1)./(Ca.*(R^2+R.*Ro+Ro^2))+ ...
                2.*alpha./Ca.*1./R^2.*(R-Ro)^2./(R.*Ro).*(1+1./lam)^3.*(2-2./lam+3./lam^2- ...
                1./lam^3+1./lam^4)+alpha./Ca.*n(i).*(n(i)+1).*(R-Ro)^2./(5.*R.*Ro.*(R^2+R.*Ro+ ...
                Ro^2)).*(10+6./lam+3./lam^2+1./lam^3));
        end
        xi = xi'; eta = eta';
    end

    epnm = epnm(:);
    epnmd = epnmd(:);
    epnmeq = epnmeq(:);
    if isscalar(epmod)
        epmod = epmod.*ones(nmodes,1);
    else
        epmod = epmod(:);
    end

    elastns = elastns(:);
    sselastns = sselastns(:);
    viscns = viscns(:);
    epinertians = 3.*n(:).*Ro^3.*(-2*U+R*Udot)./((1+n(:)).*(4+n(:)).*R^3).*epnmeq;

    epddot = -eta.*epnmd - xi.*epnm + epmod - elastns - viscns - sselastns + epinertians;

    if nargout > 1
        terms = struct();
        terms.eta = eta;
        terms.xi = xi;
        terms.ep_mod = epmod;
        terms.elastns = elastns;
        terms.sselastns = sselastns;
        terms.viscns = viscns;
        terms.epinertians = epinertians;
        terms.damping_term = -eta.*epnmd;
        terms.linear_stiffness_term = -xi.*epnm;
        terms.anisotropic_term = epmod;
        terms.elastic_equilibrium_term = -elastns;
        terms.viscous_equilibrium_term = -viscns;
        terms.second_order_elastic_equilibrium_term = -sselastns;
        terms.inertial_equilibrium_term = epinertians;
    end
end
