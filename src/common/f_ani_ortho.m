function [chiS, M1, M2, M3, M4, M5] = f_ani_ortho(nVec, mVec, Nx)
% compute_angular_projection_matrices_realSH
%
% Computes modal projection coefficients for real spherical harmonics.
%
% For expression:
%
%   Ts1
% + Ts2 cos(2 theta)
% + Ts3 cos(4 theta)
% + T1_j Y_j
% + T2_j cos(2 theta) Y_j
% + T3_j cos(4 theta) Y_j
% + T4_j sin(2 theta) dY_j/dtheta
% + T5_j sin(4 theta) dY_j/dtheta
%
% projected against test modes Y_i by multiplying by
% sin(theta) Y_i and integrating over theta, phi.
%
% Inputs:
%   nVec : vector of n values
%   mVec : vector of m values, same length as nVec
%          m > 0 means cosine real harmonic
%          m < 0 means sine real harmonic
%          m = 0 means axisymmetric harmonic
%
% Optional:
%   Nx   : number of Gauss-Legendre quadrature points in x = cos(theta)
%
% Outputs:
%   chiS : nModes x 3 matrix for the superscript-s terms
%
%          chiS(:,1) = int Y_i dOmega
%          chiS(:,2) = int cos(2theta) Y_i dOmega
%          chiS(:,3) = int cos(4theta) Y_i dOmega
%
%   M1   : coupling matrix for Y_j
%   M2   : coupling matrix for cos(2theta) Y_j
%   M3   : coupling matrix for cos(4theta) Y_j
%   M4   : coupling matrix for sin(2theta) dY_j/dtheta
%   M5   : coupling matrix for sin(4theta) dY_j/dtheta
%
% Therefore, if T1,T2,T3,T4,T5 are nModes x 1 vectors,
%
%   F = chiS(:,1)*Ts1 ...
%     + chiS(:,2)*Ts2 ...
%     + chiS(:,3)*Ts3 ...
%     + M1*T1 ...
%     + M2*T2 ...
%     + M3*T3 ...
%     + M4*T4 ...
%     + M5*T5;

    if nargin < 3 || isempty(Nx)
        Nx = max(300, 8*max(nVec(:)) + 80);
    end

    nVec = nVec(:);
    mVec = mVec(:);

    if length(nVec) ~= length(mVec)
        error('nVec and mVec must have the same length.');
    end

    nModes = length(nVec);

    if any(abs(mVec) > nVec)
        error('Every mode must satisfy |m| <= n.');
    end

    % Gauss-Legendre quadrature on x = cos(theta), x in [-1,1].
    [x, w] = gaussLegendre1D(Nx);

    % theta-dependent functions written in terms of x = cos(theta)
    cos2 = 2*x.^2 - 1;
    cos4 = 8*x.^4 - 8*x.^2 + 1;

    sinTheta = sqrt(max(0, 1 - x.^2));
    sin2 = 2*x.*sinTheta;
    sin4 = 4*x.*sinTheta.*(2*x.^2 - 1);

    % Store theta-dependent pieces of each real spherical harmonic.
    Ytheta  = zeros(Nx, nModes);
    dYtheta = zeros(Nx, nModes);

    for j = 1:nModes

        n = nVec(j);
        m = abs(mVec(j));

        % MATLAB legendre includes the Condon-Shortley phase.
        P_all = legendre(n, x.');
        Pn = P_all(m+1, :).';

        % Normalization
        logFacRatio = gammaln(n - m + 1) - gammaln(n + m + 1);
        Nnm = sqrt((2*n + 1)/(4*pi) * exp(logFacRatio));

        if m == 0
            A = Nnm;
        else
            A = sqrt(2)*Nnm;
        end

        Ytheta(:,j) = A * Pn;

        % d/dtheta P_n^m(cos theta)
        if n == 0
            dPn_dtheta = zeros(size(x));
        else
            if m <= n - 1
                Pprev_all = legendre(n-1, x.');
                Pprev = Pprev_all(m+1, :).';
            else
                Pprev = zeros(size(x));
            end

            dPn_dtheta = (n*x.*Pn - (n + m)*Pprev) ./ sinTheta;
        end

        dYtheta(:,j) = A * dPn_dtheta;
    end

    % Superscript-s terms: vectors
    chiS = zeros(nModes, 3);

    for i = 1:nModes

        phiLin = phiLinearIntegral(mVec(i));

        chiS(i,1) = phiLin * sum(w .* Ytheta(:,i));
        chiS(i,2) = phiLin * sum(w .* cos2 .* Ytheta(:,i));
        chiS(i,3) = phiLin * sum(w .* cos4 .* Ytheta(:,i));
    end

    % Coupling matrices
    M1 = zeros(nModes, nModes);
    M2 = zeros(nModes, nModes);
    M3 = zeros(nModes, nModes);
    M4 = zeros(nModes, nModes);
    M5 = zeros(nModes, nModes);

    for i = 1:nModes
        Yi = Ytheta(:,i);

        for j = 1:nModes

            phiProd = phiProductIntegral(mVec(i), mVec(j));

            if phiProd == 0
                continue
            end

            Yj  = Ytheta(:,j);
            dYj = dYtheta(:,j);

            M1(i,j) = phiProd * sum(w .* Yi .* Yj);
            M2(i,j) = phiProd * sum(w .* cos2 .* Yi .* Yj);
            M3(i,j) = phiProd * sum(w .* cos4 .* Yi .* Yj);
            M4(i,j) = phiProd * sum(w .* sin2 .* Yi .* dYj);
            M5(i,j) = phiProd * sum(w .* sin4 .* Yi .* dYj);
        end
    end

    % ------------------------------------------------------------
    % Special normalization for the monopole test mode only.
    %
    % Since Y_0^0 = 1/sqrt(4*pi),
    %
    %   int 1 * Y_0^0 dOmega = 4*pi/sqrt(4*pi) = sqrt(4*pi).
    %
    % This row-normalization makes the projection of the constant function 1
    % onto the n = 0, m = 0 test equation return 1 instead of sqrt(4*pi).
    %
    % Important: this modifies only rows where the TEST mode is n = 0, m = 0.
    % It does not modify source-mode columns.
    % ------------------------------------------------------------
    monopoleTestRows = (nVec == 0) & (mVec == 0);

    monopoleNorm = 4*pi/sqrt(4*pi);   % equivalently sqrt(4*pi)

    chiS(monopoleTestRows,:) = chiS(monopoleTestRows,:) ./ monopoleNorm;

    M1(monopoleTestRows,:) = M1(monopoleTestRows,:) ./ monopoleNorm;
    M2(monopoleTestRows,:) = M2(monopoleTestRows,:) ./ monopoleNorm;
    M3(monopoleTestRows,:) = M3(monopoleTestRows,:) ./ monopoleNorm;
    M4(monopoleTestRows,:) = M4(monopoleTestRows,:) ./ monopoleNorm;
    M5(monopoleTestRows,:) = M5(monopoleTestRows,:) ./ monopoleNorm;


    % Clean small numerical roundoff
    tol = 1e-12;

    chiS(abs(chiS) < tol) = 0;
    M1(abs(M1) < tol) = 0;
    M2(abs(M2) < tol) = 0;
    M3(abs(M3) < tol) = 0;
    M4(abs(M4) < tol) = 0;
    M5(abs(M5) < tol) = 0;
end


function val = phiLinearIntegral(m)
% Integral over phi of the real-harmonic phi factor.
%
% For m = 0, phi factor is 1.
% For m ~= 0, phi factor is cos(m phi) or sin(|m| phi), whose integral is zero.

    if m == 0
        val = 2*pi;
    else
        val = 0;
    end
end


function val = phiProductIntegral(m1, m2)
% Integral over phi of product of real-harmonic phi factors.
%
% m > 0: cos(m phi)
% m < 0: sin(|m| phi)
% m = 0: 1

    if m1 == 0 && m2 == 0
        val = 2*pi;

    elseif m1 ~= 0 && m2 ~= 0 && m1 == m2
        val = pi;

    else
        val = 0;
    end
end


function [x, w] = gaussLegendre1D(N)
% Nodes and weights for Gauss-Legendre quadrature on [-1,1].

    beta = (1:N-1).' ./ sqrt(4*(1:N-1).'.^2 - 1);
    J = diag(beta,1) + diag(beta,-1);

    [V,D] = eig(J);
    x = diag(D);

    [x, idx] = sort(x);
    V = V(:,idx);

    w = 2*(V(1,:).').^2;
end