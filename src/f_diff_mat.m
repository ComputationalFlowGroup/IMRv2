% file f_diff_mat.m
% brief contains function f_diff_mat

% brief This function returns the node positions and the differentiation
% matrix of the nodal solver f_imr_nodal for either spatial discretization.
%   [y, D] = f_diff_mat(nodes, order, tm_check, scheme)
% order    : 1 for d/dy, 2 for the second-derivative operator
% tm_check : 0 for the interior (bubble) grid on y = r/R in [0, 1], where
%            the order-2 operator is the spherical Laplacian
%            d2/dy2 + (2/y)*d/dy with the symmetry condition at y = 0;
%            1 for the exterior (medium) grid on xi in [1, -1]
% scheme   : 0 finite differences (uniform grid, second-order stencils),
%            1 Chebyshev collocation
% The finite difference operators are second-order: central
% differences in the interior, (1 +- dy/y) for the Laplacian, (-6, 6)/dy^2
% at the center, the one-sided (0.5, -2, 1.5)/dy wall row for the first
% derivative, and zero boundary rows otherwise. The Chebyshev interior
% operators are built on the doubled grid of 2*nodes-1 points on [-1, 1]
% and folded onto the half y >= 0 for even functions; the first-derivative
% row at the center is zero exactly. The Chebyshev exterior operators are
% the standard differentiation matrix on the Chebyshev points of [-1, 1]
% and its square.
function [y, D] = f_diff_mat(nodes, order, tm_check, scheme)
    
    if scheme == 0
        [y, D] = f_finite_diff(nodes, order, tm_check);
    elseif scheme == 1
        [y, D] = f_chebyshev(nodes, order, tm_check);
    else
        error('f_diff_mat: scheme must be 0 (finite difference) or 1 (Chebyshev)');
    end
    
end

% finite difference operators
function [y, dmatrix] = f_finite_diff(nodes, order, Tm_check)
    
    % coordinate creation
    N = nodes-1;
    if Tm_check == 0
        deltaY = 1/N;
    else
        deltaY = -2/N;
    end
    K = (1:N+1)';
    if Tm_check == 0
        y = (K-1)*deltaY;
    else
        y = 1 + (K-1)*deltaY;
    end
    % interior rows 2..N; spdiags stores diagonal d of an n-by-n matrix
    % in column j of the band, so A(i,i+d) = band(i+d)
    interior = (K >= 2) & (K <= N);
    sub = zeros(nodes,1);
    dia = zeros(nodes,1);
    sup = zeros(nodes,1);
    
    if order == 1
        sup(vertcat(false, interior(1:end-1))) = 0.5;
        sub(vertcat(interior(2:end), false)) = -0.5;
        dmatrix = spdiags([sub dia sup], -1:1, nodes, nodes);
        if Tm_check == 0
            dmatrix(end,end-2:end) = [0.5 -2 1.5];
        end
        dmatrix = dmatrix / deltaY;
        
    elseif order == 2
        if Tm_check == 0
            % rows i = 2..N: 1 - deltaY/y(i), -2, 1 + deltaY/y(i)
            yl = (K-1)*deltaY;
            ii = find(interior);
            sup(ii+1) = 1 + deltaY./yl(ii);
            dia(ii) = -2;
            sub(ii-1) = 1 - deltaY./yl(ii);
            % symmetry condition at the bubble center
            dia(1) = -6;
            sup(2) = 6;
        else
            sup(vertcat(false, interior(1:end-1))) = 1;
            dia(interior) = -2;
            sub(vertcat(interior(2:end), false)) = 1;
        end
        dmatrix = spdiags([sub dia sup], -1:1, nodes, nodes);
        dmatrix = dmatrix / (deltaY^2);
    else
        error('f_diff_mat: order must be 1 or 2');
    end
    
end

% Chebyshev collocation operators
function [y, D] = f_chebyshev(nodes, order, Tm_check)
    
    if order ~= 1 && order ~= 2
        error('f_diff_mat: order must be 1 or 2');
    end
    if Tm_check == 0
        % doubled grid of 2*nodes-1 points on [-1, 1]; the retained half is
        % y >= 0 ordered from the center y = 0 to the wall y = 1
        M = 2*(nodes-1);
        [x, Dfull] = f_cheb(M);
        % center index (1-based) and retained rows/columns
        ic = nodes;
        ir = (nodes:-1:1)';
        % mirrored columns x -> -x
        im = M + 2 - ir;
        x(ic) = 0;
        y = x(ir);
        if order == 1
            A = Dfull;
        else
            % spherical Laplacian d2/dx2 + (2/x)*d/dx away from the center
            % and 3*d2/dx2 at the center (even functions)
            D2 = Dfull*Dfull;
            A = D2;
            notc = true(M+1,1);
            notc(ic) = false;
            A(notc,:) = A(notc,:) + (2./x(notc)).*Dfull(notc,:);
            A(ic,:) = 3*D2(ic,:);
        end
        % fold the mirrored columns onto the retained half
        D = A(ir, ir);
        notc_half = (ir ~= ic);
        D(:, notc_half) = D(:, notc_half) + A(ir, im(notc_half));
        if order == 1
            % symmetry: f'(0) = 0 exactly
            D(1,:) = 0;
        end
    else
        % standard Chebyshev points of [-1, 1] ordered from xi = 1 to -1
        [y, D1] = f_cheb(nodes-1);
        if order == 1
            D = D1;
        else
            D = D1*D1;
        end
    end
    
end

% Chebyshev points x_j = cos(pi*j/N), j = 0..N, and differentiation matrix
function [x, D] = f_cheb(N)
    if N == 0
        x = 1;
        D = 0;
        return;
    end
    j = (0:N)';
    x = cos(pi*j/N);
    c = vertcat(2, ones(N-1,1), 2).*(-1).^j;
    X = repmat(x, 1, N+1);
    dX = X - X';
    D = (c*(1./c)')./(dX + eye(N+1));
    D = D - diag(sum(D, 2));
end
