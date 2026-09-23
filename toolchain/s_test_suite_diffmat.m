% file s_test_suite_diffmat.m
% brief contains script to test the differentiation operators

% brief This script checks the operators returned by f_diff_mat: the finite
% difference stencils (scheme 0) against their defining values, and the
% Chebyshev operators (scheme 1) for spectral accuracy on smooth test
% functions and for annihilating constants
clc;
clear;

addpath('../src/');

failed = 0;
fprintf('Checking f_diff_mat...\n');

% finite difference stencils
n = 10;
[y, D1] = f_diff_mat(n, 1, 0, 0);
[~, D2] = f_diff_mat(n, 2, 0, 0);
dy = 1/(n-1);
err = max([abs(y(2) - dy), abs(D1(5,6) - 0.5/dy), abs(D1(5,4) + 0.5/dy), ...
    abs(D1(n,n) - 1.5/dy), abs(D1(n,n-1) + 2/dy), abs(D1(n,n-2) - 0.5/dy), ...
    abs(D2(1,1) + 6/dy^2), abs(D2(1,2) - 6/dy^2), ...
    abs(D2(5,6) - (1 + dy/y(5))/dy^2), abs(D2(5,4) - (1 - dy/y(5))/dy^2), ...
    abs(D2(5,5) + 2/dy^2), full(max(abs(D1(1,:)))), full(max(abs(D2(n,:))))]);
fprintf('Finite difference interior stencils: max error = %.3e\n', err);
failed = failed + (err > 1e-12);
[x, E1] = f_diff_mat(n, 1, 1, 0);
[~, E2] = f_diff_mat(n, 2, 1, 0);
dx = -2/(n-1);
err = max([abs(x(1) - 1), abs(x(n) + 1), abs(E1(5,6) - 0.5/dx), ...
    abs(E2(5,5) + 2/dx^2), abs(E2(5,4) - 1/dx^2), ...
    full(max(abs([E1(1,:) E1(n,:) E2(1,:) E2(n,:)])))]);
fprintf('Finite difference exterior stencils: max error = %.3e\n', err);
failed = failed + (err > 1e-12);

% Chebyshev interior operators on even test functions
for n = [12 24]
    [y, D1] = f_diff_mat(n, 1, 0, 1);
    [~, D2] = f_diff_mat(n, 2, 0, 1);
    f1 = cos(pi*y/2);
    df1 = -(pi/2)*sin(pi*y/2);
    L1 = -(pi/2)^2*cos(pi*y/2) - (pi./y).*sin(pi*y/2);
    L1(1) = -3*(pi/2)^2;
    f2 = y.^2;
    err = max([max(abs(D1*f1 - df1)), max(abs(D2*f1 - L1)), ...
        max(abs(D1*f2 - 2*y)), max(abs(D2*f2 - 6)), ...
        max(abs(D1*ones(n,1))), max(abs(D2*ones(n,1)))]);
    fprintf('Chebyshev interior, %d nodes: max error = %.3e\n', n, err);
    failed = failed + (err > 1e-8);
end

% Chebyshev exterior operators
for n = [12 24]
    [x, D1] = f_diff_mat(n, 1, 1, 1);
    [~, D2] = f_diff_mat(n, 2, 1, 1);
    err = max([max(abs(D1*exp(x) - exp(x))), max(abs(D2*exp(x) - exp(x))), ...
        max(abs(D1*ones(n,1)))]);
    fprintf('Chebyshev exterior, %d nodes: max error = %.3e\n', n, err);
    failed = failed + (err > 1e-8);
end

if failed == 0
    fprintf('✅ All tests PASSED.\n');
    exit(0);
else
    fprintf('❌ %d test(s) FAILED.\n', failed);
    exit(1);
end
