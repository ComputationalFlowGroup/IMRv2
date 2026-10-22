% file jacobian_pattern.m
% brief contains function imr.numerics.jacobian_pattern

% Copyright (C) The IMR developers (see README.md)
% SPDX-License-Identifier: GPL-3.0-or-later

% brief This function builds a conservative Jacobian sparsity pattern for
% the nodal solver from its state-vector index bookkeeping.
% State layout: [R, Rdot, P, theta(ith), Tm(itm), kv(ikv), Z(iz)]. The
% pattern is a superset of the true Jacobian: the bubble pressure equation
% and the wall boundary condition feed every field row, so all field rows
% couple to (R, Rdot, P) and to the wall stencil columns, while the finite
% difference operators (scheme = 0) are tridiagonal. The Chebyshev
% operators (scheme = 1) couple every node of a field, and the wall rows
% couple every node of both fields, so the pattern is dense.
function Sp = jacobian_pattern(n, ith, itm, ikv, iz, scheme)

    if scheme == 1
        Sp = sparse(ones(n));
        return
    end

    % wall stencil columns: last three bubble points and first three
    % material points
    iwall = [ith(max(1,end-2):end), ikv(max(1,end-2):end), ...
        itm(1:min(3,end))];
    % columns every row may depend on
    icore = [1 2 3 iwall];

    rows = [];
    cols = [];

    % (R, Rdot, P) and stress rows depend on the core and the stress states
    for r = [1 2 3 iz]
        c = [icore iz];
        rows = [rows, r*ones(1,numel(c))]; %#ok<AGROW>
        cols = [cols, c]; %#ok<AGROW>
    end

    % field rows: tridiagonal within each field and between the bubble
    % temperature and vapor concentration, plus the core columns
    [rt,ct] = f_band(ith, ith);
    [rk,ck] = f_band(ikv, ikv);
    [rtk,ctk] = f_band(ith, ikv);
    [rkt,ckt] = f_band(ikv, ith);
    [rm,cm] = f_band(itm, itm);
    rows = [rows, rt, rk, rtk, rkt, rm];
    cols = [cols, ct, ck, ctk, ckt, cm];
    ifield = [ith itm ikv];
    for r = ifield
        rows = [rows, r*ones(1,numel(icore))]; %#ok<AGROW>
        cols = [cols, icore]; %#ok<AGROW>
    end
    % medium heating depends on the stress states through the dissipation
    for r = itm
        rows = [rows, r*ones(1,numel(iz))]; %#ok<AGROW>
        cols = [cols, iz]; %#ok<AGROW>
    end

    Sp = sparse(rows, cols, 1, n, n);
    Sp = spones(Sp);

    % tridiagonal band between two equally-sized index blocks
    function [r,c] = f_band(ir, ic)
        r = [];
        c = [];
        m = min(numel(ir), numel(ic));
        for k = -1:1
            idx = max(1,1-k):min(m,m-k);
            r = [r, ir(idx)]; %#ok<AGROW>
            c = [c, ic(idx+k)]; %#ok<AGROW>
        end
    end

end
