% file f_imr_one_mode_static_diagnostics.m
% brief maps one-mode static radial/modal force diagnostics
function out = f_imr_one_mode_static_diagnostics(R_values, ep_values, varargin)
%F_IMR_ONE_MODE_STATIC_DIAGNOSTICS Build one-mode static force maps.
%
%   out = f_imr_one_mode_static_diagnostics(R_values, ep_values, args...)
%
%   R_values and ep_values are dimensionless vectors. Pass [] to use
%   defaults centered on the prescribed equilibrium. args are the same
%   name/value pairs passed to f_imr_fd, with these extra diagnostic names:
%
%       'plot'             true/false, default true
%       'diagnostic_time'  time used for pressure waveform evaluation
%       'theta_points'     theta samples for surface admissibility
%       'phi_points'       phi samples for non-axisymmetric admissibility
%
%   The solver args must describe exactly one perturbation mode.

    solver_dir = fileparts(mfilename('fullpath'));
    addpath(solver_dir);
    addpath(fullfile(solver_dir,'..','common'));

    [opts, solver_args] = strip_static_args(varargin);

    [~, ~, ~, init_opts, ~, ~, ~, ~, ~, ~, ~, ~, pert_opts] = ...
        evalc('f_call_params(solver_args{:});');

    Req = init_opts(7);
    n = pert_opts.n(:);
    m = pert_opts.m(:);
    epnmeq = pert_opts.epnmeq(:);

    if isempty(m)
        m = zeros(size(n));
    end
    if length(n) ~= 1
        error('f_imr_one_mode_static_diagnostics requires exactly one mode.');
    end

    epeq = epnmeq(1);
    if nargin < 1 || isempty(R_values)
        R_values = linspace(max(0.05*Req, 0.5*Req), 1.8*Req, 121);
    end
    if nargin < 2 || isempty(ep_values)
        ep_span = max([2, 4*abs(epeq), abs(epeq) + 1.5]);
        ep_values = linspace(epeq - ep_span, epeq + ep_span, 161);
    end

    R_values = R_values(:).';
    ep_values = ep_values(:).';

    [RR, EE] = meshgrid(R_values, ep_values);
    states = [RR(:).'; zeros(1,numel(RR)); nan(1,numel(RR)); ...
        EE(:).'; zeros(1,numel(RR))];

    grid_diag = f_imr_rhs_diagnostics(states, solver_args{:}, ...
        'diagnostic_time', opts.diagnostic_time);

    out = struct();
    out.R_values = R_values;
    out.ep_values = ep_values;
    out.Req = Req;
    out.epeq = epeq;
    out.n = n;
    out.m = m;
    out.Rddot = reshape(grid_diag.Rddot, size(RR));
    out.epddot = reshape(grid_diag.epddot(1,:), size(RR));
    out.rad_mod = reshape(grid_diag.rad_mod, size(RR));
    out.ep_mod = reshape(grid_diag.ep_mod(1,:), size(RR));

    fixed_R_states = [Req.*ones(size(ep_values)); zeros(size(ep_values)); ...
        nan(size(ep_values)); ep_values; zeros(size(ep_values))];
    fixed_R_diag = f_imr_rhs_diagnostics(fixed_R_states, solver_args{:}, ...
        'diagnostic_time', opts.diagnostic_time, 'Rddot_override', 0);
    out.fixed_R = struct();
    out.fixed_R.ep_values = ep_values;
    out.fixed_R.F_eps = fixed_R_diag.epddot(1,:);
    out.fixed_R.ep_mod = fixed_R_diag.ep_mod(1,:);
    out.fixed_R.xi = fixed_R_diag.xi(1,:);
    out.fixed_R.eta = fixed_R_diag.eta(1,:);
    out.fixed_R.elastns = fixed_R_diag.elastns(1,:);
    out.fixed_R.sselastns = fixed_R_diag.sselastns(1,:);
    out.fixed_R.viscns = fixed_R_diag.viscns(1,:);
    out.fixed_R.epinertians = fixed_R_diag.epinertians(1,:);
    out.fixed_R.coupled_Rddot = fixed_R_diag.Rddot;
    out.fixed_R.modal_terms = fixed_R_diag.modal_terms;

    fixed_ep_states = [R_values; zeros(size(R_values)); nan(size(R_values)); ...
        epeq.*ones(size(R_values)); zeros(size(R_values))];
    fixed_ep_diag = f_imr_rhs_diagnostics(fixed_ep_states, solver_args{:}, ...
        'diagnostic_time', opts.diagnostic_time);
    out.fixed_ep = struct();
    out.fixed_ep.R_values = R_values;
    out.fixed_ep.F_R = fixed_ep_diag.Rddot;
    out.fixed_ep.rad_mod = fixed_ep_diag.rad_mod;
    out.fixed_ep.P = fixed_ep_diag.P;
    out.fixed_ep.S = fixed_ep_diag.S;

    out.geometry = surface_admissibility(n, m, ep_values, ...
        opts.theta_points, opts.phi_points);

    R_curves = contour_segments(R_values, ep_values, out.Rddot);
    ep_curves = contour_segments(R_values, ep_values, out.epddot);
    out.nullclines = struct();
    out.nullclines.Rddot_zero = R_curves;
    out.nullclines.epddot_zero = ep_curves;
    out.nullclines.intersections = unique_points( ...
        curve_intersections(R_curves, ep_curves), ...
        grid_cluster_tol(R_values, ep_values));
    out.nullclines.count = size(out.nullclines.intersections, 2);

    if opts.plot
        out.figure = plot_static_diagnostics(out);
    end
end

function [opts, solver_args] = strip_static_args(args)
    if mod(length(args),2) == 1
        error('Inputs after ep_values must be name/value pairs.');
    end

    opts = struct();
    opts.plot = true;
    opts.diagnostic_time = 0;
    opts.theta_points = 361;
    opts.phi_points = 361;

    keep = true(size(args));
    for k = 1:2:length(args)
        name = lower(args{k});
        switch name
            case 'plot'
                opts.plot = args{k+1};
                keep(k:k+1) = false;
            case {'diagnostic_time','time','t'}
                opts.diagnostic_time = args{k+1};
                keep(k:k+1) = false;
            case 'theta_points'
                opts.theta_points = args{k+1};
                keep(k:k+1) = false;
            case 'phi_points'
                opts.phi_points = args{k+1};
                keep(k:k+1) = false;
        end
    end
    solver_args = args(keep);
end

function geom = surface_admissibility(n, m, ep_values, ntheta, nphi)
    theta = linspace(0, pi, ntheta);
    if m == 0
        phi = 0;
    else
        phi = linspace(0, 2*pi, nphi);
    end
    [TH, PH] = ndgrid(theta, phi);
    Y = real_spherical_harmonic(n, m, TH, PH);
    Y = Y(:);

    min_shape = zeros(size(ep_values));
    for k = 1:length(ep_values)
        min_shape(k) = min(1 + ep_values(k).*Y);
    end

    geom = struct();
    geom.ep_values = ep_values;
    geom.min_shape = min_shape;
    geom.admissible = min_shape > 0;
    geom.min_over_range = min(min_shape);
end

function Y = real_spherical_harmonic(n, m, theta, phi)
    x = cos(theta(:)).';
    abs_m = abs(m);
    P_all = legendre(n, x);
    Pnm = reshape(P_all(abs_m+1,:), size(theta));
    norm_factor = sqrt((2*n + 1)/(4*pi) * ...
        exp(gammaln(n - abs_m + 1) - gammaln(n + abs_m + 1)));

    if m == 0
        Y = norm_factor.*Pnm;
    elseif m > 0
        Y = sqrt(2).*norm_factor.*Pnm.*cos(abs_m.*phi);
    else
        Y = sqrt(2).*norm_factor.*Pnm.*sin(abs_m.*phi);
    end
end

function curves = contour_segments(x, y, Z)
    C = contourc(x, y, Z, [0 0]);
    curves = {};
    idx = 1;
    while idx < size(C,2)
        npts = C(2,idx);
        if npts > 1
            curves{end+1} = C(:,(idx+1):(idx+npts)); %#ok<AGROW>
        end
        idx = idx + npts + 1;
    end
end

function points = curve_intersections(curves_a, curves_b)
    points = zeros(2,0);
    for ia = 1:length(curves_a)
        A = curves_a{ia};
        for ib = 1:length(curves_b)
            B = curves_b{ib};
            for ka = 1:(size(A,2)-1)
                p = A(:,ka);
                r = A(:,ka+1) - p;
                for kb = 1:(size(B,2)-1)
                    q = B(:,kb);
                    s = B(:,kb+1) - q;
                    den = cross2(r, s);
                    if abs(den) < eps
                        continue
                    end
                    t = cross2(q - p, s)/den;
                    u = cross2(q - p, r)/den;
                    if t >= -eps && t <= 1+eps && u >= -eps && u <= 1+eps
                        points(:,end+1) = p + t*r; %#ok<AGROW>
                    end
                end
            end
        end
    end
end

function val = cross2(a, b)
    val = a(1)*b(2) - a(2)*b(1);
end

function pts = unique_points(points, tol)
    pts = zeros(2,0);
    for k = 1:size(points,2)
        p = points(:,k);
        if isempty(pts)
            pts = p;
            continue
        end
        dist = sqrt(sum((pts - p).^2, 1));
        if all(dist > tol)
            pts(:,end+1) = p; %#ok<AGROW>
        end
    end
end

function tol = grid_cluster_tol(x, y)
    dx = min(diff(sort(unique(x))));
    dy = min(diff(sort(unique(y))));
    if isempty(dx) || isempty(dy) || isnan(dx) || isnan(dy)
        tol = 0;
    else
        tol = 2*sqrt(dx^2 + dy^2);
    end
end

function fig = plot_static_diagnostics(out)
    fig = figure('Name','IMR one-mode static diagnostics');
    tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');

    nexttile
    contour(out.R_values, out.ep_values, out.Rddot, [0 0], ...
        'Color',[0 0.25 0.85], 'LineWidth',1.5);
    hold on
    contour(out.R_values, out.ep_values, out.epddot, [0 0], ...
        'Color',[0.85 0.2 0.05], 'LineWidth',1.5);
    plot(out.Req, out.epeq, 'ko', 'MarkerFaceColor','k');
    if ~isempty(out.nullclines.intersections)
        plot(out.nullclines.intersections(1,:), ...
            out.nullclines.intersections(2,:), 'ks', 'MarkerFaceColor','y');
    end
    xlabel('R');
    ylabel('\epsilon');
    legend('Rddot = 0','epddot = 0','prescribed eq','intersections', ...
        'Location','best');
    title('Static nullclines');
    grid on

    nexttile
    plot(out.fixed_R.ep_values, out.fixed_R.F_eps, 'k-', 'LineWidth',1.5);
    hold on
    xline(out.epeq, 'k--');
    yline(0, 'k:');
    xlabel('\epsilon at R = Req');
    ylabel('epddot');
    title('Fixed-R modal force');
    grid on

    nexttile
    plot(out.fixed_ep.R_values, out.fixed_ep.F_R, 'k-', 'LineWidth',1.5);
    hold on
    xline(out.Req, 'k--');
    yline(0, 'k:');
    xlabel('R at \epsilon = epeq');
    ylabel('Rddot');
    title('Fixed-epsilon radial force');
    grid on

    nexttile
    plot(out.geometry.ep_values, out.geometry.min_shape, 'k-', ...
        'LineWidth',1.5);
    hold on
    yline(0, 'r--');
    xline(out.epeq, 'k--');
    xlabel('\epsilon');
    ylabel('min(1 + \epsilon Y)');
    title('Surface admissibility');
    grid on
end
