% file f_imr_anisotropic_projection_breakdown.m
% brief decomposes anisotropic radial/modal projection terms
function out = f_imr_anisotropic_projection_breakdown(R, epnm, varargin)
%F_IMR_ANISOTROPIC_PROJECTION_BREAKDOWN Break down anisotropic projections.
%
%   out = f_imr_anisotropic_projection_breakdown(R, epnm, solver_args...)
%
%   R and epnm are nondimensional state values. solver_args are the same
%   name/value pairs passed to f_imr_fd.

    solver_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(solver_dir,'..','common'));

    [~, ~, ~, init_opts, ~, ~, ~, ~, ~, sigma_opts, ~, ~, pert_opts] = ...
        evalc('f_call_params(varargin{:});');

    Req = init_opts(7);
    Ca = sigma_opts(6);
    ani1 = sigma_opts(19);
    ani2 = sigma_opts(20);
    n = pert_opts.n(:);
    m = pert_opts.m(:);
    if isempty(m)
        m = zeros(size(n));
    end
    epnmeq = pert_opts.epnmeq(:);
    epnm = epnm(:);

    [chiS, M1, M2, M3, M4, M5] = f_ani_ortho([0; n], [0; m]);
    [Ts1, Ts2, Ts3, T1, T2, T3, T4, T5] = ...
        f_ani_ortho_time_coeffs(n, m, R/Req, Req, Ca, ani1, ani2, ...
        epnmeq, epnm);

    parts = struct();
    parts.Ts1 = chiS(:,1)*Ts1;
    parts.Ts2 = chiS(:,2)*Ts2;
    parts.Ts3 = chiS(:,3)*Ts3;
    parts.T1 = M1*T1;
    parts.T2 = M2*T2;
    parts.T3 = M3*T3;
    parts.T4 = M4*T4;
    parts.T5 = M5*T5;

    names = fieldnames(parts);
    total = zeros(size(parts.Ts1));
    for k = 1:length(names)
        total = total + parts.(names{k});
    end

    out = struct();
    out.n = n;
    out.m = m;
    out.lambda = R/Req;
    out.epnmeq = epnmeq;
    out.epnm = epnm;
    out.epnmeq_minus_epnm_lam3 = epnmeq - epnm.*out.lambda.^3;
    out.raw = struct('Ts1',Ts1,'Ts2',Ts2,'Ts3',Ts3,'T1',T1,'T2',T2, ...
        'T3',T3,'T4',T4,'T5',T5);
    out.parts = parts;
    out.total = total;
    out.rad_mod = total(1);
    out.ep_mod = total(2:end);
end
